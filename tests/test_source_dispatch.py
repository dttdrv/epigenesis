from __future__ import annotations

from copy import deepcopy
import gzip
import json
import os
from pathlib import Path
import tempfile
import unittest

from brainc._canonical import digest
from brainc.source import (
    FASTA_PROFILE,
    GFF3_PROFILE,
    GENBANK_PROFILE,
    RAW_PROFILE,
    SourceError,
    compile_source,
    load_source_bundle,
    validate_source_descriptor,
)


ROOT = Path(__file__).resolve().parents[1]
FASTA = ROOT / "tests/data/J02482.1.fasta"
GENBANK = ROOT / "tests/data/U49845.1.gb"
GFF3 = ROOT / "tests/data/J02482.1.gff3"


def _reseal(value: dict) -> dict:
    core = {key: item for key, item in value.items() if key != "artifact_sha256"}
    return {**core, "artifact_sha256": digest(core)}


class SourceDispatchTests(unittest.TestCase):
    def test_every_closed_profile_is_deterministic_and_binds_exact_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            raw = root / "raw.dna"
            raw.write_bytes(b"ACGTRYSWKMBDHVN")
            multi = root / "multi.fasta"
            multi.write_bytes(b">a first\nACGT\n>b\nNNry\n")
            compressed = root / "multi.fasta.gz"
            compressed.write_bytes(gzip.compress(multi.read_bytes(), mtime=0))
            cases = (
                (
                    RAW_PROFILE,
                    {"sequence": raw},
                    {"record_id": "raw-1"},
                    {"sequence"},
                    ["raw-1"],
                ),
                (
                    FASTA_PROFILE,
                    {"sequence": multi},
                    {"wrapper": "identity"},
                    {"sequence"},
                    ["a", "b"],
                ),
                (
                    FASTA_PROFILE,
                    {"sequence": compressed},
                    {"wrapper": "gzip"},
                    {"sequence"},
                    ["a", "b"],
                ),
                (
                    GENBANK_PROFILE,
                    {"genbank": GENBANK},
                    {},
                    {"genbank"},
                    ["U49845.1"],
                ),
                (
                    GFF3_PROFILE,
                    {"sequence": FASTA, "annotation": GFF3},
                    {
                        "sequence_profile": FASTA_PROFILE,
                        "sequence_parameters": {"wrapper": "identity"},
                    },
                    {"sequence", "annotation"},
                    ["J02482.1"],
                ),
            )
            for profile, inputs, parameters, roles, record_ids in cases:
                with self.subTest(profile=profile, parameters=parameters):
                    first = compile_source(profile, inputs, parameters=parameters)
                    second = compile_source(profile, inputs, parameters=parameters)
                    self.assertEqual(first.to_dict(), second.to_dict())
                    descriptor = first.to_dict()
                    self.assertEqual(validate_source_descriptor(descriptor), descriptor)
                    self.assertEqual(set(first.artifacts), roles)
                    self.assertEqual(
                        [record["record_id"] for record in descriptor["source_ir"]["records"]],
                        record_ids,
                    )
                    self.assertEqual(
                        set(descriptor["source_ir"]["inputs"]),
                        set(inputs),
                    )
                    self.assertTrue(
                        all(
                            set(reference) == {"sha256", "byte_length"}
                            for reference in descriptor["source_ir"]["inputs"].values()
                        )
                    )
                    self.assertEqual(
                        {
                            role: reference["byte_length"]
                            for role, reference in descriptor["source_ir"]["inputs"].items()
                        },
                        {
                            role: len(Path(path).read_bytes())
                            for role, path in inputs.items()
                        },
                    )
                    self.assertNotIn(str(root), repr(descriptor))

    def test_gff3_supports_the_two_declared_external_sequence_routes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            raw = root / "sequence.dna"
            raw.write_bytes(b"ACGTACGT")
            annotation = root / "features.gff3"
            annotation.write_bytes(
                b"##gff-version 3\nraw\t.\tgene\t1\t8\t.\t+\t.\tID=g\n"
            )
            bundle = compile_source(
                GFF3_PROFILE,
                {"sequence": raw, "annotation": annotation},
                parameters={
                    "sequence_profile": RAW_PROFILE,
                    "sequence_parameters": {"record_id": "raw"},
                },
            )
            self.assertEqual(bundle.to_dict()["source_ir"]["records"][0]["record_id"], "raw")
            self.assertEqual(bundle.artifacts["sequence"]["inputs"]["kind"], "raw-iupac")

            compressed = root / "sequence.fasta.gz"
            compressed.write_bytes(gzip.compress(FASTA.read_bytes(), mtime=0))
            gzip_bundle = compile_source(
                GFF3_PROFILE,
                {"sequence": compressed, "annotation": GFF3},
                parameters={
                    "sequence_profile": FASTA_PROFILE,
                    "sequence_parameters": {"wrapper": "gzip"},
                },
            )
            self.assertEqual(gzip_bundle.artifacts["sequence"]["inputs"]["wrapper"], "gzip")

    def test_profile_and_wrapper_are_explicit_not_sniffed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = root / "input"
            fasta.write_bytes(b">x\nACGT\n")
            compressed = root / "compressed"
            compressed.write_bytes(gzip.compress(fasta.read_bytes(), mtime=0))
            raw = root / "raw"
            raw.write_bytes(b"ACGT")
            rejected = (
                (RAW_PROFILE, {"sequence": fasta}, {"record_id": "x"}),
                (FASTA_PROFILE, {"sequence": raw}, {"wrapper": "identity"}),
                (FASTA_PROFILE, {"sequence": compressed}, {"wrapper": "identity"}),
                (FASTA_PROFILE, {"sequence": fasta}, {"wrapper": "gzip"}),
            )
            for profile, inputs, parameters in rejected:
                with self.subTest(profile=profile, parameters=parameters):
                    with self.assertRaises(SourceError):
                        compile_source(profile, inputs, parameters=parameters)

    def test_unknown_routes_roles_and_parameters_fail_before_frontend_use(self) -> None:
        missing = Path("does-not-exist")
        cases = (
            ("auto", {"sequence": missing}, {}),
            (RAW_PROFILE, {"sequence": missing}, {}),
            (RAW_PROFILE, {"sequence": missing, "extra": missing}, {"record_id": "x"}),
            (FASTA_PROFILE, {"sequence": missing}, {"wrapper": "br"}),
            (GENBANK_PROFILE, {"genbank": missing}, {"version": 273}),
            (
                GFF3_PROFILE,
                {"sequence": missing, "annotation": missing},
                {"sequence_profile": GENBANK_PROFILE, "sequence_parameters": {}},
            ),
        )
        for profile, inputs, parameters in cases:
            with self.subTest(profile=profile, parameters=parameters):
                with self.assertRaises(SourceError) as caught:
                    compile_source(profile, inputs, parameters=parameters)
                self.assertNotIn("does-not-exist", str(caught.exception))

    def test_unhashable_route_enums_raise_source_errors_not_type_errors(self) -> None:
        missing = Path("must-not-be-read")
        attacks = (
            (
                FASTA_PROFILE,
                {"sequence": missing},
                {"wrapper": []},
            ),
            (
                GFF3_PROFILE,
                {"sequence": missing, "annotation": missing},
                {"sequence_profile": [], "sequence_parameters": {}},
            ),
        )
        for profile, inputs, parameters in attacks:
            with self.subTest(profile=profile):
                with self.assertRaisesRegex(SourceError, "SOURCE001:"):
                    compile_source(profile, inputs, parameters=parameters)

    def test_descriptor_is_closed_exact_typed_and_reference_only(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "raw"
            source.write_bytes(b"ACGTACGT")
            descriptor = compile_source(
                RAW_PROFILE,
                {"sequence": source},
                parameters={"record_id": "x"},
            ).to_dict()
            self.assertNotIn("ACGTACGT", repr(descriptor))
            self.assertNotIn("sequence_artifact", repr(descriptor))

            attacks = []
            unknown = deepcopy(descriptor)
            unknown["unknown"] = None
            attacks.append(_reseal(unknown))
            integral_float = deepcopy(descriptor)
            integral_float["version"] = 1.0
            attacks.append(_reseal(integral_float))
            boolean_length = deepcopy(descriptor)
            boolean_length["source_ir"]["inputs"]["sequence"]["byte_length"] = True
            boolean_length["source_ir_sha256"] = digest(boolean_length["source_ir"])
            attacks.append(_reseal(boolean_length))
            for attack in attacks:
                with self.subTest(attack=attack):
                    with self.assertRaises(SourceError):
                        validate_source_descriptor(attack)

    def test_source_directory_round_trip_is_closed_and_link_safe(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "raw"
            source.write_bytes(b"ACGT")
            compiled = compile_source(
                RAW_PROFILE,
                {"sequence": source},
                parameters={"record_id": "x"},
            )
            directory = root / "source-bundle"
            compiled.save(directory)
            loaded = load_source_bundle(directory)
            self.assertEqual(loaded.to_dict(), compiled.to_dict())
            self.assertEqual({item.name for item in directory.iterdir()}, {"source.json", "sequence.json"})

            (directory / "unknown.json").write_text("{}", encoding="utf-8")
            with self.assertRaises(SourceError):
                load_source_bundle(directory)
            (directory / "unknown.json").unlink()

            if os.name == "posix":
                alias = root / "source-link"
                alias.symlink_to(directory, target_is_directory=True)
                with self.assertRaises(SourceError):
                    load_source_bundle(alias)

    def test_resealed_descriptor_input_digests_cannot_detach_from_native_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            raw = root / "raw.dna"
            raw.write_bytes(b"ACGTRYSWKMBDHVN")
            multi = root / "multi.fasta"
            multi.write_bytes(b">a first\nACGT\n>b\nNNry\n")
            compressed = root / "multi.fasta.gz"
            compressed.write_bytes(gzip.compress(multi.read_bytes(), mtime=0))
            cases = (
                (RAW_PROFILE, {"sequence": raw}, {"record_id": "raw-1"}, "sequence"),
                (FASTA_PROFILE, {"sequence": multi}, {"wrapper": "identity"}, "sequence"),
                (FASTA_PROFILE, {"sequence": compressed}, {"wrapper": "gzip"}, "sequence"),
                (GENBANK_PROFILE, {"genbank": GENBANK}, {}, "genbank"),
                (
                    GFF3_PROFILE,
                    {"sequence": FASTA, "annotation": GFF3},
                    {
                        "sequence_profile": FASTA_PROFILE,
                        "sequence_parameters": {"wrapper": "identity"},
                    },
                    "annotation",
                ),
            )
            for index, (profile, inputs, parameters, attacked_role) in enumerate(cases):
                with self.subTest(profile=profile, parameters=parameters):
                    directory = root / f"bundle-{index}"
                    compile_source(profile, inputs, parameters=parameters).save(directory)
                    descriptor_path = directory / "source.json"
                    descriptor = json.loads(descriptor_path.read_text(encoding="utf-8"))
                    descriptor["source_ir"]["inputs"][attacked_role]["sha256"] = "0" * 64
                    descriptor["source_ir_sha256"] = digest(descriptor["source_ir"])
                    descriptor = _reseal(descriptor)
                    descriptor_path.write_text(
                        json.dumps(descriptor, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8",
                    )
                    with self.assertRaisesRegex(
                        SourceError,
                        "input digests do not match native source inputs",
                    ):
                        load_source_bundle(directory)

    def test_genbank_input_length_cannot_detach_from_embedded_native_source(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / "genbank-bundle"
            compile_source(
                GENBANK_PROFILE,
                {"genbank": GENBANK},
                parameters={},
            ).save(directory)
            descriptor_path = directory / "source.json"
            descriptor = json.loads(descriptor_path.read_text(encoding="utf-8"))
            descriptor["source_ir"]["inputs"]["genbank"]["byte_length"] += 1
            descriptor["source_ir_sha256"] = digest(descriptor["source_ir"])
            descriptor = _reseal(descriptor)
            descriptor_path.write_text(
                json.dumps(descriptor, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(SourceError, "input byte length"):
                load_source_bundle(directory)


if __name__ == "__main__":
    unittest.main()

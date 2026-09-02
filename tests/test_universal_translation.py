from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import tempfile
import unittest

from brainc.external_profile import (
    seal_profile_manifest,
    seal_source_descriptor,
    seal_validation_report,
)
from brainc.source import (
    EXTERNAL_PROFILE,
    FASTA_PROFILE,
    GFF3_PROFILE,
    GENBANK_PROFILE,
    PROFILES,
    RAW_PROFILE,
    REFERENCE_FASTA_PROFILE,
    compile_external_source_paths,
    compile_source,
)
from tests.test_external_profile import _profile_ir, _refget, _sha
from tests.test_insdc import _record
from tests.test_validator_development import _Fixture


RECORD_ID = "SAMP0001.1"
SEQUENCE = b"ACGTRYSWKMBDHVN"


class UniversalTranslationTests(unittest.TestCase):
    def test_all_source_routes_share_record_identity_and_reach_development_ir(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sources = self._sources(root)
            records = [source.to_dict()["source_ir"]["records"] for source, _, _ in sources]

            self.assertEqual({source.profile for source, _, _ in sources}, PROFILES)
            self.assertTrue(all(value == records[0] for value in records[1:]))
            self.assertEqual(
                records[0],
                [
                    {
                        "record_id": RECORD_ID,
                        "bases": len(SEQUENCE),
                        "sequence_sha256": _sha(SEQUENCE),
                        "refget_id": _refget(SEQUENCE),
                    }
                ],
            )
            self.assertEqual(
                len({source.to_dict()["artifact_sha256"] for source, _, _ in sources}),
                len(sources),
            )

            for index, (source, inputs, native_inputs) in enumerate(sources):
                with self.subTest(profile=source.profile):
                    fixture = _Fixture(
                        root / f"development-{index}",
                        inputs=inputs,
                        source=source,
                        native_inputs=native_inputs,
                    )
                    report = fixture.validate_paths()
                    self.assertTrue(report["valid"])
                    self.assertEqual(
                        fixture.compiled.artifacts["development_module"]["format"],
                        "brainc.development-module",
                    )
                    self.assertEqual(
                        fixture.compiled.artifacts["development_module"]["version"],
                        1,
                    )

    def _sources(self, root: Path) -> list[tuple[object, dict[str, Path], dict[str, Path]]]:
        raw = root / "sample.dna"
        raw.write_bytes(SEQUENCE.lower())
        fasta = root / "sample.fasta"
        fasta.write_bytes(b">" + RECORD_ID.encode("ascii") + b"\n" + SEQUENCE.lower() + b"\n")
        genbank = root / "sample.gb"
        genbank.write_bytes(_record("SAMP0001", 1, SEQUENCE.decode("ascii")))
        annotation = root / "sample.gff3"
        annotation.write_bytes(
            b"##gff-version 3\n"
            + RECORD_ID.encode("ascii")
            + b"\t.\tregion\t1\t15\t.\t+\t.\tID=whole\n"
        )

        external, external_inputs, external_native = self._external_source(root)
        return [
            (
                compile_source(
                    RAW_PROFILE,
                    {"sequence": raw},
                    parameters={"record_id": RECORD_ID},
                ),
                {"sequence": raw},
                {},
            ),
            (
                compile_source(
                    FASTA_PROFILE,
                    {"sequence": fasta},
                    parameters={"wrapper": "identity"},
                ),
                {"sequence": fasta},
                {},
            ),
            (
                compile_source(
                    REFERENCE_FASTA_PROFILE,
                    {"sequence": fasta},
                    parameters={"wrapper": "identity"},
                ),
                {"sequence": fasta},
                {},
            ),
            (
                compile_source(GENBANK_PROFILE, {"genbank": genbank}, parameters={}),
                {"genbank": genbank},
                {},
            ),
            (
                compile_source(
                    GFF3_PROFILE,
                    {"sequence": fasta, "annotation": annotation},
                    parameters={
                        "sequence_profile": FASTA_PROFILE,
                        "sequence_parameters": {"wrapper": "identity"},
                    },
                ),
                {"sequence": fasta, "annotation": annotation},
                {},
            ),
            (external, external_inputs, external_native),
        ]

    def _external_source(
        self, root: Path
    ) -> tuple[object, dict[str, Path], dict[str, Path]]:
        original = root / "sample.external-dna"
        original.write_bytes(RECORD_ID.encode("ascii") + b"\t" + SEQUENCE + b"\n")
        native_ir_sha256 = _sha(b"org.example.dna-record-index-ir/v1")
        native_bytes = (
            b'{"format":"org.example.dna-record-index","version":1,'
            b'"index_ir_sha256":"'
            + native_ir_sha256.encode("ascii")
            + b'"}'
        )
        native = root / "sample.external-index.json"
        native.write_bytes(native_bytes)

        profile_ir = deepcopy(_profile_ir())
        profile_ir["profile"] = {"id": "org.example.dna-lines", "version": 1}
        profile_ir["grammar"] = {
            "id": "org.example.dna-lines",
            "version": "1",
            "authority": {
                "uri": "https://example.org/dna-lines/v1",
                "sha256": _sha(b"org.example.dna-lines grammar v1"),
            },
        }
        profile_ir["inputs"][0].update(
            {"role": "source", "media_type": "text/plain"}
        )
        profile_ir["native_artifacts"][0].update(
            {
                "role": "record-index",
                "format": "org.example.dna-record-index",
                "schema_sha256": _sha(b"org.example.dna-record-index.schema/v1"),
            }
        )
        profile_ir["validator"]["command"].update(
            {
                "id": "org.example.dna-lines-validator",
                "distribution": "org.example.dna-lines-validator-dist",
            }
        )
        manifest = seal_profile_manifest(profile_ir)
        input_references = [
            {
                "role": "source",
                "sha256": _sha(original.read_bytes()),
                "byte_length": len(original.read_bytes()),
            }
        ]
        native_references = [
            {
                "role": "record-index",
                "format": "org.example.dna-record-index",
                "version": 1,
                "schema_sha256": profile_ir["native_artifacts"][0]["schema_sha256"],
                "sha256": _sha(native_bytes),
                "byte_length": len(native_bytes),
                "ir_sha256": native_ir_sha256,
            }
        ]
        records = [
            {
                "ordinal": 0,
                "input_role": "source",
                "record_id": RECORD_ID,
                "bases": len(SEQUENCE),
                "sequence_sha256": _sha(SEQUENCE),
                "refget_id": _refget(SEQUENCE),
            }
        ]
        descriptor = seal_source_descriptor(
            manifest,
            input_references=input_references,
            native_artifact_references=native_references,
            records=records,
        )
        report = seal_validation_report(
            manifest,
            descriptor,
            replayed_input_references=input_references,
            replayed_native_artifact_references=native_references,
            replayed_records=records,
        )
        source_inputs = {"source": original}
        native_inputs = {"record-index": native}
        source = compile_external_source_paths(
            manifest,
            descriptor,
            report,
            original_paths=source_inputs,
            native_artifact_paths=native_inputs,
        )
        return source, source_inputs, native_inputs


if __name__ == "__main__":
    unittest.main()

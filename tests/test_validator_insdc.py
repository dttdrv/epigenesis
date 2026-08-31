from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

from brainc.insdc import GenBankCompiler
from brainc._canonical import canonical_bytes
import brainc.validator_insdc as validator
from brainc.validator_insdc import (
    INSDCValidationError,
    validate_genbank,
    validate_genbank_paths,
    validate_genbank_report,
)
from tests.test_insdc import _record, _records


ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).parent / "data"


def _reseal(payload: dict[str, object]) -> None:
    payload["artifact_sha256"] = validator._digest(
        {key: value for key, value in payload.items() if key != "artifact_sha256"}
    )


def _rebind_collection_source(payload: dict[str, object], raw: bytes) -> None:
    collection = payload["sequence_collection"]
    source = collection["inputs"]["source"]
    source["byte_length"] = len(raw)
    source["sha256"] = hashlib.sha256(raw).hexdigest()
    source["storage"] = validator._storage(raw, "test source")
    _reseal(collection)
    payload["bio_ir"]["sequence_collection_artifact_sha256"] = collection[
        "artifact_sha256"
    ]
    payload["bio_ir_sha256"] = validator._digest(payload["bio_ir"])
    _reseal(payload)


class IndependentINSDCValidatorTests(unittest.TestCase):
    def test_replay_enforces_aggregate_feature_and_qualifier_limits(self) -> None:
        raw = _record("XV000017", 1, "ACGT") + _record(
            "XW000018", 1, "TGCA"
        )
        with mock.patch.object(validator, "MAX_FEATURES", 1):
            with self.assertRaisesRegex(INSDCValidationError, "aggregate features"):
                validator._replay(raw)
        with mock.patch.object(validator, "MAX_QUALIFIERS", 2):
            with self.assertRaisesRegex(INSDCValidationError, "aggregate qualifiers"):
                validator._replay(raw)

    def test_replay_enforces_the_combined_json_member_budget(self) -> None:
        raw = _record(
            "XX000019",
            1,
            "ACGT" * 4,
            features=[
                "     misc_feature    1..4",
                *['                     /note="x"'] * 5,
                "     misc_feature    5..8",
                *['                     /note="y"'] * 5,
            ],
        )
        replay = validator._replay(raw)
        exact_members = validator._json_member_count(replay)
        with mock.patch.object(validator, "MAX_JSON_MEMBERS", exact_members - 1):
            with self.assertRaisesRegex(INSDCValidationError, "aggregate JSON members"):
                validator._replay(raw)

    def test_canonical_number_rendering_matches_the_producer(self) -> None:
        corpus = [1e19, 1e20, 1e-6, 1e-7, 1.25, -0.0]
        for value in corpus:
            with self.subTest(value=value):
                self.assertEqual(
                    validator._canonical_bytes({"value": value}),
                    canonical_bytes({"value": value}),
                )

    def test_official_ncbi_sample_is_independently_replayed(self) -> None:
        raw = (DATA / "U49845.1.gb").read_bytes()
        artifact = GenBankCompiler().compile_bytes(raw).to_dict()

        report = validate_genbank(artifact, genbank_source=raw)

        self.assertTrue(report["valid"])
        self.assertEqual(report["version"], 2)
        self.assertEqual(report["inputs"]["genbank_source_sha256"], hashlib.sha256(raw).hexdigest())
        self.assertEqual(report["inputs"]["genbank_artifact_sha256"], artifact["artifact_sha256"])
        self.assertEqual(
            report["replay"]["sequence_collection_artifact_sha256"],
            artifact["sequence_collection"]["artifact_sha256"],
        )
        self.assertEqual(report["summary"], {
            "records": 1,
            "source_bytes": len(raw),
            "sequence_bases": 5028,
            "sequence_chunks": 1,
            "features": 6,
            "segments": 6,
            "unresolved_references": [],
        })
        validate_genbank_report(report)

    def test_multirecord_locations_qualifiers_and_collection_match(self) -> None:
        raw = _records()
        artifact = GenBankCompiler().compile_bytes(raw).to_dict()

        report = validate_genbank(artifact, genbank_source=raw)

        self.assertEqual(report["summary"]["records"], 2)
        self.assertEqual(report["summary"]["sequence_bases"], 204)
        self.assertEqual(report["summary"]["features"], 6)
        self.assertEqual(report["summary"]["segments"], 8)

    def test_padded_origin_blank_record_separator_and_final_blank_replay(self) -> None:
        first = _record("VB000014", 1, "ACGT" * 20).replace(
            b"ORIGIN\n",
            b"ORIGIN      \n",
            1,
        )
        second = _record("VC000015", 1, "TGCA" * 20)
        raw = first + b"\n" + second + b"\n"
        artifact = GenBankCompiler().compile_bytes(raw).to_dict()

        report = validate_genbank(artifact, genbank_source=raw)

        self.assertEqual(report["summary"]["records"], 2)
        self.assertEqual(report["summary"]["source_bytes"], len(raw))
        self.assertEqual(
            "".join(
                artifact["sequence_collection"]["inputs"]["source"]["storage"][
                    "chunks"
                ]
            ).encode("ascii"),
            raw,
        )
        header = artifact["bio_ir"]["records"][0]["headers"][-1]
        self.assertEqual(header["value"], "")
        self.assertNotIn("text", header["source_map"][0])

    def test_external_remote_within_and_structured_qualifier_replay(self) -> None:
        raw = _record(
            "VH000008",
            8,
            "ACGT" * 20,
            features=[
                "     misc_feature    ZZ999999.1:2..9",
                '                     /note="unresolved remote"',
                "     misc_feature    complement(2.9)",
                "                     /transl_except=(pos: 2..4,",
                "                     aa: Trp)",
            ],
        )
        artifact = GenBankCompiler().compile_bytes(raw).to_dict()

        report = validate_genbank(artifact, genbank_source=raw)

        self.assertEqual(report["summary"]["unresolved_references"], ["ZZ999999.1"])
        self.assertTrue(
            all(
                feature["feature_id"].startswith(validator.FEATURE_ID_PREFIX)
                for feature in artifact["bio_ir"]["records"][0]["features"]
            )
        )
        within = artifact["bio_ir"]["records"][0]["features"][2]["location"]
        self.assertEqual(within["ast"]["location"]["kind"], "within")
        self.assertEqual(within["segments"][0]["kind"], "uncertain-point")

        unresolved_circular = raw.replace(b"ZZ999999.1:2..9", b"ZZ999999.1:10^1", 1)
        with self.assertRaisesRegex(INSDCValidationError, "adjacent or circular"):
            validator._replay(unresolved_circular)

    def test_direct_and_fully_resealed_semantic_tampering_fail(self) -> None:
        raw = (DATA / "U49845.1.gb").read_bytes()
        artifact = GenBankCompiler().compile_bytes(raw).to_dict()
        direct = copy.deepcopy(artifact)
        direct["bio_ir"]["records"][0]["definition"] = "tampered"
        with self.assertRaisesRegex(INSDCValidationError, "artifact_sha256"):
            validate_genbank(direct, genbank_source=raw)

        resealed = copy.deepcopy(artifact)
        record = resealed["bio_ir"]["records"][0]
        record["definition"] = "tampered"
        record["record_ir_sha256"] = validator._digest(
            {key: value for key, value in record.items() if key != "record_ir_sha256"}
        )
        resealed["bio_ir_sha256"] = validator._digest(resealed["bio_ir"])
        _reseal(resealed)
        with self.assertRaisesRegex(INSDCValidationError, "independent raw GenBank replay"):
            validate_genbank(resealed, genbank_source=raw)

        for replacement in (True, 1.0):
            with self.subTest(type_confusion=type(replacement).__name__):
                confused = copy.deepcopy(artifact)
                confused["bio_ir"]["records"][0]["features"][0]["ordinal"] = replacement
                confused["bio_ir_sha256"] = validator._digest(confused["bio_ir"])
                _reseal(confused)
                with self.assertRaisesRegex(
                    INSDCValidationError,
                    "independent raw GenBank replay|unexpected JSON number",
                ):
                    validate_genbank(confused, genbank_source=raw)

    def test_blank_header_lines_fail_with_validator_errors(self) -> None:
        valid = _record("VQ000019", 1, "ACGT" * 20)
        malformed = (
            valid.replace(
                b"DBLINK      BioProject: PRJNA1\n",
                b"DBLINK      BioProject: PRJNA1\n\n",
                1,
            ),
            valid.replace(
                b"            artificial sequences.\n",
                b"            artificial sequences.\n\n",
                1,
            ),
        )
        for raw in malformed:
            with self.subTest(raw=raw[:160]):
                with self.assertRaises(INSDCValidationError):
                    validator._replay(raw)

    def test_segment_order_and_shape_match_the_genbank_grammar(self) -> None:
        valid = _record("VS000020", 1, "ACGT" * 20)
        segmented = valid.replace(
            b"KEYWORDS    .\n",
            b"KEYWORDS    .\nSEGMENT     1 of 2\n",
            1,
        )
        artifact = GenBankCompiler().compile_bytes(segmented).to_dict()
        self.assertTrue(validate_genbank(artifact, genbank_source=segmented)["valid"])

        malformed = (
            segmented.replace(b"SEGMENT     1 of 2", b"SEGMENT     0 of 2", 1),
            segmented.replace(b"SEGMENT     1 of 2", b"SEGMENT     2 of 1", 1),
            segmented.replace(
                b"SEGMENT     1 of 2",
                b"SEGMENT     1 of\n            2",
                1,
            ),
            valid.replace(
                b"            artificial sequences.\n",
                b"            artificial sequences.\nSEGMENT     1 of 2\n",
                1,
            ),
        )
        for raw in malformed:
            with self.subTest(raw=raw[:180]):
                with self.assertRaises(INSDCValidationError):
                    validator._replay(raw)

    def test_wrong_original_bytes_and_nested_collection_forgery_fail(self) -> None:
        raw = _records()
        artifact = GenBankCompiler().compile_bytes(raw).to_dict()
        with self.assertRaisesRegex(INSDCValidationError, "source storage"):
            validate_genbank(artifact, genbank_source=raw.replace(b"Generated", b"Altered__", 1))

        forged = copy.deepcopy(artifact)
        collection = forged["sequence_collection"]
        collection["refget_seqcol"]["level_2"]["lengths"][0] -= 1
        _reseal(collection)
        forged["bio_ir"]["sequence_collection_artifact_sha256"] = collection[
            "artifact_sha256"
        ]
        forged["bio_ir_sha256"] = validator._digest(forged["bio_ir"])
        _reseal(forged)
        with self.assertRaisesRegex(INSDCValidationError, "refget Sequence Collection"):
            validate_genbank(forged, genbank_source=raw)

        noncanonical = copy.deepcopy(artifact)
        collection = noncanonical["sequence_collection"]
        source_chunks = collection["inputs"]["source"]["storage"]["chunks"]
        source_chunks[:] = [source_chunks[0][:1], source_chunks[0][1:]]
        _reseal(collection)
        noncanonical["bio_ir"]["sequence_collection_artifact_sha256"] = collection[
            "artifact_sha256"
        ]
        noncanonical["bio_ir_sha256"] = validator._digest(noncanonical["bio_ir"])
        _reseal(noncanonical)
        with self.assertRaisesRegex(INSDCValidationError, "short non-final chunk"):
            validate_genbank(noncanonical, genbank_source=raw)

    def test_reference_journal_requires_exact_column_three(self) -> None:
        raw = _record("VI000009", 9, "ACGT" * 20)
        artifact = GenBankCompiler().compile_bytes(raw).to_dict()
        malformed = raw.replace(b"  JOURNAL   Submitted", b"            JOURNAL   Submitted", 1)
        forged = copy.deepcopy(artifact)
        _rebind_collection_source(forged, malformed)
        with self.assertRaisesRegex(INSDCValidationError, "exactly one JOURNAL"):
            validate_genbank(forged, genbank_source=malformed)

    def test_report_tampering_and_closed_shape_fail(self) -> None:
        raw = _record("VJ000010", 10, "ACGT" * 20)
        artifact = GenBankCompiler().compile_bytes(raw).to_dict()
        report = validate_genbank(artifact, genbank_source=raw)
        altered = copy.deepcopy(report)
        altered["summary"]["features"] += 1
        with self.assertRaisesRegex(INSDCValidationError, "report_sha256"):
            validate_genbank_report(altered)
        altered = copy.deepcopy(report)
        altered["unknown"] = True
        with self.assertRaisesRegex(INSDCValidationError, "unknown"):
            validate_genbank_report(altered)

        resealed = copy.deepcopy(report)
        oversized_reference = "Z" + "Z" * 61 + "1.1"
        self.assertEqual(len(oversized_reference), 65)
        resealed["summary"]["unresolved_references"] = [oversized_reference]
        resealed["report_sha256"] = validator._digest(
            {key: value for key, value in resealed.items() if key != "report_sha256"}
        )
        with self.assertRaisesRegex(INSDCValidationError, "at most 64"):
            validate_genbank_report(resealed)

        impossible = copy.deepcopy(report)
        impossible["summary"]["features"] = 0
        impossible["summary"]["segments"] = 0
        impossible["report_sha256"] = validator._digest(
            {key: value for key, value in impossible.items() if key != "report_sha256"}
        )
        with self.assertRaisesRegex(INSDCValidationError, "impossible cross-field"):
            validate_genbank_report(impossible)

        encoded = validator._report_bytes(report)
        with mock.patch.object(validator, "MAX_REPORT_BYTES", len(encoded) - 1):
            with self.assertRaisesRegex(INSDCValidationError, "output limit"):
                validate_genbank_report(report)

    def test_path_api_writes_a_sealed_report_and_rejects_symlinks(self) -> None:
        raw = _record("VK000011", 11, "ACGT" * 20)
        artifact = GenBankCompiler().compile_bytes(raw).to_dict()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_path = root / "source.gb"
            artifact_path = root / "artifact.json"
            report_path = root / "report.json"
            source_path.write_bytes(raw)
            artifact_path.write_text(json.dumps(artifact), encoding="utf-8")
            report = validate_genbank_paths(artifact_path, source_path, report_path=report_path)
            self.assertEqual(json.loads(report_path.read_text(encoding="utf-8")), report)
            validate_genbank_report(report)
            symlink = root / "source-link.gb"
            symlink.symlink_to(source_path)
            with self.assertRaisesRegex(INSDCValidationError, "regular non-linked"):
                validate_genbank_paths(artifact_path, symlink)
            linked_report = root / "report-link.json"
            linked_report.symlink_to(source_path)
            with self.assertRaisesRegex(INSDCValidationError, "regular non-linked"):
                validate_genbank_paths(
                    artifact_path,
                    source_path,
                    report_path=linked_report,
                )
            self.assertEqual(source_path.read_bytes(), raw)

            real_parent = root / "real-parent"
            real_parent.mkdir()
            linked_parent = root / "linked-parent"
            linked_parent.symlink_to(real_parent, target_is_directory=True)
            with self.assertRaisesRegex(INSDCValidationError, "parent.*non-linked"):
                validate_genbank_paths(
                    artifact_path,
                    source_path,
                    report_path=linked_parent / "report.json",
                )
            self.assertFalse((real_parent / "report.json").exists())

            report_path.write_bytes(b"preserve-existing-report")
            with mock.patch.object(
                validator.os,
                "replace",
                side_effect=OSError("simulated publish failure"),
            ):
                with self.assertRaisesRegex(INSDCValidationError, "cannot write"):
                    validate_genbank_paths(
                        artifact_path,
                        source_path,
                        report_path=report_path,
                    )
            self.assertEqual(report_path.read_bytes(), b"preserve-existing-report")
            self.assertEqual(
                list(root.glob(f".{report_path.name}.*")),
                [],
            )

            fallback_report = root / "portable-report.json"
            with mock.patch.object(validator.os, "supports_dir_fd", set()):
                portable = validate_genbank_paths(
                    artifact_path,
                    source_path,
                    report_path=fallback_report,
                )
            self.assertEqual(
                json.loads(fallback_report.read_text(encoding="utf-8")),
                portable,
            )
            fallback_report.write_bytes(b"portable-old-target")
            with (
                mock.patch.object(validator.os, "supports_dir_fd", set()),
                mock.patch.object(
                    validator.os,
                    "replace",
                    side_effect=OSError("portable publish failure"),
                ),
            ):
                with self.assertRaisesRegex(INSDCValidationError, "cannot write"):
                    validate_genbank_paths(
                        artifact_path,
                        source_path,
                        report_path=fallback_report,
                    )
            self.assertEqual(fallback_report.read_bytes(), b"portable-old-target")

    def test_primary_version_identity_exact_64_byte_boundary(self) -> None:
        base = _record("VM000013", 1, "ACGT" * 20)

        def with_accession(length: int) -> bytes:
            accession = "A" + "B" * (length - 2) + "1"
            return base.replace(
                b"ACCESSION   VM000013\n",
                f"ACCESSION   {accession}\n".encode("ascii"),
                1,
            ).replace(
                b"VERSION     VM000013.1\n",
                f"VERSION     {accession}.1\n".encode("ascii"),
                1,
            )

        boundary = with_accession(62)
        artifact = GenBankCompiler().compile_bytes(boundary).to_dict()
        self.assertEqual(len(artifact["bio_ir"]["records"][0]["version"]), 64)
        self.assertTrue(validate_genbank(artifact, genbank_source=boundary)["valid"])

        oversized = with_accession(63)
        forged = copy.deepcopy(artifact)
        _rebind_collection_source(forged, oversized)
        with self.assertRaisesRegex(INSDCValidationError, "VERSION exceeds 64"):
            validate_genbank(forged, genbank_source=oversized)

    def test_secondary_accessions_share_the_64_byte_identity_boundary(self) -> None:
        base = _record("VN000016", 1, "ACGT" * 20)

        def with_secondary(length: int) -> bytes:
            secondary = "A" + "B" * (length - 2) + "1"
            return base.replace(
                b"ACCESSION   VN000016\n",
                f"ACCESSION   VN000016 {secondary}\n".encode("ascii"),
                1,
            )

        boundary = with_secondary(64)
        artifact = GenBankCompiler().compile_bytes(boundary).to_dict()
        self.assertTrue(validate_genbank(artifact, genbank_source=boundary)["valid"])

        oversized = with_secondary(65)
        with self.assertRaisesRegex(Exception, "ACCESSION"):
            GenBankCompiler().compile_bytes(oversized)
        with self.assertRaisesRegex(INSDCValidationError, "ACCESSION"):
            validator._replay(oversized)

    def test_line_and_record_limits_preflight_before_parse_allocation(self) -> None:
        self.assertEqual(validator.MAX_RECORDS, 40_000)
        too_many_lines = b"\n" * (validator.MAX_LINES + 1)
        with mock.patch.object(
            validator,
            "_Line",
            side_effect=AssertionError("line objects must not be allocated"),
        ):
            with self.assertRaisesRegex(INSDCValidationError, "exceeds 1000000 lines"):
                validator._split_lines(too_many_lines)

    def test_import_isolation_and_no_fixture_identity_tokens(self) -> None:
        command = (
            "import sys; "
            f"sys.path.insert(0, {str(ROOT)!r}); "
            "import brainc.validator_insdc; "
            "forbidden={'brainc.insdc','brainc.insdc_location','brainc.sequence',"
            "'brainc.sequence_collection','brainc.sequence_collection_v2',"
            "'brainc._canonical','brainc._io'}; "
            "loaded=sorted(forbidden.intersection(sys.modules)); "
            "assert not loaded, loaded"
        )
        subprocess.run(
            [sys.executable, "-I", "-c", command],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        production = (ROOT / "brainc" / "validator_insdc.py").read_text(encoding="utf-8")
        for token in ("U49845", "AA000001", "BB000002", "J02482", "5386"):
            self.assertNotIn(token, production)
        self.assertFalse(hasattr(validator, "validate_insdc_genbank"))

    def test_isolated_cli_emits_success_and_returns_one_on_tamper(self) -> None:
        raw = _record("VL000012", 12, "ACGT" * 20)
        artifact = GenBankCompiler().compile_bytes(raw).to_dict()
        executable = ROOT / "brainc" / "validator_insdc.py"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.gb"
            compiled = root / "compiled.json"
            source.write_bytes(raw)
            compiled.write_text(json.dumps(artifact), encoding="utf-8")
            success = subprocess.run(
                [sys.executable, "-I", str(executable), str(source), str(compiled)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(success.returncode, 0, success.stderr)
            self.assertTrue(json.loads(success.stdout)["valid"])
            source.write_bytes(raw.replace(b"Generated", b"Altered__", 1))
            failure = subprocess.run(
                [sys.executable, "-I", str(executable), str(source), str(compiled)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(failure.returncode, 1)
            self.assertIn("validation failed", failure.stderr)
            self.assertEqual(failure.stdout, "")

if __name__ == "__main__":
    unittest.main()

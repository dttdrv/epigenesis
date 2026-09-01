from __future__ import annotations

import ast
import base64
from copy import deepcopy
import hashlib
from pathlib import Path
import unittest

from brainc._canonical import digest
import brainc.external_profile as external_profile
from brainc.external_profile import (
    ABI,
    RECORD_CATALOG_SCHEMA,
    VALIDATOR_PROTOCOL,
    ExternalProfileError,
    seal_profile_manifest,
    seal_source_descriptor,
    seal_validation_report,
    validate_external_evidence,
    validate_profile_manifest,
    validate_source_descriptor,
    validate_validation_report,
)


# A small syntax fixture.  The installed real-data gate is specified in the
# external-profile handoff; production code contains no record or accession.
FASTQ = (
    b"@read-1\n"
    b"GATTACA\n"
    b"+\n"
    b"IIIIIII\n"
    b"@read-2\n"
    b"ACGTN\n"
    b"+\n"
    b"IIIII\n"
)
INDEX_IR_SHA256 = hashlib.sha256(b"read-index-ir/v1").hexdigest()
NATIVE_INDEX = (
    b'{"format":"org.example.fastq-read-index","version":1,'
    b'"index_ir_sha256":"' + INDEX_IR_SHA256.encode("ascii") + b'"}'
)


def _sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _refget(raw: bytes) -> str:
    return "SQ." + base64.urlsafe_b64encode(
        hashlib.sha512(raw).digest()[:24]
    ).decode("ascii")


def _input_reference(raw: bytes) -> dict:
    return {"role": "reads", "sha256": _sha(raw), "byte_length": len(raw)}


def _native_reference(raw: bytes = NATIVE_INDEX) -> dict:
    return {
        "role": "read-index",
        "format": "org.example.fastq-read-index",
        "version": 1,
        "schema_sha256": _sha(b"org.example.fastq-read-index.schema/v1"),
        "sha256": _sha(raw),
        "byte_length": len(raw),
        "ir_sha256": INDEX_IR_SHA256,
    }


def _records() -> list[dict]:
    return [
        {
            "ordinal": 0,
            "input_role": "reads",
            "record_id": "read-1",
            "bases": 7,
            "sequence_sha256": _sha(b"GATTACA"),
            "refget_id": _refget(b"GATTACA"),
        },
        {
            "ordinal": 1,
            "input_role": "reads",
            "record_id": "read-2",
            "bases": 5,
            "sequence_sha256": _sha(b"ACGTN"),
            "refget_id": _refget(b"ACGTN"),
        },
    ]


def _profile_ir() -> dict:
    return {
        "abi": ABI,
        "profile": {"id": "org.example.fastq-sanger", "version": 1},
        "grammar": {
            "id": "ncbi.sra.fastq-sanger",
            "version": "2019-09-20/phred33",
            "authority": {
                "uri": "https://www.ncbi.nlm.nih.gov/sra/docs/submitformats",
                "sha256": _sha(b"pinned NCBI FASTQ authority bytes"),
            },
        },
        "inputs": [
            {
                "role": "reads",
                "media_type": "application/vnd.ncbi.fastq",
                "wrapper": "identity",
                "maximum_byte_length": 4096,
            }
        ],
        "native_artifacts": [
            {
                "role": "read-index",
                "media_type": "application/json",
                "format": "org.example.fastq-read-index",
                "version": 1,
                "schema_sha256": _sha(
                    b"org.example.fastq-read-index.schema/v1"
                ),
                "ir_digest_field": "index_ir_sha256",
                "maximum_byte_length": 4096,
            }
        ],
        "catalog": {"schema": RECORD_CATALOG_SCHEMA},
        "limits": {
            "maximum_total_input_bytes": 4096,
            "maximum_total_native_bytes": 4096,
            "maximum_records": 8,
            "maximum_total_bases": 1024,
            "maximum_record_id_bytes": 128,
        },
        "validator": {
            "protocol": VALIDATOR_PROTOCOL,
            "command": {
                "id": "org.example.fastq-validator",
                "distribution": "org.example.fastq-validator-dist",
                "version": "1.0.0",
                "distribution_sha256": _sha(b"validator distribution archive"),
                "executable_sha256": _sha(b"installed validator executable"),
            },
        },
    }


def _closure() -> tuple[dict, dict, dict, list[dict], list[dict], list[dict]]:
    inputs = [_input_reference(FASTQ)]
    native = [_native_reference()]
    records = _records()
    manifest = seal_profile_manifest(_profile_ir())
    descriptor = seal_source_descriptor(
        manifest,
        input_references=inputs,
        native_artifact_references=native,
        records=records,
    )
    report = seal_validation_report(
        manifest,
        descriptor,
        replayed_input_references=inputs,
        replayed_native_artifact_references=native,
        replayed_records=records,
    )
    return manifest, descriptor, report, inputs, native, records


def _reseal_manifest(value: dict) -> dict:
    value["profile_ir_sha256"] = digest(value["profile_ir"])
    value["artifact_sha256"] = digest(
        {key: member for key, member in value.items() if key != "artifact_sha256"}
    )
    return value


def _reseal_descriptor(value: dict) -> dict:
    catalog = value["frontend_ir"]["record_catalog"]
    catalog["catalog_sha256"] = digest(
        {key: member for key, member in catalog.items() if key != "catalog_sha256"}
    )
    value["frontend_ir_sha256"] = digest(value["frontend_ir"])
    value["artifact_sha256"] = digest(
        {key: member for key, member in value.items() if key != "artifact_sha256"}
    )
    return value


def _reseal_report(value: dict) -> dict:
    value["validation_ir_sha256"] = digest(value["validation_ir"])
    value["artifact_sha256"] = digest(
        {key: member for key, member in value.items() if key != "artifact_sha256"}
    )
    return value


class ExternalProfileTests(unittest.TestCase):
    def test_fastq_profile_complete_closure(self) -> None:
        manifest, descriptor, report, _, _, _ = _closure()
        observed = validate_external_evidence(
            manifest,
            descriptor,
            report,
            original_payloads={"reads": FASTQ},
            native_artifact_payloads={"read-index": NATIVE_INDEX},
        )
        self.assertEqual(observed["profile_manifest"], manifest)
        self.assertEqual(observed["source_descriptor"], descriptor)
        self.assertEqual(observed["validation_report"], report)
        self.assertEqual(
            descriptor["frontend_ir"]["profile"],
            {
                "id": "org.example.fastq-sanger",
                "version": 1,
                "manifest_sha256": manifest["artifact_sha256"],
            },
        )
        self.assertEqual(
            descriptor["frontend_ir"]["record_catalog"]["record_count"], 2
        )
        self.assertEqual(
            report["validation_ir"]["validator"],
            manifest["profile_ir"]["validator"]["command"],
        )
        self.assertEqual(_closure()[:3], (manifest, descriptor, report))

    def test_manifest_is_closed_typed_and_bounded(self) -> None:
        manifest = seal_profile_manifest(_profile_ir())

        unknown = deepcopy(manifest)
        unknown["profile_ir"]["plugin"] = "module.name"
        with self.assertRaises(ExternalProfileError):
            validate_profile_manifest(_reseal_manifest(unknown))

        for path, replacement in (
            (("version",), 1.0),
            (("profile_ir", "profile", "version"), True),
            (("profile_ir", "native_artifacts", 0, "version"), 1.0),
            (("profile_ir", "limits", "maximum_records"), True),
        ):
            attacked = deepcopy(manifest)
            cursor = attacked
            for key in path[:-1]:
                cursor = cursor[key]
            cursor[path[-1]] = replacement
            if path[0] == "profile_ir":
                attacked = _reseal_manifest(attacked)
            else:
                attacked["artifact_sha256"] = digest(
                    {
                        key: member
                        for key, member in attacked.items()
                        if key != "artifact_sha256"
                    }
                )
            with self.subTest(path=path):
                with self.assertRaises(ExternalProfileError):
                    validate_profile_manifest(attacked)

        duplicate = _profile_ir()
        duplicate["inputs"].append(deepcopy(duplicate["inputs"][0]))
        with self.assertRaisesRegex(ExternalProfileError, "duplicate role"):
            seal_profile_manifest(duplicate)

        bad_schema = _profile_ir()
        bad_schema["catalog"]["schema"] = "unknown/v1"
        with self.assertRaises(ExternalProfileError):
            seal_profile_manifest(bad_schema)

        binary_native = _profile_ir()
        binary_native["native_artifacts"][0]["media_type"] = "application/octet-stream"
        with self.assertRaisesRegex(ExternalProfileError, "application/json"):
            seal_profile_manifest(binary_native)

        non_string_key = _profile_ir()
        non_string_key[1] = "not JSON"
        with self.assertRaisesRegex(ExternalProfileError, "keys must be strings"):
            seal_profile_manifest(non_string_key)

        excessive_role = _profile_ir()
        excessive_role["inputs"][0]["maximum_byte_length"] = 4097
        with self.assertRaisesRegex(ExternalProfileError, "cumulative input"):
            seal_profile_manifest(excessive_role)

    def test_descriptor_binds_roles_schema_catalog_and_limits(self) -> None:
        manifest, descriptor, _, inputs, native, records = _closure()

        attacks: list[dict] = []
        wrong_role = deepcopy(descriptor)
        wrong_role["frontend_ir"]["inputs"][0]["role"] = "other"
        attacks.append(_reseal_descriptor(wrong_role))

        wrong_schema = deepcopy(descriptor)
        wrong_schema["frontend_ir"]["native_artifacts"][0][
            "schema_sha256"
        ] = _sha(b"other schema")
        attacks.append(_reseal_descriptor(wrong_schema))

        float_version = deepcopy(descriptor)
        float_version["frontend_ir"]["native_artifacts"][0]["version"] = 1.0
        attacks.append(_reseal_descriptor(float_version))

        wrong_ordinal = deepcopy(descriptor)
        wrong_ordinal["frontend_ir"]["record_catalog"]["records"][0][
            "ordinal"
        ] = 1
        attacks.append(_reseal_descriptor(wrong_ordinal))

        wrong_total = deepcopy(descriptor)
        wrong_total["frontend_ir"]["record_catalog"]["total_bases"] += 1
        attacks.append(_reseal_descriptor(wrong_total))

        unbound_role = deepcopy(descriptor)
        unbound_role["frontend_ir"]["record_catalog"]["records"][0][
            "input_role"
        ] = "other"
        attacks.append(_reseal_descriptor(unbound_role))

        unknown = deepcopy(descriptor)
        unknown["frontend_ir"]["unknown"] = None
        attacks.append(_reseal_descriptor(unknown))

        for index, attack in enumerate(attacks):
            with self.subTest(index=index):
                with self.assertRaises(ExternalProfileError):
                    validate_source_descriptor(manifest, attack)

        limited_profile = _profile_ir()
        limited_profile["limits"]["maximum_records"] = 1
        limited = seal_profile_manifest(limited_profile)
        with self.assertRaisesRegex(ExternalProfileError, "record ceiling"):
            seal_source_descriptor(
                limited,
                input_references=inputs,
                native_artifact_references=native,
                records=records,
            )

        short_profile = _profile_ir()
        short_profile["inputs"][0]["maximum_byte_length"] = len(FASTQ) - 1
        short_profile["limits"]["maximum_total_input_bytes"] = len(FASTQ) - 1
        short = seal_profile_manifest(short_profile)
        with self.assertRaisesRegex(ExternalProfileError, "role ceiling"):
            seal_source_descriptor(
                short,
                input_references=inputs,
                native_artifact_references=native,
                records=records,
            )

        malformed_records = deepcopy(records)
        malformed_records[0]["input_role"] = []
        with self.assertRaises(ExternalProfileError):
            seal_source_descriptor(
                manifest,
                input_references=inputs,
                native_artifact_references=native,
                records=malformed_records,
            )

    def test_frontend_output_and_replay_report_have_digest_closure(self) -> None:
        manifest, descriptor, report, inputs, native, records = _closure()

        damaged = deepcopy(descriptor)
        damaged["frontend_ir"]["record_catalog"]["records"][0][
            "sequence_sha256"
        ] = _sha(b"different sequence")
        damaged = _reseal_descriptor(damaged)
        self.assertEqual(validate_source_descriptor(manifest, damaged), damaged)
        with self.assertRaisesRegex(ExternalProfileError, "exactly bind"):
            validate_validation_report(manifest, damaged, report)

        forged_validator = deepcopy(report)
        forged_validator["validation_ir"]["validator"][
            "executable_sha256"
        ] = _sha(b"different executable")
        with self.assertRaisesRegex(ExternalProfileError, "exactly bind"):
            validate_validation_report(
                manifest, descriptor, _reseal_report(forged_validator)
            )

        unknown = deepcopy(report)
        unknown["validation_ir"]["timestamp"] = "2026-01-01T00:00:00Z"
        with self.assertRaises(ExternalProfileError):
            validate_validation_report(manifest, descriptor, _reseal_report(unknown))

        changed_records = deepcopy(records)
        changed_records[0]["sequence_sha256"] = _sha(b"not replayed")
        with self.assertRaisesRegex(ExternalProfileError, "does not reproduce"):
            seal_validation_report(
                manifest,
                descriptor,
                replayed_input_references=inputs,
                replayed_native_artifact_references=native,
                replayed_records=changed_records,
            )

        changed_native = deepcopy(native)
        changed_native[0]["ir_sha256"] = _sha(b"different native IR")
        with self.assertRaisesRegex(ExternalProfileError, "does not reproduce"):
            seal_validation_report(
                manifest,
                descriptor,
                replayed_input_references=inputs,
                replayed_native_artifact_references=changed_native,
                replayed_records=records,
            )

    def test_exact_original_and_native_payload_bytes_are_required(self) -> None:
        manifest, descriptor, report, _, _, _ = _closure()

        original_attacks = (
            {"reads": FASTQ + b"\n"},
            {"other": FASTQ},
            {"reads": bytearray(FASTQ)},
        )
        for attack in original_attacks:
            with self.subTest(attack=list(attack)):
                with self.assertRaises(ExternalProfileError):
                    validate_external_evidence(
                        manifest,
                        descriptor,
                        report,
                        original_payloads=attack,
                        native_artifact_payloads={"read-index": NATIVE_INDEX},
                    )

        for attack in (
            {"read-index": NATIVE_INDEX + b"\n"},
            {"other": NATIVE_INDEX},
        ):
            with self.subTest(attack=list(attack)):
                with self.assertRaises(ExternalProfileError):
                    validate_external_evidence(
                        manifest,
                        descriptor,
                        report,
                        original_payloads={"reads": FASTQ},
                        native_artifact_payloads=attack,
                    )

        wrong_identity = NATIVE_INDEX.replace(b'"version":1', b'"version":1.0')
        wrong_ir = NATIVE_INDEX.replace(
            INDEX_IR_SHA256.encode("ascii"), _sha(b"other IR").encode("ascii")
        )
        duplicate_identity = NATIVE_INDEX.replace(
            b'"format":', b'"format":"duplicate","format":'
        )
        for raw in (wrong_identity, wrong_ir, duplicate_identity):
            native = [_native_reference(raw)]
            descriptor = seal_source_descriptor(
                manifest,
                input_references=[_input_reference(FASTQ)],
                native_artifact_references=native,
                records=_records(),
            )
            report = seal_validation_report(
                manifest,
                descriptor,
                replayed_input_references=[_input_reference(FASTQ)],
                replayed_native_artifact_references=native,
                replayed_records=_records(),
            )
            with self.subTest(native_attack=raw[:48]):
                with self.assertRaises(ExternalProfileError):
                    validate_external_evidence(
                        manifest,
                        descriptor,
                        report,
                        original_payloads={"reads": FASTQ},
                        native_artifact_payloads={"read-index": raw},
                    )

    def test_profile_selection_never_depends_on_payload_content(self) -> None:
        profile = _profile_ir()
        profile["profile"]["id"] = "org.example.explicit-profile"
        manifest = seal_profile_manifest(profile)
        arbitrary = b"this is deliberately not FASTQ"
        inputs = [_input_reference(arbitrary)]
        descriptor = seal_source_descriptor(
            manifest,
            input_references=inputs,
            native_artifact_references=[_native_reference()],
            records=[
                {
                    "ordinal": 0,
                    "input_role": "reads",
                    "record_id": "external-record",
                    "bases": 1,
                    "sequence_sha256": _sha(b"N"),
                    "refget_id": _refget(b"N"),
                }
            ],
        )
        self.assertEqual(
            descriptor["frontend_ir"]["profile"]["id"],
            "org.example.explicit-profile",
        )

    def test_module_is_data_only(self) -> None:
        source = Path(external_profile.__file__).read_text(encoding="utf-8")
        tree = ast.parse(source)
        forbidden_import_roots = {
            "asyncio",
            "importlib",
            "multiprocessing",
            "requests",
            "socket",
            "subprocess",
            "urllib",
        }
        imported_roots: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported_roots.update(alias.name.split(".", 1)[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                imported_roots.add(node.module.split(".", 1)[0])
            elif isinstance(node, ast.Call):
                if isinstance(node.func, ast.Name):
                    self.assertNotIn(
                        node.func.id, {"__import__", "eval", "exec", "compile"}
                    )
                elif isinstance(node.func, ast.Attribute):
                    self.assertNotIn(
                        node.func.attr,
                        {"import_module", "Popen", "run", "system", "urlopen"},
                    )
        self.assertTrue(imported_roots.isdisjoint(forbidden_import_roots))
        self.assertNotIn("SRR", source)
        self.assertNotIn("training", source.lower())
        self.assertNotIn("model", source.lower())


if __name__ == "__main__":
    unittest.main()

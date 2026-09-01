from __future__ import annotations

import base64
from collections.abc import Mapping
from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import brainc.external_profile as external_profile
import brainc.external_source as external_source
from brainc.external_profile import (
    ExternalProfileError,
    seal_profile_manifest,
    seal_source_descriptor,
    seal_validation_report,
)
from brainc.external_source import build_external_source_closure
import brainc.validator_external as validator
from tests.test_external_profile import (
    FASTQ,
    INDEX_IR_SHA256,
    NATIVE_INDEX,
    _closure,
    _input_reference,
    _native_reference,
    _profile_ir,
    _records,
    _sha,
)


ROOT = Path(__file__).resolve().parents[1]


def _seal_closure(manifest: dict, descriptor: dict, report: dict) -> dict:
    closure_ir = {
        "profile_manifest": manifest,
        "source_descriptor": descriptor,
        "validation_report": report,
    }
    core = {
        "format": validator.CLOSURE_FORMAT,
        "version": validator.CLOSURE_VERSION,
        "producer": deepcopy(validator.CLOSURE_PRODUCER),
        "closure_ir": closure_ir,
        "closure_ir_sha256": validator.digest(closure_ir),
    }
    return {**core, "artifact_sha256": validator.digest(core)}


def _fixture() -> tuple[dict, dict, dict, dict]:
    manifest, descriptor, report, _, _, _ = _closure()
    closure = build_external_source_closure(
        manifest,
        descriptor,
        report,
        original_payloads={"reads": FASTQ},
        native_artifact_payloads={"read-index": NATIVE_INDEX},
    )
    return manifest, descriptor, report, closure


def _closure_for_native(raw: bytes, *, ir_sha256: str = INDEX_IR_SHA256) -> dict:
    manifest = seal_profile_manifest(_profile_ir())
    inputs = [_input_reference(FASTQ)]
    native = [_native_reference(raw)]
    native[0]["ir_sha256"] = ir_sha256
    records = _records()
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
    return _seal_closure(manifest, descriptor, report)


def _refget(sequence: bytes) -> str:
    return "SQ." + base64.urlsafe_b64encode(
        hashlib.sha512(sequence).digest()[:24]
    ).decode("ascii")


class IndependentExternalValidatorTests(unittest.TestCase):
    def test_role_mappings_reject_wrong_count_before_iteration(self) -> None:
        class WrongSizedMapping(Mapping[str, object]):
            def __getitem__(self, key: str) -> object:
                raise AssertionError("wrong-size mapping must not be inspected")

            def __iter__(self):
                raise AssertionError("wrong-size mapping must not be iterated")

            def __len__(self) -> int:
                return 2

        value = WrongSizedMapping()
        cases = (
            (
                "producer payloads",
                lambda: external_profile._payload_mapping(
                    value,
                    {"reads"},
                    "original_payloads",
                ),
                external_profile.ExternalProfileError,
            ),
            (
                "source paths",
                lambda: external_source._path_mapping(
                    value,
                    {"reads"},
                    "original_paths",
                ),
                external_source.ExternalSourceError,
            ),
            (
                "validator paths",
                lambda: validator._path_mapping(
                    value,
                    {"reads"},
                    "original_paths",
                ),
                validator.ExternalValidationError,
            ),
        )
        for label, operation, error_type in cases:
            with self.subTest(label=label), self.assertRaisesRegex(
                error_type,
                "role closure",
            ):
                operation()

    def test_external_record_ids_match_the_downstream_source_contract(self) -> None:
        profile = _profile_ir()
        profile["limits"]["maximum_record_id_bytes"] = 257
        with self.subTest(boundary="producer-limit"), self.assertRaisesRegex(
            ExternalProfileError,
            "maximum_record_id_bytes",
        ):
            seal_profile_manifest(profile)

        manifest = seal_profile_manifest(_profile_ir())
        records = _records()
        records[0]["record_id"] = "read one"
        with self.subTest(boundary="producer-whitespace"), self.assertRaisesRegex(
            ExternalProfileError,
            "record_id",
        ):
            seal_source_descriptor(
                manifest,
                input_references=[_input_reference(FASTQ)],
                native_artifact_references=[_native_reference()],
                records=records,
            )

        def coherently_reseal(attacked: dict) -> dict:
            closure_ir = attacked["closure_ir"]
            attacked_manifest = closure_ir["profile_manifest"]
            attacked_manifest["profile_ir_sha256"] = validator.digest(
                attacked_manifest["profile_ir"]
            )
            attacked_manifest["artifact_sha256"] = validator.digest(
                {
                    key: member
                    for key, member in attacked_manifest.items()
                    if key != "artifact_sha256"
                }
            )
            profile_identity = attacked_manifest["profile_ir"]["profile"]
            profile_reference = {
                "id": profile_identity["id"],
                "version": profile_identity["version"],
                "manifest_sha256": attacked_manifest["artifact_sha256"],
            }

            attacked_descriptor = closure_ir["source_descriptor"]
            frontend_ir = attacked_descriptor["frontend_ir"]
            frontend_ir["profile"] = deepcopy(profile_reference)
            catalog = frontend_ir["record_catalog"]
            catalog["catalog_sha256"] = validator.digest(
                {
                    key: member
                    for key, member in catalog.items()
                    if key != "catalog_sha256"
                }
            )
            attacked_descriptor["frontend_ir_sha256"] = validator.digest(frontend_ir)
            attacked_descriptor["artifact_sha256"] = validator.digest(
                {
                    key: member
                    for key, member in attacked_descriptor.items()
                    if key != "artifact_sha256"
                }
            )

            frontend_report = closure_ir["validation_report"]
            validation_ir = frontend_report["validation_ir"]
            validation_ir["profile"] = deepcopy(profile_reference)
            validation_ir["frontend_output_sha256"] = attacked_descriptor[
                "artifact_sha256"
            ]
            validation_ir["replay"] = {
                "inputs": deepcopy(frontend_ir["inputs"]),
                "native_artifacts": deepcopy(frontend_ir["native_artifacts"]),
                "record_catalog_sha256": catalog["catalog_sha256"],
                "frontend_ir_sha256": attacked_descriptor["frontend_ir_sha256"],
            }
            frontend_report["validation_ir_sha256"] = validator.digest(validation_ir)
            frontend_report["artifact_sha256"] = validator.digest(
                {
                    key: member
                    for key, member in frontend_report.items()
                    if key != "artifact_sha256"
                }
            )
            return _seal_closure(
                attacked_manifest,
                attacked_descriptor,
                frontend_report,
            )

        for boundary in ("validator-limit", "validator-whitespace"):
            attacked = deepcopy(_fixture()[3])
            if boundary == "validator-limit":
                attacked["closure_ir"]["profile_manifest"]["profile_ir"]["limits"][
                    "maximum_record_id_bytes"
                ] = 257
            else:
                attacked["closure_ir"]["source_descriptor"]["frontend_ir"][
                    "record_catalog"
                ]["records"][0]["record_id"] = "read one"
            attacked = coherently_reseal(attacked)
            with self.subTest(boundary=boundary), self.assertRaisesRegex(
                validator.ExternalValidationError,
                "record_id|maximum_record_id_bytes",
            ):
                validator.validate_external_source_closure(attacked)

    def test_real_fastq_paths_validate_and_report_is_exact_and_deterministic(self) -> None:
        manifest, descriptor, frontend_report, closure = _fixture()
        self.assertEqual(
            validator.validate_external_source_closure(closure), closure
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "reads.fastq"
            native = root / "read-index.json"
            closure_path = root / "external-closure.json"
            source.write_bytes(FASTQ)
            native.write_bytes(NATIVE_INDEX)
            closure_path.write_text(json.dumps(closure), encoding="utf-8")

            with mock.patch.object(
                validator, "_read_regular", wraps=validator._read_regular
            ) as reads:
                first = validator.validate_external_source(
                    closure,
                    original_paths={"reads": source},
                    native_artifact_paths={"read-index": native},
                )
            original_call = next(
                call
                for call in reads.call_args_list
                if call.kwargs["label"].startswith("original source")
            )
            self.assertIs(original_call.kwargs["capture"], False)

            second = validator.validate_external_paths(
                closure_path,
                original_paths={"reads": source},
                native_artifact_paths={"read-index": native},
            )
            self.assertEqual(first, second)
            self.assertEqual(
                validator.validate_external_report(first, closure), first
            )
            self.assertEqual(first["profile"], descriptor["frontend_ir"]["profile"])
            self.assertEqual(first["evidence"]["closure_sha256"], closure["artifact_sha256"])
            self.assertEqual(
                first["evidence"]["profile_manifest_sha256"],
                manifest["artifact_sha256"],
            )
            self.assertEqual(
                first["evidence"]["frontend_validation_sha256"],
                frontend_report["artifact_sha256"],
            )
            self.assertEqual(first["result"]["records"], 2)
            self.assertEqual(first["result"]["total_bases"], 12)
            self.assertNotIn(FASTQ, json.dumps(first).encode("utf-8"))
            self.assertNotIn(NATIVE_INDEX, json.dumps(first).encode("utf-8"))

    def test_import_isolation_loads_no_producer_compiler_or_shared_helper(self) -> None:
        command = (
            "import sys; "
            f"sys.path.insert(0, {str(ROOT)!r}); "
            "import brainc.validator_external; "
            "forbidden={'brainc.external_profile','brainc.external_source',"
            "'brainc.source','brainc.compiler','brainc.provider',"
            "'brainc._canonical','brainc._io','brainc.v2',"
            "'brainc.validator_reference'}; "
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

    def test_nested_schemas_types_and_all_seals_fail_closed(self) -> None:
        _, _, _, closure = _fixture()

        direct = deepcopy(closure)
        direct["closure_ir"]["source_descriptor"]["frontend_ir"]["record_catalog"][
            "records"
        ][0]["bases"] += 1
        with self.assertRaises(validator.ExternalValidationError):
            validator.validate_external_source_closure(direct)

        bool_version = deepcopy(closure)
        bool_version["version"] = True
        bool_version["artifact_sha256"] = validator.digest(
            {
                key: member
                for key, member in bool_version.items()
                if key != "artifact_sha256"
            }
        )
        with self.assertRaisesRegex(
            validator.ExternalValidationError, "format/version"
        ):
            validator.validate_external_source_closure(bool_version)

        unknown = deepcopy(closure)
        frontend_report = unknown["closure_ir"]["validation_report"]
        frontend_report["validation_ir"]["training_data"] = "forbidden"
        frontend_report["validation_ir_sha256"] = validator.digest(
            frontend_report["validation_ir"]
        )
        frontend_report["artifact_sha256"] = validator.digest(
            {
                key: member
                for key, member in frontend_report.items()
                if key != "artifact_sha256"
            }
        )
        unknown = _seal_closure(
            unknown["closure_ir"]["profile_manifest"],
            unknown["closure_ir"]["source_descriptor"],
            frontend_report,
        )
        with self.assertRaisesRegex(
            validator.ExternalValidationError, "unknown=.*training_data"
        ):
            validator.validate_external_source_closure(unknown)

        wrong_digest = deepcopy(closure)
        wrong_digest["closure_ir_sha256"] = "0" * 64
        wrong_digest["artifact_sha256"] = validator.digest(
            {
                key: member
                for key, member in wrong_digest.items()
                if key != "artifact_sha256"
            }
        )
        with self.assertRaisesRegex(
            validator.ExternalValidationError, "closure_ir_sha256"
        ):
            validator.validate_external_source_closure(wrong_digest)

    def test_same_profile_id_and_version_cannot_substitute_another_manifest(self) -> None:
        manifest, descriptor, report, closure = _fixture()
        other_profile = _profile_ir()
        other_profile["grammar"]["authority"]["sha256"] = _sha(
            b"different pinned grammar authority"
        )
        other_manifest = seal_profile_manifest(other_profile)
        self.assertEqual(
            other_manifest["profile_ir"]["profile"], manifest["profile_ir"]["profile"]
        )
        self.assertNotEqual(other_manifest["artifact_sha256"], manifest["artifact_sha256"])
        attacked = _seal_closure(other_manifest, descriptor, report)
        with self.assertRaisesRegex(
            validator.ExternalValidationError, "exactly bind the profile manifest"
        ):
            validator.validate_external_source_closure(attacked)
        self.assertEqual(
            validator.validate_external_source_closure(closure), closure
        )

    def test_original_roles_lengths_and_hashes_are_bound_to_streamed_files(self) -> None:
        _, _, _, closure = _fixture()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "reads.fastq"
            native = root / "read-index.json"
            source.write_bytes(FASTQ)
            native.write_bytes(NATIVE_INDEX)
            with self.assertRaisesRegex(
                validator.ExternalValidationError, "role closure"
            ):
                validator.validate_external_source(
                    closure,
                    original_paths={"wrong": source},
                    native_artifact_paths={"read-index": native},
                )
            source.write_bytes(FASTQ[:-1] + b"X")
            with self.assertRaisesRegex(
                validator.ExternalValidationError, "exact reference"
            ):
                validator.validate_external_source(
                    closure,
                    original_paths={"reads": source},
                    native_artifact_paths={"read-index": native},
                )
            source.write_bytes(FASTQ + b"X")
            with self.assertRaisesRegex(
                validator.ExternalValidationError, "byte ceiling"
            ):
                validator.validate_external_source(
                    closure,
                    original_paths={"reads": source},
                    native_artifact_paths={"read-index": native},
                )

    def test_native_json_format_version_ir_digest_duplicates_and_bytes_reject(self) -> None:
        attacks = (
            (
                b'{"format":"org.example.wrong","version":1,'
                b'"index_ir_sha256":"' + INDEX_IR_SHA256.encode("ascii") + b'"}',
                "format/version",
            ),
            (
                b'{"format":"org.example.fastq-read-index","version":1.0,'
                b'"index_ir_sha256":"' + INDEX_IR_SHA256.encode("ascii") + b'"}',
                "format/version",
            ),
            (
                b'{"format":"org.example.fastq-read-index","version":1,'
                b'"index_ir_sha256":"' + ("0" * 64).encode("ascii") + b'"}',
                "different IR digest",
            ),
            (
                b'{"format":"org.example.fastq-read-index",'
                b'"format":"org.example.fastq-read-index","version":1,'
                b'"index_ir_sha256":"' + INDEX_IR_SHA256.encode("ascii") + b'"}',
                "duplicate JSON key",
            ),
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "reads.fastq"
            native = root / "index.json"
            source.write_bytes(FASTQ)
            for raw, pattern in attacks:
                with self.subTest(pattern=pattern):
                    native.write_bytes(raw)
                    closure = _closure_for_native(raw)
                    with self.assertRaisesRegex(
                        validator.ExternalValidationError, pattern
                    ):
                        validator.validate_external_source(
                            closure,
                            original_paths={"reads": source},
                            native_artifact_paths={"read-index": native},
                        )

            _, _, _, closure = _fixture()
            native.write_bytes(NATIVE_INDEX + b" ")
            with self.assertRaisesRegex(
                validator.ExternalValidationError, "byte ceiling"
            ):
                validator.validate_external_source(
                    closure,
                    original_paths={"reads": source},
                    native_artifact_paths={"read-index": native},
                )

    @unittest.skipUnless(os.name == "posix", "link attacks require POSIX")
    def test_original_native_and_closure_paths_reject_symbolic_and_hard_links(self) -> None:
        _, _, _, closure = _fixture()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "reads.fastq"
            native = root / "index.json"
            closure_path = root / "closure.json"
            source.write_bytes(FASTQ)
            native.write_bytes(NATIVE_INDEX)
            closure_path.write_text(json.dumps(closure), encoding="utf-8")

            source_link = root / "source-link.fastq"
            source_link.symlink_to(source)
            native_link = root / "native-link.json"
            native_link.symlink_to(native)
            closure_link = root / "closure-link.json"
            closure_link.symlink_to(closure_path)
            cases = (
                (closure_path, source_link, native),
                (closure_path, source, native_link),
                (closure_link, source, native),
            )
            for closure_case, source_case, native_case in cases:
                with self.subTest(case=closure_case), self.assertRaisesRegex(
                    validator.ExternalValidationError, "single-link"
                ):
                    validator.validate_external_paths(
                        closure_case,
                        original_paths={"reads": source_case},
                        native_artifact_paths={"read-index": native_case},
                    )

        for target in ("source", "native", "closure"):
            with self.subTest(hardlink=target), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                source = root / "reads.fastq"
                native = root / "index.json"
                closure_path = root / "closure.json"
                source.write_bytes(FASTQ)
                native.write_bytes(NATIVE_INDEX)
                closure_path.write_text(json.dumps(closure), encoding="utf-8")
                selected = {"source": source, "native": native, "closure": closure_path}[
                    target
                ]
                hardlink = root / f"{target}-hardlink"
                os.link(selected, hardlink)
                with self.assertRaisesRegex(
                    validator.ExternalValidationError, "single-link"
                ):
                    validator.validate_external_paths(
                        hardlink if target == "closure" else closure_path,
                        original_paths={
                            "reads": hardlink if target == "source" else source
                        },
                        native_artifact_paths={
                            "read-index": hardlink if target == "native" else native
                        },
                    )

    def test_path_race_is_detected_after_streaming(self) -> None:
        _, _, _, closure = _fixture()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "reads.fastq"
            native = root / "index.json"
            source.write_bytes(FASTQ)
            native.write_bytes(NATIVE_INDEX)
            original_read = validator.os.read
            raced = False

            def read_and_touch(descriptor: int, size: int) -> bytes:
                nonlocal raced
                raw = original_read(descriptor, size)
                if not raced:
                    raced = True
                    metadata = source.stat()
                    os.utime(
                        source,
                        ns=(metadata.st_atime_ns, metadata.st_mtime_ns + 1_000_000),
                    )
                return raw

            with mock.patch.object(validator.os, "read", side_effect=read_and_touch):
                with self.assertRaisesRegex(
                    validator.ExternalValidationError, "changed while reading"
                ):
                    validator.validate_external_source(
                        closure,
                        original_paths={"reads": source},
                        native_artifact_paths={"read-index": native},
                    )

    def test_duplicate_deep_oversized_json_and_limits_fail_before_admission(self) -> None:
        _, _, _, closure = _fixture()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "reads.fastq"
            native = root / "index.json"
            source.write_bytes(FASTQ)
            native.write_bytes(NATIVE_INDEX)

            duplicate = root / "duplicate.json"
            duplicate.write_text('{"format":"a","format":"b"}', encoding="utf-8")
            with self.assertRaisesRegex(
                validator.ExternalValidationError, "duplicate JSON key"
            ):
                validator.validate_external_paths(
                    duplicate,
                    original_paths={"reads": source},
                    native_artifact_paths={"read-index": native},
                )

            deep = root / "deep.json"
            deep.write_bytes(b'{"x":' + b"[" * 65 + b"0" + b"]" * 65 + b"}")
            with self.assertRaisesRegex(
                validator.ExternalValidationError, "JSON depth ceiling"
            ):
                validator.validate_external_paths(
                    deep,
                    original_paths={"reads": source},
                    native_artifact_paths={"read-index": native},
                )

            closure_path = root / "closure.json"
            closure_path.write_text(json.dumps(closure), encoding="utf-8")
            with self.assertRaisesRegex(
                validator.ExternalValidationError, "byte ceiling"
            ):
                validator.validate_external_paths(
                    closure_path,
                    original_paths={"reads": source},
                    native_artifact_paths={"read-index": native},
                    maximum_closure_bytes=closure_path.stat().st_size - 1,
                )

            oversized = root / "oversized.json"
            with oversized.open("wb") as stream:
                stream.truncate(validator.MAX_CLOSURE_BYTES + 1)
            with self.assertRaisesRegex(
                validator.ExternalValidationError, "byte ceiling"
            ):
                validator.validate_external_paths(
                    oversized,
                    original_paths={"reads": source},
                    native_artifact_paths={"read-index": native},
                )

            with mock.patch.object(validator, "MAX_JSON_MEMBERS", 2):
                with self.assertRaisesRegex(
                    validator.ExternalValidationError, "JSON member ceiling"
                ):
                    validator._load_json(b'{"x":[1,2]}', "evidence")

    def test_memory_errors_are_normalized_at_every_public_boundary(self) -> None:
        _, _, _, closure = _fixture()
        with mock.patch.object(
            validator, "_validate_closure", side_effect=MemoryError("injected")
        ):
            for operation in (
                lambda: validator.validate_external_source_closure(closure),
                lambda: validator.validate_external_source(
                    closure, original_paths={}, native_artifact_paths={}
                ),
            ):
                with self.subTest(operation=operation), self.assertRaisesRegex(
                    validator.ExternalValidationError, "memory ceiling exceeded"
                ):
                    operation()
        with mock.patch.object(
            validator, "_validate_report_schema", side_effect=MemoryError("injected")
        ), self.assertRaisesRegex(
            validator.ExternalValidationError, "memory ceiling exceeded"
        ):
            validator.validate_external_report({})
        with mock.patch.object(
            validator, "_read_regular", side_effect=MemoryError("injected")
        ), self.assertRaisesRegex(
            validator.ExternalValidationError, "memory ceiling exceeded"
        ):
            validator.validate_external_paths(
                "closure.json",
                original_paths={},
                native_artifact_paths={},
            )

    def test_resealed_validation_report_cannot_change_validated_results(self) -> None:
        _, _, _, closure = _fixture()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "reads.fastq"
            native = root / "index.json"
            source.write_bytes(FASTQ)
            native.write_bytes(NATIVE_INDEX)
            report = validator.validate_external_source(
                closure,
                original_paths={"reads": source},
                native_artifact_paths={"read-index": native},
            )
        attacked = deepcopy(report)
        attacked["result"]["records"] += 1
        attacked["report_sha256"] = validator.digest(
            {
                key: member
                for key, member in attacked.items()
                if key != "report_sha256"
            }
        )
        self.assertEqual(validator.validate_external_report(attacked), attacked)
        with self.assertRaisesRegex(
            validator.ExternalValidationError, "does not bind the supplied closure"
        ):
            validator.validate_external_report(attacked, closure)

    def test_randomized_fastq_closures_replay_and_one_byte_attacks_fail(self) -> None:
        generator = random.Random(0xE91E51)
        manifest = seal_profile_manifest(_profile_ir())
        alphabet = b"ACGTN"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source_path = root / "reads.fastq"
            native_path = root / "index.json"
            for case in range(24):
                sequence = bytes(
                    alphabet[generator.randrange(len(alphabet))]
                    for _ in range(generator.randint(1, 192))
                )
                record_id = f"random-read-{case}"
                source = (
                    b"@"
                    + record_id.encode("ascii")
                    + b"\n"
                    + sequence
                    + b"\n+\n"
                    + b"I" * len(sequence)
                    + b"\n"
                )
                ir_sha256 = hashlib.sha256(
                    b"random-index-ir/v1\0" + case.to_bytes(4, "big") + sequence
                ).hexdigest()
                native = json.dumps(
                    {
                        "format": "org.example.fastq-read-index",
                        "version": 1,
                        "index_ir_sha256": ir_sha256,
                        "offsets": [0],
                    },
                    sort_keys=True,
                    separators=(",", ":"),
                ).encode("utf-8")
                input_references = [
                    {
                        "role": "reads",
                        "sha256": hashlib.sha256(source).hexdigest(),
                        "byte_length": len(source),
                    }
                ]
                native_references = [
                    {
                        "role": "read-index",
                        "format": "org.example.fastq-read-index",
                        "version": 1,
                        "schema_sha256": _sha(
                            b"org.example.fastq-read-index.schema/v1"
                        ),
                        "sha256": hashlib.sha256(native).hexdigest(),
                        "byte_length": len(native),
                        "ir_sha256": ir_sha256,
                    }
                ]
                records = [
                    {
                        "ordinal": 0,
                        "input_role": "reads",
                        "record_id": record_id,
                        "bases": len(sequence),
                        "sequence_sha256": hashlib.sha256(sequence).hexdigest(),
                        "refget_id": _refget(sequence),
                    }
                ]
                descriptor = seal_source_descriptor(
                    manifest,
                    input_references=input_references,
                    native_artifact_references=native_references,
                    records=records,
                )
                frontend_report = seal_validation_report(
                    manifest,
                    descriptor,
                    replayed_input_references=input_references,
                    replayed_native_artifact_references=native_references,
                    replayed_records=records,
                )
                closure = _seal_closure(manifest, descriptor, frontend_report)
                source_path.write_bytes(source)
                native_path.write_bytes(native)
                first = validator.validate_external_source(
                    closure,
                    original_paths={"reads": source_path},
                    native_artifact_paths={"read-index": native_path},
                )
                second = validator.validate_external_source(
                    closure,
                    original_paths={"reads": source_path},
                    native_artifact_paths={"read-index": native_path},
                )
                self.assertEqual(first, second)

                attacked = bytearray(source)
                position = generator.randrange(len(attacked))
                attacked[position] ^= 1
                source_path.write_bytes(attacked)
                with self.assertRaises(validator.ExternalValidationError):
                    validator.validate_external_source(
                        closure,
                        original_paths={"reads": source_path},
                        native_artifact_paths={"read-index": native_path},
                    )


if __name__ == "__main__":
    unittest.main()

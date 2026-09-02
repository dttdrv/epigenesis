from __future__ import annotations

from copy import deepcopy
import hashlib
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import brainc.external_runtime as external_runtime
from brainc.external_runtime import ExternalExecutionError, run_command
from brainc.external_profile import seal_profile_manifest
from brainc.external_source import (
    EXECUTABLE_VERSION,
    ExternalSourceError,
    build_executable_external_source_closure,
    source_records,
    validate_external_source_closure,
)
from brainc.validator_external import ExternalValidationError, validate_external_source
from tests.test_external_profile import FASTQ, _executable_profile_ir


ROOT = Path(__file__).resolve().parents[1]
FRONTEND = ROOT / "examples" / "external-fastq" / "frontend.py"
VALIDATOR = ROOT / "examples" / "external-fastq" / "validator.py"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _manifest(frontend: Path = FRONTEND, validator: Path = VALIDATOR) -> dict:
    profile = _executable_profile_ir()
    profile["profile"] = {"id": "org.epigenesis.fastq-sanger", "version": 1}
    profile["native_artifacts"][0].update(
        {
            "role": "fastq",
            "format": "org.epigenesis.fastq-ir",
            "schema_sha256": hashlib.sha256(
                b"org.epigenesis.fastq-ir.schema/v1"
            ).hexdigest(),
            "ir_digest_field": "fastq_ir_sha256",
        }
    )
    profile["frontend"]["command"]["executable_sha256"] = _sha(frontend)
    profile["validator"]["command"]["executable_sha256"] = _sha(validator)
    return seal_profile_manifest(profile, version=2)


class ExternalExecutionTests(unittest.TestCase):
    def _compile(self, raw: bytes, *, manifest: dict | None = None) -> dict:
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "reads.fastq"
            source.write_bytes(raw)
            return build_executable_external_source_closure(
                manifest or _manifest(),
                original_paths={"reads": source},
                frontend_executable=FRONTEND,
                validator_executable=VALIDATOR,
            )

    def test_real_fastq_frontend_and_separate_validator_build_v2_closure(self) -> None:
        closure = self._compile(FASTQ)
        self.assertEqual(closure["version"], EXECUTABLE_VERSION)
        self.assertEqual(validate_external_source_closure(closure), closure)
        self.assertEqual(
            [record["record_id"] for record in source_records(closure)],
            ["read-1", "read-2"],
        )
        self.assertEqual(
            set(closure["closure_ir"]["native_artifacts"]), {"fastq"}
        )
        native = closure["closure_ir"]["native_artifacts"]["fastq"]
        self.assertEqual(native["fastq_ir"]["quality_encoding"], "phred33")
        self.assertEqual(
            native["fastq_ir"]["records"][0],
            {
                "ordinal": 0,
                "header": "read-1 first",
                "record_id": "read-1",
                "sequence": "GATTACA",
                "quality": "IIIIIII",
            },
        )

    def test_frontend_rejects_non_dna_and_malformed_fastq(self) -> None:
        for raw in (
            b"@not-dna\nHELLO\n+\nIIIII\n",
            b"@short\nACGT\n+\nIII\n",
            b"@incomplete\nACGT\n+\n",
            b"@bare-cr\rACGT\r+\rIIII\r",
        ):
            with self.subTest(raw=raw):
                with self.assertRaisesRegex(ExternalSourceError, "external frontend"):
                    self._compile(raw)

    def test_changed_frontend_or_validator_is_rejected_before_execution(self) -> None:
        manifest = _manifest()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "reads.fastq"
            source.write_bytes(FASTQ)
            changed = root / "changed.py"
            changed.write_bytes(FRONTEND.read_bytes() + b"\n# changed\n")
            with self.assertRaisesRegex(ExternalSourceError, "SHA-256"):
                build_executable_external_source_closure(
                    manifest,
                    original_paths={"reads": source},
                    frontend_executable=changed,
                    validator_executable=VALIDATOR,
                )
            with self.assertRaisesRegex(ExternalSourceError, "SHA-256"):
                build_executable_external_source_closure(
                    manifest,
                    original_paths={"reads": source},
                    frontend_executable=FRONTEND,
                    validator_executable=changed,
                )

    def test_separately_executed_validator_rejects_a_lying_frontend(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "reads.fastq"
            source.write_bytes(FASTQ)
            liar = root / "liar.py"
            liar.write_text(
                FRONTEND.read_text(encoding="utf-8").replace(
                    "    fastq_ir = {\"quality_encoding\": \"phred33\", \"records\": native_records}\n",
                    "    records[0][\"bases\"] += 1\n"
                    "    fastq_ir = {\"quality_encoding\": \"phred33\", \"records\": native_records}\n",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ExternalSourceError, "external validator"):
                build_executable_external_source_closure(
                    _manifest(frontend=liar),
                    original_paths={"reads": source},
                    frontend_executable=liar,
                    validator_executable=VALIDATOR,
                )

    def test_validator_result_uses_exact_json_types(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "reads.fastq"
            source.write_bytes(FASTQ)
            confused = root / "validator.py"
            confused.write_text(
                VALIDATOR.read_text(encoding="utf-8").replace(
                    '        "valid": True,\n', '        "valid": 1,\n'
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ExternalSourceError, "exactly replay"):
                build_executable_external_source_closure(
                    _manifest(validator=confused),
                    original_paths={"reads": source},
                    frontend_executable=FRONTEND,
                    validator_executable=confused,
                )

    def test_embedded_native_artifact_is_part_of_the_closure(self) -> None:
        closure = self._compile(FASTQ)
        attacked = deepcopy(closure)
        attacked["closure_ir"]["native_artifacts"]["fastq"]["version"] = 2
        with self.assertRaisesRegex(ExternalSourceError, "native"):
            validate_external_source_closure(attacked)

    def test_independent_validator_reruns_pinned_program(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "reads.fastq"
            source.write_bytes(FASTQ)
            closure = build_executable_external_source_closure(
                _manifest(),
                original_paths={"reads": source},
                frontend_executable=FRONTEND,
                validator_executable=VALIDATOR,
            )
            report = validate_external_source(
                closure,
                original_paths={"reads": source},
                validator_executable=VALIDATOR,
            )
            self.assertTrue(report["valid"])

            source.write_bytes(FASTQ.replace(b"GATTACA", b"GACTACA"))
            with self.assertRaisesRegex(ExternalValidationError, "differs"):
                validate_external_source(
                    closure,
                    original_paths={"reads": source},
                    validator_executable=VALIDATOR,
                )

    def test_version_two_validation_requires_the_pinned_validator(self) -> None:
        closure = self._compile(FASTQ)
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "reads.fastq"
            source.write_bytes(FASTQ)
            with self.assertRaisesRegex(ExternalValidationError, "requires"):
                validate_external_source(
                    closure,
                    original_paths={"reads": source},
                )

    def test_process_runner_fails_closed_on_protocol_and_resource_attacks(self) -> None:
        programs = {
            "duplicate": b'import sys;sys.stdout.write(\'{"x":1,"x":2}\')\n',
            "oversized": b'import sys;sys.stdout.write("x" * 1000)\n',
            "stderr": b'import sys;sys.stderr.write("failure");raise SystemExit(7)\n',
            "timeout": b'import time;time.sleep(1)\n',
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name, raw in programs.items():
                path = root / f"{name}.py"
                path.write_bytes(raw)
                command = {
                    "runtime": "python",
                    "executable_sha256": hashlib.sha256(raw).hexdigest(),
                }
                with self.subTest(name=name):
                    if name == "oversized":
                        context = self.assertRaisesRegex(
                            ExternalExecutionError, "output byte ceiling"
                        )
                        with context:
                            run_command(
                                command,
                                path,
                                {"request": True},
                                label="test command",
                                maximum_stdout_bytes=16,
                            )
                    elif name == "timeout":
                        with mock.patch.object(
                            external_runtime, "MAX_EXECUTION_SECONDS", 0
                        ):
                            with self.assertRaisesRegex(
                                ExternalExecutionError, "exceeded 0 seconds"
                            ):
                                run_command(
                                    command,
                                    path,
                                    {"request": True},
                                    label="test command",
                                )
                    elif name == "stderr":
                        with self.assertRaisesRegex(
                            ExternalExecutionError, "status 7: failure"
                        ):
                            run_command(
                                command,
                                path,
                                {"request": True},
                                label="test command",
                            )
                    else:
                        with self.assertRaisesRegex(
                            ExternalExecutionError, "duplicate JSON key"
                        ):
                            run_command(
                                command,
                                path,
                                {"request": True},
                                label="test command",
                            )

    @unittest.skipUnless(os.name == "posix", "requires a POSIX executable")
    def test_native_runtime_executes_the_pinned_program_without_a_shell(self) -> None:
        raw = b'#!/bin/sh\nprintf \'%s\' \'{"runtime":"native"}\'\n'
        with tempfile.TemporaryDirectory() as temporary:
            executable = Path(temporary) / "frontend"
            executable.write_bytes(raw)
            result = run_command(
                {
                    "runtime": "native",
                    "executable_sha256": hashlib.sha256(raw).hexdigest(),
                },
                executable,
                {"request": True},
                label="native test command",
            )
        self.assertEqual(result, {"runtime": "native"})


if __name__ == "__main__":
    unittest.main()

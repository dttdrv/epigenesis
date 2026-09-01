"""End-to-end tests for the public compiler-v2 CLI."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from brainc._canonical import digest
from brainc.v2 import inline_storage, pack, policy_artifact, save, target_artifact
from brainc.v2._common import seal
from brainc.v2.target import DEV_DOMAIN, OP_UNIT_CREATE


ROOT = Path(__file__).parents[1]


class CompilerV2CliTests(unittest.TestCase):
    def _run(
        self, *arguments: str | Path, expected: int = 0
    ) -> subprocess.CompletedProcess[str]:
        result = subprocess.run(
            [sys.executable, "-m", "brainc", *(str(value) for value in arguments)],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(
            result.returncode,
            expected,
            msg=f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}",
        )
        return result

    def _build_chain(self, root: Path) -> dict[str, Path]:
        fasta = root / "input.fasta"
        fasta.write_text(
            ">arbitrary-record\nACGTRYSWKMBDHVNACGT\n",
            encoding="utf-8",
        )
        source = root / "source.json"
        self._run("compile-collection", fasta, "-o", source)

        target_path = root / "target.json"
        target = target_artifact(
            {
                "id": "org.example.cli-test",
                "abi_major": 1,
                "opsets": [{"domain": DEV_DOMAIN, "version": 1}],
                "unit_schemas": [
                    {
                        "id": "node",
                        "fields": [
                            {
                                "id": "state",
                                "type": {"dtype": "u64", "shape": []},
                                "unit": None,
                                "mutability": "state",
                                "numeric": {"kind": "integer"},
                            }
                        ],
                    }
                ],
                "edge_schemas": [],
                "ports": [],
                "rules": [],
            }
        )
        save(target, target_path)
        target_binding = {
            "id": target["contract"]["id"],
            "abi_major": target["contract"]["abi_major"],
            "contract_sha256": target["contract_sha256"],
        }

        outputs = [
            {
                "id": "state.count",
                "type": {"dtype": "u64", "shape": []},
                "unit": None,
                "axes": [],
                "storage": inline_storage(pack("u64", [3])),
            },
            {
                "id": "state.value",
                "type": {"dtype": "u64", "shape": []},
                "unit": None,
                "axes": [],
                "storage": inline_storage(pack("u64", [17])),
            },
        ]
        provider = {"name": "cli-test-provider", "version": "1"}
        model_identity = {
            "kind": "content-sha256",
            "value": digest(
                {
                    "format": "cli-test-provider-semantics",
                    "version": 1,
                    "meaning": "three nodes initialized to seventeen",
                }
            ),
        }
        manifest_path = root / "manifest.json"
        manifest = seal(
            {
                "format": "brainc.provider-manifest",
                "version": 2,
                "provider": provider,
                "model_identity": model_identity,
                "accepts": ["brain01.sequence-collection-ir/v1"],
                "outputs": [
                    {
                        key: output[key]
                        for key in ("id", "type", "unit", "axes")
                    }
                    for output in outputs
                ],
            }
        )
        save(manifest, manifest_path)

        request_path = root / "request.json"
        self._run(
            "make-request-v2",
            source,
            "--manifest",
            manifest_path,
            "--output-id",
            "state.count",
            "--output-id",
            "state.value",
            "-o",
            request_path,
        )
        request = json.loads(request_path.read_text(encoding="utf-8"))

        response_path = root / "response.json"
        response = seal(
            {
                "format": "brainc.prediction-response",
                "version": 2,
                "request_artifact_sha256": request["artifact_sha256"],
                "provider": provider,
                "model_identity": model_identity,
                "outputs": outputs,
            }
        )
        save(response, response_path)

        policy_path = root / "policy.json"
        policy = policy_artifact(
            "cli-test-development",
            target_binding,
            [
                {"id": "state.count", "from_output": "state.count"},
                {"id": "state.value", "from_output": "state.value"},
            ],
            [
                {
                    "id": "create.nodes",
                    "op": OP_UNIT_CREATE,
                    "version": 1,
                    "schema": "node",
                    "count": "state.count",
                    "initializers": [
                        {"field": "state", "tensor": "state.value"}
                    ],
                }
            ],
        )
        save(policy, policy_path)

        module_path = root / "module.json"
        self._run(
            "compile-v2",
            source,
            "--manifest",
            manifest_path,
            "--request",
            request_path,
            "--response",
            response_path,
            "--policy",
            policy_path,
            "--target",
            target_path,
            "-o",
            module_path,
        )
        return {
            "fasta": fasta,
            "source": source,
            "manifest": manifest_path,
            "request": request_path,
            "response": response_path,
            "policy": policy_path,
            "target": target_path,
            "module": module_path,
        }

    def test_complete_cli_round_trip_and_independent_checks(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = self._build_chain(Path(temporary))
            module = json.loads(paths["module"].read_text(encoding="utf-8"))
            self.assertEqual(module["module"]["budgets"]["units"], 3)
            self.assertEqual(module["module"]["budgets"]["operations"], 1)

            checks: list[tuple[str | Path, ...]] = [
                ("source", paths["source"]),
                ("manifest", paths["manifest"]),
                ("request", paths["request"]),
                (
                    "response",
                    paths["response"],
                    "--source",
                    paths["source"],
                ),
                ("policy", paths["policy"]),
                ("target", paths["target"]),
                (
                    "module",
                    paths["module"],
                    "--target",
                    paths["target"],
                ),
            ]
            for arguments in checks:
                with self.subTest(kind=arguments[0]):
                    result = self._run("check-v2", *arguments)
                    self.assertTrue(json.loads(result.stdout)["valid"])

            report_path = Path(temporary) / "report.json"
            result = self._run(
                "validate-v2",
                "--source-input",
                paths["fasta"],
                "--source-artifact",
                paths["source"],
                "--manifest",
                paths["manifest"],
                "--request",
                paths["request"],
                "--response",
                paths["response"],
                "--policy",
                paths["policy"],
                "--target",
                paths["target"],
                "--module",
                paths["module"],
                "--report",
                report_path,
            )
            report = json.loads(result.stdout)
            self.assertTrue(report["valid"])
            self.assertEqual(report["passed_checks"], 4)
            self.assertEqual(
                json.loads(report_path.read_text(encoding="utf-8")), report
            )

    def test_validate_v2_rejects_tampering_with_a_sealed_failure_report(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = self._build_chain(Path(temporary))
            module = json.loads(paths["module"].read_text(encoding="utf-8"))
            module["module"]["budgets"]["units"] += 1
            tampered = Path(temporary) / "tampered-module.json"
            tampered.write_text(json.dumps(module), encoding="utf-8")
            report_path = Path(temporary) / "failure-report.json"
            result = self._run(
                "validate-v2",
                "--source-input",
                paths["fasta"],
                "--source-artifact",
                paths["source"],
                "--manifest",
                paths["manifest"],
                "--request",
                paths["request"],
                "--response",
                paths["response"],
                "--policy",
                paths["policy"],
                "--target",
                paths["target"],
                "--module",
                tampered,
                "--report",
                report_path,
                expected=3,
            )
            report = json.loads(result.stdout)
            self.assertFalse(report["valid"])
            self.assertEqual(len(report["report_sha256"]), 64)
            self.assertEqual(
                json.loads(report_path.read_text(encoding="utf-8")), report
            )

    def test_check_v2_requires_binding_artifacts(self) -> None:
        response = self._run(
            "check-v2", "response", "response.json", expected=2
        )
        self.assertIn("--source", response.stderr)
        module = self._run("check-v2", "module", "module.json", expected=2)
        self.assertIn("--target", module.stderr)

    def test_compile_v2_reports_malformed_request_source_identity_cleanly(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = self._build_chain(root)
            request = json.loads(paths["request"].read_text(encoding="utf-8"))
            request["source"]["version"] = []
            malformed = root / "malformed-request.json"
            save(
                seal(
                    {
                        key: value
                        for key, value in request.items()
                        if key != "artifact_sha256"
                    }
                ),
                malformed,
            )
            output = root / "malformed-module.json"
            result = self._run(
                "compile-v2",
                paths["source"],
                "--manifest",
                paths["manifest"],
                "--request",
                malformed,
                "--response",
                paths["response"],
                "--policy",
                paths["policy"],
                "--target",
                paths["target"],
                "--output",
                output,
                expected=2,
            )
            self.assertIn("format/version is unsupported", result.stderr)
            self.assertNotIn("Traceback", result.stderr)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()

"""End-to-end tests for the public GenBank compiler and validator CLI."""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).parent / "data"


class GenBankCliTests(unittest.TestCase):
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

    def test_official_genbank_round_trip_writes_a_sealed_report(self) -> None:
        source = DATA / "U49845.1.gb"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            artifact_path = root / "compiled.json"
            report_path = root / "validation.json"

            compiled = self._run(
                "compile-genbank", source, "-o", artifact_path
            )
            self.assertEqual(compiled.stdout, "")
            artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
            self.assertEqual(artifact["version"], 2)
            self.assertEqual(artifact["bio_ir"]["records"][0]["version"], "U49845.1")

            validated = self._run(
                "validate-genbank",
                source,
                artifact_path,
                "--report",
                report_path,
            )
            report = json.loads(validated.stdout)
            self.assertTrue(report["valid"])
            self.assertEqual(report["version"], 2)
            self.assertEqual(report["summary"]["sequence_bases"], 5028)
            self.assertEqual(
                json.loads(report_path.read_text(encoding="utf-8")),
                report,
            )

            stdout_only = self._run(
                "validate-genbank", source, artifact_path
            )
            self.assertEqual(json.loads(stdout_only.stdout), report)

    def test_tamper_returns_three_and_creates_no_success_report(self) -> None:
        source = DATA / "U49845.1.gb"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            artifact_path = root / "compiled.json"
            self._run("compile-genbank", source, "-o", artifact_path)
            artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
            artifact["bio_ir"]["records"][0]["definition"] = "tampered"
            tampered_path = root / "tampered.json"
            tampered_path.write_text(json.dumps(artifact), encoding="utf-8")
            report_path = root / "must-not-exist.json"

            failure = self._run(
                "validate-genbank",
                source,
                tampered_path,
                "--report",
                report_path,
                expected=3,
            )

            self.assertEqual(failure.stdout, "")
            self.assertIn("brainc: validation failed:", failure.stderr)
            self.assertFalse(report_path.exists())

    @unittest.skipUnless(hasattr(os, "symlink"), "symlinks are unavailable")
    def test_genbank_paths_reject_symlink_ingress_and_report_targets(self) -> None:
        source = DATA / "U49845.1.gb"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            linked_source = root / "source-link.gb"
            linked_source.symlink_to(source)
            artifact_path = root / "compiled.json"

            compile_failure = self._run(
                "compile-genbank",
                linked_source,
                "-o",
                artifact_path,
                expected=2,
            )
            self.assertIn("regular non-linked file", compile_failure.stderr)
            self.assertFalse(artifact_path.exists())

            self._run("compile-genbank", source, "-o", artifact_path)
            linked_artifact = root / "artifact-link.json"
            linked_artifact.symlink_to(artifact_path)
            failed_report = root / "failed-validation.json"
            validation_failure = self._run(
                "validate-genbank",
                source,
                linked_artifact,
                "--report",
                failed_report,
                expected=3,
            )
            self.assertIn("regular non-linked file", validation_failure.stderr)
            self.assertFalse(failed_report.exists())

            preserved = root / "preserved.txt"
            preserved.write_text("preserve", encoding="utf-8")
            linked_report = root / "report-link.json"
            linked_report.symlink_to(preserved)
            report_failure = self._run(
                "validate-genbank",
                source,
                artifact_path,
                "--report",
                linked_report,
                expected=3,
            )
            self.assertIn("regular non-linked file", report_failure.stderr)
            self.assertEqual(preserved.read_text(encoding="utf-8"), "preserve")


if __name__ == "__main__":
    unittest.main()

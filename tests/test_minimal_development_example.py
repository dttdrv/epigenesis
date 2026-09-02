from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).parents[1]
EXAMPLE = ROOT / "examples" / "minimal-development"


class MinimalDevelopmentExampleTests(unittest.TestCase):
    def test_committed_inputs_compile_and_validate_with_recorded_identities(self) -> None:
        self.assertTrue(
            (EXAMPLE / "expected.json").is_file(),
            "the committed example must include expected.json",
        )
        expected = json.loads((EXAMPLE / "expected.json").read_text(encoding="utf-8"))
        interpretation = EXAMPLE / "interpretation"

        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            source = output / "source"
            development = output / "development"
            validation = output / "validation.json"

            self._run(
                "compile-source",
                "--profile",
                "fasta-reference-dna/v1",
                "--sequence",
                str(EXAMPLE / "sequence.fasta"),
                "--wrapper",
                "identity",
                "--output",
                str(source),
            )
            self._run(
                "compile-development",
                str(source),
                "--manifest",
                str(interpretation / "manifest.json"),
                "--request",
                str(interpretation / "request.json"),
                "--response",
                str(interpretation / "response.json"),
                "--policy",
                str(interpretation / "lowering-policy.json"),
                "--target",
                str(interpretation / "target-contract.json"),
                "--output",
                str(development),
            )
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "brainc.validator_development",
                    str(source),
                    str(development),
                    "--source-input",
                    f"sequence={EXAMPLE / 'sequence.fasta'}",
                    "--output",
                    str(validation),
                ],
                cwd=ROOT,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

            observed = {
                "source_artifact_sha256": self._load(source / "source.json")[
                    "artifact_sha256"
                ],
                "development_bundle_artifact_sha256": self._load(
                    development / "bundle.json"
                )["artifact_sha256"],
                "development_module_artifact_sha256": self._load(
                    development / "development_module.json"
                )["artifact_sha256"],
                "validation_report_sha256": self._load(validation)["report_sha256"],
            }
            self.assertEqual(observed, expected)
            self.assertTrue(self._load(validation)["valid"])

    def _run(self, *arguments: str) -> None:
        result = subprocess.run(
            [sys.executable, "-m", "brainc", *arguments],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    @staticmethod
    def _load(path: Path) -> dict:
        return json.loads(path.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()

"""End-to-end test for the public GenBank feature-state compiler CLI."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).parent / "data"


class INSDCGraphCliTests(unittest.TestCase):
    def _run(self, *arguments: str | Path) -> subprocess.CompletedProcess[str]:
        result = subprocess.run(
            [sys.executable, "-m", "brainc", *(str(value) for value in arguments)],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(
            result.returncode,
            0,
            msg=f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}",
        )
        return result

    def test_official_genbank_reaches_reference_only_feature_state(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source_artifact = root / "source.json"
            state_directory = root / "state"
            self._run(
                "compile-genbank",
                DATA / "U49845.1.gb",
                "-o",
                source_artifact,
            )
            result = self._run(
                "compile-insdc-graph",
                source_artifact,
                "-o",
                state_directory,
            )
            summary = json.loads(result.stdout)
            self.assertEqual(Path(summary["bundle"]), state_directory / "bundle.json")
            self.assertEqual(len(list(state_directory.glob("*.json"))), 10)
            index = json.loads(
                (state_directory / "bundle.json").read_text(encoding="utf-8")
            )
            self.assertEqual(len(index["artifacts"]), 9)
            self.assertTrue(
                all(
                    set(reference) == {"artifact_sha256"}
                    for reference in index["artifacts"].values()
                )
            )
            module = json.loads(
                (state_directory / "development_module.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(module["module"]["budgets"]["units"], 6)
            self.assertEqual(module["module"]["budgets"]["edges"], 0)


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "epigenesis_cli.py"
RELEASE_CLI = ROOT / "tools/release/epigenesis.py"


class EpigenesisCliTests(unittest.TestCase):
    def _run(self, *args: str, env: dict[str, str] | None = None, cwd: Path | None = None) -> subprocess.CompletedProcess:
        return subprocess.run(
            [sys.executable, str(CLI), *args],
            cwd=cwd or ROOT,
            env=env,
            capture_output=True,
            text=True,
        )

    def test_help_lists_all_unified_commands(self):
        res = self._run("--help")
        self.assertEqual(res.returncode, 0, res.stderr)
        for cmd in ("studio", "circuit", "scan", "predict", "evidence", "patterning", "construction", "check"):
            self.assertIn(cmd, res.stdout)

    def test_missing_repo_fails_cleanly_when_isolated(self):
        with tempfile.TemporaryDirectory() as tmp:
            isolated = Path(tmp) / "epigenesis_cli.py"
            isolated.write_text(CLI.read_text())
            env = {k: v for k, v in os.environ.items() if k != "EPIGENESIS_HOME"}
            res = subprocess.run(
                [sys.executable, str(isolated), "circuit", "--output", str(Path(tmp) / "out")],
                cwd=tmp,
                env=env,
                capture_output=True,
                text=True,
            )
            self.assertEqual(res.returncode, 1)
            self.assertIn("could not locate", res.stderr)

    def test_circuit_compare_and_scan_with_explicit_repo_and_env(self):
        ref = "ACTG" * 5 + "GACCACCCA" + "TGCA" * 5
        var = "ACTG" * 5 + "GACCACTTTCCA" + "TGCA" * 5
        with tempfile.TemporaryDirectory() as tmp:
            troot = Path(tmp)
            fasta = troot / "enhancer.fasta"
            fasta.write_text(f">enh\n{ref}\n")
            cmp_out = troot / "compare_out"
            scan_out = troot / "scan_out"

            res_cmp = self._run(
                "--repo", str(ROOT),
                "circuit", "--reference", ref, "--variant", var,
                "--target", "Olig2", "--output", str(cmp_out),
                cwd=troot,
            )
            self.assertEqual(res_cmp.returncode, 0, res_cmp.stderr)
            comp = json.loads((cmp_out / "compare.json").read_text())
            self.assertEqual(comp["edit_kind"], "insertion")
            self.assertEqual(comp["target"], "Olig2")
            self.assertEqual(comp["compiled"]["reference_length"], len(ref))
            self.assertEqual(comp["compiled"]["variant_length"], len(var))

            env = dict(os.environ, EPIGENESIS_HOME=str(ROOT))
            res_scan = self._run(
                "scan", str(fasta), "--target", "Both", "--output", str(scan_out),
                env=env, cwd=troot,
            )
            self.assertEqual(res_scan.returncode, 0, res_scan.stderr)
            scanned = json.loads((scan_out / "scan.json").read_text())
            self.assertEqual(scanned["target"], "Both")
            self.assertEqual(scanned["length"], len(ref))
            self.assertEqual(len(scanned["edits"]), 3 * len(ref))

    def test_validation_errors_exit_with_code_two(self):
        with tempfile.TemporaryDirectory() as tmp:
            dummy_in = Path(tmp) / "in.json"
            dummy_in.write_text("[]")
            dummy_out = Path(tmp) / "out.json"
            bad_cases = (
                ("circuit", "--reference", "ACGTACGTACGT"),
                ("predict", "--input", str(dummy_in), "--output", str(dummy_out), "--families", "f.json"),
                ("predict", "--input", str(dummy_in), "--output", str(dummy_out), "--activity", "--budget", "3"),
                ("predict", "--input", str(dummy_in), "--output", str(dummy_out), "--kreimer", "--device", "cpu"),
            )
            for args in bad_cases:
                with self.subTest(args=args):
                    res = self._run(*args)
                    self.assertEqual(res.returncode, 2, res.stderr)

    def test_release_launcher_new_subcommand_help(self):
        for subcmd in ("circuit", "scan", "construction"):
            with self.subTest(subcmd=subcmd):
                res = subprocess.run(
                    [sys.executable, str(RELEASE_CLI), subcmd, "--help"],
                    cwd=ROOT,
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(res.returncode, 0, res.stderr)


        res_bad = subprocess.run(
            [sys.executable, str(RELEASE_CLI), "circuit", "--reference", "ACGTACGTAC"],
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        self.assertEqual(res_bad.returncode, 2)


if __name__ == "__main__":
    unittest.main()

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
        for cmd in ("studio", "circuit", "scan", "radar", "predict", "evidence", "patterning", "construction", "check"):
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
                ("radar", "--fasta", "a.fasta", "--variants", str(dummy_in)),
                ("radar", "--benchmark", "--fasta", "a.fasta"),
                ("predict", "--input", str(dummy_in), "--output", str(dummy_out), "--families", "f.json"),
                ("predict", "--input", str(dummy_in), "--output", str(dummy_out), "--activity", "--budget", "3"),
                ("predict", "--input", str(dummy_in), "--output", str(dummy_out), "--kreimer", "--device", "cpu"),
            )
            for args in bad_cases:
                with self.subTest(args=args):
                    res = self._run(*args)
                    self.assertEqual(res.returncode, 2, res.stderr)

    def test_release_launcher_new_subcommand_help(self):
        for subcmd in ("circuit", "scan", "radar", "construction"):
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

    def test_circuit_list_presets_and_motifs(self):
        res_c = self._run("circuit", "--list-circuits")
        self.assertEqual(res_c.returncode, 0, res_c.stderr)
        presets = json.loads(res_c.stdout)
        self.assertEqual(set(presets), {"neural-tube", "mutual-toggle", "tri-lineage-switch"})

        res_m = self._run("circuit", "--list-motifs")
        self.assertEqual(res_m.returncode, 0, res_m.stderr)
        motifs = json.loads(res_m.stdout)
        self.assertEqual(motifs["hallikas_gli"], ["GLI2", "GLI3"])
        self.assertEqual(motifs["jaspar_count"], 532)

    def test_circuit_and_scan_generalized_jaspar_and_custom_specs(self):
        ref = "ACTG" * 4 + "TTCACGCATGAGTT" + "TGCA" * 2 + "GCTTTCGAGGAAT" + "ACGT" * 2
        var = "ACTG" * 4 + "TTCACAAATGAGTT" + "TGCA" * 2 + "GCTTTCGAGGAAT" + "ACGT" * 2
        with tempfile.TemporaryDirectory() as tmp:
            troot = Path(tmp)
            fasta = troot / "enh.fasta"
            fasta.write_text(f">enh\n{ref}\n")

            tri_out = troot / "tri_cmp"
            res_tri = self._run(
                "circuit", "--circuit", "tri-lineage-switch", "--target", "All",
                "--reference", ref, "--variant", var, "--output", str(tri_out),
            )
            self.assertEqual(res_tri.returncode, 0, res_tri.stderr)
            tri = json.loads((tri_out / "compare.json").read_text())
            self.assertEqual(tri["circuit"], "tri-lineage-switch")
            self.assertEqual(tri["genes"], ["Progenitor", "Excitatory", "Inhibitory"])
            self.assertEqual(tri["proteins"], ["MA0069.1", "MA0463.3"])
            self.assertEqual(tri["target"], "All")

            tog_out = troot / "tog_scan"
            res_tog = self._run(
                "scan", str(fasta), "--circuit", "mutual-toggle",
                "--protein", "MA0069.1", "--target", "GeneA", "--output", str(tog_out),
            )
            self.assertEqual(res_tog.returncode, 0, res_tog.stderr)
            tog = json.loads((tog_out / "scan.json").read_text())
            self.assertEqual(tog["circuit"], "mutual-toggle")
            self.assertEqual(tog["protein"], "MA0069.1")
            self.assertEqual(tog["target"], "GeneA")

            custom_circuit = troot / "custom_circuit.json"
            custom_circuit.write_text(json.dumps({
                "name": "cli-custom-2",
                "genes": ["X", "Y"],
                "basal": [1.0, 2.0],
                "input": [150.0, 20.0],
                "tf_binding": [10.0, 2.0],
                "repression": [[0.0, 30.0], [30.0, 0.0]],
                "cooperativity": 2.0,
                "production": [2.0, 2.0],
                "decay": [2.0, 2.0],
                "default_target": "X",
                "default_signal": 0.25,
                "default_protein": "MA0069.1",
            }))
            custom_pwm = troot / "custom_pwm.json"
            custom_pwm.write_text(json.dumps({
                "id": "CLI-PWM-5",
                "kind": "weights",
                "matrix": {
                    "A": [1.0, 0.1, 0.1, 1.0, 0.1],
                    "C": [0.1, 1.0, 0.1, 0.1, 1.0],
                    "G": [0.1, 0.1, 1.0, 0.1, 0.1],
                    "T": [0.1, 0.1, 0.1, 0.1, 0.1],
                },
            }))
            cust_out = troot / "cust_cmp"
            res_cust = self._run(
                "circuit", "--circuit", str(custom_circuit), "--protein", str(custom_pwm),
                "--reference", "TTTTACGACGGGG", "--variant", "TTTTACTTCGGGG",
                "--output", str(cust_out),
            )
            self.assertEqual(res_cust.returncode, 0, res_cust.stderr)
            cust = json.loads((cust_out / "compare.json").read_text())
            self.assertEqual(cust["circuit"], "cli-custom-2")
            self.assertEqual(cust["protein"], "CLI-PWM-5")
            self.assertEqual(cust["genes"], ["X", "Y"])

    def test_radar_variants_and_fasta_cli(self):
        ref = "ACTG" * 5 + "GACCACCCA" + "TGCA" * 5
        with tempfile.TemporaryDirectory() as tmp:
            troot = Path(tmp)
            vjson = troot / "variants.json"
            vjson.write_text(json.dumps([
                {"id": "wt", "reference": ref, "sequence": ref},
                {"id": "mut_del", "reference": ref, "sequence": "ACTG" * 5 + "GA" + "TGCA" * 5},
            ]))
            vout = troot / "vout"
            res_v = self._run(
                "radar", "--variants", str(vjson), "--steps", "7",
                "--no-sequence-models", "--output", str(vout),
            )
            self.assertEqual(res_v.returncode, 0, res_v.stderr)
            rep_v = json.loads((vout / "variants_report.json").read_text())
            self.assertEqual(rep_v["variants_tested"], 2)
            self.assertFalse(rep_v["biological_acceptance"])
            self.assertTrue((vout / "wetlab_candidates.csv").is_file())

            fasta = troot / "enh.fasta"
            fasta.write_text(f">enh\n{ref}\n")
            fout = troot / "fout"
            res_f = self._run(
                "radar", "--fasta", str(fasta), "--steps", "7", "--shuffles", "1",
                "--no-sequence-models", "--output", str(fout),
            )
            self.assertEqual(res_f.returncode, 0, res_f.stderr)
            rep_f = json.loads((fout / "radar_report.json").read_text())
            self.assertGreaterEqual(rep_f["scan"]["category_counts"]["cryptic_bifurcation"], 0)
            self.assertTrue(rep_f["causal_controls"]["all_passed"])
            self.assertTrue((fout / "index.html").is_file())


if __name__ == "__main__":
    unittest.main()

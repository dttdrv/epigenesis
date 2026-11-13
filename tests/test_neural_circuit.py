from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
CIRCUIT_PATH = ROOT / "examples/neural-circuit/circuit.py"
PATTERNING_PATH = ROOT / "examples/neural-patterning/patterning.py"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


class NeuralCircuitTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pat = _load("neural_patterning_for_circuit", PATTERNING_PATH)
        cls.cir = _load("neural_circuit", CIRCUIT_PATH)

    def test_zero_pseudocount_partition_matches_exact_hallikas_profile(self):
        core = "GACCACCCA"
        for mut in ("GAACACCCA", "GACCACACA", "TGGGTGTTC"):
            q_cir = self.cir.affinity_ratio(core, mut, "GLI2", epsilon=0.0)
            oriented = mut if mut.startswith("G") else self.pat.reverse_complement(mut)
            q_pat = self.pat.profile(oriented, "GLI2")["profile_ratio"]
            self.assertAlmostEqual(q_cir, q_pat, places=14)

    def test_gain_of_affinity_reverses_occupancy_sign_opposite_to_loss(self):
        # Embed a weak suboptimal GLI site and strengthen it to consensus (q > 1)
        weak = "AAAA" + "GAACACCCA" + "TTTT"
        strong = "AAAA" + "GACCACCCA" + "TTTT"
        q_gain = self.cir.affinity_ratio(weak, strong, "GLI2")
        self.assertGreater(q_gain, 10.0)
        model = self.pat.load_data()["model"]
        state = (0.2, 0.4, 0.1, 0.5)
        # Gain of affinity represses more below S=0.1 and activates more above S=0.1
        self.assertLess(self.pat.occupancy(model, state, 0.02, q_gain)[2], self.pat.occupancy(model, state, 0.02, 1.0)[2])
        self.assertGreater(self.pat.occupancy(model, state, 0.50, q_gain)[2], self.pat.occupancy(model, state, 0.50, 1.0)[2])
        self.assertAlmostEqual(self.pat.occupancy(model, state, 0.10, q_gain)[2], self.pat.occupancy(model, state, 0.10, 1.0)[2], places=14)

    def test_scan_detects_cryptic_bifurcation_and_neutral_flank_control(self):
        enhancer = "ACTG" * 10 + "GACCACCCA" + "TGCA" * 10
        scan = self.cir.scan_enhancer(enhancer, signal=0.20, protein="GLI2")
        self.assertGreater(scan["category_counts"]["cryptic_bifurcation"], 0)
        # Position 1 base edit far from GLI site must be neutral with zero fate flip
        first = scan["edits"][0]
        self.assertEqual(first["category"], "neutral")
        self.assertFalse(first["fate_flipped"])
        self.assertAlmostEqual(first["affinity_ratio"], 1.0, places=6)

    def test_brainc_bundle_roundtrip_and_tamper_detection(self):
        ref = "ACTG" * 5 + "GACCACCCA" + "TGCA" * 5
        var = "ACTG" * 5 + "GAACACCCA" + "TGCA" * 5
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "bundle"
            receipt = self.cir.compile_circuit_bundle(ref, var, output=out)
            loaded = self.cir.verify_circuit_bundle(out, receipt)
            self.assertLess(loaded["affinity_ratio"], 0.1)
            # Mutate FASTA without updating contract -> validation or mismatch rejects
            (out / "pair.fasta").write_text(f">reference\n{ref}\n>variant\n{ref}\n")
            with self.assertRaises(ValueError):
                self.cir.verify_circuit_bundle(out, receipt)

    def test_demo_end_to_end(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "demo"
            report = self.cir.demo(out)
            self.assertTrue(report["all_checks_passed"], report["checks"])
            self.assertGreater(report["scan"]["category_counts"]["cryptic_bifurcation"], 0)
            self.assertTrue((out / "index.html").is_file())
            self.assertTrue((out / "report.json").is_file())


if __name__ == "__main__":
    unittest.main()

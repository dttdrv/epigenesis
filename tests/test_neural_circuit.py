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

    def test_multi_target_ratios_and_validation(self):
        model = self.pat.load_data()["model"]
        state = (0.15, 0.35, 0.25, 0.45)
        self.assertEqual(self.cir.resolve_target_ratios("Olig2", 0.25), {"Olig2": 0.25})
        self.assertEqual(self.cir.resolve_target_ratios("Both", 0.25), {"Nkx2.2": 0.25, "Olig2": 0.25})
        self.assertEqual(
            self.cir.occupancy_ratios(model, state, 0.2, {"Nkx2.2": 0.4}),
            self.pat.occupancy(model, state, 0.2, 0.4),
        )
        olig_occ = self.cir.occupancy_ratios(model, state, 0.2, {"Olig2": 0.1})
        wt_occ = self.pat.occupancy(model, state, 0.2, 1.0)
        self.assertAlmostEqual(olig_occ[2], wt_occ[2], places=14)
        self.assertLess(olig_occ[1], wt_occ[1])
        with self.assertRaises(ValueError):
            self.cir.resolve_target_ratios("Pax6", 0.5)
        with self.assertRaises(ValueError):
            self.cir.resolve_target_ratios("Nkx2.2", 0.0)

    def test_scan_enhancer_supports_olig2_and_both_targets(self):
        enhancer = "ACTG" * 10 + "GACCACCCA" + "TGCA" * 10
        scan_olig = self.cir.scan_enhancer(enhancer, signal=0.20, target="Olig2")
        scan_both = self.cir.scan_enhancer(enhancer, signal=0.20, target="Both")
        self.assertEqual(scan_olig["target"], "Olig2")
        self.assertEqual(scan_both["target"], "Both")
        # Target delta tracks Olig2 when target='Olig2'
        self.assertTrue(all(e["target_delta"] == e["olig2_delta"] for e in scan_olig["edits"]))

    def test_compare_variants_mnv_insertion_deletion_and_brainc_roundtrip(self):
        ref = "ACTG" * 5 + "GACCACCCA" + "TGCA" * 5
        mnv = "ACTG" * 5 + "GAACACCAA" + "TGCA" * 5
        ins = "ACTG" * 5 + "GACCACTTTCCA" + "TGCA" * 5
        dele = "ACTG" * 5 + "GA" + "TGCA" * 5
        for kind, var in (("mnv", mnv), ("insertion", ins), ("deletion", dele)):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as tmp:
                comp = self.cir.compare_variants(ref, var, compile_dir=Path(tmp) / "bundle")
                self.assertEqual(comp["edit_kind"], kind)
                self.assertEqual(comp["compiled"]["reference_length"], len(ref))
                self.assertEqual(comp["compiled"]["variant_length"], len(var))
                self.assertLess(comp["affinity_ratio"], 0.15)
                self.assertEqual(comp["reversion_max_difference"], 0.0)
                # Tampering pair.fasta length triggers refusal
                (Path(tmp) / "bundle/pair.fasta").write_text(f">reference\n{ref}\n>variant\n{ref}\n")
                with self.assertRaises(ValueError):
                    self.cir.verify_circuit_bundle(Path(tmp) / "bundle", comp["compiled"]["contract_sha256"])

    def test_second_site_rescue_on_chr18_enhancer(self):
        reg_mod = _load("neural_regulation_for_circuit", ROOT / "examples/neural-regulation/regulation.py")
        rows, _ = reg_mod.load_data()
        chr18 = next(r.reference for r in rows if r.locus == "chr18:56475717-56475888")
        scan = self.cir.scan_enhancer(chr18, signal=0.20, protein="GLI2")
        self.assertGreater(scan["category_counts"]["cryptic_bifurcation"], 0)
        featured = scan["cryptic_switches"][0]
        rescue = self.cir._find_second_site_rescue(chr18, featured, "GLI2", self.cir.DEFAULT_EPSILON, signal=0.20)
        self.assertIsNotNone(rescue)
        self.assertTrue(rescue["fate_rescued"], rescue)
        self.assertEqual(rescue["double_fate"], "Nkx2.2")

    def test_multi_locus_benchmark(self):
        bench = self.cir.multi_locus_benchmark()
        self.assertEqual(bench["format"], "epigenesis.neural-circuit-benchmark/v1")
        self.assertEqual(bench["loci_tested"], len(self.cir.CURATED_LOCI))
        self.assertEqual(bench["total_cryptic_bifurcations"], 33)
        self.assertEqual(bench["total_measured_variants"], 85)
        rescues = [loc["second_site_rescue"] for loc in bench["loci"] if loc["second_site_rescue"]]
        self.assertTrue(any(r["fate_rescued"] for r in rescues))

    def test_bridge_to_construction_with_cx3d_and_lif(self):
        import shutil
        if shutil.which("javac") is None or shutil.which("patch") is None:
            self.skipTest("javac and patch are required for CX3D construction bridge test")
        ref = "ACTG" * 5 + "GACCACCCA" + "TGCA" * 5
        var = "ACTG" * 5 + "GAACACCCA" + "TGCA" * 5
        with tempfile.TemporaryDirectory() as tmp:
            res = self.cir.bridge_to_construction(ref, var, Path(tmp) / "bridge", dev_duration=1.5)
            self.assertFalse(res["biological_acceptance"])
            self.assertTrue(res["all_checks_passed"], res)
            self.assertTrue(res["checks"]["all_compiled_and_lif_controls_pass"])
            self.assertTrue(res["checks"]["variant_diverges_morphology"])
            self.assertTrue(res["checks"]["allele_rescue_matches_reference_exactly"])


if __name__ == "__main__":
    unittest.main()

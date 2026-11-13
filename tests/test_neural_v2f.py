from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
RADAR_PATH = ROOT / "examples/neural-v2f/radar.py"
CIRCUIT_PATH = ROOT / "examples/neural-circuit/circuit.py"
REG_PATH = ROOT / "examples/neural-regulation/regulation.py"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


class NeuralV2FRadarTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cir = _load("neural_circuit_for_v2f", CIRCUIT_PATH)
        cls.radar = _load("neural_v2f_radar", RADAR_PATH)
        cls.reg = _load("neural_reg_for_v2f", REG_PATH)
        rows, _ = cls.reg.load_data()
        cls.chr15 = next(r.reference for r in rows if r.locus == "chr15:90936956-90937127")
        cls.chr18 = next(r.reference for r in rows if r.locus == "chr18:56475717-56475888")

    def test_dinucleotide_shuffle_exact_eulerian_invariant(self):
        for seq in (self.chr15, self.chr18, "ACGTAGGTCATACG"):
            s0 = self.radar.dinucleotide_shuffle(seq, seed=17)
            s0_repeat = self.radar.dinucleotide_shuffle(seq, seed=17)
            s1 = self.radar.dinucleotide_shuffle(seq, seed=99)
            self.assertEqual(s0, s0_repeat)
            self.assertNotEqual(s0, seq)
            self.assertNotEqual(s0, s1)
            self.assertTrue(self.radar.verify_dinucleotide_counts(seq, s0))
            self.assertTrue(self.radar.verify_dinucleotide_counts(seq, s1))
        # Periodic ACGT repeat has a unique Eulerian dinucleotide walk
        periodic = "ACGTACGTACGT"
        self.assertEqual(self.radar.dinucleotide_shuffle(periodic, seed=17), periodic)
        with self.assertRaises(ValueError):
            self.radar.dinucleotide_shuffle("ACGTN")

    def test_scramble_motif_pwm_control(self):
        for motif in ("GLI2", "MA0069.1"):
            orig = self.cir.load_motif_profile(motif, epsilon=0.0)
            scram = self.radar.scramble_motif_pwm(motif, seed=42)
            scram_prof = self.cir.load_motif_profile(scram, epsilon=0.0)
            self.assertEqual(scram_prof["width"], orig["width"])
            self.assertNotEqual(scram_prof["consensus"], orig["consensus"])

    def test_regime_scan_detects_cryptic_switch_and_exact_tipping_point_shift(self):
        scan = self.radar.scan_enhancer_regime(
            self.chr15, signal_min=0.06, signal_max=0.85, signal_steps=9, include_sequence_models=False
        )
        self.assertIsNotNone(scan["wt_critical_signal"])
        self.assertGreater(scan["category_counts"]["cryptic_bifurcation"], 0)
        top = scan["cryptic_switches"][0]
        self.assertEqual(top["id"], "62:A>C")
        self.assertLess(top["max_abs_reporter_delta"], 0.02)
        self.assertGreater(top["delta_s_crit"], 0.08)

        # Verify bisected tipping point directly against independent RK4 ODE solves
        ctx = self.cir._resolve_context(None, "GLI2", "Nkx2.2", 0.20, None, None)
        s_crit = top["var_critical_signal"]
        lo_fate = self.radar._solve_at_signal(ctx, {"Nkx2.2": top["affinity_ratio"]}, max(0.0, s_crit - 0.005))[1]
        hi_fate = self.radar._solve_at_signal(ctx, {"Nkx2.2": top["affinity_ratio"]}, min(1.0, s_crit + 0.005))[1]
        self.assertNotEqual(lo_fate, hi_fate)
        self.assertEqual(hi_fate, "Nkx2.2")

    def test_dual_rescue_trans_and_regime_cis_epistasis(self):
        # 1) Trans signal rescue on chr15 62:A>C restores Nkx2.2
        with tempfile.TemporaryDirectory() as tmp:
            rep = self.radar.run_radar(
                self.chr15, Path(tmp) / "chr15", label="chr15",
                n_shuffles=2, top_k=3, include_sequence_models=False,
            )
            self.assertFalse(rep["biological_acceptance"])
            self.assertTrue(rep["causal_controls"]["all_passed"])
            trans = rep["featured_rescues"]["trans_rescue"]
            self.assertTrue(trans["rescued"])
            self.assertEqual(trans["rescued_fate"], "Nkx2.2")
            self.assertGreater(trans["signal_shift"], 0.0)
            for fname in ("radar_report.json", "wetlab_candidates.csv", "wetlab_candidates.json", "index.html"):
                self.assertTrue((Path(tmp) / "chr15" / fname).is_file(), fname)

        # 2) Cis epistatic second-site rescue across full regime on chr18
        scan18 = self.radar.scan_enhancer_regime(
            self.chr18, signal_min=0.06, signal_max=0.85, signal_steps=9, include_sequence_models=False
        )
        self.assertGreater(len(scan18["cryptic_switches"]), 0)
        ctx = self.cir._resolve_context(None, "GLI2", "Nkx2.2", 0.20, None, None)
        wt_base = {
            "wt_profile": scan18["wt_profile"],
            "wt_transitions": scan18["wt_transitions"],
            "wt_primary_critical": scan18["wt_critical_signal"],
            "wt_preferred_fate": scan18["wt_preferred_fate"],
        }
        cand18 = dict(scan18["cryptic_switches"][0], wt_critical_signal=scan18["wt_critical_signal"])
        resc18 = self.radar._rescue_candidate(self.chr18, cand18, ctx, wt_base, scan18["signal_grid"])
        self.assertIsNotNone(resc18["cis_rescue"], resc18)
        self.assertGreaterEqual(resc18["cis_rescue"]["regime_reduction"], 0.0)

    def test_arbitrary_variant_table_snv_mnv_indel_and_negative_control(self):
        ref = "ACTG" * 5 + "GACCACCCA" + "TGCA" * 5
        rows = [
            {"id": "wt_same", "reference": ref, "sequence": ref},
            {"id": "snv_core", "reference": ref, "pos": 22, "alt": "A"},
            {"id": "mnv_core", "reference": ref, "sequence": "ACTG" * 5 + "GAACACCAA" + "TGCA" * 5},
            {"id": "ins_core", "reference": ref, "sequence": "ACTG" * 5 + "GACCACTTTCCA" + "TGCA" * 5},
            {"id": "del_core", "reference": ref, "sequence": "ACTG" * 5 + "GA" + "TGCA" * 5},
        ]
        with tempfile.TemporaryDirectory() as tmp:
            vjson = Path(tmp) / "variants.json"
            vjson.write_text(json.dumps(rows))
            out = Path(tmp) / "vout"
            rep = self.radar.analyze_variant_table(vjson, out, include_sequence_models=False)
            self.assertFalse(rep["biological_acceptance"])
            self.assertEqual(rep["variants_tested"], 5)
            self.assertTrue(rep["causal_controls"]["all_passed"])
            by_id = {r["id"]: r for r in rep["variants"]}
            self.assertEqual(by_id["wt_same"]["category"], "neutral")
            self.assertEqual(by_id["wt_same"]["regime_flip_fraction"], 0.0)
            self.assertEqual(by_id["ins_core"]["edit_kind"], "insertion")
            self.assertEqual(by_id["del_core"]["edit_kind"], "deletion")
            self.assertGreater(by_id["del_core"]["regime_flip_fraction"], 0.0)
            self.assertTrue(by_id["del_core"]["signal_locked"])
            with (out / "wetlab_candidates.csv").open() as fh:
                csv_rows = list(csv.DictReader(fh))
            self.assertEqual(len(csv_rows), 5)

        # Negative control: featureless repeat has 0 cryptic bifurcations
        repeat = "ACGTACGT" * 6
        neg_scan = self.radar.scan_enhancer_regime(repeat, include_sequence_models=False)
        self.assertEqual(neg_scan["category_counts"]["cryptic_bifurcation"], 0)

    def test_generalized_mutual_toggle_and_optional_sequence_scores(self):
        p_pax = self.cir.load_motif_profile("MA0069.1")["consensus"]
        enh = "ACTG" * 4 + p_pax + "TGCA" * 4
        with tempfile.TemporaryDirectory() as tmp:
            rep = self.radar.run_radar(
                enh, Path(tmp) / "tog", label="synthetic-pax6",
                circuit="mutual-toggle", protein="MA0069.1", target="GeneA",
                signal_min=0.06, signal_max=0.60, signal_steps=7,
                n_shuffles=2, top_k=3, include_sequence_models=False,
            )
            self.assertEqual(rep["circuit"], "mutual-toggle")
            self.assertEqual(rep["protein"], "MA0069.1")
            self.assertTrue(rep["causal_controls"]["all_passed"])
            self.assertEqual(rep["null_controls"]["dinucleotide_shuffle"]["n_shuffles"], 2)

        try:
            import numpy  # noqa: F401
        except ImportError:
            pass
        else:
            scan_seq = self.radar.scan_enhancer_regime(
                self.chr15, signal_min=0.06, signal_max=0.85, signal_steps=9, include_sequence_models=True
            )
            self.assertTrue(scan_seq["sequence_models_attached"])
            self.assertEqual(scan_seq["sequence_gate_model"], "kreimer_0_72h")
            self.assertIn("0.049", scan_seq["siraj_informational_note"])
            self.assertGreaterEqual(scan_seq["doubly_silent_count"], 1)
            c62 = next(e for e in scan_seq["cryptic_switches"] if e["id"] == "62:A>C")
            self.assertTrue(c62["doubly_silent"])
            self.assertIsNotNone(c62["siraj_score"])
            self.assertLess(abs(c62["kreimer_mpra_mean"]), 0.15)
            # Confirm Siraj is informational only: find any edit with |siraj_score| >= 0.15 and |kreimer_mpra_mean| < 0.15 if present
            for e in scan_seq["edits"]:
                self.assertEqual(e["sequence_silent"], abs(e["kreimer_mpra_mean"]) < 0.15)

    def test_benchmark_radar_two_loci(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "bench"
            bench = self.radar.benchmark_radar(
                out,
                loci=self.cir.CURATED_LOCI[:2],
                signal_steps=7,
                n_shuffles=1,
                include_sequence_models=False,
            )
            self.assertFalse(bench["biological_acceptance"])
            self.assertEqual(bench["loci_tested"], 2)
            self.assertGreater(bench["total_measured_variants"], 0)
            self.assertTrue(bench["causal_controls"]["all_passed"])
            self.assertTrue((out / "radar_benchmark.json").is_file())
            self.assertTrue((out / "measured_discordance.csv").is_file())

    def test_nonspecific_floor_suppresses_scrambled_pwm_false_positives(self):
        rows, _ = self.reg.load_data()
        chr1 = next(r.reference for r in rows if r.locus == self.cir.CURATED_LOCI[4])
        for locus_seq, min_unfloored in ((self.chr15, 2), (chr1, 1)):
            scan_floored = self.radar.scan_enhancer_regime(
                locus_seq, signal_min=0.06, signal_max=0.85, signal_steps=9, include_sequence_models=False
            )
            self.assertAlmostEqual(scan_floored["background_floor"], self.cir.DEFAULT_BACKGROUND_FLOOR, places=12)
            ctx = self.cir._resolve_context(None, "GLI2", "Nkx2.2", 0.455, None, None)
            null_floored = self.radar.run_null_controls(
                locus_seq, scan_floored, ctx=ctx, signal_min=0.06, signal_max=0.85, signal_steps=9, n_shuffles=0, seed=17,
            )
            self.assertEqual(null_floored["motif_pwm_scramble"]["scrambled_cryptic_count"], 0)

            scan_unfloored = self.radar.scan_enhancer_regime(
                locus_seq, signal_min=0.06, signal_max=0.85, signal_steps=9,
                background_floor=0.0, include_sequence_models=False,
            )
            null_unfloored = self.radar.run_null_controls(
                locus_seq, scan_unfloored, ctx=ctx, signal_min=0.06, signal_max=0.85, signal_steps=9,
                n_shuffles=0, background_floor=0.0, seed=17,
            )
            self.assertGreaterEqual(null_unfloored["motif_pwm_scramble"]["scrambled_cryptic_count"], min_unfloored)

    def test_default_shuffles_fast_solver_and_null_parity(self):
        self.assertGreaterEqual(self.radar.DEFAULT_SHUFFLES, 1000)
        ctx = self.cir._resolve_context(None, "GLI2", "Nkx2.2", 0.455, None, None)
        bound = sorted(ctx["bindings"])
        for sig in (0.06, 0.25, 0.70):
            for ratio in (0.1, 1.0, 4.0):
                rmap = {g: ratio for g in bound}
                ex = self.cir.simulate_circuit(
                    [(40.0, sig)], ratios=rmap, target=ctx["target"], dt=self.radar.SCAN_DT,
                    initial=(0.0,) * len(ctx["genes"]), circuit=ctx["circuit_arg"],
                    roles=ctx["roles"], mechanisms=ctx["mechanisms"],
                )[-1][1:]
                fast = self.radar._fast_end_state(
                    ctx["model"], ctx["genes"], rmap, sig, 40.0, self.radar.SCAN_DT,
                    (0.0,) * len(ctx["genes"]), ctx["roles"], ctx["mechanisms"],
                )
                for a, b in zip(ex, fast):
                    self.assertLess(abs(a - b), 1e-6)

        scan = self.radar.scan_enhancer_regime(
            self.chr15, signal_min=0.06, signal_max=0.85, signal_steps=9, include_sequence_models=False, ctx=ctx
        )
        wt_base = {
            "wt_profile": scan["wt_profile"],
            "wt_transitions": scan["wt_transitions"],
            "wt_primary_critical": scan["wt_critical_signal"],
            "wt_preferred_fate": scan["wt_preferred_fate"],
        }
        tf_profiles = [(lbl, self.cir.load_motif_profile(spec, self.cir.DEFAULT_EPSILON)) for lbl, spec in ctx["unique_specs"]]
        gene_lbls = [self.cir._motif_spec_label(ctx["bindings"][g]["protein"]) for g in bound]
        decoy = self.radar.dinucleotide_shuffle(self.chr15, seed=17)
        exact = self.radar.scan_enhancer_regime(
            decoy, signal_min=0.06, signal_max=0.85, signal_steps=9, include_sequence_models=False, ctx=ctx
        )
        fast_null = self.radar._scan_decoy_null_stats(
            decoy, ctx=ctx, signal_grid=scan["signal_grid"], wt_base=wt_base, tf_profiles=tf_profiles,
            bound_genes=bound, gene_lbls=gene_lbls, eff_floor=self.cir.DEFAULT_BACKGROUND_FLOOR,
            duration=40.0, reporter_threshold=self.radar.DEFAULT_REPORTER_SILENT,
        )
        self.assertEqual(fast_null["cryptic_count"], exact["category_counts"]["cryptic_bifurcation"])
        self.assertEqual(fast_null["point_silent_switches"], exact["point_silent_switches"])

    def test_heldout_freeze_integrity_and_results_recorded(self):
        freeze_dir = ROOT / "examples/neural-v2f/heldout-freeze"
        sys.path.insert(0, str(freeze_dir))
        try:
            import score as heldout_score
        finally:
            if str(freeze_dir) in sys.path:
                sys.path.remove(str(freeze_dir))
        freeze = heldout_score.verify_freeze()
        self.assertTrue(freeze["sealed_before_outcomes_scored"])
        self.assertFalse(freeze["biological_acceptance"])
        res = json.loads((freeze_dir / "results.json").read_text(encoding="utf-8"))
        self.assertFalse(res["biological_acceptance"])
        self.assertTrue(res["freeze_verified"])
        self.assertEqual(res["frozen_digests"], freeze["files"])
        self.assertEqual(res["unseen_test_loci"], 118)
        self.assertEqual(res["measured_test_variants_scored"], 1885)


if __name__ == "__main__":
    unittest.main()

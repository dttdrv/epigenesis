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

    def test_jaspar_catalog_and_custom_pwm_parity_strand_symmetry_and_validation(self):
        motifs = self.cir.list_available_motifs()
        self.assertEqual(motifs["hallikas_gli"], ["GLI2", "GLI3"])
        self.assertEqual(motifs["jaspar_count"], 532)
        self.assertIn("MA0069.1", motifs["jaspar_ids"])
        self.assertIn("MA0463.3", motifs["jaspar_ids"])

        for mid in ("MA0069.1", "MA0463.3"):
            prof = self.cir.load_motif_profile(mid)
            seq = "ACGT" * 4 + prof["consensus"] + "TGCA" * 4
            rc = self.pat.reverse_complement(seq)
            z_fwd = self.cir.boltzmann_affinity(seq, mid)["partition_sum"]
            z_rc = self.cir.boltzmann_affinity(rc, mid)["partition_sum"]
            self.assertGreater(z_fwd, 1.0)
            self.assertAlmostEqual(z_fwd, z_rc, places=11)
            self.assertAlmostEqual(self.cir.affinity_ratio(seq, rc, mid), 1.0, places=12)

        custom_weights = {
            "id": "CUSTOM-TOY-6",
            "kind": "weights",
            "matrix": {
                "A": [1.0, 0.1, 0.1, 1.0, 0.1, 0.1],
                "C": [0.1, 1.0, 0.1, 0.1, 1.0, 0.1],
                "G": [0.1, 0.1, 1.0, 0.1, 0.1, 1.0],
                "T": [0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
            },
        }
        custom_log_odds = {
            "id": "CUSTOM-LO-4",
            "kind": "log_odds",
            "matrix": {
                "A": [2.0, -3.0, -3.0, 2.0],
                "C": [-3.0, 2.0, -3.0, -3.0],
                "G": [-3.0, -3.0, 2.0, -3.0],
                "T": [-3.0, -3.0, -3.0, -3.0],
            },
        }
        for spec in (custom_weights, custom_log_odds):
            prof = self.cir.load_motif_profile(spec)
            seq = "TTTT" + prof["consensus"] + "GGGG"
            alt_pos = 4 + 1
            alt_base = "T"
            alt_seq = seq[:alt_pos] + alt_base + seq[alt_pos + 1 :]
            delta = self.cir._delta_partition_single(seq, alt_pos, alt_base, prof["ratios"], prof["width"])
            z_wt = self.cir.boltzmann_affinity(seq, spec)["partition_sum"]
            z_alt = self.cir.boltzmann_affinity(alt_seq, spec)["partition_sum"]
            self.assertAlmostEqual(z_wt + delta, z_alt, places=12)
            self.assertAlmostEqual(self.cir.affinity_ratio(seq, self.pat.reverse_complement(seq), spec), 1.0, places=12)

        bad_pwms = [
            "NOT_A_MOTIF",
            {"matrix": {"A": [1.0], "C": [1.0], "G": [1.0]}},  # missing T
            {"matrix": {"A": [1.0, 1.0], "C": [1.0], "G": [1.0], "T": [1.0]}},  # unequal lengths
            {"matrix": {"A": [-0.5], "C": [1.0], "G": [1.0], "T": [1.0]}},  # negative weight
            {"matrix": {"A": [0.0], "C": [0.0], "G": [0.0], "T": [0.0]}},  # all zero column
            {"matrix": {"A": [float("nan")], "C": [1.0], "G": [1.0], "T": [1.0]}},
        ]
        for bad in bad_pwms:
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                self.cir.load_motif_profile(bad)

    def test_generalized_circuits_presets_custom_5gene_roles_and_negative_control(self):
        p_pax = self.cir.load_motif_profile("MA0069.1")["consensus"]
        p_bcl = self.cir.load_motif_profile("MA0463.3")["consensus"]

        # 1) 2-gene mutual-toggle preset
        enh_toggle = "ACTG" * 4 + p_pax + "TGCA" * 4
        scan_tog = self.cir.scan_enhancer(enh_toggle, circuit="mutual-toggle")
        self.assertEqual(scan_tog["circuit"], "mutual-toggle")
        self.assertEqual(scan_tog["genes"], ["GeneA", "GeneB"])
        self.assertEqual(scan_tog["protein"], "MA0069.1")
        self.assertTrue(any(e["fate_flipped"] for e in scan_tog["edits"]))
        flipped = next(e for e in scan_tog["edits"] if e["fate_flipped"])
        self.assertLess(flipped["gene_deltas"]["GeneA"] * flipped["gene_deltas"]["GeneB"], 0.0)

        # 2) 3-gene tristable tri-lineage-switch preserving multi-TF edge_bindings even with target='All'
        enh_tri = "ACTG" * 3 + p_pax + "TGCA" * 2 + p_bcl + "ACGT" * 3
        scan_tri = self.cir.scan_enhancer(enh_tri, circuit="tri-lineage-switch", target="All")
        self.assertEqual(scan_tri["circuit"], "tri-lineage-switch")
        self.assertEqual(scan_tri["genes"], ["Progenitor", "Excitatory", "Inhibitory"])
        self.assertEqual(scan_tri["proteins"], ["MA0069.1", "MA0463.3"])
        self.assertEqual(scan_tri["target"], "All")

        # Perturbing MA0069.1 core changes Excitatory TF ratio while MA0463.3 stays ~1
        pax_start = 12
        var_pax = enh_tri[: pax_start + 3] + ("A" if enh_tri[pax_start + 3] != "A" else "C") + enh_tri[pax_start + 4 :]
        cmp_pax = self.cir.compare_variants(enh_tri, var_pax, circuit="tri-lineage-switch", target="All")
        self.assertLess(cmp_pax["affinity_ratios"]["MA0069.1"], 0.25)
        self.assertAlmostEqual(cmp_pax["affinity_ratios"]["MA0463.3"], 1.0, delta=0.05)

        # Repressor role inverts occupancy shift direction relative to activator
        tri_spec = self.cir.circuit_preset("tri-lineage-switch")
        state = (0.2, 0.3, 0.3)
        occ_act = self.cir.occupancy_ratios(tri_spec, state, 0.35, {"Inhibitory": 0.2}, roles={"Inhibitory": "activator"})
        occ_rep = self.cir.occupancy_ratios(tri_spec, state, 0.35, {"Inhibitory": 0.2}, roles={"Inhibitory": "repressor"})
        occ_wt = self.cir.occupancy_ratios(tri_spec, state, 0.35, {"Inhibitory": 1.0}, roles={"Inhibitory": "activator"})
        self.assertLess(occ_act[2], occ_wt[2])
        self.assertGreater(occ_rep[2], occ_wt[2])

        # 3) Custom 5-gene CircuitSpec with 5x5 cooperativity matrix
        custom_5 = {
            "name": "custom-5-gene",
            "genes": ["G0", "G1", "G2", "G3", "G4"],
            "basal": [2.0, 1.5, 1.0, 1.2, 0.8],
            "input": [100.0, 60.0, 40.0, 20.0, 10.0],
            "tf_binding": [8.0, 4.0, 2.0, 1.0, 1.0],
            "repression": [
                [0.0, 18.0, 0.0, 0.0, 5.0],
                [18.0, 0.0, 12.0, 0.0, 0.0],
                [0.0, 12.0, 0.0, 10.0, 0.0],
                [0.0, 0.0, 10.0, 0.0, 8.0],
                [5.0, 0.0, 0.0, 8.0, 0.0],
            ],
            "cooperativity": [[2.0 if i != j else 0.0 for j in range(5)] for i in range(5)],
            "production": [2.0] * 5,
            "decay": [2.0] * 5,
            "default_target": "G0",
            "default_signal": 0.3,
            "default_protein": "MA0069.1",
            "fates": {"threshold": 0.45, "labels": {"G0": "StateZero"}},
        }
        traj = self.cir.simulate_circuit([(20.0, 0.3)], circuit=custom_5, ratio=0.2, target="G0")
        self.assertTrue(all(0.0 <= v <= 1.0 for row in traj for v in row[1:]))
        self.assertEqual(self.cir.fate_label((0.8, 0.1, 0.1, 0.1, 0.1), custom_5["genes"], fates=custom_5["fates"]), "StateZero")
        with self.assertRaises(ValueError):
            self.cir.validate_circuit_spec(dict(custom_5, repression=[[0.0, 1.0]]))

        # 4) Negative control: pure alternating repeat has no consensus site -> 0 fate flips in neural-tube
        repeat_ctrl = "ACGTACGT" * 6
        ctrl_scan = self.cir.scan_enhancer(repeat_ctrl, signal=0.20, protein="GLI2")
        self.assertEqual(sum(1 for e in ctrl_scan["edits"] if e["fate_flipped"]), 0)

    def test_brainc_v3_multitf_custom_circuit_and_exhaustive_tampering(self):
        p_pax = self.cir.load_motif_profile("MA0069.1")["consensus"]
        custom_pwm = {
            "id": "CUSTOM-CORE-6",
            "kind": "log_odds",
            "matrix": {
                "A": [2.0, -3.0, -3.0, 2.0, -3.0, -3.0],
                "C": [-3.0, 2.0, -3.0, -3.0, 2.0, -3.0],
                "G": [-3.0, -3.0, 2.0, -3.0, -3.0, 2.0],
                "T": [-3.0, -3.0, -3.0, -3.0, -3.0, -3.0],
            },
        }
        ref = "TTTT" + p_pax + "ACGTACGT" + "ACGACG" + "GGGG"
        var = "TTTT" + p_pax[:4] + "AAA" + p_pax[7:] + "ACGTACGT" + "ACGTCG" + "GGGG"
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "v3bundle"
            receipt = self.cir.compile_circuit_bundle(
                ref, var, protein="MA0069.1", motif_specs=["MA0069.1", custom_pwm],
                circuit="tri-lineage-switch", output=out,
            )
            verified = self.cir.verify_circuit_bundle(out, receipt)
            self.assertEqual(set(verified["affinity_ratios"]), {"MA0069.1", "CUSTOM-CORE-6"})
            self.assertTrue((out / "circuit_spec.json").is_file())
            self.assertTrue((out / "motifs.json").is_file())

            # Exhaustive tamper checks against every sealed v3 artifact
            for target_file, mutator in (
                ("pair.fasta", lambda p: p.write_text(f">reference\n{ref}\n>variant\n{ref}\n")),
                ("circuit_spec.json", lambda p: p.write_text(p.read_text().replace("tri-lineage-switch", "tampered"))),
                ("motifs.json", lambda p: p.write_text(p.read_text().replace("CUSTOM-CORE-6", "TAMPERED"))),
                ("contract.json", lambda p: p.write_text(p.read_text().replace("0.0001", "0.0002"))),
            ):
                with self.subTest(tampered=target_file):
                    saved = (out / target_file).read_bytes()
                    try:
                        mutator(out / target_file)
                        with self.assertRaises(ValueError):
                            self.cir.verify_circuit_bundle(out, receipt)
                    finally:
                        (out / target_file).write_bytes(saved)

            # Corrupting development bundle also fails verification
            dev_path = out / "development/development_module.json"
            saved_dev = dev_path.read_bytes()
            try:
                dev_path.write_text("{}")
                with self.assertRaises(ValueError):
                    self.cir.verify_circuit_bundle(out, receipt)
            finally:
                dev_path.write_bytes(saved_dev)

    def test_bridge_to_construction_generalized_mutual_toggle_and_custom_mapping(self):
        import shutil
        if shutil.which("javac") is None or shutil.which("patch") is None:
            self.skipTest("javac and patch are required for CX3D construction bridge test")
        p_pax = self.cir.load_motif_profile("MA0069.1")["consensus"]
        ref = "ACTG" * 4 + p_pax + "TGCA" * 4
        var = "ACTG" * 4 + "A" * len(p_pax) + "TGCA" * 4
        mapping = {
            "axon_speed": {"base": 60.0, "weights": {"GeneA": 40.0, "GeneB": -10.0}, "min": 40.0, "max": 120.0},
            "guidance_gain": {"base": 1.2, "weights": {"GeneA": 1.0}, "min": 0.5, "max": 3.0},
        }
        with tempfile.TemporaryDirectory() as tmp:
            res = self.cir.bridge_to_construction(
                ref, var, Path(tmp) / "bridge_toggle",
                circuit="mutual-toggle", mapping=mapping, dev_duration=1.5,
            )
            self.assertFalse(res["biological_acceptance"])
            self.assertEqual(res["circuit"], "mutual-toggle")
            self.assertEqual(res["genes"], ["GeneA", "GeneB"])
            self.assertTrue(res["all_checks_passed"], res)
            self.assertTrue(res["checks"]["allele_rescue_matches_reference_exactly"])


if __name__ == "__main__":
    unittest.main()

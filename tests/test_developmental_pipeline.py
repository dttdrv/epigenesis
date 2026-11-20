from __future__ import annotations

import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "examples/developmental-pipeline"))

import pipeline as dp  # noqa: E402
import synthesize as syn  # noqa: E402
from brainc.validator_development import DevelopmentValidationError, validate_development_paths  # noqa: E402
from brainc.validator_variant import validate_variant_paths  # noqa: E402


class DevelopmentalPipelineTests(unittest.TestCase):
    def test_multi_locus_cis_and_coding_pipeline_fast_controls_and_exports(self) -> None:
        preset = dp.builtin_preset("neural-tube-multi-locus")
        # Add a nonsense truncation in orf_pax6 (CAA/CAG -> TAA at codon 1: CAG at index 3..6 -> TAG)
        spec = preset["spec"]
        variants = list(preset["variants"]) + [
            {
                "id": "pax6_early_stop_C4T",
                "record_id": "orf_pax6",
                "start": 3,
                "end": 4,
                "ref": "C",
                "alt": "T",
            }
        ]
        with tempfile.TemporaryDirectory() as tmp:
            troot = Path(tmp)
            fasta = troot / "loci.fasta"
            dp._write_fasta(preset["records"], fasta)
            out = troot / "run_out"

            report = dp.run_pipeline(
                fasta,
                spec,
                variants,
                out,
                run_cx3d=False,
                export_bundles=True,
            )

            self.assertEqual(report["format"], dp.FORMAT)
            self.assertFalse(report["biological_acceptance"])
            self.assertTrue(report["provenance"]["variants_valid"])
            self.assertTrue(report["provenance"]["development_valid"])
            self.assertEqual(len(report["variants"]), 2)
            for v in report["variants"]:
                self.assertTrue(v["vrs_allele_id"].startswith("ga4gh:VA."))
                self.assertTrue(v["vrs_location_id"].startswith("ga4gh:SL."))

            # Both cis enhancer and coding truncation effects should be recorded
            by_rec = {loc["record_id"]: loc for loc in report["loci"]}
            self.assertLess(by_rec["enh_nkx2_2"]["effective_affinity_ratio"], 0.1)
            self.assertEqual(by_rec["orf_pax6"]["consequence"], "nonsense_truncation")
            self.assertEqual(by_rec["orf_pax6"]["coding_activity_factor"], 0.05)

            # Built-in causal & negative controls must hold exactly
            self.assertTrue(report["controls"]["reversion_exact"])
            self.assertEqual(report["controls"]["reversion_max_state_diff"], 0.0)
            self.assertTrue(report["controls"]["closed_chromatin_cis_exact"])
            self.assertEqual(report["controls"]["closed_chromatin_cis_max_state_diff"], 0.0)
            self.assertTrue(report["controls"]["passive_reporter_isolated"])
            self.assertEqual(report["controls"]["passive_reporter_isolation_max_diff"], 0.0)

            # Parameter-uncertainty ensemble sanity
            ens = report["uncertainty_ensemble"]
            self.assertEqual(ens["samples"], dp.DEFAULT_ENSEMBLE_SAMPLES)
            self.assertGreaterEqual(ens["fate_flip_probability"], 0.0)
            self.assertLessEqual(ens["fate_flip_probability"], 1.0)
            self.assertAlmostEqual(sum(ens["wt_fate_distribution"].values()), 1.0, places=9)
            self.assertAlmostEqual(sum(ens["variant_fate_distribution"].values()), 1.0, places=9)

            # Portable exports and HTML report exist
            self.assertTrue((out / "pipeline_report.json").is_file())
            self.assertTrue((out / "index.html").is_file())
            self.assertTrue((out / "rocrate" / "ro-crate-metadata.json").is_file())
            self.assertTrue((out / "bundle.omex").is_file())

            # Independent replay checks pass
            v_check = validate_variant_paths(
                out / "source",
                out / "variants" / "variants.json",
                source_inputs={"sequence": out / "reference.fasta"},
            )
            d_check = validate_development_paths(
                out / "source",
                out / "development",
                source_inputs={"sequence": out / "reference.fasta"},
            )
            self.assertTrue(v_check["valid"])
            self.assertTrue(d_check["valid"])

            # Rust cross-check passes whenever the compiled Rust validator binary is present
            rc_check = report["provenance"]["rust_cross_check"]
            if rc_check["available"]:
                self.assertTrue(rc_check["variants_valid"])
                self.assertTrue(rc_check["development_valid"])

            # Tampering with reference.fasta breaks both variant and development replay
            tampered = out / "reference.fasta"
            orig_text = tampered.read_text(encoding="utf-8")
            tampered.write_text(orig_text.replace("GACCACCCA", "GACCACCCT", 1), encoding="utf-8")
            self.assertFalse(
                validate_variant_paths(
                    out / "source",
                    out / "variants" / "variants.json",
                    source_inputs={"sequence": tampered},
                )["valid"]
            )
            with self.assertRaises(DevelopmentValidationError):
                validate_development_paths(
                    out / "source",
                    out / "development",
                    source_inputs={"sequence": tampered},
                )

    def test_chromatin_gating_monotonicity_and_tri_lineage_preset(self) -> None:
        preset = dp.builtin_preset("tri-lineage-chromatin")
        with tempfile.TemporaryDirectory() as tmp:
            troot = Path(tmp)
            fasta = troot / "tri.fasta"
            dp._write_fasta(preset["records"], fasta)

            closed_spec = json.loads(json.dumps(preset["spec"]))
            closed_spec["loci"][0]["accessibility"] = 0.0
            closed_spec["uncertainty"]["samples"] = 6
            r_closed = dp.run_pipeline(
                fasta, closed_spec, preset["variants"], troot / "closed", run_cx3d=False
            )

            open_spec = json.loads(json.dumps(preset["spec"]))
            open_spec["loci"][0]["accessibility"] = 1.0
            open_spec["uncertainty"]["samples"] = 6
            r_open = dp.run_pipeline(
                fasta, open_spec, preset["variants"], troot / "open", run_cx3d=False
            )

            self.assertEqual(r_closed["loci"][0]["effective_affinity_ratio"], 1.0)
            self.assertNotEqual(r_open["loci"][0]["effective_affinity_ratio"], 1.0)
            self.assertEqual(r_closed["nominal"]["gene_deltas"]["Excitatory"], 0.0)
            self.assertNotEqual(r_open["nominal"]["gene_deltas"]["Excitatory"], 0.0)

    def test_empirical_falsification_accepts_consistent_and_rejects_contradictory(self) -> None:
        preset = dp.builtin_preset("neural-tube-multi-locus")
        spec = json.loads(json.dumps(preset["spec"]))
        spec["uncertainty"]["samples"] = 4
        with tempfile.TemporaryDirectory() as tmp:
            troot = Path(tmp)
            fasta = troot / "loci.fasta"
            dp._write_fasta(preset["records"], fasta)

            base = dp.run_pipeline(fasta, spec, preset["variants"], troot / "base", run_cx3d=False)
            pred = base["nominal"]["gene_deltas"]

            matching_obs = [
                {"gene": "Nkx2.2", "observed_delta": pred["Nkx2.2"]},
                {"gene": "Olig2", "observed_delta": pred["Olig2"]},
            ]
            contradicting_obs = [
                {"gene": "Nkx2.2", "observed_delta": -pred["Nkx2.2"] if pred["Nkx2.2"] != 0 else 0.9},
                {"gene": "Olig2", "observed_delta": -pred["Olig2"] if pred["Olig2"] != 0 else -0.9},
            ]

            r_ok = dp.run_pipeline(
                fasta, spec, preset["variants"], troot / "ok", run_cx3d=False, observations=matching_obs
            )
            r_bad = dp.run_pipeline(
                fasta, spec, preset["variants"], troot / "bad", run_cx3d=False, observations=contradicting_obs
            )
            self.assertFalse(r_ok["falsification"]["falsified"])
            self.assertAlmostEqual(r_ok["falsification"]["rmse"], 0.0, places=9)
            self.assertTrue(r_bad["falsification"]["falsified"])
            self.assertGreater(r_bad["falsification"]["rmse"], 0.25)

    def test_overlapping_variants_rejected(self) -> None:
        preset = dp.builtin_preset("neural-tube-multi-locus")
        bad_variants = [
            {"id": "v1", "record_id": "enh_nkx2_2", "start": 10, "end": 13, "ref": "GAC", "alt": "TTT"},
            {"id": "v2", "record_id": "enh_nkx2_2", "start": 12, "end": 14, "ref": "CC", "alt": "AA"},
        ]
        with tempfile.TemporaryDirectory() as tmp:
            troot = Path(tmp)
            fasta = troot / "loci.fasta"
            dp._write_fasta(preset["records"], fasta)
            with self.assertRaisesRegex(ValueError, "overlapping"):
                dp.run_pipeline(fasta, preset["spec"], bad_variants, troot / "fail", run_cx3d=False)

    @unittest.skipUnless(shutil.which("javac") and shutil.which("java"), "JDK required for CX3D morphogenesis")
    def test_end_to_end_with_3d_cx3d_morphogenesis_and_lif_assay(self) -> None:
        preset = dp.builtin_preset("neural-tube-multi-locus")
        spec = json.loads(json.dumps(preset["spec"]))
        spec["uncertainty"]["samples"] = 4
        with tempfile.TemporaryDirectory() as tmp:
            troot = Path(tmp)
            fasta = troot / "loci.fasta"
            dp._write_fasta(preset["records"], fasta)
            out = troot / "cx3d_run"
            report = dp.run_pipeline(
                fasta,
                spec,
                preset["variants"],
                out,
                run_cx3d=True,
                cx3d_duration=1.5,
                cx3d_seed=7,
            )
            morph = report["morphogenesis"]
            self.assertTrue(morph["executed"])
            self.assertTrue(morph["allele_rescue_matches_reference"])
            self.assertTrue(morph["morphology_diverged"])
            for branch in ("reference", "variant", "allele_rescue"):
                self.assertTrue(morph["branches"][branch]["compiled_valid"])
                self.assertTrue(morph["branches"][branch]["activity"]["transmission_controls_pass"])

    def test_inverse_compiler_epistatic_rescue_and_reprogram(self) -> None:
        seq18 = syn._load_gse_locus("chr18:56475717-56475888")
        with tempfile.TemporaryDirectory() as tmp:
            troot = Path(tmp)
            rep_reprog = syn.synthesize_reprogram(
                seq18,
                "Olig2",
                troot / "reprog",
                record_id="chr18",
                ensemble_samples=6,
                run_cx3d=False,
                export_bundles=False,
            )
            self.assertEqual(rep_reprog["synthesized_edit"]["nominal_variant_fate"], "Olig2")
            self.assertTrue(rep_reprog["synthesized_edit"]["vrs_allele_id"].startswith("ga4gh:VA."))
            self.assertTrue(rep_reprog["controls"]["reversion_exact"])
            self.assertFalse(rep_reprog["biological_acceptance"])

            rep_rescue = syn.synthesize_epistatic_rescue(
                seq18,
                "165:G>T",
                troot / "rescue",
                record_id="chr18",
                ensemble_samples=6,
                run_cx3d=bool(shutil.which("javac") and shutil.which("java")),
                cx3d_duration=1.2,
                export_bundles=True,
            )
            esc = rep_rescue["synthesized_epistatic_rescue"]
            self.assertTrue(esc["primary_kept_intact"])
            self.assertEqual(esc["restored_nominal_fate"], "Nkx2.2")
            self.assertGreaterEqual(esc["ensemble_wt_fate_concordance"], 0.99)
            self.assertLess(esc["rescue_vs_wt_max_gene_delta"], 1e-5)
            self.assertEqual(len(esc["vrs_allele_ids"]), 2)
            if rep_rescue["cx3d_3d_restoration"]["executed"]:
                self.assertGreater(rep_rescue["cx3d_3d_restoration"]["neurite_morphology_recovery_pct"], 99.9)
                self.assertTrue(rep_rescue["cx3d_3d_restoration"]["all_lif_transmission_controls_pass"])


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

from brainc.validator_development import validate_development_paths
from brainc.validator_variant import validate_variant_paths


ROOT = Path(__file__).resolve().parents[1]
MPRA_PATH = ROOT / "examples" / "mpra-replay" / "mpra_replay.py"
RUST_DIR = ROOT / "validators" / "rust"


def _load_mpra():
    spec = importlib.util.spec_from_file_location("mpra_replay", MPRA_PATH)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["mpra_replay"] = mod
    spec.loader.exec_module(mod)
    return mod


class MPRAReplayTests(unittest.TestCase):
    def test_dual_parameter_gse188264_replay_and_flips(self) -> None:
        mpra = _load_mpra()
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "mpra_out"
            rep = mpra.run_mpra_replay(out, split="test", limit=150, export_bundles=True)
            self.assertEqual(rep["format"], "epigenesis.mpra-replay-comparison")
            self.assertEqual(rep["dataset"]["accession"], "GSE188264")
            self.assertEqual(rep["dataset"]["variants_evaluated"], 150)
            self.assertGreater(rep["summary"]["total_flipped_calls"], 0)
            self.assertTrue(rep["sealed_provenance"]["parameter_set_a"]["valid"])
            self.assertTrue(rep["sealed_provenance"]["parameter_set_b"]["valid"])

            # Independent Python replay of variants + both parameter bundles
            src_dir = out / "source"
            vc_path = out / "variants.json"
            fasta_path = out / "wt_loci.fasta"
            v_rep = validate_variant_paths(src_dir, vc_path)
            self.assertTrue(v_rep["valid"])
            self.assertEqual(v_rep["summary"]["variant_count"], 150)

            for sub in ("param_a_permissive_mean_fc", "param_b_strict_late_qc"):
                d_rep = validate_development_paths(
                    src_dir,
                    out / sub / "development",
                    source_inputs={"sequence": fasta_path},
                )
                self.assertTrue(d_rep["valid"], sub)

            self.assertTrue((out / "exports" / "param_a_rocrate" / "ro-crate-metadata.json").is_file())
            self.assertTrue((out / "exports" / "param_b_strict_late_qc.omex").is_file())

            # Independent Rust replay if cargo is available
            if shutil.which("cargo") and (RUST_DIR / "target" / "debug" / "epigenesis-validate-rs").is_file():
                rs_bin = RUST_DIR / "target" / "debug" / "epigenesis-validate-rs"
                rs_var = subprocess.run(
                    [str(rs_bin), "validate-variants", str(src_dir), str(vc_path)],
                    check=False,
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(rs_var.returncode, 0, rs_var.stderr + rs_var.stdout)
                self.assertTrue(json.loads(rs_var.stdout)["valid"])

                rs_dev = subprocess.run(
                    [
                        str(rs_bin),
                        "validate-development",
                        str(src_dir),
                        str(out / "param_b_strict_late_qc" / "development"),
                    ],
                    check=False,
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(rs_dev.returncode, 0, rs_dev.stderr + rs_dev.stdout)
                self.assertEqual(
                    json.loads(rs_dev.stdout)["development_bundle_artifact_sha256"],
                    rep["sealed_provenance"]["parameter_set_b"]["development_bundle_sha256"],
                )


if __name__ == "__main__":
    unittest.main()

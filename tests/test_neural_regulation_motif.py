from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

try:
    import numpy as np  # noqa: F401
except ImportError as failure:
    raise unittest.SkipTest(f"NumPy is required for motif-temporal neural regulation tests: {failure}")


ROOT = Path(__file__).resolve().parents[1]
REG_PATH = ROOT / "examples/neural-regulation/regulation.py"
MOD_PATH = ROOT / "examples/neural-regulation/motif_temporal.py"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


class NeuralRegulationMotifTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reg = _load("neural_regulation_base", REG_PATH)
        cls.mt = _load("neural_regulation_motif_temporal", MOD_PATH)

    def _synthetic_rows(self):
        # Construct families separated by real JASPAR GLI motif MA1491.3 (GGACCACCCACGTCG)
        ref = "A" * 80 + "GGACCACCCACGTCG" + "C" * 76
        mut = "A" * 80 + "GGAACACCCACGTCG" + "C" * 76
        rows = []
        for idx, (fam, split) in enumerate((("f1", "train"), ("f2", "train"), ("f3", "train"), ("f4", "test"))):
            rows.append(self.reg.Row(f"{fam}_loss", fam, fam, split, "targeted", ref, mut, (0.2,) * 7, ((1.0,) * 7,)))
            # Null / tiny distant flanking edit
            alt2 = "G" + ref[1:]
            rows.append(self.reg.Row(f"{fam}_flank", fam, fam, split, "targeted", ref, alt2, (1.0 + 0.01 * idx,) * 7, ((1.0,) * 7,)))
        return rows

    def test_fit_predict_antisymmetry_and_family_invariance(self):
        rows = [r for r in self._synthetic_rows() if r.split == "train"]
        model, coef = self.mt.fit(rows, 1.0)
        ref, alt = rows[0].reference, rows[0].sequence
        pred_fwd = self.mt.predict_rows(model, coef, [{"reference": ref, "sequence": alt}])[0]
        pred_rev = self.mt.predict_rows(model, coef, [{"reference": alt, "sequence": ref}])[0]
        pred_same = self.mt.predict_rows(model, coef, [{"reference": ref, "sequence": ref}])[0]
        self.assertLess(np.max(np.abs(pred_fwd + pred_rev)), 1e-12)
        self.assertLess(np.max(np.abs(pred_same)), 1e-12)
        # Stronger loss on GLI perturbation than zero prediction
        expected = self.reg.target(rows[0], 1.0)[0]
        self.assertLess(abs(pred_fwd[0] - expected), abs(expected))

        # Duplicating one family must leave weighted solution invariant
        dup_rows = rows + [r._replace(id=r.id + "_dup") for r in rows[:2]]
        _, coef_dup = self.mt.fit(dup_rows, 1.0)
        self.assertLess(np.max(np.abs(coef - coef_dup)), 1e-12)

    def test_shuffle_control_degrades_fit(self):
        rows = [r for r in self._synthetic_rows() if r.split == "train"]
        m_real, c_real = self.mt.fit(rows, 1.0)
        m_shuf, c_shuf = self.mt.fit(rows, 1.0, shuffle_seed=19)
        pairs = [{"reference": r.reference, "sequence": r.sequence} for r in rows]
        y = np.array([self.reg.target(r, 1.0) for r in rows])
        mse_real = float(np.mean((y - self.mt.predict_rows(m_real, c_real, pairs)) ** 2))
        mse_shuf = float(np.mean((y - self.mt.predict_rows(m_shuf, c_shuf, pairs)) ** 2))
        self.assertLess(mse_real, mse_shuf)

    def test_compile_prediction_roundtrip_and_tamper_rejection(self):
        rows = [r for r in self._synthetic_rows() if r.split == "train"]
        model, coef = self.mt.fit(rows, 1.0)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mdir = root / "model_dir"
            mdir.mkdir()
            np.save(mdir / "coefficients.npy", coef, allow_pickle=False)
            import hashlib
            model["coefficients_sha256"] = hashlib.sha256((mdir / "coefficients.npy").read_bytes()).hexdigest()
            (mdir / "model.json").write_text(json.dumps(model))
            ref_p, var_p = root / "ref.dna", root / "var.dna"
            ref_p.write_text(rows[0].reference)
            var_p.write_text(rows[0].sequence)
            out = root / "compiled"
            compiled = self.mt.compile_prediction(mdir, ref_p, var_p, out)
            reloaded = self.mt.load_prediction(out)
            self.assertEqual(compiled, reloaded)
            # Tampering with pair.fasta triggers mismatch
            (out / "pair.fasta").write_text(f">reference\n{rows[0].reference}\n>variant\n{rows[0].reference}\n")
            with self.assertRaises(ValueError):
                self.mt.load_prediction(out)


if __name__ == "__main__":
    unittest.main()

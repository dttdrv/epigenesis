from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest

try:
    import numpy as np
    import scipy  # noqa: F401
except ImportError as failure:
    raise unittest.SkipTest(f"NumPy and SciPy are required for the thermodynamic transfer tests: {failure}")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import sites  # noqa: E402
import thermo  # noqa: E402


class SiteSumTests(unittest.TestCase):
    def test_matches_brute_force_and_strand_symmetry(self):
        rng = np.random.default_rng(3)
        matrices = [rng.normal(size=(4, 5)), rng.normal(size=(4, 7))]
        seqs = ["".join(rng.choice(list("ACGT"), size=40)) for _ in range(3)]
        kappa = 0.05
        got = sites.site_sums(seqs, matrices, [kappa])[0]
        comp = str.maketrans("ACGT", "TGCA")
        for i, seq in enumerate(seqs):
            for j, m in enumerate(matrices):
                width, total = m.shape[1], 0.0
                for start in range(len(seq) - width + 1):
                    word = seq[start:start + width]
                    fwd = sum(m["ACGT".index(c), p] for p, c in enumerate(word))
                    rc = word.translate(comp)[::-1]
                    rev = sum(m["ACGT".index(c), p] for p, c in enumerate(rc))
                    odds = (np.exp(fwd) + np.exp(rev)) / 2
                    total += kappa * odds / (1 + kappa * odds)
                self.assertAlmostEqual(got[i, j], total, places=10)
        rc_seqs = [s.translate(comp)[::-1] for s in seqs]
        np.testing.assert_allclose(sites.site_sums(rc_seqs, matrices, [kappa])[0], got, rtol=1e-10)

    def test_shuffle_keeps_column_composition_and_changes_order(self):
        rng = np.random.default_rng(5)
        matrices = [rng.normal(size=(4, 9)) for _ in range(4)]
        shuffled = sites.shuffle_columns(matrices, 7)
        for before, after in zip(matrices, shuffled):
            self.assertEqual(sorted(map(tuple, before.T.round(12))), sorted(map(tuple, after.T.round(12))))
        self.assertTrue(any(not np.array_equal(b, a) for b, a in zip(matrices, shuffled)))
        again = sites.shuffle_columns(matrices, 7)
        for a, b in zip(shuffled, again):
            np.testing.assert_array_equal(a, b)


class ModelTests(unittest.TestCase):
    def test_analytic_gradient(self):
        rng = np.random.default_rng(0)
        m, n = 10, 25
        var, par = rng.random((n, m)), rng.random((n, m))
        y, w = rng.normal(size=(n, thermo.T)), np.full(n, 1 / n)
        x = rng.random(2 * m * thermo.T) * 0.5
        smooth = thermo._smooth_matrix(10.0)
        _, grad = thermo.objective(x, var, par, y, w, 0.01, smooth, m)
        numeric = np.zeros_like(x)
        for i in range(0, x.size, 7):
            e = np.zeros_like(x)
            e[i] = 1e-6
            numeric[i] = (thermo.objective(x + e, var, par, y, w, 0.01, smooth, m)[0]
                          - thermo.objective(x - e, var, par, y, w, 0.01, smooth, m)[0]) / 2e-6
            self.assertAlmostEqual(numeric[i], grad[i], places=7)

    def test_recovers_planted_nonnegative_strengths(self):
        rng = np.random.default_rng(1)
        m, n = 8, 600
        var, par = rng.random((n, m)) * 2, rng.random((n, m)) * 2
        a_true = np.zeros((m, thermo.T))
        a_true[2] = np.linspace(0.2, 1.5, thermo.T)
        r_true = np.zeros((m, thermo.T))
        r_true[5] = np.linspace(1.0, 0.1, thermo.T)
        y = thermo.predict_effects(a_true, r_true, var, par)
        fit = thermo.fit_thermo(var, par, y, np.full(n, 1 / n), 1e-8, 0.0)
        self.assertTrue((fit["A"] >= 0).all() and (fit["R"] >= 0).all())
        pred = thermo.predict_effects(fit["A"], fit["R"], var, par)
        self.assertLess(float(((pred - y) ** 2).mean()), 1e-3 * float((y ** 2).mean()))

    def test_bootstrap_interval_covers_known_gain(self):
        base = {f"g{i}": 1.0 for i in range(40)}
        cand = {f"g{i}": 0.9 + 0.01 * (i % 3) for i in range(40)}
        out = thermo.bootstrap_gain(base, cand)
        self.assertLess(out["lower"], out["gain"])
        self.assertLess(out["gain"], out["upper"])
        self.assertGreater(out["lower"], 0.0)


class FreezeTests(unittest.TestCase):
    def test_frozen_files_unchanged(self):
        freeze = json.loads((HERE / "FREEZE.json").read_text())
        for name, key in (("PLAN.md", "plan_sha256"), ("thermo.py", "thermo_sha256"), ("sites.py", "sites_sha256"),
                          ("confirm.py", "confirm_sha256"), ("selection.json", "selection_sha256")):
            self.assertEqual(thermo.sha256_file(HERE / name), freeze[key], name)


if __name__ == "__main__":
    unittest.main()

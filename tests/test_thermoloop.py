"""CPU-only scientific invariant and contract tests for Epigenesis ThermoLoop / ThermoInterlace."""

from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "studio"))

try:
    import numpy as np
    import torch
    from epigenesis_studio.thermoloop import (
        ThermoLoopModel,
        ThermoLoopV2Model,
        ThermoLoopV3Model,
        encode_sequences_onehot,
        load_verified_jaspar_bank,
        misul_rational_gate,
        reverse_complement,
    )
    HAS_TORCH = True
except ImportError:
    HAS_TORCH = False


@unittest.skipUnless(HAS_TORCH, "PyTorch not installed in stdlib test interpreter")
class TestThermoLoopInvariants(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        torch.manual_seed(1001)
        cls.wf, cls.wr, cls.mids, cls.sha = load_verified_jaspar_bank(21)
        # Compact models for fast CPU verification
        cls.loop_model = ThermoLoopModel(
            cls.wf, cls.wr, cls.mids,
            denovo_ch=16, bio_spatial_ch=16, core_dim=32, n_loops=4, n_heads=4,
            use_biophys=True, use_rational=True, use_delta_mem=False,
        )
        cls.interlace_model = ThermoLoopModel(
            cls.wf, cls.wr, cls.mids,
            denovo_ch=16, bio_spatial_ch=16, core_dim=32, n_loops=4, n_heads=4,
            use_biophys=True, use_rational=True, use_delta_mem=True,
        )
        cls.v2_loop_model = ThermoLoopV2Model(
            cls.wf, cls.wr, cls.mids,
            denovo_ch=16, bio_spatial_ch=16, core_dim=48, n_loops=4, n_heads=4,
            use_biophys=True, use_rational=True, use_delta_mem=False,
        )
        cls.v2_interlace_model = ThermoLoopV2Model(
            cls.wf, cls.wr, cls.mids,
            denovo_ch=16, bio_spatial_ch=16, core_dim=48, n_loops=4, n_heads=4,
            use_biophys=True, use_rational=True, use_delta_mem=True,
        )
        cls.v3_loop_model = ThermoLoopV3Model(
            cls.wf, cls.wr, cls.mids,
            denovo_ch=16, bio_spatial_ch=16, core_dim=48, n_loops=4, n_heads=4,
            use_biophys=True, use_rational=True, use_delta_mem=False,
        )
        cls.v3_interlace_model = ThermoLoopV3Model(
            cls.wf, cls.wr, cls.mids,
            denovo_ch=16, bio_spatial_ch=16, core_dim=48, n_loops=4, n_heads=4,
            use_biophys=True, use_rational=True, use_delta_mem=True,
        )
        for m in (
            cls.loop_model,
            cls.interlace_model,
            cls.v2_loop_model,
            cls.v2_interlace_model,
            cls.v3_loop_model,
            cls.v3_interlace_model,
        ):
            m.eval()
        # Deterministic 230-bp test sequence
        rng = np.random.default_rng(42)
        bases = np.array(list("ACGT"))
        cls.wt = "".join(bases[rng.integers(0, 4, size=230)])
        cls.mut_a = cls.wt[:80] + ("A" if cls.wt[80] != "A" else "G") + cls.wt[81:]
        cls.mut_b = cls.wt[:140] + ("C" if cls.wt[140] != "C" else "T") + cls.wt[141:]
        cls.mut_ab = (
            cls.wt[:80] + cls.mut_a[80] + cls.wt[81:140] + cls.mut_b[140] + cls.wt[141:]
        )

    def test_jaspar_bank_verified_and_nonpositive_energy(self) -> None:
        self.assertEqual(self.wf.shape, (532, 4, 21))
        self.assertEqual(len(self.mids), 532)
        self.assertEqual(len(self.sha), 64)
        self.assertLessEqual(float(self.wf.max()), 1e-6)
        self.assertLessEqual(float(self.wr.max()), 1e-6)

    def test_misul_rational_gate_strictly_bounded_and_finite(self) -> None:
        # Must stay strictly bounded within [-1/sqrt(2), +1/sqrt(2)] even for extreme inputs +-1e6
        vals = torch.tensor([-1e6, -1e3, -10.0, -1.0, 0.0, 1.0, 10.0, 1e3, 1e6], requires_grad=True)
        a = vals.unsqueeze(0).expand(9, 9)
        b = vals.unsqueeze(1).expand(9, 9)
        y = misul_rational_gate(a, b)
        self.assertTrue(torch.isfinite(y).all())
        bound = 1.0 / math.sqrt(2.0) + 1e-5
        self.assertLessEqual(float(y.detach().abs().max()), bound)
        y.sum().backward()
        self.assertTrue(torch.isfinite(vals.grad).all())

    def test_exact_variant_symmetries_thermoloop_and_interlace(self) -> None:
        for name, model in (
            ("ThermoLoop", self.loop_model),
            ("ThermoInterlace", self.interlace_model),
            ("ThermoLoopV2", self.v2_loop_model),
            ("ThermoInterlaceV2", self.v2_interlace_model),
            ("ThermoLoopV3", self.v3_loop_model),
            ("ThermoInterlaceV3", self.v3_interlace_model),
        ):
            with self.subTest(model=name):
                # 1. Exact reversion identity: delta(WT, WT) == 0.0
                rev = model.predict_variant(self.wt, self.wt)
                self.assertEqual(rev["delta_effect"], 0.0)

                # 2. Exact antisymmetry: delta(A, WT) == -delta(WT, A)
                fwd = model.predict_variant(self.wt, self.mut_a)
                bwd = model.predict_variant(self.mut_a, self.wt)
                self.assertAlmostEqual(fwd["delta_effect"] + bwd["delta_effect"], 0.0, places=6)

                # 3. Exact reverse-complement invariance: delta(RC(WT), RC(A)) == delta(WT, A)
                rc_res = model.predict_variant(reverse_complement(self.wt), reverse_complement(self.mut_a))
                self.assertAlmostEqual(fwd["delta_effect"], rc_res["delta_effect"], places=6)

                # 4. Calibrated uncertainty & equilibrium residual must be finite & positive
                self.assertGreater(fwd["uncertainty_sigma"], 0.0)
                self.assertTrue(math.isfinite(fwd["equilibrium_residual"]))
                self.assertGreaterEqual(fwd["equilibrium_residual"], 0.0)

    def test_non_additive_epistasis_and_loop_supervision_gradients(self) -> None:
        for name, model in (
            ("ThermoLoop", self.loop_model),
            ("ThermoInterlace", self.interlace_model),
            ("ThermoLoopV2", self.v2_loop_model),
            ("ThermoInterlaceV2", self.v2_interlace_model),
            ("ThermoLoopV3", self.v3_loop_model),
            ("ThermoInterlaceV3", self.v3_interlace_model),
        ):
            with self.subTest(model=name):
                epi = model.decompose_epistasis(self.wt, self.mut_a, self.mut_b, self.mut_ab)
                self.assertAlmostEqual(
                    epi["delta_ab"],
                    epi["additive_expectation"] + epi["non_additive_epistasis"],
                    places=6,
                )
                # Verify backward pass across all loop steps produces finite gradients
                model.train()
                oh = encode_sequences_onehot([self.wt, self.mut_a], torch.device("cpu"))
                out = model.forward_loops(oh, return_all_steps=True)
                self.assertGreaterEqual(len(out["step_mus"]), 2)
                loss = sum(m.pow(2).mean() for m in out["step_mus"]) + out["logvar"].pow(2).mean()
                model.zero_grad(set_to_none=True)
                loss.backward()
                for pname, p in model.named_parameters():
                    if p.requires_grad and p.grad is not None:
                        self.assertTrue(torch.isfinite(p.grad).all(), f"non-finite grad in {pname}")
                model.eval()


if __name__ == "__main__":
    unittest.main()

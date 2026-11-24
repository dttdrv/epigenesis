"""Unit tests for Epigenesis HelixLoop-Omega (Tier 1 & Tier 2) and HybridMuonCWD."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "studio"))

try:
    import numpy as np
    import torch
    import torch.nn as nn

    HAS_TORCH = True
except ImportError:
    HAS_TORCH = False


@unittest.skipUnless(HAS_TORCH, "PyTorch not installed in this Python environment")
class HelixLoopOmegaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        from epigenesis_studio.helixloop import (
            BiGatedDeltaProduct1d,
            CanonicalDNAEngram,
            Derf,
            FlashzoiHelixLoopOmegaModel,
            ForgettingDiffGatedAttention,
            HelixLoopOmegaModel,
            HybridMuonCWD,
            PairwiseTriangleEpistasis2d,
            RCPSDNAStem,
            RMSBatchNorm1d,
            StandardizedConv1d,
            encode_sequences_onehot,
            load_verified_jaspar_bank,
            reverse_complement,
        )

        cls.hl = dict(
            BiGatedDeltaProduct1d=BiGatedDeltaProduct1d,
            CanonicalDNAEngram=CanonicalDNAEngram,
            Derf=Derf,
            FlashzoiHelixLoopOmegaModel=FlashzoiHelixLoopOmegaModel,
            ForgettingDiffGatedAttention=ForgettingDiffGatedAttention,
            HelixLoopOmegaModel=HelixLoopOmegaModel,
            HybridMuonCWD=HybridMuonCWD,
            PairwiseTriangleEpistasis2d=PairwiseTriangleEpistasis2d,
            RCPSDNAStem=RCPSDNAStem,
            RMSBatchNorm1d=RMSBatchNorm1d,
            StandardizedConv1d=StandardizedConv1d,
            encode_sequences_onehot=encode_sequences_onehot,
            reverse_complement=reverse_complement,
        )
        cls.w_fwd, cls.w_rc, cls.motif_ids, _ = load_verified_jaspar_bank(21)

    def test_canonical_engram_exact_rc_equivariance(self) -> None:
        torch.manual_seed(0)
        engram = self.hl["CanonicalDNAEngram"](embed_dim=16).eval()
        s = "AGGACCGGATCAACTGCGCATGCGATCATGCTAGCATGCTAGCATGC"
        s_rc = self.hl["reverse_complement"](s)
        x = self.hl["encode_sequences_onehot"]([s, s_rc], torch.device("cpu"))
        y = engram(x)
        # First 16 channels (6-mer embeddings) + 4 symmetric shape channels must be EXACTLY RC-equivariant!
        emb_fwd = y[0, :20]
        emb_rc_flipped = torch.flip(y[1, :20], dims=[-1])
        self.assertLess((emb_fwd - emb_rc_flipped).abs().max().item(), 1e-6)
        # Last 2 channels (AT skew, GC skew) are exact antisymmetric under RC
        skew_fwd = y[0, 20:]
        skew_rc_flipped = torch.flip(y[1, 20:], dims=[-1])
        self.assertLess((skew_fwd + skew_rc_flipped).abs().max().item(), 1e-6)

    def test_deltaproduct_wy_matches_exact_sequential_householder_recurrence(self) -> None:
        torch.manual_seed(99)
        dp = self.hl["BiGatedDeltaProduct1d"](dim=32, n_heads=4, key_dim=8).eval()
        x = torch.randn(2, 12, 32)
        out_wy = dp._scan_dir(x)

        # Reference sequential n_h=2 Householder recurrence
        b, l, _ = x.shape
        h, kd, vd = dp.n_heads, dp.key_dim, dp.val_dim
        q = torch.nn.functional.linear(x, dp.q_proj.weight).view(b, l, h, kd)
        q = q * torch.rsqrt(q.pow(2).sum(dim=-1, keepdim=True) + 1e-5)
        ks = torch.nn.functional.linear(x, dp.k_proj.weight).view(b, l, 2, h, kd)
        ks = ks * torch.rsqrt(ks.pow(2).sum(dim=-1, keepdim=True) + 1e-5)
        vs = torch.nn.functional.linear(x, dp.v_proj.weight).view(b, l, 2, h, vd)
        alpha = torch.exp(-0.05 * torch.nn.functional.softplus(dp.decay_proj(x).clamp(-10.0, 10.0)))
        betas = 2.0 * torch.sigmoid(dp.beta_proj(x)).view(b, l, 2, h)
        S = torch.zeros(b, h, kd, vd)
        outs = []
        for t in range(l):
            S = S * alpha[:, t, :, None, None]
            for r in range(2):
                kr, vr, br = ks[:, t, r], vs[:, t, r], betas[:, t, r, :, None]
                pred = torch.einsum("bhk,bhkv->bhv", kr, S)
                S = S + torch.einsum("bhk,bhv->bhkv", kr, (vr - pred) * br)
            outs.append(torch.einsum("bhk,bhkv->bhv", q[:, t], S))
        out_seq = torch.stack(outs, dim=1)
        out_seq = (out_seq * torch.rsqrt(out_seq.pow(2).mean(dim=-1, keepdim=True) + 1e-5)).reshape(b, l, dp.dim)
        max_diff = (out_wy - out_seq).abs().max().item()
        self.assertLess(max_diff, 1e-5, f"WY triangular solve diverged from exact sequential Householder: {max_diff}")

    def test_helixloop_omega_single_and_multicell_invariants(self) -> None:
        torch.manual_seed(7)
        for out_dim in (1, 3):
            model = self.hl["HelixLoopOmegaModel"](
                self.w_fwd,
                self.w_rc,
                self.motif_ids,
                out_dim=out_dim,
                core_dim=128,
                n_loops=4,
                n_heads=4,
            ).eval()
            rng = np.random.default_rng(11 + out_dim)
            bases = np.array(list("ACGT"))
            ref = "".join(rng.choice(bases, size=230))
            alt = ref[:115] + ("C" if ref[115] != "C" else "G") + ref[116:]

            # Reversion identity must be bit-for-bit 0.0
            r_same = model.predict_variant(ref, ref)
            self.assertEqual(r_same["delta_effect"], 0.0)
            self.assertEqual(len(r_same["delta_per_head"]), out_dim)

            # Antisymmetry
            r_fa = model.predict_variant(ref, alt)
            r_af = model.predict_variant(alt, ref)
            self.assertAlmostEqual(r_fa["delta_effect"] + r_af["delta_effect"], 0.0, places=6)

            # Reverse-complement invariance
            rc = self.hl["reverse_complement"]
            r_rc = model.predict_variant(rc(ref), rc(alt))
            self.assertAlmostEqual(r_fa["delta_effect"], r_rc["delta_effect"], places=6)

            # Epistasis decomposition runs cleanly
            alt_b = ref[:140] + ("A" if ref[140] != "A" else "T") + ref[141:]
            alt_ab = alt[:140] + alt_b[140] + alt[141:]
            epi = model.decompose_epistasis(ref, alt, alt_b, alt_ab)
            self.assertIn("non_additive_epistasis", epi)
            self.assertTrue(np.isfinite(epi["non_additive_epistasis"]))

    def test_helixloop_omega_training_step_with_hybrid_muon_cwd(self) -> None:
        torch.manual_seed(13)
        model = self.hl["HelixLoopOmegaModel"](
            self.w_fwd,
            self.w_rc,
            self.motif_ids,
            out_dim=3,
            core_dim=128,
            n_loops=4,
            n_heads=4,
        )
        opt = self.hl["HybridMuonCWD"](model, lr_muon=0.005, lr_adam=0.003, weight_decay=0.05)
        self.assertGreater(len(opt.param_groups[0]["params"]), 0, "Expected 2D hidden matrices in Muon group")

        x = torch.zeros(4, 4, 230)
        idx = torch.randint(0, 4, (4, 230))
        x.scatter_(1, idx.unsqueeze(1), 1.0)
        ct = torch.tensor([0, 1, 2, 1], dtype=torch.long)
        y = torch.randn(4)

        out = model.forward_loops(x, return_all_steps=True, cell_idx=ct)
        self.assertEqual(out["mu"].shape, (4, 3))
        self.assertGreaterEqual(len(out["step_mus"]), 2)
        ar = torch.arange(4)
        loss = (out["mu"][ar, ct] - y).pow(2).mean()
        for sm in out["step_mus"][:-1]:
            loss = loss + 0.1 * (sm[ar, ct] - y).pow(2).mean()
        opt.zero_grad()
        loss.backward()
        opt.step()
        self.assertTrue(torch.isfinite(loss).item())
        for name, p in model.named_parameters():
            self.assertTrue(torch.isfinite(p).all().item(), f"Non-finite param after MuonCWD step: {name}")

    def test_flashzoi_hybrid_forward_backward(self) -> None:
        torch.manual_seed(17)
        dummy_fz = nn.Sequential(
            nn.Conv1d(4, 32, kernel_size=7, padding=3),
            nn.MaxPool1d(2),
            nn.Conv1d(32, 64, kernel_size=5, padding=2),
            nn.MaxPool1d(2),
            nn.Conv1d(64, 64, kernel_size=5, padding=2),
            nn.MaxPool1d(2),
        )
        model = self.hl["FlashzoiHelixLoopOmegaModel"](
            dummy_fz,
            flashzoi_out_ch=64,
            w_fwd=self.w_fwd,
            w_rc=self.w_rc,
            motif_ids=self.motif_ids,
            out_dim=3,
            core_dim=128,
            n_loops=4,
            n_heads=4,
        )
        x = torch.zeros(4, 4, 230)
        x.scatter_(1, torch.randint(0, 4, (4, 230)).unsqueeze(1), 1.0)
        ct = torch.tensor([0, 1, 2, 0], dtype=torch.long)
        out = model.forward_loops(x, return_all_steps=True, cell_idx=ct)
        self.assertEqual(out["mu"].shape, (4, 3))
        loss = out["mu"].pow(2).mean()
        loss.backward()
        self.assertTrue(torch.isfinite(loss).item())
        model.eval()
        mu_eval = model(x)
        self.assertEqual(mu_eval.shape, (4, 3))
        self.assertTrue(torch.isfinite(mu_eval).all().item())

    def test_helixloop_v4_cured_invariants_and_gradients(self) -> None:
        from epigenesis_studio.helixloop import HelixLoopV4Model

        torch.manual_seed(23)
        for out_dim in (1, 3):
            model = HelixLoopV4Model(
                self.w_fwd,
                self.w_rc,
                self.motif_ids,
                out_dim=out_dim,
                core_dim=128,
                n_loops=4,
                n_heads=4,
            ).eval()
            rng = np.random.default_rng(31 + out_dim)
            bases = np.array(list("ACGT"))
            ref = "".join(rng.choice(bases, size=230))
            alt = ref[:100] + ("T" if ref[100] != "T" else "A") + ref[101:]
            self.assertEqual(model.predict_variant(ref, ref)["delta_effect"], 0.0)
            r_fa = model.predict_variant(ref, alt)
            r_af = model.predict_variant(alt, ref)
            self.assertAlmostEqual(r_fa["delta_effect"] + r_af["delta_effect"], 0.0, places=6)
            rc = self.hl["reverse_complement"]
            r_rc = model.predict_variant(rc(ref), rc(alt))
            self.assertAlmostEqual(r_fa["delta_effect"], r_rc["delta_effect"], places=6)

            model.train()
            x = torch.zeros(4, 4, 230)
            x.scatter_(1, torch.randint(0, 4, (4, 230)).unsqueeze(1), 1.0)
            ct = torch.randint(0, out_dim, (4,)) if out_dim > 1 else None
            out = model.forward_loops(x, return_all_steps=True, cell_idx=ct)
            loss = out["mu"].pow(2).mean()
            loss.backward()
            self.assertTrue(torch.isfinite(loss).item())

    def test_bundled_helixloop_v4_engine_and_studio_receipt_replay(self) -> None:
        from epigenesis_studio import helixloop_engine, server

        self.assertTrue(helixloop_engine.available())
        ident = helixloop_engine.identity()
        self.assertEqual(ident["model"], "epigenesis-helixloop-v4-multicell")
        self.assertEqual(ident["cell_lines"], ["HepG2", "K562", "WTC11"])

        ref = (
            "AGGACCGGATCAACTGCGCATGCGATCATGCTAGCATGCTAGCATGCTGACCTGAGGAGAACTGTGCT"
            "CCGCCTTCAGAGTACCACCGAAATCTGTGCAGAGGACAACGCAGCTCCGCCCTCGCGGTGCTCTCCTG"
            "GAGAACTGTGCTCCGCCTTCAGAGTACCACCGAAATCTGTGCAGAGGACAACGCAGCTCCGCCCTCGC"
            "GGTGCTCTCATTGCGTGAACCGA"
        )
        alt = ref[:115] + ("A" if ref[115] != "A" else "G") + ref[116:]
        alt_b = ref[:135] + ("C" if ref[135] != "C" else "T") + ref[136:]
        alt_ab = alt[:135] + alt_b[135] + alt[136:]
        rc = self.hl["reverse_complement"]

        preds = helixloop_engine.predict([
            dict(id="rev", reference=ref, sequence=ref),
            dict(id="fwd", reference=ref, sequence=alt),
            dict(id="anti", reference=alt, sequence=ref),
            dict(id="rc", reference=rc(ref), sequence=rc(alt)),
        ])
        self.assertEqual(preds[0]["effect"], 0.0)
        self.assertEqual(set(preds[1]["cell_lines"]), {"HepG2", "K562", "WTC11"})
        self.assertAlmostEqual(preds[1]["effect"] + preds[2]["effect"], 0.0, places=6)
        self.assertAlmostEqual(preds[1]["effect"], preds[3]["effect"], places=6)
        self.assertGreater(preds[1]["uncertainty_sigma"], 0.0)
        self.assertGreater(len(preds[1]["top_motifs"]), 0)

        epi = helixloop_engine.decompose_epistasis(ref, alt, alt_b, alt_ab)
        self.assertTrue(np.isfinite(epi["non_additive_epistasis"]))

        receipt = server.score(dict(mode="pairs", model="helixloop", rows=[dict(id="v1", reference=ref, sequence=alt)]))
        check = server.verify(receipt)
        self.assertTrue(check["ok"], check)


if __name__ == "__main__":
    unittest.main()

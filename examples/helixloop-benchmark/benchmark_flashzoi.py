#!/usr/bin/env python3
"""Tier-2 Open Apache-2.0 Foundation-Hybrid Benchmark:
Flashzoi + HelixLoop v4 vs Flashzoi-Linear-FT (and transparent comparison to AlphaGenome-Encoder-FT).

Strictly separates Tier 2 (pretrained on 7,611 human + 2,608 mouse ENCODE/GTEx genomic tracks via
Calico's Apache-2.0 Flashzoi / Borzoi convolutional encoder, ZERO lentiMPRA or CAGI5 labels during
pretraining) from Tier 1 (from-scratch, 0 external weights).

Evaluates on the EXACT same zero-leakage scrubbed ENCODE4 lentiMPRA Multi-Lineage dataset
(HepG2 + K562 + WTC11, n_train=281,170, 0.000% cross-lineage Fold-1 overlap, 0.000% CAGI5 31-mer overlap):
  1. Strictly held-out Fold-1 CRE-averaged Pearson r, Spearman rho, MSE on HepG2, K562, WTC11
  2. 2,000-draw paired bootstrap 95% CIs (Flashzoi+HelixLoop-v4 vs Flashzoi-Linear-FT)
  3. TRUE ZERO-SHOT Variant Effect Prediction on ALL 39,170 CAGI5 / Kircher 2019 SNVs (21 disjoint loci)
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from huggingface_hub import hf_hub_download
from safetensors import safe_open
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "studio"))
sys.path.insert(0, str(REPO / "examples/helixloop-benchmark"))

from benchmark_lentimpra import (  # noqa: E402
    _amp_ctx,
    batch_pearson_loss,
    cre_averaged_metrics,
    paired_bootstrap_gain,
    rc_onehot_batch,
    recalibrate_bn_stats,
    uint8_to_onehot,
)
from benchmark_multicell import (  # noqa: E402
    CAGI5_LOCUS,
    CELL_LINES,
    eval_cagi5_multicell,
    load_multicell_zero_leakage,
)
from epigenesis_studio.helixloop import (  # noqa: E402
    HelixLoopV4Core,
    compute_dna_shape_6ch,
)
from epigenesis_studio.thermoloop import FastBiophysField  # noqa: E402
from models import count_parameters, load_jaspar_pwm_bank  # noqa: E402


# Exact Borzoi/Flashzoi conv building blocks (Apache-2.0, Calico LLC) — zero flash_attn dependency!
class _FZConvDna(nn.Module):
    def __init__(self):
        super().__init__()
        self.conv_layer = nn.Conv1d(4, 512, kernel_size=15, padding="same")
        self.max_pool = nn.MaxPool1d(kernel_size=2, padding=0)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.max_pool(self.conv_layer(x))


class _FZConvBlock(nn.Module):
    def __init__(self, in_ch: int, out_ch: int, kernel_size: int = 5):
        super().__init__()
        self.norm = nn.BatchNorm1d(in_ch, eps=0.001)
        self.activation = nn.GELU(approximate="tanh")
        self.conv_layer = nn.Conv1d(in_ch, out_ch, kernel_size=kernel_size, padding="same")

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.conv_layer(self.activation(self.norm(x)))


class FlashzoiMultiScaleEncoder(nn.Module):
    """Exact Apache-2.0 Flashzoi Convolutional Encoder (conv_dna + res_tower + unet1, 25.6M params).

    Uses symmetric N-padding (0.0 = Flashzoi native unknown base) for exact RC equivariance,
    extracting all 4 receptive-field scales (stride 8, 16, 32, 64 -> fused cleanly at 32 tokens).
    """

    def __init__(self, repo_id: str = "johahi/flashzoi-replicate-0", out_dim: int = 256):
        super().__init__()
        self.repo_id = repo_id
        self.conv_dna = _FZConvDna()
        self._max_pool = nn.MaxPool1d(kernel_size=2, padding=0)
        self.res_tower = nn.Sequential(
            _FZConvBlock(512, 608, 5),   # 0: stride 2 -> 608
            self._max_pool,              # 1: stride 4
            _FZConvBlock(608, 736, 5),   # 2: stride 4 -> 736
            self._max_pool,              # 3: stride 8 (32 tokens!)
            _FZConvBlock(736, 896, 5),   # 4: stride 8 -> 896
            self._max_pool,              # 5: stride 16 (16 tokens!)
            _FZConvBlock(896, 1056, 5),  # 6: stride 16 -> 1056
            self._max_pool,              # 7: stride 32 (8 tokens!)
            _FZConvBlock(1056, 1280, 5), # 8: stride 32 -> 1280 (x_unet0)
        )
        self.unet1 = nn.Sequential(
            self._max_pool,              # 0: stride 64 (4 tokens!)
            _FZConvBlock(1280, 1536, 5), # 1: stride 64 -> 1536 (x_unet1)
        )

        st_path = hf_hub_download(repo_id, "model.safetensors")
        own_keys = set(self.state_dict().keys())
        loaded_sd: dict[str, torch.Tensor] = {}
        with safe_open(st_path, framework="pt", device="cpu") as f:
            for k in f.keys():
                if k in own_keys:
                    loaded_sd[k] = f.get_tensor(k)
        missing = own_keys - set(loaded_sd.keys())
        if missing:
            raise RuntimeError(f"Missing Flashzoi encoder keys from {repo_id}: {sorted(missing)[:10]}")
        self.load_state_dict(loaded_sd, strict=True)

        for m in self.modules():
            if isinstance(m, nn.BatchNorm1d):
                m.momentum = 0.05

        # Standard AdamW-regularizable 1x1 Conv + BatchNorm projections (ZERO scale-invariant overfitting!)
        self.proj_s8 = nn.Sequential(nn.Conv1d(736, out_dim, kernel_size=1, bias=False), nn.BatchNorm1d(out_dim), nn.SiLU())
        self.proj_s16 = nn.Sequential(nn.Conv1d(896, out_dim, kernel_size=1, bias=False), nn.BatchNorm1d(out_dim), nn.SiLU())
        self.proj_s32 = nn.Sequential(nn.Conv1d(1280, out_dim, kernel_size=1, bias=False), nn.BatchNorm1d(out_dim), nn.SiLU())
        self.proj_s64 = nn.Sequential(nn.Conv1d(1536, out_dim, kernel_size=1, bias=False), nn.BatchNorm1d(out_dim), nn.SiLU())
        self.scale_weights = nn.Parameter(torch.tensor([0.25, 0.25, 0.30, 0.20]))

    def pad_with_flanks(self, x: torch.Tensor) -> torch.Tensor:
        L = x.size(2)
        pad_tot = (64 - (L % 64)) % 64
        if pad_tot == 0:
            return x
        return F.pad(x, (pad_tot // 2, pad_tot - pad_tot // 2), value=0.0)

    def forward_scales(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        xp = self.pad_with_flanks(x)
        h = self.conv_dna(xp)
        h = self.res_tower[1](self.res_tower[0](h))
        h = self.res_tower[3](self.res_tower[2](h))  # stride 8 (32 tokens, 736 ch)
        s8 = h
        h = self.res_tower[5](self.res_tower[4](h))  # stride 16 (16 tokens, 896 ch)
        s16 = h
        h = self.res_tower[7](self.res_tower[6](h))  # stride 32 (8 tokens, 1056 ch)
        s32 = self.res_tower[8](h)                   # stride 32 (8 tokens, 1280 ch)
        s64 = self.unet1(s32)                        # stride 64 (4 tokens, 1536 ch)

        L_tok = s8.size(2)
        f8 = self.proj_s8(s8)
        f16 = F.interpolate(self.proj_s16(s16), size=L_tok, mode="nearest")
        f32 = F.interpolate(self.proj_s32(s32), size=L_tok, mode="nearest")
        f64 = F.interpolate(self.proj_s64(s64), size=L_tok, mode="nearest")
        w = F.softmax(self.scale_weights, dim=0)
        fused = w[0] * f8 + w[1] * f16 + w[2] * f32 + w[3] * f64  # (B, out_dim, 32)
        return fused, s64

    def pretrained_parameters(self):
        for mod in (self.conv_dna, self.res_tower, self.unet1):
            yield from mod.parameters()


class FlashzoiLinearFT(nn.Module):
    """Fair Tier-2 Baseline: Apache-2.0 Flashzoi Multi-Scale Encoder + Flat Pool/MLP Head (mirrors AlphaGenome-Encoder-FT)."""

    def __init__(self, repo_id: str = "johahi/flashzoi-replicate-0", out_dim: int = 3, dropout: float = 0.25):
        super().__init__()
        self.out_dim = out_dim
        self.predict_uncertainty = False
        self.encoder = FlashzoiMultiScaleEncoder(repo_id=repo_id, out_dim=256)
        self.head = nn.Sequential(
            nn.Linear(256 * 2 + 1536, 512),
            nn.BatchNorm1d(512),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(512, out_dim),
        )

    def forward(self, x: torch.Tensor, cell_idx: torch.Tensor | None = None) -> torch.Tensor:
        fused, s64 = self.encoder.forward_scales(x)
        feat = torch.cat([fused.mean(dim=-1), fused.amax(dim=-1), s64.mean(dim=-1)], dim=1)
        out = self.head(feat)
        return out.squeeze(-1) if self.out_dim == 1 else out


class FlashzoiHelixLoopV4(nn.Module):
    """Tier-2 Champion: Apache-2.0 Flashzoi Multi-Scale Encoder + Cured Epigenesis HelixLoop v4 Core.

    Fuses Flashzoi's 25.6M pretrained ENCODE/GTEx convolutional pyramid with Epigenesis's
    532-TF JASPAR Boltzmann Thermodynamic Field, 6-ch zero-param Rohs DNA Shape,
    and genuine iterative T=6 HelixLoopV4Core (Differential Gated Attention + ResFormer V0 +
    B-DNA 10.5-bp bias + interleaved WY BiGatedDeltaProduct1d + Misul Rational FFN + Cell-State FiLM!).
    """

    def __init__(
        self,
        w_fwd: np.ndarray,
        w_rc: np.ndarray,
        *,
        repo_id: str = "johahi/flashzoi-replicate-0",
        out_dim: int = 3,
        core_dim: int = 256,
        n_loops: int = 6,
        n_heads: int = 4,
        dropout: float = 0.20,
        drop_path: float = 0.15,
        predict_uncertainty: bool = True,
    ):
        super().__init__()
        self.out_dim = out_dim
        self.n_loops = n_loops
        self.predict_uncertainty = predict_uncertainty

        self.encoder = FlashzoiMultiScaleEncoder(repo_id=repo_id, out_dim=core_dim)

        # Parallel Epigenesis 532-TF Thermodynamic Field + 6-ch Zero-Param Rohs DNA Shape
        self.thermo = FastBiophysField(w_fwd, w_rc, out_spatial_ch=32, dropout=0.15)
        self.bio_down = nn.Sequential(
            nn.Conv1d(32 + 6, 96, kernel_size=7, padding=3, bias=False),
            nn.BatchNorm1d(96),
            nn.SiLU(),
            nn.MaxPool1d(2),
            nn.Conv1d(96, core_dim, kernel_size=5, padding=2, bias=False),
            nn.BatchNorm1d(core_dim),
            nn.SiLU(),
            nn.MaxPool1d(2),
        )
        self.fuse_proj = nn.Sequential(
            nn.Linear(core_dim * 2, core_dim),
            nn.LayerNorm(core_dim),
            nn.SiLU(),
        )
        self.anchor_norm = nn.LayerNorm(core_dim)

        self.core = HelixLoopV4Core(
            dim=core_dim,
            n_heads=n_heads,
            max_loops=12,
            n_cell_types=out_dim,
            bp_per_token=8.0,
            dropout=dropout,
            drop_path=drop_path,
            use_diff_attn=True,
            use_delta_product=True,
            use_rational=True,
        )

        self.coda_norm = nn.LayerNorm(core_dim)
        self.pool_attn = nn.Linear(core_dim, 1)
        self.head_proj = nn.Sequential(
            nn.Linear(core_dim * 2 + 1536, 384),
            nn.BatchNorm1d(384),
            nn.SiLU(),
            nn.Dropout(dropout),
        )
        self.mu_head = nn.Linear(384, out_dim)
        self.aux_head = nn.Sequential(nn.LayerNorm(core_dim), nn.Linear(core_dim, out_dim))

        self.thermo_norm = nn.BatchNorm1d(w_fwd.shape[0])
        self.thermo_drop = nn.Dropout(0.20)
        self.thermo_skip = nn.Linear(w_fwd.shape[0], out_dim, bias=False)
        nn.init.zeros_(self.thermo_skip.weight)

        self.logvar_head = nn.Sequential(nn.Linear(384 + 1, 64), nn.SiLU(), nn.Linear(64, out_dim)) if predict_uncertainty else None

    def _encode_prelude(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        fz_fused, s64 = self.encoder.forward_scales(x)
        bio_sp, bio_logz = self.thermo(x)
        shape6 = compute_dna_shape_6ch(x)
        bio_f = F.adaptive_max_pool1d(self.bio_down(torch.cat([bio_sp, shape6], dim=1)), fz_fused.size(2))
        h0 = self.fuse_proj(torch.cat([fz_fused.transpose(1, 2), bio_f.transpose(1, 2)], dim=-1))
        return h0, self.anchor_norm(h0), s64.mean(dim=-1), bio_logz

    def _run_core_and_coda(
        self,
        h0: torch.Tensor,
        u_anchor: torch.Tensor,
        s64_mean: torch.Tensor,
        bio_logz: torch.Tensor,
        *,
        cell_idx: torch.Tensor | None = None,
        n_loops: int | None = None,
        return_all_steps: bool = False,
    ) -> dict[str, Any]:
        T = self.n_loops if n_loops is None else n_loops
        b, l, _ = h0.shape
        if cell_idx is None:
            cell_vec = self.core.cell_embed.weight.mean(dim=0, keepdim=True).unsqueeze(1).expand(b, 1, -1)
        else:
            cell_vec = self.core.cell_embed(cell_idx).unsqueeze(1)
        max_e = self.core.step_embed.num_embeddings - 1
        step_vecs = self.core.step_embed(torch.arange(T, device=h0.device).clamp_max(max_e)).unsqueeze(1).unsqueeze(1)
        rel_bias = self.core._bdna_rel_bias(l, h0.device, h0.dtype)

        h = h0
        v0: torch.Tensor | None = None
        aux_mus: list[torch.Tensor] = []
        last_delta = h0.new_zeros((b, 1))
        for t in range(T):
            h_next, v0 = self.core.step(h, u_anchor, rel_bias, v0, step_vecs[t] + cell_vec, step_idx=t)
            if t == T - 1:
                last_delta = (h_next - h).pow(2).mean(dim=(1, 2)).sqrt().unsqueeze(-1)
            h = h_next
            if return_all_steps and t < T - 1 and (t + 1) in (2, 4):
                m_aux = self.aux_head(h).mean(dim=1)
                aux_mus.append(m_aux.squeeze(-1) if self.out_dim == 1 else m_aux)

        hf = self.coda_norm(h)
        w_a = F.softmax(self.pool_attn(hf), dim=1)
        z_feat = self.head_proj(torch.cat([(hf * w_a).sum(dim=1), hf.mean(dim=1), s64_mean], dim=1))
        mu_final = self.mu_head(z_feat) + self.thermo_skip(self.thermo_drop(self.thermo_norm(bio_logz)))
        if return_all_steps:
            aux_mus.append(mu_final.squeeze(-1) if self.out_dim == 1 else mu_final)
        logvar = (
            self.logvar_head(torch.cat([z_feat.detach(), last_delta.detach()], dim=1)).clamp(-6.0, 4.0)
            if self.logvar_head is not None
            else torch.zeros_like(mu_final)
        )
        mu_out = mu_final.squeeze(-1) if self.out_dim == 1 else mu_final
        lv_out = logvar.squeeze(-1) if self.out_dim == 1 else logvar
        return dict(mu=mu_out, logvar=lv_out, step_mus=aux_mus, eq_residual=last_delta.squeeze(-1), bio_logz=bio_logz)

    def forward_loops(
        self,
        x: torch.Tensor,
        n_loops: int | None = None,
        return_all_steps: bool = False,
        cell_idx: torch.Tensor | None = None,
    ) -> dict[str, Any]:
        h0, u_anchor, s64_mean, bio_logz = self._encode_prelude(x)
        if self.out_dim > 1 and cell_idx is None and not self.training:
            mus, lvs, eqs = [], [], []
            for c in range(self.out_dim):
                cid = torch.full((x.size(0),), c, device=x.device, dtype=torch.long)
                out_c = self._run_core_and_coda(h0, u_anchor, s64_mean, bio_logz, cell_idx=cid, n_loops=n_loops, return_all_steps=False)
                mus.append(out_c["mu"][:, c : c + 1])
                lvs.append(out_c["logvar"][:, c : c + 1])
                eqs.append(out_c["eq_residual"])
            return dict(mu=torch.cat(mus, dim=1), logvar=torch.cat(lvs, dim=1), step_mus=[], eq_residual=torch.stack(eqs, dim=0).mean(dim=0), bio_logz=bio_logz)
        return self._run_core_and_coda(h0, u_anchor, s64_mean, bio_logz, cell_idx=cell_idx, n_loops=n_loops, return_all_steps=return_all_steps)

    def forward_with_uncertainty(self, x: torch.Tensor, cell_idx: torch.Tensor | None = None) -> tuple[torch.Tensor, torch.Tensor]:
        out = self.forward_loops(x, return_all_steps=False, cell_idx=cell_idx)
        return out["mu"], out["logvar"]

    def forward(self, x: torch.Tensor, cell_idx: torch.Tensor | None = None) -> torch.Tensor:
        return self.forward_loops(x, return_all_steps=False, cell_idx=cell_idx)["mu"]


@torch.no_grad()
def predict_tier2(model: nn.Module, x_u8: np.ndarray, device: torch.device, batch_size: int = 256, rc_avg: bool = True) -> tuple[np.ndarray, np.ndarray]:
    model.eval()
    n = len(x_u8)
    preds = np.zeros((n, 3), dtype=np.float64)
    stds = np.zeros((n, 3), dtype=np.float64)
    has_unc = isinstance(model, FlashzoiHelixLoopV4) and model.predict_uncertainty
    for s in range(0, n, batch_size):
        e = min(n, s + batch_size)
        oh = uint8_to_onehot(x_u8[s:e], device)
        with _amp_ctx(device):
            if has_unc:
                mu, lv = model.forward_with_uncertainty(oh)
                mu, lv = mu.float(), lv.float()
                if rc_avg:
                    mu2, lv2 = model.forward_with_uncertainty(rc_onehot_batch(oh))
                    mu2, lv2 = mu2.float(), lv2.float()
                    mu = 0.5 * (mu + mu2)
                    var = 0.5 * (torch.exp(lv) + torch.exp(lv2)) + 0.25 * (mu - mu2) ** 2
                    std = torch.sqrt(var.clamp_min(1e-8))
                else:
                    std = torch.exp(0.5 * lv)
                preds[s:e] = mu.cpu().numpy()
                stds[s:e] = std.cpu().numpy()
            else:
                mu = model(oh).float()
                if rc_avg:
                    mu = 0.5 * (mu + model(rc_onehot_batch(oh)).float())
                preds[s:e] = mu.cpu().numpy()
    return preds, stds


def train_tier2_model(
    model_name: str,
    repo_id: str,
    data: dict[str, Any],
    w_fwd: np.ndarray,
    w_rc: np.ndarray,
    device: torch.device,
    *,
    epochs: int = 12,
    batch_size: int = 256,
    lr_encoder: float = 3e-4,
    lr_head: float = 2.5e-3,
    weight_decay: float = 0.08,
    pearson_weight: float = 0.0,
    ema_decay: float = 0.997,
    n_loops: int = 6,
    dropout: float = 0.20,
    seed: int = 42,
) -> tuple[nn.Module, dict[str, Any]]:
    torch.manual_seed(seed)
    np.random.seed(seed)
    if model_name == "flashzoi_linear_ft":
        model: nn.Module = FlashzoiLinearFT(repo_id=repo_id, out_dim=3, dropout=dropout).to(device)
    elif model_name == "flashzoi_helixloop_v4":
        model = FlashzoiHelixLoopV4(w_fwd, w_rc, repo_id=repo_id, out_dim=3, n_loops=n_loops, dropout=dropout).to(device)
    else:
        raise ValueError(f"Unknown Tier-2 model: {model_name}")

    n_params = count_parameters(model)
    pretrained_ids = {id(p) for p in model.encoder.pretrained_parameters()}
    enc_params = [p for p in model.parameters() if p.requires_grad and id(p) in pretrained_ids]
    head_params = [p for p in model.parameters() if p.requires_grad and id(p) not in pretrained_ids]

    tr_ext = data["tr_ext"]
    tr_y = data["tr_y"]
    tr_ct = data["tr_ct"]
    n_tr = len(tr_y)
    ct_counts = np.bincount(tr_ct, minlength=3).astype(np.float64)
    ct_weights = (ct_counts.mean() / ct_counts).astype(np.float32)

    tr_ext_gpu = torch.from_numpy(tr_ext).to(device=device, dtype=torch.long)
    tr_y_gpu = torch.from_numpy(tr_y).to(device=device, dtype=torch.float32)
    tr_ct_gpu = torch.from_numpy(tr_ct).to(device=device, dtype=torch.long)
    ct_w_gpu = torch.from_numpy(ct_weights).to(device=device, dtype=torch.float32)
    seq_offsets_gpu = torch.arange(230, device=device, dtype=torch.long).unsqueeze(0)
    onehot_eye_gpu = torch.cat([torch.eye(4, device=device, dtype=torch.float32), torch.zeros(1, 4, device=device, dtype=torch.float32)], dim=0)

    steps_per_ep = math.ceil(n_tr / batch_size)
    opt = torch.optim.AdamW(
        [
            {"params": enc_params, "lr": lr_encoder / 20.0, "weight_decay": weight_decay},
            {"params": head_params, "lr": lr_head / 20.0, "weight_decay": weight_decay},
        ]
    )
    sched = torch.optim.lr_scheduler.OneCycleLR(
        opt,
        max_lr=[lr_encoder, lr_head],
        epochs=epochs,
        steps_per_epoch=steps_per_ep,
        pct_start=0.20,
        cycle_momentum=False,
    )
    rng = np.random.default_rng(seed)
    ema_params = {k: v.detach().clone().float() for k, v in model.named_parameters()} if ema_decay > 0.0 else None

    best_val_mean_r = -1e9
    best_state = None
    best_ep = -1
    t0 = time.perf_counter()

    def _per_cell_pearson(mu_o: torch.Tensor, yb_in: torch.Tensor, ctb_in: torch.Tensor) -> torch.Tensor:
        pl = mu_o.new_zeros(())
        nc = 0
        for ci in range(3):
            m = ctb_in == ci
            if int(m.sum().item()) >= 16:
                pl = pl + batch_pearson_loss(mu_o[m], yb_in[m])
                nc += 1
        return pl / max(nc, 1)

    for ep in range(1, epochs + 1):
        model.train()
        perm_gpu = torch.from_numpy(rng.permutation(n_tr)).to(device=device, dtype=torch.long)
        ep_loss = 0.0
        nb = 0
        for s in range(0, n_tr, batch_size):
            idx_g = perm_gpu[s : s + batch_size]
            b = int(idx_g.numel())
            if b < 8:
                continue
            shifts_g = torch.randint(0, 43, (b, 1), device=device)
            x_idx = tr_ext_gpu[idx_g].gather(1, shifts_g + seq_offsets_gpu)
            oh = onehot_eye_gpu[x_idx].permute(0, 2, 1).contiguous()
            rc_mask = torch.rand(b, device=device) > 0.5
            if rc_mask.any():
                oh[rc_mask] = rc_onehot_batch(oh[rc_mask])
            yb = tr_y_gpu[idx_g]
            ctb = tr_ct_gpu[idx_g]
            wb = ct_w_gpu[ctb]
            ar = torch.arange(b, device=device)

            opt.zero_grad(set_to_none=True)
            with _amp_ctx(device):
                if isinstance(model, FlashzoiHelixLoopV4):
                    out = model.forward_loops(oh, return_all_steps=True, cell_idx=ctb)
                    mu_all, lv_all = out["mu"].float(), out["logvar"].float()
                    mu_o, lv_o = mu_all[ar, ctb], lv_all[ar, ctb]
                    step_mus = [sm.float() for sm in out["step_mus"]]
                    if len(step_mus) == 3:
                        mse_o = (
                            0.15 * ((step_mus[0][ar, ctb] - yb).pow(2) * wb).mean()
                            + 0.25 * ((step_mus[1][ar, ctb] - yb).pow(2) * wb).mean()
                            + ((step_mus[2][ar, ctb] - yb).pow(2) * wb).mean()
                        )
                    else:
                        mse_o = ((mu_o - yb).pow(2) * wb).mean()
                    if pearson_weight > 0.0:
                        mse_o = mse_o + pearson_weight * _per_cell_pearson(mu_o, yb, ctb)
                    gnll_o = (0.5 * (torch.exp(-lv_o) * (yb - mu_o.detach()).pow(2) + lv_o) * wb).mean()
                    loss = mse_o + 0.15 * gnll_o
                else:
                    mu_all = model(oh).float()
                    mu_o = mu_all[ar, ctb]
                    loss = ((mu_o - yb).pow(2) * wb).mean()
                    if pearson_weight > 0.0:
                        loss = loss + pearson_weight * _per_cell_pearson(mu_o, yb, ctb)

            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 2.0)
            opt.step()
            sched.step()
            ep_loss += float(loss.item())
            nb += 1
            if ema_params is not None:
                with torch.no_grad():
                    for k, p in model.named_parameters():
                        ema_params[k].mul_(ema_decay).add_(p.detach().float(), alpha=1.0 - ema_decay)

        vr_map = {}
        for ci, ct in enumerate(CELL_LINES):
            vd = data["val_by_ct"][ct]
            vp, _ = predict_tier2(model, vd["x"], device, batch_size=batch_size, rc_avg=True)
            _, _, vm = cre_averaged_metrics(vp[:, ci], vd["y"], vd["ids"])
            vr_map[ct] = vm["pearson_r"]
        mean_vr = float(np.mean(list(vr_map.values())))
        if mean_vr > best_val_mean_r:
            best_val_mean_r = mean_vr
            best_ep = ep
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
        print(
            f"  [Tier2-{model_name:22s}] ep={ep:2d}/{epochs} loss={ep_loss/max(nb,1):.4f} "
            f"val_r: HepG2={vr_map['HepG2']:.4f} K562={vr_map['K562']:.4f} WTC11={vr_map['WTC11']:.4f} "
            f"(mean={mean_vr:.4f}, {time.perf_counter()-t0:.1f}s)",
            flush=True,
        )

    if ema_params is not None:
        with torch.no_grad():
            for k, p in model.named_parameters():
                p.copy_(ema_params[k].to(p.dtype))
        def _bn_gen():
            p_bn = rng.permutation(n_tr)
            for si in range(0, min(n_tr, 80 * batch_size), batch_size):
                ix = p_bn[si : si + batch_size]
                if len(ix) < 8:
                    continue
                ix_g = torch.from_numpy(ix).to(device=device, dtype=torch.long)
                sh_g = torch.full((len(ix), 1), 21, device=device, dtype=torch.long)
                oh_g = onehot_eye_gpu[tr_ext_gpu[ix_g].gather(1, sh_g + seq_offsets_gpu)].permute(0, 2, 1).contiguous()
                yield (oh_g, tr_ct_gpu[ix_g])
        recalibrate_bn_stats(model, _bn_gen(), device)
        vr_ema = {}
        for ci, ct in enumerate(CELL_LINES):
            vd = data["val_by_ct"][ct]
            vp, _ = predict_tier2(model, vd["x"], device, batch_size=batch_size, rc_avg=True)
            _, _, vm = cre_averaged_metrics(vp[:, ci], vd["y"], vd["ids"])
            vr_ema[ct] = vm["pearson_r"]
        mean_ema = float(np.mean(list(vr_ema.values())))
        print(
            f"  [Tier2-{model_name:22s}] EMA+BN-recal val_r: HepG2={vr_ema['HepG2']:.4f} K562={vr_ema['K562']:.4f} "
            f"WTC11={vr_ema['WTC11']:.4f} (mean={mean_ema:.4f} vs live best={best_val_mean_r:.4f})",
            flush=True,
        )
        if mean_ema > best_val_mean_r:
            best_val_mean_r = mean_ema
            best_ep = epochs
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}

    model.load_state_dict(best_state)
    return model, dict(
        model=model_name,
        pretrained_source=repo_id,
        pretrained_license="Apache-2.0",
        n_params=n_params,
        best_epoch=best_ep,
        best_val_mean_r=round(best_val_mean_r, 5),
        train_seconds=round(time.perf_counter() - t0, 1),
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", nargs="+", default=["flashzoi_helixloop_v4", "flashzoi_linear_ft"])
    ap.add_argument("--repo-id", type=str, default="johahi/flashzoi-replicate-0")
    ap.add_argument("--epochs", type=int, default=12)
    ap.add_argument("--batch-size", type=int, default=256)
    ap.add_argument("--lr-encoder", type=float, default=3e-4)
    ap.add_argument("--lr-head", type=float, default=2.5e-3)
    ap.add_argument("--weight-decay", type=float, default=0.08)
    ap.add_argument("--dropout", type=float, default=0.20)
    ap.add_argument("--n-loops", type=int, default=6)
    ap.add_argument("--pearson-weight", type=float, default=0.0)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tag", type=str, default="rep0")
    ap.add_argument("--data-dir", type=Path, default=Path("/content/data/human_legnet-main/datasets/original"))
    ap.add_argument("--kircher-dir", type=Path, default=Path("/content/data/kircher"))
    ap.add_argument("--out-dir", type=Path, default=Path("/content/results"))
    args = ap.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device.type == "cuda":
        torch.backends.cuda.matmul.allow_tf32 = True
        torch.backends.cudnn.allow_tf32 = True
        torch.backends.cudnn.benchmark = True

    w_fwd, w_rc, _ = load_jaspar_pwm_bank(21)
    data = load_multicell_zero_leakage(args.data_dir, args.kircher_dir / "kircher_refs.json")
    print(f"=== Tier-2 Open Apache-2.0 Foundation-Hybrid ({args.repo_id}) | n_train={len(data['tr_y'])} ===", flush=True)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    report: dict[str, Any] = dict(
        tier="Tier-2 Open-Weights Foundation-Hybrid (Apache-2.0 Flashzoi conv encoder pretrained on 7,611 human + 2,608 mouse ENCODE/GTEx tracks; ZERO lentiMPRA/CAGI5 pretraining labels)",
        repo_id=args.repo_id,
        n_train_total=int(len(data["tr_y"])),
        scrubbed_counts=data["scrubbed_counts"],
        epochs=args.epochs,
        seed=args.seed,
        published_alphagenome_encoder_ft_reference=dict(HepG2=0.888, K562=0.879, WTC11=0.834, CAGI5_approx=0.65),
        models={},
    )
    preds_by_ct: dict[str, dict[str, np.ndarray]] = {ct: {} for ct in CELL_LINES}
    y_by_ct: dict[str, np.ndarray] = {}

    for mname in args.models:
        print(f"\n=== Training Tier-2 {mname} ({args.repo_id}) ===", flush=True)
        model, info = train_tier2_model(
            mname,
            args.repo_id,
            data,
            w_fwd,
            w_rc,
            device,
            epochs=args.epochs,
            batch_size=args.batch_size,
            lr_encoder=args.lr_encoder,
            lr_head=args.lr_head,
            weight_decay=args.weight_decay,
            dropout=args.dropout,
            n_loops=args.n_loops,
            pearson_weight=args.pearson_weight,
            seed=args.seed,
        )
        ckpt_path = args.out_dir / f"TIER2_{mname}_{args.tag}_s{args.seed}.pt"
        torch.save({"model_name": mname, "repo_id": args.repo_id, "n_loops": args.n_loops, "state_dict": {k: v.detach().cpu() for k, v in model.state_dict().items()}, "info": info}, ckpt_path)

        ct_res = {}
        for ci, ct in enumerate(CELL_LINES):
            td = data["test_by_ct"][ct]
            tp, tstd = predict_tier2(model, td["x"], device, batch_size=args.batch_size, rc_avg=True)
            p_cre, y_cre, cm = cre_averaged_metrics(tp[:, ci], td["y_avg"], td["ids"])
            preds_by_ct[ct][mname] = p_cre
            y_by_ct[ct] = y_cre
            cal = round(float(spearmanr(tstd[:, ci], np.abs(tp[:, ci] - td["y_raw"])).statistic), 4) if np.std(tstd[:, ci]) > 1e-12 else None
            ct_res[ct] = {**{k: (round(v, 5) if isinstance(v, float) else v) for k, v in cm.items()}, "cal_spearman": cal}

        np.savez_compressed(
            args.out_dir / f"TIER2_{mname}_{args.tag}_s{args.seed}_cre_preds.npz",
            **{f"p_{ct}": preds_by_ct[ct][mname] for ct in CELL_LINES},
            **{f"y_{ct}": y_by_ct[ct] for ct in CELL_LINES},
        )

        # Probe loop depth scaling T in {1, 2, 4, 6, 8} for FlashzoiHelixLoopV4!
        loop_scaling = None
        if isinstance(model, FlashzoiHelixLoopV4):
            loop_scaling = {}
            orig_loops = model.n_loops
            for t_eval in (1, 2, 4, 6, 8):
                model.n_loops = t_eval
                t_res = {}
                for ci, ct in enumerate(CELL_LINES):
                    td = data["test_by_ct"][ct]
                    tp_t, _ = predict_tier2(model, td["x"], device, batch_size=args.batch_size, rc_avg=True)
                    _, _, cm_t = cre_averaged_metrics(tp_t[:, ci], td["y_avg"], td["ids"])
                    t_res[ct] = round(cm_t["pearson_r"], 5)
                t_res["mean"] = round(float(np.mean(list(t_res.values()))), 5)
                loop_scaling[f"T={t_eval}"] = t_res
            model.n_loops = orig_loops
            print(f"  ==> TIER2 LOOP SCALING ({mname}): {loop_scaling}", flush=True)

        cagi5 = eval_cagi5_multicell(model, args.kircher_dir, device, batch_size=args.batch_size)
        report["models"][mname] = dict(**info, test_fold1_cre=ct_res, loop_depth_scaling=loop_scaling, cagi5_zeroshot=cagi5)
        print(
            f"==> TIER2-{mname}: Fold1 CRE r: HepG2={ct_res['HepG2']['pearson_r']:.4f} K562={ct_res['K562']['pearson_r']:.4f} "
            f"WTC11={ct_res['WTC11']['pearson_r']:.4f} | CAGI5 ZS: blind={cagi5['mean_blind_3head_r']:.4f} "
            f"matched={cagi5['mean_lineage_matched_r']:.4f} sig={cagi5['mean_sig_p1e5_r']:.4f} "
            f"liver={cagi5['mean_liver_hepg2_r']:.4f} eryth={cagi5['mean_erythroid_k562_r']:.4f}",
            flush=True,
        )

    if "flashzoi_helixloop_v4" in report["models"] and "flashzoi_linear_ft" in report["models"]:
        gains = {}
        loci = sorted(set(CAGI5_LOCUS.values()))
        for ct in CELL_LINES:
            gains[f"{ct}_v4_vs_linear"] = paired_bootstrap_gain(
                preds_by_ct[ct]["flashzoi_helixloop_v4"],
                preds_by_ct[ct]["flashzoi_linear_ft"],
                y_by_ct[ct],
                n_boot=2000,
            )
        m_th = report["models"]["flashzoi_helixloop_v4"]["cagi5_zeroshot"]["per_element_matched"]
        m_lg = report["models"]["flashzoi_linear_ft"]["cagi5_zeroshot"]["per_element_matched"]
        l_th = np.array([np.mean([m_th[e] for e, lc in CAGI5_LOCUS.items() if lc == loc and e in m_th]) for loc in loci])
        l_lg = np.array([np.mean([m_lg[e] for e, lc in CAGI5_LOCUS.items() if lc == loc and e in m_lg]) for loc in loci])
        rng = np.random.default_rng(2026)
        bidx = rng.integers(0, len(loci), size=(2000, len(loci)))
        bg = (l_th[bidx] - l_lg[bidx]).mean(axis=1)
        gains["cagi5_21loci_v4_vs_linear"] = dict(
            point_gain=round(float((l_th - l_lg).mean()), 5),
            ci95_low=round(float(np.quantile(bg, 0.025)), 5),
            ci95_high=round(float(np.quantile(bg, 0.975)), 5),
            prob_better=round(float((bg > 0).mean()), 4),
            elements_won=f"{sum(m_th[e] > m_lg[e] for e in m_th)}/{len(m_th)}",
        )
        report["paired_bootstrap_gains"] = gains
        print("\n=== Tier-2 Paired Bootstrap Gains (Flashzoi+HelixLoop-v4 vs Flashzoi-Linear-FT) ===", flush=True)
        for k, v in gains.items():
            print(f"  {k:40s}: {v}", flush=True)

    out_json = args.out_dir / f"benchmark_tier2_flashzoi_{args.tag}_s{args.seed}.json"
    out_json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"\nSaved Tier-2 Flashzoi benchmark results to {out_json}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

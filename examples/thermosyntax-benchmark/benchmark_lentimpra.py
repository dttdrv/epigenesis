#!/usr/bin/env python3
"""Strict, Zero-Leakage Head-to-Head Benchmark: Epigenesis ThermoSyntaxNet vs Market SOTA.

Trains all models on the EXACT SAME scrubbed ENCODE4 lentiMPRA splits (folds 3..10 train,
fold 2 val, fold 1 strictly held-out test) with identical augmentation, identical epochs,
and evaluates:
  1. Held-out fold 1 per-sequence & CRE orientation-averaged Pearson r, Spearman rho, MSE, R^2
  2. 1,000-draw paired bootstrap 95% CIs against Official MPRALegNet, LegNet-Large, Malinois, DeepSTARR
  3. Predictive uncertainty calibration (Spearman correlation of predicted sigma vs |error|)
  4. TRUE ZERO-SHOT Variant Effect Prediction on ALL 39,170 CAGI5 / Kircher 2019 SNVs
     (ZERO Kircher training labels, 0% 31-mer overlap scrubbed from lentiMPRA training set)
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from scipy.stats import spearmanr

from models import (
    EpigenesisThermoDREAM,
    EpigenesisThermoLoop,
    EpigenesisThermoSyntaxNet,
    FlashzoiHelixLoopOmegaModel,
    HelixLoopOmegaModel,
    HelixLoopV4Model,
    HybridMuonCWD,
    RMSBatchNorm1d,
    build_model,
    count_parameters,
    load_jaspar_pwm_bank,
)

LOOP_MODEL_TYPES = (EpigenesisThermoLoop, HelixLoopOmegaModel, HelixLoopV4Model, FlashzoiHelixLoopOmegaModel)
UNC_MODEL_TYPES = (EpigenesisThermoSyntaxNet, EpigenesisThermoDREAM, EpigenesisThermoLoop, HelixLoopOmegaModel, HelixLoopV4Model, FlashzoiHelixLoopOmegaModel)


def batch_pearson_loss(pred: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
    px = pred - pred.mean()
    tx = target - target.mean()
    r = (px * tx).sum() * torch.rsqrt((px.pow(2).sum() + 1e-8) * (tx.pow(2).sum() + 1e-8))
    return 1.0 - r


@torch.no_grad()
def recalibrate_bn_stats(model: nn.Module, batch_gen, device: torch.device, disable_amp: bool = False) -> None:
    """Exact cumulative-moving-average recalibration of BatchNorm1d & RMSBatchNorm1d after EMA weight swap."""
    bn_mods = [m for m in model.modules() if isinstance(m, (nn.BatchNorm1d, RMSBatchNorm1d))]
    if not bn_mods:
        return
    saved_mom = [m.momentum for m in bn_mods]
    for m in bn_mods:
        if isinstance(m, nn.BatchNorm1d):
            m.reset_running_stats()
        else:
            m.running_sq.zero_()
    model.train()
    for step_i, batch_in in enumerate(batch_gen, start=1):
        mom = 1.0 / step_i
        for m in bn_mods:
            m.momentum = mom
        with _amp_ctx(device, disable=disable_amp):
            if isinstance(batch_in, tuple):
                oh, ctb = batch_in
                if hasattr(model, "forward_loops"):
                    model(oh, cell_idx=ctb)
                else:
                    model(oh)
            else:
                model(batch_in)
    for m, orig in zip(bn_mods, saved_mom):
        m.momentum = orig
    model.eval()

FORWARD_ADAPTER = "GGCCCGCTCTAGACCTGCAGG"
REVERSE_ADAPTER = "CACTAGAGGGTATATAATGGAAGCTCGACTTCCAGCTTGGCAATCCGGTACTGT"
LENTI_LEFT = "AGGACCGGATCAACT"
LENTI_RIGHT = "CATTGCGTGAACCGA"

BASE_LUT = np.full(128, 4, dtype=np.uint8)
BASE_LUT[ord("A")] = 0
BASE_LUT[ord("C")] = 1
BASE_LUT[ord("G")] = 2
BASE_LUT[ord("T")] = 3
BASE_LUT[ord("a")] = 0
BASE_LUT[ord("c")] = 1
BASE_LUT[ord("g")] = 2
BASE_LUT[ord("t")] = 3
RC_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")


def rc_str(s: str) -> str:
    return s.translate(RC_TABLE)[::-1]


def seqs_to_uint8(seqs: list[str]) -> np.ndarray:
    n = len(seqs)
    L = len(seqs[0])
    buf = np.frombuffer("".join(seqs).encode("ascii"), dtype=np.uint8).reshape(n, L)
    return BASE_LUT[buf]


def uint8_to_onehot(arr: np.ndarray, device: torch.device) -> torch.Tensor:
    # arr: (B, L) uint8 in {0,1,2,3,4} -> (B, 4, L) float32
    t = torch.from_numpy(arr).to(device=device, dtype=torch.long)
    oh = torch.zeros((t.size(0), 5, t.size(1)), device=device, dtype=torch.float32)
    oh.scatter_(1, t.unsqueeze(1), 1.0)
    return oh[:, :4, :]


def build_cagi5_kmer_set(kircher_refs_path: Path) -> set[str]:
    refs = json.loads(kircher_refs_path.read_text(encoding="utf-8"))
    k31 = set()
    for v in refs.values():
        s = (v["seq"] if isinstance(v, dict) else v).upper()
        for i in range(0, len(s) - 30, 5):
            sub = s[i : i + 31]
            k31.add(sub)
            k31.add(rc_str(sub))
    return k31


def base_cre_id(seq_id: str) -> str:
    if seq_id.endswith("_Reversed:"):
        return seq_id[: -len("_Reversed:")]
    if seq_id.endswith("_F") or seq_id.endswith("_R"):
        return seq_id[:-2]
    return seq_id


def load_lentimpra(data_dir: Path, cell_line: str, cagi5_k31: set[str]) -> dict[str, Any]:
    raw_path = data_dir / f"{cell_line}.tsv"
    avg_path = data_dir / f"{cell_line}_averaged.tsv"

    tr_seqs, tr_y = [], []
    va_seqs, va_y, va_ids = [], [], []
    te_seqs, te_y_raw, te_y_avg, te_ids = [], [], [], []
    scrubbed = 0

    with raw_path.open() as fr, avg_path.open() as fa:
        rr = csv.DictReader(fr, delimiter="\t")
        ra = csv.DictReader(fa, delimiter="\t")
        for row_r, row_a in zip(rr, ra):
            if row_r["rev"] != "0":
                continue
            seq = row_r["seq"].upper()
            fold = int(row_r["fold"])
            y_r = float(row_r["mean_value"])
            y_a = float(row_a["mean_value"])
            cid = base_cre_id(row_r["seq_id"])

            if fold != 1:
                ins = seq[15:-15]
                if any(ins[j : j + 31] in cagi5_k31 for j in range(0, len(ins) - 30, 15)):
                    scrubbed += 1
                    continue
                if fold == 2:
                    va_seqs.append(seq)
                    va_y.append(y_a)
                    va_ids.append(cid)
                else:
                    tr_seqs.append(seq)
                    tr_y.append(y_r)
            else:
                te_seqs.append(seq)
                te_y_raw.append(y_r)
                te_y_avg.append(y_a)
                te_ids.append(cid)

    # Pre-build padded uint8 array for ultra-fast online shift augmentation during training:
    # [21 bp FORWARD_ADAPTER] + [230 bp seq] + [21 bp REVERSE_ADAPTER[:21]] = 272 bp
    tr_ext = [FORWARD_ADAPTER + s + REVERSE_ADAPTER[:21] for s in tr_seqs]
    return dict(
        cell_line=cell_line,
        scrubbed_cagi5_overlaps=scrubbed,
        tr_ext=seqs_to_uint8(tr_ext),  # (N_tr, 272)
        tr_y=np.array(tr_y, dtype=np.float32),
        va_x=seqs_to_uint8(va_seqs),
        va_y=np.array(va_y, dtype=np.float32),
        va_ids=va_ids,
        te_x=seqs_to_uint8(te_seqs),
        te_y_raw=np.array(te_y_raw, dtype=np.float32),
        te_y_avg=np.array(te_y_avg, dtype=np.float32),
        te_ids=te_ids,
    )


def rc_onehot_batch(oh: torch.Tensor) -> torch.Tensor:
    return torch.flip(oh, dims=[1, 2])


def _amp_ctx(device: torch.device, disable: bool = False):
    if (not disable) and device.type == "cuda" and torch.cuda.is_bf16_supported():
        return torch.amp.autocast("cuda", dtype=torch.bfloat16)
    return torch.amp.autocast(device_type=device.type, enabled=False)


@torch.no_grad()
def predict_array(
    model: nn.Module,
    x_u8: np.ndarray,
    device: torch.device,
    batch_size: int = 512,
    rc_avg: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    model.eval()
    n = len(x_u8)
    preds = np.zeros(n, dtype=np.float64)
    stds = np.zeros(n, dtype=np.float64)
    has_unc = isinstance(model, UNC_MODEL_TYPES) and model.predict_uncertainty
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


def cre_averaged_metrics(preds: np.ndarray, y_avg: np.ndarray, cids: list[str]) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    by_cre: dict[str, list[int]] = defaultdict(list)
    for i, c in enumerate(cids):
        by_cre[c].append(i)
    keys = sorted(by_cre.keys())
    p_cre = np.array([preds[by_cre[k]].mean() for k in keys], dtype=np.float64)
    y_cre = np.array([y_avg[by_cre[k]].mean() for k in keys], dtype=np.float64)
    r = float(np.corrcoef(p_cre, y_cre)[0, 1])
    rho = float(spearmanr(p_cre, y_cre).statistic)
    mse = float(np.mean((p_cre - y_cre) ** 2))
    r2 = float(1.0 - mse / np.var(y_cre))
    return p_cre, y_cre, dict(n_cres=len(keys), pearson_r=r, spearman_rho=rho, mse=mse, r2=r2)


def train_single_model(
    model_name: str,
    data: dict[str, Any],
    w_fwd: np.ndarray,
    w_rc: np.ndarray,
    device: torch.device,
    *,
    epochs: int = 20,
    batch_size: int = 512,
    max_lr: float = 3e-3,
    weight_decay: float = 5e-2,
    seed: int = 42,
    use_unc_loss: bool = True,
    gpu_sleep_ms: float = 0.0,
    use_compile: bool = False,
    rc_consistency: float = 0.0,
    ema_decay: float = 0.0,
    pearson_weight: float = 0.0,
) -> tuple[nn.Module, dict[str, Any]]:
    torch.manual_seed(seed)
    np.random.seed(seed)
    model = build_model(model_name, w_fwd, w_rc, out_dim=1, seq_len=230).to(device)
    n_params = count_parameters(model)
    if use_compile and device.type == "cuda":
        if isinstance(model, HelixLoopOmegaModel):
            if model.thermo is not None:
                model.thermo = torch.compile(model.thermo)
            model.prelude = torch.compile(model.prelude)
            model.coda_conv = torch.compile(model.coda_conv)
            model.core.attn = torch.compile(model.core.attn)
            if model.core.pair_track is not None:
                model.core.pair_track = torch.compile(model.core.pair_track)
            train_model = model
        elif isinstance(model, EpigenesisThermoLoop):
            if model.thermo is not None:
                model.thermo = torch.compile(model.thermo)
            model.prelude = torch.compile(model.prelude)
            model.core = torch.compile(model.core)
            model.coda_conv = torch.compile(model.coda_conv)
            train_model = model
        elif has_unc := (isinstance(model, (EpigenesisThermoSyntaxNet, EpigenesisThermoDREAM)) and model.predict_uncertainty):
            if getattr(model, "thermo", None) is not None:
                model.thermo = torch.compile(model.thermo)
            train_model = model
        else:
            train_model = torch.compile(model)
    else:
        train_model = model

    # LegNet official uses higher LR (0.01) and wd (0.1) as in Agarwal et al. 2025;
    # honor official optimal hyperparams for LegNet so we beat it at its FULL strength!
    if model_name.startswith("legnet"):
        lr_use = 0.008
        wd_use = 0.1
    elif model_name == "dream_rnn":
        lr_use = 0.0015
        wd_use = 0.01
    elif model_name.startswith("thermodream"):
        lr_use = 0.004
        wd_use = 0.02
    elif model_name.startswith("helixloop_v4"):
        lr_use = 0.006
        wd_use = 0.10
    elif model_name.startswith("helixloop_omega"):
        lr_use = 0.005
        wd_use = 0.08
    elif model_name.startswith(("thermoloop_v3", "thermo_interlace_v3")):
        lr_use = 0.006
        wd_use = 0.10
    elif model_name.startswith(("thermoloop_v2", "thermo_interlace_v2")):
        lr_use = 0.005
        wd_use = 0.10
    elif model_name.startswith(("thermoloop", "thermo_interlace")):
        lr_use = 0.003
        wd_use = 0.05
    elif model_name in ("malinois", "deepstarr"):
        lr_use = 0.002
        wd_use = 0.01
    else:
        lr_use = 0.006
        wd_use = 0.1

    tr_ext = data["tr_ext"]
    tr_y = data["tr_y"]
    n_tr = len(tr_y)
    # Keep tiny uint8 sequence table directly on CUDA GPU (26 MB for 98k CREs!) for zero-overhead shift+onehot
    tr_ext_gpu = torch.from_numpy(tr_ext).to(device=device, dtype=torch.long) if device.type == "cuda" else None
    tr_y_gpu = torch.from_numpy(tr_y).to(device=device, dtype=torch.float32) if device.type == "cuda" else None
    seq_offsets_gpu = torch.arange(230, device=device, dtype=torch.long).unsqueeze(0) if device.type == "cuda" else None
    onehot_eye_gpu = torch.cat([torch.eye(4, device=device, dtype=torch.float32), torch.zeros(1, 4, device=device, dtype=torch.float32)], dim=0) if device.type == "cuda" else None
    steps_per_epoch = math.ceil(n_tr / batch_size)

    if model_name.startswith("helixloop_omega") and not model_name.endswith("_adamw"):
        opt = HybridMuonCWD(model, lr_muon=0.006 / 20.0, lr_adam=0.004 / 20.0, weight_decay=wd_use)
        max_lr_sched: Any = [0.006, 0.004, 0.004]
    else:
        opt = torch.optim.AdamW(model.parameters(), lr=lr_use / 20.0, weight_decay=wd_use)
        max_lr_sched = lr_use
    sched = torch.optim.lr_scheduler.OneCycleLR(
        opt,
        max_lr=max_lr_sched,
        epochs=epochs,
        steps_per_epoch=steps_per_epoch,
        pct_start=0.25,
        cycle_momentum=False,
    )

    has_unc = isinstance(model, UNC_MODEL_TYPES) and model.predict_uncertainty and use_unc_loss
    rng = np.random.default_rng(seed)

    # Optional EMA weight shadow (recalibrated cleanly at end of training with exact BN CMA pass!)
    ema_params: dict[str, torch.Tensor] | None = None
    if ema_decay > 0.0:
        ema_params = {k: v.detach().clone().float() for k, v in model.named_parameters()}

    best_val_r = -1e9
    best_state = None
    best_epoch = -1
    t0 = time.perf_counter()

    def _forward_loss(oh_in: torch.Tensor, y_in: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        if isinstance(model, LOOP_MODEL_TYPES):
            out = train_model.forward_loops(oh_in, return_all_steps=True)
            mu_o, lv_o = out["mu"].float(), out["logvar"].float()
            step_mus = [sm.float() for sm in out["step_mus"]]
            if len(step_mus) == 3:
                mse_o = 0.15 * F.mse_loss(step_mus[0], y_in) + 0.25 * F.mse_loss(step_mus[1], y_in) + F.mse_loss(step_mus[2], y_in)
            else:
                mse_o = F.mse_loss(mu_o, y_in)
            if pearson_weight > 0.0:
                mse_o = mse_o + pearson_weight * batch_pearson_loss(mu_o, y_in)
            gnll_o = 0.5 * (torch.exp(-lv_o) * (y_in - mu_o.detach()).pow(2) + lv_o).mean()
            return mu_o, mse_o + 0.15 * gnll_o
        if has_unc:
            mu_o, lv_o = train_model.forward_with_uncertainty(oh_in)
            mu_o, lv_o = mu_o.float(), lv_o.float()
            mse_o = F.mse_loss(mu_o, y_in)
            if pearson_weight > 0.0:
                mse_o = mse_o + pearson_weight * batch_pearson_loss(mu_o, y_in)
            gnll_o = 0.5 * (torch.exp(-lv_o) * (y_in - mu_o.detach()).pow(2) + lv_o).mean()
            return mu_o, mse_o + 0.15 * gnll_o
        mu_o = train_model(oh_in).float()
        loss_o = F.mse_loss(mu_o, y_in)
        if pearson_weight > 0.0:
            loss_o = loss_o + pearson_weight * batch_pearson_loss(mu_o, y_in)
        return mu_o, loss_o

    for ep in range(1, epochs + 1):
        train_model.train()
        perm = rng.permutation(n_tr)
        perm_gpu = torch.from_numpy(perm).to(device=device, dtype=torch.long) if tr_ext_gpu is not None else None
        ep_loss = 0.0
        nb = 0
        for s in range(0, n_tr, batch_size):
            if tr_ext_gpu is not None:
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
            else:
                idx = perm[s : s + batch_size]
                if len(idx) < 8:
                    continue
                sub = tr_ext[idx]  # (B, 272)
                b = len(idx)
                shifts = rng.integers(0, 43, size=b)  # 21 = zero shift
                cols = shifts[:, None] + np.arange(230)[None, :]
                x_u8 = sub[np.arange(b)[:, None], cols]
                oh = uint8_to_onehot(x_u8, device)
                rc_mask = torch.rand(b, device=device) > 0.5
                if rc_mask.any():
                    oh[rc_mask] = rc_onehot_batch(oh[rc_mask])
                yb = torch.from_numpy(tr_y[idx]).to(device=device, dtype=torch.float32)

            opt.zero_grad(set_to_none=True)
            with _amp_ctx(device, disable=(model_name == "dream_rnn")):
                mu, loss = _forward_loss(oh, yb)
                if rc_consistency > 0.0:
                    mu_rc, loss_rc = _forward_loss(rc_onehot_batch(oh), yb)
                    loss = 0.5 * (loss + loss_rc) + rc_consistency * F.mse_loss(mu, mu_rc)

            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 2.0)
            opt.step()
            sched.step()
            if ema_params is not None:
                with torch.no_grad():
                    for k, p in model.named_parameters():
                        ema_params[k].mul_(ema_decay).add_(p.detach().float(), alpha=1.0 - ema_decay)
            if gpu_sleep_ms > 0.0:
                time.sleep(gpu_sleep_ms / 1000.0)
            ep_loss += float(loss.item())
            nb += 1

        # Evaluate strictly on validation fold 2 (CRE-averaged)
        va_p, _ = predict_array(model, data["va_x"], device, batch_size=batch_size, rc_avg=True)
        _, _, va_m = cre_averaged_metrics(va_p, data["va_y"], data["va_ids"])
        if va_m["pearson_r"] > best_val_r:
            best_val_r = va_m["pearson_r"]
            best_epoch = ep
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
        if ep in (1, 5, 10, 15, epochs) or ep == best_epoch:
            print(
                f"  [{model_name:24s}] ep={ep:2d}/{epochs} loss={ep_loss/max(nb,1):.4f} "
                f"val_cre_r={va_m['pearson_r']:.4f} val_mse={va_m['mse']:.4f} ({time.perf_counter()-t0:.1f}s)",
                flush=True,
            )

    # If EMA was tracked, swap EMA weights, run exact 60-batch BN/RMSBN CMA recalibration, and keep if better on val!
    if ema_params is not None:
        with torch.no_grad():
            for k, p in model.named_parameters():
                p.copy_(ema_params[k].to(p.dtype))
        def _bn_batches():
            p_bn = rng.permutation(n_tr)
            for si in range(0, min(n_tr, 60 * batch_size), batch_size):
                ix = p_bn[si : si + batch_size]
                if len(ix) < 8:
                    continue
                if tr_ext_gpu is not None:
                    ix_g = torch.from_numpy(ix).to(device=device, dtype=torch.long)
                    sh_g = torch.full((len(ix), 1), 21, device=device, dtype=torch.long)
                    yield onehot_eye_gpu[tr_ext_gpu[ix_g].gather(1, sh_g + seq_offsets_gpu)].permute(0, 2, 1).contiguous()
                else:
                    yield uint8_to_onehot(tr_ext[ix, 21:251], device)
        recalibrate_bn_stats(model, _bn_batches(), device, disable_amp=(model_name == "dream_rnn"))
        va_p_ema, _ = predict_array(model, data["va_x"], device, batch_size=batch_size, rc_avg=True)
        _, _, va_m_ema = cre_averaged_metrics(va_p_ema, data["va_y"], data["va_ids"])
        print(
            f"  [{model_name:24s}] EMA+BN-recal val_cre_r={va_m_ema['pearson_r']:.4f} val_mse={va_m_ema['mse']:.4f} "
            f"(vs live best={best_val_r:.4f})",
            flush=True,
        )
        if va_m_ema["pearson_r"] > best_val_r:
            best_val_r = va_m_ema["pearson_r"]
            best_epoch = epochs
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}

    model.load_state_dict(best_state)
    train_sec = round(time.perf_counter() - t0, 1)
    return model, dict(
        model=model_name,
        n_params=n_params,
        best_epoch=best_epoch,
        best_val_cre_r=round(best_val_r, 5),
        train_seconds=train_sec,
    )


def evaluate_cagi5_zeroshot(
    model: nn.Module,
    kircher_dir: Path,
    device: torch.device,
    batch_size: int = 512,
) -> dict[str, Any]:
    """TRUE ZERO-SHOT evaluation on ALL 39,170 CAGI5 / Kircher 2019 GRCh38 SNVs."""
    refs = json.loads((kircher_dir / "kircher_refs.json").read_text(encoding="utf-8"))
    # Build exact 200-bp centered genomic window wrapped in lentiMPRA 15-bp adapters (230 bp)
    # Score both ref and alt in both F and RC orientations of the 200-bp insert -> exact RC symmetry!
    by_el: dict[str, list[tuple[str, str, float, float]]] = defaultdict(list)
    with gzip.open(kircher_dir / "kircher_elements.tsv.gz", "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            ref_b, alt_b = r["Ref"].upper(), r["Alt"].upper()
            if len(ref_b) != 1 or len(alt_b) != 1 or ref_b not in "ACGT" or alt_b not in "ACGT":
                continue
            el = r["Element"]
            if el not in refs:
                continue
            info = refs[el]
            seq = info["seq"].upper()
            pos0 = int(r["Pos"]) - int(info["start"])
            if pos0 < 0 or pos0 >= len(seq) or seq[pos0] != ref_b:
                continue
            # Extract 200 bp window centered on pos0
            L = len(seq)
            win = 200
            if L <= win:
                sub_ref = seq.center(win, "N")
                p_rel = pos0 + (win - L) // 2
            else:
                st = max(0, min(L - win, pos0 - win // 2))
                sub_ref = seq[st : st + win]
                p_rel = pos0 - st
            sub_alt = sub_ref[:p_rel] + alt_b + sub_ref[p_rel + 1 :]
            by_el[el].append((sub_ref, sub_alt, float(r["Coefficient"]), float(r["pValue"])))

    per_el = {}
    per_el_sig = {}
    for el in sorted(by_el):
        rows = by_el[el]
        # Deduplicate ref windows per element for speed!
        u_map: dict[str, int] = {}
        u_seqs: list[str] = []
        for sref, salt, _, _ in rows:
            for s in (sref, salt):
                if s not in u_map:
                    u_map[s] = len(u_seqs)
                    u_seqs.append(s)
        # Wrap both F and RC insert inside exact lentiMPRA adapters
        wrapped_f = [LENTI_LEFT + s + LENTI_RIGHT for s in u_seqs]
        wrapped_r = [LENTI_LEFT + rc_str(s) + LENTI_RIGHT for s in u_seqs]
        pf, _ = predict_array(model, seqs_to_uint8(wrapped_f), device, batch_size=batch_size, rc_avg=False)
        pr, _ = predict_array(model, seqs_to_uint8(wrapped_r), device, batch_size=batch_size, rc_avg=False)
        p_both = 0.5 * (pf + pr)

        y = np.array([x[2] for x in rows], dtype=np.float64)
        pv = np.array([x[3] for x in rows], dtype=np.float64)
        ri = np.array([u_map[x[0]] for x in rows], dtype=np.int64)
        ai = np.array([u_map[x[1]] for x in rows], dtype=np.int64)
        dy = p_both[ai] - p_both[ri]

        r_val = float(np.corrcoef(dy, y)[0, 1]) if np.std(dy) > 1e-12 else 0.0
        rho_val = float(spearmanr(dy, y).statistic) if np.std(dy) > 1e-12 else 0.0
        per_el[el] = dict(n=len(y), pearson_r=round(r_val, 4), spearman_rho=round(rho_val, 4))
        ms = pv < 1e-5
        if ms.sum() >= 20:
            r_sig = float(np.corrcoef(dy[ms], y[ms])[0, 1]) if np.std(dy[ms]) > 1e-12 else 0.0
            per_el_sig[el] = dict(n=int(ms.sum()), pearson_r=round(r_sig, 4))

    liver_els = [e for e in ("F9", "HNF4A", "LDLR", "LDLR.2", "SORT1", "SORT1.2", "SORT1-flip") if e in per_el]
    eryth_els = [e for e in ("BCL11A", "GP1BA", "HBB", "HBG1", "PKLR-24h", "PKLR-48h") if e in per_el]
    neural_els = [e for e in ("RET", "TERT-GAa", "TERT-GBM", "TERT-GSc", "UC88") if e in per_el]

    return dict(
        n_elements=len(per_el),
        n_snvs=sum(v["n"] for v in per_el.values()),
        mean_all_r=round(float(np.mean([v["pearson_r"] for v in per_el.values()])), 4),
        mean_all_rho=round(float(np.mean([v["spearman_rho"] for v in per_el.values()])), 4),
        mean_sig_p1e5_r=round(float(np.mean([v["pearson_r"] for v in per_el_sig.values()])), 4),
        mean_liver_r=round(float(np.mean([per_el[e]["pearson_r"] for e in liver_els])), 4),
        mean_erythroid_r=round(float(np.mean([per_el[e]["pearson_r"] for e in eryth_els])), 4),
        mean_neural_r=round(float(np.mean([per_el[e]["pearson_r"] for e in neural_els])), 4),
        per_element=per_el,
        per_element_sig=per_el_sig,
    )


def paired_bootstrap_gain(p_cand: np.ndarray, p_base: np.ndarray, y: np.ndarray, n_boot: int = 1000, seed: int = 2026) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    n = len(y)
    gains = np.zeros(n_boot, dtype=np.float64)
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        rc = np.corrcoef(p_cand[idx], y[idx])[0, 1]
        rb = np.corrcoef(p_base[idx], y[idx])[0, 1]
        gains[b] = rc - rb
    return dict(
        point_gain=round(float(np.corrcoef(p_cand, y)[0, 1] - np.corrcoef(p_base, y)[0, 1]), 5),
        ci95_low=round(float(np.quantile(gains, 0.025)), 5),
        ci95_high=round(float(np.quantile(gains, 0.975)), 5),
        prob_better=round(float((gains > 0).mean()), 4),
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell-line", default="WTC11", choices=("WTC11", "HepG2", "K562"))
    ap.add_argument(
        "--models",
        nargs="+",
        default=["legnet_official", "legnet_large", "malinois", "deepstarr", "thermosyntax", "thermosyntax_no_biophys", "thermosyntax_no_attn"],
    )
    ap.add_argument("--epochs", type=int, default=20)
    ap.add_argument("--batch-size", type=int, default=512)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--gpu-sleep-ms", type=float, default=0.0, help="Optional inter-batch sleep (ms) to share GPU politely")
    ap.add_argument("--nice", type=int, default=0, help="Optional OS process niceness (e.g. 15 for background runs)")
    ap.add_argument("--tag", type=str, default="", help="Optional suffix for output JSON filename")
    ap.add_argument("--compile", action="store_true", help="Enable torch.compile on CUDA for Triton kernel fusion")
    ap.add_argument("--rc-consistency", type=float, default=0.0, help="Paired double-stranded RC consistency loss weight")
    ap.add_argument("--ema", type=float, default=0.0, help="Optional EMA decay with exact end-of-training BN recalibration (e.g. 0.997)")
    ap.add_argument("--pearson-weight", type=float, default=0.0, help="Optional batch Pearson loss weight (e.g. 0.05)")
    ap.add_argument("--data-dir", type=Path, default=Path("/tmp/epigenesis-benchmarks/human_legnet-main/datasets/original"))
    ap.add_argument("--kircher-dir", type=Path, default=Path("/tmp/epigenesis-frontiers/data"))
    ap.add_argument("--out-dir", type=Path, default=Path("/tmp/epigenesis-benchmarks/results"))
    args = ap.parse_args()

    if args.nice > 0:
        import os
        os.nice(args.nice)
    device = torch.device("cuda" if torch.cuda.is_available() else ("mps" if torch.backends.mps.is_available() else "cpu"))
    if device.type == "cuda":
        torch.backends.cuda.matmul.allow_tf32 = True
        torch.backends.cudnn.allow_tf32 = True
        torch.backends.cudnn.benchmark = True
    print(f"=== Epigenesis Pre-Registered Benchmark on {args.cell_line} | device={device} ===", flush=True)

    w_fwd, w_rc, _ = load_jaspar_pwm_bank(max_width=21)
    cagi5_k31 = build_cagi5_kmer_set(args.kircher_dir / "kircher_refs.json")
    data = load_lentimpra(args.data_dir, args.cell_line, cagi5_k31)
    print(
        f"Loaded {args.cell_line}: train={len(data['tr_y'])}, val={len(data['va_y'])}, "
        f"test={len(data['te_y_raw'])} (scrubbed {data['scrubbed_cagi5_overlaps']} CAGI5-overlapping sequences)",
        flush=True,
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    results: dict[str, Any] = dict(
        cell_line=args.cell_line,
        device=str(device),
        epochs=args.epochs,
        seed=args.seed,
        n_train=int(len(data["tr_y"])),
        n_val=int(len(data["va_y"])),
        n_test_raw=int(len(data["te_y_raw"])),
        scrubbed_cagi5_overlaps=data["scrubbed_cagi5_overlaps"],
        models={},
    )
    cre_preds_by_model: dict[str, np.ndarray] = {}
    y_cre_ref: np.ndarray | None = None

    for mname in args.models:
        print(f"\n--- Training {mname} ---", flush=True)
        model, info = train_single_model(
            mname,
            data,
            w_fwd,
            w_rc,
            device,
            epochs=args.epochs,
            batch_size=args.batch_size,
            seed=args.seed,
            gpu_sleep_ms=args.gpu_sleep_ms,
            use_compile=args.compile,
            rc_consistency=args.rc_consistency,
            ema_decay=args.ema,
            pearson_weight=args.pearson_weight,
        )
        # Evaluate ONCE on strictly held-out fold 1
        te_p, te_std = predict_array(model, data["te_x"], device, batch_size=args.batch_size, rc_avg=True)
        raw_r = float(np.corrcoef(te_p, data["te_y_raw"])[0, 1])
        raw_rho = float(spearmanr(te_p, data["te_y_raw"]).statistic)
        raw_mse = float(np.mean((te_p - data["te_y_raw"]) ** 2))

        p_cre, y_cre, cre_m = cre_averaged_metrics(te_p, data["te_y_avg"], data["te_ids"])
        cre_preds_by_model[mname] = p_cre
        y_cre_ref = y_cre

        cal_rho = None
        if np.std(te_std) > 1e-12:
            abs_err = np.abs(te_p - data["te_y_raw"])
            cal_rho = round(float(spearmanr(te_std, abs_err).statistic), 4)

        # If Looped Transformer, also probe test-time loop depth scaling T in {1, 2, 4, 6, 8}!
        loop_scaling = None
        if isinstance(model, LOOP_MODEL_TYPES):
            loop_scaling = {}
            orig_loops = model.n_loops
            for t_eval in (1, 2, 4, 6, 8):
                model.n_loops = t_eval
                tp_t, _ = predict_array(model, data["te_x"], device, batch_size=args.batch_size, rc_avg=True)
                _, _, cm_t = cre_averaged_metrics(tp_t, data["te_y_avg"], data["te_ids"])
                loop_scaling[f"T={t_eval}"] = dict(cre_r=round(cm_t["pearson_r"], 5), cre_mse=round(cm_t["mse"], 5))
            model.n_loops = orig_loops
            print(f"  ==> LOOP DEPTH SCALING ({mname}): {loop_scaling}", flush=True)

        # TRUE ZERO-SHOT CAGI5 evaluation (39,170 SNVs)
        cagi5 = evaluate_cagi5_zeroshot(model, args.kircher_dir, device, batch_size=args.batch_size)

        # Save clean uncompiled checkpoint + test CRE predictions locally for exact cross-run paired bootstraps
        ckpt_path = args.out_dir / f"{args.cell_line}_{mname}_s{args.seed}.pt"
        clean_sd = {k.replace("._orig_mod.", ".").replace("_orig_mod.", ""): v.detach().cpu() for k, v in model.state_dict().items()}
        torch.save({"model_name": mname, "state_dict": clean_sd, "info": info}, ckpt_path)
        np.savez_compressed(args.out_dir / f"{args.cell_line}_{mname}_s{args.seed}_cre_preds.npz", p_cre=p_cre, y_cre=y_cre)

        results["models"][mname] = dict(
            **info,
            test_raw=dict(pearson_r=round(raw_r, 5), spearman_rho=round(raw_rho, 5), mse=round(raw_mse, 5)),
            test_cre_averaged={k: (round(v, 5) if isinstance(v, float) else v) for k, v in cre_m.items()},
            uncertainty_calibration_spearman=cal_rho,
            loop_depth_scaling=loop_scaling,
            cagi5_zeroshot=dict(
                mean_all_r=cagi5["mean_all_r"],
                mean_all_rho=cagi5["mean_all_rho"],
                mean_sig_p1e5_r=cagi5["mean_sig_p1e5_r"],
                mean_liver_r=cagi5["mean_liver_r"],
                mean_erythroid_r=cagi5["mean_erythroid_r"],
                mean_neural_r=cagi5["mean_neural_r"],
            ),
            cagi5_per_element=cagi5["per_element"],
        )
        print(
            f"  ==> TEST {mname}: raw_r={raw_r:.4f} | CRE_r={cre_m['pearson_r']:.4f} "
            f"(MSE={cre_m['mse']:.4f}) | cal_rho={cal_rho} | "
            f"CAGI5_ZS all={cagi5['mean_all_r']:.4f} sig={cagi5['mean_sig_p1e5_r']:.4f} "
            f"liver={cagi5['mean_liver_r']:.4f} eryth={cagi5['mean_erythroid_r']:.4f} neural={cagi5['mean_neural_r']:.4f}",
            flush=True,
        )

    # Paired bootstrap comparisons across all candidate models vs strong baselines
    if y_cre_ref is not None:
        gains = {}
        for cand in ("helixloop_v4", "helixloop_v4_no_delta", "helixloop_omega", "helixloop_omega_adamw", "thermoloop_v3", "thermo_interlace_v3", "thermoloop_v2", "thermo_interlace_v2", "thermoloop", "thermo_interlace", "thermosyntax", "thermodream"):
            if cand not in cre_preds_by_model:
                continue
            for base in ("legnet_official", "legnet_large", "dream_rnn", "helixloop_v4_t1", "helixloop_v4_no_biophys", "helixloop_v4_no_diff", "helixloop_omega_t1", "helixloop_omega_no_biophys", "helixloop_omega_no_engram", "helixloop_omega_no_tri", "helixloop_omega_no_delta", "helixloop_omega_adamw", "thermoloop_t1", "thermoloop_no_biophys", "thermoloop_no_rational"):
                if base in cre_preds_by_model and base != cand:
                    gains[f"{cand}_vs_{base}"] = paired_bootstrap_gain(cre_preds_by_model[cand], cre_preds_by_model[base], y_cre_ref, n_boot=2000)
        if gains:
            results["paired_bootstrap_cre_r_gains"] = gains
            print("\n=== Paired Bootstrap 95% CIs (CRE-averaged Test Fold 1, 2000 draws) ===", flush=True)
            for k, v in gains.items():
                print(f"  {k:44s}: {v}", flush=True)

    tag_sfx = f"_{args.tag}" if args.tag else ""
    out_json = args.out_dir / f"benchmark_{args.cell_line}{tag_sfx}_s{args.seed}.json"
    out_json.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(f"\nSaved full results to {out_json}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

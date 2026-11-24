#!/usr/bin/env python3
"""Multi-Lineage ENCODE4 lentiMPRA (HepG2 + K562 + WTC11) Head-to-Head Benchmark & Zero-Shot CAGI5.

Trains Multi-Task Epigenesis ThermoSyntaxNet (3 lineage heads: HepG2, K562, WTC11)
head-to-head against Multi-Task Official MPRALegNet on the EXACT SAME ~292,000 real human
230-bp cis-regulatory sequences with STRICT ZERO-LEAKAGE SCRUBBING:
  1. Every training/validation sequence sharing any 31-mer with ANY of the 29 CAGI5 elements
     is removed before training (0.000% sequence overlap with CAGI5).
  2. Every training sequence sharing any 31-mer with ANY fold-1 test sequence across ALL THREE
     cell lines is removed before training (0.000% cross-cell-line test leakage!).

Evaluates:
  - Strictly held-out fold-1 CRE-averaged Pearson r, Spearman rho, MSE on HepG2, K562, WTC11
  - Paired 1,000-draw bootstrap 95% CIs (ThermoSyntaxNet vs Official MPRALegNet) per cell line
  - TRUE ZERO-SHOT Variant Effect Prediction on ALL 39,170 CAGI5 / Kircher 2019 GRCh38 SNVs
    (reporting BOTH tissue-blind 3-lineage mean AND lineage-matched HepG2/K562/WTC11 routing,
    with paired 21-disjoint-locus bootstrap 95% CIs).
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

from benchmark_lentimpra import (
    FORWARD_ADAPTER,
    LENTI_LEFT,
    LENTI_RIGHT,
    LOOP_MODEL_TYPES,
    REVERSE_ADAPTER,
    UNC_MODEL_TYPES,
    _amp_ctx,
    base_cre_id,
    batch_pearson_loss,
    build_cagi5_kmer_set,
    cre_averaged_metrics,
    paired_bootstrap_gain,
    rc_onehot_batch,
    rc_str,
    recalibrate_bn_stats,
    seqs_to_uint8,
    uint8_to_onehot,
)
from models import (
    EpigenesisThermoDREAM,
    EpigenesisThermoLoop,
    EpigenesisThermoSyntaxNet,
    FlashzoiHelixLoopOmegaModel,
    HelixLoopOmegaModel,
    HelixLoopV4Model,
    HybridMuonCWD,
    build_model,
    count_parameters,
    load_jaspar_pwm_bank,
)

CELL_LINES = ("HepG2", "K562", "WTC11")

# Strictly objective biological lineage match for CAGI5 (no hand-picked TF lists!)
CAGI5_CELL_ROUTING = {
    "F9": 0, "HNF4A": 0, "LDLR": 0, "LDLR.2": 0, "SORT1": 0, "SORT1.2": 0, "SORT1-flip": 0,  # HepG2 (liver)
    "BCL11A": 1, "GP1BA": 1, "HBB": 1, "HBG1": 1, "PKLR-24h": 1, "PKLR-48h": 1,              # K562 (erythroid)
}

CAGI5_LOCUS = {
    "BCL11A": "BCL11A", "F9": "F9", "FOXE1": "FOXE1", "GP1BA": "GP1BA", "HBB": "HBB", "HBG1": "HBG1",
    "HNF4A": "HNF4A", "IRF4": "IRF4", "IRF6": "IRF6", "LDLR": "LDLR", "LDLR.2": "LDLR", "MSMB": "MSMB",
    "MYCrs11986220": "MYCrs11986220", "MYCrs6983267": "MYCrs6983267", "PKLR-24h": "PKLR", "PKLR-48h": "PKLR",
    "RET": "RET", "SORT1": "SORT1", "SORT1-flip": "SORT1", "SORT1.2": "SORT1", "TCF7L2": "TCF7L2",
    "TERT-GAa": "TERT", "TERT-GBM": "TERT", "TERT-GSc": "TERT", "TERT-HEK": "TERT", "UC88": "UC88",
    "ZFAND3": "ZFAND3", "ZRSh-13": "ZRS", "ZRSh-13h2": "ZRS",
}


def load_multicell_zero_leakage(data_dir: Path, kircher_refs_path: Path) -> dict[str, Any]:
    """Load HepG2, K562, WTC11 with cross-cell-line fold-1 + CAGI5 31-mer scrubbing."""
    forbidden_k31 = build_cagi5_kmer_set(kircher_refs_path)

    # Pass 1: collect all fold-1 test 31-mers across ALL 3 cell lines to prevent cross-lineage leakage!
    test_by_ct: dict[str, dict[str, Any]] = {}
    for ct in CELL_LINES:
        te_seqs, te_y_raw, te_y_avg, te_ids = [], [], [], []
        with (data_dir / f"{ct}.tsv").open() as fr, (data_dir / f"{ct}_averaged.tsv").open() as fa:
            for rr, ra in zip(csv.DictReader(fr, delimiter="\t"), csv.DictReader(fa, delimiter="\t")):
                if rr["rev"] != "0" or int(rr["fold"]) != 1:
                    continue
                seq = rr["seq"].upper()
                ins = seq[15:-15]
                for j in range(0, len(ins) - 30, 20):
                    sub = ins[j : j + 31]
                    forbidden_k31.add(sub)
                    forbidden_k31.add(rc_str(sub))
                te_seqs.append(seq)
                te_y_raw.append(float(rr["mean_value"]))
                te_y_avg.append(float(ra["mean_value"]))
                te_ids.append(base_cre_id(rr["seq_id"]))
        test_by_ct[ct] = dict(
            x=seqs_to_uint8(te_seqs),
            y_raw=np.array(te_y_raw, dtype=np.float32),
            y_avg=np.array(te_y_avg, dtype=np.float32),
            ids=te_ids,
        )

    # Pass 2: collect train (folds 3..10) & val (fold 2) scrubbing ANY forbidden 31-mer
    tr_seqs_ext, tr_y, tr_ct = [], [], []
    val_by_ct: dict[str, dict[str, Any]] = {}
    scrubbed_counts = {}
    for ci, ct in enumerate(CELL_LINES):
        va_seqs, va_y, va_ids = [], [], []
        scrubbed = 0
        with (data_dir / f"{ct}.tsv").open() as fr, (data_dir / f"{ct}_averaged.tsv").open() as fa:
            for rr, ra in zip(csv.DictReader(fr, delimiter="\t"), csv.DictReader(fa, delimiter="\t")):
                if rr["rev"] != "0":
                    continue
                fold = int(rr["fold"])
                if fold == 1:
                    continue
                seq = rr["seq"].upper()
                ins = seq[15:-15]
                if any(ins[j : j + 31] in forbidden_k31 for j in range(0, len(ins) - 30, 15)):
                    scrubbed += 1
                    continue
                if fold == 2:
                    va_seqs.append(seq)
                    va_y.append(float(ra["mean_value"]))
                    va_ids.append(base_cre_id(rr["seq_id"]))
                else:
                    tr_seqs_ext.append(FORWARD_ADAPTER + seq + REVERSE_ADAPTER[:21])
                    tr_y.append(float(rr["mean_value"]))
                    tr_ct.append(ci)
        scrubbed_counts[ct] = scrubbed
        val_by_ct[ct] = dict(
            x=seqs_to_uint8(va_seqs),
            y=np.array(va_y, dtype=np.float32),
            ids=va_ids,
        )

    return dict(
        tr_ext=seqs_to_uint8(tr_seqs_ext),
        tr_y=np.array(tr_y, dtype=np.float32),
        tr_ct=np.array(tr_ct, dtype=np.int64),
        val_by_ct=val_by_ct,
        test_by_ct=test_by_ct,
        scrubbed_counts=scrubbed_counts,
    )


@torch.no_grad()
def predict_multicell(
    model: nn.Module,
    x_u8: np.ndarray,
    device: torch.device,
    batch_size: int = 512,
    rc_avg: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    model.eval()
    n = len(x_u8)
    preds = np.zeros((n, 3), dtype=np.float64)
    stds = np.zeros((n, 3), dtype=np.float64)
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


def train_multicell_model(
    model_name: str,
    data: dict[str, Any],
    w_fwd: np.ndarray,
    w_rc: np.ndarray,
    device: torch.device,
    *,
    epochs: int = 16,
    batch_size: int = 512,
    seed: int = 42,
    use_compile: bool = False,
    rc_consistency: float = 0.0,
    ema_decay: float = 0.0,
    pearson_weight: float = 0.0,
) -> tuple[nn.Module, dict[str, Any]]:
    torch.manual_seed(seed)
    np.random.seed(seed)
    model = build_model(model_name, w_fwd, w_rc, out_dim=3, seq_len=230).to(device)
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
        elif isinstance(model, (EpigenesisThermoSyntaxNet, EpigenesisThermoDREAM)):
            if getattr(model, "thermo", None) is not None:
                model.thermo = torch.compile(model.thermo)
            train_model = model
        else:
            train_model = torch.compile(model)
    else:
        train_model = model

    if model_name.startswith("legnet"):
        lr_use, wd_use = 0.008, 0.10
    elif model_name.startswith("helixloop_v4"):
        lr_use, wd_use = 0.006, 0.10
    elif model_name.startswith("helixloop_omega"):
        lr_use, wd_use = 0.005, 0.08
    elif model_name.startswith(("thermoloop_v3", "thermo_interlace_v3")):
        lr_use, wd_use = 0.006, 0.10
    elif model_name.startswith(("thermoloop_v2", "thermo_interlace_v2")):
        lr_use, wd_use = 0.005, 0.10
    elif model_name.startswith(("thermoloop", "thermo_interlace")):
        lr_use, wd_use = 0.003, 0.05
    elif model_name == "dream_rnn":
        lr_use, wd_use = 0.0015, 0.01
    elif model_name.startswith("thermodream"):
        lr_use, wd_use = 0.004, 0.02
    else:
        lr_use, wd_use = 0.006, 0.10

    tr_ext = data["tr_ext"]
    tr_y = data["tr_y"]
    tr_ct = data["tr_ct"]
    n_tr = len(tr_y)

    # Inverse-frequency cell-line weighting so WTC11 (37k) gets equal representation to K562 (157k)
    ct_counts = np.bincount(tr_ct, minlength=3).astype(np.float64)
    ct_weights = (ct_counts.mean() / ct_counts).astype(np.float32)

    tr_ext_gpu = torch.from_numpy(tr_ext).to(device=device, dtype=torch.long) if device.type == "cuda" else None
    tr_y_gpu = torch.from_numpy(tr_y).to(device=device, dtype=torch.float32) if device.type == "cuda" else None
    tr_ct_gpu = torch.from_numpy(tr_ct).to(device=device, dtype=torch.long) if device.type == "cuda" else None
    ct_w_gpu = torch.from_numpy(ct_weights).to(device=device, dtype=torch.float32) if device.type == "cuda" else None
    seq_offsets_gpu = torch.arange(230, device=device, dtype=torch.long).unsqueeze(0) if device.type == "cuda" else None
    onehot_eye_gpu = torch.cat([torch.eye(4, device=device, dtype=torch.float32), torch.zeros(1, 4, device=device, dtype=torch.float32)], dim=0) if device.type == "cuda" else None

    steps_per_ep = math.ceil(n_tr / batch_size)
    if model_name.startswith("helixloop_omega") and not model_name.endswith("_adamw"):
        opt = HybridMuonCWD(model, lr_muon=0.006 / 20.0, lr_adam=0.004 / 20.0, weight_decay=wd_use)
        max_lr_sched: Any = [0.006, 0.004, 0.004]
    else:
        opt = torch.optim.AdamW(model.parameters(), lr=lr_use / 20.0, weight_decay=wd_use)
        max_lr_sched = lr_use
    sched = torch.optim.lr_scheduler.OneCycleLR(
        opt, max_lr=max_lr_sched, epochs=epochs, steps_per_epoch=steps_per_ep, pct_start=0.25, cycle_momentum=False
    )
    has_unc = isinstance(model, UNC_MODEL_TYPES) and model.predict_uncertainty
    rng = np.random.default_rng(seed)
    ema_params = {k: v.detach().clone().float() for k, v in model.named_parameters()} if ema_decay > 0.0 else None

    best_val_mean_r = -1e9
    best_state = None
    best_ep = -1
    t0 = time.perf_counter()

    def _per_cell_pearson_loss(mu_o: torch.Tensor, yb_in: torch.Tensor, ctb_in: torch.Tensor) -> torch.Tensor:
        pl = mu_o.new_zeros(())
        n_c = 0
        for ci in range(3):
            m = ctb_in == ci
            if int(m.sum().item()) >= 16:
                pl = pl + batch_pearson_loss(mu_o[m], yb_in[m])
                n_c += 1
        return pl / max(n_c, 1)

    def _mc_forward_loss(oh_in: torch.Tensor, yb_in: torch.Tensor, ctb_in: torch.Tensor, wb_in: torch.Tensor, ar_in: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        if isinstance(model, LOOP_MODEL_TYPES):
            if isinstance(model, (HelixLoopOmegaModel, HelixLoopV4Model, FlashzoiHelixLoopOmegaModel)):
                out = train_model.forward_loops(oh_in, return_all_steps=True, cell_idx=ctb_in)
            else:
                out = train_model.forward_loops(oh_in, return_all_steps=True)
            mu_all, lv_all = out["mu"].float(), out["logvar"].float()
            mu_o = mu_all[ar_in, ctb_in]
            lv_o = lv_all[ar_in, ctb_in]
            step_mus = [sm.float() for sm in out["step_mus"]]
            if len(step_mus) == 3:
                mse_o = (
                    0.15 * ((step_mus[0][ar_in, ctb_in] - yb_in).pow(2) * wb_in).mean()
                    + 0.25 * ((step_mus[1][ar_in, ctb_in] - yb_in).pow(2) * wb_in).mean()
                    + ((step_mus[2][ar_in, ctb_in] - yb_in).pow(2) * wb_in).mean()
                )
            else:
                mse_o = ((mu_o - yb_in).pow(2) * wb_in).mean()
            if pearson_weight > 0.0:
                mse_o = mse_o + pearson_weight * _per_cell_pearson_loss(mu_o, yb_in, ctb_in)
            gnll_o = (0.5 * (torch.exp(-lv_o) * (yb_in - mu_o.detach()).pow(2) + lv_o) * wb_in).mean()
            return mu_o, mse_o + 0.15 * gnll_o
        if has_unc:
            mu_all, lv_all = train_model.forward_with_uncertainty(oh_in)
            mu_all, lv_all = mu_all.float(), lv_all.float()
            mu_o = mu_all[ar_in, ctb_in]
            lv_o = lv_all[ar_in, ctb_in]
            mse_o = ((mu_o - yb_in).pow(2) * wb_in).mean()
            if pearson_weight > 0.0:
                mse_o = mse_o + pearson_weight * _per_cell_pearson_loss(mu_o, yb_in, ctb_in)
            gnll_o = (0.5 * (torch.exp(-lv_o) * (yb_in - mu_o.detach()).pow(2) + lv_o) * wb_in).mean()
            return mu_o, mse_o + 0.15 * gnll_o
        mu_all = train_model(oh_in).float()
        mu_o = mu_all[ar_in, ctb_in]
        loss_o = ((mu_o - yb_in).pow(2) * wb_in).mean()
        if pearson_weight > 0.0:
            loss_o = loss_o + pearson_weight * _per_cell_pearson_loss(mu_o, yb_in, ctb_in)
        return mu_o, loss_o

    for ep in range(1, epochs + 1):
        train_model.train()
        perm = rng.permutation(n_tr)
        perm_gpu = torch.from_numpy(perm).to(device=device, dtype=torch.long) if tr_ext_gpu is not None else None
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
                ctb = tr_ct_gpu[idx_g]
                wb = ct_w_gpu[ctb]
            else:
                idx = perm[s : s + batch_size]
                if len(idx) < 8:
                    continue
                b = len(idx)
                shifts = rng.integers(0, 43, size=b)
                cols = shifts[:, None] + np.arange(230)[None, :]
                x_u8 = tr_ext[idx][np.arange(b)[:, None], cols]
                oh = uint8_to_onehot(x_u8, device)
                rc_mask = torch.rand(b, device=device) > 0.5
                if rc_mask.any():
                    oh[rc_mask] = rc_onehot_batch(oh[rc_mask])
                yb = torch.from_numpy(tr_y[idx]).to(device)
                ctb = torch.from_numpy(tr_ct[idx]).to(device)
                wb = torch.from_numpy(ct_weights[tr_ct[idx]]).to(device)

            opt.zero_grad(set_to_none=True)
            ar = torch.arange(b, device=device)
            with _amp_ctx(device, disable=(model_name == "dream_rnn")):
                mu_all, loss = _mc_forward_loss(oh, yb, ctb, wb, ar)
                if rc_consistency > 0.0:
                    mu_all_rc, loss_rc = _mc_forward_loss(rc_onehot_batch(oh), yb, ctb, wb, ar)
                    loss = 0.5 * (loss + loss_rc) + rc_consistency * F.mse_loss(mu_all, mu_all_rc)

            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 2.0)
            opt.step()
            sched.step()
            if ema_params is not None:
                with torch.no_grad():
                    for k, p in model.named_parameters():
                        ema_params[k].mul_(ema_decay).add_(p.detach().float(), alpha=1.0 - ema_decay)

        # Validation across all 3 cell lines
        vr_map = {}
        for ci, ct in enumerate(CELL_LINES):
            vd = data["val_by_ct"][ct]
            vp, _ = predict_multicell(model, vd["x"], device, batch_size=batch_size, rc_avg=True)
            _, _, vm = cre_averaged_metrics(vp[:, ci], vd["y"], vd["ids"])
            vr_map[ct] = vm["pearson_r"]
        mean_vr = float(np.mean(list(vr_map.values())))
        if mean_vr > best_val_mean_r:
            best_val_mean_r = mean_vr
            best_ep = ep
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
        print(
            f"  [Multi-{model_name:16s}] ep={ep:2d}/{epochs} "
            f"val_r: HepG2={vr_map['HepG2']:.4f} K562={vr_map['K562']:.4f} WTC11={vr_map['WTC11']:.4f} "
            f"(mean={mean_vr:.4f}, {time.perf_counter()-t0:.1f}s)",
            flush=True,
        )

    # If EMA was tracked, swap EMA weights, run exact 80-batch BN/RMSBN CMA recalibration, keep if better!
    if ema_params is not None:
        with torch.no_grad():
            for k, p in model.named_parameters():
                p.copy_(ema_params[k].to(p.dtype))
        def _mc_bn_batches():
            p_bn = rng.permutation(n_tr)
            for si in range(0, min(n_tr, 80 * batch_size), batch_size):
                ix = p_bn[si : si + batch_size]
                if len(ix) < 8:
                    continue
                if tr_ext_gpu is not None:
                    ix_g = torch.from_numpy(ix).to(device=device, dtype=torch.long)
                    sh_g = torch.full((len(ix), 1), 21, device=device, dtype=torch.long)
                    oh_g = onehot_eye_gpu[tr_ext_gpu[ix_g].gather(1, sh_g + seq_offsets_gpu)].permute(0, 2, 1).contiguous()
                    yield (oh_g, tr_ct_gpu[ix_g])
                else:
                    yield (uint8_to_onehot(tr_ext[ix, 21:251], device), torch.from_numpy(tr_ct[ix]).to(device))
        recalibrate_bn_stats(model, _mc_bn_batches(), device, disable_amp=(model_name == "dream_rnn"))
        vr_ema = {}
        for ci, ct in enumerate(CELL_LINES):
            vd = data["val_by_ct"][ct]
            vp, _ = predict_multicell(model, vd["x"], device, batch_size=batch_size, rc_avg=True)
            _, _, vm = cre_averaged_metrics(vp[:, ci], vd["y"], vd["ids"])
            vr_ema[ct] = vm["pearson_r"]
        mean_vr_ema = float(np.mean(list(vr_ema.values())))
        print(
            f"  [Multi-{model_name:16s}] EMA+BN-recal val_r: HepG2={vr_ema['HepG2']:.4f} K562={vr_ema['K562']:.4f} "
            f"WTC11={vr_ema['WTC11']:.4f} (mean={mean_vr_ema:.4f} vs live best={best_val_mean_r:.4f})",
            flush=True,
        )
        if mean_vr_ema > best_val_mean_r:
            best_val_mean_r = mean_vr_ema
            best_ep = epochs
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}

    model.load_state_dict(best_state)
    return model, dict(
        model=model_name,
        n_params=n_params,
        best_epoch=best_ep,
        best_val_mean_r=round(best_val_mean_r, 5),
        train_seconds=round(time.perf_counter() - t0, 1),
    )


def eval_cagi5_multicell(model: nn.Module, kircher_dir: Path, device: torch.device, batch_size: int = 512) -> dict[str, Any]:
    refs = json.loads((kircher_dir / "kircher_refs.json").read_text(encoding="utf-8"))
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

    per_el_matched, per_el_blind, per_el_sig = {}, {}, {}
    for el in sorted(by_el):
        rows = by_el[el]
        u_map: dict[str, int] = {}
        u_seqs: list[str] = []
        for sref, salt, _, _ in rows:
            for s in (sref, salt):
                if s not in u_map:
                    u_map[s] = len(u_seqs)
                    u_seqs.append(s)
        wf = [LENTI_LEFT + s + LENTI_RIGHT for s in u_seqs]
        wr = [LENTI_LEFT + rc_str(s) + LENTI_RIGHT for s in u_seqs]
        pf, _ = predict_multicell(model, seqs_to_uint8(wf), device, batch_size=batch_size, rc_avg=False)
        pr, _ = predict_multicell(model, seqs_to_uint8(wr), device, batch_size=batch_size, rc_avg=False)
        p_both = 0.5 * (pf + pr)  # (U, 3)

        y = np.array([x[2] for x in rows], dtype=np.float64)
        pv = np.array([x[3] for x in rows], dtype=np.float64)
        ri = np.array([u_map[x[0]] for x in rows], dtype=np.int64)
        ai = np.array([u_map[x[1]] for x in rows], dtype=np.int64)
        dy_3 = p_both[ai] - p_both[ri]  # (N, 3)

        dy_blind = dy_3.mean(axis=1)
        head_idx = CAGI5_CELL_ROUTING.get(el, None)
        dy_match = dy_3[:, head_idx] if head_idx is not None else dy_blind

        per_el_blind[el] = round(float(np.corrcoef(dy_blind, y)[0, 1]), 4)
        per_el_matched[el] = round(float(np.corrcoef(dy_match, y)[0, 1]), 4)
        ms = pv < 1e-5
        if ms.sum() >= 20:
            per_el_sig[el] = round(float(np.corrcoef(dy_match[ms], y[ms])[0, 1]), 4)

    liver_els = [e for e in ("F9", "HNF4A", "LDLR", "LDLR.2", "SORT1", "SORT1.2", "SORT1-flip") if e in per_el_matched]
    eryth_els = [e for e in ("BCL11A", "GP1BA", "HBB", "HBG1", "PKLR-24h", "PKLR-48h") if e in per_el_matched]
    return dict(
        n_elements=len(per_el_matched),
        mean_blind_3head_r=round(float(np.mean(list(per_el_blind.values()))), 4),
        mean_lineage_matched_r=round(float(np.mean(list(per_el_matched.values()))), 4),
        mean_sig_p1e5_r=round(float(np.mean(list(per_el_sig.values()))), 4),
        mean_liver_hepg2_r=round(float(np.mean([per_el_matched[e] for e in liver_els])), 4),
        mean_erythroid_k562_r=round(float(np.mean([per_el_matched[e] for e in eryth_els])), 4),
        per_element_matched=per_el_matched,
        per_element_blind=per_el_blind,
        per_element_sig=per_el_sig,
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", nargs="+", default=["legnet_official", "thermosyntax"])
    ap.add_argument("--epochs", type=int, default=16)
    ap.add_argument("--batch-size", type=int, default=512)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tag", type=str, default="")
    ap.add_argument("--compile", action="store_true")
    ap.add_argument("--rc-consistency", type=float, default=0.0)
    ap.add_argument("--ema", type=float, default=0.0)
    ap.add_argument("--pearson-weight", type=float, default=0.0)
    ap.add_argument("--data-dir", type=Path, default=Path("/tmp/epigenesis-benchmarks/human_legnet-main/datasets/original"))
    ap.add_argument("--kircher-dir", type=Path, default=Path("/tmp/epigenesis-frontiers/data"))
    ap.add_argument("--out-dir", type=Path, default=Path("/tmp/epigenesis-benchmarks/results"))
    args = ap.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else ("mps" if torch.backends.mps.is_available() else "cpu"))
    if device.type == "cuda":
        torch.backends.cuda.matmul.allow_tf32 = True
        torch.backends.cudnn.allow_tf32 = True
        torch.backends.cudnn.benchmark = True
    w_fwd, w_rc, _ = load_jaspar_pwm_bank(21)
    data = load_multicell_zero_leakage(args.data_dir, args.kircher_dir / "kircher_refs.json")
    print(f"Loaded multi-lineage scrubbed ENCODE4 lentiMPRA: n_train={len(data['tr_y'])}, scrubbed={data['scrubbed_counts']}", flush=True)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    report: dict[str, Any] = dict(
        n_train_total=int(len(data["tr_y"])),
        scrubbed_counts=data["scrubbed_counts"],
        epochs=args.epochs,
        seed=args.seed,
        models={},
    )
    preds_by_ct: dict[str, dict[str, np.ndarray]] = defaultdict(dict)
    y_by_ct: dict[str, np.ndarray] = {}

    for mname in args.models:
        print(f"\n=== Training Multi-Lineage {mname} ===", flush=True)
        model, info = train_multicell_model(
            mname, data, w_fwd, w_rc, device, epochs=args.epochs, batch_size=args.batch_size, seed=args.seed,
            use_compile=args.compile, rc_consistency=args.rc_consistency, ema_decay=args.ema, pearson_weight=args.pearson_weight
        )
        ckpt_path = args.out_dir / f"MULTICELL_{mname}_s{args.seed}.pt"
        clean_sd = {k.replace("._orig_mod.", ".").replace("_orig_mod.", ""): v.detach().cpu() for k, v in model.state_dict().items()}
        torch.save({"model_name": mname, "state_dict": clean_sd, "info": info}, ckpt_path)

        ct_res = {}
        for ci, ct in enumerate(CELL_LINES):
            td = data["test_by_ct"][ct]
            tp, tstd = predict_multicell(model, td["x"], device, batch_size=args.batch_size, rc_avg=True)
            p_cre, y_cre, cm = cre_averaged_metrics(tp[:, ci], td["y_avg"], td["ids"])
            preds_by_ct[ct][mname] = p_cre
            y_by_ct[ct] = y_cre
            cal = round(float(spearmanr(tstd[:, ci], np.abs(tp[:, ci] - td["y_raw"])).statistic), 4) if np.std(tstd[:, ci]) > 1e-12 else None
            ct_res[ct] = {**{k: (round(v, 5) if isinstance(v, float) else v) for k, v in cm.items()}, "cal_spearman": cal}

        # Save predictions for exact cross-run paired bootstraps
        np.savez_compressed(
            args.out_dir / f"MULTICELL_{mname}_s{args.seed}_cre_preds.npz",
            **{f"p_{ct}": preds_by_ct[ct][mname] for ct in CELL_LINES},
            **{f"y_{ct}": y_by_ct[ct] for ct in CELL_LINES},
        )

        # If Looped Transformer, probe test-time loop depth scaling T in {1, 2, 4, 6, 8}
        loop_scaling = None
        if isinstance(model, LOOP_MODEL_TYPES):
            loop_scaling = {}
            orig_loops = model.n_loops
            for t_eval in (1, 2, 4, 6, 8):
                model.n_loops = t_eval
                t_res = {}
                for ci, ct in enumerate(CELL_LINES):
                    td = data["test_by_ct"][ct]
                    tp_t, _ = predict_multicell(model, td["x"], device, batch_size=args.batch_size, rc_avg=True)
                    _, _, cm_t = cre_averaged_metrics(tp_t[:, ci], td["y_avg"], td["ids"])
                    t_res[ct] = round(cm_t["pearson_r"], 5)
                t_res["mean"] = round(float(np.mean(list(t_res.values()))), 5)
                loop_scaling[f"T={t_eval}"] = t_res
            model.n_loops = orig_loops
            print(f"  ==> MULTI-LOOP SCALING ({mname}): {loop_scaling}", flush=True)

        cagi5 = eval_cagi5_multicell(model, args.kircher_dir, device, batch_size=args.batch_size)
        report["models"][mname] = dict(**info, test_fold1_cre=ct_res, loop_depth_scaling=loop_scaling, cagi5_zeroshot=cagi5)
        print(
            f"==> MULTI-{mname}: Fold1 CRE r: HepG2={ct_res['HepG2']['pearson_r']:.4f} K562={ct_res['K562']['pearson_r']:.4f} "
            f"WTC11={ct_res['WTC11']['pearson_r']:.4f} | CAGI5 ZS: blind={cagi5['mean_blind_3head_r']:.4f} "
            f"matched={cagi5['mean_lineage_matched_r']:.4f} sig={cagi5['mean_sig_p1e5_r']:.4f} "
            f"liver={cagi5['mean_liver_hepg2_r']:.4f} eryth={cagi5['mean_erythroid_k562_r']:.4f}",
            flush=True,
        )

    gains = {}
    loci = sorted(set(CAGI5_LOCUS.values()))
    for cand in ("helixloop_v4", "helixloop_v4_no_delta", "helixloop_omega", "helixloop_omega_adamw", "thermoloop_v3", "thermo_interlace_v3", "thermoloop_v2", "thermo_interlace_v2", "thermoloop", "thermo_interlace", "thermosyntax", "thermodream"):
        if cand not in report["models"]:
            continue
        for base in ("legnet_official", "dream_rnn", "legnet_large"):
            if base not in report["models"] or base == cand:
                continue
            for ct in CELL_LINES:
                gains[f"{ct}_{cand}_vs_{base}"] = paired_bootstrap_gain(preds_by_ct[ct][cand], preds_by_ct[ct][base], y_by_ct[ct], n_boot=2000)
            m_th = report["models"][cand]["cagi5_zeroshot"]["per_element_matched"]
            m_lg = report["models"][base]["cagi5_zeroshot"]["per_element_matched"]
            l_th = np.array([np.mean([m_th[e] for e, lc in CAGI5_LOCUS.items() if lc == loc and e in m_th]) for loc in loci])
            l_lg = np.array([np.mean([m_lg[e] for e, lc in CAGI5_LOCUS.items() if lc == loc and e in m_lg]) for loc in loci])
            rng = np.random.default_rng(2026)
            bidx = rng.integers(0, len(loci), size=(2000, len(loci)))
            bg = (l_th[bidx] - l_lg[bidx]).mean(axis=1)
            gains[f"cagi5_21loci_{cand}_vs_{base}"] = dict(
                point_gain=round(float((l_th - l_lg).mean()), 5),
                ci95_low=round(float(np.quantile(bg, 0.025)), 5),
                ci95_high=round(float(np.quantile(bg, 0.975)), 5),
                prob_better=round(float((bg > 0).mean()), 4),
                elements_won=f"{sum(m_th[e] > m_lg[e] for e in m_th)}/{len(m_th)}",
            )
    if gains:
        report["paired_bootstrap_gains"] = gains
        print("\n=== Multi-Lineage Paired Bootstrap Gains ===", flush=True)
        for k, v in gains.items():
            print(f"  {k:44s}: {v}", flush=True)

    tag_sfx = f"_{args.tag}" if args.tag else ""
    out_json = args.out_dir / f"benchmark_multicell{tag_sfx}_s{args.seed}.json"
    out_json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"\nSaved multi-lineage results to {out_json}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

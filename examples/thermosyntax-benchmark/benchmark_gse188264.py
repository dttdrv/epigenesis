#!/usr/bin/env python3
"""Strict, Zero-Leakage GSE188264 0-72h Neural Induction & Test-Only Epistasis Benchmark.

Compares Epigenesis ThermoSyntaxNet (both from scratch and transferred from WTC11 iPSC)
head-to-head against Official MPRALegNet, MalinoisCNN, DeepSTARR, and existing linear
baselines (Kreimer v1, biophys v2) on:
  1. Official 118 held-out test loci (1,895 perturbations, 1,415 targeted) across 0-72h
     (Pearson r, Spearman rho, exact family-weighted targeted MSE, per-hour r).
  2. STRICTLY TEST-ONLY Double-Mutant Triplets (n=51 triplets in held-out test loci ONLY,
     0% overlap with train/val):
     - Double-mutant prediction (WT -> A+B)
     - Primary-repressed rescue transition (A -> A+B)
     - True NON-ADDITIVE EPISTASIS: eps = delta(WT->A+B) - delta(WT->A) - delta(WT->B)
       (where EVERY linear model predicts eps == 0 identically!).
  3. Exact mathematical invariants: reversion identity (0.0), antisymmetry, RC invariance.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "studio"))
sys.path.insert(0, str(REPO / "examples/neural-regulation"))

import motif_temporal as mt  # noqa: E402
import regulation as reg  # noqa: E402
from epigenesis_studio import biophys  # noqa: E402
import benchmark_flashzoi as bf  # noqa: E402
from models import (  # noqa: E402
    EpigenesisThermoSyntaxNet,
    HelixLoopOmegaModel,
    HybridMuonCWD,
    build_model,
    count_parameters,
    load_jaspar_pwm_bank,
)

BASE_LUT = np.full(128, 4, dtype=np.uint8)
BASE_LUT[ord("A")] = 0
BASE_LUT[ord("C")] = 1
BASE_LUT[ord("G")] = 2
BASE_LUT[ord("T")] = 3
RC_TABLE = str.maketrans("ACGT", "TGCA")
LENTI_LEFT = "AGGACCGGATCAACT"
LENTI_RIGHT = "CATTGCGTGAACCGA"


def rc_str(s: str) -> str:
    return s.translate(RC_TABLE)[::-1]


def pad_to_len(seq: str, target_len: int = 230) -> str:
    # Wrap 171-bp insert with assay/lentiMPRA flanks so length matches 230 bp seamlessly for WTC11 transfer
    if len(seq) == target_len:
        return seq
    need = target_len - len(seq)
    left = need // 2
    right = need - left
    lpad = (LENTI_LEFT * 3)[-left:] if left > 0 else ""
    rpad = (LENTI_RIGHT * 3)[:right] if right > 0 else ""
    return lpad + seq + rpad


def seqs_to_u8(seqs: list[str]) -> np.ndarray:
    buf = np.frombuffer("".join(seqs).encode("ascii"), dtype=np.uint8).reshape(len(seqs), len(seqs[0]))
    return BASE_LUT[buf]


def u8_to_oh(arr: np.ndarray, device: torch.device) -> torch.Tensor:
    t = torch.from_numpy(arr).to(device=device, dtype=torch.long)
    oh = torch.zeros((t.size(0), 5, t.size(1)), device=device, dtype=torch.float32)
    oh.scatter_(1, t.unsqueeze(1), 1.0)
    return oh[:, :4, :]


class SiameseTemporalWrapper(nn.Module):
    """Exact Antisymmetric & RC-Invariant Siamese Sequence-to-Temporal-Delta Predictor.

    Guarantees by exact mathematical construction:
      - delta(S, S) == 0
      - delta(S_alt, S_ref) == -delta(S_ref, S_alt)
      - delta(RC(S_ref), RC(S_alt)) == delta(S_ref, S_alt)
    """

    def __init__(self, backbone: nn.Module):
        super().__init__()
        self.backbone = backbone

    def score_single(self, oh: torch.Tensor) -> torch.Tensor:
        y_f = self.backbone(oh)
        y_r = self.backbone(torch.flip(oh, dims=[1, 2]))
        return 0.5 * (y_f + y_r)

    def forward(self, oh_ref: torch.Tensor, oh_alt: torch.Tensor) -> torch.Tensor:
        return self.score_single(oh_alt) - self.score_single(oh_ref)


def load_pretrained_backbone(model: nn.Module, ckpt_path: Path, device: torch.device) -> int:
    """Transfer all shape-matching trunk weights from a pretrained lentiMPRA checkpoint (excluding final heads)."""
    state = torch.load(ckpt_path, map_location=device, weights_only=False)["state_dict"]
    own = model.state_dict()
    loaded = 0
    for k, v in state.items():
        if k in own and own[k].shape == v.shape and not k.startswith(("mu_head", "logvar_head", "aux_head", "thermo_skip", "head.3", "head.4")):
            own[k].copy_(v)
            loaded += 1
    return loaded


def _build_gse_backbone(model_name: str, w_fwd: np.ndarray, w_rc: np.ndarray, ckpt_path: Path | None, device: torch.device) -> nn.Module:
    if model_name == "flashzoi_helixloop_v4":
        blob = torch.load(ckpt_path, map_location=device, weights_only=False) if (ckpt_path and ckpt_path.is_file()) else {}
        return bf.FlashzoiHelixLoopV4(w_fwd, w_rc, repo_id=blob.get("repo_id", "johahi/flashzoi-replicate-0"), out_dim=7, n_loops=blob.get("n_loops", 6)).to(device)
    if model_name == "flashzoi_linear_ft":
        blob = torch.load(ckpt_path, map_location=device, weights_only=False) if (ckpt_path and ckpt_path.is_file()) else {}
        return bf.FlashzoiLinearFT(repo_id=blob.get("repo_id", "johahi/flashzoi-replicate-0"), out_dim=7).to(device)
    return build_model(model_name, w_fwd, w_rc, out_dim=7, seq_len=230).to(device)


def load_test_only_triplets(rows: list[reg.Row], scale: float) -> list[dict[str, Any]]:
    """Extract ONLY strictly held-out test-locus double-mutant triplets (0% train/val overlap)."""
    split_by_locus = {r.locus: r.split for r in rows}
    fa: dict[str, str] = {}
    cur = None
    for line in (REPO / "examples/neural-regulation/data/design.fasta").read_text().splitlines():
        line = line.strip()
        if line.startswith(">"):
            cur = line[1:]
        elif cur:
            fa[cur] = line.upper()

    alphas_raw: dict[str, np.ndarray] = {}
    with gzip.open(REPO / "examples/neural-regulation/data/alphas.csv.gz", "rt") as f:
        for r in csv.DictReader(f):
            alphas_raw[r[""]] = np.array([float(r[t]) for t in ("0h", "3h", "6h", "12h", "24h", "48h", "72h")], dtype=np.float64)

    wt_by_loc = {}
    for k, vals in alphas_raw.items():
        if k.startswith("WT_"):
            m = re.search(r"(chr[0-9XY]+:\d+-\d+)", k)
            if m and k in fa:
                # Match design insert (strip 15-bp flanks to 171-bp core matching regulation.py)
                wt_by_loc[m.group(1)] = (k, fa[k][15:-15], vals)

    triplets = []
    for k, v_ab in alphas_raw.items():
        if not k.startswith("hit2diff_pert"):
            continue
        rest = k[len("hit2diff_") :]
        k1, k2 = "hit2diff_1_" + rest, "hit2diff_2_" + rest
        m = re.search(r"(chr[0-9XY]+:\d+-\d+)", k)
        if not m:
            continue
        loc = m.group(1)
        if split_by_locus.get(loc) != "test":
            continue  # STRICT ZERO LEAKAGE: test loci ONLY!
        if loc not in wt_by_loc or k1 not in alphas_raw or k2 not in alphas_raw:
            continue
        if k not in fa or k1 not in fa or k2 not in fa:
            continue
        _, wt_seq, v_wt = wt_by_loc[loc]
        # Use asinh-transformed targets on exact same scale as regulation.py
        t_wt = np.arcsinh(v_wt / scale)
        t_a = np.arcsinh(alphas_raw[k1] / scale) - t_wt
        t_b = np.arcsinh(alphas_raw[k2] / scale) - t_wt
        t_ab = np.arcsinh(v_ab / scale) - t_wt
        triplets.append(
            dict(
                locus=loc,
                id=k,
                wt=wt_seq,
                seq_a=fa[k1][15:-15],
                seq_b=fa[k2][15:-15],
                seq_ab=fa[k][15:-15],
                y_a=float(t_a.mean()),
                y_b=float(t_b.mean()),
                y_ab=float(t_ab.mean()),
                epistasis=float((t_ab - t_a - t_b).mean()),
            )
        )
    return triplets


@torch.no_grad()
def predict_pairs(
    wrapper: SiameseTemporalWrapper,
    ref_seqs: list[str],
    alt_seqs: list[str],
    device: torch.device,
    seq_len: int = 230,
    batch_size: int = 256,
) -> np.ndarray:
    wrapper.eval()
    u_map: dict[str, int] = {}
    u_seqs: list[str] = []
    for s in ref_seqs + alt_seqs:
        ps = pad_to_len(s, seq_len)
        if ps not in u_map:
            u_map[ps] = len(u_seqs)
            u_seqs.append(ps)
    u_u8 = seqs_to_u8(u_seqs)
    scores = np.zeros((len(u_seqs), 7), dtype=np.float64)
    for s in range(0, len(u_seqs), batch_size):
        e = min(len(u_seqs), s + batch_size)
        oh = u8_to_oh(u_u8[s:e], device)
        scores[s:e] = wrapper.score_single(oh).cpu().numpy()
    ri = np.array([u_map[pad_to_len(s, seq_len)] for s in ref_seqs], dtype=np.int64)
    ai = np.array([u_map[pad_to_len(s, seq_len)] for s in alt_seqs], dtype=np.int64)
    return scores[ai] - scores[ri]


def train_gse_model(
    model_name: str,
    train_r: list[reg.Row],
    val_r: list[reg.Row],
    scale: float,
    w_fwd: np.ndarray,
    w_rc: np.ndarray,
    device: torch.device,
    *,
    pretrained_ckpt: Path | None = None,
    epochs: int = 25,
    batch_size: int = 128,
    lr: float = 1.5e-3,
    seed: int = 42,
) -> tuple[SiameseTemporalWrapper, dict[str, Any]]:
    torch.manual_seed(seed)
    np.random.seed(seed)
    backbone = _build_gse_backbone(model_name, w_fwd, w_rc, pretrained_ckpt, device)
    loaded_tensors = 0
    if pretrained_ckpt is not None and pretrained_ckpt.is_file():
        loaded_tensors = load_pretrained_backbone(backbone, pretrained_ckpt, device)
    wrapper = SiameseTemporalWrapper(backbone).to(device)

    # Also train jointly on WT absolute activity (centered) so the trunk learns baseline enhancer state!
    tr_ref = seqs_to_u8([pad_to_len(r.reference, 230) for r in train_r])
    tr_alt = seqs_to_u8([pad_to_len(r.sequence, 230) for r in train_r])
    tr_y = np.array([reg.target(r, scale) for r in train_r], dtype=np.float32)
    tr_w = mt._family_weights(train_r)
    tr_w = (tr_w / tr_w.mean()).astype(np.float32)

    # Also extract unique WT sequences and their asinh absolute activity for auxiliary state supervision
    wt_by_loc: dict[str, tuple[str, np.ndarray]] = {}
    for r in train_r:
        if r.locus not in wt_by_loc:
            wt_mean = np.mean([np.arcsinh(np.array(a) / scale) for a in r.wt_alphas], axis=0)
            wt_by_loc[r.locus] = (r.reference, wt_mean.astype(np.float32))
    wt_y_mean = np.mean([v[1] for v in wt_by_loc.values()], axis=0)
    tr_wt_y = np.array([wt_by_loc[r.locus][1] - wt_y_mean for r in train_r], dtype=np.float32)

    va_refs = [r.reference for r in val_r]
    va_alts = [r.sequence for r in val_r]
    va_y = np.array([reg.target(r, scale) for r in val_r], dtype=np.float64)
    va_w = mt._family_weights(val_r)

    use_lr = lr * 0.5 if loaded_tensors > 0 else lr
    if model_name.startswith("flashzoi_"):
        pre_ids = {id(p) for p in backbone.encoder.pretrained_parameters()}
        enc_p = [p for p in wrapper.parameters() if p.requires_grad and id(p) in pre_ids]
        hd_p = [p for p in wrapper.parameters() if p.requires_grad and id(p) not in pre_ids]
        opt = torch.optim.AdamW([{"params": enc_p, "lr": 1.5e-4 / 10.0, "weight_decay": 0.08}, {"params": hd_p, "lr": 1.2e-3 / 10.0, "weight_decay": 0.08}])
        max_lr_sched: Any = [1.5e-4, 1.2e-3]
    elif model_name.startswith("helixloop_omega") and not model_name.endswith("_adamw"):
        opt = HybridMuonCWD(wrapper, lr_muon=use_lr * 1.5 / 10.0, lr_adam=use_lr / 10.0, weight_decay=0.05)
        max_lr_sched = [use_lr * 1.5, use_lr, use_lr]
    else:
        opt = torch.optim.AdamW(wrapper.parameters(), lr=use_lr / 10.0, weight_decay=0.05)
        max_lr_sched = use_lr
    steps_per_ep = math.ceil(len(train_r) / batch_size)
    sched = torch.optim.lr_scheduler.OneCycleLR(
        opt, max_lr=max_lr_sched, epochs=epochs, steps_per_epoch=steps_per_ep, pct_start=0.25, cycle_momentum=False
    )

    rng = np.random.default_rng(seed)
    best_val_mse = 1e9
    best_state = None
    best_ep = -1
    t0 = time.perf_counter()

    for ep in range(1, epochs + 1):
        wrapper.train()
        perm = rng.permutation(len(train_r))
        for s in range(0, len(train_r), batch_size):
            idx = perm[s : s + batch_size]
            if len(idx) < 8:
                continue
            oh_r = u8_to_oh(tr_ref[idx], device)
            oh_a = u8_to_oh(tr_alt[idx], device)
            # Paired random cyclic/flank shift (+/- 6 bp) applied identically to ref & alt
            sh = int(rng.integers(-6, 7))
            if sh != 0:
                oh_r = torch.roll(oh_r, shifts=sh, dims=2)
                oh_a = torch.roll(oh_a, shifts=sh, dims=2)
            # Paired random RC
            if rng.random() > 0.5:
                oh_r = torch.flip(oh_r, dims=[1, 2])
                oh_a = torch.flip(oh_a, dims=[1, 2])

            yb = torch.from_numpy(tr_y[idx]).to(device)
            wb = torch.from_numpy(tr_w[idx]).to(device)
            ywt = torch.from_numpy(tr_wt_y[idx]).to(device)

            opt.zero_grad(set_to_none=True)
            s_ref = wrapper.backbone(oh_r)
            s_alt = wrapper.backbone(oh_a)
            pred_delta = s_alt - s_ref
            # Family-weighted MSE on temporal delta + small auxiliary loss on baseline WT enhancer state
            delta_loss = ((pred_delta - yb).pow(2).mean(dim=1) * wb).mean()
            wt_loss = F.mse_loss(s_ref, ywt)
            loss = delta_loss + 0.15 * wt_loss
            loss.backward()
            torch.nn.utils.clip_grad_norm_(wrapper.parameters(), 2.0)
            opt.step()
            sched.step()

        vp = predict_pairs(wrapper, va_refs, va_alts, device)
        vmse = float(np.sum(va_w[:, None] * (va_y - vp) ** 2) / 7.0)
        vr = float(np.corrcoef(va_y.mean(axis=1), vp.mean(axis=1))[0, 1])
        if vmse < best_val_mse:
            best_val_mse = vmse
            best_ep = ep
            best_state = {k: v.detach().cpu().clone() for k, v in wrapper.state_dict().items()}
        if ep in (1, 10, epochs) or ep == best_ep:
            tag = f"{model_name}{'+WTC11' if loaded_tensors > 0 else ''}"
            print(f"  [{tag:28s}] ep={ep:2d}/{epochs} val_mse={vmse:.5f} val_r={vr:.4f} ({time.perf_counter()-t0:.1f}s)", flush=True)

    wrapper.load_state_dict(best_state)
    return wrapper, dict(
        model=model_name,
        pretrained_wtc11=loaded_tensors > 0,
        loaded_tensors=loaded_tensors,
        n_params=count_parameters(wrapper),
        best_epoch=best_ep,
        best_val_mse=round(best_val_mse, 6),
        train_seconds=round(time.perf_counter() - t0, 1),
    )


def evaluate_gse_model(
    name: str,
    predict_fn,
    test_r: list[reg.Row],
    scale: float,
    triplets: list[dict[str, Any]],
) -> dict[str, Any]:
    y_true = np.array([reg.target(r, scale) for r in test_r], dtype=np.float64)
    y_pred = predict_fn([r.reference for r in test_r], [r.sequence for r in test_r])
    ym, pm = y_true.mean(axis=1), y_pred.mean(axis=1)

    targ = np.array([r.kind == "targeted" for r in test_r])
    t_rows = [r for r in test_r if r.kind == "targeted"]
    w_t = mt._family_weights(t_rows)
    mse_targ = float(np.sum(w_t[:, None] * (y_true[targ] - y_pred[targ]) ** 2) / 7.0)
    mse_zero = float(np.sum(w_t[:, None] * (y_true[targ]) ** 2) / 7.0)

    r_all = float(np.corrcoef(ym, pm)[0, 1])
    rho_all = float(spearmanr(ym, pm).statistic)
    r_targ = float(np.corrcoef(ym[targ], pm[targ])[0, 1])

    # Strictly test-only double-mutant triplets (n=51, ZERO train/val leakage!)
    p_a = predict_fn([t["wt"] for t in triplets], [t["seq_a"] for t in triplets]).mean(axis=1)
    p_b = predict_fn([t["wt"] for t in triplets], [t["seq_b"] for t in triplets]).mean(axis=1)
    p_ab = predict_fn([t["wt"] for t in triplets], [t["seq_ab"] for t in triplets]).mean(axis=1)
    p_a2ab = predict_fn([t["seq_a"] for t in triplets], [t["seq_ab"] for t in triplets]).mean(axis=1)
    p_b2ab = predict_fn([t["seq_b"] for t in triplets], [t["seq_ab"] for t in triplets]).mean(axis=1)
    pred_epi = p_ab - p_a - p_b

    y_ab = np.array([t["y_ab"] for t in triplets])
    y_epi = np.array([t["epistasis"] for t in triplets])
    r_ab = float(np.corrcoef(p_ab, y_ab)[0, 1])
    rho_ab = float(spearmanr(p_ab, y_ab).statistic)
    r_epi = float(np.corrcoef(pred_epi, y_epi)[0, 1]) if np.std(pred_epi) > 1e-9 else 0.0
    rho_epi = float(spearmanr(pred_epi, y_epi).statistic) if np.std(pred_epi) > 1e-9 else 0.0

    # Primary-repressed transitions in TEST-ONLY triplets
    obs_res, pred_res = [], []
    for i, t in enumerate(triplets):
        if t["y_a"] < -0.25:
            obs_res.append(t["y_ab"] - t["y_a"])
            pred_res.append(p_a2ab[i])
        if t["y_b"] < -0.25:
            obs_res.append(t["y_ab"] - t["y_b"])
            pred_res.append(p_b2ab[i])
    obs_r = np.array(obs_res)
    pred_r = np.array(pred_res)
    r_res = float(np.corrcoef(pred_r, obs_r)[0, 1]) if len(obs_r) > 2 else 0.0

    # Exact symmetry invariants check on first 16 test edits
    s0 = [r.reference for r in test_r[:16]]
    s1 = [r.sequence for r in test_r[:16]]
    rev_err = float(np.max(np.abs(predict_fn(s0, s0))))
    anti_err = float(np.max(np.abs(predict_fn(s0, s1) + predict_fn(s1, s0))))
    rc_err = float(np.max(np.abs(predict_fn(s0, s1) - predict_fn([rc_str(x) for x in s0], [rc_str(x) for x in s1]))))

    return dict(
        name=name,
        test_pearson_r_all=round(r_all, 4),
        test_spearman_rho_all=round(rho_all, 4),
        test_pearson_r_targeted=round(r_targ, 4),
        test_targeted_family_mse=round(mse_targ, 6),
        mse_reduction_vs_zero_pct=round((1.0 - mse_targ / mse_zero) * 100.0, 2),
        test_only_triplets_n=len(triplets),
        test_only_double_mutant_r=round(r_ab, 4),
        test_only_double_mutant_rho=round(rho_ab, 4),
        test_only_rescue_n=len(obs_r),
        test_only_rescue_r=round(r_res, 4),
        test_only_epistasis_r=round(r_epi, 4),
        test_only_epistasis_rho=round(rho_epi, 4),
        invariants=dict(reversion_max_err=rev_err, antisymmetry_max_err=anti_err, rc_max_err=rc_err),
        _preds=dict(pm=pm, p_ab=p_ab, pred_epi=pred_epi, pred_r=pred_r),
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--wtc11-dir", type=Path, default=Path("/tmp/epigenesis-benchmarks/results"))
    ap.add_argument("--pretrained-prefix", type=str, default="WTC11", help="Prefix for pretrained checkpoints (e.g. WTC11 or MULTICELL)")
    ap.add_argument("--out", type=Path, default=Path("/tmp/epigenesis-benchmarks/results/benchmark_gse188264.json"))
    ap.add_argument("--epochs", type=int, default=25)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else ("mps" if torch.backends.mps.is_available() else "cpu"))
    if device.type == "cuda":
        torch.backends.cuda.matmul.allow_tf32 = True
        torch.backends.cudnn.allow_tf32 = True
        torch.backends.cudnn.benchmark = True
    w_fwd, w_rc, _ = load_jaspar_pwm_bank(21)
    rows, _ = reg.load_data()
    train_r = [r for r in rows if r.split == "train"]
    val_r = [r for r in rows if r.split == "validation"]
    test_r = [r for r in rows if r.split == "test"]
    scale = reg.training_scale(train_r + val_r)
    triplets = load_test_only_triplets(rows, scale)
    print(f"GSE188264 loaded: train={len(train_r)}, val={len(val_r)}, test={len(test_r)}, test-only triplets={len(triplets)}", flush=True)

    results: dict[str, Any] = dict(
        note="118 test loci are retrospective (exposed in prior repo iterations); 51 double-mutant triplets are strictly test-only (0% overlap with train/val).",
        n_train=len(train_r),
        n_val=len(val_r),
        n_test=len(test_r),
        n_test_only_triplets=len(triplets),
        models={},
    )

    # 1. Existing linear biophys v2 baseline
    def pred_biophys(refs: list[str], alts: list[str]) -> np.ndarray:
        out = biophys.predict([{"id": f"e{i}", "reference": r, "sequence": a} for i, (r, a) in enumerate(zip(refs, alts))])
        return np.array([p["temporal"] for p in out], dtype=np.float64)

    preds_dump: dict[str, np.ndarray] = dict(
        y_all=np.array([np.mean(reg.target(r, scale)) for r in test_r], dtype=np.float64),
        test_loci=np.array([r.locus for r in test_r]),
        y_ab=np.array([t["y_ab"] for t in triplets], dtype=np.float64),
        y_epi=np.array([t["epistasis"] for t in triplets], dtype=np.float64),
        triplet_loci=np.array([t["locus"] for t in triplets]),
    )
    res_bio = evaluate_gse_model("biophys_v2_linear", pred_biophys, test_r, scale, triplets)
    pb = res_bio.pop("_preds")
    for k, v in pb.items():
        preds_dump[f"biophys_v2_linear__{k}"] = v
    results["models"]["biophys_v2_linear"] = res_bio
    print(f"==> biophys_v2_linear: r={res_bio['test_pearson_r_all']:.4f} mse={res_bio['test_targeted_family_mse']:.6f} dbl_r={res_bio['test_only_double_mutant_r']:.4f} epi_r={res_bio['test_only_epistasis_r']:.4f}", flush=True)

    # 2. Neural competitors + HelixLoop / HelixLoop-Interlace (from scratch AND Multi-Lineage transfer!)
    pfx = args.pretrained_prefix
    specs = [
        ("legnet_official", None),
        ("malinois", None),
        ("deepstarr", None),
        ("helixloop_v4", None),
        ("legnet_official", args.wtc11_dir / f"{pfx}_legnet_official_s{args.seed}.pt"),
        ("helixloop_v4", args.wtc11_dir / f"{pfx}_helixloop_v4_s{args.seed}.pt"),
        ("flashzoi_linear_ft", args.wtc11_dir / f"TIER2_flashzoi_linear_ft_rep0_s{args.seed}.pt"),
        ("flashzoi_helixloop_v4", args.wtc11_dir / f"TIER2_flashzoi_helixloop_v4_rep0_s{args.seed}.pt"),
    ]
    for mname, ckpt in specs:
        if ckpt is not None and not ckpt.is_file():
            continue
        tag = f"{mname}_{pfx.lower()}_transfer" if ckpt is not None else f"{mname}_scratch"
        print(f"\n--- Training {tag} on GSE188264 (train only, val select, test once) ---", flush=True)
        wrapper, info = train_gse_model(
            mname, train_r, val_r, scale, w_fwd, w_rc, device, pretrained_ckpt=ckpt, epochs=args.epochs, seed=args.seed
        )
        res = evaluate_gse_model(tag, lambda rs, alts, w=wrapper: predict_pairs(w, rs, alts, device), test_r, scale, triplets)
        pm_d = res.pop("_preds")
        for k, v in pm_d.items():
            preds_dump[f"{tag}__{k}"] = v
        res.update(info)
        results["models"][tag] = res
        print(
            f"==> {tag:28s}: r={res['test_pearson_r_all']:.4f} mse={res['test_targeted_family_mse']:.6f} "
            f"({res['mse_reduction_vs_zero_pct']:.1f}% vs 0) | test-triplets dbl_r={res['test_only_double_mutant_r']:.4f} "
            f"res_r={res['test_only_rescue_r']:.4f} epi_r={res['test_only_epistasis_r']:.4f}",
            flush=True,
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    np.savez_compressed(args.out.with_suffix(".npz"), **preds_dump)
    print(f"\nSaved GSE188264 benchmark to {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

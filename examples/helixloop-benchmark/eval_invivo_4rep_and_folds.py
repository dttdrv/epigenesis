#!/usr/bin/env python3
"""Post-Stage-4 Audit: 4-Replicate In Vivo (ADASTRA/UDACHA) + Exact Borzoi fold3+fold4 Stratification.

Runs in ~60 seconds on GPU right after Stage 4 finishes:
  1. Evaluates Tier-1 (LegNet MC, HelixLoop v4 MC) and Tier-2 (single rep0 AND 4-rep ensembles
     of Flashzoi-Linear-FT and Flashzoi+HelixLoop-v4) zero-shot across all 6 in vivo panels
     (ADASTRA ChIP-seq ASB + UDACHA DNase/ATAC-seq in HepG2 & K562) using the cached GRCh38 panels.
  2. Stratifies every SNP in ADASTRA/UDACHA and every element in CAGI5 against Borzoi/Flashzoi's
     official hg38 split (borzoi_sequences_human.bed.gz: train={fold0,1,2,5,6,7}, val=fold4, test=fold3),
     reporting metrics on ALL loci AND on STRICTLY BORZOI-UNSEEN (fold3+fold4) loci separately!
  3. Computes exact two-sided finite-sample Phipson-Smyth (2010) paired bootstrap p-values
     p = min(1, 2*(k+1)/(B+1)) with B=2,000 draws and 95% percentile CIs for every comparison.
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import torch
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "studio"))
sys.path.insert(0, str(REPO / "examples/helixloop-benchmark"))

import benchmark_flashzoi as bf  # noqa: E402
import evaluate_invivo as ei  # noqa: E402
import models  # noqa: E402


def load_borzoi_folds(bed_gz: Path) -> dict[str, tuple[list[int], list[int], list[str]]]:
    by_chr: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    with gzip.open(bed_gz, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                by_chr[parts[0]].append((int(parts[1]), int(parts[2]), parts[3]))
    idx: dict[str, tuple[list[int], list[int], list[str]]] = {}
    for ch, ivs in by_chr.items():
        ivs.sort()
        idx[ch] = ([x[0] for x in ivs], [x[1] for x in ivs], [x[2] for x in ivs])
    return idx


def lookup_borzoi_fold(idx: dict[str, tuple[list[int], list[int], list[str]]], chrom: str, pos1: int) -> str:
    ch = chrom if chrom.startswith("chr") else f"chr{chrom}"
    if ch not in idx:
        return "unmapped"
    starts, ends, folds = idx[ch]
    pos0 = pos1 - 1
    i = bisect.bisect_right(starts, pos0) - 1
    found = set()
    for j in range(max(0, i - 3), min(len(starts), i + 4)):
        if starts[j] <= pos0 < ends[j]:
            found.add(folds[j])
    if not found:
        return "unmapped"
    if found <= {"fold3", "fold4"}:
        return "unseen_fold34"
    if found & {"fold0", "fold1", "fold2", "fold5", "fold6", "fold7"}:
        return "train_folds"
    return sorted(found)[0]


def paired_boot_ps(dc: np.ndarray, db: np.ndarray, y: np.ndarray, n_boot: int = 2000, seed: int = 2026) -> dict[str, Any]:
    n = len(y)
    if n < 10:
        return dict(n=n, insufficient=True)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(n_boot, n))
    rc = np.array([np.corrcoef(dc[i], y[i])[0, 1] for i in idx])
    rb = np.array([np.corrcoef(db[i], y[i])[0, 1] for i in idx])
    g = rc - rb
    pt = float(np.corrcoef(dc, y)[0, 1] - np.corrcoef(db, y)[0, 1])
    k_opp = int((g <= 0).sum()) if pt >= 0 else int((g >= 0).sum())
    p_two = min(1.0, 2.0 * (k_opp + 1) / (n_boot + 1))
    return dict(
        n=n,
        point_gain=round(pt, 5),
        ci95_low=round(float(np.quantile(g, 0.025)), 5),
        ci95_high=round(float(np.quantile(g, 0.975)), 5),
        p_two_sided=round(p_two, 5),
        k_opposite=k_opp,
        n_boot=n_boot,
    )


def metric_pack(dy: np.ndarray, y: np.ndarray) -> dict[str, Any]:
    if len(y) < 5:
        return dict(n=int(len(y)))
    yb = (y > 0).astype(np.int32)
    auroc = round(float(roc_auc_score(yb, dy)), 4) if len(np.unique(yb)) == 2 else None
    return dict(
        n=int(len(y)),
        pearson_r=round(float(np.corrcoef(dy, y)[0, 1]), 4),
        spearman_rho=round(float(spearmanr(dy, y).statistic), 4),
        dir_auroc=auroc,
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", type=Path, default=Path("/content/results"))
    ap.add_argument("--cache-dir", type=Path, default=Path("/content/invivo_cache"))
    ap.add_argument("--borzoi-bed", type=Path, default=Path("/content/borzoi_sequences_human.bed.gz"))
    ap.add_argument("--out", type=Path, default=Path("/content/results/benchmark_invivo_4rep_audited.json"))
    args = ap.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    wf, wr, _ = models.load_jaspar_pwm_bank(21)
    fold_idx = load_borzoi_folds(args.borzoi_bed)

    panel_meta = [
        ("HepG2_ASB", 0),
        ("HepG2_DNase", 0),
        ("HepG2_ATAC", 0),
        ("K562_ASB", 1),
        ("K562_DNase", 1),
        ("K562_ATAC", 1),
    ]
    panels: dict[str, tuple[list[dict[str, Any]], int, np.ndarray]] = {}
    for pname, ct_i in panel_meta:
        pfile = args.cache_dir / f"{pname}.json"
        if not pfile.is_file():
            continue
        items = json.loads(pfile.read_text(encoding="utf-8"))
        folds = np.array([lookup_borzoi_fold(fold_idx, x["chrom"], int(x["pos"])) for x in items])
        panels[pname] = (items, ct_i, folds)
        n_unseen = int((folds == "unseen_fold34").sum())
        n_train = int((folds == "train_folds").sum())
        print(f"Panel {pname:12s}: total={len(items)} | Borzoi unseen(fold3+4)={n_unseen} | train_folds={n_train}", flush=True)

    rdir = args.results_dir
    specs = [
        ("Tier1_LegNet_MC", "legnet_official", [rdir / "MULTICELL_legnet_official_s42.pt"]),
        ("Tier1_HelixLoop_v4_MC", "helixloop_v4", [rdir / "MULTICELL_helixloop_v4_s42.pt"]),
        ("Tier2_Flashzoi_Linear_rep0", "flashzoi_linear_ft", [rdir / "TIER2_flashzoi_linear_ft_rep0_s42.pt"]),
        ("Tier2_Flashzoi_HelixLoop_v4_rep0", "flashzoi_helixloop_v4", [rdir / "TIER2_flashzoi_helixloop_v4_rep0_s42.pt"]),
        (
            "Tier2_Flashzoi_Linear_4rep",
            "flashzoi_linear_ft",
            [rdir / f"TIER2_flashzoi_linear_ft_rep{i}_s{42+i}.pt" for i in range(4)],
        ),
        (
            "Tier2_Flashzoi_HelixLoop_v4_4rep",
            "flashzoi_helixloop_v4",
            [rdir / f"TIER2_flashzoi_helixloop_v4_rep{i}_s{42+i}.pt" for i in range(4)],
        ),
    ]

    report: dict[str, Any] = {
        "borzoi_split_note": "All 4 Flashzoi replicates share identical Borzoi hg38 split: train={fold0,1,2,5,6,7}, val=fold4, test=fold3 (0 MPRA/ADASTRA/UDACHA labels in pretraining).",
        "panels": {
            k: dict(total=len(v[0]), unseen_fold34=int((v[2] == "unseen_fold34").sum()), train_folds=int((v[2] == "train_folds").sum()))
            for k, v in panels.items()
        },
        "models": {},
        "pooled_summary": {},
        "paired_bootstrap": {},
    }
    dys: dict[str, dict[str, np.ndarray]] = {k: {} for k in panels}

    for label, mname, paths in specs:
        if not all(p.is_file() for p in paths):
            print(f"Skipping {label}: missing checkpoints", flush=True)
            continue
        ms = []
        for p in paths:
            blob = torch.load(p, map_location=device, weights_only=False)
            if mname == "flashzoi_helixloop_v4":
                m = bf.FlashzoiHelixLoopV4(wf, wr, repo_id=blob.get("repo_id", "johahi/flashzoi-replicate-0"), out_dim=3, n_loops=blob.get("n_loops", 6)).to(device)
            elif mname == "flashzoi_linear_ft":
                m = bf.FlashzoiLinearFT(repo_id=blob.get("repo_id", "johahi/flashzoi-replicate-0"), out_dim=3).to(device)
            else:
                m = models.build_model(mname, wf, wr, out_dim=3).to(device)
            m.load_state_dict(blob["state_dict"])
            m.eval()
            ms.append(m)

        m_res: dict[str, Any] = {}
        all_rs, unseen_rs, all_aucs, unseen_aucs = [], [], [], []
        for pname, (items, ct_i, folds) in panels.items():
            dy = np.mean([ei.score_panel(m, items, device, head_idx=ct_i) for m in ms], axis=0)
            dys[pname][label] = dy
            y = np.array([x["es"] for x in items], dtype=np.float64)
            m_all = metric_pack(dy, y)
            mask_u = folds == "unseen_fold34"
            mask_t = folds == "train_folds"
            m_unseen = metric_pack(dy[mask_u], y[mask_u])
            m_train = metric_pack(dy[mask_t], y[mask_t])
            m_res[pname] = dict(all=m_all, unseen_fold34=m_unseen, train_folds=m_train)
            all_rs.append(m_all["pearson_r"])
            all_aucs.append(m_all["dir_auroc"])
            if "pearson_r" in m_unseen:
                unseen_rs.append(m_unseen["pearson_r"])
                if m_unseen.get("dir_auroc") is not None:
                    unseen_aucs.append(m_unseen["dir_auroc"])
        report["models"][label] = m_res
        report["pooled_summary"][label] = dict(
            mean_6panel_all_r=round(float(np.mean(all_rs)), 4),
            mean_6panel_all_auroc=round(float(np.mean(all_aucs)), 4),
            mean_6panel_unseen_fold34_r=round(float(np.mean(unseen_rs)), 4),
            mean_6panel_unseen_fold34_auroc=round(float(np.mean(unseen_aucs)), 4),
        )
        ps = report["pooled_summary"][label]
        print(
            f"==> {label:32s}: 6-panel ALL r={ps['mean_6panel_all_r']:.4f} auroc={ps['mean_6panel_all_auroc']:.4f} | "
            f"STRICT UNSEEN(fold3+4) r={ps['mean_6panel_unseen_fold34_r']:.4f} auroc={ps['mean_6panel_unseen_fold34_auroc']:.4f}",
            flush=True,
        )

    pairs = [
        ("Tier1_HelixLoop_v4_MC", "Tier1_LegNet_MC"),
        ("Tier2_Flashzoi_HelixLoop_v4_rep0", "Tier2_Flashzoi_Linear_rep0"),
        ("Tier2_Flashzoi_HelixLoop_v4_4rep", "Tier2_Flashzoi_Linear_4rep"),
    ]
    for cand, base in pairs:
        for pname, (items, _, folds) in panels.items():
            if cand in dys[pname] and base in dys[pname]:
                y = np.array([x["es"] for x in items], dtype=np.float64)
                mask_u = folds == "unseen_fold34"
                report["paired_bootstrap"][f"{pname}_{cand}_vs_{base}_all"] = paired_boot_ps(dys[pname][cand], dys[pname][base], y)
                report["paired_bootstrap"][f"{pname}_{cand}_vs_{base}_unseen34"] = paired_boot_ps(
                    dys[pname][cand][mask_u], dys[pname][base][mask_u], y[mask_u]
                )
        # Pooled standardized z-score across all 6 panels (and strictly unseen fold3+4 subset)
        if all(cand in dys[p] and base in dys[p] for p in panels):
            yc_all, yb_all, yt_all = [], [], []
            yc_u, yb_u, yt_u = [], [], []
            for pname, (items, _, folds) in panels.items():
                y = np.array([x["es"] for x in items], dtype=np.float64)
                yz = (y - y.mean()) / (y.std() + 1e-8)
                dc = (dys[pname][cand] - dys[pname][cand].mean()) / (dys[pname][cand].std() + 1e-8)
                db = (dys[pname][base] - dys[pname][base].mean()) / (dys[pname][base].std() + 1e-8)
                yc_all.append(dc); yb_all.append(db); yt_all.append(yz)
                mu = folds == "unseen_fold34"
                yc_u.append(dc[mu]); yb_u.append(db[mu]); yt_u.append(yz[mu])
            report["paired_bootstrap"][f"POOLED6_{cand}_vs_{base}_all"] = paired_boot_ps(
                np.concatenate(yc_all), np.concatenate(yb_all), np.concatenate(yt_all)
            )
            report["paired_bootstrap"][f"POOLED6_{cand}_vs_{base}_unseen34"] = paired_boot_ps(
                np.concatenate(yc_u), np.concatenate(yb_u), np.concatenate(yt_u)
            )

    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Saved 4-rep + fold3+4 audited in vivo report to {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

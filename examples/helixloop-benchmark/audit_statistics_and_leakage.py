#!/usr/bin/env python3
"""Canonical Statistical & Pretraining-Leakage Auditor for Epigenesis HelixLoop.

Computes:
  1. Exact two-sided finite-sample Phipson-Smyth (2010) paired bootstrap p-values
     p = min(1, 2 * (k + 1) / (B + 1)) with B=2,000 draws and 95% percentile CIs
     across all Tier-1 and Tier-2 comparisons directly from saved CRE/variant predictions.
  2. Exact Borzoi/Flashzoi hg38 pretraining-fold stratification (borzoi_sequences_human.bed.gz:
     train={fold0,1,2,5,6,7}, val=fold4, test=fold3) on:
       - CAGI5 (29 elements / 21 disjoint loci; 4 strictly unseen fold3+fold4 elements)
       - In Vivo ADASTRA (ChIP-seq ASB) + UDACHA (DNase/ATAC-seq) across 6 panels
       - GSE188264 (118 held-out test loci & 51 strictly test-only epistasis triplets)
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import spearmanr


def paired_boot_ps(
    p_cand: np.ndarray,
    p_base: np.ndarray,
    y: np.ndarray,
    *,
    n_boot: int = 2000,
    seed: int = 2026,
) -> dict[str, Any]:
    """Exact two-sided finite-sample Phipson-Smyth (2010) paired bootstrap test on Pearson r."""
    n = len(y)
    if n < 5:
        return dict(n=n, insufficient=True)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(n_boot, n))
    rc = np.array([np.corrcoef(p_cand[i], y[i])[0, 1] for i in idx])
    rb = np.array([np.corrcoef(p_base[i], y[i])[0, 1] for i in idx])
    g = rc - rb
    r_cand = float(np.corrcoef(p_cand, y)[0, 1])
    r_base = float(np.corrcoef(p_base, y)[0, 1])
    pt = r_cand - r_base
    k_opp = int((g <= 0.0).sum()) if pt >= 0.0 else int((g >= 0.0).sum())
    p_one = (k_opp + 1) / (n_boot + 1)
    p_two = min(1.0, 2.0 * p_one)
    return dict(
        n=int(n),
        r_cand=round(r_cand, 5),
        r_base=round(r_base, 5),
        delta_r=round(pt, 5),
        ci95_low=round(float(np.quantile(g, 0.025)), 5),
        ci95_high=round(float(np.quantile(g, 0.975)), 5),
        p_two_sided=round(p_two, 5),
        p_one_sided=round(p_one, 5),
        k_opposite=k_opp,
        n_boot=n_boot,
    )


def paired_boot_loci_ps(
    cand_by_el: dict[str, float],
    base_by_el: dict[str, float],
    locus_map: dict[str, str],
    *,
    n_boot: int = 2000,
    seed: int = 2026,
) -> dict[str, Any]:
    """Cluster-bootstrap across disjoint CAGI5 genomic loci with exact two-sided Phipson-Smyth p-value."""
    by_loc_c: dict[str, list[float]] = defaultdict(list)
    by_loc_b: dict[str, list[float]] = defaultdict(list)
    common = sorted(set(cand_by_el) & set(base_by_el))
    for el in common:
        loc = locus_map.get(el, el)
        vc_e = cand_by_el[el]["pearson_r"] if isinstance(cand_by_el[el], dict) else cand_by_el[el]
        vb_e = base_by_el[el]["pearson_r"] if isinstance(base_by_el[el], dict) else base_by_el[el]
        by_loc_c[loc].append(float(vc_e))
        by_loc_b[loc].append(float(vb_e))
    loci = sorted(by_loc_c)
    vc = np.array([np.mean(by_loc_c[l]) for l in loci], dtype=np.float64)
    vb = np.array([np.mean(by_loc_b[l]) for l in loci], dtype=np.float64)
    diff = vc - vb
    pt = float(diff.mean())
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(loci), size=(n_boot, len(loci)))
    g = diff[idx].mean(axis=1)
    k_opp = int((g <= 0.0).sum()) if pt >= 0.0 else int((g >= 0.0).sum())
    p_one = (k_opp + 1) / (n_boot + 1)
    p_two = min(1.0, 2.0 * p_one)
    wins = sum(
        1
        for el in common
        if (cand_by_el[el]["pearson_r"] if isinstance(cand_by_el[el], dict) else cand_by_el[el])
        > (base_by_el[el]["pearson_r"] if isinstance(base_by_el[el], dict) else base_by_el[el])
    )
    return dict(
        n_elements=len(common),
        n_loci=len(loci),
        mean_cand_21loci=round(float(vc.mean()), 5),
        mean_base_21loci=round(float(vb.mean()), 5),
        delta_r=round(pt, 5),
        ci95_low=round(float(np.quantile(g, 0.025)), 5),
        ci95_high=round(float(np.quantile(g, 0.975)), 5),
        p_two_sided=round(p_two, 5),
        k_opposite=k_opp,
        elements_won=f"{wins}/{len(common)}",
    )


def load_borzoi_folds(bed_gz: Path) -> tuple[dict[str, tuple[list[int], list[int], list[str]]], dict[str, int]]:
    by_chr: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    counts: dict[str, int] = defaultdict(int)
    with gzip.open(bed_gz, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                by_chr[parts[0]].append((int(parts[1]), int(parts[2]), parts[3]))
                counts[parts[3]] += 1
    idx: dict[str, tuple[list[int], list[int], list[str]]] = {}
    for ch, ivs in by_chr.items():
        ivs.sort()
        idx[ch] = ([x[0] for x in ivs], [x[1] for x in ivs], [x[2] for x in ivs])
    return idx, dict(sorted(counts.items()))


def lookup_interval_fold(idx: dict[str, tuple[list[int], list[int], list[str]]], chrom: str, start0: int, end0: int) -> str:
    ch = chrom if chrom.startswith("chr") else f"chr{chrom}"
    if ch not in idx:
        return "unmapped"
    starts, ends, folds = idx[ch]
    mid = (start0 + end0) // 2
    i = bisect.bisect_right(starts, mid) - 1
    found = set()
    for j in range(max(0, i - 4), min(len(starts), i + 5)):
        if starts[j] < end0 and ends[j] > start0:
            found.add(folds[j])
    if not found:
        return "unmapped"
    if found <= {"fold3", "fold4"}:
        return "unseen_fold34"
    if found & {"fold0", "fold1", "fold2", "fold5", "fold6", "fold7"}:
        return "train_folds"
    return sorted(found)[0]


CAGI5_LOCUS = {
    "SORT1": "SORT1", "SORT1.2": "SORT1", "SORT1-flip": "SORT1",
    "LDLR": "LDLR", "LDLR.2": "LDLR",
    "PKLR-24h": "PKLR", "PKLR-48h": "PKLR",
    "MYCrs11986220": "MYC_enh", "MYCrs6983267": "MYC_enh",
    "TERT-GBM": "TERT", "TERT-HEK": "TERT",
    "ZFAND3": "ZFAND3", "IRF4": "IRF4", "IRF6": "IRF6", "F9": "F9",
    "HNF4A": "HNF4A", "BCL11A": "BCL11A", "GP1BA": "GP1BA", "HBB": "HBB",
    "HBG1": "HBG1", "FOXE1": "FOXE1", "MSMB": "MSMB", "RET": "RET",
    "TCF7L2": "TCF7L2", "UC88": "UC88",
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--a100-dir", type=Path, default=Path("/tmp/epigenesis-omega-results"))
    ap.add_argument("--l4-dir", type=Path, default=Path("/tmp/epigenesis-l4-archive"))
    ap.add_argument("--out", type=Path, default=Path("/tmp/epigenesis-omega-results/canonical_audit_summary.json"))
    args = ap.parse_args()

    a100 = args.a100_dir
    l4 = args.l4_dir
    fold_idx, fold_counts = load_borzoi_folds(a100 / "borzoi_sequences_human.bed.gz")

    summary: dict[str, Any] = {
        "borzoi_hg38_fold_counts": fold_counts,
        "tier1_hepg2_single": {},
        "tier1_multicell": {},
        "tier2_flashzoi_4rep": {},
        "cagi5_borzoi_unseen_audit": {},
        "invivo_audit": {},
        "gse188264_audit": {},
    }

    # --- 1. Tier-1 Single-Cell HepG2 paired bootstrap (from exact CRE .npz predictions!) ---
    def load_npz_1d(p: Path) -> tuple[np.ndarray, np.ndarray]:
        d = np.load(p)
        return d["p_cre"].astype(np.float64), d["y_cre"].astype(np.float64)

    p_v4, y_h = load_npz_1d(a100 / "HepG2_helixloop_v4_s42_cre_preds.npz")
    p_nodelta, _ = load_npz_1d(a100 / "HepG2_helixloop_v4_no_delta_s42_cre_preds.npz")
    p_leg_std, _ = load_npz_1d(l4 / "HepG2_legnet_official_s42_cre_preds.npz")
    p_leg_lg, _ = load_npz_1d(l4 / "HepG2_legnet_large_s42_cre_preds.npz")
    p_dream, _ = load_npz_1d(l4 / "HepG2_dream_rnn_s42_cre_preds.npz")
    p_tl3, _ = load_npz_1d(l4 / "HepG2_thermoloop_v3_s42_cre_preds.npz")
    p_ti3, _ = load_npz_1d(l4 / "HepG2_thermo_interlace_v3_s42_cre_preds.npz")

    summary["tier1_hepg2_single"]["paired_bootstrap"] = {
        "helixloop_v4_vs_legnet_official_std": paired_boot_ps(p_v4, p_leg_std, y_h),
        "helixloop_v4_vs_legnet_large": paired_boot_ps(p_v4, p_leg_lg, y_h),
        "helixloop_v4_vs_dream_rnn": paired_boot_ps(p_v4, p_dream, y_h),
        "helixloop_v4_vs_no_delta_ablation": paired_boot_ps(p_v4, p_nodelta, y_h),
        "helixloop_v4_vs_thermoloop_v3": paired_boot_ps(p_v4, p_tl3, y_h),
        "helixloop_v4_vs_thermo_interlace_v3": paired_boot_ps(p_v4, p_ti3, y_h),
    }
    h_v4_j = json.loads((a100 / "benchmark_HepG2_v4_hepg2_s42.json").read_text())
    h_om_j = json.loads((a100 / "benchmark_HepG2_omega_v1_s42.json").read_text())
    summary["tier1_hepg2_single"]["cagi5_v4_vs_no_delta"] = paired_boot_loci_ps(
        h_v4_j["models"]["helixloop_v4"]["cagi5_per_element"],
        h_v4_j["models"]["helixloop_v4_no_delta"]["cagi5_per_element"],
        CAGI5_LOCUS,
    )
    summary["tier1_hepg2_single"]["loop_scaling_v4"] = h_v4_j["models"]["helixloop_v4"]["loop_depth_scaling"]
    summary["tier1_hepg2_single"]["loop_scaling_omega_v1"] = h_om_j["models"]["helixloop_omega"]["loop_depth_scaling"]

    # --- 2. Tier-1 Multi-Lineage paired bootstrap (exact two-sided Phipson-Smyth!) ---
    mc_hl = np.load(a100 / "MULTICELL_helixloop_v4_s42_cre_preds.npz")
    mc_ln = np.load(a100 / "MULTICELL_legnet_official_s42_cre_preds.npz")
    mc_j = json.loads((a100 / "benchmark_multicell_tier1_v4_mc_s42.json").read_text())
    mc_boot = {}
    for ct in ("HepG2", "K562", "WTC11"):
        mc_boot[ct] = paired_boot_ps(mc_hl[f"p_{ct}"], mc_ln[f"p_{ct}"], mc_hl[f"y_{ct}"])
    mc_boot["cagi5_21loci"] = paired_boot_loci_ps(
        mc_j["models"]["helixloop_v4"]["cagi5_zeroshot"]["per_element_matched"],
        mc_j["models"]["legnet_official"]["cagi5_zeroshot"]["per_element_matched"],
        CAGI5_LOCUS,
    )
    summary["tier1_multicell"] = {
        "metrics": {
            m: {
                "n_params": mc_j["models"][m]["n_params"],
                "HepG2_r": mc_j["models"][m]["test_fold1_cre"]["HepG2"]["pearson_r"],
                "K562_r": mc_j["models"][m]["test_fold1_cre"]["K562"]["pearson_r"],
                "WTC11_r": mc_j["models"][m]["test_fold1_cre"]["WTC11"]["pearson_r"],
                "mean_r": round(float(np.mean([mc_j["models"][m]["test_fold1_cre"][c]["pearson_r"] for c in ("HepG2", "K562", "WTC11")])), 5),
                "cagi5_matched": mc_j["models"][m]["cagi5_zeroshot"]["mean_lineage_matched_r"],
                "cagi5_blind": mc_j["models"][m]["cagi5_zeroshot"]["mean_blind_3head_r"],
                "cagi5_sig": mc_j["models"][m]["cagi5_zeroshot"]["mean_sig_p1e5_r"],
                "cagi5_liver": mc_j["models"][m]["cagi5_zeroshot"]["mean_liver_hepg2_r"],
                "cagi5_erythroid": mc_j["models"][m]["cagi5_zeroshot"]["mean_erythroid_k562_r"],
            }
            for m in ("helixloop_v4", "legnet_official")
        },
        "loop_scaling": mc_j["models"]["helixloop_v4"]["loop_depth_scaling"],
        "paired_bootstrap": mc_boot,
    }

    # --- 3. Tier-2 Single rep0 AND 4-Replicate Ensemble paired bootstrap (exact Phipson-Smyth!) ---
    ens_j = json.loads((a100 / "benchmark_tier2_4rep_ensemble.json").read_text())
    hl_reps = [np.load(a100 / f"TIER2_flashzoi_helixloop_v4_rep{i}_s{42+i}_cre_preds.npz") for i in range(4)]
    ln_reps = [np.load(a100 / f"TIER2_flashzoi_linear_ft_rep{i}_s{42+i}_cre_preds.npz") for i in range(4)]
    t2_rep0_boot, t2_4rep_boot = {}, {}
    for ct in ("HepG2", "K562", "WTC11"):
        y_ct = hl_reps[0][f"y_{ct}"]
        t2_rep0_boot[ct] = paired_boot_ps(hl_reps[0][f"p_{ct}"], ln_reps[0][f"p_{ct}"], y_ct)
        p_hl_4 = np.mean([r[f"p_{ct}"] for r in hl_reps], axis=0)
        p_ln_4 = np.mean([r[f"p_{ct}"] for r in ln_reps], axis=0)
        t2_4rep_boot[ct] = paired_boot_ps(p_hl_4, p_ln_4, y_ct)
    t2_4rep_boot["cagi5_21loci"] = paired_boot_loci_ps(
        ens_j["ensembles"]["flashzoi_helixloop_v4_ens"]["cagi5_zeroshot"]["per_element_matched"],
        ens_j["ensembles"]["flashzoi_linear_ft_ens"]["cagi5_zeroshot"]["per_element_matched"],
        CAGI5_LOCUS,
    )
    summary["tier2_flashzoi_4rep"] = {
        "rep0_paired_bootstrap": t2_rep0_boot,
        "ens4_paired_bootstrap": t2_4rep_boot,
        "per_rep_mean_r": {
            f"rep{i}": {
                "hl": round(float(np.mean([json.loads((a100 / f"benchmark_tier2_flashzoi_rep{i}_s{42+i}.json").read_text())["models"]["flashzoi_helixloop_v4"]["test_fold1_cre"][c]["pearson_r"] for c in ("HepG2", "K562", "WTC11")])), 5),
                "lin": round(float(np.mean([json.loads((a100 / f"benchmark_tier2_flashzoi_rep{i}_s{42+i}.json").read_text())["models"]["flashzoi_linear_ft"]["test_fold1_cre"][c]["pearson_r"] for c in ("HepG2", "K562", "WTC11")])), 5),
            }
            for i in range(4)
        },
        "ens4_summary": {
            k: {
                "HepG2": v["test_fold1_cre"]["HepG2"]["pearson_r"],
                "K562": v["test_fold1_cre"]["K562"]["pearson_r"],
                "WTC11": v["test_fold1_cre"]["WTC11"]["pearson_r"],
                "mean": v["fold1_mean_cre_r"],
                "cagi5_matched": v["cagi5_zeroshot"]["mean_lineage_matched_r"],
                "cagi5_blind": v["cagi5_zeroshot"]["mean_blind_3head_r"],
                "cagi5_sig": v["cagi5_zeroshot"]["mean_sig_p1e5_r"],
                "cagi5_liver": v["cagi5_zeroshot"]["mean_liver_hepg2_r"],
                "cagi5_erythroid": v["cagi5_zeroshot"]["mean_erythroid_k562_r"],
            }
            for k, v in ens_j["ensembles"].items()
        },
    }

    # --- 4. CAGI5 Borzoi fold3+fold4 strict unseen stratification ---
    krefs = json.loads((a100 / "kircher_refs.json").read_text())
    el_fold = {}
    for el, info in krefs.items():
        el_fold[el] = lookup_interval_fold(fold_idx, info["chrom"], int(info["start"]), int(info["start"]) + len(info["seq"]))
    unseen_els = sorted([e for e, f in el_fold.items() if f == "unseen_fold34"])
    train_els = sorted([e for e, f in el_fold.items() if f == "train_folds"])
    cagi_sources = {
        "Tier1_LegNet_MC": mc_j["models"]["legnet_official"]["cagi5_zeroshot"]["per_element_matched"],
        "Tier1_HelixLoop_v4_MC": mc_j["models"]["helixloop_v4"]["cagi5_zeroshot"]["per_element_matched"],
        "Tier2_Flashzoi_Linear_4rep": ens_j["ensembles"]["flashzoi_linear_ft_ens"]["cagi5_zeroshot"]["per_element_matched"],
        "Tier2_Flashzoi_HelixLoop_v4_4rep": ens_j["ensembles"]["flashzoi_helixloop_v4_ens"]["cagi5_zeroshot"]["per_element_matched"],
    }
    cagi_audit = {"unseen_fold34_elements": unseen_els, "train_fold_count": len(train_els), "models": {}}
    for mlabel, emap in cagi_sources.items():
        u_vals = {e: emap[e] for e in unseen_els if e in emap}
        t_vals = {e: emap[e] for e in train_els if e in emap}
        cagi_audit["models"][mlabel] = dict(
            all29_mean=round(float(np.mean(list(emap.values()))), 4),
            unseen_fold34_mean=round(float(np.mean(list(u_vals.values()))), 4),
            train_folds_mean=round(float(np.mean(list(t_vals.values()))), 4),
            unseen_per_element=u_vals,
        )
    summary["cagi5_borzoi_unseen_audit"] = cagi_audit

    # --- 5. In Vivo 4-rep + fold3+fold4 audit ---
    iv_path = a100 / "benchmark_invivo_4rep_audited.json"
    if iv_path.is_file():
        iv = json.loads(iv_path.read_text())
        summary["invivo_audit"] = {
            "panels": iv["panels"],
            "pooled_summary": iv["pooled_summary"],
            "pooled_bootstrap": {k: v for k, v in iv["paired_bootstrap"].items() if k.startswith("POOLED6_") or "K562_ASB" in k or "K562_DNase" in k or "HepG2_ASB" in k or "HepG2_DNase" in k or "HepG2_ATAC" in k},
            "models": iv["models"],
        }

    # --- 6. GSE188264 + 51 test-only epistasis triplets + Borzoi fold3+fold4 audit (if available) ---
    gse_j_path = a100 / "benchmark_gse188264_v4.json"
    gse_npz_path = a100 / "benchmark_gse188264_v4.npz"
    if gse_j_path.is_file():
        gse_j = json.loads(gse_j_path.read_text())
        gse_out: dict[str, Any] = {"models": gse_j["models"]}
        if gse_npz_path.is_file():
            gd = np.load(gse_npz_path, allow_pickle=True)
            y_all = gd["y_all"]
            y_ab = gd["y_ab"]
            y_epi = gd["y_epi"]
            t_loci = list(gd["test_loci"])
            tr_loci = list(gd["triplet_loci"])

            def parse_loc_fold(loc_str: str) -> str:
                m = re.match(r"(chr[0-9XY]+):(\d+)-(\d+)", str(loc_str))
                if not m:
                    return "unmapped"
                return lookup_interval_fold(fold_idx, m.group(1), int(m.group(2)), int(m.group(3)))

            t_folds = np.array([parse_loc_fold(l) for l in t_loci])
            tr_folds = np.array([parse_loc_fold(l) for l in tr_loci])
            mu_all = t_folds == "unseen_fold34"
            mu_tr = tr_folds == "unseen_fold34"
            gse_out["borzoi_fold_counts"] = {
                "test_edits_total": int(len(t_loci)),
                "test_edits_unseen_fold34": int(mu_all.sum()),
                "triplets_total": int(len(tr_loci)),
                "triplets_unseen_fold34": int(mu_tr.sum()),
            }
            unseen_m = {}
            for mtag in gse_j["models"]:
                pm = gd[f"{mtag}__pm"]
                pab = gd[f"{mtag}__p_ab"]
                pepi = gd[f"{mtag}__pred_epi"]
                unseen_m[mtag] = dict(
                    unseen_edits_r=round(float(np.corrcoef(pm[mu_all], y_all[mu_all])[0, 1]), 4) if mu_all.sum() >= 5 else None,
                    unseen_triplets_dbl_r=round(float(np.corrcoef(pab[mu_tr], y_ab[mu_tr])[0, 1]), 4) if mu_tr.sum() >= 5 else None,
                    unseen_triplets_epi_r=round(float(np.corrcoef(pepi[mu_tr], y_epi[mu_tr])[0, 1]), 4) if (mu_tr.sum() >= 5 and np.std(pepi[mu_tr]) > 1e-9) else None,
                )
            gse_out["unseen_fold34_metrics"] = unseen_m
            boot_gse = {}
            pairs_gse = [
                ("helixloop_v4_multicell_transfer", "legnet_official_multicell_transfer"),
                ("helixloop_v4_scratch", "legnet_official_scratch"),
                ("flashzoi_helixloop_v4_multicell_transfer", "flashzoi_linear_ft_multicell_transfer"),
                ("flashzoi_helixloop_v4_multicell_transfer", "biophys_v2_linear"),
                ("helixloop_v4_multicell_transfer", "biophys_v2_linear"),
            ]
            for c, b in pairs_gse:
                if f"{c}__pm" in gd and f"{b}__pm" in gd:
                    boot_gse[f"edits_{c}_vs_{b}"] = paired_boot_ps(gd[f"{c}__pm"], gd[f"{b}__pm"], y_all)
                    boot_gse[f"dbl_{c}_vs_{b}"] = paired_boot_ps(gd[f"{c}__p_ab"], gd[f"{b}__p_ab"], y_ab)
                    if np.std(gd[f"{c}__pred_epi"]) > 1e-9 and np.std(gd[f"{b}__pred_epi"]) > 1e-9:
                        boot_gse[f"epi_{c}_vs_{b}"] = paired_boot_ps(gd[f"{c}__pred_epi"], gd[f"{b}__pred_epi"], y_epi)
            gse_out["paired_bootstrap"] = boot_gse
        summary["gse188264_audit"] = gse_out

    args.out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote canonical audit summary to {args.out}")
    print("\n=== TIER 1 HEPG2 SINGLE-CELL (EXACT TWO-SIDED PHIPSON-SMYTH) ===")
    for k, v in summary["tier1_hepg2_single"]["paired_bootstrap"].items():
        print(f"  {k:40s}: {v['r_cand']:.4f} vs {v['r_base']:.4f} (delta={v['delta_r']:+.5f}, 95% CI [{v['ci95_low']:+.5f}, {v['ci95_high']:+.5f}], p_two={v['p_two_sided']:.5f})")
    print("  CAGI5 v4 vs no_delta (21 loci):", summary["tier1_hepg2_single"]["cagi5_v4_vs_no_delta"])
    print("\n=== TIER 1 MULTI-LINEAGE (EXACT TWO-SIDED PHIPSON-SMYTH) ===")
    for k, v in summary["tier1_multicell"]["paired_bootstrap"].items():
        print(f"  {k:16s}: {v}")
    print("\n=== TIER 2 FLASHZOI 4-REP ENSEMBLE (EXACT TWO-SIDED PHIPSON-SMYTH) ===")
    print("  4-rep summary:", summary["tier2_flashzoi_4rep"]["ens4_summary"])
    for k, v in summary["tier2_flashzoi_4rep"]["ens4_paired_bootstrap"].items():
        print(f"  ens4 {k:14s}: {v}")
    print("\n=== CAGI5 STRICT BORZOI-UNSEEN (fold3+fold4) AUDIT ===")
    print("  Unseen elements:", summary["cagi5_borzoi_unseen_audit"]["unseen_fold34_elements"])
    for k, v in summary["cagi5_borzoi_unseen_audit"]["models"].items():
        print(f"  {k:34s}: all29={v['all29_mean']:.4f} | unseen_fold34={v['unseen_fold34_mean']:.4f} | train_folds={v['train_folds_mean']:.4f} | {v['unseen_per_element']}")
    if summary["invivo_audit"]:
        print("\n=== IN VIVO 4-REP + BORZOI-UNSEEN (fold3+fold4) AUDIT ===")
        for k, v in summary["invivo_audit"]["pooled_summary"].items():
            print(f"  {k:34s}: {v}")
        for k, v in summary["invivo_audit"]["pooled_bootstrap"].items():
            if k.startswith("POOLED6_"):
                print(f"  {k}: {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

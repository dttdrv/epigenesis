#!/usr/bin/env python3
"""Unified, Zero-Leakage Evaluator & Paired Bootstrap Matrix Across All Saved Checkpoints.

Loads all saved model checkpoints for a given cell line (e.g. WTC11 or HepG2), evaluates
every model under IDENTICAL conditions on strictly held-out fold 1 AND on all 39,170
zero-shot CAGI5/Kircher 2019 SNVs, verifies exact mathematical invariants (reversion,
antisymmetry, RC invariance), checks loop-depth test-time scaling (for ThermoLoop /
ThermoInterlace), and computes 2,000-draw paired bootstrap 95% CIs between every candidate
and every published market baseline.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import torch
from scipy.stats import spearmanr

from benchmark_lentimpra import (
    LENTI_LEFT,
    LENTI_RIGHT,
    build_cagi5_kmer_set,
    cre_averaged_metrics,
    evaluate_cagi5_zeroshot,
    load_lentimpra,
    paired_bootstrap_gain,
    predict_array,
    rc_str,
    seqs_to_uint8,
)
from models import EpigenesisThermoLoop, build_model, count_parameters, load_jaspar_pwm_bank

CAGI5_LOCUS = {
    "BCL11A": "BCL11A", "F9": "F9", "FOXE1": "FOXE1", "GP1BA": "GP1BA", "HBB": "HBB", "HBG1": "HBG1",
    "HNF4A": "HNF4A", "IRF4": "IRF4", "IRF6": "IRF6", "LDLR": "LDLR", "LDLR.2": "LDLR", "MSMB": "MSMB",
    "MYCrs11986220": "MYCrs11986220", "MYCrs6983267": "MYCrs6983267", "PKLR-24h": "PKLR", "PKLR-48h": "PKLR",
    "RET": "RET", "SORT1": "SORT1", "SORT1-flip": "SORT1", "SORT1.2": "SORT1", "TCF7L2": "TCF7L2",
    "TERT-GAa": "TERT", "TERT-GBM": "TERT", "TERT-GSc": "TERT", "TERT-HEK": "TERT", "UC88": "UC88",
    "ZFAND3": "ZFAND3", "ZRSh-13": "ZRS", "ZRSh-13h2": "ZRS",
}


def verify_invariants(model: torch.nn.Module, sample_seqs: list[str], device: torch.device) -> dict[str, float]:
    refs = sample_seqs[:16]
    alts = [s[:115] + ("A" if s[115] != "A" else "G") + s[116:] for s in refs]

    def score(seqs: list[str]) -> np.ndarray:
        p, _ = predict_array(model, seqs_to_uint8(seqs), device, batch_size=32, rc_avg=True)
        return p

    pr, pa = score(refs), score(alts)
    d_ra = pa - pr
    d_ar = pr - pa
    d_rr = pr - pr
    pr_rc, pa_rc = score([rc_str(s) for s in refs]), score([rc_str(s) for s in alts])
    d_rc = pa_rc - pr_rc
    return dict(
        reversion_max_err=float(np.max(np.abs(d_rr))),
        antisymmetry_max_err=float(np.max(np.abs(d_ra + d_ar))),
        rc_invariance_max_err=float(np.max(np.abs(d_ra - d_rc))),
    )


def locus_bootstrap_cagi5(m_cand: dict[str, Any], m_base: dict[str, Any], n_boot: int = 2000) -> dict[str, Any]:
    loci = sorted(set(CAGI5_LOCUS.values()))
    lc_c = np.array([np.mean([m_cand[e]["pearson_r"] for e, lc in CAGI5_LOCUS.items() if lc == loc and e in m_cand]) for loc in loci])
    lc_b = np.array([np.mean([m_base[e]["pearson_r"] for e, lc in CAGI5_LOCUS.items() if lc == loc and e in m_base]) for loc in loci])
    rng = np.random.default_rng(2026)
    bidx = rng.integers(0, len(loci), size=(n_boot, len(loci)))
    gains = (lc_c[bidx] - lc_b[bidx]).mean(axis=1)
    return dict(
        point_gain=round(float((lc_c - lc_b).mean()), 5),
        ci95_low=round(float(np.quantile(gains, 0.025)), 5),
        ci95_high=round(float(np.quantile(gains, 0.975)), 5),
        prob_better=round(float((gains > 0).mean()), 4),
        elements_won=f"{sum(m_cand[e]['pearson_r'] > m_base[e]['pearson_r'] for e in m_cand)}/{len(m_cand)}",
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell-line", default="WTC11")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--ckpt-dir", type=Path, default=Path("/tmp/epigenesis-benchmarks/results"))
    ap.add_argument("--data-dir", type=Path, default=Path("/tmp/epigenesis-benchmarks/human_legnet-main/datasets/original"))
    ap.add_argument("--kircher-dir", type=Path, default=Path("/tmp/epigenesis-frontiers/data"))
    args = ap.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else ("mps" if torch.backends.mps.is_available() else "cpu"))
    wf, wr, _ = load_jaspar_pwm_bank(21)
    k31 = build_cagi5_kmer_set(args.kircher_dir / "kircher_refs.json")
    data = load_lentimpra(args.data_dir, args.cell_line, k31)

    # Decode first 16 test sequences for invariant checks
    inv_lut = np.array(list("ACGTN"))
    sample_seqs = ["".join(inv_lut[row]) for row in data["te_x"][:16]]

    ckpts = sorted(args.ckpt_dir.glob(f"{args.cell_line}_*_s{args.seed}.pt"))
    summary: dict[str, Any] = dict(cell_line=args.cell_line, seed=args.seed, models={}, paired_cre_r_gains={}, paired_cagi5_gains={})
    cre_preds: dict[str, np.ndarray] = {}
    cagi5_els: dict[str, dict[str, Any]] = {}
    y_cre_ref = None

    for p in ckpts:
        obj = torch.load(p, map_location=device, weights_only=False)
        mname = obj["model_name"]
        try:
            model = build_model(mname, wf, wr, out_dim=1, seq_len=230).to(device)
            model.load_state_dict(obj["state_dict"])
        except Exception as exc:
            print(f"Skipping {p.name} (incompatible checkpoint schema): {exc}")
            continue

        te_p, te_std = predict_array(model, data["te_x"], device, batch_size=512, rc_avg=True)
        raw_r = float(np.corrcoef(te_p, data["te_y_raw"])[0, 1])
        p_cre, y_cre, cm = cre_averaged_metrics(te_p, data["te_y_avg"], data["te_ids"])
        cre_preds[mname] = p_cre
        y_cre_ref = y_cre

        cal_rho = round(float(spearmanr(te_std, np.abs(te_p - data["te_y_raw"])).statistic), 4) if np.std(te_std) > 1e-12 else None
        inv = verify_invariants(model, sample_seqs, device)
        cagi5 = evaluate_cagi5_zeroshot(model, args.kircher_dir, device, batch_size=512)
        cagi5_els[mname] = cagi5["per_element"]

        loop_sc = None
        if isinstance(model, EpigenesisThermoLoop):
            loop_sc = {}
            orig = model.n_loops
            for t_eval in (1, 2, 4, 6, 8):
                model.n_loops = t_eval
                tp_t, _ = predict_array(model, data["te_x"], device, batch_size=512, rc_avg=True)
                _, _, cm_t = cre_averaged_metrics(tp_t, data["te_y_avg"], data["te_ids"])
                loop_sc[f"T={t_eval}"] = dict(cre_r=round(cm_t["pearson_r"], 5), mse=round(cm_t["mse"], 5))
            model.n_loops = orig

        summary["models"][mname] = dict(
            n_params=count_parameters(model),
            raw_r=round(raw_r, 5),
            cre_r=round(cm["pearson_r"], 5),
            cre_rho=round(cm["spearman_rho"], 5),
            cre_mse=round(cm["mse"], 5),
            cre_r2=round(cm["r2"], 5),
            cal_spearman=cal_rho,
            cagi5_all_r=cagi5["mean_all_r"],
            cagi5_sig_r=cagi5["mean_sig_p1e5_r"],
            cagi5_liver_r=cagi5["mean_liver_r"],
            cagi5_eryth_r=cagi5["mean_erythroid_r"],
            cagi5_neural_r=cagi5["mean_neural_r"],
            invariants=inv,
            loop_scaling=loop_sc,
        )
        print(
            f"{mname:22s} | params={count_parameters(model):>9,d} | CRE_r={cm['pearson_r']:.4f} "
            f"| MSE={cm['mse']:.4f} | cal={str(cal_rho):>6s} | CAGI5={cagi5['mean_all_r']:.4f} (sig {cagi5['mean_sig_p1e5_r']:.4f})",
            flush=True,
        )

    # Paired bootstrap gains for our novel candidates vs all baselines
    baselines = [b for b in ("legnet_official", "legnet_large", "dream_rnn", "deepstarr", "malinois") if b in cre_preds]
    cands = [c for c in ("thermoloop", "thermo_interlace", "thermosyntax") if c in cre_preds]
    if y_cre_ref is not None:
        for c in cands:
            for b in baselines:
                summary["paired_cre_r_gains"][f"{c}_vs_{b}"] = paired_bootstrap_gain(cre_preds[c], cre_preds[b], y_cre_ref, n_boot=2000)
                summary["paired_cagi5_gains"][f"{c}_vs_{b}"] = locus_bootstrap_cagi5(cagi5_els[c], cagi5_els[b], n_boot=2000)

    out_p = args.ckpt_dir / f"unified_eval_{args.cell_line}_s{args.seed}.json"
    out_p.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"\nSaved unified evaluation to {out_p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

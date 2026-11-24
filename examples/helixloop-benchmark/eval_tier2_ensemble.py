#!/usr/bin/env python3
"""Evaluate multi-replicate Tier-2 Open Apache-2.0 Flashzoi+HelixLoop-v4 Ensemble (rep0..3)
on strictly held-out ENCODE4 Fold 1 (HepG2, K562, WTC11) and TRUE ZERO-SHOT CAGI5 (39,170 SNVs).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import torch
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "studio"))
sys.path.insert(0, str(REPO / "examples/helixloop-benchmark"))

import benchmark_flashzoi as bf  # noqa: E402
from benchmark_lentimpra import LENTI_LEFT, LENTI_RIGHT, cre_averaged_metrics, paired_bootstrap_gain, rc_str, seqs_to_uint8  # noqa: E402
from benchmark_multicell import CAGI5_CELL_ROUTING, CAGI5_LOCUS, CELL_LINES, load_multicell_zero_leakage  # noqa: E402
from models import load_jaspar_pwm_bank  # noqa: E402


def load_tier2_ckpt(path: Path, wf: np.ndarray, wr: np.ndarray, device: torch.device) -> torch.nn.Module:
    blob = torch.load(path, map_location=device, weights_only=False)
    mname = blob["model_name"]
    repo_id = blob.get("repo_id", "johahi/flashzoi-replicate-0")
    n_loops = blob.get("n_loops", 6)
    if mname == "flashzoi_helixloop_v4":
        m = bf.FlashzoiHelixLoopV4(wf, wr, repo_id=repo_id, out_dim=3, n_loops=n_loops).to(device)
    elif mname == "flashzoi_linear_ft":
        m = bf.FlashzoiLinearFT(repo_id=repo_id, out_dim=3).to(device)
    else:
        raise ValueError(f"Unsupported checkpoint model_name: {mname}")
    m.load_state_dict(blob["state_dict"])
    m.eval()
    return m


def eval_cagi5_ensemble(models_list: list[torch.nn.Module], kircher_dir: Path, device: torch.device, batch_size: int = 512) -> dict[str, Any]:
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
            L, win = len(seq), 200
            if L <= win:
                sub_ref = seq.center(win, "N")
                p_rel = pos0 + (win - L) // 2
            else:
                st = max(0, min(L - win, pos0 - win // 2))
                sub_ref = seq[st : st + win]
                p_rel = pos0 - st
            by_el[el].append((sub_ref, sub_ref[:p_rel] + alt_b + sub_ref[p_rel + 1 :], float(r["Coefficient"]), float(r["pValue"])))

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
        wf_u8 = seqs_to_uint8([LENTI_LEFT + s + LENTI_RIGHT for s in u_seqs])
        wr_u8 = seqs_to_uint8([LENTI_LEFT + rc_str(s) + LENTI_RIGHT for s in u_seqs])
        p_sum = np.zeros((len(u_seqs), 3), dtype=np.float64)
        for m in models_list:
            pf, _ = bf.predict_tier2(m, wf_u8, device, batch_size=batch_size, rc_avg=False)
            pr, _ = bf.predict_tier2(m, wr_u8, device, batch_size=batch_size, rc_avg=False)
            p_sum += 0.5 * (pf + pr)
        p_both = p_sum / len(models_list)

        y = np.array([x[2] for x in rows], dtype=np.float64)
        pv = np.array([x[3] for x in rows], dtype=np.float64)
        ri = np.array([u_map[x[0]] for x in rows], dtype=np.int64)
        ai = np.array([u_map[x[1]] for x in rows], dtype=np.int64)
        dy_3 = p_both[ai] - p_both[ri]
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
    ap.add_argument("--ckpts", nargs="+", type=Path, required=True)
    ap.add_argument("--linear-ckpts", nargs="*", type=Path, default=[])
    ap.add_argument("--data-dir", type=Path, default=Path("/content/data/human_legnet-main/datasets/original"))
    ap.add_argument("--kircher-dir", type=Path, default=Path("/content/data/kircher"))
    ap.add_argument("--out", type=Path, default=Path("/content/results/benchmark_tier2_ensemble.json"))
    args = ap.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    wf, wr, _ = load_jaspar_pwm_bank(21)
    data = load_multicell_zero_leakage(args.data_dir, args.kircher_dir / "kircher_refs.json")

    def eval_group(paths: list[Path], label: str) -> tuple[dict[str, Any], dict[str, np.ndarray], dict[str, np.ndarray]]:
        ms = [load_tier2_ckpt(p, wf, wr, device) for p in paths]
        ct_res = {}
        preds_cre, ys_cre = {}, {}
        for ci, ct in enumerate(CELL_LINES):
            td = data["test_by_ct"][ct]
            p_sum = np.zeros((len(td["x"]), 3), dtype=np.float64)
            for m in ms:
                tp, _ = bf.predict_tier2(m, td["x"], device, batch_size=512, rc_avg=True)
                p_sum += tp
            p_avg = p_sum / len(ms)
            p_cre, y_cre, cm = cre_averaged_metrics(p_avg[:, ci], td["y_avg"], td["ids"])
            preds_cre[ct] = p_cre
            ys_cre[ct] = y_cre
            ct_res[ct] = {k: (round(v, 5) if isinstance(v, float) else v) for k, v in cm.items()}
        mean_r = round(float(np.mean([ct_res[ct]["pearson_r"] for ct in CELL_LINES])), 5)
        cagi5 = eval_cagi5_ensemble(ms, args.kircher_dir, device)
        print(
            f"==> ENSEMBLE {label} (k={len(paths)}): Fold1 CRE r: "
            f"HepG2={ct_res['HepG2']['pearson_r']:.4f} K562={ct_res['K562']['pearson_r']:.4f} "
            f"WTC11={ct_res['WTC11']['pearson_r']:.4f} (mean={mean_r:.4f}) | "
            f"CAGI5 ZS: blind={cagi5['mean_blind_3head_r']:.4f} matched={cagi5['mean_lineage_matched_r']:.4f} "
            f"sig={cagi5['mean_sig_p1e5_r']:.4f} liver={cagi5['mean_liver_hepg2_r']:.4f} eryth={cagi5['mean_erythroid_k562_r']:.4f}",
            flush=True,
        )
        return dict(checkpoints=[str(p.name) for p in paths], n_models=len(paths), fold1_mean_cre_r=mean_r, test_fold1_cre=ct_res, cagi5_zeroshot=cagi5), preds_cre, ys_cre

    report: dict[str, Any] = dict(
        published_alphagenome_encoder_ft_reference=dict(HepG2=0.888, K562=0.879, WTC11=0.834, mean=0.867),
        ensembles={},
    )
    res_hl, p_hl, y_ref = eval_group(args.ckpts, "flashzoi_helixloop_v4_ens")
    report["ensembles"]["flashzoi_helixloop_v4_ens"] = res_hl
    if args.linear_ckpts:
        res_lin, p_lin, _ = eval_group(args.linear_ckpts, "flashzoi_linear_ft_ens")
        report["ensembles"]["flashzoi_linear_ft_ens"] = res_lin
        gains = {}
        for ct in CELL_LINES:
            gains[f"{ct}_hl_ens_vs_lin_ens"] = paired_bootstrap_gain(p_hl[ct], p_lin[ct], y_ref[ct], n_boot=2000)
        report["paired_bootstrap_gains"] = gains
        print("Paired bootstrap gains:", gains, flush=True)

    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Saved ensemble evaluation to {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""True Zero-Shot In Vivo Human Genetics Benchmark: ADASTRA (ChIP-seq ASB) & UDACHA (DNase/ATAC-seq).

Evaluates lentiMPRA-trained models ZERO-SHOT on completely independent in vivo human allele-specific
regulatory variants from ADASTRA (TF binding) and UDACHA (DNase-seq & ATAC-seq chromatin accessibility):
  - Fetches & caches exact 200-bp GRCh38/hg38 reference windows via UCSC REST API (verifying ref base).
  - Evaluates exact RC-averaged delta = f(alt) - f(ref) against experimental allelic imbalance
    (alt_es - ref_es) using Pearson r, Spearman rho, directional AUROC, and 2,000-draw paired bootstrap CIs.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np
import torch
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "examples/helixloop-benchmark"))
import benchmark_flashzoi as bf  # noqa: E402
import benchmark_lentimpra as bl  # noqa: E402
import benchmark_multicell as bm  # noqa: E402
import models  # noqa: E402


def _fetch_one_ucsc(item: tuple[str, int, str, str, float]) -> dict[str, Any] | None:
    chrom, pos1, ref, alt, es = item
    # 1-based pos1 -> 0-based [pos1-101, pos1+99) gives 200 bp with SNP at index 100
    start0 = max(0, pos1 - 101)
    end0 = start0 + 200
    url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={chrom};start={start0};end={end0}"
    try:
        with urllib.request.urlopen(url, timeout=15) as resp:
            d = json.loads(resp.read().decode("utf-8"))
        seq = d.get("dna", "").upper()
        rel = (pos1 - 1) - start0
        if len(seq) != 200 or rel < 0 or rel >= 200 or seq[rel] != ref or any(c not in "ACGT" for c in seq):
            return None
        alt_seq = seq[:rel] + alt + seq[rel + 1 :]
        return dict(chrom=chrom, pos=pos1, ref=ref, alt=alt, es=es, ref_seq=seq, alt_seq=alt_seq)
    except Exception:
        return None


def load_or_build_invivo_panel(
    tsv_path: Path,
    kind: str,
    cache_path: Path,
    max_variants: int = 400,
) -> list[dict[str, Any]]:
    if cache_path.is_file():
        return json.loads(cache_path.read_text(encoding="utf-8"))

    cands: list[tuple[float, str, int, str, str, float]] = []
    with tsv_path.open() as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            chrom = row.get("#chr") or row.get("chr") or ""
            if not chrom.startswith("chr"):
                continue
            ref = row.get("ref", "").upper()
            alt = row.get("alt", "").upper()
            if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
                continue
            try:
                if kind == "asb":
                    pos1 = int(row["pos"])
                    fdr = min(float(row["fdrp_bh_ref"]), float(row["fdrp_bh_alt"]))
                    if fdr >= 0.01:
                        continue
                    es = float(row["es_mean_alt"]) - float(row["es_mean_ref"])
                else:
                    pos1 = int(row["end"])  # 0-based start, 1-based end
                    fdr = float(row.get("fdr_comb_pval", "1.0"))
                    if fdr >= 0.05:
                        continue
                    es = float(row["alt_comb_es"]) - float(row["ref_comb_es"])
            except (KeyError, ValueError):
                continue
            if abs(es) < 0.2:
                continue
            cands.append((fdr, chrom, pos1, ref, alt, es))

    # Balance positive (alt-preferring) and negative (ref-preferring) top-confidence SNPs
    pos_c = sorted([c for c in cands if c[5] > 0], key=lambda x: x[0])[: max_variants // 2 + 40]
    neg_c = sorted([c for c in cands if c[5] < 0], key=lambda x: x[0])[: max_variants // 2 + 40]
    selected = [(c[1], c[2], c[3], c[4], c[5]) for c in pos_c + neg_c]

    panel: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=12) as ex:
        for res in ex.map(_fetch_one_ucsc, selected):
            if res is not None:
                panel.append(res)
            if len(panel) >= max_variants:
                break

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(panel), encoding="utf-8")
    return panel


@torch.no_grad()
def score_panel(
    model: torch.nn.Module,
    panel: list[dict[str, Any]],
    device: torch.device,
    *,
    head_idx: int | None = None,
    batch_size: int = 512,
) -> np.ndarray:
    useqs: list[str] = []
    umap: dict[str, int] = {}
    for x in panel:
        for s in (x["ref_seq"], x["alt_seq"]):
            if s not in umap:
                umap[s] = len(useqs)
                useqs.append(s)
    wf = [bl.LENTI_LEFT + s + bl.LENTI_RIGHT for s in useqs]
    wr = [bl.LENTI_LEFT + bl.rc_str(s) + bl.LENTI_RIGHT for s in useqs]
    if head_idx is not None:
        pf, _ = bm.predict_multicell(model, bl.seqs_to_uint8(wf), device, batch_size=batch_size, rc_avg=False)
        pr, _ = bm.predict_multicell(model, bl.seqs_to_uint8(wr), device, batch_size=batch_size, rc_avg=False)
        pb = 0.5 * (pf[:, head_idx] + pr[:, head_idx])
    else:
        pf, _ = bl.predict_array(model, bl.seqs_to_uint8(wf), device, batch_size=batch_size, rc_avg=False)
        pr, _ = bl.predict_array(model, bl.seqs_to_uint8(wr), device, batch_size=batch_size, rc_avg=False)
        pb = 0.5 * (pf + pr)
    ri = np.array([umap[x["ref_seq"]] for x in panel], dtype=np.int64)
    ai = np.array([umap[x["alt_seq"]] for x in panel], dtype=np.int64)
    return pb[ai] - pb[ri]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", type=Path, default=Path("/content/results"))
    ap.add_argument("--datasets-dir", type=Path, default=Path("/tmp/epigenesis-benchmarks/human_legnet-main/datasets"))
    ap.add_argument("--cache-dir", type=Path, default=Path("/content/invivo_cache"))
    ap.add_argument("--checkpoints", nargs="+", required=True, help="List of label=path:model_name:out_dim:head_idx_hepg2:head_idx_k562")
    ap.add_argument("--max-variants", type=int, default=400)
    ap.add_argument("--out", type=Path, default=Path("/content/results/invivo_zeroshot.json"))
    args = ap.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    wf, wr, _ = models.load_jaspar_pwm_bank(21)

    panels_spec = [
        ("HepG2_ASB", args.datasets_dir / "asb/HepG2__hepatoblastoma_.asb.tsv", "asb", 0),
        ("HepG2_DNase", args.datasets_dir / "coverage/HepG2__hepatoblastoma_.dnase.tsv", "cov", 0),
        ("HepG2_ATAC", args.datasets_dir / "coverage/HepG2__hepatoblastoma_.atac.tsv", "cov", 0),
        ("K562_ASB", args.datasets_dir / "asb/K562__myelogenous_leukemia_.asb.tsv", "asb", 1),
        ("K562_DNase", args.datasets_dir / "coverage/K562__myelogenous_leukemia_.dnase.tsv", "cov", 1),
        ("K562_ATAC", args.datasets_dir / "coverage/K562__myelogenous_leukemia_.atac.tsv", "cov", 1),
    ]
    panels: dict[str, tuple[list[dict[str, Any]], int]] = {}
    for pname, tsv, kind, ct_i in panels_spec:
        if tsv.is_file():
            p = load_or_build_invivo_panel(tsv, kind, args.cache_dir / f"{pname}.json", max_variants=args.max_variants)
            print(f"Loaded in vivo panel {pname:12s}: n={len(p)} verified GRCh38 SNPs", flush=True)
            panels[pname] = (p, ct_i)

    report: dict[str, Any] = {"panels": {k: len(v[0]) for k, v in panels.items()}, "models": {}}
    dys_by_panel: dict[str, dict[str, np.ndarray]] = {k: {} for k in panels}

    for spec in args.checkpoints:
        label, rest = spec.split("=", 1)
        ckpt_s, mname, od_s = rest.split(":")
        out_dim = int(od_s)
        ckpt_paths = [p.strip() for p in ckpt_s.split("+") if p.strip()]
        loaded_models = []
        for cp in ckpt_paths:
            blob = torch.load(cp, map_location=device, weights_only=False)
            sd = blob["state_dict"]
            if mname == "flashzoi_helixloop_v4":
                repo_id = blob.get("repo_id", "johahi/flashzoi-replicate-0")
                m = bf.FlashzoiHelixLoopV4(wf, wr, repo_id=repo_id, out_dim=out_dim, n_loops=blob.get("n_loops", 6)).to(device)
            elif mname == "flashzoi_linear_ft":
                repo_id = blob.get("repo_id", "johahi/flashzoi-replicate-0")
                m = bf.FlashzoiLinearFT(repo_id=repo_id, out_dim=out_dim).to(device)
            else:
                m = models.build_model(mname, wf, wr, out_dim=out_dim).to(device)
            m.load_state_dict(sd)
            m.eval()
            loaded_models.append(m)

        m_res = {}
        for pname, (panel, ct_i) in panels.items():
            if out_dim == 1 and ct_i != 0:
                continue  # HepG2 single-cell model only evaluated on HepG2 panels
            h_idx = ct_i if out_dim == 3 else None
            dy = np.mean([score_panel(m, panel, device, head_idx=h_idx) for m in loaded_models], axis=0)
            dys_by_panel[pname][label] = dy
            y = np.array([x["es"] for x in panel], dtype=np.float64)
            yb = (y > 0).astype(np.int32)
            m_res[pname] = dict(
                n=len(panel),
                pearson_r=round(float(np.corrcoef(dy, y)[0, 1]), 4),
                spearman_rho=round(float(spearmanr(dy, y).statistic), 4),
                dir_auroc=round(float(roc_auc_score(yb, dy)), 4),
            )
        report["models"][label] = m_res
        print(f"==> {label:32s}: " + " | ".join(f"{k}: r={v['pearson_r']:.3f} auroc={v['dir_auroc']:.3f}" for k, v in m_res.items()), flush=True)

    # Exact two-sided Phipson-Smyth paired bootstrap p-values (B=2000) across all 6 in vivo panels!
    pairs_to_test = [
        ("Tier1_HelixLoop_v4_MC", "Tier1_LegNet_MC"),
        ("Tier2_Flashzoi_HelixLoop_v4_rep0", "Tier2_Flashzoi_Linear_rep0"),
        ("Tier2_Flashzoi_HelixLoop_v4_4rep", "Tier2_Flashzoi_Linear_4rep"),
    ]
    boot_gains = {}
    rng = np.random.default_rng(2026)
    for cand, base in pairs_to_test:
        for pname, (panel, _) in panels.items():
            if cand in dys_by_panel[pname] and base in dys_by_panel[pname]:
                y = np.array([x["es"] for x in panel], dtype=np.float64)
                dc, db = dys_by_panel[pname][cand], dys_by_panel[pname][base]
                n = len(y)
                idx = rng.integers(0, n, size=(2000, n))
                rc = np.array([np.corrcoef(dc[i], y[i])[0, 1] for i in idx])
                rb = np.array([np.corrcoef(db[i], y[i])[0, 1] for i in idx])
                g = rc - rb
                pt = float(np.corrcoef(dc, y)[0, 1] - np.corrcoef(db, y)[0, 1])
                k_opp = int((g <= 0).sum()) if pt >= 0 else int((g >= 0).sum())
                p_two = min(1.0, 2.0 * (k_opp + 1) / 2001.0)
                boot_gains[f"{pname}_{cand}_vs_{base}"] = dict(
                    point_gain=round(pt, 5),
                    ci95_low=round(float(np.quantile(g, 0.025)), 5),
                    ci95_high=round(float(np.quantile(g, 0.975)), 5),
                    p_two_sided=round(p_two, 5),
                )
    if boot_gains:
        report["paired_bootstrap_invivo"] = boot_gains

    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Saved zero-shot in vivo genetics results to {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

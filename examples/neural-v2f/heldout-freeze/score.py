"""Single-shot scorer for frozen V2F predictions on the 118 unseen GSE188264 test loci.

Verifies FREEZE.json SHA-256 digests before reading measured r.alpha labels and scoring once.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import statistics
import sys

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "examples/neural-v2f"))
sys.path.insert(0, str(REPO / "examples/neural-circuit"))
sys.path.insert(0, str(REPO / "examples/neural-regulation"))

import circuit as nc
import radar as nr
import regulation as reg


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _pearson(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 2:
        return None
    mx = statistics.fmean(xs)
    my = statistics.fmean(ys)
    dx = [x - mx for x in xs]
    dy = [y - my for y in ys]
    vx = math.fsum(a * a for a in dx)
    vy = math.fsum(b * b for b in dy)
    if vx <= 0.0 or vy <= 0.0:
        return None
    return math.fsum(a * b for a, b in zip(dx, dy)) / math.sqrt(vx * vy)


def _ranks(vals: list[float]) -> list[float]:
    indexed = sorted(enumerate(vals), key=lambda kv: kv[1])
    out = [0.0] * len(vals)
    i = 0
    while i < len(indexed):
        j = i
        while j + 1 < len(indexed) and indexed[j + 1][1] == indexed[i][1]:
            j += 1
        avg_rank = 0.5 * (i + j) + 1.0
        for k in range(i, j + 1):
            out[indexed[k][0]] = avg_rank
        i = j + 1
    return out


def _spearman(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 2:
        return None
    return _pearson(_ranks(xs), _ranks(ys))


def verify_freeze() -> dict:
    freeze_path = HERE / "FREEZE.json"
    if not freeze_path.is_file():
        raise FileNotFoundError("FREEZE.json is required before scoring")
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    if not freeze.get("sealed_before_outcomes_scored"):
        raise ValueError("FREEZE.json must declare sealed_before_outcomes_scored=true")
    files = freeze.get("files")
    if not isinstance(files, dict) or not files:
        raise ValueError("FREEZE.json must map frozen file digests")
    for rel, expected in files.items():
        path = HERE / rel
        if not path.is_file():
            raise FileNotFoundError(f"missing frozen file: {rel}")
        actual = _sha256(path)
        if actual != expected:
            raise ValueError(f"digest mismatch for {rel}: expected {expected}, got {actual}")
    return freeze


def main() -> int:
    freeze = verify_freeze()
    preds = json.loads((HERE / "predictions.json").read_text(encoding="utf-8"))

    rows, _ = reg.load_data()
    curated = set(nc.CURATED_LOCI)
    measured_by_id = {
        r.id: r
        for r in rows
        if r.split == "test" and r.locus not in curated and r.kind != "wild-type" and r.sequence != r.reference
    }

    loci = preds["loci"]
    variants = preds["variants"]
    loci_with_cryptic = [loc for loc in loci if loc["category_counts"]["cryptic_bifurcation"] > 0]
    total_snv_cryptic = sum(loc["category_counts"]["cryptic_bifurcation"] for loc in loci)
    total_snv_doubly_silent = sum(loc["doubly_silent_count"] for loc in loci)

    # Check H1: null control on top unseen cryptic locus (if any)
    top_unseen = None
    top_unseen_null = None
    if loci_with_cryptic:
        top_loc = max(
            loci_with_cryptic,
            key=lambda x: (x["doubly_silent_count"], x["category_counts"]["cryptic_bifurcation"], x["wt_peak_weight"]),
        )
        top_unseen = top_loc["locus"]
        ref_seq = next(r.reference for r in rows if r.locus == top_unseen)
        ctx = nc._resolve_context(
            None, None, None, 0.5 * (nr.DEFAULT_SIGNAL_MIN + nr.DEFAULT_SIGNAL_MAX), None, None
        )
        eff_floor = nc._resolve_background_floor(None, nc.DEFAULT_EPSILON)
        top_scan = nr.scan_enhancer_regime(
            ref_seq,
            signal_min=nr.DEFAULT_SIGNAL_MIN,
            signal_max=nr.DEFAULT_SIGNAL_MAX,
            signal_steps=nr.DEFAULT_SIGNAL_STEPS,
            epsilon=nc.DEFAULT_EPSILON,
            background_floor=eff_floor,
            include_sequence_models=False,
            ctx=ctx,
        )
        nulls = nr.run_null_controls(
            ref_seq,
            top_scan,
            ctx=ctx,
            signal_min=nr.DEFAULT_SIGNAL_MIN,
            signal_max=nr.DEFAULT_SIGNAL_MAX,
            signal_steps=nr.DEFAULT_SIGNAL_STEPS,
            n_shuffles=100,
            epsilon=nc.DEFAULT_EPSILON,
            background_floor=eff_floor,
            seed=17,
        )
        top_unseen_null = {
            "locus": top_unseen,
            "true_cryptic_count": top_loc["category_counts"]["cryptic_bifurcation"],
            "n_shuffles": 100,
            "dinucleotide_null_mean": nulls["dinucleotide_shuffle"]["null_cryptic_mean"],
            "dinucleotide_null_max": nulls["dinucleotide_shuffle"]["null_cryptic_max"],
            "empirical_p_ge": nulls["dinucleotide_shuffle"]["empirical_p_ge"],
            "scrambled_pwm_cryptic": nulls["motif_pwm_scramble"]["scrambled_cryptic_count"],
        }

    # Check H2 & H3: score every predicted test variant against measured r.alpha
    scored_rows = []
    by_cat: dict[str, list[dict]] = {
        "cryptic_bifurcation": [],
        "concordant_switch": [],
        "canalized_buffered": [],
        "neutral": [],
    }
    for v in variants:
        r = measured_by_id.get(v["id"])
        if r is None:
            raise ValueError(f"predicted variant {v['id']} absent from test split")
        alphas = list(r.alpha)
        meas_mean = statistics.fmean(alphas)
        meas_72h = alphas[-1]
        meas_max_abs = max(abs(a) for a in alphas)
        rec = {
            **v,
            "measured_alpha_mean": round(meas_mean, 6),
            "measured_alpha_72h": round(meas_72h, 6),
            "measured_alpha_max_abs": round(meas_max_abs, 6),
            "empirically_silent_015": abs(meas_mean) < 0.15,
        }
        scored_rows.append(rec)
        by_cat[v["category"]].append(rec)

    abs_rep = [r["max_abs_reporter_delta"] for r in scored_rows]
    abs_meas = [abs(r["measured_alpha_mean"]) for r in scored_rows]
    abs_kreimer = [abs(r["kreimer_mpra_mean"]) for r in scored_rows if r["kreimer_mpra_mean"] is not None]
    meas_for_kreimer = [abs(r["measured_alpha_mean"]) for r in scored_rows if r["kreimer_mpra_mean"] is not None]
    signed_kreimer = [r["kreimer_mpra_mean"] for r in scored_rows if r["kreimer_mpra_mean"] is not None]
    signed_meas = [r["measured_alpha_mean"] for r in scored_rows if r["kreimer_mpra_mean"] is not None]

    cat_summary = {}
    for cat, items in by_cat.items():
        if items:
            cat_summary[cat] = {
                "count": len(items),
                "mean_abs_measured_alpha": round(statistics.fmean(abs(x["measured_alpha_mean"]) for x in items), 6),
                "median_abs_measured_alpha": round(statistics.median(abs(x["measured_alpha_mean"]) for x in items), 6),
                "empirically_silent_fraction": round(sum(x["empirically_silent_015"] for x in items) / len(items), 6),
            }
        else:
            cat_summary[cat] = {
                "count": 0,
                "mean_abs_measured_alpha": None,
                "median_abs_measured_alpha": None,
                "empirically_silent_fraction": None,
            }

    results = {
        "format": "epigenesis.neural-v2f-heldout-results/v1",
        "biological_acceptance": False,
        "freeze_verified": True,
        "frozen_digests": freeze["files"],
        "unseen_test_loci": len(loci),
        "unseen_loci_with_snv_cryptic": len(loci_with_cryptic),
        "total_snv_cryptic_across_unseen_loci": total_snv_cryptic,
        "total_snv_doubly_silent_across_unseen_loci": total_snv_doubly_silent,
        "top_unseen_locus_null_control": top_unseen_null,
        "measured_test_variants_scored": len(scored_rows),
        "measured_category_summary": cat_summary,
        "correlations": {
            "pearson_abs_reporter_vs_abs_alpha": round(_pearson(abs_rep, abs_meas), 6) if _pearson(abs_rep, abs_meas) is not None else None,
            "spearman_abs_reporter_vs_abs_alpha": round(_spearman(abs_rep, abs_meas), 6) if _spearman(abs_rep, abs_meas) is not None else None,
            "pearson_kreimer_vs_measured_alpha": round(_pearson(signed_kreimer, signed_meas), 6) if _pearson(signed_kreimer, signed_meas) is not None else None,
            "spearman_abs_kreimer_vs_abs_alpha": round(_spearman(abs_kreimer, meas_for_kreimer), 6) if _spearman(abs_kreimer, meas_for_kreimer) is not None else None,
        },
        "measured_cryptic_variants": by_cat["cryptic_bifurcation"],
    }
    (HERE / "results.json").write_text(
        json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(json.dumps({k: v for k, v in results.items() if k != "measured_cryptic_variants"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

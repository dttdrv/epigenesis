"""Generate frozen pre-scoring V2F predictions on the 118 unseen GSE188264 test loci.

Reads sequence and metadata only (never touches r.alpha or r.wt_alphas) and writes
predictions.json for cryptographic hashing before single-shot scoring.
"""

from __future__ import annotations

import json
from pathlib import Path
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


def main() -> int:
    rows, _ = reg.load_data()
    curated = set(nc.CURATED_LOCI)
    test_rows = [r for r in rows if r.split == "test" and r.locus not in curated]
    by_locus: dict[str, list] = {}
    for r in test_rows:
        by_locus.setdefault(r.locus, []).append(r)

    loci = sorted(by_locus)
    ctx = nc._resolve_context(
        None, None, None, 0.5 * (nr.DEFAULT_SIGNAL_MIN + nr.DEFAULT_SIGNAL_MAX), None, None
    )
    eff_floor = nc._resolve_background_floor(None, nc.DEFAULT_EPSILON)
    signal_grid = nr._signal_grid(
        nr.DEFAULT_SIGNAL_MIN, nr.DEFAULT_SIGNAL_MAX, nr.DEFAULT_SIGNAL_STEPS
    )
    bindings = ctx["bindings"]
    bound_genes = sorted(bindings)

    locus_predictions = []
    variant_predictions = []
    cryptic_pairs = []
    variant_pairs = []

    for locus in loci:
        lrows = by_locus[locus]
        ref = lrows[0].reference
        scan = nr.scan_enhancer_regime(
            ref,
            signal_min=nr.DEFAULT_SIGNAL_MIN,
            signal_max=nr.DEFAULT_SIGNAL_MAX,
            signal_steps=nr.DEFAULT_SIGNAL_STEPS,
            epsilon=nc.DEFAULT_EPSILON,
            background_floor=eff_floor,
            include_sequence_models=False,
            ctx=ctx,
        )
        for sw in scan["cryptic_switches"]:
            cryptic_pairs.append({
                "id": f"{locus}::{sw['id']}",
                "reference": ref,
                "sequence": sw["variant"],
            })
        locus_predictions.append({
            "locus": locus,
            "group": lrows[0].group,
            "length": len(ref),
            "wt_critical_signal": scan["wt_critical_signal"],
            "wt_peak_weight": round(scan["wt_partition"]["peak_site"]["weight"], 6),
            "wt_partition_sum": round(scan["wt_partition"]["partition_sum"], 6),
            "category_counts": scan["category_counts"],
            "doubly_silent_count": 0,
            "destabilized_count": scan["destabilized_count"],
            "cryptic_switches": [
                {
                    "id": sw["id"],
                    "affinity_ratio": sw["affinity_ratio"],
                    "max_abs_reporter_delta": sw["max_abs_reporter_delta"],
                    "max_abs_target_delta": sw["max_abs_target_delta"],
                    "regime_flip_fraction": sw["regime_flip_fraction"],
                    "delta_s_crit": sw["delta_s_crit"],
                }
                for sw in scan["cryptic_switches"]
            ],
            "top_cryptic": None,
        })

        wt_base = {
            "wt_profile": scan["wt_profile"],
            "wt_transitions": scan["wt_transitions"],
            "wt_primary_critical": scan["wt_critical_signal"],
            "wt_preferred_fate": scan["wt_preferred_fate"],
        }
        cache: dict[tuple, dict] = {}
        seen_ids: set[str] = set()
        for r in lrows:
            if r.kind == "wild-type" or r.sequence == ref or r.id in seen_ids:
                continue
            seen_ids.add(r.id)
            q_by_tf = nc.multi_tf_affinity(
                ref, r.sequence, dict(ctx["unique_specs"]), nc.DEFAULT_EPSILON, background_floor=eff_floor
            )["ratios"]
            ratio_map = {g: q_by_tf[nc._motif_spec_label(bindings[g]["protein"])] for g in bound_genes}
            reg_eval = nr._evaluate_regime_for_ratios(
                ctx, ratio_map, wt_base["wt_profile"], wt_base["wt_primary_critical"],
                wt_base["wt_preferred_fate"], signal_grid, cache, bound_genes=bound_genes,
            )
            variant_pairs.append({"id": r.id, "reference": ref, "sequence": r.sequence})
            variant_predictions.append({
                "id": r.id,
                "locus": locus,
                "group": r.group,
                "kind": r.kind,
                "edit_kind": nc._edit_kind(ref, r.sequence),
                "affinity_ratio": round(q_by_tf[ctx["primary_label"]], 6),
                "max_abs_reporter_delta": reg_eval["max_abs_reporter_delta"],
                "regime_flip_fraction": reg_eval["regime_flip_fraction"],
                "max_abs_target_delta": reg_eval["max_abs_target_delta"],
                "delta_s_crit": reg_eval["delta_s_crit"],
                "signal_locked": reg_eval["signal_locked"],
                "category": reg_eval["category"],
                "siraj_score": None,
                "kreimer_mpra_mean": None,
                "kreimer_mpra_72h": None,
                "sequence_silent": None,
                "doubly_silent": False,
            })

    all_pairs = cryptic_pairs + variant_pairs
    s_sc, k_sc = nr._optional_sequence_scores(all_pairs)
    for loc in locus_predictions:
        switches = loc.pop("cryptic_switches")
        scored_sw = []
        for sw in switches:
            cid = f"{loc['locus']}::{sw['id']}"
            s_val = round(s_sc[cid], 6) if s_sc is not None else None
            k_mean = round(k_sc[cid]["mean"], 6) if k_sc is not None else None
            d_sil = bool(k_mean is not None and abs(k_mean) < nr.DEFAULT_SEQUENCE_SILENT)
            scored_sw.append({**sw, "siraj_score": s_val, "kreimer_mpra_mean": k_mean, "doubly_silent": d_sil})
        loc["doubly_silent_count"] = sum(1 for x in scored_sw if x["doubly_silent"])
        if scored_sw:
            scored_sw.sort(key=lambda e: (not e["doubly_silent"], -e["regime_flip_fraction"], -e["max_abs_target_delta"], e["max_abs_reporter_delta"], e["id"]))
            loc["top_cryptic"] = scored_sw[0]

    for item in variant_predictions:
        if s_sc is not None:
            item["siraj_score"] = round(s_sc[item["id"]], 6)
        if k_sc is not None:
            item["kreimer_mpra_mean"] = round(k_sc[item["id"]]["mean"], 6)
            item["kreimer_mpra_72h"] = round(k_sc[item["id"]]["h72"], 6)
            seq_silent = abs(item["kreimer_mpra_mean"]) < nr.DEFAULT_SEQUENCE_SILENT
            item["sequence_silent"] = seq_silent
            item["doubly_silent"] = bool(item["category"] == "cryptic_bifurcation" and seq_silent)

    payload = {
        "format": "epigenesis.neural-v2f-heldout-predictions/v1",
        "biological_acceptance": False,
        "sealed_before_outcomes_scored": True,
        "curated_loci_excluded": sorted(curated),
        "test_loci_count": len(locus_predictions),
        "test_variants_count": len(variant_predictions),
        "background_floor": eff_floor,
        "signal_grid": signal_grid,
        "sequence_models_attached": k_sc is not None,
        "loci": locus_predictions,
        "variants": variant_predictions,
    }
    (HERE / "predictions.json").write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(f"Wrote {len(locus_predictions)} unseen loci and {len(variant_predictions)} variants to predictions.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

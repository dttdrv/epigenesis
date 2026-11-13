"""Variant-to-Function (V2F) Tipping-Point & Rescue Radar for non-coding regulatory DNA.

Identifies "passive-reporter silent, endogenous-circuit tipping" non-coding variants
across full developmental signaling regimes (SNVs, MNVs, insertions, deletions),
quantifies morphogen tipping-point shifts (Delta S_crit) via exact bisection,
searches for dual rescues (second-site cis epistasis and trans morphogen dose shifts),
evaluates statistical null controls (Altschul-Erickson dinucleotide shuffles and
column-scrambled motif PWMs), optionally attaches passive Kreimer 0-72h and Siraj
sequence predictors when NumPy is available, and seals top candidates in verified
brainc Development Modules.

All biological development claims remain conditional (biological_acceptance = False).
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import html
import json
import math
from pathlib import Path
import random
import statistics
import sys

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "examples/neural-circuit"))
sys.path.insert(0, str(REPO / "examples/neural-patterning"))

from brainc._io import load_json_object, read_regular_file
import circuit as nc
import patterning as pat

DEFAULT_SIGNAL_MIN = 0.06
DEFAULT_SIGNAL_MAX = 0.85
DEFAULT_SIGNAL_STEPS = 9
DEFAULT_REPORTER_SILENT = 0.02
DEFAULT_SEQUENCE_SILENT = 0.15
DEFAULT_SHUFFLES = 10
DEFAULT_TOP_K = 10
BISECTION_STEPS = 12
SCAN_DT = 0.25


def _signal_grid(signal_min: float, signal_max: float, steps: int) -> list[float]:
    if not (0.0 <= signal_min < signal_max <= 1.0):
        raise ValueError("signal bounds must satisfy 0.0 <= signal_min < signal_max <= 1.0")
    if steps < 3:
        raise ValueError("signal steps must be >= 3")
    return [
        round(signal_min + (signal_max - signal_min) * i / (steps - 1), 6)
        for i in range(steps)
    ]


def dinucleotide_shuffle(sequence: str, seed: int = 0) -> str:
    """Exact Altschul-Erickson Eulerian dinucleotide-preserving DNA shuffle in pure stdlib."""
    seq = sequence.upper()
    if set(seq) - set("ACGT"):
        raise ValueError("dinucleotide_shuffle requires uppercase ACGT characters")
    if len(seq) <= 2:
        return seq

    adj_orig: dict[str, list[str]] = {b: [] for b in "ACGT"}
    for i in range(len(seq) - 1):
        adj_orig[seq[i]].append(seq[i + 1])

    start_base = seq[0]
    end_base = seq[-1]
    active = [b for b in "ACGT" if adj_orig[b] or b == end_base]
    rng = random.Random(seed)

    # Rejection-sample a spanning arborescence routed into end_base (Altschul & Erickson 1985)
    for _ in range(512):
        adj = {b: list(adj_orig[b]) for b in "ACGT"}
        last_edge: dict[str, str] = {}
        for b in active:
            if b != end_base and adj[b]:
                idx = rng.randrange(len(adj[b]))
                last_edge[b] = adj[b].pop(idx)

        # Verify every active vertex != end_base reaches end_base via last_edge without cycles
        reachable = True
        for b in active:
            if b == end_base or not adj_orig[b]:
                continue
            cur = b
            seen = set()
            while cur != end_base:
                if cur in seen or cur not in last_edge:
                    reachable = False
                    break
                seen.add(cur)
                cur = last_edge[cur]
            if not reachable:
                break

        if not reachable:
            continue

        for b in active:
            rng.shuffle(adj[b])
            if b in last_edge:
                adj[b].append(last_edge[b])

        out = [start_base]
        cur = start_base
        pos_map = {b: 0 for b in "ACGT"}
        while pos_map[cur] < len(adj[cur]):
            nxt = adj[cur][pos_map[cur]]
            pos_map[cur] += 1
            out.append(nxt)
            cur = nxt
        shuffled = "".join(out)
        if len(shuffled) == len(seq):
            return shuffled

    raise RuntimeError("failed to construct Eulerian dinucleotide shuffle")


def verify_dinucleotide_counts(a: str, b: str) -> bool:
    return (
        len(a) == len(b)
        and Counter(a) == Counter(b)
        and Counter(zip(a, a[1:])) == Counter(zip(b, b[1:]))
    )


def scramble_motif_pwm(protein: str | dict, seed: int = 0) -> dict:
    """Column-permuted PWM negative control preserving column composition but destroying motif word order."""
    prof = nc.load_motif_profile(protein, epsilon=0.0)
    width = prof["width"]
    order = list(range(width))
    if width >= 2:
        rng = random.Random(seed)
        cols = [tuple(round(prof["ratios"][b][i], 10) for b in "ACGT") for i in order]
        if len(set(cols)) > 1:
            for _ in range(16):
                rng.shuffle(order)
                if order != list(range(width)):
                    break
    matrix = {b: [float(prof["ratios"][b][i]) for i in order] for b in "ACGT"}
    return {
        "id": f"{prof['id']}-scrambled",
        "kind": "weights",
        "matrix": matrix,
    }


def _solve_at_signal(
    ctx: dict,
    ratio_map: dict[str, float],
    signal: float,
    *,
    duration: float = 40.0,
    dt: float = SCAN_DT,
    initial: tuple[float, ...] | None = None,
) -> tuple[tuple[float, ...], str, float]:
    model = ctx["model"]
    genes = ctx["genes"]
    eff_target = ctx["target"]
    rep_idx = nc._target_gene_index(genes, eff_target)
    init_t = initial if initial is None else tuple(initial)
    solve_cache = ctx.setdefault("_solve_cache", {})
    ckey = None
    if initial is None:
        ckey = (
            tuple((g, round(ratio_map[g], 6)) for g in sorted(ratio_map)),
            round(float(signal), 8),
            round(float(duration), 6),
            round(float(dt), 6),
        )
        hit = solve_cache.get(ckey)
        if hit is not None:
            return hit

    end_state = nc.simulate_circuit(
        [(duration, float(signal))],
        ratios=ratio_map,
        target=eff_target,
        dt=dt,
        initial=(0.0,) * len(genes) if init_t is None else init_t,
        circuit=ctx["circuit_arg"],
        roles=ctx["roles"],
        mechanisms=ctx["mechanisms"],
    )[-1][1:]
    fate = nc.fate_label(end_state, genes, fates=model.get("fates"))
    rep = nc.occupancy_ratios(
        model, end_state, float(signal), ratio_map, roles=ctx["roles"], mechanisms=ctx["mechanisms"]
    )[rep_idx]
    res = (end_state, fate, float(rep))
    if ckey is not None:
        solve_cache[ckey] = res
    return res


def find_critical_transitions(
    ctx: dict,
    ratio_map: dict[str, float],
    *,
    duration: float = 40.0,
    coarse_steps: int = 21,
    seed_points: list[tuple[float, str]] | None = None,
) -> list[dict]:
    """Locate cell-fate transition thresholds along S in [0, 1] via grid + bisection."""
    if seed_points:
        pts: dict[float, str] = {round(float(s), 6): f for s, f in seed_points}
        for edge in (0.0, 1.0):
            if edge not in pts:
                pts[edge] = _solve_at_signal(ctx, ratio_map, edge, duration=duration)[1]
        grid = sorted(pts)
        fates = [pts[s] for s in grid]
    else:
        grid = [round(i / (coarse_steps - 1), 6) for i in range(coarse_steps)]
        fates = [_solve_at_signal(ctx, ratio_map, s, duration=duration)[1] for s in grid]
    transitions = []
    for idx in range(len(grid) - 1):
        if fates[idx] != fates[idx + 1]:
            lo_s, hi_s = grid[idx], grid[idx + 1]
            lo_f, hi_f = fates[idx], fates[idx + 1]
            for _ in range(BISECTION_STEPS):
                mid = 0.5 * (lo_s + hi_s)
                mid_f = _solve_at_signal(ctx, ratio_map, mid, duration=duration)[1]
                if mid_f == lo_f:
                    lo_s = mid
                else:
                    hi_s = mid
            transitions.append({
                "from_fate": lo_f,
                "to_fate": hi_f,
                "signal": round(0.5 * (lo_s + hi_s), 6),
                "bracket": [lo_s, hi_s],
            })
    return transitions


def _choose_primary_transition(
    transitions: list[dict],
    preferred_fate: str | None = None,
    *,
    allow_from: bool = False,
) -> dict | None:
    if not transitions:
        return None
    if preferred_fate is not None:
        for tr in transitions:
            if tr["to_fate"] == preferred_fate:
                return tr
        if allow_from:
            for tr in transitions:
                if tr["from_fate"] == preferred_fate:
                    return tr
        return None
    return transitions[-1]


def _optional_sequence_scores(pairs: list[dict]) -> tuple[dict | None, dict | None]:
    """Evaluate bundled Kreimer 0-72h and Siraj Huber-linear sequence predictors if NumPy is installed."""
    if not pairs:
        return None, None
    try:
        import numpy as np  # noqa: F401
    except ImportError:
        return None, None

    siraj_scores = None
    kreimer_scores = None
    try:
        model_dir = str(REPO / "studio/epigenesis_studio/model")
        studio_dir = str(REPO / "studio")
        if studio_dir not in sys.path:
            sys.path.insert(0, studio_dir)
        if model_dir not in sys.path:
            sys.path.insert(0, model_dir)
        from epigenesis_studio.model import predict as siraj_predict
        from epigenesis_studio import kreimer as kreimer_mod

        siraj_scores = {r["id"]: float(r["effect"]) for r in siraj_predict.predict(pairs)}
        kreimer_scores = {
            r["id"]: {"mean": float(r["effect"]), "h72": float(r["temporal"][-1])}
            for r in kreimer_mod.predict(pairs)
        }
    except (ImportError, ValueError, OSError, KeyError):
        return siraj_scores, kreimer_scores
    return siraj_scores, kreimer_scores


def find_trans_signal_rescue(
    ctx: dict,
    ratio_map: dict[str, float],
    fail_signal: float,
    desired_fate: str,
    *,
    duration: float = 40.0,
) -> dict:
    """Find minimal morphogen signal dose shift Delta S_rescue restoring desired WT fate."""
    grid = [round(i / 20, 4) for i in range(21)]
    grid_fates = [_solve_at_signal(ctx, ratio_map, s, duration=duration)[1] for s in grid]
    matching = [s for s, f in zip(grid, grid_fates) if f == desired_fate]
    if not matching:
        return {
            "rescued": False,
            "desired_fate": desired_fate,
            "fail_signal": fail_signal,
            "rescue_signal": None,
            "signal_shift": None,
            "reason": "no signal in [0, 1] restores wild-type fate",
        }

    nearest = min(matching, key=lambda s: (abs(s - fail_signal), s))
    # Bracket toward fail_signal for minimal threshold shift
    step = 0.05
    if nearest > fail_signal:
        lo_s = max(0.0, nearest - step)
        hi_s = nearest
    elif nearest < fail_signal:
        lo_s = nearest
        hi_s = min(1.0, nearest + step)
    else:
        lo_s = hi_s = nearest

    lo_f = _solve_at_signal(ctx, ratio_map, lo_s, duration=duration)[1]
    hi_f = _solve_at_signal(ctx, ratio_map, hi_s, duration=duration)[1]
    if lo_f == desired_fate and hi_f != desired_fate:
        good_s, bad_s = lo_s, hi_s
    elif hi_f == desired_fate and lo_f != desired_fate:
        good_s, bad_s = hi_s, lo_s
    else:
        good_s, bad_s = nearest, nearest

    if good_s != bad_s:
        for _ in range(12):
            mid = 0.5 * (good_s + bad_s)
            if _solve_at_signal(ctx, ratio_map, mid, duration=duration)[1] == desired_fate:
                good_s = mid
            else:
                bad_s = mid

    # Step 1e-3 past boundary on the rescuing side so steady state is squarely inside desired attractor
    eps = 1e-3 if good_s >= bad_s else -1e-3
    safe_s = min(1.0, max(0.0, round(good_s + eps, 6)))
    safe_fate = _solve_at_signal(ctx, ratio_map, safe_s, duration=duration)[1]
    if safe_fate != desired_fate:
        safe_s = round(good_s, 6)
        safe_fate = _solve_at_signal(ctx, ratio_map, safe_s, duration=duration)[1]

    return {
        "rescued": safe_fate == desired_fate,
        "desired_fate": desired_fate,
        "fail_signal": round(fail_signal, 6),
        "rescue_signal": safe_s,
        "signal_shift": round(safe_s - fail_signal, 6),
        "rescued_fate": safe_fate,
    }


def _evaluate_regime_for_ratios(
    ctx: dict,
    ratio_map: dict[str, float],
    wt_profile: list[dict],
    wt_primary_crit: float | None,
    wt_preferred_fate: str | None,
    signal_grid: list[float],
    cache: dict[tuple, dict],
    *,
    bound_genes: list[str],
    duration: float = 40.0,
    reporter_threshold: float = DEFAULT_REPORTER_SILENT,
    compute_transitions: bool = True,
    quantize_for_null: bool = False,
) -> dict:
    model = ctx["model"]
    genes = ctx["genes"]
    eff_target = ctx["target"]
    rep_idx = nc._target_gene_index(genes, eff_target)
    if quantize_for_null:
        def _qbin(q: float) -> float:
            step = 0.025
            b = round(math.log10(max(q, 1e-12)) / step) * step
            return round(10.0 ** b, 6)
        q_ratios = {g: _qbin(ratio_map[g]) for g in bound_genes}
        key = tuple(q_ratios[g] for g in bound_genes)
        cached = cache.get(key)
        if cached is not None:
            return cached
        eval_ratios = q_ratios
    else:
        key = tuple(round(ratio_map[g], 6) for g in bound_genes)
        cached = cache.get(key)
        if cached is not None:
            return cached
        eval_ratios = ratio_map

    per_signal = []
    for w_row in wt_profile:
        s = w_row["signal"]
        end_state, mut_fate, _ = _solve_at_signal(ctx, eval_ratios, s, duration=duration)
        wt_end = tuple(w_row["end_state"])
        # Passive reporter occupancy is assayed on WT endogenous TF state with mutant cis-ratio
        mut_rep = nc.occupancy_ratios(
            model, wt_end, s, eval_ratios, roles=ctx["roles"], mechanisms=ctx["mechanisms"]
        )[rep_idx]
        d_rep = float(mut_rep - w_row["reporter"])
        d_target = float(end_state[rep_idx] - wt_end[rep_idx])
        per_signal.append({
            "signal": s,
            "wt_fate": w_row["fate"],
            "mut_fate": mut_fate,
            "flipped": mut_fate != w_row["fate"],
            "d_rep": d_rep,
            "d_target": d_target,
            "mut_target": float(end_state[rep_idx]),
            "wt_target": float(wt_end[rep_idx]),
            "mut_state": [float(v) for v in end_state],
        })

    max_abs_rep = max(abs(r["d_rep"]) for r in per_signal)
    max_abs_target = max(abs(r["d_target"]) for r in per_signal)
    flipped_rows = [r for r in per_signal if r["flipped"]]

    # If using quantized null bins and this bin flips fate, refine on exact ratio_map
    if quantize_for_null and flipped_rows:
        return _evaluate_regime_for_ratios(
            ctx, ratio_map, wt_profile, wt_primary_crit, wt_preferred_fate, signal_grid,
            ctx.setdefault("_null_exact_cache", {}),
            bound_genes=bound_genes, duration=duration, reporter_threshold=reporter_threshold,
            compute_transitions=compute_transitions, quantize_for_null=False,
        )

    flip_fraction = len(flipped_rows) / len(signal_grid)
    point_silent_flip = any(abs(r["d_rep"]) < reporter_threshold for r in flipped_rows)
    worst_row = max(
        flipped_rows if flipped_rows else per_signal,
        key=lambda r: (abs(r["d_target"]), -abs(r["d_rep"])),
    )

    passive_silent = max_abs_rep < reporter_threshold
    if passive_silent and flipped_rows:
        category = "cryptic_bifurcation"
    elif (not passive_silent) and flipped_rows:
        category = "concordant_switch"
    elif (not passive_silent) and not flipped_rows:
        category = "canalized_buffered"
    else:
        category = "neutral"

    var_crit = None
    delta_scrit = None
    signal_locked = False
    if compute_transitions:
        if flipped_rows or max_abs_target >= 0.05 or any(abs(q - 1.0) >= 0.02 for q in eval_ratios.values()):
            seed_pts = [(r["signal"], r["mut_fate"]) for r in per_signal]
            var_transitions = find_critical_transitions(ctx, eval_ratios, duration=duration, seed_points=seed_pts)
            tr = _choose_primary_transition(var_transitions, wt_preferred_fate)
            if tr is not None:
                var_crit = tr["signal"]
                if wt_primary_crit is not None:
                    delta_scrit = round(var_crit - wt_primary_crit, 6)
            elif wt_primary_crit is not None:
                signal_locked = True
        else:
            var_crit = wt_primary_crit
            delta_scrit = 0.0 if wt_primary_crit is not None else None

    summary = {
        "max_abs_reporter_delta": round(max_abs_rep, 8),
        "max_abs_target_delta": round(max_abs_target, 8),
        "amplification_ratio": round(max_abs_target / max(max_abs_rep, 1e-12), 4),
        "passive_silent": passive_silent,
        "point_silent_flip": point_silent_flip,
        "regime_flip_fraction": round(flip_fraction, 6),
        "flipped_signal_count": len(flipped_rows),
        "flipped_signals": [r["signal"] for r in flipped_rows],
        "worst_signal": worst_row["signal"],
        "worst_wt_fate": worst_row["wt_fate"],
        "worst_mut_fate": worst_row["mut_fate"],
        "worst_d_rep": round(worst_row["d_rep"], 8),
        "worst_d_target": round(worst_row["d_target"], 8),
        "category": category,
        "var_critical_signal": var_crit,
        "delta_s_crit": delta_scrit,
        "signal_locked": signal_locked,
        "per_signal": per_signal,
        "_ratio_key": key,
    }
    cache[key] = summary
    return summary


def _build_wt_baseline(ctx: dict, signal_grid: list[float], bound_genes: list[str], *, duration: float = 40.0) -> dict:
    genes = ctx["genes"]
    eff_target = ctx["target"]
    rep_idx = nc._target_gene_index(genes, eff_target)
    wt_ratios = dict.fromkeys(bound_genes, 1.0)
    wt_profile = []
    for s in signal_grid:
        end_state, fate, rep = _solve_at_signal(ctx, wt_ratios, s, duration=duration)
        wt_profile.append({
            "signal": s,
            "fate": fate,
            "reporter": rep,
            "target_level": float(end_state[rep_idx]),
            "end_state": [float(v) for v in end_state],
        })
    wt_transitions = find_critical_transitions(ctx, wt_ratios, duration=duration)
    preferred_fate = eff_target if (isinstance(eff_target, str) and eff_target in genes) else None
    wt_primary = _choose_primary_transition(wt_transitions, preferred_fate, allow_from=True)
    return {
        "wt_ratios": wt_ratios,
        "wt_profile": wt_profile,
        "wt_transitions": wt_transitions,
        "wt_primary_critical": wt_primary["signal"] if wt_primary else None,
        "wt_preferred_fate": wt_primary["to_fate"] if wt_primary else preferred_fate,
    }


def scan_enhancer_regime(
    reference: str,
    *,
    circuit: str | dict | Path | None = None,
    protein: str | dict | Path | None = None,
    target: str | list[str] | tuple[str, ...] | None = None,
    signal_min: float = DEFAULT_SIGNAL_MIN,
    signal_max: float = DEFAULT_SIGNAL_MAX,
    signal_steps: int = DEFAULT_SIGNAL_STEPS,
    epsilon: float = nc.DEFAULT_EPSILON,
    background_floor: float | None = None,
    duration: float = 40.0,
    reporter_threshold: float = DEFAULT_REPORTER_SILENT,
    sequence_threshold: float = DEFAULT_SEQUENCE_SILENT,
    include_sequence_models: bool = True,
    compute_transitions: bool = True,
    quantize_for_null: bool = False,
    ctx: dict | None = None,
) -> dict:
    """Scan all 3*L SNVs across a full signaling regime [signal_min, signal_max]."""
    eff_floor = nc._resolve_background_floor(background_floor, epsilon)
    if ctx is None:
        ctx = nc._resolve_context(circuit, protein, target, 0.5 * (signal_min + signal_max), None, None)
    signal_grid = _signal_grid(signal_min, signal_max, signal_steps)
    bindings = ctx["bindings"]
    bound_genes = sorted(bindings)

    wt_base = _build_wt_baseline(ctx, signal_grid, bound_genes, duration=duration)
    wt_profile = wt_base["wt_profile"]

    tf_tables: dict[str, dict] = {}
    for lbl, spec in ctx["unique_specs"]:
        prof = nc.load_motif_profile(spec, epsilon)
        wt_part = nc.boltzmann_affinity(reference, spec, epsilon)
        tf_tables[lbl] = {"width": prof["width"], "ratios": prof["ratios"], "wt_part": wt_part, "z_wt": wt_part["partition_sum"]}
    primary_label = ctx["primary_label"]

    cache: dict[tuple, dict] = ctx.setdefault("_null_qbin_cache", {}) if quantize_for_null else {}
    edits = []
    pairs = []
    for pos, base in enumerate(reference):
        for alt in "ACGT":
            if alt == base:
                continue
            eid = f"{pos + 1}:{base}>{alt}"
            variant = reference[:pos] + alt + reference[pos + 1 :]
            q_by_tf = {}
            for lbl, tbl in tf_tables.items():
                z_mut = tbl["z_wt"] + nc._delta_partition_single(reference, pos, alt, tbl["ratios"], tbl["width"])
                q_by_tf[lbl] = (z_mut + eff_floor) / (tbl["z_wt"] + eff_floor)
            ratio_map = {g: q_by_tf[nc._motif_spec_label(bindings[g]["protein"])] for g in bound_genes}
            reg_eval = _evaluate_regime_for_ratios(
                ctx, ratio_map, wt_profile, wt_base["wt_primary_critical"],
                wt_base["wt_preferred_fate"], signal_grid, cache,
                bound_genes=bound_genes, duration=duration, reporter_threshold=reporter_threshold,
                compute_transitions=compute_transitions,
                quantize_for_null=quantize_for_null,
            )
            pairs.append({"id": eid, "reference": reference, "sequence": variant})
            entry = {
                "id": eid,
                "pos": pos,
                "ref": base,
                "alt": alt,
                "edit_kind": "snv",
                "variant": variant,
                "affinity_ratio": q_by_tf[primary_label],
                "affinity_ratios": q_by_tf,
                "max_abs_reporter_delta": reg_eval["max_abs_reporter_delta"],
                "max_abs_target_delta": reg_eval["max_abs_target_delta"],
                "amplification_ratio": reg_eval["amplification_ratio"],
                "passive_silent": reg_eval["passive_silent"],
                "point_silent_flip": reg_eval["point_silent_flip"],
                "regime_flip_fraction": reg_eval["regime_flip_fraction"],
                "flipped_signal_count": reg_eval["flipped_signal_count"],
                "flipped_signals": reg_eval["flipped_signals"],
                "worst_signal": reg_eval["worst_signal"],
                "worst_wt_fate": reg_eval["worst_wt_fate"],
                "worst_mut_fate": reg_eval["worst_mut_fate"],
                "worst_d_rep": reg_eval["worst_d_rep"],
                "worst_d_target": reg_eval["worst_d_target"],
                "category": reg_eval["category"],
                "var_critical_signal": reg_eval["var_critical_signal"],
                "delta_s_crit": reg_eval["delta_s_crit"],
                "signal_locked": reg_eval["signal_locked"],
                "siraj_score": None,
                "kreimer_mpra_mean": None,
                "kreimer_mpra_72h": None,
                "sequence_silent": None,
                "doubly_silent": False,
            }
            edits.append(entry)

    siraj_scores, kreimer_scores = (None, None)
    if include_sequence_models:
        siraj_scores, kreimer_scores = _optional_sequence_scores(pairs)

    for e in edits:
        if siraj_scores is not None:
            e["siraj_score"] = round(siraj_scores[e["id"]], 6)
        if kreimer_scores is not None:
            e["kreimer_mpra_mean"] = round(kreimer_scores[e["id"]]["mean"], 6)
            e["kreimer_mpra_72h"] = round(kreimer_scores[e["id"]]["h72"], 6)
        if siraj_scores is not None and kreimer_scores is not None:
            seq_silent = abs(e["siraj_score"]) < sequence_threshold and abs(e["kreimer_mpra_mean"]) < sequence_threshold
            e["sequence_silent"] = seq_silent
            e["doubly_silent"] = bool(e["category"] == "cryptic_bifurcation" and seq_silent)

    cryptic = sorted(
        [e for e in edits if e["category"] == "cryptic_bifurcation"],
        key=lambda e: (
            not e["doubly_silent"],
            -e["regime_flip_fraction"],
            -e["max_abs_target_delta"],
            e["max_abs_reporter_delta"],
            e["id"],
        ),
    )
    doubly_silent_list = [e for e in cryptic if e["doubly_silent"]]
    counts = {
        cat: sum(e["category"] == cat for e in edits)
        for cat in ("cryptic_bifurcation", "concordant_switch", "canalized_buffered", "neutral")
    }
    destabilized = sum(
        1 for e in edits
        if e["signal_locked"] or (e["delta_s_crit"] is not None and abs(e["delta_s_crit"]) >= 0.01)
    )
    return {
        "circuit": ctx["circuit_name"],
        "genes": ctx["genes"],
        "protein": primary_label,
        "proteins": [lbl for lbl, _ in ctx["unique_specs"]],
        "target": ctx["target"],
        "epsilon": epsilon,
        "background_floor": eff_floor,
        "length": len(reference),
        "signal_grid": signal_grid,
        "wt_profile": wt_profile,
        "wt_transitions": wt_base["wt_transitions"],
        "wt_critical_signal": wt_base["wt_primary_critical"],
        "wt_preferred_fate": wt_base["wt_preferred_fate"],
        "wt_partition": tf_tables[primary_label]["wt_part"],
        "unique_regime_solves": len(cache),
        "category_counts": counts,
        "point_silent_switches": sum(1 for e in edits if e["point_silent_flip"]),
        "destabilized_count": destabilized,
        "doubly_silent_count": len(doubly_silent_list),
        "sequence_models_attached": siraj_scores is not None and kreimer_scores is not None,
        "cryptic_switches": cryptic,
        "edits": edits,
    }


def _rescue_candidate(
    reference: str,
    candidate: dict,
    ctx: dict,
    wt_base: dict,
    signal_grid: list[float],
    *,
    epsilon: float = nc.DEFAULT_EPSILON,
    background_floor: float | None = None,
    duration: float = 40.0,
) -> dict:
    """Compute both trans signal rescue and regime-evaluated cis epistatic rescue for a candidate."""
    eff_floor = nc._resolve_background_floor(background_floor, epsilon)
    bindings = ctx["bindings"]
    bound_genes = sorted(bindings)
    fail_sig = candidate["worst_signal"]
    desired_fate = candidate["worst_wt_fate"]

    # Variant ratios
    if candidate.get("variant"):
        q_by_tf = nc.multi_tf_affinity(
            reference, candidate["variant"], dict(ctx["unique_specs"]), epsilon, background_floor=eff_floor
        )["ratios"]
    else:
        q_by_tf = candidate["affinity_ratios"]
    var_ratio_map = {g: q_by_tf[nc._motif_spec_label(bindings[g]["protein"])] for g in bound_genes}

    trans_rescue = find_trans_signal_rescue(ctx, var_ratio_map, fail_sig, desired_fate, duration=duration)

    cis_rescue = None
    if candidate.get("edit_kind") == "snv" and "pos" in candidate and "alt" in candidate:
        raw_cis = nc.find_second_site_rescue(
            reference,
            candidate,
            protein=ctx["primary_protein"],
            epsilon=epsilon,
            background_floor=eff_floor,
            signal=fail_sig,
            target=ctx["target"],
            duration=duration,
            circuit=ctx["circuit_arg"],
            roles=ctx["roles"],
            mechanisms=ctx["mechanisms"],
        )
        if raw_cis is not None:
            dbl_seq = raw_cis["double_variant"]
            dbl_q_by_tf = nc.multi_tf_affinity(
                reference, dbl_seq, dict(ctx["unique_specs"]), epsilon, background_floor=eff_floor
            )["ratios"]
            dbl_ratio_map = {g: dbl_q_by_tf[nc._motif_spec_label(bindings[g]["protein"])] for g in bound_genes}
            dbl_eval = _evaluate_regime_for_ratios(
                ctx, dbl_ratio_map, wt_base["wt_profile"], wt_base["wt_primary_critical"],
                wt_base["wt_preferred_fate"], signal_grid, {},
                bound_genes=bound_genes, duration=duration,
            )
            before_frac = candidate["regime_flip_fraction"]
            after_frac = dbl_eval["regime_flip_fraction"]
            cis_rescue = {
                "second_site_id": raw_cis["id"],
                "second_site_pos": raw_cis["pos"],
                "double_affinity_ratio": round(raw_cis["affinity_ratio"], 6),
                "single_point_rescued": bool(raw_cis["fate_rescued"]),
                "regime_flip_fraction_before": before_frac,
                "regime_flip_fraction_after": after_frac,
                "regime_reduction": round(before_frac - after_frac, 6),
                "full_regime_rescued": after_frac == 0.0,
                "double_critical_signal": dbl_eval["var_critical_signal"],
                "double_delta_s_crit": dbl_eval["delta_s_crit"],
                "double_variant": dbl_seq,
                "double_profile": [
                    {"signal": r["signal"], "mut_target": r["mut_target"], "mut_fate": r["mut_fate"]}
                    for r in dbl_eval["per_signal"]
                ],
            }

    return {"trans_rescue": trans_rescue, "cis_rescue": cis_rescue}


def run_null_controls(
    reference: str,
    true_scan: dict,
    *,
    ctx: dict,
    signal_min: float,
    signal_max: float,
    signal_steps: int,
    n_shuffles: int = DEFAULT_SHUFFLES,
    epsilon: float = nc.DEFAULT_EPSILON,
    background_floor: float | None = None,
    duration: float = 40.0,
    reporter_threshold: float = DEFAULT_REPORTER_SILENT,
    seed: int = 17,
) -> dict:
    """Falsification harness: compare true enhancer against dinucleotide-shuffled decoys and scrambled PWM."""
    eff_floor = nc._resolve_background_floor(background_floor, epsilon)
    decoy_records = []
    for idx in range(n_shuffles):
        decoy_seq = dinucleotide_shuffle(reference, seed=seed + idx)
        if not verify_dinucleotide_counts(reference, decoy_seq):
            raise AssertionError("dinucleotide shuffle invariant violated")
        d_scan = scan_enhancer_regime(
            decoy_seq,
            signal_min=signal_min,
            signal_max=signal_max,
            signal_steps=signal_steps,
            epsilon=epsilon,
            background_floor=eff_floor,
            duration=duration,
            reporter_threshold=reporter_threshold,
            include_sequence_models=False,
            compute_transitions=False,
            quantize_for_null=True,
            ctx=ctx,
        )
        decoy_records.append({
            "decoy_index": idx,
            "seed": seed + idx,
            "cryptic_count": d_scan["category_counts"]["cryptic_bifurcation"],
            "point_silent_switches": d_scan["point_silent_switches"],
            "concordant_count": d_scan["category_counts"]["concordant_switch"],
            "regime_flip_edits": sum(1 for e in d_scan["edits"] if e["regime_flip_fraction"] > 0),
            "wt_partition_sum": round(d_scan["wt_partition"]["partition_sum"], 6),
            "wt_peak_weight": round(d_scan["wt_partition"]["peak_site"]["weight"], 6),
        })

    true_cryptic = true_scan["category_counts"]["cryptic_bifurcation"]
    true_peak = true_scan["wt_partition"]["peak_site"]["weight"]
    true_z = true_scan["wt_partition"]["partition_sum"]
    if decoy_records:
        c_vals = [r["cryptic_count"] for r in decoy_records]
        p_vals = [r["wt_peak_weight"] for r in decoy_records]
        z_vals = [r["wt_partition_sum"] for r in decoy_records]
        null_mean = statistics.mean(c_vals)
        null_std = statistics.pstdev(c_vals) if len(c_vals) > 1 else 0.0
        null_max = max(c_vals)
        emp_p = (1 + sum(c >= true_cryptic for c in c_vals)) / (1 + len(c_vals))
        fold = round(true_cryptic / null_mean, 4) if null_mean > 0 else (None if true_cryptic == 0 else float("inf"))
        peak_mean = statistics.mean(p_vals)
        peak_ratio = round(true_peak / peak_mean, 4) if peak_mean > 0 else None
        z_mean = statistics.mean(z_vals)
        z_ratio = round(true_z / z_mean, 4) if z_mean > 0 else None
        dino_summary = {
            "n_shuffles": n_shuffles,
            "dinucleotide_invariant_verified": True,
            "true_cryptic_count": true_cryptic,
            "null_cryptic_mean": round(null_mean, 4),
            "null_cryptic_std": round(null_std, 4),
            "null_cryptic_max": null_max,
            "empirical_p_ge": round(emp_p, 6),
            "fold_enrichment_over_null": None if fold == float("inf") else fold,
            "enrichment_infinite_zero_null": fold == float("inf"),
            "true_peak_weight": round(true_peak, 6),
            "null_peak_weight_mean": round(peak_mean, 6),
            "peak_weight_ratio": peak_ratio,
            "true_partition_sum": round(true_z, 6),
            "null_partition_sum_mean": round(z_mean, 6),
            "partition_sum_ratio": z_ratio,
            "true_exceeds_null_mean": true_cryptic > null_mean,
            "decoys": decoy_records,
        }
    else:
        dino_summary = {"n_shuffles": 0, "skipped": True}

    # Column-scrambled primary PWM control
    scrambled_spec = scramble_motif_pwm(ctx["primary_protein"], seed=seed + 1000)
    scram_ctx = nc._resolve_context(
        ctx["circuit_arg"], scrambled_spec, ctx["target"], 0.5 * (signal_min + signal_max),
        ctx["roles"], ctx["mechanisms"],
    )
    scram_ctx["_solve_cache"] = ctx.setdefault("_solve_cache", {})
    scram_ctx["_null_qbin_cache"] = ctx.setdefault("_null_qbin_cache", {})
    scram_scan = scan_enhancer_regime(
        reference,
        signal_min=signal_min,
        signal_max=signal_max,
        signal_steps=signal_steps,
        epsilon=epsilon,
        background_floor=eff_floor,
        duration=duration,
        reporter_threshold=reporter_threshold,
        include_sequence_models=False,
        compute_transitions=False,
        quantize_for_null=True,
        ctx=scram_ctx,
    )
    scram_cryptic = scram_scan["category_counts"]["cryptic_bifurcation"]
    scram_z = scram_scan["wt_partition"]["partition_sum"]
    motif_scramble_summary = {
        "scrambled_protein": scrambled_spec["id"],
        "background_floor": eff_floor,
        "true_cryptic_count": true_cryptic,
        "scrambled_cryptic_count": scram_cryptic,
        "true_regime_flip_edits": sum(1 for e in true_scan["edits"] if e["regime_flip_fraction"] > 0),
        "scrambled_regime_flip_edits": sum(1 for e in scram_scan["edits"] if e["regime_flip_fraction"] > 0),
        "true_peak_weight": round(true_peak, 6),
        "scrambled_peak_weight": round(scram_scan["wt_partition"]["peak_site"]["weight"], 6),
        "true_partition_sum": round(true_z, 6),
        "scrambled_partition_sum": round(scram_z, 6),
        "partition_sum_ratio": round(true_z / max(scram_z, 1e-12), 4),
        "true_exceeds_scrambled_pwm": true_cryptic > scram_cryptic,
        "biophysical_floor_note": (
            f"Partition ratios use q = (Z_mut + Z_ns) / (Z_ref + Z_ns) with non-specific floor Z_ns={eff_floor:g} "
            "to damp ratio blow-ups when baseline Z_ref << 1 in featureless/scrambled DNA."
        ),
    }
    return {
        "dinucleotide_shuffle": dino_summary,
        "motif_pwm_scramble": motif_scramble_summary,
    }


def verify_causal_controls(
    reference: str,
    top_candidate: dict | None,
    top_rescue: dict | None,
    *,
    ctx: dict,
    signal_grid: list[float],
    output: Path,
    epsilon: float = nc.DEFAULT_EPSILON,
    background_floor: float | None = None,
    duration: float = 40.0,
) -> dict:
    """Verify exact allele reversion, passive-reporter isolation, and sealed brainc compilation."""
    if top_candidate is None:
        return {"skipped": True, "reason": "no variant available for causal compilation check", "all_passed": True}

    eff_floor = nc._resolve_background_floor(background_floor, epsilon)
    bindings = ctx["bindings"]
    bound_genes = sorted(bindings)
    variant = top_candidate["variant"]
    motif_list = [spec for _, spec in ctx["unique_specs"]]

    compiled_dir = output / "compiled"
    compiled_dir.mkdir(parents=True, exist_ok=True)
    top_receipt = nc.compile_circuit_bundle(
        reference, variant, protein=ctx["primary_protein"], epsilon=epsilon,
        background_floor=eff_floor, output=compiled_dir / "top_candidate", motif_specs=motif_list, circuit=ctx["circuit_arg"],
    )
    top_compiled = nc.verify_circuit_bundle(compiled_dir / "top_candidate", top_receipt)

    cis_compiled = None
    cis_seq = (top_rescue or {}).get("cis_rescue", {}) or {}
    if cis_seq.get("double_variant"):
        cis_receipt = nc.compile_circuit_bundle(
            reference, cis_seq["double_variant"], protein=ctx["primary_protein"], epsilon=epsilon,
            background_floor=eff_floor, output=compiled_dir / "cis_rescue", motif_specs=motif_list, circuit=ctx["circuit_arg"],
        )
        cis_compiled = nc.verify_circuit_bundle(compiled_dir / "cis_rescue", cis_receipt)

    # Exact reversion and reporter isolation checks across the signal grid
    rev_q = nc.multi_tf_affinity(
        reference, reference, dict(ctx["unique_specs"]), epsilon, background_floor=eff_floor
    )["ratios"]
    rev_map = {g: rev_q[nc._motif_spec_label(bindings[g]["protein"])] for g in bound_genes}
    wt_map = dict.fromkeys(bound_genes, 1.0)

    var_q = nc.multi_tf_affinity(
        reference, variant, dict(ctx["unique_specs"]), epsilon, background_floor=eff_floor
    )["ratios"]
    var_map = {g: var_q[nc._motif_spec_label(bindings[g]["protein"])] for g in bound_genes}

    rev_max_diff = 0.0
    rep_max_diff = 0.0
    for s in signal_grid:
        wt_state = nc.simulate_circuit(
            [(duration, s)], ratios=wt_map, target=ctx["target"],
            circuit=ctx["circuit_arg"], roles=ctx["roles"], mechanisms=ctx["mechanisms"],
        )[-1][1:]
        rev_state = nc.simulate_circuit(
            [(duration, s)], ratios=rev_map, target=ctx["target"],
            circuit=ctx["circuit_arg"], roles=ctx["roles"], mechanisms=ctx["mechanisms"],
        )[-1][1:]
        rep_assay = nc.run_circuit_assay(
            [(duration, s)], ratio=var_map, mode="reporter", target=ctx["target"],
            circuit=ctx["circuit_arg"], roles=ctx["roles"], mechanisms=ctx["mechanisms"],
        )
        rep_state = rep_assay["trajectory"][-1][1:]
        rev_max_diff = max(rev_max_diff, max(abs(a - b) for a, b in zip(wt_state, rev_state)))
        rep_max_diff = max(rep_max_diff, max(abs(a - b) for a, b in zip(wt_state, rep_state)))

    checks = [
        {"name": "exact_allele_reversion_zero_diff", "passed": rev_max_diff == 0.0, "evidence": {"max_state_diff": rev_max_diff}},
        {"name": "passive_reporter_isolation", "passed": rep_max_diff == 0.0, "evidence": {"max_state_diff": rep_max_diff}},
        {"name": "top_candidate_brainc_bundle_verified", "passed": bool(top_compiled.get("bundle_sha256")), "evidence": {"contract_sha256": top_compiled["contract_sha256"]}},
    ]
    if cis_compiled is not None:
        checks.append({
            "name": "cis_rescue_brainc_bundle_verified",
            "passed": bool(cis_compiled.get("bundle_sha256")),
            "evidence": {"contract_sha256": cis_compiled["contract_sha256"]},
        })

    return {
        "top_compiled": top_compiled,
        "cis_compiled": cis_compiled,
        "checks": checks,
        "all_passed": all(c["passed"] for c in checks),
    }


def _candidate_row_for_export(rank: int, cand: dict, rescue_bundle: dict | None) -> dict:
    tr = (rescue_bundle or {}).get("trans_rescue") or {}
    cr = (rescue_bundle or {}).get("cis_rescue") or {}
    return {
        "rank": rank,
        "id": cand["id"],
        "edit_kind": cand.get("edit_kind", "snv"),
        "pos": cand.get("pos"),
        "ref": cand.get("ref"),
        "alt": cand.get("alt"),
        "category": cand["category"],
        "doubly_silent": cand.get("doubly_silent", False),
        "affinity_ratio": round(cand["affinity_ratio"], 6),
        "wt_critical_signal": cand.get("wt_critical_signal"),
        "var_critical_signal": cand.get("var_critical_signal"),
        "delta_s_crit": cand.get("delta_s_crit"),
        "signal_locked": cand.get("signal_locked", False),
        "max_abs_reporter_delta": cand["max_abs_reporter_delta"],
        "regime_flip_fraction": cand["regime_flip_fraction"],
        "max_abs_target_delta": cand["max_abs_target_delta"],
        "amplification_ratio": cand["amplification_ratio"],
        "siraj_score": cand.get("siraj_score"),
        "kreimer_mpra_mean": cand.get("kreimer_mpra_mean"),
        "kreimer_mpra_72h": cand.get("kreimer_mpra_72h"),
        "cis_rescue_id": cr.get("second_site_id"),
        "cis_rescue_double_ratio": cr.get("double_affinity_ratio"),
        "cis_rescue_full_regime": cr.get("full_regime_rescued"),
        "cis_rescue_regime_reduction": cr.get("regime_reduction"),
        "trans_rescue_rescued": tr.get("rescued"),
        "trans_rescue_signal": tr.get("rescue_signal"),
        "trans_rescue_signal_shift": tr.get("signal_shift"),
    }


def _write_candidates_csv(path: Path, rows: list[dict]) -> None:
    fields = [
        "rank", "id", "edit_kind", "pos", "ref", "alt", "category", "doubly_silent",
        "affinity_ratio", "wt_critical_signal", "var_critical_signal", "delta_s_crit", "signal_locked",
        "max_abs_reporter_delta", "regime_flip_fraction", "max_abs_target_delta", "amplification_ratio",
        "siraj_score", "kreimer_mpra_mean", "kreimer_mpra_72h",
        "cis_rescue_id", "cis_rescue_double_ratio", "cis_rescue_full_regime", "cis_rescue_regime_reduction",
        "trans_rescue_rescued", "trans_rescue_signal", "trans_rescue_signal_shift",
    ]
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k) for k in fields})


def _render_radar_html(report: dict) -> str:
    scan = report["scan"]
    edits = scan["edits"]
    target_gene = str(report["target"])
    protein = str(report["protein"])
    circuit_name = str(report["circuit"])
    null_dino = report["null_controls"]["dinucleotide_shuffle"]
    null_scram = report["null_controls"]["motif_pwm_scramble"]
    top_rows = report["wetlab_candidates"]
    featured = report.get("featured_candidate")

    w, h, pad = 540, 290, 48
    max_rep = max([0.05] + [e["max_abs_reporter_delta"] for e in edits])

    # SVG 1: Regime Tipping-Point Curve
    grid = scan["signal_grid"]
    wt_curve = [(w_r["signal"], w_r["target_level"]) for w_r in scan["wt_profile"]]

    def _to_pts(curve_vals: list[tuple[float, float]]) -> str:
        s_lo, s_hi = grid[0], grid[-1]
        span = max(s_hi - s_lo, 1e-6)
        return " ".join(
            f"{pad + (w - 2 * pad) * (s - s_lo) / span:.1f},{h - pad - (h - 2 * pad) * min(1.0, max(0.0, v)):.1f}"
            for s, v in curve_vals
        )

    curve_svg = f'<polyline fill="none" stroke="#111827" stroke-width="2.5" points="{_to_pts(wt_curve)}"/>'
    if report.get("featured_variant_profile"):
        mut_curve = [(r["signal"], r["mut_target"]) for r in report["featured_variant_profile"]]
        curve_svg += f'<polyline fill="none" stroke="#dc2626" stroke-width="2.3" stroke-dasharray="6 4" points="{_to_pts(mut_curve)}"/>'
    cis_prof = ((report.get("featured_rescues") or {}).get("cis_rescue") or {}).get("double_profile")
    if cis_prof:
        cis_curve = [(r["signal"], r["mut_target"]) for r in cis_prof]
        curve_svg += f'<polyline fill="none" stroke="#059669" stroke-width="2.2" stroke-dasharray="2 3" points="{_to_pts(cis_curve)}"/>'

    # SVG 2: Passive Reporter vs Endogenous Regime Shift scatter
    dots = []
    for e in edits:
        cx = pad + (w - 2 * pad) * min(1.0, e["max_abs_reporter_delta"] / max_rep)
        cy = h - pad - (h - 2 * pad) * min(1.0, e["max_abs_target_delta"])
        if e.get("doubly_silent"):
            col = "#7c3aed"
            rad = 4.5
        elif e["category"] == "cryptic_bifurcation":
            col = "#dc2626"
            rad = 3.8
        elif e["category"] == "concordant_switch":
            col = "#d97706"
            rad = 3.0
        elif e["category"] == "canalized_buffered":
            col = "#2563eb"
            rad = 3.0
        else:
            col = "#9ca3af"
            rad = 2.4
        if featured and e["id"] == featured["id"]:
            rad = 5.8
        title = html.escape(
            f"{e['id']} | q={e['affinity_ratio']:.4f} | max|dRep|={e['max_abs_reporter_delta']:.4f} "
            f"| max|dTgt|={e['max_abs_target_delta']:.3f} | flipFrac={e['regime_flip_fraction']:.2f} ({e['category']})"
        )
        dots.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{rad}" fill="{col}" opacity="0.82"><title>{title}</title></circle>')

    # Plain-English interpretation summary
    n_cryptic = scan["category_counts"]["cryptic_bifurcation"]
    n_doubly = scan["doubly_silent_count"]
    null_note = (
        f"Tested {null_dino['n_shuffles']} matched dinucleotide-shuffled controls: true enhancer has "
        f"{null_dino['true_cryptic_count']} cryptic variants vs null mean {null_dino['null_cryptic_mean']} "
        f"(range up to {null_dino['null_cryptic_max']}, empirical p={null_dino['empirical_p_ge']}); "
        f"column-scrambled motif control yielded {null_scram['scrambled_cryptic_count']} cryptic variants."
        if not null_dino.get("skipped")
        else "Dinucleotide shuffle null controls were skipped (--shuffles 0)."
    )

    row_chunks = []
    for r in top_rows:
        if r.get("signal_locked"):
            dscrit_str = "LOCKED"
        elif r.get("delta_s_crit") is not None:
            dscrit_str = f"{r['delta_s_crit']:+.4f}"
        else:
            dscrit_str = "n/a"
        siraj_str = "n/a" if r["siraj_score"] is None else f"{r['siraj_score']:+.3f}"
        kreimer_str = "n/a" if r["kreimer_mpra_mean"] is None else f"{r['kreimer_mpra_mean']:+.3f}"
        cis_tag = ""
        if r.get("cis_rescue_full_regime"):
            cis_tag = " (full)"
        elif r.get("cis_rescue_regime_reduction"):
            cis_tag = f" (-{r['cis_rescue_regime_reduction']:.2f})"
        trans_str = "none" if not r.get("trans_rescue_rescued") else f"dS={r['trans_rescue_signal_shift']:+.3f}"
        silent_tag = " [2x-silent]" if r.get("doubly_silent") else ""
        row_chunks.append(
            f"<tr><td>{r['rank']}</td><td><code>{html.escape(str(r['id']))}</code></td>"
            f"<td>{html.escape(str(r['category']))}{silent_tag}</td>"
            f"<td>{r['affinity_ratio']:.4f}</td>"
            f"<td>{dscrit_str}</td>"
            f"<td>{r['max_abs_reporter_delta']:.4f}</td>"
            f"<td>{r['regime_flip_fraction']:.2f}</td>"
            f"<td>{r['max_abs_target_delta']:.3f}</td>"
            f"<td>{siraj_str} / {kreimer_str}</td>"
            f"<td>{html.escape(str(r['cis_rescue_id'] or 'none'))}{cis_tag}</td>"
            f"<td>{trans_str}</td></tr>"
        )
    cand_rows_html = "".join(row_chunks)
    if not cand_rows_html:
        cand_rows_html = '<tr><td colspan="11">No cryptic or fate-flipping edits found in this regime.</td></tr>'

    checks_html = "".join(
        f"<li><strong>{'PASS' if c['passed'] else 'FAIL'}</strong> &mdash; {html.escape(c['name'])}: "
        f"<code>{html.escape(json.dumps(c['evidence']))}</code></li>"
        for c in report["causal_controls"].get("checks", [])
    )

    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>Epigenesis · V2F Tipping-Point &amp; Rescue Radar</title>
<style>
body{{font-family:system-ui,-apple-system,sans-serif;max-width:1180px;margin:24px auto;padding:0 20px;color:#111;line-height:1.45}}
.grid{{display:grid;grid-template-columns:1fr 1fr;gap:20px;margin:16px 0}}
.callout{{background:#f8fafc;border-left:4px solid #0f172a;padding:12px 16px;margin:12px 0}}
svg{{width:100%;border:1px solid #e5e7eb;background:#fff}}
table{{width:100%;border-collapse:collapse;font-size:13px;margin:12px 0}}
th,td{{padding:7px 8px;border-bottom:1px solid #e5e7eb;text-align:left}}
th{{background:#f8fafc}}
code{{background:#f1f5f9;padding:1px 4px;border-radius:3px}}
</style></head><body>
<h1>Variant-to-Function (V2F) Tipping-Point &amp; Rescue Radar</h1>
<div class="callout">
<strong>Plain-English Summary:</strong> Standard lab reporter assays test an enhancer without coupled downstream feedback, often at a single condition.
Here we sweep the entire developmental signaling window <code>[{grid[0]:.2f}, {grid[-1]:.2f}]</code> for <strong>{html.escape(report['label'])}</strong>
({scan['length']} bp, circuit <code>{html.escape(circuit_name)}</code>, TF <code>{html.escape(protein)}</code>, target <code>{html.escape(target_gene)}</code>).
We detected <strong>{n_cryptic}</strong> passive-reporter-silent edits that still flip endogenous cell fate
{f'(including <strong>{n_doubly}</strong> that also score near zero in Siraj/Kreimer sequence predictors)' if scan['sequence_models_attached'] else ''}.
<br><strong>Null-Control Reality Check:</strong> {html.escape(null_note)}
<br><em>Conditionality:</em> Mechanistic assumptions only; biological embryo/clinical effects are NOT established (<code>biological_acceptance = false</code>).
</div>
<div class="grid">
<div><h3>Morphogen Tipping-Point Curve (Solid=WT, Red=Top Variant, Green=Cis Rescue)</h3>
<svg viewBox="0 0 {w} {h}"><line x1="{pad}" y1="{h-pad}" x2="{w-pad}" y2="{h-pad}" stroke="#111"/>
<line x1="{pad}" y1="{pad}" x2="{pad}" y2="{h-pad}" stroke="#111"/>
<text x="{w/2}" y="{h-10}" text-anchor="middle" font-size="12">Developmental Signal S in [{grid[0]:.2f}, {grid[-1]:.2f}]</text>
<text x="16" y="{h/2}" transform="rotate(-90 16 {h/2})" text-anchor="middle" font-size="12">{html.escape(target_gene)} Steady Level</text>
{curve_svg}</svg></div>
<div><h3>Passive Reporter vs. Endogenous Regime Impact</h3>
<svg viewBox="0 0 {w} {h}"><line x1="{pad}" y1="{h-pad}" x2="{w-pad}" y2="{h-pad}" stroke="#111"/>
<line x1="{pad}" y1="{pad}" x2="{pad}" y2="{h-pad}" stroke="#111"/>
<text x="{w/2}" y="{h-10}" text-anchor="middle" font-size="12">Max |Passive Reporter Delta| across Regime</text>
<text x="16" y="{h/2}" transform="rotate(-90 16 {h/2})" text-anchor="middle" font-size="12">Max |Endogenous {html.escape(target_gene)} Delta|</text>
{''.join(dots)}</svg></div>
</div>
<h3>Prioritized Experimental Candidates &amp; Rescues</h3>
<table><thead><tr>
<th>#</th><th>Edit</th><th>Class</th><th>q</th><th>&Delta;S<sub>crit</sub></th>
<th>Max |&Delta;Rep|</th><th>Flip Frac</th><th>Max |&Delta;Tgt|</th>
<th>Siraj / Kreimer</th><th>Cis Rescue</th><th>Trans Rescue</th>
</tr></thead><tbody>{cand_rows_html}</tbody></table>
<h3>Causal &amp; Sealed <code>brainc</code> Verification Controls</h3>
<ul>{checks_html}</ul>
</body></html>"""


def run_radar(
    sequence: str,
    output: Path,
    *,
    label: str = "custom-enhancer",
    circuit: str | dict | Path | None = None,
    protein: str | dict | Path | None = None,
    target: str | list[str] | tuple[str, ...] | None = None,
    signal_min: float = DEFAULT_SIGNAL_MIN,
    signal_max: float = DEFAULT_SIGNAL_MAX,
    signal_steps: int = DEFAULT_SIGNAL_STEPS,
    n_shuffles: int = DEFAULT_SHUFFLES,
    top_k: int = DEFAULT_TOP_K,
    epsilon: float = nc.DEFAULT_EPSILON,
    background_floor: float | None = None,
    duration: float = 40.0,
    include_sequence_models: bool = True,
    seed: int = 17,
) -> dict:
    """Run full V2F tipping-point scan, rescues, null shuffles, brainc verification, and exports."""
    output.mkdir(parents=True, exist_ok=True)
    eff_floor = nc._resolve_background_floor(background_floor, epsilon)
    ctx = nc._resolve_context(circuit, protein, target, 0.5 * (signal_min + signal_max), None, None)
    scan = scan_enhancer_regime(
        sequence,
        circuit=circuit,
        protein=protein,
        target=target,
        signal_min=signal_min,
        signal_max=signal_max,
        signal_steps=signal_steps,
        epsilon=epsilon,
        background_floor=eff_floor,
        duration=duration,
        include_sequence_models=include_sequence_models,
        ctx=ctx,
    )
    signal_grid = scan["signal_grid"]
    wt_base = {
        "wt_profile": scan["wt_profile"],
        "wt_transitions": scan["wt_transitions"],
        "wt_primary_critical": scan["wt_critical_signal"],
        "wt_preferred_fate": scan["wt_preferred_fate"],
    }

    # Pick top candidates for full dual rescue evaluation
    pool = scan["cryptic_switches"] if scan["cryptic_switches"] else sorted(
        scan["edits"], key=lambda e: (-e["regime_flip_fraction"], -e["max_abs_target_delta"], e["id"])
    )
    selected = pool[:top_k]
    rescues_by_id: dict[str, dict] = {}
    export_rows = []
    for rank, cand in enumerate(selected, start=1):
        cand_with_wt = dict(cand, wt_critical_signal=scan["wt_critical_signal"])
        rescues = _rescue_candidate(
            sequence, cand_with_wt, ctx, wt_base, signal_grid,
            epsilon=epsilon, background_floor=eff_floor, duration=duration,
        )
        rescues_by_id[cand["id"]] = rescues
        export_rows.append(_candidate_row_for_export(rank, cand_with_wt, rescues))

    featured = selected[0] if selected else None
    featured_rescues = rescues_by_id.get(featured["id"]) if featured else None
    featured_profile = None
    if featured is not None:
        bindings = ctx["bindings"]
        var_q = featured["affinity_ratios"]
        var_map = {g: var_q[nc._motif_spec_label(bindings[g]["protein"])] for g in sorted(bindings)}
        featured_profile = _evaluate_regime_for_ratios(
            ctx, var_map, scan["wt_profile"], scan["wt_critical_signal"],
            scan["wt_preferred_fate"], signal_grid, {}, bound_genes=sorted(bindings), duration=duration,
        )["per_signal"]

    nulls = run_null_controls(
        sequence,
        scan,
        ctx=ctx,
        signal_min=signal_min,
        signal_max=signal_max,
        signal_steps=signal_steps,
        n_shuffles=n_shuffles,
        epsilon=epsilon,
        background_floor=eff_floor,
        duration=duration,
        seed=seed,
    )
    causal = verify_causal_controls(
        sequence,
        featured,
        featured_rescues,
        ctx=ctx,
        signal_grid=signal_grid,
        output=output,
        epsilon=epsilon,
        background_floor=eff_floor,
        duration=duration,
    )

    # Strip full variant strings from bulk report for concise clean JSON
    compact_scan = {
        **{k: v for k, v in scan.items() if k not in ("edits", "cryptic_switches")},
        "cryptic_switches": [{k: v for k, v in e.items() if k != "variant"} for e in scan["cryptic_switches"]],
        "edits": [{k: v for k, v in e.items() if k != "variant"} for e in scan["edits"]],
    }
    compact_rescues = None
    if featured_rescues is not None:
        compact_rescues = {
            "trans_rescue": featured_rescues["trans_rescue"],
            "cis_rescue": (
                {k: v for k, v in featured_rescues["cis_rescue"].items() if k != "double_variant"}
                if featured_rescues["cis_rescue"]
                else None
            ),
        }

    report = {
        "format": "epigenesis.neural-v2f-radar/v1",
        "biological_acceptance": False,
        "scope": (
            "Computational Variant-to-Function tipping-point, rescue, and null-control radar. "
            "Parameters and rates are uncalibrated model assumptions; biological phenotype is NOT established."
        ),
        "label": label,
        "circuit": scan["circuit"],
        "genes": scan["genes"],
        "protein": scan["protein"],
        "proteins": scan["proteins"],
        "target": scan["target"],
        "epsilon": epsilon,
        "background_floor": eff_floor,
        "wt_critical_signal": scan["wt_critical_signal"],
        "featured_candidate": {k: v for k, v in featured.items() if k != "variant"} if featured else None,
        "featured_variant_profile": featured_profile,
        "featured_rescues": compact_rescues,
        "wetlab_candidates": export_rows,
        "null_controls": nulls,
        "causal_controls": causal,
        "scan": compact_scan,
    }
    (output / "radar_report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (output / "wetlab_candidates.json").write_text(json.dumps(export_rows, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    _write_candidates_csv(output / "wetlab_candidates.csv", export_rows)
    (output / "index.html").write_text(_render_radar_html(report), encoding="utf-8")
    return report


def _load_variants_file(path: Path) -> list[dict]:
    """Load targeted variant table from JSON, CSV, or TSV."""
    raw = read_regular_file(path, maximum_bytes=4 * 1024 * 1024).decode("utf-8")
    suf = path.suffix.lower()
    if suf == ".json":
        data = json.loads(raw)
        if not isinstance(data, list) or not data:
            raise ValueError("variants JSON must be a non-empty list")
        return data
    delim = "\t" if suf == ".tsv" else ","
    reader = csv.DictReader(raw.splitlines(), delimiter=delim)
    rows = [dict(r) for r in reader]
    if not rows:
        raise ValueError("variants CSV/TSV contains no rows")
    return rows


def analyze_variant_table(
    variants_path: Path,
    output: Path,
    *,
    circuit: str | dict | Path | None = None,
    protein: str | dict | Path | None = None,
    target: str | list[str] | tuple[str, ...] | None = None,
    signal_min: float = DEFAULT_SIGNAL_MIN,
    signal_max: float = DEFAULT_SIGNAL_MAX,
    signal_steps: int = DEFAULT_SIGNAL_STEPS,
    epsilon: float = nc.DEFAULT_EPSILON,
    background_floor: float | None = None,
    duration: float = 40.0,
    include_sequence_models: bool = True,
) -> dict:
    """Evaluate arbitrary SNV/MNV/insertion/deletion variants from a JSON/CSV/TSV table across the full signal regime."""
    output.mkdir(parents=True, exist_ok=True)
    eff_floor = nc._resolve_background_floor(background_floor, epsilon)
    rows = _load_variants_file(variants_path)
    ctx = nc._resolve_context(circuit, protein, target, 0.5 * (signal_min + signal_max), None, None)
    signal_grid = _signal_grid(signal_min, signal_max, signal_steps)
    bindings = ctx["bindings"]
    bound_genes = sorted(bindings)
    wt_base = _build_wt_baseline(ctx, signal_grid, bound_genes, duration=duration)
    cache: dict[tuple, dict] = {}

    parsed = []
    seq_pairs = []
    for idx, row in enumerate(rows):
        eid = str(row.get("id") or f"var-{idx + 1}")
        ref = str(row.get("reference") or "").strip().upper()
        alt_seq = row.get("sequence") or row.get("variant")
        if not ref:
            raise ValueError(f"variant row {eid} requires 'reference' DNA")
        if alt_seq:
            var = str(alt_seq).strip().upper()
        elif row.get("pos") is not None and row.get("alt"):
            p0 = int(row["pos"]) - 1
            alt_b = str(row["alt"]).strip().upper()
            if not 0 <= p0 < len(ref) or len(alt_b) != 1 or alt_b not in "ACGT":
                raise ValueError(f"variant row {eid} has invalid 1-based SNV pos/alt")
            var = ref[:p0] + alt_b + ref[p0 + 1 :]
        else:
            raise ValueError(f"variant row {eid} requires 'sequence'/'variant' or 1-based 'pos'+'alt'")
        if set(ref) - set("ACGT") or set(var) - set("ACGT"):
            raise ValueError(f"variant row {eid} contains non-ACGT bases")

        q_by_tf = nc.multi_tf_affinity(
            ref, var, dict(ctx["unique_specs"]), epsilon, background_floor=eff_floor
        )["ratios"]
        ratio_map = {g: q_by_tf[nc._motif_spec_label(bindings[g]["protein"])] for g in bound_genes}
        reg_eval = _evaluate_regime_for_ratios(
            ctx, ratio_map, wt_base["wt_profile"], wt_base["wt_primary_critical"],
            wt_base["wt_preferred_fate"], signal_grid, cache, bound_genes=bound_genes, duration=duration,
        )
        kind = nc._edit_kind(ref, var)
        pos_val = None
        ref_b = None
        alt_b = None
        if kind == "snv":
            pos_val = next(i for i, (a, b) in enumerate(zip(ref, var)) if a != b)
            ref_b = ref[pos_val]
            alt_b = var[pos_val]
        seq_pairs.append({"id": eid, "reference": ref, "sequence": var})
        parsed.append({
            "id": eid,
            "edit_kind": kind,
            "pos": pos_val,
            "ref": ref_b,
            "alt": alt_b,
            "reference": ref,
            "variant": var,
            "affinity_ratio": q_by_tf[ctx["primary_label"]],
            "affinity_ratios": q_by_tf,
            "wt_critical_signal": wt_base["wt_primary_critical"],
            **{k: reg_eval[k] for k in (
                "max_abs_reporter_delta", "max_abs_target_delta", "amplification_ratio",
                "passive_silent", "regime_flip_fraction", "flipped_signal_count", "flipped_signals",
                "worst_signal", "worst_wt_fate", "worst_mut_fate", "worst_d_rep", "worst_d_target",
                "category", "var_critical_signal", "delta_s_crit", "signal_locked",
            )},
            "siraj_score": None,
            "kreimer_mpra_mean": None,
            "kreimer_mpra_72h": None,
            "sequence_silent": None,
            "doubly_silent": False,
        })

    siraj_scores, kreimer_scores = (None, None)
    if include_sequence_models:
        siraj_scores, kreimer_scores = _optional_sequence_scores(seq_pairs)

    export_rows = []
    for idx, item in enumerate(parsed, start=1):
        if siraj_scores is not None:
            item["siraj_score"] = round(siraj_scores[item["id"]], 6)
        if kreimer_scores is not None:
            item["kreimer_mpra_mean"] = round(kreimer_scores[item["id"]]["mean"], 6)
            item["kreimer_mpra_72h"] = round(kreimer_scores[item["id"]]["h72"], 6)
        if siraj_scores is not None and kreimer_scores is not None:
            seq_silent = abs(item["siraj_score"]) < DEFAULT_SEQUENCE_SILENT and abs(item["kreimer_mpra_mean"]) < DEFAULT_SEQUENCE_SILENT
            item["sequence_silent"] = seq_silent
            item["doubly_silent"] = bool(item["category"] == "cryptic_bifurcation" and seq_silent)
        rescues = _rescue_candidate(
            item["reference"], item, ctx, wt_base, signal_grid,
            epsilon=epsilon, background_floor=eff_floor, duration=duration,
        )
        item["rescues"] = {
            "trans_rescue": rescues["trans_rescue"],
            "cis_rescue": (
                {k: v for k, v in rescues["cis_rescue"].items() if k != "double_variant"}
                if rescues["cis_rescue"] else None
            ),
        }
        export_rows.append(_candidate_row_for_export(idx, item, rescues))

    top_cand = max(parsed, key=lambda e: (e["category"] == "cryptic_bifurcation", e["regime_flip_fraction"], e["max_abs_target_delta"]))
    top_rescue = _rescue_candidate(
        top_cand["reference"], top_cand, ctx, wt_base, signal_grid,
        epsilon=epsilon, background_floor=eff_floor, duration=duration,
    )
    causal = verify_causal_controls(
        top_cand["reference"], top_cand, top_rescue, ctx=ctx,
        signal_grid=signal_grid, output=output, epsilon=epsilon,
        background_floor=eff_floor, duration=duration,
    )

    report = {
        "format": "epigenesis.neural-v2f-variants/v1",
        "biological_acceptance": False,
        "scope": "Targeted arbitrary-length variant table analysis (SNVs, MNVs, indels); biological phenotype is NOT established.",
        "circuit": ctx["circuit_name"],
        "genes": ctx["genes"],
        "protein": ctx["primary_label"],
        "target": ctx["target"],
        "epsilon": epsilon,
        "background_floor": eff_floor,
        "signal_grid": signal_grid,
        "wt_critical_signal": wt_base["wt_primary_critical"],
        "variants_tested": len(parsed),
        "cryptic_count": sum(p["category"] == "cryptic_bifurcation" for p in parsed),
        "doubly_silent_count": sum(bool(p["doubly_silent"]) for p in parsed),
        "variants": [{k: v for k, v in p.items() if k not in ("reference", "variant")} for p in parsed],
        "wetlab_candidates": export_rows,
        "causal_controls": causal,
    }
    (output / "variants_report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (output / "wetlab_candidates.json").write_text(json.dumps(export_rows, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    _write_candidates_csv(output / "wetlab_candidates.csv", export_rows)
    return report


def benchmark_radar(
    output: Path,
    *,
    loci: tuple[str, ...] | list[str] = nc.CURATED_LOCI,
    circuit: str | dict | Path | None = None,
    protein: str | dict | Path | None = None,
    target: str | list[str] | tuple[str, ...] | None = None,
    signal_min: float = DEFAULT_SIGNAL_MIN,
    signal_max: float = DEFAULT_SIGNAL_MAX,
    signal_steps: int = DEFAULT_SIGNAL_STEPS,
    n_shuffles: int = 5,
    epsilon: float = nc.DEFAULT_EPSILON,
    background_floor: float | None = None,
    duration: float = 40.0,
    include_sequence_models: bool = True,
    seed: int = 17,
) -> dict:
    """Multi-locus GSE188264 benchmark + real measured perturbation discordance + dinucleotide null audit."""
    sys.path.insert(0, str(REPO / "examples/neural-regulation"))
    import regulation as reg

    rows, _ = reg.load_data()
    by_locus: dict[str, list] = {}
    for r in rows:
        by_locus.setdefault(r.locus, []).append(r)

    output.mkdir(parents=True, exist_ok=True)
    eff_floor = nc._resolve_background_floor(background_floor, epsilon)
    ctx = nc._resolve_context(circuit, protein, target, 0.5 * (signal_min + signal_max), None, None)
    signal_grid = _signal_grid(signal_min, signal_max, signal_steps)
    bindings = ctx["bindings"]
    bound_genes = sorted(bindings)

    locus_reports = []
    measured_records = []
    all_wetlab_rows = []
    for l_idx, locus in enumerate(loci):
        locus_rows = by_locus.get(locus)
        if not locus_rows:
            raise ValueError(f"unknown GSE188264 locus: {locus}")
        ref = locus_rows[0].reference
        scan = scan_enhancer_regime(
            ref, signal_min=signal_min, signal_max=signal_max, signal_steps=signal_steps,
            epsilon=epsilon, background_floor=eff_floor, duration=duration,
            include_sequence_models=include_sequence_models, ctx=ctx,
        )
        wt_base = {
            "wt_profile": scan["wt_profile"],
            "wt_transitions": scan["wt_transitions"],
            "wt_primary_critical": scan["wt_critical_signal"],
            "wt_preferred_fate": scan["wt_preferred_fate"],
        }
        nulls = run_null_controls(
            ref, scan, ctx=ctx, signal_min=signal_min, signal_max=signal_max, signal_steps=signal_steps,
            n_shuffles=n_shuffles, epsilon=epsilon, background_floor=eff_floor,
            duration=duration, seed=seed + 100 * l_idx,
        )

        # Evaluate real measured sequence perturbations at this locus
        seen_seq: set[str] = set()
        loc_measured = []
        meas_pairs = []
        for r in locus_rows:
            if r.kind == "wild-type" or r.sequence == ref or r.sequence in seen_seq:
                continue
            seen_seq.add(r.sequence)
            q_by_tf = nc.multi_tf_affinity(
                ref, r.sequence, dict(ctx["unique_specs"]), epsilon, background_floor=eff_floor
            )["ratios"]
            ratio_map = {g: q_by_tf[nc._motif_spec_label(bindings[g]["protein"])] for g in bound_genes}
            reg_eval = _evaluate_regime_for_ratios(
                ctx, ratio_map, wt_base["wt_profile"], wt_base["wt_primary_critical"],
                wt_base["wt_preferred_fate"], signal_grid, {}, bound_genes=bound_genes, duration=duration,
            )
            alphas = list(r.alpha)
            meas_pairs.append({"id": r.id, "reference": ref, "sequence": r.sequence})
            loc_measured.append({
                "locus": locus,
                "id": r.id,
                "split": r.split,
                "kind": r.kind,
                "edit_kind": nc._edit_kind(ref, r.sequence),
                "affinity_ratio": round(q_by_tf[ctx["primary_label"]], 6),
                "measured_alpha_mean": round(sum(alphas) / len(alphas), 6),
                "measured_alpha_72h": round(alphas[-1], 6),
                "measured_alpha_max_abs": round(max(abs(v) for v in alphas), 6),
                "max_abs_reporter_delta": reg_eval["max_abs_reporter_delta"],
                "regime_flip_fraction": reg_eval["regime_flip_fraction"],
                "max_abs_target_delta": reg_eval["max_abs_target_delta"],
                "delta_s_crit": reg_eval["delta_s_crit"],
                "signal_locked": reg_eval["signal_locked"],
                "category": reg_eval["category"],
                "siraj_score": None,
                "kreimer_mpra_mean": None,
            })

        if include_sequence_models and meas_pairs:
            s_sc, k_sc = _optional_sequence_scores(meas_pairs)
            for m in loc_measured:
                if s_sc is not None:
                    m["siraj_score"] = round(s_sc[m["id"]], 6)
                if k_sc is not None:
                    m["kreimer_mpra_mean"] = round(k_sc[m["id"]]["mean"], 6)
        measured_records.extend(loc_measured)

        top = scan["cryptic_switches"][0] if scan["cryptic_switches"] else None
        top_rescue = None
        if top is not None:
            top_with_wt = dict(top, wt_critical_signal=scan["wt_critical_signal"])
            top_rescue = _rescue_candidate(
                ref, top_with_wt, ctx, wt_base, signal_grid,
                epsilon=epsilon, background_floor=eff_floor, duration=duration,
            )
            all_wetlab_rows.append(dict(_candidate_row_for_export(len(all_wetlab_rows) + 1, top_with_wt, top_rescue), locus=locus))

        locus_reports.append({
            "locus": locus,
            "length": len(ref),
            "wt_critical_signal": scan["wt_critical_signal"],
            "wt_peak_weight": round(scan["wt_partition"]["peak_site"]["weight"], 6),
            "category_counts": scan["category_counts"],
            "doubly_silent_count": scan["doubly_silent_count"],
            "destabilized_count": scan["destabilized_count"],
            "null_dinucleotide_mean": nulls["dinucleotide_shuffle"].get("null_cryptic_mean"),
            "null_dinucleotide_p": nulls["dinucleotide_shuffle"].get("empirical_p_ge"),
            "scrambled_pwm_cryptic": nulls["motif_pwm_scramble"]["scrambled_cryptic_count"],
            "top_cryptic": {k: v for k, v in top.items() if k != "variant"} if top else None,
            "top_rescue": (
                {
                    "cis_id": (top_rescue["cis_rescue"] or {}).get("second_site_id"),
                    "cis_full_regime": (top_rescue["cis_rescue"] or {}).get("full_regime_rescued"),
                    "trans_shift": top_rescue["trans_rescue"].get("signal_shift"),
                }
                if top_rescue
                else None
            ),
            "measured_variants_tested": len(loc_measured),
            "measured_cryptic_switches": sum(m["category"] == "cryptic_bifurcation" for m in loc_measured),
            "measured_any_switches": sum(m["regime_flip_fraction"] > 0 for m in loc_measured),
        })

    # Verify one sealed brainc bundle from the first locus with a cryptic candidate (or first locus)
    first_ref = by_locus[loci[0]][0].reference
    first_scan = scan_enhancer_regime(
        first_ref, signal_min=signal_min, signal_max=signal_max, signal_steps=signal_steps,
        epsilon=epsilon, background_floor=eff_floor, duration=duration, include_sequence_models=False, ctx=ctx,
    )
    rep_cand = first_scan["cryptic_switches"][0] if first_scan["cryptic_switches"] else first_scan["edits"][0]
    causal = verify_causal_controls(
        first_ref, rep_cand, None, ctx=ctx, signal_grid=signal_grid, output=output,
        epsilon=epsilon, background_floor=eff_floor, duration=duration,
    )

    summary = {
        "format": "epigenesis.neural-v2f-benchmark/v1",
        "biological_acceptance": False,
        "scope": "Multi-locus GSE188264 tipping-point benchmark, measured perturbation discordance, and dinucleotide null controls; biological acceptance is FALSE.",
        "circuit": ctx["circuit_name"],
        "protein": ctx["primary_label"],
        "target": ctx["target"],
        "epsilon": epsilon,
        "background_floor": eff_floor,
        "signal_grid": signal_grid,
        "loci_tested": len(locus_reports),
        "total_cryptic_bifurcations": sum(r["category_counts"]["cryptic_bifurcation"] for r in locus_reports),
        "total_doubly_silent": sum(r["doubly_silent_count"] for r in locus_reports),
        "total_measured_variants": len(measured_records),
        "total_measured_switches": sum(m["regime_flip_fraction"] > 0 for m in measured_records),
        "total_measured_cryptic": sum(m["category"] == "cryptic_bifurcation" for m in measured_records),
        "loci": locus_reports,
        "measured_variants": measured_records,
        "causal_controls": causal,
    }
    (output / "radar_benchmark.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    _write_candidates_csv(output / "wetlab_candidates.csv", all_wetlab_rows)

    m_fields = [
        "locus", "id", "split", "kind", "edit_kind", "affinity_ratio",
        "measured_alpha_mean", "measured_alpha_72h", "measured_alpha_max_abs",
        "max_abs_reporter_delta", "regime_flip_fraction", "max_abs_target_delta",
        "delta_s_crit", "signal_locked", "category", "siraj_score", "kreimer_mpra_mean",
    ]
    with (output / "measured_discordance.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=m_fields)
        w.writeheader()
        for m in measured_records:
            w.writerow({k: m.get(k) for k in m_fields})
    return summary


def _add_regime_args(sp: argparse.ArgumentParser) -> None:
    sp.add_argument("--circuit", help="circuit preset ('neural-tube', 'mutual-toggle', 'tri-lineage-switch') or JSON spec")
    sp.add_argument("--protein", help="TF motif: GLI2, GLI3, JASPAR ID (e.g. MA0069.1), or custom PWM JSON")
    sp.add_argument("--target", help="target gene in chosen circuit (e.g. Nkx2.2, Olig2, GeneA, Excitatory, Both, All)")
    sp.add_argument("--signal-min", type=float, default=DEFAULT_SIGNAL_MIN)
    sp.add_argument("--signal-max", type=float, default=DEFAULT_SIGNAL_MAX)
    sp.add_argument("--steps", type=int, default=DEFAULT_SIGNAL_STEPS)
    sp.add_argument("--background-floor", type=float, default=None, help=f"non-specific partition floor Z_ns (default {nc.DEFAULT_BACKGROUND_FLOOR:g} when epsilon>0)")
    sp.add_argument("--no-sequence-models", action="store_true", help="skip optional Siraj/Kreimer sequence predictor scoring")
    sp.add_argument("--output", type=Path, required=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    d = sub.add_parser("demo", help="run V2F tipping-point & rescue radar on a curated GSE188264 human enhancer")
    d.add_argument("--locus", default=nc.CURATED_LOCI[0])
    d.add_argument("--shuffles", type=int, default=DEFAULT_SHUFFLES)
    d.add_argument("--top-k", type=int, default=DEFAULT_TOP_K)
    _add_regime_args(d)

    s = sub.add_parser("scan", help="run V2F tipping-point & rescue radar on any FASTA enhancer")
    s.add_argument("--fasta", type=Path, required=True)
    s.add_argument("--shuffles", type=int, default=DEFAULT_SHUFFLES)
    s.add_argument("--top-k", type=int, default=DEFAULT_TOP_K)
    _add_regime_args(s)

    v = sub.add_parser("variants", help="evaluate arbitrary SNV/MNV/indel variant table across the signaling regime")
    v.add_argument("--variants", type=Path, required=True, help="JSON, CSV, or TSV table of variants")
    _add_regime_args(v)

    b = sub.add_parser("benchmark", help="benchmark across curated GSE188264 loci + real measured CRE perturbations")
    b.add_argument("--shuffles", type=int, default=5)
    _add_regime_args(b)

    args = parser.parse_args()
    try:
        seq_on = not args.no_sequence_models
        if args.command == "demo":
            sys.path.insert(0, str(REPO / "examples/neural-regulation"))
            import regulation as reg
            rows, _ = reg.load_data()
            ref = next((r.reference for r in rows if r.locus == args.locus), None)
            if ref is None:
                raise ValueError(f"unknown GSE188264 locus: {args.locus}")
            rep = run_radar(
                ref, args.output, label=args.locus, circuit=args.circuit, protein=args.protein, target=args.target,
                signal_min=args.signal_min, signal_max=args.signal_max, signal_steps=args.steps,
                n_shuffles=args.shuffles, top_k=args.top_k, background_floor=args.background_floor,
                include_sequence_models=seq_on,
            )
            print(json.dumps({
                "output": str(args.output),
                "label": rep["label"],
                "circuit": rep["circuit"],
                "protein": rep["protein"],
                "background_floor": rep["background_floor"],
                "wt_critical_signal": rep["wt_critical_signal"],
                "category_counts": rep["scan"]["category_counts"],
                "doubly_silent_count": rep["scan"]["doubly_silent_count"],
                "featured_candidate": rep["featured_candidate"],
                "featured_rescues": rep["featured_rescues"],
                "null_controls": {
                    "dinucleotide_shuffle": {k: v for k, v in rep["null_controls"]["dinucleotide_shuffle"].items() if k != "decoys"},
                    "motif_pwm_scramble": rep["null_controls"]["motif_pwm_scramble"],
                },
                "all_checks_passed": rep["causal_controls"]["all_passed"],
            }, indent=2))
            return 0 if rep["causal_controls"]["all_passed"] else 1

        if args.command == "scan":
            ref = nc._read_dna_arg(str(args.fasta))
            rep = run_radar(
                ref, args.output, label=args.fasta.name, circuit=args.circuit, protein=args.protein, target=args.target,
                signal_min=args.signal_min, signal_max=args.signal_max, signal_steps=args.steps,
                n_shuffles=args.shuffles, top_k=args.top_k, background_floor=args.background_floor,
                include_sequence_models=seq_on,
            )
            print(json.dumps({
                "output": str(args.output),
                "label": rep["label"],
                "background_floor": rep["background_floor"],
                "category_counts": rep["scan"]["category_counts"],
                "doubly_silent_count": rep["scan"]["doubly_silent_count"],
                "featured_candidate": rep["featured_candidate"],
                "all_checks_passed": rep["causal_controls"]["all_passed"],
            }, indent=2))
            return 0 if rep["causal_controls"]["all_passed"] else 1

        if args.command == "variants":
            rep = analyze_variant_table(
                args.variants, args.output, circuit=args.circuit, protein=args.protein, target=args.target,
                signal_min=args.signal_min, signal_max=args.signal_max, signal_steps=args.steps,
                background_floor=args.background_floor, include_sequence_models=seq_on,
            )
            print(json.dumps({
                "output": str(args.output),
                "background_floor": rep["background_floor"],
                "variants_tested": rep["variants_tested"],
                "cryptic_count": rep["cryptic_count"],
                "doubly_silent_count": rep["doubly_silent_count"],
                "all_checks_passed": rep["causal_controls"]["all_passed"],
            }, indent=2))
            return 0 if rep["causal_controls"]["all_passed"] else 1

        bench = benchmark_radar(
            args.output, circuit=args.circuit, protein=args.protein, target=args.target,
            signal_min=args.signal_min, signal_max=args.signal_max, signal_steps=args.steps,
            n_shuffles=args.shuffles, background_floor=args.background_floor, include_sequence_models=seq_on,
        )
        print(json.dumps({k: v for k, v in bench.items() if k != "measured_variants"}, indent=2))
        return 0 if bench["causal_controls"]["all_passed"] else 1
    except (ValueError, OSError, RuntimeError, json.JSONDecodeError) as err:
        print(f"neural-v2f failed: {err}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

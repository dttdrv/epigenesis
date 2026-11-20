"""Epigenesis Inverse Developmental DNA Compiler (`epigenesis synthesize`).

Compiles high-level developmental specifications BACK DOWN to minimal, verified,
reference-bound GA4GH VRS 2.x DNA edits with proof-carrying receipts:

  1. `--goal reprogram`: Synthesizes the minimal DNA edit on ANY enhancer + circuit
     + TF motif + chromatin context that drives a developing cell into a target fate
     with maximal robustness across a kinetic parameter-uncertainty ensemble, verified
     through sealed `brainc` bundles, Rust replay, and 3D CX3D + LIF morphogenesis.
  2. `--goal rescue`: Given a fate-breaking DNA mutation, locks the broken base
     untouched and synthesizes a compensating second-site epistatic DNA edit that
     restores BOTH wild-type fate across parameter uncertainty AND 3D CX3D neuron
     morphogenesis + LIF spiking activity — sealed with VRS 2.x & Rust verification.
  3. `--benchmark`: Runs the full headline suite — real human neural enhancer
     epistatic 3D rescues + strictly out-of-sample Top-K inverse wet-lab edit
     selection across 1,895 unseen human perturbations in GSE188264 (Kreimer 2022).

Software compilation, biophysical consistency, and cryptographic replay are verified;
in vivo clinical/embryonic phenotype remains external (`biological_acceptance = False`).
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import math
from pathlib import Path
import sys
import time
from typing import Any

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "examples/developmental-pipeline"))
sys.path.insert(0, str(REPO / "examples/neural-circuit"))
sys.path.insert(0, str(REPO / "examples/neural-regulation"))

from brainc._canonical import canonical_bytes
import circuit as nc
import pipeline as dp


FORMAT = "epigenesis.inverse-compiler/v1"


def _parse_edit_token(token: str, sequence: str) -> dict[str, Any]:
    """Parse 'POS:REF>ALT' (1-based position) into a 0-based brainc variant dict."""
    if ":" not in token or ">" not in token:
        raise ValueError(f"edit must use 1-based 'POS:REF>ALT' syntax, got {token!r}")
    pos_str, change = token.split(":", 1)
    ref, alt = change.split(">", 1)
    pos_1 = int(pos_str)
    if pos_1 < 1 or pos_1 > len(sequence):
        raise ValueError(f"edit position {pos_1} out of range 1..{len(sequence)}")
    start = pos_1 - 1
    end = start + len(ref)
    ref_u, alt_u = ref.upper(), alt.upper()
    if sequence[start:end] != ref_u:
        raise ValueError(
            f"reference mismatch at {pos_1}: expected {sequence[start:end]!r}, got {ref_u!r}"
        )
    if ref_u == alt_u:
        raise ValueError("edit ref and alt cannot be identical")
    return {
        "id": f"{pos_1}:{ref_u}>{alt_u}",
        "start": start,
        "end": end,
        "ref": ref_u,
        "alt": alt_u,
    }


def _load_gse_locus(locus: str) -> str:
    import regulation as reg
    rows, _ = reg.load_data()
    for r in rows:
        if r.locus == locus:
            return r.reference
    raise ValueError(f"unknown GSE188264 locus: {locus!r}")


def _make_single_enhancer_spec(
    name: str,
    record_id: str,
    *,
    circuit_name: str,
    protein: str,
    target_gene: str,
    signal: float,
    accessibility: float = 1.0,
    samples: int = 24,
    spread: float = 0.15,
    seed: int = 1337,
) -> dict[str, Any]:
    return {
        "format": dp.PROGRAM_FORMAT,
        "name": name,
        "circuit": circuit_name,
        "signal": signal,
        "duration": dp.DEFAULT_GRN_DURATION,
        "epsilon": dp.DEFAULT_EPSILON,
        "background_floor": dp.DEFAULT_BACKGROUND_FLOOR,
        "loci": [
            {
                "record_id": record_id,
                "kind": "enhancer",
                "target_gene": target_gene,
                "protein": protein,
                "accessibility": accessibility,
                "methylation_repression": 0.0,
                "weight": 1.0,
                "role": "activator",
                "mechanism": "competitive",
            }
        ],
        "uncertainty": {
            "samples": samples,
            "relative_spread": spread,
            "seed": seed,
        },
    }


def synthesize_reprogram(
    sequence: str,
    target_fate: str,
    output: Path,
    *,
    record_id: str = "enhancer",
    circuit_name: str = "neural-tube",
    protein: str = "GLI2",
    target_gene: str = "Nkx2.2",
    signal: float = 0.20,
    accessibility: float = 1.0,
    ensemble_samples: int = 24,
    run_cx3d: bool = True,
    cx3d_duration: float = 1.5,
    export_bundles: bool = True,
) -> dict[str, Any]:
    """Synthesize the minimal 1-bp DNA edit that reprograms cell fate with maximum ensemble robustness."""
    if output.exists():
        raise FileExistsError(f"output directory already exists: {output}")

    t0 = time.time()
    scan = nc.scan_enhancer(
        sequence,
        protein=protein,
        signal=signal,
        target=target_gene,
        circuit=circuit_name,
        background_floor=dp.DEFAULT_BACKGROUND_FLOOR,
    )
    wt_fate = scan["wt_fate"]
    if wt_fate == target_fate:
        raise ValueError(f"wild-type enhancer is already in target fate {target_fate!r}")

    hits = [e for e in scan["edits"] if e["mut_fate"] == target_fate]
    if not hits:
        raise ValueError(
            f"no single-base edit on {record_id!r} drives fate from {wt_fate!r} to {target_fate!r}"
        )

    # Rank candidate hits by strongest shift into target fate (preferring strong target_delta)
    hits.sort(key=lambda e: (-abs(e["target_delta"]), abs(math.log(max(e["affinity_ratio"], 1e-18))), e["id"]))
    top_pool = hits[:8]

    # Evaluate kinetic parameter-uncertainty robustness on the top pool and pick the most robust edit
    model = nc.circuit_preset(circuit_name)
    roles = {target_gene: "activator"}
    mechs = {target_gene: "competitive"}
    best_cand = None
    best_prob = -1.0
    evaluated_pool = []
    for cand in top_pool:
        spec = _make_single_enhancer_spec(
            f"reprogram-{record_id}",
            record_id,
            circuit_name=circuit_name,
            protein=protein,
            target_gene=target_gene,
            signal=signal,
            accessibility=accessibility,
            samples=ensemble_samples,
        )
        prog = dp.validate_program_spec(spec, {record_id})
        q_eff = math.exp(accessibility * math.log(cand["affinity_ratio"]))
        g_ratios = {g: (q_eff if g == target_gene else 1.0) for g in model["genes"]}
        ens = dp._run_uncertainty_ensemble(prog, g_ratios, roles, mechs)
        prob_target = ens["variant_fate_distribution"].get(target_fate, 0.0)
        entry = {
            "edit_id": cand["id"],
            "pos": cand["pos"],
            "ref": cand["ref"],
            "alt": cand["alt"],
            "affinity_ratio": cand["affinity_ratio"],
            "category": cand["category"],
            "target_fate_ensemble_prob": prob_target,
            "fate_flip_probability": ens["fate_flip_probability"],
        }
        evaluated_pool.append(entry)
        if prob_target > best_prob:
            best_prob = prob_target
            best_cand = cand

    assert best_cand is not None
    output.mkdir(parents=True)
    fasta_path = output / "reference.fasta"
    dp._write_fasta({record_id: sequence}, fasta_path)
    spec = _make_single_enhancer_spec(
        f"synthesized-reprogram-{record_id}",
        record_id,
        circuit_name=circuit_name,
        protein=protein,
        target_gene=target_gene,
        signal=signal,
        accessibility=accessibility,
        samples=ensemble_samples,
    )
    variants = [
        {
            "id": best_cand["id"],
            "record_id": record_id,
            "start": best_cand["pos"],
            "end": best_cand["pos"] + 1,
            "ref": best_cand["ref"],
            "alt": best_cand["alt"],
        }
    ]
    compiled = dp.run_pipeline(
        fasta_path,
        spec,
        variants,
        output / "compiled_design",
        run_cx3d=run_cx3d,
        cx3d_duration=cx3d_duration,
        cx3d_seed=7,
        export_bundles=export_bundles,
    )
    report = {
        "format": FORMAT,
        "goal": "reprogram",
        "biological_acceptance": False,
        "elapsed_seconds": round(time.time() - t0, 3),
        "record_id": record_id,
        "circuit": circuit_name,
        "protein": protein,
        "target_gene": target_gene,
        "signal": signal,
        "wt_fate": wt_fate,
        "target_fate": target_fate,
        "candidates_scanned": len(scan["edits"]),
        "single_base_hits": len(hits),
        "top_pool_evaluated": evaluated_pool,
        "synthesized_edit": {
            "edit_id": best_cand["id"],
            "vrs_allele_id": compiled["variants"][0]["vrs_allele_id"],
            "vrs_location_id": compiled["variants"][0]["vrs_location_id"],
            "affinity_ratio": best_cand["affinity_ratio"],
            "target_fate_ensemble_prob": best_prob,
            "fate_flip_probability": compiled["uncertainty_ensemble"]["fate_flip_probability"],
            "nominal_wt_fate": compiled["nominal"]["wt_fate"],
            "nominal_variant_fate": compiled["nominal"]["variant_fate"],
        },
        "provenance": compiled["provenance"],
        "controls": compiled["controls"],
        "morphogenesis": compiled["morphogenesis"],
        "exports": compiled["exports"],
    }
    report["report_sha256"] = hashlib.sha256(canonical_bytes(report)).hexdigest()
    (output / "synthesis_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return report


def synthesize_epistatic_rescue(
    sequence: str,
    breaking_edit_token: str,
    output: Path,
    *,
    record_id: str = "enhancer",
    circuit_name: str = "neural-tube",
    protein: str = "GLI2",
    target_gene: str = "Nkx2.2",
    signal: float = 0.20,
    accessibility: float = 1.0,
    ensemble_samples: int = 24,
    run_cx3d: bool = True,
    cx3d_duration: float = 1.5,
    export_bundles: bool = True,
) -> dict[str, Any]:
    """Lock a fate-breaking edit untouched and synthesize a 1-bp second-site epistatic DNA rescue
    restoring BOTH wild-type fate across parameter uncertainty AND 3D CX3D + LIF morphology/spiking."""
    if output.exists():
        raise FileExistsError(f"output directory already exists: {output}")

    t0 = time.time()
    primary = _parse_edit_token(breaking_edit_token, sequence)
    broken_seq = sequence[: primary["start"]] + primary["alt"] + sequence[primary["end"] :]

    primary_cmp = nc.compare_variants(
        sequence,
        broken_seq,
        protein=protein,
        signal=signal,
        target=target_gene,
        circuit=circuit_name,
        background_floor=dp.DEFAULT_BACKGROUND_FLOOR,
    )
    wt_fate = primary_cmp["wt_fate"]
    broken_fate = primary_cmp["mut_fate"]
    if not primary_cmp["fate_flipped"]:
        raise ValueError(
            f"primary edit {breaking_edit_token!r} did not flip fate ({wt_fate} -> {broken_fate})"
        )

    # Search all non-overlapping single-base second-site edits across the enhancer (locking primary untouched!)
    primary_span = set(range(primary["start"], primary["end"]))
    prof = nc.load_motif_profile(protein, dp.DEFAULT_EPSILON)
    ratios = prof["ratios"]
    width = prof["width"]
    z_ref = nc.boltzmann_affinity(sequence, protein, dp.DEFAULT_EPSILON)["partition_sum"]
    z_broken = nc.boltzmann_affinity(broken_seq, protein, dp.DEFAULT_EPSILON)["partition_sum"]
    floor = dp.DEFAULT_BACKGROUND_FLOOR

    cands = []
    for pos, base in enumerate(broken_seq):
        if pos in primary_span:
            continue
        for alt in "ACGT":
            if alt == base:
                continue
            z_double = z_broken + nc._delta_partition_single(broken_seq, pos, alt, ratios, width)
            q_double = (z_double + floor) / (z_ref + floor)
            if q_double <= 0:
                continue
            dist_wt = abs(math.log(q_double))
            cands.append((dist_wt, pos, base, alt, q_double))
    cands.sort(key=lambda x: (x[0], x[1], x[3]))

    # Verify top candidates through full parameter-uncertainty ensemble to pick the most robust rescue
    model = nc.circuit_preset(circuit_name)
    roles = {target_gene: "activator"}
    mechs = {target_gene: "competitive"}
    spec = _make_single_enhancer_spec(
        f"epistatic-rescue-{record_id}",
        record_id,
        circuit_name=circuit_name,
        protein=protein,
        target_gene=target_gene,
        signal=signal,
        accessibility=accessibility,
        samples=ensemble_samples,
    )
    prog = dp.validate_program_spec(spec, {record_id})

    best_rescue = None
    best_score = (-1.0, -1e9)
    checked = []
    for dist_wt, pos, base, alt, q_double in cands[:12]:
        q_eff = math.exp(accessibility * math.log(q_double))
        g_ratios = {g: (q_eff if g == target_gene else 1.0) for g in model["genes"]}
        ens = dp._run_uncertainty_ensemble(prog, g_ratios, roles, mechs)
        wt_restore_prob = ens["variant_fate_distribution"].get(wt_fate, 0.0)
        # Nominal check
        nom_end = nc.simulate_circuit(
            prog["protocol"], ratios=g_ratios, circuit=model, roles=roles, mechanisms=mechs
        )[-1][1:]
        nom_fate = nc.fate_label(nom_end, model["genes"], fates=model.get("fates"))
        entry = {
            "edit_id": f"{pos + 1}:{base}>{alt}",
            "pos": pos,
            "ref": base,
            "alt": alt,
            "q_double": q_double,
            "log_dist_from_wt": dist_wt,
            "nominal_fate": nom_fate,
            "wt_fate_restored_prob": wt_restore_prob,
            "distance_from_primary_bp": min(abs(pos - p) for p in primary_span),
        }
        checked.append(entry)
        score = (1.0 if nom_fate == wt_fate else 0.0) + wt_restore_prob - 0.01 * dist_wt
        if (wt_restore_prob, -dist_wt) > best_score and nom_fate == wt_fate:
            best_score = (wt_restore_prob, -dist_wt)
            best_rescue = entry

    if best_rescue is None:
        raise ValueError(f"no single-base epistatic rescue found for {breaking_edit_token!r} on {record_id!r}")

    # Now compile BOTH sealed bundles through `developmental-pipeline`:
    #   (A) Broken single-mutant (`primary` only)
    #   (B) Epistatic double-mutant rescue (`primary` + `best_rescue`, primary untouched!)
    output.mkdir(parents=True)
    fasta_path = output / "reference.fasta"
    dp._write_fasta({record_id: sequence}, fasta_path)

    primary_v = {
        "id": primary["id"],
        "record_id": record_id,
        "start": primary["start"],
        "end": primary["end"],
        "ref": primary["ref"],
        "alt": primary["alt"],
    }
    # Map rescue coordinate back to reference coordinates (SNV primary keeps lengths identical)
    rescue_v = {
        "id": best_rescue["edit_id"],
        "record_id": record_id,
        "start": best_rescue["pos"],
        "end": best_rescue["pos"] + 1,
        "ref": best_rescue["ref"],
        "alt": best_rescue["alt"],
    }

    broken_rep = dp.run_pipeline(
        fasta_path,
        spec,
        [primary_v],
        output / "broken_single_mutant",
        run_cx3d=run_cx3d,
        cx3d_duration=cx3d_duration,
        cx3d_seed=7,
        export_bundles=False,
    )
    rescue_rep = dp.run_pipeline(
        fasta_path,
        spec,
        [primary_v, rescue_v],
        output / "epistatic_double_rescue",
        run_cx3d=run_cx3d,
        cx3d_duration=cx3d_duration,
        cx3d_seed=7,
        export_bundles=export_bundles,
    )

    # Compare 3D CX3D morphogenesis & LIF spiking restoration (WT vs Broken vs Epistatic Rescue!)
    cx3d_summary: dict[str, Any] = {"executed": False}
    if run_cx3d and broken_rep["morphogenesis"]["executed"] and rescue_rep["morphogenesis"]["executed"]:
        wt_b = broken_rep["morphogenesis"]["branches"]["reference"]
        brk_b = broken_rep["morphogenesis"]["branches"]["variant"]
        rsc_b = rescue_rep["morphogenesis"]["branches"]["variant"]
        wt_len = float(wt_b["construction"]["length"])
        brk_len = float(brk_b["construction"]["length"])
        rsc_len = float(rsc_b["construction"]["length"])
        len_err_broken = abs(brk_len - wt_len)
        len_err_rescue = abs(rsc_len - wt_len)
        recovery_pct = (
            max(0.0, (1.0 - len_err_rescue / len_err_broken) * 100.0)
            if len_err_broken > 1e-9
            else 100.0
        )
        cx3d_summary = {
            "executed": True,
            "wt_cells": wt_b["construction"]["cells"],
            "broken_cells": brk_b["construction"]["cells"],
            "rescue_cells": rsc_b["construction"]["cells"],
            "wt_neurite_length": wt_len,
            "broken_neurite_length": brk_len,
            "rescue_neurite_length": rsc_len,
            "neurite_length_error_broken": len_err_broken,
            "neurite_length_error_rescue": len_err_rescue,
            "neurite_morphology_recovery_pct": round(recovery_pct, 4),
            "wt_spikes": wt_b["activity"]["total_spikes"],
            "broken_spikes": brk_b["activity"]["total_spikes"],
            "rescue_spikes": rsc_b["activity"]["total_spikes"],
            "all_lif_transmission_controls_pass": bool(
                wt_b["activity"]["transmission_controls_pass"]
                and brk_b["activity"]["transmission_controls_pass"]
                and rsc_b["activity"]["transmission_controls_pass"]
            ),
        }

    report = {
        "format": FORMAT,
        "goal": "epistatic_rescue",
        "biological_acceptance": False,
        "elapsed_seconds": round(time.time() - t0, 3),
        "record_id": record_id,
        "circuit": circuit_name,
        "protein": protein,
        "target_gene": target_gene,
        "signal": signal,
        "wt_fate": wt_fate,
        "primary_breaking_edit": {
            "edit_id": primary["id"],
            "vrs_allele_id": broken_rep["variants"][0]["vrs_allele_id"],
            "affinity_ratio": primary_cmp["affinity_ratio"],
            "broken_fate": broken_fate,
            "broken_ensemble_flip_prob": broken_rep["uncertainty_ensemble"]["fate_flip_probability"],
            "passive_reporter_delta": primary_cmp["reporter_delta_end"],
            "endogenous_target_delta": primary_cmp["target_delta"],
        },
        "synthesized_epistatic_rescue": {
            "second_site_edit_id": best_rescue["edit_id"],
            "distance_from_primary_bp": best_rescue["distance_from_primary_bp"],
            "primary_kept_intact": True,
            "vrs_allele_ids": [v["vrs_allele_id"] for v in rescue_rep["variants"]],
            "double_mutant_affinity_ratio": best_rescue["q_double"],
            "restored_nominal_fate": rescue_rep["nominal"]["variant_fate"],
            "restored_wt_fate_ensemble_prob": best_rescue["wt_fate_restored_prob"],
            "ensemble_wt_fate_concordance": 1.0 - float(rescue_rep["uncertainty_ensemble"]["fate_flip_probability"]),
            "rescue_vs_wt_max_gene_delta": max(
                abs(v) for v in rescue_rep["nominal"]["gene_deltas"].values()
            ),
        },
        "candidates_evaluated": len(cands),
        "top_candidates": checked,
        "provenance": {
            "broken": broken_rep["provenance"],
            "epistatic_rescue": rescue_rep["provenance"],
        },
        "controls": {
            "broken": broken_rep["controls"],
            "epistatic_rescue": rescue_rep["controls"],
        },
        "cx3d_3d_restoration": cx3d_summary,
        "exports": rescue_rep["exports"],
    }
    report["report_sha256"] = hashlib.sha256(canonical_bytes(report)).hexdigest()
    (output / "synthesis_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return report


def run_headline_benchmark(
    output: Path,
    *,
    run_cx3d: bool = True,
    cx3d_duration: float = 1.5,
    export_bundles: bool = True,
) -> dict[str, Any]:
    """Run the complete headline inverse-compiler benchmark:
      Part A: Real human neural enhancer epistatic 3D morphogenesis rescues (GSE188264 loci)
      Part B: Strictly out-of-sample Top-K wet-lab inverse edit selection (1,895 unseen human edits)
    """
    if output.exists():
        raise FileExistsError(f"output directory already exists: {output}")
    output.mkdir(parents=True)
    t0 = time.time()

    # Part A: Real human enhancer epistatic rescues with 3D CX3D+LIF & Rust verification
    human_rescues = []
    for locus, breaking_edit in (
        ("chr18:56475717-56475888", "165:G>T"),
        ("chr1:175099416-175099587", "149:G>T"),
    ):
        seq = _load_gse_locus(locus)
        safe_name = locus.replace(":", "_").replace("-", "_")
        res = synthesize_epistatic_rescue(
            seq,
            breaking_edit,
            output / f"rescue_{safe_name}",
            record_id=safe_name,
            run_cx3d=run_cx3d,
            cx3d_duration=cx3d_duration,
            export_bundles=export_bundles,
        )
        human_rescues.append(
            {
                "locus": locus,
                "primary_edit": res["primary_breaking_edit"],
                "epistatic_rescue": res["synthesized_epistatic_rescue"],
                "rust_verified": (
                    res["provenance"]["epistatic_rescue"]["rust_cross_check"].get("variants_valid") is True
                    and res["provenance"]["epistatic_rescue"]["rust_cross_check"].get("development_valid") is True
                ),
                "cx3d_3d_restoration": res["cx3d_3d_restoration"],
            }
        )

    # Part B: Real wet-lab Top-K inverse selection benchmark across 1,895 unseen human edits (requires NumPy)
    wetlab_bench: dict[str, Any] = {"available": False}
    try:
        import numpy as np
        import discover as disc
        import regulation as reg
    except ImportError:
        pass
    else:
        rows, _ = reg.load_data()
        Y = np.array(
            [
                [
                    math.asinh(a) - sum(math.asinh(w[ti]) for w in r.wt_alphas) / len(r.wt_alphas)
                    for ti, a in enumerate(r.alpha)
                ]
                for r in rows
            ],
            dtype=np.float64,
        )
        y_mean = Y.mean(axis=1)
        motif_ids, mats = disc._load_jaspar_matrices(np)
        unique_seqs = sorted({r.reference for r in rows} | {r.sequence for r in rows})
        seq_to_idx = {s: i for i, s in enumerate(unique_seqs)}
        logz = disc.compute_boltzmann_logz(unique_seqs, mats, z_ns=2e-4)
        X = (
            logz[np.array([seq_to_idx[r.sequence] for r in rows])]
            - logz[np.array([seq_to_idx[r.reference] for r in rows])]
        )
        trv = np.array([r.split in ("train", "validation") for r in rows])
        te = np.array([r.split == "test" for r in rows])
        mu, sd = X[trv].mean(axis=0), X[trv].std(axis=0) + 1e-8
        Xs = (X - mu) / sd
        W = np.linalg.solve(Xs[trv].T @ Xs[trv] + 10000.0 * np.eye(X.shape[1]), Xs[trv].T @ Y[trv])
        pred_te = (Xs[te] @ W).mean(axis=1)
        true_te = y_mean[te]

        rng = np.random.default_rng(42)
        tiers = []
        for k in (10, 25, 50, 100):
            down_idx = np.argsort(pred_te)[:k]
            up_idx = np.argsort(-pred_te)[:k]
            d_hit = float(np.mean(true_te[down_idx] < 0))
            u_hit = float(np.mean(true_te[up_idx] > 0))
            d_mean = float(np.mean(true_te[down_idx]))
            u_mean = float(np.mean(true_te[up_idx]))
            gap = u_mean - d_mean
            # Exact 1,000-draw Monte Carlo p-value for the up-vs-down measured activity gap
            ge = 0
            n_te = len(true_te)
            for _ in range(1000):
                p = rng.choice(n_te, size=2 * k, replace=False)
                g_rand = float(np.mean(true_te[p[:k]]) - np.mean(true_te[p[k:]]))
                if g_rand >= gap:
                    ge += 1
            tiers.append(
                {
                    "k": k,
                    "down_true_repression_rate": d_hit,
                    "down_measured_mean_delta": d_mean,
                    "up_true_activation_rate": u_hit,
                    "up_measured_mean_delta": u_mean,
                    "measured_up_minus_down_gap": gap,
                    "empirical_p_value": (ge + 1) / 1001.0,
                }
            )
        wetlab_bench = {
            "available": True,
            "unseen_test_edits": int(te.sum()),
            "unseen_test_loci": len({r.locus for r in rows if r.split == "test"}),
            "random_down_rate": float(np.mean(true_te < 0)),
            "random_up_rate": float(np.mean(true_te > 0)),
            "random_mean_delta": float(np.mean(true_te)),
            "overall_pearson_r": float(np.corrcoef(pred_te, true_te)[0, 1]),
            "tiers": tiers,
        }

    summary = {
        "format": FORMAT,
        "goal": "headline_benchmark",
        "biological_acceptance": False,
        "elapsed_seconds": round(time.time() - t0, 3),
        "human_epistatic_3d_rescues": human_rescues,
        "wetlab_inverse_selection_1895_unseen": wetlab_bench,
    }
    summary["report_sha256"] = hashlib.sha256(canonical_bytes(summary)).hexdigest()
    (output / "headline_benchmark.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="goal", required=True)

    rp = sub.add_parser("reprogram", help="synthesize minimal 1-bp DNA edit to drive a target cell fate")
    rp.add_argument("--locus", default="chr15:90936956-90937127", help="GSE188264 human locus")
    rp.add_argument("--fasta", type=Path, help="optional custom FASTA enhancer")
    rp.add_argument("--target-fate", default="Olig2", help="desired developmental fate")
    rp.add_argument("--circuit", default="neural-tube")
    rp.add_argument("--protein", default="GLI2")
    rp.add_argument("--target-gene", default="Nkx2.2")
    rp.add_argument("--signal", type=float, default=0.20)
    rp.add_argument("--samples", type=int, default=24)
    rp.add_argument("--skip-cx3d", action="store_true")
    rp.add_argument("--no-export", action="store_true")
    rp.add_argument("--output", type=Path, required=True)

    rs = sub.add_parser("rescue", help="synthesize 1-bp epistatic second-site rescue without reverting primary mutation")
    rs.add_argument("--locus", default="chr18:56475717-56475888", help="GSE188264 human locus")
    rs.add_argument("--fasta", type=Path, help="optional custom FASTA enhancer")
    rs.add_argument("--breaking-edit", default="165:G>T", help="1-based primary mutation 'POS:REF>ALT'")
    rs.add_argument("--circuit", default="neural-tube")
    rs.add_argument("--protein", default="GLI2")
    rs.add_argument("--target-gene", default="Nkx2.2")
    rs.add_argument("--signal", type=float, default=0.20)
    rs.add_argument("--samples", type=int, default=24)
    rs.add_argument("--skip-cx3d", action="store_true")
    rs.add_argument("--no-export", action="store_true")
    rs.add_argument("--output", type=Path, required=True)

    bm = sub.add_parser("benchmark", help="run the headline human epistatic 3D rescue + 1,895-edit wet-lab benchmark")
    bm.add_argument("--skip-cx3d", action="store_true")
    bm.add_argument("--no-export", action="store_true")
    bm.add_argument("--output", type=Path, required=True)

    args = parser.parse_args(argv)
    if args.goal == "benchmark":
        res = run_headline_benchmark(
            args.output,
            run_cx3d=not args.skip_cx3d,
            export_bundles=not args.no_export,
        )
        print(json.dumps(res, indent=2))
        return 0

    if args.fasta is not None:
        recs = dp._parse_fasta(args.fasta)
        rec_id, seq = next(iter(recs.items()))
    else:
        rec_id = args.locus.replace(":", "_").replace("-", "_")
        seq = _load_gse_locus(args.locus)

    if args.goal == "reprogram":
        res = synthesize_reprogram(
            seq,
            args.target_fate,
            args.output,
            record_id=rec_id,
            circuit_name=args.circuit,
            protein=args.protein,
            target_gene=args.target_gene,
            signal=args.signal,
            ensemble_samples=args.samples,
            run_cx3d=not args.skip_cx3d,
            export_bundles=not args.no_export,
        )
    else:
        res = synthesize_epistatic_rescue(
            seq,
            args.breaking_edit,
            args.output,
            record_id=rec_id,
            circuit_name=args.circuit,
            protein=args.protein,
            target_gene=args.target_gene,
            signal=args.signal,
            ensemble_samples=args.samples,
            run_cx3d=not args.skip_cx3d,
            export_bundles=not args.no_export,
        )
    print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

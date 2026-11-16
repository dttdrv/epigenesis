"""The single frozen test scoring for PLAN.md. Refuses to run twice or on altered inputs."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import platform
import statistics
import time

import numpy as np
import scipy
from scipy.stats import mannwhitneyu

import motif_temporal as mt
import regulation as reg
import sites
import thermo as th

LINEAR_PENALTIES = (0.01, 0.1, 1.0, 10.0, 100.0)


def _verify_freeze(selection_path: Path) -> dict:
    freeze_path = th.HERE / "FREEZE.json"
    freeze = json.loads(freeze_path.read_text())
    expected = {"PLAN.md": freeze["plan_sha256"], "thermo.py": freeze["thermo_sha256"],
                "sites.py": freeze["sites_sha256"], "confirm.py": freeze["confirm_sha256"]}
    for name, digest in expected.items():
        if th.sha256_file(th.HERE / name) != digest:
            raise ValueError(f"frozen file changed after freeze: {name}")
    if th.sha256_file(selection_path) != freeze["selection_sha256"]:
        raise ValueError("selection.json differs from the frozen selection")
    return freeze


def _predictions(pred: np.ndarray, rows: list[reg.Row]) -> dict[str, list[float]]:
    return {r.id: [float(v) for v in pred[i]] for i, r in enumerate(rows)}


def _center_of_mass(a: np.ndarray) -> np.ndarray:
    hours = np.array(reg.TIMES, dtype=np.float64)
    total = a.sum(axis=1)
    com = np.full(a.shape[0], np.nan)
    ok = total > 1e-9
    com[ok] = (a[ok] * hours).sum(axis=1) / total[ok]
    return com


def run_confirm(output: Path, selection_path: Path) -> dict:
    freeze = _verify_freeze(selection_path)
    selection = json.loads(selection_path.read_text())
    chosen = selection["chosen"]
    kappa, lam, rho = chosen["kappa"], chosen["lambda"], chosen["rho"]
    output.mkdir(parents=True)
    started = time.time()

    rows, qc = reg.load_data()
    training = [r for r in rows if r.split == "train"]
    validation_t = [r for r in rows if r.split == "validation" and r.kind == "targeted"]
    fitting = [r for r in rows if r.split != "test"]
    testing = [r for r in rows if r.split == "test"]
    primary = [r for r in testing if r.kind == "targeted"]
    sel_scale = reg.training_scale(training)
    scale = reg.training_scale(fitting)
    ids, matrices, names = th.load_motifs()

    # ---- baselines, refitted by the repository's own code (gate G0)
    model, coefficients = mt.fit(fitting, scale)
    man, mats, _ = mt._load_jaspar()
    mt_pred = mt.predict_rows(model, coefficients, testing, mats)
    k2_model, _ = mt._fit_kmer_baseline(fitting, 2, reg.PENALTIES, scale, validation_t, training, sel_scale)
    k5_model, _ = mt._fit_kmer_baseline(fitting, 5, reg.PENALTIES, scale, validation_t, training, sel_scale)
    k2_pred = np.array([reg.predict(k2_model, r.reference, r.sequence) for r in testing])
    k5_pred = np.array([reg.predict(k5_model, r.reference, r.sequence) for r in testing])
    zero_pred = np.zeros((len(testing), th.T))
    is_primary = np.array([r.kind == "targeted" for r in testing])

    def family_scores(pred):
        return th.row_losses(pred[is_primary], primary, scale)

    base = {"zero_effect": family_scores(zero_pred), "dinucleotide": family_scores(k2_pred),
            "k5": family_scores(k5_pred), "motif_temporal": family_scores(mt_pred)}
    base_mse = {k: statistics.mean(v.values()) for k, v in base.items()}
    g0 = {k: {"computed": base_mse[k], "published": th.PUBLISHED[k], "pass": abs(base_mse[k] - th.PUBLISHED[k]) < 1e-5}
          for k in th.PUBLISHED}
    print("baselines", {k: round(v, 6) for k, v in base_mse.items()}, flush=True)

    # ---- the thermodynamic model on real PWMs
    index, table = th.theta_table(rows, matrices, (kappa,))
    var, par, y, w = th.design(fitting, index, table[0], scale)
    tvar, tpar, _, _ = th.design(testing, index, table[0], scale)
    fitted = th.fit_thermo(var, par, y, w, lam, rho)
    thermo_pred = th.predict_effects(fitted["A"], fitted["R"], tvar, tpar)
    real = family_scores(thermo_pred)
    real_mse = statistics.mean(real.values())
    print("thermo test MSE", round(real_mse, 6), "converged", fitted["converged"], fitted["iterations"], flush=True)

    # ---- linear same-feature control (report-only)
    trvar, trpar, try_, trw = th.design(training, index, table[0], sel_scale)
    vvar, vpar, _, _ = th.design(validation_t, index, table[0], sel_scale)
    lin_choices = []
    for penalty in LINEAR_PENALTIES:
        coef = th.fit_linear(trvar, trpar, try_, trw, penalty)
        lin_choices.append((th.mean_loss((vvar - vpar) @ coef, validation_t, sel_scale), penalty))
    lin_penalty = min(lin_choices)[1]
    lin_coef = th.fit_linear(var, par, y, w, lin_penalty)
    lin_pred = (tvar - tpar) @ lin_coef
    linear = family_scores(lin_pred)

    # ---- PWM-shuffle nulls (decisive)
    nulls = []
    for seed in th.NULL_SEEDS:
        shuffled = sites.shuffle_columns(matrices, seed)
        nidx, ntab = th.theta_table(rows, shuffled, (kappa,))
        nvar, npar, ny, nw = th.design(fitting, nidx, ntab[0], scale)
        ntvar, ntpar, _, _ = th.design(testing, nidx, ntab[0], scale)
        nfit = th.fit_thermo(nvar, npar, ny, nw, lam, rho)
        nscores = family_scores(th.predict_effects(nfit["A"], nfit["R"], ntvar, ntpar))
        nulls.append({"seed": seed, "mse": statistics.mean(nscores.values()), "converged": nfit["converged"],
                      "scores": nscores})
        print(f"null seed {seed}: {nulls[-1]['mse']:.6f}", flush=True)
    null_avg = {g: statistics.mean(n["scores"][g] for n in nulls) for g in real}

    # ---- gates
    gains = {name: th.bootstrap_gain(scores, real) for name, scores in
             (("vs_zero_effect", base["zero_effect"]), ("vs_dinucleotide", base["dinucleotide"]), ("vs_k5", base["k5"]),
              ("vs_motif_temporal", base["motif_temporal"]), ("vs_linear_same_features", linear),
              ("vs_pwm_shuffle_null_average", null_avg))}
    beaten = sum(real_mse < n["mse"] for n in nulls)
    gates = {
        "G0_integrity": all(v["pass"] for v in g0.values()),
        "G1_vs_zero": gains["vs_zero_effect"]["lower"] > 0,
        "G2_vs_k5": gains["vs_k5"]["lower"] > 0,
        "G3a_vs_null_average": gains["vs_pwm_shuffle_null_average"]["lower"] > 0,
        "G3b_beats_all_nulls": beaten == len(nulls),
        "G4_noninferior_motif_temporal": gains["vs_motif_temporal"]["gain"] >= -0.005,
    }
    supported = all(gates.values())
    stretch = {"G5_superior_to_motif_temporal": gains["vs_motif_temporal"]["lower"] > 0}

    # ---- exploratory: dynamics sanity (G6)
    a, r = fitted["A"], fitted["R"]
    com = _center_of_mass(a)
    early = [com[i] for i, n in enumerate(names) if n in th.EARLY_SET and np.isfinite(com[i])]
    neural = [com[i] for i, n in enumerate(names) if n in th.NEURAL_SET and np.isfinite(com[i])]
    g6 = {"early_n": len(early), "neural_n": len(neural), "early_median_hours": float(np.median(early)) if early else None,
          "neural_median_hours": float(np.median(neural)) if neural else None}
    if len(early) >= 3 and len(neural) >= 3:
        test = mannwhitneyu(early, neural, alternative="less")
        g6.update({"mannwhitney_u": float(test.statistic), "p_one_sided_early_earlier": float(test.pvalue),
                   "pass": bool(test.pvalue < 0.05)})
    top = []
    for t, hour in enumerate(reg.TIMES):
        order = np.argsort(a[:, t])[::-1][:5]
        rorder = np.argsort(r[:, t])[::-1][:5]
        top.append({"hour": hour, "top_activators": [(names[i], ids[i], float(a[i, t])) for i in order],
                    "top_quenchers": [(names[i], ids[i], float(r[i, t])) for i in rorder]})

    categories = {}
    for kind in ("targeted", "random", "scrambled"):
        sel = np.array([x.kind == kind for x in testing])
        sub = [x for x in testing if x.kind == kind]
        categories[kind] = {name: th.mean_loss(pred[sel], sub, scale) for name, pred in
                            (("thermo", thermo_pred), ("motif_temporal", mt_pred), ("k5", k5_pred), ("zero_effect", zero_pred))}

    # ---- artifacts
    np.savez_compressed(output / "strengths.npz", A=a, R=r, motif_ids=np.array(ids), times_hours=np.array(reg.TIMES))
    with (output / "tf-strengths.csv").open("x", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("motif_id", "name", "kind", *[f"{h}h" for h in reg.TIMES]))
        for i, (mid, name) in enumerate(zip(ids, names)):
            writer.writerow((mid, name, "activator", *[f"{v:.6g}" for v in a[i]]))
            writer.writerow((mid, name, "quencher", *[f"{v:.6g}" for v in r[i]]))
    with (output / "test-predictions.csv").open("x", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("id", "group", "kind", "hour", "observed_effect", "thermo", "motif_temporal", "k5"))
        for i, row in enumerate(testing):
            obs = reg.target(row, scale)
            for h, hour in enumerate(reg.TIMES):
                writer.writerow((row.id, row.group, row.kind, hour, f"{obs[h]:.6g}", f"{thermo_pred[i, h]:.6g}",
                                 f"{mt_pred[i, h]:.6g}", f"{k5_pred[i, h]:.6g}"))
    report = {
        "status": "evaluation-complete", "plan_sha256": freeze["plan_sha256"], "selection": chosen,
        "scope": "retrospective held-out-family reporter effects; test families previously examined; one assay",
        "test_groups": len(real), "test_targeted_rows": len(primary), "final_scale": scale,
        "baseline_mse": base_mse, "g0_integrity": g0, "thermo_mse": real_mse,
        "thermo_fit": {k: fitted[k] for k in ("loss", "iterations", "converged", "message")},
        "linear_same_features": {"penalty": lin_penalty, "mse": statistics.mean(linear.values())},
        "pwm_shuffle_nulls": [{k: v for k, v in n.items() if k != "scores"} for n in nulls],
        "null_mean_mse": statistics.mean(n["mse"] for n in nulls),
        "null_min_mse": min(n["mse"] for n in nulls), "nulls_beaten_by_real": beaten,
        "gains_95_paired_family_bootstrap": gains, "gates": gates, "stretch": stretch,
        "thermodynamic_transfer_supported": supported,
        "g6_dynamics_exploratory": g6, "top_motifs_by_hour": top, "categories": categories,
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "runtime_seconds": time.time() - started,
        "artifacts": {n: th.sha256_file(output / n) for n in ("strengths.npz", "tf-strengths.csv", "test-predictions.csv")},
    }
    (output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"gates": gates, "supported": supported, "gains": gains, "g6": g6}, indent=2))
    return report

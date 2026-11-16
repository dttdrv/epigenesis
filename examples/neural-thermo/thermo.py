"""Thermodynamic transfer test on the GSE188264 neural-induction reporter data.

Implements the model and protocol frozen in PLAN.md. Commands:

  select   fit the frozen grid on train, score validation only, write selection.json
  confirm  refit on train+validation, score the test families once (refuses to rerun)

Outputs are retrospective held-out-family reporter effects, not endogenous
expression or cell fates.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import itertools
import json
import math
from pathlib import Path
import random
import statistics
import sys
import time

import numpy as np
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
for extra in (REPO, HERE, REPO / "examples/neural-regulation", REPO / "studio/epigenesis_studio/model"):
    sys.path.insert(0, str(extra))

import motif_temporal as mt  # noqa: E402
import regulation as reg  # noqa: E402
import sites  # noqa: E402

LN2 = math.log(2.0)
T = len(reg.TIMES)
KAPPAS = (1e-4, 1e-3, 1e-2, 1e-1)
LAMBDAS = (1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
RHOS = (0.0, 10.0, 100.0)
MAXITER = 500
FTOL = 1e-10
BOOTSTRAPS = 10_000
BOOTSTRAP_SEED = 20261116
NULL_SEEDS = tuple(range(1, 21))
PUBLISHED = {"zero_effect": 0.360870, "dinucleotide": 0.319358, "k5": 0.285185, "motif_temporal": 0.252079}
EARLY_SET = ("POU2F1::SOX2", "KLF2", "KLF17", "TFAP2C", "ESRRB", "TBX3", "FOXD3", "SOX15")
NEURAL_SET = ("PAX6", "ZIC1", "ZIC4", "SOX1", "SOX21", "SOX4", "NEUROG1", "NEUROG2", "ASCL1", "POU3F2", "POU3F1",
              "EMX2", "FOXG1", "LHX5", "SIX3", "OLIG2", "HES5", "RFX4", "NR2F1", "NFIA", "NFIX")


# ---------------------------------------------------------------- features

def load_motifs() -> tuple[list[str], list[np.ndarray], list[str]]:
    manifest, matrices, _ = mt._load_jaspar()
    names = json.loads((HERE / "data/motif_names.json").read_text())
    return manifest["motif_ids"], matrices, [names[i] for i in manifest["motif_ids"]]


def theta_table(rows: list[reg.Row], matrices: list[np.ndarray], kappas) -> tuple[dict[str, int], np.ndarray]:
    """Theta for every distinct parent and variant sequence: shape (len(kappas), n_seq, n_motif)."""
    unique = sorted({row.sequence for row in rows} | {row.reference for row in rows})
    return {s: i for i, s in enumerate(unique)}, sites.site_sums(unique, matrices, kappas)


def design(rows: list[reg.Row], index: dict[str, int], table: np.ndarray, scale: float):
    var = table[np.array([index[r.sequence] for r in rows])]
    par = table[np.array([index[r.reference] for r in rows])]
    y = np.array([reg.target(r, scale) for r in rows], dtype=np.float64)
    counts = Counter(r.group for r in rows)
    w = np.array([1.0 / (len(counts) * counts[r.group]) for r in rows], dtype=np.float64)
    return var, par, y, w


# ------------------------------------------------------------------- model

def _smooth_matrix(rho: float) -> np.ndarray:
    d = np.zeros((T - 1, T))
    d[np.arange(T - 1), np.arange(T - 1)] = -1.0
    d[np.arange(T - 1), np.arange(1, T)] = 1.0
    return rho * (d.T @ d)


def _unpack(x: np.ndarray, m: int) -> tuple[np.ndarray, np.ndarray]:
    half = m * T
    return x[:half].reshape(m, T), x[half:].reshape(m, T)


def predict_effects(a: np.ndarray, r: np.ndarray, var: np.ndarray, par: np.ndarray) -> np.ndarray:
    def f(theta):
        return (np.log1p(theta @ a) - np.log1p(theta @ r)) / LN2
    return f(var) - f(par)


def objective(x, var, par, y, w, lam, smooth, m):
    a, r = _unpack(x, m)
    ua_v, ua_p = var @ a, par @ a
    ur_v, ur_p = var @ r, par @ r
    pred = (np.log1p(ua_v) - np.log1p(ur_v) - np.log1p(ua_p) + np.log1p(ur_p)) / LN2
    err = pred - y
    loss = float((w[:, None] * err * err).sum() / T)
    g = (2.0 / T) * w[:, None] * err / LN2
    grad_a = var.T @ (g / (1 + ua_v)) - par.T @ (g / (1 + ua_p))
    grad_r = -(var.T @ (g / (1 + ur_v)) - par.T @ (g / (1 + ur_p)))
    sa, sr = a @ smooth, r @ smooth
    loss += lam * float((a * a).sum() + (r * r).sum()) + lam * float((a * sa).sum() + (r * sr).sum())
    grad_a = grad_a + 2 * lam * a + 2 * lam * sa
    grad_r = grad_r + 2 * lam * r + 2 * lam * sr
    return loss, np.concatenate((grad_a.ravel(), grad_r.ravel()))


def fit_thermo(var, par, y, w, lam: float, rho: float, init: np.ndarray | None = None) -> dict:
    m = var.shape[1]
    smooth = _smooth_matrix(rho)
    x0 = np.full(2 * m * T, 0.01) if init is None else init
    res = minimize(objective, x0, args=(var, par, y, w, lam, smooth, m), jac=True, method="L-BFGS-B",
                   bounds=[(0.0, None)] * x0.size, options={"maxiter": MAXITER, "ftol": FTOL, "maxcor": 20})
    a, r = _unpack(res.x, m)
    return {"A": a, "R": r, "x": res.x, "loss": float(res.fun), "iterations": int(res.nit),
            "converged": bool(res.success), "message": str(res.message)}


def fit_linear(var, par, y, w, penalty: float) -> np.ndarray:
    x = var - par
    xw = x * w[:, None]
    gram = xw.T @ x
    gram[np.diag_indices_from(gram)] += penalty
    return np.linalg.solve(gram, xw.T @ y)


# -------------------------------------------------------------- evaluation

def row_losses(pred: np.ndarray, rows: list[reg.Row], scale: float) -> dict[str, float]:
    return reg.family_losses(rows, {r.id: [float(v) for v in pred[i]] for i, r in enumerate(rows)}, scale)


def mean_loss(pred: np.ndarray, rows: list[reg.Row], scale: float) -> float:
    return statistics.mean(row_losses(pred, rows, scale).values())


def bootstrap_gain(baseline: dict[str, float], candidate: dict[str, float]) -> dict:
    keys = sorted(baseline)
    if keys != sorted(candidate):
        raise ValueError("paired family scores must share keys")
    d = np.array([baseline[k] - candidate[k] for k in keys])
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    means = np.sort(d[rng.integers(0, len(d), size=(BOOTSTRAPS, len(d)))].mean(axis=1))
    return {"gain": float(d.mean()), "lower": float(means[int(0.025 * BOOTSTRAPS)]),
            "upper": float(means[int(0.975 * BOOTSTRAPS)]), "families": len(d)}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ------------------------------------------------------------------ select

def cmd_select(output: Path) -> dict:
    output.mkdir(parents=True)
    rows, _ = reg.load_data()
    training = [r for r in rows if r.split == "train"]
    validation = [r for r in rows if r.split == "validation" and r.kind == "targeted"]
    ids, matrices, names = load_motifs()
    scale = reg.training_scale(training)
    t0 = time.time()
    index, table = theta_table(training + validation, matrices, KAPPAS)
    print(f"features {table.shape} in {time.time() - t0:.1f}s", flush=True)
    results = []
    for k, kappa in enumerate(KAPPAS):
        var, par, y, w = design(training, index, table[k], scale)
        vvar, vpar, _, _ = design(validation, index, table[k], scale)
        for rho in RHOS:
            init = None
            for lam in sorted(LAMBDAS, reverse=True):
                t1 = time.time()
                fitted = fit_thermo(var, par, y, w, lam, rho, init)
                init = fitted["x"]  # warm start along the decreasing-lambda path
                val = mean_loss(predict_effects(fitted["A"], fitted["R"], vvar, vpar), validation, scale)
                results.append({"kappa": kappa, "lambda": lam, "rho": rho, "validation_mse": val,
                                "train_objective": fitted["loss"], "iterations": fitted["iterations"],
                                "converged": fitted["converged"]})
                print(f"kappa={kappa:g} rho={rho:g} lambda={lam:g} val={val:.6f} "
                      f"obj={fitted['loss']:.5f} it={fitted['iterations']} {time.time() - t1:.1f}s", flush=True)
    chosen = min(results, key=lambda c: (c["validation_mse"], -c["lambda"], c["kappa"]))
    zero = mean_loss(np.zeros((len(validation), T)), validation, scale)
    report = {"selection_scale": scale, "validation_zero_mse": zero, "grid": results, "chosen": chosen,
              "validation_families": len({r.group for r in validation}),
              "implementation_sha256": sha256_file(Path(__file__))}
    (output / "selection.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(chosen, indent=2))
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("select", "confirm"):
        sub = commands.add_parser(name)
        sub.add_argument("--output", type=Path, required=True)
        if name == "confirm":
            sub.add_argument("--selection", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "select":
        cmd_select(args.output)
        return 0
    from confirm import run_confirm  # noqa: PLC0415  (frozen protocol lives in confirm.py)
    run_confirm(args.output, args.selection)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

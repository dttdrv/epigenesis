"""Motif-plus-temporal-contrast sequence predictor for GSE188264 neural induction.

Reuses the family splits, asinh transform, family weighting, and bootstrap
evaluation from regulation.py alongside bundled JASPAR 2024 PWMs and vectorized
feature extractors from studio/epigenesis_studio/model/features.py.

Selected purely on the train -> validation split (no test labels used during
tuning). Outputs are retrospective held-out-family reporter activity differences
at [0, 3, 6, 12, 24, 48, 72]h, not endogenous gene expression or cell fates.
"""

from __future__ import annotations

import argparse
import base64
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import random
import statistics
import struct
import sys

import numpy as np

# coding-minimalism-rung: reuse (reuses repo brainc, regulation.py, and studio features/JASPAR assets)
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(REPO / "studio/epigenesis_studio/model"))

from brainc._io import load_json_object, read_regular_file
from brainc.development_bundle import compile_development
from brainc.source import FASTA_PROFILE, compile_source
from brainc.v2 import inline_storage, make_development_request, pack, policy_artifact, save, target_artifact
from brainc.v2._common import seal
from brainc.v2.target import DEV_DOMAIN, OP_UNIT_CREATE
from brainc.validator_development import load_development_bundle_directory, validate_development_paths
import features
import regulation as reg


MODEL_FORMAT = "epigenesis.neural-regulation.motif-temporal-ridge/v1"
TARGET = "org.epigenesis.example.neural-regulation.motif-temporal"
DEFAULT_MOTIF_SCALE = 1.0
DEFAULT_WORD_PENALTY = 0.3
DEFAULT_MOTIF_PENALTY = 3.0
DEFAULT_TEMPORAL_MULTIPLIER = 10.0


def _load_jaspar() -> tuple[dict, list[np.ndarray], str]:
    model_dir = REPO / "studio/epigenesis_studio/model"
    manifest_bytes = read_regular_file(model_dir / "model.json", maximum_bytes=1024 * 1024)
    manifest = load_json_object(manifest_bytes, "bundled JASPAR manifest")
    for name, expected in manifest["sha256"].items():
        if Path(name).name != name or hashlib.sha256((model_dir / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"bundled model asset digest mismatch: {name}")
    with np.load(model_dir / "weights.npz", allow_pickle=False) as saved:
        matrices = [saved[f"motif_{i}"] for i in range(len(manifest["motif_ids"]))]
    return manifest, matrices, hashlib.sha256(manifest_bytes).hexdigest()


def _temporal_basis(n_times: int = len(reg.TIMES)) -> np.ndarray:
    q = np.linalg.qr(np.column_stack([np.ones(n_times), np.eye(n_times)[:, 1:]]))[0]
    if q[:, 0].sum() < 0:
        q[:, 0] = -q[:, 0]
    return q


def _family_weights(rows: list[reg.Row]) -> np.ndarray:
    counts = Counter(row.group for row in rows)
    if not counts:
        raise ValueError("training rows are empty")
    n_groups = len(counts)
    return np.array([1.0 / (n_groups * counts[row.group]) for row in rows], dtype=np.float64)


def _solve_weighted_ridge(x: np.ndarray, y: np.ndarray, w: np.ndarray, penalties: np.ndarray) -> np.ndarray:
    xw = x * w[:, None]
    gram = xw.T @ x
    gram[np.diag_indices_from(gram)] += penalties
    return np.linalg.solve(gram, xw.T @ y)


def extract_features(rows: list[reg.Row] | list[dict], matrices: list[np.ndarray], motif_scale: float) -> tuple[np.ndarray, np.ndarray]:
    if not rows:
        raise ValueError("feature extraction requires at least one sequence pair")
    refs = [row["reference"] if isinstance(row, dict) else row.reference for row in rows]
    alts = [row["sequence"] if isinstance(row, dict) else row.sequence for row in rows]
    pairs = [{"reference": ref, "sequence": alt} for ref, alt in zip(refs, alts)]
    k5 = features.whole_counts(alts, 5) - features.whole_counts(refs, 5)
    md = features.deltas(pairs, matrices, motif_scale)
    return k5, md


def fit(
    rows: list[reg.Row],
    scale: float,
    *,
    motif_scale: float = DEFAULT_MOTIF_SCALE,
    word_penalty: float = DEFAULT_WORD_PENALTY,
    motif_penalty: float = DEFAULT_MOTIF_PENALTY,
    temporal_multiplier: float = DEFAULT_TEMPORAL_MULTIPLIER,
    shuffle_seed: int | None = None,
) -> tuple[dict, np.ndarray]:
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError("training scale must be finite and positive")
    for name, val in (("motif_scale", motif_scale), ("word_penalty", word_penalty),
                      ("motif_penalty", motif_penalty), ("temporal_multiplier", temporal_multiplier)):
        if not math.isfinite(val) or val <= 0:
            raise ValueError(f"{name} must be finite and positive")

    manifest, matrices, manifest_sha = _load_jaspar()
    k5, md = extract_features(rows, matrices, motif_scale)
    x = np.hstack([k5, md])
    y = np.array([reg.target(row, scale) for row in rows], dtype=np.float64)
    if shuffle_seed is not None:
        perm = list(range(len(y)))
        random.Random(shuffle_seed).shuffle(perm)
        y = y[perm]

    w = _family_weights(rows)
    q = _temporal_basis(len(reg.TIMES))
    y_q = y @ q

    base_penalties = np.concatenate([
        np.full(k5.shape[1], word_penalty, dtype=np.float64),
        np.full(md.shape[1], motif_penalty, dtype=np.float64),
    ])
    coef_mean = _solve_weighted_ridge(x, y_q[:, :1], w, base_penalties)
    coef_dev = _solve_weighted_ridge(x, y_q[:, 1:], w, base_penalties * temporal_multiplier)
    coefficients = np.hstack([coef_mean, coef_dev]) @ q.T
    if not np.isfinite(coefficients).all():
        raise ValueError("fitted coefficients contain nonfinite values")

    model = {
        "format": MODEL_FORMAT,
        "k": 5,
        "times_hours": reg.TIMES,
        "scale": float(scale),
        "motif_scale": float(motif_scale),
        "word_penalty": float(word_penalty),
        "motif_penalty": float(motif_penalty),
        "temporal_multiplier": float(temporal_multiplier),
        "motif_ids": manifest["motif_ids"],
        "motif_manifest_sha256": manifest_sha,
        "training_rows": len(rows),
        "training_groups": sorted({row.group for row in rows}),
        "provenance": {},
    }
    return model, coefficients


def predict_rows(model: dict, coefficients: np.ndarray, rows: list[reg.Row] | list[dict], matrices: list[np.ndarray] | None = None) -> np.ndarray:
    if model.get("format") != MODEL_FORMAT or model.get("times_hours") != reg.TIMES:
        raise ValueError("unsupported motif-temporal model format")
    if matrices is None:
        manifest, matrices, manifest_sha = _load_jaspar()
        if manifest_sha != model.get("motif_manifest_sha256") or manifest["motif_ids"] != model.get("motif_ids"):
            raise ValueError("motif manifest mismatch")
    expected_rows = len(features.vocabulary(5)) + len(model["motif_ids"])
    if coefficients.shape != (expected_rows, len(reg.TIMES)) or not np.isfinite(coefficients).all():
        raise ValueError("invalid coefficient matrix shape or values")
    k5, md = extract_features(rows, matrices, model["motif_scale"])
    pred = np.hstack([k5, md]) @ coefficients
    if not np.isfinite(pred).all():
        raise ValueError("predictions must be finite")
    return pred


def explain_rows(
    model: dict,
    coefficients: np.ndarray,
    rows: list[reg.Row] | list[dict],
    matrices: list[np.ndarray] | None = None,
    *,
    top_n: int = 5,
) -> list[list[dict]]:
    """Return top JASPAR motif contributors (by mean temporal effect magnitude) per edit row."""
    if model.get("format") != MODEL_FORMAT or model.get("times_hours") != reg.TIMES:
        raise ValueError("unsupported motif-temporal model format")
    if matrices is None:
        manifest, matrices, manifest_sha = _load_jaspar()
        if manifest_sha != model.get("motif_manifest_sha256") or manifest["motif_ids"] != model.get("motif_ids"):
            raise ValueError("motif manifest mismatch")
    n_k5 = len(features.vocabulary(5))
    m_coefs = coefficients[n_k5:, :]
    _, md = extract_features(rows, matrices, model["motif_scale"])
    motif_ids = model["motif_ids"]
    mean_c = m_coefs.mean(axis=1)
    h72_c = m_coefs[:, -1]
    out = []
    for idx in range(len(rows)):
        contribs = md[idx] * mean_c
        nz = [j for j in range(len(motif_ids)) if contribs[j] != 0.0]
        nz.sort(key=lambda j: (-abs(float(contribs[j])), motif_ids[j]))
        out.append([
            {
                "id": motif_ids[j],
                "delta": float(md[idx, j]),
                "effect": float(contribs[j]),
                "h72_effect": float(md[idx, j] * h72_c[j]),
            }
            for j in nz[:top_n]
        ])
    return out


def _fit_kmer_baseline(rows: list[reg.Row], k: int, penalties: tuple[float, ...], scale: float,
                       val_rows: list[reg.Row], sel_rows: list[reg.Row], sel_scale: float) -> tuple[dict, float]:
    sufficient = reg.moments(sel_rows, k, sel_scale)
    choices = []
    for penalty in penalties:
        candidate = reg.fit(sel_rows, k, penalty, sel_scale, sufficient=sufficient)
        preds = {row.id: reg.predict(candidate, row.reference, row.sequence) for row in val_rows}
        mse = statistics.mean(reg.family_losses(val_rows, preds, sel_scale).values())
        choices.append((mse, -penalty, penalty))
    best_penalty = min(choices)[2]
    return reg.fit(rows, k, best_penalty, scale), best_penalty


def _target(n_times: int = len(reg.TIMES)) -> dict:
    return target_artifact({
        "id": TARGET, "abi_major": 1, "opsets": [{"domain": DEV_DOMAIN, "version": 1}],
        "unit_schemas": [{"id": "reporter-effect", "fields": [{"id": "effect", "type": {"dtype": "f64", "shape": [n_times]},
          "unit": None, "mutability": "constant", "numeric": {"kind": "float"}}]}],
        "edge_schemas": [], "rules": [], "ports": [],
    })


def _clean_dna_file(path: Path, min_bases: int, max_bases: int = 16384) -> str:
    raw = read_regular_file(path, maximum_bytes=max_bases * 2 + 1024).decode("ascii")
    lines = [line.strip() for line in raw.splitlines() if not line.startswith(">")]
    seq = "".join("".join(lines).split()).upper()
    if not min_bases <= len(seq) <= max_bases or set(seq) - set("ACGT"):
        raise ValueError(f"prediction requires {min_bases}..{max_bases} unambiguous ACGT bases")
    return seq


def compile_prediction(model_dir: Path, reference_path: Path, variant_path: Path, output: Path) -> list[float]:
    model_bytes = read_regular_file(model_dir / "model.json", maximum_bytes=2 * 1024 * 1024)
    model = load_json_object(model_bytes, "motif-temporal model")
    coef_bytes = read_regular_file(model_dir / "coefficients.npy", maximum_bytes=2 * 1024 * 1024)
    if hashlib.sha256(coef_bytes).hexdigest() != model.get("coefficients_sha256"):
        raise ValueError("coefficients.npy identity mismatch")
    coefficients = np.load(model_dir / "coefficients.npy", allow_pickle=False)
    min_bases = int(model.get("k", 5))
    reference = _clean_dna_file(reference_path, min_bases)
    variant = _clean_dna_file(variant_path, min_bases)
    effect = [float(x) for x in predict_rows(model, coefficients, [{"reference": reference, "sequence": variant}])[0]]

    output.mkdir()
    (output / "model.json").write_bytes(model_bytes)
    (output / "coefficients.npy").write_bytes(coef_bytes)
    pair = output / "pair.fasta"
    pair.write_text(f">reference\n{reference}\n>variant\n{variant}\n")
    source = compile_source(FASTA_PROFILE, {"sequence": pair}, parameters={"wrapper": "identity"})
    source.save(output / "source")
    tensors = [
        {"id": name, "type": {"dtype": dtype, "shape": shape}, "unit": None,
         "axes": [None] * len(shape), "storage": inline_storage(pack(dtype, values))}
        for name, dtype, shape, values in (("count", "u64", [], [1]), ("effect", "f64", [len(effect)], effect))
    ]
    identity = {"kind": "content-sha256", "value": hashlib.sha256(model_bytes + coef_bytes).hexdigest()}
    manifest = seal({
        "format": "brainc.provider-manifest", "version": 2,
        "provider": {"name": TARGET, "version": "1"}, "model_identity": identity,
        "accepts": [f"brainc.source-descriptor/v1;profile={FASTA_PROFILE}"],
        "outputs": [{k: t[k] for k in ("id", "type", "unit", "axes")} for t in tensors],
    })
    save(manifest, output / "manifest.json")
    request = make_development_request(source, output / "manifest.json", [t["id"] for t in tensors])
    save(request, output / "request.json")
    save(seal({
        "format": "brainc.prediction-response", "version": 2,
        "request_artifact_sha256": request["artifact_sha256"],
        "provider": manifest["provider"], "model_identity": identity, "outputs": tensors,
    }), output / "response.json")
    target_contract = _target(len(effect))
    save(target_contract, output / "target.json")
    policy = policy_artifact(
        TARGET, {"id": TARGET, "abi_major": 1, "contract_sha256": target_contract["contract_sha256"]},
        [{"id": t["id"], "from_output": t["id"]} for t in tensors],
        [{"id": "create.effect", "op": OP_UNIT_CREATE, "version": 1, "schema": "reporter-effect",
          "count": "count", "initializers": [{"field": "effect", "tensor": "effect"}]}],
    )
    save(policy, output / "policy.json")
    compile_development(source, output / "manifest.json", output / "request.json",
                        output / "response.json", output / "policy.json", output / "target.json").save(output / "development")
    validation = validate_development_paths(output / "source", output / "development", source_inputs={"sequence": pair})
    if not validation["valid"]:
        raise ValueError(f"development validation failed: {validation}")
    (output / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    return load_prediction(output)


def load_prediction(directory: Path) -> list[float]:
    validation = validate_development_paths(directory / "source", directory / "development",
                                            source_inputs={"sequence": directory / "pair.fasta"})
    if not validation["valid"]:
        raise ValueError(f"prediction artifact validation failed: {validation}")
    _, artifacts = load_development_bundle_directory(directory / "development")
    model_bytes = read_regular_file(directory / "model.json", maximum_bytes=2 * 1024 * 1024)
    model = load_json_object(model_bytes, "motif-temporal model")
    n_times = len(model.get("times_hours", reg.TIMES))
    coef_bytes = read_regular_file(directory / "coefficients.npy", maximum_bytes=2 * 1024 * 1024)
    expected_id = hashlib.sha256(model_bytes + coef_bytes).hexdigest()
    if artifacts["provider_manifest"]["model_identity"]["value"] != expected_id:
        raise ValueError("prediction model identity mismatch")
    if artifacts["target_contract"] != _target(n_times):
        raise ValueError("unsupported prediction target")
    lines = read_regular_file(directory / "pair.fasta", maximum_bytes=65536).decode("ascii").splitlines()
    if len(lines) != 4 or lines[0] != ">reference" or lines[2] != ">variant":
        raise ValueError("unsupported prediction source layout")
    module = artifacts["development_module"]["module"]
    tensors = {t["id"]: base64.b64decode(t["storage"]["data"], validate=True) for t in module["tensors"]}
    if struct.unpack("<Q", tensors["count"]) != (1,):
        raise ValueError("expected one perturbation effect")
    effect = list(struct.unpack(f"<{n_times}d", tensors["effect"]))
    if hashlib.sha256(coef_bytes).hexdigest() != model.get("coefficients_sha256"):
        raise ValueError("coefficients.npy identity mismatch")
    coefficients = np.load(directory / "coefficients.npy", allow_pickle=False)
    recomputed = [float(x) for x in predict_rows(model, coefficients, [{"reference": lines[1], "sequence": lines[3]}])[0]]
    if any(not math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12) for a, b in zip(effect, recomputed)):
        raise ValueError("compiled effect disagrees with recomputed prediction")
    return effect


def benchmark(output: Path, data: Path | None = None) -> dict:
    output.mkdir()
    rows, qc = reg.load_data(data)
    training = [row for row in rows if row.split == "train"]
    validation = [row for row in rows if row.split == "validation" and row.kind == "targeted"]
    fitting = [row for row in rows if row.split != "test"]
    testing = [row for row in rows if row.split == "test"]
    primary = [row for row in testing if row.kind == "targeted"]

    selection_scale = reg.training_scale(training)
    scale = reg.training_scale(fitting)
    manifest, matrices, _ = _load_jaspar()

    # Fixed hyperparameters selected on validation only (prior to test evaluation).
    model, coefficients = fit(fitting, scale)
    np.save(output / "coefficients.npy", coefficients, allow_pickle=False)
    coef_sha = hashlib.sha256((output / "coefficients.npy").read_bytes()).hexdigest()
    model["coefficients_sha256"] = coef_sha
    model["provenance"] = {
        "sources": qc["source_provenance"],
        "selection_scale": selection_scale,
        "implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (output / "model.json").write_text(json.dumps(model, indent=2, allow_nan=False) + "\n")

    k2_model, k2_penalty = _fit_kmer_baseline(fitting, 2, reg.PENALTIES, scale, validation, training, selection_scale)
    k5_model, k5_penalty = _fit_kmer_baseline(fitting, 5, reg.PENALTIES, scale, validation, training, selection_scale)

    pred_mat = predict_rows(model, coefficients, testing, matrices)
    predicted = {row.id: [float(x) for x in pred_mat[i]] for i, row in enumerate(testing)}
    k2_preds = {row.id: reg.predict(k2_model, row.reference, row.sequence) for row in testing}
    k5_preds = {row.id: reg.predict(k5_model, row.reference, row.sequence) for row in testing}
    zero = {row.id: [0.0] * len(reg.TIMES) for row in testing}
    averaged = {key: [statistics.mean(vals)] * len(reg.TIMES) for key, vals in predicted.items()}

    scores = reg.family_losses(primary, predicted, scale)
    intervals = {
        name: reg.improvement_interval(reg.family_losses(primary, other, scale), scores)
        for name, other in (("zero_effect", zero), ("dinucleotide", k2_preds), ("k5_baseline", k5_preds))
    }
    temporal_ci = reg.improvement_interval(reg.family_losses(primary, averaged, scale), scores)

    permutations = []
    for seed in (19, 73, 211):
        shuf_model, shuf_coef = fit(fitting, scale, shuffle_seed=seed)
        shuf_mat = predict_rows(shuf_model, shuf_coef, primary, matrices)
        shuf_preds = {row.id: [float(x) for x in shuf_mat[i]] for i, row in enumerate(primary)}
        permutations.append({"seed": seed, "mse": statistics.mean(reg.family_losses(primary, shuf_preds, scale).values())})

    # Inspect largest learned motif coefficients (mean & temporal 72h-0h slope)
    n_k5 = len(features.vocabulary(5))
    motif_coef = coefficients[n_k5:, :]
    motif_mean = motif_coef.mean(axis=1)
    motif_slope = motif_coef[:, -1] - motif_coef[:, 0]
    top_motifs = [
        {"motif_id": manifest["motif_ids"][idx],
         "mean_effect": round(float(motif_mean[idx]), 6),
         "late_minus_early": round(float(motif_slope[idx]), 6)}
        for idx in np.argsort(np.linalg.norm(motif_coef, axis=1))[::-1][:10]
    ]

    with (output / "test-predictions.csv").open("x", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("id", "locus", "group", "kind", "hour", "observed_effect",
                         "motif_temporal_effect", "k5_effect", "dinucleotide_effect"))
        for row in testing:
            observed = reg.target(row, scale)
            for h, hour in enumerate(reg.TIMES):
                writer.writerow((row.id, row.locus, row.group, row.kind, hour, observed[h],
                                 predicted[row.id][h], k5_preds[row.id][h], k2_preds[row.id][h]))

    # Verify brainc compilation on the first targeted held-out variant
    sample = primary[0]
    ref_path = output / "sample_ref.dna"
    var_path = output / "sample_var.dna"
    ref_path.write_text(sample.reference)
    var_path.write_text(sample.sequence)
    compiled_sample = compile_prediction(output, ref_path, var_path, output / "sample_prediction")

    report = {
        "status": "evaluation-complete",
        "model_format": MODEL_FORMAT,
        "scope": "retrospective held-out-family reporter effects; hyperparameters frozen on validation after October 2026 k=5 release",
        "predictive_acceptance_vs_zero_and_dinucleotide": all(intervals[k][0] > 0 for k in ("zero_effect", "dinucleotide")),
        "predictive_acceptance_vs_k5": intervals["k5_baseline"][0] > 0,
        "biological_development_accepted": False,
        "hyperparameters": {
            "motif_scale": DEFAULT_MOTIF_SCALE,
            "word_penalty": DEFAULT_WORD_PENALTY,
            "motif_penalty": DEFAULT_MOTIF_PENALTY,
            "temporal_multiplier": DEFAULT_TEMPORAL_MULTIPLIER,
            "k2_penalty": k2_penalty,
            "k5_penalty": k5_penalty,
        },
        "test_groups": len(scores),
        "test_targeted_rows": len(primary),
        "test_targeted_mse": {
            "motif_temporal": statistics.mean(scores.values()),
            "k5_baseline": statistics.mean(reg.family_losses(primary, k5_preds, scale).values()),
            "dinucleotide": statistics.mean(reg.family_losses(primary, k2_preds, scale).values()),
            "zero_effect": statistics.mean(reg.family_losses(primary, zero, scale).values()),
        },
        "improvement_intervals_95": intervals,
        "temporal_improvement_interval_95": temporal_ci,
        "temporal_contrast_mse": {
            name: statistics.mean(reg.family_losses(primary, vals, scale, temporal=True).values())
            for name, vals in (("motif_temporal", predicted), ("k5_baseline", k5_preds),
                               ("dinucleotide", k2_preds), ("time_averaged", averaged))
        },
        "label_permutation_diagnostics": permutations,
        "categories": {
            kind: {
                name: statistics.mean(reg.family_losses([row for row in testing if row.kind == kind], vals, scale).values())
                for name, vals in (("motif_temporal", predicted), ("k5_baseline", k5_preds),
                                   ("dinucleotide", k2_preds), ("zero_effect", zero))
            }
            for kind in ("targeted", "random", "scrambled")
        },
        "top_motifs": top_motifs,
        "brainc_sample_replay_match": max(abs(a - b) for a, b in zip(compiled_sample, predicted[sample.id])) < 1e-12,
    }
    (output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    trial = commands.add_parser("benchmark")
    trial.add_argument("--output", type=Path, required=True)
    trial.add_argument("--data", type=Path)
    pred = commands.add_parser("predict")
    for name in ("model-dir", "reference", "variant", "output"):
        pred.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == "benchmark":
            print(json.dumps(benchmark(args.output, args.data), indent=2))
        else:
            compiled = compile_prediction(args.model_dir, args.reference, args.variant, args.output)
            print(json.dumps({"times_hours": reg.TIMES, "predicted_effect": compiled}))
    except (ValueError, OSError, OverflowError) as failure:
        print(f"motif-temporal experiment failed: {failure}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

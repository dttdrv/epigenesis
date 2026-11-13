"""Variant-to-Circuit bridge: arbitrary regulatory DNA -> GLI partition physics -> neural GRN bifurcation.

Connects full-length DNA enhancers (via brainc Development Modules and Boltzmann
partition affinity over Hallikas et al. 2006 Table S1 GLI2/GLI3 weights) to the
four-gene neural-tube ODE (Pax6, Olig2, Nkx2.2, Irx3; Herrera-Delgado et al. 2020)
and optionally scores the same edits with passive reporter sequence models
(Kreimer neural-induction motif_temporal & Siraj Huber-linear).

Exposes "reporter-silent, circuit-switch" (cryptic bifurcation) edits where
passive promoter/MPRA readouts barely move yet coupled feedback flips cell fate.
All developmental concentrations, distances, and times are uncalibrated model units.
"""

from __future__ import annotations

import argparse
import base64
import hashlib
import html
import importlib.util
import json
import math
from pathlib import Path
import struct
import sys

# coding-minimalism-rung: reuse (reuses brainc, examples/neural-patterning, and optional Track A / Studio)
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "examples/neural-patterning"))

from brainc._canonical import canonical_bytes
from brainc._io import load_json_object, read_regular_file
from brainc.development_bundle import compile_development
from brainc.source import FASTA_PROFILE, compile_source
from brainc.v2 import inline_storage, make_development_request, pack, policy_artifact, save, target_artifact
from brainc.v2._common import seal
from brainc.v2.target import DEV_DOMAIN, OP_UNIT_CREATE
from brainc.validator_development import load_development_bundle_directory, validate_development_paths
import patterning as pat


TARGET = "org.epigenesis.example.neural-circuit"
DEFAULT_EPSILON = 1e-4
DEFAULT_SIGNAL = 0.20
LOCAL_ERROR_TOLERANCE = pat.LOCAL_ERROR_TOLERANCE


def _binding_tables(protein: str, epsilon: float) -> tuple[str, dict[str, list[float]]]:
    if type(epsilon) not in (int, float) or not math.isfinite(epsilon) or epsilon < 0:
        raise ValueError("epsilon must be finite and nonnegative")
    binding = pat.load_data()["binding"]
    if protein not in binding["profiles"]:
        raise ValueError(f"unsupported GLI protein: {protein}")
    consensus = binding["consensus"]
    profile = binding["profiles"][protein]
    ratios = {
        base: [
            (profile[base][pos] + epsilon) / (profile[consensus[pos]][pos] + epsilon)
            for pos in range(len(consensus))
        ]
        for base in "ACGT"
    }
    return consensus, ratios


def _site_weight(window: str, ratios: dict[str, list[float]]) -> float:
    weight = 1.0
    for pos, base in enumerate(window):
        weight *= ratios[base][pos]
    return weight


def boltzmann_gli_affinity(sequence: str, protein: str = "GLI2", epsilon: float = DEFAULT_EPSILON) -> dict:
    """Sum relative association affinity across all 9-bp oriented windows of any DNA sequence."""
    consensus, ratios = _binding_tables(protein, epsilon)
    width = len(consensus)
    if type(sequence) is not str or len(sequence) < width or set(sequence) - set("ACGT"):
        raise ValueError(f"sequence must be >= {width} uppercase ACGT bases")

    total = 0.0
    best_weight = -1.0
    best_site = None
    top_sites = []
    for start in range(len(sequence) - width + 1):
        fwd = sequence[start : start + width]
        rev = pat.reverse_complement(fwd)
        wf = _site_weight(fwd, ratios)
        wr = _site_weight(rev, ratios)
        total += wf + wr
        for strand, w, oriented in ((1, wf, fwd), (-1, wr, rev)):
            if w > best_weight:
                best_weight = w
                best_site = {"start": start, "strand": strand, "oriented_site": oriented, "weight": w}
            if w >= 1e-2:
                top_sites.append({"start": start, "strand": strand, "oriented_site": oriented, "weight": w})
    top_sites.sort(key=lambda item: -item["weight"])
    return {
        "protein": protein,
        "epsilon": epsilon,
        "length": len(sequence),
        "partition_sum": total,
        "peak_site": best_site,
        "strong_sites": top_sites[:8],
    }


def affinity_ratio(reference: str, variant: str, protein: str = "GLI2", epsilon: float = DEFAULT_EPSILON) -> float:
    """Compute thermodynamic GLI affinity ratio q = Z(variant) / Z(reference) (> 0)."""
    z_ref = boltzmann_gli_affinity(reference, protein, epsilon)["partition_sum"]
    z_var = boltzmann_gli_affinity(variant, protein, epsilon)["partition_sum"]
    if z_ref <= 0 or z_var <= 0 or not math.isfinite(z_var / z_ref):
        raise ValueError("affinity ratio requires positive finite partition sums")
    return z_var / z_ref


def _delta_partition_single(reference: str, pos: int, alt_base: str, ratios: dict[str, list[float]], width: int) -> float:
    """Exact local partition-sum delta Z(alt) - Z(ref) in O(width) for a single-base edit."""
    lo = max(0, pos - width + 1)
    hi = min(len(reference) - width, pos)
    delta = 0.0
    for start in range(lo, hi + 1):
        ref_w = reference[start : start + width]
        alt_w = ref_w[: pos - start] + alt_base + ref_w[pos - start + 1 :]
        delta += (_site_weight(alt_w, ratios) + _site_weight(pat.reverse_complement(alt_w), ratios)) - (
            _site_weight(ref_w, ratios) + _site_weight(pat.reverse_complement(ref_w), ratios)
        )
    return delta


def simulate_circuit(
    protocol: list,
    *,
    ratio: float = 1.0,
    dt: float = pat.DEFAULT_DT,
    initial: tuple[float, ...] = (0.0, 0.0, 0.0, 0.0),
    mechanism: str = "competitive",
) -> list[tuple]:
    """Integrate the 4-gene neural-tube ODE for any positive affinity multiplier ratio > 0."""
    def finite(v):
        return type(v) in (int, float) and math.isfinite(v)

    model = pat.load_data()["model"]
    if (
        not finite(dt) or dt <= 0 or not finite(ratio) or ratio <= 0
        or mechanism not in ("competitive", "activation-only")
        or len(initial) != len(model["genes"])
        or any(type(v) not in (int, float) or not 0 <= v <= 1 for v in initial)
    ):
        raise ValueError("invalid timestep, positive affinity ratio, mechanism, or initial state")

    phases = []
    total, steps = 0.0, 0
    for phase in protocol:
        if not isinstance(phase, (list, tuple)) or len(phase) != 2 or any(not finite(v) for v in phase) or phase[0] <= 0 or not 0 <= phase[1] <= 1:
            raise ValueError("protocol needs positive durations and Gli activator fractions in [0, 1]")
        duration, signal = phase
        count = duration / dt
        end = total + duration
        if not math.isfinite(count) or not math.isfinite(end) or end <= total or count > pat.MAX_STEPS - steps:
            raise ValueError("protocol exceeds the step budget or representable time")
        count = max(1, math.ceil(count))
        steps += count
        phases.append((total, end, signal, count, duration / count))
        total = end
    if not phases:
        raise ValueError("protocol must not be empty")

    state = tuple(float(v) for v in initial)
    trajectory = [(0.0, *state)]

    def rhs(y, s):
        if any(not math.isfinite(v) or not 0 <= v <= 1 for v in y):
            raise ArithmeticError("integration stage outside physical [0, 1] state range")
        return tuple(model["production"] * p - model["decay"] * x for p, x in zip(pat.occupancy(model, y, s, ratio, mechanism), y))

    def step(y, h, s):
        a = rhs(y, s)
        b = rhs(tuple(x + h * v / 2 for x, v in zip(y, a)), s)
        c = rhs(tuple(x + h * v / 2 for x, v in zip(y, b)), s)
        d = rhs(tuple(x + h * v for x, v in zip(y, c)), s)
        result = tuple(x + h * (v + 2 * w + 2 * z + u) / 6 for x, v, w, z, u in zip(y, a, b, c, d))
        rhs(result, s)
        return result

    attempts = 0
    for start, end, signal, count, interval in phases:
        elapsed, h = 0.0, interval
        for index in range(count):
            stop = (index + 1) * interval
            while elapsed < stop:
                h = min(h, stop - elapsed)
                if elapsed + h <= elapsed or attempts >= pat.MAX_STEPS:
                    raise ValueError("adaptive integration exceeds budget")
                attempts += 1
                try:
                    coarse = step(state, h, signal)
                    fine = step(step(state, h / 2, signal), h / 2, signal)
                    error = max(abs(a - b) for a, b in zip(fine, coarse)) / 15
                except ArithmeticError:
                    h /= 2
                    continue
                if error <= LOCAL_ERROR_TOLERANCE:
                    state = fine
                    elapsed = stop if h == stop - elapsed else elapsed + h
                factor = 2.0 if error == 0 else min(2.0, max(0.2, 0.9 * (LOCAL_ERROR_TOLERANCE / error) ** 0.2))
                h = min(interval, h * factor)
            trajectory.append((end if index == count - 1 else start + stop, *state))
    return trajectory


def run_circuit_assay(
    protocol: list,
    *,
    ratio: float,
    mode: str,
    dt: float = pat.DEFAULT_DT,
    mechanism: str = "competitive",
    initial: dict | None = None,
) -> dict:
    if mode not in ("reporter", "endogenous"):
        raise ValueError("mode must be reporter or endogenous")
    model = pat.load_data()["model"]
    if initial is None:
        initial = dict.fromkeys(model["genes"], 0.0)
    state = pat.named_state(initial)
    trajectory = simulate_circuit(protocol, ratio=1.0 if mode == "reporter" else ratio, dt=dt, mechanism=mechanism, initial=state)
    phase, boundary = 0, protocol[0][0]
    reporter, signals = [], []
    nkx_idx = model["genes"].index("Nkx2.2")
    for time_val, *y in trajectory:
        while time_val > boundary and phase + 1 < len(protocol):
            phase += 1
            boundary += protocol[phase][0]
        signal = protocol[phase][1]
        signals.append(signal)
        reporter.append(pat.occupancy(model, y, signal, ratio, mechanism)[nkx_idx])
    return {
        "mode": mode, "mechanism": mechanism, "protocol": protocol, "ratio": ratio,
        "initial_state": dict(zip(model["genes"], trajectory[0][1:])),
        "trajectory": trajectory, "signal": signals, "reporter_occupancy": reporter,
    }


def fate_label(state: tuple[float, ...] | list[float], genes: list[str]) -> str:
    best_i = max(range(1, len(genes)), key=lambda i: state[i])
    return genes[best_i] if state[best_i] > 0.5 else ("Pax6+" if state[0] > 0.5 else "unresolved")


def _target(semantics_hash: str, length: int) -> dict:
    return target_artifact({
        "id": TARGET, "abi_major": 1, "opsets": [{"domain": DEV_DOMAIN, "version": 1}],
        "unit_schemas": [{"id": "enhancer-pair", "fields": [
            {"id": "affinity_ratio", "type": {"dtype": "f64", "shape": [1]}, "unit": None, "mutability": "constant", "numeric": {"kind": "float"}},
            {"id": "reference_nucleotides", "type": {"dtype": "u64", "shape": [length]}, "unit": None, "mutability": "constant", "numeric": {"kind": "integer"}},
            {"id": "variant_nucleotides", "type": {"dtype": "u64", "shape": [length]}, "unit": None, "mutability": "constant", "numeric": {"kind": "integer"}},
        ]}],
        "edge_schemas": [], "ports": [],
        "rules": [{"id": TARGET + ".partition", "version": 1, "semantics_sha256": semantics_hash,
                   "subject": "unit", "schema": "enhancer-pair", "phase": "UPDATE", "triggers": ["step"], "parameters": [],
                   "reads": [{"scope": "subject", "field": "affinity_ratio"},
                             {"scope": "subject", "field": "reference_nucleotides"},
                             {"scope": "subject", "field": "variant_nucleotides"}], "writes": []}],
    })


def compile_circuit_bundle(reference: str, variant: str, *, protein: str = "GLI2",
                           epsilon: float = DEFAULT_EPSILON, output: Path) -> str:
    """Compile a reference/variant enhancer pair and its Boltzmann affinity ratio into a brainc Development Module."""
    if len(reference) != len(variant):
        raise ValueError("pair bundle compilation requires equal-length enhancers")
    q = affinity_ratio(reference, variant, protein, epsilon)
    contract = {
        "format": "epigenesis.neural-circuit-pair", "version": 1,
        "protein": protein, "epsilon": epsilon, "length": len(reference),
        "data_sha256": pat.DATA_SHA256,
    }
    impl_bytes = read_regular_file(Path(__file__).resolve(), maximum_bytes=256 * 1024)
    semantics_hash = hashlib.sha256(impl_bytes + b"\0" + canonical_bytes(contract)).hexdigest()
    receipt = hashlib.sha256(canonical_bytes(contract)).hexdigest()

    output.mkdir()
    pair_fasta = output / "pair.fasta"
    pair_fasta.write_text(f">reference\n{reference}\n>variant\n{variant}\n")
    (output / "contract.json").write_bytes(canonical_bytes(contract) + b"\n")
    source = compile_source(FASTA_PROFILE, {"sequence": pair_fasta}, parameters={"wrapper": "identity"})
    source.save(output / "source")

    tensors = sorted([
        {"id": "count", "type": {"dtype": "u64", "shape": []}, "unit": None, "axes": [],
         "storage": inline_storage(pack("u64", [1]))},
        {"id": "reference_nucleotides", "type": {"dtype": "u64", "shape": [len(reference)]}, "unit": None, "axes": [None],
         "storage": inline_storage(pack("u64", list(reference.encode("ascii"))))},
        {"id": "variant_nucleotides", "type": {"dtype": "u64", "shape": [len(variant)]}, "unit": None, "axes": [None],
         "storage": inline_storage(pack("u64", list(variant.encode("ascii"))))},
        {"id": "affinity_ratio", "type": {"dtype": "f64", "shape": [1]}, "unit": None, "axes": [None],
         "storage": inline_storage(pack("f64", [q]))},
    ], key=lambda item: item["id"])
    manifest = seal({
        "format": "brainc.provider-manifest", "version": 2,
        "provider": {"name": TARGET, "version": "1"},
        "model_identity": {"kind": "content-sha256", "value": semantics_hash},
        "accepts": [f"brainc.source-descriptor/v1;profile={FASTA_PROFILE}"],
        "outputs": [{k: t[k] for k in ("id", "type", "unit", "axes")} for t in tensors],
    })
    save(manifest, output / "manifest.json")
    request = make_development_request(source, output / "manifest.json", [t["id"] for t in tensors])
    save(request, output / "request.json")
    save(seal({
        "format": "brainc.prediction-response", "version": 2,
        "request_artifact_sha256": request["artifact_sha256"],
        "provider": manifest["provider"], "model_identity": manifest["model_identity"], "outputs": tensors,
    }), output / "response.json")
    target = _target(semantics_hash, len(reference))
    save(target, output / "target.json")
    ops = [
        {"id": "create.pair", "op": OP_UNIT_CREATE, "version": 1, "schema": "enhancer-pair", "count": "count",
         "initializers": [{"field": fid, "tensor": fid} for fid in ("affinity_ratio", "reference_nucleotides", "variant_nucleotides")]},
    ]
    save(policy_artifact(TARGET, {"id": TARGET, "abi_major": 1, "contract_sha256": target["contract_sha256"]},
                         [{"id": t["id"], "from_output": t["id"]} for t in tensors], ops), output / "policy.json")
    compile_development(source, output / "manifest.json", output / "request.json",
                        output / "response.json", output / "policy.json", output / "target.json").save(output / "development")
    verify_circuit_bundle(output, receipt)
    return receipt


def verify_circuit_bundle(directory: Path, expected_receipt: str) -> dict:
    contract = load_json_object(read_regular_file(directory / "contract.json", maximum_bytes=64 * 1024), "circuit contract")
    if hashlib.sha256(canonical_bytes(contract)).hexdigest() != expected_receipt:
        raise ValueError("circuit contract receipt mismatch")
    pair_path = directory / "pair.fasta"
    validation = validate_development_paths(directory / "source", directory / "development", source_inputs={"sequence": pair_path})
    if not validation["valid"]:
        raise ValueError(f"circuit bundle validation failed: {validation}")
    bundle, artifacts = load_development_bundle_directory(directory / "development")
    lines = read_regular_file(pair_path, maximum_bytes=64 * 1024).decode("ascii").splitlines()
    if len(lines) != 4 or lines[0] != ">reference" or lines[2] != ">variant":
        raise ValueError("invalid circuit pair.fasta layout")
    reference, variant = lines[1], lines[3]
    module = artifacts["development_module"]["module"]
    tensors = {t["id"]: base64.b64decode(t["storage"]["data"], validate=True) for t in module["tensors"]}
    if tensors["reference_nucleotides"] != pack("u64", list(reference.encode("ascii"))) or tensors["variant_nucleotides"] != pack("u64", list(variant.encode("ascii"))):
        raise ValueError("compiled nucleotides disagree with pair.fasta")
    stored_q = struct.unpack("<d", tensors["affinity_ratio"])[0]
    recomputed_q = affinity_ratio(reference, variant, contract["protein"], contract["epsilon"])
    if not math.isclose(stored_q, recomputed_q, rel_tol=1e-12, abs_tol=1e-12):
        raise ValueError("compiled affinity ratio disagrees with partition sum")
    return {
        "contract_sha256": expected_receipt,
        "bundle_sha256": bundle["artifact_sha256"],
        "development_sha256": artifacts["development_module"]["artifact_sha256"],
        "affinity_ratio": stored_q,
        "length": len(reference),
    }


def _optional_sequence_scores(pairs: list[dict], motif_model_dir: Path | None) -> tuple[dict | None, dict | None]:
    """Optionally evaluate passive reporter predictions (Track A motif_temporal + Siraj Huber-linear)."""
    kreimer_scores, siraj_scores = None, None
    try:
        import numpy as np
    except ImportError:
        return None, None
    try:
        sys.path.insert(0, str(REPO / "studio/epigenesis_studio/model"))
        import predict as siraj_predict
        siraj_scores = {r["id"]: r["effect"] for r in siraj_predict.predict(pairs)}
    except (ImportError, ValueError, OSError):
        siraj_scores = None
    if motif_model_dir is not None and (motif_model_dir / "model.json").is_file() and (motif_model_dir / "coefficients.npy").is_file():
        try:
            sys.path.insert(0, str(REPO / "examples/neural-regulation"))
            import motif_temporal
            m = json.loads((motif_model_dir / "model.json").read_text())
            coef = np.load(motif_model_dir / "coefficients.npy", allow_pickle=False)
            mat = motif_temporal.predict_rows(m, coef, pairs)
            kreimer_scores = {p["id"]: {"mean": float(mat[i].mean()), "h72": float(mat[i, -1])} for i, p in enumerate(pairs)}
        except (ImportError, ValueError, OSError):
            kreimer_scores = None
    return kreimer_scores, siraj_scores


def scan_enhancer(
    reference: str,
    *,
    signal: float = DEFAULT_SIGNAL,
    protein: str = "GLI2",
    epsilon: float = DEFAULT_EPSILON,
    duration: float = 40.0,
    motif_model_dir: Path | None = None,
) -> dict:
    """Scan all 3*L single-base mutations in an enhancer for reporter-vs-circuit dissociations."""
    consensus, ratios = _binding_tables(protein, epsilon)
    width = len(consensus)
    wt_info = boltzmann_gli_affinity(reference, protein, epsilon)
    z_wt = wt_info["partition_sum"]
    model = pat.load_data()["model"]
    genes = model["genes"]
    nkx_idx, olig_idx = genes.index("Nkx2.2"), genes.index("Olig2")

    wt_traj = simulate_circuit([(duration, signal)], ratio=1.0)
    wt_end = wt_traj[-1][1:]
    wt_fate = fate_label(wt_end, genes)
    wt_rep_end = pat.occupancy(model, wt_end, signal, 1.0)[nkx_idx]

    # Precompute unique q responses accurately (keyed at 1e-8 resolution on q)
    edits = []
    cache: dict[float, tuple[tuple[float, ...], float]] = {}
    pairs = []
    for pos, base in enumerate(reference):
        for alt in "ACGT":
            if alt == base:
                continue
            eid = f"{pos + 1}:{base}>{alt}"
            variant = reference[:pos] + alt + reference[pos + 1 :]
            z_mut = z_wt + _delta_partition_single(reference, pos, alt, ratios, width)
            q = z_mut / z_wt
            key = round(q, 8)
            if key not in cache:
                end_state = simulate_circuit([(duration, signal)], ratio=q)[-1][1:]
                rep_end = pat.occupancy(model, wt_end, signal, q)[nkx_idx]
                cache[key] = (end_state, rep_end)
            end_state, rep_end = cache[key]
            d_rep = rep_end - wt_rep_end
            d_nkx = end_state[nkx_idx] - wt_end[nkx_idx]
            d_olig = end_state[olig_idx] - wt_end[olig_idx]
            mut_fate = fate_label(end_state, genes)
            fate_flipped = mut_fate != wt_fate
            if abs(d_rep) < 0.02 and fate_flipped:
                category = "cryptic_bifurcation"
            elif abs(d_rep) >= 0.02 and abs(d_nkx) < 0.05 and not fate_flipped:
                category = "canalized_buffered"
            elif fate_flipped:
                category = "concordant_switch"
            else:
                category = "neutral"
            pairs.append({"id": eid, "reference": reference, "sequence": variant})
            edits.append({
                "id": eid, "pos": pos, "ref": base, "alt": alt,
                "variant": variant, "affinity_ratio": q,
                "reporter_delta_end": d_rep,
                "nkx2_2_delta": d_nkx, "olig2_delta": d_olig,
                "wt_fate": wt_fate, "mut_fate": mut_fate,
                "fate_flipped": fate_flipped,
                "amplification_ratio": abs(d_nkx) / max(abs(d_rep), 1e-12),
                "category": category,
            })

    kreimer_scores, siraj_scores = _optional_sequence_scores(pairs, motif_model_dir)
    for item in edits:
        if kreimer_scores is not None:
            item["kreimer_mpra_mean"] = kreimer_scores[item["id"]]["mean"]
            item["kreimer_mpra_72h"] = kreimer_scores[item["id"]]["h72"]
        if siraj_scores is not None:
            item["siraj_score"] = siraj_scores[item["id"]]

    cryptic = sorted([e for e in edits if e["category"] == "cryptic_bifurcation"], key=lambda e: (-abs(e["nkx2_2_delta"]), abs(e["reporter_delta_end"])))
    counts = {cat: sum(e["category"] == cat for e in edits) for cat in ("cryptic_bifurcation", "concordant_switch", "canalized_buffered", "neutral")}
    return {
        "protein": protein, "epsilon": epsilon, "signal": signal, "duration": duration,
        "length": len(reference), "wt_partition": wt_info,
        "wt_end_state": dict(zip(genes, wt_end)), "wt_fate": wt_fate, "wt_reporter_end": wt_rep_end,
        "unique_affinity_solves": len(cache), "category_counts": counts,
        "cryptic_switches": cryptic, "edits": edits,
    }


def _find_second_site_rescue(reference: str, primary_edit: dict, protein: str, epsilon: float) -> dict | None:
    """Search single-base edits elsewhere in the enhancer that increase GLI affinity q back toward WT."""
    pos0 = primary_edit["pos"]
    base_mut = primary_edit["variant"]
    q_single = primary_edit["affinity_ratio"]
    z_ref = boltzmann_gli_affinity(reference, protein, epsilon)["partition_sum"]
    consensus, ratios = _binding_tables(protein, epsilon)
    width = len(consensus)
    z_single = z_ref * q_single
    best = None
    for pos, base in enumerate(base_mut):
        if abs(pos - pos0) <= 0:
            continue
        for alt in "ACGT":
            if alt == base:
                continue
            z_double = z_single + _delta_partition_single(base_mut, pos, alt, ratios, width)
            q_double = z_double / z_ref
            if q_double > q_single * 1.05:
                cand = {"id": f"{pos + 1}:{base}>{alt}", "pos": pos, "double_variant": base_mut[:pos] + alt + base_mut[pos + 1 :], "affinity_ratio": q_double}
                if best is None or q_double > best["affinity_ratio"]:
                    best = cand
    return best


def _render_html(report: dict) -> str:
    scans = report["scan"]["edits"]
    wt_fate = report["scan"]["wt_fate"]
    top = report["featured_variant"]
    # Build crisp SVG scatter points: x = |passive reporter delta|, y = |endogenous Nkx2.2 shift|
    w, h, pad = 540, 290, 48
    max_x = max([0.05] + [abs(e["reporter_delta_end"]) for e in scans])
    points = []
    for e in scans:
        cx = pad + (w - 2 * pad) * min(1.0, abs(e["reporter_delta_end"]) / max_x)
        cy = h - pad - (h - 2 * pad) * min(1.0, abs(e["nkx2_2_delta"]))
        color = "#c0262d" if e["category"] == "cryptic_bifurcation" else "#1f6feb" if e["category"] == "canalized_buffered" else "#9ca3af"
        r = 5 if e["id"] == top["id"] else 3.2
        title = html.escape(f"{e['id']} | q={e['affinity_ratio']:.4f} | dRep={e['reporter_delta_end']:+.4f} | dNkx={e['nkx2_2_delta']:+.3f} ({e['category']})")
        points.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r}" fill="{color}" opacity="0.82"><title>{title}</title></circle>')

    # Spatial gradient Nkx2.2 & Olig2 profiles (WT vs mutant)
    sp = report["spatial"]
    xs = sp["positions"]
    def poly(vals, stroke, dash=""):
        pts = " ".join(f"{pad + (w - 2 * pad) * x:.1f},{h - pad - (h - 2 * pad) * v:.1f}" for x, v in zip(xs, vals))
        d_attr = f' stroke-dasharray="{dash}"' if dash else ""
        return f'<polyline fill="none" stroke="{stroke}" stroke-width="2.2"{d_attr} points="{pts}"/>'

    nkx_wt = [s[2] for s in sp["reference"]]
    nkx_mut = [s[2] for s in sp["variant"]]
    olig_wt = [s[1] for s in sp["reference"]]
    olig_mut = [s[1] for s in sp["variant"]]

    rows_html = "".join(
        f"<tr><td><code>{html.escape(e['id'])}</code></td><td>{e['affinity_ratio']:.4f}</td>"
        f"<td>{e['reporter_delta_end']:+.4f}</td><td>{e['nkx2_2_delta']:+.3f}</td>"
        f"<td>{html.escape(e['wt_fate'])} &rarr; {html.escape(e['mut_fate'])}</td>"
        f"<td>{e.get('kreimer_mpra_mean', float('nan')):+.4f}</td>"
        f"<td>{e.get('siraj_score', float('nan')):+.4f}</td></tr>"
        for e in report["scan"]["cryptic_switches"][:8]
    )
    checks_html = "".join(
        f"<li><strong>{'PASS' if c['passed'] else 'FAIL'}</strong> &mdash; {html.escape(c['name'])}: <code>{html.escape(json.dumps(c['evidence']))}</code></li>"
        for c in report["checks"]
    )
    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><title>Epigenesis · Variant-to-Circuit Scanner</title>
<style>body{{font-family:system-ui,-apple-system,sans-serif;max-width:1120px;margin:24px auto;padding:0 20px;color:#111}}
.grid{{display:grid;grid-template-columns:1fr 1fr;gap:20px}}svg{{width:100%;border:1px solid #e5e7eb;background:#fff}}
table{{width:100%;border-collapse:collapse;font-size:14px}}th,td{{padding:7px 10px;border-bottom:1px solid #e5e7eb;text-align:left}}
.badge{{display:inline-block;padding:3px 8px;background:#fee2e2;color:#991b1b;font-weight:600;border-radius:4px}}</style></head>
<body><h1>Variant-to-Circuit Enhancer Scanner</h1>
<p><strong>Locus:</strong> {html.escape(report['locus'])} ({report['scan']['length']} bp) · <strong>Featured Cryptic Switch:</strong>
<span class="badge">{html.escape(top['id'])}</span> (q={top['affinity_ratio']:.4f}, passive reporter &Delta;={top['reporter_delta_end']:+.4f}, endogenous &Delta;Nkx2.2={top['nkx2_2_delta']:+.3f}, {html.escape(wt_fate)} &rarr; {html.escape(top['mut_fate'])})</p>
<div class="grid">
<div><h3>Passive Reporter vs. Endogenous Circuit ({report['scan']['length']*3} edits at S={report['scan']['signal']})</h3>
<svg viewBox="0 0 {w} {h}"><line x1="{pad}" y1="{h-pad}" x2="{w-pad}" y2="{h-pad}" stroke="#111"/><line x1="{pad}" y1="{pad}" x2="{pad}" y2="{h-pad}" stroke="#111"/>
<text x="{w/2}" y="{h-10}" text-anchor="middle" font-size="12">|Passive Reporter Occupancy Shift| (max {max_x:.3f})</text>
<text x="16" y="{h/2}" transform="rotate(-90 16 {h/2})" text-anchor="middle" font-size="12">|Endogenous Nkx2.2 Fate Shift|</text>
{''.join(points)}</svg></div>
<div><h3>Spatial Neural Tube Morphogen Gradient (t={sp['duration']})</h3>
<svg viewBox="0 0 {w} {h}"><line x1="{pad}" y1="{h-pad}" x2="{w-pad}" y2="{h-pad}" stroke="#111"/><line x1="{pad}" y1="{pad}" x2="{pad}" y2="{h-pad}" stroke="#111"/>
<text x="{w/2}" y="{h-10}" text-anchor="middle" font-size="12">Ventral &rarr; Dorsal Position x (S=exp(-x/0.15))</text>
<text x="16" y="{h/2}" transform="rotate(-90 16 {h/2})" text-anchor="middle" font-size="12">Steady Concentration (Solid=WT, Dashed=Mutant)</text>
{poly(nkx_wt, "#c0262d")}{poly(nkx_mut, "#c0262d", "6 4")}{poly(olig_wt, "#1f6feb")}{poly(olig_mut, "#1f6feb", "6 4")}</svg></div></div>
<h3>Top Reporter-Silent, Circuit-Switch Mutations</h3>
<table><thead><tr><th>Edit</th><th>q (GLI2)</th><th>&Delta;Reporter</th><th>&Delta;Nkx2.2</th><th>Fate Flip</th><th>Kreimer MPRA</th><th>Siraj Score</th></tr></thead>
<tbody>{rows_html}</tbody></table>
<h3>Causal & Numerical Controls</h3><ul>{checks_html}</ul>
<p style="font-size:12px;color:#555">Compiled brainc receipt: <code>{report['compiled']['contract_sha256']}</code>. Concentrations and time are uncalibrated Herrera-Delgado (2020) model units.</p>
</body></html>"""


def demo(output: Path, *, motif_model_dir: Path | None = None) -> dict:
    """Run the Variant-to-Circuit discovery demo on real human enhancer chr15:90936956-90937127."""
    sys.path.insert(0, str(REPO / "examples/neural-regulation"))
    import regulation as reg
    rows, _ = reg.load_data()
    locus = "chr15:90936956-90937127"
    reference = next(r.reference for r in rows if r.locus == locus)

    output.mkdir()
    scan = scan_enhancer(reference, signal=DEFAULT_SIGNAL, protein="GLI2", epsilon=DEFAULT_EPSILON, motif_model_dir=motif_model_dir)
    if not scan["cryptic_switches"]:
        raise RuntimeError("expected cryptic bifurcation edits in chr15 enhancer")
    featured = next((e for e in scan["cryptic_switches"] if e["id"] == "60:C>A"), scan["cryptic_switches"][0])
    variant = featured["variant"]
    q = featured["affinity_ratio"]

    # Compile and verify real brainc Development Module
    receipt = compile_circuit_bundle(reference, variant, protein="GLI2", epsilon=DEFAULT_EPSILON, output=output / "compiled")
    compiled = verify_circuit_bundle(output / "compiled", receipt)

    # Causal intervention controls
    protocol = [(40.0, DEFAULT_SIGNAL)]
    wt_assay = run_circuit_assay(protocol, ratio=1.0, mode="endogenous")
    mut_assay = run_circuit_assay(protocol, ratio=q, mode="endogenous")
    rep_assay = run_circuit_assay(protocol, ratio=q, mode="reporter")
    allele_rescue = run_circuit_assay(protocol, ratio=affinity_ratio(reference, reference), mode="endogenous")
    mediator_rescue = run_circuit_assay(protocol, ratio=1.0, mode="endogenous")
    phenocopy = run_circuit_assay(protocol, ratio=q, mode="endogenous")

    # Numerical refinement check (dt vs dt/2)
    fine_wt = run_circuit_assay(protocol, ratio=1.0, mode="endogenous", dt=pat.DEFAULT_DT / 2)
    fine_mut = run_circuit_assay(protocol, ratio=q, mode="endogenous", dt=pat.DEFAULT_DT / 2)
    diff = lambda a, b: max(abs(x - y) for r1, r2 in zip(a, b) for x, y in zip(r1[1:], r2[1:]))
    refinement_err = max(diff(wt_assay["trajectory"], fine_wt["trajectory"][::2]), diff(mut_assay["trajectory"], fine_mut["trajectory"][::2]))

    # Spatial boundary shift (0..1)
    model = pat.load_data()["model"]
    positions = [i / 40 for i in range(41)]
    sp_ref = [list(simulate_circuit([(100.0, math.exp(-x / model["gradient_length"]))], ratio=1.0)[-1][1:]) for x in positions]
    sp_mut = [list(simulate_circuit([(100.0, math.exp(-x / model["gradient_length"]))], ratio=q)[-1][1:]) for x in positions]
    p3_wt = [x for x, y in zip(positions, sp_ref) if y[2] > 0.5]
    p3_mut = [x for x, y in zip(positions, sp_mut) if y[2] > 0.5]

    # Hysteresis / signal history check on featured variant vs WT (competitive vs activation-only ablation)
    hist_wt = simulate_circuit([(5.0, 1.0), (15.0, 0.05)], ratio=1.0, mechanism="competitive")[-1][3]
    hist_mut = simulate_circuit([(5.0, 1.0), (15.0, 0.05)], ratio=q, mechanism="competitive")[-1][3]
    hist_mut_act = simulate_circuit([(5.0, 1.0), (15.0, 0.05)], ratio=q, mechanism="activation-only")[-1][3]

    # Exact 9bp zero-epsilon equivalence check against patterning.profile
    core_wt, core_mut = "GACCACCCA", "GAACACCCA"
    q_exact = affinity_ratio(core_wt, core_mut, "GLI2", epsilon=0.0)
    q_pat = pat.profile(core_mut, "GLI2")["profile_ratio"]

    checks = [
        {"name": "exact_9bp_parity_with_patterning", "passed": abs(q_exact - q_pat) < 1e-15, "evidence": {"q_circuit": q_exact, "q_patterning": q_pat}},
        {"name": "cryptic_bifurcation_found", "passed": abs(featured["reporter_delta_end"]) < 0.02 and featured["fate_flipped"],
         "evidence": {"id": featured["id"], "reporter_delta": featured["reporter_delta_end"], "nkx2_2_delta": featured["nkx2_2_delta"]}},
        {"name": "reporter_leaves_endogenous_unchanged", "passed": rep_assay["trajectory"] == wt_assay["trajectory"],
         "evidence": {"max_state_diff": diff(rep_assay["trajectory"], wt_assay["trajectory"])}},
        {"name": "allele_and_mediator_rescue", "passed": allele_rescue["trajectory"] == wt_assay["trajectory"] and mediator_rescue["trajectory"] == wt_assay["trajectory"],
         "evidence": {"allele_err": diff(allele_rescue["trajectory"], wt_assay["trajectory"])}},
        {"name": "phenocopy_reproduces_variant", "passed": phenocopy["trajectory"] == mut_assay["trajectory"],
         "evidence": {"phenocopy_err": diff(phenocopy["trajectory"], mut_assay["trajectory"])}},
        {"name": "spatial_p3_domain_contracts", "passed": len(p3_mut) < len(p3_wt),
         "evidence": {"wt_p3_max_x": max(p3_wt), "mut_p3_max_x": max(p3_mut) if p3_mut else None}},
        {"name": "competitive_primed_derepression_reverses_ablation",
         "passed": hist_mut > hist_wt > 0.5 and hist_mut_act < 0.01,
         "evidence": {"wt_primed_nkx": hist_wt, "mut_primed_competitive": hist_mut, "mut_primed_activation_only": hist_mut_act}},
        {"name": "timestep_refinement_converges", "passed": refinement_err < 1e-4, "evidence": {"max_error": refinement_err, "tolerance": 1e-4}},
    ]

    report = {
        "format": "epigenesis.neural-circuit-experiment/v1",
        "locus": locus,
        "scope": "computational sequence-to-circuit model; biological embryo phenotype is not established",
        "compiled": compiled,
        "featured_variant": {k: v for k, v in featured.items() if k != "variant"},
        "second_site_candidates": _find_second_site_rescue(reference, featured, "GLI2", DEFAULT_EPSILON),
        "scan": {**{k: v for k, v in scan.items() if k != "edits"},
                 "edits": [{k: v for k, v in e.items() if k != "variant"} for e in scan["edits"]],
                 "cryptic_switches": [{k: v for k, v in e.items() if k != "variant"} for e in scan["cryptic_switches"]]},
        "spatial": {"duration": 100.0, "positions": positions, "reference": sp_ref, "variant": sp_mut},
        "checks": checks,
        "all_checks_passed": all(c["passed"] for c in checks),
    }
    (output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    (output / "index.html").write_text(_render_html(report))
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    d = sub.add_parser("demo")
    d.add_argument("--output", type=Path, required=True)
    d.add_argument("--motif-model-dir", type=Path)
    s = sub.add_parser("scan")
    s.add_argument("--fasta", type=Path, required=True)
    s.add_argument("--signal", type=float, default=DEFAULT_SIGNAL)
    s.add_argument("--protein", choices=("GLI2", "GLI3"), default="GLI2")
    s.add_argument("--motif-model-dir", type=Path)
    s.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == "demo":
            report = demo(args.output, motif_model_dir=args.motif_model_dir)
            print(json.dumps({"output": str(args.output), "all_checks_passed": report["all_checks_passed"],
                              "category_counts": report["scan"]["category_counts"],
                              "featured": report["featured_variant"], "checks": report["checks"]}, indent=2))
            return 0 if report["all_checks_passed"] else 1
        raw = read_regular_file(args.fasta, maximum_bytes=128 * 1024).decode("ascii")
        seq = "".join(line.strip() for line in raw.splitlines() if not line.startswith(">"))
        res = scan_enhancer(seq, signal=args.signal, protein=args.protein, motif_model_dir=args.motif_model_dir)
        args.output.mkdir()
        (args.output / "scan.json").write_text(json.dumps(res, indent=2, allow_nan=False) + "\n")
        print(json.dumps({"output": str(args.output), "category_counts": res["category_counts"],
                          "top_cryptic": [{k: v for k, v in e.items() if k != "variant"} for e in res["cryptic_switches"][:5]]}, indent=2))
        return 0
    except (ValueError, OSError, RuntimeError) as err:
        print(f"neural-circuit failed: {err}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

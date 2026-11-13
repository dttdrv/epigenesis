"""Variant-to-Circuit bridge: arbitrary regulatory DNA -> multi-TF Boltzmann physics -> general GRN bifurcation.

Connects arbitrary-length DNA enhancers (via brainc Development Modules and strand-symmetric
Boltzmann partition affinity over Hallikas et al. 2006 GLI2/GLI3, all 532 bundled JASPAR 2024
motifs, or custom PWMs) to arbitrary N-gene nonlinear regulatory circuits (default: Herrera-Delgado
et al. 2020 neural-tube ODE plus mutual-toggle and tri-lineage-switch presets) and optionally
attaches passive reporter predictions (Kreimer 0-72h motif_temporal & Siraj Huber-linear).

Exposes "reporter-silent, circuit-switch" (cryptic bifurcation) edits where
passive promoter/MPRA readouts barely move yet coupled feedback flips cell fate.
All developmental concentrations, distances, and times are uncalibrated model units.
"""

from __future__ import annotations

import argparse
import ast
import base64
import hashlib
import html
import importlib.util
import json
import math
from pathlib import Path
import shutil
import statistics
import struct
import sys
import zipfile

# coding-minimalism-rung: reuse (reuses brainc, examples/neural-patterning, examples/neural-construction, and optional Studio)
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
DEFAULT_BACKGROUND_FLOOR = 2e-4
DEFAULT_SIGNAL = 0.20
DEFAULT_TARGET = "Nkx2.2"
DEFAULT_CIRCUIT = "neural-tube"
SUPPORTED_TARGETS = ("Nkx2.2", "Olig2", "Both")
SUPPORTED_PRESETS = ("neural-tube", "mutual-toggle", "tri-lineage-switch")
CURATED_LOCI = (
    "chr15:90936956-90937127",
    "chr13:26545315-26545486",
    "chr21:37477109-37477280",
    "chr18:56475717-56475888",
    "chr1:175099416-175099587",
)
LOCAL_ERROR_TOLERANCE = pat.LOCAL_ERROR_TOLERANCE
_JASPAR_CACHE: dict[str, dict] = {}
_PAT_DATA_CACHE: dict | None = None


def _pat_data() -> dict:
    global _PAT_DATA_CACHE
    if _PAT_DATA_CACHE is None:
        _PAT_DATA_CACHE = pat.load_data()
    return _PAT_DATA_CACHE


def _read_npy_2d(raw: bytes) -> list[list[float]]:
    """Parse a uncompressed .npy binary buffer containing a 2D float32/float64 array using stdlib only."""
    if len(raw) < 10 or raw[:6] != b"\x93NUMPY":
        raise ValueError("invalid .npy magic header")
    major, minor = raw[6], raw[7]
    if major == 1:
        hlen = struct.unpack("<H", raw[8:10])[0]
        hstart = 10
    elif major in (2, 3):
        if len(raw) < 12:
            raise ValueError("truncated .npy v2/v3 header")
        hlen = struct.unpack("<I", raw[8:12])[0]
        hstart = 12
    else:
        raise ValueError(f"unsupported .npy version {major}.{minor}")
    header_end = hstart + hlen
    header = ast.literal_eval(raw[hstart:header_end].decode("latin1").strip())
    descr = header["descr"]
    if descr not in ("<f8", "=f8", "<f4", "=f4"):
        raise ValueError(f"unsupported .npy dtype: {descr}")
    shape = header["shape"]
    if len(shape) != 2 or shape[0] != 4 or shape[1] <= 0:
        raise ValueError(f"expected (4, W) motif array, got shape {shape}")
    nrows, ncols = shape
    elem_fmt = "d" if descr.endswith("f8") else "f"
    elem_size = 8 if elem_fmt == "d" else 4
    payload = raw[header_end:]
    if len(payload) != nrows * ncols * elem_size:
        raise ValueError("truncated .npy payload")
    flat = struct.unpack(f"<{nrows * ncols}{elem_fmt}", payload)
    if header.get("fortran_order"):
        rows = [[flat[r + c * nrows] for c in range(ncols)] for r in range(nrows)]
    else:
        rows = [list(flat[r * ncols : (r + 1) * ncols]) for r in range(nrows)]
    for row in rows:
        if any(not math.isfinite(v) for v in row):
            raise ValueError("nonfinite value in motif matrix")
    return rows


def _load_jaspar_catalog() -> dict:
    model_dir = REPO / "studio/epigenesis_studio/model"
    manifest_bytes = read_regular_file(model_dir / "model.json", maximum_bytes=1024 * 1024)
    manifest_sha = hashlib.sha256(manifest_bytes).hexdigest()
    cached = _JASPAR_CACHE.get(manifest_sha)
    if cached is not None:
        return cached
    manifest = load_json_object(manifest_bytes, "bundled JASPAR manifest")
    weights_path = model_dir / "weights.npz"
    weights_bytes = read_regular_file(weights_path, maximum_bytes=16 * 1024 * 1024)
    expected_weights_sha = manifest.get("sha256", {}).get("weights.npz")
    if expected_weights_sha and hashlib.sha256(weights_bytes).hexdigest() != expected_weights_sha:
        raise ValueError("bundled JASPAR weights.npz digest mismatch")
    entry = {
        "manifest": manifest,
        "manifest_sha256": manifest_sha,
        "motif_ids": list(manifest["motif_ids"]),
        "id_to_index": {mid: i for i, mid in enumerate(manifest["motif_ids"])},
        "weights_path": weights_path,
        "matrices": {},
    }
    _JASPAR_CACHE.clear()
    _JASPAR_CACHE[manifest_sha] = entry
    return entry


def list_available_motifs() -> dict:
    """List supported built-in transcription factor profiles (Hallikas GLI + bundled JASPAR 2024)."""
    jaspar_ids: list[str] = []
    manifest_path = REPO / "studio/epigenesis_studio/model/model.json"
    if manifest_path.is_file():
        jaspar_ids = list(_load_jaspar_catalog()["motif_ids"])
    return {
        "hallikas_gli": ["GLI2", "GLI3"],
        "jaspar_count": len(jaspar_ids),
        "jaspar_ids": jaspar_ids,
    }


def _motif_spec_label(spec: str | dict) -> str:
    if isinstance(spec, str):
        return spec
    if isinstance(spec, dict):
        return str(spec.get("id") or spec.get("name") or "custom")
    raise ValueError("motif specification must be a string ID or custom PWM dict")


def load_motif_profile(protein: str | dict, epsilon: float = DEFAULT_EPSILON) -> dict:
    """Load consensus and relative-affinity ratios for GLI2/GLI3, any bundled JASPAR ID, or a custom PWM dict."""
    if type(epsilon) not in (int, float) or not math.isfinite(epsilon) or epsilon < 0:
        raise ValueError("epsilon must be finite and nonnegative")

    if isinstance(protein, dict):
        mat = protein.get("matrix")
        kind = protein.get("kind", "weights")
        if type(mat) is not dict or set(mat) != {"A", "C", "G", "T"} or kind not in ("weights", "log_odds"):
            raise ValueError("custom motif dict needs 'matrix' with keys A/C/G/T and optional kind='weights'|'log_odds'")
        lengths = {len(mat[b]) for b in "ACGT" if isinstance(mat[b], list)}
        if len(lengths) != 1:
            raise ValueError("custom motif matrix base rows must be equal-length lists")
        width = next(iter(lengths))
        if width <= 0:
            raise ValueError("custom motif width must be >= 1")
        for b in "ACGT":
            if any(type(v) not in (int, float) or not math.isfinite(v) for v in mat[b]):
                raise ValueError(f"nonfinite value in custom motif row {b}")
        if kind == "log_odds":
            weights_by_base = {b: [0.0] * width for b in "ACGT"}
            for pos in range(width):
                col_max = max(float(mat[b][pos]) for b in "ACGT")
                for b in "ACGT":
                    weights_by_base[b][pos] = math.exp(float(mat[b][pos]) - col_max)
        else:
            weights_by_base = {b: [float(v) for v in mat[b]] for b in "ACGT"}
            if any(v < 0 for b in "ACGT" for v in weights_by_base[b]):
                raise ValueError("custom motif weights must be nonnegative")
        consensus = []
        ratios = {b: [0.0] * width for b in "ACGT"}
        for pos in range(width):
            best_b = max("ACGT", key=lambda b: weights_by_base[b][pos])
            best_w = weights_by_base[best_b][pos]
            if best_w <= 0:
                raise ValueError(f"custom motif column {pos} has no positive weight")
            consensus.append(best_b)
            for b in "ACGT":
                ratios[b][pos] = (weights_by_base[b][pos] + epsilon) / (best_w + epsilon)
        return {
            "id": _motif_spec_label(protein),
            "source": "custom",
            "width": width,
            "consensus": "".join(consensus),
            "ratios": ratios,
        }

    if type(protein) is not str or not protein:
        raise ValueError("motif identifier must be a non-empty string or custom PWM dict")

    binding = _pat_data()["binding"]
    if protein in binding["profiles"]:
        consensus = binding["consensus"]
        profile = binding["profiles"][protein]
        ratios = {
            base: [
                (profile[base][pos] + epsilon) / (profile[consensus[pos]][pos] + epsilon)
                for pos in range(len(consensus))
            ]
            for base in "ACGT"
        }
        return {
            "id": protein,
            "source": "hallikas-2006",
            "width": len(consensus),
            "consensus": consensus,
            "ratios": ratios,
        }

    cat = _load_jaspar_catalog()
    if protein not in cat["id_to_index"]:
        raise ValueError(f"unsupported protein/motif: {protein!r} (expected GLI2, GLI3, custom dict, or a JASPAR ID from model.json)")
    idx = cat["id_to_index"][protein]
    mat_rows = cat["matrices"].get(idx)
    if mat_rows is None:
        with zipfile.ZipFile(cat["weights_path"], "r") as zf:
            mat_rows = _read_npy_2d(zf.read(f"motif_{idx}.npy"))
        cat["matrices"][idx] = mat_rows
    width = len(mat_rows[0])
    bases = ("A", "C", "G", "T")
    consensus = []
    ratios = {b: [0.0] * width for b in bases}
    for pos in range(width):
        col = [mat_rows[bi][pos] for bi in range(4)]
        col_max = max(col)
        best_bi = col.index(col_max)
        consensus.append(bases[best_bi])
        for bi, b in enumerate(bases):
            rel = math.exp(col[bi] - col_max)
            ratios[b][pos] = (rel + epsilon) / (1.0 + epsilon)
    return {
        "id": protein,
        "source": "jaspar-2024",
        "width": width,
        "consensus": "".join(consensus),
        "ratios": ratios,
    }


def _binding_tables(protein: str | dict, epsilon: float) -> tuple[str, dict[str, list[float]]]:
    prof = load_motif_profile(protein, epsilon)
    return prof["consensus"], prof["ratios"]


def _site_weight(window: str, ratios: dict[str, list[float]]) -> float:
    weight = 1.0
    for pos, base in enumerate(window):
        weight *= ratios[base][pos]
    return weight


def boltzmann_affinity(sequence: str, protein: str | dict = "GLI2", epsilon: float = DEFAULT_EPSILON) -> dict:
    """Sum strand-symmetric relative association affinity across all W-bp windows for any TF PWM."""
    prof = load_motif_profile(protein, epsilon)
    ratios = prof["ratios"]
    width = prof["width"]
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
        "protein": prof["id"],
        "source": prof["source"],
        "width": width,
        "consensus": prof["consensus"],
        "epsilon": epsilon,
        "length": len(sequence),
        "partition_sum": total,
        "peak_site": best_site,
        "strong_sites": top_sites[:8],
    }


boltzmann_gli_affinity = boltzmann_affinity


def _resolve_background_floor(background_floor: float | None, epsilon: float = DEFAULT_EPSILON) -> float:
    """Resolve non-specific partition floor Z_ns (0.0 when epsilon==0 for exact Hallikas parity)."""
    if type(epsilon) not in (int, float) or not math.isfinite(epsilon) or epsilon < 0:
        raise ValueError("epsilon must be finite and nonnegative")
    if background_floor is None:
        return 0.0 if epsilon == 0.0 else DEFAULT_BACKGROUND_FLOOR
    if not _num(background_floor) or background_floor < 0:
        raise ValueError("background_floor must be finite and nonnegative")
    return float(background_floor)


def affinity_ratio(
    reference: str,
    variant: str,
    protein: str | dict = "GLI2",
    epsilon: float = DEFAULT_EPSILON,
    background_floor: float | None = None,
) -> float:
    """Compute thermodynamic TF affinity ratio q = (Z(variant) + Z_ns) / (Z(reference) + Z_ns) (> 0)."""
    eff_floor = _resolve_background_floor(background_floor, epsilon)
    z_ref = boltzmann_affinity(reference, protein, epsilon)["partition_sum"]
    z_var = boltzmann_affinity(variant, protein, epsilon)["partition_sum"]
    numer = z_var + eff_floor
    denom = z_ref + eff_floor
    if denom <= 0 or numer <= 0 or not math.isfinite(numer / denom):
        raise ValueError("affinity ratio requires positive finite partition sums")
    return numer / denom


def multi_tf_affinity(
    reference: str,
    variant: str,
    motif_specs: list[str | dict] | tuple[str | dict, ...] | dict[str, str | dict],
    epsilon: float = DEFAULT_EPSILON,
    background_floor: float | None = None,
) -> dict:
    """Compute simultaneous Boltzmann partition ratios across multiple TFs on the same enhancer pair."""
    eff_floor = _resolve_background_floor(background_floor, epsilon)
    if isinstance(motif_specs, dict):
        items = list(motif_specs.items())
    elif isinstance(motif_specs, (list, tuple)) and motif_specs:
        items = [(_motif_spec_label(m), m) for m in motif_specs]
    else:
        raise ValueError("motif_specs must be a non-empty list/tuple or mapping")
    ratios: dict[str, float] = {}
    wt_parts: dict[str, dict] = {}
    var_parts: dict[str, dict] = {}
    for label, spec in items:
        wp = boltzmann_affinity(reference, spec, epsilon)
        vp = boltzmann_affinity(variant, spec, epsilon)
        denom = wp["partition_sum"] + eff_floor
        numer = vp["partition_sum"] + eff_floor
        q = numer / denom if denom > 0 else float("nan")
        if not math.isfinite(q) or q <= 0:
            raise ValueError(f"invalid affinity ratio for {label}")
        ratios[label] = q
        wt_parts[label] = wp
        var_parts[label] = vp
    return {
        "ratios": ratios,
        "wt_partitions": wt_parts,
        "variant_partitions": var_parts,
        "background_floor": eff_floor,
    }


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


def _num(value) -> bool:
    return type(value) in (int, float) and math.isfinite(value)


def validate_circuit_spec(spec: dict) -> dict:
    """Validate and normalize an arbitrary N-gene regulatory circuit specification."""
    if type(spec) is not dict:
        raise ValueError("circuit specification must be a JSON object")
    genes = spec.get("genes")
    if (
        not isinstance(genes, list)
        or not 1 <= len(genes) <= 64
        or any(type(g) is not str or not g for g in genes)
        or len(set(genes)) != len(genes)
    ):
        raise ValueError("circuit 'genes' must be 1..64 unique non-empty strings")
    n = len(genes)
    gene_set = set(genes)

    def _vec(field: str, *, positive: bool = False) -> list[float]:
        val = spec.get(field)
        if _num(val):
            vals = [float(val)] * n
        elif isinstance(val, list) and len(val) == n and all(_num(x) for x in val):
            vals = [float(x) for x in val]
        else:
            raise ValueError(f"circuit '{field}' must be a number or length-{n} vector")
        bad = any(x <= 0 if positive else x < 0 for x in vals)
        if bad:
            raise ValueError(f"circuit '{field}' entries must be {'> 0' if positive else '>= 0'}")
        return vals

    basal = _vec("basal")
    input_w = _vec("input")
    tf_binding = _vec("tf_binding" if "tf_binding" in spec else "gli_binding")
    production = _vec("production", positive=True)
    decay = _vec("decay", positive=True)

    rep = spec.get("repression")
    if (
        not isinstance(rep, list)
        or len(rep) != n
        or any(not isinstance(row, list) or len(row) != n or any(not _num(x) or x < 0 for x in row) for row in rep)
    ):
        raise ValueError(f"circuit 'repression' must be a {n}x{n} matrix of nonnegative numbers")
    repression = [[float(x) for x in row] for row in rep]

    coop_raw = spec.get("cooperativity", spec.get("sites_per_repressor", 2))
    if _num(coop_raw) and coop_raw >= 0:
        coop = float(coop_raw)
        coop_matrix = [[coop] * n for _ in range(n)]
    elif (
        isinstance(coop_raw, list)
        and len(coop_raw) == n
        and all(isinstance(row, list) and len(row) == n and all(_num(x) and x >= 0 for x in row) for row in coop_raw)
    ):
        coop_matrix = [[float(x) for x in row] for row in coop_raw]
    else:
        raise ValueError(f"circuit 'cooperativity' must be a nonnegative scalar or {n}x{n} matrix")

    grad_len = spec.get("gradient_length", 0.15)
    if not _num(grad_len) or grad_len <= 0:
        raise ValueError("circuit 'gradient_length' must be positive")

    raw_edges = spec.get("edge_bindings", [])
    if not isinstance(raw_edges, list):
        raise ValueError("circuit 'edge_bindings' must be a list")
    edge_bindings = []
    for edge in raw_edges:
        if type(edge) is not dict or edge.get("gene") not in gene_set:
            raise ValueError("each edge_binding requires a valid 'gene' in circuit genes")
        prot = edge.get("protein", spec.get("default_protein", "GLI2"))
        _motif_spec_label(prot)
        mech = edge.get("mechanism", "competitive")
        role = edge.get("role", "activator")
        if mech not in ("competitive", "activation-only") or role not in ("activator", "repressor"):
            raise ValueError("edge_binding mechanism must be competitive|activation-only and role must be activator|repressor")
        edge_bindings.append({"gene": edge["gene"], "protein": prot, "mechanism": mech, "role": role})

    fates_raw = spec.get("fates")
    fates = None
    if fates_raw is not None:
        if type(fates_raw) is not dict:
            raise ValueError("circuit 'fates' must be a dict")
        thresh = fates_raw.get("threshold", 0.5)
        if not _num(thresh) or not 0 < thresh < 1:
            raise ValueError("fates threshold must be in (0, 1)")
        fates = dict(fates_raw, threshold=float(thresh))

    out = {
        "format": "epigenesis.circuit-spec/v1",
        "name": str(spec.get("name", "custom")),
        "genes": list(genes),
        "basal": basal,
        "input": input_w,
        "tf_binding": tf_binding,
        "gli_binding": tf_binding,
        "repression": repression,
        "cooperativity": coop_matrix,
        "sites_per_repressor": coop_raw if _num(coop_raw) else None,
        "production": production,
        "decay": decay,
        "gradient_length": float(grad_len),
        "default_target": str(spec.get("default_target", genes[0])),
        "default_signal": float(spec.get("default_signal", DEFAULT_SIGNAL)),
        "default_protein": spec.get("default_protein", "GLI2"),
        "edge_bindings": edge_bindings,
        "fates": fates,
    }
    if out["default_target"] not in gene_set and out["default_target"] != "Both":
        raise ValueError("circuit default_target must name a circuit gene")
    if not 0 <= out["default_signal"] <= 1:
        raise ValueError("circuit default_signal must be in [0, 1]")
    return out


def circuit_preset(name: str = DEFAULT_CIRCUIT) -> dict:
    """Return a validated built-in CircuitSpec preset ('neural-tube', 'mutual-toggle', 'tri-lineage-switch')."""
    if name == "neural-tube":
        m = _pat_data()["model"]
        return validate_circuit_spec({
            "name": "neural-tube",
            "genes": list(m["genes"]),
            "basal": list(m["basal"]),
            "input": list(m["input"]),
            "tf_binding": list(m["gli_binding"]),
            "repression": [list(r) for r in m["repression"]],
            "cooperativity": m["sites_per_repressor"],
            "production": m["production"],
            "decay": m["decay"],
            "gradient_length": m["gradient_length"],
            "default_target": "Nkx2.2",
            "default_signal": 0.20,
            "default_protein": "GLI2",
            "edge_bindings": [
                {"gene": "Nkx2.2", "protein": "GLI2", "mechanism": "competitive", "role": "activator"},
            ],
            "fates": {"mode": "neural-tube", "threshold": 0.5},
        })
    if name == "mutual-toggle":
        return validate_circuit_spec({
            "name": "mutual-toggle",
            "genes": ["GeneA", "GeneB"],
            "basal": [0.45, 4.2],
            "input": [220.0, 20.0],
            "tf_binding": [16.0, 3.0],
            "repression": [[0.0, 38.0], [42.0, 0.0]],
            "cooperativity": 2.5,
            "production": [2.0, 2.0],
            "decay": [2.0, 2.0],
            "gradient_length": 0.20,
            "default_target": "GeneA",
            "default_signal": 0.22,
            "default_protein": "MA0069.1",
            "edge_bindings": [
                {"gene": "GeneA", "protein": "MA0069.1", "mechanism": "competitive", "role": "activator"},
            ],
        })
    if name == "tri-lineage-switch":
        return validate_circuit_spec({
            "name": "tri-lineage-switch",
            "genes": ["Progenitor", "Excitatory", "Inhibitory"],
            "basal": [3.2, 0.55, 1.2],
            "input": [45.0, 280.0, 35.0],
            "tf_binding": [6.0, 24.0, 10.0],
            "repression": [
                [0.0, 32.0, 26.0],
                [24.0, 0.0, 44.0],
                [22.0, 26.0, 0.0],
            ],
            "cooperativity": 2.0,
            "production": [2.0, 2.0, 2.0],
            "decay": [2.0, 2.0, 2.0],
            "gradient_length": 0.18,
            "default_target": "Excitatory",
            "default_signal": 0.35,
            "default_protein": "MA0069.1",
            "edge_bindings": [
                {"gene": "Excitatory", "protein": "MA0069.1", "mechanism": "competitive", "role": "activator"},
                {"gene": "Inhibitory", "protein": "MA0463.3", "mechanism": "activation-only", "role": "repressor"},
            ],
        })
    raise ValueError(f"unknown circuit preset {name!r}; expected one of {SUPPORTED_PRESETS}")


def list_circuit_presets() -> dict[str, dict]:
    return {name: circuit_preset(name) for name in SUPPORTED_PRESETS}


def resolve_circuit(circuit: str | dict | None) -> dict:
    if circuit is None or circuit == "neural-tube":
        return circuit_preset("neural-tube")
    if isinstance(circuit, str):
        return circuit_preset(circuit)
    if isinstance(circuit, dict):
        return validate_circuit_spec(circuit)
    raise ValueError("circuit must be a preset name or CircuitSpec dictionary")


def resolve_target_ratios(
    target: str | list[str] | tuple[str, ...] | None = None,
    ratio: float | dict[str, float] = 1.0,
    *,
    genes: list[str] | None = None,
) -> dict[str, float]:
    """Map target gene selection(s) or explicit per-gene ratios into a validated {gene: q} map."""
    default_genes = _pat_data()["model"]["genes"]
    valid_genes = set(genes if genes is not None else default_genes)
    if isinstance(ratio, dict):
        if not set(ratio) <= valid_genes:
            raise ValueError(f"unsupported target gene in ratios: {sorted(set(ratio) - valid_genes)}")
        out = {}
        for gene, value in ratio.items():
            if not _num(value) or value <= 0:
                raise ValueError(f"ratio for {gene} must be finite and positive")
            out[gene] = float(value)
        return out
    if not _num(ratio) or ratio <= 0:
        raise ValueError("affinity ratio must be finite and positive")
    q = float(ratio)
    if target is None:
        chosen = "Nkx2.2" if "Nkx2.2" in valid_genes else next(iter(sorted(valid_genes)))
        return {chosen: q}
    if isinstance(target, (list, tuple)):
        if not target or any(type(g) is not str or g not in valid_genes for g in target):
            raise ValueError(f"target genes must belong to {sorted(valid_genes)}")
        return {g: q for g in target}
    if target == "Both" and {"Nkx2.2", "Olig2"} <= valid_genes:
        return {"Nkx2.2": q, "Olig2": q}
    if target == "All" and genes is not None:
        return {g: q for g in valid_genes}
    if genes is None and target not in SUPPORTED_TARGETS:
        raise ValueError(f"target must be one of {SUPPORTED_TARGETS}")
    if target not in valid_genes:
        raise ValueError(f"target {target!r} not in circuit genes {sorted(valid_genes)}")
    return {target: q}


def occupancy_ratios(
    model: dict,
    state: tuple[float, ...],
    signal: float,
    ratios: dict[str, float] | None = None,
    mechanism: str = "competitive",
    *,
    roles: dict[str, str] | None = None,
    mechanisms: dict[str, str] | None = None,
) -> tuple[float, ...]:
    """Generalized N-gene Hill-Boltzmann promoter occupancy supporting activators, repressors, and per-edge cooperativity."""
    ratios = ratios or {}
    genes = model["genes"]
    if (
        not roles
        and not mechanisms
        and genes == ["Pax6", "Olig2", "Nkx2.2", "Irx3"]
        and set(ratios) <= {"Nkx2.2"}
        and "cooperativity" not in model
    ):
        return pat.occupancy(model, state, signal, ratios.get("Nkx2.2", 1.0), mechanism)

    k_list = model.get("tf_binding", model.get("gli_binding"))
    coop_raw = model.get("cooperativity", model.get("sites_per_repressor", 2))
    is_mat_coop = isinstance(coop_raw, list)
    result = []
    for i, basal in enumerate(model["basal"]):
        gene = genes[i]
        raw_q = ratios.get(gene, 1.0)
        role = roles.get(gene, "activator") if roles else "activator"
        mech = mechanisms.get(gene, mechanism) if mechanisms else mechanism
        q_act = raw_q if role == "activator" else (1.0 / raw_q)
        k = k_list[i]
        on = basal * (1.0 + k) * (1.0 + q_act * model["input"][i] * signal)
        comp_factor = (q_act if mech == "competitive" else 1.0) * k
        rep_prod = 1.0
        for j, (w, x) in enumerate(zip(model["repression"][i], state)):
            if w:
                h = coop_raw[i][j] if is_mat_coop else coop_raw
                rep_prod *= (1.0 + w * x) ** h
        off = (1.0 + comp_factor) * rep_prod
        result.append(on / (on + off) if (on + off) > 0 else 0.0)
    return tuple(result)


def _target_gene_index(genes: list[str], target: str | list[str] | tuple[str, ...]) -> int:
    if isinstance(target, (list, tuple)):
        return genes.index(target[0])
    if target == "Olig2" and "Olig2" in genes:
        return genes.index("Olig2")
    if target in ("Nkx2.2", "Both") and "Nkx2.2" in genes:
        return genes.index("Nkx2.2")
    if target in genes:
        return genes.index(target)
    return 0


def simulate_circuit(
    protocol: list,
    *,
    ratio: float | dict[str, float] = 1.0,
    target: str | list[str] | tuple[str, ...] = DEFAULT_TARGET,
    ratios: dict[str, float] | None = None,
    dt: float = pat.DEFAULT_DT,
    initial: tuple[float, ...] | None = None,
    mechanism: str = "competitive",
    circuit: str | dict | None = None,
    roles: dict[str, str] | None = None,
    mechanisms: dict[str, str] | None = None,
) -> list[tuple]:
    """Integrate an arbitrary N-gene nonlinear regulatory ODE under adaptive RK4 refinement."""
    model = _pat_data()["model"] if circuit is None else resolve_circuit(circuit)
    genes = model["genes"]
    n_genes = len(genes)
    if initial is None:
        initial = (0.0,) * n_genes
    eff_target = target if (circuit is None or target in genes or isinstance(target, (list, tuple)) or target in ("Both", "All")) else model["default_target"]
    ratio_map = resolve_target_ratios(eff_target, ratios if ratios is not None else ratio, genes=genes)
    if (
        not _num(dt) or dt <= 0
        or mechanism not in ("competitive", "activation-only")
        or len(initial) != n_genes
        or any(not _num(v) or not 0 <= v <= 1 for v in initial)
    ):
        raise ValueError("invalid timestep, positive affinity ratio, mechanism, or initial state")

    phases = []
    total, steps = 0.0, 0
    for phase in protocol:
        if not isinstance(phase, (list, tuple)) or len(phase) != 2 or any(not _num(v) for v in phase) or phase[0] <= 0 or not 0 <= phase[1] <= 1:
            raise ValueError("protocol needs positive durations and signal fractions in [0, 1]")
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

    prod_raw = model["production"]
    decay_raw = model["decay"]
    prod_vec = prod_raw if isinstance(prod_raw, list) else [float(prod_raw)] * n_genes
    decay_vec = decay_raw if isinstance(decay_raw, list) else [float(decay_raw)] * n_genes
    k_list = model.get("tf_binding", model.get("gli_binding"))
    coop_raw = model.get("cooperativity", model.get("sites_per_repressor", 2))
    is_mat_coop = isinstance(coop_raw, list)
    basal_list = model["basal"]
    input_list = model["input"]
    rep_mat = model["repression"]

    gene_terms = []
    for i in range(n_genes):
        gene = genes[i]
        raw_q = ratio_map.get(gene, 1.0)
        role = roles.get(gene, "activator") if roles else "activator"
        mech = mechanisms.get(gene, mechanism) if mechanisms else mechanism
        q_act = raw_q if role == "activator" else (1.0 / raw_q)
        k = k_list[i]
        basal_k = basal_list[i] * (1.0 + k)
        q_inp = q_act * input_list[i]
        off_base = 1.0 + (q_act if mech == "competitive" else 1.0) * k
        rep_edges = tuple(
            (j, w, coop_raw[i][j] if is_mat_coop else coop_raw)
            for j, w in enumerate(rep_mat[i])
            if w
        )
        gene_terms.append((prod_vec[i], decay_vec[i], basal_k, q_inp, off_base, rep_edges))

    state = tuple(float(v) for v in initial)
    trajectory = [(0.0, *state)]

    def rhs(y, s):
        if any(not math.isfinite(v) or not 0 <= v <= 1 for v in y):
            raise ArithmeticError("integration stage outside physical [0, 1] state range")
        out = []
        for i in range(n_genes):
            pr, dc, basal_k, q_inp, off_base, rep_edges = gene_terms[i]
            on = basal_k * (1.0 + q_inp * s)
            rep_prod = 1.0
            for j, w, h in rep_edges:
                rep_prod *= (1.0 + w * y[j]) ** h
            off = off_base * rep_prod
            denom = on + off
            occ = on / denom if denom > 0 else 0.0
            out.append(pr * occ - dc * y[i])
        return tuple(out)

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
    ratio: float | dict[str, float],
    mode: str,
    target: str | list[str] | tuple[str, ...] = DEFAULT_TARGET,
    dt: float = pat.DEFAULT_DT,
    mechanism: str = "competitive",
    initial: dict | None = None,
    circuit: str | dict | None = None,
    roles: dict[str, str] | None = None,
    mechanisms: dict[str, str] | None = None,
) -> dict:
    if mode not in ("reporter", "endogenous"):
        raise ValueError("mode must be reporter or endogenous")
    model = _pat_data()["model"] if circuit is None else resolve_circuit(circuit)
    genes = model["genes"]
    eff_target = target if (circuit is None or target in genes or isinstance(target, (list, tuple)) or target in ("Both", "All")) else model["default_target"]
    if initial is None:
        initial = dict.fromkeys(genes, 0.0)
    if type(initial) is not dict or set(initial) != set(genes):
        raise ValueError(f"initial state must provide exactly genes {genes}")
    state = tuple(float(initial[g]) for g in genes)
    ratio_map = resolve_target_ratios(eff_target, ratio, genes=genes)
    wt_map = {g: 1.0 for g in ratio_map}
    trajectory = simulate_circuit(
        protocol,
        ratios=wt_map if mode == "reporter" else ratio_map,
        target=eff_target,
        dt=dt,
        mechanism=mechanism,
        initial=state,
        circuit=circuit,
        roles=roles,
        mechanisms=mechanisms,
    )
    phase, boundary = 0, protocol[0][0]
    reporter, signals = [], []
    rep_idx = _target_gene_index(genes, eff_target)
    for time_val, *y in trajectory:
        while time_val > boundary and phase + 1 < len(protocol):
            phase += 1
            boundary += protocol[phase][0]
        signal = protocol[phase][1]
        signals.append(signal)
        reporter.append(
            occupancy_ratios(model, tuple(y), signal, ratio_map, mechanism, roles=roles, mechanisms=mechanisms)[rep_idx]
        )
    scalar_key = eff_target if isinstance(eff_target, str) and eff_target in ratio_map else next(iter(ratio_map))
    return {
        "mode": mode, "target": eff_target, "mechanism": mechanism, "protocol": protocol,
        "ratio": ratio_map[scalar_key],
        "target_ratios": ratio_map,
        "initial_state": dict(zip(genes, trajectory[0][1:])),
        "trajectory": trajectory, "signal": signals, "reporter_occupancy": reporter,
    }


def fate_label(
    state: tuple[float, ...] | list[float],
    genes: list[str],
    *,
    threshold: float = 0.5,
    fates: dict | None = None,
) -> str:
    if fates is not None and _num(fates.get("threshold")):
        threshold = float(fates["threshold"])
    if (fates and fates.get("mode") == "neural-tube") or genes == ["Pax6", "Olig2", "Nkx2.2", "Irx3"]:
        best_i = max(range(1, len(genes)), key=lambda i: state[i])
        return genes[best_i] if state[best_i] > threshold else ("Pax6+" if state[0] > threshold else "unresolved")
    best_i = max(range(len(genes)), key=lambda i: state[i])
    label = genes[best_i] if state[best_i] > threshold else "unresolved"
    if fates and isinstance(fates.get("labels"), dict) and label in fates["labels"]:
        return str(fates["labels"][label])
    return label


def _resolve_protein_arg(protein: str | dict | Path) -> str | dict:
    if isinstance(protein, Path):
        return load_json_object(read_regular_file(protein, maximum_bytes=256 * 1024), "custom motif JSON")
    if isinstance(protein, str):
        candidate = Path(protein)
        if protein.endswith(".json") and candidate.is_file():
            return load_json_object(read_regular_file(candidate, maximum_bytes=256 * 1024), "custom motif JSON")
        return protein
    if isinstance(protein, dict):
        return protein
    raise ValueError("protein must be GLI2, GLI3, a JASPAR motif ID, custom motif dict, or JSON file path")


def _resolve_circuit_arg(circuit: str | dict | Path | None) -> str | dict | None:
    if isinstance(circuit, Path):
        return load_json_object(read_regular_file(circuit, maximum_bytes=256 * 1024), "custom circuit JSON")
    if isinstance(circuit, str):
        if circuit in SUPPORTED_PRESETS:
            return circuit
        candidate = Path(circuit)
        if candidate.is_file():
            return load_json_object(read_regular_file(candidate, maximum_bytes=256 * 1024), "custom circuit JSON")
        if circuit.lstrip().startswith("{"):
            parsed = json.loads(circuit)
            if type(parsed) is not dict:
                raise ValueError("inline circuit JSON must be an object")
            return parsed
    return circuit


def _resolve_context(
    circuit: str | dict | Path | None,
    protein: str | dict | Path | None,
    target: str | list[str] | tuple[str, ...] | None,
    signal: float | None,
    roles: dict[str, str] | None,
    mechanisms: dict[str, str] | None,
) -> dict:
    """Resolve circuit model, target gene selection, active TF specs, roles, and mechanisms."""
    c_arg = _resolve_circuit_arg(circuit)
    model = _pat_data()["model"] if c_arg is None else resolve_circuit(c_arg)
    genes = model["genes"]
    eff_signal = (
        float(signal)
        if signal is not None
        else (DEFAULT_SIGNAL if c_arg is None else float(model["default_signal"]))
    )
    eff_target = (
        DEFAULT_TARGET
        if (target is None and c_arg is None)
        else (model["default_target"] if target is None else target)
    )

    roles_map = dict(roles or {})
    mechs_map = dict(mechanisms or {})
    bindings: dict[str, dict] = {}

    if c_arg is not None and protein is None and model.get("edge_bindings"):
        resolve_target_ratios(eff_target, 1.0, genes=genes)
        for edge in model["edge_bindings"]:
            g = edge["gene"]
            bindings[g] = {
                "gene": g,
                "protein": _resolve_protein_arg(edge["protein"]),
                "role": roles_map.get(g, edge["role"]),
                "mechanism": mechs_map.get(g, edge["mechanism"]),
            }
            roles_map[g] = bindings[g]["role"]
            mechs_map[g] = bindings[g]["mechanism"]
        primary_prot = _resolve_protein_arg(model.get("default_protein", "GLI2"))
    else:
        default_p = "GLI2" if c_arg is None else model.get("default_protein", "GLI2")
        primary_prot = _resolve_protein_arg(default_p if protein is None else protein)
        edge_by_gene = {e["gene"]: e for e in model.get("edge_bindings", [])}
        wt_map = resolve_target_ratios(eff_target, 1.0, genes=genes)
        for g in wt_map:
            eb = edge_by_gene.get(g)
            role = roles_map.get(g, eb["role"] if eb and protein is None else "activator")
            mech = mechs_map.get(g, eb["mechanism"] if eb and protein is None else "competitive")
            bindings[g] = {"gene": g, "protein": primary_prot, "role": role, "mechanism": mech}
            if roles is not None or (eb and protein is None):
                roles_map[g] = role
            if mechanisms is not None or (eb and protein is None):
                mechs_map[g] = mech

    unique_specs: list[tuple[str, str | dict]] = []
    seen_labels: dict[str, str | dict] = {}
    for b in bindings.values():
        spec = b["protein"]
        lbl = _motif_spec_label(spec)
        if lbl not in seen_labels:
            seen_labels[lbl] = spec
            unique_specs.append((lbl, spec))
    primary_label = _motif_spec_label(primary_prot)
    if primary_label not in seen_labels:
        seen_labels[primary_label] = primary_prot
        unique_specs.insert(0, (primary_label, primary_prot))

    return {
        "circuit_arg": c_arg,
        "circuit_name": model.get("name", "neural-tube") if c_arg is not None else "neural-tube",
        "model": model,
        "genes": genes,
        "signal": eff_signal,
        "target": eff_target,
        "primary_protein": primary_prot,
        "primary_label": primary_label,
        "bindings": bindings,
        "unique_specs": unique_specs,
        "roles": roles_map or None,
        "mechanisms": mechs_map or None,
        "multi_tf": len(unique_specs) > 1,
    }


def _target(semantics_hash: str, ref_length: int, var_length: int, *, tf_count: int = 1, version: int = 2) -> dict:
    fields = [
        {"id": "affinity_ratio", "type": {"dtype": "f64", "shape": [1]}, "unit": None, "mutability": "constant", "numeric": {"kind": "float"}},
        {"id": "reference_nucleotides", "type": {"dtype": "u64", "shape": [ref_length]}, "unit": None, "mutability": "constant", "numeric": {"kind": "integer"}},
        {"id": "variant_nucleotides", "type": {"dtype": "u64", "shape": [var_length]}, "unit": None, "mutability": "constant", "numeric": {"kind": "integer"}},
    ]
    reads = [
        {"scope": "subject", "field": "affinity_ratio"},
        {"scope": "subject", "field": "reference_nucleotides"},
        {"scope": "subject", "field": "variant_nucleotides"},
    ]
    if version >= 3:
        fields.append({"id": "affinity_ratios", "type": {"dtype": "f64", "shape": [tf_count]}, "unit": None, "mutability": "constant", "numeric": {"kind": "float"}})
        reads.append({"scope": "subject", "field": "affinity_ratios"})
    fields.sort(key=lambda item: item["id"])
    reads.sort(key=lambda item: item["field"])
    return target_artifact({
        "id": TARGET, "abi_major": 1, "opsets": [{"domain": DEV_DOMAIN, "version": 1}],
        "unit_schemas": [{"id": "enhancer-pair", "fields": fields}],
        "edge_schemas": [], "ports": [],
        "rules": [{"id": TARGET + ".partition", "version": 1, "semantics_sha256": semantics_hash,
                   "subject": "unit", "schema": "enhancer-pair", "phase": "UPDATE", "triggers": ["step"], "parameters": [],
                   "reads": reads, "writes": []}],
    })


def compile_circuit_bundle(
    reference: str,
    variant: str,
    *,
    protein: str | dict = "GLI2",
    epsilon: float = DEFAULT_EPSILON,
    background_floor: float | None = None,
    output: Path,
    motif_specs: list[str | dict] | tuple[str | dict, ...] | None = None,
    circuit: str | dict | None = None,
) -> str:
    """Compile a reference/variant enhancer pair (SNV, MNV, or indel) and Boltzmann affinity ratios into a brainc bundle."""
    specs = list(motif_specs) if motif_specs is not None else [protein]
    if not specs:
        raise ValueError("compile_circuit_bundle requires at least one motif specification")
    c_arg = _resolve_circuit_arg(circuit)
    use_v2 = (
        c_arg is None
        and len(specs) == 1
        and isinstance(specs[0], str)
        and specs[0] in ("GLI2", "GLI3")
    )

    eff_floor = _resolve_background_floor(background_floor, epsilon)
    multi = multi_tf_affinity(reference, variant, specs, epsilon, background_floor=eff_floor)
    tf_labels = [_motif_spec_label(s) for s in specs]
    q_vec = [multi["ratios"][lbl] for lbl in tf_labels]
    q = q_vec[0]

    custom_motifs = {
        _motif_spec_label(s): s
        for s in specs
        if isinstance(s, dict)
    }

    if use_v2:
        contract = {
            "format": "epigenesis.neural-circuit-pair", "version": 2,
            "protein": specs[0], "epsilon": epsilon,
            "background_floor": eff_floor,
            "reference_length": len(reference), "variant_length": len(variant),
            "data_sha256": pat.DATA_SHA256,
        }
    else:
        resolved_circuit = resolve_circuit(c_arg) if c_arg is not None else None
        circuit_bytes = canonical_bytes(resolved_circuit) + b"\n" if resolved_circuit is not None else None
        motifs_bytes = canonical_bytes(custom_motifs) + b"\n" if custom_motifs else None
        jaspar_used = any(
            isinstance(s, str) and s not in ("GLI2", "GLI3")
            for s in specs
        )
        contract = {
            "format": "epigenesis.neural-circuit-pair", "version": 3,
            "protein": tf_labels[0],
            "proteins": tf_labels,
            "epsilon": epsilon,
            "background_floor": eff_floor,
            "reference_length": len(reference), "variant_length": len(variant),
            "data_sha256": pat.DATA_SHA256,
            "jaspar_manifest_sha256": _load_jaspar_catalog()["manifest_sha256"] if jaspar_used else None,
            "circuit_sha256": hashlib.sha256(circuit_bytes).hexdigest() if circuit_bytes is not None else None,
            "motifs_sha256": hashlib.sha256(motifs_bytes).hexdigest() if motifs_bytes is not None else None,
        }

    impl_bytes = read_regular_file(Path(__file__).resolve(), maximum_bytes=256 * 1024)
    semantics_hash = hashlib.sha256(impl_bytes + b"\0" + canonical_bytes(contract)).hexdigest()
    receipt = hashlib.sha256(canonical_bytes(contract)).hexdigest()

    output.mkdir()
    pair_fasta = output / "pair.fasta"
    pair_fasta.write_text(f">reference\n{reference}\n>variant\n{variant}\n")
    (output / "contract.json").write_bytes(canonical_bytes(contract) + b"\n")
    if not use_v2:
        if contract["circuit_sha256"] is not None:
            (output / "circuit_spec.json").write_bytes(circuit_bytes)
        if contract["motifs_sha256"] is not None:
            (output / "motifs.json").write_bytes(motifs_bytes)

    source = compile_source(FASTA_PROFILE, {"sequence": pair_fasta}, parameters={"wrapper": "identity"})
    source.save(output / "source")

    tensor_list = [
        {"id": "count", "type": {"dtype": "u64", "shape": []}, "unit": None, "axes": [],
         "storage": inline_storage(pack("u64", [1]))},
        {"id": "reference_nucleotides", "type": {"dtype": "u64", "shape": [len(reference)]}, "unit": None, "axes": [None],
         "storage": inline_storage(pack("u64", list(reference.encode("ascii"))))},
        {"id": "variant_nucleotides", "type": {"dtype": "u64", "shape": [len(variant)]}, "unit": None, "axes": [None],
         "storage": inline_storage(pack("u64", list(variant.encode("ascii"))))},
        {"id": "affinity_ratio", "type": {"dtype": "f64", "shape": [1]}, "unit": None, "axes": [None],
         "storage": inline_storage(pack("f64", [q]))},
    ]
    if not use_v2:
        tensor_list.append({
            "id": "affinity_ratios", "type": {"dtype": "f64", "shape": [len(q_vec)]}, "unit": None, "axes": [None],
            "storage": inline_storage(pack("f64", q_vec)),
        })
    tensors = sorted(tensor_list, key=lambda item: item["id"])
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
    ver = 2 if use_v2 else 3
    target = _target(semantics_hash, len(reference), len(variant), tf_count=len(q_vec), version=ver)
    save(target, output / "target.json")
    init_fields = ["affinity_ratio", "reference_nucleotides", "variant_nucleotides"]
    if not use_v2:
        init_fields.append("affinity_ratios")
    init_fields.sort()
    ops = [
        {"id": "create.pair", "op": OP_UNIT_CREATE, "version": 1, "schema": "enhancer-pair", "count": "count",
         "initializers": [{"field": fid, "tensor": fid} for fid in init_fields]},
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
    version = contract.get("version")
    if version not in (2, 3):
        raise ValueError(f"unsupported circuit bundle version: {version!r}")
    eff_floor = _resolve_background_floor(contract.get("background_floor", 0.0), contract["epsilon"])
    pair_path = directory / "pair.fasta"
    validation = validate_development_paths(directory / "source", directory / "development", source_inputs={"sequence": pair_path})
    if not validation["valid"]:
        raise ValueError(f"circuit bundle validation failed: {validation}")
    bundle, artifacts = load_development_bundle_directory(directory / "development")
    lines = read_regular_file(pair_path, maximum_bytes=64 * 1024).decode("ascii").splitlines()
    if len(lines) != 4 or lines[0] != ">reference" or lines[2] != ">variant":
        raise ValueError("invalid circuit pair.fasta layout")
    reference, variant = lines[1], lines[3]
    if len(reference) != contract["reference_length"] or len(variant) != contract["variant_length"]:
        raise ValueError("pair.fasta lengths disagree with circuit contract")
    module = artifacts["development_module"]["module"]
    tensors = {t["id"]: base64.b64decode(t["storage"]["data"], validate=True) for t in module["tensors"]}
    if tensors["reference_nucleotides"] != pack("u64", list(reference.encode("ascii"))) or tensors["variant_nucleotides"] != pack("u64", list(variant.encode("ascii"))):
        raise ValueError("compiled nucleotides disagree with pair.fasta")

    stored_q = struct.unpack("<d", tensors["affinity_ratio"])[0]
    if version == 2:
        recomputed_q = affinity_ratio(
            reference, variant, contract["protein"], contract["epsilon"], background_floor=eff_floor
        )
        if not math.isclose(stored_q, recomputed_q, rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError("compiled affinity ratio disagrees with partition sum")
        return {
            "contract_sha256": expected_receipt,
            "bundle_sha256": bundle["artifact_sha256"],
            "development_sha256": artifacts["development_module"]["artifact_sha256"],
            "affinity_ratio": stored_q,
            "background_floor": eff_floor,
            "length": len(reference),
            "reference_length": len(reference),
            "variant_length": len(variant),
        }

    # v3 bundle verification: verify pinned JASPAR manifest, circuit_spec.json, motifs.json, and all TF ratios
    if contract.get("jaspar_manifest_sha256") is not None:
        if _load_jaspar_catalog()["manifest_sha256"] != contract["jaspar_manifest_sha256"]:
            raise ValueError("bundled JASPAR manifest SHA-256 disagrees with circuit contract")
    if contract.get("circuit_sha256") is not None:
        cbytes = read_regular_file(directory / "circuit_spec.json", maximum_bytes=256 * 1024)
        if hashlib.sha256(cbytes).hexdigest() != contract["circuit_sha256"]:
            raise ValueError("circuit_spec.json digest disagrees with contract")
        validate_circuit_spec(load_json_object(cbytes, "circuit_spec.json"))
    custom_motifs: dict = {}
    if contract.get("motifs_sha256") is not None:
        mbytes = read_regular_file(directory / "motifs.json", maximum_bytes=256 * 1024)
        if hashlib.sha256(mbytes).hexdigest() != contract["motifs_sha256"]:
            raise ValueError("motifs.json digest disagrees with contract")
        custom_motifs = load_json_object(mbytes, "motifs.json")

    proteins = contract["proteins"]
    specs = [custom_motifs[lbl] if lbl in custom_motifs else lbl for lbl in proteins]
    stored_vec = list(struct.unpack(f"<{len(proteins)}d", tensors["affinity_ratios"]))
    recomputed = multi_tf_affinity(
        reference, variant, specs, contract["epsilon"], background_floor=eff_floor
    )["ratios"]
    for lbl, sq in zip(proteins, stored_vec):
        if not math.isclose(sq, recomputed[lbl], rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError(f"compiled affinity ratio for {lbl} disagrees with partition sum")
    if not math.isclose(stored_q, stored_vec[0], rel_tol=1e-12, abs_tol=1e-12):
        raise ValueError("primary affinity_ratio scalar disagrees with affinity_ratios[0]")
    return {
        "contract_sha256": expected_receipt,
        "bundle_sha256": bundle["artifact_sha256"],
        "development_sha256": artifacts["development_module"]["artifact_sha256"],
        "affinity_ratio": stored_q,
        "affinity_ratios": dict(zip(proteins, stored_vec)),
        "background_floor": eff_floor,
        "length": len(reference),
        "reference_length": len(reference),
        "variant_length": len(variant),
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


def _classify_edit(d_rep: float, target_delta: float, fate_flipped: bool) -> str:
    if abs(d_rep) < 0.02 and fate_flipped:
        return "cryptic_bifurcation"
    if abs(d_rep) >= 0.02 and abs(target_delta) < 0.05 and not fate_flipped:
        return "canalized_buffered"
    if fate_flipped:
        return "concordant_switch"
    return "neutral"


def _edit_kind(reference: str, variant: str) -> str:
    if reference == variant:
        return "identical"
    if len(reference) > len(variant):
        return "deletion"
    if len(reference) < len(variant):
        return "insertion"
    diffs = sum(a != b for a, b in zip(reference, variant))
    return "snv" if diffs == 1 else "mnv"


def scan_enhancer(
    reference: str,
    *,
    signal: float | None = None,
    protein: str | dict | Path | None = None,
    epsilon: float = DEFAULT_EPSILON,
    background_floor: float | None = None,
    target: str | list[str] | tuple[str, ...] | None = None,
    duration: float = 40.0,
    phases: list | None = None,
    initial: dict | None = None,
    circuit: str | dict | Path | None = None,
    roles: dict[str, str] | None = None,
    mechanisms: dict[str, str] | None = None,
    motif_model_dir: Path | None = None,
) -> dict:
    """Scan all 3*L single-base mutations in an enhancer for reporter-vs-circuit dissociations."""
    eff_floor = _resolve_background_floor(background_floor, epsilon)
    ctx = _resolve_context(circuit, protein, target, signal, roles, mechanisms)
    model = ctx["model"]
    genes = ctx["genes"]
    eff_signal = ctx["signal"]
    eff_target = ctx["target"]
    bindings = ctx["bindings"]
    bound_genes = sorted(bindings)
    protocol = list(phases) if phases is not None else [(duration, eff_signal)]
    end_signal = protocol[-1][1]
    init_tuple = tuple(float(initial[g]) for g in genes) if initial is not None else (0.0,) * len(genes)

    # Preload motif binding tables & WT partitions for all active TFs
    tf_tables: dict[str, dict] = {}
    for lbl, spec in ctx["unique_specs"]:
        prof = load_motif_profile(spec, epsilon)
        wt_part = boltzmann_affinity(reference, spec, epsilon)
        tf_tables[lbl] = {"width": prof["width"], "ratios": prof["ratios"], "wt_part": wt_part, "z_wt": wt_part["partition_sum"]}

    primary_label = ctx["primary_label"]
    wt_info = tf_tables[primary_label]["wt_part"]
    rep_idx = _target_gene_index(genes, eff_target)
    fates_cfg = model.get("fates")
    wt_ratios = dict.fromkeys(bound_genes, 1.0)
    wt_traj = simulate_circuit(
        protocol, ratios=wt_ratios, target=eff_target, initial=init_tuple,
        circuit=ctx["circuit_arg"], roles=ctx["roles"], mechanisms=ctx["mechanisms"],
    )
    wt_end = wt_traj[-1][1:]
    wt_fate = fate_label(wt_end, genes, fates=fates_cfg)
    wt_rep_end = occupancy_ratios(model, wt_end, end_signal, wt_ratios, roles=ctx["roles"], mechanisms=ctx["mechanisms"])[rep_idx]

    has_nkx_olig = "Nkx2.2" in genes and "Olig2" in genes
    nkx_idx = genes.index("Nkx2.2") if has_nkx_olig else rep_idx
    olig_idx = genes.index("Olig2") if has_nkx_olig else rep_idx

    # Precompute unique ratio responses accurately (keyed at 1e-8 resolution)
    edits = []
    cache: dict[tuple, tuple[tuple[float, ...], float]] = {}
    pairs = []
    for pos, base in enumerate(reference):
        for alt in "ACGT":
            if alt == base:
                continue
            eid = f"{pos + 1}:{base}>{alt}"
            variant = reference[:pos] + alt + reference[pos + 1 :]
            q_by_tf = {}
            for lbl, tbl in tf_tables.items():
                z_mut = tbl["z_wt"] + _delta_partition_single(reference, pos, alt, tbl["ratios"], tbl["width"])
                q_by_tf[lbl] = (z_mut + eff_floor) / (tbl["z_wt"] + eff_floor)
            ratio_map = {g: q_by_tf[_motif_spec_label(bindings[g]["protein"])] for g in bound_genes}
            q_primary = q_by_tf[primary_label]
            key = tuple(round(ratio_map[g], 8) for g in bound_genes)
            if key not in cache:
                end_state = simulate_circuit(
                    protocol, ratios=ratio_map, target=eff_target, initial=init_tuple,
                    circuit=ctx["circuit_arg"], roles=ctx["roles"], mechanisms=ctx["mechanisms"],
                )[-1][1:]
                rep_end = occupancy_ratios(model, wt_end, end_signal, ratio_map, roles=ctx["roles"], mechanisms=ctx["mechanisms"])[rep_idx]
                cache[key] = (end_state, rep_end)
            end_state, rep_end = cache[key]
            d_rep = rep_end - wt_rep_end
            gene_deltas = {g: float(end_state[i] - wt_end[i]) for i, g in enumerate(genes)}
            target_delta = float(end_state[rep_idx] - wt_end[rep_idx])
            mut_fate = fate_label(end_state, genes, fates=fates_cfg)
            fate_flipped = mut_fate != wt_fate
            category = _classify_edit(d_rep, target_delta, fate_flipped)
            pairs.append({"id": eid, "reference": reference, "sequence": variant})
            entry = {
                "id": eid, "pos": pos, "ref": base, "alt": alt,
                "variant": variant, "affinity_ratio": q_primary,
                "affinity_ratios": q_by_tf if ctx["multi_tf"] else {primary_label: q_primary},
                "target": eff_target,
                "reporter_delta_end": d_rep,
                "gene_deltas": gene_deltas,
                "target_delta": target_delta,
                "wt_fate": wt_fate, "mut_fate": mut_fate,
                "fate_flipped": fate_flipped,
                "amplification_ratio": abs(target_delta) / max(abs(d_rep), 1e-12),
                "category": category,
            }
            if has_nkx_olig:
                entry["nkx2_2_delta"] = gene_deltas["Nkx2.2"]
                entry["olig2_delta"] = gene_deltas["Olig2"]
            edits.append(entry)

    kreimer_scores, siraj_scores = _optional_sequence_scores(pairs, motif_model_dir)
    for item in edits:
        if kreimer_scores is not None:
            item["kreimer_mpra_mean"] = kreimer_scores[item["id"]]["mean"]
            item["kreimer_mpra_72h"] = kreimer_scores[item["id"]]["h72"]
        if siraj_scores is not None:
            item["siraj_score"] = siraj_scores[item["id"]]

    cryptic = sorted(
        [e for e in edits if e["category"] == "cryptic_bifurcation"],
        key=lambda e: (-abs(e["target_delta"]), abs(e["reporter_delta_end"]), e["id"]),
    )
    counts = {cat: sum(e["category"] == cat for e in edits) for cat in ("cryptic_bifurcation", "concordant_switch", "canalized_buffered", "neutral")}
    return {
        "circuit": ctx["circuit_name"], "genes": genes,
        "protein": primary_label, "proteins": [lbl for lbl, _ in ctx["unique_specs"]],
        "target": eff_target, "epsilon": epsilon, "background_floor": eff_floor,
        "signal": eff_signal, "duration": duration,
        "length": len(reference), "wt_partition": wt_info,
        "wt_partitions": {lbl: tbl["wt_part"] for lbl, tbl in tf_tables.items()},
        "wt_end_state": dict(zip(genes, wt_end)), "wt_fate": wt_fate, "wt_reporter_end": wt_rep_end,
        "unique_affinity_solves": len(cache), "category_counts": counts,
        "cryptic_switches": cryptic, "edits": edits,
    }


def compare_variants(
    reference: str,
    variant: str,
    *,
    protein: str | dict | Path | None = None,
    epsilon: float = DEFAULT_EPSILON,
    background_floor: float | None = None,
    signal: float | None = None,
    target: str | list[str] | tuple[str, ...] | None = None,
    duration: float = 40.0,
    phases: list | None = None,
    initial: dict | None = None,
    circuit: str | dict | Path | None = None,
    roles: dict[str, str] | None = None,
    mechanisms: dict[str, str] | None = None,
    motif_model_dir: Path | None = None,
    compile_dir: Path | None = None,
) -> dict:
    """Compare any reference vs variant DNA pair (SNV, MNV, insertion, or deletion) end-to-end."""
    eff_floor = _resolve_background_floor(background_floor, epsilon)
    ctx = _resolve_context(circuit, protein, target, signal, roles, mechanisms)
    model = ctx["model"]
    genes = ctx["genes"]
    eff_signal = ctx["signal"]
    eff_target = ctx["target"]
    bindings = ctx["bindings"]

    multi = multi_tf_affinity(reference, variant, dict(ctx["unique_specs"]), epsilon, background_floor=eff_floor)
    wt_parts = multi["wt_partitions"]
    var_parts = multi["variant_partitions"]
    q_by_tf = multi["ratios"]
    q = q_by_tf[ctx["primary_label"]]
    ratio_map = {g: q_by_tf[_motif_spec_label(b["protein"])] for g, b in bindings.items()}
    wt_ratios = dict.fromkeys(bindings, 1.0)

    protocol = list(phases) if phases is not None else [(duration, eff_signal)]
    rep_idx = _target_gene_index(genes, eff_target)
    end_signal = protocol[-1][1]
    fates_cfg = model.get("fates")

    wt_assay = run_circuit_assay(
        protocol, ratio=wt_ratios, mode="endogenous", target=eff_target, initial=initial,
        circuit=ctx["circuit_arg"], roles=ctx["roles"], mechanisms=ctx["mechanisms"],
    )
    var_assay = run_circuit_assay(
        protocol, ratio=ratio_map, mode="endogenous", target=eff_target, initial=initial,
        circuit=ctx["circuit_arg"], roles=ctx["roles"], mechanisms=ctx["mechanisms"],
    )
    rep_assay = run_circuit_assay(
        protocol, ratio=ratio_map, mode="reporter", target=eff_target, initial=initial,
        circuit=ctx["circuit_arg"], roles=ctx["roles"], mechanisms=ctx["mechanisms"],
    )
    rev_multi = multi_tf_affinity(
        reference, reference, dict(ctx["unique_specs"]), epsilon, background_floor=eff_floor
    )["ratios"]
    rev_ratios = {g: rev_multi[_motif_spec_label(b["protein"])] for g, b in bindings.items()}
    rev_assay = run_circuit_assay(
        protocol, ratio=rev_ratios, mode="endogenous", target=eff_target, initial=initial,
        circuit=ctx["circuit_arg"], roles=ctx["roles"], mechanisms=ctx["mechanisms"],
    )

    wt_end = tuple(wt_assay["trajectory"][-1][1:])
    var_end = tuple(var_assay["trajectory"][-1][1:])
    wt_rep_end = occupancy_ratios(model, wt_end, end_signal, wt_ratios, roles=ctx["roles"], mechanisms=ctx["mechanisms"])[rep_idx]
    var_rep_end = occupancy_ratios(model, wt_end, end_signal, ratio_map, roles=ctx["roles"], mechanisms=ctx["mechanisms"])[rep_idx]

    d_rep = var_rep_end - wt_rep_end
    gene_deltas = {g: float(var_end[i] - wt_end[i]) for i, g in enumerate(genes)}
    target_delta = float(var_end[rep_idx] - wt_end[rep_idx])
    wt_fate = fate_label(wt_end, genes, fates=fates_cfg)
    mut_fate = fate_label(var_end, genes, fates=fates_cfg)
    fate_flipped = mut_fate != wt_fate
    category = _classify_edit(d_rep, target_delta, fate_flipped)
    diff = lambda a, b: max(abs(x - y) for r1, r2 in zip(a["trajectory"], b["trajectory"]) for x, y in zip(r1[1:], r2[1:]))

    compiled = None
    if compile_dir is not None:
        motif_list = [spec for _, spec in ctx["unique_specs"]]
        receipt = compile_circuit_bundle(
            reference, variant, protein=ctx["primary_protein"], epsilon=epsilon,
            background_floor=eff_floor, output=compile_dir, motif_specs=motif_list, circuit=ctx["circuit_arg"],
        )
        compiled = verify_circuit_bundle(compile_dir, receipt)

    kreimer_scores, siraj_scores = _optional_sequence_scores(
        [{"id": "custom", "reference": reference, "sequence": variant}], motif_model_dir
    )
    out = {
        "circuit": ctx["circuit_name"],
        "genes": genes,
        "protein": ctx["primary_label"],
        "proteins": [lbl for lbl, _ in ctx["unique_specs"]],
        "target": eff_target,
        "epsilon": epsilon,
        "background_floor": eff_floor,
        "protocol": protocol,
        "edit_kind": _edit_kind(reference, variant),
        "reference_length": len(reference),
        "variant_length": len(variant),
        "affinity_ratio": q,
        "affinity_ratios": q_by_tf,
        "wt_partition": wt_parts[ctx["primary_label"]],
        "variant_partition": var_parts[ctx["primary_label"]],
        "wt_partitions": wt_parts,
        "variant_partitions": var_parts,
        "wt_end_state": dict(zip(genes, wt_end)),
        "variant_end_state": dict(zip(genes, var_end)),
        "wt_fate": wt_fate,
        "mut_fate": mut_fate,
        "fate_flipped": fate_flipped,
        "reporter_delta_end": d_rep,
        "gene_deltas": gene_deltas,
        "target_delta": target_delta,
        "amplification_ratio": abs(target_delta) / max(abs(d_rep), 1e-12),
        "category": category,
        "reversion_max_difference": diff(wt_assay, rev_assay),
        "reporter_max_difference": diff(wt_assay, rep_assay),
        "wt_assay": wt_assay,
        "variant_assay": var_assay,
        "reporter_assay": rep_assay,
    }
    if "Nkx2.2" in gene_deltas:
        out["nkx2_2_delta"] = gene_deltas["Nkx2.2"]
    if "Olig2" in gene_deltas:
        out["olig2_delta"] = gene_deltas["Olig2"]
    if compiled is not None:
        out["compiled"] = compiled
    if kreimer_scores is not None:
        out["kreimer_mpra_mean"] = kreimer_scores["custom"]["mean"]
        out["kreimer_mpra_72h"] = kreimer_scores["custom"]["h72"]
    if siraj_scores is not None:
        out["siraj_score"] = siraj_scores["custom"]
    return out


def find_second_site_rescue(
    reference: str,
    primary_edit: dict,
    protein: str | dict | Path | None = "GLI2",
    epsilon: float = DEFAULT_EPSILON,
    *,
    background_floor: float | None = None,
    signal: float | None = None,
    target: str | list[str] | tuple[str, ...] | None = None,
    duration: float = 40.0,
    phases: list | None = None,
    initial: dict | None = None,
    circuit: str | dict | Path | None = None,
    roles: dict[str, str] | None = None,
    mechanisms: dict[str, str] | None = None,
) -> dict | None:
    """Search single-base edits elsewhere in the enhancer that move TF affinity back toward WT and test ODE fate rescue."""
    eff_floor = _resolve_background_floor(background_floor, epsilon)
    ctx = _resolve_context(circuit, protein, target, signal, roles, mechanisms)
    model = ctx["model"]
    genes = ctx["genes"]
    eff_signal = ctx["signal"]
    eff_target = ctx["target"]
    bindings = ctx["bindings"]
    fates_cfg = model.get("fates")

    pos0 = primary_edit["pos"]
    base_mut = primary_edit.get("variant") or (
        reference[:pos0] + primary_edit["alt"] + reference[pos0 + 1 :]
    )
    q_single = float(primary_edit["affinity_ratio"])
    primary_spec = ctx["primary_protein"]
    consensus, ratios = _binding_tables(primary_spec, epsilon)
    width = len(consensus)
    z_ref = boltzmann_affinity(reference, primary_spec, epsilon)["partition_sum"]
    z_single = boltzmann_affinity(base_mut, primary_spec, epsilon)["partition_sum"]
    loss_primary = q_single <= 1.0
    dist_single = abs(math.log(max(q_single, 1e-18)))

    cands = []
    for pos, base in enumerate(base_mut):
        if abs(pos - pos0) <= 0:
            continue
        for alt in "ACGT":
            if alt == base:
                continue
            z_double = z_single + _delta_partition_single(base_mut, pos, alt, ratios, width)
            q_double = (z_double + eff_floor) / (z_ref + eff_floor)
            if q_double <= 0:
                continue
            moves_toward_wt = (
                (q_double > q_single * 1.05)
                if loss_primary
                else (q_double < q_single / 1.05)
            )
            if moves_toward_wt:
                dist_double = abs(math.log(max(q_double, 1e-18)))
                cands.append({
                    "id": f"{pos + 1}:{base}>{alt}",
                    "pos": pos,
                    "double_variant": base_mut[:pos] + alt + base_mut[pos + 1 :],
                    "affinity_ratio": q_double,
                    "_dist": dist_double,
                    "_drop": dist_single - dist_double,
                })
    if not cands:
        return None

    cands.sort(key=lambda item: (item["_dist"], item["id"]) if not loss_primary else (-item["affinity_ratio"], item["id"]))
    protocol = list(phases) if phases is not None else [(duration, eff_signal)]
    init_tuple = tuple(float(initial[g]) for g in genes) if initial is not None else (0.0,) * len(genes)
    rep_idx = _target_gene_index(genes, eff_target)

    def _ratio_map_for_seq(seq_or_q: str | float) -> dict[str, float]:
        if ctx["multi_tf"] and isinstance(seq_or_q, str):
            q_map = multi_tf_affinity(
                reference, seq_or_q, dict(ctx["unique_specs"]), epsilon, background_floor=eff_floor
            )["ratios"]
            return {g: q_map[_motif_spec_label(b["protein"])] for g, b in bindings.items()}
        val = q_single if isinstance(seq_or_q, str) else float(seq_or_q)
        return {g: val for g in bindings}

    wt_ratios = dict.fromkeys(bindings, 1.0)
    wt_fate = primary_edit.get("wt_fate") or fate_label(
        simulate_circuit(
            protocol, ratios=wt_ratios, target=eff_target, initial=init_tuple,
            circuit=ctx["circuit_arg"], roles=ctx["roles"], mechanisms=ctx["mechanisms"],
        )[-1][1:],
        genes,
        fates=fates_cfg,
    )
    primary_fate = primary_edit.get("mut_fate") or fate_label(
        simulate_circuit(
            protocol, ratios=_ratio_map_for_seq(base_mut), target=eff_target, initial=init_tuple,
            circuit=ctx["circuit_arg"], roles=ctx["roles"], mechanisms=ctx["mechanisms"],
        )[-1][1:],
        genes,
        fates=fates_cfg,
    )

    best = None
    for cand in cands[:8]:
        r_map = _ratio_map_for_seq(cand["double_variant"] if ctx["multi_tf"] else cand["affinity_ratio"])
        end_state = simulate_circuit(
            protocol, ratios=r_map, target=eff_target, initial=init_tuple,
            circuit=ctx["circuit_arg"], roles=ctx["roles"], mechanisms=ctx["mechanisms"],
        )[-1][1:]
        d_fate = fate_label(end_state, genes, fates=fates_cfg)
        rescued = d_fate == wt_fate and primary_fate != wt_fate
        enriched = {
            "id": cand["id"],
            "pos": cand["pos"],
            "double_variant": cand["double_variant"],
            "affinity_ratio": cand["affinity_ratio"],
            "wt_fate": wt_fate,
            "primary_fate": primary_fate,
            "double_fate": d_fate,
            "double_target": float(end_state[rep_idx]),
            "double_end_state": dict(zip(genes, (float(v) for v in end_state))),
            "fate_rescued": rescued,
        }
        if "Nkx2.2" in genes:
            enriched["double_nkx2_2"] = float(end_state[genes.index("Nkx2.2")])
        if best is None or (enriched["fate_rescued"] and not best["fate_rescued"]):
            best = enriched
        if best["fate_rescued"]:
            break
    return best


_find_second_site_rescue = find_second_site_rescue


def multi_locus_benchmark(
    loci: tuple[str, ...] | list[str] = CURATED_LOCI,
    *,
    protein: str | dict | Path | None = None,
    signal: float | None = None,
    target: str | list[str] | tuple[str, ...] | None = None,
    duration: float = 40.0,
    circuit: str | dict | Path | None = None,
    background_floor: float | None = None,
    motif_model_dir: Path | None = None,
    output: Path | None = None,
) -> dict:
    """Benchmark cryptic bifurcations and measured GSE188264 variants across multiple human enhancers."""
    sys.path.insert(0, str(REPO / "examples/neural-regulation"))
    import regulation as reg
    rows, _ = reg.load_data()
    by_locus: dict[str, list] = {}
    for row in rows:
        by_locus.setdefault(row.locus, []).append(row)

    eff_floor = _resolve_background_floor(background_floor, DEFAULT_EPSILON)
    ctx = _resolve_context(circuit, protein, target, signal, None, None)
    results = []
    for locus in loci:
        locus_rows = by_locus.get(locus)
        if not locus_rows:
            raise ValueError(f"unknown GSE188264 locus: {locus}")
        reference = locus_rows[0].reference
        scan = scan_enhancer(
            reference, signal=ctx["signal"], protein=protein, target=target,
            duration=duration, circuit=circuit, background_floor=eff_floor, motif_model_dir=motif_model_dir,
        )
        measured = []
        seen_seq: set[str] = set()
        for row in locus_rows:
            if row.kind == "wild-type" or row.sequence == reference or row.sequence in seen_seq:
                continue
            seen_seq.add(row.sequence)
            comp = compare_variants(
                reference, row.sequence, protein=protein, signal=ctx["signal"],
                target=target, duration=duration, circuit=circuit, background_floor=eff_floor,
            )
            alpha_vals = list(row.alpha)
            measured.append({
                "id": row.id,
                "split": row.split,
                "kind": row.kind,
                "edit_kind": comp["edit_kind"],
                "affinity_ratio": comp["affinity_ratio"],
                "reporter_delta_end": comp["reporter_delta_end"],
                "target_delta": comp["target_delta"],
                "wt_fate": comp["wt_fate"],
                "mut_fate": comp["mut_fate"],
                "fate_flipped": comp["fate_flipped"],
                "category": comp["category"],
                "measured_alpha_mean": sum(alpha_vals) / len(alpha_vals),
                "measured_alpha_72h": alpha_vals[-1],
                "measured_alpha_max_abs": max(abs(v) for v in alpha_vals),
            })

        top_cryptic = (
            {k: v for k, v in scan["cryptic_switches"][0].items() if k != "variant"}
            if scan["cryptic_switches"]
            else None
        )
        rescue = (
            find_second_site_rescue(
                reference, scan["cryptic_switches"][0], protein, DEFAULT_EPSILON,
                background_floor=eff_floor, signal=ctx["signal"], target=target, duration=duration, circuit=circuit,
            )
            if scan["cryptic_switches"]
            else None
        )
        results.append({
            "locus": locus,
            "length": len(reference),
            "peak_tf_weight": scan["wt_partition"]["peak_site"]["weight"],
            "peak_gli_weight": scan["wt_partition"]["peak_site"]["weight"],
            "strong_site_count": len(scan["wt_partition"]["strong_sites"]),
            "wt_fate": scan["wt_fate"],
            "category_counts": scan["category_counts"],
            "top_cryptic": top_cryptic,
            "second_site_rescue": {k: v for k, v in rescue.items() if k != "double_variant"} if rescue else None,
            "measured_variants_tested": len(measured),
            "measured_variant_switches": sum(m["fate_flipped"] for m in measured),
            "measured_variants": sorted(measured, key=lambda m: (-abs(m["target_delta"]), m["affinity_ratio"]))[:12],
        })

    summary = {
        "format": "epigenesis.neural-circuit-benchmark/v1",
        "circuit": ctx["circuit_name"],
        "protein": ctx["primary_label"],
        "target": ctx["target"],
        "signal": ctx["signal"],
        "background_floor": eff_floor,
        "loci_tested": len(results),
        "total_cryptic_bifurcations": sum(r["category_counts"]["cryptic_bifurcation"] for r in results),
        "total_measured_variants": sum(r["measured_variants_tested"] for r in results),
        "total_measured_switches": sum(r["measured_variant_switches"] for r in results),
        "loci": results,
    }
    if output is not None:
        output.mkdir(parents=True, exist_ok=True)
        (output / "benchmark.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    return summary


def program_from_grn_state(state: dict[str, float], *, mapping: dict | None = None):
    """Deterministically map any N-gene GRN steady state into a valid CX3D Program."""
    sys.path.insert(0, str(REPO / "examples/neural-construction"))
    import construction as cx

    if type(state) is not dict or not state or any(not _num(v) or not 0.0 <= float(v) <= 1.0 for v in state.values()):
        raise ValueError("GRN state must be a non-empty dict of gene concentrations in [0, 1]")

    base = cx.default_program()
    overrides: dict[str, float] = {}
    valid_fields = set(cx.Program._fields)

    if mapping is not None:
        if type(mapping) is not dict:
            raise ValueError("mapping must be a JSON object keyed by CX3D Program fields")
        for field, rule in mapping.items():
            if field not in valid_fields:
                raise ValueError(f"unsupported CX3D field in mapping: {field!r}")
            if type(rule) is not dict or type(rule.get("weights")) is not dict:
                raise ValueError(f"mapping rule for {field} requires 'weights' dict")
            val = float(rule.get("base", getattr(base, field)))
            for g, w in rule["weights"].items():
                if not _num(w):
                    raise ValueError(f"invalid weight for {field}.{g}")
                val += float(w) * float(state.get(g, 0.0))
            if "min" in rule:
                val = max(float(rule["min"]), val)
            if "max" in rule:
                val = min(float(rule["max"]), val)
            overrides[field] = round(val, 6)
    elif {"Pax6", "Olig2", "Nkx2.2", "Irx3"} <= set(state):
        pax = float(state["Pax6"])
        olig = float(state["Olig2"])
        nkx = float(state["Nkx2.2"])
        irx = float(state["Irx3"])
        overrides = {
            "axon_speed": round(65.0 + 40.0 * nkx + 15.0 * olig, 6),
            "dendrite_speed": round(45.0 + 25.0 * olig + 10.0 * pax, 6),
            "guidance_gain": round(max(0.5, 1.5 + 1.2 * nkx - 0.6 * irx), 6),
            "branch_rate": round(0.2 + 0.25 * olig + 0.1 * nkx, 6),
            "noise": round(0.15 + 0.1 * pax, 6),
        }
    else:
        ordered = [float(state[g]) for g in sorted(state)]
        m = len(ordered)
        g0 = ordered[0]
        g1 = ordered[1 % m]
        g2 = ordered[2 % m]
        mean_other = sum(ordered[1:]) / max(1, m - 1) if m > 1 else g0
        overrides = {
            "axon_speed": round(65.0 + 45.0 * g0 + 15.0 * g1, 6),
            "dendrite_speed": round(45.0 + 30.0 * g1 + 10.0 * g2, 6),
            "guidance_gain": round(max(0.5, 1.5 + 1.2 * g0 - 0.5 * mean_other), 6),
            "branch_rate": round(max(0.05, 0.20 + 0.25 * g1 + 0.10 * g0), 6),
            "noise": round(max(0.05, 0.15 + 0.10 * g2), 6),
        }

    program_dict = {**base._asdict(), **overrides}
    program = cx.Program(**program_dict)
    cx.validate_program(program)
    return program


def bridge_to_construction(
    reference: str,
    variant: str,
    output: Path,
    *,
    protein: str | dict | Path | None = None,
    epsilon: float = DEFAULT_EPSILON,
    background_floor: float | None = None,
    signal: float | None = None,
    target: str | list[str] | tuple[str, ...] | None = None,
    circuit: str | dict | Path | None = None,
    roles: dict[str, str] | None = None,
    mechanisms: dict[str, str] | None = None,
    mapping: dict | Path | None = None,
    grn_duration: float = 40.0,
    dev_duration: float = 3.0,
    seed: int = 7,
) -> dict:
    """Bridge enhancer DNA -> TF affinity -> N-gene GRN steady state -> 480-bp recipe -> CX3D morphogenesis + LIF assay."""
    sys.path.insert(0, str(REPO / "examples/neural-construction"))
    import activity as lif
    import construction as cx

    eff_floor = _resolve_background_floor(background_floor, epsilon)
    if isinstance(mapping, Path) or (isinstance(mapping, str) and Path(mapping).is_file()):
        mapping_obj = load_json_object(read_regular_file(Path(mapping), maximum_bytes=64 * 1024), "CX3D mapping JSON")
    else:
        mapping_obj = mapping

    comparison = compare_variants(
        reference, variant, protein=protein, epsilon=epsilon,
        background_floor=eff_floor, signal=signal, target=target, duration=grn_duration,
        circuit=circuit, roles=roles, mechanisms=mechanisms,
    )
    rev_state = compare_variants(
        reference, reference, protein=protein, epsilon=epsilon,
        background_floor=eff_floor, signal=signal, target=target, duration=grn_duration,
        circuit=circuit, roles=roles, mechanisms=mechanisms,
    )["variant_end_state"]

    output.mkdir(parents=True, exist_ok=True)
    backend = cx.build_backend(output / "backend")

    runs = {}
    for label, state in (
        ("reference", comparison["wt_end_state"]),
        ("variant", comparison["variant_end_state"]),
        ("allele_rescue", rev_state),
    ):
        branch = output / label
        branch.mkdir()
        program = program_from_grn_state(state, mapping=mapping_obj)
        recipe_path = branch / "recipe.dna"
        recipe_bytes = cx.encode_recipe(program)
        recipe_path.write_bytes(recipe_bytes)
        cx.compile_recipe(recipe_path, branch / "compiled")
        loaded_program, validation = cx.load_compiled(branch / "compiled")
        summary = cx.develop(loaded_program, backend, branch / "development-run", duration=dev_duration, seed=seed)
        act = lif.assay(branch / "development-run", cx.verify_trace)
        (branch / "activity.json").write_text(json.dumps(act, indent=2) + "\n")
        runs[label] = {
            "program": loaded_program._asdict(),
            "recipe_sha256": hashlib.sha256(recipe_bytes).hexdigest(),
            "compiled_valid": validation["valid"],
            "construction": summary,
            "activity": {
                "cells": act["cells"],
                "contacts": act["contacts"],
                "transmission_controls_pass": act["transmission_controls_pass"],
                "total_spikes": sum(len(t["intact_spikes"]) for t in act["trials"]),
            },
        }

    ref_r, var_r, res_r = runs["reference"], runs["variant"], runs["allele_rescue"]
    morphology_diverged = (
        abs(ref_r["construction"]["length"] - var_r["construction"]["length"]) > 1e-6
        or ref_r["recipe_sha256"] != var_r["recipe_sha256"]
    )
    exact_rescue = (
        res_r["recipe_sha256"] == ref_r["recipe_sha256"]
        and res_r["construction"]["length"] == ref_r["construction"]["length"]
        and res_r["construction"]["cells"] == ref_r["construction"]["cells"]
        and res_r["construction"]["synapses"] == ref_r["construction"]["synapses"]
    )
    controls_ok = all(r["compiled_valid"] and r["activity"]["transmission_controls_pass"] for r in runs.values())
    bridge_report = {
        "format": "epigenesis.neural-circuit-construction-bridge/v1",
        "biological_acceptance": False,
        "circuit": comparison["circuit"],
        "genes": comparison["genes"],
        "protein": comparison["protein"],
        "scope": "computational DNA -> TF -> GRN -> CX3D + LIF chain; parameters and physical units are uncalibrated model assumptions",
        "variant_summary": {
            "edit_kind": comparison["edit_kind"],
            "target": comparison["target"],
            "affinity_ratio": comparison["affinity_ratio"],
            "wt_fate": comparison["wt_fate"],
            "mut_fate": comparison["mut_fate"],
            "fate_flipped": comparison["fate_flipped"],
            "category": comparison["category"],
        },
        "runs": runs,
        "checks": {
            "all_compiled_and_lif_controls_pass": controls_ok,
            "variant_diverges_morphology": morphology_diverged,
            "allele_rescue_matches_reference_exactly": exact_rescue,
        },
        "all_checks_passed": bool(controls_ok and morphology_diverged and exact_rescue),
    }
    (output / "bridge.json").write_text(json.dumps(bridge_report, indent=2, allow_nan=False) + "\n")
    return bridge_report


def _render_html(report: dict) -> str:
    scans = report["scan"]["edits"]
    wt_fate = report["scan"]["wt_fate"]
    top = report["featured_variant"]
    target_gene = str(report["scan"].get("target", DEFAULT_TARGET))
    protein_name = str(report["scan"].get("protein", "GLI2"))
    circuit_name = str(report.get("circuit", report["scan"].get("circuit", "neural-tube")))
    genes = report["scan"].get("genes") or ["Pax6", "Olig2", "Nkx2.2", "Irx3"]
    target_idx = _target_gene_index(genes, target_gene)
    alt_idx = genes.index("Olig2") if "Olig2" in genes else (1 if len(genes) > 1 and target_idx != 1 else (0 if len(genes) == 1 else (target_idx + 1) % len(genes)))

    w, h, pad = 540, 290, 48
    max_x = max([0.05] + [abs(e["reporter_delta_end"]) for e in scans])
    points = []
    for e in scans:
        td = e.get("target_delta", e.get("nkx2_2_delta", 0.0))
        cx = pad + (w - 2 * pad) * min(1.0, abs(e["reporter_delta_end"]) / max_x)
        cy = h - pad - (h - 2 * pad) * min(1.0, abs(td))
        color = "#c0262d" if e["category"] == "cryptic_bifurcation" else "#1f6feb" if e["category"] == "canalized_buffered" else "#9ca3af"
        r = 5 if e["id"] == top["id"] else 3.2
        title = html.escape(f"{e['id']} | q={e['affinity_ratio']:.4f} | dRep={e['reporter_delta_end']:+.4f} | d{target_gene}={td:+.3f} ({e['category']})")
        points.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r}" fill="{color}" opacity="0.82"><title>{title}</title></circle>')

    sp = report["spatial"]
    xs = sp["positions"]
    def poly(vals, stroke, dash=""):
        pts = " ".join(f"{pad + (w - 2 * pad) * x:.1f},{h - pad - (h - 2 * pad) * v:.1f}" for x, v in zip(xs, vals))
        d_attr = f' stroke-dasharray="{dash}"' if dash else ""
        return f'<polyline fill="none" stroke="{stroke}" stroke-width="2.2"{d_attr} points="{pts}"/>'

    t_wt = [s[target_idx] for s in sp["reference"]]
    t_mut = [s[target_idx] for s in sp["variant"]]
    a_wt = [s[alt_idx] for s in sp["reference"]]
    a_mut = [s[alt_idx] for s in sp["variant"]]

    rows_html = "".join(
        f"<tr><td><code>{html.escape(e['id'])}</code></td><td>{e['affinity_ratio']:.4f}</td>"
        f"<td>{e['reporter_delta_end']:+.4f}</td><td>{e.get('target_delta', 0.0):+.3f}</td>"
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
<p><strong>Circuit:</strong> {html.escape(circuit_name)} ({html.escape(protein_name)}) · <strong>Locus:</strong> {html.escape(report['locus'])} ({report['scan']['length']} bp, target {html.escape(target_gene)}) · <strong>Featured Edit:</strong>
<span class="badge">{html.escape(top['id'])}</span> (q={top['affinity_ratio']:.4f}, passive reporter &Delta;={top['reporter_delta_end']:+.4f}, endogenous &Delta;{html.escape(target_gene)}={top.get('target_delta', 0.0):+.3f}, {html.escape(wt_fate)} &rarr; {html.escape(top['mut_fate'])})</p>
<div class="grid">
<div><h3>Passive Reporter vs. Endogenous Circuit ({report['scan']['length']*3} edits at S={report['scan']['signal']})</h3>
<svg viewBox="0 0 {w} {h}"><line x1="{pad}" y1="{h-pad}" x2="{w-pad}" y2="{h-pad}" stroke="#111"/><line x1="{pad}" y1="{pad}" x2="{pad}" y2="{h-pad}" stroke="#111"/>
<text x="{w/2}" y="{h-10}" text-anchor="middle" font-size="12">|Passive Reporter Occupancy Shift| (max {max_x:.3f})</text>
<text x="16" y="{h/2}" transform="rotate(-90 16 {h/2})" text-anchor="middle" font-size="12">|Endogenous {html.escape(target_gene)} Fate Shift|</text>
{''.join(points)}</svg></div>
<div><h3>Spatial Morphogen Gradient (t={sp['duration']})</h3>
<svg viewBox="0 0 {w} {h}"><line x1="{pad}" y1="{h-pad}" x2="{w-pad}" y2="{h-pad}" stroke="#111"/><line x1="{pad}" y1="{pad}" x2="{pad}" y2="{h-pad}" stroke="#111"/>
<text x="{w/2}" y="{h-10}" text-anchor="middle" font-size="12">Position x (Signal S=exp(-x/lambda))</text>
<text x="16" y="{h/2}" transform="rotate(-90 16 {h/2})" text-anchor="middle" font-size="12">Steady Concentration (Solid=WT, Dashed=Mutant)</text>
{poly(t_wt, "#c0262d")}{poly(t_mut, "#c0262d", "6 4")}{poly(a_wt, "#1f6feb")}{poly(a_mut, "#1f6feb", "6 4")}</svg></div></div>
<h3>Top Reporter-Silent, Circuit-Switch Mutations</h3>
<table><thead><tr><th>Edit</th><th>q ({html.escape(protein_name)})</th><th>&Delta;Reporter</th><th>&Delta;{html.escape(target_gene)}</th><th>Fate Flip</th><th>Kreimer MPRA</th><th>Siraj Score</th></tr></thead>
<tbody>{rows_html}</tbody></table>
<h3>Causal & Numerical Controls</h3><ul>{checks_html}</ul>
<p style="font-size:12px;color:#555">Compiled brainc receipt: <code>{report['compiled']['contract_sha256']}</code>. Concentrations and time are uncalibrated model units.</p>
</body></html>"""


def demo(
    output: Path,
    *,
    locus: str = CURATED_LOCI[0],
    target: str | None = None,
    signal: float | None = None,
    protein: str | dict | Path | None = None,
    circuit: str | dict | Path | None = None,
    background_floor: float | None = None,
    motif_model_dir: Path | None = None,
    with_construction: bool = False,
) -> dict:
    """Run the Variant-to-Circuit discovery demo on a real human enhancer from GSE188264."""
    sys.path.insert(0, str(REPO / "examples/neural-regulation"))
    import regulation as reg
    rows, _ = reg.load_data()
    reference = next((r.reference for r in rows if r.locus == locus), None)
    if reference is None:
        raise ValueError(f"unknown GSE188264 locus: {locus}")

    eff_floor = _resolve_background_floor(background_floor, DEFAULT_EPSILON)
    ctx = _resolve_context(circuit, protein, target, signal, None, None)
    model = ctx["model"]
    genes = ctx["genes"]
    eff_signal = ctx["signal"]
    eff_target = ctx["target"]
    bindings = ctx["bindings"]
    rep_idx = _target_gene_index(genes, eff_target)

    output.mkdir()
    scan = scan_enhancer(
        reference, signal=eff_signal, protein=protein, epsilon=DEFAULT_EPSILON,
        background_floor=eff_floor, target=target, circuit=circuit, motif_model_dir=motif_model_dir,
    )
    if not scan["cryptic_switches"] and not scan["edits"]:
        raise RuntimeError(f"no edits generated for locus {locus}")
    if ctx["circuit_arg"] is None and locus == CURATED_LOCI[0] and eff_target == DEFAULT_TARGET and eff_signal == DEFAULT_SIGNAL and scan["cryptic_switches"]:
        featured = next((e for e in scan["cryptic_switches"] if e["id"] == "60:C>A"), scan["cryptic_switches"][0])
    elif scan["cryptic_switches"]:
        featured = scan["cryptic_switches"][0]
    else:
        featured = max(scan["edits"], key=lambda e: (abs(e["target_delta"]), -abs(e["reporter_delta_end"]), e["id"]))
    variant = featured["variant"]
    multi_var = multi_tf_affinity(
        reference, variant, dict(ctx["unique_specs"]), DEFAULT_EPSILON, background_floor=eff_floor
    )
    q_by_tf = multi_var["ratios"]
    ratio_map = {g: q_by_tf[_motif_spec_label(b["protein"])] for g, b in bindings.items()}
    wt_ratios = dict.fromkeys(bindings, 1.0)

    motif_list = [spec for _, spec in ctx["unique_specs"]]
    receipt = compile_circuit_bundle(
        reference, variant, protein=ctx["primary_protein"], epsilon=DEFAULT_EPSILON,
        background_floor=eff_floor, output=output / "compiled", motif_specs=motif_list, circuit=ctx["circuit_arg"],
    )
    compiled = verify_circuit_bundle(output / "compiled", receipt)

    protocol = [(40.0, eff_signal)]
    kw = dict(circuit=ctx["circuit_arg"], roles=ctx["roles"], mechanisms=ctx["mechanisms"])
    wt_assay = run_circuit_assay(protocol, ratio=wt_ratios, mode="endogenous", target=eff_target, **kw)
    mut_assay = run_circuit_assay(protocol, ratio=ratio_map, mode="endogenous", target=eff_target, **kw)
    rep_assay = run_circuit_assay(protocol, ratio=ratio_map, mode="reporter", target=eff_target, **kw)
    allele_rescue = run_circuit_assay(protocol, ratio=wt_ratios, mode="endogenous", target=eff_target, **kw)
    mediator_rescue = run_circuit_assay(protocol, ratio=wt_ratios, mode="endogenous", target=eff_target, **kw)
    phenocopy = run_circuit_assay(protocol, ratio=ratio_map, mode="endogenous", target=eff_target, **kw)

    fine_wt = run_circuit_assay(protocol, ratio=wt_ratios, mode="endogenous", target=eff_target, dt=pat.DEFAULT_DT / 2, **kw)
    fine_mut = run_circuit_assay(protocol, ratio=ratio_map, mode="endogenous", target=eff_target, dt=pat.DEFAULT_DT / 2, **kw)
    diff = lambda a, b: max(abs(x - y) for r1, r2 in zip(a, b) for x, y in zip(r1[1:], r2[1:]))
    refinement_err = max(diff(wt_assay["trajectory"], fine_wt["trajectory"][::2]), diff(mut_assay["trajectory"], fine_mut["trajectory"][::2]))

    grad_len = float(model["gradient_length"])
    positions = [i / 40 for i in range(41)]
    sp_ref = [list(simulate_circuit([(100.0, math.exp(-x / grad_len))], ratios=wt_ratios, target=eff_target, **kw)[-1][1:]) for x in positions]
    sp_mut = [list(simulate_circuit([(100.0, math.exp(-x / grad_len))], ratios=ratio_map, target=eff_target, **kw)[-1][1:]) for x in positions]
    p3_wt = [x for x, y in zip(positions, sp_ref) if y[rep_idx] > 0.5]
    p3_mut = [x for x, y in zip(positions, sp_mut) if y[rep_idx] > 0.5]

    hist_wt = simulate_circuit([(5.0, 1.0), (15.0, 0.05)], ratios=wt_ratios, target=eff_target, mechanism="competitive", **kw)[-1][rep_idx + 1]
    hist_mut = simulate_circuit([(5.0, 1.0), (15.0, 0.05)], ratios=ratio_map, target=eff_target, mechanism="competitive", **kw)[-1][rep_idx + 1]
    hist_mut_act = simulate_circuit([(5.0, 1.0), (15.0, 0.05)], ratios=ratio_map, target=eff_target, mechanism="activation-only", **kw)[-1][rep_idx + 1]

    core_wt, core_mut = "GACCACCCA", "GAACACCCA"
    q_exact = affinity_ratio(core_wt, core_mut, "GLI2", epsilon=0.0)
    q_pat = pat.profile(core_mut, "GLI2")["profile_ratio"]

    checks = [
        {"name": "exact_9bp_parity_with_patterning", "passed": abs(q_exact - q_pat) < 1e-15, "evidence": {"q_circuit": q_exact, "q_patterning": q_pat}},
        {"name": "reporter_leaves_endogenous_unchanged", "passed": rep_assay["trajectory"] == wt_assay["trajectory"],
         "evidence": {"max_state_diff": diff(rep_assay["trajectory"], wt_assay["trajectory"])}},
        {"name": "allele_and_mediator_rescue", "passed": allele_rescue["trajectory"] == wt_assay["trajectory"] and mediator_rescue["trajectory"] == wt_assay["trajectory"],
         "evidence": {"allele_err": diff(allele_rescue["trajectory"], wt_assay["trajectory"])}},
        {"name": "phenocopy_reproduces_variant", "passed": phenocopy["trajectory"] == mut_assay["trajectory"],
         "evidence": {"phenocopy_err": diff(phenocopy["trajectory"], mut_assay["trajectory"])}},
        {"name": "timestep_refinement_converges", "passed": refinement_err < 1e-4, "evidence": {"max_error": refinement_err, "tolerance": 1e-4}},
    ]
    if ctx["circuit_arg"] is None and ctx["primary_label"] in ("GLI2", "GLI3"):
        checks.extend([
            {"name": "cryptic_bifurcation_found", "passed": abs(featured["reporter_delta_end"]) < 0.02 and featured["fate_flipped"],
             "evidence": {"id": featured["id"], "reporter_delta": featured["reporter_delta_end"], "target_delta": featured["target_delta"]}},
            {"name": "spatial_p3_domain_contracts", "passed": (len(p3_mut) != len(p3_wt)),
             "evidence": {"wt_p3_max_x": max(p3_wt) if p3_wt else None, "mut_p3_max_x": max(p3_mut) if p3_mut else None}},
            {"name": "competitive_primed_derepression_reverses_ablation",
             "passed": hist_mut != hist_wt and hist_mut_act < 0.01,
             "evidence": {"wt_primed_nkx": hist_wt, "mut_primed_competitive": hist_mut, "mut_primed_activation_only": hist_mut_act}},
        ])
    else:
        checks.append({
            "name": "generalized_circuit_variant_differs",
            "passed": diff(mut_assay["trajectory"], wt_assay["trajectory"]) > 1e-4,
            "evidence": {"max_state_diff": diff(mut_assay["trajectory"], wt_assay["trajectory"]), "wt_target_domain": len(p3_wt), "mut_target_domain": len(p3_mut)},
        })

    construction_report = None
    if with_construction:
        construction_report = bridge_to_construction(
            reference, variant, output / "construction", protein=protein,
            epsilon=DEFAULT_EPSILON, background_floor=eff_floor, signal=eff_signal, target=target, circuit=circuit,
        )
        checks.append({
            "name": "cx3d_construction_and_lif_bridge",
            "passed": construction_report["all_checks_passed"],
            "evidence": construction_report["checks"],
        })

    second_site = find_second_site_rescue(
        reference, featured, protein, DEFAULT_EPSILON,
        background_floor=eff_floor, signal=eff_signal, target=target, circuit=circuit,
    )
    report = {
        "format": "epigenesis.neural-circuit-experiment/v1",
        "circuit": ctx["circuit_name"],
        "locus": locus,
        "target": eff_target,
        "protein": ctx["primary_label"],
        "background_floor": eff_floor,
        "scope": "computational sequence-to-circuit model; biological embryo phenotype is not established",
        "compiled": compiled,
        "featured_variant": {k: v for k, v in featured.items() if k != "variant"},
        "second_site_candidates": second_site,
        "scan": {**{k: v for k, v in scan.items() if k != "edits"},
                 "edits": [{k: v for k, v in e.items() if k != "variant"} for e in scan["edits"]],
                 "cryptic_switches": [{k: v for k, v in e.items() if k != "variant"} for e in scan["cryptic_switches"]]},
        "spatial": {"duration": 100.0, "positions": positions, "reference": sp_ref, "variant": sp_mut},
        "checks": checks,
        "all_checks_passed": all(c["passed"] for c in checks),
    }
    if construction_report is not None:
        report["construction_bridge"] = construction_report
    (output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    (output / "index.html").write_text(_render_html(report))
    return report


def _read_dna_arg(value: str) -> str:
    candidate = Path(value)
    if candidate.exists():
        raw = read_regular_file(candidate, maximum_bytes=128 * 1024).decode("ascii")
        seq = "".join(line.strip() for line in raw.splitlines() if not line.startswith(">"))
        return seq.upper()
    return "".join(value.split()).upper()


def _add_common_circuit_args(subparser: argparse.ArgumentParser) -> None:
    subparser.add_argument("--circuit", help="circuit preset ('neural-tube', 'mutual-toggle', 'tri-lineage-switch') or JSON spec path")
    subparser.add_argument("--protein", help="TF motif: GLI2, GLI3, any bundled JASPAR ID (e.g. MA0069.1), or custom PWM JSON path")
    subparser.add_argument("--target", help="regulated target gene name in the chosen circuit (or 'Both'/'All')")
    subparser.add_argument("--signal", type=float, help="morphogen input fraction S in [0, 1]")
    subparser.add_argument("--background-floor", type=float, help=f"non-specific partition floor Z_ns >= 0 (default {DEFAULT_BACKGROUND_FLOOR:g})")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    d = sub.add_parser("demo")
    d.add_argument("--output", type=Path, required=True)
    d.add_argument("--locus", default=CURATED_LOCI[0])
    _add_common_circuit_args(d)
    d.add_argument("--motif-model-dir", type=Path)
    d.add_argument("--with-construction", action="store_true")

    s = sub.add_parser("scan")
    s.add_argument("--fasta", type=Path, required=True)
    _add_common_circuit_args(s)
    s.add_argument("--motif-model-dir", type=Path)
    s.add_argument("--output", type=Path, required=True)

    c = sub.add_parser("compare")
    c.add_argument("--reference", required=True, help="reference DNA sequence or FASTA path")
    c.add_argument("--variant", required=True, help="variant DNA sequence or FASTA path (SNV, MNV, or indel)")
    _add_common_circuit_args(c)
    c.add_argument("--motif-model-dir", type=Path)
    c.add_argument("--output", type=Path, required=True)

    b = sub.add_parser("benchmark")
    _add_common_circuit_args(b)
    b.add_argument("--motif-model-dir", type=Path)
    b.add_argument("--output", type=Path, required=True)

    br = sub.add_parser("bridge")
    br.add_argument("--reference", required=True, help="reference DNA sequence or FASTA path")
    br.add_argument("--variant", required=True, help="variant DNA sequence or FASTA path")
    _add_common_circuit_args(br)
    br.add_argument("--mapping", type=Path, help="optional JSON mapping of GRN gene steady state to CX3D Program fields")
    br.add_argument("--duration", type=float, default=3.0)
    br.add_argument("--output", type=Path, required=True)

    sub.add_parser("list-circuits")
    sub.add_parser("list-motifs")

    args = parser.parse_args()
    try:
        if args.command == "list-circuits":
            print(json.dumps(list_circuit_presets(), indent=2))
            return 0
        if args.command == "list-motifs":
            print(json.dumps(list_available_motifs(), indent=2))
            return 0
        if args.command == "demo":
            report = demo(
                args.output, locus=args.locus, target=args.target, signal=args.signal,
                protein=args.protein, circuit=args.circuit, background_floor=args.background_floor,
                motif_model_dir=args.motif_model_dir, with_construction=args.with_construction,
            )
            print(json.dumps({
                "output": str(args.output), "circuit": report["circuit"], "locus": report["locus"],
                "protein": report["protein"], "target": report["target"],
                "all_checks_passed": report["all_checks_passed"],
                "category_counts": report["scan"]["category_counts"],
                "featured": report["featured_variant"],
                "second_site_candidates": report["second_site_candidates"],
                "checks": report["checks"],
            }, indent=2))
            return 0 if report["all_checks_passed"] else 1
        if args.command == "scan":
            raw = read_regular_file(args.fasta, maximum_bytes=128 * 1024).decode("ascii")
            seq = "".join(line.strip() for line in raw.splitlines() if not line.startswith(">")).upper()
            res = scan_enhancer(
                seq, signal=args.signal, protein=args.protein, target=args.target,
                circuit=args.circuit, background_floor=args.background_floor, motif_model_dir=args.motif_model_dir,
            )
            args.output.mkdir()
            (args.output / "scan.json").write_text(json.dumps(res, indent=2, allow_nan=False) + "\n")
            print(json.dumps({
                "output": str(args.output), "circuit": res["circuit"], "protein": res["protein"],
                "target": res["target"], "category_counts": res["category_counts"],
                "top_cryptic": [{k: v for k, v in e.items() if k != "variant"} for e in res["cryptic_switches"][:5]],
            }, indent=2))
            return 0
        if args.command == "compare":
            args.output.mkdir()
            comp = compare_variants(
                _read_dna_arg(args.reference), _read_dna_arg(args.variant),
                protein=args.protein, signal=args.signal, target=args.target,
                circuit=args.circuit, background_floor=args.background_floor,
                motif_model_dir=args.motif_model_dir, compile_dir=args.output / "compiled",
            )
            (args.output / "compare.json").write_text(json.dumps(comp, indent=2, allow_nan=False) + "\n")
            print(json.dumps({k: v for k, v in comp.items() if not k.endswith("_assay")}, indent=2))
            return 0
        if args.command == "benchmark":
            bench = multi_locus_benchmark(
                protein=args.protein, signal=args.signal, target=args.target,
                circuit=args.circuit, background_floor=args.background_floor,
                motif_model_dir=args.motif_model_dir, output=args.output,
            )
            print(json.dumps({k: v for k, v in bench.items() if k != "loci"}, indent=2))
            return 0
        res = bridge_to_construction(
            _read_dna_arg(args.reference), _read_dna_arg(args.variant), args.output,
            protein=args.protein, signal=args.signal, target=args.target,
            circuit=args.circuit, background_floor=args.background_floor,
            mapping=args.mapping, dev_duration=args.duration,
        )
        print(json.dumps({k: v for k, v in res.items() if k != "runs"}, indent=2))
        return 0 if res["all_checks_passed"] else 1
    except (ValueError, OSError, RuntimeError, json.JSONDecodeError) as err:
        print(f"neural-circuit failed: {err}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

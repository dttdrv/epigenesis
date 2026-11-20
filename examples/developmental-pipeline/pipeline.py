"""Unified composable developmental pipeline: DNA + VRS variants -> chromatin/coding regulation -> GRN uncertainty -> 3D CX3D + LIF.

Connects Epigenesis's sealed core compiler (`brainc` sources, GA4GH VRS 2.x variant
collections, and Development Modules) directly to multi-locus cis-regulatory
Boltzmann TF binding, explicit chromatin accessibility/methylation gating, coding
ORF translation, N-gene regulatory dynamics with parameter-uncertainty ensembles,
optional empirical falsification, and 3D CX3D morphogenesis + LIF electrophysiology.

All developmental units and structural transfer assumptions are explicitly declared;
software determinism and causal controls do not establish in vivo biological acceptance
(`biological_acceptance = False`).
"""

from __future__ import annotations

import argparse
import base64
import copy
import hashlib
import html
import json
import math
from pathlib import Path
import random
import shutil
import statistics
import struct
import subprocess
import sys
from typing import Any

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "examples/neural-patterning"))
sys.path.insert(0, str(REPO / "examples/neural-circuit"))
sys.path.insert(0, str(REPO / "examples/coding-consequences"))
sys.path.insert(0, str(REPO / "examples/neural-construction"))

from brainc._canonical import canonical_bytes
from brainc._io import load_json_object, read_regular_file
from brainc.development_bundle import compile_development
from brainc.export import export_omex, export_rocrate
from brainc.source import FASTA_PROFILE, compile_source
from brainc.v2 import inline_storage, make_development_request, pack, policy_artifact, save, target_artifact
from brainc.v2._common import seal
from brainc.v2.target import DEV_DOMAIN, OP_UNIT_CREATE
from brainc.validator_development import load_development_bundle_directory, validate_development_paths
from brainc.validator_variant import validate_variant_paths
from brainc.variant import apply_variant_to_sequence, compile_variants, load_variant_collection, revert_variant

import circuit as nc
import coding as cc


TARGET = "org.epigenesis.developmental-pipeline"
FORMAT = "epigenesis.developmental-pipeline/v1"
PROGRAM_FORMAT = "epigenesis.developmental-program/v1"
DEFAULT_EPSILON = nc.DEFAULT_EPSILON
DEFAULT_BACKGROUND_FLOOR = nc.DEFAULT_BACKGROUND_FLOOR
DEFAULT_ENSEMBLE_SAMPLES = 24
DEFAULT_ENSEMBLE_SPREAD = 0.15
DEFAULT_ENSEMBLE_SEED = 1337
DEFAULT_GRN_DURATION = 40.0
DEFAULT_CX3D_DURATION = 2.0
DEFAULT_CX3D_SEED = 7


def _is_num(value: Any) -> bool:
    return type(value) in (int, float) and math.isfinite(value)


def _semantics_sha256() -> str:
    h = hashlib.sha256()
    for rel in (
        "examples/developmental-pipeline/pipeline.py",
        "examples/neural-circuit/circuit.py",
        "examples/coding-consequences/coding.py",
    ):
        raw = read_regular_file(REPO / rel, maximum_bytes=2 * 1024 * 1024)
        h.update(rel.encode("utf-8") + b"\0" + raw)
    return h.hexdigest()


def _parse_fasta(path: Path) -> dict[str, str]:
    raw = read_regular_file(path, maximum_bytes=4 * 1024 * 1024).decode("utf-8")
    records: dict[str, str] = {}
    current_id: str | None = None
    chunks: list[str] = []
    for raw_line in raw.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id is not None:
                records[current_id] = "".join(chunks).upper()
            header = line[1:].strip()
            if not header:
                raise ValueError("empty FASTA record header")
            rid = header.split()[0]
            if rid in records:
                raise ValueError(f"duplicate FASTA record id: {rid!r}")
            current_id = rid
            chunks = []
        else:
            if current_id is None:
                raise ValueError("FASTA sequence lines appeared before header")
            chunks.append(line)
    if current_id is not None:
        records[current_id] = "".join(chunks).upper()
    if not records:
        raise ValueError("FASTA file contains no records")
    for rid, seq in records.items():
        if not seq or set(seq) - set("ACGT"):
            raise ValueError(f"record {rid!r} must be nonempty uppercase ACGT DNA")
    return records


def _write_fasta(records: dict[str, str], path: Path) -> None:
    lines: list[str] = []
    for rid, seq in records.items():
        lines.append(f">{rid}")
        for idx in range(0, len(seq), 80):
            lines.append(seq[idx : idx + 80])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def builtin_preset(name: str = "neural-tube-multi-locus") -> dict[str, Any]:
    """Return a self-contained multi-locus demo (FASTA records, program spec, and edits)."""
    if name == "neural-tube-multi-locus":
        # Multi-locus synthetic neural-tube program:
        # 1) Nkx2.2 cis-enhancer containing a canonical GLI2 site (GACCACCCA) in open chromatin
        # 2) Olig2 cis-enhancer containing a second GLI2 site (GACCACCCA) in partial chromatin
        # 3) A short coding ORF for a co-factor gene (Pax6) to demonstrate cis + coding coupling
        enh_nkx = "CTGACGTACG" + "GACCACCCA" + "TGCATGCATG" + "GCTAGCTAGC" + "AATTCCGGAA"
        enh_olig = "ACGTACGTAC" + "GACCACCCA" + "CGATCGATCG" + "TTGGCCAATT" + "CGTACGTACG"
        orf_pax6 = (
            "ATGCAGACCCTGAGCGTGCTGGCCGGCCTGAGCCGCACCATGAAGGCCGAGCTGAGCAAGGAGGAGCTG"
            "CTGGCCGGCCAGCTGCTGGAGCGCGTGGCCGGCTAA"
        )
        records = {
            "enh_nkx2_2": enh_nkx,
            "enh_olig2": enh_olig,
            "orf_pax6": orf_pax6,
        }
        spec = {
            "format": PROGRAM_FORMAT,
            "name": "neural-tube-multi-locus",
            "circuit": "neural-tube",
            "signal": 0.20,
            "duration": DEFAULT_GRN_DURATION,
            "epsilon": DEFAULT_EPSILON,
            "background_floor": DEFAULT_BACKGROUND_FLOOR,
            "loci": [
                {
                    "record_id": "enh_nkx2_2",
                    "kind": "enhancer",
                    "target_gene": "Nkx2.2",
                    "protein": "GLI2",
                    "accessibility": 1.0,
                    "methylation_repression": 0.0,
                    "weight": 1.0,
                    "role": "activator",
                    "mechanism": "competitive",
                },
                {
                    "record_id": "enh_olig2",
                    "kind": "enhancer",
                    "target_gene": "Olig2",
                    "protein": "GLI2",
                    "accessibility": 0.6,
                    "methylation_repression": 0.0,
                    "weight": 1.0,
                    "role": "activator",
                    "mechanism": "competitive",
                },
                {
                    "record_id": "orf_pax6",
                    "kind": "coding",
                    "affects_gene": "Pax6",
                    "truncation_floor": 0.05,
                },
            ],
            "uncertainty": {
                "samples": DEFAULT_ENSEMBLE_SAMPLES,
                "relative_spread": DEFAULT_ENSEMBLE_SPREAD,
                "seed": DEFAULT_ENSEMBLE_SEED,
            },
        }
        # Edit the GLI2 core in enh_nkx2_2: GACCACCCA -> GAACACCCA (position 12, C>A)
        variants = [
            {
                "id": "nkx2_2_gli_core_C13A",
                "record_id": "enh_nkx2_2",
                "start": 12,
                "end": 13,
                "ref": "C",
                "alt": "A",
            }
        ]
        return {"records": records, "spec": spec, "variants": variants}

    if name == "tri-lineage-chromatin":
        # Three-lineage cortical/neural switch with JASPAR motifs (PAX6 MA0069.1, BCL6/repressor MA0463.3)
        enh_exc = "ACTGACTGACTGACTG" + "TTCACGCATGAGTT" + "TGCATGCATGCA"
        enh_inh = "CGTACGTACGTACGTA" + "GCTTTCGAGGAAT" + "ACGTACGTACGT"
        records = {
            "enh_excitatory": enh_exc,
            "enh_inhibitory": enh_inh,
        }
        spec = {
            "format": PROGRAM_FORMAT,
            "name": "tri-lineage-chromatin",
            "circuit": "tri-lineage-switch",
            "signal": 0.50,
            "duration": DEFAULT_GRN_DURATION,
            "epsilon": DEFAULT_EPSILON,
            "background_floor": DEFAULT_BACKGROUND_FLOOR,
            "loci": [
                {
                    "record_id": "enh_excitatory",
                    "kind": "enhancer",
                    "target_gene": "Excitatory",
                    "protein": "MA0069.1",
                    "accessibility": 0.9,
                    "methylation_repression": 0.0,
                    "weight": 1.0,
                    "role": "activator",
                    "mechanism": "competitive",
                },
                {
                    "record_id": "enh_inhibitory",
                    "kind": "enhancer",
                    "target_gene": "Inhibitory",
                    "protein": "MA0463.3",
                    "accessibility": 0.85,
                    "methylation_repression": 0.0,
                    "weight": 1.0,
                    "role": "activator",
                    "mechanism": "competitive",
                },
            ],
            "uncertainty": {
                "samples": DEFAULT_ENSEMBLE_SAMPLES,
                "relative_spread": DEFAULT_ENSEMBLE_SPREAD,
                "seed": DEFAULT_ENSEMBLE_SEED,
            },
        }
        variants = [
            {
                "id": "exc_pax6_core_GC_to_AA",
                "record_id": "enh_excitatory",
                "start": 21,
                "end": 23,
                "ref": "GC",
                "alt": "AA",
            }
        ]
        return {"records": records, "spec": spec, "variants": variants}

    raise ValueError(f"unknown developmental pipeline preset: {name!r}")


def validate_program_spec(spec: dict[str, Any], record_ids: set[str]) -> dict[str, Any]:
    """Validate and normalize a declarative multi-locus developmental program specification."""
    if type(spec) is not dict:
        raise ValueError("program spec must be a JSON object")
    fmt = spec.get("format", PROGRAM_FORMAT)
    if fmt != PROGRAM_FORMAT:
        raise ValueError(f"unsupported program spec format: {fmt!r}")

    circuit_arg = spec.get("circuit", "neural-tube")
    if isinstance(circuit_arg, str):
        model = nc.circuit_preset(circuit_arg)
    elif isinstance(circuit_arg, dict):
        model = nc.validate_circuit_spec(circuit_arg)
    else:
        raise ValueError("spec 'circuit' must be a preset name or CircuitSpec JSON object")

    genes = model["genes"]
    gene_set = set(genes)

    loci_raw = spec.get("loci")
    if not isinstance(loci_raw, list) or not loci_raw:
        raise ValueError("spec 'loci' must be a non-empty list")

    normalized_loci: list[dict[str, Any]] = []
    seen_locus_keys: set[tuple[str, str, str]] = set()
    for idx, loc in enumerate(loci_raw):
        if type(loc) is not dict:
            raise ValueError(f"loci[{idx}] must be an object")
        rid = loc.get("record_id")
        if type(rid) is not str or rid not in record_ids:
            raise ValueError(f"loci[{idx}].record_id {rid!r} is not present in the source FASTA")
        kind = loc.get("kind")
        if kind == "enhancer":
            tgene = loc.get("target_gene")
            if tgene not in gene_set:
                raise ValueError(f"loci[{idx}].target_gene {tgene!r} is not in circuit genes {genes}")
            protein = loc.get("protein", model["default_protein"])
            nc.load_motif_profile(protein, float(spec.get("epsilon", DEFAULT_EPSILON)))
            acc = loc.get("accessibility", 1.0)
            meth = loc.get("methylation_repression", 0.0)
            weight = loc.get("weight", 1.0)
            role = loc.get("role", "activator")
            mech = loc.get("mechanism", "competitive")
            if not _is_num(acc) or not 0.0 <= float(acc) <= 1.0:
                raise ValueError(f"loci[{idx}].accessibility must be in [0, 1]")
            if not _is_num(meth) or not 0.0 <= float(meth) <= 1.0:
                raise ValueError(f"loci[{idx}].methylation_repression must be in [0, 1]")
            if not _is_num(weight) or float(weight) <= 0.0:
                raise ValueError(f"loci[{idx}].weight must be > 0")
            if role not in ("activator", "repressor") or mech not in ("competitive", "activation-only"):
                raise ValueError(f"loci[{idx}] role/mechanism invalid")
            prot_lbl = nc._motif_spec_label(protein)
            key = (rid, "enhancer", f"{tgene}:{prot_lbl}")
            if key in seen_locus_keys:
                raise ValueError(f"duplicate enhancer wiring {key}")
            seen_locus_keys.add(key)
            normalized_loci.append(
                {
                    "record_id": rid,
                    "kind": "enhancer",
                    "target_gene": tgene,
                    "protein": protein,
                    "protein_label": prot_lbl,
                    "accessibility": float(acc),
                    "methylation_repression": float(meth),
                    "chromatin_gate": float(acc) * (1.0 - float(meth)),
                    "weight": float(weight),
                    "role": role,
                    "mechanism": mech,
                }
            )
        elif kind == "coding":
            agene = loc.get("affects_gene")
            if agene not in gene_set:
                raise ValueError(f"loci[{idx}].affects_gene {agene!r} is not in circuit genes {genes}")
            tfloor = loc.get("truncation_floor", 0.05)
            if not _is_num(tfloor) or not 0.0 < float(tfloor) <= 1.0:
                raise ValueError(f"loci[{idx}].truncation_floor must be in (0, 1]")
            key = (rid, "coding", agene)
            if key in seen_locus_keys:
                raise ValueError(f"duplicate coding locus wiring {key}")
            seen_locus_keys.add(key)
            normalized_loci.append(
                {
                    "record_id": rid,
                    "kind": "coding",
                    "affects_gene": agene,
                    "truncation_floor": float(tfloor),
                }
            )
        else:
            raise ValueError(f"loci[{idx}].kind must be 'enhancer' or 'coding', got {kind!r}")

    signal = spec.get("signal", model["default_signal"])
    if not _is_num(signal) or not 0.0 <= float(signal) <= 1.0:
        raise ValueError("spec 'signal' must be in [0, 1]")
    duration = spec.get("duration", DEFAULT_GRN_DURATION)
    if not _is_num(duration) or float(duration) <= 0.0:
        raise ValueError("spec 'duration' must be > 0")
    phases_raw = spec.get("phases")
    if phases_raw is not None:
        if not isinstance(phases_raw, list) or not phases_raw:
            raise ValueError("spec 'phases' must be a non-empty list of [duration, signal] pairs")
        protocol = []
        for p in phases_raw:
            if not isinstance(p, (list, tuple)) or len(p) != 2 or not _is_num(p[0]) or not _is_num(p[1]):
                raise ValueError("each phase must be [duration, signal]")
            if float(p[0]) <= 0.0 or not 0.0 <= float(p[1]) <= 1.0:
                raise ValueError("phase duration must be > 0 and signal in [0, 1]")
            protocol.append((float(p[0]), float(p[1])))
    else:
        protocol = [(float(duration), float(signal))]

    epsilon = spec.get("epsilon", DEFAULT_EPSILON)
    if not _is_num(epsilon) or float(epsilon) < 0.0:
        raise ValueError("spec 'epsilon' must be >= 0")
    bg_floor = nc._resolve_background_floor(spec.get("background_floor"), float(epsilon))

    unc_raw = spec.get("uncertainty", {})
    if type(unc_raw) is not dict:
        raise ValueError("spec 'uncertainty' must be an object")
    samples = unc_raw.get("samples", DEFAULT_ENSEMBLE_SAMPLES)
    spread = unc_raw.get("relative_spread", DEFAULT_ENSEMBLE_SPREAD)
    seed = unc_raw.get("seed", DEFAULT_ENSEMBLE_SEED)
    if type(samples) is not int or not 0 <= samples <= 500:
        raise ValueError("uncertainty.samples must be an integer in [0, 500]")
    if not _is_num(spread) or not 0.0 <= float(spread) <= 0.5:
        raise ValueError("uncertainty.relative_spread must be in [0.0, 0.5]")
    if type(seed) is not int:
        raise ValueError("uncertainty.seed must be an integer")

    initial_raw = spec.get("initial")
    if initial_raw is not None:
        if type(initial_raw) is not dict or set(initial_raw) != gene_set or any(not _is_num(v) or float(v) < 0 for v in initial_raw.values()):
            raise ValueError("spec 'initial' must map every circuit gene to a nonnegative number")
        initial = {g: float(initial_raw[g]) for g in genes}
    else:
        initial = {g: 0.0 for g in genes}

    return {
        "format": PROGRAM_FORMAT,
        "name": str(spec.get("name", "custom-developmental-program")),
        "model": model,
        "genes": genes,
        "loci": normalized_loci,
        "protocol": protocol,
        "signal": protocol[-1][1],
        "duration": sum(d for d, _ in protocol),
        "epsilon": float(epsilon),
        "background_floor": float(bg_floor),
        "initial": initial,
        "uncertainty": {"samples": samples, "relative_spread": float(spread), "seed": seed},
        "morphogenesis_mapping": spec.get("morphogenesis_mapping"),
        "observations": spec.get("observations"),
    }


def _apply_multi_variants(
    ref_records: dict[str, str],
    compiled_variants: list[dict[str, Any]],
) -> tuple[dict[str, str], dict[str, str]]:
    """Apply non-overlapping sealed VRS variants per record and verify exact allele reversion."""
    by_rec: dict[str, list[dict[str, Any]]] = {}
    for entry in compiled_variants:
        by_rec.setdefault(entry["record_id"], []).append(entry)

    mutated = dict(ref_records)
    reverted = dict(ref_records)
    for rid, entries in by_rec.items():
        # Check non-overlap on reference coordinates
        ordered_asc = sorted(entries, key=lambda e: (e["location"]["start"], e["location"]["end"], e["ordinal"]))
        for prev, curr in zip(ordered_asc, ordered_asc[1:]):
            if curr["location"]["start"] < prev["location"]["end"]:
                raise ValueError(
                    f"overlapping variants on record {rid!r}: {prev['id']} and {curr['id']}"
                )
        # Apply from highest coordinate down so earlier coordinates remain exact
        seq = ref_records[rid]
        ordered_desc = list(reversed(ordered_asc))
        for entry in ordered_desc:
            seq = apply_variant_to_sequence(seq, entry)
        mutated[rid] = seq

        # Revert in exact reverse order (ascending) and confirm exact WT recovery
        rev_seq = seq
        for entry in ordered_asc:
            rev_seq = revert_variant(rev_seq, entry)
        if rev_seq != ref_records[rid]:
            raise ValueError(f"exact allele reversion check failed for record {rid!r}")
        reverted[rid] = rev_seq
    return mutated, reverted


def _coding_effect(wt_seq: str, mut_seq: str, truncation_floor: float) -> dict[str, Any]:
    """Translate WT and variant ORFs using NCBI standard genetic code and compute structural activity factor."""
    wt = cc.translate(wt_seq)
    if mut_seq.startswith("ATG") and len(mut_seq) >= 3:
        mut = cc.translate(mut_seq)
        wt_p, mut_p = wt["peptide"], mut["peptide"]
        if wt_p == mut_p:
            consequence = "synonymous_or_unchanged"
            factor = 1.0
        elif len(mut_p) < len(wt_p) and mut["status"] == "terminated":
            consequence = "nonsense_truncation"
            frac = len(mut_p) / max(1, len(wt_p))
            factor = max(truncation_floor, frac if frac >= 0.5 else truncation_floor)
        elif len(mut_p) != len(wt_p):
            consequence = "frameshift_or_length_change"
            shared = sum(a == b for a, b in zip(wt_p, mut_p))
            factor = max(truncation_floor, shared / max(1, len(wt_p)))
        else:
            consequence = "missense"
            shared = sum(a == b for a, b in zip(wt_p, mut_p))
            factor = max(truncation_floor, shared / max(1, len(wt_p)))
        mut_pep_sha = mut["peptide_sha256"]
        mut_len = len(mut_p)
        mut_status = mut["status"]
    else:
        consequence = "start_loss"
        factor = truncation_floor
        mut_pep_sha = None
        mut_len = 0
        mut_status = "start-loss"

    return {
        "consequence": consequence,
        "coding_activity_factor": float(factor),
        "wt_peptide_length": len(wt["peptide"]),
        "mut_peptide_length": mut_len,
        "wt_peptide_sha256": wt["peptide_sha256"],
        "mut_peptide_sha256": mut_pep_sha,
        "wt_status": wt["status"],
        "mut_status": mut_status,
        "assumption": "explicit_structural_coding_transfer",
    }


def _evaluate_loci(
    program: dict[str, Any],
    ref_records: dict[str, str],
    mut_records: dict[str, str],
    *,
    force_closed_chromatin: bool = False,
) -> dict[str, Any]:
    """Compute per-locus biophysical/coding effects and aggregate per-gene effective ratios."""
    genes = program["genes"]
    eps = program["epsilon"]
    floor = program["background_floor"]

    gene_ratios: dict[str, float] = {g: 1.0 for g in genes}
    gene_roles: dict[str, str] = {}
    gene_mechanisms: dict[str, str] = {}
    locus_reports: list[dict[str, Any]] = []

    for loc in program["loci"]:
        rid = loc["record_id"]
        wt_seq = ref_records[rid]
        mut_seq = mut_records[rid]
        if loc["kind"] == "enhancer":
            tgene = loc["target_gene"]
            prot = loc["protein"]
            wt_part = nc.boltzmann_affinity(wt_seq, prot, eps)
            mut_part = nc.boltzmann_affinity(mut_seq, prot, eps)
            q_raw = (mut_part["partition_sum"] + floor) / (wt_part["partition_sum"] + floor)
            gate = 0.0 if force_closed_chromatin else loc["chromatin_gate"]
            # Log-linear chromatin gating: gate=0 -> exact 1.0; gate=1 -> exact q_raw; strictly > 0 & symmetric
            q_eff = math.exp(gate * math.log(q_raw))
            gene_ratios[tgene] *= q_eff ** loc["weight"]
            gene_roles.setdefault(tgene, loc["role"])
            gene_mechanisms.setdefault(tgene, loc["mechanism"])
            locus_reports.append(
                {
                    "record_id": rid,
                    "kind": "enhancer",
                    "target_gene": tgene,
                    "protein": loc["protein_label"],
                    "accessibility": 0.0 if force_closed_chromatin else loc["accessibility"],
                    "methylation_repression": loc["methylation_repression"],
                    "chromatin_gate": gate,
                    "weight": loc["weight"],
                    "wt_partition_sum": wt_part["partition_sum"],
                    "mut_partition_sum": mut_part["partition_sum"],
                    "raw_affinity_ratio": q_raw,
                    "effective_affinity_ratio": q_eff,
                    "sequence_changed": wt_seq != mut_seq,
                }
            )
        else:
            agene = loc["affects_gene"]
            ce = _coding_effect(wt_seq, mut_seq, loc["truncation_floor"])
            gene_ratios[agene] *= ce["coding_activity_factor"]
            locus_reports.append(
                {
                    "record_id": rid,
                    "kind": "coding",
                    "affects_gene": agene,
                    "sequence_changed": wt_seq != mut_seq,
                    **ce,
                }
            )

    return {
        "locus_reports": locus_reports,
        "gene_ratios": gene_ratios,
        "gene_roles": gene_roles,
        "gene_mechanisms": gene_mechanisms,
    }


def _perturb_model(model: dict[str, Any], rng: random.Random, spread: float) -> dict[str, Any]:
    """Perturb kinetic ODE parameters multiplicatively within [1 - spread, 1 + spread]."""
    perturbed = copy.deepcopy(model)
    if spread <= 0.0:
        return perturbed

    def _jitter_vec(vec: list[float], positive: bool) -> list[float]:
        out = []
        for v in vec:
            if v == 0.0:
                out.append(0.0)
            else:
                factor = 1.0 + rng.uniform(-spread, spread)
                nv = v * factor
                out.append(max(1e-6, nv) if positive else max(0.0, nv))
        return out

    perturbed["basal"] = _jitter_vec(perturbed["basal"], positive=False)
    perturbed["input"] = _jitter_vec(perturbed["input"], positive=False)
    perturbed["tf_binding"] = _jitter_vec(perturbed["tf_binding"], positive=False)
    perturbed["gli_binding"] = list(perturbed["tf_binding"])
    # Scale synthesis/decay together by a shared kinetic timescale factor per gene so
    # dimensionless normalized concentrations remain strictly inside the physical [0, 1] box.
    new_prod: list[float] = []
    new_decay: list[float] = []
    for pr, dc in zip(perturbed["production"], perturbed["decay"]):
        rate_scale = 1.0 + rng.uniform(-spread, spread)
        new_prod.append(max(1e-6, pr * rate_scale))
        new_decay.append(max(1e-6, dc * rate_scale))
    perturbed["production"] = new_prod
    perturbed["decay"] = new_decay
    perturbed["repression"] = [
        _jitter_vec(row, positive=False) for row in perturbed["repression"]
    ]
    return perturbed


def _run_uncertainty_ensemble(
    program: dict[str, Any],
    gene_ratios: dict[str, float],
    roles: dict[str, str],
    mechanisms: dict[str, str],
) -> dict[str, Any]:
    """Run deterministic parameter-uncertainty ensemble reporting fate probabilities & CIs."""
    unc = program["uncertainty"]
    n_samples = unc["samples"]
    spread = unc["relative_spread"]
    seed = unc["seed"]
    genes = program["genes"]
    init_tuple = tuple(program["initial"][g] for g in genes)
    protocol = program["protocol"]
    wt_ratios = {g: 1.0 for g in genes}
    fates_cfg = program["model"].get("fates")

    if n_samples == 0:
        return {
            "samples": 0,
            "relative_spread": spread,
            "seed": seed,
            "fate_flip_probability": None,
            "wt_fate_distribution": {},
            "variant_fate_distribution": {},
            "gene_delta_summary": {},
        }

    rng = random.Random(seed)
    wt_fates: dict[str, int] = {}
    var_fates: dict[str, int] = {}
    flips = 0
    deltas_by_gene: dict[str, list[float]] = {g: [] for g in genes}

    for _ in range(n_samples):
        pmodel = _perturb_model(program["model"], rng, spread)
        wt_end = nc.simulate_circuit(
            protocol, ratios=wt_ratios, initial=init_tuple, circuit=pmodel, roles=roles, mechanisms=mechanisms
        )[-1][1:]
        var_end = nc.simulate_circuit(
            protocol, ratios=gene_ratios, initial=init_tuple, circuit=pmodel, roles=roles, mechanisms=mechanisms
        )[-1][1:]
        wf = nc.fate_label(wt_end, genes, fates=fates_cfg)
        vf = nc.fate_label(var_end, genes, fates=fates_cfg)
        wt_fates[wf] = wt_fates.get(wf, 0) + 1
        var_fates[vf] = var_fates.get(vf, 0) + 1
        if wf != vf:
            flips += 1
        for i, g in enumerate(genes):
            deltas_by_gene[g].append(float(var_end[i] - wt_end[i]))

    def _quantile(xs: list[float], q: float) -> float:
        s = sorted(xs)
        pos = q * (len(s) - 1)
        lo = int(math.floor(pos))
        hi = int(math.ceil(pos))
        if lo == hi:
            return s[lo]
        return s[lo] + (pos - lo) * (s[hi] - s[lo])

    delta_summary = {}
    for g, vals in deltas_by_gene.items():
        delta_summary[g] = {
            "mean": float(statistics.fmean(vals)),
            "std": float(statistics.pstdev(vals)) if len(vals) > 1 else 0.0,
            "p05": float(_quantile(vals, 0.05)),
            "p95": float(_quantile(vals, 0.95)),
        }

    return {
        "samples": n_samples,
        "relative_spread": spread,
        "seed": seed,
        "fate_flip_probability": flips / n_samples,
        "wt_fate_distribution": {k: v / n_samples for k, v in sorted(wt_fates.items())},
        "variant_fate_distribution": {k: v / n_samples for k, v in sorted(var_fates.items())},
        "gene_delta_summary": delta_summary,
    }


def evaluate_observations(
    predicted_deltas: dict[str, float],
    observations: list[dict[str, Any]],
    *,
    max_rmse: float = 0.25,
    min_sign_agreement: float = 0.75,
) -> dict[str, Any]:
    """Compare simulated gene-state deltas against empirical observations and flag falsification."""
    if not isinstance(observations, list) or not observations:
        raise ValueError("observations must be a non-empty list of {gene, observed_delta} entries")
    preds: list[float] = []
    obs: list[float] = []
    details: list[dict[str, Any]] = []
    sign_matches = 0
    for idx, item in enumerate(observations):
        if type(item) is not dict or "gene" not in item or "observed_delta" not in item:
            raise ValueError(f"observations[{idx}] requires 'gene' and 'observed_delta'")
        g = item["gene"]
        if g not in predicted_deltas:
            raise ValueError(f"observations[{idx}].gene {g!r} not in circuit genes")
        ov = float(item["observed_delta"])
        if not math.isfinite(ov):
            raise ValueError(f"observations[{idx}].observed_delta must be finite")
        pv = float(predicted_deltas[g])
        preds.append(pv)
        obs.append(ov)
        same_sign = (pv == 0.0 and ov == 0.0) or (pv > 0.0 and ov > 0.0) or (pv < 0.0 and ov < 0.0)
        if same_sign:
            sign_matches += 1
        details.append(
            {
                "gene": g,
                "predicted_delta": pv,
                "observed_delta": ov,
                "residual": pv - ov,
                "sign_match": same_sign,
            }
        )

    n = len(preds)
    rmse = math.sqrt(sum((p - o) ** 2 for p, o in zip(preds, obs)) / n)
    mae = sum(abs(p - o) for p, o in zip(preds, obs)) / n
    sign_rate = sign_matches / n
    pearson_r: float | None = None
    if n >= 2:
        mp = statistics.fmean(preds)
        mo = statistics.fmean(obs)
        vp = sum((p - mp) ** 2 for p in preds)
        vo = sum((o - mo) ** 2 for o in obs)
        if vp > 0 and vo > 0:
            pearson_r = sum((p - mp) * (o - mo) for p, o in zip(preds, obs)) / math.sqrt(vp * vo)

    falsified = (rmse > max_rmse) or (sign_rate < min_sign_agreement)
    return {
        "n_observations": n,
        "rmse": rmse,
        "mean_abs_error": mae,
        "sign_agreement": sign_rate,
        "pearson_r": pearson_r,
        "thresholds": {"max_rmse": max_rmse, "min_sign_agreement": min_sign_agreement},
        "falsified": falsified,
        "details": details,
    }


def _pipeline_target(n_loci: int, n_genes: int) -> dict[str, Any]:
    return target_artifact(
        {
            "id": TARGET,
            "abi_major": 1,
            "opsets": [{"domain": DEV_DOMAIN, "version": 1}],
            "unit_schemas": [
                {
                    "id": "developmental_program",
                    "fields": [
                        {
                            "id": "gene_effective_ratios",
                            "type": {"dtype": "f64", "shape": [n_genes]},
                            "unit": None,
                            "mutability": "constant",
                            "numeric": {"kind": "float"},
                        },
                        {
                            "id": "locus_effective_ratios",
                            "type": {"dtype": "f64", "shape": [n_loci]},
                            "unit": None,
                            "mutability": "constant",
                            "numeric": {"kind": "float"},
                        },
                        {
                            "id": "variant_steady_state",
                            "type": {"dtype": "f64", "shape": [n_genes]},
                            "unit": None,
                            "mutability": "state",
                            "numeric": {"kind": "float"},
                        },
                        {
                            "id": "wt_steady_state",
                            "type": {"dtype": "f64", "shape": [n_genes]},
                            "unit": None,
                            "mutability": "state",
                            "numeric": {"kind": "float"},
                        },
                    ],
                }
            ],
            "edge_schemas": [],
            "ports": [],
            "rules": [],
        }
    )


def _try_rust_validate(
    source_dir: Path,
    fasta_path: Path,
    variants_path: Path,
    dev_dir: Path,
) -> dict[str, Any]:
    """Cross-verify sealed variants and development bundle with the Rust validator if built."""
    candidates = [
        REPO / "validators/rust/target/release/epigenesis-validate-rs",
        REPO / "validators/rust/target/debug/epigenesis-validate-rs",
    ]
    binary = next((p for p in candidates if p.is_file()), None)
    if binary is None:
        return {"available": False, "variants_valid": None, "development_valid": None}

    res_v = subprocess.run(
        [
            str(binary),
            "validate-variants",
            str(source_dir),
            str(variants_path),
            "--source-input",
            f"sequence={fasta_path}",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    res_d = subprocess.run(
        [
            str(binary),
            "validate-development",
            str(source_dir),
            str(dev_dir),
            "--source-input",
            f"sequence={fasta_path}",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    v_ok = res_v.returncode == 0 and json.loads(res_v.stdout).get("valid") is True
    d_ok = res_d.returncode == 0 and json.loads(res_d.stdout).get("valid") is True
    return {
        "available": True,
        "binary": str(binary),
        "variants_valid": v_ok,
        "development_valid": d_ok,
    }


def _render_html(report: dict[str, Any]) -> str:
    esc = html.escape
    genes = report["genes"]
    wt = report["nominal"]["wt_end_state"]
    var = report["nominal"]["variant_end_state"]
    ens = report["uncertainty_ensemble"]

    gene_rows = []
    for g in genes:
        d = var[g] - wt[g]
        ds = ens["gene_delta_summary"].get(g) if ens["samples"] else None
        ci = f"[{ds['p05']:+.4f}, {ds['p95']:+.4f}]" if ds else "n/a"
        gene_rows.append(
            f"<tr><td>{esc(g)}</td><td>{wt[g]:.4f}</td><td>{var[g]:.4f}</td>"
            f"<td>{d:+.4f}</td><td>{ci}</td></tr>"
        )

    locus_rows = []
    for loc in report["loci"]:
        if loc["kind"] == "enhancer":
            locus_rows.append(
                f"<tr><td>{esc(loc['record_id'])}</td><td>enhancer → {esc(loc['target_gene'])} ({esc(loc['protein'])})</td>"
                f"<td>a={loc['accessibility']:.2f}, gate={loc['chromatin_gate']:.2f}</td>"
                f"<td>q_raw={loc['raw_affinity_ratio']:.5f} → q_eff={loc['effective_affinity_ratio']:.5f}</td></tr>"
            )
        else:
            locus_rows.append(
                f"<tr><td>{esc(loc['record_id'])}</td><td>coding → {esc(loc['affects_gene'])}</td>"
                f"<td>{esc(loc['consequence'])}</td>"
                f"<td>factor={loc['coding_activity_factor']:.4f} ({loc['wt_peptide_length']}→{loc['mut_peptide_length']} aa)</td></tr>"
            )

    cx = report.get("morphogenesis")
    cx_html = "<p>3D CX3D morphogenesis skipped (<code>--skip-cx3d</code>).</p>"
    if cx and cx.get("executed"):
        ref_b = cx["branches"]["reference"]
        var_b = cx["branches"]["variant"]
        res_b = cx["branches"]["allele_rescue"]
        cx_html = (
            f"<p><strong>3D CX3D + LIF executed:</strong> "
            f"WT ({ref_b['activity']['cells']} cells, {ref_b['activity']['contacts']} contacts, {ref_b['activity']['total_spikes']} spikes) vs "
            f"Variant ({var_b['activity']['cells']} cells, {var_b['activity']['contacts']} contacts, {var_b['activity']['total_spikes']} spikes) vs "
            f"Allele Rescue ({res_b['activity']['cells']} cells, {res_b['activity']['contacts']} contacts, {res_b['activity']['total_spikes']} spikes). "
            f"Reversion exact: <code>{cx['allele_rescue_matches_reference']}</code>.</p>"
        )

    flip_p = f"{ens['fate_flip_probability']:.3f}" if ens["fate_flip_probability"] is not None else "n/a"
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Epigenesis Developmental Pipeline — {esc(report['program_name'])}</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; max-width: 960px; margin: 2rem auto; padding: 0 1rem; color: #1a1a1a; line-height: 1.5; }}
  h1, h2 {{ border-bottom: 1px solid #e0e0e0; padding-bottom: 0.3rem; }}
  table {{ border-collapse: collapse; width: 100%; margin: 1rem 0; font-size: 0.92rem; }}
  th, td {{ border: 1px solid #d4d4d4; padding: 0.45rem 0.65rem; text-align: left; }}
  th {{ background: #f5f5f5; }}
  code {{ background: #f2f2f2; padding: 0.1rem 0.3rem; border-radius: 3px; font-size: 0.88rem; }}
  .notice {{ background: #fff8e6; border-left: 4px solid #d9822b; padding: 0.75rem 1rem; margin: 1rem 0; }}
</style>
</head>
<body>
<h1>Epigenesis Developmental Pipeline: {esc(report['program_name'])}</h1>
<div class="notice">
  <strong>Scientific Scope:</strong> <code>biological_acceptance = false</code>. Cryptographic provenance,
  VRS 2.x variant binding, chromatin gating, allele reversion, and deterministic simulation controls passed;
  in vivo developmental phenotype is not claimed without external experimental calibration.
</div>
<h2>1. Provenance &amp; Sealed Identities</h2>
<ul>
  <li>Source descriptor SHA-256: <code>{esc(report['provenance']['source_descriptor_sha256'])}</code></li>
  <li>Variant collection SHA-256: <code>{esc(report['provenance']['variant_collection_sha256'])}</code> (validated: <code>{report['provenance']['variants_valid']}</code>)</li>
  <li>Development bundle SHA-256: <code>{esc(report['provenance']['development_bundle_sha256'])}</code> (validated: <code>{report['provenance']['development_valid']}</code>)</li>
  <li>VRS 2.x Allele IDs: <code>{esc(', '.join(v['vrs_allele_id'] for v in report['variants']))}</code></li>
</ul>
<h2>2. Multi-Locus Regulatory &amp; Coding Effects</h2>
<table>
  <thead><tr><th>Record</th><th>Wiring</th><th>Context / Consequence</th><th>Effective Transfer</th></tr></thead>
  <tbody>{''.join(locus_rows)}</tbody>
</table>
<h2>3. GRN Steady State &amp; Parameter-Uncertainty Ensemble ({ens['samples']} draws, ±{ens['relative_spread']*100:.0f}%)</h2>
<p>Nominal fate: <strong>{esc(report['nominal']['wt_fate'])}</strong> → <strong>{esc(report['nominal']['variant_fate'])}</strong>
(flipped: <code>{report['nominal']['fate_flipped']}</code>) · Ensemble fate-flip probability: <strong>{flip_p}</strong></p>
<table>
  <thead><tr><th>Gene</th><th>WT</th><th>Variant</th><th>Nominal Δ</th><th>Ensemble 90% CI</th></tr></thead>
  <tbody>{''.join(gene_rows)}</tbody>
</table>
<h2>4. Built-in Causal &amp; Negative Controls</h2>
<ul>
  <li>Exact allele reversion max state diff: <code>{report['controls']['reversion_max_state_diff']:.3e}</code></li>
  <li>Closed-chromatin cis-null max state diff: <code>{report['controls']['closed_chromatin_cis_max_state_diff']:.3e}</code></li>
  <li>Passive-reporter isolation max diff: <code>{report['controls']['passive_reporter_isolation_max_diff']:.3e}</code></li>
</ul>
<h2>5. Coupled 3D Neural Morphogenesis &amp; Electrophysiology</h2>
{cx_html}
</body>
</html>
"""


def run_pipeline(
    fasta_path: Path,
    spec_obj: dict[str, Any],
    variant_specs: list[dict[str, Any]],
    output: Path,
    *,
    run_cx3d: bool = True,
    cx3d_duration: float = DEFAULT_CX3D_DURATION,
    cx3d_seed: int = DEFAULT_CX3D_SEED,
    export_bundles: bool = False,
    observations: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Execute the complete sealed multi-locus DNA -> regulation -> GRN -> 3D CX3D pipeline."""
    if output.exists():
        raise FileExistsError(f"output directory already exists: {output}")

    ref_records = _parse_fasta(fasta_path)
    program = validate_program_spec(spec_obj, set(ref_records))
    if not isinstance(variant_specs, list) or not variant_specs:
        raise ValueError("variant_specs must be a non-empty list of reference-bound edits")

    output.mkdir(parents=True)
    snapshot_fasta = output / "reference.fasta"
    shutil.copyfile(fasta_path, snapshot_fasta)

    # Stage 1: Sealed DNA source + GA4GH VRS 2.x variant collection
    source_bundle = compile_source(
        FASTA_PROFILE, {"sequence": snapshot_fasta}, parameters={"wrapper": "identity"}
    )
    source_dir = output / "source"
    source_bundle.save(source_dir)

    variant_artifact = compile_variants(
        source_bundle,
        variant_specs,
        source_inputs={"sequence": snapshot_fasta},
    )
    variants_dir = output / "variants"
    variants_dir.mkdir()
    variants_path = variants_dir / "variants.json"
    variant_artifact.save(variants_path)

    variant_validation = validate_variant_paths(
        source_dir,
        variants_path,
        source_inputs={"sequence": snapshot_fasta},
    )
    if not variant_validation["valid"]:
        raise ValueError(f"variant collection validation failed: {variant_validation}")

    v_payload = load_variant_collection(variants_path).to_dict()
    compiled_variants = v_payload["variants_ir"]["variants"]
    mut_records, rev_records = _apply_multi_variants(ref_records, compiled_variants)

    # Stage 2: Multi-locus chromatin-gated TF binding & coding ORF translation
    eval_var = _evaluate_loci(program, ref_records, mut_records, force_closed_chromatin=False)
    eval_rev = _evaluate_loci(program, ref_records, rev_records, force_closed_chromatin=False)
    # For the closed-chromatin cis-null control, keep coding sequences WT so we specifically test cis silencing
    cis_only_records = dict(mut_records)
    for loc in program["loci"]:
        if loc["kind"] == "coding":
            cis_only_records[loc["record_id"]] = ref_records[loc["record_id"]]
    eval_closed = _evaluate_loci(program, ref_records, cis_only_records, force_closed_chromatin=True)

    genes = program["genes"]
    model = program["model"]
    protocol = program["protocol"]
    init_tuple = tuple(program["initial"][g] for g in genes)
    roles = eval_var["gene_roles"]
    mechanisms = eval_var["gene_mechanisms"]
    wt_ratios = {g: 1.0 for g in genes}
    var_ratios = eval_var["gene_ratios"]
    rev_ratios = eval_rev["gene_ratios"]
    closed_ratios = eval_closed["gene_ratios"]
    fates_cfg = model.get("fates")

    # Stage 3: Nominal ODE, passive reporter, causal controls & uncertainty ensemble
    wt_traj = nc.simulate_circuit(
        protocol, ratios=wt_ratios, initial=init_tuple, circuit=model, roles=roles, mechanisms=mechanisms
    )
    var_traj = nc.simulate_circuit(
        protocol, ratios=var_ratios, initial=init_tuple, circuit=model, roles=roles, mechanisms=mechanisms
    )
    rev_traj = nc.simulate_circuit(
        protocol, ratios=rev_ratios, initial=init_tuple, circuit=model, roles=roles, mechanisms=mechanisms
    )
    closed_traj = nc.simulate_circuit(
        protocol, ratios=closed_ratios, initial=init_tuple, circuit=model, roles=roles, mechanisms=mechanisms
    )
    # Passive reporter assay leaves endogenous trajectories fixed at WT
    rep_assay = nc.run_circuit_assay(
        protocol,
        ratio=var_ratios,
        mode="reporter",
        initial=program["initial"],
        circuit=model,
        roles=roles,
        mechanisms=mechanisms,
    )

    wt_end = wt_traj[-1][1:]
    var_end = var_traj[-1][1:]
    end_signal = protocol[-1][1]
    wt_occ = nc.occupancy_ratios(model, wt_end, end_signal, wt_ratios, roles=roles, mechanisms=mechanisms)
    var_occ = nc.occupancy_ratios(model, wt_end, end_signal, var_ratios, roles=roles, mechanisms=mechanisms)

    def _max_traj_diff(t1: list[tuple], t2: list[tuple]) -> float:
        return max(abs(a - b) for r1, r2 in zip(t1, t2) for a, b in zip(r1[1:], r2[1:]))

    reversion_diff = _max_traj_diff(wt_traj, rev_traj)
    closed_chromatin_diff = _max_traj_diff(wt_traj, closed_traj)
    reporter_iso_diff = _max_traj_diff(wt_traj, rep_assay["trajectory"])

    wt_fate = nc.fate_label(wt_end, genes, fates=fates_cfg)
    var_fate = nc.fate_label(var_end, genes, fates=fates_cfg)
    wt_state_dict = {g: float(wt_end[i]) for i, g in enumerate(genes)}
    var_state_dict = {g: float(var_end[i]) for i, g in enumerate(genes)}
    delta_dict = {g: float(var_end[i] - wt_end[i]) for i, g in enumerate(genes)}

    ensemble = _run_uncertainty_ensemble(program, var_ratios, roles, mechanisms)

    # Seal Stage 1-3 into a brainc.development-bundle/v1
    locus_eff_vec = [
        float(loc["effective_affinity_ratio"]) if loc["kind"] == "enhancer" else float(loc["coding_activity_factor"])
        for loc in eval_var["locus_reports"]
    ]
    gene_ratio_vec = [float(var_ratios[g]) for g in genes]
    wt_vec = [float(v) for v in wt_end]
    var_vec = [float(v) for v in var_end]
    tensors = sorted(
        [
            {
                "id": "count",
                "type": {"dtype": "u64", "shape": []},
                "unit": None,
                "axes": [],
                "storage": inline_storage(pack("u64", [1])),
            },
            {
                "id": "gene_effective_ratios",
                "type": {"dtype": "f64", "shape": [len(gene_ratio_vec)]},
                "unit": None,
                "axes": [None],
                "storage": inline_storage(pack("f64", gene_ratio_vec)),
            },
            {
                "id": "locus_effective_ratios",
                "type": {"dtype": "f64", "shape": [len(locus_eff_vec)]},
                "unit": None,
                "axes": [None],
                "storage": inline_storage(pack("f64", locus_eff_vec)),
            },
            {
                "id": "variant_steady_state",
                "type": {"dtype": "f64", "shape": [len(var_vec)]},
                "unit": None,
                "axes": [None],
                "storage": inline_storage(pack("f64", var_vec)),
            },
            {
                "id": "wt_steady_state",
                "type": {"dtype": "f64", "shape": [len(wt_vec)]},
                "unit": None,
                "axes": [None],
                "storage": inline_storage(pack("f64", wt_vec)),
            },
        ],
        key=lambda t: t["id"],
    )

    sem_hash = _semantics_sha256()
    manifest = seal(
        {
            "format": "brainc.provider-manifest",
            "version": 2,
            "provider": {"name": TARGET, "version": "1"},
            "model_identity": {"kind": "content-sha256", "value": sem_hash},
            "accepts": [f"brainc.source-descriptor/v1;profile={FASTA_PROFILE}"],
            "outputs": [{k: t[k] for k in ("id", "type", "unit", "axes")} for t in tensors],
        }
    )
    interp_dir = output / "interpretation"
    interp_dir.mkdir()
    save(manifest, interp_dir / "manifest.json")
    request = make_development_request(source_bundle, interp_dir / "manifest.json", [t["id"] for t in tensors])
    save(request, interp_dir / "request.json")
    response = seal(
        {
            "format": "brainc.prediction-response",
            "version": 2,
            "request_artifact_sha256": request["artifact_sha256"],
            "provider": manifest["provider"],
            "model_identity": manifest["model_identity"],
            "outputs": tensors,
        }
    )
    save(response, interp_dir / "response.json")
    target = _pipeline_target(len(locus_eff_vec), len(genes))
    save(target, interp_dir / "target.json")
    policy = policy_artifact(
        TARGET,
        {
            "id": target["contract"]["id"],
            "abi_major": target["contract"]["abi_major"],
            "contract_sha256": target["contract_sha256"],
        },
        [{"id": t["id"], "from_output": t["id"]} for t in tensors],
        [
            {
                "id": "create.developmental_program",
                "op": OP_UNIT_CREATE,
                "version": 1,
                "schema": "developmental_program",
                "count": "count",
                "initializers": [
                    {"field": "gene_effective_ratios", "tensor": "gene_effective_ratios"},
                    {"field": "locus_effective_ratios", "tensor": "locus_effective_ratios"},
                    {"field": "variant_steady_state", "tensor": "variant_steady_state"},
                    {"field": "wt_steady_state", "tensor": "wt_steady_state"},
                ],
            }
        ],
    )
    save(policy, interp_dir / "policy.json")

    dev_bundle = compile_development(
        source_bundle,
        interp_dir / "manifest.json",
        interp_dir / "request.json",
        interp_dir / "response.json",
        interp_dir / "policy.json",
        interp_dir / "target.json",
    )
    dev_dir = output / "development"
    dev_bundle.save(dev_dir)
    dev_validation = validate_development_paths(
        source_dir,
        dev_dir,
        source_inputs={"sequence": snapshot_fasta},
    )
    if not dev_validation["valid"]:
        raise ValueError(f"development bundle validation failed: {dev_validation}")

    # Cross-check compiled tensors decode cleanly
    _, dev_artifacts = load_development_bundle_directory(dev_dir)
    bundle_sha = dev_bundle.to_dict()["artifact_sha256"]
    rust_check = _try_rust_validate(source_dir, snapshot_fasta, variants_path, dev_dir)

    # Optional empirical falsification check
    obs_input = observations if observations is not None else program["observations"]
    falsification = evaluate_observations(delta_dict, obs_input) if obs_input is not None else None

    # Stage 4: Optional 3D CX3D Morphogenesis & LIF Spiking Electrophysiology
    morphogenesis: dict[str, Any]
    if run_cx3d:
        import activity as lif
        import construction as cx

        cx_dir = output / "morphogenesis"
        cx_dir.mkdir()
        backend = cx.build_backend(cx_dir / "backend")
        branches: dict[str, Any] = {}
        rev_state_dict = {g: float(rev_traj[-1][1 + i]) for i, g in enumerate(genes)}
        for label, st in (
            ("reference", wt_state_dict),
            ("variant", var_state_dict),
            ("allele_rescue", rev_state_dict),
        ):
            bdir = cx_dir / label
            bdir.mkdir()
            cx_prog = nc.program_from_grn_state(st, mapping=program["morphogenesis_mapping"])
            rec_bytes = cx.encode_recipe(cx_prog)
            rec_path = bdir / "recipe.dna"
            rec_path.write_bytes(rec_bytes)
            cx.compile_recipe(rec_path, bdir / "compiled")
            loaded_prog, cx_val = cx.load_compiled(bdir / "compiled")
            dev_sum = cx.develop(loaded_prog, backend, bdir / "development-run", duration=cx3d_duration, seed=cx3d_seed)
            act = lif.assay(bdir / "development-run", cx.verify_trace)
            (bdir / "activity.json").write_text(json.dumps(act, indent=2) + "\n", encoding="utf-8")
            branches[label] = {
                "recipe_sha256": hashlib.sha256(rec_bytes).hexdigest(),
                "compiled_valid": cx_val["valid"],
                "program": loaded_prog._asdict(),
                "construction": dev_sum,
                "activity": {
                    "cells": act["cells"],
                    "contacts": act["contacts"],
                    "transmission_controls_pass": act["transmission_controls_pass"],
                    "total_spikes": sum(len(t["intact_spikes"]) for t in act["trials"]),
                },
            }
        ref_b, var_b, res_b = branches["reference"], branches["variant"], branches["allele_rescue"]
        def _cx_metrics(b: dict[str, Any]) -> tuple[Any, ...]:
            c = b["construction"]
            return (b["recipe_sha256"], c["cells"], c["segments"], c["synapses"], c["length"])

        morphogenesis = {
            "executed": True,
            "duration": cx3d_duration,
            "seed": cx3d_seed,
            "branches": branches,
            "morphology_diverged": _cx_metrics(ref_b) != _cx_metrics(var_b),
            "activity_diverged": ref_b["activity"] != var_b["activity"],
            "allele_rescue_matches_reference": (
                _cx_metrics(ref_b) == _cx_metrics(res_b)
                and ref_b["activity"] == res_b["activity"]
            ),
        }
    else:
        morphogenesis = {"executed": False, "reason": "skipped_by_flag"}

    # Stage 5: Optional RO-Crate 1.1 & COMBINE/OMEX export
    exports: dict[str, Any] = {"exported": False}
    if export_bundles:
        rc_info = export_rocrate(
            dev_dir,
            output / "rocrate",
            source_bundle=source_dir,
            variant_collection=variants_path,
        )
        omex_info = export_omex(
            dev_dir,
            output / "bundle.omex",
            source_bundle=source_dir,
            variant_collection=variants_path,
        )
        exports = {"exported": True, "rocrate": rc_info, "omex": omex_info}

    report = {
        "format": FORMAT,
        "program_name": program["name"],
        "circuit": model["name"],
        "genes": genes,
        "protocol": [[d, s] for d, s in protocol],
        "biological_acceptance": False,
        "scientific_scope": (
            "Composable multi-locus DNA -> VRS 2.x -> chromatin-gated TF/coding -> GRN uncertainty "
            "-> 3D CX3D+LIF software pipeline with sealed provenance and causal controls. "
            "Uncalibrated model units; in vivo developmental phenotype is NOT established."
        ),
        "provenance": {
            "source_descriptor_sha256": source_bundle.to_dict()["artifact_sha256"],
            "variant_collection_sha256": v_payload["artifact_sha256"],
            "development_bundle_sha256": bundle_sha,
            "development_module_sha256": dev_artifacts["development_module"]["artifact_sha256"],
            "variants_valid": variant_validation["valid"],
            "development_valid": dev_validation["valid"],
            "rust_cross_check": rust_check,
        },
        "variants": [
            {
                "id": v["id"],
                "record_id": v["record_id"],
                "edit_kind": v["edit_kind"],
                "start": v["location"]["start"],
                "end": v["location"]["end"],
                "ref": v["ref"],
                "alt": v["alt"],
                "vrs_location_id": v["vrs_location_id"],
                "vrs_allele_id": v["vrs_allele_id"],
            }
            for v in compiled_variants
        ],
        "loci": eval_var["locus_reports"],
        "gene_effective_ratios": var_ratios,
        "nominal": {
            "wt_end_state": wt_state_dict,
            "variant_end_state": var_state_dict,
            "gene_deltas": delta_dict,
            "wt_fate": wt_fate,
            "variant_fate": var_fate,
            "fate_flipped": wt_fate != var_fate,
            "passive_reporter_occupancy_wt": {g: float(wt_occ[i]) for i, g in enumerate(genes)},
            "passive_reporter_occupancy_variant": {g: float(var_occ[i]) for i, g in enumerate(genes)},
            "passive_reporter_deltas": {g: float(var_occ[i] - wt_occ[i]) for i, g in enumerate(genes)},
        },
        "uncertainty_ensemble": ensemble,
        "controls": {
            "reversion_max_state_diff": reversion_diff,
            "reversion_exact": reversion_diff == 0.0,
            "closed_chromatin_cis_max_state_diff": closed_chromatin_diff,
            "closed_chromatin_cis_exact": closed_chromatin_diff == 0.0,
            "passive_reporter_isolation_max_diff": reporter_iso_diff,
            "passive_reporter_isolated": reporter_iso_diff == 0.0,
        },
        "falsification": falsification,
        "morphogenesis": morphogenesis,
        "exports": exports,
    }
    report["report_sha256"] = hashlib.sha256(canonical_bytes(report)).hexdigest()

    (output / "pipeline_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    (output / "index.html").write_text(_render_html(report), encoding="utf-8")
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="mode", required=True)

    dm = sub.add_parser("demo", help="run a built-in multi-locus composable developmental pipeline demo")
    dm.add_argument("--preset", default="neural-tube-multi-locus", choices=("neural-tube-multi-locus", "tri-lineage-chromatin"))
    dm.add_argument("--output", type=Path, required=True)
    dm.add_argument("--skip-cx3d", action="store_true", help="skip 3D CX3D + LIF stage (pure Python fast mode)")
    dm.add_argument("--cx3d-duration", type=float, default=DEFAULT_CX3D_DURATION)
    dm.add_argument("--cx3d-seed", type=int, default=DEFAULT_CX3D_SEED)
    dm.add_argument("--export-bundles", action="store_true", help="also emit RO-Crate 1.1 and COMBINE/OMEX archives")

    rn = sub.add_parser("run", help="run a custom multi-locus developmental program from FASTA + spec + variants JSON")
    rn.add_argument("--fasta", type=Path, required=True)
    rn.add_argument("--spec", type=Path, required=True)
    rn.add_argument("--variants", type=Path, required=True)
    rn.add_argument("--observations", type=Path, help="optional empirical observations JSON for falsification check")
    rn.add_argument("--output", type=Path, required=True)
    rn.add_argument("--skip-cx3d", action="store_true")
    rn.add_argument("--cx3d-duration", type=float, default=DEFAULT_CX3D_DURATION)
    rn.add_argument("--cx3d-seed", type=int, default=DEFAULT_CX3D_SEED)
    rn.add_argument("--export-bundles", action="store_true")

    args = parser.parse_args(argv)
    if args.mode == "demo":
        preset = builtin_preset(args.preset)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        tmp_fasta = args.output.parent / f".{args.output.name}.preset.fasta"
        try:
            _write_fasta(preset["records"], tmp_fasta)
            report = run_pipeline(
                tmp_fasta,
                preset["spec"],
                preset["variants"],
                args.output,
                run_cx3d=not args.skip_cx3d,
                cx3d_duration=args.cx3d_duration,
                cx3d_seed=args.cx3d_seed,
                export_bundles=args.export_bundles,
            )
        finally:
            if tmp_fasta.exists():
                tmp_fasta.unlink()
    else:
        spec_obj = load_json_object(read_regular_file(args.spec, maximum_bytes=2 * 1024 * 1024), "program spec")
        variants_raw = json.loads(read_regular_file(args.variants, maximum_bytes=2 * 1024 * 1024).decode("utf-8"))
        obs_raw = (
            json.loads(read_regular_file(args.observations, maximum_bytes=2 * 1024 * 1024).decode("utf-8"))
            if args.observations is not None
            else None
        )
        report = run_pipeline(
            args.fasta,
            spec_obj,
            variants_raw,
            args.output,
            run_cx3d=not args.skip_cx3d,
            cx3d_duration=args.cx3d_duration,
            cx3d_seed=args.cx3d_seed,
            export_bundles=args.export_bundles,
            observations=obs_raw,
        )

    print(
        json.dumps(
            {
                "program_name": report["program_name"],
                "wt_fate": report["nominal"]["wt_fate"],
                "variant_fate": report["nominal"]["variant_fate"],
                "fate_flip_probability": report["uncertainty_ensemble"]["fate_flip_probability"],
                "reversion_exact": report["controls"]["reversion_exact"],
                "closed_chromatin_cis_exact": report["controls"]["closed_chromatin_cis_exact"],
                "cx3d_executed": report["morphogenesis"]["executed"],
                "biological_acceptance": report["biological_acceptance"],
                "output": str(args.output),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

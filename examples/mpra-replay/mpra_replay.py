"""Replayable MPRA variant-effect analysis on Kreimer et al. 2022 (GSE188264).

Demonstrates the Epigenesis provenance pipeline for variant-effect claims:
1. Compiles verified GSE188264 locus references under fasta-dna/v1 (brainc.source-descriptor/v1).
2. Binds tested variant alleles to exact reference loci as GA4GH VRS 2.x identifiers
   (brainc.variant-collection/v1), verified by validate_variant_paths.
3. Compiles two defensible MPRA calling parameter sets (Permissive Mean-Delta vs.
   Strict Late-QC) into separate sealed Development Modules (brainc.development-bundle/v1),
   each independently replayed via validate_development_paths.
4. Reports exact VRS 2.x variant calls that flip between parameter sets, optionally
   exporting deterministic RO-Crate 1.1 / COMBINE OMEX provenance bundles.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys
from typing import Any

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "examples/neural-regulation"))

from brainc._canonical import canonical_bytes, digest
from brainc._io import load_json_object, read_regular_file
from brainc.development_bundle import compile_development
from brainc.export import export_omex, export_rocrate
from brainc.source import FASTA_PROFILE, compile_source
from brainc.v2 import inline_storage, make_development_request, pack, policy_artifact, save, target_artifact
from brainc.v2._common import seal
from brainc.v2.target import DEV_DOMAIN, OP_UNIT_CREATE
from brainc.validator_development import load_development_bundle_directory, validate_development_paths
from brainc.validator_variant import validate_variant_paths
from brainc.variant import compile_variants
import regulation as reg


TARGET_ID = "org.epigenesis.example.mpra-replay"

PARAMETER_SETS: dict[str, dict[str, Any]] = {
    "permissive_mean_fc": {
        "name": "permissive_mean_fc",
        "description": "All-timepoint mean |asinh log2 FC| >= 0.50 with zero WT floor",
        "timepoints_hours": [0, 3, 6, 12, 24, 48, 72],
        "time_indices": [0, 1, 2, 3, 4, 5, 6],
        "min_wt_mean_alpha": 0.0,
        "effect_threshold": 0.50,
        "require_sign_consistency": False,
    },
    "strict_late_qc": {
        "name": "strict_late_qc",
        "description": "Late-window (24-72h) mean |asinh log2 FC| >= 1.00, sign-consistent, WT floor >= 1.00",
        "timepoints_hours": [24, 48, 72],
        "time_indices": [4, 5, 6],
        "min_wt_mean_alpha": 1.00,
        "effect_threshold": 1.00,
        "require_sign_consistency": True,
    },
}

CALL_LABELS = {0: "silent", 1: "active", 2: "qc_filtered"}


def _minimal_edit_span(reference: str, variant: str) -> tuple[int, int, str, str]:
    start = 0
    limit = min(len(reference), len(variant))
    while start < limit and reference[start] == variant[start]:
        start += 1
    suffix = 0
    while (
        suffix < (len(reference) - start)
        and suffix < (len(variant) - start)
        and reference[len(reference) - 1 - suffix] == variant[len(variant) - 1 - suffix]
    ):
        suffix += 1
    end = len(reference) - suffix
    var_end = len(variant) - suffix
    return start, end, reference[start:end], variant[start:var_end]


def _target_contract(model_digest: str, count: int) -> dict[str, Any]:
    return target_artifact(
        {
            "id": TARGET_ID,
            "abi_major": 1,
            "opsets": [{"domain": DEV_DOMAIN, "version": 1}],
            "unit_schemas": [
                {
                    "id": "mpra-call-summary",
                    "fields": sorted(
                        [
                            {
                                "id": "variant-count",
                                "type": {"dtype": "u64", "shape": []},
                                "unit": None,
                                "mutability": "constant",
                                "numeric": {"kind": "integer"},
                            },
                            {
                                "id": "effect-scores",
                                "type": {"dtype": "f64", "shape": [count]},
                                "unit": None,
                                "mutability": "constant",
                                "numeric": {"kind": "float"},
                            },
                            {
                                "id": "call-codes",
                                "type": {"dtype": "u64", "shape": [count]},
                                "unit": None,
                                "mutability": "constant",
                                "numeric": {"kind": "integer"},
                            },
                        ],
                        key=lambda item: item["id"],
                    ),
                }
            ],
            "edge_schemas": [],
            "ports": [],
            "rules": [],
        }
    )


def _evaluate_parameter_set(
    rows: list[reg.Row],
    scale: float,
    params: dict[str, Any],
) -> tuple[list[float], list[int], dict[str, int]]:
    indices = params["time_indices"]
    min_wt = params["min_wt_mean_alpha"]
    threshold = params["effect_threshold"]
    consistent = params["require_sign_consistency"]

    effects: list[float] = []
    calls: list[int] = []
    counts = {"silent": 0, "active": 0, "qc_filtered": 0}
    for row in rows:
        deltas = reg.target(row, scale)
        window_deltas = [deltas[i] for i in indices]
        mean_delta = statistics.mean(window_deltas)
        wt_window_means = [
            statistics.mean(wt[i] for wt in row.wt_alphas)
            for i in indices
        ]
        if min(wt_window_means) < min_wt:
            code = 2  # qc_filtered
        else:
            sign_ok = True
            if consistent:
                signs = {1 if d > 0 else (-1 if d < 0 else 0) for d in window_deltas}
                sign_ok = len(signs) == 1 and 0 not in signs
            code = 1 if (abs(mean_delta) >= threshold and sign_ok) else 0
        effects.append(mean_delta)
        calls.append(code)
        counts[CALL_LABELS[code]] += 1
    return effects, calls, counts


def _compile_parameter_bundle(
    source_dir: Path,
    fasta_path: Path,
    output_dir: Path,
    params: dict[str, Any],
    variant_artifact_sha256: str,
    provenance_sha256: str,
    effects: list[float],
    calls: list[int],
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    count = len(calls)
    param_record = {
        "format": "epigenesis.mpra-parameter-set",
        "version": 1,
        "parameters": params,
        "variant_collection_sha256": variant_artifact_sha256,
        "source_provenance_sha256": provenance_sha256,
    }
    param_bytes = canonical_bytes(param_record) + b"\n"
    (output_dir / "parameters.json").write_bytes(param_bytes)
    model_digest = hashlib.sha256(param_bytes).hexdigest()

    tensors = sorted(
        [
            {
                "id": "unit_count",
                "type": {"dtype": "u64", "shape": []},
                "unit": None,
                "axes": [],
                "storage": inline_storage(pack("u64", [1])),
            },
            {
                "id": "variant-count",
                "type": {"dtype": "u64", "shape": []},
                "unit": None,
                "axes": [],
                "storage": inline_storage(pack("u64", [count])),
            },
            {
                "id": "effect-scores",
                "type": {"dtype": "f64", "shape": [count]},
                "unit": None,
                "axes": [None],
                "storage": inline_storage(pack("f64", effects)),
            },
            {
                "id": "call-codes",
                "type": {"dtype": "u64", "shape": [count]},
                "unit": None,
                "axes": [None],
                "storage": inline_storage(pack("u64", calls)),
            },
        ],
        key=lambda item: item["id"],
    )
    manifest = seal(
        {
            "format": "brainc.provider-manifest",
            "version": 2,
            "provider": {"name": TARGET_ID, "version": "1"},
            "model_identity": {"kind": "content-sha256", "value": model_digest},
            "accepts": [f"brainc.source-descriptor/v1;profile={FASTA_PROFILE}"],
            "outputs": [{k: t[k] for k in ("id", "type", "unit", "axes")} for t in tensors],
        }
    )
    save(manifest, output_dir / "manifest.json")
    request = make_development_request(source_dir, output_dir / "manifest.json", [t["id"] for t in tensors])
    save(request, output_dir / "request.json")
    save(
        seal(
            {
                "format": "brainc.prediction-response",
                "version": 2,
                "request_artifact_sha256": request["artifact_sha256"],
                "provider": manifest["provider"],
                "model_identity": manifest["model_identity"],
                "outputs": tensors,
            }
        ),
        output_dir / "response.json",
    )
    target = _target_contract(model_digest, count)
    save(target, output_dir / "target.json")
    ops = [
        {
            "id": "create.mpra",
            "op": OP_UNIT_CREATE,
            "version": 1,
            "schema": "mpra-call-summary",
            "count": "unit_count",
            "initializers": [
                {"field": fid, "tensor": fid}
                for fid in ("call-codes", "effect-scores", "variant-count")
            ],
        }
    ]
    save(
        policy_artifact(
            TARGET_ID,
            {"id": TARGET_ID, "abi_major": 1, "contract_sha256": target["contract_sha256"]},
            [{"id": t["id"], "from_output": t["id"]} for t in tensors],
            ops,
        ),
        output_dir / "policy.json",
    )
    dev_dir = output_dir / "development"
    dev_bundle = compile_development(
        source_dir,
        output_dir / "manifest.json",
        output_dir / "request.json",
        output_dir / "response.json",
        output_dir / "policy.json",
        output_dir / "target.json",
    )
    dev_bundle.save(dev_dir)
    val = validate_development_paths(source_dir, dev_dir, source_inputs={"sequence": fasta_path})
    if not val["valid"]:
        raise RuntimeError(f"development validation failed for {params['name']}: {val}")
    return {
        "parameter_set": params["name"],
        "parameters_sha256": model_digest,
        "development_bundle_sha256": dev_bundle.to_dict()["artifact_sha256"],
        "development_module_sha256": val["result"]["module_artifact_sha256"],
        "validation_report_sha256": val["report_sha256"],
        "valid": True,
    }


def run_mpra_replay(
    output: Path,
    *,
    split: str = "test",
    limit: int | None = None,
    export_bundles: bool = True,
) -> dict[str, Any]:
    """Run replayable dual-parameter MPRA analysis on verified GSE188264 rows."""
    output.mkdir(parents=True, exist_ok=True)
    all_rows, qc = reg.load_data()
    train_rows = [r for r in all_rows if r.split == "train"]
    scale = reg.training_scale(train_rows)

    selected = all_rows if split == "all" else [r for r in all_rows if r.split == split]
    if limit is not None:
        selected = selected[:limit]

    # Filter out exact no-op duplicates if any exist so every variant binds a non-noop VRS 2.x Allele
    eligible: list[reg.Row] = []
    skipped_noop = 0
    for r in selected:
        if r.reference == r.sequence:
            skipped_noop += 1
        else:
            eligible.append(r)
    if not eligible:
        raise ValueError("no non-noop variants selected for MPRA replay")

    loci_seqs = sorted({(r.locus, r.reference) for r in eligible}, key=lambda x: x[0])
    fasta_path = output / "wt_loci.fasta"
    fasta_path.write_text(
        "".join(f">{locus}\n{seq}\n" for locus, seq in loci_seqs),
        encoding="utf-8",
    )
    source_dir = output / "source"
    source_bundle = compile_source(
        FASTA_PROFILE,
        {"sequence": fasta_path},
        parameters={"wrapper": "identity"},
    )
    source_bundle.save(source_dir)

    variant_specs = []
    for r in eligible:
        start, end, ref_s, alt_s = _minimal_edit_span(r.reference, r.sequence)
        variant_specs.append(
            {
                "id": r.id,
                "record_id": r.locus,
                "start": start,
                "end": end,
                "ref": ref_s,
                "alt": alt_s,
            }
        )
    vc = compile_variants(source_dir, variant_specs)
    vc_path = output / "variants.json"
    vc.save(vc_path)
    vc_val = validate_variant_paths(source_dir, vc_path)
    if not vc_val["valid"]:
        raise RuntimeError(f"MPRA variant collection failed validation: {vc_val}")

    v_entries = vc.to_dict()["variants_ir"]["variants"]
    provenance_sha = hashlib.sha256(
        canonical_bytes(qc["source_provenance"])
    ).hexdigest()

    eff_a, calls_a, counts_a = _evaluate_parameter_set(eligible, scale, PARAMETER_SETS["permissive_mean_fc"])
    eff_b, calls_b, counts_b = _evaluate_parameter_set(eligible, scale, PARAMETER_SETS["strict_late_qc"])

    seal_a = _compile_parameter_bundle(
        source_dir,
        fasta_path,
        output / "param_a_permissive_mean_fc",
        PARAMETER_SETS["permissive_mean_fc"],
        vc.artifact_sha256,
        provenance_sha,
        eff_a,
        calls_a,
    )
    seal_b = _compile_parameter_bundle(
        source_dir,
        fasta_path,
        output / "param_b_strict_late_qc",
        PARAMETER_SETS["strict_late_qc"],
        vc.artifact_sha256,
        provenance_sha,
        eff_b,
        calls_b,
    )

    flipped_calls = []
    transition_counts: dict[str, int] = {}
    for row, entry, ea, ca, eb, cb in zip(eligible, v_entries, eff_a, calls_a, eff_b, calls_b):
        la, lb = CALL_LABELS[ca], CALL_LABELS[cb]
        key = f"{la}->{lb}"
        transition_counts[key] = transition_counts.get(key, 0) + 1
        if ca != cb:
            flipped_calls.append(
                {
                    "id": row.id,
                    "locus": row.locus,
                    "kind": row.kind,
                    "vrs_allele_id": entry["vrs_allele_id"],
                    "vrs_location_id": entry["vrs_location_id"],
                    "edit_kind": entry["edit_kind"],
                    "permissive_call": la,
                    "permissive_effect": round(ea, 6),
                    "strict_call": lb,
                    "strict_effect": round(eb, 6),
                }
            )

    exports = {}
    if export_bundles:
        (output / "exports").mkdir(parents=True, exist_ok=True)
        exports["param_a_rocrate"] = export_rocrate(
            output / "param_a_permissive_mean_fc" / "development",
            output / "exports" / "param_a_rocrate",
            source_bundle=source_dir,
            variant_collection=vc_path,
        )["rocrate_metadata_sha256"]
        exports["param_b_omex"] = export_omex(
            output / "param_b_strict_late_qc" / "development",
            output / "exports" / "param_b_strict_late_qc.omex",
            source_bundle=source_dir,
            variant_collection=vc_path,
        )["archive_sha256"]

    report_core = {
        "format": "epigenesis.mpra-replay-comparison",
        "version": 1,
        "dataset": {
            "study": qc["source_provenance"]["study"],
            "accession": qc["source_provenance"]["accession"],
            "doi": qc["source_provenance"]["doi"],
            "split": split,
            "loci_count": len(loci_seqs),
            "variants_evaluated": len(eligible),
            "skipped_noop": skipped_noop,
            "training_scale": round(scale, 9),
            "provenance_sha256": provenance_sha,
        },
        "sealed_provenance": {
            "source_artifact_sha256": source_bundle.to_dict()["artifact_sha256"],
            "variant_collection_sha256": vc.artifact_sha256,
            "variant_validation_report_sha256": vc_val["report_sha256"],
            "parameter_set_a": seal_a,
            "parameter_set_b": seal_b,
            "exports": exports,
        },
        "summary": {
            "counts_permissive_mean_fc": counts_a,
            "counts_strict_late_qc": counts_b,
            "total_flipped_calls": len(flipped_calls),
            "flip_fraction": round(len(flipped_calls) / len(eligible), 6),
            "transition_matrix": dict(sorted(transition_counts.items())),
        },
        "flipped_examples_top25": flipped_calls[:25],
    }
    report = {**report_core, "report_sha256": digest(report_core)}
    (output / "mpra_comparison.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    (output / "flipped_variants.json").write_text(
        json.dumps(flipped_calls, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--split", choices=("test", "validation", "train", "all"), default="test")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--no-export", action="store_true")
    parser.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args()
    rep = run_mpra_replay(
        args.output,
        split=args.split,
        limit=args.limit,
        export_bundles=not args.no_export,
    )
    print(json.dumps({k: v for k, v in rep.items() if k != "flipped_examples_top25"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

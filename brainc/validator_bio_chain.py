"""Standalone raw-DNA-to-feature-graph validation chain.

This module snapshots every external input once, independently replays raw
sequence and GFF3 into the supplied biological artifacts, independently
replays those artifacts into the supplied feature-graph bundle, and emits one
cryptographically sealed report.  It imports only the two independent
validators and the Python standard library.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import stat
import sys
import tempfile
from typing import Any

import brainc.validator_bio as validator_bio
import brainc.validator_bio_graph as validator_bio_graph


REPORT_FORMAT = "brainc.bio.feature-graph-chain-validation-report"
REPORT_VERSION = 1
VALIDATOR = {
    "name": "brainc-independent-bio-chain-validator",
    "version": "0.1.0",
    "strategy": "two-stage-independent-replay",
}
CHECKS = [
    "raw-sequence-to-sequence-collection-replay",
    "raw-gff3-to-bio-ir-replay",
    "bio-ir-to-feature-graph-replay",
    "cross-stage-artifact-bindings",
    "cross-stage-cardinalities",
    "canonical-report-seals",
]

READ_CHUNK_BYTES = 64 * 1024
_SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
_GRAPH_REPORT_CHECKS = [
    "sequence-collection-identity",
    "bio-ir-structure-and-binding",
    "published-normative-backend-specification",
    "complete-normalized-bio-ir-sidecar-and-dictionaries",
    "numeric-tensor-encoding",
    "explicit-relationship-cardinality",
    "target-policy-provider-chain",
    "development-module-linkage-and-budgets",
    "compilation-record-and-bundle-bindings",
]


class BioChainValidationError(ValueError):
    """An input snapshot or cross-stage chain invariant failed."""


def _fail(detail: str) -> BioChainValidationError:
    return BioChainValidationError(f"BIOCHAINVAL001: {detail}")


def _keys(value: Any, expected: set[str], label: str) -> dict[str, Any]:
    if type(value) is not dict:
        raise _fail(f"{label} must be an object")
    missing = sorted(expected - value.keys())
    extra = sorted(value.keys() - expected)
    if missing or extra:
        raise _fail(
            f"{label} keys invalid; missing={missing or 'none'}, "
            f"unknown={extra or 'none'}"
        )
    return value


def _integer(value: Any, label: str) -> int:
    if type(value) is not int or value < 0 or value > 2**53 - 1:
        raise _fail(f"{label} must be a nonnegative safe integer")
    return value


def _sha256(value: Any, label: str) -> str:
    if type(value) is not str or _SHA256_RE.fullmatch(value) is None:
        raise _fail(f"{label} must be a lowercase SHA-256 digest")
    return value


def _read_regular(path: str | Path, label: str, maximum_bytes: int) -> bytes:
    """Read one stable regular-file snapshot without following a final link."""

    source = Path(path)
    descriptor = -1
    try:
        inspected = source.lstat()
        if stat.S_ISLNK(inspected.st_mode) or not stat.S_ISREG(inspected.st_mode):
            raise _fail(f"{label} must be a regular non-linked file")
        if inspected.st_size > maximum_bytes:
            raise _fail(f"{label} exceeds {maximum_bytes} bytes")
        flags = os.O_RDONLY
        flags |= getattr(os, "O_CLOEXEC", 0)
        flags |= getattr(os, "O_NOFOLLOW", 0)
        flags |= getattr(os, "O_NONBLOCK", 0)
        descriptor = os.open(source, flags)
        opened = os.fstat(descriptor)
        if not stat.S_ISREG(opened.st_mode):
            raise _fail(f"{label} must be a regular non-linked file")
        if (opened.st_dev, opened.st_ino) != (inspected.st_dev, inspected.st_ino):
            raise _fail(f"{label} changed while being opened")
        if opened.st_size > maximum_bytes:
            raise _fail(f"{label} exceeds {maximum_bytes} bytes")

        chunks: list[bytes] = []
        total = 0
        while True:
            request = min(READ_CHUNK_BYTES, maximum_bytes + 1 - total)
            chunk = os.read(descriptor, request)
            if not chunk:
                break
            chunks.append(chunk)
            total += len(chunk)
            if total > maximum_bytes:
                raise _fail(f"{label} exceeds {maximum_bytes} bytes")
        finished = os.fstat(descriptor)
        stable = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
        if any(getattr(opened, field) != getattr(finished, field) for field in stable):
            raise _fail(f"{label} changed while being read")
        if total != opened.st_size:
            raise _fail(f"{label} did not match its inspected length")
        return b"".join(chunks)
    except BioChainValidationError:
        raise
    except OSError as failure:
        raise _fail(f"cannot read {label}: {failure}") from failure
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise _fail(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def _reject_constant(token: str) -> None:
    raise _fail(f"non-finite JSON number {token!r}")


def _json_object(raw: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(
            raw.decode("utf-8"),
            object_pairs_hook=_reject_duplicates,
            parse_constant=_reject_constant,
        )
    except BioChainValidationError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError, RecursionError) as failure:
        raise _fail(f"invalid {label} JSON: {failure}") from failure
    if type(value) is not dict:
        raise _fail(f"{label} JSON must be an object")
    return value


def _validate_graph_report(report: Any) -> dict[str, Any]:
    root = _keys(
        report,
        {
            "format",
            "version",
            "stage",
            "valid",
            "inputs",
            "result",
            "checks",
            "report_sha256",
        },
        "feature-graph stage report",
    )
    if (
        root["format"] != "brainc.bio.feature-graph-validation-report"
        or type(root["version"]) is not int
        or root["version"] != 1
        or root["stage"] != "bio-ir-feature-graph-bridge"
        or root["valid"] is not True
        or root["checks"] != _GRAPH_REPORT_CHECKS
    ):
        raise _fail("feature-graph stage report identity or checks are invalid")
    inputs = _keys(
        root["inputs"],
        {
            "sequence_artifact_sha256",
            "bio_ir_artifact_sha256",
            "bundle_artifact_sha256",
            "backend_spec_sha256",
        },
        "feature-graph stage inputs",
    )
    for name, value in inputs.items():
        _sha256(value, f"feature-graph stage inputs.{name}")
    result = _keys(
        root["result"],
        {
            "operations",
            "tensor_bytes",
            "units",
            "edges",
            "attachments",
            "relationship_edges",
            "semantics_sha256",
            "module_sha256",
        },
        "feature-graph stage result",
    )
    for name in (
        "operations",
        "tensor_bytes",
        "units",
        "edges",
        "attachments",
        "relationship_edges",
    ):
        _integer(result[name], f"feature-graph stage result.{name}")
    _sha256(result["semantics_sha256"], "feature-graph stage result.semantics_sha256")
    _sha256(result["module_sha256"], "feature-graph stage result.module_sha256")
    expected = validator_bio_graph.digest(
        {key: value for key, value in root.items() if key != "report_sha256"}
    )
    if root["report_sha256"] != expected:
        raise _fail("feature-graph stage report seal is invalid")
    return root


def _expected_result(
    biological: dict[str, Any], graph: dict[str, Any]
) -> dict[str, Any]:
    bio_summary = biological["summary"]
    graph_result = graph["result"]
    return {
        "records": bio_summary["records"],
        "sequence_bases": bio_summary["sequence_bases"],
        "features": bio_summary["features"],
        "segments": bio_summary["segments"],
        "relationships": bio_summary["relationships"],
        "operations": graph_result["operations"],
        "tensor_bytes": graph_result["tensor_bytes"],
        "units": graph_result["units"],
        "edges": graph_result["edges"],
        "module_sha256": graph_result["module_sha256"],
    }


def _check_stage_bindings(
    biological: dict[str, Any], graph: dict[str, Any]
) -> None:
    bio_inputs = biological["inputs"]
    graph_inputs = graph["inputs"]
    if (
        bio_inputs["sequence_collection_artifact_sha256"]
        != graph_inputs["sequence_artifact_sha256"]
    ):
        raise _fail("sequence collection binding differs between validation stages")
    if bio_inputs["bio_artifact_sha256"] != graph_inputs["bio_ir_artifact_sha256"]:
        raise _fail("BioIR binding differs between validation stages")

    bio_summary = biological["summary"]
    graph_result = graph["result"]
    if bio_summary["features"] != graph_result["units"]:
        raise _fail("validated feature count differs from emitted unit count")
    if not (
        bio_summary["relationships"]
        == graph_result["edges"]
        == graph_result["relationship_edges"]
    ):
        raise _fail("validated relationship count differs from emitted edge count")


def validate_bio_feature_graph_chain(
    *,
    fasta_source: bytes,
    gff3_source: bytes,
    sequence_collection: dict[str, Any],
    bio_artifact: dict[str, Any],
    feature_graph_bundle: dict[str, Any],
    source_metadata: bytes | dict[str, Any] | None = None,
    record_id: str | None = None,
) -> dict[str, Any]:
    """Validate one in-memory raw-source-to-feature-graph transaction."""

    if type(fasta_source) is not bytes or type(gff3_source) is not bytes:
        raise _fail("raw sequence and GFF3 snapshots must be bytes")
    biological = validator_bio.validate_bio_chain(
        bio_artifact,
        sequence_collection,
        fasta_source=fasta_source,
        gff3_source=gff3_source,
        source_metadata=source_metadata,
        record_id=record_id,
    )
    graph = validator_bio_graph.validate_feature_graph_bundle(
        feature_graph_bundle,
        sequence_collection,
        bio_artifact,
    )
    validator_bio.validate_bio_report(biological)
    _validate_graph_report(graph)
    _check_stage_bindings(biological, graph)
    core = {
        "format": REPORT_FORMAT,
        "version": REPORT_VERSION,
        "validator": dict(VALIDATOR),
        "valid": True,
        "inputs": {
            "fasta_sha256": biological["inputs"]["fasta_sha256"],
            "gff3_sha256": biological["inputs"]["gff3_sha256"],
            "sequence_collection_artifact_sha256": biological["inputs"][
                "sequence_collection_artifact_sha256"
            ],
            "bio_artifact_sha256": biological["inputs"]["bio_artifact_sha256"],
            "bundle_artifact_sha256": graph["inputs"]["bundle_artifact_sha256"],
        },
        "checks": list(CHECKS),
        "stages": {
            "source_to_bio_ir": biological,
            "bio_ir_to_feature_graph": graph,
        },
        "result": _expected_result(biological, graph),
    }
    report = {**core, "report_sha256": validator_bio_graph.digest(core)}
    return validate_chain_report(report)


def validate_chain_report(report: dict[str, Any]) -> dict[str, Any]:
    """Validate a saved combined report's closed schema, stages, and seal."""

    root = _keys(
        report,
        {
            "format",
            "version",
            "validator",
            "valid",
            "inputs",
            "checks",
            "stages",
            "result",
            "report_sha256",
        },
        "combined validation report",
    )
    if (
        root["format"] != REPORT_FORMAT
        or type(root["version"]) is not int
        or root["version"] != REPORT_VERSION
        or root["validator"] != VALIDATOR
        or root["valid"] is not True
        or root["checks"] != CHECKS
    ):
        raise _fail("combined validation report identity or checks are invalid")
    inputs = _keys(
        root["inputs"],
        {
            "fasta_sha256",
            "gff3_sha256",
            "sequence_collection_artifact_sha256",
            "bio_artifact_sha256",
            "bundle_artifact_sha256",
        },
        "combined validation inputs",
    )
    for name, value in inputs.items():
        _sha256(value, f"combined validation inputs.{name}")
    stages = _keys(
        root["stages"],
        {"source_to_bio_ir", "bio_ir_to_feature_graph"},
        "combined validation stages",
    )
    biological = validator_bio.validate_bio_report(stages["source_to_bio_ir"])
    graph = _validate_graph_report(stages["bio_ir_to_feature_graph"])
    _check_stage_bindings(biological, graph)
    expected_inputs = {
        "fasta_sha256": biological["inputs"]["fasta_sha256"],
        "gff3_sha256": biological["inputs"]["gff3_sha256"],
        "sequence_collection_artifact_sha256": biological["inputs"][
            "sequence_collection_artifact_sha256"
        ],
        "bio_artifact_sha256": biological["inputs"]["bio_artifact_sha256"],
        "bundle_artifact_sha256": graph["inputs"]["bundle_artifact_sha256"],
    }
    if inputs != expected_inputs:
        raise _fail("combined report inputs do not match its sealed stage reports")
    result = _keys(
        root["result"],
        {
            "records",
            "sequence_bases",
            "features",
            "segments",
            "relationships",
            "operations",
            "tensor_bytes",
            "units",
            "edges",
            "module_sha256",
        },
        "combined validation result",
    )
    for name in result.keys() - {"module_sha256"}:
        _integer(result[name], f"combined validation result.{name}")
    _sha256(result["module_sha256"], "combined validation result.module_sha256")
    if result != _expected_result(biological, graph):
        raise _fail("combined result does not match its sealed stage reports")
    expected_sha256 = validator_bio_graph.digest(
        {key: value for key, value in root.items() if key != "report_sha256"}
    )
    if root["report_sha256"] != expected_sha256:
        raise _fail("combined validation report seal is invalid")
    return root


def validate_bio_feature_graph_paths(
    *,
    sequence_input: str | Path,
    gff3_source: str | Path,
    sequence_collection: str | Path,
    bio_artifact: str | Path,
    feature_graph_bundle: str | Path,
    source_metadata: str | Path | None = None,
    record_id: str | None = None,
) -> dict[str, Any]:
    """Snapshot all path inputs once, then validate the complete chain."""

    specifications = [
        ("fasta", sequence_input, "sequence source", validator_bio.MAX_SOURCE_INPUT_BYTES),
        ("gff3", gff3_source, "GFF3 source", validator_bio.MAX_GFF3_BYTES),
        (
            "sequence_collection",
            sequence_collection,
            "sequence collection artifact",
            validator_bio_graph.MAX_SOURCE_JSON_BYTES,
        ),
        ("bio", bio_artifact, "BioIR artifact", validator_bio_graph.MAX_ARTIFACT_BYTES),
        (
            "bundle",
            feature_graph_bundle,
            "feature-graph bundle",
            validator_bio_graph.MAX_ARTIFACT_BYTES,
        ),
    ]
    if source_metadata is not None:
        specifications.append(
            ("metadata", source_metadata, "source metadata", validator_bio.MAX_TEXT_BYTES)
        )
    snapshots = {
        name: _read_regular(path, label, maximum)
        for name, path, label, maximum in specifications
    }
    return validate_bio_feature_graph_chain(
        fasta_source=snapshots["fasta"],
        gff3_source=snapshots["gff3"],
        sequence_collection=_json_object(
            snapshots["sequence_collection"], "sequence collection artifact"
        ),
        bio_artifact=_json_object(snapshots["bio"], "BioIR artifact"),
        feature_graph_bundle=_json_object(
            snapshots["bundle"], "feature-graph bundle"
        ),
        source_metadata=snapshots.get("metadata"),
        record_id=record_id,
    )


def _render(report: dict[str, Any]) -> bytes:
    validate_chain_report(report)
    try:
        return (
            json.dumps(
                report,
                indent=2,
                sort_keys=True,
                ensure_ascii=False,
                allow_nan=False,
            )
            + "\n"
        ).encode("utf-8")
    except (TypeError, ValueError, UnicodeError) as failure:
        raise _fail(f"validation report is not I-JSON: {failure}") from failure


def save_report(report: dict[str, Any], path: str | Path) -> None:
    """Atomically save a validated deterministic report."""

    rendered = _render(report)
    destination = Path(path)
    descriptor = -1
    temporary: str | None = None
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        try:
            existing = destination.lstat()
        except FileNotFoundError:
            existing = None
        if existing is not None and (
            stat.S_ISLNK(existing.st_mode) or not stat.S_ISREG(existing.st_mode)
        ):
            raise _fail("validation report path must be a regular non-linked file")
        descriptor, temporary = tempfile.mkstemp(
            prefix=f".{destination.name}.", dir=destination.parent
        )
        with os.fdopen(descriptor, "wb") as stream:
            descriptor = -1
            stream.write(rendered)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, destination)
        temporary = None
    except BioChainValidationError:
        raise
    except OSError as failure:
        raise _fail(f"cannot write validation report: {failure}") from failure
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        if temporary is not None:
            try:
                os.unlink(temporary)
            except OSError:
                pass


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sequence_input", type=Path)
    parser.add_argument("gff3_source", type=Path)
    parser.add_argument("sequence_collection", type=Path)
    parser.add_argument("bio_artifact", type=Path)
    parser.add_argument("feature_graph_bundle", type=Path)
    parser.add_argument("--source-metadata", type=Path)
    parser.add_argument("--record-id")
    parser.add_argument("--report", type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the complete standalone validator from paths."""

    arguments = _parser().parse_args(argv)
    try:
        report = validate_bio_feature_graph_paths(
            sequence_input=arguments.sequence_input,
            gff3_source=arguments.gff3_source,
            sequence_collection=arguments.sequence_collection,
            bio_artifact=arguments.bio_artifact,
            feature_graph_bundle=arguments.feature_graph_bundle,
            source_metadata=arguments.source_metadata,
            record_id=arguments.record_id,
        )
        if arguments.report is not None:
            save_report(report, arguments.report)
        sys.stdout.buffer.write(_render(report))
        return 0
    except (
        BioChainValidationError,
        validator_bio.BioValidationError,
        validator_bio_graph.BioGraphValidationError,
    ) as failure:
        print(f"validation failed: {failure}", file=sys.stderr)
        return 1


__all__ = [
    "BioChainValidationError",
    "main",
    "save_report",
    "validate_bio_feature_graph_chain",
    "validate_bio_feature_graph_paths",
    "validate_chain_report",
]


if __name__ == "__main__":
    raise SystemExit(main())

"""Explicit lowering from normalized GFF3 BioIR to a development graph.

This backend is deliberately mechanical.  A BioIR feature becomes one unit and
an explicit ``Parent`` or ``Derives_from`` relationship becomes one edge.  The
numeric runtime state is accompanied by a content-bound semantics artifact that
preserves the complete BioIR and every numeric-code dictionary.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
import errno
import json
import os
from pathlib import Path
import secrets
import stat
import tempfile
from typing import Any

from ._canonical import ContractError, artifact_digest, canonical_bytes, digest
from .bio import (
    GFF3Artifact,
    GFF3Error,
    validate_gff3_artifact,
)
from .v2 import (
    compile_module,
    inline_storage,
    make_request,
    pack,
    policy_artifact,
    save,
    target_artifact,
)
from .v2._common import V2Error, pretty_bytes, seal
from .v2.compiler import COMPILER as MODULE_COMPILER
from .v2.limits import MAX_JSON_BYTES as MAX_GRAPH_JSON_BYTES
from .v2.provider import load_source
from .v2.target import DEV_DOMAIN, OP_EDGE_CREATE, OP_UNIT_CREATE


BACKEND_ID = "io.github.dttdrv.epigenesis.bio.feature-graph"
BACKEND_VERSION = 1
SEMANTICS_FORMAT = "brainc.bio.feature-graph-semantics"
RECORD_FORMAT = "brainc.bio.feature-graph-compilation"
BUNDLE_FORMAT = "brainc.bio.feature-graph-bundle"

RELATIONSHIP_KINDS = ("derives-from", "parent")
STRANDS = ("+", "-", ".", "?")
MAX_U64 = 2**64 - 1
_HAS_DIRECTORY_DESCRIPTOR = (
    os.name == "posix"
    and all(
        operation in getattr(os, "supports_dir_fd", set())
        for operation in (
            os.open,
            os.unlink,
            os.rename,
            os.stat,
            os.mkdir,
            os.rmdir,
        )
    )
)
_directory_fsync = os.fsync

FEATURE_FIELD_SPECS: tuple[dict[str, Any], ...] = (
    {
        "id": "attribute_count",
        "dtype": "u64",
        "tensor_id": "t.feature.attribute_count",
        "formula": {
            "op": "sum",
            "over": "feature.segments",
            "value": {"op": "length", "path": "segment.attributes"},
        },
    },
    {
        "id": "declared_id",
        "dtype": "bool",
        "tensor_id": "t.feature.declared_id",
        "formula": {"op": "is-not-null", "path": "feature.declared_id"},
    },
    {
        "id": "interval_count",
        "dtype": "u64",
        "tensor_id": "t.feature.interval_count",
        "formula": {
            "op": "sum",
            "over": "feature.segments",
            "value": {"op": "length", "path": "segment.intervals"},
        },
    },
    {
        "id": "location_max_end",
        "dtype": "u64",
        "tensor_id": "t.feature.location_max_end",
        "formula": {
            "op": "max",
            "over": "feature.segments[].intervals[]",
            "value": {"path": "interval.end"},
        },
    },
    {
        "id": "location_min_start",
        "dtype": "u64",
        "tensor_id": "t.feature.location_min_start",
        "formula": {
            "op": "min",
            "over": "feature.segments[].intervals[]",
            "value": {"path": "interval.start"},
        },
    },
    {
        "id": "phase_mask",
        "dtype": "u64",
        "tensor_id": "t.feature.phase_mask",
        "formula": {
            "op": "sum-distinct",
            "over": "feature.segments",
            "where": {"op": "is-not-null", "path": "segment.phase"},
            "value": {
                "op": "left-shift",
                "left": 1,
                "right": {"path": "segment.phase"},
            },
        },
    },
    {
        "id": "segment_bases",
        "dtype": "u64",
        "tensor_id": "t.feature.segment_bases",
        "formula": {
            "op": "sum",
            "over": "feature.segments[].intervals[]",
            "value": {
                "op": "subtract",
                "left": {"path": "interval.end"},
                "right": {"path": "interval.start"},
            },
        },
    },
    {
        "id": "segment_count",
        "dtype": "u64",
        "tensor_id": "t.feature.segment_count",
        "formula": {"op": "length", "path": "feature.segments"},
    },
    {
        "id": "seqid_code",
        "dtype": "u64",
        "tensor_id": "t.feature.seqid_code",
        "formula": {
            "op": "dictionary-code",
            "dictionary": "seqids",
            "key": {"path": "feature.seqid"},
        },
    },
    {
        "id": "source_code",
        "dtype": "u64",
        "tensor_id": "t.feature.source_code",
        "formula": {
            "op": "dictionary-code",
            "dictionary": "sources",
            "key": {"path": "feature.source"},
        },
    },
    {
        "id": "strand_code",
        "dtype": "u64",
        "tensor_id": "t.feature.strand_code",
        "formula": {
            "op": "dictionary-code",
            "dictionary": "strands",
            "key": {"path": "feature.strand"},
        },
    },
    {
        "id": "type_code",
        "dtype": "u64",
        "tensor_id": "t.feature.type_code",
        "formula": {
            "op": "dictionary-code",
            "dictionary": "types",
            "key": {"path": "feature.type"},
        },
    },
    {
        "id": "wraps_origin",
        "dtype": "bool",
        "tensor_id": "t.feature.wraps_origin",
        "formula": {
            "op": "any",
            "over": "feature.segments",
            "predicate": {
                "op": "greater-than",
                "left": {"op": "length", "path": "segment.intervals"},
                "right": 1,
            },
        },
    },
)
FEATURE_FIELDS: tuple[tuple[str, str], ...] = tuple(
    (field["id"], field["dtype"]) for field in FEATURE_FIELD_SPECS
)

TARGET_CONTRACT_SPEC = {
    "id": BACKEND_ID,
    "abi_major": 1,
    "opsets": [{"domain": DEV_DOMAIN, "version": 1}],
    "unit_schemas": [
        {
            "id": "feature",
            "fields": [
                {
                    "id": field["id"],
                    "type": {"dtype": field["dtype"], "shape": []},
                    "unit": None,
                    "mutability": "constant",
                    "numeric": {
                        "kind": "boolean" if field["dtype"] == "bool" else "integer"
                    },
                }
                for field in FEATURE_FIELD_SPECS
            ],
        }
    ],
    "edge_schemas": [
        {
            "id": f"relationship.{kind}",
            "fields": [
                {
                    "id": "kind_code",
                    "type": {"dtype": "u64", "shape": []},
                    "unit": None,
                    "mutability": "constant",
                    "numeric": {"kind": "integer"},
                }
            ],
        }
        for kind in RELATIONSHIP_KINDS
    ],
    "ports": [],
    "rules": [],
}

OPERATION_SPECS: tuple[dict[str, Any], ...] = (
    {
        "id": "create.features",
        "op": OP_UNIT_CREATE,
        "version": 1,
        "schema": "feature",
        "count": "t.feature.count",
        "initializers": [
            {"field": field["id"], "tensor": field["tensor_id"]}
            for field in FEATURE_FIELD_SPECS
        ],
    },
    *(
        {
            "id": f"create.relationship.{kind}",
            "op": OP_EDGE_CREATE,
            "version": 1,
            "schema": f"relationship.{kind}",
            "sources": "create.features",
            "targets": "create.features",
            "pairs": f"t.relationship.{kind}.pairs",
            "initializers": [
                {
                    "field": "kind_code",
                    "tensor": f"t.relationship.{kind}.kind_code",
                }
            ],
        }
        for kind in RELATIONSHIP_KINDS
    ),
)

BIO_GRAPH_PRODUCER = {
    "name": "brainc-bio-feature-graph",
    "version": "0.6.0",
    "passes": [
        "validate-bound-gff3-bio-ir",
        "construct-complete-normalized-bio-ir-sidecar",
        "encode-feature-tensors",
        "encode-explicit-relationship-tensors",
        "compile-development-module",
        "bind-compilation-record",
    ],
}

BACKEND_SPEC = {
    "format": "brainc.bio.feature-graph-backend-spec",
    "version": 1,
    "id": BACKEND_ID,
    "backend_version": BACKEND_VERSION,
    "input": {
        "sequence": {"format": "brain01.sequence-collection-ir", "version": 1},
        "bio_ir": {
            "format": "brainc.bio.gff3-ir",
            "version": 1,
            "coordinate_system": "0-based-half-open",
        },
        "binding": "exact-artifact-and-ir-sha256",
    },
    "semantics_sidecar": {
        "bio_ir": "exact-deep-copy-of-normalized-bio-ir",
        "unit_order": "bio-ir-feature-array-order",
        "edge_order": {
            "kind_order": list(RELATIONSHIP_KINDS),
            "within_kind": "bio-ir-relationship-array-order",
            "entry_fields": [
                "edge_index",
                "bio_ir_relationship_index",
                "source_unit_index",
                "target_unit_index",
                "source",
                "kind",
                "target",
            ],
            "endpoints": "zero-based-indices-into-unit-order",
            "indices": {
                "edge_index": "zero-based-position-in-emitted-edge-order",
                "bio_ir_relationship_index": (
                    "zero-based-position-in-bio-ir-relationships"
                ),
            },
            "copied_fields": {
                "source": "relationship.source",
                "kind": "relationship.kind",
                "target": "relationship.target",
            },
        },
        "dictionary_wire": {
            "container": "object-keyed-by-dictionary-id",
            "value": "array",
            "entry": {
                "code": "zero-based-array-position",
                "value": "derived-dictionary-value",
            },
            "entry_order": "ascending-code",
        },
    },
    "dictionaries": [
        {
            "id": "relationship_kinds",
            "values": {"op": "literal", "items": list(RELATIONSHIP_KINDS)},
            "codes": "zero-based-array-position",
        },
        {
            "id": "seqids",
            "values": {
                "op": "distinct-sort",
                "path": "bio_ir.features[].seqid",
                "order": "ascending-utf8-bytes",
            },
            "codes": "zero-based-array-position",
        },
        {
            "id": "sources",
            "values": {
                "op": "prefix-distinct-sort",
                "prefix": [None],
                "path": "bio_ir.features[].source",
                "exclude": [None],
                "order": "ascending-utf8-bytes",
            },
            "codes": "zero-based-array-position",
        },
        {
            "id": "strands",
            "values": {"op": "literal", "items": list(STRANDS)},
            "codes": "zero-based-array-position",
        },
        {
            "id": "types",
            "values": {
                "op": "distinct-sort",
                "path": "bio_ir.features[].type",
                "order": "ascending-utf8-bytes",
            },
            "codes": "zero-based-array-position",
        },
    ],
    "tensors": {
        "encoding": {
            "byte_order": "little-endian",
            "bool_bytes": {"false": 0, "true": 1},
            "storage": "inline-base64-with-byte-length-and-sha256",
            "unit": None,
            "axes": {
                "scalar": [],
                "ranked": "one-null-entry-per-tensor-dimension",
            },
            "output_order": "ascending-tensor-id",
        },
        "feature_count": {
            "id": "t.feature.count",
            "dtype": "u64",
            "shape": [],
            "formula": {"op": "length", "path": "bio_ir.features"},
        },
        "feature_fields": [
            {**copy.deepcopy(field), "shape": ["feature-count"]}
            for field in FEATURE_FIELD_SPECS
        ],
        "relationships": [
            {
                "kind": kind,
                "kind_code": code,
                "selection": {
                    "op": "filter",
                    "path": "semantics.edge_order",
                    "where": {"field": "kind", "equals": kind},
                },
                "relationship_count": {
                    "op": "length",
                    "of": "selection",
                },
                "kind_tensor": {
                    "id": f"t.relationship.{kind}.kind_code",
                    "dtype": "u64",
                    "shape": [],
                    "value": code,
                },
                "pairs_tensor": {
                    "id": f"t.relationship.{kind}.pairs",
                    "dtype": "u64",
                    "shape": ["relationship-count", 2],
                    "row": ["source_unit_index", "target_unit_index"],
                },
            }
            for code, kind in enumerate(RELATIONSHIP_KINDS)
        ],
    },
    "target_contract": copy.deepcopy(TARGET_CONTRACT_SPEC),
    "operations": copy.deepcopy(list(OPERATION_SPECS)),
    "development_abi": {
        "domain": DEV_DOMAIN,
        "version": 1,
        "initializer_semantics": {
            "applies_to": [OP_UNIT_CREATE, OP_EDGE_CREATE],
            "field_order": "exact-target-schema-field-order",
            "field_coverage": "every-field-exactly-once",
            "created_count": {
                OP_UNIT_CREATE: "value-of-unitless-u64-count-scalar",
                OP_EDGE_CREATE: "row-count-of-unitless-u64-pairs-tensor",
            },
            "allowed_shapes": [
                {
                    "shape": "field-shape",
                    "meaning": "broadcast-one-value-to-every-created-item",
                },
                {
                    "shape": "[created-count]+field-shape",
                    "meaning": "one-value-for-each-created-item",
                },
            ],
        },
    },
    "artifact_profile": {
        "backend_spec": {
            "format": "brainc.bio.feature-graph-backend-spec",
            "version": 1,
        },
        "backend_semantics": {"format": SEMANTICS_FORMAT, "version": 1},
        "provider_manifest": {
            "format": "brainc.provider-manifest",
            "version": 2,
            "provider": {"name": BACKEND_ID, "version": str(BACKEND_VERSION)},
            "model_identity": {
                "kind": "content-sha256",
                "value": "backend-spec-artifact-sha256",
            },
            "accepts": ["brain01.sequence-collection-ir/v1"],
        },
        "prediction_request": {
            "format": "brainc.prediction-request",
            "version": 2,
        },
        "prediction_response": {
            "format": "brainc.prediction-response",
            "version": 2,
        },
        "target_contract": {"format": "brainc.target-contract", "version": 1},
        "lowering_policy": {
            "format": "brainc.lowering-policy",
            "version": 2,
            "id": f"{BACKEND_ID}.lowering",
        },
        "development_module": {
            "format": "brainc.development-module",
            "version": 1,
            "producer": copy.deepcopy(MODULE_COMPILER),
        },
        "compilation_record": {
            "format": RECORD_FORMAT,
            "version": 1,
            "producer": copy.deepcopy(BIO_GRAPH_PRODUCER),
        },
        "bundle": {"format": BUNDLE_FORMAT, "version": 1},
    },
    "capabilities": {
        "ontology_interpretation": False,
        "inference": False,
        "learning": False,
        "rules": False,
        "ports": False,
    },
}
_BACKEND_SPEC_CANONICAL = canonical_bytes(BACKEND_SPEC)
BACKEND_SPEC_SHA256 = digest(BACKEND_SPEC)
PRODUCER = copy.deepcopy(BIO_GRAPH_PRODUCER)


def backend_spec_artifact() -> dict[str, Any]:
    """Return a fresh sealed copy of the immutable normative backend core."""

    core = json.loads(_BACKEND_SPEC_CANONICAL.decode("utf-8"))
    return seal(core)


def _require_backend_implementation_integrity(
    backend_spec: dict[str, Any],
) -> None:
    """Reject execution if public implementation constants drift from the ABI."""

    try:
        stored = artifact_digest(backend_spec)
    except ContractError as failure:
        raise _fail(f"normative backend specification is invalid: {failure}") from failure
    if stored != BACKEND_SPEC_SHA256:
        raise _fail("normative backend specification digest changed")
    definitions = {
        item["id"]: item for item in backend_spec["dictionaries"]
    }
    feature_fields = tuple(
        {key: copy.deepcopy(value) for key, value in item.items() if key != "shape"}
        for item in backend_spec["tensors"]["feature_fields"]
    )
    profile = backend_spec["artifact_profile"]
    expected_constants = (
        backend_spec["id"] == BACKEND_ID,
        backend_spec["backend_version"] == BACKEND_VERSION,
        tuple(definitions["relationship_kinds"]["values"]["items"])
        == RELATIONSHIP_KINDS,
        tuple(definitions["strands"]["values"]["items"]) == STRANDS,
        feature_fields == FEATURE_FIELD_SPECS,
        backend_spec["target_contract"] == TARGET_CONTRACT_SPEC,
        tuple(backend_spec["operations"]) == OPERATION_SPECS,
        profile["backend_semantics"]
        == {"format": SEMANTICS_FORMAT, "version": 1},
        profile["compilation_record"]["format"] == RECORD_FORMAT,
        profile["bundle"] == {"format": BUNDLE_FORMAT, "version": 1},
        profile["development_module"]["producer"] == MODULE_COMPILER,
        profile["compilation_record"]["producer"] == PRODUCER,
    )
    if not all(expected_constants):
        raise _fail("backend implementation constants drifted from the normative spec")

class FeatureGraphError(ContractError):
    """A BioIR feature-graph compilation violates the closed backend contract."""


def _fail(detail: str) -> FeatureGraphError:
    return FeatureGraphError(f"BIOGRAPH001: {detail}")


def _preflight_wire(value: dict[str, Any], label: str) -> None:
    """Require an artifact to fit the exact JSON writer used by this backend."""

    try:
        pretty_bytes(value)
    except V2Error as failure:
        raise _fail(
            f"{label} is outside the {MAX_GRAPH_JSON_BYTES}-byte graph wire contract: "
            f"{failure}"
        ) from failure


def _save_staged_artifact(
    artifact: dict[str, Any],
    role: str,
    staging: Path,
    directory_descriptor: int,
) -> None:
    """Write one graph artifact without resolving a substituted POSIX path."""

    if directory_descriptor < 0:
        save(artifact, staging / f"{role}.json")
        return
    raw = pretty_bytes(artifact)
    final_name = f"{role}.json"
    temporary_name = f".{final_name}.{secrets.token_hex(12)}"
    descriptor = -1
    temporary_exists = False
    try:
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
        flags |= getattr(os, "O_CLOEXEC", 0)
        flags |= getattr(os, "O_NOFOLLOW", 0)
        descriptor = os.open(
            temporary_name,
            flags,
            0o600,
            dir_fd=directory_descriptor,
        )
        temporary_exists = True
        view = memoryview(raw)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise OSError("artifact write made no progress")
            view = view[written:]
        os.fsync(descriptor)
        os.close(descriptor)
        descriptor = -1
        try:
            os.stat(
                final_name,
                dir_fd=directory_descriptor,
                follow_symlinks=False,
            )
        except FileNotFoundError:
            pass
        else:
            raise OSError(f"staged artifact already exists: {final_name}")
        os.rename(
            temporary_name,
            final_name,
            src_dir_fd=directory_descriptor,
            dst_dir_fd=directory_descriptor,
        )
        temporary_exists = False
    except V2Error:
        raise
    except OSError as failure:
        raise V2Error(f"cannot write staged {role}: {failure}") from failure
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        if temporary_exists:
            try:
                os.unlink(temporary_name, dir_fd=directory_descriptor)
            except OSError:
                pass


def _bio_payload(value: GFF3Artifact | dict[str, Any]) -> dict[str, Any]:
    if isinstance(value, GFF3Artifact):
        return value.to_dict()
    if type(value) is dict:
        return copy.deepcopy(value)
    raise _fail("BioIR input must be a GFF3Artifact or artifact object")


def _source_binding(source: dict[str, Any]) -> dict[str, Any]:
    if (source.get("format"), source.get("version")) != (
        "brain01.sequence-collection-ir",
        1,
    ):
        raise _fail("the feature-graph backend requires Sequence Collection IR v1")
    return {
        "format": source["format"],
        "version": source["version"],
        "artifact_sha256": source["artifact_sha256"],
        "ir_sha256": source["collection_ir_sha256"],
    }


def _bio_binding(bio: dict[str, Any]) -> dict[str, Any]:
    return {
        "format": bio["format"],
        "version": bio["version"],
        "artifact_sha256": bio["artifact_sha256"],
        "bio_ir_sha256": bio["bio_ir_sha256"],
    }


def _indexed(values: list[Any]) -> list[dict[str, Any]]:
    return [{"code": index, "value": value} for index, value in enumerate(values)]


def _utf8_sorted(values: set[str]) -> list[str]:
    return sorted(values, key=lambda value: value.encode("utf-8"))


def _semantics_artifact(
    source: dict[str, Any], bio: dict[str, Any]
) -> dict[str, Any]:
    bio_ir = bio["bio_ir"]
    features = bio_ir["features"]
    seqids = _utf8_sorted({feature["seqid"] for feature in features})
    sources = [None] + _utf8_sorted(
        {feature["source"] for feature in features if feature["source"] is not None}
    )
    types = _utf8_sorted({feature["type"] for feature in features})
    feature_index = {
        feature["entity_id"]: index for index, feature in enumerate(features)
    }
    edge_order: list[dict[str, Any]] = []
    next_edge = 0
    for kind in RELATIONSHIP_KINDS:
        for relationship_index, relationship in enumerate(bio_ir["relationships"]):
            if relationship["kind"] != kind:
                continue
            edge_order.append(
                {
                    "edge_index": next_edge,
                    "bio_ir_relationship_index": relationship_index,
                    "source_unit_index": feature_index[relationship["source"]],
                    "target_unit_index": feature_index[relationship["target"]],
                    "source": relationship["source"],
                    "kind": relationship["kind"],
                    "target": relationship["target"],
                }
            )
            next_edge += 1
    semantics = {
        "dictionaries": {
            "relationship_kinds": _indexed(list(RELATIONSHIP_KINDS)),
            "seqids": _indexed(seqids),
            "sources": _indexed(sources),
            "strands": _indexed(list(STRANDS)),
            "types": _indexed(types),
        },
        "unit_order": [feature["entity_id"] for feature in features],
        "edge_order": edge_order,
        "bio_ir": copy.deepcopy(bio_ir),
    }
    return seal(
        {
            "format": SEMANTICS_FORMAT,
            "version": 1,
            "backend": {
                "id": BACKEND_ID,
                "version": BACKEND_VERSION,
                "spec_sha256": BACKEND_SPEC_SHA256,
            },
            "inputs": {
                "sequence": _source_binding(source),
                "bio_ir": _bio_binding(bio),
            },
            "semantics": semantics,
            "semantics_sha256": digest(semantics),
        }
    )


def _dictionary_codes(
    semantics: dict[str, Any], name: str
) -> dict[Any, int]:
    return {
        entry["value"]: entry["code"]
        for entry in semantics["semantics"]["dictionaries"][name]
    }


def _checked_u64(value: int, label: str) -> int:
    if type(value) is not int or value < 0 or value > MAX_U64:
        raise _fail(f"{label} is outside the u64 range")
    return value


def _feature_values(
    semantics: dict[str, Any],
) -> dict[str, list[bool | int]]:
    features = semantics["semantics"]["bio_ir"]["features"]
    seqid_codes = _dictionary_codes(semantics, "seqids")
    source_codes = _dictionary_codes(semantics, "sources")
    strand_codes = _dictionary_codes(semantics, "strands")
    type_codes = _dictionary_codes(semantics, "types")
    values: dict[str, list[bool | int]] = {field: [] for field, _ in FEATURE_FIELDS}
    for feature_index, feature in enumerate(features):
        intervals = [
            interval
            for segment in feature["segments"]
            for interval in segment["intervals"]
        ]
        segment_bases = sum(
            interval["end"] - interval["start"] for interval in intervals
        )
        phases = {
            segment["phase"]
            for segment in feature["segments"]
            if segment["phase"] is not None
        }
        phase_mask = sum(1 << phase for phase in phases)
        attribute_count = sum(
            len(segment["attributes"]) for segment in feature["segments"]
        )
        numeric = {
            "attribute_count": attribute_count,
            "interval_count": len(intervals),
            "location_max_end": max(interval["end"] for interval in intervals),
            "location_min_start": min(interval["start"] for interval in intervals),
            "phase_mask": phase_mask,
            "segment_bases": segment_bases,
            "segment_count": len(feature["segments"]),
            "seqid_code": seqid_codes[feature["seqid"]],
            "source_code": source_codes[feature["source"]],
            "strand_code": strand_codes[feature["strand"]],
            "type_code": type_codes[feature["type"]],
        }
        for field, value in numeric.items():
            values[field].append(
                _checked_u64(value, f"feature {feature_index} field {field}")
            )
        values["declared_id"].append(feature["declared_id"] is not None)
        values["wraps_origin"].append(
            any(len(segment["intervals"]) > 1 for segment in feature["segments"])
        )
    return values


def _tensor(
    identifier: str,
    dtype: str,
    shape: list[int],
    values: list[bool | int],
) -> dict[str, Any]:
    return {
        "id": identifier,
        "type": {"dtype": dtype, "shape": shape},
        "unit": None,
        "axes": [None] * len(shape),
        "storage": inline_storage(pack(dtype, values)),
    }


def _outputs(semantics: dict[str, Any]) -> list[dict[str, Any]]:
    features = semantics["semantics"]["bio_ir"]["features"]
    feature_values = _feature_values(semantics)
    result = [_tensor("t.feature.count", "u64", [], [len(features)])]
    for field, dtype in FEATURE_FIELDS:
        result.append(
            _tensor(
                f"t.feature.{field}",
                dtype,
                [len(features)],
                feature_values[field],
            )
        )
    kind_codes = _dictionary_codes(semantics, "relationship_kinds")
    edge_order = semantics["semantics"]["edge_order"]
    for kind in RELATIONSHIP_KINDS:
        edges = [edge for edge in edge_order if edge["kind"] == kind]
        pairs = [
            endpoint
            for edge in edges
            for endpoint in (edge["source_unit_index"], edge["target_unit_index"])
        ]
        stem = f"t.relationship.{kind}"
        result.append(_tensor(f"{stem}.kind_code", "u64", [], [kind_codes[kind]]))
        result.append(_tensor(f"{stem}.pairs", "u64", [len(edges), 2], pairs))
    return sorted(result, key=lambda value: value["id"])


def _target() -> dict[str, Any]:
    return target_artifact(copy.deepcopy(TARGET_CONTRACT_SPEC))


def _operations() -> list[dict[str, Any]]:
    return copy.deepcopy(list(OPERATION_SPECS))


def _manifest(outputs: list[dict[str, Any]]) -> dict[str, Any]:
    return seal(
        {
            "format": "brainc.provider-manifest",
            "version": 2,
            "provider": {"name": BACKEND_ID, "version": str(BACKEND_VERSION)},
            "model_identity": {
                "kind": "content-sha256",
                "value": BACKEND_SPEC_SHA256,
            },
            "accepts": ["brain01.sequence-collection-ir/v1"],
            "outputs": [
                {
                    key: copy.deepcopy(output[key])
                    for key in ("id", "type", "unit", "axes")
                }
                for output in outputs
            ],
        }
    )


def _response(
    request: dict[str, Any], manifest: dict[str, Any], outputs: list[dict[str, Any]]
) -> dict[str, Any]:
    return seal(
        {
            "format": "brainc.prediction-response",
            "version": 2,
            "request_artifact_sha256": request["artifact_sha256"],
            "provider": copy.deepcopy(manifest["provider"]),
            "model_identity": copy.deepcopy(manifest["model_identity"]),
            "outputs": copy.deepcopy(outputs),
        }
    )


def _record(
    source: dict[str, Any],
    bio: dict[str, Any],
    semantics: dict[str, Any],
    manifest: dict[str, Any],
    request: dict[str, Any],
    response: dict[str, Any],
    target: dict[str, Any],
    policy: dict[str, Any],
    module: dict[str, Any],
    backend_spec: dict[str, Any],
) -> dict[str, Any]:
    return seal(
        {
            "format": RECORD_FORMAT,
            "version": 1,
            "producer": copy.deepcopy(PRODUCER),
            "backend": {
                "id": BACKEND_ID,
                "version": BACKEND_VERSION,
                "spec_sha256": BACKEND_SPEC_SHA256,
                "semantics_sha256": semantics["semantics_sha256"],
            },
            "inputs": {
                "sequence": _source_binding(source),
                "bio_ir": _bio_binding(bio),
            },
            "artifacts": {
                "backend_spec": {
                    "artifact_sha256": backend_spec["artifact_sha256"]
                },
                "backend_semantics": {
                    "artifact_sha256": semantics["artifact_sha256"]
                },
                "development_module": {
                    "artifact_sha256": module["artifact_sha256"]
                },
                "lowering_policy": {
                    "artifact_sha256": policy["artifact_sha256"]
                },
                "prediction_request": {
                    "artifact_sha256": request["artifact_sha256"]
                },
                "prediction_response": {
                    "artifact_sha256": response["artifact_sha256"]
                },
                "provider_manifest": {
                    "artifact_sha256": manifest["artifact_sha256"]
                },
                "target_contract": {
                    "artifact_sha256": target["artifact_sha256"]
                },
            },
            "result": {
                "module_sha256": module["module_sha256"],
                "operations": module["module"]["budgets"]["operations"],
                "tensor_bytes": module["module"]["budgets"]["tensor_bytes"],
                "units": module["module"]["budgets"]["units"],
                "edges": module["module"]["budgets"]["edges"],
                "attachments": 0,
                "ports": 0,
                "inference": False,
                "learning": False,
            },
        }
    )


@dataclass(frozen=True, init=False)
class FeatureGraphBundle:
    """All independently hash-bound artifacts for one feature-graph lowering."""

    _artifact_wires: tuple[tuple[str, bytes], ...] = field(repr=False)

    def __init__(self) -> None:
        raise TypeError(
            "FeatureGraphBundle instances are created by compile_feature_graph"
        )

    @classmethod
    def _from_artifacts(
        cls,
        backend_spec: dict[str, Any],
        semantics: dict[str, Any],
        provider_manifest: dict[str, Any],
        prediction_request: dict[str, Any],
        prediction_response: dict[str, Any],
        target_contract: dict[str, Any],
        lowering_policy: dict[str, Any],
        development_module: dict[str, Any],
        compilation_record: dict[str, Any],
    ) -> FeatureGraphBundle:
        artifacts = (
            ("backend_spec", backend_spec),
            ("backend_semantics", semantics),
            ("compilation_record", compilation_record),
            ("development_module", development_module),
            ("lowering_policy", lowering_policy),
            ("prediction_request", prediction_request),
            ("prediction_response", prediction_response),
            ("provider_manifest", provider_manifest),
            ("target_contract", target_contract),
        )
        try:
            wires = tuple(
                (role, canonical_bytes(copy.deepcopy(artifact)))
                for role, artifact in artifacts
            )
        except (ContractError, TypeError, ValueError) as failure:
            raise _fail(f"cannot freeze feature-graph bundle: {failure}") from failure
        bundle = object.__new__(cls)
        object.__setattr__(bundle, "_artifact_wires", wires)
        bundle._validated_artifacts()
        return bundle

    def _artifact(self, role: str) -> dict[str, Any]:
        for candidate, wire in self._artifact_wires:
            if candidate == role:
                value = json.loads(wire)
                if type(value) is not dict:
                    raise _fail(f"invalid {role} child artifact")
                return value
        raise _fail(f"missing {role} child artifact")

    @property
    def backend_spec(self) -> dict[str, Any]:
        return self._artifact("backend_spec")

    @property
    def semantics(self) -> dict[str, Any]:
        return self._artifact("backend_semantics")

    @property
    def provider_manifest(self) -> dict[str, Any]:
        return self._artifact("provider_manifest")

    @property
    def prediction_request(self) -> dict[str, Any]:
        return self._artifact("prediction_request")

    @property
    def prediction_response(self) -> dict[str, Any]:
        return self._artifact("prediction_response")

    @property
    def target_contract(self) -> dict[str, Any]:
        return self._artifact("target_contract")

    @property
    def lowering_policy(self) -> dict[str, Any]:
        return self._artifact("lowering_policy")

    @property
    def development_module(self) -> dict[str, Any]:
        return self._artifact("development_module")

    @property
    def compilation_record(self) -> dict[str, Any]:
        return self._artifact("compilation_record")

    def _validated_artifacts(self) -> dict[str, dict[str, Any]]:
        artifacts = {
            role: self._artifact(role) for role, _ in self._artifact_wires
        }
        for role, artifact in artifacts.items():
            try:
                artifact_digest(artifact)
            except (ContractError, KeyError, TypeError) as failure:
                raise _fail(f"invalid {role} child artifact: {failure}") from failure
        if artifacts["backend_spec"] != backend_spec_artifact():
            raise _fail("backend specification child is not the normative artifact")
        record = artifacts["compilation_record"]
        bindings = record.get("artifacts")
        bound_roles = set(artifacts) - {"compilation_record"}
        if type(bindings) is not dict or set(bindings) != bound_roles:
            raise _fail("compilation record child bindings are incomplete")
        for role in sorted(bound_roles):
            if bindings[role] != {
                "artifact_sha256": artifacts[role]["artifact_sha256"]
            }:
                raise _fail(f"compilation record does not bind child {role}")
        return artifacts

    def to_dict(self) -> dict[str, Any]:
        core = {
            "format": BUNDLE_FORMAT,
            "version": 1,
            "artifacts": self._validated_artifacts(),
        }
        payload = seal(core)
        _preflight_wire(payload, "feature-graph bundle")
        return payload

    def save(self, directory: str | Path) -> dict[str, Path]:
        destination = Path(directory)
        bundle_payload = self.to_dict()
        destination.parent.mkdir(parents=True, exist_ok=True)
        try:
            existing = destination.lstat()
        except FileNotFoundError:
            existing = None
        except OSError as failure:
            raise _fail(f"cannot inspect output directory: {failure}") from failure
        if existing is not None:
            kind = "symbolic link" if stat.S_ISLNK(existing.st_mode) else "existing path"
            raise _fail(
                f"feature-graph output directory must be absent, not a {kind}"
            )
        artifacts = {
            **bundle_payload["artifacts"],
            "bundle": bundle_payload,
        }
        parent_descriptor = -1
        parent_identity: tuple[int, int] | None = None
        if _HAS_DIRECTORY_DESCRIPTOR:
            try:
                parent_inspected = destination.parent.lstat()
                if stat.S_ISLNK(parent_inspected.st_mode) or not stat.S_ISDIR(
                    parent_inspected.st_mode
                ):
                    raise _fail(
                        "feature-graph output parent must be a non-linked directory"
                    )
                parent_flags = os.O_RDONLY
                parent_flags |= getattr(os, "O_CLOEXEC", 0)
                parent_flags |= getattr(os, "O_DIRECTORY", 0)
                parent_flags |= getattr(os, "O_NOFOLLOW", 0)
                parent_descriptor = os.open(destination.parent, parent_flags)
                parent_opened = os.fstat(parent_descriptor)
                parent_identity = (parent_opened.st_dev, parent_opened.st_ino)
                if (
                    not stat.S_ISDIR(parent_opened.st_mode)
                    or parent_identity
                    != (parent_inspected.st_dev, parent_inspected.st_ino)
                ):
                    raise _fail(
                        "feature-graph output parent changed while opening"
                    )
            except (FeatureGraphError, OSError) as failure:
                if parent_descriptor >= 0:
                    os.close(parent_descriptor)
                if isinstance(failure, FeatureGraphError):
                    raise
                raise _fail(
                    f"cannot open feature-graph output parent: {failure}"
                ) from failure
        try:
            staging = Path(
                tempfile.mkdtemp(
                    prefix=f".{destination.name}.", dir=destination.parent
                )
            )
        except OSError as failure:
            if parent_descriptor >= 0:
                os.close(parent_descriptor)
            raise _fail(
                f"cannot create feature-graph staging directory: {failure}"
            ) from failure
        staging_metadata = staging.lstat()
        staging_identity = (staging_metadata.st_dev, staging_metadata.st_ino)
        staging_descriptor = -1
        published = False
        staging_moved = False
        reservation_created = False
        reservation_identity: tuple[int, int] | None = None

        if _HAS_DIRECTORY_DESCRIPTOR:
            directory_flags = os.O_RDONLY
            directory_flags |= getattr(os, "O_CLOEXEC", 0)
            directory_flags |= getattr(os, "O_DIRECTORY", 0)
            directory_flags |= getattr(os, "O_NOFOLLOW", 0)
            try:
                staging_descriptor = os.open(staging, directory_flags)
                opened = os.fstat(staging_descriptor)
                if (
                    not stat.S_ISDIR(opened.st_mode)
                    or (opened.st_dev, opened.st_ino) != staging_identity
                ):
                    raise _fail(
                        "feature-graph staging directory changed while opening"
                    )
                parent_entry = os.stat(
                    staging.name,
                    dir_fd=parent_descriptor,
                    follow_symlinks=False,
                )
                if (
                    not stat.S_ISDIR(parent_entry.st_mode)
                    or (parent_entry.st_dev, parent_entry.st_ino)
                    != staging_identity
                ):
                    raise _fail(
                        "feature-graph staging directory is outside output parent"
                    )
            except (FeatureGraphError, OSError) as failure:
                if staging_descriptor >= 0:
                    os.close(staging_descriptor)
                    staging_descriptor = -1
                if parent_descriptor >= 0:
                    os.close(parent_descriptor)
                    parent_descriptor = -1
                try:
                    current = staging.lstat()
                    if (
                        stat.S_ISDIR(current.st_mode)
                        and not stat.S_ISLNK(current.st_mode)
                        and (current.st_dev, current.st_ino) == staging_identity
                    ):
                        staging.rmdir()
                except OSError:
                    pass
                if isinstance(failure, FeatureGraphError):
                    raise
                raise _fail(
                    f"cannot open feature-graph staging directory: {failure}"
                ) from failure

        def require_staging_identity() -> None:
            current = staging.lstat()
            if (
                stat.S_ISLNK(current.st_mode)
                or not stat.S_ISDIR(current.st_mode)
                or (current.st_dev, current.st_ino) != staging_identity
            ):
                raise _fail("feature-graph staging directory changed identity")

        try:
            for role, artifact in artifacts.items():
                require_staging_identity()
                try:
                    _save_staged_artifact(
                        artifact, role, staging, staging_descriptor
                    )
                except V2Error as failure:
                    raise _fail(f"cannot save {role}: {failure}") from failure
                require_staging_identity()
            if staging_descriptor >= 0:
                try:
                    _directory_fsync(staging_descriptor)
                except OSError as failure:
                    unsupported = {
                        errno.EBADF,
                        errno.EINVAL,
                        getattr(errno, "ENOTSUP", errno.EINVAL),
                        getattr(errno, "EOPNOTSUPP", errno.EINVAL),
                    }
                    if failure.errno not in unsupported:
                        raise
            try:
                if parent_descriptor >= 0:
                    os.mkdir(
                        destination.name,
                        0o700,
                        dir_fd=parent_descriptor,
                    )
                    reserved = os.stat(
                        destination.name,
                        dir_fd=parent_descriptor,
                        follow_symlinks=False,
                    )
                else:
                    # Windows directory renames are no-clobber: renaming a
                    # directory to an existing name fails.  It also cannot
                    # replace an empty directory, so do not create the POSIX
                    # reservation used by the descriptor-capable path.
                    reserved = None
            except FileExistsError as failure:
                raise _fail(
                    "feature-graph output path appeared during staging"
                ) from failure
            if parent_descriptor >= 0:
                assert reserved is not None
                reservation_created = True
                reservation_identity = (reserved.st_dev, reserved.st_ino)
                if not stat.S_ISDIR(reserved.st_mode) or stat.S_ISLNK(
                    reserved.st_mode
                ):
                    raise _fail(
                        "feature-graph output reservation is not a directory"
                    )
                os.replace(
                    staging.name,
                    destination.name,
                    src_dir_fd=parent_descriptor,
                    dst_dir_fd=parent_descriptor,
                )
                published_entry = os.stat(
                    destination.name,
                    dir_fd=parent_descriptor,
                    follow_symlinks=False,
                )
            else:
                try:
                    os.rename(staging, destination)
                except FileExistsError as failure:
                    raise _fail(
                        "feature-graph output path appeared during staging"
                    ) from failure
                published_entry = destination.lstat()
            staging_moved = True
            reservation_created = False
            if (
                not stat.S_ISDIR(published_entry.st_mode)
                or stat.S_ISLNK(published_entry.st_mode)
                or (published_entry.st_dev, published_entry.st_ino)
                != staging_identity
            ):
                raise _fail(
                    "feature-graph staging identity changed during publication"
                )
            if parent_descriptor >= 0:
                parent_finished = destination.parent.lstat()
                if (
                    stat.S_ISLNK(parent_finished.st_mode)
                    or not stat.S_ISDIR(parent_finished.st_mode)
                    or (parent_finished.st_dev, parent_finished.st_ino)
                    != parent_identity
                ):
                    raise _fail(
                        "feature-graph output parent changed during publication"
                    )
                try:
                    _directory_fsync(parent_descriptor)
                except OSError as failure:
                    unsupported = {
                        errno.EBADF,
                        errno.EINVAL,
                        getattr(errno, "ENOTSUP", errno.EINVAL),
                        getattr(errno, "EOPNOTSUPP", errno.EINVAL),
                    }
                    if failure.errno not in unsupported:
                        raise
            published = True
        except FeatureGraphError:
            raise
        except OSError as failure:
            raise _fail(f"cannot publish feature-graph directory: {failure}") from failure
        finally:
            if not published and staging_descriptor >= 0:
                # Address the original staging directory through its still-open
                # descriptor.  A hostile pathname substitution cannot redirect
                # cleanup into another directory.
                for role in artifacts:
                    try:
                        os.unlink(f"{role}.json", dir_fd=staging_descriptor)
                    except OSError:
                        pass
            elif not published:
                # Portable fallback for hosts without descriptor-relative
                # directory operations.  The staging directory is private;
                # recheck its identity around every exact known-name unlink.
                for role in artifacts:
                    try:
                        require_staging_identity()
                        child = staging / f"{role}.json"
                        child_metadata = child.lstat()
                        if not (
                            stat.S_ISREG(child_metadata.st_mode)
                            or stat.S_ISLNK(child_metadata.st_mode)
                        ):
                            break
                        child.unlink()
                        require_staging_identity()
                    except FileNotFoundError:
                        continue
                    except (FeatureGraphError, OSError):
                        break
            if staging_descriptor >= 0:
                os.close(staging_descriptor)
            if not published and not staging_moved:
                try:
                    current = staging.lstat()
                    if (
                        stat.S_ISDIR(current.st_mode)
                        and not stat.S_ISLNK(current.st_mode)
                        and (current.st_dev, current.st_ino) == staging_identity
                    ):
                        staging.rmdir()
                except OSError:
                    pass
            if reservation_created and reservation_identity is not None:
                try:
                    if parent_descriptor >= 0:
                        current = os.stat(
                            destination.name,
                            dir_fd=parent_descriptor,
                            follow_symlinks=False,
                        )
                        if (
                            stat.S_ISDIR(current.st_mode)
                            and (current.st_dev, current.st_ino)
                            == reservation_identity
                        ):
                            os.rmdir(
                                destination.name,
                                dir_fd=parent_descriptor,
                            )
                    else:
                        current = destination.lstat()
                        if (
                            stat.S_ISDIR(current.st_mode)
                            and not stat.S_ISLNK(current.st_mode)
                            and (current.st_dev, current.st_ino)
                            == reservation_identity
                        ):
                            destination.rmdir()
                except OSError:
                    pass
            if parent_descriptor >= 0:
                os.close(parent_descriptor)
        return {
            role: destination / f"{role}.json" for role in artifacts
        }


def compile_feature_graph(
    source_path: str | Path,
    bio_artifact: GFF3Artifact | dict[str, Any],
    *,
    gff3_source: bytes | None = None,
) -> FeatureGraphBundle:
    """Compile validated GFF3 BioIR to an explicit structural graph module."""

    try:
        backend_spec = backend_spec_artifact()
        _require_backend_implementation_integrity(backend_spec)
        source, _ = load_source(source_path)
        _source_binding(source)
        bio = _bio_payload(bio_artifact)
        validate_gff3_artifact(bio, source, gff3_source=gff3_source)
        _preflight_wire(bio, "BioIR input")
        semantics = _semantics_artifact(source, bio)
        _preflight_wire(semantics, "backend semantics")
        outputs = _outputs(semantics)
        manifest = _manifest(outputs)
        target = _target()
        target_binding = {
            "id": target["contract"]["id"],
            "abi_major": target["contract"]["abi_major"],
            "contract_sha256": target["contract_sha256"],
        }
        with tempfile.TemporaryDirectory(prefix="brainc-bio-graph-") as temporary:
            root = Path(temporary)
            source_snapshot = root / "source.json"
            manifest_path = root / "manifest.json"
            target_path = root / "target.json"
            save(source, source_snapshot)
            save(manifest, manifest_path)
            save(target, target_path)
            request = make_request(
                source_snapshot,
                manifest_path,
                [output["id"] for output in outputs],
            )
            response = _response(request, manifest, outputs)
            policy = policy_artifact(
                f"{BACKEND_ID}.lowering",
                target_binding,
                [
                    {"id": output["id"], "from_output": output["id"]}
                    for output in outputs
                ],
                _operations(),
            )
            request_path = root / "request.json"
            response_path = root / "response.json"
            policy_path = root / "policy.json"
            save(request, request_path)
            save(response, response_path)
            save(policy, policy_path)
            module = compile_module(
                source_snapshot,
                manifest_path,
                request_path,
                response_path,
                policy_path,
                target_path,
            )
        record = _record(
            source,
            bio,
            semantics,
            manifest,
            request,
            response,
            target,
            policy,
            module,
            backend_spec,
        )
        bundle = FeatureGraphBundle._from_artifacts(
            backend_spec,
            semantics,
            manifest,
            request,
            response,
            target,
            policy,
            module,
            record,
        )
        for role, artifact in (
            ("backend specification", backend_spec),
            ("provider manifest", manifest),
            ("prediction request", request),
            ("prediction response", response),
            ("target contract", target),
            ("lowering policy", policy),
            ("development module", module),
            ("compilation record", record),
        ):
            _preflight_wire(artifact, role)
        _require_backend_implementation_integrity(backend_spec)
        bundle.to_dict()
        return bundle
    except FeatureGraphError:
        raise
    except (GFF3Error, V2Error, OSError) as failure:
        raise _fail(str(failure)) from failure


def _bundle_payload(
    value: FeatureGraphBundle | dict[str, Any]
) -> dict[str, Any]:
    if isinstance(value, FeatureGraphBundle):
        return value.to_dict()
    if type(value) is dict:
        return copy.deepcopy(value)
    raise _fail("bundle must be a FeatureGraphBundle or bundle artifact object")


def validate_feature_graph_bundle(
    value: FeatureGraphBundle | dict[str, Any],
    source_path: str | Path,
    bio_artifact: GFF3Artifact | dict[str, Any],
    *,
    gff3_source: bytes | None = None,
) -> dict[str, Any]:
    """Replay the backend and require exact equality, including every digest."""

    supplied = _bundle_payload(value)
    try:
        artifact_digest(supplied)
    except ContractError as failure:
        raise _fail(f"invalid bundle artifact digest: {failure}") from failure
    artifacts = supplied.get("artifacts")
    if type(artifacts) is not dict:
        raise _fail("bundle artifacts must be an object")
    for role, artifact in artifacts.items():
        if type(artifact) is not dict:
            raise _fail(f"bundle artifact {role!r} must be an object")
        try:
            artifact_digest(artifact)
        except ContractError as failure:
            raise _fail(f"invalid {role} artifact digest: {failure}") from failure
    backend_spec = artifacts.get("backend_spec")
    if type(backend_spec) is not dict:
        raise _fail("bundle is missing the normative backend specification")
    _require_backend_implementation_integrity(backend_spec)
    expected = compile_feature_graph(
        source_path, bio_artifact, gff3_source=gff3_source
    ).to_dict()
    if supplied != expected:
        raise _fail("bundle does not match deterministic lowering replay")
    return supplied


__all__ = [
    "BACKEND_ID",
    "BACKEND_SPEC",
    "BACKEND_SPEC_SHA256",
    "BACKEND_VERSION",
    "FeatureGraphBundle",
    "FeatureGraphError",
    "backend_spec_artifact",
    "compile_feature_graph",
    "validate_feature_graph_bundle",
]

"""Deterministic GenBank BioIR lowering to a feature-state module.

The adapter is mechanical: every BioIR feature becomes one unit, no edges are
inferred, and the original GenBank artifact remains an external content-bound
source.  The returned bundle owns nine independently sealed child artifacts;
its serialized index contains references only.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
import json
from pathlib import Path
import tempfile
from typing import Any, Iterator

from ._canonical import ContractError, artifact_digest, canonical_bytes, digest
from ._publish import (
    HAS_LINUX_DIRECTORY_PUBLICATION,
    PublicationError,
    publish_directory,
)
from .insdc import (
    AUTHORITY_MANIFEST_SHA256 as GENBANK_AUTHORITY_MANIFEST_SHA256,
    FORMAT as GENBANK_FORMAT,
    PROFILE as GENBANK_PROFILE,
    VERSION as GENBANK_VERSION,
    GenBankArtifact,
    GenBankError,
    load_genbank_artifact,
    validate_genbank_artifact,
)
from .v2 import (
    inline_storage,
    pack,
    policy_artifact,
    save,
    target_artifact,
)
from .v2._common import V2Error, pretty_bytes, seal
from .v2.compiler import (
    COMPILER as MODULE_COMPILER,
    _compile_module_from_validated_genbank_source,
)
from .v2.limits import MAX_JSON_BYTES as MAX_CHILD_BYTES
from .v2.provider import _make_request_from_validated_genbank_source
from .v2.target import DEV_DOMAIN, OP_UNIT_CREATE


BACKEND_ID = "io.github.dttdrv.epigenesis.bio.insdc-graph"
BACKEND_VERSION = 1
BACKEND_SPEC_FORMAT = "brainc.bio.insdc-graph-backend-spec"
SEMANTICS_FORMAT = "brainc.bio.insdc-graph-semantics"
RECORD_FORMAT = "brainc.bio.insdc-graph-compilation"
BUNDLE_FORMAT = "brainc.bio.insdc-graph-bundle"

MAX_U64 = 2**64 - 1
_HAS_DIRECTORY_DESCRIPTOR = HAS_LINUX_DIRECTORY_PUBLICATION
KIND_BITS = (("between", 2), ("interval", 1), ("uncertain-point", 4))
ORIENTATION_BITS = ((-1, 2), (1, 1))
_KIND_BITS_CANONICAL = canonical_bytes(
    [{"key": kind, "bit": bit} for kind, bit in KIND_BITS]
)
_ORIENTATION_BITS_CANONICAL = canonical_bytes(
    [
        {"key": orientation, "bit": bit}
        for orientation, bit in ORIENTATION_BITS
    ]
)
CHILD_ROLES = (
    "backend_spec",
    "backend_semantics",
    "compilation_record",
    "development_module",
    "lowering_policy",
    "prediction_request",
    "prediction_response",
    "provider_manifest",
    "target_contract",
)
BUNDLE_FILENAME = "bundle.json"
CHILD_FILENAMES = tuple((role, f"{role}.json") for role in CHILD_ROLES)


FEATURE_FIELD_SPECS: tuple[dict[str, Any], ...] = (
    {
        "id": "fuzzy_boundary_count",
        "dtype": "u64",
        "tensor_id": "t.feature.fuzzy_boundary_count",
        "formula": {
            "op": "count-non-null",
            "over": "feature.location.segments",
            "paths": ["segment.start_fuzz", "segment.end_fuzz"],
            "excludes": "within-uncertainty-represented-by-uncertain-point-kind",
        },
    },
    {
        "id": "key_code",
        "dtype": "u64",
        "tensor_id": "t.feature.key_code",
        "formula": {
            "op": "dictionary-code",
            "dictionary": "keys",
            "key": {"path": "feature.key"},
        },
    },
    {
        "id": "kind_mask",
        "dtype": "u64",
        "tensor_id": "t.feature.kind_mask",
        "formula": {
            "op": "bitwise-or",
            "over": "feature.location.segments",
            "key": {"path": "segment.kind"},
            "mapping": [
                {"key": kind, "bit": bit} for kind, bit in KIND_BITS
            ],
            "losses": ["multiplicity", "order", "segment-association"],
        },
    },
    {
        "id": "ordinal",
        "dtype": "u64",
        "tensor_id": "t.feature.ordinal",
        "formula": {
            "op": "identity",
            "path": "feature.ordinal",
            "base": 1,
        },
    },
    {
        "id": "orientation_mask",
        "dtype": "u64",
        "tensor_id": "t.feature.orientation_mask",
        "formula": {
            "op": "bitwise-or",
            "over": "feature.location.segments",
            "key": {"path": "segment.orientation"},
            "mapping": [
                {"key": orientation, "bit": bit}
                for orientation, bit in ORIENTATION_BITS
            ],
            "losses": ["multiplicity", "order", "segment-association"],
        },
    },
    {
        "id": "qualifier_count",
        "dtype": "u64",
        "tensor_id": "t.feature.qualifier_count",
        "formula": {"op": "length", "path": "feature.qualifiers"},
    },
    {
        "id": "record_code",
        "dtype": "u64",
        "tensor_id": "t.feature.record_code",
        "formula": {
            "op": "dictionary-code",
            "dictionary": "records",
            "key": {"path": "record.record_id"},
        },
    },
    {
        "id": "remote_segment_count",
        "dtype": "u64",
        "tensor_id": "t.feature.remote_segment_count",
        "formula": {
            "op": "count",
            "over": "feature.location.segments",
            "where": {"op": "is-not-null", "path": "segment.reference"},
        },
    },
    {
        "id": "segment_count",
        "dtype": "u64",
        "tensor_id": "t.feature.segment_count",
        "formula": {
            "op": "length",
            "path": "feature.location.segments",
        },
    },
    {
        "id": "segment_extent_sum",
        "dtype": "u64",
        "tensor_id": "t.feature.segment_extent_sum",
        "formula": {
            "op": "sum",
            "over": "feature.location.segments",
            "value": {
                "op": "subtract",
                "left": {"path": "segment.end"},
                "right": {"path": "segment.start"},
            },
            "coordinate_system": "normalized-0-based-half-open",
            "between_contribution": 0,
            "uncertain_point_contribution": "normalized-envelope-width",
            "overlap_policy": "count-each-segment-extent-with-multiplicity",
        },
    },
    {
        "id": "unresolved_segment_count",
        "dtype": "u64",
        "tensor_id": "t.feature.unresolved_segment_count",
        "formula": {
            "op": "count",
            "over": "feature.location.segments",
            "where": {
                "op": "equals",
                "path": "segment.bounds_status",
                "value": "unresolved",
            },
            "meaning": "reference-resolution-and-bounds-status",
        },
    },
)
FEATURE_FIELDS = tuple(field["id"] for field in FEATURE_FIELD_SPECS)
_FEATURE_FIELD_SPECS_CANONICAL = canonical_bytes(list(FEATURE_FIELD_SPECS))

TARGET_CONTRACT_SPEC = {
    "id": BACKEND_ID,
    "abi_major": 1,
    "opsets": [{"domain": DEV_DOMAIN, "version": 1}],
    "unit_schemas": [
        {
            "id": "feature",
            "fields": [
                {
                    "id": field,
                    "type": {"dtype": "u64", "shape": []},
                    "unit": None,
                    "mutability": "constant",
                    "numeric": {"kind": "integer"},
                }
                for field in FEATURE_FIELDS
            ],
        }
    ],
    "edge_schemas": [],
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
)

PRODUCER = {
    "name": "brainc-insdc-feature-state",
    "version": "0.1.0",
    "passes": [
        "validate-genbank-v2-artifact",
        "bind-whole-source-and-bio-ir",
        "encode-structural-feature-state",
        "compile-development-module",
        "emit-reference-only-bundle-index",
    ],
}

BACKEND_SPEC = {
    "format": BACKEND_SPEC_FORMAT,
    "version": 1,
    "id": BACKEND_ID,
    "backend_version": BACKEND_VERSION,
    "input": {
        "format": GENBANK_FORMAT,
        "version": GENBANK_VERSION,
        "profile": GENBANK_PROFILE,
        "authority_manifest_sha256": GENBANK_AUTHORITY_MANIFEST_SHA256,
        "coordinate_system": "0-based-half-open",
        "binding": "whole-artifact-sha256-and-bio-ir-sha256",
    },
    "projection": {
        "unit": "one-per-ordered-insdc-feature-table-entry",
        "columns": {
            "count": len(FEATURE_FIELD_SPECS),
            "kind": "structural-summary",
        },
        "semantic_authority": {
            "artifact": "bound-genbank-bio-ir",
            "role": "source-dependency",
        },
        "dictionary_codes": "bundle-local-and-may-renumber-with-membership",
        "masks": "discard-multiplicity-order-and-segment-association",
        "fuzzy_boundary_count": "excludes-within-uncertainty",
        "unresolved_segment_count": "reference-resolution-and-bounds-status",
        "edges": "none-inferred-or-emitted",
    },
    "semantics_sidecar": {
        "format": SEMANTICS_FORMAT,
        "version": 1,
        "bio_ir_storage": "reference-only",
        "source_binding_fields": [
            "format",
            "version",
            "artifact_sha256",
            "bio_ir_sha256",
        ],
        "dictionary_wire": {
            "container": "object-keyed-by-dictionary-id",
            "entry": {"code": "zero-based-array-position", "value": "text"},
            "entry_order": "ascending-code",
        },
        "unit_order": "bio-ir-record-array-order-then-one-based-feature-ordinal",
        "unit_identity": "feature.feature_id",
        "edge_order": "empty",
        "counts": ["units", "edges"],
    },
    "dictionaries": [
        {
            "id": "keys",
            "values": {
                "op": "distinct-sort",
                "path": "bio_ir.records[].features[].key",
                "order": "ascending-utf8-bytes",
            },
            "codes": "zero-based-array-position",
            "stability": "bundle-local; may-renumber-when-membership-changes",
        },
        {
            "id": "records",
            "values": {
                "op": "distinct-sort",
                "path": "bio_ir.records[].record_id",
                "order": "ascending-utf8-bytes",
            },
            "codes": "zero-based-array-position",
            "stability": "bundle-local; may-renumber-when-membership-changes",
        },
    ],
    "tensors": {
        "encoding": {
            "byte_order": "little-endian",
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
            "formula": {"op": "length", "path": "semantics.unit_order"},
        },
        "feature_fields": [
            {**copy.deepcopy(field), "shape": ["feature-count"]}
            for field in FEATURE_FIELD_SPECS
        ],
    },
    "target_contract": copy.deepcopy(TARGET_CONTRACT_SPEC),
    "operations": copy.deepcopy(list(OPERATION_SPECS)),
    "development_abi": {
        "domain": DEV_DOMAIN,
        "version": 1,
        "unit_order": "unit-local-id-equals-semantics-unit-order-index",
        "initializer_semantics": {
            "field_order": "exact-target-schema-field-order",
            "field_coverage": "every-field-exactly-once",
            "created_count": "value-of-unitless-u64-count-scalar",
            "field_shape": "[created-count]",
        },
    },
    "artifact_profile": {
        "backend_spec": {"format": BACKEND_SPEC_FORMAT, "version": 1},
        "backend_semantics": {"format": SEMANTICS_FORMAT, "version": 1},
        "provider_manifest": {
            "format": "brainc.provider-manifest",
            "version": 2,
            "provider": {"name": BACKEND_ID, "version": str(BACKEND_VERSION)},
            "model_identity": {
                "kind": "content-sha256",
                "value": "backend-spec-artifact-sha256",
            },
            "accepts": [f"{GENBANK_FORMAT}/v{GENBANK_VERSION}"],
        },
        "prediction_request": {"format": "brainc.prediction-request", "version": 2},
        "prediction_response": {"format": "brainc.prediction-response", "version": 2},
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
            "producer": copy.deepcopy(PRODUCER),
        },
        "bundle": {"format": BUNDLE_FORMAT, "version": 1},
    },
    "capabilities": {
        "edges": False,
        "inference": False,
        "learning": False,
        "ports": False,
        "rules": False,
    },
}
_BACKEND_SPEC_CANONICAL = canonical_bytes(BACKEND_SPEC)
BACKEND_SPEC_SHA256 = digest(BACKEND_SPEC)
_PRODUCER_CANONICAL = canonical_bytes(PRODUCER)


class INSDCGraphError(ContractError):
    """A GenBank feature-state lowering violates its closed ABI."""


def _fail(detail: str) -> INSDCGraphError:
    return INSDCGraphError(f"INSDCGRAPH001: {detail}")


def _fresh_producer() -> dict[str, Any]:
    return json.loads(_PRODUCER_CANONICAL)


def _fresh_field_specs() -> list[dict[str, Any]]:
    return json.loads(_FEATURE_FIELD_SPECS_CANONICAL)


def _fresh_mask(canonical: bytes) -> dict[Any, int]:
    return {entry["key"]: entry["bit"] for entry in json.loads(canonical)}


def backend_spec_artifact() -> dict[str, Any]:
    """Return a fresh sealed copy of the immutable normative backend spec."""

    return seal(json.loads(_BACKEND_SPEC_CANONICAL))


def _require_backend_integrity(artifact: dict[str, Any]) -> None:
    try:
        stored = artifact_digest(artifact)
    except ContractError as failure:
        raise _fail(f"backend specification is invalid: {failure}") from failure
    if stored != BACKEND_SPEC_SHA256:
        raise _fail("backend specification content identity changed")
    if artifact != backend_spec_artifact():
        raise _fail("backend specification differs from its normative content")


def _preflight_child(artifact: dict[str, Any], role: str) -> None:
    try:
        raw = pretty_bytes(artifact)
    except V2Error as failure:
        raise _fail(f"{role} is outside the child wire contract: {failure}") from failure
    if len(raw) > MAX_CHILD_BYTES:
        raise _fail(f"{role} exceeds {MAX_CHILD_BYTES} bytes")


def _source_payload(
    value: str | Path | GenBankArtifact | dict[str, Any],
) -> dict[str, Any]:
    validated = False
    if isinstance(value, GenBankArtifact):
        payload = value.to_dict()
    elif type(value) is dict:
        payload = copy.deepcopy(value)
    elif isinstance(value, (str, Path)):
        payload = load_genbank_artifact(value).to_dict()
        validated = True
    else:
        raise _fail("source must be a GenBank artifact object or artifact path")
    if not validated:
        validate_genbank_artifact(payload)
    return payload


def _source_binding(source: dict[str, Any]) -> dict[str, Any]:
    return {
        "format": source["format"],
        "version": source["version"],
        "artifact_sha256": source["artifact_sha256"],
        "bio_ir_sha256": source["bio_ir_sha256"],
    }


def _generic_source_binding(source: dict[str, Any]) -> dict[str, Any]:
    binding = _source_binding(source)
    return {
        "format": binding["format"],
        "version": binding["version"],
        "artifact_sha256": binding["artifact_sha256"],
        "ir_sha256": binding["bio_ir_sha256"],
    }


def _utf8_sorted(values: set[str]) -> list[str]:
    return sorted(values, key=lambda value: value.encode("utf-8"))


def _indexed(values: list[str]) -> list[dict[str, Any]]:
    return [{"code": index, "value": value} for index, value in enumerate(values)]


def _features(
    source: dict[str, Any],
) -> Iterator[tuple[int, dict[str, Any], dict[str, Any]]]:
    seen_ids: set[str] = set()
    for record_index, record in enumerate(source["bio_ir"]["records"]):
        for feature_index, feature in enumerate(record["features"]):
            if feature["ordinal"] != feature_index + 1:
                raise _fail(
                    f"record {record_index} feature ordinals are not contiguous and one-based"
                )
            feature_id = feature["feature_id"]
            if feature_id in seen_ids:
                raise _fail(f"duplicate feature identity {feature_id!r}")
            seen_ids.add(feature_id)
            yield record_index, record, feature


def _semantics_artifact(source: dict[str, Any]) -> dict[str, Any]:
    rows = list(_features(source))
    keys = _utf8_sorted({feature["key"] for _, _, feature in rows})
    records = _utf8_sorted(
        {record["record_id"] for record in source["bio_ir"]["records"]}
    )
    semantic_core = {
        "source": _source_binding(source),
        "dictionaries": {
            "keys": _indexed(keys),
            "records": _indexed(records),
        },
        "unit_order": [feature["feature_id"] for _, _, feature in rows],
        "edge_order": [],
        "counts": {"units": len(rows), "edges": 0},
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
            **semantic_core,
            "semantics_sha256": digest(semantic_core),
        }
    )


def _dictionary_codes(
    semantics: dict[str, Any], name: str
) -> dict[str, int]:
    return {
        entry["value"]: entry["code"]
        for entry in semantics["dictionaries"][name]
    }


def _checked_u64(value: int, label: str) -> int:
    if type(value) is not int or not 0 <= value <= MAX_U64:
        raise _fail(f"{label} is outside the u64 range")
    return value


def _mask(values: Iterator[Any], mapping: dict[Any, int], label: str) -> int:
    result = 0
    for value in values:
        try:
            result |= mapping[value]
        except KeyError as failure:
            raise _fail(f"{label} has unsupported value {value!r}") from failure
    return result


def _feature_values(
    source: dict[str, Any], semantics: dict[str, Any]
) -> dict[str, list[int]]:
    key_codes = _dictionary_codes(semantics, "keys")
    record_codes = _dictionary_codes(semantics, "records")
    kind_bits = _fresh_mask(_KIND_BITS_CANONICAL)
    orientation_bits = _fresh_mask(_ORIENTATION_BITS_CANONICAL)
    values = {field: [] for field in FEATURE_FIELDS}
    for unit_index, (_, record, feature) in enumerate(_features(source)):
        segments = feature["location"]["segments"]
        row = {
            "fuzzy_boundary_count": sum(
                segment[boundary] is not None
                for segment in segments
                for boundary in ("start_fuzz", "end_fuzz")
            ),
            "key_code": key_codes[feature["key"]],
            "kind_mask": _mask(
                (segment["kind"] for segment in segments),
                kind_bits,
                f"unit {unit_index} segment kind",
            ),
            "orientation_mask": _mask(
                (segment["orientation"] for segment in segments),
                orientation_bits,
                f"unit {unit_index} segment orientation",
            ),
            "ordinal": feature["ordinal"],
            "qualifier_count": len(feature["qualifiers"]),
            "record_code": record_codes[record["record_id"]],
            "remote_segment_count": sum(
                segment["reference"] is not None for segment in segments
            ),
            "segment_count": len(segments),
            "segment_extent_sum": sum(
                segment["end"] - segment["start"] for segment in segments
            ),
            "unresolved_segment_count": sum(
                segment["bounds_status"] == "unresolved" for segment in segments
            ),
        }
        for field in FEATURE_FIELDS:
            values[field].append(
                _checked_u64(row[field], f"unit {unit_index} field {field}")
            )
    return values


def _tensor(identifier: str, shape: list[int], values: list[int]) -> dict[str, Any]:
    return {
        "id": identifier,
        "type": {"dtype": "u64", "shape": shape},
        "unit": None,
        "axes": [None] * len(shape),
        "storage": inline_storage(pack("u64", values)),
    }


def _outputs(
    source: dict[str, Any], semantics: dict[str, Any]
) -> list[dict[str, Any]]:
    count = semantics["counts"]["units"]
    values = _feature_values(source, semantics)
    outputs = [_tensor("t.feature.count", [], [count])]
    outputs.extend(
        _tensor(field["tensor_id"], [count], values[field["id"]])
        for field in _fresh_field_specs()
    )
    return sorted(outputs, key=lambda output: output["id"])


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
            "accepts": [f"{GENBANK_FORMAT}/v{GENBANK_VERSION}"],
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
    request: dict[str, Any],
    manifest: dict[str, Any],
    outputs: list[dict[str, Any]],
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


def _compilation_record(
    source: dict[str, Any], artifacts: dict[str, dict[str, Any]]
) -> dict[str, Any]:
    semantics = artifacts["backend_semantics"]
    module = artifacts["development_module"]
    budgets = module["module"]["budgets"]
    return seal(
        {
            "format": RECORD_FORMAT,
            "version": 1,
            "producer": _fresh_producer(),
            "backend": {
                "id": BACKEND_ID,
                "version": BACKEND_VERSION,
                "spec_sha256": BACKEND_SPEC_SHA256,
                "semantics_sha256": semantics["semantics_sha256"],
            },
            "source": _source_binding(source),
            "artifacts": {
                role: {"artifact_sha256": artifact["artifact_sha256"]}
                for role, artifact in artifacts.items()
            },
            "result": {
                "module_sha256": module["module_sha256"],
                "operations": budgets["operations"],
                "tensor_bytes": budgets["tensor_bytes"],
                "units": budgets["units"],
                "edges": budgets["edges"],
                "attachments": budgets["attachments"],
                "ports": 0,
                "inference": False,
                "learning": False,
            },
        }
    )


@dataclass(frozen=True, init=False)
class INSDCGraphBundle:
    """Nine frozen children plus a reference-only bundle index."""

    _artifact_wires: tuple[tuple[str, bytes], ...] = field(repr=False)

    def __init__(self) -> None:
        raise TypeError("INSDCGraphBundle instances are created by compile_insdc_graph")

    @classmethod
    def _from_artifacts(
        cls, artifacts: dict[str, dict[str, Any]]
    ) -> INSDCGraphBundle:
        if tuple(artifacts) != CHILD_ROLES:
            raise _fail("bundle child roles or order are invalid")
        wires = tuple(
            (role, canonical_bytes(copy.deepcopy(artifact)))
            for role, artifact in artifacts.items()
        )
        bundle = object.__new__(cls)
        object.__setattr__(bundle, "_artifact_wires", wires)
        bundle._validated_artifacts()
        return bundle

    def artifact(self, role: str) -> dict[str, Any]:
        for candidate, wire in self._artifact_wires:
            if candidate == role:
                value = json.loads(wire)
                if type(value) is not dict:
                    raise _fail(f"invalid {role} child artifact")
                return value
        raise _fail(f"unknown child role {role!r}")

    @property
    def artifacts(self) -> dict[str, dict[str, Any]]:
        return {role: self.artifact(role) for role in CHILD_ROLES}

    @property
    def backend_spec(self) -> dict[str, Any]:
        return self.artifact("backend_spec")

    @property
    def semantics(self) -> dict[str, Any]:
        return self.artifact("backend_semantics")

    @property
    def compilation_record(self) -> dict[str, Any]:
        return self.artifact("compilation_record")

    @property
    def development_module(self) -> dict[str, Any]:
        return self.artifact("development_module")

    @property
    def lowering_policy(self) -> dict[str, Any]:
        return self.artifact("lowering_policy")

    @property
    def prediction_request(self) -> dict[str, Any]:
        return self.artifact("prediction_request")

    @property
    def prediction_response(self) -> dict[str, Any]:
        return self.artifact("prediction_response")

    @property
    def provider_manifest(self) -> dict[str, Any]:
        return self.artifact("provider_manifest")

    @property
    def target_contract(self) -> dict[str, Any]:
        return self.artifact("target_contract")

    def _validated_artifacts(self) -> dict[str, dict[str, Any]]:
        artifacts = self.artifacts
        for role, artifact in artifacts.items():
            try:
                artifact_digest(artifact)
            except ContractError as failure:
                raise _fail(f"invalid {role} child digest: {failure}") from failure
            _preflight_child(artifact, role)
        _require_backend_integrity(artifacts["backend_spec"])
        record = artifacts["compilation_record"]
        bound = set(CHILD_ROLES) - {"compilation_record"}
        if set(record["artifacts"]) != bound:
            raise _fail("compilation record child references are incomplete")
        for role in sorted(bound):
            if record["artifacts"][role] != {
                "artifact_sha256": artifacts[role]["artifact_sha256"]
            }:
                raise _fail(f"compilation record does not bind {role}")
        semantics = artifacts["backend_semantics"]
        if record["source"] != semantics["source"]:
            raise _fail("compilation record and semantics source bindings differ")
        expected_generic = {
            "format": record["source"]["format"],
            "version": record["source"]["version"],
            "artifact_sha256": record["source"]["artifact_sha256"],
            "ir_sha256": record["source"]["bio_ir_sha256"],
        }
        if artifacts["prediction_request"]["source"] != expected_generic:
            raise _fail("prediction request does not bind the GenBank BioIR")
        if artifacts["development_module"]["sources"]["sequence"] != expected_generic:
            raise _fail("development module does not bind the GenBank BioIR")
        return artifacts

    def to_dict(self) -> dict[str, Any]:
        """Return the reference-only index manifest, not a standalone bundle."""

        artifacts = self._validated_artifacts()
        record = artifacts["compilation_record"]
        index = seal(
            {
                "format": BUNDLE_FORMAT,
                "version": 1,
                "backend": {
                    "id": BACKEND_ID,
                    "version": BACKEND_VERSION,
                    "spec_sha256": BACKEND_SPEC_SHA256,
                },
                "source": copy.deepcopy(record["source"]),
                "artifacts": {
                    role: {
                        "artifact_sha256": artifacts[role]["artifact_sha256"]
                    }
                    for role in CHILD_ROLES
                },
            }
        )
        _preflight_child(index, "bundle index")
        return index

    def save(self, directory: str | Path) -> dict[str, Path]:
        """Publish the reference-only index and nine children as one directory."""

        artifacts = {
            **{role: self.artifact(role) for role in CHILD_ROLES},
            "bundle": self.to_dict(),
        }
        filenames = {
            **dict(CHILD_FILENAMES),
            "bundle": BUNDLE_FILENAME,
        }
        try:
            paths = publish_directory(
                directory,
                {
                    filenames[role]: artifact
                    for role, artifact in artifacts.items()
                },
                pretty_bytes,
            )
        except PublicationError as failure:
            state = " after commit" if failure.committed else ""
            raise _fail(
                f"cannot publish bundle directory{state}: {failure.strerror}"
            ) from failure
        except V2Error as failure:
            raise _fail(f"cannot encode bundle child: {failure}") from failure
        return {
            "bundle": paths[BUNDLE_FILENAME],
            **{
                role: paths[filename]
                for role, filename in CHILD_FILENAMES
            },
        }


def compile_insdc_graph(
    source_artifact: str | Path | GenBankArtifact | dict[str, Any],
) -> INSDCGraphBundle:
    """Lower a validated GenBank v2 artifact to structural feature state."""

    try:
        source = _source_payload(source_artifact)
        records = {
            member["record_id"]: member["sequence"]["bases"]
            for member in source["sequence_collection"]["members"]
        }
        backend_spec = backend_spec_artifact()
        _require_backend_integrity(backend_spec)
        semantics = _semantics_artifact(source)
        outputs = _outputs(source, semantics)
        manifest = _manifest(outputs)
        target = target_artifact(copy.deepcopy(backend_spec["target_contract"]))
        target_binding = {
            "id": target["contract"]["id"],
            "abi_major": target["contract"]["abi_major"],
            "contract_sha256": target["contract_sha256"],
        }
        with tempfile.TemporaryDirectory(prefix="brainc-insdc-graph-") as temporary:
            root = Path(temporary)
            manifest_path = root / "manifest.json"
            target_path = root / "target.json"
            save(manifest, manifest_path)
            save(target, target_path)
            request = _make_request_from_validated_genbank_source(
                source,
                records,
                manifest_path,
                [output["id"] for output in outputs],
            )
            if request["source"] != _generic_source_binding(source):
                raise _fail("generic source ABI did not bind the GenBank BioIR")
            response = _response(request, manifest, outputs)
            policy = policy_artifact(
                f"{BACKEND_ID}.lowering",
                target_binding,
                [
                    {"id": output["id"], "from_output": output["id"]}
                    for output in outputs
                ],
                copy.deepcopy(backend_spec["operations"]),
            )
            request_path = root / "request.json"
            response_path = root / "response.json"
            policy_path = root / "policy.json"
            save(request, request_path)
            save(response, response_path)
            save(policy, policy_path)
            module = _compile_module_from_validated_genbank_source(
                source,
                records,
                manifest_path,
                request_path,
                response_path,
                policy_path,
                target_path,
            )
        non_record = {
            "backend_spec": backend_spec,
            "backend_semantics": semantics,
            "development_module": module,
            "lowering_policy": policy,
            "prediction_request": request,
            "prediction_response": response,
            "provider_manifest": manifest,
            "target_contract": target,
        }
        record = _compilation_record(source, non_record)
        artifacts = {
            "backend_spec": backend_spec,
            "backend_semantics": semantics,
            "compilation_record": record,
            "development_module": module,
            "lowering_policy": policy,
            "prediction_request": request,
            "prediction_response": response,
            "provider_manifest": manifest,
            "target_contract": target,
        }
        bundle = INSDCGraphBundle._from_artifacts(artifacts)
        bundle.to_dict()
        _require_backend_integrity(backend_spec)
        return bundle
    except INSDCGraphError:
        raise
    except (GenBankError, V2Error, OSError) as failure:
        raise _fail(str(failure)) from failure


__all__ = [
    "BACKEND_ID",
    "BACKEND_SPEC",
    "BACKEND_SPEC_SHA256",
    "BACKEND_VERSION",
    "BUNDLE_FILENAME",
    "BUNDLE_FORMAT",
    "CHILD_FILENAMES",
    "CHILD_ROLES",
    "FEATURE_FIELD_SPECS",
    "INSDCGraphBundle",
    "INSDCGraphError",
    "RECORD_FORMAT",
    "SEMANTICS_FORMAT",
    "backend_spec_artifact",
    "compile_insdc_graph",
]

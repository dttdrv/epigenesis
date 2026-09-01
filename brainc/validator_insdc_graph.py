"""Independent validator for the deterministic GenBank feature-state adapter.

Only Python's standard library and the two independent compiler validators are
allowed here.  Producer modules and the producer's canonicalization helpers are
deliberately outside this validator's trust boundary.
"""

from __future__ import annotations

import argparse
import base64
import errno
import functools
import hashlib
import json
import math
import os
from pathlib import Path
import re
import secrets
import stat
import struct
import sys
from typing import Any


SAFE_INTEGER = 2**53 - 1
MAX_U64 = 2**64 - 1
MAX_JSON_DEPTH = 64
MAX_JSON_MEMBERS = 1_000_000
MAX_TEXT_BYTES = 1024 * 1024
MAX_CHILD_BYTES = 16 * 1024 * 1024
MAX_SOURCE_ARTIFACT_BYTES = 64 * 1024 * 1024
MAX_GENBANK_BYTES = 16 * 1024 * 1024
MAX_REPORT_BYTES = 16 * 1024 * 1024
SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")

BACKEND_ID = "io.github.dttdrv.epigenesis.bio.insdc-graph"
BACKEND_VERSION = 1
BACKEND_SPEC_FORMAT = "brainc.bio.insdc-graph-backend-spec"
SEMANTICS_FORMAT = "brainc.bio.insdc-graph-semantics"
COMPILATION_FORMAT = "brainc.bio.insdc-graph-compilation"
BUNDLE_FORMAT = "brainc.bio.insdc-graph-bundle"
BACKEND_SPEC_SHA256 = "a8f2494dddd259f7f41604001f5dcae9e19b0e491eac7c13aec298f2f99a5da4"
DEV_DOMAIN = "io.github.dttdrv.epigenesis.dev"
OP_UNIT_CREATE = f"{DEV_DOMAIN}.unit.create"

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
BUNDLE_FILENAMES = frozenset(
    {BUNDLE_FILENAME, *(filename for _, filename in CHILD_FILENAMES)}
)
REPORT_CHECKS = (
    "genbank-source-artifact-binding",
    "independent-raw-genbank-replay",
    "normative-backend-specification-commitment",
    "utf8-dictionaries-and-record-ordinal-unit-order",
    "eleven-field-u64-feature-state-replay",
    "little-endian-inline-tensor-replay",
    "zero-edge-target-and-lowering-contract",
    "generic-module-and-compilation-reference-replay",
    "reference-only-bundle-index",
)

_HAS_DIRECTORY_DESCRIPTOR = (
    os.name == "posix"
    and os.listdir in getattr(os, "supports_fd", set())
    and all(
        operation in getattr(os, "supports_dir_fd", set())
        for operation in (os.link, os.open, os.stat, os.unlink)
    )
    and getattr(os, "O_DIRECTORY", 0) != 0
    and getattr(os, "O_NOFOLLOW", 0) != 0
)

MODULE_PRODUCER = {
    "name": "brainc",
    "version": "0.5.0",
    "passes": [
        "validate-sequence-source",
        "bind-provider-tensor-contract",
        "validate-target-contract",
        "validate-lowering-policy",
        "type-check-tensors",
        "validate-development-operations",
        "compute-exact-budgets",
        "emit-development-module",
    ],
}

GRAPH_PRODUCER = {
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

FEATURE_FIELDS = (
    "fuzzy_boundary_count",
    "key_code",
    "kind_mask",
    "ordinal",
    "orientation_mask",
    "qualifier_count",
    "record_code",
    "remote_segment_count",
    "segment_count",
    "segment_extent_sum",
    "unresolved_segment_count",
)
KIND_BITS = {"interval": 1, "between": 2, "uncertain-point": 4}
ORIENTATION_BITS = {1: 1, -1: 2}


class INSDCGraphValidationError(ValueError):
    """The adapter bundle differs from independent deterministic replay."""


def _fail(detail: str) -> INSDCGraphValidationError:
    return INSDCGraphValidationError(f"INSDCGRAPHVAL001: {detail}")


def _bounded_validation(function: Any) -> Any:
    @functools.wraps(function)
    def guarded(*args: Any, **kwargs: Any) -> Any:
        try:
            return function(*args, **kwargs)
        except INSDCGraphValidationError:
            raise
        except MemoryError as failure:
            raise _fail("validation memory ceiling exceeded") from failure

    return guarded


def _jcs_string(value: str) -> str:
    if any(0xD800 <= ord(character) <= 0xDFFF for character in value):
        raise ValueError("lone Unicode surrogate")
    return json.dumps(value, ensure_ascii=False, allow_nan=False)


def _jcs_number(value: float) -> str:
    if not math.isfinite(value):
        raise ValueError("non-finite number")
    if value == 0:
        return "0"
    sign = "-" if value < 0 else ""
    mantissa, separator, exponent_text = repr(abs(value)).lower().partition("e")
    exponent = int(exponent_text) if separator else 0
    digits = mantissa.replace(".", "")
    point = (mantissa.find(".") if "." in mantissa else len(mantissa)) + exponent
    while len(digits) > 1 and digits[0] == "0":
        digits = digits[1:]
        point -= 1
    if 1e-6 <= abs(value) < 1e21:
        if point <= 0:
            rendered = "0." + "0" * -point + digits
        elif point >= len(digits):
            rendered = digits + "0" * (point - len(digits))
        else:
            rendered = digits[:point] + "." + digits[point:]
        if "." in rendered:
            rendered = rendered.rstrip("0").rstrip(".")
        return sign + rendered
    digits = digits.rstrip("0")
    normalized_exponent = point - 1
    rendered = digits[0] + ("." + digits[1:] if len(digits) > 1 else "")
    return sign + rendered + ("e+" if normalized_exponent >= 0 else "e") + str(
        normalized_exponent
    )


def _jcs(value: Any) -> str:
    if value is None:
        return "null"
    if type(value) is bool:
        return "true" if value else "false"
    if type(value) is str:
        return _jcs_string(value)
    if type(value) is int:
        if not -SAFE_INTEGER <= value <= SAFE_INTEGER:
            raise ValueError("unsafe integer")
        return str(value)
    if type(value) is float:
        return _jcs_number(value)
    if type(value) is list:
        return "[" + ",".join(_jcs(item) for item in value) + "]"
    if type(value) is dict:
        if any(type(key) is not str for key in value):
            raise TypeError("non-string object key")
        keys = sorted(value, key=lambda key: key.encode("utf-16be"))
        return "{" + ",".join(
            _jcs_string(key) + ":" + _jcs(value[key]) for key in keys
        ) + "}"
    raise TypeError(f"unsupported value type {type(value).__name__}")


def canonical_bytes(value: Any) -> bytes:
    """Encode one I-JSON value with this validator's RFC 8785 implementation."""

    try:
        return _jcs(value).encode("utf-8")
    except MemoryError as failure:
        raise _fail(
            "canonical JSON exceeds the validation memory ceiling"
        ) from failure
    except (TypeError, ValueError, UnicodeError, RecursionError) as failure:
        raise _fail(f"value is not RFC 8785 canonical JSON: {failure}") from failure


def digest(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def _seal(core: dict[str, Any]) -> dict[str, Any]:
    return {**core, "artifact_sha256": digest(core)}


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


def _integer(
    value: Any, label: str, *, minimum: int = 0, maximum: int = SAFE_INTEGER
) -> int:
    if type(value) is not int or not minimum <= value <= maximum:
        raise _fail(f"{label} must be an integer in [{minimum}, {maximum}]")
    return value


def _u64(value: Any, label: str) -> int:
    if type(value) is not int or not 0 <= value <= MAX_U64:
        raise _fail(f"{label} must be an integer in [0, {MAX_U64}]")
    return value


def _text(value: Any, label: str) -> str:
    if type(value) is not str or not value:
        raise _fail(f"{label} must be nonempty text")
    if any(0xD800 <= ord(character) <= 0xDFFF for character in value):
        raise _fail(f"{label} contains a lone Unicode surrogate")
    if len(value.encode("utf-8")) > MAX_TEXT_BYTES:
        raise _fail(f"{label} exceeds the text byte ceiling")
    return value


def _sha(value: Any, label: str) -> str:
    if type(value) is not str or SHA256_RE.fullmatch(value) is None:
        raise _fail(f"{label} must be a lowercase SHA-256 digest")
    return value


def _json_string_size(value: str, label: str) -> int:
    try:
        utf8_bytes = len(value.encode("utf-8"))
    except UnicodeEncodeError as failure:
        raise _fail(f"{label} contains a lone Unicode surrogate") from failure
    if utf8_bytes > MAX_TEXT_BYTES:
        raise _fail(f"{label} exceeds the text byte ceiling")
    return len(
        json.dumps(value, ensure_ascii=False, allow_nan=False).encode("utf-8")
    )


def _validate_tree(
    value: Any,
    label: str,
    *,
    allow_float: bool = False,
    maximum_bytes: int | None = None,
) -> int:
    members = 0
    total_bytes = 0

    def add(size: int) -> None:
        nonlocal total_bytes
        total_bytes += size
        if maximum_bytes is not None and total_bytes > maximum_bytes:
            raise _fail(f"{label} exceeds {maximum_bytes} canonical JSON bytes")

    def visit(current: Any, depth: int) -> None:
        nonlocal members
        if depth > MAX_JSON_DEPTH:
            raise _fail(f"{label} exceeds JSON depth {MAX_JSON_DEPTH}")
        if current is None:
            add(4)
        elif type(current) is bool:
            add(4 if current else 5)
        elif type(current) is int:
            if not -SAFE_INTEGER <= current <= SAFE_INTEGER:
                raise _fail(f"{label} contains an unsafe JSON integer")
            add(len(str(current)))
        elif type(current) is float:
            if not math.isfinite(current):
                raise _fail(f"{label} contains a non-finite number")
            if not allow_float:
                raise _fail(f"{label} contains a floating-point JSON number")
            add(len(_jcs_number(current)))
        elif type(current) is str:
            add(_json_string_size(current, label))
        elif type(current) is list:
            members += len(current)
            if members > MAX_JSON_MEMBERS:
                raise _fail(f"{label} exceeds JSON member ceiling {MAX_JSON_MEMBERS}")
            add(2 + max(0, len(current) - 1))
            for item in current:
                visit(item, depth + 1)
        elif type(current) is dict:
            members += len(current)
            if members > MAX_JSON_MEMBERS:
                raise _fail(f"{label} exceeds JSON member ceiling {MAX_JSON_MEMBERS}")
            add(2 + max(0, len(current) - 1) + len(current))
            for key, item in current.items():
                if type(key) is not str:
                    raise _fail(f"{label} contains a non-string key")
                if depth + 1 > MAX_JSON_DEPTH:
                    raise _fail(f"{label} exceeds JSON depth {MAX_JSON_DEPTH}")
                add(_json_string_size(key, label))
                visit(item, depth + 1)
        else:
            raise _fail(f"{label} contains unsupported type {type(current).__name__}")

    try:
        visit(value, 0)
    except MemoryError as failure:
        raise _fail(f"{label} exceeds the validation memory ceiling") from failure
    return total_bytes


def _artifact(value: Any, label: str) -> dict[str, Any]:
    if type(value) is not dict:
        raise _fail(f"{label} must be an object")
    try:
        _validate_tree(value, label, maximum_bytes=MAX_CHILD_BYTES)
        stored = _sha(value.get("artifact_sha256"), f"{label}.artifact_sha256")
        expected = digest(
            {key: item for key, item in value.items() if key != "artifact_sha256"}
        )
    except MemoryError as failure:
        raise _fail(f"{label} exceeds the validation memory ceiling") from failure
    if stored != expected:
        raise _fail(f"{label}.artifact_sha256 does not match canonical content")
    return value


def _source_artifact(value: Any) -> dict[str, Any]:
    if type(value) is not dict:
        raise _fail("GenBank source must be an object")
    try:
        _validate_tree(
            value,
            "GenBank source",
            allow_float=True,
            maximum_bytes=MAX_SOURCE_ARTIFACT_BYTES,
        )
        source = _keys(
            value,
            {
                "format",
                "version",
                "profile",
                "authority",
                "compiler",
                "sequence_collection",
                "bio_ir",
                "bio_ir_sha256",
                "artifact_sha256",
            },
            "GenBank source",
        )
        if (
            source["format"] != "brainc.bio.insdc-genbank-ir"
            or type(source["version"]) is not int
            or source["version"] != 2
        ):
            raise _fail("unsupported GenBank source artifact identity")
        stored_bio = _sha(source["bio_ir_sha256"], "GenBank source.bio_ir_sha256")
        if stored_bio != digest(source["bio_ir"]):
            raise _fail("GenBank source.bio_ir_sha256 does not match BioIR")
        stored_artifact = _sha(
            source["artifact_sha256"], "GenBank source.artifact_sha256"
        )
        if stored_artifact != digest(
            {key: item for key, item in source.items() if key != "artifact_sha256"}
        ):
            raise _fail("GenBank source.artifact_sha256 does not match canonical content")
        return source
    except MemoryError as failure:
        raise _fail("GenBank source exceeds the validation memory ceiling") from failure


def _backend_spec(value: Any) -> dict[str, Any]:
    specification = _artifact(value, "backend specification")
    if (
        specification.get("format") != BACKEND_SPEC_FORMAT
        or type(specification.get("version")) is not int
        or specification["version"] != 1
        or specification.get("id") != BACKEND_ID
        or type(specification.get("backend_version")) is not int
        or specification["backend_version"] != BACKEND_VERSION
    ):
        raise _fail("unsupported backend specification identity")
    if specification["artifact_sha256"] != BACKEND_SPEC_SHA256:
        raise _fail(
            "backend specification does not match the independent normative commitment"
        )
    return specification


def _dictionary(values: set[str]) -> list[dict[str, Any]]:
    ordered = sorted(values, key=lambda value: value.encode("utf-8"))
    return [{"code": code, "value": value} for code, value in enumerate(ordered)]


def _codes(entries: list[dict[str, Any]]) -> dict[str, int]:
    return {entry["value"]: entry["code"] for entry in entries}


def _checked_add(total: int, value: int, label: str) -> int:
    _integer(value, label)
    result = total + value
    if result > MAX_U64:
        raise _fail(f"{label} sum exceeds u64")
    return result


def _feature_state(bio_ir: dict[str, Any]) -> dict[str, Any]:
    """Reconstruct dictionaries, ordered unit ids, and all 11 u64 columns."""

    records = bio_ir.get("records")
    if type(records) is not list:
        raise _fail("BioIR records must be an array")
    record_ids: set[str] = set()
    keys: set[str] = set()
    feature_ids: set[str] = set()
    ordered: list[tuple[str, dict[str, Any]]] = []
    for record_index, record in enumerate(records):
        if type(record) is not dict:
            raise _fail(f"BioIR record {record_index} must be an object")
        record_id = _text(record.get("record_id"), f"record {record_index}.record_id")
        if record_id in record_ids:
            raise _fail("BioIR record ids must be unique")
        record_ids.add(record_id)
        features = record.get("features")
        if type(features) is not list:
            raise _fail(f"record {record_index}.features must be an array")
        for feature_index, feature in enumerate(features):
            if type(feature) is not dict:
                raise _fail(
                    f"record {record_index}.features[{feature_index}] must be an object"
                )
            ordinal = _integer(
                feature.get("ordinal"),
                f"record {record_index}.features[{feature_index}].ordinal",
                minimum=1,
            )
            if ordinal != feature_index + 1:
                raise _fail("feature unit order must be record-array order then ordinal")
            feature_id = _text(
                feature.get("feature_id"),
                f"record {record_index}.features[{feature_index}].feature_id",
            )
            if feature_id in feature_ids:
                raise _fail("feature identities must be globally unique")
            feature_ids.add(feature_id)
            keys.add(_text(feature.get("key"), "feature key"))
            ordered.append((record_id, feature))

    dictionaries = {
        "keys": _dictionary(keys),
        "records": _dictionary(record_ids),
    }
    key_codes = _codes(dictionaries["keys"])
    record_codes = _codes(dictionaries["records"])
    values: dict[str, list[int]] = {field: [] for field in FEATURE_FIELDS}
    ordered_feature_ids: list[str] = []

    for unit_index, (record_id, feature) in enumerate(ordered):
        ordered_feature_ids.append(feature["feature_id"])
        qualifiers = feature.get("qualifiers")
        location = feature.get("location")
        if type(qualifiers) is not list or type(location) is not dict:
            raise _fail(f"unit {unit_index} qualifiers/location are invalid")
        segments = location.get("segments")
        if type(segments) is not list or not segments:
            raise _fail(f"unit {unit_index} segments must be a nonempty array")

        fuzzy_boundaries = 0
        kind_mask = 0
        orientation_mask = 0
        remote_segments = 0
        segment_extent_sum = 0
        unresolved_segments = 0
        for segment_index, segment in enumerate(segments):
            if type(segment) is not dict:
                raise _fail(f"unit {unit_index} segment {segment_index} must be an object")
            kind = segment.get("kind")
            orientation = segment.get("orientation")
            if kind not in KIND_BITS:
                raise _fail(f"unit {unit_index} segment {segment_index} kind is unsupported")
            if orientation not in ORIENTATION_BITS:
                raise _fail(
                    f"unit {unit_index} segment {segment_index} orientation is unsupported"
                )
            kind_mask |= KIND_BITS[kind]
            orientation_mask |= ORIENTATION_BITS[orientation]
            fuzzy_boundaries += int(segment.get("start_fuzz") is not None)
            fuzzy_boundaries += int(segment.get("end_fuzz") is not None)
            remote_segments += int(segment.get("reference") is not None)
            unresolved_segments += int(segment.get("bounds_status") == "unresolved")
            start = _integer(segment.get("start"), "segment start")
            end = _integer(segment.get("end"), "segment end")
            if end < start:
                raise _fail(f"unit {unit_index} segment {segment_index} has reversed bounds")
            extent = end - start
            segment_extent_sum = _checked_add(
                segment_extent_sum, extent, "segment_extent_sum"
            )

        row = {
            "fuzzy_boundary_count": fuzzy_boundaries,
            "key_code": key_codes[feature["key"]],
            "kind_mask": kind_mask,
            "orientation_mask": orientation_mask,
            "ordinal": feature["ordinal"],
            "qualifier_count": len(qualifiers),
            "record_code": record_codes[record_id],
            "remote_segment_count": remote_segments,
            "segment_extent_sum": segment_extent_sum,
            "segment_count": len(segments),
            "unresolved_segment_count": unresolved_segments,
        }
        for field in FEATURE_FIELDS:
            values[field].append(_u64(row[field], f"unit {unit_index}.{field}"))

    return {
        "dictionaries": dictionaries,
        "feature_ids": ordered_feature_ids,
        "values": values,
    }


def _pack_u64(values: list[int]) -> bytes:
    for index, value in enumerate(values):
        _u64(value, f"u64 value {index}")
    return struct.pack("<" + "Q" * len(values), *values)


def _storage(raw: bytes) -> dict[str, Any]:
    return {
        "kind": "inline-base64",
        "data": base64.b64encode(raw).decode("ascii"),
        "byte_length": len(raw),
        "sha256": hashlib.sha256(raw).hexdigest(),
    }


def _tensor(identifier: str, values: list[int]) -> dict[str, Any]:
    return {
        "id": identifier,
        "type": {"dtype": "u64", "shape": [len(values)]},
        "unit": None,
        "axes": [None],
        "storage": _storage(_pack_u64(values)),
    }


def _scalar_tensor(identifier: str, value: int) -> dict[str, Any]:
    return {
        "id": identifier,
        "type": {"dtype": "u64", "shape": []},
        "unit": None,
        "axes": [],
        "storage": _storage(_pack_u64([value])),
    }


def _source_binding(genbank: dict[str, Any]) -> dict[str, Any]:
    if (
        genbank.get("format") != "brainc.bio.insdc-genbank-ir"
        or type(genbank.get("version")) is not int
        or genbank["version"] != 2
    ):
        raise _fail("unsupported GenBank source artifact identity")
    return {
        "format": genbank["format"],
        "version": genbank["version"],
        "artifact_sha256": _sha(
            genbank.get("artifact_sha256"), "GenBank source artifact_sha256"
        ),
        "bio_ir_sha256": _sha(
            genbank.get("bio_ir_sha256"), "GenBank source bio_ir_sha256"
        ),
    }


def _generic_source_binding(genbank: dict[str, Any]) -> dict[str, Any]:
    binding = _source_binding(genbank)
    return {
        "format": binding["format"],
        "version": binding["version"],
        "artifact_sha256": binding["artifact_sha256"],
        "ir_sha256": binding["bio_ir_sha256"],
    }


def _expected_semantics(
    genbank: dict[str, Any], backend_spec_sha256: str
) -> tuple[dict[str, Any], dict[str, list[int]]]:
    state = _feature_state(genbank["bio_ir"])
    body = {
        "source": _source_binding(genbank),
        "dictionaries": state["dictionaries"],
        "unit_order": state["feature_ids"],
        "edge_order": [],
        "counts": {"units": len(state["feature_ids"]), "edges": 0},
    }
    return (
        _seal(
            {
                "format": SEMANTICS_FORMAT,
                "version": 1,
                "backend": {
                    "id": BACKEND_ID,
                    "version": BACKEND_VERSION,
                    "spec_sha256": backend_spec_sha256,
                },
                **body,
                "semantics_sha256": digest(body),
            }
        ),
        state["values"],
    )


def _outputs(
    semantics: dict[str, Any], values: dict[str, list[int]]
) -> list[dict[str, Any]]:
    units = semantics["counts"]["units"]
    if any(len(values[field]) != units for field in FEATURE_FIELDS):
        raise _fail("feature tensor columns do not match semantic unit count")
    result = [_scalar_tensor("t.feature.count", units)]
    result.extend(_tensor(f"t.feature.{field}", values[field]) for field in FEATURE_FIELDS)
    return sorted(result, key=lambda output: output["id"])


def _field(identifier: str) -> dict[str, Any]:
    return {
        "id": identifier,
        "type": {"dtype": "u64", "shape": []},
        "unit": None,
        "mutability": "constant",
        "numeric": {"kind": "integer"},
    }


def _target() -> dict[str, Any]:
    contract = {
        "id": BACKEND_ID,
        "abi_major": 1,
        "opsets": [{"domain": DEV_DOMAIN, "version": 1}],
        "unit_schemas": [
            {"id": "feature", "fields": [_field(field) for field in FEATURE_FIELDS]}
        ],
        "edge_schemas": [],
        "ports": [],
        "rules": [],
    }
    return _seal(
        {
            "format": "brainc.target-contract",
            "version": 1,
            "contract": contract,
            "contract_sha256": digest(contract),
        }
    )


def _operations() -> list[dict[str, Any]]:
    return [
        {
            "id": "create.features",
            "op": OP_UNIT_CREATE,
            "version": 1,
            "schema": "feature",
            "count": "t.feature.count",
            "initializers": [
                {"field": field, "tensor": f"t.feature.{field}"}
                for field in FEATURE_FIELDS
            ],
        }
    ]


def _output_contract(output: dict[str, Any]) -> dict[str, Any]:
    return {
        key: output[key]
        for key in ("id", "type", "unit", "axes")
    }


def _expected_generic_chain(
    genbank: dict[str, Any],
    outputs: list[dict[str, Any]],
    backend_spec_sha256: str,
) -> dict[str, dict[str, Any]]:
    contracts = [_output_contract(output) for output in outputs]
    manifest = _seal(
        {
            "format": "brainc.provider-manifest",
            "version": 2,
            "provider": {"name": BACKEND_ID, "version": "1"},
            "model_identity": {
                "kind": "content-sha256",
                "value": backend_spec_sha256,
            },
            "accepts": ["brainc.bio.insdc-genbank-ir/v2"],
            "outputs": contracts,
        }
    )
    request = _seal(
        {
            "format": "brainc.prediction-request",
            "version": 2,
            "source": _generic_source_binding(genbank),
            "provider_manifest_sha256": manifest["artifact_sha256"],
            "requested_outputs": contracts,
        }
    )
    response = _seal(
        {
            "format": "brainc.prediction-response",
            "version": 2,
            "request_artifact_sha256": request["artifact_sha256"],
            "provider": manifest["provider"],
            "model_identity": manifest["model_identity"],
            "outputs": outputs,
        }
    )
    target = _target()
    target_binding = {
        "id": BACKEND_ID,
        "abi_major": 1,
        "contract_sha256": target["contract_sha256"],
    }
    operations = _operations()
    policy = _seal(
        {
            "format": "brainc.lowering-policy",
            "version": 2,
            "id": f"{BACKEND_ID}.lowering",
            "target": target_binding,
            "values": [
                {"id": output["id"], "from_output": output["id"]}
                for output in outputs
            ],
            "operations": operations,
        }
    )
    tensors = [
        {
            **output,
            "lineage": {
                "response_output": output["id"],
                "lowering_value": output["id"],
            },
        }
        for output in outputs
    ]
    count_raw = base64.b64decode(
        next(output for output in outputs if output["id"] == "t.feature.count")[
            "storage"
        ]["data"],
        validate=True,
    )
    budgets = {
        "operations": 1,
        "tensor_bytes": sum(output["storage"]["byte_length"] for output in outputs),
        "units": struct.unpack("<Q", count_raw)[0],
        "edges": 0,
        "attachments": 0,
    }
    body = {
        "target": target_binding,
        "requirements": [{"domain": DEV_DOMAIN, "version": 1}],
        "budgets": budgets,
        "tensors": tensors,
        "entrypoint": {"name": "develop", "operations": operations},
    }
    module = _seal(
        {
            "format": "brainc.development-module",
            "version": 1,
            "producer": MODULE_PRODUCER,
            "sources": {
                "sequence": _generic_source_binding(genbank),
                "provider_manifest": {"artifact_sha256": manifest["artifact_sha256"]},
                "prediction_request": {"artifact_sha256": request["artifact_sha256"]},
                "prediction_response": {"artifact_sha256": response["artifact_sha256"]},
                "lowering_policy": {"artifact_sha256": policy["artifact_sha256"]},
                "target_contract": {"artifact_sha256": target["artifact_sha256"]},
            },
            "module": body,
            "module_sha256": digest(body),
        }
    )
    return {
        "development_module": module,
        "lowering_policy": policy,
        "prediction_request": request,
        "prediction_response": response,
        "provider_manifest": manifest,
        "target_contract": target,
    }


def _expected_compilation_record(
    source: dict[str, Any], artifacts: dict[str, dict[str, Any]]
) -> dict[str, Any]:
    semantics = artifacts["backend_semantics"]
    module = artifacts["development_module"]
    budgets = module["module"]["budgets"]
    return _seal(
        {
            "format": COMPILATION_FORMAT,
            "version": 1,
            "producer": GRAPH_PRODUCER,
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


def _expected_artifacts(
    source: dict[str, Any], backend_spec: dict[str, Any]
) -> dict[str, dict[str, Any]]:
    semantics, values = _expected_semantics(source, BACKEND_SPEC_SHA256)
    outputs = _outputs(semantics, values)
    generic = _expected_generic_chain(source, outputs, BACKEND_SPEC_SHA256)
    non_record = {
        "backend_spec": backend_spec,
        "backend_semantics": semantics,
        **generic,
    }
    record = _expected_compilation_record(source, non_record)
    return {
        "backend_spec": backend_spec,
        "backend_semantics": semantics,
        "compilation_record": record,
        **generic,
    }


def _expected_bundle_index(
    source: dict[str, Any], artifacts: dict[str, dict[str, Any]]
) -> dict[str, Any]:
    return _seal(
        {
            "format": BUNDLE_FORMAT,
            "version": 1,
            "backend": {
                "id": BACKEND_ID,
                "version": BACKEND_VERSION,
                "spec_sha256": BACKEND_SPEC_SHA256,
            },
            "source": _source_binding(source),
            "artifacts": {
                role: {"artifact_sha256": artifacts[role]["artifact_sha256"]}
                for role in CHILD_ROLES
            },
        }
    )


@_bounded_validation
def validate_insdc_graph_bundle(
    bundle_index: dict[str, Any],
    artifacts_by_role: dict[str, dict[str, Any]],
    genbank_artifact: dict[str, Any],
    *,
    genbank_source: bytes | None = None,
) -> dict[str, Any]:
    """Replay one reference-only GenBank feature-state bundle independently."""

    source = _source_artifact(genbank_artifact)
    if genbank_source is not None and (
        type(genbank_source) is not bytes or len(genbank_source) > MAX_GENBANK_BYTES
    ):
        raise _fail(
            f"GenBank source evidence must be bytes at most {MAX_GENBANK_BYTES} bytes"
        )
    try:
        from . import validator_v2

        validator_v2.validate_source_artifact(
            source,
            raw_source=genbank_source,
        )
    except (OSError, ValueError) as failure:
        raise _fail(f"independent GenBank source replay failed: {failure}") from failure

    if (
        type(artifacts_by_role) is not dict
        or len(artifacts_by_role) != len(CHILD_ROLES)
        or any(role not in artifacts_by_role for role in CHILD_ROLES)
    ):
        raise _fail("child artifact roles are incomplete or unknown")
    supplied_artifacts = {
        role: _artifact(artifacts_by_role[role], f"child artifact {role}")
        for role in CHILD_ROLES
    }
    try:
        validator_v2.validate_manifest_artifact(
            supplied_artifacts["provider_manifest"]
        )
        validator_v2.validate_request_artifact(
            supplied_artifacts["prediction_request"]
        )
        validator_v2.validate_response_artifact(
            supplied_artifacts["prediction_response"],
            sequence_sha256=source["artifact_sha256"],
        )
        validator_v2.validate_policy_artifact(
            supplied_artifacts["lowering_policy"]
        )
        validator_v2.validate_target_artifact(
            supplied_artifacts["target_contract"]
        )
        validator_v2.validate_module_artifact(
            supplied_artifacts["development_module"],
            supplied_artifacts["target_contract"],
        )
    except (OSError, ValueError) as failure:
        raise _fail(f"independent generic-chain validation failed: {failure}") from failure
    backend_spec = _backend_spec(supplied_artifacts["backend_spec"])
    expected_artifacts = _expected_artifacts(source, backend_spec)
    for role in CHILD_ROLES:
        actual, expected = supplied_artifacts[role], expected_artifacts[role]
        if not _exact_json_equal(actual, expected):
            raise _fail(
                f"{role} differs from independent adapter replay at "
                + _first_difference(actual, expected)
            )

    supplied_index = _artifact(bundle_index, "bundle index")
    expected_index = _expected_bundle_index(source, expected_artifacts)
    if not _exact_json_equal(supplied_index, expected_index):
        raise _fail(
            "bundle index differs from independent reference replay at "
            + _first_difference(supplied_index, expected_index)
        )

    semantics = expected_artifacts["backend_semantics"]
    module = expected_artifacts["development_module"]
    core = {
        "format": "brainc.bio.insdc-graph-validation-report",
        "version": 1,
        "stage": "genbank-feature-state-adapter",
        "valid": True,
        "inputs": {
            "genbank_artifact_sha256": source["artifact_sha256"],
            "bio_ir_sha256": source["bio_ir_sha256"],
            "bundle_artifact_sha256": expected_index["artifact_sha256"],
            "backend_spec_sha256": BACKEND_SPEC_SHA256,
        },
        "result": {
            **module["module"]["budgets"],
            "semantics_sha256": semantics["semantics_sha256"],
            "module_sha256": module["module_sha256"],
        },
        "checks": list(REPORT_CHECKS),
    }
    return {**core, "report_sha256": digest(core)}


def _exact_json_equal(left: Any, right: Any) -> bool:
    pending = [(left, right)]
    while pending:
        left_item, right_item = pending.pop()
        if type(left_item) is not type(right_item):
            return False
        if type(left_item) is dict:
            if left_item.keys() != right_item.keys():
                return False
            pending.extend(
                (left_item[key], right_item[key]) for key in left_item
            )
        elif type(left_item) is list:
            if len(left_item) != len(right_item):
                return False
            pending.extend(zip(left_item, right_item))
        elif left_item != right_item:
            return False
    return True


def _first_difference(actual: Any, expected: Any, path: str = "$") -> str:
    if type(actual) is not type(expected):
        return f"{path} type {type(actual).__name__} != {type(expected).__name__}"
    if type(actual) is dict:
        actual_keys, expected_keys = set(actual), set(expected)
        if actual_keys != expected_keys:
            return (
                f"{path} keys differ; missing={sorted(expected_keys - actual_keys)}, "
                f"unknown={sorted(actual_keys - expected_keys)}"
            )
        for key in sorted(
            actual,
            key=lambda item: (item.endswith("sha256"), item.encode("utf-8")),
        ):
            left, right = actual[key], expected[key]
            if not _exact_json_equal(left, right):
                return _first_difference(left, right, f"{path}.{key}")
    elif type(actual) is list:
        if len(actual) != len(expected):
            return f"{path} length {len(actual)} != {len(expected)}"
        for index, (left, right) in enumerate(zip(actual, expected)):
            if not _exact_json_equal(left, right):
                return _first_difference(left, right, f"{path}[{index}]")
    elif actual != expected:
        return f"{path} value differs"
    return path


def _read_descriptor(
    descriptor: int,
    opened: os.stat_result,
    label: str,
    maximum: int,
) -> bytes:
    if not stat.S_ISREG(opened.st_mode):
        raise _fail(f"{label} must be a regular non-linked file")
    if opened.st_size > maximum:
        raise _fail(f"{label} exceeds {maximum} bytes")
    chunks: list[bytes] = []
    total = 0
    while True:
        chunk = os.read(descriptor, min(64 * 1024, maximum + 1 - total))
        if not chunk:
            break
        chunks.append(chunk)
        total += len(chunk)
        if total > maximum:
            raise _fail(f"{label} exceeds {maximum} bytes")
    finished = os.fstat(descriptor)
    stable = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
    if any(getattr(opened, field) != getattr(finished, field) for field in stable):
        raise _fail(f"{label} changed while being read")
    if total != opened.st_size:
        raise _fail(f"{label} length changed while being read")
    return b"".join(chunks)


def _read_regular(path: str | Path, label: str, maximum: int) -> bytes:
    descriptor = -1
    try:
        source = Path(path)
        inspected = source.lstat()
        if stat.S_ISLNK(inspected.st_mode) or not stat.S_ISREG(inspected.st_mode):
            raise _fail(f"{label} must be a regular non-linked file")
        flags = os.O_RDONLY
        flags |= getattr(os, "O_CLOEXEC", 0)
        flags |= getattr(os, "O_NOFOLLOW", 0)
        flags |= getattr(os, "O_NONBLOCK", 0)
        descriptor = os.open(source, flags)
        opened = os.fstat(descriptor)
        if (opened.st_dev, opened.st_ino) != (inspected.st_dev, inspected.st_ino):
            raise _fail(f"{label} changed while being opened")
        raw = _read_descriptor(descriptor, opened, label, maximum)
        attached = source.lstat()
        stable = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
        if stat.S_ISLNK(attached.st_mode) or any(
            getattr(opened, field) != getattr(attached, field) for field in stable
        ):
            raise _fail(f"{label} path changed while being read")
        return raw
    except INSDCGraphValidationError:
        raise
    except (OSError, TypeError, ValueError) as failure:
        raise _fail(f"cannot read {label}: {failure}") from failure
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _closed_names(names: set[str]) -> None:
    missing = sorted(BUNDLE_FILENAMES - names)
    unknown = sorted(names - BUNDLE_FILENAMES)
    if missing or unknown:
        raise _fail(
            "bundle directory is not closed; "
            f"missing={missing or 'none'}, unknown={unknown or 'none'}"
        )


def _read_bundle_directory(directory: str | Path) -> dict[str, bytes]:
    """Snapshot the exact ten-file bundle convention through one directory fd."""

    descriptor = -1
    try:
        if not _HAS_DIRECTORY_DESCRIPTOR:
            raise _fail("secure bundle directory input is unsupported on this platform")
        source = Path(directory)
        inspected = source.lstat()
        if stat.S_ISLNK(inspected.st_mode) or not stat.S_ISDIR(inspected.st_mode):
            raise _fail("bundle path must be a regular non-linked directory")
        flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW
        flags |= getattr(os, "O_CLOEXEC", 0)
        descriptor = os.open(source, flags)
        opened = os.fstat(descriptor)
        identity = (opened.st_dev, opened.st_ino)
        if not stat.S_ISDIR(opened.st_mode) or identity != (
            inspected.st_dev,
            inspected.st_ino,
        ):
            raise _fail("bundle directory changed while being opened")
        _closed_names(set(os.listdir(descriptor)))

        snapshots: dict[str, bytes] = {}
        identities: dict[str, os.stat_result] = {}
        file_flags = os.O_RDONLY | os.O_NOFOLLOW
        file_flags |= getattr(os, "O_CLOEXEC", 0)
        file_flags |= getattr(os, "O_NONBLOCK", 0)
        for filename in sorted(BUNDLE_FILENAMES):
            inspected_child = os.stat(
                filename,
                dir_fd=descriptor,
                follow_symlinks=False,
            )
            if stat.S_ISLNK(inspected_child.st_mode) or not stat.S_ISREG(
                inspected_child.st_mode
            ):
                raise _fail(
                    f"bundle artifact {filename} must be a regular non-linked file"
                )
            identities[filename] = inspected_child
            child = os.open(filename, file_flags, dir_fd=descriptor)
            try:
                opened_child = os.fstat(child)
                if (opened_child.st_dev, opened_child.st_ino) != (
                    inspected_child.st_dev,
                    inspected_child.st_ino,
                ):
                    raise _fail(f"bundle artifact {filename} changed while opening")
                snapshots[filename] = _read_descriptor(
                    child,
                    opened_child,
                    f"bundle artifact {filename}",
                    MAX_CHILD_BYTES,
                )
            finally:
                os.close(child)

        _closed_names(set(os.listdir(descriptor)))
        stable = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
        for filename, inspected_child in identities.items():
            attached_child = os.stat(
                filename,
                dir_fd=descriptor,
                follow_symlinks=False,
            )
            if stat.S_ISLNK(attached_child.st_mode) or any(
                getattr(inspected_child, field) != getattr(attached_child, field)
                for field in stable
            ):
                raise _fail(f"bundle artifact {filename} path changed while reading")
        attached = source.lstat()
        finished = os.fstat(descriptor)
        directory_stable = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
        if (
            stat.S_ISLNK(attached.st_mode)
            or not stat.S_ISDIR(attached.st_mode)
            or (attached.st_dev, attached.st_ino) != identity
            or (finished.st_dev, finished.st_ino) != identity
            or any(
                getattr(opened, field) != getattr(finished, field)
                for field in directory_stable
            )
        ):
            raise _fail("bundle directory changed while being read")
        return snapshots
    except INSDCGraphValidationError:
        raise
    except (OSError, TypeError, ValueError) as failure:
        raise _fail(f"cannot read bundle directory: {failure}") from failure
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _json_object(raw: bytes, label: str) -> dict[str, Any]:
    def pairs(values: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in values:
            if key in result:
                raise _fail(f"{label} contains duplicate JSON key {key!r}")
            result[key] = value
        return result

    def constant(token: str) -> None:
        raise _fail(f"{label} contains non-finite number {token}")

    try:
        value = json.loads(
            raw.decode("utf-8"),
            object_pairs_hook=pairs,
            parse_constant=constant,
        )
    except INSDCGraphValidationError:
        raise
    except MemoryError as failure:
        raise _fail(f"invalid {label} JSON: memory ceiling exceeded") from failure
    except (
        UnicodeDecodeError,
        json.JSONDecodeError,
        ValueError,
        RecursionError,
    ) as failure:
        raise _fail(f"invalid {label} JSON: {failure}") from failure
    if type(value) is not dict:
        raise _fail(f"{label} must contain a JSON object")
    return value


@_bounded_validation
def validate_insdc_graph_report(report: dict[str, Any]) -> dict[str, Any]:
    """Validate the exact sealed success-report contract."""

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
        "validation report",
    )
    _validate_tree(root, "validation report")
    if (
        root["format"] != "brainc.bio.insdc-graph-validation-report"
        or type(root["version"]) is not int
        or root["version"] != 1
        or root["stage"] != "genbank-feature-state-adapter"
        or type(root["valid"]) is not bool
        or root["valid"] is not True
    ):
        raise _fail("validation report identity is invalid")
    inputs = _keys(
        root["inputs"],
        {
            "genbank_artifact_sha256",
            "bio_ir_sha256",
            "bundle_artifact_sha256",
            "backend_spec_sha256",
        },
        "validation report.inputs",
    )
    for name, value in inputs.items():
        _sha(value, f"validation report.inputs.{name}")
    result = _keys(
        root["result"],
        {
            "operations",
            "tensor_bytes",
            "units",
            "edges",
            "attachments",
            "semantics_sha256",
            "module_sha256",
        },
        "validation report.result",
    )
    for name in ("operations", "tensor_bytes", "units", "edges", "attachments"):
        _integer(result[name], f"validation report.result.{name}")
    for name in ("semantics_sha256", "module_sha256"):
        _sha(result[name], f"validation report.result.{name}")
    if root["checks"] != list(REPORT_CHECKS):
        raise _fail("validation report checks differ from the closed contract")
    claimed = _sha(root["report_sha256"], "validation report.report_sha256")
    if claimed != digest(
        {key: value for key, value in root.items() if key != "report_sha256"}
    ):
        raise _fail("validation report seal is invalid")
    return root


@_bounded_validation
def _report_bytes(report: dict[str, Any]) -> bytes:
    validate_insdc_graph_report(report)
    try:
        raw = (
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
    if len(raw) > MAX_REPORT_BYTES:
        raise _fail(f"validation report exceeds {MAX_REPORT_BYTES} bytes")
    return raw


def _atomic_write_report(path: str | Path, raw: bytes) -> None:
    destination = Path(path)
    parent_descriptor = -1
    output_descriptor = -1
    temporary_name: str | None = None
    temporary_identity: tuple[int, int] | None = None
    backup_name: str | None = None
    backup_identity: tuple[int, int] | None = None
    existing_identity: tuple[int, int] | None = None
    destination_was_absent: bool | None = None
    succeeded = False
    try:
        if destination.name in {"", ".", ".."}:
            raise _fail("validation report path must name a file")
        if not _HAS_DIRECTORY_DESCRIPTOR:
            raise _fail("secure validation report output is unsupported on this platform")
        destination.parent.mkdir(parents=True, exist_ok=True)
        inspected_parent = destination.parent.lstat()
        if stat.S_ISLNK(inspected_parent.st_mode) or not stat.S_ISDIR(
            inspected_parent.st_mode
        ):
            raise _fail("validation report parent must be a non-linked directory")
        parent_flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW
        parent_flags |= getattr(os, "O_CLOEXEC", 0)
        parent_descriptor = os.open(destination.parent, parent_flags)
        opened_parent = os.fstat(parent_descriptor)
        parent_identity = (opened_parent.st_dev, opened_parent.st_ino)
        if not stat.S_ISDIR(opened_parent.st_mode) or parent_identity != (
            inspected_parent.st_dev,
            inspected_parent.st_ino,
        ):
            raise _fail("validation report parent changed while opening")
        try:
            existing = os.stat(
                destination.name,
                dir_fd=parent_descriptor,
                follow_symlinks=False,
            )
        except FileNotFoundError:
            existing = None
        destination_was_absent = existing is None
        if existing is not None and not stat.S_ISREG(existing.st_mode):
            raise _fail(
                "validation report path must be absent or a regular non-linked file"
            )
        if existing is not None:
            existing_identity = (existing.st_dev, existing.st_ino)

        temporary_name = f".{destination.name}.{secrets.token_hex(12)}"
        output_flags = os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW
        output_flags |= getattr(os, "O_CLOEXEC", 0)
        output_descriptor = os.open(
            temporary_name,
            output_flags,
            0o600,
            dir_fd=parent_descriptor,
        )
        temporary_metadata = os.fstat(output_descriptor)
        temporary_identity = (temporary_metadata.st_dev, temporary_metadata.st_ino)
        view = memoryview(raw)
        while view:
            written = os.write(output_descriptor, view)
            if written <= 0:
                raise OSError("validation report write made no progress")
            view = view[written:]
        os.fsync(output_descriptor)
        finished = os.fstat(output_descriptor)
        if finished.st_size != len(raw):
            raise _fail("validation report output has the wrong size")

        try:
            current = os.stat(
                destination.name,
                dir_fd=parent_descriptor,
                follow_symlinks=False,
            )
        except FileNotFoundError:
            current = None
        if current is not None and not stat.S_ISREG(current.st_mode):
            raise _fail(
                "validation report path must be absent or a regular non-linked file"
            )
        current_identity = (
            None if current is None else (current.st_dev, current.st_ino)
        )
        if current_identity != existing_identity:
            raise _fail("validation report path changed before publishing")
        if existing_identity is not None:
            backup_name = f".{destination.name}.previous.{secrets.token_hex(12)}"
            os.link(
                destination.name,
                backup_name,
                src_dir_fd=parent_descriptor,
                dst_dir_fd=parent_descriptor,
                follow_symlinks=False,
            )
            backup = os.stat(
                backup_name,
                dir_fd=parent_descriptor,
                follow_symlinks=False,
            )
            backup_identity = (backup.st_dev, backup.st_ino)
            if not stat.S_ISREG(backup.st_mode) or backup_identity != existing_identity:
                raise _fail("validation report backup changed while being created")
        attached_temporary = os.stat(
            temporary_name,
            dir_fd=parent_descriptor,
            follow_symlinks=False,
        )
        if (
            not stat.S_ISREG(attached_temporary.st_mode)
            or (attached_temporary.st_dev, attached_temporary.st_ino)
            != temporary_identity
        ):
            raise _fail("validation report temporary changed before publishing")
        os.replace(
            temporary_name,
            destination.name,
            src_dir_fd=parent_descriptor,
            dst_dir_fd=parent_descriptor,
        )
        published = os.stat(
            destination.name,
            dir_fd=parent_descriptor,
            follow_symlinks=False,
        )
        if (
            not stat.S_ISREG(published.st_mode)
            or (published.st_dev, published.st_ino) != temporary_identity
        ):
            raise _fail("validation report path changed while publishing")
        attached_parent = destination.parent.lstat()
        if (
            stat.S_ISLNK(attached_parent.st_mode)
            or not stat.S_ISDIR(attached_parent.st_mode)
            or (attached_parent.st_dev, attached_parent.st_ino) != parent_identity
        ):
            raise _fail("validation report parent changed while publishing")
        try:
            os.fsync(parent_descriptor)
        except OSError as failure:
            unsupported = {
                errno.EBADF,
                errno.EINVAL,
                getattr(errno, "ENOTSUP", errno.EINVAL),
                getattr(errno, "EOPNOTSUPP", errno.EINVAL),
            }
            if failure.errno not in unsupported:
                raise
        final = os.stat(
            destination.name,
            dir_fd=parent_descriptor,
            follow_symlinks=False,
        )
        final_parent = destination.parent.lstat()
        if (
            not stat.S_ISREG(final.st_mode)
            or (final.st_dev, final.st_ino) != temporary_identity
            or stat.S_ISLNK(final_parent.st_mode)
            or not stat.S_ISDIR(final_parent.st_mode)
            or (final_parent.st_dev, final_parent.st_ino) != parent_identity
        ):
            raise _fail("validation report path changed before completion")
        if backup_name is not None:
            backup = os.stat(
                backup_name,
                dir_fd=parent_descriptor,
                follow_symlinks=False,
            )
            if (
                not stat.S_ISREG(backup.st_mode)
                or (backup.st_dev, backup.st_ino) != backup_identity
            ):
                raise _fail("validation report backup changed before cleanup")
        os.lseek(output_descriptor, 0, os.SEEK_SET)
        chunks: list[bytes] = []
        total = 0
        while True:
            chunk = os.read(
                output_descriptor,
                min(64 * 1024, len(raw) + 1 - total),
            )
            if not chunk:
                break
            chunks.append(chunk)
            total += len(chunk)
            if total > len(raw):
                raise _fail("validation report content changed during publication")
        if b"".join(chunks) != raw:
            raise _fail("validation report content changed during publication")
        verified = os.fstat(output_descriptor)
        final = os.stat(
            destination.name,
            dir_fd=parent_descriptor,
            follow_symlinks=False,
        )
        stable = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
        if not stat.S_ISREG(final.st_mode) or any(
            getattr(verified, field) != getattr(final, field) for field in stable
        ):
            raise _fail("validation report path changed during cleanup")
        if backup_name is not None:
            os.unlink(backup_name, dir_fd=parent_descriptor)
            backup_name = None
        succeeded = True
    except INSDCGraphValidationError:
        raise
    except (OSError, TypeError, ValueError) as failure:
        raise _fail(f"cannot write validation report: {failure}") from failure
    finally:
        if not succeeded and parent_descriptor >= 0:
            try:
                try:
                    current = os.stat(
                        destination.name,
                        dir_fd=parent_descriptor,
                        follow_symlinks=False,
                    )
                    current_identity = (current.st_dev, current.st_ino)
                except FileNotFoundError:
                    current_identity = None
                if existing_identity is not None and backup_name is not None:
                    try:
                        backup = os.stat(
                            backup_name,
                            dir_fd=parent_descriptor,
                            follow_symlinks=False,
                        )
                    except FileNotFoundError:
                        backup = None
                    if current_identity == existing_identity:
                        os.unlink(backup_name, dir_fd=parent_descriptor)
                        backup_name = None
                    elif backup is not None and (
                        stat.S_ISREG(backup.st_mode)
                        and (backup.st_dev, backup.st_ino) == backup_identity
                    ):
                        os.replace(
                            backup_name,
                            destination.name,
                            src_dir_fd=parent_descriptor,
                            dst_dir_fd=parent_descriptor,
                        )
                        backup_name = None
                elif destination_was_absent is True and current_identity is not None:
                    try:
                        os.unlink(destination.name, dir_fd=parent_descriptor)
                    except IsADirectoryError:
                        pass
                if temporary_identity is not None:
                    for name in os.listdir(parent_descriptor):
                        try:
                            candidate = os.stat(
                                name,
                                dir_fd=parent_descriptor,
                                follow_symlinks=False,
                            )
                            if (
                                stat.S_ISREG(candidate.st_mode)
                                and (candidate.st_dev, candidate.st_ino)
                                == temporary_identity
                            ):
                                os.unlink(name, dir_fd=parent_descriptor)
                        except OSError:
                            pass
            except OSError:
                pass
        if output_descriptor >= 0:
            os.close(output_descriptor)
        if temporary_name is not None and parent_descriptor >= 0:
            try:
                os.unlink(temporary_name, dir_fd=parent_descriptor)
            except OSError:
                pass
        if parent_descriptor >= 0:
            os.close(parent_descriptor)


@_bounded_validation
def load_insdc_graph_bundle_directory(
    path: str | Path,
) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    """Load the exact closed index-and-nine-children directory convention."""

    raw = _read_bundle_directory(path)
    index = _json_object(raw[BUNDLE_FILENAME], "bundle index")
    artifacts = {
        role: _json_object(raw[filename], f"child artifact {role}")
        for role, filename in CHILD_FILENAMES
    }
    return index, artifacts


@_bounded_validation
def validate_insdc_graph_paths(
    genbank_path: str | Path,
    genbank_artifact_path: str | Path,
    bundle_directory: str | Path,
    *,
    report_path: str | Path | None = None,
) -> dict[str, Any]:
    """Snapshot and validate one closed bundle against its original GenBank."""

    index, artifacts = load_insdc_graph_bundle_directory(bundle_directory)
    source_raw = _read_regular(
        genbank_artifact_path,
        "compiled GenBank artifact",
        MAX_SOURCE_ARTIFACT_BYTES,
    )
    genbank_raw = _read_regular(
        genbank_path,
        "original GenBank source",
        MAX_GENBANK_BYTES,
    )
    source = _json_object(source_raw, "compiled GenBank artifact")
    report = validate_insdc_graph_bundle(
        index,
        artifacts,
        source,
        genbank_source=genbank_raw,
    )
    validate_insdc_graph_report(report)
    if report_path is not None:
        _atomic_write_report(report_path, _report_bytes(report))
    return report


def main(argv: list[str] | None = None) -> int:
    """Run the independent GenBank feature-state validator from paths."""

    parser = argparse.ArgumentParser(
        description="Independently validate a GenBank feature-state bundle"
    )
    parser.add_argument("genbank", help="original GenBank flat-file source")
    parser.add_argument("artifact", help="compiled GenBank v2 artifact JSON")
    parser.add_argument("bundle", help="closed graph bundle directory")
    parser.add_argument(
        "-o",
        "--output",
        "--report",
        dest="output",
        help="atomically write the sealed validation report",
    )
    arguments = parser.parse_args(argv)
    try:
        report = validate_insdc_graph_paths(
            arguments.genbank,
            arguments.artifact,
            arguments.bundle,
            report_path=arguments.output,
        )
    except INSDCGraphValidationError as failure:
        parser.exit(1, f"validation failed: {failure}\n")
    if arguments.output is None:
        sys.stdout.buffer.write(_report_bytes(report))
    return 0


__all__ = [
    "BACKEND_SPEC_SHA256",
    "BUNDLE_FILENAME",
    "CHILD_FILENAMES",
    "CHILD_ROLES",
    "INSDCGraphValidationError",
    "canonical_bytes",
    "digest",
    "load_insdc_graph_bundle_directory",
    "validate_insdc_graph_bundle",
    "validate_insdc_graph_paths",
    "validate_insdc_graph_report",
]


if __name__ == "__main__":
    raise SystemExit(main())

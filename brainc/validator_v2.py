"""Independent standard-library validator for the tensor compiler chain.

This module intentionally imports no :mod:`brainc` implementation code.  It
parses every wire artifact again, implements RFC 8785 canonicalization again,
recomputes all content identities, replays the compiler-owned lowering, and
checks the emitted development module against that replay.  External provider
execution and the biological truth of provider outputs remain external facts.
"""

from __future__ import annotations

import argparse
import base64
import binascii
from dataclasses import dataclass
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import stat
import struct
import tempfile
from typing import Any
import zlib

from .validator_insdc import INSDCValidationError, validate_genbank


SAFE_INTEGER = 2**53 - 1
SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
IUPAC = frozenset("ACGTRYSWKMBDHVN")
DTYPE_BYTES = {"bool": 1, "i64": 8, "u64": 8, "f64": 8}
DEV_DOMAIN = "io.github.dttdrv.epigenesis.dev"
DEV_VERSION = 1
OP_UNIT_CREATE = f"{DEV_DOMAIN}.unit.create"
OP_EDGE_CREATE = f"{DEV_DOMAIN}.edge.create"
OP_RULE_ATTACH = f"{DEV_DOMAIN}.rule.attach"
OP_PORT_BIND = f"{DEV_DOMAIN}.port.bind"
SUPPORTED_OPERATIONS = {
    OP_UNIT_CREATE,
    OP_EDGE_CREATE,
    OP_RULE_ATTACH,
    OP_PORT_BIND,
}
PHASES = ("INPUT", "DELIVER", "UPDATE", "EMIT", "LEARN", "OUTPUT", "SEAL")
SCOPES = ("subject", "source", "target", "event")
SOURCE_ROLES = (
    "sequence",
    "provider_manifest",
    "prediction_request",
    "prediction_response",
    "lowering_policy",
    "target_contract",
)
SEQUENCE_COMPILER = {
    "name": "brainc-dna",
    "version": "0.3.0",
    "passes": [
        "parse-fasta",
        "parse-context",
        "validate-iupac",
        "resolve-reference",
        "canonicalize-sequence",
        "emit-sequence-ir",
    ],
}
COLLECTION_COMPILER = {"name": "brainc-dna-collection", "version": "0.2.0"}
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
GENBANK_SOURCE = ("brainc.bio.insdc-genbank-ir", 2)
SOURCE_IR_FIELDS = {
    ("brain01.sequence-ir", 2): "ir_sha256",
    ("brain01.sequence-collection-ir", 1): "collection_ir_sha256",
    GENBANK_SOURCE: "bio_ir_sha256",
}
_DOMAIN_SEGMENT = r"[a-z0-9](?:[a-z0-9-]*[a-z0-9])?"
DOMAIN_RE = re.compile(rf"(?:{_DOMAIN_SEGMENT}\.)+{_DOMAIN_SEGMENT}\Z")


class ValidationError(ValueError):
    """An independently evaluated compiler invariant failed."""


@dataclass(frozen=True)
class ValidationLimits:
    """Hard ceilings enforced before allocating compiler-described state."""

    input_bytes: int = 16 * 1024 * 1024
    source_artifact_bytes: int = 64 * 1024 * 1024
    decompressed_bytes: int = 64 * 1024 * 1024
    json_depth: int = 64
    json_members: int = 1_000_000
    string_bytes: int = 1 * 1024 * 1024
    identifier_bytes: int = 256
    tensor_rank: int = 16
    tensor_bytes: int = 64 * 1024 * 1024
    total_tensor_bytes: int = 256 * 1024 * 1024
    tensors: int = 100_000
    operations: int = 100_000
    source_records: int = 100_000
    units: int = 10_000_000
    edges: int = 10_000_000
    attachments: int = 100_000
    ports: int = 100_000
    schemas: int = 100_000
    rules: int = 100_000
    fields_per_schema: int = 100_000

    def validate(self) -> None:
        if type(self) is not ValidationLimits:
            raise ValidationError("limits must be an exact ValidationLimits value")
        maxima = {
            "input_bytes": 16 * 1024 * 1024,
            "source_artifact_bytes": 64 * 1024 * 1024,
            "decompressed_bytes": 64 * 1024 * 1024,
            "json_depth": 64,
            "json_members": 1_000_000,
            "string_bytes": 1 * 1024 * 1024,
            "identifier_bytes": 256,
            "tensor_rank": 16,
            "tensor_bytes": 64 * 1024 * 1024,
            "total_tensor_bytes": 256 * 1024 * 1024,
            "tensors": 100_000,
            "operations": 100_000,
            "source_records": 100_000,
            "units": 10_000_000,
            "edges": 10_000_000,
            "attachments": 100_000,
            "ports": 100_000,
            "schemas": 100_000,
            "rules": 100_000,
            "fields_per_schema": 100_000,
        }
        for name, maximum in maxima.items():
            value = getattr(self, name)
            if type(value) is not int or value < 0 or value > maximum:
                raise ValidationError(
                    f"limits.{name} must be an integer in [0, {maximum}]"
                )


DEFAULT_LIMITS = ValidationLimits()


def _reject_constant(value: str) -> None:
    raise ValidationError(f"non-finite JSON number is not allowed: {value}")


def _parse_integer(value: str) -> int:
    parsed = int(value)
    if not -SAFE_INTEGER <= parsed <= SAFE_INTEGER:
        raise ValidationError("JSON integer exceeds the I-JSON safe range")
    return parsed


def _parse_float(value: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValidationError("non-finite JSON number is not allowed")
    return parsed


def _no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValidationError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _validate_tree(value: Any, label: str, limits: ValidationLimits) -> None:
    members = 0
    stack: list[tuple[Any, int]] = [(value, 0)]
    while stack:
        current, depth = stack.pop()
        if depth > limits.json_depth:
            raise ValidationError(f"{label} exceeds JSON depth limit {limits.json_depth}")
        if current is None or type(current) is bool:
            continue
        if type(current) is int:
            if not -SAFE_INTEGER <= current <= SAFE_INTEGER:
                raise ValidationError(f"{label} contains an unsafe JSON integer")
            continue
        if type(current) is float:
            if not math.isfinite(current):
                raise ValidationError(f"{label} contains a non-finite JSON number")
            continue
        if type(current) is str:
            if any(0xD800 <= ord(character) <= 0xDFFF for character in current):
                raise ValidationError(f"{label} contains a lone Unicode surrogate")
            if len(current.encode("utf-8")) > limits.string_bytes:
                raise ValidationError(
                    f"{label} exceeds JSON string byte limit {limits.string_bytes}"
                )
            continue
        if type(current) is list:
            members += len(current)
            stack.extend((item, depth + 1) for item in reversed(current))
        elif type(current) is dict:
            members += len(current)
            for key, item in reversed(list(current.items())):
                if type(key) is not str:
                    raise ValidationError(f"{label} contains a non-string object key")
                stack.append((key, depth + 1))
                stack.append((item, depth + 1))
        else:
            raise ValidationError(
                f"{label} contains unsupported type {type(current).__name__}"
            )
        if members > limits.json_members:
            raise ValidationError(
                f"{label} exceeds JSON member limit {limits.json_members}"
            )


def _loads_bounded(
    raw: bytes,
    label: str,
    limits: ValidationLimits,
    maximum: int,
) -> dict[str, Any]:
    limits.validate()
    if type(raw) is not bytes:
        raise ValidationError(f"{label} JSON input must be bytes")
    if len(raw) > maximum:
        raise ValidationError(f"{label} exceeds JSON byte limit {maximum}")
    try:
        value = json.loads(
            raw.decode("utf-8"),
            object_pairs_hook=_no_duplicates,
            parse_constant=_reject_constant,
            parse_int=_parse_integer,
            parse_float=_parse_float,
        )
    except ValidationError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError, RecursionError) as failure:
        raise ValidationError(f"invalid {label} JSON: {failure}") from failure
    if type(value) is not dict:
        raise ValidationError(f"{label} must be a JSON object")
    _validate_tree(value, label, limits)
    return value


def loads(
    raw: bytes, label: str = "artifact", *, limits: ValidationLimits = DEFAULT_LIMITS
) -> dict[str, Any]:
    """Decode one closed I-JSON object under the generic artifact ceiling."""

    return _loads_bounded(raw, label, limits, limits.input_bytes)


def _read_regular(path: str | Path, label: str, maximum: int) -> bytes:
    source = Path(path)
    descriptor = -1
    try:
        metadata = source.lstat()
        if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISREG(metadata.st_mode):
            raise ValidationError(f"{label} must be a regular non-linked file")
        if metadata.st_size > maximum:
            raise ValidationError(f"{label} exceeds byte limit {maximum}")
        descriptor = os.open(source, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
        opened = os.fstat(descriptor)
        if (
            not stat.S_ISREG(opened.st_mode)
            or (opened.st_dev, opened.st_ino) != (metadata.st_dev, metadata.st_ino)
            or opened.st_size > maximum
        ):
            raise ValidationError(f"{label} changed or is not a safe regular file")
        chunks: list[bytes] = []
        remaining = opened.st_size
        while remaining:
            chunk = os.read(descriptor, min(remaining, 1024 * 1024))
            if not chunk:
                raise ValidationError(f"{label} ended before its inspected length")
            chunks.append(chunk)
            remaining -= len(chunk)
        if os.read(descriptor, 1):
            raise ValidationError(f"{label} grew while being read")
        finished = os.fstat(descriptor)
        stable_fields = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
        if any(
            getattr(opened, field) != getattr(finished, field)
            for field in stable_fields
        ):
            raise ValidationError(f"{label} changed while being read")
        return b"".join(chunks)
    except ValidationError:
        raise
    except OSError as failure:
        raise ValidationError(f"cannot read {label}: {failure}") from failure
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def load(
    path: str | Path,
    label: str = "artifact",
    *,
    limits: ValidationLimits = DEFAULT_LIMITS,
) -> tuple[dict[str, Any], bytes]:
    raw = _read_regular(path, label, limits.input_bytes)
    return loads(raw, label, limits=limits), raw


def _load_source_artifact(
    path: str | Path,
    label: str,
    limits: ValidationLimits,
) -> tuple[dict[str, Any], bytes]:
    raw = _read_regular(path, label, limits.source_artifact_bytes)
    value = _loads_bounded(raw, label, limits, limits.source_artifact_bytes)
    if (
        (value.get("format"), value.get("version")) != GENBANK_SOURCE
        and len(raw) > limits.input_bytes
    ):
        raise ValidationError(
            f"{label} exceeds JSON byte limit {limits.input_bytes}"
        )
    return value, raw


def _jcs_string(value: str) -> str:
    if any(0xD800 <= ord(character) <= 0xDFFF for character in value):
        raise ValueError("lone Unicode surrogate is not allowed")
    return json.dumps(value, ensure_ascii=False, allow_nan=False)


def _jcs_number(value: float) -> str:
    if not math.isfinite(value):
        raise ValueError("non-finite number is not allowed")
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
    rendered = digits[0] + (("." + digits[1:]) if len(digits) > 1 else "")
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
            raise ValueError(f"integer exceeds I-JSON safe range: {value}")
        return str(value)
    if type(value) is float:
        return _jcs_number(value)
    if type(value) is list:
        return "[" + ",".join(_jcs(item) for item in value) + "]"
    if type(value) is dict:
        if any(type(key) is not str for key in value):
            raise TypeError("JSON object keys must be strings")
        ordered = sorted(value, key=lambda key: key.encode("utf-16be"))
        return "{" + ",".join(
            _jcs_string(key) + ":" + _jcs(value[key]) for key in ordered
        ) + "}"
    raise TypeError(f"unsupported JSON value type: {type(value).__name__}")


def canonical_bytes(
    value: Any, *, limits: ValidationLimits = DEFAULT_LIMITS
) -> bytes:
    """Serialize an I-JSON value with an independent RFC 8785 implementation."""

    try:
        _validate_tree(value, "canonical value", limits)
        return _jcs(value).encode("utf-8")
    except ValidationError:
        raise
    except (TypeError, ValueError, UnicodeError) as failure:
        raise ValidationError(f"value is not RFC 8785 canonical JSON: {failure}") from failure


def digest(value: Any, *, limits: ValidationLimits = DEFAULT_LIMITS) -> str:
    return hashlib.sha256(canonical_bytes(value, limits=limits)).hexdigest()


def bytes_digest(value: bytes) -> str:
    if type(value) is not bytes:
        raise ValidationError("byte digest input must be bytes")
    return hashlib.sha256(value).hexdigest()


def _keys(value: Any, expected: set[str], label: str) -> dict[str, Any]:
    if type(value) is not dict:
        raise ValidationError(f"{label} must be an object")
    missing = sorted(expected - value.keys())
    extra = sorted(value.keys() - expected)
    if missing or extra:
        raise ValidationError(
            f"{label} keys invalid; missing={missing or 'none'}, unknown={extra or 'none'}"
        )
    return value


def _integer(
    value: Any,
    label: str,
    *,
    minimum: int = 0,
    maximum: int = SAFE_INTEGER,
) -> int:
    if type(value) is not int or not minimum <= value <= maximum:
        raise ValidationError(
            f"{label} must be an integer in [{minimum}, {maximum}]"
        )
    return value


def _text(value: Any, label: str, *, allow_empty: bool = False) -> str:
    if type(value) is not str or (not allow_empty and not value) or value != value.strip():
        raise ValidationError(
            f"{label} must be a {'possibly empty ' if allow_empty else ''}trimmed string"
        )
    return value


def _identifier(value: Any, label: str, limits: ValidationLimits) -> str:
    parsed = _text(value, label)
    if len(parsed.encode("utf-8")) > limits.identifier_bytes:
        raise ValidationError(f"{label} exceeds identifier byte ceiling")
    return parsed


def _domain(value: Any, label: str, limits: ValidationLimits) -> str:
    parsed = _identifier(value, label, limits)
    if DOMAIN_RE.fullmatch(parsed) is None:
        raise ValidationError(
            f"{label} must be a lowercase dotted namespaced identifier"
        )
    return parsed


def _optional_unit(value: Any, label: str, limits: ValidationLimits) -> str | None:
    if value is None:
        return None
    return _domain(value, label, limits)


def _sha256(value: Any, label: str) -> str:
    if type(value) is not str or SHA256_RE.fullmatch(value) is None:
        raise ValidationError(f"{label} must be a lowercase SHA-256 digest")
    return value


def _artifact(value: dict[str, Any], label: str, limits: ValidationLimits) -> str:
    stored = _sha256(value.get("artifact_sha256"), f"{label}.artifact_sha256")
    expected = digest(
        {key: item for key, item in value.items() if key != "artifact_sha256"},
        limits=limits,
    )
    if stored != expected:
        raise ValidationError(f"{label}.artifact_sha256 does not match canonical content")
    return stored


def _sorted_unique(values: list[str], label: str) -> None:
    if values != sorted(values):
        raise ValidationError(f"{label} must be sorted")
    if len(values) != len(set(values)):
        raise ValidationError(f"{label} contains duplicate ids")


def _provider(value: Any, label: str, limits: ValidationLimits) -> dict[str, Any]:
    parsed = _keys(value, {"name", "version"}, label)
    _identifier(parsed["name"], f"{label}.name", limits)
    _identifier(parsed["version"], f"{label}.version", limits)
    return parsed


def _model_identity(value: Any, label: str, limits: ValidationLimits) -> dict[str, Any]:
    parsed = _keys(value, {"kind", "value"}, label)
    if parsed["kind"] == "content-sha256":
        _sha256(parsed["value"], f"{label}.value")
    elif parsed["kind"] == "opaque":
        _text(parsed["value"], f"{label}.value")
    else:
        raise ValidationError(f"{label}.kind is unsupported")
    return parsed


def _shape(value: Any, label: str, limits: ValidationLimits) -> tuple[int, ...]:
    if type(value) is not list:
        raise ValidationError(f"{label} must be an array")
    if len(value) > limits.tensor_rank:
        raise ValidationError(f"{label} exceeds tensor rank ceiling")
    result: list[int] = []
    product = 1
    for index, dimension in enumerate(value):
        parsed = _integer(dimension, f"{label}[{index}]")
        product *= parsed
        if product > SAFE_INTEGER:
            raise ValidationError(f"{label} product exceeds I-JSON safe range")
        result.append(parsed)
    return tuple(result)


def _tensor_type(value: Any, label: str, limits: ValidationLimits) -> dict[str, Any]:
    parsed = _keys(value, {"dtype", "shape"}, label)
    dtype = _identifier(parsed["dtype"], f"{label}.dtype", limits)
    if dtype not in DTYPE_BYTES:
        raise ValidationError(f"{label}.dtype is unsupported")
    shape = _shape(parsed["shape"], f"{label}.shape", limits)
    elements = math.prod(shape)
    byte_length = elements * DTYPE_BYTES[dtype]
    if byte_length > limits.tensor_bytes:
        raise ValidationError(f"{label} exceeds per-tensor byte ceiling")
    return {"dtype": dtype, "shape": list(shape)}


def _type_bytes(value: dict[str, Any]) -> int:
    return math.prod(value["shape"]) * DTYPE_BYTES[value["dtype"]]


def _axis(
    value: Any,
    label: str,
    dimension: int,
    limits: ValidationLimits,
    *,
    sequence_sha256: str | None,
    records: dict[str, int] | None,
) -> dict[str, Any] | None:
    if value is None:
        return None
    parsed = _keys(
        value,
        {
            "kind",
            "source_artifact_sha256",
            "record_id",
            "start",
            "step",
            "span",
            "coordinate_system",
        },
        label,
    )
    if (
        parsed["kind"] != "source-coordinate"
        or parsed["coordinate_system"] != "0-based-half-open"
    ):
        raise ValidationError(f"{label} has unsupported source-coordinate semantics")
    source = _sha256(parsed["source_artifact_sha256"], f"{label}.source_artifact_sha256")
    if sequence_sha256 is not None and source != sequence_sha256:
        raise ValidationError(f"{label} does not bind the supplied sequence source")
    record_id = _identifier(parsed["record_id"], f"{label}.record_id", limits)
    start = _integer(parsed["start"], f"{label}.start")
    step = _integer(parsed["step"], f"{label}.step", minimum=1)
    span = _integer(parsed["span"], f"{label}.span", minimum=1)
    if records is not None:
        if record_id not in records:
            raise ValidationError(f"{label} references unknown sequence record {record_id!r}")
        if start + span > records[record_id]:
            raise ValidationError(f"{label} starts outside its sequence record")
    if dimension and start + (dimension - 1) * step + span > SAFE_INTEGER:
        raise ValidationError(f"{label} endpoint exceeds I-JSON safe range")
    if (
        records is not None
        and dimension
        and start + (dimension - 1) * step + span > records[record_id]
    ):
        raise ValidationError(f"{label} coordinate exceeds its sequence record")
    return parsed


def _axes(
    value: Any,
    tensor_type: dict[str, Any],
    label: str,
    limits: ValidationLimits,
    *,
    sequence_sha256: str | None,
    records: dict[str, int] | None = None,
) -> list[dict[str, Any] | None]:
    shape = tensor_type["shape"]
    if type(value) is not list or len(value) != len(shape):
        raise ValidationError(f"{label} must have one entry per tensor dimension")
    return [
        _axis(
            item,
            f"{label}[{index}]",
            shape[index],
            limits,
            sequence_sha256=sequence_sha256,
            records=records,
        )
        for index, item in enumerate(value)
    ]


def _read_blob(
    root: str | Path | None,
    blob_sha256: str,
    expected_length: int,
    label: str,
    limits: ValidationLimits,
) -> bytes:
    if root is None:
        raise ValidationError(f"{label} requires an explicit blob root")
    if expected_length > limits.tensor_bytes:
        raise ValidationError(f"{label} exceeds the tensor byte ceiling")
    if os.name != "posix":
        raise ValidationError(
            f"{label} external blob storage requires a POSIX host; use inline-base64"
        )
    root_path = Path(root)
    root_descriptor = -1
    descriptor = -1
    try:
        root_metadata = root_path.lstat()
        if stat.S_ISLNK(root_metadata.st_mode) or not stat.S_ISDIR(root_metadata.st_mode):
            raise ValidationError("blob root must be a non-linked directory")
        root_descriptor = os.open(
            root_path,
            os.O_RDONLY
            | getattr(os, "O_DIRECTORY", 0)
            | getattr(os, "O_NOFOLLOW", 0),
        )
        opened_root = os.fstat(root_descriptor)
        if (
            not stat.S_ISDIR(opened_root.st_mode)
            or (opened_root.st_dev, opened_root.st_ino)
            != (root_metadata.st_dev, root_metadata.st_ino)
        ):
            raise ValidationError("blob root changed or is not a safe directory")
        before = os.stat(blob_sha256, dir_fd=root_descriptor, follow_symlinks=False)
        if (
            stat.S_ISLNK(before.st_mode)
            or not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
        ):
            raise ValidationError(f"{label} must be one regular non-linked file")
        if before.st_size != expected_length:
            raise ValidationError(f"{label} length does not match storage descriptor")
        descriptor = os.open(
            blob_sha256,
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=root_descriptor,
        )
        after = os.fstat(descriptor)
        if (
            not stat.S_ISREG(after.st_mode)
            or after.st_nlink != 1
            or (after.st_dev, after.st_ino) != (before.st_dev, before.st_ino)
            or after.st_size != expected_length
        ):
            raise ValidationError(f"{label} changed or is not a safe regular file")
        chunks: list[bytes] = []
        remaining = expected_length
        while remaining:
            chunk = os.read(descriptor, min(remaining, 1024 * 1024))
            if not chunk:
                raise ValidationError(f"{label} ended before its declared length")
            chunks.append(chunk)
            remaining -= len(chunk)
        if os.read(descriptor, 1):
            raise ValidationError(f"{label} exceeds its declared length")
        raw = b"".join(chunks)
    except ValidationError:
        raise
    except OSError as failure:
        raise ValidationError(f"cannot read {label}: {failure}") from failure
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        if root_descriptor >= 0:
            os.close(root_descriptor)
    if bytes_digest(raw) != blob_sha256:
        raise ValidationError(f"{label} digest mismatch")
    return raw


def _validate_tensor_elements(raw: bytes, dtype: str, label: str) -> None:
    if dtype == "bool":
        if any(value not in (0, 1) for value in raw):
            raise ValidationError(f"{label} contains an invalid bool byte")
        return
    if dtype != "f64":
        return
    for offset in range(0, len(raw), 8):
        bits = int.from_bytes(raw[offset : offset + 8], "little")
        if bits & 0x7FF0000000000000 == 0x7FF0000000000000:
            raise ValidationError(f"{label} contains non-finite f64")
        if bits == 0x8000000000000000:
            raise ValidationError(f"{label} contains negative-zero f64")


def _storage(
    value: Any,
    expected_length: int,
    label: str,
    limits: ValidationLimits,
    *,
    blob_root: str | Path | None,
) -> tuple[dict[str, Any], bytes]:
    if type(value) is not dict or "kind" not in value:
        raise ValidationError(f"{label} must be a storage descriptor")
    if expected_length > limits.tensor_bytes:
        raise ValidationError(f"{label} exceeds per-tensor byte ceiling")
    kind = value["kind"]
    if kind == "inline-base64":
        parsed = _keys(value, {"kind", "data", "byte_length", "sha256"}, label)
        data = parsed["data"]
        if type(data) is not str or any(character.isspace() for character in data):
            raise ValidationError(f"{label}.data must be strict padded base64")
        encoded_length = ((expected_length + 2) // 3) * 4
        if len(data) != encoded_length:
            raise ValidationError(f"{label}.data length does not match tensor type")
        if _integer(parsed["byte_length"], f"{label}.byte_length") != expected_length:
            raise ValidationError(f"{label}.byte_length does not match tensor type")
        try:
            raw = base64.b64decode(data, validate=True)
        except (ValueError, binascii.Error) as failure:
            raise ValidationError(f"{label}.data is invalid base64") from failure
        if base64.b64encode(raw).decode("ascii") != data:
            raise ValidationError(f"{label}.data is not canonical padded base64")
        stored = _sha256(parsed["sha256"], f"{label}.sha256")
        if len(raw) != expected_length or bytes_digest(raw) != stored:
            raise ValidationError(f"{label} bytes do not match descriptor")
    elif kind == "sha256-blob":
        parsed = _keys(value, {"kind", "byte_length", "sha256"}, label)
        if _integer(parsed["byte_length"], f"{label}.byte_length") != expected_length:
            raise ValidationError(f"{label}.byte_length does not match tensor type")
        stored = _sha256(parsed["sha256"], f"{label}.sha256")
        raw = _read_blob(blob_root, stored, expected_length, label, limits)
    else:
        raise ValidationError(f"{label}.kind is unsupported")
    return parsed, raw


def _output_contract(
    value: Any,
    label: str,
    limits: ValidationLimits,
    *,
    sequence_sha256: str | None,
) -> dict[str, Any]:
    parsed = _keys(value, {"id", "type", "unit", "axes"}, label)
    _identifier(parsed["id"], f"{label}.id", limits)
    tensor_type = _tensor_type(parsed["type"], f"{label}.type", limits)
    _optional_unit(parsed["unit"], f"{label}.unit", limits)
    _axes(
        parsed["axes"],
        tensor_type,
        f"{label}.axes",
        limits,
        sequence_sha256=sequence_sha256,
    )
    return parsed


def _tensor_payload(
    value: Any,
    label: str,
    limits: ValidationLimits,
    *,
    sequence_sha256: str,
    records: dict[str, int] | None,
    blob_root: str | Path | None,
    lineage: bool,
) -> dict[str, Any]:
    expected = {"id", "type", "unit", "axes", "storage"}
    if lineage:
        expected.add("lineage")
    parsed = _keys(value, expected, label)
    tensor_id = _identifier(parsed["id"], f"{label}.id", limits)
    tensor_type = _tensor_type(parsed["type"], f"{label}.type", limits)
    unit = _optional_unit(parsed["unit"], f"{label}.unit", limits)
    axes = _axes(
        parsed["axes"],
        tensor_type,
        f"{label}.axes",
        limits,
        sequence_sha256=sequence_sha256,
        records=records,
    )
    storage, raw = _storage(
        parsed["storage"],
        _type_bytes(tensor_type),
        f"{label}.storage",
        limits,
        blob_root=blob_root,
    )
    _validate_tensor_elements(raw, tensor_type["dtype"], label)
    result = {
        "id": tensor_id,
        "type": tensor_type,
        "unit": unit,
        "axes": axes,
        "storage": storage,
        "raw": raw,
    }
    if lineage:
        parsed_lineage = _keys(
            parsed["lineage"], {"response_output", "lowering_value"}, f"{label}.lineage"
        )
        response_output = _identifier(
            parsed_lineage["response_output"], f"{label}.lineage.response_output", limits
        )
        lowering_value = _identifier(
            parsed_lineage["lowering_value"], f"{label}.lineage.lowering_value", limits
        )
        if lowering_value != tensor_id:
            raise ValidationError(f"{label}.lineage.lowering_value must equal tensor id")
        result["response_output"] = response_output
        result["lowering_value"] = lowering_value
    return result


def _source_binding(value: dict[str, Any]) -> dict[str, Any]:
    identity = (value["format"], value["version"])
    try:
        ir_field = SOURCE_IR_FIELDS[identity]
    except KeyError as failure:
        raise ValidationError("sequence source format/version is unsupported") from failure
    return {
        "format": value["format"],
        "version": value["version"],
        "artifact_sha256": value["artifact_sha256"],
        "ir_sha256": value[ir_field],
    }


def _source_reference(value: Any, label: str) -> tuple[dict[str, Any], str]:
    source = _keys(
        value,
        {"format", "version", "artifact_sha256", "ir_sha256"},
        label,
    )
    if (
        type(source["format"]) is not str
        or type(source["version"]) is not int
        or (source["format"], source["version"]) not in SOURCE_IR_FIELDS
    ):
        raise ValidationError(f"{label} format/version is unsupported")
    source_sha = _sha256(source["artifact_sha256"], f"{label}.artifact_sha256")
    _sha256(source["ir_sha256"], f"{label}.ir_sha256")
    return source, source_sha


def validate_manifest_artifact(
    value: dict[str, Any], *, limits: ValidationLimits = DEFAULT_LIMITS
) -> dict[str, Any]:
    top = _keys(
        value,
        {
            "format",
            "version",
            "provider",
            "model_identity",
            "accepts",
            "outputs",
            "artifact_sha256",
        },
        "provider manifest",
    )
    if top["format"] != "brainc.provider-manifest" or _integer(
        top["version"], "provider manifest.version", minimum=1
    ) != 2:
        raise ValidationError("unsupported provider manifest format/version")
    _provider(top["provider"], "provider manifest.provider", limits)
    _model_identity(top["model_identity"], "provider manifest.model_identity", limits)
    accepts = top["accepts"]
    if type(accepts) is not list or not accepts:
        raise ValidationError("provider manifest.accepts must be a nonempty array")
    normalized_accepts = [
        _identifier(item, f"provider manifest.accepts[{index}]", limits)
        for index, item in enumerate(accepts)
    ]
    if normalized_accepts != sorted(normalized_accepts) or len(normalized_accepts) != len(
        set(normalized_accepts)
    ):
        raise ValidationError("provider manifest.accepts must be sorted and unique")
    outputs = top["outputs"]
    if type(outputs) is not list or not outputs:
        raise ValidationError("provider manifest.outputs must be a nonempty array")
    if len(outputs) > limits.tensors:
        raise ValidationError("provider manifest.outputs exceeds tensor-count ceiling")
    parsed_outputs = [
        _output_contract(
            item,
            f"provider manifest.outputs[{index}]",
            limits,
            sequence_sha256=None,
        )
        for index, item in enumerate(outputs)
    ]
    _sorted_unique([item["id"] for item in parsed_outputs], "provider manifest.outputs")
    _artifact(top, "provider manifest", limits)
    return top


def validate_request_artifact(
    value: dict[str, Any], *, limits: ValidationLimits = DEFAULT_LIMITS
) -> dict[str, Any]:
    top = _keys(
        value,
        {
            "format",
            "version",
            "source",
            "provider_manifest_sha256",
            "requested_outputs",
            "artifact_sha256",
        },
        "prediction request",
    )
    if top["format"] != "brainc.prediction-request" or _integer(
        top["version"], "prediction request.version", minimum=1
    ) != 2:
        raise ValidationError("unsupported prediction request format/version")
    _, source_sha = _source_reference(top["source"], "prediction request.source")
    _sha256(top["provider_manifest_sha256"], "prediction request provider manifest")
    outputs = top["requested_outputs"]
    if type(outputs) is not list or not outputs:
        raise ValidationError("prediction request.requested_outputs must be nonempty")
    parsed_outputs = [
        _output_contract(
            item,
            f"prediction request.requested_outputs[{index}]",
            limits,
            sequence_sha256=source_sha,
        )
        for index, item in enumerate(outputs)
    ]
    _sorted_unique(
        [item["id"] for item in parsed_outputs], "prediction request.requested_outputs"
    )
    _artifact(top, "prediction request", limits)
    return top


def validate_response_artifact(
    value: dict[str, Any],
    *,
    sequence_sha256: str,
    records: dict[str, int] | None = None,
    blob_root: str | Path | None = None,
    limits: ValidationLimits = DEFAULT_LIMITS,
) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    top = _keys(
        value,
        {
            "format",
            "version",
            "request_artifact_sha256",
            "provider",
            "model_identity",
            "outputs",
            "artifact_sha256",
        },
        "prediction response",
    )
    if top["format"] != "brainc.prediction-response" or _integer(
        top["version"], "prediction response.version", minimum=1
    ) != 2:
        raise ValidationError("unsupported prediction response format/version")
    _sha256(top["request_artifact_sha256"], "prediction response request")
    _provider(top["provider"], "prediction response.provider", limits)
    _model_identity(top["model_identity"], "prediction response.model_identity", limits)
    outputs = top["outputs"]
    if type(outputs) is not list or not outputs:
        raise ValidationError("prediction response.outputs must be nonempty")
    if len(outputs) > limits.tensors:
        raise ValidationError("prediction response.outputs exceeds tensor-count ceiling")
    parsed = [
        _tensor_payload(
            item,
            f"prediction response.outputs[{index}]",
            limits,
            sequence_sha256=sequence_sha256,
            records=records,
            blob_root=blob_root,
            lineage=False,
        )
        for index, item in enumerate(outputs)
    ]
    _sorted_unique([item["id"] for item in parsed], "prediction response.outputs")
    total = sum(len(item["raw"]) for item in parsed)
    if total > limits.total_tensor_bytes:
        raise ValidationError("prediction response exceeds total tensor byte ceiling")
    _artifact(top, "prediction response", limits)
    return top, {item["id"]: item for item in parsed}


def _parse_fasta(raw: bytes) -> tuple[str, str, str, list[dict[str, int]]]:
    try:
        lines = raw.decode("utf-8").splitlines()
    except UnicodeDecodeError as failure:
        raise ValidationError("FASTA must be UTF-8") from failure
    if not lines or not lines[0].startswith(">"):
        raise ValidationError("FASTA must begin with a defline")
    head = lines[0][1:].split(maxsplit=1)
    if not head:
        raise ValidationError("FASTA record id is empty")
    record_id, description = head[0], head[1] if len(head) == 2 else ""
    parts: list[str] = []
    segments: list[dict[str, int]] = []
    offset = 0
    for line_number, line in enumerate(lines[1:], 2):
        if line.startswith(">"):
            raise ValidationError("single-record FASTA contains another defline")
        if not line:
            continue
        if any(character.isspace() or character.upper() not in IUPAC for character in line):
            raise ValidationError(f"invalid IUPAC sequence at FASTA line {line_number}")
        parts.append(line)
        segments.append(
            {"line": line_number, "normalized_start": offset, "length": len(line)}
        )
        offset += len(line)
    sequence = "".join(parts)
    if not sequence:
        raise ValidationError("FASTA sequence is empty")
    return record_id, description, sequence, segments


def _sequence_context(
    raw: bytes,
    record_id: str,
    sequence_length: int,
    limits: ValidationLimits,
) -> tuple[Any, Any]:
    value = loads(raw, "sequence context", limits=limits)
    _keys(
        value,
        {"format", "version", "record_id", "reference", "provenance"},
        "sequence context",
    )
    if (
        value["format"] != "brain01.sequence-context"
        or type(value["version"]) is not int
        or value["version"] != 1
        or value["record_id"] != record_id
    ):
        raise ValidationError("sequence context identity mismatch")
    reference = value["reference"]
    if reference is not None:
        reference = _keys(
            reference,
            {
                "assembly",
                "contig",
                "start",
                "end",
                "coordinate_system",
                "orientation",
                "aliases",
            },
            "sequence context.reference",
        )
        _text(reference["assembly"], "sequence context.reference.assembly")
        _text(reference["contig"], "sequence context.reference.contig")
        start = _integer(reference["start"], "sequence context.reference.start")
        end = _integer(reference["end"], "sequence context.reference.end")
        if end - start != sequence_length:
            raise ValidationError("sequence context reference span mismatch")
        if (
            reference["coordinate_system"] != "0-based-half-open"
            or reference["orientation"] not in ("forward", "reverse")
        ):
            raise ValidationError("sequence context coordinate convention mismatch")
        aliases = reference["aliases"]
        if type(aliases) is not list:
            raise ValidationError("sequence context aliases must be an array")
        normalized = [
            _text(item, f"sequence context.reference.aliases[{index}]")
            for index, item in enumerate(aliases)
        ]
        if len(normalized) != len(set(normalized)):
            raise ValidationError("sequence context aliases are duplicated")
    provenance = value["provenance"]
    if type(provenance) is not list:
        raise ValidationError("sequence context provenance must be an array")
    seen: set[str] = set()
    for index, raw_source in enumerate(provenance):
        if type(raw_source) is not dict:
            raise ValidationError(f"sequence context.provenance[{index}] must be an object")
        required = {"id", "kind", "uri", "version"}
        allowed = required | {"sha256"}
        missing = required - raw_source.keys()
        extra = raw_source.keys() - allowed
        if missing or extra:
            raise ValidationError(
                f"sequence context.provenance[{index}] keys are invalid"
            )
        source_id = _text(raw_source["id"], f"sequence context.provenance[{index}].id")
        if source_id in seen:
            raise ValidationError("sequence context provenance ids are duplicated")
        seen.add(source_id)
        for field in ("kind", "uri", "version"):
            _text(raw_source[field], f"sequence context.provenance[{index}].{field}")
        if "sha256" in raw_source:
            _sha256(raw_source["sha256"], f"sequence context.provenance[{index}].sha256")
    return reference, provenance


def _sequence_payload(
    raw: bytes,
    *,
    context_raw: bytes | None,
    limits: ValidationLimits,
) -> dict[str, Any]:
    record_id, description, sequence, segments = _parse_fasta(raw)
    reference: Any = None
    provenance: Any = []
    context_sha: str | None = None
    if context_raw is not None:
        context_sha = bytes_digest(context_raw)
        reference, provenance = _sequence_context(
            context_raw, record_id, len(sequence), limits
        )
    canonical = sequence.upper()
    sequence_ir = {
        "record_id": record_id,
        "description": description,
        "sequence": sequence,
        "sequence_sha256": bytes_digest(sequence.encode("ascii")),
        "canonical_sha256": bytes_digest(canonical.encode("ascii")),
        "refget_id": "SQ."
        + base64.urlsafe_b64encode(
            hashlib.sha512(canonical.encode("ascii")).digest()[:24]
        )
        .decode("ascii")
        .rstrip("="),
        "reference": reference,
        "provenance": provenance,
    }
    core = {
        "format": "brain01.sequence-ir",
        "version": 2,
        "compiler": SEQUENCE_COMPILER,
        "inputs": {
            "fasta_sha256": bytes_digest(raw),
            "context_sha256": context_sha,
        },
        "sequence_ir": sequence_ir,
        "ir_sha256": digest(sequence_ir, limits=limits),
        "source_map": {"sequence_segments": segments},
    }
    return {**core, "artifact_sha256": digest(core, limits=limits)}


def _unwrap_source(raw: bytes, limits: ValidationLimits) -> tuple[bytes, str | None]:
    if not raw.startswith(b"\x1f\x8b"):
        if len(raw) > limits.decompressed_bytes:
            raise ValidationError("logical source exceeds decompressed byte ceiling")
        return raw, None
    try:
        with gzip.GzipFile(fileobj=io.BytesIO(raw), mode="rb") as stream:
            logical = stream.read(limits.decompressed_bytes + 1)
    except (gzip.BadGzipFile, EOFError, OSError, zlib.error) as failure:
        raise ValidationError(f"malformed gzip source: {failure}") from failure
    if len(logical) > limits.decompressed_bytes:
        raise ValidationError("gzip source exceeds decompressed byte ceiling")
    return logical, bytes_digest(raw)


def _sha512t24u(value: bytes) -> str:
    return base64.urlsafe_b64encode(hashlib.sha512(value).digest()[:24]).decode("ascii")


def _collection_payload(
    raw: bytes,
    *,
    record_id: str | None,
    limits: ValidationLimits,
) -> dict[str, Any]:
    logical, compressed_sha = _unwrap_source(raw, limits)
    if logical.startswith(b">"):
        if record_id is not None:
            raise ValidationError("record_id is only valid for raw IUPAC collection input")
        kind = "fasta"
        lines = logical.splitlines(keepends=True)
        if not lines or not lines[0].startswith(b">"):
            raise ValidationError("FASTA collection must begin with a defline")
        starts = [index for index, line in enumerate(lines) if line.startswith(b">")]
        if len(starts) > limits.source_records:
            raise ValidationError("FASTA collection exceeds source record ceiling")
        records = [
            (
                start,
                b"".join(
                    lines[
                        start : (
                            starts[position + 1]
                            if position + 1 < len(starts)
                            else len(lines)
                        )
                    ]
                ),
            )
            for position, start in enumerate(starts)
        ]
    else:
        if record_id is None:
            raise ValidationError("record_id is required for raw IUPAC collection input")
        _text(record_id, "raw record_id")
        if any(character.isspace() for character in record_id):
            raise ValidationError("raw record_id must not contain whitespace")
        try:
            sequence = logical.decode("ascii")
        except UnicodeDecodeError as failure:
            raise ValidationError("raw sequence must be ASCII") from failure
        if not sequence or any(
            character.isspace() or character.upper() not in IUPAC for character in sequence
        ):
            raise ValidationError(
                "raw sequence must be nonempty whitespace-free IUPAC DNA"
            )
        kind = "raw-iupac"
        records = [(None, f">{record_id}\n{sequence}\n".encode("utf-8"))]
    members: list[dict[str, Any]] = []
    seen: set[str] = set()
    for start, record_raw in records:
        sequence_artifact = _sequence_payload(
            record_raw, context_raw=None, limits=limits
        )
        member_id = sequence_artifact["sequence_ir"]["record_id"]
        if member_id in seen:
            raise ValidationError(f"duplicate FASTA record identifier {member_id!r}")
        seen.add(member_id)
        local_segments = sequence_artifact["source_map"]["sequence_segments"]
        outer_segments = (
            [
                {
                    "line": 1,
                    "normalized_start": 0,
                    "length": len(sequence_artifact["sequence_ir"]["sequence"]),
                }
            ]
            if start is None
            else [{**segment, "line": segment["line"] + start} for segment in local_segments]
        )
        members.append(
            {
                "record_id": member_id,
                "ir_sha256": sequence_artifact["ir_sha256"],
                "artifact_sha256": sequence_artifact["artifact_sha256"],
                "input_source_map": {"sequence_segments": outer_segments},
                "sequence_artifact": sequence_artifact,
            }
        )
    collection_ir = {
        "members": [
            {"record_id": member["record_id"], "ir_sha256": member["ir_sha256"]}
            for member in members
        ]
    }
    level_2 = {
        "lengths": [
            len(member["sequence_artifact"]["sequence_ir"]["sequence"])
            for member in members
        ],
        "names": [member["record_id"] for member in members],
        "sequences": [
            member["sequence_artifact"]["sequence_ir"]["refget_id"]
            for member in members
        ],
    }
    level_1 = {
        name: _sha512t24u(canonical_bytes(value, limits=limits))
        for name, value in level_2.items()
    }
    refget_seqcol = {
        "version": "1.0.0",
        "digest": _sha512t24u(
            canonical_bytes(
                {name: level_1[name] for name in ("names", "sequences")},
                limits=limits,
            )
        ),
        "level_1": level_1,
        "level_2": level_2,
    }
    core = {
        "format": "brain01.sequence-collection-ir",
        "version": 1,
        "compiler": COLLECTION_COMPILER,
        "inputs": {
            "kind": kind,
            "wrapper": "gzip" if compressed_sha is not None else None,
            "raw_sha256": bytes_digest(raw),
            "logical_sha256": bytes_digest(logical),
            "compressed_sha256": compressed_sha,
        },
        "collection_ir": collection_ir,
        "collection_ir_sha256": digest(collection_ir, limits=limits),
        "refget_seqcol": refget_seqcol,
        "members": members,
    }
    return {**core, "artifact_sha256": digest(core, limits=limits)}


def _sequence_segments(
    value: Any,
    sequence_length: int,
    label: str,
    limits: ValidationLimits,
    *,
    first_line: int,
) -> list[tuple[int, int, int]]:
    source_map = _keys(value, {"sequence_segments"}, label)
    raw_segments = source_map["sequence_segments"]
    if type(raw_segments) is not list or not raw_segments:
        raise ValidationError(f"{label}.sequence_segments must be a nonempty array")
    if len(raw_segments) > limits.json_members:
        raise ValidationError(f"{label}.sequence_segments exceeds JSON member ceiling")
    segments: list[tuple[int, int, int]] = []
    expected_start = 0
    previous_line = first_line - 1
    for index, raw_segment in enumerate(raw_segments):
        segment_label = f"{label}.sequence_segments[{index}]"
        segment = _keys(
            raw_segment, {"line", "normalized_start", "length"}, segment_label
        )
        line = _integer(
            segment["line"],
            f"{segment_label}.line",
            minimum=first_line,
            maximum=limits.json_members + 1,
        )
        start = _integer(segment["normalized_start"], f"{segment_label}.normalized_start")
        length = _integer(segment["length"], f"{segment_label}.length", minimum=1)
        if line <= previous_line or start != expected_start:
            raise ValidationError(f"{segment_label} is not a continuous ordered mapping")
        expected_start += length
        if expected_start > sequence_length:
            raise ValidationError(f"{label} covers more bases than its sequence")
        previous_line = line
        segments.append((line, start, length))
    if expected_start != sequence_length:
        raise ValidationError(f"{label} does not cover the complete sequence")
    return segments


def _sequence_reference(
    value: Any, sequence_length: int, label: str, limits: ValidationLimits
) -> None:
    if value is None:
        return
    reference = _keys(
        value,
        {
            "assembly",
            "contig",
            "start",
            "end",
            "coordinate_system",
            "orientation",
            "aliases",
        },
        label,
    )
    _text(reference["assembly"], f"{label}.assembly")
    _text(reference["contig"], f"{label}.contig")
    start = _integer(reference["start"], f"{label}.start")
    end = _integer(reference["end"], f"{label}.end", minimum=1)
    if end <= start or end - start != sequence_length:
        raise ValidationError(f"{label} interval must exactly span the sequence")
    if reference["coordinate_system"] != "0-based-half-open":
        raise ValidationError(f"{label}.coordinate_system is unsupported")
    if reference["orientation"] not in {"forward", "reverse"}:
        raise ValidationError(f"{label}.orientation is unsupported")
    aliases = reference["aliases"]
    if type(aliases) is not list:
        raise ValidationError(f"{label}.aliases must be an array")
    if len(aliases) > limits.json_members:
        raise ValidationError(f"{label}.aliases exceeds JSON member ceiling")
    normalized = [
        _text(alias, f"{label}.aliases[{index}]")
        for index, alias in enumerate(aliases)
    ]
    if len(normalized) != len(set(normalized)):
        raise ValidationError(f"{label}.aliases contains duplicates")


def _sequence_provenance(
    value: Any, label: str, limits: ValidationLimits
) -> None:
    if type(value) is not list:
        raise ValidationError(f"{label} must be an array")
    if len(value) > limits.json_members:
        raise ValidationError(f"{label} exceeds JSON member ceiling")
    seen: set[str] = set()
    for index, raw_source in enumerate(value):
        source_label = f"{label}[{index}]"
        if type(raw_source) is not dict:
            raise ValidationError(f"{source_label} must be an object")
        required = {"id", "kind", "uri", "version"}
        optional = {"sha256"}
        missing = sorted(required - raw_source.keys())
        extra = sorted(raw_source.keys() - required - optional)
        if missing or extra:
            raise ValidationError(
                f"{source_label} keys invalid; missing={missing or 'none'}, "
                f"unknown={extra or 'none'}"
            )
        source_id = _identifier(raw_source["id"], f"{source_label}.id", limits)
        if source_id in seen:
            raise ValidationError(f"{label} contains duplicate ids")
        seen.add(source_id)
        for field in ("kind", "uri", "version"):
            _text(raw_source[field], f"{source_label}.{field}")
        if "sha256" in raw_source:
            _sha256(raw_source["sha256"], f"{source_label}.sha256")


def _validate_sequence_source(
    value: Any, label: str, limits: ValidationLimits
) -> tuple[str, int, list[tuple[int, int, int]]]:
    source = _keys(
        value,
        {
            "format",
            "version",
            "compiler",
            "inputs",
            "sequence_ir",
            "ir_sha256",
            "source_map",
            "artifact_sha256",
        },
        label,
    )
    if source["format"] != "brain01.sequence-ir" or type(source["version"]) is not int or source["version"] != 2:
        raise ValidationError(f"{label} format/version is unsupported")
    compiler = _keys(source["compiler"], {"name", "version", "passes"}, f"{label}.compiler")
    if compiler != SEQUENCE_COMPILER:
        raise ValidationError(f"{label}.compiler identity or pass set is unsupported")
    inputs = _keys(
        source["inputs"], {"fasta_sha256", "context_sha256"}, f"{label}.inputs"
    )
    _sha256(inputs["fasta_sha256"], f"{label}.inputs.fasta_sha256")
    context_sha = inputs["context_sha256"]
    if context_sha is not None:
        _sha256(context_sha, f"{label}.inputs.context_sha256")

    ir = _keys(
        source["sequence_ir"],
        {
            "record_id",
            "description",
            "sequence",
            "sequence_sha256",
            "canonical_sha256",
            "refget_id",
            "reference",
            "provenance",
        },
        f"{label}.sequence_ir",
    )
    record_id = _identifier(ir["record_id"], f"{label}.sequence_ir.record_id", limits)
    if any(character.isspace() for character in record_id):
        raise ValidationError(f"{label}.sequence_ir.record_id must be whitespace-free")
    description = ir["description"]
    if type(description) is not str or (description and description.splitlines() != [description]):
        raise ValidationError(f"{label}.sequence_ir.description must be one text line")
    sequence = ir["sequence"]
    if type(sequence) is not str or not sequence or any(
        ord(base) > 127 or base.upper() not in IUPAC for base in sequence
    ):
        raise ValidationError(f"{label}.sequence_ir.sequence must be nonempty ASCII IUPAC DNA")
    if len(sequence) > limits.string_bytes:
        raise ValidationError(f"{label}.sequence_ir.sequence exceeds JSON string byte ceiling")
    if len(sequence) > limits.decompressed_bytes:
        raise ValidationError(f"{label}.sequence_ir.sequence exceeds logical DNA byte ceiling")
    sequence_bytes = sequence.encode("ascii")
    canonical_bytes_value = sequence.upper().encode("ascii")
    if _sha256(ir["sequence_sha256"], f"{label}.sequence_ir.sequence_sha256") != bytes_digest(sequence_bytes):
        raise ValidationError(f"{label}.sequence_ir.sequence_sha256 does not match sequence")
    if _sha256(ir["canonical_sha256"], f"{label}.sequence_ir.canonical_sha256") != bytes_digest(canonical_bytes_value):
        raise ValidationError(f"{label}.sequence_ir.canonical_sha256 does not match canonical sequence")
    expected_refget = "SQ." + base64.urlsafe_b64encode(
        hashlib.sha512(canonical_bytes_value).digest()[:24]
    ).decode("ascii").rstrip("=")
    if ir["refget_id"] != expected_refget:
        raise ValidationError(f"{label}.sequence_ir.refget_id does not match canonical sequence")
    _sequence_reference(ir["reference"], len(sequence), f"{label}.sequence_ir.reference", limits)
    _sequence_provenance(ir["provenance"], f"{label}.sequence_ir.provenance", limits)
    if context_sha is None and (ir["reference"] is not None or ir["provenance"]):
        raise ValidationError(f"{label} has context assertions without a context input")
    segments = _sequence_segments(
        source["source_map"],
        len(sequence),
        f"{label}.source_map",
        limits,
        first_line=2,
    )
    stored_ir = _sha256(source["ir_sha256"], f"{label}.ir_sha256")
    if stored_ir != digest(ir, limits=limits):
        raise ValidationError(f"{label}.ir_sha256 does not match Sequence IR")
    _artifact(source, label, limits)
    return record_id, len(sequence), segments


def _validate_collection_refget(
    value: Any,
    members: list[tuple[str, int, str]],
    label: str,
    limits: ValidationLimits,
) -> None:
    refget = _keys(value, {"version", "digest", "level_1", "level_2"}, label)
    if refget["version"] != "1.0.0":
        raise ValidationError(f"{label}.version is unsupported")
    level_1 = _keys(
        refget["level_1"], {"lengths", "names", "sequences"}, f"{label}.level_1"
    )
    for name in ("lengths", "names", "sequences"):
        _text(level_1[name], f"{label}.level_1.{name}")
    level_2 = _keys(
        refget["level_2"], {"lengths", "names", "sequences"}, f"{label}.level_2"
    )
    for name in ("lengths", "names", "sequences"):
        if type(level_2[name]) is not list:
            raise ValidationError(f"{label}.level_2.{name} must be an array")
        if len(level_2[name]) != len(members):
            raise ValidationError(f"{label}.level_2.{name} length does not match members")
    for index, item in enumerate(level_2["lengths"]):
        _integer(item, f"{label}.level_2.lengths[{index}]", minimum=1)
    for index, item in enumerate(level_2["names"]):
        _identifier(item, f"{label}.level_2.names[{index}]", limits)
    for index, item in enumerate(level_2["sequences"]):
        _text(item, f"{label}.level_2.sequences[{index}]")

    expected_level_2 = {
        "lengths": [length for _, length, _ in members],
        "names": [record_id for record_id, _, _ in members],
        "sequences": [refget_id for _, _, refget_id in members],
    }
    expected_level_1 = {
        name: _sha512t24u(canonical_bytes(items, limits=limits))
        for name, items in expected_level_2.items()
    }
    expected = {
        "version": "1.0.0",
        "digest": _sha512t24u(
            canonical_bytes(
                {name: expected_level_1[name] for name in ("names", "sequences")},
                limits=limits,
            )
        ),
        "level_1": expected_level_1,
        "level_2": expected_level_2,
    }
    if canonical_bytes(refget, limits=limits) != canonical_bytes(expected, limits=limits):
        raise ValidationError(f"{label} does not match collection members")


def _validate_collection_source(
    value: Any, label: str, limits: ValidationLimits
) -> dict[str, int]:
    source = _keys(
        value,
        {
            "format",
            "version",
            "compiler",
            "inputs",
            "collection_ir",
            "collection_ir_sha256",
            "refget_seqcol",
            "members",
            "artifact_sha256",
        },
        label,
    )
    if source["format"] != "brain01.sequence-collection-ir" or type(source["version"]) is not int or source["version"] != 1:
        raise ValidationError(f"{label} format/version is unsupported")
    compiler = _keys(source["compiler"], {"name", "version"}, f"{label}.compiler")
    if compiler != COLLECTION_COMPILER:
        raise ValidationError(f"{label}.compiler identity is unsupported")
    inputs = _keys(
        source["inputs"],
        {"kind", "wrapper", "raw_sha256", "logical_sha256", "compressed_sha256"},
        f"{label}.inputs",
    )
    if inputs["kind"] not in {"fasta", "raw-iupac"}:
        raise ValidationError(f"{label}.inputs.kind is unsupported")
    raw_sha = _sha256(inputs["raw_sha256"], f"{label}.inputs.raw_sha256")
    logical_sha = _sha256(inputs["logical_sha256"], f"{label}.inputs.logical_sha256")
    compressed_sha = inputs["compressed_sha256"]
    if inputs["wrapper"] is None:
        if compressed_sha is not None or raw_sha != logical_sha:
            raise ValidationError(f"{label}.inputs unwrapped digests are inconsistent")
    elif inputs["wrapper"] == "gzip":
        if compressed_sha is None or _sha256(compressed_sha, f"{label}.inputs.compressed_sha256") != raw_sha:
            raise ValidationError(f"{label}.inputs gzip digests are inconsistent")
    else:
        raise ValidationError(f"{label}.inputs.wrapper is unsupported")

    raw_members = source["members"]
    if type(raw_members) is not list or not raw_members:
        raise ValidationError(f"{label}.members must be a nonempty array")
    if len(raw_members) > limits.source_records:
        raise ValidationError(f"{label}.members exceeds source record ceiling")
    records: dict[str, int] = {}
    total_bases = 0
    member_identities: list[tuple[str, int, str]] = []
    expected_collection_members: list[dict[str, str]] = []
    previous_outer_line = 0
    for index, raw_member in enumerate(raw_members):
        member_label = f"{label}.members[{index}]"
        member = _keys(
            raw_member,
            {
                "record_id",
                "ir_sha256",
                "artifact_sha256",
                "input_source_map",
                "sequence_artifact",
            },
            member_label,
        )
        record_id, length, local_segments = _validate_sequence_source(
            member["sequence_artifact"], f"{member_label}.sequence_artifact", limits
        )
        nested = member["sequence_artifact"]
        nested_ir = nested["sequence_ir"]
        if nested["inputs"]["context_sha256"] is not None or nested_ir["reference"] is not None or nested_ir["provenance"]:
            raise ValidationError(f"{member_label} contains context not emitted by collection ingress")
        if member["record_id"] != record_id:
            raise ValidationError(f"{member_label}.record_id does not match nested Sequence IR")
        if _sha256(member["ir_sha256"], f"{member_label}.ir_sha256") != nested["ir_sha256"]:
            raise ValidationError(f"{member_label}.ir_sha256 does not match nested Sequence IR")
        if _sha256(member["artifact_sha256"], f"{member_label}.artifact_sha256") != nested["artifact_sha256"]:
            raise ValidationError(f"{member_label}.artifact_sha256 does not match nested Sequence IR")
        if record_id in records:
            raise ValidationError(f"{label}.members contains duplicate record ids")
        records[record_id] = length
        total_bases += length
        if total_bases > limits.decompressed_bytes:
            raise ValidationError(f"{label}.members exceeds logical DNA byte ceiling")
        outer_segments = _sequence_segments(
            member["input_source_map"],
            length,
            f"{member_label}.input_source_map",
            limits,
            first_line=1,
        )
        if inputs["kind"] == "raw-iupac":
            if len(raw_members) != 1 or outer_segments != [(1, 0, length)]:
                raise ValidationError(f"{label} raw IUPAC mapping must cover one line and one member")
            if nested_ir["description"] != "":
                raise ValidationError(f"{member_label} raw IUPAC member cannot have a description")
            expected_fasta = f">{record_id}\n{nested_ir['sequence']}\n".encode("utf-8")
            if nested["inputs"]["fasta_sha256"] != bytes_digest(expected_fasta):
                raise ValidationError(f"{member_label} synthetic FASTA digest is inconsistent")
            if logical_sha != bytes_digest(nested_ir["sequence"].encode("ascii")):
                raise ValidationError(f"{label}.inputs.logical_sha256 does not match raw IUPAC DNA")
        else:
            if len(outer_segments) != len(local_segments) or any(
                outer[1:] != local[1:]
                for outer, local in zip(outer_segments, local_segments)
            ):
                raise ValidationError(f"{member_label} local/global source maps disagree")
            shifts = {
                outer[0] - local[0]
                for outer, local in zip(outer_segments, local_segments)
            }
            if len(shifts) != 1:
                raise ValidationError(f"{member_label} global source-map shift is inconsistent")
            shift = next(iter(shifts))
            if shift < 0 or (index == 0 and shift != 0):
                raise ValidationError(f"{member_label} global source-map shift is invalid")
            if index and shift + 1 <= previous_outer_line:
                raise ValidationError(f"{member_label} overlaps the preceding FASTA record")
            previous_outer_line = outer_segments[-1][0]
        refget_id = nested_ir["refget_id"]
        member_identities.append((record_id, length, refget_id))
        expected_collection_members.append(
            {"record_id": record_id, "ir_sha256": nested["ir_sha256"]}
        )

    collection_ir = _keys(source["collection_ir"], {"members"}, f"{label}.collection_ir")
    collection_members = collection_ir["members"]
    if type(collection_members) is not list or len(collection_members) != len(raw_members):
        raise ValidationError(f"{label}.collection_ir.members does not match members")
    for index, item in enumerate(collection_members):
        _keys(item, {"record_id", "ir_sha256"}, f"{label}.collection_ir.members[{index}]")
    if canonical_bytes(collection_members, limits=limits) != canonical_bytes(expected_collection_members, limits=limits):
        raise ValidationError(f"{label}.collection_ir.members does not bind nested Sequence IR")
    stored_ir = _sha256(source["collection_ir_sha256"], f"{label}.collection_ir_sha256")
    if stored_ir != digest(collection_ir, limits=limits):
        raise ValidationError(f"{label}.collection_ir_sha256 does not match collection IR")
    _validate_collection_refget(
        source["refget_seqcol"], member_identities, f"{label}.refget_seqcol", limits
    )
    _artifact(source, label, limits)
    return records


def _embedded_genbank_source(value: dict[str, Any], limits: ValidationLimits) -> bytes:
    try:
        storage = value["sequence_collection"]["inputs"]["source"]["storage"]
        chunks = storage["chunks"]
    except (KeyError, TypeError) as failure:
        raise ValidationError(
            "GenBank source artifact does not contain bounded source storage"
        ) from failure
    if type(chunks) is not list or not chunks:
        raise ValidationError("GenBank source storage chunks must be a nonempty array")
    raw = bytearray()
    for index, chunk in enumerate(chunks):
        if type(chunk) is not str:
            raise ValidationError(f"GenBank source storage chunk {index} must be text")
        try:
            encoded = chunk.encode("ascii")
        except UnicodeEncodeError as failure:
            raise ValidationError(
                f"GenBank source storage chunk {index} must be ASCII"
            ) from failure
        if len(raw) + len(encoded) > limits.input_bytes:
            raise ValidationError("raw GenBank source exceeds input byte ceiling")
        raw.extend(encoded)
    if not raw:
        raise ValidationError("raw GenBank source must not be empty")
    return bytes(raw)


def _validated_genbank_records(
    value: dict[str, Any], limits: ValidationLimits
) -> dict[str, int]:
    members = value["sequence_collection"]["members"]
    if len(members) > limits.source_records:
        raise ValidationError("sequence source exceeds source record ceiling")
    records: dict[str, int] = {}
    for index, member in enumerate(members):
        record_id = _identifier(
            member["record_id"],
            f"sequence source.sequence_collection.members[{index}].record_id",
            limits,
        )
        if record_id in records:
            raise ValidationError("sequence source contains duplicate record ids")
        records[record_id] = _integer(
            member["sequence"]["bases"],
            f"sequence source.sequence_collection.members[{index}].sequence.bases",
            minimum=1,
        )
    if not records:
        raise ValidationError("sequence source must contain at least one record")
    return records


def _validate_source_structure(
    value: Any,
    limits: ValidationLimits,
    *,
    genbank_source: bytes | None = None,
) -> dict[str, int]:
    limits.validate()
    _validate_tree(value, "sequence source", limits)
    if type(value) is not dict:
        raise ValidationError("sequence source must be an object")
    identity = (value.get("format"), value.get("version"))
    if identity == ("brain01.sequence-ir", 2):
        if limits.source_records < 1:
            raise ValidationError("sequence source exceeds source record ceiling")
        record_id, length, _ = _validate_sequence_source(value, "sequence source", limits)
        return {record_id: length}
    if identity == ("brain01.sequence-collection-ir", 1):
        return _validate_collection_source(value, "sequence source", limits)
    if identity == GENBANK_SOURCE:
        raw = (
            _embedded_genbank_source(value, limits)
            if genbank_source is None
            else genbank_source
        )
        try:
            validate_genbank(value, genbank_source=raw)
        except INSDCValidationError as failure:
            raise ValidationError(f"invalid GenBank source artifact: {failure}") from failure
        return _validated_genbank_records(value, limits)
    raise ValidationError("sequence source format/version is unsupported")


def _validate_source_and_records(
    value: dict[str, Any],
    *,
    raw_source: bytes | None,
    context_raw: bytes | None,
    record_id: str | None,
    limits: ValidationLimits,
) -> dict[str, int]:
    identity = (
        (value.get("format"), value.get("version"))
        if type(value) is dict
        else (None, None)
    )
    if raw_source is not None:
        if type(raw_source) is not bytes:
            raise ValidationError("raw sequence source must be bytes")
        if len(raw_source) > limits.input_bytes:
            raise ValidationError("raw sequence source exceeds input byte ceiling")
    if context_raw is not None and type(context_raw) is not bytes:
        raise ValidationError("sequence context must be bytes")
    if context_raw is not None and len(context_raw) > limits.input_bytes:
        raise ValidationError("sequence context exceeds input byte ceiling")
    if identity == GENBANK_SOURCE:
        if context_raw is not None or record_id is not None:
            raise ValidationError("GenBank source replay does not accept context/record_id")
        return _validate_source_structure(
            value,
            limits,
            genbank_source=raw_source,
        )

    records = _validate_source_structure(value, limits)
    if raw_source is not None:
        if identity == ("brain01.sequence-ir", 2) and record_id is not None:
            raise ValidationError("record_id is only valid for collection source replay")
        if identity == ("brain01.sequence-collection-ir", 1) and context_raw is not None:
            raise ValidationError("sequence collection replay does not accept a context file")
        expected = (
            _sequence_payload(raw_source, context_raw=context_raw, limits=limits)
            if identity == ("brain01.sequence-ir", 2)
            else _collection_payload(raw_source, record_id=record_id, limits=limits)
        )
        if canonical_bytes(value, limits=limits) != canonical_bytes(expected, limits=limits):
            raise ValidationError("sequence source does not replay from supplied DNA")
    elif context_raw is not None or record_id is not None:
        raise ValidationError("context/record_id require raw sequence source bytes")
    return records


def validate_source_artifact(
    value: dict[str, Any],
    *,
    raw_source: bytes | None = None,
    context_raw: bytes | None = None,
    record_id: str | None = None,
    limits: ValidationLimits = DEFAULT_LIMITS,
) -> dict[str, Any]:
    """Validate a closed source artifact and optionally replay its source bytes."""

    _validate_source_and_records(
        value,
        raw_source=raw_source,
        context_raw=context_raw,
        record_id=record_id,
        limits=limits,
    )
    return value


def source_record_lengths(
    value: dict[str, Any], *, limits: ValidationLimits = DEFAULT_LIMITS
) -> dict[str, int]:
    """Validate a source artifact and return a fresh exact record-length map."""

    return dict(_validate_source_structure(value, limits))


def _numeric_contract(value: Any, dtype: str, label: str) -> dict[str, Any]:
    if type(value) is not dict or "kind" not in value:
        raise ValidationError(f"{label} must be a numeric contract")
    kind = value["kind"]
    if kind in {"boolean", "integer", "float"}:
        parsed = _keys(value, {"kind"}, label)
        if kind == "boolean" and dtype != "bool":
            raise ValidationError(f"{label} boolean semantics require bool")
        if kind == "integer" and dtype not in {"i64", "u64"}:
            raise ValidationError(f"{label} integer semantics require i64 or u64")
        if kind == "float" and dtype != "f64":
            raise ValidationError(f"{label} float semantics require f64")
        return parsed
    if kind == "binary-fixed":
        parsed = _keys(
            value, {"kind", "fraction_bits", "rounding", "overflow"}, label
        )
        if dtype != "i64":
            raise ValidationError(f"{label} binary-fixed semantics require i64")
        _integer(parsed["fraction_bits"], f"{label}.fraction_bits", maximum=62)
        if parsed["rounding"] != "nearest-even" or parsed["overflow"] != "saturate":
            raise ValidationError(f"{label} has unsupported binary-fixed semantics")
        return parsed
    raise ValidationError(f"{label}.kind is unsupported")


def _field_contract(
    value: Any, label: str, limits: ValidationLimits
) -> dict[str, Any]:
    parsed = _keys(value, {"id", "type", "unit", "mutability", "numeric"}, label)
    field_id = _identifier(parsed["id"], f"{label}.id", limits)
    tensor_type = _tensor_type(parsed["type"], f"{label}.type", limits)
    unit = _optional_unit(parsed["unit"], f"{label}.unit", limits)
    if parsed["mutability"] not in {"constant", "state"}:
        raise ValidationError(f"{label}.mutability must be constant or state")
    numeric = _numeric_contract(parsed["numeric"], tensor_type["dtype"], f"{label}.numeric")
    return {
        "id": field_id,
        "type": tensor_type,
        "unit": unit,
        "mutability": parsed["mutability"],
        "numeric": numeric,
    }


def _schema_contract(
    value: Any, label: str, limits: ValidationLimits
) -> dict[str, Any]:
    parsed = _keys(value, {"id", "fields"}, label)
    schema_id = _identifier(parsed["id"], f"{label}.id", limits)
    fields = parsed["fields"]
    if type(fields) is not list or not fields:
        raise ValidationError(f"{label}.fields must be a nonempty array")
    if len(fields) > limits.fields_per_schema:
        raise ValidationError(f"{label}.fields exceeds field ceiling")
    normalized = [
        _field_contract(item, f"{label}.fields[{index}]", limits)
        for index, item in enumerate(fields)
    ]
    _sorted_unique([item["id"] for item in normalized], f"{label}.fields")
    return {
        "id": schema_id,
        "fields": {item["id"]: item for item in normalized},
        "field_order": [item["id"] for item in normalized],
    }


def _port_contract(
    value: Any, label: str, limits: ValidationLimits
) -> dict[str, Any]:
    parsed = _keys(value, {"id", "direction", "type", "unit", "codec"}, label)
    port_id = _identifier(parsed["id"], f"{label}.id", limits)
    if parsed["direction"] not in {"input", "output"}:
        raise ValidationError(f"{label}.direction must be input or output")
    tensor_type = _tensor_type(parsed["type"], f"{label}.type", limits)
    unit = _optional_unit(parsed["unit"], f"{label}.unit", limits)
    codec = _keys(parsed["codec"], {"id", "version", "semantics_sha256"}, f"{label}.codec")
    _domain(codec["id"], f"{label}.codec.id", limits)
    _integer(codec["version"], f"{label}.codec.version", minimum=1)
    _sha256(codec["semantics_sha256"], f"{label}.codec.semantics_sha256")
    return {
        "id": port_id,
        "direction": parsed["direction"],
        "type": tensor_type,
        "unit": unit,
        "codec": codec,
    }


def _access_contract(
    value: Any, label: str, limits: ValidationLimits
) -> dict[str, str]:
    parsed = _keys(value, {"scope", "field"}, label)
    if parsed["scope"] not in SCOPES:
        raise ValidationError(f"{label}.scope is unsupported")
    return {
        "scope": parsed["scope"],
        "field": _identifier(parsed["field"], f"{label}.field", limits),
    }


def _parameter_contract(
    value: Any, label: str, limits: ValidationLimits
) -> dict[str, Any]:
    parsed = _keys(value, {"id", "type", "unit"}, label)
    return {
        "id": _identifier(parsed["id"], f"{label}.id", limits),
        "type": _tensor_type(parsed["type"], f"{label}.type", limits),
        "unit": _optional_unit(parsed["unit"], f"{label}.unit", limits),
    }


def _rule_contract(
    value: Any, label: str, limits: ValidationLimits
) -> dict[str, Any]:
    parsed = _keys(
        value,
        {
            "id",
            "version",
            "semantics_sha256",
            "subject",
            "schema",
            "phase",
            "triggers",
            "parameters",
            "reads",
            "writes",
        },
        label,
    )
    rule_id = _domain(parsed["id"], f"{label}.id", limits)
    version = _integer(parsed["version"], f"{label}.version", minimum=1)
    semantics = _sha256(parsed["semantics_sha256"], f"{label}.semantics_sha256")
    if parsed["subject"] not in {"unit", "edge"}:
        raise ValidationError(f"{label}.subject must be unit or edge")
    schema = _identifier(parsed["schema"], f"{label}.schema", limits)
    if parsed["phase"] not in PHASES:
        raise ValidationError(f"{label}.phase is unsupported")
    triggers = parsed["triggers"]
    if type(triggers) is not list or not triggers:
        raise ValidationError(f"{label}.triggers must be a nonempty array")
    normalized_triggers = [
        _identifier(item, f"{label}.triggers[{index}]", limits)
        for index, item in enumerate(triggers)
    ]
    if normalized_triggers != sorted(normalized_triggers) or len(
        normalized_triggers
    ) != len(set(normalized_triggers)):
        raise ValidationError(f"{label}.triggers must be sorted and unique")
    if type(parsed["parameters"]) is not list:
        raise ValidationError(f"{label}.parameters must be an array")
    parameters = [
        _parameter_contract(item, f"{label}.parameters[{index}]", limits)
        for index, item in enumerate(parsed["parameters"])
    ]
    _sorted_unique([item["id"] for item in parameters], f"{label}.parameters")
    accesses: dict[str, list[dict[str, str]]] = {}
    for name in ("reads", "writes"):
        raw_accesses = parsed[name]
        if type(raw_accesses) is not list:
            raise ValidationError(f"{label}.{name} must be an array")
        normalized = [
            _access_contract(item, f"{label}.{name}[{index}]", limits)
            for index, item in enumerate(raw_accesses)
        ]
        tuples = [(item["scope"], item["field"]) for item in normalized]
        if tuples != sorted(tuples) or len(tuples) != len(set(tuples)):
            raise ValidationError(f"{label}.{name} must be sorted and unique")
        accesses[name] = normalized
    if any(item["scope"] != "subject" for item in accesses["writes"]):
        raise ValidationError(f"{label}.writes must be subject-scoped")
    if parsed["subject"] == "unit" and any(
        item["scope"] in {"source", "target"} for item in accesses["reads"]
    ):
        raise ValidationError(f"{label} unit rule cannot read source/target scope")
    return {
        "id": rule_id,
        "version": version,
        "semantics_sha256": semantics,
        "subject": parsed["subject"],
        "schema": schema,
        "phase": parsed["phase"],
        "triggers": normalized_triggers,
        "parameters": parameters,
        "reads": accesses["reads"],
        "writes": accesses["writes"],
    }


def _require_schema_field(
    schema: dict[str, Any], field_id: str, label: str
) -> dict[str, Any]:
    try:
        return schema["fields"][field_id]
    except KeyError as failure:
        raise ValidationError(f"{label} references unknown field {field_id!r}") from failure


def validate_target_artifact(
    value: dict[str, Any], *, limits: ValidationLimits = DEFAULT_LIMITS
) -> dict[str, Any]:
    """Validate a generic target ABI v1 contract and return independent metadata."""

    top = _keys(
        value,
        {"format", "version", "contract", "contract_sha256", "artifact_sha256"},
        "target contract",
    )
    if top["format"] != "brainc.target-contract" or _integer(
        top["version"], "target contract.version", minimum=1
    ) != 1:
        raise ValidationError("unsupported target contract format/version")
    body = _keys(
        top["contract"],
        {"id", "abi_major", "opsets", "unit_schemas", "edge_schemas", "ports", "rules"},
        "target contract body",
    )
    stored_contract = _sha256(top["contract_sha256"], "target contract.contract_sha256")
    if stored_contract != digest(body, limits=limits):
        raise ValidationError("target contract.contract_sha256 does not match contract")
    artifact_sha = _artifact(top, "target contract", limits)
    target_id = _domain(body["id"], "target contract.id", limits)
    abi_major = _integer(body["abi_major"], "target contract.abi_major", minimum=1)
    if abi_major != 1:
        raise ValidationError("validator supports target ABI major 1 only")
    opsets_raw = body["opsets"]
    if type(opsets_raw) is not list or not opsets_raw:
        raise ValidationError("target contract.opsets must be a nonempty array")
    if len(opsets_raw) > limits.operations:
        raise ValidationError("target contract.opsets exceeds capability ceiling")
    opsets: list[tuple[str, int]] = []
    for index, item in enumerate(opsets_raw):
        parsed = _keys(item, {"domain", "version"}, f"target contract.opsets[{index}]")
        opsets.append(
            (
                _domain(parsed["domain"], f"target contract.opsets[{index}].domain", limits),
                _integer(
                    parsed["version"],
                    f"target contract.opsets[{index}].version",
                    minimum=1,
                ),
            )
        )
    if opsets != sorted(opsets) or len({item[0] for item in opsets}) != len(opsets):
        raise ValidationError("target contract.opsets must be sorted with unique domains")
    collections: dict[str, dict[str, dict[str, Any]]] = {}
    for name, parser, allow_empty, ceiling in (
        ("unit_schemas", _schema_contract, False, limits.schemas),
        ("edge_schemas", _schema_contract, True, limits.schemas),
        ("ports", _port_contract, True, limits.ports),
        ("rules", _rule_contract, True, limits.rules),
    ):
        raw = body[name]
        if type(raw) is not list or (not allow_empty and not raw):
            qualifier = "nonempty " if not allow_empty else ""
            raise ValidationError(f"target contract.{name} must be a {qualifier}array")
        if len(raw) > ceiling:
            raise ValidationError(f"target contract.{name} exceeds host ceiling")
        normalized = [
            parser(item, f"target contract.{name}[{index}]", limits)
            for index, item in enumerate(raw)
        ]
        _sorted_unique([item["id"] for item in normalized], f"target contract.{name}")
        collections[name] = {item["id"]: item for item in normalized}
    unit_fields = {
        field_id
        for schema in collections["unit_schemas"].values()
        for field_id in schema["fields"]
    }
    for rule in collections["rules"].values():
        schemas = (
            collections["unit_schemas"]
            if rule["subject"] == "unit"
            else collections["edge_schemas"]
        )
        if rule["schema"] not in schemas:
            raise ValidationError(
                f"target rule {rule['id']!r} references unknown subject schema"
            )
        subject_schema = schemas[rule["schema"]]
        for access in rule["reads"] + rule["writes"]:
            if access["scope"] == "event":
                continue
            if access["scope"] == "subject":
                field = _require_schema_field(
                    subject_schema,
                    access["field"],
                    f"target rule {rule['id']!r}",
                )
                if access in rule["writes"] and field["mutability"] != "state":
                    raise ValidationError(
                        f"target rule {rule['id']!r} writes a constant field"
                    )
            elif access["field"] not in unit_fields:
                raise ValidationError(
                    f"target rule {rule['id']!r} {access['scope']} field does not resolve"
                )
    return {
        "artifact": top,
        "id": target_id,
        "abi_major": abi_major,
        "contract_sha256": stored_contract,
        "artifact_sha256": artifact_sha,
        "opsets": tuple(opsets),
        **collections,
    }


def _binding_list(
    value: Any,
    label: str,
    limits: ValidationLimits,
    *,
    first: str,
    second: str,
) -> list[dict[str, str]]:
    if type(value) is not list:
        raise ValidationError(f"{label} must be an array")
    result: list[dict[str, str]] = []
    for index, item in enumerate(value):
        parsed = _keys(item, {first, second}, f"{label}[{index}]")
        result.append(
            {
                first: _identifier(parsed[first], f"{label}[{index}].{first}", limits),
                second: _identifier(
                    parsed[second], f"{label}[{index}].{second}", limits
                ),
            }
        )
    _sorted_unique([item[first] for item in result], label)
    return result


def _operation_shape(
    value: Any, label: str, limits: ValidationLimits
) -> dict[str, Any]:
    if type(value) is not dict or "op" not in value:
        raise ValidationError(f"{label} must be an operation object")
    op = value["op"]
    if op not in SUPPORTED_OPERATIONS:
        raise ValidationError(f"{label}.op is unknown or unsupported")
    if op == OP_UNIT_CREATE:
        parsed = _keys(
            value, {"id", "op", "version", "schema", "count", "initializers"}, label
        )
        _identifier(parsed["schema"], f"{label}.schema", limits)
        _identifier(parsed["count"], f"{label}.count", limits)
        _binding_list(
            parsed["initializers"],
            f"{label}.initializers",
            limits,
            first="field",
            second="tensor",
        )
    elif op == OP_EDGE_CREATE:
        parsed = _keys(
            value,
            {
                "id",
                "op",
                "version",
                "schema",
                "sources",
                "targets",
                "pairs",
                "initializers",
            },
            label,
        )
        for name in ("schema", "sources", "targets", "pairs"):
            _identifier(parsed[name], f"{label}.{name}", limits)
        _binding_list(
            parsed["initializers"],
            f"{label}.initializers",
            limits,
            first="field",
            second="tensor",
        )
    elif op == OP_RULE_ATTACH:
        parsed = _keys(
            value,
            {"id", "op", "version", "rule", "rule_version", "subject", "parameters"},
            label,
        )
        _domain(parsed["rule"], f"{label}.rule", limits)
        _integer(parsed["rule_version"], f"{label}.rule_version", minimum=1)
        _identifier(parsed["subject"], f"{label}.subject", limits)
        _binding_list(
            parsed["parameters"],
            f"{label}.parameters",
            limits,
            first="id",
            second="tensor",
        )
    else:
        parsed = _keys(
            value,
            {"id", "op", "version", "port", "subject", "field", "indices"},
            label,
        )
        for name in ("port", "subject", "field", "indices"):
            _identifier(parsed[name], f"{label}.{name}", limits)
    _identifier(parsed["id"], f"{label}.id", limits)
    if _integer(parsed["version"], f"{label}.version", minimum=1) != 1:
        raise ValidationError(f"{label}.version is unsupported")
    return parsed


def _policy_info(
    value: dict[str, Any], limits: ValidationLimits
) -> dict[str, Any]:
    top = _keys(
        value,
        {
            "format",
            "version",
            "id",
            "target",
            "values",
            "operations",
            "artifact_sha256",
        },
        "lowering policy",
    )
    if top["format"] != "brainc.lowering-policy" or _integer(
        top["version"], "lowering policy.version", minimum=1
    ) != 2:
        raise ValidationError("unsupported lowering policy format/version")
    policy_id = _identifier(top["id"], "lowering policy.id", limits)
    target = _keys(
        top["target"], {"id", "abi_major", "contract_sha256"}, "lowering policy.target"
    )
    _domain(target["id"], "lowering policy.target.id", limits)
    if _integer(target["abi_major"], "lowering policy.target.abi_major", minimum=1) != 1:
        raise ValidationError("lowering policy target ABI must be 1")
    _sha256(target["contract_sha256"], "lowering policy.target.contract_sha256")
    values = top["values"]
    if type(values) is not list or not values:
        raise ValidationError("lowering policy.values must be a nonempty array")
    if len(values) > limits.tensors:
        raise ValidationError("lowering policy.values exceeds tensor-count ceiling")
    parsed_values: list[dict[str, str]] = []
    for index, item in enumerate(values):
        parsed = _keys(item, {"id", "from_output"}, f"lowering policy.values[{index}]")
        parsed_values.append(
            {
                "id": _identifier(
                    parsed["id"], f"lowering policy.values[{index}].id", limits
                ),
                "from_output": _identifier(
                    parsed["from_output"],
                    f"lowering policy.values[{index}].from_output",
                    limits,
                ),
            }
        )
    _sorted_unique([item["id"] for item in parsed_values], "lowering policy.values")
    if len({item["from_output"] for item in parsed_values}) != len(parsed_values):
        raise ValidationError("lowering policy.values maps one output more than once")
    operations = top["operations"]
    if type(operations) is not list:
        raise ValidationError("lowering policy.operations must be an array")
    if len(operations) > limits.operations:
        raise ValidationError("lowering policy.operations exceeds operation ceiling")
    parsed_operations = [
        _operation_shape(item, f"lowering policy.operations[{index}]", limits)
        for index, item in enumerate(operations)
    ]
    operation_ids = [item["id"] for item in parsed_operations]
    if len(operation_ids) != len(set(operation_ids)):
        raise ValidationError("lowering policy operation ids are duplicated")
    artifact_sha = _artifact(top, "lowering policy", limits)
    return {
        "artifact": top,
        "id": policy_id,
        "target": target,
        "values": parsed_values,
        "operations": parsed_operations,
        "artifact_sha256": artifact_sha,
    }


def validate_policy_artifact(
    value: dict[str, Any], *, limits: ValidationLimits = DEFAULT_LIMITS
) -> dict[str, Any]:
    _policy_info(value, limits)
    return value


def _require_tensor(
    tensors: dict[str, dict[str, Any]], tensor_id: str, label: str
) -> dict[str, Any]:
    try:
        return tensors[tensor_id]
    except KeyError as failure:
        raise ValidationError(f"{label} references unknown tensor {tensor_id!r}") from failure


def _u64_scalar(tensor: dict[str, Any], label: str) -> int:
    if tensor["type"] != {"dtype": "u64", "shape": []} or tensor["unit"] is not None:
        raise ValidationError(f"{label} must reference a unitless u64 scalar")
    value = struct.unpack("<Q", tensor["raw"])[0]
    if value > SAFE_INTEGER:
        raise ValidationError(f"{label} exceeds ABI v1 safe-integer range")
    return value


def _u64_vector(tensor: dict[str, Any], label: str) -> tuple[int, ...]:
    if (
        tensor["type"]["dtype"] != "u64"
        or len(tensor["type"]["shape"]) != 1
        or tensor["unit"] is not None
    ):
        raise ValidationError(f"{label} must reference a unitless u64 vector")
    values = tuple(item[0] for item in struct.iter_unpack("<Q", tensor["raw"]))
    if any(value > SAFE_INTEGER for value in values):
        raise ValidationError(f"{label} contains an ABI-unsafe integer")
    return values


def _u64_pairs(tensor: dict[str, Any], label: str) -> tuple[tuple[int, int], ...]:
    if (
        tensor["type"]["dtype"] != "u64"
        or len(tensor["type"]["shape"]) != 2
        or tensor["type"]["shape"][1] != 2
        or tensor["unit"] is not None
    ):
        raise ValidationError(f"{label} must reference a unitless u64[E,2] tensor")
    flat = tuple(item[0] for item in struct.iter_unpack("<Q", tensor["raw"]))
    if any(value > SAFE_INTEGER for value in flat):
        raise ValidationError(f"{label} contains an ABI-unsafe integer")
    return tuple((flat[index], flat[index + 1]) for index in range(0, len(flat), 2))


def _initializer_bindings(
    value: Any,
    label: str,
    limits: ValidationLimits,
) -> list[dict[str, str]]:
    return _binding_list(
        value, label, limits, first="field", second="tensor"
    )


def _validate_initializers(
    raw: Any,
    schema: dict[str, Any],
    count: int,
    tensors: dict[str, dict[str, Any]],
    label: str,
    limits: ValidationLimits,
) -> list[dict[str, str]]:
    bindings = _initializer_bindings(raw, label, limits)
    if [item["field"] for item in bindings] != schema["field_order"]:
        raise ValidationError(f"{label} must initialize every schema field exactly once")
    for index, binding in enumerate(bindings):
        field = _require_schema_field(schema, binding["field"], f"{label}[{index}]")
        tensor = _require_tensor(tensors, binding["tensor"], f"{label}[{index}]")
        if tensor["type"]["dtype"] != field["type"]["dtype"] or tensor["unit"] != field["unit"]:
            raise ValidationError(f"{label}[{index}] dtype/unit does not match field")
        suffix = field["type"]["shape"]
        if tensor["type"]["shape"] not in (suffix, [count] + suffix):
            raise ValidationError(f"{label}[{index}] uses an unsupported initializer shape")
    return bindings


def _evaluate_operations(
    operations: list[dict[str, Any]],
    tensors: dict[str, dict[str, Any]],
    target: dict[str, Any],
    limits: ValidationLimits,
) -> dict[str, int]:
    by_id: dict[str, dict[str, Any]] = {}
    next_unit = 0
    next_edge = 0
    attachment_count = 0
    seen_edges: set[tuple[str, int, int]] = set()
    attached: set[tuple[str, str]] = set()
    writes: set[tuple[str, str, str]] = set()
    bound_ports: set[str] = set()
    input_writes: set[tuple[int, str]] = set()
    for index, raw in enumerate(operations):
        label = f"development operation[{index}]"
        parsed = _operation_shape(raw, label, limits)
        operation_id = parsed["id"]
        if operation_id in by_id:
            raise ValidationError(f"duplicate development operation id {operation_id!r}")
        op = parsed["op"]
        if op == OP_UNIT_CREATE:
            schema_id = parsed["schema"]
            if schema_id not in target["unit_schemas"]:
                raise ValidationError(f"{label} references unknown unit schema")
            schema = target["unit_schemas"][schema_id]
            count = _u64_scalar(_require_tensor(tensors, parsed["count"], label), f"{label}.count")
            _validate_initializers(
                parsed["initializers"], schema, count, tensors, f"{label}.initializers", limits
            )
            if next_unit + count > limits.units or next_unit + count > SAFE_INTEGER:
                raise ValidationError(f"{label} exceeds unit ceiling")
            operation = {
                "kind": "unit",
                "id": operation_id,
                "schema": schema_id,
                "count": count,
                "first_id": next_unit,
            }
            next_unit += count
        elif op == OP_EDGE_CREATE:
            schema_id = parsed["schema"]
            if schema_id not in target["edge_schemas"]:
                raise ValidationError(f"{label} references unknown edge schema")
            source_set = by_id.get(parsed["sources"])
            target_set = by_id.get(parsed["targets"])
            if (
                source_set is None
                or source_set.get("kind") != "unit"
                or target_set is None
                or target_set.get("kind") != "unit"
            ):
                raise ValidationError(
                    f"{label} endpoints must reference earlier unit.create operations"
                )
            pairs = _u64_pairs(
                _require_tensor(tensors, parsed["pairs"], label), f"{label}.pairs"
            )
            for row, (source_local, target_local) in enumerate(pairs):
                if source_local >= source_set["count"] or target_local >= target_set["count"]:
                    raise ValidationError(f"{label}.pairs[{row}] endpoint is out of range")
                edge_key = (
                    schema_id,
                    source_set["first_id"] + source_local,
                    target_set["first_id"] + target_local,
                )
                if edge_key in seen_edges:
                    raise ValidationError(f"{label}.pairs[{row}] duplicates a typed edge")
                seen_edges.add(edge_key)
            _validate_initializers(
                parsed["initializers"],
                target["edge_schemas"][schema_id],
                len(pairs),
                tensors,
                f"{label}.initializers",
                limits,
            )
            if next_edge + len(pairs) > limits.edges or next_edge + len(pairs) > SAFE_INTEGER:
                raise ValidationError(f"{label} exceeds edge ceiling")
            operation = {
                "kind": "edge",
                "id": operation_id,
                "schema": schema_id,
                "sources": parsed["sources"],
                "targets": parsed["targets"],
                "count": len(pairs),
                "first_id": next_edge,
            }
            next_edge += len(pairs)
        elif op == OP_RULE_ATTACH:
            if parsed["rule"] not in target["rules"]:
                raise ValidationError(f"{label} references unknown target rule")
            rule = target["rules"][parsed["rule"]]
            if parsed["rule_version"] != rule["version"]:
                raise ValidationError(f"{label}.rule_version does not match target rule")
            subject = by_id.get(parsed["subject"])
            expected_kind = "unit" if rule["subject"] == "unit" else "edge"
            if (
                subject is None
                or subject.get("kind") != expected_kind
                or subject.get("schema") != rule["schema"]
            ):
                raise ValidationError(f"{label}.subject does not match rule subject/schema")
            parameters = _binding_list(
                parsed["parameters"],
                f"{label}.parameters",
                limits,
                first="id",
                second="tensor",
            )
            if [item["id"] for item in parameters] != [
                item["id"] for item in rule["parameters"]
            ]:
                raise ValidationError(f"{label} must bind every rule parameter exactly once")
            for binding, contract in zip(parameters, rule["parameters"]):
                tensor = _require_tensor(tensors, binding["tensor"], label)
                if tensor["type"] != contract["type"] or tensor["unit"] != contract["unit"]:
                    raise ValidationError(f"{label} parameter {binding['id']!r} type mismatch")
            attach_key = (parsed["rule"], parsed["subject"])
            if attach_key in attached:
                raise ValidationError(f"{label} duplicates a rule attachment")
            attached.add(attach_key)
            for access in rule["writes"]:
                key = (parsed["subject"], rule["phase"], access["field"])
                if key in writes:
                    raise ValidationError(f"{label} creates a same-phase rule write conflict")
                writes.add(key)
            if subject["kind"] == "edge":
                source_schema = target["unit_schemas"][
                    by_id[subject["sources"]]["schema"]
                ]
                target_schema = target["unit_schemas"][
                    by_id[subject["targets"]]["schema"]
                ]
                for access in rule["reads"]:
                    if access["scope"] == "source":
                        _require_schema_field(source_schema, access["field"], label)
                    elif access["scope"] == "target":
                        _require_schema_field(target_schema, access["field"], label)
            attachment_count += 1
            if attachment_count > limits.attachments:
                raise ValidationError(f"{label} exceeds attachment ceiling")
            operation = {"kind": "attachment", "id": operation_id}
        else:
            if parsed["port"] not in target["ports"]:
                raise ValidationError(f"{label} references unknown target port")
            port = target["ports"][parsed["port"]]
            if parsed["port"] in bound_ports:
                raise ValidationError(f"{label} binds a target port more than once")
            subject = by_id.get(parsed["subject"])
            if subject is None or subject.get("kind") != "unit":
                raise ValidationError(f"{label}.subject must reference an earlier unit set")
            schema = target["unit_schemas"][subject["schema"]]
            field = _require_schema_field(schema, parsed["field"], label)
            if port["direction"] == "input" and field["mutability"] != "state":
                raise ValidationError(f"{label} input port cannot bind a constant field")
            indices = _u64_vector(
                _require_tensor(tensors, parsed["indices"], label), f"{label}.indices"
            )
            if len(indices) != len(set(indices)):
                raise ValidationError(f"{label}.indices must be unique")
            if any(item >= subject["count"] for item in indices):
                raise ValidationError(f"{label}.indices contains an out-of-range index")
            expected_type = {
                "dtype": field["type"]["dtype"],
                "shape": [len(indices)] + field["type"]["shape"],
            }
            if port["type"] != expected_type or port["unit"] != field["unit"]:
                raise ValidationError(f"{label} port type/unit does not match field")
            if port["direction"] == "input":
                for item in indices:
                    key = (subject["first_id"] + item, parsed["field"])
                    if key in input_writes:
                        raise ValidationError(f"{label} overlaps another input binding")
                    input_writes.add(key)
            bound_ports.add(parsed["port"])
            if len(bound_ports) > limits.ports:
                raise ValidationError(f"{label} exceeds port ceiling")
            operation = {"kind": "port", "id": operation_id}
        by_id[operation_id] = operation
    if bound_ports != set(target["ports"]):
        raise ValidationError("development operations must bind every target port exactly once")
    return {
        "operations": len(operations),
        "tensor_bytes": sum(len(item["raw"]) for item in tensors.values()),
        "units": next_unit,
        "edges": next_edge,
        "attachments": attachment_count,
    }


def _module_sources(value: Any) -> tuple[str, dict[str, Any]]:
    parsed = _keys(value, set(SOURCE_ROLES), "development module.sources")
    _, sequence_sha = _source_reference(
        parsed["sequence"], "development module.sources.sequence"
    )
    for role in SOURCE_ROLES[1:]:
        binding = _keys(
            parsed[role], {"artifact_sha256"}, f"development module.sources.{role}"
        )
        _sha256(
            binding["artifact_sha256"],
            f"development module.sources.{role}.artifact_sha256",
        )
    return sequence_sha, parsed


def _requirements(
    value: Any, label: str, limits: ValidationLimits
) -> tuple[tuple[str, int], ...]:
    if type(value) is not list:
        raise ValidationError(f"{label} must be an array")
    result: list[tuple[str, int]] = []
    for index, item in enumerate(value):
        parsed = _keys(item, {"domain", "version"}, f"{label}[{index}]")
        result.append(
            (
                _domain(parsed["domain"], f"{label}[{index}].domain", limits),
                _integer(parsed["version"], f"{label}[{index}].version", minimum=1),
            )
        )
    if result != sorted(result) or len({item[0] for item in result}) != len(result):
        raise ValidationError(f"{label} must be sorted with unique domains")
    return tuple(result)


def _budgets(value: Any, label: str) -> dict[str, int]:
    parsed = _keys(
        value,
        {"operations", "tensor_bytes", "units", "edges", "attachments"},
        label,
    )
    return {name: _integer(parsed[name], f"{label}.{name}") for name in parsed}


def _module_info(
    value: dict[str, Any],
    target_info: dict[str, Any],
    *,
    blob_root: str | Path | None,
    limits: ValidationLimits,
) -> dict[str, Any]:
    top = _keys(
        value,
        {
            "format",
            "version",
            "producer",
            "sources",
            "module",
            "module_sha256",
            "artifact_sha256",
        },
        "development module",
    )
    if top["format"] != "brainc.development-module" or _integer(
        top["version"], "development module.version", minimum=1
    ) != 1:
        raise ValidationError("unsupported development module format/version")
    producer = _keys(
        top["producer"], {"name", "version", "passes"}, "development module.producer"
    )
    _identifier(producer["name"], "development module.producer.name", limits)
    _identifier(producer["version"], "development module.producer.version", limits)
    if type(producer["passes"]) is not list or len(producer["passes"]) > limits.operations:
        raise ValidationError("development module.producer.passes is invalid")
    for index, compiler_pass in enumerate(producer["passes"]):
        _identifier(
            compiler_pass, f"development module.producer.passes[{index}]", limits
        )
    sequence_sha, sources = _module_sources(top["sources"])
    body = _keys(
        top["module"],
        {"target", "requirements", "budgets", "tensors", "entrypoint"},
        "development module body",
    )
    stored_module = _sha256(top["module_sha256"], "development module.module_sha256")
    if stored_module != digest(body, limits=limits):
        raise ValidationError("development module.module_sha256 does not match module")
    artifact_sha = _artifact(top, "development module", limits)
    binding = _keys(
        body["target"],
        {"id", "abi_major", "contract_sha256"},
        "development module.target",
    )
    _domain(binding["id"], "development module.target.id", limits)
    _integer(binding["abi_major"], "development module.target.abi_major", minimum=1)
    _sha256(binding["contract_sha256"], "development module.target.contract_sha256")
    expected_binding = {
        "id": target_info["id"],
        "abi_major": target_info["abi_major"],
        "contract_sha256": target_info["contract_sha256"],
    }
    if binding != expected_binding:
        raise ValidationError("development module target binding does not match target")
    if sources["target_contract"]["artifact_sha256"] != target_info["artifact_sha256"]:
        raise ValidationError("development module does not bind supplied target artifact")
    requirements = _requirements(
        body["requirements"], "development module.requirements", limits
    )
    declared = _budgets(body["budgets"], "development module.budgets")
    ceilings = {
        "operations": limits.operations,
        "tensor_bytes": limits.total_tensor_bytes,
        "units": limits.units,
        "edges": limits.edges,
        "attachments": limits.attachments,
    }
    for name, maximum in ceilings.items():
        if declared[name] > maximum:
            raise ValidationError(f"development module budget {name} exceeds ceiling")
    tensors_raw = body["tensors"]
    if type(tensors_raw) is not list:
        raise ValidationError("development module.tensors must be an array")
    if len(tensors_raw) > limits.tensors:
        raise ValidationError("development module exceeds tensor-count ceiling")
    tensors_list = [
        _tensor_payload(
            item,
            f"development module.tensors[{index}]",
            limits,
            sequence_sha256=sequence_sha,
            records=None,
            blob_root=blob_root,
            lineage=True,
        )
        for index, item in enumerate(tensors_raw)
    ]
    _sorted_unique([item["id"] for item in tensors_list], "development module.tensors")
    response_outputs = [item["response_output"] for item in tensors_list]
    if len(response_outputs) != len(set(response_outputs)):
        raise ValidationError("development module tensor response lineages are duplicated")
    total_tensor_bytes = sum(len(item["raw"]) for item in tensors_list)
    if total_tensor_bytes > limits.total_tensor_bytes:
        raise ValidationError("development module exceeds total tensor byte ceiling")
    tensors = {item["id"]: item for item in tensors_list}
    entrypoint = _keys(
        body["entrypoint"], {"name", "operations"}, "development module.entrypoint"
    )
    if entrypoint["name"] != "develop":
        raise ValidationError("development module entrypoint must be develop")
    if type(entrypoint["operations"]) is not list:
        raise ValidationError("development module operations must be an array")
    exact = _evaluate_operations(entrypoint["operations"], tensors, target_info, limits)
    used_requirements = ((DEV_DOMAIN, DEV_VERSION),) if entrypoint["operations"] else ()
    if requirements != used_requirements:
        raise ValidationError("development module requirements do not equal used opsets")
    if any(requirement not in target_info["opsets"] for requirement in requirements):
        raise ValidationError("development module requires an opset absent from target")
    if declared != exact:
        raise ValidationError("development module budgets do not equal recomputed totals")
    return {
        "artifact": top,
        "target": target_info,
        "sources": sources,
        "requirements": requirements,
        "budgets": exact,
        "tensors": tensors,
        "module_sha256": stored_module,
        "artifact_sha256": artifact_sha,
    }


def validate_module_artifact(
    value: dict[str, Any],
    target: dict[str, Any],
    *,
    blob_root: str | Path | None = None,
    limits: ValidationLimits = DEFAULT_LIMITS,
) -> dict[str, Any]:
    """Validate a module against a raw, closed target-contract artifact."""

    target_info = validate_target_artifact(target, limits=limits)
    _module_info(value, target_info, blob_root=blob_root, limits=limits)
    return value


def _response_contract(output: dict[str, Any]) -> dict[str, Any]:
    return {
        "id": output["id"],
        "type": output["type"],
        "unit": output["unit"],
        "axes": output["axes"],
    }


def _chain_replay(
    source: dict[str, Any],
    manifest: dict[str, Any],
    request: dict[str, Any],
    response: dict[str, Any],
    response_tensors: dict[str, dict[str, Any]],
    policy: dict[str, Any],
    target: dict[str, Any],
    module: dict[str, Any],
    *,
    blob_root: str | Path | None,
    limits: ValidationLimits,
) -> dict[str, Any]:
    source_binding = _source_binding(source)
    if request["source"] != source_binding:
        raise ValidationError("prediction request is not bound to supplied sequence source")
    if request["provider_manifest_sha256"] != manifest["artifact_sha256"]:
        raise ValidationError("prediction request is not bound to provider manifest")
    source_tag = f'{source["format"]}/v{source["version"]}'
    if source_tag not in manifest["accepts"]:
        raise ValidationError("provider manifest does not declare source format support")
    declared_outputs = {item["id"]: item for item in manifest["outputs"]}
    for requested in request["requested_outputs"]:
        if (
            requested["id"] not in declared_outputs
            or requested != declared_outputs[requested["id"]]
        ):
            raise ValidationError("prediction request output is not declared by manifest")
    if response["request_artifact_sha256"] != request["artifact_sha256"]:
        raise ValidationError("prediction response is not bound to supplied request")
    if (
        response["provider"] != manifest["provider"]
        or response["model_identity"] != manifest["model_identity"]
    ):
        raise ValidationError("prediction response provider/model differs from manifest")
    response_contracts = [_response_contract(item) for item in response["outputs"]]
    if response_contracts != request["requested_outputs"]:
        raise ValidationError("prediction response contract/order differs from request")
    target_info = validate_target_artifact(target, limits=limits)
    policy_info = _policy_info(policy, limits)
    target_binding = {
        "id": target_info["id"],
        "abi_major": target_info["abi_major"],
        "contract_sha256": target_info["contract_sha256"],
    }
    if policy_info["target"] != target_binding:
        raise ValidationError("lowering policy is not bound to supplied target")
    if {item["from_output"] for item in policy_info["values"]} != set(response_tensors):
        raise ValidationError(
            "lowering policy values must cover every response output exactly once"
        )
    tensors: dict[str, dict[str, Any]] = {}
    expected_tensors: list[dict[str, Any]] = []
    response_json = {item["id"]: item for item in response["outputs"]}
    for binding in policy_info["values"]:
        output_id = binding["from_output"]
        if output_id not in response_tensors:
            raise ValidationError(
                f"lowering policy references unavailable provider output {output_id!r}"
            )
        source_tensor = response_tensors[output_id]
        tensor_id = binding["id"]
        tensor = {
            **source_tensor,
            "id": tensor_id,
            "response_output": output_id,
            "lowering_value": tensor_id,
        }
        tensors[tensor_id] = tensor
        source_json = response_json[output_id]
        expected_tensors.append(
            {
                "id": tensor_id,
                "type": source_json["type"],
                "unit": source_json["unit"],
                "axes": source_json["axes"],
                "storage": source_json["storage"],
                "lineage": {
                    "response_output": output_id,
                    "lowering_value": tensor_id,
                },
            }
        )
    expected_tensors.sort(key=lambda item: item["id"])
    exact_budgets = _evaluate_operations(
        policy_info["operations"], tensors, target_info, limits
    )
    expected_requirements = (
        [{"domain": DEV_DOMAIN, "version": DEV_VERSION}]
        if policy_info["operations"]
        else []
    )
    expected_body = {
        "target": target_binding,
        "requirements": expected_requirements,
        "budgets": exact_budgets,
        "tensors": expected_tensors,
        "entrypoint": {
            "name": "develop",
            "operations": policy["operations"],
        },
    }
    expected_sources = {
        "sequence": source_binding,
        "provider_manifest": {"artifact_sha256": manifest["artifact_sha256"]},
        "prediction_request": {"artifact_sha256": request["artifact_sha256"]},
        "prediction_response": {"artifact_sha256": response["artifact_sha256"]},
        "lowering_policy": {"artifact_sha256": policy["artifact_sha256"]},
        "target_contract": {"artifact_sha256": target["artifact_sha256"]},
    }
    module_info = _module_info(
        module, target_info, blob_root=blob_root, limits=limits
    )
    if module["producer"] != MODULE_PRODUCER:
        raise ValidationError("development module producer identity differs from compiler v0.5")
    if module["sources"] != expected_sources:
        raise ValidationError("development module source bindings differ from replay")
    if canonical_bytes(module["module"], limits=limits) != canonical_bytes(
        expected_body, limits=limits
    ):
        raise ValidationError("development module differs from independent lowering replay")
    return module_info


def validate_chain(
    *,
    fasta: str | Path,
    sequence: str | Path,
    manifest: str | Path,
    request: str | Path,
    response: str | Path,
    policy: str | Path,
    target: str | Path,
    module: str | Path,
    blob_root: str | Path | None = None,
    context: str | Path | None = None,
    record_id: str | None = None,
    limits: ValidationLimits = DEFAULT_LIMITS,
) -> dict[str, Any]:
    """Independently replay and validate the complete v2 compiler chain.

    The returned report is sealed even on failure.  A caller must check
    ``report["valid"]``; invalid inputs are not reported as successful partial
    validations.
    """

    limits.validate()
    named: dict[str, str | Path] = {
        "source_input": fasta,
        "source_artifact": sequence,
        "manifest": manifest,
        "request": request,
        "response": response,
        "policy": policy,
        "target": target,
        "module": module,
    }
    if blob_root is not None:
        named["blob_root"] = blob_root
    if context is not None:
        named["context"] = context
    checks: list[dict[str, Any]] = []
    hashes: dict[str, str] = {}
    failure_detail: str | None = None
    try:
        raw_source = _read_regular(fasta, "DNA source", limits.input_bytes)
        hashes["source_input"] = bytes_digest(raw_source)
        context_raw = (
            None
            if context is None
            else _read_regular(context, "sequence context", limits.input_bytes)
        )
        if context_raw is not None:
            hashes["context"] = bytes_digest(context_raw)
        source_value, raw = _load_source_artifact(
            sequence, "sequence artifact", limits
        )
        hashes["source_artifact"] = bytes_digest(raw)
        records = _validate_source_and_records(
            source_value,
            raw_source=raw_source,
            context_raw=context_raw,
            record_id=record_id,
            limits=limits,
        )
        checks.append({"name": "source-replay", "passed": True})

        manifest_value, raw = load(manifest, "provider manifest", limits=limits)
        hashes["manifest"] = bytes_digest(raw)
        validate_manifest_artifact(manifest_value, limits=limits)
        request_value, raw = load(request, "prediction request", limits=limits)
        hashes["request"] = bytes_digest(raw)
        validate_request_artifact(request_value, limits=limits)
        response_value, raw = load(response, "prediction response", limits=limits)
        hashes["response"] = bytes_digest(raw)
        response_value, response_tensors = validate_response_artifact(
            response_value,
            sequence_sha256=source_value["artifact_sha256"],
            records=records,
            blob_root=blob_root,
            limits=limits,
        )
        checks.append(
            {
                "name": "tensor-provider-contract-and-binding",
                "passed": True,
                "scope": "identity, tensor contracts, storage bytes, and bindings",
            }
        )

        policy_value, raw = load(policy, "lowering policy", limits=limits)
        hashes["policy"] = bytes_digest(raw)
        validate_policy_artifact(policy_value, limits=limits)
        target_value, raw = load(target, "target contract", limits=limits)
        hashes["target"] = bytes_digest(raw)
        validate_target_artifact(target_value, limits=limits)
        checks.append({"name": "target-and-lowering-policy", "passed": True})

        module_value, raw = load(module, "development module", limits=limits)
        hashes["module"] = bytes_digest(raw)
        _chain_replay(
            source_value,
            manifest_value,
            request_value,
            response_value,
            response_tensors,
            policy_value,
            target_value,
            module_value,
            blob_root=blob_root,
            limits=limits,
        )
        checks.append(
            {
                "name": "development-lowering-and-emission-replay",
                "passed": True,
            }
        )
    except (OSError, ValidationError) as failure:
        failure_detail = str(failure)
        checks.append(
            {"name": "compiler-v2-chain", "passed": False, "detail": failure_detail}
        )
    valid = failure_detail is None and len(checks) == 4 and all(
        item["passed"] for item in checks
    )
    core = {
        "format": "brainc.validation-report",
        "version": 2,
        "stage": "tensor-compiler-chain",
        "inputs": {key: str(value) for key, value in sorted(named.items())},
        "input_sha256": dict(sorted(hashes.items())),
        "valid": valid,
        "passed_checks": sum(bool(item["passed"]) for item in checks),
        "total_checks": 4,
        "checks": checks,
        "validated_scope": (
            "DNA replay, RFC 8785 identities, closed provider/request/response/policy/target "
            "contracts, tensor bytes and axes, all ABI-v1 development operations, exact "
            "budgets, bindings, and independently reconstructed module emission"
        ),
    }
    return {**core, "report_sha256": digest(core, limits=limits)}


def save_report(report: dict[str, Any], path: str | Path) -> None:
    """Atomically write a report without following an output symlink."""

    if type(report) is not dict:
        raise ValidationError("validation report must be an object")
    try:
        rendered = (
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
        raise ValidationError(f"validation report is not I-JSON: {failure}") from failure
    destination = Path(path)
    temporary: str | None = None
    descriptor = -1
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.is_symlink():
            raise ValidationError("validation report path must not be a symbolic link")
        descriptor, temporary = tempfile.mkstemp(
            prefix=f".{destination.name}.", dir=destination.parent
        )
        with os.fdopen(descriptor, "wb") as stream:
            descriptor = -1
            stream.write(rendered)
            stream.flush()
            os.fsync(stream.fileno())
        if destination.is_symlink():
            raise ValidationError("validation report path became a symbolic link")
        os.replace(temporary, destination)
        temporary = None
    except ValidationError:
        raise
    except OSError as failure:
        raise ValidationError(f"cannot write validation report: {failure}") from failure
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        if temporary is not None:
            try:
                os.unlink(temporary)
            except OSError:
                pass


__all__ = [
    "DEFAULT_LIMITS",
    "MODULE_PRODUCER",
    "ValidationError",
    "ValidationLimits",
    "canonical_bytes",
    "digest",
    "load",
    "loads",
    "main",
    "save_report",
    "source_record_lengths",
    "validate_chain",
    "validate_manifest_artifact",
    "validate_module_artifact",
    "validate_policy_artifact",
    "validate_request_artifact",
    "validate_response_artifact",
    "validate_source_artifact",
    "validate_target_artifact",
]


def main(argv: list[str] | None = None) -> int:
    """Standalone validator entry point used by ``python -m`` and the CLI."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_input", type=Path)
    parser.add_argument("sequence", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("request", type=Path)
    parser.add_argument("response", type=Path)
    parser.add_argument("policy", type=Path)
    parser.add_argument("target", type=Path)
    parser.add_argument("module", type=Path)
    parser.add_argument("--blob-root", type=Path)
    parser.add_argument("--context", type=Path)
    parser.add_argument("--record-id")
    parser.add_argument("--report", type=Path)
    arguments = parser.parse_args(argv)
    report = validate_chain(
        fasta=arguments.source_input,
        sequence=arguments.sequence,
        manifest=arguments.manifest,
        request=arguments.request,
        response=arguments.response,
        policy=arguments.policy,
        target=arguments.target,
        module=arguments.module,
        blob_root=arguments.blob_root,
        context=arguments.context,
        record_id=arguments.record_id,
    )
    if arguments.report is not None:
        save_report(report, arguments.report)
    print(json.dumps(report, sort_keys=True, ensure_ascii=False, allow_nan=False))
    return 0 if report["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

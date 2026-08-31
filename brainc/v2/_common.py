"""Closed-contract and bounded-I/O helpers for brainc v2."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from brainc._canonical import (
    SAFE_INTEGER,
    ContractError,
    artifact_digest,
    canonical_bytes,
    digest,
    keys,
    sha256,
)
from brainc._io import BoundedIOError, atomic_write_file, read_regular_file

from .limits import (
    MAX_IDENTIFIER_BYTES,
    MAX_JSON_BYTES,
    MAX_JSON_DEPTH,
    MAX_JSON_MEMBERS,
    MAX_SOURCE_JSON_BYTES,
    MAX_STRING_BYTES,
)


# Canonicalization is part of the v2 trust boundary, not an implementation
# detail.  Expose its one contract-error type under the v2 API name so callers
# can catch every malformed artifact uniformly, including failures raised by
# shared RFC 8785 helpers.
V2Error = ContractError


def integer(value: Any, label: str, *, minimum: int = 0, maximum: int = SAFE_INTEGER) -> int:
    if type(value) is not int or value < minimum or value > maximum:
        raise V2Error(f"{label} must be an integer in [{minimum}, {maximum}]")
    return value


def text(value: Any, label: str, *, identifier: bool = False) -> str:
    if type(value) is not str or not value or value != value.strip():
        raise V2Error(f"{label} must be a nonempty trimmed string")
    if any(0xD800 <= ord(character) <= 0xDFFF for character in value):
        raise V2Error(f"{label} contains a lone Unicode surrogate")
    ceiling = MAX_IDENTIFIER_BYTES if identifier else MAX_STRING_BYTES
    if len(value.encode("utf-8")) > ceiling:
        raise V2Error(f"{label} exceeds its UTF-8 byte ceiling")
    return value


def exact_version(value: Any, wanted: int, label: str) -> int:
    if type(value) is not int or value != wanted:
        raise V2Error(f"{label} must be exactly {wanted}")
    return value


def checked_product(shape: Any, label: str) -> int:
    if type(shape) is not list:
        raise V2Error(f"{label} must be an array")
    result = 1
    for index, dimension in enumerate(shape):
        result *= integer(dimension, f"{label}[{index}]")
        if result > SAFE_INTEGER:
            raise V2Error(f"{label} product exceeds the I-JSON safe range")
    return result


def sorted_unique(values: list[str], label: str) -> None:
    if values != sorted(values):
        raise V2Error(f"{label} must be sorted")
    if len(values) != len(set(values)):
        raise V2Error(f"{label} contains duplicates")


def _reject_constant(value: str) -> None:
    raise V2Error(f"non-finite JSON number is not allowed: {value}")


def _parse_integer(value: str) -> int:
    parsed = int(value)
    if not -SAFE_INTEGER <= parsed <= SAFE_INTEGER:
        raise V2Error("JSON integer exceeds the I-JSON safe range")
    return parsed


def _parse_float(value: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed):
        raise V2Error("non-finite JSON number is not allowed")
    return parsed


def _no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise V2Error(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def validate_tree(value: Any, label: str) -> None:
    members = 0
    stack: list[tuple[Any, int]] = [(value, 0)]
    while stack:
        current, depth = stack.pop()
        if depth > MAX_JSON_DEPTH:
            raise V2Error(f"{label} exceeds JSON depth limit {MAX_JSON_DEPTH}")
        if current is None or type(current) is bool:
            continue
        if type(current) is int:
            if not -SAFE_INTEGER <= current <= SAFE_INTEGER:
                raise V2Error(f"{label} contains an unsafe JSON integer")
            continue
        if type(current) is float:
            if not math.isfinite(current):
                raise V2Error(f"{label} contains a non-finite JSON number")
            continue
        if type(current) is str:
            if any(0xD800 <= ord(character) <= 0xDFFF for character in current):
                raise V2Error(f"{label} contains a lone Unicode surrogate")
            if len(current.encode("utf-8")) > MAX_STRING_BYTES:
                raise V2Error(f"{label} exceeds JSON string byte limit {MAX_STRING_BYTES}")
            continue
        if type(current) is list:
            members += len(current)
            stack.extend((item, depth + 1) for item in reversed(current))
        elif type(current) is dict:
            members += len(current)
            for key, item in reversed(list(current.items())):
                if type(key) is not str:
                    raise V2Error(f"{label} contains a non-string object key")
                stack.append((key, depth + 1))
                stack.append((item, depth + 1))
        else:
            raise V2Error(f"{label} contains unsupported type {type(current).__name__}")
        if members > MAX_JSON_MEMBERS:
            raise V2Error(f"{label} exceeds JSON member limit {MAX_JSON_MEMBERS}")


def _json_byte_limit(value: int) -> int:
    if type(value) is not int or not 1 <= value <= MAX_SOURCE_JSON_BYTES:
        raise V2Error(
            f"JSON byte limit must be an integer in [1, {MAX_SOURCE_JSON_BYTES}]"
        )
    return value


def loads(
    raw: bytes,
    label: str,
    *,
    maximum_bytes: int = MAX_JSON_BYTES,
) -> dict[str, Any]:
    if type(raw) is not bytes:
        raise V2Error(f"{label} JSON input must be bytes")
    maximum_bytes = _json_byte_limit(maximum_bytes)
    if len(raw) > maximum_bytes:
        raise V2Error(f"{label} exceeds JSON byte limit {maximum_bytes}")
    try:
        value = json.loads(
            raw.decode("utf-8"),
            object_pairs_hook=_no_duplicates,
            parse_constant=_reject_constant,
            parse_int=_parse_integer,
            parse_float=_parse_float,
        )
    except V2Error:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError, RecursionError) as failure:
        raise V2Error(f"invalid {label} JSON: {failure}") from failure
    if type(value) is not dict:
        raise V2Error(f"{label} must be a JSON object")
    validate_tree(value, label)
    return value


def load(
    path: str | Path,
    label: str,
    *,
    maximum_bytes: int = MAX_JSON_BYTES,
) -> tuple[dict[str, Any], bytes]:
    maximum_bytes = _json_byte_limit(maximum_bytes)
    try:
        raw = read_regular_file(
            path,
            maximum_bytes=maximum_bytes,
            label=label,
        )
    except BoundedIOError as failure:
        raise V2Error(str(failure)) from failure
    return loads(raw, label, maximum_bytes=maximum_bytes), raw


def _check_artifact(value: dict[str, Any], label: str) -> dict[str, Any]:
    if "artifact_sha256" not in value:
        raise V2Error(f"{label} is missing artifact_sha256")
    try:
        artifact_digest(value)
    except (ContractError, KeyError) as failure:
        raise V2Error(str(failure)) from failure
    return value


def checked_artifact(
    path: str | Path,
    label: str,
    *,
    maximum_bytes: int = MAX_JSON_BYTES,
) -> dict[str, Any]:
    value, _ = load(path, label, maximum_bytes=maximum_bytes)
    return _check_artifact(value, label)


def pretty_bytes(payload: dict[str, Any]) -> bytes:
    validate_tree(payload, "output artifact")
    try:
        raw = (
            json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False)
            + "\n"
        ).encode("utf-8")
        canonical_bytes(payload)
    except (TypeError, ValueError, UnicodeError, ContractError) as failure:
        raise V2Error(f"output is not canonicalizable I-JSON: {failure}") from failure
    if len(raw) > MAX_JSON_BYTES:
        raise V2Error(f"serialized output exceeds JSON byte limit {MAX_JSON_BYTES}")
    return raw


def save(payload: dict[str, Any], path: str | Path) -> None:
    raw = pretty_bytes(payload)
    try:
        atomic_write_file(
            path,
            raw,
            maximum_bytes=MAX_JSON_BYTES,
            label="v2 artifact output",
        )
    except BoundedIOError as failure:
        raise V2Error(str(failure)) from failure


def seal(core: dict[str, Any]) -> dict[str, Any]:
    validate_tree(core, "artifact")
    return {**core, "artifact_sha256": digest(core)}


__all__ = [
    "SAFE_INTEGER",
    "V2Error",
    "artifact_digest",
    "checked_artifact",
    "checked_product",
    "digest",
    "exact_version",
    "integer",
    "keys",
    "load",
    "pretty_bytes",
    "save",
    "seal",
    "sha256",
    "sorted_unique",
    "text",
    "validate_tree",
]

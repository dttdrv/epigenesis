"""Strict deterministic JSON helpers shared by the compiler implementation."""

from __future__ import annotations

from collections.abc import Iterator
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any

from ._io import (
    BoundedIOError,
    MAX_JSON_BYTES,
    MAX_JSON_DEPTH,
    MAX_JSON_MEMBERS,
    MAX_STRING_BYTES,
    atomic_write_file,
    load_json_object,
    pretty_json_bytes,
    read_regular_file,
)


SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
SAFE_INTEGER = 2**53 - 1


class ContractError(ValueError):
    """An input violates a closed compiler artifact contract."""


def reject_constant(value: str) -> None:
    raise ContractError(f"non-finite JSON number is not allowed: {value}")


def no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ContractError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def loads(raw: bytes, label: str) -> dict[str, Any]:
    try:
        return load_json_object(
            raw,
            label,
            maximum_bytes=MAX_JSON_BYTES,
            maximum_depth=MAX_JSON_DEPTH,
            maximum_members=MAX_JSON_MEMBERS,
            maximum_string_bytes=MAX_STRING_BYTES,
        )
    except BoundedIOError as failure:
        raise ContractError(f"invalid {label} JSON: {failure}") from failure


def load(path: str | Path, label: str) -> tuple[dict[str, Any], bytes]:
    try:
        raw = read_regular_file(path, maximum_bytes=MAX_JSON_BYTES, label=label)
    except BoundedIOError as failure:
        raise ContractError(f"cannot read {label}: {failure}") from failure
    return loads(raw, label), raw


def canonical_bytes(value: Any) -> bytes:
    """Serialize an I-JSON value with RFC 8785 JSON canonicalization."""
    try:
        return "".join(_jcs_chunks(value)).encode("utf-8")
    except (TypeError, ValueError, UnicodeError, RecursionError) as failure:
        raise ContractError(f"value is not RFC 8785 canonical JSON: {failure}") from failure


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
    return sign + rendered + ("e+" if normalized_exponent >= 0 else "e") + str(normalized_exponent)


def _jcs_chunks(value: Any) -> Iterator[str]:
    """Yield RFC 8785 text without constructing every enclosing subtree."""

    if value is None:
        yield "null"
    elif type(value) is bool:
        yield "true" if value else "false"
    elif type(value) is str:
        yield _jcs_string(value)
    elif type(value) is int:
        if not -SAFE_INTEGER <= value <= SAFE_INTEGER:
            raise ValueError(f"integer exceeds I-JSON safe range: {value}")
        yield str(value)
    elif type(value) is float:
        yield _jcs_number(value)
    elif type(value) is list:
        yield "["
        for index, item in enumerate(value):
            if index:
                yield ","
            yield from _jcs_chunks(item)
        yield "]"
    elif type(value) is dict:
        if any(type(key) is not str for key in value):
            raise TypeError("JSON object keys must be strings")
        yield "{"
        for index, key in enumerate(sorted(value, key=lambda item: item.encode("utf-16be"))):
            if index:
                yield ","
            yield _jcs_string(key)
            yield ":"
            yield from _jcs_chunks(value[key])
        yield "}"
    else:
        raise TypeError(f"unsupported JSON value type: {type(value).__name__}")


def digest(value: Any) -> str:
    hasher = hashlib.sha256()
    try:
        pending: list[str] = []
        pending_characters = 0
        for chunk in _jcs_chunks(value):
            pending.append(chunk)
            pending_characters += len(chunk)
            if pending_characters >= 64 * 1024:
                hasher.update("".join(pending).encode("utf-8"))
                pending.clear()
                pending_characters = 0
        if pending:
            hasher.update("".join(pending).encode("utf-8"))
    except (TypeError, ValueError, UnicodeError, RecursionError) as failure:
        raise ContractError(f"value is not RFC 8785 canonical JSON: {failure}") from failure
    return hasher.hexdigest()


def bytes_digest(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def keys(value: Any, expected: set[str], label: str) -> dict[str, Any]:
    if type(value) is not dict:
        raise ContractError(f"{label} must be an object")
    missing = sorted(expected - value.keys())
    extra = sorted(value.keys() - expected)
    if missing or extra:
        raise ContractError(f"{label} keys invalid; missing={missing or 'none'}, unknown={extra or 'none'}")
    return value


def text(value: Any, label: str, *, allow_empty: bool = False) -> str:
    if type(value) is not str or (not allow_empty and not value) or value != value.strip():
        raise ContractError(f"{label} must be a {'possibly empty ' if allow_empty else ''}trimmed string")
    return value


def sha256(value: Any, label: str) -> str:
    if type(value) is not str or SHA256_RE.fullmatch(value) is None:
        raise ContractError(f"{label} must be a lowercase SHA-256 digest")
    return value


def number(value: Any, label: str) -> int | float:
    if type(value) is int:
        if not -SAFE_INTEGER <= value <= SAFE_INTEGER:
            raise ContractError(f"{label} exceeds the I-JSON safe integer range")
        return value
    if type(value) is not float or not math.isfinite(value):
        raise ContractError(f"{label} must be a finite number")
    return value


def artifact_digest(payload: dict[str, Any], field: str = "artifact_sha256") -> str:
    stored = sha256(payload[field], field)
    expected = digest({key: value for key, value in payload.items() if key != field})
    if stored != expected:
        raise ContractError(f"{field} does not match canonical artifact content")
    return stored


def save_artifact(payload: dict[str, Any], path: str | Path) -> None:
    try:
        raw = pretty_json_bytes(
            payload,
            ensure_ascii=False,
            maximum_bytes=MAX_JSON_BYTES,
            maximum_depth=MAX_JSON_DEPTH,
            maximum_members=MAX_JSON_MEMBERS,
            maximum_string_bytes=MAX_STRING_BYTES,
        )
        atomic_write_file(
            path,
            raw,
            maximum_bytes=MAX_JSON_BYTES,
            label="artifact output",
        )
    except BoundedIOError as failure:
        raise ContractError(f"cannot save artifact: {failure}") from failure

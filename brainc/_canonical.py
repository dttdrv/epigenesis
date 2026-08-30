"""Strict deterministic JSON helpers shared by the compiler implementation."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any


SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")


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
        value = json.loads(
            raw.decode("utf-8"),
            object_pairs_hook=no_duplicates,
            parse_constant=reject_constant,
        )
    except ContractError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as failure:
        raise ContractError(f"invalid {label} JSON: {failure}") from failure
    if type(value) is not dict:
        raise ContractError(f"{label} must be a JSON object")
    return value


def load(path: str | Path, label: str) -> tuple[dict[str, Any], bytes]:
    try:
        raw = Path(path).read_bytes()
    except OSError as failure:
        raise ContractError(f"cannot read {label}: {failure}") from failure
    return loads(raw, label), raw


def canonical_bytes(value: Any) -> bytes:
    try:
        return json.dumps(
            value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False
        ).encode("utf-8")
    except (TypeError, ValueError) as failure:
        raise ContractError(f"value is not deterministic JSON: {failure}") from failure


def digest(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


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
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ContractError(f"{label} must be a finite number")
    return value


def artifact_digest(payload: dict[str, Any], field: str = "artifact_sha256") -> str:
    stored = sha256(payload[field], field)
    expected = digest({key: value for key, value in payload.items() if key != field})
    if stored != expected:
        raise ContractError(f"{field} does not match canonical artifact content")
    return stored


def save_artifact(payload: dict[str, Any], path: str | Path) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )

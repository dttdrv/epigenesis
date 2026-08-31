"""Tensor wire types, canonical byte packing, and storage validation."""

from __future__ import annotations

import base64
import binascii
from dataclasses import dataclass
import hashlib
import math
import os
from pathlib import Path
import re
import stat
import struct
from typing import Any, Iterable

from brainc._canonical import SAFE_INTEGER, bytes_digest, keys, sha256

from ._common import V2Error, checked_product, integer, text
from .limits import MAX_TENSOR_BYTES, MAX_TENSOR_RANK


DTYPE_BYTES = {"bool": 1, "i64": 8, "u64": 8, "f64": 8}
_STRUCT_CODES = {"i64": "q", "u64": "Q", "f64": "d"}
_DOMAIN_SEGMENT = r"[a-z0-9](?:[a-z0-9-]*[a-z0-9])?"
_DOMAIN_RE = re.compile(rf"(?:{_DOMAIN_SEGMENT}\.)+{_DOMAIN_SEGMENT}\Z")


@dataclass(frozen=True)
class TensorType:
    dtype: str
    shape: tuple[int, ...]

    @property
    def elements(self) -> int:
        result = 1
        for dimension in self.shape:
            result *= dimension
        return result

    @property
    def byte_length(self) -> int:
        return self.elements * DTYPE_BYTES[self.dtype]


@dataclass(frozen=True)
class ParsedTensor:
    wire: dict[str, Any]
    id: str
    type: TensorType
    unit: str | None
    axes: tuple[dict[str, Any] | None, ...]
    raw: bytes


def parse_type(value: Any, label: str) -> TensorType:
    item = keys(value, {"dtype", "shape"}, label)
    dtype = text(item["dtype"], f"{label}.dtype", identifier=True)
    if dtype not in DTYPE_BYTES:
        raise V2Error(f"{label}.dtype is unsupported: {dtype!r}")
    shape = item["shape"]
    if type(shape) is not list:
        raise V2Error(f"{label}.shape must be an array")
    if len(shape) > MAX_TENSOR_RANK:
        raise V2Error(f"{label}.shape exceeds rank ceiling {MAX_TENSOR_RANK}")
    elements = checked_product(shape, f"{label}.shape")
    byte_length = elements * DTYPE_BYTES[dtype]
    if byte_length > SAFE_INTEGER or byte_length > MAX_TENSOR_BYTES:
        raise V2Error(f"{label} exceeds tensor byte ceiling {MAX_TENSOR_BYTES}")
    return TensorType(dtype, tuple(shape))


def parse_unit(value: Any, label: str) -> str | None:
    if value is None:
        return None
    parsed = text(value, label, identifier=True)
    if _DOMAIN_RE.fullmatch(parsed) is None:
        raise V2Error(f"{label} must be a lowercase dotted namespaced identifier")
    return parsed


def parse_axis(
    value: Any,
    label: str,
    *,
    source_artifact_sha256: str | None = None,
    records: dict[str, int] | None = None,
) -> dict[str, Any] | None:
    if value is None:
        return None
    item = keys(
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
    if item["kind"] != "source-coordinate" or item["coordinate_system"] != "0-based-half-open":
        raise V2Error(f"{label} has unsupported coordinate semantics")
    source_sha = sha256(item["source_artifact_sha256"], f"{label}.source_artifact_sha256")
    if source_artifact_sha256 is not None and source_sha != source_artifact_sha256:
        raise V2Error(f"{label} does not bind the supplied sequence source")
    record_id = text(item["record_id"], f"{label}.record_id", identifier=True)
    start = integer(item["start"], f"{label}.start")
    step = integer(item["step"], f"{label}.step", minimum=1)
    span = integer(item["span"], f"{label}.span", minimum=1)
    if records is not None:
        if record_id not in records:
            raise V2Error(f"{label} references unknown sequence record {record_id!r}")
        if start + span > records[record_id]:
            raise V2Error(f"{label} starts outside its sequence record")
    return dict(item)


def parse_axes(
    value: Any,
    tensor_type: TensorType,
    label: str,
    *,
    source_artifact_sha256: str | None = None,
    records: dict[str, int] | None = None,
) -> tuple[dict[str, Any] | None, ...]:
    if type(value) is not list or len(value) != len(tensor_type.shape):
        raise V2Error(f"{label} must be an array matching tensor rank")
    result = tuple(
        parse_axis(
            axis,
            f"{label}[{index}]",
            source_artifact_sha256=source_artifact_sha256,
            records=records,
        )
        for index, axis in enumerate(value)
    )
    for dimension, axis in zip(tensor_type.shape, result):
        if axis is None or dimension == 0:
            continue
        endpoint = axis["start"] + (dimension - 1) * axis["step"] + axis["span"]
        if endpoint > SAFE_INTEGER:
            raise V2Error(f"{label} coordinate exceeds the I-JSON safe range")
        if records is not None and endpoint > records[axis["record_id"]]:
            raise V2Error(f"{label} coordinate exceeds its sequence record")
    return result


def pack(dtype: str, values: Iterable[Any]) -> bytes:
    """Pack exact scalar values in the ABI's little-endian representation."""
    if type(dtype) is not str or dtype not in DTYPE_BYTES:
        raise V2Error(f"unsupported tensor dtype: {dtype!r}")
    if isinstance(values, (str, bytes, bytearray, memoryview)):
        raise V2Error("tensor values must be an iterable of typed scalars")
    try:
        iterator = iter(values)
    except TypeError as failure:
        raise V2Error("tensor values must be iterable") from failure
    maximum_items = MAX_TENSOR_BYTES // DTYPE_BYTES[dtype]
    items: list[Any] = []
    for value in iterator:
        if len(items) >= maximum_items:
            raise V2Error(f"packed tensor exceeds byte ceiling {MAX_TENSOR_BYTES}")
        items.append(value)
    if dtype == "bool":
        if any(type(value) is not bool for value in items):
            raise V2Error("bool tensors require exact bool values")
        return bytes(1 if value else 0 for value in items)
    if dtype in {"i64", "u64"}:
        lower, upper = (-(2**63), 2**63 - 1) if dtype == "i64" else (0, 2**64 - 1)
        if any(type(value) is not int or value < lower or value > upper for value in items):
            raise V2Error(f"{dtype} tensor value is outside [{lower}, {upper}]")
    else:
        if any(type(value) is not float or not math.isfinite(value) or value == 0.0 and math.copysign(1.0, value) < 0 for value in items):
            raise V2Error("f64 tensors require finite float values and forbid negative zero")
    try:
        return struct.pack("<" + _STRUCT_CODES[dtype] * len(items), *items)
    except (struct.error, OverflowError) as failure:
        raise V2Error(f"cannot pack {dtype} tensor: {failure}") from failure


def inline_storage(raw: bytes) -> dict[str, Any]:
    if type(raw) is not bytes:
        raise V2Error("inline tensor storage requires exact bytes")
    if len(raw) > MAX_TENSOR_BYTES:
        raise V2Error(f"inline tensor exceeds byte ceiling {MAX_TENSOR_BYTES}")
    return {
        "kind": "inline-base64",
        "data": base64.b64encode(raw).decode("ascii"),
        "byte_length": len(raw),
        "sha256": bytes_digest(raw),
    }


def _read_blob(root: str | Path | None, blob_sha256: str, expected_length: int, label: str) -> bytes:
    if root is None:
        raise V2Error(f"{label} requires an explicit blob root")
    if os.name != "posix":
        raise V2Error(
            f"{label} external blob storage requires a POSIX host; use inline-base64"
        )
    root_path = Path(root)
    root_descriptor = -1
    descriptor = -1
    try:
        root_metadata = root_path.lstat()
        if stat.S_ISLNK(root_metadata.st_mode) or not stat.S_ISDIR(root_metadata.st_mode):
            raise V2Error("blob root must be a non-linked directory")
        root_descriptor = os.open(
            root_path,
            os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
        )
        opened_root = os.fstat(root_descriptor)
        if (
            not stat.S_ISDIR(opened_root.st_mode)
            or (opened_root.st_dev, opened_root.st_ino)
            != (root_metadata.st_dev, root_metadata.st_ino)
        ):
            raise V2Error("blob root changed or is not a safe directory")
        before = os.stat(blob_sha256, dir_fd=root_descriptor, follow_symlinks=False)
        if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode) or before.st_nlink != 1:
            raise V2Error(f"{label} must be one regular non-linked file")
        if before.st_size != expected_length:
            raise V2Error(f"{label} length does not match storage descriptor")
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
            raise V2Error(f"{label} changed or is not a safe regular file")
        chunks: list[bytes] = []
        remaining = expected_length
        while remaining:
            chunk = os.read(descriptor, min(remaining, 1024 * 1024))
            if not chunk:
                raise V2Error(f"{label} ended before its declared length")
            chunks.append(chunk)
            remaining -= len(chunk)
        if os.read(descriptor, 1):
            raise V2Error(f"{label} exceeds its declared length")
        raw = b"".join(chunks)
    except V2Error:
        raise
    except OSError as failure:
        raise V2Error(f"cannot read {label}: {failure}") from failure
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        if root_descriptor >= 0:
            os.close(root_descriptor)
    if hashlib.sha256(raw).hexdigest() != blob_sha256:
        raise V2Error(f"{label} digest mismatch")
    return raw


def parse_storage(
    value: Any,
    expected_length: int,
    label: str,
    *,
    blob_root: str | Path | None = None,
) -> bytes:
    if type(value) is not dict or "kind" not in value:
        raise V2Error(f"{label} must be a storage descriptor")
    if expected_length > MAX_TENSOR_BYTES:
        raise V2Error(f"{label} exceeds tensor byte ceiling")
    if value["kind"] == "inline-base64":
        item = keys(value, {"kind", "data", "byte_length", "sha256"}, label)
        declared = integer(item["byte_length"], f"{label}.byte_length")
        if declared != expected_length:
            raise V2Error(f"{label}.byte_length does not match tensor type")
        data = item["data"]
        if type(data) is not str or any(character.isspace() for character in data):
            raise V2Error(f"{label}.data must be strict padded base64")
        if len(data) != ((expected_length + 2) // 3) * 4:
            raise V2Error(f"{label}.data length does not match tensor type")
        try:
            raw = base64.b64decode(data, validate=True)
        except (ValueError, binascii.Error) as failure:
            raise V2Error(f"{label}.data is invalid base64") from failure
        if base64.b64encode(raw).decode("ascii") != data:
            raise V2Error(f"{label}.data is not canonical padded base64")
        stored_sha = sha256(item["sha256"], f"{label}.sha256")
        if len(raw) != expected_length or bytes_digest(raw) != stored_sha:
            raise V2Error(f"{label} bytes do not match descriptor")
        return raw
    if value["kind"] == "sha256-blob":
        item = keys(value, {"kind", "byte_length", "sha256"}, label)
        declared = integer(item["byte_length"], f"{label}.byte_length")
        if declared != expected_length:
            raise V2Error(f"{label}.byte_length does not match tensor type")
        stored_sha = sha256(item["sha256"], f"{label}.sha256")
        return _read_blob(blob_root, stored_sha, expected_length, label)
    raise V2Error(f"{label}.kind is unsupported")


def validate_elements(raw: bytes, dtype: str, label: str) -> None:
    if dtype == "bool":
        if any(value not in (0, 1) for value in raw):
            raise V2Error(f"{label} contains an invalid bool byte")
        return
    if dtype != "f64":
        return
    for offset in range(0, len(raw), 8):
        bits = int.from_bytes(raw[offset : offset + 8], "little")
        if bits & 0x7FF0000000000000 == 0x7FF0000000000000:
            raise V2Error(f"{label} contains a non-finite f64")
        if bits == 0x8000000000000000:
            raise V2Error(f"{label} contains negative-zero f64")


def parse_tensor_output(
    value: Any,
    label: str,
    *,
    blob_root: str | Path | None,
    source_artifact_sha256: str,
    records: dict[str, int],
) -> ParsedTensor:
    item = keys(value, {"id", "type", "unit", "axes", "storage"}, label)
    tensor_type = parse_type(item["type"], f"{label}.type")
    axes = parse_axes(
        item["axes"],
        tensor_type,
        f"{label}.axes",
        source_artifact_sha256=source_artifact_sha256,
        records=records,
    )
    raw = parse_storage(
        item["storage"], tensor_type.byte_length, f"{label}.storage", blob_root=blob_root
    )
    validate_elements(raw, tensor_type.dtype, label)
    return ParsedTensor(
        wire=dict(item),
        id=text(item["id"], f"{label}.id", identifier=True),
        type=tensor_type,
        unit=parse_unit(item["unit"], f"{label}.unit"),
        axes=axes,
        raw=raw,
    )


def unpack_u64_scalar(tensor: ParsedTensor, label: str) -> int:
    if tensor.type != TensorType("u64", ()) or tensor.unit is not None:
        raise V2Error(f"{label} must reference a unitless u64 scalar")
    value = struct.unpack("<Q", tensor.raw)[0]
    if value > SAFE_INTEGER:
        raise V2Error(f"{label} exceeds the ABI safe-integer range")
    return value


def unpack_u64_vector(tensor: ParsedTensor, label: str) -> tuple[int, ...]:
    if tensor.type.dtype != "u64" or len(tensor.type.shape) != 1 or tensor.unit is not None:
        raise V2Error(f"{label} must reference a unitless u64 vector")
    values = tuple(value[0] for value in struct.iter_unpack("<Q", tensor.raw))
    if any(value > SAFE_INTEGER for value in values):
        raise V2Error(f"{label} contains an ABI-unsafe integer")
    return values


def unpack_u64_pairs(tensor: ParsedTensor, label: str) -> tuple[tuple[int, int], ...]:
    if (
        tensor.type.dtype != "u64"
        or len(tensor.type.shape) != 2
        or tensor.type.shape[1] != 2
        or tensor.unit is not None
    ):
        raise V2Error(f"{label} must reference a unitless u64[E,2] tensor")
    flat = tuple(value[0] for value in struct.iter_unpack("<Q", tensor.raw))
    if any(value > SAFE_INTEGER for value in flat):
        raise V2Error(f"{label} contains an ABI-unsafe integer")
    return tuple((flat[index], flat[index + 1]) for index in range(0, len(flat), 2))

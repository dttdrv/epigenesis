"""Finite, race-resistant byte and JSON I/O used at compiler boundaries."""

from __future__ import annotations

import errno
import gzip
import io
import json
import math
import os
from pathlib import Path
import secrets
import stat
import tempfile
from typing import Any
import zlib


MAX_INPUT_BYTES = 16 * 1024 * 1024
MAX_DECOMPRESSED_BYTES = 64 * 1024 * 1024
MAX_JSON_BYTES = 16 * 1024 * 1024
MAX_JSON_DEPTH = 64
MAX_JSON_MEMBERS = 1_000_000
MAX_STRING_BYTES = 1 * 1024 * 1024
MAX_IDENTIFIER_BYTES = 256
MAX_SOURCE_RECORDS = 100_000
SAFE_INTEGER = 2**53 - 1
_CHUNK_BYTES = 64 * 1024


class BoundedIOError(ValueError):
    """Input or output exceeds the compiler's finite I/O contract."""


def _ceiling(value: int, label: str) -> int:
    if type(value) is not int or value < 0:
        raise BoundedIOError(f"{label} ceiling must be a non-negative integer")
    return value


def read_regular_file(
    path: str | Path,
    *,
    maximum_bytes: int = MAX_INPUT_BYTES,
    label: str = "input",
) -> bytes:
    """Read one bounded regular file without following a final-path symlink.

    The pathname is inspected before opening, the descriptor identity is checked
    after opening, and descriptor metadata is checked again after the bounded
    read.  Replacement of the pathname cannot redirect an already-open file.
    """

    limit = _ceiling(maximum_bytes, "byte")
    source = Path(path)
    descriptor = -1
    try:
        inspected = source.lstat()
        if stat.S_ISLNK(inspected.st_mode) or not stat.S_ISREG(inspected.st_mode):
            raise BoundedIOError(f"{label} must be a regular non-linked file")
        if inspected.st_size > limit:
            raise BoundedIOError(f"{label} exceeds byte limit {limit}")

        flags = os.O_RDONLY
        flags |= getattr(os, "O_BINARY", 0)
        flags |= getattr(os, "O_CLOEXEC", 0)
        flags |= getattr(os, "O_NOFOLLOW", 0)
        # Avoid blocking if a non-regular file is swapped into the pathname on
        # a platform where O_NOFOLLOW is unavailable.
        flags |= getattr(os, "O_NONBLOCK", 0)
        descriptor = os.open(source, flags)
        opened = os.fstat(descriptor)
        if not stat.S_ISREG(opened.st_mode):
            raise BoundedIOError(f"{label} must be a regular non-linked file")
        if (opened.st_dev, opened.st_ino) != (inspected.st_dev, inspected.st_ino):
            raise BoundedIOError(f"{label} changed while being opened")
        if opened.st_size > limit:
            raise BoundedIOError(f"{label} exceeds byte limit {limit}")

        chunks: list[bytes] = []
        total = 0
        while True:
            chunk = os.read(descriptor, min(_CHUNK_BYTES, limit + 1 - total))
            if not chunk:
                break
            chunks.append(chunk)
            total += len(chunk)
            if total > limit:
                raise BoundedIOError(f"{label} exceeds byte limit {limit}")

        finished = os.fstat(descriptor)
        stable_fields = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
        if any(getattr(opened, field) != getattr(finished, field) for field in stable_fields):
            raise BoundedIOError(f"{label} changed while being read")
        if total != opened.st_size:
            raise BoundedIOError(f"{label} did not match its inspected length")
        return b"".join(chunks)
    except BoundedIOError:
        raise
    except OSError as failure:
        raise BoundedIOError(f"cannot read {label}: {failure}") from failure
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def atomic_write_file(
    path: str | Path,
    raw: bytes,
    *,
    maximum_bytes: int = MAX_JSON_BYTES,
    label: str = "output",
) -> None:
    """Atomically replace a bounded regular target without following a symlink."""

    limit = _ceiling(maximum_bytes, "output byte")
    if type(raw) is not bytes:
        raise BoundedIOError(f"{label} must be bytes")
    if len(raw) > limit:
        raise BoundedIOError(f"{label} exceeds byte limit {limit}")
    destination = Path(path)
    descriptor = -1
    parent_descriptor = -1
    temporary: str | None = None
    temporary_name: str | None = None
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        descriptor_relative = (
            os.name == "posix"
            and all(
                operation in getattr(os, "supports_dir_fd", set())
                for operation in (os.open, os.stat, os.unlink)
            )
        )
        if descriptor_relative:
            parent_inspected = destination.parent.lstat()
            if stat.S_ISLNK(parent_inspected.st_mode) or not stat.S_ISDIR(
                parent_inspected.st_mode
            ):
                raise BoundedIOError(
                    f"{label} parent must be a regular non-linked directory"
                )
            parent_flags = os.O_RDONLY
            parent_flags |= getattr(os, "O_CLOEXEC", 0)
            parent_flags |= getattr(os, "O_DIRECTORY", 0)
            parent_flags |= getattr(os, "O_NOFOLLOW", 0)
            parent_descriptor = os.open(destination.parent, parent_flags)
            parent_opened = os.fstat(parent_descriptor)
            if (
                not stat.S_ISDIR(parent_opened.st_mode)
                or (parent_opened.st_dev, parent_opened.st_ino)
                != (parent_inspected.st_dev, parent_inspected.st_ino)
            ):
                raise BoundedIOError(f"{label} parent changed while opening")
            try:
                existing = os.stat(
                    destination.name,
                    dir_fd=parent_descriptor,
                    follow_symlinks=False,
                )
            except FileNotFoundError:
                existing = None
            if existing is not None and (
                stat.S_ISLNK(existing.st_mode) or not stat.S_ISREG(existing.st_mode)
            ):
                raise BoundedIOError(
                    f"{label} path must be absent or a regular non-linked file"
                )

            temporary_name = (
                f".{destination.name}.{secrets.token_hex(12)}"
            )
            flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
            flags |= getattr(os, "O_CLOEXEC", 0)
            flags |= getattr(os, "O_NOFOLLOW", 0)
            descriptor = os.open(
                temporary_name,
                flags,
                0o600,
                dir_fd=parent_descriptor,
            )
            temporary_metadata = os.fstat(descriptor)
            with os.fdopen(descriptor, "wb") as stream:
                descriptor = -1
                stream.write(raw)
                stream.flush()
                os.fsync(stream.fileno())
            try:
                current = os.stat(
                    destination.name,
                    dir_fd=parent_descriptor,
                    follow_symlinks=False,
                )
            except FileNotFoundError:
                current = None
            if current is not None and (
                stat.S_ISLNK(current.st_mode) or not stat.S_ISREG(current.st_mode)
            ):
                raise BoundedIOError(
                    f"{label} path must be absent or a regular non-linked file"
                )
            os.replace(
                temporary_name,
                destination.name,
                src_dir_fd=parent_descriptor,
                dst_dir_fd=parent_descriptor,
            )
            temporary_name = None
            published = os.stat(
                destination.name,
                dir_fd=parent_descriptor,
                follow_symlinks=False,
            )
            if (
                not stat.S_ISREG(published.st_mode)
                or (published.st_dev, published.st_ino)
                != (temporary_metadata.st_dev, temporary_metadata.st_ino)
            ):
                raise BoundedIOError(f"{label} path changed while publishing")
            parent_finished = destination.parent.lstat()
            if (
                stat.S_ISLNK(parent_finished.st_mode)
                or not stat.S_ISDIR(parent_finished.st_mode)
                or (parent_finished.st_dev, parent_finished.st_ino)
                != (parent_opened.st_dev, parent_opened.st_ino)
            ):
                raise BoundedIOError(f"{label} parent changed while publishing")
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
            return

        try:
            existing = destination.lstat()
        except FileNotFoundError:
            existing = None
        if existing is not None and (
            stat.S_ISLNK(existing.st_mode) or not stat.S_ISREG(existing.st_mode)
        ):
            raise BoundedIOError(
                f"{label} path must be absent or a regular non-linked file"
            )

        descriptor, temporary = tempfile.mkstemp(
            prefix=f".{destination.name}.", dir=destination.parent
        )
        with os.fdopen(descriptor, "wb") as stream:
            descriptor = -1
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())

        # Recheck an output that appeared or changed during the temporary write.
        try:
            current = destination.lstat()
        except FileNotFoundError:
            current = None
        if current is not None and (
            stat.S_ISLNK(current.st_mode) or not stat.S_ISREG(current.st_mode)
        ):
            raise BoundedIOError(
                f"{label} path must be absent or a regular non-linked file"
            )
        os.replace(temporary, destination)
        temporary = None
    except BoundedIOError:
        raise
    except OSError as failure:
        raise BoundedIOError(f"cannot write {label}: {failure}") from failure
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        if temporary_name is not None and parent_descriptor >= 0:
            try:
                os.unlink(temporary_name, dir_fd=parent_descriptor)
            except OSError:
                pass
        if parent_descriptor >= 0:
            os.close(parent_descriptor)
        if temporary is not None:
            try:
                os.unlink(temporary)
            except OSError:
                pass


def unwrap_gzip(
    raw: bytes,
    *,
    maximum_input_bytes: int = MAX_INPUT_BYTES,
    maximum_decompressed_bytes: int = MAX_DECOMPRESSED_BYTES,
) -> tuple[bytes, bool]:
    """Return logical bytes and whether gzip was present, with capped streaming."""

    input_limit = _ceiling(maximum_input_bytes, "input byte")
    logical_limit = _ceiling(maximum_decompressed_bytes, "decompressed byte")
    if type(raw) is not bytes:
        raise BoundedIOError("source must be bytes")
    if len(raw) > input_limit:
        raise BoundedIOError(f"source exceeds byte limit {input_limit}")
    if not raw.startswith(b"\x1f\x8b"):
        return raw, False

    chunks: list[bytes] = []
    total = 0
    try:
        with gzip.GzipFile(fileobj=io.BytesIO(raw), mode="rb") as stream:
            while True:
                chunk = stream.read(min(_CHUNK_BYTES, logical_limit + 1 - total))
                if not chunk:
                    break
                chunks.append(chunk)
                total += len(chunk)
                if total > logical_limit:
                    raise BoundedIOError(
                        f"gzip source exceeds decompressed byte limit {logical_limit}"
                    )
    except BoundedIOError:
        raise
    except (gzip.BadGzipFile, EOFError, OSError, zlib.error) as failure:
        raise BoundedIOError(f"malformed gzip source: {failure}") from failure
    return b"".join(chunks), True


def _reject_constant(value: str) -> None:
    raise BoundedIOError(f"non-finite JSON number is not allowed: {value}")


def _parse_integer(value: str) -> int:
    parsed = int(value)
    if not -SAFE_INTEGER <= parsed <= SAFE_INTEGER:
        raise BoundedIOError("JSON integer exceeds the I-JSON safe range")
    return parsed


def _parse_float(value: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed):
        raise BoundedIOError("non-finite JSON number is not allowed")
    return parsed


def _no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise BoundedIOError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def validate_json_tree(
    value: Any,
    label: str,
    *,
    maximum_depth: int = MAX_JSON_DEPTH,
    maximum_members: int = MAX_JSON_MEMBERS,
    maximum_string_bytes: int = MAX_STRING_BYTES,
) -> None:
    """Validate the same finite I-JSON tree contract consumed by compiler v2."""

    depth_limit = _ceiling(maximum_depth, "JSON depth")
    member_limit = _ceiling(maximum_members, "JSON member")
    string_limit = _ceiling(maximum_string_bytes, "JSON string byte")
    members = 0
    stack: list[tuple[Any, int]] = [(value, 0)]
    while stack:
        current, depth = stack.pop()
        if depth > depth_limit:
            raise BoundedIOError(f"{label} exceeds JSON depth limit {depth_limit}")
        if current is None or type(current) is bool:
            continue
        if type(current) is int:
            if not -SAFE_INTEGER <= current <= SAFE_INTEGER:
                raise BoundedIOError(f"{label} contains an unsafe JSON integer")
            continue
        if type(current) is float:
            if not math.isfinite(current):
                raise BoundedIOError(f"{label} contains a non-finite JSON number")
            continue
        if type(current) is str:
            if any(0xD800 <= ord(character) <= 0xDFFF for character in current):
                raise BoundedIOError(f"{label} contains a lone Unicode surrogate")
            if len(current.encode("utf-8")) > string_limit:
                raise BoundedIOError(
                    f"{label} exceeds JSON string byte limit {string_limit}"
                )
            continue
        if type(current) is list:
            members += len(current)
            stack.extend((item, depth + 1) for item in reversed(current))
        elif type(current) is dict:
            members += len(current)
            for key, item in reversed(list(current.items())):
                if type(key) is not str:
                    raise BoundedIOError(f"{label} contains a non-string object key")
                stack.append((key, depth + 1))
                stack.append((item, depth + 1))
        else:
            raise BoundedIOError(
                f"{label} contains unsupported type {type(current).__name__}"
            )
        if members > member_limit:
            raise BoundedIOError(f"{label} exceeds JSON member limit {member_limit}")


def load_json_object(
    raw: bytes,
    label: str,
    *,
    maximum_bytes: int = MAX_JSON_BYTES,
    maximum_depth: int = MAX_JSON_DEPTH,
    maximum_members: int = MAX_JSON_MEMBERS,
    maximum_string_bytes: int = MAX_STRING_BYTES,
) -> dict[str, Any]:
    """Decode one bounded, duplicate-safe I-JSON object."""

    byte_limit = _ceiling(maximum_bytes, "JSON byte")
    if type(raw) is not bytes:
        raise BoundedIOError(f"{label} JSON input must be bytes")
    if len(raw) > byte_limit:
        raise BoundedIOError(f"{label} exceeds JSON byte limit {byte_limit}")
    try:
        value = json.loads(
            raw.decode("utf-8"),
            object_pairs_hook=_no_duplicates,
            parse_constant=_reject_constant,
            parse_int=_parse_integer,
            parse_float=_parse_float,
        )
    except BoundedIOError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError, RecursionError) as failure:
        raise BoundedIOError(f"invalid {label} JSON: {failure}") from failure
    if type(value) is not dict:
        raise BoundedIOError(f"{label} must be a JSON object")
    validate_json_tree(
        value,
        label,
        maximum_depth=maximum_depth,
        maximum_members=maximum_members,
        maximum_string_bytes=maximum_string_bytes,
    )
    return value


def pretty_json_bytes(
    payload: dict[str, Any],
    *,
    ensure_ascii: bool,
    maximum_bytes: int = MAX_JSON_BYTES,
    maximum_depth: int = MAX_JSON_DEPTH,
    maximum_members: int = MAX_JSON_MEMBERS,
    maximum_string_bytes: int = MAX_STRING_BYTES,
) -> bytes:
    """Serialize exactly as the source producers do, after v2 admission checks."""

    byte_limit = _ceiling(maximum_bytes, "JSON byte")
    validate_json_tree(
        payload,
        "output artifact",
        maximum_depth=maximum_depth,
        maximum_members=maximum_members,
        maximum_string_bytes=maximum_string_bytes,
    )
    try:
        raw = (
            json.dumps(
                payload,
                indent=2,
                sort_keys=True,
                ensure_ascii=ensure_ascii,
                allow_nan=False,
            )
            + "\n"
        ).encode("utf-8")
    except (TypeError, ValueError, UnicodeError) as failure:
        raise BoundedIOError(f"output is not I-JSON: {failure}") from failure
    if len(raw) > byte_limit:
        raise BoundedIOError(f"serialized output exceeds JSON byte limit {byte_limit}")
    return raw


def compact_json_bytes(
    payload: dict[str, Any],
    *,
    ensure_ascii: bool,
    maximum_bytes: int = MAX_JSON_BYTES,
    maximum_depth: int = MAX_JSON_DEPTH,
    maximum_members: int = MAX_JSON_MEMBERS,
    maximum_string_bytes: int = MAX_STRING_BYTES,
) -> bytes:
    """Serialize bounded deterministic JSON without presentation whitespace."""

    byte_limit = _ceiling(maximum_bytes, "JSON byte")
    validate_json_tree(
        payload,
        "output artifact",
        maximum_depth=maximum_depth,
        maximum_members=maximum_members,
        maximum_string_bytes=maximum_string_bytes,
    )
    try:
        raw = (
            json.dumps(
                payload,
                separators=(",", ":"),
                sort_keys=True,
                ensure_ascii=ensure_ascii,
                allow_nan=False,
            )
            + "\n"
        ).encode("utf-8")
    except (TypeError, ValueError, UnicodeError) as failure:
        raise BoundedIOError(f"output is not I-JSON: {failure}") from failure
    if len(raw) > byte_limit:
        raise BoundedIOError(f"serialized output exceeds JSON byte limit {byte_limit}")
    return raw


__all__ = [
    "atomic_write_file",
    "BoundedIOError",
    "compact_json_bytes",
    "MAX_DECOMPRESSED_BYTES",
    "MAX_IDENTIFIER_BYTES",
    "MAX_INPUT_BYTES",
    "MAX_JSON_BYTES",
    "MAX_JSON_DEPTH",
    "MAX_JSON_MEMBERS",
    "MAX_SOURCE_RECORDS",
    "MAX_STRING_BYTES",
    "load_json_object",
    "pretty_json_bytes",
    "read_regular_file",
    "unwrap_gzip",
    "validate_json_tree",
]

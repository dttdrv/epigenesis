"""Independent replay validator for streaming reference FASTA catalogs.

This module intentionally imports only the Python standard library.  The
producer, compiler, and their canonicalization and I/O helpers are outside its
trust boundary.
"""

from __future__ import annotations

import base64
import copy
import functools
import hashlib
import json
import math
import os
from pathlib import Path
import re
import stat
from typing import Any, BinaryIO, Iterator
import zlib


SAFE_INTEGER = 2**53 - 1
MAX_IDENTIFIER_BYTES = 256
MAX_STRING_BYTES = 1024 * 1024
MAX_ARTIFACT_BYTES = 256 * 1024 * 1024
MAX_JSON_DEPTH = 64
MAX_JSON_MEMBERS = 1_000_000
MAX_RECORD_BASES = 2_147_483_647
MAX_RECORDS = 100_000
DEFAULT_MAXIMUM_GZIP_MEMBERS = 100_000
MAXIMUM_GZIP_MEMBERS = 1_000_000
DEFAULT_MAX_INPUT_BYTES = 64 * 1024**3
DEFAULT_MAX_LOGICAL_BYTES = 64 * 1024**3
CHUNK_BYTES = 64 * 1024

FORMAT = "brainc.reference-sequence-catalog"
VERSION = 1
PROFILE = "fasta-reference-dna/v1"
PRODUCER = {
    "name": "brainc-reference-fasta",
    "version": "1.0.0",
    "passes": [
        "stream-fasta",
        "validate-iupac-dna",
        "digest-reference-sequences",
        "emit-reference-catalog",
    ],
}
REPORT_FORMAT = "brainc.reference-sequence-catalog-validation-report"
REPORT_CHECKS = (
    "closed-reference-catalog-schema",
    "canonical-catalog-and-artifact-seals",
    "stable-single-link-source-path",
    "exact-physical-source-identity",
    "exact-logical-fasta-identity",
    "independent-streaming-fasta-replay",
    "sequence-digests-and-refget-identifiers",
    "sequence-collection-identity",
)

_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_REFGET = re.compile(r"SQ\.[A-Za-z0-9_-]{32}\Z")
_GA4GH_DIGEST = re.compile(r"[A-Za-z0-9_-]{32}\Z")
_INVALID_SEQUENCE = re.compile(rb"[^ACGTRYSWKMBDHVNacgtryswkmbdhvn\r\n]")
_UPPER = bytes.maketrans(
    b"acgtryswkmbdhvn",
    b"ACGTRYSWKMBDHVN",
)


class ReferenceValidationError(ValueError):
    """A reference catalog differs from independent source replay."""


def _fail(detail: str) -> ReferenceValidationError:
    return ReferenceValidationError(f"REFVAL001: {detail}")


def _bounded(function: Any) -> Any:
    @functools.wraps(function)
    def guarded(*args: Any, **kwargs: Any) -> Any:
        try:
            return function(*args, **kwargs)
        except ReferenceValidationError:
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


def _jcs_chunks(value: Any) -> Iterator[str]:
    if value is None:
        yield "null"
    elif type(value) is bool:
        yield "true" if value else "false"
    elif type(value) is str:
        yield _jcs_string(value)
    elif type(value) is int:
        if not -SAFE_INTEGER <= value <= SAFE_INTEGER:
            raise ValueError("unsafe integer")
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
            raise TypeError("non-string object key")
        yield "{"
        keys = sorted(value, key=lambda key: key.encode("utf-16be"))
        for index, key in enumerate(keys):
            if index:
                yield ","
            yield _jcs_string(key)
            yield ":"
            yield from _jcs_chunks(value[key])
        yield "}"
    else:
        raise TypeError(f"unsupported value type {type(value).__name__}")


def _canonical_hash(value: Any, hasher: Any) -> bytes:
    try:
        for chunk in _jcs_chunks(value):
            hasher.update(chunk.encode("utf-8"))
    except ReferenceValidationError:
        raise
    except MemoryError as failure:
        raise _fail("canonical JSON exceeds the validation memory ceiling") from failure
    except (TypeError, ValueError, UnicodeError, RecursionError) as failure:
        raise _fail(f"value is not RFC 8785 canonical JSON: {failure}") from failure
    return hasher.digest()


def digest(value: Any) -> str:
    """Return the independently computed RFC 8785 SHA-256 digest."""

    return _canonical_hash(value, hashlib.sha256()).hex()


def _sha512t24u(value: Any) -> str:
    raw = _canonical_hash(value, hashlib.sha512())[:24]
    return base64.urlsafe_b64encode(raw).decode("ascii")


def _keys(value: Any, expected: set[str], label: str) -> dict[str, Any]:
    if type(value) is not dict:
        raise _fail(f"{label} must be an object")
    if any(type(key) is not str for key in value):
        raise _fail(f"{label} keys must be strings")
    missing = sorted(expected - value.keys())
    extra = sorted(value.keys() - expected)
    if missing or extra:
        raise _fail(
            f"{label} keys invalid; missing={missing or 'none'}, "
            f"unknown={extra or 'none'}"
        )
    return value


def _integer(
    value: Any,
    label: str,
    *,
    minimum: int = 1,
    maximum: int = SAFE_INTEGER,
) -> int:
    if type(value) is not int or not minimum <= value <= maximum:
        raise _fail(f"{label} must be an integer in [{minimum}, {maximum}]")
    return value


def _sha(value: Any, label: str) -> str:
    if type(value) is not str or _SHA256.fullmatch(value) is None:
        raise _fail(f"{label} must be a lowercase SHA-256 digest")
    return value


def _sequence_collection(records: list[dict[str, Any]]) -> dict[str, Any]:
    level_2 = {
        "lengths": [record["bases"] for record in records],
        "names": [record["record_id"] for record in records],
        "sequences": [record["refget_id"] for record in records],
    }
    level_1 = {
        name: _sha512t24u(value)
        for name, value in level_2.items()
    }
    inherent = {name: level_1[name] for name in ("names", "sequences")}
    return {
        "version": "1.0.0",
        "digest": _sha512t24u(inherent),
        "level_1": level_1,
        "level_2": level_2,
    }


def _validate_sequence_collection(
    value: Any,
    records: list[dict[str, Any]],
) -> None:
    root = _keys(
        value,
        {"version", "digest", "level_1", "level_2"},
        "sequence_catalog.refget_seqcol",
    )
    if type(root["version"]) is not str or root["version"] != "1.0.0":
        raise _fail("sequence_catalog.refget_seqcol.version is unsupported")
    if (
        type(root["digest"]) is not str
        or _GA4GH_DIGEST.fullmatch(root["digest"]) is None
    ):
        raise _fail("sequence_catalog.refget_seqcol.digest is invalid")
    level_1 = _keys(
        root["level_1"],
        {"lengths", "names", "sequences"},
        "sequence_catalog.refget_seqcol.level_1",
    )
    for name in ("lengths", "names", "sequences"):
        if (
            type(level_1[name]) is not str
            or _GA4GH_DIGEST.fullmatch(level_1[name]) is None
        ):
            raise _fail(
                f"sequence_catalog.refget_seqcol.level_1.{name} is invalid"
            )
    level_2 = _keys(
        root["level_2"],
        {"lengths", "names", "sequences"},
        "sequence_catalog.refget_seqcol.level_2",
    )
    if type(level_2["lengths"]) is not list or any(
        type(item) is not int for item in level_2["lengths"]
    ):
        raise _fail("sequence_catalog.refget_seqcol.level_2.lengths is invalid")
    if type(level_2["names"]) is not list or any(
        type(item) is not str for item in level_2["names"]
    ):
        raise _fail("sequence_catalog.refget_seqcol.level_2.names is invalid")
    if type(level_2["sequences"]) is not list or any(
        type(item) is not str for item in level_2["sequences"]
    ):
        raise _fail("sequence_catalog.refget_seqcol.level_2.sequences is invalid")
    if root != _sequence_collection(records):
        raise _fail("sequence_catalog.refget_seqcol does not match records")


@_bounded
def validate_reference_catalog(payload: dict[str, Any]) -> dict[str, Any]:
    """Validate the closed catalog schema and both canonical seals."""

    root = _keys(
        payload,
        {
            "format",
            "version",
            "producer",
            "inputs",
            "sequence_catalog",
            "catalog_sha256",
            "artifact_sha256",
        },
        "reference sequence catalog",
    )
    if (
        type(root["format"]) is not str
        or root["format"] != FORMAT
        or type(root["version"]) is not int
        or root["version"] != VERSION
    ):
        raise _fail("unsupported reference sequence catalog identity")
    producer = _keys(
        root["producer"],
        {"name", "version", "passes"},
        "producer",
    )
    if producer != PRODUCER:
        raise _fail("reference sequence catalog producer is unsupported")
    inputs = _keys(root["inputs"], {"profile", "wrapper", "source", "logical"}, "inputs")
    if type(inputs["profile"]) is not str or inputs["profile"] != PROFILE:
        raise _fail("reference sequence catalog profile is unsupported")
    if type(inputs["wrapper"]) is not str or inputs["wrapper"] not in {
        "identity",
        "gzip",
    }:
        raise _fail("reference sequence catalog wrapper is unsupported")
    for name in ("source", "logical"):
        reference = _keys(inputs[name], {"sha256", "byte_length"}, f"inputs.{name}")
        _sha(reference["sha256"], f"inputs.{name}.sha256")
        _integer(reference["byte_length"], f"inputs.{name}.byte_length")
    if inputs["wrapper"] == "identity" and inputs["source"] != inputs["logical"]:
        raise _fail("identity wrapper source and logical identities must match")

    catalog = _keys(
        root["sequence_catalog"],
        {"records", "total_bases", "refget_seqcol"},
        "sequence_catalog",
    )
    records = catalog["records"]
    if type(records) is not list or not 1 <= len(records) <= MAX_RECORDS:
        raise _fail(f"sequence_catalog.records must contain 1..{MAX_RECORDS} records")
    seen: set[str] = set()
    for index, value in enumerate(records):
        record = _keys(
            value,
            {"record_id", "bases", "sequence_sha256", "refget_id"},
            f"sequence_catalog.records[{index}]",
        )
        record_id = record["record_id"]
        try:
            encoded = record_id.encode("utf-8") if type(record_id) is str else b""
        except UnicodeEncodeError as failure:
            raise _fail(
                f"sequence_catalog.records[{index}].record_id is invalid"
            ) from failure
        if (
            type(record_id) is not str
            or not record_id
            or record_id != record_id.strip()
            or any(character.isspace() for character in record_id)
            or any(
                ord(character) < 0x20 or 0x7F <= ord(character) <= 0x9F
                for character in record_id
            )
            or len(encoded) > MAX_IDENTIFIER_BYTES
        ):
            raise _fail(f"sequence_catalog.records[{index}].record_id is invalid")
        if record_id in seen:
            raise _fail("sequence_catalog.records contains duplicate identifiers")
        seen.add(record_id)
        _integer(
            record["bases"],
            f"sequence_catalog.records[{index}].bases",
            maximum=MAX_RECORD_BASES,
        )
        _sha(
            record["sequence_sha256"],
            f"sequence_catalog.records[{index}].sequence_sha256",
        )
        if (
            type(record["refget_id"]) is not str
            or _REFGET.fullmatch(record["refget_id"]) is None
        ):
            raise _fail(f"sequence_catalog.records[{index}].refget_id is invalid")

    total_bases = sum(record["bases"] for record in records)
    _integer(catalog["total_bases"], "sequence_catalog.total_bases")
    if total_bases > SAFE_INTEGER or catalog["total_bases"] != total_bases:
        raise _fail("sequence_catalog.total_bases does not match records")
    _validate_sequence_collection(catalog["refget_seqcol"], records)
    if _sha(root["catalog_sha256"], "catalog_sha256") != digest(catalog):
        raise _fail("catalog_sha256 does not match sequence_catalog")
    claimed = _sha(root["artifact_sha256"], "artifact_sha256")
    if claimed != digest(
        {key: value for key, value in root.items() if key != "artifact_sha256"}
    ):
        raise _fail("artifact_sha256 does not match canonical artifact content")
    return copy.deepcopy(root)


class _TrackingReader:
    __slots__ = ("stream", "maximum", "count", "hasher")

    def __init__(self, stream: BinaryIO, maximum: int) -> None:
        self.stream = stream
        self.maximum = maximum
        self.count = 0
        self.hasher = hashlib.sha256()

    def read(self, size: int = -1) -> bytes:
        raw = self.stream.read(size)
        self.count += len(raw)
        if self.count > self.maximum:
            raise _fail(f"physical source exceeds byte ceiling {self.maximum}")
        self.hasher.update(raw)
        return raw


class _StrictGzipReader:
    """Stream only complete RFC 1952 members, with no padding or junk."""

    __slots__ = (
        "source",
        "maximum_members",
        "members",
        "pending",
        "inflater",
        "finished",
    )

    def __init__(self, source: _TrackingReader, maximum_members: int) -> None:
        self.source = source
        self.maximum_members = maximum_members
        self.members = 0
        self.pending = b""
        self.inflater: Any = None
        self.finished = False

    def _more(self) -> bool:
        raw = self.source.read(CHUNK_BYTES)
        if raw:
            self.pending += raw
            return True
        return False

    def _start_member(self) -> None:
        while len(self.pending) < 2 and self._more():
            pass
        if not self.pending:
            if self.members == 0:
                raise _fail("gzip FASTA contains no members")
            self.finished = True
            return
        if len(self.pending) < 2 or self.pending[:2] != b"\x1f\x8b":
            raise _fail("gzip FASTA contains bytes outside complete members")
        if self.members >= self.maximum_members:
            raise _fail(f"gzip FASTA exceeds member ceiling {self.maximum_members}")
        self.members += 1
        self.inflater = zlib.decompressobj(wbits=16 + zlib.MAX_WBITS)

    def read(self, size: int = -1) -> bytes:
        if self.finished:
            return b""
        if size is None or size < 0:
            size = CHUNK_BYTES
        if size == 0:
            return b""
        output = bytearray()
        while len(output) < size and not self.finished:
            if self.inflater is None:
                self._start_member()
                if self.finished:
                    break
            if not self.pending and not self._more():
                raise _fail("gzip FASTA contains an incomplete member")
            assert self.inflater is not None
            try:
                produced = self.inflater.decompress(
                    self.pending,
                    size - len(output),
                )
            except zlib.error as failure:
                raise _fail(f"invalid gzip FASTA member: {failure}") from failure
            output.extend(produced)
            if self.inflater.eof:
                self.pending = self.inflater.unused_data
                self.inflater = None
            else:
                self.pending = self.inflater.unconsumed_tail
                if not produced and not self.pending and not self._more():
                    raise _fail("gzip FASTA contains an incomplete member")
        return bytes(output)


class _FastaReplay:
    __slots__ = (
        "maximum_logical_bytes",
        "logical_count",
        "logical_hasher",
        "header",
        "at_line_start",
        "record_id",
        "record_bases",
        "record_sha256",
        "record_sha512",
        "total_bases",
        "records",
        "seen",
    )

    def __init__(self, maximum_logical_bytes: int) -> None:
        self.maximum_logical_bytes = maximum_logical_bytes
        self.logical_count = 0
        self.logical_hasher = hashlib.sha256()
        self.header: bytearray | None = None
        self.at_line_start = True
        self.record_id: str | None = None
        self.record_bases = 0
        self.record_sha256 = hashlib.sha256()
        self.record_sha512 = hashlib.sha512()
        self.total_bases = 0
        self.records: list[dict[str, Any]] = []
        self.seen: set[str] = set()

    def _finish_record(self) -> None:
        if self.record_id is None:
            return
        if self.record_bases == 0:
            raise _fail(f"FASTA record {self.record_id!r} has no sequence")
        self.total_bases += self.record_bases
        if self.total_bases > SAFE_INTEGER:
            raise _fail("FASTA total base count exceeds I-JSON safe integer range")
        refget = "SQ." + base64.urlsafe_b64encode(
            self.record_sha512.digest()[:24]
        ).decode("ascii")
        self.records.append(
            {
                "record_id": self.record_id,
                "bases": self.record_bases,
                "sequence_sha256": self.record_sha256.hexdigest(),
                "refget_id": refget,
            }
        )
        self.record_id = None

    def _start_header(self) -> None:
        self._finish_record()
        self.header = bytearray()
        self.at_line_start = False

    def _append_header(self, raw: bytes) -> None:
        assert self.header is not None
        if len(self.header) + len(raw) > MAX_STRING_BYTES:
            raise _fail(f"FASTA defline exceeds byte ceiling {MAX_STRING_BYTES}")
        self.header.extend(raw)

    def _finish_header(self) -> None:
        assert self.header is not None
        try:
            defline = bytes(self.header).decode("utf-8")
        except UnicodeDecodeError as failure:
            raise _fail("FASTA defline must be UTF-8") from failure
        self.header = None
        if not defline.strip():
            raise _fail("FASTA record identifier is empty")
        fields = defline.split(maxsplit=1)
        record_id = fields[0]
        description = fields[1] if len(fields) == 2 else ""
        if len(record_id.encode("utf-8")) > MAX_IDENTIFIER_BYTES:
            raise _fail(
                f"FASTA record identifier exceeds {MAX_IDENTIFIER_BYTES} UTF-8 bytes"
            )
        if any(character.isspace() for character in record_id) or any(
            ord(character) < 0x20 or 0x7F <= ord(character) <= 0x9F
            for character in record_id
        ):
            raise _fail("FASTA record identifier is invalid")
        if len(description.encode("utf-8")) > MAX_STRING_BYTES or any(
            ord(character) < 0x20 or 0x7F <= ord(character) <= 0x9F
            for character in description
        ):
            raise _fail("FASTA description is invalid")
        if record_id in self.seen:
            raise _fail(f"duplicate FASTA record identifier: {record_id}")
        if len(self.records) >= MAX_RECORDS:
            raise _fail(f"FASTA exceeds record ceiling {MAX_RECORDS}")
        self.seen.add(record_id)
        self.record_id = record_id
        self.record_bases = 0
        self.record_sha256 = hashlib.sha256()
        self.record_sha512 = hashlib.sha512()
        self.at_line_start = True

    def _sequence(self, raw: bytes) -> None:
        if self.record_id is None:
            raise _fail("FASTA must begin with a '>' defline")
        invalid = _INVALID_SEQUENCE.search(raw)
        if invalid is not None:
            raise _fail(f"invalid IUPAC DNA byte 0x{raw[invalid.start()]:02x}")
        canonical = raw.translate(_UPPER, b"\r\n")
        bases = self.record_bases + len(canonical)
        if bases > MAX_RECORD_BASES:
            raise _fail(f"FASTA record exceeds base ceiling {MAX_RECORD_BASES}")
        self.record_bases = bases
        self.record_sha256.update(canonical)
        self.record_sha512.update(canonical)
        if raw:
            self.at_line_start = raw.endswith((b"\r", b"\n"))

    @staticmethod
    def _line_break(raw: bytes, start: int) -> tuple[int, int] | None:
        carriage = raw.find(b"\r", start)
        newline = raw.find(b"\n", start)
        if carriage < 0 and newline < 0:
            return None
        if newline < 0 or 0 <= carriage < newline:
            width = 2 if carriage + 1 < len(raw) and raw[carriage + 1] == 0x0A else 1
            return carriage, width
        return newline, 1

    @staticmethod
    def _header_marker(raw: bytes, start: int) -> int:
        after_lf = raw.find(b"\n>", start)
        after_cr = raw.find(b"\r>", start)
        if after_lf < 0:
            return after_cr
        if after_cr < 0:
            return after_lf
        return min(after_lf, after_cr)

    def _parse(self, raw: bytes) -> None:
        position = 0
        while position < len(raw):
            if self.header is not None:
                boundary = self._line_break(raw, position)
                if boundary is None:
                    self._append_header(raw[position:])
                    return
                end, width = boundary
                self._append_header(raw[position:end])
                self._finish_header()
                position = end + width
                continue
            if self.at_line_start and raw[position] == 0x3E:
                self._start_header()
                position += 1
                continue
            marker = self._header_marker(raw, position)
            if marker < 0:
                self._sequence(raw[position:])
                return
            self._sequence(raw[position : marker + 1])
            self._start_header()
            position = marker + 2

    def scan(self, stream: Any) -> tuple[list[dict[str, Any]], int, str, int]:
        carry = b""
        while True:
            chunk = stream.read(CHUNK_BYTES)
            if not chunk:
                break
            self.logical_count += len(chunk)
            if self.logical_count > self.maximum_logical_bytes:
                raise _fail(
                    f"logical FASTA exceeds byte ceiling {self.maximum_logical_bytes}"
                )
            self.logical_hasher.update(chunk)
            parse = carry + chunk
            carry = b"\r" if parse.endswith(b"\r") else b""
            if carry:
                parse = parse[:-1]
            if parse:
                self._parse(parse)
        if carry:
            self._parse(carry)
        if self.header is not None:
            self._finish_header()
        self._finish_record()
        if not self.records:
            raise _fail("FASTA contains no records")
        return (
            self.records,
            self.total_bases,
            self.logical_hasher.hexdigest(),
            self.logical_count,
        )


def _limit(value: Any, label: str) -> int:
    return _integer(value, label, maximum=SAFE_INTEGER)


def _open_regular(
    path: str | Path,
    *,
    label: str,
    maximum_bytes: int,
) -> tuple[BinaryIO, os.stat_result, Path]:
    descriptor = -1
    try:
        if type(path) is not str and not isinstance(path, Path):
            raise _fail(f"{label} path must be str or Path")
        source = Path(path)
        inspected = source.lstat()
        if (
            stat.S_ISLNK(inspected.st_mode)
            or not stat.S_ISREG(inspected.st_mode)
            or inspected.st_nlink != 1
        ):
            raise _fail(f"{label} must be a regular non-linked single-link file")
        if inspected.st_size > maximum_bytes:
            raise _fail(f"{label} exceeds byte ceiling {maximum_bytes}")
        flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
        flags |= getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
        descriptor = os.open(source, flags)
        opened = os.fstat(descriptor)
        if (
            not stat.S_ISREG(opened.st_mode)
            or opened.st_nlink != 1
            or (opened.st_dev, opened.st_ino) != (inspected.st_dev, inspected.st_ino)
        ):
            raise _fail(f"{label} changed while opening")
        if opened.st_size > maximum_bytes:
            raise _fail(f"{label} exceeds byte ceiling {maximum_bytes}")
        stream = os.fdopen(descriptor, "rb", buffering=0)
        descriptor = -1
        return stream, opened, source
    except ReferenceValidationError:
        raise
    except (OSError, TypeError, ValueError) as failure:
        raise _fail(f"cannot open {label}: {failure}") from failure
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _finish_regular(
    source: Path,
    stream: BinaryIO,
    opened: os.stat_result,
    byte_count: int,
    label: str,
) -> None:
    try:
        finished = os.fstat(stream.fileno())
        current = source.lstat()
    except OSError as failure:
        raise _fail(f"{label} changed while reading: {failure}") from failure
    fields = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns", "st_nlink")
    if any(getattr(opened, field) != getattr(finished, field) for field in fields):
        raise _fail(f"{label} changed while reading")
    if (
        stat.S_ISLNK(current.st_mode)
        or not stat.S_ISREG(current.st_mode)
        or current.st_nlink != 1
        or any(getattr(opened, field) != getattr(current, field) for field in fields)
    ):
        raise _fail(f"{label} path changed while reading")
    if byte_count != opened.st_size:
        raise _fail(f"{label} byte count changed while reading")


def _read_regular(
    path: str | Path,
    *,
    label: str,
    maximum_bytes: int,
) -> bytes:
    stream, opened, source = _open_regular(
        path,
        label=label,
        maximum_bytes=maximum_bytes,
    )
    chunks: list[bytes] = []
    count = 0
    try:
        while True:
            raw = stream.read(min(CHUNK_BYTES, maximum_bytes + 1 - count))
            if not raw:
                break
            chunks.append(raw)
            count += len(raw)
            if count > maximum_bytes:
                raise _fail(f"{label} exceeds byte ceiling {maximum_bytes}")
        _finish_regular(source, stream, opened, count, label)
        return b"".join(chunks)
    except ReferenceValidationError:
        raise
    except OSError as failure:
        raise _fail(f"cannot read {label}: {failure}") from failure
    finally:
        stream.close()


def _json_preflight(raw: bytes, label: str) -> None:
    """Bound JSON shape before the standard-library decoder allocates it."""

    stack: list[list[Any]] = []
    members = 0
    in_string = False
    escaped = False
    string_start = 0

    def add_member() -> None:
        nonlocal members
        members += 1
        if members > MAX_JSON_MEMBERS:
            raise _fail(f"{label} exceeds JSON member ceiling {MAX_JSON_MEMBERS}")

    for index, byte in enumerate(raw):
        if in_string:
            if escaped:
                escaped = False
            elif byte == 0x5C:
                escaped = True
            elif byte == 0x22:
                in_string = False
            if index - string_start > MAX_STRING_BYTES:
                raise _fail(
                    f"{label} exceeds JSON string byte ceiling {MAX_STRING_BYTES}"
                )
            continue
        if byte in b" \t\r\n":
            continue
        if byte == 0x22:
            if stack and stack[-1] == ["array", True]:
                add_member()
                stack[-1][1] = False
            in_string = True
            string_start = index + 1
        elif byte in (0x7B, 0x5B):
            if stack and stack[-1] == ["array", True]:
                add_member()
                stack[-1][1] = False
            stack.append(["object" if byte == 0x7B else "array", byte == 0x5B])
            if len(stack) - 1 > MAX_JSON_DEPTH:
                raise _fail(f"{label} exceeds JSON depth ceiling {MAX_JSON_DEPTH}")
        elif byte == 0x3A and stack and stack[-1][0] == "object":
            add_member()
        elif byte == 0x2C and stack and stack[-1][0] == "array":
            stack[-1][1] = True
        elif byte in (0x7D, 0x5D):
            if stack:
                stack.pop()
        elif stack and stack[-1] == ["array", True]:
            add_member()
            stack[-1][1] = False


def _validate_json_tree(value: Any, label: str) -> None:
    members = 0
    pending: list[tuple[Any, int]] = [(value, 0)]
    while pending:
        current, depth = pending.pop()
        if depth > MAX_JSON_DEPTH:
            raise _fail(f"{label} exceeds JSON depth ceiling {MAX_JSON_DEPTH}")
        if current is None or type(current) is bool:
            continue
        if type(current) is int:
            if not -SAFE_INTEGER <= current <= SAFE_INTEGER:
                raise _fail(f"{label} contains an unsafe JSON integer")
            continue
        if type(current) is float:
            if not math.isfinite(current):
                raise _fail(f"{label} contains a non-finite JSON number")
            continue
        if type(current) is str:
            try:
                encoded = current.encode("utf-8")
            except UnicodeEncodeError as failure:
                raise _fail(f"{label} contains a lone Unicode surrogate") from failure
            if len(encoded) > MAX_STRING_BYTES:
                raise _fail(
                    f"{label} exceeds JSON string byte ceiling {MAX_STRING_BYTES}"
                )
            continue
        if type(current) is list:
            members += len(current)
            pending.extend((item, depth + 1) for item in current)
        elif type(current) is dict:
            members += len(current)
            for key, item in current.items():
                if type(key) is not str:
                    raise _fail(f"{label} contains a non-string JSON key")
                pending.append((key, depth + 1))
                pending.append((item, depth + 1))
        else:
            raise _fail(f"{label} contains unsupported JSON value")
        if members > MAX_JSON_MEMBERS:
            raise _fail(f"{label} exceeds JSON member ceiling {MAX_JSON_MEMBERS}")


def _load_json(raw: bytes, label: str) -> dict[str, Any]:
    def pairs(items: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in items:
            if key in result:
                raise _fail(f"{label} contains duplicate JSON key {key!r}")
            result[key] = value
        return result

    def constant(token: str) -> None:
        raise _fail(f"{label} contains non-finite number {token}")

    try:
        _json_preflight(raw, label)
        value = json.loads(
            raw.decode("utf-8"),
            object_pairs_hook=pairs,
            parse_constant=constant,
        )
    except ReferenceValidationError:
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
    _validate_json_tree(value, label)
    return value


def _observed_artifact(
    path: str | Path,
    wrapper: str,
    *,
    maximum_source_bytes: int,
    maximum_logical_bytes: int,
    maximum_gzip_members: int,
) -> dict[str, Any]:
    stream, opened, source = _open_regular(
        path,
        label="reference FASTA",
        maximum_bytes=maximum_source_bytes,
    )
    physical = _TrackingReader(stream, maximum_source_bytes)
    try:
        logical: Any = physical
        if wrapper == "gzip":
            logical = _StrictGzipReader(physical, maximum_gzip_members)
        records, total_bases, logical_sha256, logical_bytes = _FastaReplay(
            maximum_logical_bytes
        ).scan(logical)
        while physical.read(CHUNK_BYTES):
            pass
        _finish_regular(source, stream, opened, physical.count, "reference FASTA")
    except ReferenceValidationError:
        raise
    except (EOFError, OSError, zlib.error) as failure:
        raise _fail(f"invalid {wrapper} FASTA source: {failure}") from failure
    finally:
        stream.close()

    catalog = {
        "records": records,
        "total_bases": total_bases,
        "refget_seqcol": _sequence_collection(records),
    }
    core = {
        "format": FORMAT,
        "version": VERSION,
        "producer": {
            "name": PRODUCER["name"],
            "version": PRODUCER["version"],
            "passes": list(PRODUCER["passes"]),
        },
        "inputs": {
            "profile": PROFILE,
            "wrapper": wrapper,
            "source": {
                "sha256": physical.hasher.hexdigest(),
                "byte_length": physical.count,
            },
            "logical": {
                "sha256": logical_sha256,
                "byte_length": logical_bytes,
            },
        },
        "sequence_catalog": catalog,
        "catalog_sha256": digest(catalog),
    }
    return {**core, "artifact_sha256": digest(core)}


def _validation_report(artifact: dict[str, Any]) -> dict[str, Any]:
    inputs = artifact["inputs"]
    catalog = artifact["sequence_catalog"]
    core = {
        "format": REPORT_FORMAT,
        "version": 1,
        "stage": "reference-fasta-admission",
        "valid": True,
        "inputs": {
            "artifact_sha256": artifact["artifact_sha256"],
            "source_sha256": inputs["source"]["sha256"],
            "logical_sha256": inputs["logical"]["sha256"],
        },
        "result": {
            "records": len(catalog["records"]),
            "total_bases": catalog["total_bases"],
            "source_bytes": inputs["source"]["byte_length"],
            "logical_bytes": inputs["logical"]["byte_length"],
            "catalog_sha256": artifact["catalog_sha256"],
            "refget_seqcol_digest": catalog["refget_seqcol"]["digest"],
        },
        "checks": list(REPORT_CHECKS),
    }
    return {**core, "report_sha256": digest(core)}


@_bounded
def validate_reference_report(report: dict[str, Any]) -> dict[str, Any]:
    """Validate the exact deterministic success-report contract and seal."""

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
    if (
        type(root["format"]) is not str
        or root["format"] != REPORT_FORMAT
        or type(root["version"]) is not int
        or root["version"] != 1
        or type(root["stage"]) is not str
        or root["stage"] != "reference-fasta-admission"
        or type(root["valid"]) is not bool
        or root["valid"] is not True
    ):
        raise _fail("validation report identity is invalid")
    inputs = _keys(
        root["inputs"],
        {"artifact_sha256", "source_sha256", "logical_sha256"},
        "validation report.inputs",
    )
    for name, value in inputs.items():
        _sha(value, f"validation report.inputs.{name}")
    result = _keys(
        root["result"],
        {
            "records",
            "total_bases",
            "source_bytes",
            "logical_bytes",
            "catalog_sha256",
            "refget_seqcol_digest",
        },
        "validation report.result",
    )
    _integer(
        result["records"],
        "validation report.result.records",
        maximum=MAX_RECORDS,
    )
    for name in ("total_bases", "source_bytes", "logical_bytes"):
        _integer(result[name], f"validation report.result.{name}")
    _sha(result["catalog_sha256"], "validation report.result.catalog_sha256")
    if (
        type(result["refget_seqcol_digest"]) is not str
        or _GA4GH_DIGEST.fullmatch(result["refget_seqcol_digest"]) is None
    ):
        raise _fail("validation report.result.refget_seqcol_digest is invalid")
    if root["checks"] != list(REPORT_CHECKS):
        raise _fail("validation report checks differ from the closed contract")
    if _sha(root["report_sha256"], "validation report.report_sha256") != digest(
        {key: value for key, value in root.items() if key != "report_sha256"}
    ):
        raise _fail("validation report seal is invalid")
    return copy.deepcopy(root)


@_bounded
def validate_reference_fasta(
    source_path: str | Path,
    payload: dict[str, Any],
    *,
    maximum_source_bytes: int = DEFAULT_MAX_INPUT_BYTES,
    maximum_logical_bytes: int = DEFAULT_MAX_LOGICAL_BYTES,
    maximum_gzip_members: int = DEFAULT_MAXIMUM_GZIP_MEMBERS,
) -> dict[str, Any]:
    """Replay one source path and return a sealed deterministic report."""

    source_limit = _limit(maximum_source_bytes, "maximum_source_bytes")
    logical_limit = _limit(maximum_logical_bytes, "maximum_logical_bytes")
    gzip_member_limit = _integer(
        maximum_gzip_members,
        "maximum_gzip_members",
        maximum=MAXIMUM_GZIP_MEMBERS,
    )
    expected = validate_reference_catalog(payload)
    observed = _observed_artifact(
        source_path,
        expected["inputs"]["wrapper"],
        maximum_source_bytes=source_limit,
        maximum_logical_bytes=logical_limit,
        maximum_gzip_members=gzip_member_limit,
    )
    if observed != expected:
        raise _fail("reference FASTA replay does not match catalog")
    report = _validation_report(observed)
    validate_reference_report(report)
    return report


@_bounded
def validate_reference_paths(
    source_path: str | Path,
    artifact_path: str | Path,
    *,
    maximum_source_bytes: int = DEFAULT_MAX_INPUT_BYTES,
    maximum_logical_bytes: int = DEFAULT_MAX_LOGICAL_BYTES,
    maximum_artifact_bytes: int = MAX_ARTIFACT_BYTES,
    maximum_gzip_members: int = DEFAULT_MAXIMUM_GZIP_MEMBERS,
) -> dict[str, Any]:
    """Load a safe catalog path, replay its source path, and return a report."""

    artifact_limit = _integer(
        maximum_artifact_bytes,
        "maximum_artifact_bytes",
        maximum=MAX_ARTIFACT_BYTES,
    )
    raw = _read_regular(
        artifact_path,
        label="reference catalog",
        maximum_bytes=artifact_limit,
    )
    payload = _load_json(raw, "reference catalog")
    return validate_reference_fasta(
        source_path,
        payload,
        maximum_source_bytes=maximum_source_bytes,
        maximum_logical_bytes=maximum_logical_bytes,
        maximum_gzip_members=maximum_gzip_members,
    )


__all__ = [
    "DEFAULT_MAX_INPUT_BYTES",
    "DEFAULT_MAX_LOGICAL_BYTES",
    "DEFAULT_MAXIMUM_GZIP_MEMBERS",
    "FORMAT",
    "MAX_ARTIFACT_BYTES",
    "MAX_RECORD_BASES",
    "MAX_RECORDS",
    "MAXIMUM_GZIP_MEMBERS",
    "PROFILE",
    "PRODUCER",
    "REPORT_CHECKS",
    "REPORT_FORMAT",
    "ReferenceValidationError",
    "VERSION",
    "digest",
    "validate_reference_catalog",
    "validate_reference_fasta",
    "validate_reference_paths",
    "validate_reference_report",
]

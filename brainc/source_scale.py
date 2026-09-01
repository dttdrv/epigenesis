"""Streaming, reference-only admission for chromosome and genome FASTA."""

from __future__ import annotations

import base64
import copy
from dataclasses import dataclass
import functools
import gzip
import hashlib
import os
from pathlib import Path
import re
import stat
from typing import Any, BinaryIO
import zlib

from ._canonical import ContractError, SAFE_INTEGER, artifact_digest, canonical_bytes, digest
from ._io import MAX_IDENTIFIER_BYTES, MAX_SOURCE_RECORDS, MAX_STRING_BYTES


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

# NCBI's documented per-sequence technical ceiling. The collection ceiling
# remains an I-JSON safe integer and does not constrain an ordinary genome.
MAX_RECORD_BASES = 2_147_483_647
MAX_RECORDS = MAX_SOURCE_RECORDS
MAX_GZIP_MEMBERS = 1_000_000
DEFAULT_GZIP_MEMBERS = 100_000
# A 64 GiB default admits ordinary whole-genome FASTA with ample headroom over
# the 3.1 Gb human reference while bounding accidental or adversarial expansion.
# Exceptional assemblies may raise either execution budget explicitly.
DEFAULT_MAX_INPUT_BYTES = 64 * 1024**3
DEFAULT_MAX_LOGICAL_BYTES = 64 * 1024**3
_CHUNK_BYTES = 64 * 1024
_UPPER = bytes.maketrans(
    b"acgtryswkmbdhvn",
    b"ACGTRYSWKMBDHVN",
)
_INVALID_SEQUENCE = re.compile(rb"[^ACGTRYSWKMBDHVNacgtryswkmbdhvn\r\n]")
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_REFGET = re.compile(r"SQ\.[A-Za-z0-9_-]{32}\Z")
_GA4GH_DIGEST = re.compile(r"[A-Za-z0-9_-]{32}\Z")


class SourceScaleError(ContractError):
    """The reference FASTA or its compact catalog is invalid."""


def _fail(detail: str) -> SourceScaleError:
    return SourceScaleError(f"SCALE001: {detail}")


def _bounded(function: Any) -> Any:
    @functools.wraps(function)
    def guarded(*args: Any, **kwargs: Any) -> Any:
        try:
            return function(*args, **kwargs)
        except SourceScaleError:
            raise
        except MemoryError as failure:
            raise _fail("source-scale memory ceiling exceeded") from failure

    return guarded


@dataclass(frozen=True, slots=True)
class ScaleLimits:
    """Caller resource budget; accepted inputs have budget-independent output."""

    maximum_input_bytes: int = DEFAULT_MAX_INPUT_BYTES
    maximum_logical_bytes: int = DEFAULT_MAX_LOGICAL_BYTES
    maximum_records: int = MAX_RECORDS
    maximum_record_bases: int = MAX_RECORD_BASES
    maximum_header_bytes: int = MAX_STRING_BYTES
    maximum_gzip_members: int = DEFAULT_GZIP_MEMBERS

    def __post_init__(self) -> None:
        bounds = (
            ("maximum_input_bytes", self.maximum_input_bytes, SAFE_INTEGER),
            ("maximum_logical_bytes", self.maximum_logical_bytes, SAFE_INTEGER),
            ("maximum_records", self.maximum_records, MAX_RECORDS),
            ("maximum_record_bases", self.maximum_record_bases, MAX_RECORD_BASES),
            ("maximum_header_bytes", self.maximum_header_bytes, MAX_STRING_BYTES),
            ("maximum_gzip_members", self.maximum_gzip_members, MAX_GZIP_MEMBERS),
        )
        for name, value, maximum in bounds:
            if type(value) is not int or not 1 <= value <= maximum:
                raise _fail(f"{name} must be an integer in [1, {maximum}]")


DEFAULT_LIMITS = ScaleLimits()


class _HashingReader:
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
            raise _fail(f"physical source exceeds byte budget {self.maximum}")
        self.hasher.update(raw)
        return raw


class _StrictGzipReader:
    """Incrementally decode a closed, bounded concatenated gzip stream."""

    __slots__ = (
        "source",
        "pending",
        "decoder",
        "members",
        "finished",
        "maximum_members",
    )

    def __init__(self, source: _HashingReader, maximum_members: int) -> None:
        self.source = source
        self.pending = b""
        self.decoder: zlib.Decompress | None = None
        self.members = 0
        self.finished = False
        self.maximum_members = maximum_members

    def _start_member(self) -> None:
        while len(self.pending) < 2:
            chunk = self.source.read(_CHUNK_BYTES)
            if not chunk:
                if self.pending:
                    raise _fail("gzip FASTA contains trailing or truncated member bytes")
                if self.members == 0:
                    raise _fail("gzip FASTA contains no gzip member")
                self.finished = True
                return
            self.pending += chunk
        if not self.pending.startswith(b"\x1f\x8b"):
            raise _fail("gzip FASTA contains bytes outside a gzip member")
        self.members += 1
        if self.members > self.maximum_members:
            raise _fail(
                f"gzip FASTA exceeds member ceiling {self.maximum_members}"
            )
        self.decoder = zlib.decompressobj(16 + zlib.MAX_WBITS)

    def read(self, size: int = -1) -> bytes:
        target = _CHUNK_BYTES if size < 0 else size
        if target == 0 or self.finished:
            return b""
        output = bytearray()
        while len(output) < target and not self.finished:
            if self.decoder is None:
                self._start_member()
                if self.finished:
                    break
            assert self.decoder is not None
            if not self.pending:
                self.pending = self.source.read(_CHUNK_BYTES)
                if not self.pending:
                    raise _fail("gzip FASTA member is truncated")
            try:
                produced = self.decoder.decompress(
                    self.pending,
                    target - len(output),
                )
            except zlib.error as failure:
                raise _fail(f"invalid gzip FASTA source: {failure}") from failure
            output.extend(produced)
            if self.decoder.eof:
                self.pending = self.decoder.unused_data
                self.decoder = None
            else:
                self.pending = self.decoder.unconsumed_tail
        return bytes(output)


def _sha512t24u(value: bytes) -> str:
    return base64.urlsafe_b64encode(hashlib.sha512(value).digest()[:24]).decode("ascii")


def _refget(hasher: Any) -> str:
    return "SQ." + base64.urlsafe_b64encode(hasher.digest()[:24]).decode("ascii")


def _seqcol(records: list[dict[str, Any]]) -> dict[str, Any]:
    level_2 = {
        "lengths": [record["bases"] for record in records],
        "names": [record["record_id"] for record in records],
        "sequences": [record["refget_id"] for record in records],
    }
    level_1 = {
        name: _sha512t24u(canonical_bytes(value))
        for name, value in level_2.items()
    }
    inherent = {name: level_1[name] for name in ("names", "sequences")}
    return {
        "version": "1.0.0",
        "digest": _sha512t24u(canonical_bytes(inherent)),
        "level_1": level_1,
        "level_2": level_2,
    }


class _FastaScanner:
    __slots__ = (
        "limits",
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

    def __init__(self, limits: ScaleLimits) -> None:
        self.limits = limits
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

    def _start_header(self) -> None:
        self._finish_record()
        self.header = bytearray()
        self.at_line_start = False

    def _append_header(self, raw: bytes) -> None:
        assert self.header is not None
        if len(self.header) + len(raw) > self.limits.maximum_header_bytes:
            raise _fail(
                f"FASTA defline exceeds byte budget {self.limits.maximum_header_bytes}"
            )
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
        first = defline.split(maxsplit=1)
        record_id = first[0]
        description = first[1] if len(first) == 2 else ""
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
        if len(self.records) >= self.limits.maximum_records:
            raise _fail(f"FASTA exceeds record budget {self.limits.maximum_records}")
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
        new_bases = self.record_bases + len(canonical)
        if new_bases > self.limits.maximum_record_bases:
            raise _fail(
                f"FASTA record exceeds base budget {self.limits.maximum_record_bases}"
            )
        self.record_bases = new_bases
        self.record_sha256.update(canonical)
        self.record_sha512.update(canonical)
        if raw:
            self.at_line_start = raw.endswith((b"\r", b"\n"))

    def _finish_record(self) -> None:
        if self.record_id is None:
            return
        if self.record_bases == 0:
            raise _fail(f"FASTA record {self.record_id!r} has no sequence")
        self.total_bases += self.record_bases
        if self.total_bases > SAFE_INTEGER:
            raise _fail("FASTA total base count exceeds I-JSON safe integer range")
        self.records.append(
            {
                "record_id": self.record_id,
                "bases": self.record_bases,
                "sequence_sha256": self.record_sha256.hexdigest(),
                "refget_id": _refget(self.record_sha512),
            }
        )
        self.record_id = None

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
            chunk = stream.read(_CHUNK_BYTES)
            if not chunk:
                break
            self.logical_count += len(chunk)
            if self.logical_count > self.limits.maximum_logical_bytes:
                raise _fail(
                    f"logical FASTA exceeds byte budget {self.limits.maximum_logical_bytes}"
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


def _open_regular(path: str | Path, maximum_bytes: int) -> tuple[BinaryIO, os.stat_result]:
    descriptor = -1
    try:
        source = Path(path)
        inspected = source.lstat()
        if (
            stat.S_ISLNK(inspected.st_mode)
            or not stat.S_ISREG(inspected.st_mode)
            or inspected.st_nlink != 1
        ):
            raise _fail("reference FASTA must be one regular non-linked file")
        if inspected.st_size > maximum_bytes:
            raise _fail(f"physical source exceeds byte budget {maximum_bytes}")
        flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
        flags |= getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
        descriptor = os.open(source, flags)
        opened = os.fstat(descriptor)
        if (
            not stat.S_ISREG(opened.st_mode)
            or opened.st_nlink != 1
            or (opened.st_dev, opened.st_ino) != (inspected.st_dev, inspected.st_ino)
        ):
            raise _fail("reference FASTA changed while opening")
        if opened.st_size > maximum_bytes:
            raise _fail(f"physical source exceeds byte budget {maximum_bytes}")
        stream = os.fdopen(descriptor, "rb", buffering=0)
        descriptor = -1
        return stream, opened
    except SourceScaleError:
        raise
    except (OSError, TypeError, ValueError) as failure:
        raise _fail(f"cannot open reference FASTA: {failure}") from failure
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _stable(
    path: str | Path,
    opened: os.stat_result,
    finished: os.stat_result,
    byte_count: int,
) -> None:
    fields = ("st_dev", "st_ino", "st_nlink", "st_size", "st_mtime_ns", "st_ctime_ns")
    if any(getattr(opened, field) != getattr(finished, field) for field in fields):
        raise _fail("reference FASTA changed while reading")
    if byte_count != opened.st_size:
        raise _fail("reference FASTA byte count changed while reading")
    attached = Path(path).lstat()
    if stat.S_ISLNK(attached.st_mode) or any(
        getattr(opened, field) != getattr(attached, field) for field in fields
    ):
        raise _fail("reference FASTA path changed while reading")


@_bounded
def compile_reference_fasta(
    path: str | Path,
    *,
    wrapper: str,
    limits: ScaleLimits = DEFAULT_LIMITS,
) -> dict[str, Any]:
    """Compile explicitly wrapped FASTA without materializing sequence strings."""

    if type(wrapper) is not str or wrapper not in {"identity", "gzip"}:
        raise _fail("wrapper must be identity or gzip")
    if type(limits) is not ScaleLimits:
        raise _fail("limits must be ScaleLimits")
    stream, opened = _open_regular(path, limits.maximum_input_bytes)
    physical = _HashingReader(stream, limits.maximum_input_bytes)
    try:
        logical: Any = physical
        if wrapper == "gzip":
            logical = _StrictGzipReader(physical, limits.maximum_gzip_members)
        records, total_bases, logical_sha256, logical_bytes = _FastaScanner(limits).scan(logical)
        while physical.read(_CHUNK_BYTES):
            pass
        _stable(path, opened, os.fstat(stream.fileno()), physical.count)
    except SourceScaleError:
        raise
    except (EOFError, gzip.BadGzipFile, OSError, zlib.error) as failure:
        raise _fail(f"invalid {wrapper} FASTA source: {failure}") from failure
    finally:
        stream.close()

    catalog = {
        "records": records,
        "total_bases": total_bases,
        "refget_seqcol": _seqcol(records),
    }
    core = {
        "format": FORMAT,
        "version": VERSION,
        "producer": copy.deepcopy(PRODUCER),
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
    artifact = {**core, "artifact_sha256": digest(core)}
    return validate_reference_catalog(artifact)


def _keys(value: Any, expected: set[str], label: str) -> dict[str, Any]:
    if type(value) is not dict:
        raise _fail(f"{label} must be an object")
    if any(type(key) is not str for key in value):
        raise _fail(f"{label} keys must be strings")
    missing = sorted(expected - value.keys())
    extra = sorted(value.keys() - expected)
    if missing or extra:
        raise _fail(
            f"{label} keys invalid; missing={missing or 'none'}, unknown={extra or 'none'}"
        )
    return value


def _integer(value: Any, label: str, *, maximum: int = SAFE_INTEGER) -> int:
    if type(value) is not int or not 1 <= value <= maximum:
        raise _fail(f"{label} must be an integer in [1, {maximum}]")
    return value


def _sha(value: Any, label: str) -> str:
    if type(value) is not str or _SHA256.fullmatch(value) is None:
        raise _fail(f"{label} must be a lowercase SHA-256 digest")
    return value


def _validate_seqcol(value: Any, records: list[dict[str, Any]]) -> None:
    seqcol = _keys(value, {"version", "digest", "level_1", "level_2"}, "refget_seqcol")
    if seqcol["version"] != "1.0.0" or type(seqcol["version"]) is not str:
        raise _fail("refget_seqcol.version is unsupported")
    if type(seqcol["digest"]) is not str or _GA4GH_DIGEST.fullmatch(seqcol["digest"]) is None:
        raise _fail("refget_seqcol.digest is invalid")
    level_1 = _keys(seqcol["level_1"], {"lengths", "names", "sequences"}, "refget_seqcol.level_1")
    if any(
        type(level_1[name]) is not str or _GA4GH_DIGEST.fullmatch(level_1[name]) is None
        for name in ("lengths", "names", "sequences")
    ):
        raise _fail("refget_seqcol.level_1 contains an invalid digest")
    level_2 = _keys(seqcol["level_2"], {"lengths", "names", "sequences"}, "refget_seqcol.level_2")
    if (
        type(level_2["lengths"]) is not list
        or any(type(value) is not int for value in level_2["lengths"])
        or type(level_2["names"]) is not list
        or any(type(value) is not str for value in level_2["names"])
        or type(level_2["sequences"]) is not list
        or any(type(value) is not str for value in level_2["sequences"])
    ):
        raise _fail("refget_seqcol.level_2 has invalid JSON types")
    if seqcol != _seqcol(records):
        raise _fail("sequence_catalog.refget_seqcol does not match records")


@_bounded
def validate_reference_catalog(payload: dict[str, Any]) -> dict[str, Any]:
    """Validate a compact catalog without trusting its derived fields."""

    item = _keys(
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
    if item["format"] != FORMAT or type(item["version"]) is not int or item["version"] != VERSION:
        raise _fail("unsupported reference sequence catalog format/version")
    if type(item["producer"]) is not dict or item["producer"] != PRODUCER:
        raise _fail("reference sequence catalog producer is unsupported")
    inputs = _keys(item["inputs"], {"profile", "wrapper", "source", "logical"}, "inputs")
    if inputs["profile"] != PROFILE or inputs["wrapper"] not in {"identity", "gzip"}:
        raise _fail("reference sequence catalog profile/wrapper is unsupported")
    for name in ("source", "logical"):
        reference = _keys(inputs[name], {"sha256", "byte_length"}, f"inputs.{name}")
        _sha(reference["sha256"], f"inputs.{name}.sha256")
        _integer(reference["byte_length"], f"inputs.{name}.byte_length")
    if inputs["wrapper"] == "identity" and inputs["source"] != inputs["logical"]:
        raise _fail("identity FASTA physical and logical inputs must be identical")

    catalog = _keys(
        item["sequence_catalog"],
        {"records", "total_bases", "refget_seqcol"},
        "sequence_catalog",
    )
    records = catalog["records"]
    if type(records) is not list or not 1 <= len(records) <= MAX_RECORDS:
        raise _fail(f"sequence_catalog.records must contain 1..{MAX_RECORDS} records")
    seen: set[str] = set()
    normalized: list[dict[str, Any]] = []
    for index, value in enumerate(records):
        record = _keys(
            value,
            {"record_id", "bases", "sequence_sha256", "refget_id"},
            f"sequence_catalog.records[{index}]",
        )
        record_id = record["record_id"]
        try:
            record_id_bytes = record_id.encode("utf-8") if type(record_id) is str else b""
        except UnicodeEncodeError as failure:
            raise _fail(f"sequence_catalog.records[{index}].record_id is invalid") from failure
        if (
            type(record_id) is not str
            or not record_id
            or record_id != record_id.strip()
            or any(character.isspace() for character in record_id)
            or any(
                ord(character) < 0x20 or 0x7F <= ord(character) <= 0x9F
                for character in record_id
            )
            or len(record_id_bytes) > MAX_IDENTIFIER_BYTES
        ):
            raise _fail(f"sequence_catalog.records[{index}].record_id is invalid")
        if record_id in seen:
            raise _fail("sequence_catalog.records contains duplicate identifiers")
        seen.add(record_id)
        _integer(record["bases"], f"sequence_catalog.records[{index}].bases", maximum=MAX_RECORD_BASES)
        _sha(record["sequence_sha256"], f"sequence_catalog.records[{index}].sequence_sha256")
        if type(record["refget_id"]) is not str or _REFGET.fullmatch(record["refget_id"]) is None:
            raise _fail(f"sequence_catalog.records[{index}].refget_id is invalid")
        normalized.append(record)
    total_bases = sum(record["bases"] for record in normalized)
    _integer(catalog["total_bases"], "sequence_catalog.total_bases")
    if total_bases > SAFE_INTEGER or catalog["total_bases"] != total_bases:
        raise _fail("sequence_catalog.total_bases does not match records")
    _validate_seqcol(catalog["refget_seqcol"], normalized)
    if _sha(item["catalog_sha256"], "catalog_sha256") != digest(catalog):
        raise _fail("catalog_sha256 does not match sequence_catalog")
    try:
        artifact_digest(item)
    except (ContractError, KeyError) as failure:
        raise _fail(f"reference sequence catalog digest is invalid: {failure}") from failure
    return copy.deepcopy(item)


@_bounded
def replay_reference_fasta(
    path: str | Path,
    payload: dict[str, Any],
    *,
    limits: ScaleLimits = DEFAULT_LIMITS,
) -> dict[str, Any]:
    """Recompile original bytes and require the exact catalog artifact."""

    expected = validate_reference_catalog(payload)
    observed = compile_reference_fasta(
        path,
        wrapper=expected["inputs"]["wrapper"],
        limits=limits,
    )
    if observed != expected:
        raise _fail("reference FASTA replay does not match catalog")
    return expected


__all__ = [
    "DEFAULT_LIMITS",
    "DEFAULT_MAX_INPUT_BYTES",
    "DEFAULT_MAX_LOGICAL_BYTES",
    "DEFAULT_GZIP_MEMBERS",
    "FORMAT",
    "MAX_RECORD_BASES",
    "MAX_RECORDS",
    "MAX_GZIP_MEMBERS",
    "PROFILE",
    "PRODUCER",
    "ScaleLimits",
    "SourceScaleError",
    "VERSION",
    "compile_reference_fasta",
    "replay_reference_fasta",
    "validate_reference_catalog",
]

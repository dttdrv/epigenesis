"""Compact, content-bound sequence collection ABI for large GenBank records.

Version 2 keeps exact source bytes and normalized DNA in canonical 1 MiB ASCII
chunks.  The legacy Sequence IR and Sequence Collection v1 modules are not
imported or reinterpreted.
"""

from __future__ import annotations

import base64
import copy
from dataclasses import dataclass
import hashlib
from pathlib import Path
import re
from typing import Any, Iterable

from ._canonical import ContractError, artifact_digest, canonical_bytes, digest
from ._io import (
    BoundedIOError,
    atomic_write_file,
    compact_json_bytes,
    load_json_object,
    read_regular_file,
    validate_json_tree,
)


FORMAT = "brain01.sequence-collection-ir"
VERSION = 2
GENBANK_PROFILE = "genbank-273-traditional-dna-physical-structural/v2"
COMPILER = {
    "name": "brainc-dna-collection",
    "version": "0.3.0",
    "passes": [
        "validate-source-binding",
        "validate-iupac",
        "normalize-uppercase",
        "chunk-inline-ascii",
        "emit-sequence-collection-v2",
    ],
}

CHUNK_KIND = "chunked-inline-ascii"
CHUNK_VERSION = 1
CHUNK_BYTES = 1024 * 1024
MAX_SOURCE_BYTES = 16 * 1024 * 1024
MAX_TOTAL_BASES = 16 * 1024 * 1024
MAX_ARTIFACT_BYTES = 64 * 1024 * 1024
MAX_JSON_DEPTH = 64
MAX_JSON_MEMBERS = 1_000_000
MAX_RECORDS = 40_000
MAX_RECORD_ID_BYTES = 256
MAX_LINES = 1_000_000

_SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
_REFGET_RE = re.compile(r"SQ\.[A-Za-z0-9_-]{32}\Z")


class SequenceCollectionV2Error(ContractError):
    """A source or artifact violates Sequence Collection v2."""


def _fail(code: str, detail: str) -> SequenceCollectionV2Error:
    return SequenceCollectionV2Error(f"{code}: {detail}")


@dataclass(frozen=True, slots=True)
class SequenceMemberInput:
    """One normalized record and the inclusive raw ORIGIN line range."""

    record_id: str
    sequence: str
    line_start: int
    line_end: int


def _keys(value: Any, expected: set[str], label: str) -> dict[str, Any]:
    if type(value) is not dict:
        raise _fail("DNAC201", f"{label} must be an object")
    missing = sorted(expected - value.keys())
    extra = sorted(value.keys() - expected)
    if missing or extra:
        raise _fail(
            "DNAC201",
            f"{label} keys invalid; missing={missing or 'none'}, unknown={extra or 'none'}",
        )
    return value


def _integer(value: Any, label: str, *, minimum: int = 0, maximum: int = 2**53 - 1) -> int:
    if type(value) is not int or not minimum <= value <= maximum:
        raise _fail("DNAC202", f"{label} must be an integer in [{minimum}, {maximum}]")
    return value


def _sha256(value: Any, label: str) -> str:
    if type(value) is not str or _SHA256_RE.fullmatch(value) is None:
        raise _fail("DNAC203", f"{label} must be a lowercase SHA-256 digest")
    return value


def _record_id(value: Any, label: str) -> str:
    if type(value) is not str or not value or value != value.strip():
        raise _fail("DNAC204", f"{label} must be nonempty trimmed text")
    try:
        encoded = value.encode("ascii")
    except UnicodeEncodeError as failure:
        raise _fail("DNAC204", f"{label} must be ASCII") from failure
    if len(encoded) > MAX_RECORD_ID_BYTES or any(octet < 0x21 or octet > 0x7E for octet in encoded):
        raise _fail(
            "DNAC204",
            f"{label} must be at most {MAX_RECORD_ID_BYTES} printable non-space ASCII bytes",
        )
    return value


def _storage(raw: bytes, label: str) -> dict[str, Any]:
    if type(raw) is not bytes or not raw:
        raise _fail("DNAC205", f"{label} must be nonempty bytes")
    chunks: list[str] = []
    for start in range(0, len(raw), CHUNK_BYTES):
        try:
            chunks.append(raw[start : start + CHUNK_BYTES].decode("ascii"))
        except UnicodeDecodeError as failure:
            raise _fail("DNAC205", f"{label} must be ASCII") from failure
    return {
        "kind": CHUNK_KIND,
        "version": CHUNK_VERSION,
        "chunk_bytes": CHUNK_BYTES,
        "chunks": chunks,
    }


def _storage_bytes(value: Any, label: str, *, maximum: int) -> bytes:
    root = _keys(value, {"kind", "version", "chunk_bytes", "chunks"}, label)
    if (
        root["kind"] != CHUNK_KIND
        or type(root["version"]) is not int
        or root["version"] != CHUNK_VERSION
        or type(root["chunk_bytes"]) is not int
        or root["chunk_bytes"] != CHUNK_BYTES
    ):
        raise _fail("DNAC205", f"{label} has an unsupported chunk contract")
    chunks = root["chunks"]
    maximum_chunks = (maximum + CHUNK_BYTES - 1) // CHUNK_BYTES
    if type(chunks) is not list or not 1 <= len(chunks) <= maximum_chunks:
        raise _fail("DNAC205", f"{label}.chunks must be a bounded nonempty array")
    encoded_chunks: list[bytes] = []
    total = 0
    for index, chunk in enumerate(chunks):
        if type(chunk) is not str:
            raise _fail("DNAC205", f"{label}.chunks[{index}] must be text")
        try:
            encoded = chunk.encode("ascii")
        except UnicodeEncodeError as failure:
            raise _fail("DNAC205", f"{label}.chunks[{index}] must be ASCII") from failure
        if index < len(chunks) - 1:
            if len(encoded) != CHUNK_BYTES:
                raise _fail("DNAC205", f"{label} has a short non-final chunk")
        elif not 1 <= len(encoded) <= CHUNK_BYTES:
            raise _fail("DNAC205", f"{label} has an empty or oversized final chunk")
        total += len(encoded)
        if total > maximum:
            raise _fail("DNAC205", f"{label} exceeds {maximum} bytes")
        encoded_chunks.append(encoded)
    return b"".join(encoded_chunks)


def _sha512t24u(raw: bytes) -> str:
    return base64.urlsafe_b64encode(hashlib.sha512(raw).digest()[:24]).decode("ascii").rstrip("=")


def _refget(sequence: bytes) -> str:
    return "SQ." + _sha512t24u(sequence)


def _refget_seqcol(members: list[dict[str, Any]]) -> dict[str, Any]:
    level_2 = {
        "lengths": [member["sequence"]["bases"] for member in members],
        "names": [member["record_id"] for member in members],
        "sequences": [member["sequence"]["refget_id"] for member in members],
    }
    try:
        level_1 = {
            name: _sha512t24u(canonical_bytes(value))
            for name, value in level_2.items()
        }
        inherent = {name: level_1[name] for name in ("names", "sequences")}
        collection_digest = _sha512t24u(canonical_bytes(inherent))
    except ContractError as failure:
        raise _fail("DNAC206", f"cannot compute refget Sequence Collection identity: {failure}") from failure
    return {
        "version": "1.0.0",
        "digest": collection_digest,
        "level_1": level_1,
        "level_2": level_2,
    }


def _source_lines(raw: bytes) -> list[bytes]:
    physical_lines = raw.count(b"\n") if raw.endswith(b"\n") else raw.count(b"\n") + 1
    if not 1 <= physical_lines <= MAX_LINES:
        raise _fail("DNAC207", f"source must contain 1..{MAX_LINES} physical lines")
    carriage = raw.find(b"\r")
    while carriage >= 0:
        if carriage + 1 >= len(raw) or raw[carriage + 1] != 0x0A:
            raise _fail("DNAC207", "source contains a bare carriage return")
        carriage = raw.find(b"\r", carriage + 2)
    lines = raw.split(b"\n")
    if lines and lines[-1] == b"":
        lines.pop()
    if not lines or len(lines) > MAX_LINES:
        raise _fail("DNAC207", f"source must contain 1..{MAX_LINES} physical lines")
    result: list[bytes] = []
    for index, line in enumerate(lines, start=1):
        if line.endswith(b"\r"):
            line = line[:-1]
        if b"\r" in line:
            raise _fail("DNAC207", f"source line {index} contains a bare carriage return")
        if any(octet < 0x20 or octet > 0x7E for octet in line):
            raise _fail("DNAC207", f"source line {index} is not printable ASCII")
        result.append(line)
    return result


def _origin_from_lines(lines: list[bytes], line_start: int, line_end: int) -> bytes:
    if not 1 <= line_start <= line_end <= len(lines):
        raise _fail("DNAC208", "GenBank ORIGIN source map is outside the exact source")
    parts: list[bytes] = []
    lengths: list[int] = []
    normalized = 0
    for line_number in range(line_start, line_end + 1):
        line = lines[line_number - 1]
        if len(line) < 11 or line[9:10] != b" ":
            raise _fail("DNAC208", f"source line {line_number} is not an ORIGIN sequence line")
        if line[:9] != f"{normalized + 1:>9}".encode("ascii"):
            raise _fail("DNAC208", f"source line {line_number} has an invalid ORIGIN index")
        groups = line[10:].split(b" ")
        if (
            not 1 <= len(groups) <= 6
            or any(not group for group in groups)
            or any(len(group) != 10 for group in groups[:-1])
            or not 1 <= len(groups[-1]) <= 10
        ):
            raise _fail("DNAC208", f"source line {line_number} has invalid ORIGIN grouping")
        length = sum(len(group) for group in groups)
        if length > 60 or any(octet not in b"acgtryswkmbdhvn" for group in groups for octet in group):
            raise _fail("DNAC208", f"source line {line_number} is not lowercase IUPAC DNA")
        parts.append(b"".join(groups).upper())
        normalized += length
        lengths.append(length)
    if any(length != 60 for length in lengths[:-1]):
        raise _fail("DNAC208", "every non-final mapped ORIGIN line must contain 60 bases")
    return b"".join(parts)


def _record_binding(
    lines: list[bytes],
    line_start: int,
    line_end: int,
) -> tuple[str, int, int, int]:
    """Return VERSION, LOCUS length, record start, and terminator line."""

    if line_start < 2:
        raise _fail("DNAC208", "mapped sequence lines must immediately follow ORIGIN")
    origin = lines[line_start - 2]
    if origin not in {b"ORIGIN", b"ORIGIN      "}:
        if not origin.startswith(b"ORIGIN      "):
            raise _fail("DNAC208", "mapped sequence lines must immediately follow ORIGIN")
        origin_value = origin[12:]
        if (
            not origin_value
            or origin_value != origin_value.strip()
            or not origin_value.endswith(b".")
        ):
            raise _fail("DNAC208", "nonblank ORIGIN text must be trimmed and end in a period")
    terminator_line = line_end + 1
    if terminator_line > len(lines) or lines[terminator_line - 1] != b"//":
        raise _fail("DNAC208", "mapped sequence lines must end immediately before //")

    record_start = line_start - 1
    while record_start > 0 and not lines[record_start - 1].startswith(b"LOCUS       "):
        if lines[record_start - 1] == b"//":
            break
        record_start -= 1
    if record_start < 1 or not lines[record_start - 1].startswith(b"LOCUS       "):
        raise _fail("DNAC208", "mapped ORIGIN has no owning LOCUS record")

    try:
        locus_tokens = lines[record_start - 1][12:].decode("ascii").split()
    except UnicodeDecodeError as failure:
        raise _fail("DNAC208", "LOCUS line is not ASCII") from failure
    if len(locus_tokens) not in {6, 7}:
        raise _fail("DNAC208", "LOCUS line has an invalid token count")
    length_token = locus_tokens[1]
    if (
        not length_token.isascii()
        or not length_token.isdecimal()
        or length_token.startswith("0")
        or len(length_token) > len(str(MAX_TOTAL_BASES))
    ):
        raise _fail("DNAC208", "LOCUS length is not a positive ASCII decimal")
    declared_length = int(length_token)
    if declared_length > MAX_TOTAL_BASES:
        raise _fail("DNAC208", f"LOCUS length exceeds {MAX_TOTAL_BASES} bases")
    if locus_tokens[2] != "bp" or re.fullmatch(r"(?:(?:ss|ds|ms)-)?DNA", locus_tokens[3]) is None:
        raise _fail("DNAC208", "LOCUS must declare base-pair DNA")
    if len(locus_tokens) == 7 and locus_tokens[4] not in {"linear", "circular"}:
        raise _fail("DNAC208", "LOCUS topology is invalid")

    versions: list[str] = []
    for line in lines[record_start: line_start - 1]:
        if line.startswith(b"VERSION     "):
            try:
                value = line[12:].decode("ascii").split()
            except UnicodeDecodeError as failure:
                raise _fail("DNAC208", "VERSION line is not ASCII") from failure
            if not value:
                raise _fail("DNAC208", "VERSION line is empty")
            versions.append(value[0])
    if len(versions) != 1:
        raise _fail("DNAC208", "each mapped record must contain exactly one VERSION")
    return versions[0], declared_length, record_start, terminator_line


def _complete_record_envelope(
    lines: list[bytes],
    starts: list[int],
    ends: list[int],
) -> None:
    if not starts or len(starts) != len(ends) or starts[0] != 1:
        raise _fail("DNAC208", "source must begin with the first mapped LOCUS record")
    for index, (start, end) in enumerate(zip(starts, ends)):
        if index and start not in {ends[index - 1] + 1, ends[index - 1] + 2}:
            raise _fail("DNAC208", "source records are not contiguous")
        if index and start == ends[index - 1] + 2 and lines[ends[index - 1]] != b"":
            raise _fail("DNAC208", "only an exact blank record separator is accepted")
    trailing = lines[ends[-1] :]
    if trailing not in ([], [b""]):
        raise _fail("DNAC208", "source contains trailing data outside the mapped records")


def _validate_refget_seqcol(value: Any, members: list[dict[str, Any]]) -> None:
    root = _keys(value, {"version", "digest", "level_1", "level_2"}, "refget_seqcol")
    if root["version"] != "1.0.0" or type(root["digest"]) is not str or re.fullmatch(
        r"[A-Za-z0-9_-]{32}", root["digest"]
    ) is None:
        raise _fail("DNAC217", "refget Sequence Collection root identity is invalid")
    level_1 = _keys(root["level_1"], {"lengths", "names", "sequences"}, "refget_seqcol.level_1")
    if any(
        type(level_1[name]) is not str
        or re.fullmatch(r"[A-Za-z0-9_-]{32}", level_1[name]) is None
        for name in ("lengths", "names", "sequences")
    ):
        raise _fail("DNAC217", "refget Sequence Collection level_1 identities are invalid")
    level_2 = _keys(root["level_2"], {"lengths", "names", "sequences"}, "refget_seqcol.level_2")
    lengths, names, sequences = level_2["lengths"], level_2["names"], level_2["sequences"]
    if (
        type(lengths) is not list
        or type(names) is not list
        or type(sequences) is not list
        or not all(type(item) is int and 1 <= item <= MAX_TOTAL_BASES for item in lengths)
        or not all(type(item) is str for item in names)
        or not all(type(item) is str and _REFGET_RE.fullmatch(item) is not None for item in sequences)
    ):
        raise _fail("DNAC217", "refget Sequence Collection level_2 values are invalid")
    expected = _refget_seqcol(members)
    try:
        same = canonical_bytes(root) == canonical_bytes(expected)
    except ContractError as failure:
        raise _fail("DNAC217", f"refget Sequence Collection is not canonical JSON: {failure}") from failure
    if not same:
        raise _fail("DNAC217", "refget Sequence Collection identity does not match members")


def _artifact_bytes(payload: dict[str, Any]) -> bytes:
    try:
        return compact_json_bytes(
            payload,
            ensure_ascii=False,
            maximum_bytes=MAX_ARTIFACT_BYTES,
            maximum_depth=MAX_JSON_DEPTH,
            maximum_members=MAX_JSON_MEMBERS,
            maximum_string_bytes=CHUNK_BYTES,
        )
    except BoundedIOError as failure:
        raise _fail("DNAC209", f"Sequence Collection v2 is outside output limits: {failure}") from failure


def _build(source: bytes, profile: str, records: Iterable[SequenceMemberInput]) -> dict[str, Any]:
    if type(source) is not bytes or not source or len(source) > MAX_SOURCE_BYTES:
        raise _fail("DNAC210", f"source must be 1..{MAX_SOURCE_BYTES} bytes")
    try:
        source.decode("ascii")
    except UnicodeDecodeError as failure:
        raise _fail("DNAC210", "source must be ASCII") from failure
    if type(profile) is not str or profile != GENBANK_PROFILE:
        raise _fail("DNAC210", f"profile must equal {GENBANK_PROFILE!r}")
    if isinstance(records, (str, bytes, bytearray)):
        raise _fail("DNAC210", "records must be an iterable of SequenceMemberInput values")

    lines = _source_lines(source)
    emitted: list[dict[str, Any]] = []
    seen: set[str] = set()
    total_bases = 0
    previous_line_end = 0
    previous_record_end = 0
    bound_record_starts: list[int] = []
    bound_record_ends: list[int] = []
    try:
        iterator = iter(records)
    except TypeError as failure:
        raise _fail("DNAC210", "records must be iterable") from failure
    for item in iterator:
        if len(emitted) >= MAX_RECORDS:
            raise _fail("DNAC210", f"records exceed {MAX_RECORDS}")
        if type(item) is not SequenceMemberInput:
            raise _fail("DNAC210", "every record must be SequenceMemberInput")
        record_id = _record_id(item.record_id, "record_id")
        if record_id in seen:
            raise _fail("DNAC210", f"duplicate record_id {record_id!r}")
        seen.add(record_id)
        line_start = _integer(item.line_start, "line_start", minimum=1, maximum=MAX_LINES)
        line_end = _integer(item.line_end, "line_end", minimum=line_start, maximum=MAX_LINES)
        if line_start <= previous_line_end:
            raise _fail("DNAC210", "record ORIGIN source maps must be ordered and disjoint")
        previous_line_end = line_end
        if type(item.sequence) is not str or not item.sequence:
            raise _fail("DNAC210", f"record {record_id!r} sequence must be nonempty text")
        try:
            sequence = item.sequence.upper().encode("ascii")
        except UnicodeEncodeError as failure:
            raise _fail("DNAC210", f"record {record_id!r} sequence must be ASCII") from failure
        if any(octet not in b"ACGTRYSWKMBDHVN" for octet in sequence):
            raise _fail("DNAC210", f"record {record_id!r} is not IUPAC DNA")
        total_bases += len(sequence)
        if total_bases > MAX_TOTAL_BASES:
            raise _fail("DNAC210", f"normalized DNA exceeds {MAX_TOTAL_BASES} bases")
        replayed = _origin_from_lines(lines, line_start, line_end)
        if replayed != sequence:
            raise _fail("DNAC210", f"record {record_id!r} does not replay from its ORIGIN source map")
        source_record_id, declared_length, record_start, record_end = _record_binding(
            lines,
            line_start,
            line_end,
        )
        if source_record_id != record_id or declared_length != len(sequence):
            raise _fail(
                "DNAC210",
                f"record {record_id!r} differs from its VERSION or LOCUS length",
            )
        if record_start <= previous_record_end:
            raise _fail("DNAC210", "GenBank records must be mapped in exact source order")
        previous_record_end = record_end
        bound_record_starts.append(record_start)
        bound_record_ends.append(record_end)
        emitted.append(
            {
                "record_id": record_id,
                "source_map": {
                    "kind": "genbank-origin-lines",
                    "version": 1,
                    "line_start": line_start,
                    "line_end": line_end,
                },
                "sequence": {
                    "alphabet": "IUPAC-DNA",
                    "normalization": "uppercase",
                    "bases": len(sequence),
                    "sha256": hashlib.sha256(sequence).hexdigest(),
                    "refget_id": _refget(sequence),
                    "storage": _storage(sequence, f"record {record_id!r} sequence"),
                },
            }
        )
    if not emitted:
        raise _fail("DNAC210", "records must be nonempty")
    source_record_starts = [
        line_number
        for line_number, line in enumerate(lines, start=1)
        if line.startswith(b"LOCUS       ")
    ]
    if source_record_starts != bound_record_starts:
        raise _fail("DNAC210", "every source LOCUS record must be mapped exactly once")
    _complete_record_envelope(lines, bound_record_starts, bound_record_ends)

    core = {
        "format": FORMAT,
        "version": VERSION,
        "compiler": copy.deepcopy(COMPILER),
        "inputs": {
            "kind": "genbank-flatfile",
            "profile": profile,
            "source": {
                "encoding": "ascii",
                "byte_length": len(source),
                "sha256": hashlib.sha256(source).hexdigest(),
                "storage": _storage(source, "GenBank source"),
            },
        },
        "members": emitted,
        "refget_seqcol": _refget_seqcol(emitted),
    }
    payload = {**core, "artifact_sha256": digest(core)}
    _artifact_bytes(payload)
    return payload


def validate_sequence_collection_v2(payload: dict[str, Any]) -> dict[str, Any]:
    """Validate the complete closed v2 artifact and exact raw-to-sequence maps."""

    try:
        validate_json_tree(
            payload,
            "Sequence Collection v2",
            maximum_depth=MAX_JSON_DEPTH,
            maximum_members=MAX_JSON_MEMBERS,
            maximum_string_bytes=CHUNK_BYTES,
        )
    except BoundedIOError as failure:
        raise _fail("DNAC211", f"artifact tree is outside contract limits: {failure}") from failure
    root = _keys(
        payload,
        {"format", "version", "compiler", "inputs", "members", "refget_seqcol", "artifact_sha256"},
        "Sequence Collection v2",
    )
    if root["format"] != FORMAT or type(root["version"]) is not int or root["version"] != VERSION:
        raise _fail("DNAC211", "unsupported Sequence Collection format or version")
    if root["compiler"] != COMPILER:
        raise _fail("DNAC211", "unsupported Sequence Collection compiler identity")
    try:
        artifact_digest(root)
    except ContractError as failure:
        raise _fail("DNAC212", str(failure)) from failure

    inputs = _keys(root["inputs"], {"kind", "profile", "source"}, "inputs")
    if inputs["kind"] != "genbank-flatfile" or inputs["profile"] != GENBANK_PROFILE:
        raise _fail("DNAC213", "input kind or profile is invalid")
    source_value = _keys(
        inputs["source"],
        {"encoding", "byte_length", "sha256", "storage"},
        "inputs.source",
    )
    if source_value["encoding"] != "ascii":
        raise _fail("DNAC213", "inputs.source.encoding must be ascii")
    byte_length = _integer(
        source_value["byte_length"],
        "inputs.source.byte_length",
        minimum=1,
        maximum=MAX_SOURCE_BYTES,
    )
    source_sha = _sha256(source_value["sha256"], "inputs.source.sha256")
    source = _storage_bytes(source_value["storage"], "inputs.source.storage", maximum=MAX_SOURCE_BYTES)
    if len(source) != byte_length or hashlib.sha256(source).hexdigest() != source_sha:
        raise _fail("DNAC213", "exact source storage differs from its length or SHA-256")
    lines = _source_lines(source)

    members_value = root["members"]
    if type(members_value) is not list or not 1 <= len(members_value) <= MAX_RECORDS:
        raise _fail("DNAC214", "members must be a bounded nonempty array")
    normalized_members: list[dict[str, Any]] = []
    seen: set[str] = set()
    total_bases = 0
    previous_line_end = 0
    previous_record_end = 0
    bound_record_starts: list[int] = []
    bound_record_ends: list[int] = []
    for index, member_value in enumerate(members_value):
        member = _keys(member_value, {"record_id", "source_map", "sequence"}, f"members[{index}]")
        record_id = _record_id(member["record_id"], f"members[{index}].record_id")
        if record_id in seen:
            raise _fail("DNAC214", f"duplicate record_id {record_id!r}")
        seen.add(record_id)
        source_map = _keys(
            member["source_map"],
            {"kind", "version", "line_start", "line_end"},
            f"members[{index}].source_map",
        )
        if (
            source_map["kind"] != "genbank-origin-lines"
            or type(source_map["version"]) is not int
            or source_map["version"] != 1
        ):
            raise _fail("DNAC214", "member source map contract is invalid")
        line_start = _integer(source_map["line_start"], "line_start", minimum=1, maximum=MAX_LINES)
        line_end = _integer(source_map["line_end"], "line_end", minimum=line_start, maximum=MAX_LINES)
        if line_start <= previous_line_end:
            raise _fail("DNAC214", "member source maps must be ordered and disjoint")
        previous_line_end = line_end

        sequence_value = _keys(
            member["sequence"],
            {"alphabet", "normalization", "bases", "sha256", "refget_id", "storage"},
            f"members[{index}].sequence",
        )
        if sequence_value["alphabet"] != "IUPAC-DNA" or sequence_value["normalization"] != "uppercase":
            raise _fail("DNAC215", "sequence alphabet or normalization is invalid")
        bases = _integer(
            sequence_value["bases"],
            f"members[{index}].sequence.bases",
            minimum=1,
            maximum=MAX_TOTAL_BASES,
        )
        sequence_sha = _sha256(sequence_value["sha256"], f"members[{index}].sequence.sha256")
        refget_id = sequence_value["refget_id"]
        if type(refget_id) is not str or _REFGET_RE.fullmatch(refget_id) is None:
            raise _fail("DNAC215", "sequence refget_id is invalid")
        sequence = _storage_bytes(
            sequence_value["storage"],
            f"members[{index}].sequence.storage",
            maximum=MAX_TOTAL_BASES,
        )
        if any(octet not in b"ACGTRYSWKMBDHVN" for octet in sequence):
            raise _fail("DNAC215", "normalized sequence storage is not uppercase IUPAC DNA")
        if (
            len(sequence) != bases
            or hashlib.sha256(sequence).hexdigest() != sequence_sha
            or _refget(sequence) != refget_id
        ):
            raise _fail("DNAC215", "sequence storage differs from its identities")
        total_bases += bases
        if total_bases > MAX_TOTAL_BASES:
            raise _fail("DNAC215", f"collection exceeds {MAX_TOTAL_BASES} bases")
        if _origin_from_lines(lines, line_start, line_end) != sequence:
            raise _fail("DNAC216", f"record {record_id!r} does not replay from exact source lines")
        source_record_id, declared_length, record_start, record_end = _record_binding(
            lines,
            line_start,
            line_end,
        )
        if source_record_id != record_id or declared_length != bases:
            raise _fail("DNAC216", "member identity or length differs from its source record")
        if record_start <= previous_record_end:
            raise _fail("DNAC216", "source records are not mapped in exact order")
        previous_record_end = record_end
        bound_record_starts.append(record_start)
        bound_record_ends.append(record_end)
        normalized_members.append(member)

    source_record_starts = [
        line_number
        for line_number, line in enumerate(lines, start=1)
        if line.startswith(b"LOCUS       ")
    ]
    if source_record_starts != bound_record_starts:
        raise _fail("DNAC216", "every source LOCUS record must be mapped exactly once")
    _complete_record_envelope(lines, bound_record_starts, bound_record_ends)

    _validate_refget_seqcol(root["refget_seqcol"], normalized_members)
    _artifact_bytes(root)
    return root


@dataclass(frozen=True, slots=True)
class SequenceCollectionV2Artifact:
    """Immutable-by-interface Sequence Collection v2 artifact."""

    _payload: dict[str, Any]

    @property
    def digest(self) -> str:
        return self._payload["artifact_sha256"]

    def to_dict(self) -> dict[str, Any]:
        return copy.deepcopy(self._payload)

    def source_bytes(self) -> bytes:
        return _storage_bytes(
            self._payload["inputs"]["source"]["storage"],
            "inputs.source.storage",
            maximum=MAX_SOURCE_BYTES,
        )

    def save(self, path: str | Path) -> None:
        raw = _artifact_bytes(self._payload)
        try:
            atomic_write_file(
                path,
                raw,
                maximum_bytes=MAX_ARTIFACT_BYTES,
                label="Sequence Collection v2 output",
            )
        except BoundedIOError as failure:
            raise _fail("DNAC218", f"cannot save Sequence Collection v2: {failure}") from failure


class SequenceCollectionV2Compiler:
    """Build the compact collection from exact GenBank bytes and parsed records."""

    def compile_genbank(
        self,
        source: bytes,
        records: Iterable[SequenceMemberInput],
        *,
        profile: str = GENBANK_PROFILE,
    ) -> SequenceCollectionV2Artifact:
        payload = _build(source, profile, records)
        return SequenceCollectionV2Artifact(payload)


def as_sequence_collection_v2(payload: dict[str, Any]) -> SequenceCollectionV2Artifact:
    validate_sequence_collection_v2(payload)
    return SequenceCollectionV2Artifact(copy.deepcopy(payload))


def load_sequence_collection_v2(path: str | Path) -> SequenceCollectionV2Artifact:
    try:
        raw = read_regular_file(
            path,
            maximum_bytes=MAX_ARTIFACT_BYTES,
            label="Sequence Collection v2 artifact",
        )
        payload = load_json_object(
            raw,
            "Sequence Collection v2 artifact",
            maximum_bytes=MAX_ARTIFACT_BYTES,
            maximum_depth=MAX_JSON_DEPTH,
            maximum_members=MAX_JSON_MEMBERS,
            maximum_string_bytes=CHUNK_BYTES,
        )
    except BoundedIOError as failure:
        raise _fail("DNAC219", f"cannot load Sequence Collection v2: {failure}") from failure
    return as_sequence_collection_v2(payload)


__all__ = [
    "CHUNK_BYTES",
    "COMPILER",
    "FORMAT",
    "GENBANK_PROFILE",
    "MAX_ARTIFACT_BYTES",
    "MAX_SOURCE_BYTES",
    "MAX_TOTAL_BASES",
    "SequenceCollectionV2Artifact",
    "SequenceCollectionV2Compiler",
    "SequenceCollectionV2Error",
    "SequenceMemberInput",
    "VERSION",
    "as_sequence_collection_v2",
    "load_sequence_collection_v2",
    "validate_sequence_collection_v2",
]

"""Closed GenBank flat-file frontend for physical DNA records.

The profile is deliberately named and versioned: it parses traditional
GenBank 273.0 records under INSDC Feature Table 11.4, preserves their source
spelling and maps, and emits deterministic structural IR.  Sequence and
feature meaning is never obtained from a network or an accession-specific
table.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
from datetime import date as calendar_date
import re
from pathlib import Path
from typing import Any

from ._canonical import ContractError, SAFE_INTEGER, artifact_digest, digest
from ._io import (
    BoundedIOError,
    MAX_JSON_MEMBERS,
    MAX_STRING_BYTES,
    atomic_write_file,
    compact_json_bytes,
    load_json_object,
    read_regular_file,
    validate_json_tree,
)
from .insdc_location import (
    Between,
    Complement,
    Endpoint,
    INSDCLocationError,
    Join,
    Location,
    MAX_ACCESSION_BYTES,
    NormalizedSegment,
    Order,
    Point,
    Remote,
    SequenceContext,
    Span,
    Within,
    normalize_segments,
    parse_location,
)
from .sequence_collection_v2 import (
    MAX_ARTIFACT_BYTES,
    MAX_RECORDS as COLLECTION_MAX_RECORDS,
    MAX_SOURCE_BYTES,
    MAX_TOTAL_BASES,
    SequenceCollectionV2Artifact,
    SequenceCollectionV2Compiler,
    SequenceCollectionV2Error,
    SequenceMemberInput,
    as_sequence_collection_v2,
)


FORMAT = "brainc.bio.insdc-genbank-ir"
VERSION = 2
PROFILE = "genbank-273-traditional-dna-physical-structural/v2"
FEATURE_ID_PREFIX = "io.github.dttdrv.epigenesis.insdc-feature.sha256."
AUTHORITY_MANIFEST_SHA256 = "6dc3658a5f7d6a774ce9cd08ff6dc5327cbef8a9951c7679fd9a34650f912bf9"
AUTHORITY = {
    "genbank_flatfile": "NCBI GenBank 273.0",
    "insdc_feature_table": "INSDC Feature Table Definition 11.4",
    "authority_manifest_sha256": AUTHORITY_MANIFEST_SHA256,
}
COMPILER = {
    "name": "brainc-insdc-genbank",
    "version": "0.8.1",
    "passes": [
        "split-record-stream",
        "parse-genbank-headers",
        "parse-feature-table",
        "parse-origin-sequence",
        "resolve-record-local-remotes",
        "normalize-locations",
        "compile-sequence-collection-v2",
        "emit-insdc-genbank-ir",
        "seal-artifact",
    ],
}

MAX_GENBANK_BYTES = MAX_SOURCE_BYTES
MAX_OUTPUT_BYTES = MAX_ARTIFACT_BYTES
MAX_LINE_BYTES = MAX_STRING_BYTES
MAX_LINES = 1_000_000
MAX_RECORDS = COLLECTION_MAX_RECORDS
MAX_FEATURES = 20_000
MAX_QUALIFIERS = 100_000
MAX_ACCESSIONS = 10_000
MAX_FEATURE_LINE_BYTES = 80
MAX_JSON_DEPTH = 64
MAX_TEXT_BYTES = MAX_STRING_BYTES
_AUTHORITY_MANIFEST_BYTES = 32 * 1024
_AUTHORITY_MANIFEST_PATH = (
    Path(__file__).with_name("standards")
    / "genbank-273-insdc-ft-11.4.authority.json"
)

_ACCESSION_RE = re.compile(r"[A-Z][A-Z0-9_]*[0-9]\Z")
_VERSION_RE = re.compile(r"([A-Z][A-Z0-9_]*[0-9])\.([1-9][0-9]*)\Z")
_MOLECULE_RE = re.compile(r"(?:(ss|ds|ms)-)?DNA\Z")
_DIVISION_RE = re.compile(r"[A-Z]{3}\Z")
_DATE_RE = re.compile(r"[0-9]{2}-[A-Z]{3}-[0-9]{4}\Z")
_MONTHS = {
    "JAN": 1,
    "FEB": 2,
    "MAR": 3,
    "APR": 4,
    "MAY": 5,
    "JUN": 6,
    "JUL": 7,
    "AUG": 8,
    "SEP": 9,
    "OCT": 10,
    "NOV": 11,
    "DEC": 12,
}
_COMPONENT_NAME_RE = re.compile(r"[A-Za-z0-9_*'\-]+\Z")
_IUPAC_LOWER = frozenset("acgtryswkmbdhvn")
_PRE_KEYWORDS_FIELDS = frozenset({"DBLINK", "PROJECT"})
_POST_ORGANISM_FIELDS = {
    "REFERENCE": 0,
    "COMMENT": 1,
    "PRIMARY": 2,
}


class GenBankError(ContractError):
    """A source or artifact violates the closed GenBank profile."""


def _fail(code: str, detail: str) -> GenBankError:
    return GenBankError(f"{code}: {detail}")


@dataclass(frozen=True, slots=True)
class _Line:
    number: int
    text: str
    eol: str


@dataclass(frozen=True, slots=True)
class _Header:
    name: str
    value: str
    chunks: tuple[dict[str, Any], ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "value": self.value,
            "source_map": [_position(chunk) for chunk in self.chunks],
        }


@dataclass(frozen=True, slots=True)
class _RawFeature:
    ordinal: int
    key: str
    location_text: str
    location_chunks: tuple[dict[str, Any], ...]
    qualifiers: tuple[dict[str, Any], ...]
    line_start: int
    line_end: int


@dataclass(slots=True)
class _QualifierBuilder:
    parts: list[str]
    chunks: list[dict[str, Any]]
    line: int
    name: str
    complete: bool
    continuation_mode: str


@dataclass(frozen=True, slots=True)
class _ParsedRecord:
    line_start: int
    line_end: int
    origin_line_start: int
    origin_line_end: int
    locus: dict[str, Any]
    definition: str
    accessions: tuple[str, ...]
    version: str
    legacy_gi: int | None
    keywords: str
    source: str
    organism: str
    lineage: str
    headers: tuple[_Header, ...]
    features: tuple[_RawFeature, ...]
    sequence: str


def _chunk(line: _Line, column: int, text: str) -> dict[str, Any]:
    if not 1 <= column <= len(line.text) + 1:
        raise _fail("INSDC001", f"line {line.number} has an invalid source column")
    return {
        "line": line.number,
        "column_start": column,
        "column_end_exclusive": column + len(text),
        "text": text,
    }


def _position(chunk: dict[str, Any]) -> dict[str, int]:
    """Project an internal parse chunk to its compact positional source map."""

    return {
        "line": chunk["line"],
        "column_start": chunk["column_start"],
        "column_end_exclusive": chunk["column_end_exclusive"],
    }


def _split_lines(raw: bytes) -> tuple[_Line, ...]:
    if type(raw) is not bytes:
        raise _fail("INSDC001", "GenBank source must be bytes")
    if not raw:
        raise _fail("INSDC002", "GenBank source is empty")
    if len(raw) > MAX_GENBANK_BYTES:
        raise _fail("INSDC003", f"GenBank source exceeds {MAX_GENBANK_BYTES} bytes")
    physical_lines = raw.count(b"\n") + (0 if raw.endswith(b"\n") else 1)
    if physical_lines > MAX_LINES:
        raise _fail("INSDC003", f"GenBank source exceeds {MAX_LINES} lines")
    lines: list[_Line] = []
    start = 0
    number = 1
    while start < len(raw):
        newline = raw.find(b"\n", start)
        if newline < 0:
            body = raw[start:]
            eol = "none"
            start = len(raw)
        else:
            if newline > start and raw[newline - 1] == 0x0D:
                body = raw[start : newline - 1]
                eol = "CRLF"
            else:
                body = raw[start:newline]
                eol = "LF"
            start = newline + 1
        if len(body) > MAX_LINE_BYTES:
            raise _fail("INSDC004", f"line {number} exceeds {MAX_LINE_BYTES} bytes")
        if b"\r" in body:
            raise _fail("INSDC005", f"line {number} contains a bare carriage return")
        if any(octet < 0x20 or octet > 0x7E for octet in body):
            raise _fail("INSDC005", f"line {number} is not printable ASCII (tabs are not accepted)")
        try:
            text = body.decode("ascii")
        except UnicodeDecodeError as failure:  # guarded above, kept as a closed boundary
            raise _fail("INSDC005", f"line {number} is not ASCII") from failure
        lines.append(_Line(number, text, eol))
        if len(lines) > MAX_LINES:
            raise _fail("INSDC003", f"GenBank source exceeds {MAX_LINES} lines")
        number += 1
    if not lines:
        raise _fail("INSDC002", "GenBank source is empty")
    for line in lines[:-1]:
        if line.eol == "none":
            raise _fail("INSDC007", f"line {line.number} is not terminated")
    return tuple(lines)


def _split_records(lines: tuple[_Line, ...]) -> tuple[tuple[_Line, ...], ...]:
    records: list[tuple[_Line, ...]] = []
    start = 0
    while start < len(lines):
        if not lines[start].text.startswith("LOCUS       "):
            raise _fail("INSDC008", f"record at line {lines[start].number} must begin with LOCUS")
        end = start
        while end < len(lines) and lines[end].text != "//":
            end += 1
        if end == len(lines):
            raise _fail("INSDC009", f"record at line {lines[start].number} has no exact // terminator")
        records.append(lines[start : end + 1])
        if len(records) > MAX_RECORDS:
            raise _fail("INSDC003", f"GenBank source exceeds {MAX_RECORDS} records")
        start = end + 1
        if start < len(lines) and lines[start].text == "":
            start += 1
    return tuple(records)


def _positive_decimal(token: str, label: str, line: int) -> int:
    if not token or not token.isascii() or not token.isdecimal() or (len(token) > 1 and token[0] == "0"):
        raise _fail("INSDC010", f"line {line} has an invalid {label}")
    value = 0
    for character in token:
        digit_value = ord(character) - ord("0")
        if value > (SAFE_INTEGER - digit_value) // 10:
            raise _fail("INSDC010", f"line {line} {label} exceeds the I-JSON safe range")
        value = value * 10 + digit_value
    if value < 1:
        raise _fail("INSDC010", f"line {line} {label} must be positive")
    return value


def _parse_locus(line: _Line) -> dict[str, Any]:
    if not line.text.startswith("LOCUS       "):
        raise _fail("INSDC011", f"line {line.number} is not a traditional LOCUS line")
    tokens = line.text[12:].split()
    if len(tokens) not in {6, 7}:
        raise _fail("INSDC011", f"line {line.number} LOCUS must contain six or seven tokens")
    if len(tokens) == 6:
        name, length_token, unit, molecule_token, division, date = tokens
        topology = "linear"
        topology_declared = False
    else:
        name, length_token, unit, molecule_token, topology, division, date = tokens
        topology_declared = True
    if len(name) > 64 or any(character.isspace() for character in name):
        raise _fail("INSDC011", f"line {line.number} has an invalid locus name")
    length = _positive_decimal(length_token, "LOCUS length", line.number)
    if length > MAX_TOTAL_BASES:
        raise _fail("INSDC003", f"line {line.number} record length exceeds {MAX_TOTAL_BASES} bases")
    if unit != "bp":
        raise _fail("INSDC012", f"line {line.number} is not a base-pair record")
    molecule_match = _MOLECULE_RE.fullmatch(molecule_token)
    if molecule_match is None:
        raise _fail("INSDC012", f"line {line.number} is not a DNA record")
    if topology not in {"linear", "circular"}:
        raise _fail("INSDC011", f"line {line.number} topology must be linear or circular")
    if _DIVISION_RE.fullmatch(division) is None:
        raise _fail("INSDC011", f"line {line.number} has an invalid division token")
    if _DATE_RE.fullmatch(date) is None:
        raise _fail("INSDC011", f"line {line.number} has an invalid date token")
    try:
        day_token, month_token, year_token = date.split("-")
        calendar_date(int(year_token), _MONTHS[month_token], int(day_token))
    except (KeyError, ValueError) as failure:
        raise _fail("INSDC011", f"line {line.number} has an invalid calendar date") from failure
    return {
        "name": name,
        "length": length,
        "unit": unit,
        "molecule": "DNA",
        "strandedness": molecule_match.group(1),
        "topology": topology,
        "topology_declared": topology_declared,
        "division": division,
        "date": date,
        "source_map": _chunk(line, 1, line.text),
    }


def _field_prefix(name: str) -> str:
    return f"{name:<12}"


def _parse_field(lines: tuple[_Line, ...], index: int, name: str) -> tuple[_Header, int]:
    if index >= len(lines) or not lines[index].text.startswith(_field_prefix(name)):
        actual = "end of record" if index >= len(lines) else repr(lines[index].text[:12].rstrip())
        raise _fail("INSDC013", f"expected {name} at line {lines[min(index, len(lines) - 1)].number}; found {actual}")
    first = lines[index]
    raw_parts = [first.text[12:]]
    chunks = [_chunk(first, 13, first.text[12:])]
    index += 1
    while index < len(lines) and lines[index].text.startswith(" " * 12):
        continuation = lines[index]
        raw_parts.append(continuation.text[12:])
        chunks.append(_chunk(continuation, 13, continuation.text[12:]))
        index += 1
    parts = [part.strip() for part in raw_parts]
    if any(not part for part in parts):
        raise _fail("INSDC014", f"{name} has an empty chunk")
    return _Header(name, " ".join(parts), tuple(chunks)), index


def _optional_block(lines: tuple[_Line, ...], index: int, name: str) -> tuple[_Header, int]:
    first = lines[index]
    if not first.text.startswith(_field_prefix(name)):
        raise _fail("INSDC015", f"line {first.number} does not begin {name}")
    block = [first]
    index += 1
    while index < len(lines) and lines[index].text.startswith(" "):
        block.append(lines[index])
        index += 1
    values: list[str] = []
    chunks: list[dict[str, Any]] = []
    for position, line in enumerate(block):
        if position == 0:
            text = line.text[12:]
            column = 13
        else:
            leading = len(line.text) - len(line.text.lstrip(" "))
            text = line.text[leading:]
            column = leading + 1
        if not text.strip() and name != "COMMENT":
            raise _fail("INSDC014", f"{name} has an empty chunk at line {line.number}")
        values.append(text.strip())
        chunks.append(_chunk(line, column, text))
    return _Header(name, " ".join(value for value in values if value), tuple(chunks)), index


def _organism_value(header: _Header) -> tuple[str, str]:
    parts = [chunk["text"].strip() for chunk in header.chunks]
    if len(parts) < 2:
        raise _fail("INSDC029", "ORGANISM must contain a taxonomic lineage")
    lineage_start = 1
    if len(parts[0]) == 68:
        lineage_start = len(parts) - 1
        for index, part in enumerate(parts[1:], start=1):
            if ";" in part and len([token for token in part.split(";") if token.strip()]) >= 2:
                lineage_start = index
                break
    organism = " ".join(parts[:lineage_start])
    lineage = " ".join(parts[lineage_start:])
    if not organism or not lineage:
        raise _fail("INSDC029", "ORGANISM name and lineage must be nonempty")
    return organism, lineage


def _parse_origin_header(line: _Line) -> _Header:
    if line.text in {"ORIGIN", _field_prefix("ORIGIN")}:
        return _Header("ORIGIN", "", (_chunk(line, 1, line.text),))
    prefix = _field_prefix("ORIGIN")
    if not line.text.startswith(prefix):
        raise _fail("INSDC022", f"line {line.number} is not an ORIGIN header")
    value = line.text[12:]
    if not value or value != value.strip() or not value.endswith("."):
        raise _fail(
            "INSDC022",
            f"line {line.number} nonblank ORIGIN text must be a trimmed value ending in a period",
        )
    return _Header("ORIGIN", value, (_chunk(line, 1, line.text),))


def _is_origin_header(line: _Line) -> bool:
    return line.text == "ORIGIN" or line.text.startswith(_field_prefix("ORIGIN"))


def _component_name(value: str, maximum: int, label: str, line: int) -> str:
    has_letter = any("A" <= character <= "Z" or "a" <= character <= "z" for character in value)
    if (
        len(value) > maximum
        or _COMPONENT_NAME_RE.fullmatch(value) is None
        or not has_letter
    ):
        raise _fail("INSDC016", f"line {line} {label} is not a {maximum}-character INSDC component name")
    return value


def _quoted_value(parts: list[str], name: str, line: int) -> tuple[bool, str | None]:
    """Return a complete normalized value from fixed-column source chunks."""

    joiner = "" if name == "translation" else " "
    source = joiner.join(part.strip() for part in parts)
    if not source.startswith('"'):
        if (
            not source
            or '"' in source
            or any(not 0x20 <= ord(character) <= 0x7E for character in source)
        ):
            raise _fail("INSDC017", f"qualifier /{name} at line {line} has an invalid unquoted value")
        return True, source
    index = 1
    value: list[str] = []
    while index < len(source):
        character = source[index]
        if character == '"':
            if index + 1 < len(source) and source[index + 1] == '"':
                value.append('"')
                index += 2
                continue
            if index != len(source) - 1:
                raise _fail("INSDC017", f"qualifier /{name} at line {line} has text after its closing quote")
            return True, "".join(value)
        if not 0x20 <= ord(character) <= 0x7E:
            raise _fail("INSDC017", f"qualifier /{name} at line {line} is not printable ASCII")
        value.append(character)
        index += 1
    return False, None


def _parse_qualifier_parts(parts: list[str], line: int) -> tuple[bool, str, str | None]:
    first = parts[0]
    if not first.startswith("/"):
        raise _fail("INSDC018", f"line {line} qualifier must begin with '/'")
    body = first[1:]
    if "=" not in body:
        name = _component_name(body, 20, "qualifier name", line)
        if len(parts) != 1:
            raise _fail("INSDC018", f"valueless qualifier /{name} cannot continue")
        return True, name, None
    name, initial = body.split("=", 1)
    name = _component_name(name, 20, "qualifier name", line)
    complete, value = _quoted_value([initial, *parts[1:]], name, line)
    return complete, name, value


def _finish_qualifier(parts: list[str], chunks: list[dict[str, Any]], line: int) -> dict[str, Any]:
    complete, name, value = _parse_qualifier_parts(parts, line)
    if not complete:
        raise _fail("INSDC017", f"qualifier /{name} at line {line} has no closing quote")
    return {"name": name, "value": value, "chunks": copy.deepcopy(chunks)}


def _continued_quoted_value_complete(part: str, name: str, line: int) -> bool:
    """Scan one continuation fragment while a quoted value remains open."""

    source = part.strip()
    index = 0
    while index < len(source):
        character = source[index]
        if character == '"':
            if index + 1 < len(source) and source[index + 1] == '"':
                index += 2
                continue
            if index != len(source) - 1:
                raise _fail(
                    "INSDC017",
                    f"qualifier /{name} at line {line} has text after its closing quote",
                )
            return True
        if not 0x20 <= ord(character) <= 0x7E:
            raise _fail(
                "INSDC017",
                f"qualifier /{name} at line {line} is not printable ASCII",
            )
        index += 1
    return False


def _parse_features(lines: tuple[_Line, ...], index: int) -> tuple[tuple[_RawFeature, ...], int]:
    features: list[_RawFeature] = []
    qualifier_count = 0
    while index < len(lines) and not _is_origin_header(lines[index]):
        line = lines[index]
        if line.text.startswith("CONTIG"):
            raise _fail("INSDC019", f"line {line.number} is an assembled CONTIG record, not physical ORIGIN DNA")
        if len(line.text) > MAX_FEATURE_LINE_BYTES:
            raise _fail("INSDC020", f"feature-table line {line.number} extends past column {MAX_FEATURE_LINE_BYTES}")
        if len(line.text) < 22 or line.text[:5] != " " * 5 or line.text[20] != " ":
            raise _fail("INSDC020", f"line {line.number} is not a fixed-column feature descriptor")
        key_field = line.text[5:20]
        key = key_field.rstrip(" ")
        if not key or key_field != key.ljust(15):
            raise _fail("INSDC020", f"line {line.number} has a misaligned feature key")
        key = _component_name(key, 15, "feature key", line.number)
        location_parts = [line.text[21:].strip()]
        if not location_parts[0]:
            raise _fail("INSDC021", f"feature at line {line.number} has no location")
        location_chunks = [_chunk(line, 22, line.text[21:])]
        line_start = line.number
        line_end = line.number
        index += 1

        qualifier_builders: list[_QualifierBuilder] = []
        current: _QualifierBuilder | None = None
        while index < len(lines):
            continuation = lines[index]
            if _is_origin_header(continuation):
                break
            if len(continuation.text) >= 22 and continuation.text[:21] == " " * 21:
                if len(continuation.text) > MAX_FEATURE_LINE_BYTES:
                    raise _fail("INSDC020", f"feature-table line {continuation.number} extends past column {MAX_FEATURE_LINE_BYTES}")
                body = continuation.text[21:]
                if not body:
                    raise _fail("INSDC020", f"feature-table line {continuation.number} has no column-22 content")
                if current is None and not qualifier_builders and not body.startswith("/"):
                    location_parts.append(body.strip())
                    location_chunks.append(_chunk(continuation, 22, body))
                else:
                    if body.startswith("/") and (current is None or current.complete):
                        complete, name, _ = _parse_qualifier_parts(
                            [body],
                            continuation.number,
                        )
                        current = _QualifierBuilder(
                            parts=[body],
                            chunks=[_chunk(continuation, 22, body)],
                            line=continuation.number,
                            name=name,
                            complete=complete,
                            continuation_mode=(
                                "valueless"
                                if "=" not in body[1:]
                                else "quoted"
                                if body.split("=", 1)[1].startswith('"')
                                else "unquoted"
                            ),
                        )
                        qualifier_builders.append(current)
                    elif current is not None and not current.complete:
                        current.parts.append(body)
                        current.chunks.append(_chunk(continuation, 22, body))
                        current.complete = _continued_quoted_value_complete(
                            body,
                            current.name,
                            current.line,
                        )
                    elif (
                        current is not None
                        and current.complete
                        and current.continuation_mode == "unquoted"
                    ):
                        current.parts.append(body)
                        current.chunks.append(_chunk(continuation, 22, body))
                    else:
                        raise _fail("INSDC018", f"line {continuation.number} is an unexpected qualifier continuation")
                line_end = continuation.number
                index += 1
                continue
            if continuation.text.startswith(" " * 5) and len(continuation.text) >= 21:
                break
            raise _fail("INSDC020", f"line {continuation.number} is not a feature, continuation, or ORIGIN")

        qualifiers: list[dict[str, Any]] = []
        for builder in qualifier_builders:
            qualifiers.append(
                _finish_qualifier(builder.parts, builder.chunks, builder.line)
            )
        qualifier_count += len(qualifiers)
        if qualifier_count > MAX_QUALIFIERS:
            raise _fail("INSDC003", f"record stream exceeds {MAX_QUALIFIERS} qualifiers")
        features.append(
            _RawFeature(
                ordinal=len(features) + 1,
                key=key,
                location_text="".join(location_parts),
                location_chunks=tuple(location_chunks),
                qualifiers=tuple(qualifiers),
                line_start=line_start,
                line_end=line_end,
            )
        )
        if len(features) > MAX_FEATURES:
            raise _fail("INSDC003", f"record stream exceeds {MAX_FEATURES} features")
    if index >= len(lines) or not _is_origin_header(lines[index]):
        line_number = lines[-1].number
        raise _fail("INSDC022", f"record ending at line {line_number} has no ORIGIN header")
    return tuple(features), index


def _parse_origin(
    lines: tuple[_Line, ...],
    index: int,
    expected_length: int,
) -> tuple[str, int, int, int]:
    first_index = index
    sequence_parts: list[str] = []
    line_lengths: list[int] = []
    normalized = 0
    while index < len(lines) and lines[index].text != "//":
        line = lines[index]
        if len(line.text) < 11 or line.text[9] != " ":
            raise _fail("INSDC023", f"line {line.number} is not an NCBI ORIGIN sequence line")
        index_token = line.text[:9]
        wanted_index = normalized + 1
        if index_token != f"{wanted_index:>9}":
            raise _fail("INSDC023", f"line {line.number} ORIGIN index must be {wanted_index}")
        payload = line.text[10:]
        groups = payload.split(" ")
        if not 1 <= len(groups) <= 6 or any(not group for group in groups):
            raise _fail("INSDC023", f"line {line.number} has invalid ORIGIN grouping")
        if any(len(group) != 10 for group in groups[:-1]) or not 1 <= len(groups[-1]) <= 10:
            raise _fail("INSDC023", f"line {line.number} ORIGIN groups must contain ten bases except the last")
        line_length = sum(map(len, groups))
        if line_length > 60:
            raise _fail("INSDC023", f"line {line.number} contains more than 60 bases")
        for group in groups:
            if any(base not in _IUPAC_LOWER for base in group):
                raise _fail("INSDC024", f"line {line.number} ORIGIN sequence must be lowercase IUPAC DNA")
            sequence_parts.append(group.upper())
            normalized += len(group)
        line_lengths.append(line_length)
        index += 1
    if not line_lengths:
        raise _fail("INSDC025", "ORIGIN must contain at least one physical sequence line")
    if any(length != 60 for length in line_lengths[:-1]):
        raise _fail("INSDC023", "every non-final ORIGIN line must contain exactly 60 bases")
    sequence = "".join(sequence_parts)
    if len(sequence) != expected_length:
        raise _fail("INSDC025", f"LOCUS declares {expected_length} bases but ORIGIN contains {len(sequence)}")
    return sequence, lines[first_index].number, lines[index - 1].number, index


def _accessions(value: str, line: int) -> tuple[str, ...]:
    tokens = value.split()
    if not tokens or len(tokens) > MAX_ACCESSIONS:
        raise _fail("INSDC026", f"line {line} has an invalid ACCESSION list")
    if len(tokens) != len(set(tokens)) or any(
        _ACCESSION_RE.fullmatch(token) is None
        or len(token.encode("ascii")) > MAX_ACCESSION_BYTES
        for token in tokens
    ):
        raise _fail("INSDC026", f"line {line} ACCESSION values must be unique uppercase accessions")
    return tuple(tokens)


def _record(record_lines: tuple[_Line, ...]) -> _ParsedRecord:
    if record_lines[-1].text != "//":
        raise _fail("INSDC009", f"record at line {record_lines[0].number} has no exact terminator")
    index = 0
    locus = _parse_locus(record_lines[index])
    index += 1
    headers: list[_Header] = []
    definition, index = _parse_field(record_lines, index, "DEFINITION")
    headers.append(definition)
    accession_header, index = _parse_field(record_lines, index, "ACCESSION")
    headers.append(accession_header)
    accessions = _accessions(accession_header.value, accession_header.chunks[0]["line"])
    version_header, index = _parse_field(record_lines, index, "VERSION")
    headers.append(version_header)
    version_tokens = version_header.value.split()
    if (
        len(version_tokens) not in {1, 2}
        or _VERSION_RE.fullmatch(version_tokens[0]) is None
        or (len(version_tokens) == 2 and not version_tokens[1].startswith("GI:"))
    ):
        raise _fail(
            "INSDC027",
            f"VERSION at line {version_header.chunks[0]['line']} must contain accession.version and optional GI:decimal",
        )
    version = version_tokens[0]
    if len(version.encode("ascii")) > MAX_ACCESSION_BYTES:
        raise _fail(
            "INSDC027",
            f"VERSION accession identity exceeds {MAX_ACCESSION_BYTES} ASCII bytes",
        )
    legacy_gi = None
    if len(version_tokens) == 2:
        legacy_gi = _positive_decimal(
            version_tokens[1][3:],
            "legacy VERSION GI",
            version_header.chunks[0]["line"],
        )
    version_match = _VERSION_RE.fullmatch(version)
    assert version_match is not None
    if version_match.group(1) != accessions[0]:
        raise _fail("INSDC027", "VERSION base accession must equal the primary ACCESSION")

    seen_pre: set[str] = set()
    while index < len(record_lines):
        name = record_lines[index].text[:12].rstrip()
        if name not in _PRE_KEYWORDS_FIELDS:
            break
        if name in seen_pre:
            raise _fail("INSDC028", f"duplicate {name} field at line {record_lines[index].number}")
        optional, index = _optional_block(record_lines, index, name)
        headers.append(optional)
        seen_pre.add(name)

    keywords, index = _parse_field(record_lines, index, "KEYWORDS")
    headers.append(keywords)
    if index < len(record_lines) and record_lines[index].text.startswith(
        _field_prefix("SEGMENT")
    ):
        segment, index = _parse_field(record_lines, index, "SEGMENT")
        segment_tokens = segment.value.split()
        if len(segment.chunks) != 1 or len(segment_tokens) != 3 or segment_tokens[1] != "of":
            raise _fail("INSDC031", "SEGMENT must be one line containing 'n of m'")
        segment_number = _positive_decimal(
            segment_tokens[0],
            "SEGMENT number",
            segment.chunks[0]["line"],
        )
        segment_total = _positive_decimal(
            segment_tokens[2],
            "SEGMENT total",
            segment.chunks[0]["line"],
        )
        if segment_total < 2 or segment_number > segment_total:
            raise _fail("INSDC031", "SEGMENT must satisfy 1 <= n <= m and m >= 2")
        headers.append(segment)
    source, index = _parse_field(record_lines, index, "SOURCE")
    headers.append(source)
    organism, index = _parse_field(record_lines, index, "  ORGANISM")
    headers.append(organism)
    organism_name, lineage = _organism_value(organism)

    last_rank = -1
    seen_singletons: set[str] = set()
    reference_count = 0
    while index < len(record_lines) and record_lines[index].text != "FEATURES             Location/Qualifiers":
        line = record_lines[index]
        if line.text.startswith("CONTIG"):
            raise _fail("INSDC019", f"line {line.number} is a CONTIG record, not physical ORIGIN DNA")
        if not line.text or line.text.startswith(" "):
            raise _fail("INSDC030", f"orphaned indented field at line {line.number}")
        name = line.text[:12].rstrip()
        if name not in _POST_ORGANISM_FIELDS:
            raise _fail("INSDC030", f"unknown top-level field {name!r} at line {line.number}")
        rank = _POST_ORGANISM_FIELDS[name]
        if rank < last_rank:
            raise _fail("INSDC031", f"misordered {name} field at line {line.number}")
        if name != "REFERENCE" and name in seen_singletons:
            raise _fail("INSDC028", f"duplicate {name} field at line {line.number}")
        optional, index = _optional_block(record_lines, index, name)
        if name == "REFERENCE":
            reference_count += 1
            reference_token = optional.chunks[0]["text"].split(maxsplit=1)[0]
            reference_number = _positive_decimal(
                reference_token,
                "REFERENCE number",
                optional.chunks[0]["line"],
            )
            if reference_number != reference_count:
                raise _fail(
                    "INSDC031",
                    f"REFERENCE at line {optional.chunks[0]['line']} must be numbered {reference_count}",
                )
            journals = [
                chunk
                for chunk in optional.chunks[1:]
                if chunk["column_start"] == 3
                and chunk["text"].startswith("JOURNAL   ")
            ]
            if len(journals) != 1:
                raise _fail(
                    "INSDC031",
                    f"REFERENCE at line {optional.chunks[0]['line']} must contain exactly one JOURNAL subkeyword",
                )
        headers.append(optional)
        seen_singletons.add(name)
        last_rank = rank
    if index >= len(record_lines) or record_lines[index].text != "FEATURES             Location/Qualifiers":
        raise _fail("INSDC032", f"record at line {record_lines[0].number} has no exact FEATURES header")
    if reference_count == 0:
        raise _fail("INSDC031", f"record at line {record_lines[0].number} must contain a REFERENCE")
    headers.append(
        _Header(
            "FEATURES",
            "Location/Qualifiers",
            (_chunk(record_lines[index], 1, record_lines[index].text),),
        )
    )
    index += 1
    features, index = _parse_features(record_lines, index)
    origin_header = _parse_origin_header(record_lines[index])
    headers.append(origin_header)
    sequence, origin_line_start, origin_line_end, index = _parse_origin(
        record_lines,
        index + 1,
        locus["length"],
    )
    if index >= len(record_lines) or record_lines[index].text != "//" or index != len(record_lines) - 1:
        raise _fail("INSDC009", f"record at line {record_lines[0].number} has an invalid terminator")
    return _ParsedRecord(
        line_start=record_lines[0].number,
        line_end=record_lines[-1].number,
        origin_line_start=origin_line_start,
        origin_line_end=origin_line_end,
        locus=locus,
        definition=definition.value,
        accessions=accessions,
        version=version,
        legacy_gi=legacy_gi,
        keywords=keywords.value,
        source=source.value,
        organism=organism_name,
        lineage=lineage,
        headers=tuple(headers),
        features=features,
        sequence=sequence,
    )


def _ast(node: Location) -> dict[str, Any]:
    if type(node) is Point:
        return {"kind": "point", "position": node.position}
    if type(node) is Span:
        return {
            "kind": "span",
            "start": {"position": node.start.position, "fuzz": node.start.fuzz},
            "end": {"position": node.end.position, "fuzz": node.end.fuzz},
        }
    if type(node) is Between:
        return {"kind": "between", "left": node.left, "right": node.right}
    if type(node) is Within:
        return {"kind": "within", "left": node.left, "right": node.right}
    if type(node) is Remote:
        return {"kind": "remote", "reference": node.reference, "location": _ast(node.location)}
    if type(node) is Complement:
        return {"kind": "complement", "location": _ast(node.location)}
    if type(node) is Join:
        return {"kind": "join", "children": [_ast(child) for child in node.children]}
    if type(node) is Order:
        return {"kind": "order", "children": [_ast(child) for child in node.children]}
    raise _fail("INSDC033", "location parser returned an unsupported AST node")


def _segment(segment: NormalizedSegment) -> dict[str, Any]:
    return {
        "reference": segment.reference,
        "start": segment.start,
        "end": segment.end,
        "orientation": segment.orientation,
        "kind": segment.kind,
        "start_fuzz": segment.start_fuzz,
        "end_fuzz": segment.end_fuzz,
        "bounds_status": segment.bounds_status,
    }


def _source_semantics(record: _ParsedRecord, features: list[dict[str, Any]]) -> None:
    sources = [feature for feature in features if feature["key"] == "source"]
    if not sources:
        raise _fail("INSDC034", f"record {record.version!r} must contain a source feature")

    coverage: list[tuple[int, int]] = []
    present_mol_types: list[str] = []
    missing_mol_type = False
    for source in sources:
        for segment in source["location"]["segments"]:
            if segment["reference"] is not None or segment["kind"] != "interval":
                raise _fail(
                    "INSDC034",
                    f"record {record.version!r} source locations must be local base spans",
                )
            coverage.append((segment["start"], segment["end"]))
        qualifiers = source["qualifiers"]
        organisms = [
            entry["value"] for entry in qualifiers if entry["name"] == "organism"
        ]
        mol_types = [
            entry["value"] for entry in qualifiers if entry["name"] == "mol_type"
        ]
        if len(organisms) != 1 or not organisms[0]:
            raise _fail(
                "INSDC035",
                f"record {record.version!r} each source /organism must occur once with a value",
            )
        if len(mol_types) > 1 or (mol_types and not mol_types[0]):
            raise _fail(
                "INSDC035",
                f"record {record.version!r} each source /mol_type may occur once and must have a value",
            )
        if mol_types:
            present_mol_types.append(mol_types[0])
        else:
            missing_mol_type = True

    covered_until = 0
    for start, end in sorted(coverage):
        if start > covered_until:
            raise _fail(
                "INSDC034",
                f"record {record.version!r} source features leave bases uncovered",
            )
        covered_until = max(covered_until, end)
    if covered_until != record.locus["length"]:
        raise _fail(
            "INSDC034",
            f"record {record.version!r} source features do not span the complete sequence",
        )
    if present_mol_types and (
        missing_mol_type or len(set(present_mol_types)) != 1
    ):
        raise _fail(
            "INSDC035",
            f"record {record.version!r} source /mol_type values must be uniformly present and equal",
        )


def _json_member_count(value: Any) -> int:
    members = 0
    stack = [value]
    while stack:
        current = stack.pop()
        if type(current) is dict:
            members += len(current)
            stack.extend(current.values())
        elif type(current) is list:
            members += len(current)
            stack.extend(current)
    return members


def _emit_record(
    record: _ParsedRecord,
    contexts: dict[str, SequenceContext],
    *,
    maximum_members: int,
) -> tuple[dict[str, Any], int]:
    local_context = contexts[record.version]
    locus = copy.deepcopy(record.locus)
    locus["source_map"] = _position(locus["source_map"])
    emitted_features: list[dict[str, Any]] = []
    core = {
        "record_id": record.version,
        "accessions": list(record.accessions),
        "version": record.version,
        "legacy_gi": record.legacy_gi,
        "locus": locus,
        "definition": record.definition,
        "keywords": record.keywords,
        "source": record.source,
        "organism": {"name": record.organism, "lineage": record.lineage},
        "headers": [header.to_dict() for header in record.headers],
        "features": emitted_features,
        "relationships": [],
        "source_map": {"line_start": record.line_start, "line_end": record.line_end},
    }
    record_members = _json_member_count(
        {**core, "record_ir_sha256": "0" * 64}
    )
    if record_members > maximum_members:
        raise _fail("INSDC039", "record metadata exceeds the remaining JSON member budget")
    for raw_feature in record.features:
        try:
            location = parse_location(
                raw_feature.location_text,
                context=local_context,
                references=contexts,
            )
            segments = normalize_segments(location, context=local_context, references=contexts)
        except INSDCLocationError as failure:
            raise _fail(
                "INSDC036",
                f"feature {raw_feature.ordinal} in {record.version} has an invalid location: {failure}",
            ) from failure
        identity_core = {
            "record_id": record.version,
            "ordinal": raw_feature.ordinal,
            "key": raw_feature.key,
            "location": raw_feature.location_text,
            "qualifiers": [
                {"name": qualifier["name"], "value": qualifier["value"]}
                for qualifier in raw_feature.qualifiers
            ],
        }
        emitted = {
            "feature_id": FEATURE_ID_PREFIX + digest(identity_core),
            "ordinal": raw_feature.ordinal,
            "key": raw_feature.key,
            "location": {
                "text": raw_feature.location_text,
                "ast": _ast(location),
                "segments": [_segment(segment) for segment in segments],
                "unresolved_references": sorted(
                    {
                        segment.reference
                        for segment in segments
                        if segment.bounds_status == "unresolved"
                        and segment.reference is not None
                    }
                ),
                "source_map": [
                    _position(chunk) for chunk in raw_feature.location_chunks
                ],
            },
            "qualifiers": [
                {
                    "name": qualifier["name"],
                    "value": qualifier["value"],
                    "source_map": [
                        _position(chunk) for chunk in qualifier["chunks"]
                    ],
                }
                for qualifier in raw_feature.qualifiers
            ],
            "source_map": {
                "line_start": raw_feature.line_start,
                "line_end": raw_feature.line_end,
            },
        }
        contribution = 1 + _json_member_count(emitted)
        if record_members + contribution > maximum_members:
            raise _fail(
                "INSDC039",
                f"source stream exceeds {MAX_JSON_MEMBERS} aggregate JSON members",
            )
        emitted_features.append(emitted)
        record_members += contribution
    _source_semantics(record, emitted_features)
    return {**core, "record_ir_sha256": digest(core)}, record_members


def _compile_parts(raw: bytes) -> tuple[SequenceCollectionV2Artifact, dict[str, Any]]:
    lines = _split_lines(raw)
    parsed = tuple(_record(record) for record in _split_records(lines))
    feature_count = sum(len(record.features) for record in parsed)
    qualifier_count = sum(
        len(feature.qualifiers)
        for record in parsed
        for feature in record.features
    )
    if feature_count > MAX_FEATURES:
        raise _fail(
            "INSDC003",
            f"source stream exceeds {MAX_FEATURES} aggregate features",
        )
    if qualifier_count > MAX_QUALIFIERS:
        raise _fail(
            "INSDC003",
            f"source stream exceeds {MAX_QUALIFIERS} aggregate qualifiers",
        )
    versions = [record.version for record in parsed]
    primary_accessions = [record.accessions[0] for record in parsed]
    if len(versions) != len(set(versions)):
        raise _fail("INSDC037", "record VERSION identities must be unique within the input")
    if len(primary_accessions) != len(set(primary_accessions)):
        raise _fail("INSDC037", "primary ACCESSION identities must be unique within the input")
    contexts = {
        record.version: SequenceContext(record.locus["length"], record.locus["topology"])
        for record in parsed
    }
    try:
        collection = SequenceCollectionV2Compiler().compile_genbank(
            raw,
            (
                SequenceMemberInput(
                    record_id=record.version,
                    sequence=record.sequence,
                    line_start=record.origin_line_start,
                    line_end=record.origin_line_end,
                )
                for record in parsed
            ),
            profile=PROFILE,
        )
    except SequenceCollectionV2Error as failure:
        raise _fail("INSDC038", f"cannot bind parsed ORIGIN sequences: {failure}") from failure
    collection_payload = collection.to_dict()
    records: list[dict[str, Any]] = []
    bio_ir = {
        "profile": PROFILE,
        "authority": copy.deepcopy(AUTHORITY),
        "sequence_collection_artifact_sha256": collection.digest,
        "coordinate_system": {
            "source": "1-based-closed",
            "normalized": "0-based-half-open",
        },
        "records": records,
    }
    member_skeleton = {
        "format": FORMAT,
        "version": VERSION,
        "profile": PROFILE,
        "authority": copy.deepcopy(AUTHORITY),
        "compiler": copy.deepcopy(COMPILER),
        "sequence_collection": collection_payload,
        "bio_ir": bio_ir,
        "bio_ir_sha256": "0" * 64,
        "artifact_sha256": "0" * 64,
    }
    used_members = _json_member_count(member_skeleton)
    if used_members > MAX_JSON_MEMBERS:
        raise _fail(
            "INSDC039",
            "collection and record metadata exceed the aggregate JSON member budget",
        )
    for record in parsed:
        remaining = MAX_JSON_MEMBERS - used_members - 1
        if remaining < 0:
            raise _fail(
                "INSDC039",
                f"source stream exceeds {MAX_JSON_MEMBERS} aggregate JSON members",
            )
        emitted, record_members = _emit_record(
            record,
            contexts,
            maximum_members=remaining,
        )
        records.append(emitted)
        used_members += 1 + record_members
    return collection, bio_ir


def _artifact_bytes(payload: dict[str, Any]) -> bytes:
    try:
        return compact_json_bytes(
            payload,
            ensure_ascii=False,
            maximum_bytes=MAX_OUTPUT_BYTES,
            maximum_depth=MAX_JSON_DEPTH,
            maximum_members=MAX_JSON_MEMBERS,
            maximum_string_bytes=MAX_TEXT_BYTES,
        )
    except BoundedIOError as failure:
        raise _fail("INSDC039", f"artifact is outside output limits: {failure}") from failure


def _seal(collection: SequenceCollectionV2Artifact, bio_ir: dict[str, Any]) -> dict[str, Any]:
    core = {
        "format": FORMAT,
        "version": VERSION,
        "profile": PROFILE,
        "authority": copy.deepcopy(AUTHORITY),
        "compiler": copy.deepcopy(COMPILER),
        "sequence_collection": collection.to_dict(),
        "bio_ir": bio_ir,
        "bio_ir_sha256": digest(bio_ir),
    }
    payload = {**core, "artifact_sha256": digest(core)}
    _artifact_bytes(payload)
    return payload


@dataclass(frozen=True, slots=True)
class GenBankArtifact:
    """Immutable-by-interface compiled GenBank artifact."""

    _payload: dict[str, Any]

    @property
    def digest(self) -> str:
        return self._payload["artifact_sha256"]

    @property
    def sequence_collection(self) -> SequenceCollectionV2Artifact:
        return as_sequence_collection_v2(copy.deepcopy(self._payload["sequence_collection"]))

    def to_dict(self) -> dict[str, Any]:
        return copy.deepcopy(self._payload)

    def save(self, path: str | Path) -> None:
        raw = _artifact_bytes(self._payload)
        try:
            atomic_write_file(path, raw, maximum_bytes=MAX_OUTPUT_BYTES, label="GenBank artifact output")
        except BoundedIOError as failure:
            raise _fail("INSDC040", f"cannot save GenBank artifact: {failure}") from failure


class GenBankCompiler:
    """Compile the named GenBank 273.0 physical-DNA profile."""

    def compile_bytes(self, raw: bytes) -> GenBankArtifact:
        _load_authority_manifest()
        collection, bio_ir = _compile_parts(raw)
        return GenBankArtifact(_seal(collection, bio_ir))

    def compile_file(self, path: str | Path) -> GenBankArtifact:
        try:
            raw = read_regular_file(path, maximum_bytes=MAX_GENBANK_BYTES, label="GenBank source")
        except BoundedIOError as failure:
            raise _fail("INSDC041", f"cannot read GenBank source: {failure}") from failure
        return self.compile_bytes(raw)


def _exact_keys(value: Any, expected: set[str], label: str) -> dict[str, Any]:
    if type(value) is not dict:
        raise _fail("INSDC042", f"{label} must be an object")
    missing = sorted(expected - value.keys())
    extra = sorted(value.keys() - expected)
    if missing or extra:
        raise _fail("INSDC042", f"{label} keys invalid; missing={missing or 'none'}, unknown={extra or 'none'}")
    return value


def _exact_json_equal(left: Any, right: Any) -> bool:
    """Compare replay trees without Python's bool/int or int/float coercion."""

    pending = [(left, right)]
    while pending:
        left_item, right_item = pending.pop()
        if type(left_item) is not type(right_item):
            return False
        if type(left_item) is dict:
            if left_item.keys() != right_item.keys():
                return False
            pending.extend((left_item[key], right_item[key]) for key in left_item)
        elif type(left_item) is list:
            if len(left_item) != len(right_item):
                return False
            pending.extend(zip(left_item, right_item))
        elif left_item != right_item:
            return False
    return True


def _load_authority_manifest() -> dict[str, Any]:
    """Load and verify the exact offline standards snapshot bound by this profile."""

    try:
        raw = read_regular_file(
            _AUTHORITY_MANIFEST_PATH,
            maximum_bytes=_AUTHORITY_MANIFEST_BYTES,
            label="GenBank authority manifest",
        )
        manifest = load_json_object(
            raw,
            "GenBank authority manifest",
            maximum_bytes=_AUTHORITY_MANIFEST_BYTES,
            maximum_depth=16,
            maximum_members=256,
            maximum_string_bytes=4096,
        )
    except BoundedIOError as failure:
        raise _fail("INSDC051", f"cannot load the authority manifest: {failure}") from failure
    root = _exact_keys(
        manifest,
        {"format", "version", "profile", "sources", "artifact_sha256"},
        "authority manifest",
    )
    claimed = root["artifact_sha256"]
    core = {key: value for key, value in root.items() if key != "artifact_sha256"}
    if (
        claimed != AUTHORITY_MANIFEST_SHA256
        or digest(core) != AUTHORITY_MANIFEST_SHA256
    ):
        raise _fail("INSDC051", "authority manifest does not match the profile content pin")
    if (
        root["format"] != "brainc.bio.authority-manifest"
        or type(root["version"]) is not int
        or root["version"] != 1
        or root["profile"] != PROFILE
        or type(root["sources"]) is not dict
    ):
        raise _fail("INSDC051", "authority manifest identity is invalid")
    return copy.deepcopy(root)


def validate_genbank_artifact(
    payload: dict[str, Any],
    *,
    genbank_source: bytes | None = None,
) -> dict[str, Any]:
    """Validate the complete artifact, replaying its preserved GenBank bytes."""

    _load_authority_manifest()
    try:
        validate_json_tree(
            payload,
            "GenBank artifact",
            maximum_depth=MAX_JSON_DEPTH,
            maximum_members=MAX_JSON_MEMBERS,
            maximum_string_bytes=MAX_TEXT_BYTES,
        )
    except BoundedIOError as failure:
        raise _fail(
            "INSDC042",
            f"GenBank artifact is outside JSON limits: {failure}",
        ) from failure

    root = _exact_keys(
        payload,
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
        "GenBank artifact",
    )
    if root["format"] != FORMAT or type(root["version"]) is not int or root["version"] != VERSION:
        raise _fail("INSDC044", "unsupported GenBank artifact format or version")
    if root["profile"] != PROFILE or root["authority"] != AUTHORITY or root["compiler"] != COMPILER:
        raise _fail("INSDC044", "unsupported GenBank profile, authority, or compiler identity")
    try:
        artifact_digest(root)
    except ContractError as failure:
        raise _fail("INSDC045", str(failure)) from failure
    bio_ir = _exact_keys(
        root["bio_ir"],
        {
            "profile",
            "authority",
            "sequence_collection_artifact_sha256",
            "coordinate_system",
            "records",
        },
        "bio_ir",
    )
    if root["bio_ir_sha256"] != digest(bio_ir):
        raise _fail("INSDC047", "bio_ir_sha256 does not match Bio IR")
    if bio_ir["profile"] != PROFILE or bio_ir["authority"] != AUTHORITY:
        raise _fail("INSDC047", "Bio IR profile or authority is invalid")
    if bio_ir["coordinate_system"] != {
        "source": "1-based-closed",
        "normalized": "0-based-half-open",
    }:
        raise _fail("INSDC047", "Bio IR coordinate declaration is invalid")
    try:
        supplied_collection = as_sequence_collection_v2(root["sequence_collection"])
    except (SequenceCollectionV2Error, TypeError, ValueError) as failure:
        raise _fail("INSDC048", f"embedded Sequence Collection is invalid: {failure}") from failure
    if bio_ir["sequence_collection_artifact_sha256"] != supplied_collection.digest:
        raise _fail("INSDC048", "Bio IR does not bind the embedded Sequence Collection")
    reconstructed = supplied_collection.source_bytes()
    if genbank_source is not None:
        if type(genbank_source) is not bytes or len(genbank_source) > MAX_GENBANK_BYTES:
            raise _fail("INSDC046", "supplied GenBank source is not bounded bytes")
        if genbank_source != reconstructed:
            raise _fail("INSDC046", "supplied GenBank source does not match embedded source")
    replay_collection, replay_ir = _compile_parts(reconstructed)
    if not _exact_json_equal(replay_ir, bio_ir):
        raise _fail("INSDC049", "Bio IR does not match GenBank frontend replay")
    if not _exact_json_equal(
        replay_collection.to_dict(),
        supplied_collection.to_dict(),
    ):
        raise _fail("INSDC048", "embedded Sequence Collection does not match ORIGIN replay")
    _artifact_bytes(root)
    return root


def load_genbank_artifact(
    path: str | Path,
    *,
    genbank_source: bytes | None = None,
) -> GenBankArtifact:
    try:
        raw = read_regular_file(path, maximum_bytes=MAX_OUTPUT_BYTES, label="GenBank artifact")
        payload = load_json_object(
            raw,
            "GenBank artifact",
            maximum_bytes=MAX_OUTPUT_BYTES,
            maximum_depth=MAX_JSON_DEPTH,
            maximum_members=MAX_JSON_MEMBERS,
            maximum_string_bytes=MAX_TEXT_BYTES,
        )
    except BoundedIOError as failure:
        raise _fail("INSDC050", f"cannot load GenBank artifact: {failure}") from failure
    validate_genbank_artifact(payload, genbank_source=genbank_source)
    return GenBankArtifact(payload)


__all__ = [
    "AUTHORITY",
    "AUTHORITY_MANIFEST_SHA256",
    "COMPILER",
    "FORMAT",
    "PROFILE",
    "VERSION",
    "GenBankArtifact",
    "GenBankCompiler",
    "GenBankError",
    "load_genbank_artifact",
    "validate_genbank_artifact",
]

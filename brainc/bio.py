"""A bounded, target-neutral GFF3 external-sequence profile.

The frontend deliberately compiles syntax and relationships, not biological
ontology claims.  Feature types and application-defined lowercase attributes
remain opaque until a separately versioned ontology/policy stage interprets
them.  Sequence authority is the separately bound Sequence Collection IR, so
the optional GFF3 ``##FASTA`` payload form is rejected by this profile.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
import hashlib
import json
from pathlib import Path
import re
from typing import Any

from ._canonical import ContractError, SAFE_INTEGER, artifact_digest, canonical_bytes, digest
from ._io import BoundedIOError, atomic_write_file, read_regular_file
from .sequence_collection import (
    FORMAT as COLLECTION_FORMAT,
    VERSION as COLLECTION_VERSION,
    SequenceCollectionArtifact,
    SequenceCollectionError,
    _as_collection,
)


FORMAT = "brainc.bio.gff3-ir"
VERSION = 1
COMPILER = {
    "name": "brainc-bio-gff3",
    "version": "0.6.0",
    "passes": [
        "bind-sequence-collection",
        "parse-gff3",
        "normalize-coordinates",
        "resolve-feature-identities",
        "resolve-relationships",
        "emit-gff3-ir",
    ],
}
SOURCE_COORDINATES = "1-based-closed"
NORMALIZED_COORDINATES = "0-based-half-open"

# These ceilings are contract, not tuning hints.  They bound parser and emitted
# artifact work independently of any particular genome or annotation source.
MAX_GFF3_BYTES = 64 * 1024 * 1024
MAX_OUTPUT_BYTES = 256 * 1024 * 1024
MAX_LINE_BYTES = 1024 * 1024
MAX_LINES = 1_000_000
MAX_FEATURE_ROWS = 250_000
MAX_SEQUENCE_RECORDS = 100_000
MAX_DIRECTIVES = 100_000
MAX_ATTRIBUTES_PER_ROW = 256
MAX_VALUES_PER_ATTRIBUTE = 10_000
MAX_RELATIONSHIPS = 1_000_000
MAX_TEXT_BYTES = 1024 * 1024
MAX_JSON_DEPTH = 32
MAX_JSON_MEMBERS = 5_000_000

_GFF_VERSION_RE = re.compile(r"3(?:\.[0-9]+){0,2}\Z")
_SCORE_RE = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z")
_POSITIVE_DECIMAL_RE = re.compile(r"[0-9]+\Z")
_DIRECTIVE_RE = re.compile(r"[a-z][a-z0-9-]*\Z")
_ATTRIBUTE_TAG_RE = re.compile(r"[A-Za-z][A-Za-z0-9_.:-]*\Z")
_SAFE_SEQID = frozenset("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.:^*$@!+_?-|")
_HEX = frozenset("0123456789abcdefABCDEF")
_GENERAL_ESCAPE_BYTES = frozenset(range(0x20)) | frozenset({0x25, 0x7F})
_ATTRIBUTE_ESCAPE_BYTES = _GENERAL_ESCAPE_BYTES | frozenset(map(ord, ";=&,"))
_RESERVED_ATTRIBUTES = frozenset(
    {
        "ID",
        "Name",
        "Alias",
        "Parent",
        "Target",
        "Gap",
        "Derives_from",
        "Note",
        "Dbxref",
        "Ontology_term",
        "Is_circular",
    }
)
_CANONICAL_MULTI_VALUE_ATTRIBUTES = frozenset(
    {"Parent", "Alias", "Note", "Dbxref", "Ontology_term", "Derives_from"}
)


class GFF3Error(ContractError):
    """A source or artifact violates the closed GFF3 frontend contract."""


def _fail(code: str, detail: str) -> GFF3Error:
    return GFF3Error(f"{code}: {detail}")


def _keys(value: Any, expected: set[str], label: str) -> dict[str, Any]:
    if type(value) is not dict:
        raise _fail("BIO001", f"{label} must be an object")
    missing = sorted(expected - value.keys())
    extra = sorted(value.keys() - expected)
    if missing or extra:
        raise _fail(
            "BIO002",
            f"{label} keys invalid; missing={missing or 'none'}, unknown={extra or 'none'}",
        )
    return value


def _integer(value: Any, label: str, *, minimum: int = 0, maximum: int = SAFE_INTEGER) -> int:
    if type(value) is not int or value < minimum or value > maximum:
        raise _fail("BIO003", f"{label} must be an integer in [{minimum}, {maximum}]")
    return value


def _text(value: Any, label: str, *, empty: bool = False) -> str:
    if type(value) is not str or (not empty and not value):
        raise _fail("BIO004", f"{label} must be {'possibly empty ' if empty else 'nonempty '}text")
    if any(0xD800 <= ord(character) <= 0xDFFF for character in value):
        raise _fail("BIO004", f"{label} contains a lone Unicode surrogate")
    if len(value.encode("utf-8")) > MAX_TEXT_BYTES:
        raise _fail("BIO004", f"{label} exceeds {MAX_TEXT_BYTES} UTF-8 bytes")
    return value


def _sorted_text(values: list[str], label: str) -> None:
    wanted = sorted(values, key=lambda item: item.encode("utf-8"))
    if values != wanted:
        raise _fail("BIO005", f"{label} must be sorted by UTF-8 bytes")
    if len(values) != len(set(values)):
        raise _fail("BIO005", f"{label} contains duplicate values")


def _validate_json_tree(value: Any) -> None:
    members = 0
    stack: list[tuple[Any, int]] = [(value, 0)]
    while stack:
        current, depth = stack.pop()
        if depth > MAX_JSON_DEPTH:
            raise _fail("BIO065", f"artifact exceeds JSON depth {MAX_JSON_DEPTH}")
        if current is None or type(current) is bool:
            continue
        if type(current) is int:
            if not -SAFE_INTEGER <= current <= SAFE_INTEGER:
                raise _fail("BIO065", "artifact contains an unsafe JSON integer")
            continue
        if type(current) is float:
            raise _fail("BIO065", "artifact contains an unexpected JSON number")
        if type(current) is str:
            _text(current, "artifact string", empty=True)
            continue
        if type(current) is list:
            members += len(current)
            stack.extend((item, depth + 1) for item in reversed(current))
        elif type(current) is dict:
            members += len(current)
            for key, item in reversed(list(current.items())):
                if type(key) is not str:
                    raise _fail("BIO065", "artifact contains a non-string object key")
                stack.append((key, depth + 1))
                stack.append((item, depth + 1))
        else:
            raise _fail("BIO065", f"artifact contains unsupported type {type(current).__name__}")
        if members > MAX_JSON_MEMBERS:
            raise _fail("BIO065", f"artifact exceeds {MAX_JSON_MEMBERS} JSON members")


def _read_regular(path: str | Path, maximum: int, label: str) -> bytes:
    try:
        return read_regular_file(path, maximum_bytes=maximum, label=label)
    except BoundedIOError as failure:
        raise _fail("BIO041", str(failure)) from failure


def _atomic_write(path: str | Path, raw: bytes) -> None:
    try:
        atomic_write_file(
            path,
            raw,
            maximum_bytes=MAX_OUTPUT_BYTES,
            label="GFF3 artifact output",
        )
    except BoundedIOError as failure:
        raise _fail("BIO040", str(failure)) from failure


def _pretty_wire_size(payload: dict[str, Any]) -> int:
    """Count the exact UTF-8 bytes emitted by :meth:`GFF3Artifact.save`."""

    total = 1  # The writer terminates JSON with one newline.
    try:
        encoder = json.JSONEncoder(
            indent=2,
            sort_keys=True,
            ensure_ascii=False,
            allow_nan=False,
        )
        for chunk in encoder.iterencode(payload):
            total += len(chunk.encode("utf-8"))
            if total > MAX_OUTPUT_BYTES:
                return total
    except (TypeError, ValueError, UnicodeError) as failure:
        raise _fail("BIO040", f"artifact cannot be serialized: {failure}") from failure
    return total


def _preflight_pretty_wire(payload: dict[str, Any]) -> None:
    if _pretty_wire_size(payload) > MAX_OUTPUT_BYTES:
        raise _fail(
            "BIO040",
            f"serialized artifact exceeds {MAX_OUTPUT_BYTES} bytes",
        )


def _percent_decode(
    value: str,
    label: str,
    *,
    allowed_escape_bytes: frozenset[int] | None,
) -> str:
    encoded = bytearray()
    index = 0
    while index < len(value):
        character = value[index]
        if character == "%":
            if index + 2 >= len(value) or value[index + 1] not in _HEX or value[index + 2] not in _HEX:
                raise _fail("BIO006", f"{label} contains a malformed percent escape")
            octet = int(value[index + 1 : index + 3], 16)
            if allowed_escape_bytes is not None and octet not in allowed_escape_bytes:
                raise _fail("BIO006", f"{label} percent-encodes a character that must be literal")
            encoded.append(octet)
            index += 3
            continue
        try:
            encoded.extend(character.encode("utf-8"))
        except UnicodeEncodeError as failure:
            raise _fail("BIO006", f"{label} is not valid UTF-8 text") from failure
        index += 1
    try:
        decoded = encoded.decode("utf-8")
    except UnicodeDecodeError as failure:
        raise _fail("BIO006", f"{label} percent escapes do not decode as UTF-8") from failure
    return _text(decoded, label)


def _decode_seqid(value: str, label: str) -> str:
    if not value or value == ".":
        raise _fail("BIO007", f"{label} must identify a sequence record")
    index = 0
    while index < len(value):
        if value[index] == "%":
            if index + 2 >= len(value) or value[index + 1] not in _HEX or value[index + 2] not in _HEX:
                raise _fail("BIO006", f"{label} contains a malformed percent escape")
            octet = int(value[index + 1 : index + 3], 16)
            if octet < 128 and chr(octet) in _SAFE_SEQID:
                raise _fail("BIO006", f"{label} redundantly percent-encodes a permitted character")
            index += 3
        elif value[index] not in _SAFE_SEQID:
            raise _fail("BIO007", f"{label} contains an unescaped character {value[index]!r}")
        else:
            index += 1
    decoded = _percent_decode(value, label, allowed_escape_bytes=None)
    if any(
        character.isspace()
        or ord(character) < 0x20
        or 0x7F <= ord(character) <= 0x9F
        for character in decoded
    ):
        raise _fail("BIO007", f"{label} decodes to whitespace or a control character")
    return decoded


def _collection(value: SequenceCollectionArtifact | dict[str, Any]) -> SequenceCollectionArtifact:
    try:
        if isinstance(value, SequenceCollectionArtifact):
            return _as_collection(value.to_dict())
        if type(value) is dict:
            return _as_collection(value)
    except SequenceCollectionError as failure:
        raise _fail("BIO008", f"invalid sequence collection: {failure}") from failure
    raise _fail("BIO008", "sequence collection must be a validated artifact or artifact object")


def _record_bindings(collection: SequenceCollectionArtifact) -> list[dict[str, Any]]:
    if len(collection.members) > MAX_SEQUENCE_RECORDS:
        raise _fail("BIO009", f"sequence collection exceeds {MAX_SEQUENCE_RECORDS} records")
    bindings = [
        {
            "seqid": member.artifact.record_id,
            "length": len(member.artifact.sequence),
            "sequence_sha256": member.artifact.sequence_sha256,
            "sequence_artifact_sha256": member.artifact.digest,
        }
        for member in collection.members
    ]
    return sorted(bindings, key=lambda item: item["seqid"].encode("utf-8"))


def _parse_score(value: str, line: int) -> str | None:
    if value == ".":
        return None
    if _SCORE_RE.fullmatch(value) is None:
        raise _fail("BIO010", f"line {line} score is not a decimal literal")
    try:
        score = Decimal(value)
    except InvalidOperation as failure:
        raise _fail("BIO010", f"line {line} score is invalid") from failure
    if not score.is_finite():
        raise _fail("BIO010", f"line {line} score must be finite")
    _text(value, f"line {line} score")
    return value


def _parse_coordinate(value: str, line: int, column: str) -> int:
    if _POSITIVE_DECIMAL_RE.fullmatch(value) is None:
        raise _fail(
            "BIO033",
            f"line {line} {column} must be an ASCII positive-decimal integer",
        )
    try:
        parsed = int(value)
    except ValueError as failure:
        raise _fail("BIO033", f"line {line} {column} is too large to parse") from failure
    return _integer(parsed, f"line {line} {column}", minimum=1)


def _parse_attributes(value: str, line: int) -> list[dict[str, Any]]:
    if value == ".":
        return []
    if not value:
        raise _fail("BIO011", f"line {line} attributes cannot be empty")
    fields = value.split(";")
    if fields[-1] == "":
        fields.pop()
    if not fields or any(field == "" for field in fields):
        raise _fail("BIO011", f"line {line} contains an empty attribute")
    if len(fields) > MAX_ATTRIBUTES_PER_ROW:
        raise _fail("BIO011", f"line {line} exceeds {MAX_ATTRIBUTES_PER_ROW} attributes")
    attributes: list[dict[str, Any]] = []
    seen: set[str] = set()
    for field in fields:
        if field.count("=") != 1:
            raise _fail("BIO012", f"line {line} attribute must contain exactly one unescaped '='")
        raw_tag, raw_values = field.split("=", 1)
        tag = _percent_decode(
            raw_tag,
            f"line {line} attribute tag",
            allowed_escape_bytes=_ATTRIBUTE_ESCAPE_BYTES,
        )
        if _ATTRIBUTE_TAG_RE.fullmatch(tag) is None:
            raise _fail("BIO012", f"line {line} attribute tag {tag!r} is not a GFF3 name")
        if tag in seen:
            raise _fail("BIO013", f"line {line} repeats attribute {tag!r}")
        seen.add(tag)
        if tag[0].isupper() and tag not in _RESERVED_ATTRIBUTES:
            raise _fail("BIO014", f"line {line} uses unsupported reserved attribute {tag!r}")
        value_parts = raw_values.split(",")
        if not value_parts or any(part == "" for part in value_parts):
            raise _fail("BIO015", f"line {line} attribute {tag!r} has an empty value")
        if len(value_parts) > MAX_VALUES_PER_ATTRIBUTE:
            raise _fail("BIO015", f"line {line} attribute {tag!r} has too many values")
        if any("&" in part for part in value_parts):
            raise _fail("BIO015", f"line {line} attribute {tag!r} contains unescaped '&'")
        values = [
            _percent_decode(
                part,
                f"line {line} attribute {tag!r}",
                allowed_escape_bytes=_ATTRIBUTE_ESCAPE_BYTES,
            )
            for part in value_parts
        ]
        if tag in _CANONICAL_MULTI_VALUE_ATTRIBUTES:
            values.sort(key=lambda item: item.encode("utf-8"))
            if len(values) != len(set(values)):
                raise _fail("BIO015", f"line {line} attribute {tag!r} repeats a value")
        if tag in {"ID", "Is_circular", "Target", "Gap"} and len(values) != 1:
            raise _fail("BIO015", f"line {line} attribute {tag!r} requires one value")
        if tag == "Is_circular" and values[0] not in {"true", "false"}:
            raise _fail("BIO016", f"line {line} Is_circular must be 'true' or 'false'")
        attributes.append({"tag": tag, "values": values})
    attributes.sort(key=lambda item: item["tag"].encode("utf-8"))
    return attributes


def _attribute(attributes: list[dict[str, Any]], tag: str) -> list[str]:
    for attribute in attributes:
        if attribute["tag"] == tag:
            return attribute["values"]
    return []


def _intervals(start: int, end: int, length: int, circular: bool, label: str) -> list[dict[str, int]]:
    if start > length:
        raise _fail("BIO017", f"{label} starts beyond sequence length {length}")
    span = end - start + 1
    if end <= length:
        return [{"start": start - 1, "end": end}]
    if not circular:
        raise _fail("BIO018", f"{label} crosses the origin without Is_circular=true")
    if span > length:
        raise _fail("BIO018", f"{label} spans more than one circular sequence traversal")
    first_start = start - 1
    overflow = first_start + span - length
    return [{"start": first_start, "end": length}, {"start": 0, "end": overflow}]


def _row_signature(row: dict[str, Any]) -> dict[str, Any]:
    return {
        "seqid": row["seqid"],
        "source": row["source"],
        "type": row["type"],
        "strand": row["strand"],
        "segment": row["segment"],
    }


def _anonymous_id(row: dict[str, Any]) -> str:
    return "anon:" + digest(_row_signature(row))


def _cycle_check(relationships: list[dict[str, str]]) -> None:
    graph: dict[str, set[str]] = {}
    for relationship in relationships:
        graph.setdefault(relationship["source"], set()).add(relationship["target"])
    state: dict[str, int] = {}
    for node in sorted(graph, key=lambda item: item.encode("utf-8")):
        if state.get(node) == 2:
            continue
        state[node] = 1
        stack: list[tuple[str, list[str], int]] = [
            (node, sorted(graph.get(node, ()), key=lambda item: item.encode("utf-8")), 0)
        ]
        while stack:
            current, targets, index = stack[-1]
            if index == len(targets):
                state[current] = 2
                stack.pop()
                continue
            target = targets[index]
            stack[-1] = (current, targets, index + 1)
            marker = state.get(target, 0)
            if marker == 1:
                raise _fail("BIO019", f"feature relationship cycle includes {target!r}")
            if marker == 2:
                continue
            state[target] = 1
            stack.append(
                (
                    target,
                    sorted(graph.get(target, ()), key=lambda item: item.encode("utf-8")),
                    0,
                )
            )


def _parse(raw: bytes, collection: SequenceCollectionArtifact) -> dict[str, Any]:
    if type(raw) is not bytes:
        raise _fail("BIO020", "GFF3 source must be bytes")
    if len(raw) > MAX_GFF3_BYTES:
        raise _fail("BIO020", f"GFF3 source exceeds {MAX_GFF3_BYTES} bytes")
    try:
        source = raw.decode("utf-8")
    except UnicodeDecodeError as failure:
        raise _fail("BIO021", "GFF3 source must be strict UTF-8") from failure
    if any(
        (ord(character) < 0x20 and character not in "\t\n\r")
        or 0x7F <= ord(character) <= 0x9F
        for character in source
    ):
        raise _fail("BIO021", "GFF3 source contains an unescaped control character")
    lines = source.split("\n")
    lines = [line[:-1] if line.endswith("\r") else line for line in lines]
    if any("\r" in line for line in lines):
        raise _fail("BIO021", "GFF3 source contains a carriage return outside CRLF")
    if not lines:
        raise _fail("BIO022", "GFF3 source is empty")
    if len(lines) > MAX_LINES:
        raise _fail("BIO022", f"GFF3 source exceeds {MAX_LINES} lines")
    if not lines[0].startswith("##gff-version "):
        raise _fail("BIO023", "##gff-version must be the first physical line")
    version = lines[0][len("##gff-version ") :]
    if _GFF_VERSION_RE.fullmatch(version) is None:
        raise _fail("BIO023", "GFF version must be 3, 3.x, or 3.x.y")

    records = _record_bindings(collection)
    lengths = {record["seqid"]: record["length"] for record in records}
    sequence_regions: list[dict[str, Any]] = []
    region_by_seqid: dict[str, tuple[int, int]] = {}
    directives: list[dict[str, Any]] = []
    rows: list[dict[str, Any]] = []

    for line_number, line in enumerate(lines[1:], start=2):
        if len(line.encode("utf-8")) > MAX_LINE_BYTES:
            raise _fail("BIO024", f"line {line_number} exceeds {MAX_LINE_BYTES} bytes")
        if not line:
            continue
        if line == "##FASTA" or line.startswith(">"):
            raise _fail(
                "BIO066",
                "embedded FASTA is outside the external-sequence GFF3 profile",
            )
        if line.startswith("###"):
            if line != "###":
                raise _fail("BIO025", f"line {line_number} has content after feature boundary")
            continue
        if line.startswith("##"):
            body = line[2:]
            name, separator, value = body.partition(" ")
            if not name or (separator and (not value or value != value.strip())):
                raise _fail("BIO026", f"line {line_number} has malformed directive syntax")
            if name == "gff-version":
                raise _fail("BIO023", "GFF3 source contains more than one version directive")
            if name == "sequence-region":
                parts = value.split()
                if len(parts) != 3:
                    raise _fail("BIO027", f"line {line_number} sequence-region needs seqid start end")
                seqid = _decode_seqid(parts[0], f"line {line_number} sequence-region seqid")
                try:
                    start = _parse_coordinate(parts[1], line_number, "sequence-region start")
                    end = _parse_coordinate(parts[2], line_number, "sequence-region end")
                except GFF3Error as failure:
                    raise _fail("BIO027", str(failure)) from failure
                if end < start:
                    raise _fail("BIO027", f"line {line_number} sequence-region end precedes start")
                if seqid not in lengths:
                    raise _fail("BIO028", f"line {line_number} references unknown sequence {seqid!r}")
                if end > lengths[seqid]:
                    raise _fail("BIO028", f"line {line_number} sequence-region exceeds sequence length")
                if seqid in region_by_seqid:
                    raise _fail("BIO027", f"duplicate sequence-region for {seqid!r}")
                region_by_seqid[seqid] = (start, end)
                sequence_regions.append(
                    {
                        "line": line_number,
                        "seqid": seqid,
                        "raw_start": start,
                        "raw_end": end,
                        "start": start - 1,
                        "end": end,
                    }
                )
                continue
            if _DIRECTIVE_RE.fullmatch(name) is None:
                raise _fail("BIO029", f"line {line_number} uses unsupported reserved directive {name!r}")
            if len(directives) >= MAX_DIRECTIVES:
                raise _fail("BIO029", f"source exceeds {MAX_DIRECTIVES} application directives")
            directives.append({"line": line_number, "name": name, "value": _text(value, f"line {line_number} directive value", empty=True)})
            continue
        if line.startswith("#"):
            continue
        if len(rows) >= MAX_FEATURE_ROWS:
            raise _fail("BIO030", f"source exceeds {MAX_FEATURE_ROWS} feature rows")
        columns = line.split("\t")
        if len(columns) != 9:
            raise _fail("BIO031", f"line {line_number} must contain exactly nine tab-delimited columns")
        raw_seqid, raw_source, raw_type, raw_start, raw_end, raw_score, strand, raw_phase, raw_attributes = columns
        seqid = _decode_seqid(raw_seqid, f"line {line_number} seqid")
        if seqid not in lengths:
            raise _fail("BIO028", f"line {line_number} references unknown sequence {seqid!r}")
        source_value = (
            None
            if raw_source == "."
            else _percent_decode(
                raw_source,
                f"line {line_number} source",
                allowed_escape_bytes=_GENERAL_ESCAPE_BYTES,
            )
        )
        if raw_type == ".":
            raise _fail("BIO032", f"line {line_number} type cannot be '.'")
        feature_type = _percent_decode(
            raw_type,
            f"line {line_number} type",
            allowed_escape_bytes=_GENERAL_ESCAPE_BYTES,
        )
        start = _parse_coordinate(raw_start, line_number, "start")
        end = _parse_coordinate(raw_end, line_number, "end")
        if end < start:
            raise _fail("BIO033", f"line {line_number} end precedes start")
        if strand not in {"+", "-", ".", "?"}:
            raise _fail("BIO034", f"line {line_number} strand must be +, -, ., or ?")
        if feature_type == "CDS":
            if raw_phase not in {"0", "1", "2"}:
                raise _fail("BIO035", f"line {line_number} CDS phase must be 0, 1, or 2")
            phase: int | None = int(raw_phase)
        else:
            if raw_phase != ".":
                raise _fail("BIO035", f"line {line_number} non-CDS phase must be '.'")
            phase = None
        attributes = _parse_attributes(raw_attributes, line_number)
        rows.append(
            {
                "seqid": seqid,
                "source": source_value,
                "type": feature_type,
                "strand": strand,
                "raw_start": start,
                "raw_end": end,
                "score": _parse_score(raw_score, line_number),
                "phase": phase,
                "attributes": attributes,
                "line": line_number,
            }
        )

    circular_records = {
        row["seqid"]
        for row in rows
        if _attribute(row["attributes"], "Is_circular") == ["true"]
        and row["raw_start"] == 1
        and row["raw_end"] == lengths[row["seqid"]]
    }
    for row in rows:
        intervals = _intervals(
            row["raw_start"],
            row["raw_end"],
            lengths[row["seqid"]],
            row["seqid"] in circular_records,
            f"line {row['line']}",
        )
        region = region_by_seqid.get(row["seqid"])
        if region is not None:
            region_start, region_end = region
            for interval in intervals:
                if interval["start"] < region_start - 1 or interval["end"] > region_end:
                    raise _fail("BIO036", f"line {row['line']} lies outside its sequence-region")
        row["segment"] = {
            "line": row.pop("line"),
            "raw_start": row.pop("raw_start"),
            "raw_end": row.pop("raw_end"),
            "intervals": intervals,
            "score": row.pop("score"),
            "phase": row.pop("phase"),
            "attributes": row.pop("attributes"),
        }

    entities: dict[str, dict[str, Any]] = {}
    for row in rows:
        declared_values = _attribute(row["segment"]["attributes"], "ID")
        declared_id = declared_values[0] if declared_values else None
        entity_id = "id:" + declared_id if declared_id is not None else _anonymous_id(row)
        entity = entities.get(entity_id)
        if entity is None:
            entities[entity_id] = {
                "entity_id": entity_id,
                "declared_id": declared_id,
                "seqid": row["seqid"],
                "source": row["source"],
                "type": row["type"],
                "strand": row["strand"],
                "segments": [row["segment"]],
            }
            continue
        if declared_id is None:
            raise _fail("BIO037", f"anonymous feature identity collision at line {row['segment']['line']}")
        common = ("seqid", "source", "type", "strand")
        if any(entity[key] != row[key] for key in common):
            raise _fail("BIO038", f"discontinuous ID {declared_id!r} changes core feature fields")
        if any(segment["line"] == row["segment"]["line"] for segment in entity["segments"]):
            raise _fail("BIO038", f"discontinuous ID {declared_id!r} repeats a source line")
        entity["segments"].append(row["segment"])

    features = sorted(entities.values(), key=lambda item: item["entity_id"].encode("utf-8"))
    for entity in features:
        entity["segments"].sort(key=lambda item: item["line"])

    known_ids = set(entities)
    relationship_set: set[tuple[str, str, str]] = set()
    for entity in features:
        for segment in entity["segments"]:
            for tag, kind in (("Parent", "parent"), ("Derives_from", "derives-from")):
                for target_id in _attribute(segment["attributes"], tag):
                    target = "id:" + target_id
                    if target not in known_ids:
                        raise _fail("BIO039", f"line {segment['line']} has dangling {tag} {target_id!r}")
                    relationship_set.add((entity["entity_id"], kind, target))
                    if len(relationship_set) > MAX_RELATIONSHIPS:
                        raise _fail("BIO039", f"source exceeds {MAX_RELATIONSHIPS} relationships")
    relationships = [
        {"source": source, "kind": kind, "target": target}
        for source, kind, target in sorted(
            relationship_set,
            key=lambda item: tuple(part.encode("utf-8") for part in item),
        )
    ]
    _cycle_check(relationships)
    return {
        "gff3_version": version,
        "coordinate_system": {"source": SOURCE_COORDINATES, "normalized": NORMALIZED_COORDINATES},
        "records": records,
        "sequence_regions": sequence_regions,
        "directives": directives,
        "features": features,
        "relationships": relationships,
    }


def _input_binding(collection: SequenceCollectionArtifact, raw: bytes) -> dict[str, Any]:
    collection_payload = collection.to_dict()
    return {
        "gff3_sha256": hashlib.sha256(raw).hexdigest(),
        "sequence_collection": {
            "format": COLLECTION_FORMAT,
            "version": COLLECTION_VERSION,
            "artifact_sha256": collection_payload["artifact_sha256"],
            "collection_ir_sha256": collection_payload["collection_ir_sha256"],
        },
    }


def _seal(inputs: dict[str, Any], bio_ir: dict[str, Any]) -> dict[str, Any]:
    core = {
        "format": FORMAT,
        "version": VERSION,
        "compiler": copy.deepcopy(COMPILER),
        "inputs": inputs,
        "bio_ir": bio_ir,
        "bio_ir_sha256": digest(bio_ir),
    }
    payload = {**core, "artifact_sha256": digest(core)}
    if len(canonical_bytes(payload)) > MAX_OUTPUT_BYTES:
        raise _fail("BIO040", f"emitted artifact exceeds {MAX_OUTPUT_BYTES} canonical bytes")
    _preflight_pretty_wire(payload)
    return payload


@dataclass(frozen=True)
class GFF3Artifact:
    """An immutable-by-interface compiled GFF3 artifact."""

    _payload: dict[str, Any]

    @property
    def digest(self) -> str:
        return self._payload["artifact_sha256"]

    def to_dict(self) -> dict[str, Any]:
        return copy.deepcopy(self._payload)

    def save(self, path: str | Path) -> None:
        _preflight_pretty_wire(self._payload)
        raw = (
            json.dumps(self._payload, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False)
            + "\n"
        ).encode("utf-8")
        if len(raw) > MAX_OUTPUT_BYTES:
            raise _fail("BIO040", f"serialized artifact exceeds {MAX_OUTPUT_BYTES} bytes")
        _atomic_write(path, raw)


class GFF3Compiler:
    """Compile the bounded GFF3 profile against a validated sequence collection."""

    def compile_bytes(
        self,
        raw: bytes,
        sequence_collection: SequenceCollectionArtifact | dict[str, Any],
    ) -> GFF3Artifact:
        collection = _collection(sequence_collection)
        bio_ir = _parse(raw, collection)
        payload = _seal(_input_binding(collection, raw), bio_ir)
        validate_gff3_artifact(payload, collection, gff3_source=raw)
        return GFF3Artifact(payload)

    def compile_file(
        self,
        path: str | Path,
        sequence_collection: SequenceCollectionArtifact | dict[str, Any],
    ) -> GFF3Artifact:
        return self.compile_bytes(
            _read_regular(path, MAX_GFF3_BYTES, "GFF3 source"), sequence_collection
        )


def _validate_attribute_array(value: Any, label: str) -> list[dict[str, Any]]:
    if type(value) is not list or len(value) > MAX_ATTRIBUTES_PER_ROW:
        raise _fail("BIO042", f"{label} must be a bounded array")
    result: list[dict[str, Any]] = []
    tags: list[str] = []
    for index, item_value in enumerate(value):
        item = _keys(item_value, {"tag", "values"}, f"{label}[{index}]")
        tag = _text(item["tag"], f"{label}[{index}].tag")
        if _ATTRIBUTE_TAG_RE.fullmatch(tag) is None:
            raise _fail("BIO012", f"{label}[{index}] tag {tag!r} is not a GFF3 name")
        if tag[0].isupper() and tag not in _RESERVED_ATTRIBUTES:
            raise _fail("BIO014", f"{label}[{index}] uses unsupported reserved attribute {tag!r}")
        values = item["values"]
        if type(values) is not list or not values or len(values) > MAX_VALUES_PER_ATTRIBUTE:
            raise _fail("BIO042", f"{label}[{index}].values must be a bounded nonempty array")
        checked = [_text(entry, f"{label}[{index}].values") for entry in values]
        if tag in _CANONICAL_MULTI_VALUE_ATTRIBUTES:
            _sorted_text(checked, f"{label}[{index}].values")
        if tag in {"ID", "Is_circular", "Target", "Gap"} and len(checked) != 1:
            raise _fail("BIO015", f"attribute {tag!r} requires one value")
        if tag == "Is_circular" and checked[0] not in {"true", "false"}:
            raise _fail("BIO016", "Is_circular must be 'true' or 'false'")
        tags.append(tag)
        result.append({"tag": tag, "values": checked})
    _sorted_text(tags, f"{label} tags")
    return result


def validate_gff3_artifact(
    payload: dict[str, Any],
    sequence_collection: SequenceCollectionArtifact | dict[str, Any],
    *,
    gff3_source: bytes | None = None,
) -> dict[str, Any]:
    """Validate the full closed artifact and its external collection binding.

    Supplying the collection is intentional: a digest reference alone cannot
    establish that the referenced biological input is actually present.
    """

    _validate_json_tree(payload)
    collection = _collection(sequence_collection)
    root = _keys(
        payload,
        {"format", "version", "compiler", "inputs", "bio_ir", "bio_ir_sha256", "artifact_sha256"},
        "GFF3 artifact",
    )
    if root["format"] != FORMAT or type(root["version"]) is not int or root["version"] != VERSION:
        raise _fail("BIO043", "unsupported GFF3 artifact format or version")
    if root["compiler"] != COMPILER:
        raise _fail("BIO044", "unsupported GFF3 compiler identity")
    try:
        artifact_digest(root)
    except ContractError as failure:
        raise _fail("BIO045", str(failure)) from failure
    inputs = _keys(root["inputs"], {"gff3_sha256", "sequence_collection"}, "inputs")
    if type(inputs["gff3_sha256"]) is not str or re.fullmatch(r"[0-9a-f]{64}", inputs["gff3_sha256"]) is None:
        raise _fail("BIO046", "inputs.gff3_sha256 must be lowercase SHA-256")
    replayed_bio_ir: dict[str, Any] | None = None
    if gff3_source is not None:
        if type(gff3_source) is not bytes or len(gff3_source) > MAX_GFF3_BYTES:
            raise _fail("BIO046", "supplied GFF3 source is not bounded bytes")
        if hashlib.sha256(gff3_source).hexdigest() != inputs["gff3_sha256"]:
            raise _fail("BIO046", "inputs.gff3_sha256 does not match supplied GFF3 source")
        replayed_bio_ir = _parse(gff3_source, collection)
    binding = _keys(
        inputs["sequence_collection"],
        {"format", "version", "artifact_sha256", "collection_ir_sha256"},
        "inputs.sequence_collection",
    )
    expected_collection = collection.to_dict()
    if binding != {
        "format": COLLECTION_FORMAT,
        "version": COLLECTION_VERSION,
        "artifact_sha256": expected_collection["artifact_sha256"],
        "collection_ir_sha256": expected_collection["collection_ir_sha256"],
    }:
        raise _fail("BIO047", "sequence collection binding does not match supplied artifact")
    bio_ir = _keys(
        root["bio_ir"],
        {"gff3_version", "coordinate_system", "records", "sequence_regions", "directives", "features", "relationships"},
        "bio_ir",
    )
    if type(root["bio_ir_sha256"]) is not str or root["bio_ir_sha256"] != digest(bio_ir):
        raise _fail("BIO048", "bio_ir_sha256 does not match canonical Bio IR")
    if replayed_bio_ir is not None and bio_ir != replayed_bio_ir:
        raise _fail("BIO048", "Bio IR does not match independent GFF3 frontend replay")
    if type(bio_ir["gff3_version"]) is not str or _GFF_VERSION_RE.fullmatch(bio_ir["gff3_version"]) is None:
        raise _fail("BIO049", "bio_ir.gff3_version is not GFF3")
    coordinates = _keys(bio_ir["coordinate_system"], {"source", "normalized"}, "coordinate_system")
    if coordinates != {"source": SOURCE_COORDINATES, "normalized": NORMALIZED_COORDINATES}:
        raise _fail("BIO050", "coordinate system declaration is invalid")
    if bio_ir["records"] != _record_bindings(collection):
        raise _fail("BIO051", "Bio IR record bindings do not match sequence collection")
    lengths = {record["seqid"]: record["length"] for record in bio_ir["records"]}

    regions = bio_ir["sequence_regions"]
    if type(regions) is not list:
        raise _fail("BIO052", "sequence_regions must be an array")
    region_map: dict[str, tuple[int, int]] = {}
    occupied_source_lines = {1}
    last_region_line = 0
    for index, region_value in enumerate(regions):
        region = _keys(
            region_value,
            {"line", "seqid", "raw_start", "raw_end", "start", "end"},
            f"sequence_regions[{index}]",
        )
        line = _integer(
            region["line"], f"sequence_regions[{index}].line", minimum=2, maximum=MAX_LINES
        )
        if line <= last_region_line:
            raise _fail("BIO052", "sequence_regions must be in source-line order")
        if line in occupied_source_lines:
            raise _fail("BIO052", "sequence-region source lines must be unique")
        occupied_source_lines.add(line)
        last_region_line = line
        seqid = _text(region["seqid"], f"sequence_regions[{index}].seqid")
        if seqid not in lengths or seqid in region_map:
            raise _fail("BIO052", f"invalid or duplicate sequence-region {seqid!r}")
        raw_start = _integer(region["raw_start"], f"sequence_regions[{index}].raw_start", minimum=1)
        raw_end = _integer(region["raw_end"], f"sequence_regions[{index}].raw_end", minimum=1)
        if raw_end < raw_start or raw_end > lengths[seqid]:
            raise _fail("BIO052", f"sequence-region {seqid!r} has invalid bounds")
        if region["start"] != raw_start - 1 or region["end"] != raw_end:
            raise _fail("BIO052", f"sequence-region {seqid!r} has invalid normalization")
        region_map[seqid] = (raw_start, raw_end)

    directives = bio_ir["directives"]
    if type(directives) is not list or len(directives) > MAX_DIRECTIVES:
        raise _fail("BIO053", "directives must be a bounded array")
    last_directive_line = 0
    for index, directive_value in enumerate(directives):
        directive = _keys(directive_value, {"line", "name", "value"}, f"directives[{index}]")
        line = _integer(
            directive["line"], f"directives[{index}].line", minimum=2, maximum=MAX_LINES
        )
        if line <= last_directive_line:
            raise _fail("BIO053", "directives must be in source-line order")
        if line in occupied_source_lines:
            raise _fail("BIO053", "directive source lines must be unique")
        occupied_source_lines.add(line)
        last_directive_line = line
        name = _text(directive["name"], f"directives[{index}].name")
        if _DIRECTIVE_RE.fullmatch(name) is None or name in {"gff-version", "sequence-region"}:
            raise _fail("BIO053", f"unsupported opaque directive {name!r}")
        _text(directive["value"], f"directives[{index}].value", empty=True)

    features = bio_ir["features"]
    if type(features) is not list or len(features) > MAX_FEATURE_ROWS:
        raise _fail("BIO054", "features must be a bounded array")
    feature_ids: list[str] = []
    feature_map: dict[str, dict[str, Any]] = {}
    source_lines: set[int] = set()
    total_segments = 0
    for feature_index, feature_value in enumerate(features):
        feature = _keys(
            feature_value,
            {"entity_id", "declared_id", "seqid", "source", "type", "strand", "segments"},
            f"features[{feature_index}]",
        )
        entity_id = _text(feature["entity_id"], f"features[{feature_index}].entity_id")
        declared_id = feature["declared_id"]
        if declared_id is not None:
            declared_id = _text(declared_id, f"features[{feature_index}].declared_id")
            if entity_id != "id:" + declared_id:
                raise _fail("BIO055", f"feature {entity_id!r} has inconsistent declared ID")
        elif not entity_id.startswith("anon:"):
            raise _fail("BIO055", f"anonymous feature {entity_id!r} has invalid identity")
        if entity_id in feature_map:
            raise _fail("BIO055", f"duplicate feature identity {entity_id!r}")
        feature_ids.append(entity_id)
        feature_map[entity_id] = feature
        seqid = _text(feature["seqid"], f"features[{feature_index}].seqid")
        if seqid not in lengths:
            raise _fail("BIO056", f"feature {entity_id!r} references unknown sequence")
        if feature["source"] is not None:
            _text(feature["source"], f"features[{feature_index}].source")
        feature_type = _text(feature["type"], f"features[{feature_index}].type")
        if feature["strand"] not in {"+", "-", ".", "?"}:
            raise _fail("BIO057", f"feature {entity_id!r} has invalid strand")
        segments = feature["segments"]
        if type(segments) is not list or not segments:
            raise _fail("BIO058", f"feature {entity_id!r} needs at least one segment")
        if declared_id is None and len(segments) != 1:
            raise _fail("BIO058", f"anonymous feature {entity_id!r} cannot be discontinuous")
        last_line = 0
        for segment_index, segment_value in enumerate(segments):
            segment = _keys(
                segment_value,
                {"line", "raw_start", "raw_end", "intervals", "score", "phase", "attributes"},
                f"features[{feature_index}].segments[{segment_index}]",
            )
            line = _integer(segment["line"], "feature segment line", minimum=2, maximum=MAX_LINES)
            if line <= last_line or line in source_lines or line in occupied_source_lines:
                raise _fail("BIO058", "feature segment source lines must be unique and ordered")
            last_line = line
            source_lines.add(line)
            total_segments += 1
            if total_segments > MAX_FEATURE_ROWS:
                raise _fail("BIO054", "artifact exceeds feature-row ceiling")
            raw_start = _integer(segment["raw_start"], "feature raw_start", minimum=1)
            raw_end = _integer(segment["raw_end"], "feature raw_end", minimum=1)
            if raw_end < raw_start:
                raise _fail("BIO059", "feature raw_end precedes raw_start")
            if segment["score"] is not None:
                if type(segment["score"]) is not str:
                    raise _fail("BIO010", "feature score must be a decimal string or null")
                _parse_score(segment["score"], line)
            if feature_type == "CDS":
                if type(segment["phase"]) is not int or segment["phase"] not in {0, 1, 2}:
                    raise _fail("BIO035", "CDS phase must be 0, 1, or 2")
            elif segment["phase"] is not None:
                raise _fail("BIO035", "non-CDS phase must be null")
            attributes = _validate_attribute_array(segment["attributes"], "feature attributes")
            ids = _attribute(attributes, "ID")
            if (ids[0] if ids else None) != declared_id:
                raise _fail("BIO055", f"feature {entity_id!r} segment ID does not match entity")
    _sorted_text(feature_ids, "feature entity IDs")

    circular_records = {
        feature["seqid"]
        for feature in features
        for segment in feature["segments"]
        if _attribute(segment["attributes"], "Is_circular") == ["true"]
        and segment["raw_start"] == 1
        and segment["raw_end"] == lengths[feature["seqid"]]
    }
    for feature in features:
        for segment in feature["segments"]:
            wanted = _intervals(
                segment["raw_start"],
                segment["raw_end"],
                lengths[feature["seqid"]],
                feature["seqid"] in circular_records,
                f"line {segment['line']}",
            )
            if segment["intervals"] != wanted:
                raise _fail("BIO060", f"line {segment['line']} normalized intervals are invalid")
            region = region_map.get(feature["seqid"])
            if region is not None and any(
                interval["start"] < region[0] - 1 or interval["end"] > region[1]
                for interval in wanted
            ):
                raise _fail("BIO036", f"line {segment['line']} lies outside its sequence-region")
            if feature["declared_id"] is None:
                row = {
                    "seqid": feature["seqid"],
                    "source": feature["source"],
                    "type": feature["type"],
                    "strand": feature["strand"],
                    "segment": segment,
                }
                if feature["entity_id"] != _anonymous_id(row):
                    raise _fail("BIO061", f"line {segment['line']} anonymous identity is not deterministic")

    relationships = bio_ir["relationships"]
    if type(relationships) is not list or len(relationships) > MAX_RELATIONSHIPS:
        raise _fail("BIO062", "relationships must be a bounded array")
    actual: list[tuple[str, str, str]] = []
    for index, relationship_value in enumerate(relationships):
        relationship = _keys(relationship_value, {"source", "kind", "target"}, f"relationships[{index}]")
        source = _text(relationship["source"], f"relationships[{index}].source")
        target = _text(relationship["target"], f"relationships[{index}].target")
        kind = relationship["kind"]
        if kind not in {"parent", "derives-from"}:
            raise _fail("BIO062", f"relationship {index} has invalid kind")
        if source not in feature_map or target not in feature_map:
            raise _fail("BIO062", f"relationship {index} has a dangling endpoint")
        actual.append((source, kind, target))
    if actual != sorted(actual, key=lambda item: tuple(part.encode("utf-8") for part in item)) or len(actual) != len(set(actual)):
        raise _fail("BIO062", "relationships must be sorted and unique")
    expected_relationships: set[tuple[str, str, str]] = set()
    for feature in features:
        for segment in feature["segments"]:
            for tag, kind in (("Parent", "parent"), ("Derives_from", "derives-from")):
                for target_id in _attribute(segment["attributes"], tag):
                    expected_relationships.add((feature["entity_id"], kind, "id:" + target_id))
    if actual != sorted(expected_relationships, key=lambda item: tuple(part.encode("utf-8") for part in item)):
        raise _fail("BIO063", "relationships do not exactly match GFF3 attributes")
    _cycle_check(relationships)
    if len(canonical_bytes(root)) > MAX_OUTPUT_BYTES:
        raise _fail("BIO040", f"artifact exceeds {MAX_OUTPUT_BYTES} canonical bytes")
    _preflight_pretty_wire(root)
    return root


def load_gff3_artifact(
    path: str | Path,
    sequence_collection: SequenceCollectionArtifact | dict[str, Any],
    *,
    gff3_source: bytes | None = None,
) -> GFF3Artifact:
    try:
        raw = _read_regular(path, MAX_OUTPUT_BYTES, "GFF3 artifact")
        payload = json.loads(raw.decode("utf-8"), object_pairs_hook=_reject_duplicates, parse_constant=_reject_constant)
    except GFF3Error:
        raise
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, ValueError, RecursionError) as failure:
        raise _fail("BIO064", f"cannot load GFF3 artifact: {failure}") from failure
    if type(payload) is not dict:
        raise _fail("BIO064", "GFF3 artifact must be a JSON object")
    validate_gff3_artifact(payload, sequence_collection, gff3_source=gff3_source)
    return GFF3Artifact(payload)


def _reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise _fail("BIO064", f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def _reject_constant(value: str) -> None:
    raise _fail("BIO064", f"non-finite JSON number {value!r}")


__all__ = [
    "FORMAT",
    "VERSION",
    "GFF3Artifact",
    "GFF3Compiler",
    "GFF3Error",
    "load_gff3_artifact",
    "validate_gff3_artifact",
]

"""Independent raw-source validator for the GenBank frontend.

The validator is deliberately self-contained.  It imports only Python's
standard library, parses the original GenBank bytes with its own grammar,
reconstructs the complete GenBank IR and embedded Sequence Collection, and
compares their RFC 8785 canonical forms.  No compiler implementation is
imported or executed during validation.
"""

from __future__ import annotations

import argparse
import base64
from collections.abc import Iterator
from dataclasses import dataclass
from datetime import date as calendar_date
import errno
import hashlib
import json
import math
import os
from pathlib import Path
import re
import secrets
import stat
import sys
import tempfile
from typing import Any


SAFE_INTEGER = 2**53 - 1

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

COLLECTION_FORMAT = "brain01.sequence-collection-ir"
COLLECTION_VERSION = 2
COLLECTION_COMPILER = {
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

REPORT_FORMAT = "brainc.bio.insdc-genbank-validation-report"
REPORT_VERSION = 2
VALIDATOR = {
    "name": "brainc-independent-insdc-genbank-validator",
    "version": "0.2.0",
    "strategy": "independent-genbank-and-sequence-collection-v2-replay",
}
CHECKS = [
    "bounded-original-genbank-bytes",
    "independent-genbank-record-grammar",
    "independent-insdc-location-replay",
    "exact-source-storage-and-origin-line-binding",
    "independent-sequence-collection-v2-replay",
    "sequence-sha-refget-and-seqcol-identities",
    "feature-identity-and-coordinate-normalization",
    "position-only-source-map-replay",
    "canonical-bio-ir-and-artifact-digests",
]

MAX_GENBANK_BYTES = 16 * 1024 * 1024
MAX_ARTIFACT_BYTES = 64 * 1024 * 1024
MAX_REPORT_BYTES = 16 * 1024 * 1024
MAX_LINE_BYTES = 1024 * 1024
MAX_LINES = 1_000_000
MAX_RECORDS = 40_000
MAX_TOTAL_BASES = 16 * 1024 * 1024
MAX_FEATURES = 20_000
MAX_QUALIFIERS = 100_000
MAX_ACCESSIONS = 10_000
MAX_RECORD_REFERENCE_BYTES = 64
MAX_FEATURE_LINE_BYTES = 80
MAX_TEXT_BYTES = 1024 * 1024
MAX_JSON_DEPTH = 64
MAX_JSON_MEMBERS = 1_000_000
MAX_LOCATION_BYTES = 1024 * 1024
MAX_LOCATION_DEPTH = 32
MAX_LOCATION_NODES = 300_005
MAX_LOCATION_CHILDREN = 100_000
MAX_AUTHORITY_MANIFEST_BYTES = 32 * 1024
_AUTHORITY_MANIFEST_PATH = (
    Path(__file__).with_name("standards")
    / "genbank-273-insdc-ft-11.4.authority.json"
)

_SHA_RE = re.compile(r"[0-9a-f]{64}\Z")
_ACCESSION_RE = re.compile(r"[A-Z][A-Z0-9_]*[0-9]\Z")
_VERSION_RE = re.compile(r"([A-Z][A-Z0-9_]*[0-9])\.([1-9][0-9]*)\Z")
_REMOTE_RE = re.compile(r"([A-Z][A-Z0-9_]*[0-9])\.([1-9][0-9]*):")
_MOLECULE_RE = re.compile(r"(?:(ss|ds|ms)-)?DNA\Z")
_DIVISION_RE = re.compile(r"[A-Z]{3}\Z")
_DATE_RE = re.compile(r"[0-9]{2}-[A-Z]{3}-[0-9]{4}\Z")
_COMPONENT_RE = re.compile(r"[A-Za-z0-9_*'\-]+\Z")
_MONTHS = {
    "JAN": 1, "FEB": 2, "MAR": 3, "APR": 4, "MAY": 5, "JUN": 6,
    "JUL": 7, "AUG": 8, "SEP": 9, "OCT": 10, "NOV": 11, "DEC": 12,
}
_IUPAC_LOWER = frozenset("acgtryswkmbdhvn")
_PRE_FIELDS = frozenset({"DBLINK", "PROJECT"})
_POST_FIELDS = {"REFERENCE": 0, "COMMENT": 1, "PRIMARY": 2}


class INSDCValidationError(ValueError):
    """The evidence, artifact, report, or path boundary is invalid."""


def _fail(code: str, detail: str) -> INSDCValidationError:
    return INSDCValidationError(f"{code}: {detail}")


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
    rendered = digits[0] + (("." + digits[1:]) if len(digits) > 1 else "")
    return sign + rendered + ("e+" if normalized_exponent >= 0 else "e") + str(normalized_exponent)


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
        for index, key in enumerate(sorted(value, key=lambda item: item.encode("utf-16be"))):
            if index:
                yield ","
            yield _jcs_string(key)
            yield ":"
            yield from _jcs_chunks(value[key])
        yield "}"
    else:
        raise TypeError(f"unsupported value type {type(value).__name__}")


def _canonical_bytes(value: Any) -> bytes:
    try:
        return "".join(_jcs_chunks(value)).encode("utf-8")
    except (TypeError, ValueError, UnicodeError, RecursionError) as failure:
        raise _fail("GIV001", f"value is not RFC 8785 canonical JSON: {failure}") from failure


def _digest(value: Any) -> str:
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
        raise _fail("GIV001", f"value is not RFC 8785 canonical JSON: {failure}") from failure
    return hasher.hexdigest()


def _canonical_size(value: Any, maximum: int) -> int:
    total = 0
    try:
        for chunk in _jcs_chunks(value):
            total += len(chunk.encode("utf-8"))
            if total > maximum:
                return total
    except (TypeError, ValueError, UnicodeError, RecursionError) as failure:
        raise _fail("GIV001", f"value is not RFC 8785 canonical JSON: {failure}") from failure
    return total


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


def _sha512t24u(value: bytes) -> str:
    return base64.urlsafe_b64encode(hashlib.sha512(value).digest()[:24]).decode("ascii")


def _refget(sequence: str) -> str:
    raw = hashlib.sha512(sequence.upper().encode("ascii")).digest()[:24]
    return "SQ." + base64.urlsafe_b64encode(raw).decode("ascii").rstrip("=")


def _validate_tree(value: Any, label: str) -> None:
    members = 0
    stack: list[tuple[Any, int]] = [(value, 0)]
    while stack:
        current, depth = stack.pop()
        if depth > MAX_JSON_DEPTH:
            raise _fail("GIV002", f"{label} exceeds JSON depth {MAX_JSON_DEPTH}")
        if current is None or type(current) is bool:
            continue
        if type(current) is int:
            if not -SAFE_INTEGER <= current <= SAFE_INTEGER:
                raise _fail("GIV002", f"{label} contains an unsafe integer")
            continue
        if type(current) is float:
            raise _fail("GIV002", f"{label} contains an unexpected JSON number")
        if type(current) is str:
            try:
                encoded = current.encode("utf-8")
            except UnicodeEncodeError as failure:
                raise _fail("GIV002", f"{label} contains invalid Unicode") from failure
            if len(encoded) > MAX_TEXT_BYTES:
                raise _fail("GIV002", f"{label} contains an oversized string")
            continue
        if type(current) is list:
            members += len(current)
            stack.extend((item, depth + 1) for item in reversed(current))
        elif type(current) is dict:
            members += len(current)
            for key, item in reversed(list(current.items())):
                if type(key) is not str:
                    raise _fail("GIV002", f"{label} contains a non-string key")
                stack.append((key, depth + 1))
                stack.append((item, depth + 1))
        else:
            raise _fail("GIV002", f"{label} contains unsupported {type(current).__name__}")
        if members > MAX_JSON_MEMBERS:
            raise _fail("GIV002", f"{label} exceeds {MAX_JSON_MEMBERS} JSON members")


def _exact_keys(value: Any, expected: set[str], label: str) -> dict[str, Any]:
    if type(value) is not dict:
        raise _fail("GIV003", f"{label} must be an object")
    missing = sorted(expected - value.keys())
    extra = sorted(value.keys() - expected)
    if missing or extra:
        raise _fail("GIV003", f"{label} keys invalid; missing={missing or 'none'}, unknown={extra or 'none'}")
    return value


def _assert_sha(value: Any, label: str) -> str:
    if type(value) is not str or _SHA_RE.fullmatch(value) is None:
        raise _fail("GIV004", f"{label} must be lowercase SHA-256")
    return value


def _integer(
    value: Any,
    label: str,
    minimum: int = 0,
    maximum: int = SAFE_INTEGER,
) -> int:
    if type(value) is not int or not minimum <= value <= maximum:
        raise _fail("GIV004", f"{label} must be an integer in [{minimum}, {maximum}]")
    return value


def _positive_decimal(token: str, label: str, line: int) -> int:
    if not token or not token.isascii() or not token.isdecimal() or (len(token) > 1 and token[0] == "0"):
        raise _fail("GIV005", f"line {line} has an invalid {label}")
    value = 0
    for character in token:
        digit = ord(character) - 48
        if value > (SAFE_INTEGER - digit) // 10:
            raise _fail("GIV005", f"line {line} {label} exceeds the I-JSON safe range")
        value = value * 10 + digit
    if value < 1:
        raise _fail("GIV005", f"line {line} {label} must be positive")
    return value


@dataclass(frozen=True, slots=True)
class _Line:
    number: int
    text: str
    eol: str


def _chunk(line: _Line, column: int, text: str) -> dict[str, Any]:
    if not 1 <= column <= len(line.text) + 1:
        raise _fail("GIV006", f"line {line.number} has an invalid source column")
    return {"line": line.number, "column_start": column, "column_end_exclusive": column + len(text), "text": text}


def _position(chunk: dict[str, Any]) -> dict[str, int]:
    return {
        "line": chunk["line"],
        "column_start": chunk["column_start"],
        "column_end_exclusive": chunk["column_end_exclusive"],
    }


def _split_lines(raw: bytes) -> tuple[_Line, ...]:
    if type(raw) is not bytes or not raw:
        raise _fail("GIV007", "GenBank source must be nonempty bytes")
    if len(raw) > MAX_GENBANK_BYTES:
        raise _fail("GIV007", f"GenBank source exceeds {MAX_GENBANK_BYTES} bytes")
    physical_lines = raw.count(b"\n") + (0 if raw.endswith(b"\n") else 1)
    if physical_lines > MAX_LINES:
        raise _fail("GIV007", f"source exceeds {MAX_LINES} lines")
    lines: list[_Line] = []
    start = 0
    number = 1
    while start < len(raw):
        newline = raw.find(b"\n", start)
        if newline < 0:
            body, eol, start = raw[start:], "none", len(raw)
        elif newline > start and raw[newline - 1] == 13:
            body, eol, start = raw[start:newline - 1], "CRLF", newline + 1
        else:
            body, eol, start = raw[start:newline], "LF", newline + 1
        if len(body) > MAX_LINE_BYTES or b"\r" in body:
            raise _fail("GIV007", f"line {number} exceeds limits or contains bare CR")
        if any(octet < 0x20 or octet > 0x7E for octet in body):
            raise _fail("GIV007", f"line {number} must be printable ASCII")
        lines.append(_Line(number, body.decode("ascii"), eol))
        if len(lines) > MAX_LINES:
            raise _fail("GIV007", f"source exceeds {MAX_LINES} lines")
        number += 1
    return tuple(lines)


def _split_records(lines: tuple[_Line, ...]) -> tuple[tuple[_Line, ...], ...]:
    records: list[tuple[_Line, ...]] = []
    start = 0
    while start < len(lines):
        if not lines[start].text.startswith("LOCUS       "):
            raise _fail("GIV008", f"record at line {lines[start].number} must begin LOCUS")
        end = start
        while end < len(lines) and lines[end].text != "//":
            end += 1
        if end == len(lines):
            raise _fail("GIV008", f"record at line {lines[start].number} has no exact terminator")
        records.append(lines[start:end + 1])
        if len(records) > MAX_RECORDS:
            raise _fail("GIV008", f"source exceeds {MAX_RECORDS} records")
        start = end + 1
        if start < len(lines) and lines[start].text == "":
            start += 1
    return tuple(records)


def _reference(value: str) -> str:
    if type(value) is not str or len(value) > MAX_RECORD_REFERENCE_BYTES:
        raise _fail(
            "GIV009",
            f"remote reference must be at most {MAX_RECORD_REFERENCE_BYTES} ASCII bytes",
        )
    match = _VERSION_RE.fullmatch(value)
    if match is None:
        raise _fail("GIV009", "remote reference must be an exact uppercase accession.version")
    _positive_decimal(match.group(2), "accession version", 0)
    return value


def _bound(node: dict[str, Any], context: dict[str, Any] | None) -> None:
    kind = node["kind"]
    if kind == "point":
        coordinates = (node["position"],)
    elif kind == "span":
        coordinates = (node["start"]["position"], node["end"]["position"])
    else:
        coordinates = (node["left"], node["right"])
    if context is not None and any(value > context["length"] for value in coordinates):
        raise _fail("GIV009", "location coordinate exceeds explicit sequence length")
    if kind == "between":
        left, right = coordinates
        adjacent = left < SAFE_INTEGER and right == left + 1
        circular = (
            right == 1
            and context is not None
            and context["topology"] == "circular"
            and left == context["length"]
        )
        if not adjacent and not circular:
            raise _fail("GIV009", "between-base coordinates are not adjacent or circular length^1")
    if kind == "within" and coordinates[0] >= coordinates[1]:
        raise _fail("GIV009", "within-range coordinates must increase")


class _LocationParser:
    """Second implementation of the closed structural location grammar."""

    def __init__(self, source: str, local: dict[str, Any], references: dict[str, dict[str, Any]]) -> None:
        if type(source) is not str or not source or len(source) > MAX_LOCATION_BYTES:
            raise _fail("GIV009", "location must be bounded nonempty text")
        try:
            source.encode("ascii")
        except UnicodeEncodeError as failure:
            raise _fail("GIV009", "location must be ASCII") from failure
        if any(character.isspace() and character != " " for character in source):
            raise _fail("GIV009", "location may contain ASCII spaces only")
        if re.search(r"[A-Za-z0-9_*'\-] +[A-Za-z0-9_*'\-]", source):
            raise _fail("GIV009", "spaces split a location component")
        if re.search(r"(?:complement|join|order) +\(", source):
            raise _fail("GIV009", "spaces split a location operator")
        self.source = source.replace(" ", "")
        self.local = local
        self.references = references
        self.position = 0
        self.nodes = 0
        self.children = 0

    def parse(self) -> dict[str, Any]:
        if not self.source:
            raise _fail("GIV009", "location contains no descriptor")
        result = self._location(1, False)
        if self.position != len(self.source):
            raise _fail("GIV009", f"unexpected location token at byte {self.position + 1}")
        return result

    def _reserve(self, depth: int) -> None:
        if depth > MAX_LOCATION_DEPTH:
            raise _fail("GIV009", f"location exceeds depth {MAX_LOCATION_DEPTH}")
        self.nodes += 1
        if self.nodes > MAX_LOCATION_NODES:
            raise _fail("GIV009", f"location exceeds {MAX_LOCATION_NODES} nodes")

    def _location(self, depth: int, compound: bool) -> dict[str, Any]:
        if self.source.startswith("complement(", self.position):
            self._reserve(depth)
            self.position += 11
            child = self._location(depth + 1, compound)
            self._literal(")")
            if child["kind"] == "complement" or (compound and child["kind"] in {"join", "order"}):
                raise _fail("GIV009", "invalid nested complement or compound")
            return {"kind": "complement", "location": child}
        if self.source.startswith("join(", self.position):
            return self._compound(depth, "join", compound)
        if self.source.startswith("order(", self.position):
            return self._compound(depth, "order", compound)
        if self.position < len(self.source) and "A" <= self.source[self.position] <= "Z":
            return self._remote(depth)
        return self._atomic(depth, self.local)

    def _compound(self, depth: int, operator: str, inside: bool) -> dict[str, Any]:
        if inside:
            raise _fail("GIV009", "compound locations cannot be nested or mixed")
        self._reserve(depth)
        self.position += len(operator) + 1
        children: list[dict[str, Any]] = []
        while True:
            self.children += 1
            if self.children > MAX_LOCATION_CHILDREN:
                raise _fail("GIV009", f"location exceeds {MAX_LOCATION_CHILDREN} children")
            children.append(self._location(depth + 1, True))
            if len(children) > MAX_LOCATION_CHILDREN:
                raise _fail("GIV009", f"{operator} has too many children")
            if self._take(","):
                continue
            self._literal(")")
            break
        if len(children) < 2:
            raise _fail("GIV009", f"{operator} requires two children")
        return {"kind": operator, "children": children}

    def _remote(self, depth: int) -> dict[str, Any]:
        self._reserve(depth)
        match = _REMOTE_RE.match(self.source, self.position)
        if match is None:
            raise _fail("GIV009", "remote location requires accession.version")
        reference = _reference(match.group(1) + "." + match.group(2))
        self.position = match.end()
        child = self._atomic(depth + 1, self.references.get(reference))
        return {"kind": "remote", "reference": reference, "location": child}

    def _atomic(self, depth: int, context: dict[str, Any] | None) -> dict[str, Any]:
        self._reserve(depth)
        left_fuzz = self._fuzz()
        left = self._coordinate("coordinate")
        if self._take(".."):
            right_fuzz = self._fuzz()
            right = self._coordinate("coordinate")
            if left_fuzz not in {None, "<"} or right_fuzz not in {None, ">"} or left > right:
                raise _fail("GIV009", "invalid fuzzy span")
            node = {
                "kind": "span",
                "start": {"position": left, "fuzz": left_fuzz},
                "end": {"position": right, "fuzz": right_fuzz},
            }
        elif self._take("^"):
            if left_fuzz is not None or self._fuzz() is not None:
                raise _fail("GIV009", "between-base coordinates cannot be fuzzy")
            node = {"kind": "between", "left": left, "right": self._coordinate("between coordinate")}
        elif self._take("."):
            if left_fuzz is not None:
                raise _fail("GIV009", "within-range coordinates cannot be fuzzy")
            node = {"kind": "within", "left": left, "right": self._coordinate("within coordinate")}
        else:
            if left_fuzz is not None:
                raise _fail("GIV009", "fuzz is accepted only on span endpoints")
            node = {"kind": "point", "position": left}
        _bound(node, context)
        return node

    def _coordinate(self, label: str) -> int:
        start = self.position
        if start >= len(self.source) or not self.source[start].isdigit() or self.source[start] == "0":
            raise _fail("GIV009", f"invalid {label} at byte {start + 1}")
        value = 0
        while self.position < len(self.source) and self.source[self.position].isdigit():
            digit = ord(self.source[self.position]) - 48
            if value > (SAFE_INTEGER - digit) // 10:
                raise _fail("GIV009", f"{label} exceeds safe integer range")
            value = value * 10 + digit
            self.position += 1
        return value

    def _fuzz(self) -> str | None:
        if self._take("<"):
            return "<"
        if self._take(">"):
            return ">"
        return None

    def _take(self, token: str) -> bool:
        if self.source.startswith(token, self.position):
            self.position += len(token)
            return True
        return False

    def _literal(self, token: str) -> None:
        if not self._take(token):
            raise _fail("GIV009", f"expected {token!r} at byte {self.position + 1}")


def _segments(
    node: dict[str, Any],
    local: dict[str, Any],
    references: dict[str, dict[str, Any]],
    active_reference: str | None = None,
    active_context: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    if active_context is None and active_reference is None:
        active_context = local
    kind = node["kind"]
    if kind == "remote":
        return _segments(node["location"], local, references, node["reference"], references.get(node["reference"]))
    if kind in {"point", "span", "between", "within"}:
        _bound(node, active_context)
        status = "verified" if active_context is not None else "unresolved"
        if kind == "point":
            start, end, segment_kind, start_fuzz, end_fuzz = node["position"] - 1, node["position"], "interval", None, None
        elif kind == "span":
            start, end, segment_kind = node["start"]["position"] - 1, node["end"]["position"], "interval"
            start_fuzz, end_fuzz = node["start"]["fuzz"], node["end"]["fuzz"]
        elif kind == "between":
            start = 0 if node["right"] == 1 and (active_context is None or node["left"] == active_context["length"]) else node["left"]
            end, segment_kind, start_fuzz, end_fuzz = start, "between", None, None
        else:
            start, end, segment_kind, start_fuzz, end_fuzz = node["left"] - 1, node["right"], "uncertain-point", None, None
        return [{
            "reference": active_reference,
            "start": start,
            "end": end,
            "orientation": 1,
            "kind": segment_kind,
            "start_fuzz": start_fuzz,
            "end_fuzz": end_fuzz,
            "bounds_status": status,
        }]
    if kind == "complement":
        children = _segments(node["location"], local, references, active_reference, active_context)
        return [{**segment, "orientation": -segment["orientation"]} for segment in reversed(children)]
    if kind in {"join", "order"}:
        result: list[dict[str, Any]] = []
        for child in node["children"]:
            result.extend(_segments(child, local, references, active_reference, active_context))
        return result
    raise _fail("GIV009", "unsupported location AST node")


def _field_prefix(name: str) -> str:
    return f"{name:<12}"


def _header(name: str, value: str, chunks: list[dict[str, Any]]) -> dict[str, Any]:
    return {"name": name, "value": value, "chunks": chunks}


def _parse_field(lines: tuple[_Line, ...], index: int, name: str) -> tuple[dict[str, Any], int]:
    if index >= len(lines) or not lines[index].text.startswith(_field_prefix(name)):
        raise _fail("GIV010", f"expected {name} at source line boundary")
    first = lines[index]
    parts = [first.text[12:]]
    chunks = [_chunk(first, 13, first.text[12:])]
    index += 1
    while index < len(lines) and lines[index].text.startswith(" " * 12):
        continuation = lines[index]
        parts.append(continuation.text[12:])
        chunks.append(_chunk(continuation, 13, continuation.text[12:]))
        index += 1
    stripped = [part.strip() for part in parts]
    if any(not part for part in stripped):
        raise _fail("GIV010", f"{name} contains an empty chunk")
    return _header(name, " ".join(stripped), chunks), index


def _optional_block(lines: tuple[_Line, ...], index: int, name: str) -> tuple[dict[str, Any], int]:
    first = lines[index]
    if not first.text.startswith(_field_prefix(name)):
        raise _fail("GIV010", f"line {first.number} does not begin {name}")
    block = [first]
    index += 1
    while index < len(lines) and lines[index].text.startswith(" "):
        block.append(lines[index])
        index += 1
    values: list[str] = []
    chunks: list[dict[str, Any]] = []
    for position, line in enumerate(block):
        if position == 0:
            text, column = line.text[12:], 13
        else:
            leading = len(line.text) - len(line.text.lstrip(" "))
            text, column = line.text[leading:], leading + 1
        if not text.strip() and name != "COMMENT":
            raise _fail("GIV010", f"{name} contains an empty chunk")
        values.append(text.strip())
        chunks.append(_chunk(line, column, text))
    return _header(name, " ".join(value for value in values if value), chunks), index


def _parse_locus(line: _Line) -> dict[str, Any]:
    if not line.text.startswith("LOCUS       "):
        raise _fail("GIV011", f"line {line.number} is not a traditional LOCUS")
    tokens = line.text[12:].split()
    if len(tokens) not in {6, 7}:
        raise _fail("GIV011", f"line {line.number} LOCUS must have six or seven tokens")
    if len(tokens) == 6:
        name, length_token, unit, molecule_token, division, date = tokens
        topology, declared = "linear", False
    else:
        name, length_token, unit, molecule_token, topology, division, date = tokens
        declared = True
    length = _positive_decimal(length_token, "LOCUS length", line.number)
    if length > MAX_TOTAL_BASES:
        raise _fail("GIV011", f"line {line.number} record length exceeds {MAX_TOTAL_BASES} bases")
    molecule = _MOLECULE_RE.fullmatch(molecule_token)
    if (
        not name or len(name) > 64 or any(character.isspace() for character in name)
        or unit != "bp" or molecule is None or topology not in {"linear", "circular"}
        or _DIVISION_RE.fullmatch(division) is None or _DATE_RE.fullmatch(date) is None
    ):
        raise _fail("GIV011", f"line {line.number} contains invalid LOCUS fields")
    try:
        day, month, year = date.split("-")
        calendar_date(int(year), _MONTHS[month], int(day))
    except (KeyError, ValueError) as failure:
        raise _fail("GIV011", f"line {line.number} contains an invalid calendar date") from failure
    return {
        "name": name,
        "length": length,
        "unit": unit,
        "molecule": "DNA",
        "strandedness": molecule.group(1),
        "topology": topology,
        "topology_declared": declared,
        "division": division,
        "date": date,
        "source_map": _chunk(line, 1, line.text),
    }


def _component(value: str, maximum: int, label: str, line: int) -> str:
    if (
        not value or len(value) > maximum or _COMPONENT_RE.fullmatch(value) is None
        or not any(character.isalpha() and character.isascii() for character in value)
    ):
        raise _fail("GIV012", f"line {line} has an invalid {label}")
    return value


def _qualifier_value(parts: list[str], name: str, line: int) -> tuple[bool, str | None]:
    source = ("" if name == "translation" else " ").join(part.strip() for part in parts)
    if not source.startswith('"'):
        if not source or '"' in source or any(not 0x20 <= ord(character) <= 0x7E for character in source):
            raise _fail("GIV012", f"qualifier /{name} at line {line} has invalid unquoted text")
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
                raise _fail("GIV012", f"qualifier /{name} has trailing text")
            return True, "".join(value)
        value.append(character)
        index += 1
    return False, None


def _qualifier_parts(parts: list[str], line: int) -> tuple[bool, str, str | None]:
    if not parts[0].startswith("/"):
        raise _fail("GIV012", f"line {line} qualifier lacks slash")
    body = parts[0][1:]
    if "=" not in body:
        name = _component(body, 20, "qualifier name", line)
        if len(parts) != 1:
            raise _fail("GIV012", f"valueless qualifier /{name} cannot continue")
        return True, name, None
    name, initial = body.split("=", 1)
    name = _component(name, 20, "qualifier name", line)
    complete, value = _qualifier_value([initial, *parts[1:]], name, line)
    return complete, name, value


def _quoted_continuation_complete(part: str, name: str, line: int) -> bool:
    source = part.strip()
    index = 0
    while index < len(source):
        if source[index] == '"':
            if index + 1 < len(source) and source[index + 1] == '"':
                index += 2
                continue
            if index != len(source) - 1:
                raise _fail("GIV012", f"qualifier /{name} has trailing text")
            return True
        index += 1
    return False


def _parse_features(lines: tuple[_Line, ...], index: int) -> tuple[list[dict[str, Any]], int]:
    features: list[dict[str, Any]] = []
    qualifiers_seen = 0
    while index < len(lines) and not (lines[index].text == "ORIGIN" or lines[index].text.startswith(_field_prefix("ORIGIN"))):
        line = lines[index]
        if line.text.startswith("CONTIG"):
            raise _fail("GIV013", f"line {line.number} is a CONTIG record")
        if len(line.text) > MAX_FEATURE_LINE_BYTES or len(line.text) < 22 or line.text[:5] != " " * 5 or line.text[20] != " ":
            raise _fail("GIV013", f"line {line.number} is not a fixed-column feature")
        key_field = line.text[5:20]
        key = key_field.rstrip(" ")
        if not key or key_field != key.ljust(15):
            raise _fail("GIV013", f"line {line.number} has a misaligned feature key")
        key = _component(key, 15, "feature key", line.number)
        location_parts = [line.text[21:].strip()]
        if not location_parts[0]:
            raise _fail("GIV013", f"line {line.number} has no location")
        location_chunks = [_chunk(line, 22, line.text[21:])]
        line_start = line.number
        line_end = line.number
        index += 1
        builders: list[dict[str, Any]] = []
        current: dict[str, Any] | None = None
        while index < len(lines):
            continuation = lines[index]
            if continuation.text == "ORIGIN" or continuation.text.startswith(_field_prefix("ORIGIN")):
                break
            if len(continuation.text) >= 22 and continuation.text[:21] == " " * 21:
                if len(continuation.text) > MAX_FEATURE_LINE_BYTES:
                    raise _fail("GIV013", f"line {continuation.number} exceeds feature width")
                body = continuation.text[21:]
                if not body:
                    raise _fail("GIV013", f"line {continuation.number} has empty continuation")
                source_chunk = _chunk(continuation, 22, body)
                if current is None and not builders and not body.startswith("/"):
                    location_parts.append(body.strip())
                    location_chunks.append(source_chunk)
                elif body.startswith("/") and (current is None or current["complete"]):
                    complete, name, _ = _qualifier_parts([body], continuation.number)
                    remainder = body.split("=", 1)[1] if "=" in body[1:] else None
                    mode = "valueless" if remainder is None else "quoted" if remainder.startswith('"') else "unquoted"
                    current = {"parts": [body], "chunks": [source_chunk], "line": continuation.number, "name": name, "complete": complete, "mode": mode}
                    builders.append(current)
                elif current is not None and not current["complete"]:
                    current["parts"].append(body)
                    current["chunks"].append(source_chunk)
                    current["complete"] = _quoted_continuation_complete(body, current["name"], current["line"])
                elif current is not None and current["mode"] == "unquoted":
                    current["parts"].append(body)
                    current["chunks"].append(source_chunk)
                else:
                    raise _fail("GIV013", f"line {continuation.number} has an unexpected continuation")
                line_end = continuation.number
                index += 1
                continue
            if continuation.text.startswith(" " * 5) and len(continuation.text) >= 21:
                break
            raise _fail("GIV013", f"line {continuation.number} is outside feature grammar")
        qualifiers: list[dict[str, Any]] = []
        for builder in builders:
            complete, name, value = _qualifier_parts(builder["parts"], builder["line"])
            if not complete:
                raise _fail("GIV013", f"qualifier /{name} has no closing quote")
            qualifiers.append({"name": name, "value": value, "chunks": builder["chunks"]})
        qualifiers_seen += len(qualifiers)
        if qualifiers_seen > MAX_QUALIFIERS:
            raise _fail("GIV013", "too many qualifiers")
        features.append({
            "ordinal": len(features) + 1,
            "key": key,
            "location_text": "".join(location_parts),
            "location_chunks": location_chunks,
            "qualifiers": qualifiers,
            "line_start": line_start,
            "line_end": line_end,
        })
        if len(features) > MAX_FEATURES:
            raise _fail("GIV013", "too many features")
    return features, index


def _organism_value(header: dict[str, Any]) -> tuple[str, str]:
    parts = [chunk["text"].strip() for chunk in header["chunks"]]
    if len(parts) < 2:
        raise _fail("GIV014", "ORGANISM must contain a taxonomic lineage")
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
        raise _fail("GIV014", "ORGANISM name and lineage must be nonempty")
    return organism, lineage


def _origin_header(line: _Line) -> dict[str, Any]:
    if line.text in {"ORIGIN", _field_prefix("ORIGIN")}:
        return _header("ORIGIN", "", [_chunk(line, 1, line.text)])
    if not line.text.startswith(_field_prefix("ORIGIN")):
        raise _fail("GIV015", f"line {line.number} is not ORIGIN")
    value = line.text[12:]
    if not value or value != value.strip() or not value.endswith("."):
        raise _fail("GIV015", f"line {line.number} has an invalid nonblank ORIGIN value")
    return _header("ORIGIN", value, [_chunk(line, 1, line.text)])


def _parse_origin(
    lines: tuple[_Line, ...],
    index: int,
    expected: int,
) -> tuple[str, int, int, int]:
    first_index = index
    parts: list[str] = []
    lengths: list[int] = []
    normalized = 0
    while index < len(lines) and lines[index].text != "//":
        line = lines[index]
        if len(line.text) < 11 or line.text[9] != " " or line.text[:9] != f"{normalized + 1:>9}":
            raise _fail("GIV015", f"line {line.number} has invalid ORIGIN indexing")
        groups = line.text[10:].split(" ")
        if (
            not 1 <= len(groups) <= 6 or any(not group for group in groups)
            or any(len(group) != 10 for group in groups[:-1])
            or not 1 <= len(groups[-1]) <= 10
        ):
            raise _fail("GIV015", f"line {line.number} has invalid ORIGIN groups")
        line_length = sum(len(group) for group in groups)
        if line_length > 60:
            raise _fail("GIV015", f"line {line.number} has too many bases")
        for group in groups:
            if any(base not in _IUPAC_LOWER for base in group):
                raise _fail("GIV015", f"line {line.number} is not lowercase IUPAC DNA")
            parts.append(group.upper())
            normalized += len(group)
        lengths.append(line_length)
        index += 1
    if not lengths or any(length != 60 for length in lengths[:-1]):
        raise _fail("GIV015", "ORIGIN must contain full non-final sequence lines")
    sequence = "".join(parts)
    if len(sequence) != expected:
        raise _fail("GIV015", f"LOCUS length {expected} differs from ORIGIN length {len(sequence)}")
    return sequence, lines[first_index].number, lines[index - 1].number, index


def _accessions(value: str, line: int) -> list[str]:
    tokens = value.split()
    if (
        not tokens or len(tokens) > MAX_ACCESSIONS or len(tokens) != len(set(tokens))
        or any(
            _ACCESSION_RE.fullmatch(token) is None
            or len(token.encode("ascii")) > MAX_RECORD_REFERENCE_BYTES
            for token in tokens
        )
    ):
        raise _fail("GIV016", f"line {line} has an invalid ACCESSION list")
    return tokens


def _parse_record(lines: tuple[_Line, ...]) -> dict[str, Any]:
    locus = _parse_locus(lines[0])
    index = 1
    headers: list[dict[str, Any]] = []
    definition, index = _parse_field(lines, index, "DEFINITION")
    accession_header, index = _parse_field(lines, index, "ACCESSION")
    version_header, index = _parse_field(lines, index, "VERSION")
    headers.extend((definition, accession_header, version_header))
    accessions = _accessions(accession_header["value"], accession_header["chunks"][0]["line"])
    version_tokens = version_header["value"].split()
    if (
        len(version_tokens) not in {1, 2} or _VERSION_RE.fullmatch(version_tokens[0]) is None
        or (len(version_tokens) == 2 and not version_tokens[1].startswith("GI:"))
    ):
        raise _fail("GIV016", "VERSION must contain accession.version and optional GI")
    version = version_tokens[0]
    if len(version.encode("ascii")) > MAX_RECORD_REFERENCE_BYTES:
        raise _fail(
            "GIV016",
            f"primary VERSION exceeds {MAX_RECORD_REFERENCE_BYTES} ASCII bytes",
        )
    legacy_gi = None if len(version_tokens) == 1 else _positive_decimal(
        version_tokens[1][3:], "legacy VERSION GI", version_header["chunks"][0]["line"]
    )
    match = _VERSION_RE.fullmatch(version)
    assert match is not None
    if match.group(1) != accessions[0]:
        raise _fail("GIV016", "VERSION base accession differs from primary ACCESSION")

    seen_pre: set[str] = set()
    while index < len(lines):
        name = lines[index].text[:12].rstrip()
        if name not in _PRE_FIELDS:
            break
        if name in seen_pre:
            raise _fail("GIV017", f"duplicate {name}")
        optional, index = _optional_block(lines, index, name)
        headers.append(optional)
        seen_pre.add(name)

    keywords, index = _parse_field(lines, index, "KEYWORDS")
    headers.append(keywords)
    if index < len(lines) and lines[index].text.startswith(_field_prefix("SEGMENT")):
        segment, index = _parse_field(lines, index, "SEGMENT")
        segment_tokens = segment["value"].split()
        if len(segment["chunks"]) != 1 or len(segment_tokens) != 3 or segment_tokens[1] != "of":
            raise _fail("GIV017", "SEGMENT must be one line containing 'n of m'")
        segment_number = _positive_decimal(
            segment_tokens[0],
            "SEGMENT number",
            segment["chunks"][0]["line"],
        )
        segment_total = _positive_decimal(
            segment_tokens[2],
            "SEGMENT total",
            segment["chunks"][0]["line"],
        )
        if segment_total < 2 or segment_number > segment_total:
            raise _fail("GIV017", "SEGMENT must satisfy 1 <= n <= m and m >= 2")
        headers.append(segment)
    source, index = _parse_field(lines, index, "SOURCE")
    organism, index = _parse_field(lines, index, "  ORGANISM")
    headers.extend((source, organism))
    organism_name, lineage = _organism_value(organism)

    last_rank = -1
    singletons: set[str] = set()
    references = 0
    while index < len(lines) and lines[index].text != "FEATURES             Location/Qualifiers":
        line = lines[index]
        if not line.text or line.text.startswith("CONTIG") or line.text.startswith(" "):
            raise _fail("GIV017", f"line {line.number} is outside accepted top-level grammar")
        name = line.text[:12].rstrip()
        if name not in _POST_FIELDS:
            raise _fail("GIV017", f"unknown top-level field {name!r}")
        rank = _POST_FIELDS[name]
        if rank < last_rank or (name != "REFERENCE" and name in singletons):
            raise _fail("GIV017", f"misordered or duplicate {name}")
        optional, index = _optional_block(lines, index, name)
        if name == "REFERENCE":
            references += 1
            token = optional["chunks"][0]["text"].split(maxsplit=1)[0]
            if _positive_decimal(token, "REFERENCE number", optional["chunks"][0]["line"]) != references:
                raise _fail("GIV017", "REFERENCE numbers must be consecutive")
            journals = [
                chunk for chunk in optional["chunks"][1:]
                if chunk["column_start"] == 3 and chunk["text"].startswith("JOURNAL   ")
            ]
            if len(journals) != 1:
                raise _fail("GIV017", "REFERENCE must contain exactly one JOURNAL")
        headers.append(optional)
        singletons.add(name)
        last_rank = rank
    if index >= len(lines) or lines[index].text != "FEATURES             Location/Qualifiers" or references == 0:
        raise _fail("GIV017", "record requires REFERENCES and exact FEATURES header")
    headers.append(_header("FEATURES", "Location/Qualifiers", [_chunk(lines[index], 1, lines[index].text)]))
    index += 1
    features, index = _parse_features(lines, index)
    if index >= len(lines):
        raise _fail("GIV015", "record has no ORIGIN")
    headers.append(_origin_header(lines[index]))
    sequence, origin_line_start, origin_line_end, index = _parse_origin(
        lines,
        index + 1,
        locus["length"],
    )
    if index != len(lines) - 1 or lines[index].text != "//":
        raise _fail("GIV008", "record has invalid terminator placement")
    return {
        "line_start": lines[0].number,
        "line_end": lines[-1].number,
        "origin_line_start": origin_line_start,
        "origin_line_end": origin_line_end,
        "locus": locus,
        "definition": definition["value"],
        "accessions": accessions,
        "version": version,
        "legacy_gi": legacy_gi,
        "keywords": keywords["value"],
        "source": source["value"],
        "organism": organism_name,
        "lineage": lineage,
        "headers": headers,
        "features": features,
        "sequence": sequence,
    }


def _storage(raw: bytes, label: str) -> dict[str, Any]:
    if type(raw) is not bytes or not raw:
        raise _fail("GIV030", f"{label} must be nonempty bytes")
    chunks: list[str] = []
    for start in range(0, len(raw), CHUNK_BYTES):
        try:
            chunks.append(raw[start:start + CHUNK_BYTES].decode("ascii"))
        except UnicodeDecodeError as failure:
            raise _fail("GIV030", f"{label} must be ASCII") from failure
    return {
        "kind": CHUNK_KIND,
        "version": CHUNK_VERSION,
        "chunk_bytes": CHUNK_BYTES,
        "chunks": chunks,
    }


def _refget_seqcol(members: list[dict[str, Any]]) -> dict[str, Any]:
    level_2 = {
        "lengths": [member["sequence"]["bases"] for member in members],
        "names": [member["record_id"] for member in members],
        "sequences": [member["sequence"]["refget_id"] for member in members],
    }
    level_1 = {
        name: _sha512t24u(_canonical_bytes(value))
        for name, value in level_2.items()
    }
    inherent = {name: level_1[name] for name in ("names", "sequences")}
    return {
        "version": "1.0.0",
        "digest": _sha512t24u(_canonical_bytes(inherent)),
        "level_1": level_1,
        "level_2": level_2,
    }


def _sequence_collection(raw: bytes, records: list[dict[str, Any]]) -> dict[str, Any]:
    if len(records) > MAX_RECORDS:
        raise _fail("GIV030", f"collection exceeds {MAX_RECORDS} records")
    members: list[dict[str, Any]] = []
    total_bases = 0
    previous_line_end = 0
    for record in records:
        sequence = record["sequence"].upper().encode("ascii")
        total_bases += len(sequence)
        if total_bases > MAX_TOTAL_BASES:
            raise _fail("GIV030", f"collection exceeds {MAX_TOTAL_BASES} normalized bases")
        if record["origin_line_start"] <= previous_line_end:
            raise _fail("GIV030", "ORIGIN source maps must be ordered and disjoint")
        previous_line_end = record["origin_line_end"]
        members.append({
            "record_id": record["version"],
            "source_map": {
                "kind": "genbank-origin-lines",
                "version": 1,
                "line_start": record["origin_line_start"],
                "line_end": record["origin_line_end"],
            },
            "sequence": {
                "alphabet": "IUPAC-DNA",
                "normalization": "uppercase",
                "bases": len(sequence),
                "sha256": hashlib.sha256(sequence).hexdigest(),
                "refget_id": _refget(sequence.decode("ascii")),
                "storage": _storage(sequence, f"record {record['version']!r} sequence"),
            },
        })
    core = {
        "format": COLLECTION_FORMAT,
        "version": COLLECTION_VERSION,
        "compiler": {**COLLECTION_COMPILER, "passes": list(COLLECTION_COMPILER["passes"])},
        "inputs": {
            "kind": "genbank-flatfile",
            "profile": PROFILE,
            "source": {
                "encoding": "ascii",
                "byte_length": len(raw),
                "sha256": hashlib.sha256(raw).hexdigest(),
                "storage": _storage(raw, "GenBank source"),
            },
        },
        "members": members,
        "refget_seqcol": _refget_seqcol(members),
    }
    return {**core, "artifact_sha256": _digest(core)}


def _storage_bytes(value: Any, label: str, *, maximum: int) -> bytes:
    root = _exact_keys(value, {"kind", "version", "chunk_bytes", "chunks"}, label)
    if (
        root["kind"] != CHUNK_KIND
        or type(root["version"]) is not int
        or root["version"] != CHUNK_VERSION
        or type(root["chunk_bytes"]) is not int
        or root["chunk_bytes"] != CHUNK_BYTES
    ):
        raise _fail("GIV031", f"{label} has an unsupported chunk contract")
    chunks = root["chunks"]
    maximum_chunks = (maximum + CHUNK_BYTES - 1) // CHUNK_BYTES
    if type(chunks) is not list or not 1 <= len(chunks) <= maximum_chunks:
        raise _fail("GIV031", f"{label}.chunks must be a bounded nonempty array")
    encoded: list[bytes] = []
    total = 0
    for index, chunk in enumerate(chunks):
        if type(chunk) is not str:
            raise _fail("GIV031", f"{label}.chunks[{index}] must be text")
        try:
            raw = chunk.encode("ascii")
        except UnicodeEncodeError as failure:
            raise _fail("GIV031", f"{label}.chunks[{index}] must be ASCII") from failure
        if index < len(chunks) - 1 and len(raw) != CHUNK_BYTES:
            raise _fail("GIV031", f"{label} has a short non-final chunk")
        if index == len(chunks) - 1 and not 1 <= len(raw) <= CHUNK_BYTES:
            raise _fail("GIV031", f"{label} has an empty or oversized final chunk")
        total += len(raw)
        if total > maximum:
            raise _fail("GIV031", f"{label} exceeds {maximum} bytes")
        encoded.append(raw)
    return b"".join(encoded)


def _validate_collection_v2(
    value: Any,
    *,
    genbank_source: bytes,
) -> dict[str, Any]:
    root = _exact_keys(
        value,
        {"format", "version", "compiler", "inputs", "members", "refget_seqcol", "artifact_sha256"},
        "Sequence Collection v2",
    )
    if (
        root["format"] != COLLECTION_FORMAT
        or type(root["version"]) is not int
        or root["version"] != COLLECTION_VERSION
        or root["compiler"] != COLLECTION_COMPILER
    ):
        raise _fail("GIV032", "unsupported Sequence Collection v2 identity")
    _assert_seal(root, "artifact_sha256", "Sequence Collection v2")
    inputs = _exact_keys(root["inputs"], {"kind", "profile", "source"}, "collection inputs")
    if inputs["kind"] != "genbank-flatfile" or inputs["profile"] != PROFILE:
        raise _fail("GIV032", "collection input kind or profile is invalid")
    source_value = _exact_keys(
        inputs["source"],
        {"encoding", "byte_length", "sha256", "storage"},
        "collection source",
    )
    if source_value["encoding"] != "ascii":
        raise _fail("GIV032", "collection source encoding must be ascii")
    byte_length = _integer(source_value["byte_length"], "collection source byte_length", 1)
    if byte_length > MAX_GENBANK_BYTES:
        raise _fail("GIV032", "collection source exceeds the GenBank byte limit")
    source = _storage_bytes(source_value["storage"], "collection source storage", maximum=MAX_GENBANK_BYTES)
    if (
        len(source) != byte_length
        or _assert_sha(source_value["sha256"], "collection source sha256")
        != hashlib.sha256(source).hexdigest()
        or source != genbank_source
    ):
        raise _fail("GIV032", "collection source storage differs from raw GenBank evidence")

    members_value = root["members"]
    if type(members_value) is not list or not 1 <= len(members_value) <= MAX_RECORDS:
        raise _fail("GIV033", "collection members must be a bounded nonempty array")
    normalized: list[dict[str, Any]] = []
    seen: set[str] = set()
    total_bases = 0
    previous_line_end = 0
    for index, member_value in enumerate(members_value):
        member = _exact_keys(member_value, {"record_id", "source_map", "sequence"}, f"collection member {index}")
        record_id = _reference(member["record_id"])
        if record_id in seen:
            raise _fail("GIV033", f"duplicate collection record_id {record_id!r}")
        seen.add(record_id)
        source_map = _exact_keys(
            member["source_map"],
            {"kind", "version", "line_start", "line_end"},
            f"collection member {index} source map",
        )
        if (
            source_map["kind"] != "genbank-origin-lines"
            or type(source_map["version"]) is not int
            or source_map["version"] != 1
        ):
            raise _fail("GIV033", "collection member source map identity is invalid")
        line_start = _integer(source_map["line_start"], "member line_start", 1)
        line_end = _integer(source_map["line_end"], "member line_end", line_start)
        if line_end > MAX_LINES or line_start <= previous_line_end:
            raise _fail("GIV033", "collection member source maps must be ordered and disjoint")
        previous_line_end = line_end
        sequence_value = _exact_keys(
            member["sequence"],
            {"alphabet", "normalization", "bases", "sha256", "refget_id", "storage"},
            f"collection member {index} sequence",
        )
        if sequence_value["alphabet"] != "IUPAC-DNA" or sequence_value["normalization"] != "uppercase":
            raise _fail("GIV034", "collection sequence alphabet or normalization is invalid")
        bases = _integer(sequence_value["bases"], "member sequence bases", 1)
        if bases > MAX_TOTAL_BASES:
            raise _fail("GIV034", "collection member exceeds the base limit")
        sequence = _storage_bytes(
            sequence_value["storage"],
            f"collection member {index} sequence storage",
            maximum=MAX_TOTAL_BASES,
        )
        if any(octet not in b"ACGTRYSWKMBDHVN" for octet in sequence):
            raise _fail("GIV034", "collection sequence storage is not uppercase IUPAC DNA")
        if (
            len(sequence) != bases
            or _assert_sha(sequence_value["sha256"], "member sequence sha256")
            != hashlib.sha256(sequence).hexdigest()
            or sequence_value["refget_id"] != _refget(sequence.decode("ascii"))
        ):
            raise _fail("GIV034", "collection sequence differs from its length or identities")
        total_bases += bases
        if total_bases > MAX_TOTAL_BASES:
            raise _fail("GIV034", "collection exceeds the total base limit")
        normalized.append(member)
    if root["refget_seqcol"] != _refget_seqcol(normalized):
        raise _fail("GIV035", "refget Sequence Collection identity does not match ordered members")
    return root


def _source_semantics(record: dict[str, Any], features: list[dict[str, Any]]) -> None:
    sources = [feature for feature in features if feature["key"] == "source"]
    if not sources:
        raise _fail("GIV018", f"record {record['version']} lacks source feature")
    coverage: list[tuple[int, int]] = []
    present_mol_types: list[str] = []
    missing_mol_type = False
    for source in sources:
        for segment in source["location"]["segments"]:
            if segment["reference"] is not None or segment["kind"] != "interval":
                raise _fail("GIV018", "source locations must be local base spans")
            coverage.append((segment["start"], segment["end"]))
        organisms = [entry["value"] for entry in source["qualifiers"] if entry["name"] == "organism"]
        mol_types = [entry["value"] for entry in source["qualifiers"] if entry["name"] == "mol_type"]
        if len(organisms) != 1 or not organisms[0]:
            raise _fail("GIV018", "each source feature requires one /organism")
        if len(mol_types) > 1 or (mol_types and not mol_types[0]):
            raise _fail("GIV018", "source /mol_type may occur once and must have a value")
        if mol_types:
            present_mol_types.append(mol_types[0])
        else:
            missing_mol_type = True
    covered = 0
    for start, end in sorted(coverage):
        if start > covered:
            raise _fail("GIV018", "source features leave bases uncovered")
        covered = max(covered, end)
    if covered != record["locus"]["length"]:
        raise _fail("GIV018", "source features do not span the sequence")
    if present_mol_types and (missing_mol_type or len(set(present_mol_types)) != 1):
        raise _fail("GIV018", "source /mol_type values are not uniformly present and equal")


def _emit_record(
    record: dict[str, Any],
    contexts: dict[str, dict[str, Any]],
    *,
    maximum_members: int,
) -> tuple[dict[str, Any], int]:
    local = contexts[record["version"]]
    locus = dict(record["locus"])
    locus["source_map"] = _position(locus["source_map"])
    features: list[dict[str, Any]] = []
    core = {
        "record_id": record["version"],
        "accessions": record["accessions"],
        "version": record["version"],
        "legacy_gi": record["legacy_gi"],
        "locus": locus,
        "definition": record["definition"],
        "keywords": record["keywords"],
        "source": record["source"],
        "organism": {"name": record["organism"], "lineage": record["lineage"]},
        "headers": [
            {
                "name": header["name"],
                "value": header["value"],
                "source_map": [_position(chunk) for chunk in header["chunks"]],
            }
            for header in record["headers"]
        ],
        "features": features,
        "relationships": [],
        "source_map": {"line_start": record["line_start"], "line_end": record["line_end"]},
    }
    record_members = _json_member_count(
        {**core, "record_ir_sha256": "0" * 64}
    )
    if record_members > maximum_members:
        raise _fail("GIV002", "record metadata exceeds the remaining JSON member budget")
    for raw in record["features"]:
        ast = _LocationParser(raw["location_text"], local, contexts).parse()
        segments = _segments(ast, local, contexts)
        identity = {
            "record_id": record["version"],
            "ordinal": raw["ordinal"],
            "key": raw["key"],
            "location": raw["location_text"],
            "qualifiers": [
                {"name": qualifier["name"], "value": qualifier["value"]}
                for qualifier in raw["qualifiers"]
            ],
        }
        emitted = {
            "feature_id": FEATURE_ID_PREFIX + _digest(identity),
            "ordinal": raw["ordinal"],
            "key": raw["key"],
            "location": {
                "text": raw["location_text"],
                "ast": ast,
                "segments": segments,
                "unresolved_references": sorted({
                    segment["reference"] for segment in segments
                    if segment["bounds_status"] == "unresolved" and segment["reference"] is not None
                }),
                "source_map": [_position(chunk) for chunk in raw["location_chunks"]],
            },
            "qualifiers": [
                {
                    "name": qualifier["name"],
                    "value": qualifier["value"],
                    "source_map": [_position(chunk) for chunk in qualifier["chunks"]],
                }
                for qualifier in raw["qualifiers"]
            ],
            "source_map": {"line_start": raw["line_start"], "line_end": raw["line_end"]},
        }
        contribution = 1 + _json_member_count(emitted)
        if record_members + contribution > maximum_members:
            raise _fail(
                "GIV002",
                f"source stream exceeds {MAX_JSON_MEMBERS} aggregate JSON members",
            )
        features.append(emitted)
        record_members += contribution
    _source_semantics(record, features)
    return {**core, "record_ir_sha256": _digest(core)}, record_members


def _replay(raw: bytes) -> dict[str, Any]:
    records = [_parse_record(record) for record in _split_records(_split_lines(raw))]
    feature_count = sum(len(record["features"]) for record in records)
    qualifier_count = sum(
        len(feature["qualifiers"])
        for record in records
        for feature in record["features"]
    )
    if feature_count > MAX_FEATURES:
        raise _fail("GIV013", f"source stream exceeds {MAX_FEATURES} aggregate features")
    if qualifier_count > MAX_QUALIFIERS:
        raise _fail("GIV013", f"source stream exceeds {MAX_QUALIFIERS} aggregate qualifiers")
    versions = [record["version"] for record in records]
    primaries = [record["accessions"][0] for record in records]
    if len(versions) != len(set(versions)) or len(primaries) != len(set(primaries)):
        raise _fail("GIV019", "record identities must be unique")
    contexts = {
        record["version"]: {"length": record["locus"]["length"], "topology": record["locus"]["topology"]}
        for record in records
    }
    collection = _sequence_collection(raw, records)
    emitted: list[dict[str, Any]] = []
    bio_ir = {
        "profile": PROFILE,
        "authority": dict(AUTHORITY),
        "sequence_collection_artifact_sha256": collection["artifact_sha256"],
        "coordinate_system": {"source": "1-based-closed", "normalized": "0-based-half-open"},
        "records": emitted,
    }
    member_skeleton = {
        "format": FORMAT,
        "version": VERSION,
        "profile": PROFILE,
        "authority": dict(AUTHORITY),
        "compiler": {**COMPILER, "passes": list(COMPILER["passes"])},
        "sequence_collection": collection,
        "bio_ir": bio_ir,
        "bio_ir_sha256": "0" * 64,
        "artifact_sha256": "0" * 64,
    }
    used_members = _json_member_count(member_skeleton)
    if used_members > MAX_JSON_MEMBERS:
        raise _fail(
            "GIV002",
            "collection and record metadata exceed the aggregate JSON member budget",
        )
    for record in records:
        remaining = MAX_JSON_MEMBERS - used_members - 1
        if remaining < 0:
            raise _fail(
                "GIV002",
                f"source stream exceeds {MAX_JSON_MEMBERS} aggregate JSON members",
            )
        replayed, record_members = _emit_record(
            record,
            contexts,
            maximum_members=remaining,
        )
        emitted.append(replayed)
        used_members += 1 + record_members
    core = {
        "format": FORMAT,
        "version": VERSION,
        "profile": PROFILE,
        "authority": dict(AUTHORITY),
        "compiler": {**COMPILER, "passes": list(COMPILER["passes"])},
        "sequence_collection": collection,
        "bio_ir": bio_ir,
        "bio_ir_sha256": _digest(bio_ir),
    }
    return {**core, "artifact_sha256": _digest(core)}


def _assert_seal(value: dict[str, Any], field: str, label: str) -> None:
    stored = _assert_sha(value.get(field), f"{label}.{field}")
    wanted = _digest({key: item for key, item in value.items() if key != field})
    if stored != wanted:
        raise _fail("GIV020", f"{label}.{field} does not match canonical content")


def validate_genbank(
    artifact: dict[str, Any],
    *,
    genbank_source: bytes,
) -> dict[str, Any]:
    """Independently replay original GenBank bytes into a sealed report."""

    _load_authority_manifest()
    if type(genbank_source) is not bytes or not genbank_source or len(genbank_source) > MAX_GENBANK_BYTES:
        raise _fail("GIV021", "GenBank evidence must be bounded nonempty bytes")
    _validate_tree(artifact, "GenBank artifact")
    if _canonical_size(artifact, MAX_ARTIFACT_BYTES) > MAX_ARTIFACT_BYTES:
        raise _fail("GIV021", "GenBank artifact exceeds canonical byte limit")
    root = _exact_keys(
        artifact,
        {
            "format", "version", "profile", "authority", "compiler",
            "sequence_collection", "bio_ir", "bio_ir_sha256", "artifact_sha256",
        },
        "GenBank artifact",
    )
    if (
        root["format"] != FORMAT or type(root["version"]) is not int or root["version"] != VERSION
        or root["profile"] != PROFILE or root["authority"] != AUTHORITY or root["compiler"] != COMPILER
    ):
        raise _fail("GIV022", "unsupported GenBank artifact identity")
    _assert_seal(root, "artifact_sha256", "GenBank artifact")
    source_sha = hashlib.sha256(genbank_source).hexdigest()
    collection = _validate_collection_v2(
        root["sequence_collection"],
        genbank_source=genbank_source,
    )
    bio_ir = _exact_keys(
        root["bio_ir"],
        {
            "profile",
            "authority",
            "sequence_collection_artifact_sha256",
            "coordinate_system",
            "records",
        },
        "Bio IR",
    )
    if _assert_sha(root["bio_ir_sha256"], "bio_ir_sha256") != _digest(bio_ir):
        raise _fail("GIV020", "bio_ir_sha256 does not match canonical Bio IR")
    if (
        bio_ir["profile"] != PROFILE
        or bio_ir["authority"] != AUTHORITY
        or bio_ir["coordinate_system"]
        != {"source": "1-based-closed", "normalized": "0-based-half-open"}
        or bio_ir["sequence_collection_artifact_sha256"]
        != collection["artifact_sha256"]
    ):
        raise _fail("GIV023", "Bio IR identity or Sequence Collection binding is invalid")

    expected = _replay(genbank_source)
    if not _exact_json_equal(root, expected):
        raise _fail("GIV024", "artifact does not match independent raw GenBank replay")

    collection = expected["sequence_collection"]
    records = expected["bio_ir"]["records"]
    features = [feature for record in records for feature in record["features"]]
    segments = [segment for feature in features for segment in feature["location"]["segments"]]
    unresolved = sorted({
        reference
        for feature in features
        for reference in feature["location"]["unresolved_references"]
    })
    report_core = {
        "format": REPORT_FORMAT,
        "version": REPORT_VERSION,
        "validator": dict(VALIDATOR),
        "inputs": {
            "genbank_source_sha256": source_sha,
            "genbank_artifact_sha256": expected["artifact_sha256"],
            "sequence_collection_artifact_sha256": collection["artifact_sha256"],
        },
        "checks": list(CHECKS),
        "replay": {
            "source_sha256": collection["inputs"]["source"]["sha256"],
            "sequence_collection_artifact_sha256": collection["artifact_sha256"],
            "refget_seqcol_digest": collection["refget_seqcol"]["digest"],
            "bio_ir_sha256": expected["bio_ir_sha256"],
            "genbank_artifact_sha256": expected["artifact_sha256"],
        },
        "summary": {
            "records": len(records),
            "source_bytes": len(genbank_source),
            "sequence_bases": sum(member["sequence"]["bases"] for member in collection["members"]),
            "sequence_chunks": sum(len(member["sequence"]["storage"]["chunks"]) for member in collection["members"]),
            "features": len(features),
            "segments": len(segments),
            "unresolved_references": unresolved,
        },
        "valid": True,
    }
    report = {**report_core, "report_sha256": _digest(report_core)}
    validate_genbank_report(report)
    return report


def validate_genbank_report(report: dict[str, Any]) -> dict[str, Any]:
    """Validate the closed report schema and its canonical seal."""

    _validate_tree(report, "GenBank validation report")
    root = _exact_keys(
        report,
        {"format", "version", "validator", "inputs", "checks", "replay", "summary", "valid", "report_sha256"},
        "GenBank validation report",
    )
    if (
        root["format"] != REPORT_FORMAT or type(root["version"]) is not int or root["version"] != REPORT_VERSION
        or root["validator"] != VALIDATOR or root["checks"] != CHECKS or root["valid"] is not True
    ):
        raise _fail("GIV025", "unsupported or unsuccessful validation report")
    inputs = _exact_keys(
        root["inputs"],
        {
            "genbank_source_sha256",
            "genbank_artifact_sha256",
            "sequence_collection_artifact_sha256",
        },
        "report inputs",
    )
    replay = _exact_keys(
        root["replay"],
        {
            "source_sha256",
            "sequence_collection_artifact_sha256",
            "refget_seqcol_digest",
            "bio_ir_sha256",
            "genbank_artifact_sha256",
        },
        "report replay",
    )
    for name, value in inputs.items():
        _assert_sha(value, f"report input digest {name}")
    for name in (
        "source_sha256",
        "sequence_collection_artifact_sha256",
        "bio_ir_sha256",
        "genbank_artifact_sha256",
    ):
        _assert_sha(replay[name], f"report replay digest {name}")
    if (
        type(replay["refget_seqcol_digest"]) is not str
        or re.fullmatch(r"[A-Za-z0-9_-]{32}", replay["refget_seqcol_digest"]) is None
    ):
        raise _fail("GIV025", "report replay refget_seqcol_digest is invalid")
    if (
        replay["source_sha256"] != inputs["genbank_source_sha256"]
        or replay["sequence_collection_artifact_sha256"]
        != inputs["sequence_collection_artifact_sha256"]
        or replay["genbank_artifact_sha256"]
        != inputs["genbank_artifact_sha256"]
    ):
        raise _fail("GIV025", "report replay bindings are inconsistent")
    summary = _exact_keys(
        root["summary"],
        {
            "records",
            "source_bytes",
            "sequence_bases",
            "sequence_chunks",
            "features",
            "segments",
            "unresolved_references",
        },
        "report summary",
    )
    _integer(summary["records"], "summary.records", 1, MAX_RECORDS)
    _integer(summary["source_bytes"], "summary.source_bytes", 1, MAX_GENBANK_BYTES)
    _integer(summary["sequence_bases"], "summary.sequence_bases", 1, MAX_TOTAL_BASES)
    _integer(summary["sequence_chunks"], "summary.sequence_chunks", 1, MAX_TOTAL_BASES)
    _integer(summary["features"], "summary.features", 0, MAX_FEATURES)
    _integer(summary["segments"], "summary.segments", 0, MAX_JSON_MEMBERS)
    if type(summary["unresolved_references"]) is not list:
        raise _fail("GIV025", "summary.unresolved_references must be an array")
    references = summary["unresolved_references"]
    if any(type(value) is not str for value in references):
        raise _fail("GIV025", "summary.unresolved_references must contain text")
    for value in references:
        _reference(value)
    if references != sorted(set(references)):
        raise _fail("GIV025", "summary.unresolved_references must be sorted unique accession.versions")
    _assert_seal(root, "report_sha256", "GenBank validation report")
    if (
        summary["sequence_bases"] < summary["records"]
        or summary["sequence_chunks"] < summary["records"]
        or summary["features"] < summary["records"]
        or summary["segments"] < summary["features"]
        or len(references) > summary["segments"]
    ):
        raise _fail("GIV025", "report summary contains impossible cross-field counts")
    _report_bytes(root)
    return root


def _read_regular(path: str | Path, label: str, maximum: int) -> bytes:
    source = Path(path)
    descriptor = -1
    try:
        inspected = source.lstat()
        if stat.S_ISLNK(inspected.st_mode) or not stat.S_ISREG(inspected.st_mode):
            raise _fail("GIV026", f"{label} must be a regular non-linked file")
        if inspected.st_size > maximum:
            raise _fail("GIV026", f"{label} exceeds {maximum} bytes")
        flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
        descriptor = os.open(source, flags)
        opened = os.fstat(descriptor)
        if not stat.S_ISREG(opened.st_mode) or (opened.st_dev, opened.st_ino) != (inspected.st_dev, inspected.st_ino):
            raise _fail("GIV026", f"{label} changed while being opened")
        chunks: list[bytes] = []
        total = 0
        while True:
            chunk = os.read(descriptor, min(64 * 1024, maximum + 1 - total))
            if not chunk:
                break
            chunks.append(chunk)
            total += len(chunk)
            if total > maximum:
                raise _fail("GIV026", f"{label} exceeds {maximum} bytes")
        finished = os.fstat(descriptor)
        for field in ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns"):
            if getattr(opened, field) != getattr(finished, field):
                raise _fail("GIV026", f"{label} changed while being read")
        if total != opened.st_size:
            raise _fail("GIV026", f"{label} length changed while being read")
        return b"".join(chunks)
    except INSDCValidationError:
        raise
    except OSError as failure:
        raise _fail("GIV026", f"cannot read {label}: {failure}") from failure
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _duplicate_safe_json(raw: bytes, label: str) -> dict[str, Any]:
    def pairs(values: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in values:
            if key in result:
                raise _fail("GIV027", f"{label} contains duplicate key {key!r}")
            result[key] = value
        return result

    def constant(token: str) -> None:
        raise _fail("GIV027", f"{label} contains non-finite number {token}")

    try:
        value = json.loads(raw.decode("utf-8"), object_pairs_hook=pairs, parse_constant=constant)
    except INSDCValidationError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError, RecursionError) as failure:
        raise _fail("GIV027", f"invalid {label} JSON: {failure}") from failure
    if type(value) is not dict:
        raise _fail("GIV027", f"{label} must contain a JSON object")
    return value


def _load_authority_manifest() -> dict[str, Any]:
    """Independently verify the offline standards evidence bound to this ABI."""

    raw = _read_regular(
        _AUTHORITY_MANIFEST_PATH,
        "GenBank authority manifest",
        MAX_AUTHORITY_MANIFEST_BYTES,
    )
    manifest = _duplicate_safe_json(raw, "GenBank authority manifest")
    _validate_tree(manifest, "GenBank authority manifest")
    root = _exact_keys(
        manifest,
        {"format", "version", "profile", "sources", "artifact_sha256"},
        "GenBank authority manifest",
    )
    claimed = _assert_sha(root["artifact_sha256"], "authority manifest artifact_sha256")
    if claimed != AUTHORITY_MANIFEST_SHA256 or _digest(
        {key: value for key, value in root.items() if key != "artifact_sha256"}
    ) != AUTHORITY_MANIFEST_SHA256:
        raise _fail("GIV029", "authority manifest differs from the compiler ABI content pin")
    if (
        root["format"] != "brainc.bio.authority-manifest"
        or type(root["version"]) is not int
        or root["version"] != 1
        or root["profile"] != PROFILE
        or type(root["sources"]) is not dict
    ):
        raise _fail("GIV029", "authority manifest identity is invalid")
    return root


def _report_bytes(report: dict[str, Any]) -> bytes:
    try:
        raw = (json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")
    except (TypeError, ValueError, UnicodeError) as failure:
        raise _fail("GIV028", f"report is not serializable JSON: {failure}") from failure
    if len(raw) > MAX_REPORT_BYTES:
        raise _fail("GIV028", "report exceeds output limit")
    return raw


def _atomic_write(path: str | Path, raw: bytes) -> None:
    target = Path(path)
    descriptor = -1
    parent_descriptor = -1
    temporary_name: str | None = None
    try:
        if target.name in {"", ".", ".."}:
            raise _fail("GIV028", "report path must name a file")
        target.parent.mkdir(parents=True, exist_ok=True)
        descriptor_relative = (
            os.name == "posix"
            and getattr(os, "O_DIRECTORY", 0) != 0
            and getattr(os, "O_NOFOLLOW", 0) != 0
            and all(
                operation in getattr(os, "supports_dir_fd", set())
                for operation in (os.open, os.stat, os.unlink)
            )
        )
        if not descriptor_relative:
            _portable_atomic_write(target, raw)
            return

        parent_inspected = target.parent.lstat()
        if stat.S_ISLNK(parent_inspected.st_mode) or not stat.S_ISDIR(
            parent_inspected.st_mode
        ):
            raise _fail(
                "GIV028",
                "report parent must be a regular non-linked directory",
            )
        parent_flags = os.O_RDONLY | os.O_CLOEXEC | os.O_DIRECTORY | os.O_NOFOLLOW
        parent_descriptor = os.open(target.parent, parent_flags)
        parent_opened = os.fstat(parent_descriptor)
        parent_identity = (parent_opened.st_dev, parent_opened.st_ino)
        if (
            not stat.S_ISDIR(parent_opened.st_mode)
            or parent_identity
            != (parent_inspected.st_dev, parent_inspected.st_ino)
        ):
            raise _fail("GIV028", "report parent changed while being opened")

        try:
            existing = os.stat(
                target.name,
                dir_fd=parent_descriptor,
                follow_symlinks=False,
            )
        except FileNotFoundError:
            existing = None
        if existing is not None and (
            stat.S_ISLNK(existing.st_mode) or not stat.S_ISREG(existing.st_mode)
        ):
            raise _fail(
                "GIV028",
                "report path must be absent or a regular non-linked file",
            )

        temporary_name = f".{target.name}.{secrets.token_hex(12)}"
        temporary_flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC | os.O_NOFOLLOW
        descriptor = os.open(
            temporary_name,
            temporary_flags,
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
                target.name,
                dir_fd=parent_descriptor,
                follow_symlinks=False,
            )
        except FileNotFoundError:
            current = None
        if current is not None and (
            stat.S_ISLNK(current.st_mode) or not stat.S_ISREG(current.st_mode)
        ):
            raise _fail(
                "GIV028",
                "report path must be absent or a regular non-linked file",
            )

        os.replace(
            temporary_name,
            target.name,
            src_dir_fd=parent_descriptor,
            dst_dir_fd=parent_descriptor,
        )
        temporary_name = None
        published = os.stat(
            target.name,
            dir_fd=parent_descriptor,
            follow_symlinks=False,
        )
        if (
            not stat.S_ISREG(published.st_mode)
            or (published.st_dev, published.st_ino)
            != (temporary_metadata.st_dev, temporary_metadata.st_ino)
        ):
            raise _fail("GIV028", "report path changed while being published")

        parent_finished = target.parent.lstat()
        if (
            stat.S_ISLNK(parent_finished.st_mode)
            or not stat.S_ISDIR(parent_finished.st_mode)
            or (parent_finished.st_dev, parent_finished.st_ino) != parent_identity
        ):
            raise _fail("GIV028", "report parent changed while being published")
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
    except INSDCValidationError:
        raise
    except (OSError, TypeError, NotImplementedError) as failure:
        raise _fail("GIV028", f"cannot write validation report: {failure}") from failure
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


def _portable_atomic_write(target: Path, raw: bytes) -> None:
    """Portable same-parent replacement under the documented trusted-parent boundary."""

    descriptor = -1
    temporary: str | None = None
    try:
        parent_inspected = target.parent.lstat()
        parent_identity = (parent_inspected.st_dev, parent_inspected.st_ino)
        if stat.S_ISLNK(parent_inspected.st_mode) or not stat.S_ISDIR(
            parent_inspected.st_mode
        ):
            raise _fail(
                "GIV028",
                "report parent must be a regular non-linked directory",
            )
        try:
            existing = target.lstat()
        except FileNotFoundError:
            existing = None
        if existing is not None and (
            stat.S_ISLNK(existing.st_mode) or not stat.S_ISREG(existing.st_mode)
        ):
            raise _fail(
                "GIV028",
                "report path must be absent or a regular non-linked file",
            )
        descriptor, temporary = tempfile.mkstemp(
            prefix=f".{target.name}.",
            dir=target.parent,
        )
        temporary_metadata = os.fstat(descriptor)
        with os.fdopen(descriptor, "wb") as stream:
            descriptor = -1
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        parent_current = target.parent.lstat()
        if (
            stat.S_ISLNK(parent_current.st_mode)
            or not stat.S_ISDIR(parent_current.st_mode)
            or (parent_current.st_dev, parent_current.st_ino) != parent_identity
        ):
            raise _fail("GIV028", "report parent changed while being published")
        try:
            current = target.lstat()
        except FileNotFoundError:
            current = None
        if current is not None and (
            stat.S_ISLNK(current.st_mode) or not stat.S_ISREG(current.st_mode)
        ):
            raise _fail(
                "GIV028",
                "report path must be absent or a regular non-linked file",
            )
        os.replace(temporary, target)
        temporary = None
        published = target.lstat()
        if (
            not stat.S_ISREG(published.st_mode)
            or (published.st_dev, published.st_ino)
            != (temporary_metadata.st_dev, temporary_metadata.st_ino)
        ):
            raise _fail("GIV028", "report path changed while being published")
        parent_finished = target.parent.lstat()
        if (
            stat.S_ISLNK(parent_finished.st_mode)
            or not stat.S_ISDIR(parent_finished.st_mode)
            or (parent_finished.st_dev, parent_finished.st_ino) != parent_identity
        ):
            raise _fail("GIV028", "report parent changed while being published")
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        if temporary is not None:
            try:
                os.unlink(temporary)
            except OSError:
                pass


def validate_genbank_paths(
    artifact_path: str | Path,
    genbank_path: str | Path,
    *,
    report_path: str | Path | None = None,
) -> dict[str, Any]:
    """Read stable path snapshots, validate them, and optionally save a report."""

    artifact = _duplicate_safe_json(_read_regular(artifact_path, "GenBank artifact", MAX_ARTIFACT_BYTES), "GenBank artifact")
    source = _read_regular(genbank_path, "GenBank source", MAX_GENBANK_BYTES)
    report = validate_genbank(artifact, genbank_source=source)
    if report_path is not None:
        _atomic_write(report_path, _report_bytes(report))
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Independently validate a compiled GenBank artifact against original bytes")
    parser.add_argument("genbank", help="original GenBank flat-file source")
    parser.add_argument("artifact", help="compiled GenBank artifact JSON")
    parser.add_argument(
        "-o",
        "--output",
        "--report",
        dest="output",
        help="write the sealed validation report",
    )
    arguments = parser.parse_args(argv)
    try:
        report = validate_genbank_paths(arguments.artifact, arguments.genbank, report_path=arguments.output)
    except INSDCValidationError as failure:
        parser.exit(1, f"validation failed: {failure}\n")
    if arguments.output is None:
        sys.stdout.buffer.write(_report_bytes(report))
    return 0


__all__ = [
    "CHECKS",
    "INSDCValidationError",
    "REPORT_FORMAT",
    "REPORT_VERSION",
    "VALIDATOR",
    "validate_genbank",
    "validate_genbank_paths",
    "validate_genbank_report",
]


if __name__ == "__main__":
    raise SystemExit(main())

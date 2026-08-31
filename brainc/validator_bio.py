"""Independent source-replay validator for the biological GFF3 frontend.

This module intentionally imports no compiler implementation.  It reconstructs
both the Sequence Collection artifact and the GFF3 BioIR from the supplied raw
sources, compares the reconstructed artifacts byte-for-byte under RFC 8785
canonicalization, and emits a sealed validation report.
"""

from __future__ import annotations

import base64
from decimal import Decimal, InvalidOperation
import gzip
import hashlib
import io
import json
import math
import re
from typing import Any


SAFE_INTEGER = 2**53 - 1

COLLECTION_FORMAT = "brain01.sequence-collection-ir"
COLLECTION_VERSION = 1
COLLECTION_COMPILER = {"name": "brainc-dna-collection", "version": "0.2.0"}
SEQUENCE_FORMAT = "brain01.sequence-ir"
SEQUENCE_VERSION = 2
SEQUENCE_COMPILER = {
    "name": "brainc-dna",
    "version": "0.3.0",
    "passes": [
        "parse-fasta",
        "parse-context",
        "validate-iupac",
        "resolve-reference",
        "canonicalize-sequence",
        "emit-sequence-ir",
    ],
}

BIO_FORMAT = "brainc.bio.gff3-ir"
BIO_VERSION = 1
BIO_COMPILER = {
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
REPORT_FORMAT = "brainc.bio.validation-report"
REPORT_VERSION = 1
VALIDATOR = {
    "name": "brainc-independent-bio-validator",
    "version": "0.1.0",
    "strategy": "independent-source-replay",
}

SOURCE_COORDINATES = "1-based-closed"
NORMALIZED_COORDINATES = "0-based-half-open"
IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")

MAX_SOURCE_INPUT_BYTES = 16 * 1024 * 1024
MAX_LOGICAL_DNA_BYTES = 64 * 1024 * 1024
MAX_SOURCE_ARTIFACT_BYTES = 16 * 1024 * 1024
MAX_GFF3_BYTES = 64 * 1024 * 1024
MAX_BIO_ARTIFACT_BYTES = 256 * 1024 * 1024
MAX_LINE_BYTES = 1024 * 1024
MAX_LINES = 1_000_000
MAX_SEQUENCE_LINE_BREAKS = 1_000_000
MAX_IDENTIFIER_BYTES = 256
MAX_SEQUENCE_RECORDS = 100_000
MAX_FEATURE_ROWS = 250_000
MAX_DIRECTIVES = 100_000
MAX_ATTRIBUTES_PER_ROW = 256
MAX_VALUES_PER_ATTRIBUTE = 10_000
MAX_RELATIONSHIPS = 1_000_000
MAX_TEXT_BYTES = 1024 * 1024
MAX_JSON_DEPTH = 32
MAX_JSON_MEMBERS = 5_000_000

_SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
_GFF_VERSION_RE = re.compile(r"3(?:\.[0-9]+){0,2}\Z")
_SCORE_RE = re.compile(
    r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z"
)
_POSITIVE_COORDINATE_RE = re.compile(r"[0-9]+\Z")
_DIRECTIVE_RE = re.compile(r"[a-z][a-z0-9-]*\Z")
_ATTRIBUTE_TAG_RE = re.compile(r"[A-Za-z][A-Za-z0-9_.:-]*\Z")
_SAFE_SEQID = frozenset(
    "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.:^*$@!+_?-|"
)
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
_NORMALIZED_MULTI_ATTRIBUTES = frozenset(
    {"Parent", "Alias", "Note", "Dbxref", "Ontology_term", "Derives_from"}
)
_CHECKS = [
    "bounded-closed-artifacts",
    "sequence-collection-source-replay",
    "sequence-artifact-source-maps",
    "sequence-and-refget-digests",
    "gff3-profile-source-replay",
    "coordinate-and-circular-normalization",
    "feature-identity-and-discontinuous-segments",
    "relationship-resolution-and-acyclicity",
    "external-input-bindings",
    "canonical-artifact-digests",
]


class BioValidationError(ValueError):
    """Raw evidence or an artifact fails independent BioIR validation."""


def _fail(code: str, detail: str) -> BioValidationError:
    return BioValidationError(f"{code}: {detail}")


def _keys(value: Any, expected: set[str], label: str) -> dict[str, Any]:
    if type(value) is not dict:
        raise _fail("BIV001", f"{label} must be an object")
    missing = sorted(expected - value.keys())
    extra = sorted(value.keys() - expected)
    if missing or extra:
        raise _fail(
            "BIV002",
            f"{label} keys invalid; missing={missing or 'none'}, unknown={extra or 'none'}",
        )
    return value


def _integer(
    value: Any,
    label: str,
    *,
    minimum: int = 0,
    maximum: int = SAFE_INTEGER,
) -> int:
    if type(value) is not int or value < minimum or value > maximum:
        raise _fail("BIV003", f"{label} must be an integer in [{minimum}, {maximum}]")
    return value


def _text(
    value: Any,
    label: str,
    *,
    empty: bool = False,
    maximum: int = MAX_TEXT_BYTES,
) -> str:
    if type(value) is not str or (not empty and not value):
        raise _fail("BIV004", f"{label} must be {'possibly empty ' if empty else ''}text")
    if any(0xD800 <= ord(character) <= 0xDFFF for character in value):
        raise _fail("BIV004", f"{label} contains a lone Unicode surrogate")
    if len(value.encode("utf-8")) > maximum:
        raise _fail("BIV004", f"{label} exceeds {maximum} UTF-8 bytes")
    return value


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _sha512t24u(value: bytes) -> str:
    return base64.urlsafe_b64encode(hashlib.sha512(value).digest()[:24]).decode("ascii")


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
    return sign + rendered + ("e+" if normalized_exponent >= 0 else "e") + str(
        normalized_exponent
    )


def _jcs(value: Any) -> str:
    if value is None:
        return "null"
    if type(value) is bool:
        return "true" if value else "false"
    if type(value) is str:
        return _jcs_string(value)
    if type(value) is int:
        if not -SAFE_INTEGER <= value <= SAFE_INTEGER:
            raise ValueError("unsafe integer")
        return str(value)
    if type(value) is float:
        return _jcs_number(value)
    if type(value) is list:
        return "[" + ",".join(_jcs(item) for item in value) + "]"
    if type(value) is dict:
        if any(type(key) is not str for key in value):
            raise TypeError("non-string object key")
        keys = sorted(value, key=lambda key: key.encode("utf-16be"))
        return "{" + ",".join(_jcs_string(key) + ":" + _jcs(value[key]) for key in keys) + "}"
    raise TypeError(f"unsupported value type {type(value).__name__}")


def _canonical_bytes(value: Any) -> bytes:
    try:
        return _jcs(value).encode("utf-8")
    except (TypeError, ValueError, UnicodeError) as failure:
        raise _fail("BIV005", f"value is not RFC 8785 canonical JSON: {failure}") from failure


def _digest(value: Any) -> str:
    return _sha256_bytes(_canonical_bytes(value))


def _validate_tree(value: Any, label: str, *, string_limit: int) -> None:
    members = 0
    stack: list[tuple[Any, int]] = [(value, 0)]
    while stack:
        current, depth = stack.pop()
        if depth > MAX_JSON_DEPTH:
            raise _fail("BIV006", f"{label} exceeds JSON depth {MAX_JSON_DEPTH}")
        if current is None or type(current) is bool:
            continue
        if type(current) is int:
            if not -SAFE_INTEGER <= current <= SAFE_INTEGER:
                raise _fail("BIV006", f"{label} contains an unsafe integer")
            continue
        if type(current) is float:
            raise _fail("BIV006", f"{label} contains an unexpected JSON number")
        if type(current) is str:
            _text(current, f"{label} string", empty=True, maximum=string_limit)
            continue
        if type(current) is list:
            members += len(current)
            stack.extend((item, depth + 1) for item in reversed(current))
        elif type(current) is dict:
            members += len(current)
            for key, item in reversed(list(current.items())):
                if type(key) is not str:
                    raise _fail("BIV006", f"{label} contains a non-string object key")
                stack.append((key, depth + 1))
                stack.append((item, depth + 1))
        else:
            raise _fail("BIV006", f"{label} contains unsupported {type(current).__name__}")
        if members > MAX_JSON_MEMBERS:
            raise _fail("BIV006", f"{label} exceeds {MAX_JSON_MEMBERS} JSON members")


def _validate_wire_size(value: dict[str, Any], label: str, maximum: int) -> None:
    """Bound the exact indented UTF-8 file form used by source producers."""

    total = 1
    try:
        encoder = json.JSONEncoder(
            indent=2,
            sort_keys=True,
            ensure_ascii=False,
            allow_nan=False,
        )
        for chunk in encoder.iterencode(value):
            total += len(chunk.encode("utf-8"))
            if total > maximum:
                raise _fail("BIV027", f"{label} exceeds its wire byte limit")
    except BioValidationError:
        raise
    except (TypeError, ValueError, UnicodeError) as failure:
        raise _fail("BIV027", f"{label} is not serializable I-JSON: {failure}") from failure


def _assert_digest(payload: dict[str, Any], field: str, label: str) -> None:
    stored = payload.get(field)
    if type(stored) is not str or _SHA256_RE.fullmatch(stored) is None:
        raise _fail("BIV007", f"{label}.{field} must be lowercase SHA-256")
    wanted = _digest({key: value for key, value in payload.items() if key != field})
    if stored != wanted:
        raise _fail("BIV007", f"{label}.{field} does not match canonical content")


def _canonical_equal(actual: Any, expected: Any, label: str) -> None:
    if _canonical_bytes(actual) != _canonical_bytes(expected):
        raise _fail("BIV008", f"{label} does not match independent source replay")


def _decompress_fasta(raw: bytes) -> tuple[bytes, str | None]:
    if type(raw) is not bytes or len(raw) > MAX_SOURCE_INPUT_BYTES:
        raise _fail(
            "BIV009",
            f"sequence source must be at most {MAX_SOURCE_INPUT_BYTES} bytes",
        )
    if not raw.startswith(b"\x1f\x8b"):
        return raw, None
    try:
        with gzip.GzipFile(fileobj=io.BytesIO(raw)) as stream:
            logical = stream.read(MAX_LOGICAL_DNA_BYTES + 1)
    except (EOFError, OSError) as failure:
        raise _fail("BIV009", f"malformed gzip FASTA: {failure}") from failure
    if len(logical) > MAX_LOGICAL_DNA_BYTES:
        raise _fail(
            "BIV009",
            f"logical sequence source exceeds {MAX_LOGICAL_DNA_BYTES} bytes",
        )
    return logical, _sha256_bytes(raw)


def _refget_id(sequence: str) -> str:
    raw = hashlib.sha512(sequence.upper().encode("ascii")).digest()[:24]
    return "SQ." + base64.urlsafe_b64encode(raw).decode("ascii").rstrip("=")


def _parse_fasta_record(raw: bytes) -> dict[str, Any]:
    try:
        source = raw.decode("utf-8")
    except UnicodeDecodeError as failure:
        raise _fail("BIV010", "FASTA must be UTF-8") from failure
    line_breaks = source.count("\n") + source.count("\r") - source.count("\r\n")
    line_breaks += sum(
        source.count(separator)
        for separator in ("\v", "\f", "\x1c", "\x1d", "\x1e", "\x85", "\u2028", "\u2029")
    )
    if line_breaks > MAX_SEQUENCE_LINE_BREAKS:
        raise _fail(
            "BIV010",
            f"FASTA exceeds {MAX_SEQUENCE_LINE_BREAKS} physical line breaks",
        )
    lines = source.splitlines()
    if not lines or not lines[0].startswith(">"):
        raise _fail("BIV010", "FASTA record must begin with a defline")
    defline = lines[0][1:]
    if not defline.strip():
        raise _fail("BIV010", "FASTA record identifier is empty")
    first = defline.split(maxsplit=1)
    record_id = first[0]
    description = first[1] if len(first) == 2 else ""
    if len(record_id.encode("utf-8")) > MAX_IDENTIFIER_BYTES:
        raise _fail(
            "BIV010",
            f"FASTA record identifier exceeds {MAX_IDENTIFIER_BYTES} UTF-8 bytes",
        )
    if len(description.encode("utf-8")) > MAX_TEXT_BYTES:
        raise _fail(
            "BIV010",
            f"FASTA description exceeds {MAX_TEXT_BYTES} UTF-8 bytes",
        )
    if any(character.isspace() for character in record_id):
        raise _fail("BIV010", "FASTA record identifier contains whitespace")
    if any(
        ord(character) < 0x20 or 0x7F <= ord(character) <= 0x9F
        for character in record_id
    ):
        raise _fail("BIV010", "FASTA record identifier contains a control character")
    if any(
        ord(character) < 0x20 or 0x7F <= ord(character) <= 0x9F
        for character in description
    ):
        raise _fail("BIV010", "FASTA description contains a control character")
    sequence_parts: list[str] = []
    segments: list[dict[str, int]] = []
    offset = 0
    for line_number, line in enumerate(lines[1:], start=2):
        if line.startswith(">"):
            raise _fail("BIV010", "FASTA member contains a second defline")
        if not line:
            continue
        if len(line.encode("utf-8")) > MAX_LINE_BYTES:
            raise _fail("BIV010", f"FASTA line {line_number} exceeds its byte limit")
        for column, character in enumerate(line, start=1):
            if character.isspace() or character.upper() not in IUPAC_DNA:
                raise _fail(
                    "BIV010",
                    f"invalid IUPAC DNA at line {line_number}, column {column}",
                )
        sequence_parts.append(line)
        segments.append({"line": line_number, "normalized_start": offset, "length": len(line)})
        offset += len(line)
    sequence = "".join(sequence_parts)
    if not sequence:
        raise _fail("BIV010", "FASTA sequence is empty")
    if len(sequence.encode("utf-8")) > MAX_TEXT_BYTES:
        raise _fail("BIV010", f"FASTA sequence exceeds {MAX_TEXT_BYTES} bytes")
    return {
        "record_id": record_id,
        "description": description,
        "sequence": sequence,
        "segments": segments,
    }


def _sequence_artifact(record_raw: bytes, record: dict[str, Any]) -> dict[str, Any]:
    sequence = record["sequence"]
    sequence_ir = {
        "record_id": record["record_id"],
        "description": record["description"],
        "sequence": sequence,
        "sequence_sha256": _sha256_bytes(sequence.encode("ascii")),
        "canonical_sha256": _sha256_bytes(sequence.upper().encode("ascii")),
        "refget_id": _refget_id(sequence),
        "reference": None,
        "provenance": [],
    }
    core = {
        "format": SEQUENCE_FORMAT,
        "version": SEQUENCE_VERSION,
        "compiler": SEQUENCE_COMPILER,
        "inputs": {"fasta_sha256": _sha256_bytes(record_raw), "context_sha256": None},
        "sequence_ir": sequence_ir,
        "ir_sha256": _digest(sequence_ir),
        "source_map": {"sequence_segments": record["segments"]},
    }
    return {**core, "artifact_sha256": _digest(core)}


def _replay_collection(
    raw: bytes, record_id: str | None = None
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    logical, compressed_sha256 = _decompress_fasta(raw)
    line_breaks = logical.count(b"\n") + logical.count(b"\r") - logical.count(b"\r\n")
    line_breaks += sum(
        logical.count(separator)
        for separator in (b"\v", b"\f", b"\x1c", b"\x1d", b"\x1e", b"\x85")
    )
    if line_breaks > MAX_SEQUENCE_LINE_BREAKS:
        raise _fail(
            "BIV011",
            f"FASTA exceeds {MAX_SEQUENCE_LINE_BREAKS} physical line breaks",
        )
    lines = logical.splitlines(keepends=True)
    raw_input = not logical.startswith(b">")
    if raw_input:
        if type(record_id) is not str or not record_id or any(
            character.isspace() for character in record_id
        ):
            raise _fail(
                "BIV011",
                "raw IUPAC collection replay requires a nonempty whitespace-free record_id",
            )
        if any(
            ord(character) < 0x20 or 0x7F <= ord(character) <= 0x9F
            for character in record_id
        ):
            raise _fail("BIV011", "raw record_id contains a control character")
        try:
            identifier_bytes = len(record_id.encode("utf-8"))
        except UnicodeEncodeError as failure:
            raise _fail("BIV011", "raw record_id must be valid UTF-8") from failure
        if identifier_bytes > MAX_IDENTIFIER_BYTES:
            raise _fail(
                "BIV011",
                f"raw record_id exceeds {MAX_IDENTIFIER_BYTES} UTF-8 bytes",
            )
        try:
            sequence = logical.decode("ascii")
        except UnicodeDecodeError as failure:
            raise _fail("BIV011", "raw sequence must be ASCII IUPAC DNA") from failure
        if not sequence or any(
            character.isspace() or character.upper() not in IUPAC_DNA
            for character in sequence
        ):
            raise _fail(
                "BIV011", "raw sequence must be nonempty whitespace-free IUPAC DNA"
            )
        if len(sequence) > MAX_TEXT_BYTES:
            raise _fail("BIV011", f"raw sequence exceeds {MAX_TEXT_BYTES} bytes")
        try:
            encoded_record = f">{record_id}\n{sequence}\n".encode("utf-8")
        except UnicodeEncodeError as failure:
            raise _fail("BIV011", "raw record_id must be valid UTF-8") from failure
        source_records: list[tuple[int | None, bytes]] = [(None, encoded_record)]
        input_kind = "raw-iupac"
    else:
        if record_id is not None:
            raise _fail("BIV011", "record_id is only valid for raw IUPAC replay")
        starts = [index for index, line in enumerate(lines) if line.startswith(b">")]
        if not starts or len(starts) > MAX_SEQUENCE_RECORDS:
            raise _fail("BIV011", f"FASTA exceeds {MAX_SEQUENCE_RECORDS} records")
        source_records = [
            (
                start,
                b"".join(
                    lines[
                        start : (
                            starts[position + 1]
                            if position + 1 < len(starts)
                            else len(lines)
                        )
                    ]
                ),
            )
            for position, start in enumerate(starts)
        ]
        input_kind = "fasta"
    members: list[dict[str, Any]] = []
    records: list[dict[str, Any]] = []
    seen: set[str] = set()
    for start, record_raw in source_records:
        record = _parse_fasta_record(record_raw)
        record_id = record["record_id"]
        if record_id in seen:
            raise _fail("BIV011", f"duplicate FASTA record {record_id!r}")
        seen.add(record_id)
        sequence_artifact = _sequence_artifact(record_raw, record)
        outer_segments = (
            [
                {
                    "line": 1,
                    "normalized_start": 0,
                    "length": len(record["sequence"]),
                }
            ]
            if start is None
            else [
                {
                    "line": segment["line"] + start,
                    "normalized_start": segment["normalized_start"],
                    "length": segment["length"],
                }
                for segment in record["segments"]
            ]
        )
        members.append(
            {
                "record_id": record_id,
                "ir_sha256": sequence_artifact["ir_sha256"],
                "artifact_sha256": sequence_artifact["artifact_sha256"],
                "input_source_map": {"sequence_segments": outer_segments},
                "sequence_artifact": sequence_artifact,
            }
        )
        records.append(
            {
                "record_id": record_id,
                "sequence": record["sequence"],
                "sequence_sha256": sequence_artifact["sequence_ir"]["sequence_sha256"],
                "sequence_artifact_sha256": sequence_artifact["artifact_sha256"],
                "ir_sha256": sequence_artifact["ir_sha256"],
                "refget_id": sequence_artifact["sequence_ir"]["refget_id"],
            }
        )
    collection_ir = {
        "members": [
            {"record_id": member["record_id"], "ir_sha256": member["ir_sha256"]}
            for member in members
        ]
    }
    level_2 = {
        "lengths": [len(record["sequence"]) for record in records],
        "names": [record["record_id"] for record in records],
        "sequences": [record["refget_id"] for record in records],
    }
    level_1 = {name: _sha512t24u(_canonical_bytes(value)) for name, value in level_2.items()}
    inherent = {name: level_1[name] for name in ("names", "sequences")}
    core = {
        "format": COLLECTION_FORMAT,
        "version": COLLECTION_VERSION,
        "compiler": COLLECTION_COMPILER,
        "inputs": {
            "kind": input_kind,
            "wrapper": "gzip" if compressed_sha256 is not None else None,
            "raw_sha256": _sha256_bytes(raw),
            "logical_sha256": _sha256_bytes(logical),
            "compressed_sha256": compressed_sha256,
        },
        "collection_ir": collection_ir,
        "collection_ir_sha256": _digest(collection_ir),
        "refget_seqcol": {
            "version": "1.0.0",
            "digest": _sha512t24u(_canonical_bytes(inherent)),
            "level_1": level_1,
            "level_2": level_2,
        },
        "members": members,
    }
    return {**core, "artifact_sha256": _digest(core)}, records


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
            if (
                index + 2 >= len(value)
                or value[index + 1] not in _HEX
                or value[index + 2] not in _HEX
            ):
                raise _fail("BIV012", f"{label} contains a malformed percent escape")
            octet = int(value[index + 1 : index + 3], 16)
            if allowed_escape_bytes is not None and octet not in allowed_escape_bytes:
                raise _fail("BIV012", f"{label} percent-encodes a character that must be literal")
            encoded.append(octet)
            index += 3
            continue
        if ord(character) < 0x20 or ord(character) == 0x7F:
            raise _fail("BIV012", f"{label} contains a literal control character")
        encoded.extend(character.encode("utf-8"))
        index += 1
    try:
        return _text(encoded.decode("utf-8"), label)
    except UnicodeDecodeError as failure:
        raise _fail("BIV012", f"{label} percent escapes are not UTF-8") from failure


def _decode_seqid(value: str, label: str) -> str:
    if not value or value == ".":
        raise _fail("BIV013", f"{label} must identify a sequence")
    index = 0
    while index < len(value):
        if value[index] == "%":
            if (
                index + 2 >= len(value)
                or value[index + 1] not in _HEX
                or value[index + 2] not in _HEX
            ):
                raise _fail("BIV012", f"{label} contains a malformed percent escape")
            octet = int(value[index + 1 : index + 3], 16)
            if octet < 128 and chr(octet) in _SAFE_SEQID:
                raise _fail("BIV012", f"{label} redundantly percent-encodes a permitted character")
            if octet <= 0x20 or octet == 0x7F:
                raise _fail("BIV012", f"{label} percent-encodes forbidden whitespace/control")
            index += 3
        elif value[index] not in _SAFE_SEQID:
            raise _fail("BIV013", f"{label} contains an unescaped character")
        else:
            index += 1
    decoded = _percent_decode(value, label, allowed_escape_bytes=None)
    if any(
        character.isspace()
        or ord(character) < 0x20
        or 0x7F <= ord(character) <= 0x9F
        for character in decoded
    ):
        raise _fail("BIV013", f"{label} decodes to whitespace or a control character")
    return decoded


def _parse_score(value: str, line: int) -> str | None:
    if value == ".":
        return None
    if _SCORE_RE.fullmatch(value) is None:
        raise _fail("BIV014", f"line {line} score is not a decimal literal")
    try:
        score = Decimal(value)
    except InvalidOperation as failure:
        raise _fail("BIV014", f"line {line} score is invalid") from failure
    if not score.is_finite():
        raise _fail("BIV014", f"line {line} score is not finite")
    return _text(value, f"line {line} score")


def _read_coordinate(token: str, line: int, field: str) -> int:
    if _POSITIVE_COORDINATE_RE.fullmatch(token) is None:
        raise _fail(
            "BIV022",
            f"line {line} {field} is not an ASCII positive-decimal integer",
        )
    try:
        coordinate = int(token, 10)
    except ValueError as failure:
        raise _fail("BIV022", f"line {line} {field} is too large to parse") from failure
    return _integer(coordinate, f"line {line} {field}", minimum=1)


def _parse_attributes(value: str, line: int) -> list[dict[str, Any]]:
    if value == ".":
        return []
    if not value:
        raise _fail("BIV015", f"line {line} attributes cannot be empty")
    fields = value.split(";")
    if fields[-1] == "":
        fields.pop()
    if not fields or any(field == "" for field in fields):
        raise _fail("BIV015", f"line {line} contains an empty attribute")
    if len(fields) > MAX_ATTRIBUTES_PER_ROW:
        raise _fail("BIV015", f"line {line} has too many attributes")
    attributes: list[dict[str, Any]] = []
    seen: set[str] = set()
    for field in fields:
        if field.count("=") != 1:
            raise _fail("BIV015", f"line {line} attribute must contain one '='")
        raw_tag, raw_values = field.split("=", 1)
        tag = _percent_decode(
            raw_tag,
            f"line {line} attribute tag",
            allowed_escape_bytes=_ATTRIBUTE_ESCAPE_BYTES,
        )
        if _ATTRIBUTE_TAG_RE.fullmatch(tag) is None:
            raise _fail("BIV015", f"line {line} attribute tag {tag!r} is not a GFF3 name")
        if tag in seen:
            raise _fail("BIV015", f"line {line} repeats attribute {tag!r}")
        seen.add(tag)
        if tag[0].isupper() and tag not in _RESERVED_ATTRIBUTES:
            raise _fail("BIV015", f"line {line} uses unsupported reserved attribute {tag!r}")
        parts = raw_values.split(",")
        if not parts or any(part == "" for part in parts):
            raise _fail("BIV015", f"line {line} attribute {tag!r} has an empty value")
        if len(parts) > MAX_VALUES_PER_ATTRIBUTE:
            raise _fail("BIV015", f"line {line} attribute {tag!r} has too many values")
        if any("&" in part for part in parts):
            raise _fail("BIV015", f"line {line} attribute {tag!r} contains unescaped '&'")
        values = [
            _percent_decode(
                part,
                f"line {line} attribute {tag!r}",
                allowed_escape_bytes=_ATTRIBUTE_ESCAPE_BYTES,
            )
            for part in parts
        ]
        if tag in _NORMALIZED_MULTI_ATTRIBUTES:
            values.sort(key=lambda item: item.encode("utf-8"))
            if len(values) != len(set(values)):
                raise _fail("BIV015", f"line {line} attribute {tag!r} repeats a value")
        if tag in {"ID", "Is_circular", "Target", "Gap"} and len(values) != 1:
            raise _fail("BIV015", f"line {line} attribute {tag!r} needs one value")
        if tag == "Is_circular" and values[0] not in {"true", "false"}:
            raise _fail("BIV015", f"line {line} Is_circular must be true or false")
        attributes.append({"tag": tag, "values": values})
    attributes.sort(key=lambda item: item["tag"].encode("utf-8"))
    return attributes


def _attribute(attributes: list[dict[str, Any]], tag: str) -> list[str]:
    for attribute in attributes:
        if attribute["tag"] == tag:
            return attribute["values"]
    return []


def _intervals(
    start: int,
    end: int,
    length: int,
    circular: bool,
    label: str,
) -> list[dict[str, int]]:
    if start > length:
        raise _fail("BIV016", f"{label} starts beyond sequence length {length}")
    span = end - start + 1
    if end <= length:
        return [{"start": start - 1, "end": end}]
    if not circular:
        raise _fail("BIV016", f"{label} crosses a non-circular sequence origin")
    if span > length:
        raise _fail("BIV016", f"{label} spans more than one circular traversal")
    first_start = start - 1
    overflow = first_start + span - length
    return [{"start": first_start, "end": length}, {"start": 0, "end": overflow}]


def _anonymous_id(row: dict[str, Any]) -> str:
    signature = {
        "seqid": row["seqid"],
        "source": row["source"],
        "type": row["type"],
        "strand": row["strand"],
        "segment": row["segment"],
    }
    return "anon:" + _digest(signature)


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
                raise _fail("BIV017", f"relationship cycle includes {target!r}")
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


def _bio_records(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    result = [
        {
            "seqid": record["record_id"],
            "length": len(record["sequence"]),
            "sequence_sha256": record["sequence_sha256"],
            "sequence_artifact_sha256": record["sequence_artifact_sha256"],
        }
        for record in records
    ]
    return sorted(result, key=lambda item: item["seqid"].encode("utf-8"))


def _replay_gff3(raw: bytes, records: list[dict[str, Any]]) -> dict[str, Any]:
    if type(raw) is not bytes or len(raw) > MAX_GFF3_BYTES:
        raise _fail("BIV018", f"GFF3 source must be at most {MAX_GFF3_BYTES} bytes")
    try:
        source = raw.decode("utf-8")
    except UnicodeDecodeError as failure:
        raise _fail("BIV018", "GFF3 source must be strict UTF-8") from failure
    if any(
        (ord(character) < 0x20 and character not in "\t\n\r")
        or 0x7F <= ord(character) <= 0x9F
        for character in source
    ):
        raise _fail("BIV018", "GFF3 source contains an unescaped control character")
    lines = source.split("\n")
    lines = [line[:-1] if line.endswith("\r") else line for line in lines]
    if any("\r" in line for line in lines):
        raise _fail("BIV018", "GFF3 source contains carriage return outside CRLF")
    if not lines or len(lines) > MAX_LINES:
        raise _fail("BIV018", "GFF3 source has an invalid line count")
    if not lines[0].startswith("##gff-version "):
        raise _fail("BIV019", "GFF version must be the first physical line")
    version = lines[0][len("##gff-version ") :]
    if _GFF_VERSION_RE.fullmatch(version) is None:
        raise _fail("BIV019", "unsupported GFF3 version declaration")

    bindings = _bio_records(records)
    lengths = {record["seqid"]: record["length"] for record in bindings}
    sequence_regions: list[dict[str, Any]] = []
    region_by_seqid: dict[str, tuple[int, int]] = {}
    directives: list[dict[str, Any]] = []
    rows: list[dict[str, Any]] = []

    for line_number, line in enumerate(lines[1:], start=2):
        if len(line.encode("utf-8")) > MAX_LINE_BYTES:
            raise _fail("BIV020", f"line {line_number} exceeds its byte limit")
        if not line:
            continue
        if line == "##FASTA" or line.startswith(">"):
            raise _fail("BIV020", "embedded FASTA is outside this external-sequence profile")
        if line.startswith("###"):
            if line != "###":
                raise _fail("BIV020", f"line {line_number} has content after '###'")
            continue
        if line.startswith("##"):
            body = line[2:]
            name, separator, value = body.partition(" ")
            if not name or (separator and (not value or value != value.strip())):
                raise _fail("BIV021", f"line {line_number} has malformed directive syntax")
            if name == "gff-version":
                raise _fail("BIV019", "GFF3 source repeats its version directive")
            if name == "sequence-region":
                parts = value.split()
                if len(parts) != 3:
                    raise _fail("BIV021", f"line {line_number} has malformed sequence-region")
                seqid = _decode_seqid(parts[0], f"line {line_number} sequence-region seqid")
                start = _read_coordinate(parts[1], line_number, "sequence-region start")
                end = _read_coordinate(parts[2], line_number, "sequence-region end")
                if end < start:
                    raise _fail("BIV021", "sequence-region end precedes start")
                if seqid not in lengths or end > lengths[seqid]:
                    raise _fail("BIV021", f"sequence-region {seqid!r} is outside its sequence")
                if seqid in region_by_seqid:
                    raise _fail("BIV021", f"duplicate sequence-region {seqid!r}")
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
                raise _fail("BIV021", f"unsupported reserved directive {name!r}")
            if len(directives) >= MAX_DIRECTIVES:
                raise _fail("BIV021", "too many application directives")
            directives.append(
                {
                    "line": line_number,
                    "name": name,
                    "value": _text(value, f"line {line_number} directive", empty=True),
                }
            )
            continue
        if line.startswith("#"):
            continue
        if len(rows) >= MAX_FEATURE_ROWS:
            raise _fail("BIV022", "too many feature rows")
        columns = line.split("\t")
        if len(columns) != 9:
            raise _fail("BIV022", f"line {line_number} does not have nine columns")
        (
            raw_seqid,
            raw_source,
            raw_type,
            raw_start,
            raw_end,
            raw_score,
            strand,
            raw_phase,
            raw_attributes,
        ) = columns
        seqid = _decode_seqid(raw_seqid, f"line {line_number} seqid")
        if seqid not in lengths:
            raise _fail("BIV022", f"line {line_number} refers to unknown sequence {seqid!r}")
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
            raise _fail("BIV022", f"line {line_number} feature type is absent")
        feature_type = _percent_decode(
            raw_type,
            f"line {line_number} type",
            allowed_escape_bytes=_GENERAL_ESCAPE_BYTES,
        )
        start = _read_coordinate(raw_start, line_number, "start")
        end = _read_coordinate(raw_end, line_number, "end")
        if end < start:
            raise _fail("BIV022", f"line {line_number} end precedes start")
        if strand not in {"+", "-", ".", "?"}:
            raise _fail("BIV022", f"line {line_number} has an invalid strand")
        if feature_type == "CDS":
            if raw_phase not in {"0", "1", "2"}:
                raise _fail("BIV022", f"line {line_number} CDS phase is invalid")
            phase: int | None = int(raw_phase)
        else:
            if raw_phase != ".":
                raise _fail("BIV022", f"line {line_number} non-CDS phase is present")
            phase = None
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
                "attributes": _parse_attributes(raw_attributes, line_number),
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
        if region is not None and any(
            interval["start"] < region[0] - 1 or interval["end"] > region[1]
            for interval in intervals
        ):
            raise _fail("BIV023", f"line {row['line']} lies outside sequence-region")
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
            raise _fail("BIV024", "anonymous feature identity collision")
        if any(entity[key] != row[key] for key in ("seqid", "source", "type", "strand")):
            raise _fail("BIV024", f"discontinuous ID {declared_id!r} changes core fields")
        entity["segments"].append(row["segment"])

    features = sorted(entities.values(), key=lambda item: item["entity_id"].encode("utf-8"))
    for feature in features:
        feature["segments"].sort(key=lambda item: item["line"])

    known_ids = set(entities)
    relationship_set: set[tuple[str, str, str]] = set()
    for feature in features:
        for segment in feature["segments"]:
            for tag, kind in (("Parent", "parent"), ("Derives_from", "derives-from")):
                for target_id in _attribute(segment["attributes"], tag):
                    target = "id:" + target_id
                    if target not in known_ids:
                        raise _fail("BIV025", f"line {segment['line']} has dangling {tag}")
                    relationship_set.add((feature["entity_id"], kind, target))
                    if len(relationship_set) > MAX_RELATIONSHIPS:
                        raise _fail("BIV025", "too many feature relationships")
    relationships = [
        {"source": source_id, "kind": kind, "target": target}
        for source_id, kind, target in sorted(
            relationship_set,
            key=lambda item: tuple(part.encode("utf-8") for part in item),
        )
    ]
    _cycle_check(relationships)
    return {
        "gff3_version": version,
        "coordinate_system": {
            "source": SOURCE_COORDINATES,
            "normalized": NORMALIZED_COORDINATES,
        },
        "records": bindings,
        "sequence_regions": sequence_regions,
        "directives": directives,
        "features": features,
        "relationships": relationships,
    }


def _reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise _fail("BIV026", f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def _reject_constant(value: str) -> None:
    raise _fail("BIV026", f"non-finite JSON number {value!r}")


def _source_metadata(
    value: bytes | dict[str, Any] | None,
    fasta_raw: bytes,
    records: list[dict[str, Any]],
) -> str | None:
    if value is None:
        return None
    if type(value) is bytes:
        if len(value) > MAX_TEXT_BYTES:
            raise _fail("BIV026", "source metadata exceeds its byte limit")
        try:
            metadata = json.loads(
                value.decode("utf-8"),
                object_pairs_hook=_reject_duplicates,
                parse_constant=_reject_constant,
            )
        except BioValidationError:
            raise
        except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as failure:
            raise _fail("BIV026", f"invalid source metadata JSON: {failure}") from failure
        evidence_sha256 = _sha256_bytes(value)
    elif type(value) is dict:
        metadata = value
        _validate_tree(metadata, "source metadata", string_limit=MAX_TEXT_BYTES)
        evidence_sha256 = _sha256_bytes(_canonical_bytes(value))
    else:
        raise _fail("BIV026", "source metadata must be bytes, an object, or null")
    _validate_tree(metadata, "source metadata", string_limit=MAX_TEXT_BYTES)
    metadata = _keys(
        metadata,
        {
            "accession",
            "bases",
            "retrieved",
            "source",
            "fasta_sha256",
            "normalized_sequence_sha256",
        },
        "source metadata",
    )
    accession = _text(metadata["accession"], "source metadata accession")
    _text(metadata["retrieved"], "source metadata retrieved")
    _text(metadata["source"], "source metadata source")
    bases = _integer(metadata["bases"], "source metadata bases", minimum=1)
    if type(metadata["fasta_sha256"]) is not str or _SHA256_RE.fullmatch(
        metadata["fasta_sha256"]
    ) is None:
        raise _fail("BIV026", "source metadata fasta_sha256 is invalid")
    if type(metadata["normalized_sequence_sha256"]) is not str or _SHA256_RE.fullmatch(
        metadata["normalized_sequence_sha256"]
    ) is None:
        raise _fail("BIV026", "source metadata normalized_sequence_sha256 is invalid")
    matching = [record for record in records if record["record_id"] == accession]
    if len(matching) != 1:
        raise _fail("BIV026", "source metadata accession does not identify one FASTA member")
    record = matching[0]
    if bases != len(record["sequence"]):
        raise _fail("BIV026", "source metadata base count does not match FASTA")
    if metadata["fasta_sha256"] != _sha256_bytes(fasta_raw):
        raise _fail("BIV026", "source metadata FASTA digest does not match raw evidence")
    if metadata["normalized_sequence_sha256"] != _sha256_bytes(
        record["sequence"].upper().encode("ascii")
    ):
        raise _fail("BIV026", "source metadata normalized sequence digest does not match FASTA")
    return evidence_sha256


def _expected_bio_artifact(
    collection: dict[str, Any],
    records: list[dict[str, Any]],
    gff3_source: bytes,
) -> dict[str, Any]:
    bio_ir = _replay_gff3(gff3_source, records)
    inputs = {
        "gff3_sha256": _sha256_bytes(gff3_source),
        "sequence_collection": {
            "format": COLLECTION_FORMAT,
            "version": COLLECTION_VERSION,
            "artifact_sha256": collection["artifact_sha256"],
            "collection_ir_sha256": collection["collection_ir_sha256"],
        },
    }
    core = {
        "format": BIO_FORMAT,
        "version": BIO_VERSION,
        "compiler": BIO_COMPILER,
        "inputs": inputs,
        "bio_ir": bio_ir,
        "bio_ir_sha256": _digest(bio_ir),
    }
    return {**core, "artifact_sha256": _digest(core)}


def validate_bio_chain(
    bio_artifact: dict[str, Any],
    sequence_collection: dict[str, Any],
    *,
    fasta_source: bytes,
    gff3_source: bytes,
    source_metadata: bytes | dict[str, Any] | None = None,
    record_id: str | None = None,
) -> dict[str, Any]:
    """Replay Sequence Collection source bytes and GFF3 into a sealed report.

    The raw sources are mandatory because internally consistent hashes alone do
    not establish that an artifact is an honest compilation of external input.
    ``record_id`` is required only when ``fasta_source`` contains raw IUPAC DNA
    rather than FASTA. Optional source metadata is evidence only; it never
    changes compilation.
    """

    _validate_tree(
        sequence_collection,
        "sequence collection",
        string_limit=MAX_TEXT_BYTES,
    )
    _validate_tree(bio_artifact, "BioIR artifact", string_limit=MAX_TEXT_BYTES)
    if len(_canonical_bytes(sequence_collection)) > MAX_SOURCE_ARTIFACT_BYTES:
        raise _fail("BIV027", "sequence collection artifact exceeds its canonical byte limit")
    if len(_canonical_bytes(bio_artifact)) > MAX_BIO_ARTIFACT_BYTES:
        raise _fail("BIV027", "BioIR artifact exceeds its canonical byte limit")
    _validate_wire_size(
        sequence_collection,
        "sequence collection artifact",
        MAX_SOURCE_ARTIFACT_BYTES,
    )
    _validate_wire_size(
        bio_artifact,
        "BioIR artifact",
        MAX_BIO_ARTIFACT_BYTES,
    )

    collection_root = _keys(
        sequence_collection,
        {
            "format",
            "version",
            "compiler",
            "inputs",
            "collection_ir",
            "collection_ir_sha256",
            "refget_seqcol",
            "members",
            "artifact_sha256",
        },
        "sequence collection",
    )
    if (
        collection_root["format"] != COLLECTION_FORMAT
        or type(collection_root["version"]) is not int
        or collection_root["version"] != COLLECTION_VERSION
    ):
        raise _fail("BIV028", "unsupported sequence collection ABI")
    if collection_root["compiler"] != COLLECTION_COMPILER:
        raise _fail("BIV028", "unsupported sequence collection compiler identity")
    _assert_digest(collection_root, "artifact_sha256", "sequence collection")
    if collection_root.get("collection_ir_sha256") != _digest(collection_root["collection_ir"]):
        raise _fail("BIV028", "collection_ir_sha256 does not match collection_ir")

    bio_root = _keys(
        bio_artifact,
        {
            "format",
            "version",
            "compiler",
            "inputs",
            "bio_ir",
            "bio_ir_sha256",
            "artifact_sha256",
        },
        "BioIR artifact",
    )
    if (
        bio_root["format"] != BIO_FORMAT
        or type(bio_root["version"]) is not int
        or bio_root["version"] != BIO_VERSION
    ):
        raise _fail("BIV029", "unsupported BioIR ABI")
    if bio_root["compiler"] != BIO_COMPILER:
        raise _fail("BIV029", "unsupported BioIR compiler identity")
    _assert_digest(bio_root, "artifact_sha256", "BioIR artifact")
    if bio_root.get("bio_ir_sha256") != _digest(bio_root["bio_ir"]):
        raise _fail("BIV029", "bio_ir_sha256 does not match BioIR")

    expected_collection, records = _replay_collection(fasta_source, record_id)
    _canonical_equal(collection_root, expected_collection, "sequence collection")
    metadata_sha256 = _source_metadata(source_metadata, fasta_source, records)
    expected_bio = _expected_bio_artifact(expected_collection, records, gff3_source)

    inputs = _keys(bio_root["inputs"], {"gff3_sha256", "sequence_collection"}, "BioIR inputs")
    if inputs["gff3_sha256"] != _sha256_bytes(gff3_source):
        raise _fail("BIV030", "BioIR GFF3 digest does not match raw source")
    binding = _keys(
        inputs["sequence_collection"],
        {"format", "version", "artifact_sha256", "collection_ir_sha256"},
        "BioIR sequence binding",
    )
    if binding != expected_bio["inputs"]["sequence_collection"]:
        raise _fail("BIV030", "BioIR sequence collection binding is invalid")
    _canonical_equal(bio_root, expected_bio, "BioIR artifact")

    bio_ir = expected_bio["bio_ir"]
    segments = sum(len(feature["segments"]) for feature in bio_ir["features"])
    record_lengths = {record["seqid"]: record["length"] for record in bio_ir["records"]}
    circular_records = sorted(
        {
            feature["seqid"]
            for feature in bio_ir["features"]
            for segment in feature["segments"]
            if _attribute(segment["attributes"], "Is_circular") == ["true"]
            and segment["raw_start"] == 1
            and segment["raw_end"] == record_lengths[feature["seqid"]]
        },
        key=lambda item: item.encode("utf-8"),
    )
    report_core = {
        "format": REPORT_FORMAT,
        "version": REPORT_VERSION,
        "validator": dict(VALIDATOR),
        "inputs": {
            "sequence_collection_artifact_sha256": expected_collection["artifact_sha256"],
            "bio_artifact_sha256": expected_bio["artifact_sha256"],
            "fasta_sha256": _sha256_bytes(fasta_source),
            "gff3_sha256": _sha256_bytes(gff3_source),
            "source_metadata_sha256": metadata_sha256,
        },
        "checks": list(_CHECKS),
        "replay": {
            "sequence_collection_sha256": expected_collection["artifact_sha256"],
            "collection_ir_sha256": expected_collection["collection_ir_sha256"],
            "bio_ir_sha256": expected_bio["bio_ir_sha256"],
            "bio_artifact_sha256": expected_bio["artifact_sha256"],
        },
        "summary": {
            "records": len(records),
            "sequence_bases": sum(len(record["sequence"]) for record in records),
            "features": len(bio_ir["features"]),
            "segments": segments,
            "relationships": len(bio_ir["relationships"]),
            "circular_records": circular_records,
            "source_metadata_verified": metadata_sha256 is not None,
        },
        "valid": True,
    }
    report = {**report_core, "report_sha256": _digest(report_core)}
    validate_bio_report(report)
    return report


def _report_sha(value: Any, label: str, *, nullable: bool = False) -> str | None:
    if nullable and value is None:
        return None
    if type(value) is not str or _SHA256_RE.fullmatch(value) is None:
        raise _fail("BIV031", f"{label} must be lowercase SHA-256")
    return value


def validate_bio_report(report: dict[str, Any]) -> dict[str, Any]:
    """Validate the closed shape and cryptographic seal of a report."""

    _validate_tree(report, "validation report", string_limit=MAX_TEXT_BYTES)
    root = _keys(
        report,
        {
            "format",
            "version",
            "validator",
            "inputs",
            "checks",
            "replay",
            "summary",
            "valid",
            "report_sha256",
        },
        "validation report",
    )
    if (
        root["format"] != REPORT_FORMAT
        or type(root["version"]) is not int
        or root["version"] != REPORT_VERSION
    ):
        raise _fail("BIV031", "unsupported validation report ABI")
    if root["validator"] != VALIDATOR or root["checks"] != _CHECKS or root["valid"] is not True:
        raise _fail("BIV031", "validation report identity, checks, or result is invalid")
    inputs = _keys(
        root["inputs"],
        {
            "sequence_collection_artifact_sha256",
            "bio_artifact_sha256",
            "fasta_sha256",
            "gff3_sha256",
            "source_metadata_sha256",
        },
        "validation report inputs",
    )
    for name in (
        "sequence_collection_artifact_sha256",
        "bio_artifact_sha256",
        "fasta_sha256",
        "gff3_sha256",
    ):
        _report_sha(inputs[name], f"validation report inputs.{name}")
    _report_sha(
        inputs["source_metadata_sha256"],
        "validation report inputs.source_metadata_sha256",
        nullable=True,
    )
    replay = _keys(
        root["replay"],
        {
            "sequence_collection_sha256",
            "collection_ir_sha256",
            "bio_ir_sha256",
            "bio_artifact_sha256",
        },
        "validation report replay",
    )
    for name, value in replay.items():
        _report_sha(value, f"validation report replay.{name}")
    if replay["sequence_collection_sha256"] != inputs["sequence_collection_artifact_sha256"]:
        raise _fail("BIV031", "validation report collection replay digest is inconsistent")
    if replay["bio_artifact_sha256"] != inputs["bio_artifact_sha256"]:
        raise _fail("BIV031", "validation report BioIR replay digest is inconsistent")
    summary = _keys(
        root["summary"],
        {
            "records",
            "sequence_bases",
            "features",
            "segments",
            "relationships",
            "circular_records",
            "source_metadata_verified",
        },
        "validation report summary",
    )
    _integer(summary["records"], "summary.records", minimum=1)
    _integer(summary["sequence_bases"], "summary.sequence_bases", minimum=1)
    _integer(summary["features"], "summary.features")
    _integer(summary["segments"], "summary.segments")
    _integer(summary["relationships"], "summary.relationships")
    if type(summary["circular_records"]) is not list:
        raise _fail("BIV031", "summary.circular_records must be an array")
    circular = [_text(value, "summary.circular_records value") for value in summary["circular_records"]]
    if circular != sorted(circular, key=lambda item: item.encode("utf-8")) or len(circular) != len(
        set(circular)
    ):
        raise _fail("BIV031", "summary.circular_records must be sorted and unique")
    if type(summary["source_metadata_verified"]) is not bool:
        raise _fail("BIV031", "summary.source_metadata_verified must be boolean")
    _assert_digest(root, "report_sha256", "validation report")
    return root


__all__ = [
    "BioValidationError",
    "validate_bio_chain",
    "validate_bio_report",
]

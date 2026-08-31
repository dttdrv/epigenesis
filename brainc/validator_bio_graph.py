"""Independent validator for the GFF3 BioIR feature-graph bridge.

The validator deliberately has no dependency on the compiler package.  It
validates the supplied Sequence Collection IR and BioIR, reconstructs every
backend dictionary and tensor, reconstructs the complete v2 development
chain, and then requires byte-identity at the RFC 8785 content boundary.
"""

from __future__ import annotations

import base64
from decimal import Decimal, InvalidOperation
import hashlib
import json
import math
import os
from pathlib import Path
import re
import stat
import struct
from typing import Any


SAFE_INTEGER = 2**53 - 1
MAX_U64 = 2**64 - 1
# The bridge consumes the v2 wire form, whose serialized source, BioIR, and
# bundle artifacts are each capped at 16 MiB.  The standalone BioIR frontend's
# larger internal ceiling does not apply after admission to this lowering ABI.
MAX_SOURCE_INPUT_BYTES = 16 * 1024 * 1024
MAX_SOURCE_JSON_BYTES = MAX_SOURCE_INPUT_BYTES
MAX_ARTIFACT_BYTES = MAX_SOURCE_JSON_BYTES
MAX_LOGICAL_DNA_BYTES = 64 * 1024 * 1024
MAX_SOURCE_STRING_BYTES = 1 * 1024 * 1024
MAX_SOURCE_IDENTIFIER_BYTES = 256
MAX_SOURCE_RECORDS = 100_000
MAX_SOURCE_LINE = 1_000_001
MAX_JSON_DEPTH = 64
MAX_JSON_MEMBERS = 1_000_000
MAX_SOURCE_JSON_MEMBERS = MAX_JSON_MEMBERS
MAX_TEXT_BYTES = 1024 * 1024
MAX_FEATURE_ROWS = 250_000
MAX_RELATIONSHIPS = 1_000_000
READ_CHUNK_BYTES = 64 * 1024
SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
SCORE_RE = re.compile(
    r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z"
)
ATTRIBUTE_RE = re.compile(r"[A-Za-z][A-Za-z0-9_.:-]*\Z")
DIRECTIVE_RE = re.compile(r"[a-z][a-z0-9-]*\Z")
IUPAC = frozenset("ACGTRYSWKMBDHVN")

BACKEND_ID = "io.github.dttdrv.epigenesis.bio.feature-graph"
BACKEND_VERSION = 1
SEMANTICS_FORMAT = "brainc.bio.feature-graph-semantics"
RECORD_FORMAT = "brainc.bio.feature-graph-compilation"
BUNDLE_FORMAT = "brainc.bio.feature-graph-bundle"
DEV_DOMAIN = "io.github.dttdrv.epigenesis.dev"
OP_UNIT_CREATE = f"{DEV_DOMAIN}.unit.create"
OP_EDGE_CREATE = f"{DEV_DOMAIN}.edge.create"

# Literal rather than imported or recomputed from compiler data: the independent
# validator commits to the normative backend specification published by v0.6
# while separately implementing every formula and lowering rule below.
BACKEND_SPEC_SHA256 = "7b1287e4289c3b8e6328fed05dea22eb0b416fee0c2ab5b6a14635e5ccaec5e0"

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
COLLECTION_COMPILER = {"name": "brainc-dna-collection", "version": "0.2.0"}
MODULE_PRODUCER = {
    "name": "brainc",
    "version": "0.5.0",
    "passes": [
        "validate-sequence-source",
        "bind-provider-tensor-contract",
        "validate-target-contract",
        "validate-lowering-policy",
        "type-check-tensors",
        "validate-development-operations",
        "compute-exact-budgets",
        "emit-development-module",
    ],
}
BIO_GRAPH_PRODUCER = {
    "name": "brainc-bio-feature-graph",
    "version": "0.6.0",
    "passes": [
        "validate-bound-gff3-bio-ir",
        "construct-complete-normalized-bio-ir-sidecar",
        "encode-feature-tensors",
        "encode-explicit-relationship-tensors",
        "compile-development-module",
        "bind-compilation-record",
    ],
}

RELATIONSHIP_KINDS = ("derives-from", "parent")
STRANDS = ("+", "-", ".", "?")
FEATURE_FIELDS: tuple[tuple[str, str], ...] = (
    ("attribute_count", "u64"),
    ("declared_id", "bool"),
    ("interval_count", "u64"),
    ("location_max_end", "u64"),
    ("location_min_start", "u64"),
    ("phase_mask", "u64"),
    ("segment_bases", "u64"),
    ("segment_count", "u64"),
    ("seqid_code", "u64"),
    ("source_code", "u64"),
    ("strand_code", "u64"),
    ("type_code", "u64"),
    ("wraps_origin", "bool"),
)
RESERVED_ATTRIBUTES = frozenset(
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
CANONICAL_MULTI_VALUE_ATTRIBUTES = frozenset(
    {"Parent", "Alias", "Note", "Dbxref", "Ontology_term", "Derives_from"}
)


class BioGraphValidationError(ValueError):
    """An independently checked feature-graph invariant failed."""


def _fail(detail: str) -> BioGraphValidationError:
    return BioGraphValidationError(f"BIOGRAPHVAL001: {detail}")


def _keys(value: Any, expected: set[str], label: str) -> dict[str, Any]:
    if type(value) is not dict:
        raise _fail(f"{label} must be an object")
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
    minimum: int = 0,
    maximum: int = SAFE_INTEGER,
) -> int:
    if type(value) is not int or value < minimum or value > maximum:
        raise _fail(f"{label} must be an integer in [{minimum}, {maximum}]")
    return value


def _text(value: Any, label: str, *, empty: bool = False) -> str:
    if type(value) is not str or (not empty and not value):
        raise _fail(f"{label} must be {'possibly empty ' if empty else 'nonempty '}text")
    if any(0xD800 <= ord(character) <= 0xDFFF for character in value):
        raise _fail(f"{label} contains a lone Unicode surrogate")
    if len(value.encode("utf-8")) > MAX_TEXT_BYTES:
        raise _fail(f"{label} exceeds the text byte ceiling")
    return value


def _sha(value: Any, label: str) -> str:
    if type(value) is not str or SHA256_RE.fullmatch(value) is None:
        raise _fail(f"{label} must be a lowercase SHA-256 digest")
    return value


def _utf8_sorted(values: list[str], label: str, *, unique: bool = True) -> None:
    if values != sorted(values, key=lambda item: item.encode("utf-8")):
        raise _fail(f"{label} must be sorted by UTF-8 bytes")
    if unique and len(values) != len(set(values)):
        raise _fail(f"{label} contains duplicates")


def _validate_tree(
    value: Any,
    label: str,
    *,
    maximum_depth: int = MAX_JSON_DEPTH,
    maximum_members: int = MAX_JSON_MEMBERS,
    maximum_string_bytes: int = MAX_TEXT_BYTES,
) -> None:
    members = 0
    stack: list[tuple[Any, int]] = [(value, 0)]
    while stack:
        current, depth = stack.pop()
        if depth > maximum_depth:
            raise _fail(f"{label} exceeds JSON depth {maximum_depth}")
        if current is None or type(current) is bool:
            continue
        if type(current) is int:
            if not -SAFE_INTEGER <= current <= SAFE_INTEGER:
                raise _fail(f"{label} contains an unsafe JSON integer")
            continue
        if type(current) is float:
            if not math.isfinite(current):
                raise _fail(f"{label} contains a non-finite number")
            continue
        if type(current) is str:
            if any(0xD800 <= ord(character) <= 0xDFFF for character in current):
                raise _fail(f"{label} contains a lone Unicode surrogate")
            if len(current.encode("utf-8")) > maximum_string_bytes:
                raise _fail(
                    f"{label} exceeds JSON string byte ceiling "
                    f"{maximum_string_bytes}"
                )
            continue
        if type(current) is list:
            members += len(current)
            stack.extend((item, depth + 1) for item in reversed(current))
        elif type(current) is dict:
            members += len(current)
            for key, item in reversed(list(current.items())):
                if type(key) is not str:
                    raise _fail(f"{label} contains a non-string key")
                stack.append((key, depth + 1))
                stack.append((item, depth + 1))
        else:
            raise _fail(f"{label} contains unsupported type {type(current).__name__}")
        if members > maximum_members:
            raise _fail(f"{label} exceeds JSON member ceiling {maximum_members}")


def _validate_wire_size(
    value: dict[str, Any], label: str, maximum_bytes: int
) -> None:
    try:
        encoder = json.JSONEncoder(
            indent=2,
            sort_keys=True,
            ensure_ascii=False,
            allow_nan=False,
        )
        total = 1  # Producers terminate their JSON wire form with one newline.
        for chunk in encoder.iterencode(value):
            total += len(chunk.encode("utf-8"))
            if total > maximum_bytes:
                raise _fail(f"{label} wire form exceeds {maximum_bytes} bytes")
    except BioGraphValidationError:
        raise
    except (TypeError, ValueError, UnicodeError) as failure:
        raise _fail(f"{label} is not serializable I-JSON: {failure}") from failure


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
        ordered = sorted(value, key=lambda key: key.encode("utf-16be"))
        return "{" + ",".join(
            _jcs_string(key) + ":" + _jcs(value[key]) for key in ordered
        ) + "}"
    raise TypeError(f"unsupported JSON type {type(value).__name__}")


def canonical_bytes(value: Any) -> bytes:
    """Return independent RFC 8785 bytes for a bounded I-JSON value."""

    try:
        _validate_tree(value, "canonical value")
        return _jcs(value).encode("utf-8")
    except BioGraphValidationError:
        raise
    except (TypeError, ValueError, UnicodeError) as failure:
        raise _fail(f"value is not RFC 8785 canonical JSON: {failure}") from failure


def digest(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def _seal(core: dict[str, Any]) -> dict[str, Any]:
    return {**core, "artifact_sha256": digest(core)}


def _artifact(value: Any, label: str) -> dict[str, Any]:
    if type(value) is not dict:
        raise _fail(f"{label} must be an object")
    _validate_tree(value, label)
    _validate_wire_size(value, label, MAX_ARTIFACT_BYTES)
    stored = _sha(value.get("artifact_sha256"), f"{label}.artifact_sha256")
    expected = digest({key: item for key, item in value.items() if key != "artifact_sha256"})
    if stored != expected:
        raise _fail(f"{label}.artifact_sha256 does not match canonical content")
    return value


def _backend_spec(value: Any) -> dict[str, Any]:
    """Require the published normative mapping committed by this validator."""

    item = _keys(
        _artifact(value, "backend specification"),
        {
            "format",
            "version",
            "id",
            "backend_version",
            "input",
            "semantics_sidecar",
            "dictionaries",
            "tensors",
            "target_contract",
            "operations",
            "development_abi",
            "artifact_profile",
            "capabilities",
            "artifact_sha256",
        },
        "backend specification",
    )
    if (
        item["format"],
        item["version"],
        item["id"],
        item["backend_version"],
    ) != (
        "brainc.bio.feature-graph-backend-spec",
        1,
        BACKEND_ID,
        BACKEND_VERSION,
    ):
        raise _fail("unsupported backend specification identity")
    if item["artifact_sha256"] != BACKEND_SPEC_SHA256:
        raise _fail(
            "backend specification does not match the independent normative commitment"
        )
    return item


def _no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise _fail(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def _read_regular(path: str | Path, label: str, maximum_bytes: int) -> bytes:
    source = Path(path)
    descriptor = -1
    try:
        inspected = source.lstat()
        if stat.S_ISLNK(inspected.st_mode) or not stat.S_ISREG(inspected.st_mode):
            raise _fail(f"{label} path must be a regular non-linked file")
        if inspected.st_size > maximum_bytes:
            raise _fail(f"{label} path exceeds {maximum_bytes} bytes")
        flags = os.O_RDONLY
        flags |= getattr(os, "O_CLOEXEC", 0)
        flags |= getattr(os, "O_NOFOLLOW", 0)
        flags |= getattr(os, "O_NONBLOCK", 0)
        descriptor = os.open(source, flags)
        opened = os.fstat(descriptor)
        if not stat.S_ISREG(opened.st_mode):
            raise _fail(f"{label} path must be a regular non-linked file")
        if (opened.st_dev, opened.st_ino) != (inspected.st_dev, inspected.st_ino):
            raise _fail(f"{label} path changed while being opened")
        if opened.st_size > maximum_bytes:
            raise _fail(f"{label} path exceeds {maximum_bytes} bytes")
        chunks: list[bytes] = []
        total = 0
        while True:
            request = min(READ_CHUNK_BYTES, maximum_bytes + 1 - total)
            chunk = os.read(descriptor, request)
            if not chunk:
                break
            chunks.append(chunk)
            total += len(chunk)
            if total > maximum_bytes:
                raise _fail(f"{label} path exceeds {maximum_bytes} bytes")
        finished = os.fstat(descriptor)
        stable = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
        if any(getattr(opened, field) != getattr(finished, field) for field in stable):
            raise _fail(f"{label} path changed while being read")
        if total != opened.st_size:
            raise _fail(f"{label} path did not match its inspected length")
        return b"".join(chunks)
    except BioGraphValidationError:
        raise
    except OSError as failure:
        raise _fail(f"cannot read {label}: {failure}") from failure
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _as_object(
    value: dict[str, Any] | str | Path,
    label: str,
    *,
    maximum_bytes: int = MAX_ARTIFACT_BYTES,
) -> dict[str, Any]:
    if type(value) is dict:
        return value
    if not isinstance(value, (str, Path)):
        raise _fail(f"{label} must be an object or JSON path")
    try:
        raw = _read_regular(value, label, maximum_bytes)
        payload = json.loads(
            raw.decode("utf-8"),
            object_pairs_hook=_no_duplicates,
            parse_constant=lambda token: (_ for _ in ()).throw(
                _fail(f"non-finite JSON number {token!r}")
            ),
        )
    except BioGraphValidationError:
        raise
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, ValueError) as failure:
        raise _fail(f"cannot load {label}: {failure}") from failure
    if type(payload) is not dict:
        raise _fail(f"{label} JSON must be an object")
    return payload


def _refget_id(sequence: str) -> str:
    raw = hashlib.sha512(sequence.upper().encode("ascii")).digest()[:24]
    return "SQ." + base64.urlsafe_b64encode(raw).decode("ascii").rstrip("=")


def _sha512t24u(value: Any) -> str:
    raw = hashlib.sha512(canonical_bytes(value)).digest()[:24]
    return base64.urlsafe_b64encode(raw).decode("ascii")


def _source_identifier(value: Any, label: str) -> str:
    identifier = _text(value, label)
    if identifier != identifier.strip() or any(
        character.isspace() for character in identifier
    ):
        raise _fail(f"{label} must be trimmed and whitespace-free")
    if len(identifier.encode("utf-8")) > MAX_SOURCE_IDENTIFIER_BYTES:
        raise _fail(
            f"{label} exceeds {MAX_SOURCE_IDENTIFIER_BYTES} UTF-8 bytes"
        )
    if any(
        ord(character) < 0x20 or 0x7F <= ord(character) <= 0x9F
        for character in identifier
    ):
        raise _fail(f"{label} contains a control character")
    return identifier


def _validate_reference(value: Any, length: int, label: str) -> None:
    if value is None:
        return
    item = _keys(
        value,
        {
            "assembly",
            "contig",
            "start",
            "end",
            "coordinate_system",
            "orientation",
            "aliases",
        },
        label,
    )
    _text(item["assembly"], f"{label}.assembly")
    _text(item["contig"], f"{label}.contig")
    start = _integer(item["start"], f"{label}.start")
    end = _integer(item["end"], f"{label}.end", minimum=1)
    if end <= start or end - start != length:
        raise _fail(f"{label} span must equal sequence length")
    if item["coordinate_system"] != "0-based-half-open":
        raise _fail(f"{label} has an unsupported coordinate system")
    if item["orientation"] not in {"forward", "reverse"}:
        raise _fail(f"{label} has an unsupported orientation")
    aliases = item["aliases"]
    if type(aliases) is not list:
        raise _fail(f"{label}.aliases must be an array")
    checked = [_text(alias, f"{label}.aliases") for alias in aliases]
    if len(checked) != len(set(checked)):
        raise _fail(f"{label}.aliases contains duplicates")


def _validate_provenance(value: Any, label: str) -> None:
    if type(value) is not list:
        raise _fail(f"{label} must be an array")
    identifiers: set[str] = set()
    for index, raw in enumerate(value):
        item_label = f"{label}[{index}]"
        if type(raw) is not dict:
            raise _fail(f"{item_label} must be an object")
        expected = {"id", "kind", "uri", "version"}
        if "sha256" in raw:
            expected.add("sha256")
        item = _keys(raw, expected, item_label)
        identifier = _text(item["id"], f"{item_label}.id")
        if identifier in identifiers:
            raise _fail(f"{label} contains duplicate ids")
        identifiers.add(identifier)
        for field in ("kind", "uri", "version"):
            _text(item[field], f"{item_label}.{field}")
        if "sha256" in item:
            _sha(item["sha256"], f"{item_label}.sha256")


def _validate_segments(
    value: Any,
    sequence_length: int,
    label: str,
    *,
    minimum_line: int = 1,
) -> list[tuple[int, int, int]]:
    body = _keys(value, {"sequence_segments"}, label)
    segments = body["sequence_segments"]
    if type(segments) is not list or not segments:
        raise _fail(f"{label}.sequence_segments must be a nonempty array")
    offset = 0
    last_line = minimum_line - 1
    result: list[tuple[int, int, int]] = []
    for index, raw in enumerate(segments):
        item = _keys(
            raw,
            {"line", "normalized_start", "length"},
            f"{label}.sequence_segments[{index}]",
        )
        line = _integer(
            item["line"],
            f"{label}.line",
            minimum=minimum_line,
            maximum=MAX_SOURCE_LINE,
        )
        start = _integer(item["normalized_start"], f"{label}.normalized_start")
        length = _integer(item["length"], f"{label}.length", minimum=1)
        if line <= last_line or start != offset:
            raise _fail(f"{label}.sequence_segments are discontinuous")
        offset += length
        last_line = line
        result.append((line, start, length))
    if offset != sequence_length:
        raise _fail(f"{label} does not cover the sequence")
    return result


def _validate_sequence_artifact(value: Any, label: str) -> dict[str, Any]:
    item = _keys(
        _artifact(value, label),
        {
            "format",
            "version",
            "compiler",
            "inputs",
            "sequence_ir",
            "ir_sha256",
            "source_map",
            "artifact_sha256",
        },
        label,
    )
    if (item["format"], item["version"], item["compiler"]) != (
        "brain01.sequence-ir",
        2,
        SEQUENCE_COMPILER,
    ):
        raise _fail(f"{label} has an unsupported compiler contract")
    inputs = _keys(item["inputs"], {"fasta_sha256", "context_sha256"}, f"{label}.inputs")
    _sha(inputs["fasta_sha256"], f"{label}.inputs.fasta_sha256")
    if inputs["context_sha256"] is not None:
        _sha(inputs["context_sha256"], f"{label}.inputs.context_sha256")
    ir = _keys(
        item["sequence_ir"],
        {
            "record_id",
            "description",
            "sequence",
            "sequence_sha256",
            "canonical_sha256",
            "refget_id",
            "reference",
            "provenance",
        },
        f"{label}.sequence_ir",
    )
    record_id = _source_identifier(
        ir["record_id"], f"{label}.sequence_ir.record_id"
    )
    description = _text(
        ir["description"], f"{label}.sequence_ir.description", empty=True
    )
    if description and description.splitlines() != [description]:
        raise _fail(f"{label}.sequence_ir.description must be one text line")
    if any(
        ord(character) < 0x20 or 0x7F <= ord(character) <= 0x9F
        for character in description
    ):
        raise _fail(f"{label}.sequence_ir.description contains a control character")
    sequence = _text(ir["sequence"], f"{label}.sequence_ir.sequence")
    try:
        sequence_bytes = sequence.encode("ascii")
    except UnicodeEncodeError as failure:
        raise _fail(f"{label}.sequence_ir.sequence must be ASCII") from failure
    if any(base.upper() not in IUPAC for base in sequence):
        raise _fail(f"{label}.sequence_ir.sequence is not IUPAC DNA")
    if len(sequence_bytes) > MAX_SOURCE_STRING_BYTES:
        raise _fail(
            f"{label}.sequence_ir.sequence exceeds {MAX_SOURCE_STRING_BYTES} bytes"
        )
    if len(sequence_bytes) > MAX_LOGICAL_DNA_BYTES:
        raise _fail(f"{label}.sequence_ir.sequence exceeds the logical DNA ceiling")
    if ir["sequence_sha256"] != hashlib.sha256(sequence_bytes).hexdigest():
        raise _fail(f"{label}.sequence_ir.sequence_sha256 mismatch")
    if ir["canonical_sha256"] != hashlib.sha256(sequence.upper().encode("ascii")).hexdigest():
        raise _fail(f"{label}.sequence_ir.canonical_sha256 mismatch")
    if ir["refget_id"] != _refget_id(sequence):
        raise _fail(f"{label}.sequence_ir.refget_id mismatch")
    _validate_reference(ir["reference"], len(sequence), f"{label}.sequence_ir.reference")
    _validate_provenance(ir["provenance"], f"{label}.sequence_ir.provenance")
    if inputs["context_sha256"] is None and (
        ir["reference"] is not None or ir["provenance"]
    ):
        raise _fail(f"{label} has context assertions without a context input")
    _validate_segments(
        item["source_map"],
        len(sequence),
        f"{label}.source_map",
        minimum_line=2,
    )
    if item["ir_sha256"] != digest(ir):
        raise _fail(f"{label}.ir_sha256 mismatch")
    return item


def _validate_collection(value: Any) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    _validate_tree(
        value,
        "sequence collection",
        maximum_depth=MAX_JSON_DEPTH,
        maximum_members=MAX_SOURCE_JSON_MEMBERS,
        maximum_string_bytes=MAX_SOURCE_STRING_BYTES,
    )
    item = _keys(
        _artifact(value, "sequence collection"),
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
    if (item["format"], item["version"], item["compiler"]) != (
        "brain01.sequence-collection-ir",
        1,
        COLLECTION_COMPILER,
    ):
        raise _fail("unsupported sequence collection contract")
    inputs = _keys(
        item["inputs"],
        {"kind", "wrapper", "raw_sha256", "logical_sha256", "compressed_sha256"},
        "sequence collection.inputs",
    )
    if inputs["kind"] not in {"raw-iupac", "fasta"}:
        raise _fail("sequence collection input kind is unsupported")
    raw_sha = _sha(inputs["raw_sha256"], "sequence collection.inputs.raw_sha256")
    logical_sha = _sha(inputs["logical_sha256"], "sequence collection.inputs.logical_sha256")
    if inputs["wrapper"] is None:
        if inputs["compressed_sha256"] is not None or raw_sha != logical_sha:
            raise _fail("unwrapped collection digest fields are inconsistent")
    elif inputs["wrapper"] == "gzip":
        if _sha(inputs["compressed_sha256"], "compressed_sha256") != raw_sha:
            raise _fail("compressed collection digest fields are inconsistent")
    else:
        raise _fail("sequence collection wrapper is unsupported")
    members = item["members"]
    if (
        type(members) is not list
        or not members
        or len(members) > MAX_SOURCE_RECORDS
    ):
        raise _fail("sequence collection members are outside limits")
    record_ids: list[str] = []
    bindings: list[dict[str, Any]] = []
    collection_members: list[dict[str, str]] = []
    total_bases = 0
    previous_outer_line = 0
    for index, raw in enumerate(members):
        label = f"sequence collection.members[{index}]"
        member = _keys(
            raw,
            {
                "record_id",
                "ir_sha256",
                "artifact_sha256",
                "input_source_map",
                "sequence_artifact",
            },
            label,
        )
        sequence = _validate_sequence_artifact(member["sequence_artifact"], f"{label}.sequence_artifact")
        ir = sequence["sequence_ir"]
        if (
            sequence["inputs"]["context_sha256"] is not None
            or ir["reference"] is not None
            or ir["provenance"]
        ):
            raise _fail(f"{label} contains context not emitted by collection ingress")
        record_id = _source_identifier(member["record_id"], f"{label}.record_id")
        if record_id != ir["record_id"]:
            raise _fail(f"{label}.record_id differs from embedded Sequence IR")
        if member["ir_sha256"] != sequence["ir_sha256"]:
            raise _fail(f"{label}.ir_sha256 differs from embedded Sequence IR")
        if member["artifact_sha256"] != sequence["artifact_sha256"]:
            raise _fail(f"{label}.artifact_sha256 differs from embedded Sequence IR")
        total_bases += len(ir["sequence"])
        if total_bases > MAX_LOGICAL_DNA_BYTES:
            raise _fail(
                f"sequence collection exceeds {MAX_LOGICAL_DNA_BYTES} logical DNA bytes"
            )
        outer_segments = _validate_segments(
            member["input_source_map"],
            len(ir["sequence"]),
            f"{label}.input_source_map",
        )
        local_segments = [
            (segment["line"], segment["normalized_start"], segment["length"])
            for segment in sequence["source_map"]["sequence_segments"]
        ]
        if inputs["kind"] == "raw-iupac":
            if len(members) != 1 or outer_segments != [(1, 0, len(ir["sequence"]))]:
                raise _fail("raw-IUPAC collection must map one member to source line 1")
            if ir["description"] != "":
                raise _fail("raw-IUPAC collection member cannot have a description")
            expected_fasta = f">{record_id}\n{ir['sequence']}\n".encode("utf-8")
            if sequence["inputs"]["fasta_sha256"] != hashlib.sha256(
                expected_fasta
            ).hexdigest():
                raise _fail("raw-IUPAC synthetic FASTA digest is inconsistent")
            if logical_sha != hashlib.sha256(
                ir["sequence"].encode("ascii")
            ).hexdigest():
                raise _fail("raw-IUPAC logical digest does not match sequence")
        else:
            if len(outer_segments) != len(local_segments) or any(
                outer[1:] != local[1:]
                for outer, local in zip(outer_segments, local_segments)
            ):
                raise _fail("FASTA collection local and outer source maps differ")
            shifts = {
                outer[0] - local[0]
                for outer, local in zip(outer_segments, local_segments)
            }
            if len(shifts) != 1:
                raise _fail("FASTA collection source-map line shifts are inconsistent")
            shift = next(iter(shifts))
            if shift < 0 or (index == 0 and shift != 0):
                raise _fail("FASTA collection source-map line shift is invalid")
            if index and shift + 1 <= previous_outer_line:
                raise _fail("FASTA collection member overlaps the preceding record")
            previous_outer_line = outer_segments[-1][0]
        record_ids.append(record_id)
        collection_members.append({"record_id": record_id, "ir_sha256": sequence["ir_sha256"]})
        bindings.append(
            {
                "seqid": record_id,
                "length": len(ir["sequence"]),
                "sequence_sha256": ir["sequence_sha256"],
                "sequence_artifact_sha256": sequence["artifact_sha256"],
            }
        )
    if len(record_ids) != len(set(record_ids)):
        raise _fail("sequence collection contains duplicate record ids")
    expected_ir = {"members": collection_members}
    if item["collection_ir"] != expected_ir or item["collection_ir_sha256"] != digest(expected_ir):
        raise _fail("sequence collection IR does not match its members")
    level_2 = {
        "lengths": [binding["length"] for binding in bindings],
        "names": record_ids,
        "sequences": [member["sequence_artifact"]["sequence_ir"]["refget_id"] for member in members],
    }
    level_1 = {name: _sha512t24u(values) for name, values in level_2.items()}
    expected_refget = {
        "version": "1.0.0",
        "digest": _sha512t24u({name: level_1[name] for name in ("names", "sequences")}),
        "level_1": level_1,
        "level_2": level_2,
    }
    if item["refget_seqcol"] != expected_refget:
        raise _fail("sequence collection refget identity mismatch")
    return item, sorted(bindings, key=lambda binding: binding["seqid"].encode("utf-8"))


def _attribute(attributes: list[dict[str, Any]], tag: str) -> list[str]:
    for attribute in attributes:
        if attribute["tag"] == tag:
            return attribute["values"]
    return []


def _validate_attributes(value: Any, label: str) -> list[dict[str, Any]]:
    if type(value) is not list or len(value) > 256:
        raise _fail(f"{label} must be a bounded array")
    tags: list[str] = []
    for index, raw in enumerate(value):
        item = _keys(raw, {"tag", "values"}, f"{label}[{index}]")
        tag = _text(item["tag"], f"{label}[{index}].tag")
        if ATTRIBUTE_RE.fullmatch(tag) is None:
            raise _fail(f"{label}[{index}].tag is invalid")
        if tag[0].isupper() and tag not in RESERVED_ATTRIBUTES:
            raise _fail(f"{label}[{index}] uses an unsupported reserved attribute")
        values = item["values"]
        if type(values) is not list or not values or len(values) > 10_000:
            raise _fail(f"{label}[{index}].values must be a bounded nonempty array")
        checked = [_text(entry, f"{label}[{index}].values") for entry in values]
        if tag in CANONICAL_MULTI_VALUE_ATTRIBUTES:
            _utf8_sorted(checked, f"{label}[{index}].values")
        if tag in {"ID", "Is_circular", "Target", "Gap"} and len(checked) != 1:
            raise _fail(f"{label}[{index}] requires exactly one value")
        if tag == "Is_circular" and checked[0] not in {"true", "false"}:
            raise _fail(f"{label}[{index}] has an invalid Is_circular value")
        tags.append(tag)
    _utf8_sorted(tags, f"{label} tags")
    return value


def _normalized_intervals(start: int, end: int, length: int, circular: bool) -> list[dict[str, int]]:
    if start > length:
        raise _fail("feature starts beyond its sequence")
    span = end - start + 1
    if end <= length:
        return [{"start": start - 1, "end": end}]
    if not circular or span > length:
        raise _fail("feature crosses a non-circular boundary or more than one traversal")
    first_start = start - 1
    return [
        {"start": first_start, "end": length},
        {"start": 0, "end": first_start + span - length},
    ]


def _cycle_check(relationships: list[dict[str, str]]) -> None:
    graph: dict[str, list[str]] = {}
    for relationship in relationships:
        graph.setdefault(relationship["source"], []).append(relationship["target"])
    state: dict[str, int] = {}
    for root in sorted(graph, key=lambda item: item.encode("utf-8")):
        if state.get(root) == 2:
            continue
        state[root] = 1
        stack: list[tuple[str, int]] = [(root, 0)]
        while stack:
            node, index = stack[-1]
            targets = graph.get(node, [])
            if index == len(targets):
                state[node] = 2
                stack.pop()
                continue
            target = targets[index]
            stack[-1] = (node, index + 1)
            if state.get(target) == 1:
                raise _fail("BioIR relationship graph contains a cycle")
            if state.get(target) == 2:
                continue
            state[target] = 1
            stack.append((target, 0))


def _validate_bio(value: Any, source: dict[str, Any], records: list[dict[str, Any]]) -> dict[str, Any]:
    item = _keys(
        _artifact(value, "BioIR artifact"),
        {"format", "version", "compiler", "inputs", "bio_ir", "bio_ir_sha256", "artifact_sha256"},
        "BioIR artifact",
    )
    if (item["format"], item["version"], item["compiler"]) != (
        "brainc.bio.gff3-ir",
        1,
        BIO_COMPILER,
    ):
        raise _fail("unsupported BioIR compiler contract")
    inputs = _keys(item["inputs"], {"gff3_sha256", "sequence_collection"}, "BioIR inputs")
    _sha(inputs["gff3_sha256"], "BioIR inputs.gff3_sha256")
    source_binding = _keys(
        inputs["sequence_collection"],
        {"format", "version", "artifact_sha256", "collection_ir_sha256"},
        "BioIR sequence binding",
    )
    expected_source_binding = {
        "format": "brain01.sequence-collection-ir",
        "version": 1,
        "artifact_sha256": source["artifact_sha256"],
        "collection_ir_sha256": source["collection_ir_sha256"],
    }
    if source_binding != expected_source_binding:
        raise _fail("BioIR does not bind the supplied sequence collection")
    ir = _keys(
        item["bio_ir"],
        {
            "gff3_version",
            "coordinate_system",
            "records",
            "sequence_regions",
            "directives",
            "features",
            "relationships",
        },
        "BioIR body",
    )
    if item["bio_ir_sha256"] != digest(ir):
        raise _fail("BioIR body digest mismatch")
    if type(ir["gff3_version"]) is not str or re.fullmatch(r"3(?:\.[0-9]+){0,2}", ir["gff3_version"]) is None:
        raise _fail("BioIR GFF3 version is invalid")
    if ir["coordinate_system"] != {"source": "1-based-closed", "normalized": "0-based-half-open"}:
        raise _fail("BioIR coordinate declaration is invalid")
    if ir["records"] != records:
        raise _fail("BioIR record bindings differ from the supplied sequence collection")
    lengths = {record["seqid"]: record["length"] for record in records}
    regions = ir["sequence_regions"]
    if type(regions) is not list:
        raise _fail("BioIR sequence_regions must be an array")
    region_map: dict[str, tuple[int, int]] = {}
    occupied_lines = {1}
    last_line = 0
    for index, raw in enumerate(regions):
        region = _keys(raw, {"line", "seqid", "raw_start", "raw_end", "start", "end"}, f"sequence_regions[{index}]")
        line = _integer(region["line"], f"sequence_regions[{index}].line", minimum=2, maximum=1_000_000)
        if line <= last_line or line in occupied_lines:
            raise _fail("BioIR sequence regions are not in unique source-line order")
        last_line = line
        occupied_lines.add(line)
        seqid = _text(region["seqid"], f"sequence_regions[{index}].seqid")
        start = _integer(region["raw_start"], f"sequence_regions[{index}].raw_start", minimum=1)
        end = _integer(region["raw_end"], f"sequence_regions[{index}].raw_end", minimum=1)
        if seqid not in lengths or seqid in region_map or end < start or end > lengths[seqid]:
            raise _fail("BioIR sequence region is invalid")
        if (region["start"], region["end"]) != (start - 1, end):
            raise _fail("BioIR sequence region normalization is invalid")
        region_map[seqid] = (start - 1, end)
    directives = ir["directives"]
    if type(directives) is not list or len(directives) > 100_000:
        raise _fail("BioIR directives are outside limits")
    last_line = 0
    for index, raw in enumerate(directives):
        directive = _keys(raw, {"line", "name", "value"}, f"directives[{index}]")
        line = _integer(directive["line"], f"directives[{index}].line", minimum=2, maximum=1_000_000)
        if line <= last_line or line in occupied_lines:
            raise _fail("BioIR directives are not in unique source-line order")
        last_line = line
        occupied_lines.add(line)
        name = _text(directive["name"], f"directives[{index}].name")
        if DIRECTIVE_RE.fullmatch(name) is None or name in {"gff-version", "sequence-region"}:
            raise _fail("BioIR opaque directive name is invalid")
        _text(directive["value"], f"directives[{index}].value", empty=True)
    features = ir["features"]
    if type(features) is not list or len(features) > MAX_FEATURE_ROWS:
        raise _fail("BioIR features are outside limits")
    feature_ids: list[str] = []
    feature_map: dict[str, dict[str, Any]] = {}
    feature_lines: set[int] = set()
    total_segments = 0
    for feature_index, raw in enumerate(features):
        label = f"features[{feature_index}]"
        feature = _keys(raw, {"entity_id", "declared_id", "seqid", "source", "type", "strand", "segments"}, label)
        entity_id = _text(feature["entity_id"], f"{label}.entity_id")
        declared_id = feature["declared_id"]
        if declared_id is None:
            if not entity_id.startswith("anon:"):
                raise _fail(f"{label} anonymous identity is invalid")
        elif entity_id != "id:" + _text(declared_id, f"{label}.declared_id"):
            raise _fail(f"{label} declared identity is inconsistent")
        if entity_id in feature_map:
            raise _fail("BioIR contains duplicate feature identities")
        feature_ids.append(entity_id)
        feature_map[entity_id] = feature
        seqid = _text(feature["seqid"], f"{label}.seqid")
        if seqid not in lengths:
            raise _fail(f"{label} references an unknown sequence")
        if feature["source"] is not None:
            _text(feature["source"], f"{label}.source")
        feature_type = _text(feature["type"], f"{label}.type")
        if feature["strand"] not in STRANDS:
            raise _fail(f"{label}.strand is invalid")
        segments = feature["segments"]
        if type(segments) is not list or not segments or (declared_id is None and len(segments) != 1):
            raise _fail(f"{label}.segments are invalid")
        segment_lines: list[int] = []
        for segment_index, raw_segment in enumerate(segments):
            total_segments += 1
            if total_segments > MAX_FEATURE_ROWS:
                raise _fail(
                    f"BioIR exceeds {MAX_FEATURE_ROWS} total feature segments"
                )
            segment_label = f"{label}.segments[{segment_index}]"
            segment = _keys(
                raw_segment,
                {"line", "raw_start", "raw_end", "intervals", "score", "phase", "attributes"},
                segment_label,
            )
            line = _integer(segment["line"], f"{segment_label}.line", minimum=2, maximum=1_000_000)
            if line in occupied_lines or line in feature_lines:
                raise _fail("BioIR feature source lines are not unique")
            feature_lines.add(line)
            segment_lines.append(line)
            start = _integer(segment["raw_start"], f"{segment_label}.raw_start", minimum=1)
            end = _integer(segment["raw_end"], f"{segment_label}.raw_end", minimum=1)
            if end < start:
                raise _fail(f"{segment_label}.raw_end precedes raw_start")
            score = segment["score"]
            if score is not None:
                if type(score) is not str or SCORE_RE.fullmatch(score) is None:
                    raise _fail(f"{segment_label}.score is invalid")
                try:
                    if not Decimal(score).is_finite():
                        raise _fail(f"{segment_label}.score is non-finite")
                except InvalidOperation as failure:
                    raise _fail(f"{segment_label}.score is invalid") from failure
            if feature_type == "CDS":
                if segment["phase"] not in {0, 1, 2} or type(segment["phase"]) is not int:
                    raise _fail(f"{segment_label}.phase is invalid for CDS")
            elif segment["phase"] is not None:
                raise _fail(f"{segment_label}.phase must be null outside CDS")
            attributes = _validate_attributes(segment["attributes"], f"{segment_label}.attributes")
            ids = _attribute(attributes, "ID")
            if (ids[0] if ids else None) != declared_id:
                raise _fail(f"{segment_label} ID attribute differs from feature identity")
        if segment_lines != sorted(segment_lines):
            raise _fail(f"{label}.segments are not in source-line order")
    _utf8_sorted(feature_ids, "BioIR feature identities")
    circular = {
        feature["seqid"]
        for feature in features
        for segment in feature["segments"]
        if _attribute(segment["attributes"], "Is_circular") == ["true"]
        and segment["raw_start"] == 1
        and segment["raw_end"] == lengths[feature["seqid"]]
    }
    for feature in features:
        for segment in feature["segments"]:
            expected = _normalized_intervals(
                segment["raw_start"], segment["raw_end"], lengths[feature["seqid"]], feature["seqid"] in circular
            )
            if segment["intervals"] != expected:
                raise _fail("BioIR normalized intervals differ from source coordinates")
            region = region_map.get(feature["seqid"])
            if region is not None and any(interval["start"] < region[0] or interval["end"] > region[1] for interval in expected):
                raise _fail("BioIR feature lies outside its sequence region")
            if feature["declared_id"] is None:
                signature = {
                    "seqid": feature["seqid"],
                    "source": feature["source"],
                    "type": feature["type"],
                    "strand": feature["strand"],
                    "segment": segment,
                }
                if feature["entity_id"] != "anon:" + digest(signature):
                    raise _fail("BioIR anonymous identity is not deterministic")
    relationships = ir["relationships"]
    if type(relationships) is not list or len(relationships) > MAX_RELATIONSHIPS:
        raise _fail("BioIR relationships are outside limits")
    actual: list[tuple[str, str, str]] = []
    for index, raw in enumerate(relationships):
        item_label = f"relationships[{index}]"
        relationship = _keys(raw, {"source", "kind", "target"}, item_label)
        source_id = _text(relationship["source"], f"{item_label}.source")
        target_id = _text(relationship["target"], f"{item_label}.target")
        if relationship["kind"] not in RELATIONSHIP_KINDS:
            raise _fail(f"{item_label}.kind is unsupported")
        if source_id not in feature_map or target_id not in feature_map:
            raise _fail(f"{item_label} contains a dangling endpoint")
        actual.append((source_id, relationship["kind"], target_id))
    wanted_order = sorted(actual, key=lambda row: tuple(item.encode("utf-8") for item in row))
    if actual != wanted_order or len(actual) != len(set(actual)):
        raise _fail("BioIR relationships are not sorted and unique")
    expected_relationships = {
        (feature["entity_id"], kind, "id:" + target)
        for feature in features
        for segment in feature["segments"]
        for tag, kind in (("Parent", "parent"), ("Derives_from", "derives-from"))
        for target in _attribute(segment["attributes"], tag)
    }
    if actual != sorted(expected_relationships, key=lambda row: tuple(item.encode("utf-8") for item in row)):
        raise _fail("BioIR relationships do not exactly match feature attributes")
    _cycle_check(relationships)
    return item


def _source_binding(source: dict[str, Any]) -> dict[str, Any]:
    return {
        "format": source["format"],
        "version": source["version"],
        "artifact_sha256": source["artifact_sha256"],
        "ir_sha256": source["collection_ir_sha256"],
    }


def _bio_binding(bio: dict[str, Any]) -> dict[str, Any]:
    return {
        "format": bio["format"],
        "version": bio["version"],
        "artifact_sha256": bio["artifact_sha256"],
        "bio_ir_sha256": bio["bio_ir_sha256"],
    }


def _indexed(values: list[Any]) -> list[dict[str, Any]]:
    return [{"code": index, "value": value} for index, value in enumerate(values)]


def _expected_semantics(source: dict[str, Any], bio: dict[str, Any]) -> dict[str, Any]:
    ir = bio["bio_ir"]
    features = ir["features"]
    seqids = sorted({feature["seqid"] for feature in features}, key=lambda value: value.encode("utf-8"))
    sources = [None] + sorted(
        {feature["source"] for feature in features if feature["source"] is not None},
        key=lambda value: value.encode("utf-8"),
    )
    types = sorted({feature["type"] for feature in features}, key=lambda value: value.encode("utf-8"))
    feature_index = {feature["entity_id"]: index for index, feature in enumerate(features)}
    edge_order: list[dict[str, Any]] = []
    for kind in RELATIONSHIP_KINDS:
        for relationship_index, relationship in enumerate(ir["relationships"]):
            if relationship["kind"] == kind:
                edge_order.append(
                    {
                        "edge_index": len(edge_order),
                        "bio_ir_relationship_index": relationship_index,
                        "source_unit_index": feature_index[relationship["source"]],
                        "target_unit_index": feature_index[relationship["target"]],
                        "source": relationship["source"],
                        "kind": kind,
                        "target": relationship["target"],
                    }
                )
    semantics = {
        "dictionaries": {
            "relationship_kinds": _indexed(list(RELATIONSHIP_KINDS)),
            "seqids": _indexed(seqids),
            "sources": _indexed(sources),
            "strands": _indexed(list(STRANDS)),
            "types": _indexed(types),
        },
        "unit_order": [feature["entity_id"] for feature in features],
        "edge_order": edge_order,
        "bio_ir": ir,
    }
    backend_spec_sha = BACKEND_SPEC_SHA256
    return _seal(
        {
            "format": SEMANTICS_FORMAT,
            "version": 1,
            "backend": {"id": BACKEND_ID, "version": 1, "spec_sha256": backend_spec_sha},
            "inputs": {"sequence": _source_binding(source), "bio_ir": _bio_binding(bio)},
            "semantics": semantics,
            "semantics_sha256": digest(semantics),
        }
    )


def _dictionary_codes(semantics: dict[str, Any], name: str) -> dict[Any, int]:
    return {entry["value"]: entry["code"] for entry in semantics["semantics"]["dictionaries"][name]}


def _feature_values(semantics: dict[str, Any]) -> dict[str, list[bool | int]]:
    features = semantics["semantics"]["bio_ir"]["features"]
    seqids = _dictionary_codes(semantics, "seqids")
    sources = _dictionary_codes(semantics, "sources")
    strands = _dictionary_codes(semantics, "strands")
    types = _dictionary_codes(semantics, "types")
    result: dict[str, list[bool | int]] = {name: [] for name, _ in FEATURE_FIELDS}
    for feature in features:
        intervals = [interval for segment in feature["segments"] for interval in segment["intervals"]]
        phases = {segment["phase"] for segment in feature["segments"] if segment["phase"] is not None}
        numeric = {
            "attribute_count": sum(len(segment["attributes"]) for segment in feature["segments"]),
            "interval_count": len(intervals),
            "location_max_end": max(interval["end"] for interval in intervals),
            "location_min_start": min(interval["start"] for interval in intervals),
            "phase_mask": sum(1 << phase for phase in phases),
            "segment_bases": sum(interval["end"] - interval["start"] for interval in intervals),
            "segment_count": len(feature["segments"]),
            "seqid_code": seqids[feature["seqid"]],
            "source_code": sources[feature["source"]],
            "strand_code": strands[feature["strand"]],
            "type_code": types[feature["type"]],
        }
        for name, value in numeric.items():
            if type(value) is not int or not 0 <= value <= MAX_U64:
                raise _fail(f"feature field {name} exceeds u64")
            result[name].append(value)
        result["declared_id"].append(feature["declared_id"] is not None)
        result["wraps_origin"].append(any(len(segment["intervals"]) > 1 for segment in feature["segments"]))
    return result


def _pack(dtype: str, values: list[bool | int]) -> bytes:
    if dtype == "bool":
        if any(type(value) is not bool for value in values):
            raise _fail("boolean tensor contains a non-boolean value")
        return bytes(1 if value else 0 for value in values)
    if dtype == "u64":
        if any(type(value) is not int or not 0 <= value <= MAX_U64 for value in values):
            raise _fail("u64 tensor contains an out-of-range value")
        return struct.pack("<" + "Q" * len(values), *values)
    raise _fail(f"unsupported bridge tensor dtype {dtype!r}")


def _storage(raw: bytes) -> dict[str, Any]:
    return {
        "kind": "inline-base64",
        "data": base64.b64encode(raw).decode("ascii"),
        "byte_length": len(raw),
        "sha256": hashlib.sha256(raw).hexdigest(),
    }


def _tensor(identifier: str, dtype: str, shape: list[int], values: list[bool | int]) -> dict[str, Any]:
    return {
        "id": identifier,
        "type": {"dtype": dtype, "shape": shape},
        "unit": None,
        "axes": [None] * len(shape),
        "storage": _storage(_pack(dtype, values)),
    }


def _outputs(semantics: dict[str, Any]) -> list[dict[str, Any]]:
    features = semantics["semantics"]["bio_ir"]["features"]
    values = _feature_values(semantics)
    outputs = [_tensor("t.feature.count", "u64", [], [len(features)])]
    for name, dtype in FEATURE_FIELDS:
        outputs.append(_tensor(f"t.feature.{name}", dtype, [len(features)], values[name]))
    kind_codes = _dictionary_codes(semantics, "relationship_kinds")
    edge_order = semantics["semantics"]["edge_order"]
    for kind in RELATIONSHIP_KINDS:
        edges = [edge for edge in edge_order if edge["kind"] == kind]
        pairs = [endpoint for edge in edges for endpoint in (edge["source_unit_index"], edge["target_unit_index"])]
        stem = f"t.relationship.{kind}"
        outputs.append(_tensor(f"{stem}.kind_code", "u64", [], [kind_codes[kind]]))
        outputs.append(_tensor(f"{stem}.pairs", "u64", [len(edges), 2], pairs))
    return sorted(outputs, key=lambda output: output["id"])


def _field(identifier: str, dtype: str) -> dict[str, Any]:
    return {
        "id": identifier,
        "type": {"dtype": dtype, "shape": []},
        "unit": None,
        "mutability": "constant",
        "numeric": {"kind": "boolean" if dtype == "bool" else "integer"},
    }


def _target() -> dict[str, Any]:
    contract = {
        "id": BACKEND_ID,
        "abi_major": 1,
        "opsets": [{"domain": DEV_DOMAIN, "version": 1}],
        "unit_schemas": [{"id": "feature", "fields": [_field(name, dtype) for name, dtype in FEATURE_FIELDS]}],
        "edge_schemas": [
            {"id": f"relationship.{kind}", "fields": [_field("kind_code", "u64")]}
            for kind in RELATIONSHIP_KINDS
        ],
        "ports": [],
        "rules": [],
    }
    return _seal(
        {
            "format": "brainc.target-contract",
            "version": 1,
            "contract": contract,
            "contract_sha256": digest(contract),
        }
    )


def _operations() -> list[dict[str, Any]]:
    operations: list[dict[str, Any]] = [
        {
            "id": "create.features",
            "op": OP_UNIT_CREATE,
            "version": 1,
            "schema": "feature",
            "count": "t.feature.count",
            "initializers": [
                {"field": name, "tensor": f"t.feature.{name}"} for name, _ in FEATURE_FIELDS
            ],
        }
    ]
    for kind in RELATIONSHIP_KINDS:
        stem = f"relationship.{kind}"
        operations.append(
            {
                "id": f"create.{stem}",
                "op": OP_EDGE_CREATE,
                "version": 1,
                "schema": stem,
                "sources": "create.features",
                "targets": "create.features",
                "pairs": f"t.{stem}.pairs",
                "initializers": [{"field": "kind_code", "tensor": f"t.{stem}.kind_code"}],
            }
        )
    return operations


def _contract(output: dict[str, Any]) -> dict[str, Any]:
    return {key: output[key] for key in ("id", "type", "unit", "axes")}


def _expected_chain(
    source: dict[str, Any],
    bio: dict[str, Any],
    backend_spec: dict[str, Any],
) -> dict[str, dict[str, Any]]:
    semantics = _expected_semantics(source, bio)
    outputs = _outputs(semantics)
    backend_spec_sha = BACKEND_SPEC_SHA256
    manifest = _seal(
        {
            "format": "brainc.provider-manifest",
            "version": 2,
            "provider": {"name": BACKEND_ID, "version": "1"},
            "model_identity": {"kind": "content-sha256", "value": backend_spec_sha},
            "accepts": ["brain01.sequence-collection-ir/v1"],
            "outputs": [_contract(output) for output in outputs],
        }
    )
    request = _seal(
        {
            "format": "brainc.prediction-request",
            "version": 2,
            "source": _source_binding(source),
            "provider_manifest_sha256": manifest["artifact_sha256"],
            "requested_outputs": [_contract(output) for output in outputs],
        }
    )
    response = _seal(
        {
            "format": "brainc.prediction-response",
            "version": 2,
            "request_artifact_sha256": request["artifact_sha256"],
            "provider": manifest["provider"],
            "model_identity": manifest["model_identity"],
            "outputs": outputs,
        }
    )
    target = _target()
    target_binding = {
        "id": BACKEND_ID,
        "abi_major": 1,
        "contract_sha256": target["contract_sha256"],
    }
    operations = _operations()
    policy = _seal(
        {
            "format": "brainc.lowering-policy",
            "version": 2,
            "id": f"{BACKEND_ID}.lowering",
            "target": target_binding,
            "values": [{"id": output["id"], "from_output": output["id"]} for output in outputs],
            "operations": operations,
        }
    )
    emitted_tensors = [
        {
            **output,
            "lineage": {"response_output": output["id"], "lowering_value": output["id"]},
        }
        for output in outputs
    ]
    budgets = {
        "operations": len(operations),
        "tensor_bytes": sum(output["storage"]["byte_length"] for output in outputs),
        "units": len(bio["bio_ir"]["features"]),
        "edges": len(bio["bio_ir"]["relationships"]),
        "attachments": 0,
    }
    body = {
        "target": target_binding,
        "requirements": [{"domain": DEV_DOMAIN, "version": 1}],
        "budgets": budgets,
        "tensors": emitted_tensors,
        "entrypoint": {"name": "develop", "operations": operations},
    }
    module = _seal(
        {
            "format": "brainc.development-module",
            "version": 1,
            "producer": MODULE_PRODUCER,
            "sources": {
                "sequence": _source_binding(source),
                "provider_manifest": {"artifact_sha256": manifest["artifact_sha256"]},
                "prediction_request": {"artifact_sha256": request["artifact_sha256"]},
                "prediction_response": {"artifact_sha256": response["artifact_sha256"]},
                "lowering_policy": {"artifact_sha256": policy["artifact_sha256"]},
                "target_contract": {"artifact_sha256": target["artifact_sha256"]},
            },
            "module": body,
            "module_sha256": digest(body),
        }
    )
    record = _seal(
        {
            "format": RECORD_FORMAT,
            "version": 1,
            "producer": BIO_GRAPH_PRODUCER,
            "backend": {
                "id": BACKEND_ID,
                "version": 1,
                "spec_sha256": backend_spec_sha,
                "semantics_sha256": semantics["semantics_sha256"],
            },
            "inputs": {"sequence": _source_binding(source), "bio_ir": _bio_binding(bio)},
            "artifacts": {
                "backend_spec": {
                    "artifact_sha256": backend_spec["artifact_sha256"]
                },
                "backend_semantics": {"artifact_sha256": semantics["artifact_sha256"]},
                "development_module": {"artifact_sha256": module["artifact_sha256"]},
                "lowering_policy": {"artifact_sha256": policy["artifact_sha256"]},
                "prediction_request": {"artifact_sha256": request["artifact_sha256"]},
                "prediction_response": {"artifact_sha256": response["artifact_sha256"]},
                "provider_manifest": {"artifact_sha256": manifest["artifact_sha256"]},
                "target_contract": {"artifact_sha256": target["artifact_sha256"]},
            },
            "result": {
                "module_sha256": module["module_sha256"],
                "operations": budgets["operations"],
                "tensor_bytes": budgets["tensor_bytes"],
                "units": budgets["units"],
                "edges": budgets["edges"],
                "attachments": 0,
                "ports": 0,
                "inference": False,
                "learning": False,
            },
        }
    )
    return {
        "backend_spec": backend_spec,
        "backend_semantics": semantics,
        "compilation_record": record,
        "development_module": module,
        "lowering_policy": policy,
        "prediction_request": request,
        "prediction_response": response,
        "provider_manifest": manifest,
        "target_contract": target,
    }


def _first_difference(actual: Any, expected: Any, path: str = "$") -> str:
    if type(actual) is not type(expected):
        return f"{path} type {type(actual).__name__} != {type(expected).__name__}"
    if type(actual) is dict:
        actual_keys, expected_keys = set(actual), set(expected)
        if actual_keys != expected_keys:
            return f"{path} keys differ; missing={sorted(expected_keys-actual_keys)}, unknown={sorted(actual_keys-expected_keys)}"
        for key in sorted(actual, key=lambda item: item.encode("utf-8")):
            if actual[key] != expected[key]:
                return _first_difference(actual[key], expected[key], f"{path}.{key}")
    elif type(actual) is list:
        if len(actual) != len(expected):
            return f"{path} length {len(actual)} != {len(expected)}"
        for index, (left, right) in enumerate(zip(actual, expected)):
            if left != right:
                return _first_difference(left, right, f"{path}[{index}]")
    elif actual != expected:
        return f"{path} value differs"
    return path


def validate_feature_graph_bundle(
    bundle: dict[str, Any] | str | Path,
    sequence_collection: dict[str, Any] | str | Path,
    bio_ir: dict[str, Any] | str | Path,
) -> dict[str, Any]:
    """Independently replay and validate one complete BioIR bridge bundle.

    The supplied sequence and BioIR are external validation inputs, not values
    trusted from inside the bundle.  A successful call returns a sealed report;
    any contract or replay mismatch raises :class:`BioGraphValidationError`.
    """

    source_value, records = _validate_collection(
        _as_object(
            sequence_collection,
            "sequence collection",
            maximum_bytes=MAX_SOURCE_JSON_BYTES,
        )
    )
    bio_value = _validate_bio(
        _as_object(
            bio_ir,
            "BioIR artifact",
            maximum_bytes=MAX_ARTIFACT_BYTES,
        ),
        source_value,
        records,
    )
    supplied = _keys(
        _artifact(
            _as_object(
                bundle,
                "feature-graph bundle",
                maximum_bytes=MAX_ARTIFACT_BYTES,
            ),
            "feature-graph bundle",
        ),
        {"format", "version", "artifacts", "artifact_sha256"},
        "feature-graph bundle",
    )
    if (supplied["format"], supplied["version"]) != (BUNDLE_FORMAT, 1):
        raise _fail("unsupported feature-graph bundle format or version")
    supplied_artifacts = _keys(
        supplied["artifacts"],
        {
            "backend_spec",
            "backend_semantics",
            "compilation_record",
            "development_module",
            "lowering_policy",
            "prediction_request",
            "prediction_response",
            "provider_manifest",
            "target_contract",
        },
        "feature-graph bundle.artifacts",
    )
    for role, artifact in supplied_artifacts.items():
        _artifact(artifact, f"feature-graph bundle.artifacts.{role}")
    backend_spec = _backend_spec(supplied_artifacts["backend_spec"])
    expected_artifacts = _expected_chain(source_value, bio_value, backend_spec)
    expected_bundle = _seal(
        {"format": BUNDLE_FORMAT, "version": 1, "artifacts": expected_artifacts}
    )
    if supplied != expected_bundle:
        raise _fail(
            "bundle differs from independent BioIR-to-module replay at "
            + _first_difference(supplied, expected_bundle)
        )
    budgets = expected_artifacts["development_module"]["module"]["budgets"]
    core = {
        "format": "brainc.bio.feature-graph-validation-report",
        "version": 1,
        "stage": "bio-ir-feature-graph-bridge",
        "valid": True,
        "inputs": {
            "sequence_artifact_sha256": source_value["artifact_sha256"],
            "bio_ir_artifact_sha256": bio_value["artifact_sha256"],
            "bundle_artifact_sha256": supplied["artifact_sha256"],
            "backend_spec_sha256": backend_spec["artifact_sha256"],
        },
        "result": {
            **budgets,
            "relationship_edges": len(bio_value["bio_ir"]["relationships"]),
            "semantics_sha256": expected_artifacts["backend_semantics"]["semantics_sha256"],
            "module_sha256": expected_artifacts["development_module"]["module_sha256"],
        },
        "checks": [
            "sequence-collection-identity",
            "bio-ir-structure-and-binding",
            "published-normative-backend-specification",
            "complete-normalized-bio-ir-sidecar-and-dictionaries",
            "numeric-tensor-encoding",
            "explicit-relationship-cardinality",
            "target-policy-provider-chain",
            "development-module-linkage-and-budgets",
            "compilation-record-and-bundle-bindings",
        ],
    }
    return {**core, "report_sha256": digest(core)}


__all__ = [
    "BioGraphValidationError",
    "canonical_bytes",
    "digest",
    "validate_feature_graph_bundle",
]

"""Deterministic compilation of one ordinary FASTA record into Sequence IR.

Biological assertions are inputs, not compiler constants.  FASTA carries the
sequence and record label; an optional context document carries reference and
provenance assertions whose origin remains visible in the emitted artifact.
"""

from __future__ import annotations

from dataclasses import dataclass
import base64
import hashlib
from pathlib import Path
import re
from typing import Any

from ._canonical import ContractError, canonical_bytes
from ._io import (
    BoundedIOError,
    MAX_IDENTIFIER_BYTES,
    MAX_INPUT_BYTES,
    MAX_JSON_BYTES,
    MAX_JSON_DEPTH,
    MAX_JSON_MEMBERS,
    MAX_STRING_BYTES,
    atomic_write_file,
    load_json_object,
    pretty_json_bytes,
    read_regular_file,
)


SEQUENCE_COMPILER_VERSION = "0.3.0"
SEQUENCE_PASSES = (
    "parse-fasta",
    "parse-context",
    "validate-iupac",
    "resolve-reference",
    "canonicalize-sequence",
    "emit-sequence-ir",
)
IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}\Z")
COORDINATE_SYSTEM = "0-based-half-open"
MAX_SOURCE_LINE_BREAKS = MAX_JSON_MEMBERS


class SequenceCompilerError(ValueError):
    """A deterministic sequence compilation failure."""


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _canonical_bytes(value: Any) -> bytes:
    try:
        return canonical_bytes(value)
    except ContractError as failure:
        raise SequenceCompilerError(f"DNA037: {failure}") from failure


def _canonical_digest(value: Any) -> str:
    return _sha256_bytes(_canonical_bytes(value))


def _load_json_bytes(raw: bytes, label: str) -> dict[str, Any]:
    try:
        return load_json_object(
            raw,
            label,
            maximum_bytes=MAX_JSON_BYTES,
            maximum_depth=MAX_JSON_DEPTH,
            maximum_members=MAX_JSON_MEMBERS,
            maximum_string_bytes=MAX_STRING_BYTES,
        )
    except BoundedIOError as failure:
        raise SequenceCompilerError(f"DNA020: invalid {label} JSON: {failure}") from failure


def _artifact_bytes(payload: dict[str, Any]) -> bytes:
    try:
        return pretty_json_bytes(
            payload,
            ensure_ascii=True,
            maximum_bytes=MAX_JSON_BYTES,
            maximum_depth=MAX_JSON_DEPTH,
            maximum_members=MAX_JSON_MEMBERS,
            maximum_string_bytes=MAX_STRING_BYTES,
        )
    except BoundedIOError as failure:
        raise SequenceCompilerError(f"DNA038: sequence artifact is outside compiler limits: {failure}") from failure


def _require_keys(value: dict[str, Any], required: set[str], optional: set[str], label: str) -> None:
    missing = sorted(required - value.keys())
    extra = sorted(value.keys() - required - optional)
    if missing:
        raise SequenceCompilerError(f"DNA021: {label} is missing keys: {', '.join(missing)}")
    if extra:
        raise SequenceCompilerError(f"DNA022: {label} has unknown keys: {', '.join(extra)}")


def _require_text(value: Any, label: str) -> str:
    if type(value) is not str or not value or value != value.strip():
        raise SequenceCompilerError(f"DNA023: {label} must be a non-empty trimmed string")
    return value


def _require_record_id(value: Any, label: str) -> str:
    record_id = _require_text(value, label)
    if any(character.isspace() for character in record_id):
        raise SequenceCompilerError(f"DNA023: {label} must be whitespace-free")
    if any(
        ord(character) < 0x20 or 0x7F <= ord(character) <= 0x9F
        for character in record_id
    ):
        raise SequenceCompilerError(f"DNA023: {label} contains a control character")
    if len(record_id.encode("utf-8")) > MAX_IDENTIFIER_BYTES:
        raise SequenceCompilerError(
            f"DNA023: {label} exceeds {MAX_IDENTIFIER_BYTES} UTF-8 bytes"
        )
    return record_id


def _refget_id(sequence: str) -> str:
    digest = hashlib.sha512(sequence.upper().encode("ascii")).digest()[:24]
    return "SQ." + base64.urlsafe_b64encode(digest).decode("ascii").rstrip("=")


@dataclass(frozen=True)
class SequenceSegment:
    line: int
    normalized_start: int
    length: int

    def to_dict(self) -> dict[str, int]:
        return {"line": self.line, "normalized_start": self.normalized_start, "length": self.length}


@dataclass(frozen=True)
class FastaRecord:
    record_id: str
    description: str
    sequence: str
    segments: tuple[SequenceSegment, ...]


@dataclass(frozen=True)
class ReferenceInterval:
    assembly: str
    contig: str
    start: int
    end: int
    coordinate_system: str = COORDINATE_SYSTEM
    orientation: str = "forward"
    aliases: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "assembly": self.assembly,
            "contig": self.contig,
            "start": self.start,
            "end": self.end,
            "coordinate_system": self.coordinate_system,
            "orientation": self.orientation,
            "aliases": list(self.aliases),
        }


@dataclass(frozen=True)
class ProvenanceSource:
    source_id: str
    kind: str
    uri: str
    version: str
    sha256: str | None = None

    def to_dict(self) -> dict[str, Any]:
        result: dict[str, Any] = {
            "id": self.source_id,
            "kind": self.kind,
            "uri": self.uri,
            "version": self.version,
        }
        if self.sha256 is not None:
            result["sha256"] = self.sha256
        return result


@dataclass(frozen=True)
class SequenceContext:
    record_id: str
    reference: ReferenceInterval | None
    provenance: tuple[ProvenanceSource, ...]


@dataclass(frozen=True)
class SequenceArtifact:
    fasta_sha256: str
    context_sha256: str | None
    record_id: str
    description: str
    sequence: str
    reference: ReferenceInterval | None
    provenance: tuple[ProvenanceSource, ...]
    source_map: tuple[SequenceSegment, ...]

    @property
    def sequence_sha256(self) -> str:
        return _sha256_bytes(self.sequence.encode("ascii"))

    @property
    def canonical_sha256(self) -> str:
        return _sha256_bytes(self.sequence.upper().encode("ascii"))

    @property
    def refget_id(self) -> str:
        return _refget_id(self.sequence)

    @property
    def stats(self) -> dict[str, int | float]:
        canonical = self.sequence.upper()
        gc_bases = canonical.count("G") + canonical.count("C")
        return {
            "bases": len(self.sequence),
            "gc_bases": gc_bases,
            "gc_fraction": round(gc_bases / len(self.sequence), 12),
            "masked_bases": sum(base.islower() for base in self.sequence),
            "ambiguous_bases": sum(base not in "ACGT" for base in canonical),
        }

    def _sequence_ir(self) -> dict[str, Any]:
        return {
            "record_id": self.record_id,
            "description": self.description,
            "sequence": self.sequence,
            "sequence_sha256": self.sequence_sha256,
            "canonical_sha256": self.canonical_sha256,
            "refget_id": self.refget_id,
            "reference": None if self.reference is None else self.reference.to_dict(),
            "provenance": [source.to_dict() for source in self.provenance],
        }

    def _core_dict(self) -> dict[str, Any]:
        sequence_ir = self._sequence_ir()
        return {
            "format": "brain01.sequence-ir",
            "version": 2,
            "compiler": {
                "name": "brainc-dna",
                "version": SEQUENCE_COMPILER_VERSION,
                "passes": list(SEQUENCE_PASSES),
            },
            "inputs": {"fasta_sha256": self.fasta_sha256, "context_sha256": self.context_sha256},
            "sequence_ir": sequence_ir,
            "ir_sha256": _canonical_digest(sequence_ir),
            "source_map": {"sequence_segments": [segment.to_dict() for segment in self.source_map]},
        }

    @property
    def digest(self) -> str:
        return _canonical_digest(self._core_dict())

    def to_dict(self) -> dict[str, Any]:
        payload = self._core_dict()
        payload["artifact_sha256"] = self.digest
        return payload

    def save(self, path: str | Path) -> None:
        raw = _artifact_bytes(self.to_dict())
        try:
            atomic_write_file(
                path,
                raw,
                maximum_bytes=MAX_JSON_BYTES,
                label="sequence artifact output",
            )
        except BoundedIOError as failure:
            raise SequenceCompilerError(f"DNA038: cannot save sequence artifact: {failure}") from failure


def _parse_fasta(raw: bytes) -> FastaRecord:
    try:
        source = raw.decode("utf-8")
    except UnicodeDecodeError as failure:
        raise SequenceCompilerError("DNA001: FASTA must be UTF-8 text") from failure
    line_breaks = source.count("\n") + source.count("\r") - source.count("\r\n")
    line_breaks += sum(
        source.count(separator)
        for separator in ("\v", "\f", "\x1c", "\x1d", "\x1e", "\x85", "\u2028", "\u2029")
    )
    if line_breaks > MAX_SOURCE_LINE_BREAKS:
        raise SequenceCompilerError(
            f"DNA008: FASTA exceeds {MAX_SOURCE_LINE_BREAKS} physical line breaks"
        )
    lines = source.splitlines()
    if not lines or not lines[0].startswith(">"):
        raise SequenceCompilerError("DNA002: FASTA must begin with a '>' defline")
    defline = lines[0][1:]
    if not defline.strip():
        raise SequenceCompilerError("DNA003: FASTA record identifier is empty")
    first = defline.split(maxsplit=1)
    record_id = first[0]
    description = first[1] if len(first) == 2 else ""
    if len(record_id.encode("utf-8")) > MAX_IDENTIFIER_BYTES:
        raise SequenceCompilerError(
            f"DNA009: FASTA record identifier exceeds {MAX_IDENTIFIER_BYTES} UTF-8 bytes"
        )
    if len(description.encode("utf-8")) > MAX_STRING_BYTES:
        raise SequenceCompilerError(
            f"DNA010: FASTA description exceeds {MAX_STRING_BYTES} UTF-8 bytes"
        )
    if any(character.isspace() for character in record_id):
        raise SequenceCompilerError("DNA003: FASTA record identifier contains whitespace")
    if any(
        ord(character) < 0x20 or 0x7F <= ord(character) <= 0x9F
        for character in record_id
    ):
        raise SequenceCompilerError("DNA003: FASTA record identifier contains a control character")
    if any(
        ord(character) < 0x20 or 0x7F <= ord(character) <= 0x9F
        for character in description
    ):
        raise SequenceCompilerError("DNA010: FASTA description contains a control character")
    sequence_parts: list[str] = []
    segments: list[SequenceSegment] = []
    offset = 0
    for line_number, line in enumerate(lines[1:], start=2):
        if line.startswith(">"):
            raise SequenceCompilerError(f"DNA004: exactly one FASTA record is required (line {line_number})")
        if not line:
            continue
        for column, character in enumerate(line, start=1):
            if character.isspace():
                raise SequenceCompilerError(
                    f"DNA005: whitespace inside sequence at line {line_number}, column {column}"
                )
            if character.upper() not in IUPAC_DNA:
                raise SequenceCompilerError(
                    f"DNA006: invalid IUPAC DNA symbol {character!r} at line {line_number}, column {column}"
                )
        sequence_parts.append(line)
        segments.append(SequenceSegment(line_number, offset, len(line)))
        offset += len(line)
    sequence = "".join(sequence_parts)
    if not sequence:
        raise SequenceCompilerError("DNA007: FASTA sequence is empty")
    if len(sequence.encode("ascii")) > MAX_STRING_BYTES:
        raise SequenceCompilerError(
            f"DNA011: FASTA sequence exceeds emitted JSON string limit {MAX_STRING_BYTES}"
        )
    return FastaRecord(record_id, description, sequence, tuple(segments))


def _parse_context(raw: bytes) -> SequenceContext:
    value = _load_json_bytes(raw, "sequence context")
    _require_keys(value, {"format", "version", "record_id", "reference", "provenance"}, set(), "context")
    if (
        value["format"] != "brain01.sequence-context"
        or type(value["version"]) is not int
        or value["version"] != 1
    ):
        raise SequenceCompilerError("DNA024: unsupported sequence context format or version")
    record_id = _require_record_id(value["record_id"], "context.record_id")

    reference_value = value["reference"]
    reference: ReferenceInterval | None
    if reference_value is None:
        reference = None
    elif isinstance(reference_value, dict):
        _require_keys(
            reference_value,
            {"assembly", "contig", "start", "end", "coordinate_system", "orientation", "aliases"},
            set(),
            "context.reference",
        )
        assembly = _require_text(reference_value["assembly"], "context.reference.assembly")
        contig = _require_text(reference_value["contig"], "context.reference.contig")
        start = reference_value["start"]
        end = reference_value["end"]
        if isinstance(start, bool) or isinstance(end, bool) or not isinstance(start, int) or not isinstance(end, int):
            raise SequenceCompilerError("DNA025: reference start and end must be integers")
        if start < 0 or end <= start:
            raise SequenceCompilerError("DNA026: invalid reference interval")
        if reference_value["coordinate_system"] != COORDINATE_SYSTEM:
            raise SequenceCompilerError(f"DNA027: coordinate_system must be {COORDINATE_SYSTEM!r}")
        orientation = reference_value["orientation"]
        if not isinstance(orientation, str) or orientation not in {"forward", "reverse"}:
            raise SequenceCompilerError("DNA028: orientation must be 'forward' or 'reverse'")
        aliases_value = reference_value["aliases"]
        if not isinstance(aliases_value, list) or any(
            not isinstance(alias, str) or not alias or alias != alias.strip() for alias in aliases_value
        ):
            raise SequenceCompilerError("DNA029: reference aliases must be non-empty trimmed strings")
        if len(set(aliases_value)) != len(aliases_value):
            raise SequenceCompilerError("DNA029: reference aliases must be unique")
        reference = ReferenceInterval(assembly, contig, start, end, COORDINATE_SYSTEM, orientation, tuple(aliases_value))
    else:
        raise SequenceCompilerError("DNA030: context.reference must be an object or null")

    provenance_value = value["provenance"]
    if not isinstance(provenance_value, list):
        raise SequenceCompilerError("DNA031: context.provenance must be an array")
    provenance: list[ProvenanceSource] = []
    seen_ids: set[str] = set()
    for index, item in enumerate(provenance_value):
        if not isinstance(item, dict):
            raise SequenceCompilerError(f"DNA032: provenance item {index} must be an object")
        _require_keys(item, {"id", "kind", "uri", "version"}, {"sha256"}, f"provenance item {index}")
        source_id = _require_text(item["id"], f"provenance[{index}].id")
        if source_id in seen_ids:
            raise SequenceCompilerError(f"DNA033: duplicate provenance id: {source_id}")
        seen_ids.add(source_id)
        sha256 = item.get("sha256")
        if "sha256" in item and (not isinstance(sha256, str) or SHA256_PATTERN.fullmatch(sha256) is None):
            raise SequenceCompilerError(f"DNA034: provenance[{index}].sha256 must be lowercase SHA-256")
        provenance.append(
            ProvenanceSource(
                source_id,
                _require_text(item["kind"], f"provenance[{index}].kind"),
                _require_text(item["uri"], f"provenance[{index}].uri"),
                _require_text(item["version"], f"provenance[{index}].version"),
                sha256,
            )
        )
    return SequenceContext(record_id, reference, tuple(provenance))


class SequenceCompiler:
    """Compile one standard FASTA record and optional explicit context."""

    def compile_file(self, path: str | Path, context_path: str | Path | None = None) -> SequenceArtifact:
        try:
            fasta_raw = read_regular_file(path, maximum_bytes=MAX_INPUT_BYTES, label="FASTA source")
            context_raw = (
                None
                if context_path is None
                else read_regular_file(
                    context_path,
                    maximum_bytes=MAX_JSON_BYTES,
                    label="sequence context",
                )
            )
        except BoundedIOError as failure:
            raise SequenceCompilerError(f"DNA039: unsafe or oversized input file: {failure}") from failure
        return self.compile_bytes(fasta_raw, context_raw)

    def compile_text(self, source: str, source_name: str = "<memory>") -> SequenceArtifact:
        del source_name  # Filenames are deliberately excluded from semantic identity.
        return self.compile_bytes(source.encode("utf-8"), None)

    def compile_bytes(self, fasta_raw: bytes, context_raw: bytes | None = None) -> SequenceArtifact:
        if type(fasta_raw) is not bytes or len(fasta_raw) > MAX_INPUT_BYTES:
            raise SequenceCompilerError(
                f"DNA012: FASTA source must be bytes at most {MAX_INPUT_BYTES} bytes"
            )
        if context_raw is not None and (
            type(context_raw) is not bytes or len(context_raw) > MAX_JSON_BYTES
        ):
            raise SequenceCompilerError(
                f"DNA020: sequence context must be bytes at most {MAX_JSON_BYTES} bytes"
            )
        record = _parse_fasta(fasta_raw)
        context = None if context_raw is None else _parse_context(context_raw)
        if context is not None and context.record_id != record.record_id:
            raise SequenceCompilerError("DNA035: context record_id does not match FASTA record identifier")
        reference = None if context is None else context.reference
        if reference is not None and reference.end - reference.start != len(record.sequence):
            raise SequenceCompilerError(
                f"DNA036: reference span {reference.end - reference.start} does not equal sequence length "
                f"{len(record.sequence)}"
            )
        artifact = SequenceArtifact(
            fasta_sha256=_sha256_bytes(fasta_raw),
            context_sha256=None if context_raw is None else _sha256_bytes(context_raw),
            record_id=record.record_id,
            description=record.description,
            sequence=record.sequence,
            reference=reference,
            provenance=() if context is None else context.provenance,
            source_map=record.segments,
        )
        _artifact_bytes(artifact.to_dict())
        return artifact


def _as_artifact(payload: dict[str, Any]) -> SequenceArtifact:
    _require_keys(
        payload,
        {"format", "version", "compiler", "inputs", "sequence_ir", "ir_sha256", "source_map", "artifact_sha256"},
        set(),
        "sequence artifact",
    )
    if (
        payload["format"] != "brain01.sequence-ir"
        or type(payload["version"]) is not int
        or payload["version"] != 2
    ):
        raise SequenceCompilerError("DNA040: unsupported sequence artifact format or version")
    compiler = payload["compiler"]
    if not isinstance(compiler, dict):
        raise SequenceCompilerError("DNA041: compiler must be an object")
    _require_keys(compiler, {"name", "version", "passes"}, set(), "compiler")
    if compiler != {"name": "brainc-dna", "version": SEQUENCE_COMPILER_VERSION, "passes": list(SEQUENCE_PASSES)}:
        raise SequenceCompilerError("DNA042: unsupported compiler identity or pass set")
    inputs = payload["inputs"]
    if not isinstance(inputs, dict):
        raise SequenceCompilerError("DNA043: inputs must be an object")
    _require_keys(inputs, {"fasta_sha256", "context_sha256"}, set(), "inputs")
    if not isinstance(inputs["fasta_sha256"], str) or SHA256_PATTERN.fullmatch(inputs["fasta_sha256"]) is None:
        raise SequenceCompilerError("DNA044: inputs.fasta_sha256 is invalid")
    if inputs["context_sha256"] is not None and (
        not isinstance(inputs["context_sha256"], str) or SHA256_PATTERN.fullmatch(inputs["context_sha256"]) is None
    ):
        raise SequenceCompilerError("DNA044: inputs.context_sha256 is invalid")
    ir = payload["sequence_ir"]
    if not isinstance(ir, dict):
        raise SequenceCompilerError("DNA045: sequence_ir must be an object")
    _require_keys(
        ir,
        {
            "record_id", "description", "sequence", "sequence_sha256", "canonical_sha256",
            "refget_id", "reference", "provenance",
        },
        set(),
        "sequence_ir",
    )
    record_id = _require_record_id(ir["record_id"], "sequence_ir.record_id")
    if not isinstance(ir["description"], str) or ir["description"].splitlines() not in (
        [],
        [ir["description"]],
    ):
        raise SequenceCompilerError("DNA046: sequence_ir.description must be one text line")
    sequence = ir["sequence"]
    if not isinstance(sequence, str) or not sequence or any(base.upper() not in IUPAC_DNA for base in sequence):
        raise SequenceCompilerError("DNA047: sequence_ir.sequence is not valid IUPAC DNA")

    reference: ReferenceInterval | None = None
    if ir["reference"] is not None:
        context_payload = {
            "format": "brain01.sequence-context", "version": 1, "record_id": record_id,
            "reference": ir["reference"], "provenance": ir["provenance"],
        }
        context = _parse_context(_canonical_bytes(context_payload))
        reference = context.reference
        provenance = context.provenance
    else:
        context_payload = {
            "format": "brain01.sequence-context", "version": 1, "record_id": record_id,
            "reference": None, "provenance": ir["provenance"],
        }
        provenance = _parse_context(_canonical_bytes(context_payload)).provenance
    if reference is not None and reference.end - reference.start != len(sequence):
        raise SequenceCompilerError("DNA048: artifact reference span does not equal sequence length")

    source_map = payload["source_map"]
    if not isinstance(source_map, dict):
        raise SequenceCompilerError("DNA049: source_map must be an object")
    _require_keys(source_map, {"sequence_segments"}, set(), "source_map")
    segments_value = source_map["sequence_segments"]
    if not isinstance(segments_value, list) or not segments_value:
        raise SequenceCompilerError("DNA050: source_map.sequence_segments must be a non-empty array")
    segments: list[SequenceSegment] = []
    expected_offset = 0
    last_line = 1
    for index, item in enumerate(segments_value):
        if not isinstance(item, dict):
            raise SequenceCompilerError(f"DNA051: sequence segment {index} must be an object")
        _require_keys(item, {"line", "normalized_start", "length"}, set(), f"sequence segment {index}")
        line, start, length = item["line"], item["normalized_start"], item["length"]
        if any(isinstance(value, bool) or not isinstance(value, int) for value in (line, start, length)):
            raise SequenceCompilerError(f"DNA052: sequence segment {index} values must be integers")
        if (
            line <= last_line
            or line > MAX_SOURCE_LINE_BREAKS + 1
            or start != expected_offset
            or length <= 0
        ):
            raise SequenceCompilerError(f"DNA053: sequence segment {index} is discontinuous")
        segments.append(SequenceSegment(line, start, length))
        expected_offset += length
        last_line = line
    if expected_offset != len(sequence):
        raise SequenceCompilerError("DNA054: source map length does not equal sequence length")

    artifact = SequenceArtifact(
        inputs["fasta_sha256"], inputs["context_sha256"], record_id, ir["description"], sequence,
        reference, provenance, tuple(segments),
    )
    if artifact.context_sha256 is None and (
        artifact.reference is not None or artifact.provenance
    ):
        raise SequenceCompilerError(
            "DNA058: sequence context assertions require inputs.context_sha256"
        )
    expected_sequence_fields = artifact._sequence_ir()
    if ir != expected_sequence_fields:
        raise SequenceCompilerError("DNA055: derived sequence IR fields do not recompute")
    if payload["ir_sha256"] != _canonical_digest(ir):
        raise SequenceCompilerError("DNA056: sequence IR digest mismatch")
    if payload["artifact_sha256"] != artifact.digest:
        raise SequenceCompilerError("DNA057: sequence artifact digest mismatch")
    return artifact


def load_sequence_artifact(path: str | Path) -> SequenceArtifact:
    try:
        raw = read_regular_file(path, maximum_bytes=MAX_JSON_BYTES, label="sequence artifact")
    except BoundedIOError as failure:
        raise SequenceCompilerError(f"DNA039: unsafe or oversized sequence artifact: {failure}") from failure
    payload = _load_json_bytes(raw, "sequence artifact")
    return _as_artifact(payload)

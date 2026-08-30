"""Raw-DNA and multi-FASTA ingress for collections of Sequence IR artifacts."""

from __future__ import annotations

from dataclasses import dataclass
import gzip
import hashlib
import json
from pathlib import Path
import re
from typing import Any
import zlib

from .sequence import (
    SequenceArtifact,
    SequenceCompiler,
    SequenceCompilerError,
    SequenceSegment,
    _as_artifact as _as_sequence_artifact,
)


FORMAT = "brain01.sequence-collection-ir"
VERSION = 1
COMPILER = {"name": "brainc-dna-collection", "version": "0.1.0"}
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}\Z")


class SequenceCollectionError(SequenceCompilerError):
    """A deterministic raw or multi-FASTA compilation failure."""


def _sha256(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _canonical_bytes(value: Any) -> bytes:
    try:
        return json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        ).encode("utf-8")
    except (TypeError, ValueError) as failure:
        raise SequenceCollectionError(f"DNAC001: value is not canonical JSON: {failure}") from failure


def _digest(value: Any) -> str:
    return _sha256(_canonical_bytes(value))


def _reject_constant(value: str) -> None:
    raise SequenceCollectionError(f"DNAC002: non-finite JSON number is not allowed: {value}")


def _reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise SequenceCollectionError(f"DNAC003: duplicate JSON key: {key}")
        result[key] = value
    return result


def _load_json(raw: bytes) -> dict[str, Any]:
    try:
        value = json.loads(
            raw.decode("utf-8"),
            object_pairs_hook=_reject_duplicates,
            parse_constant=_reject_constant,
        )
    except SequenceCollectionError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as failure:
        raise SequenceCollectionError(f"DNAC004: invalid collection JSON: {failure}") from failure
    if not isinstance(value, dict):
        raise SequenceCollectionError("DNAC004: collection artifact must be an object")
    return value


def _keys(value: dict[str, Any], expected: set[str], label: str) -> None:
    missing = sorted(expected - value.keys())
    extra = sorted(value.keys() - expected)
    if missing or extra:
        details = []
        if missing:
            details.append(f"missing {', '.join(missing)}")
        if extra:
            details.append(f"unknown {', '.join(extra)}")
        raise SequenceCollectionError(f"DNAC005: {label} has invalid keys: {'; '.join(details)}")


def _mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise SequenceCollectionError(f"DNAC006: {label} must be an object")
    return value


def _valid_digest(value: Any, label: str) -> str:
    if not isinstance(value, str) or SHA256_PATTERN.fullmatch(value) is None:
        raise SequenceCollectionError(f"DNAC007: {label} must be a lowercase SHA-256 digest")
    return value


def _unwrap(raw: bytes) -> tuple[bytes, str | None]:
    if not isinstance(raw, bytes):
        raise SequenceCollectionError("DNAC008: source must be bytes")
    if not raw.startswith(b"\x1f\x8b"):
        return raw, None
    try:
        return gzip.decompress(raw), _sha256(raw)
    except (gzip.BadGzipFile, EOFError, OSError, zlib.error) as failure:
        raise SequenceCollectionError(f"DNAC009: malformed gzip source: {failure}") from failure


def _source_segments(value: Any, sequence_length: int, label: str) -> tuple[SequenceSegment, ...]:
    source_map = _mapping(value, label)
    _keys(source_map, {"sequence_segments"}, label)
    items = source_map["sequence_segments"]
    if not isinstance(items, list) or not items:
        raise SequenceCollectionError(f"DNAC010: {label}.sequence_segments must be a non-empty array")
    result: list[SequenceSegment] = []
    offset = 0
    last_line = 0
    for index, item_value in enumerate(items):
        item = _mapping(item_value, f"{label}.sequence_segments[{index}]")
        _keys(item, {"line", "normalized_start", "length"}, f"{label}.sequence_segments[{index}]")
        line, start, length = item["line"], item["normalized_start"], item["length"]
        if any(type(number) is not int for number in (line, start, length)):
            raise SequenceCollectionError(f"DNAC011: {label} segment values must be integers")
        if line <= last_line or start != offset or length <= 0:
            raise SequenceCollectionError(f"DNAC012: {label} segments are discontinuous")
        result.append(SequenceSegment(line, start, length))
        offset += length
        last_line = line
    if offset != sequence_length:
        raise SequenceCollectionError(f"DNAC013: {label} length does not equal sequence length")
    return tuple(result)


@dataclass(frozen=True)
class SequenceCollectionMember:
    artifact: SequenceArtifact
    input_source_map: tuple[SequenceSegment, ...]

    def to_dict(self) -> dict[str, Any]:
        sequence_payload = self.artifact.to_dict()
        return {
            "record_id": self.artifact.record_id,
            "ir_sha256": sequence_payload["ir_sha256"],
            "artifact_sha256": sequence_payload["artifact_sha256"],
            "input_source_map": {
                "sequence_segments": [segment.to_dict() for segment in self.input_source_map]
            },
            "sequence_artifact": sequence_payload,
        }


@dataclass(frozen=True)
class SequenceCollectionArtifact:
    input_kind: str
    raw_sha256: str
    logical_sha256: str
    compressed_sha256: str | None
    members: tuple[SequenceCollectionMember, ...]

    def _collection_ir(self) -> dict[str, Any]:
        return {
            "members": [
                {"record_id": member.artifact.record_id, "ir_sha256": member.artifact.to_dict()["ir_sha256"]}
                for member in self.members
            ]
        }

    def _core_dict(self) -> dict[str, Any]:
        collection_ir = self._collection_ir()
        return {
            "format": FORMAT,
            "version": VERSION,
            "compiler": dict(COMPILER),
            "inputs": {
                "kind": self.input_kind,
                "wrapper": "gzip" if self.compressed_sha256 is not None else None,
                "raw_sha256": self.raw_sha256,
                "logical_sha256": self.logical_sha256,
                "compressed_sha256": self.compressed_sha256,
            },
            "collection_ir": collection_ir,
            "collection_ir_sha256": _digest(collection_ir),
            "members": [member.to_dict() for member in self.members],
        }

    @property
    def digest(self) -> str:
        return _digest(self._core_dict())

    def to_dict(self) -> dict[str, Any]:
        value = self._core_dict()
        value["artifact_sha256"] = self.digest
        return value

    def save(self, path: str | Path) -> None:
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(
            (
                json.dumps(self.to_dict(), indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False)
                + "\n"
            ).encode("utf-8")
        )


def _collection(
    kind: str,
    raw: bytes,
    logical: bytes,
    compressed_sha256: str | None,
    members: tuple[SequenceCollectionMember, ...],
) -> SequenceCollectionArtifact:
    artifact = SequenceCollectionArtifact(
        input_kind=kind,
        raw_sha256=_sha256(raw),
        logical_sha256=_sha256(logical),
        compressed_sha256=compressed_sha256,
        members=members,
    )
    _validate_artifact(artifact)
    return artifact


class SequenceCollectionCompiler:
    """Compile explicit raw IUPAC DNA or a FASTA record collection."""

    def compile_raw(self, sequence: str | bytes, record_id: str) -> SequenceCollectionArtifact:
        if not isinstance(record_id, str) or not record_id or any(character.isspace() for character in record_id):
            raise SequenceCollectionError("DNAC014: raw input record_id must be a non-empty whitespace-free string")
        if isinstance(sequence, str):
            try:
                raw = sequence.encode("utf-8")
            except UnicodeEncodeError as failure:
                raise SequenceCollectionError("DNAC016: raw sequence must be ASCII IUPAC DNA") from failure
        elif isinstance(sequence, bytes):
            raw = sequence
        else:
            raise SequenceCollectionError("DNAC015: raw sequence must be text or bytes")
        logical, compressed_sha256 = _unwrap(raw)
        try:
            text = logical.decode("ascii")
        except UnicodeDecodeError as failure:
            raise SequenceCollectionError("DNAC016: raw sequence must be ASCII IUPAC DNA") from failure
        if not text or any(character.isspace() for character in text):
            raise SequenceCollectionError("DNAC017: raw sequence must be non-empty and contain no whitespace")
        try:
            fasta = f">{record_id}\n{text}\n".encode("utf-8")
        except UnicodeEncodeError as failure:
            raise SequenceCollectionError("DNAC014: raw input record_id must be valid UTF-8 text") from failure
        try:
            member_artifact = SequenceCompiler().compile_bytes(fasta)
        except SequenceCompilerError as failure:
            raise SequenceCollectionError(f"DNAC018: invalid raw sequence: {failure}") from failure
        source_map = (SequenceSegment(line=1, normalized_start=0, length=len(text)),)
        member = SequenceCollectionMember(member_artifact, source_map)
        return _collection("raw-iupac", raw, logical, compressed_sha256, (member,))

    def compile_fasta_bytes(self, raw: bytes) -> SequenceCollectionArtifact:
        logical, compressed_sha256 = _unwrap(raw)
        lines = logical.splitlines(keepends=True)
        if not lines or not lines[0].startswith(b">"):
            raise SequenceCollectionError("DNAC019: FASTA collection must begin with a '>' defline")
        starts = [index for index, line in enumerate(lines) if line.startswith(b">")]
        members: list[SequenceCollectionMember] = []
        seen: set[str] = set()
        for position, start in enumerate(starts):
            end = starts[position + 1] if position + 1 < len(starts) else len(lines)
            record_raw = b"".join(lines[start:end])
            try:
                local = SequenceCompiler().compile_bytes(record_raw)
            except SequenceCompilerError as failure:
                raise SequenceCollectionError(
                    f"DNAC020: invalid FASTA record beginning at line {start + 1}: {failure}"
                ) from failure
            if local.record_id in seen:
                raise SequenceCollectionError(f"DNAC021: duplicate FASTA record identifier: {local.record_id}")
            seen.add(local.record_id)
            shifted = tuple(
                SequenceSegment(segment.line + start, segment.normalized_start, segment.length)
                for segment in local.source_map
            )
            artifact = SequenceArtifact(
                local.fasta_sha256,
                local.context_sha256,
                local.record_id,
                local.description,
                local.sequence,
                local.reference,
                local.provenance,
                shifted,
            )
            members.append(SequenceCollectionMember(artifact, shifted))
        return _collection("fasta", raw, logical, compressed_sha256, tuple(members))

    def compile_file(self, path: str | Path) -> SequenceCollectionArtifact:
        return self.compile_fasta_bytes(Path(path).read_bytes())


def _validate_artifact(artifact: SequenceCollectionArtifact) -> None:
    if artifact.input_kind not in {"raw-iupac", "fasta"}:
        raise SequenceCollectionError("DNAC022: unsupported collection input kind")
    _valid_digest(artifact.raw_sha256, "inputs.raw_sha256")
    _valid_digest(artifact.logical_sha256, "inputs.logical_sha256")
    if artifact.compressed_sha256 is None:
        if artifact.raw_sha256 != artifact.logical_sha256:
            raise SequenceCollectionError("DNAC023: unwrapped raw and logical digests must match")
    else:
        _valid_digest(artifact.compressed_sha256, "inputs.compressed_sha256")
        if artifact.compressed_sha256 != artifact.raw_sha256:
            raise SequenceCollectionError("DNAC024: compressed and raw digests must match")
    if not artifact.members:
        raise SequenceCollectionError("DNAC025: collection must contain at least one member")
    ids: set[str] = set()
    for index, member in enumerate(artifact.members):
        record_id = member.artifact.record_id
        if record_id in ids:
            raise SequenceCollectionError(f"DNAC026: duplicate member record_id: {record_id}")
        ids.add(record_id)
        _as_sequence_artifact(member.artifact.to_dict())
        segments = _source_segments(
            {"sequence_segments": [segment.to_dict() for segment in member.input_source_map]},
            len(member.artifact.sequence),
            f"members[{index}].input_source_map",
        )
        if artifact.input_kind == "raw-iupac":
            if len(artifact.members) != 1 or segments != (
                SequenceSegment(1, 0, len(member.artifact.sequence)),
            ):
                raise SequenceCollectionError("DNAC027: raw input must map one member to source line 1")
        elif segments != member.artifact.source_map:
            raise SequenceCollectionError("DNAC028: FASTA input source map must match its member artifact")


def _as_collection(payload: dict[str, Any]) -> SequenceCollectionArtifact:
    _keys(
        payload,
        {
            "format", "version", "compiler", "inputs", "collection_ir",
            "collection_ir_sha256", "members", "artifact_sha256",
        },
        "collection artifact",
    )
    if payload["format"] != FORMAT or type(payload["version"]) is not int or payload["version"] != VERSION:
        raise SequenceCollectionError("DNAC029: unsupported collection format or version")
    if _canonical_bytes(payload["compiler"]) != _canonical_bytes(COMPILER):
        raise SequenceCollectionError("DNAC030: unsupported collection compiler")
    inputs = _mapping(payload["inputs"], "inputs")
    _keys(
        inputs,
        {"kind", "wrapper", "raw_sha256", "logical_sha256", "compressed_sha256"},
        "inputs",
    )
    kind = inputs["kind"]
    if not isinstance(kind, str) or kind not in {"raw-iupac", "fasta"}:
        raise SequenceCollectionError("DNAC022: unsupported collection input kind")
    compressed = inputs["compressed_sha256"]
    if inputs["wrapper"] is None:
        if compressed is not None:
            raise SequenceCollectionError("DNAC031: unwrapped input cannot have a compressed digest")
    elif inputs["wrapper"] != "gzip" or compressed is None:
        raise SequenceCollectionError("DNAC031: wrapper must be null or gzip with a compressed digest")
    members_value = payload["members"]
    if not isinstance(members_value, list) or not members_value:
        raise SequenceCollectionError("DNAC025: collection must contain at least one member")
    members: list[SequenceCollectionMember] = []
    for index, member_value in enumerate(members_value):
        member_payload = _mapping(member_value, f"members[{index}]")
        _keys(
            member_payload,
            {"record_id", "ir_sha256", "artifact_sha256", "input_source_map", "sequence_artifact"},
            f"members[{index}]",
        )
        sequence_payload = _mapping(member_payload["sequence_artifact"], f"members[{index}].sequence_artifact")
        try:
            sequence_artifact = _as_sequence_artifact(sequence_payload)
        except SequenceCompilerError as failure:
            raise SequenceCollectionError(f"DNAC032: invalid member Sequence IR: {failure}") from failure
        if not isinstance(member_payload["record_id"], str) or member_payload["record_id"] != sequence_artifact.record_id:
            raise SequenceCollectionError(f"DNAC033: members[{index}].record_id does not match Sequence IR")
        if _valid_digest(member_payload["ir_sha256"], f"members[{index}].ir_sha256") != sequence_payload["ir_sha256"]:
            raise SequenceCollectionError(f"DNAC034: members[{index}].ir_sha256 does not match Sequence IR")
        if _valid_digest(member_payload["artifact_sha256"], f"members[{index}].artifact_sha256") != sequence_artifact.digest:
            raise SequenceCollectionError(f"DNAC035: members[{index}].artifact_sha256 does not match Sequence IR")
        source_map = _source_segments(
            member_payload["input_source_map"],
            len(sequence_artifact.sequence),
            f"members[{index}].input_source_map",
        )
        members.append(SequenceCollectionMember(sequence_artifact, source_map))
    artifact = SequenceCollectionArtifact(
        input_kind=kind,
        raw_sha256=_valid_digest(inputs["raw_sha256"], "inputs.raw_sha256"),
        logical_sha256=_valid_digest(inputs["logical_sha256"], "inputs.logical_sha256"),
        compressed_sha256=None if compressed is None else _valid_digest(compressed, "inputs.compressed_sha256"),
        members=tuple(members),
    )
    _validate_artifact(artifact)
    expected_ir = artifact._collection_ir()
    collection_ir = _mapping(payload["collection_ir"], "collection_ir")
    if _canonical_bytes(collection_ir) != _canonical_bytes(expected_ir):
        raise SequenceCollectionError("DNAC036: collection IR does not match its members")
    if _valid_digest(payload["collection_ir_sha256"], "collection_ir_sha256") != _digest(expected_ir):
        raise SequenceCollectionError("DNAC037: collection IR digest mismatch")
    if _valid_digest(payload["artifact_sha256"], "artifact_sha256") != artifact.digest:
        raise SequenceCollectionError("DNAC038: collection artifact digest mismatch")
    return artifact


def load_sequence_collection(path: str | Path) -> SequenceCollectionArtifact:
    return _as_collection(_load_json(Path(path).read_bytes()))


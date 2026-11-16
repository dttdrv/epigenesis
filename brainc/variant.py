"""Sealed reference-bound variant collections with GA4GH VRS 2.x identifiers."""

from __future__ import annotations

import base64
from collections.abc import Mapping, Sequence
import copy
import hashlib
import os
from pathlib import Path
import re
from typing import Any
import zlib

from ._canonical import (
    ContractError,
    SAFE_INTEGER,
    canonical_bytes,
    digest,
    save_artifact,
)
from ._io import (
    BoundedIOError,
    MAX_DECOMPRESSED_BYTES,
    MAX_IDENTIFIER_BYTES,
    MAX_INPUT_BYTES,
    MAX_JSON_BYTES,
    MAX_JSON_DEPTH,
    MAX_JSON_MEMBERS,
    MAX_STRING_BYTES,
    load_json_object,
    pretty_json_bytes,
    read_regular_file,
)
from .source import (
    EXTERNAL_PROFILE,
    FASTA_PROFILE,
    GENBANK_PROFILE,
    GFF3_PROFILE,
    RAW_PROFILE,
    REFERENCE_FASTA_PROFILE,
    SourceBundle,
    SourceError,
    load_source_bundle,
    validate_source_descriptor,
)


FORMAT = "brainc.variant-collection"
VERSION = 1
VRS_VERSION = "2.x"
PRODUCER = {
    "name": "brainc-variant",
    "version": "1.0.0",
    "passes": [
        "validate-source-bundle",
        "extract-reference-sequences",
        "normalize-and-verify-variants",
        "compute-vrs-v2-identifiers",
        "seal-variant-collection",
    ],
}

IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")
MAX_VARIANTS = 100_000
MAX_VARIANT_BASES = 1_000_000
_SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
_REFGET_RE = re.compile(r"SQ\.[A-Za-z0-9_-]{32}\Z")
_VRS_DIGEST_RE = re.compile(r"[A-Za-z0-9_-]{32}\Z")
_GZIP_HEADER = b"\x1f\x8b"


class VariantError(ContractError):
    """A variant input or sealed variant collection violates the contract."""


def _fail(detail: str) -> VariantError:
    return VariantError(f"VARIANT001: {detail}")


def _sha512t24u(raw: bytes) -> str:
    return base64.urlsafe_b64encode(hashlib.sha512(raw).digest()[:24]).decode("ascii").rstrip("=")


def _refget_sequence_id(sequence: str) -> str:
    return "SQ." + _sha512t24u(sequence.encode("ascii"))


def _keys(value: Any, expected: set[str], label: str) -> dict[str, Any]:
    if type(value) is not dict:
        raise _fail(f"{label} must be a JSON object")
    missing = sorted(expected - value.keys())
    extra = sorted(value.keys() - expected)
    if missing or extra:
        raise _fail(
            f"{label} keys invalid; missing={missing or 'none'}, unknown={extra or 'none'}"
        )
    return value


def _exact_equal(left: Any, right: Any) -> bool:
    if type(left) is not type(right):
        return False
    if type(left) is dict:
        return left.keys() == right.keys() and all(_exact_equal(left[k], right[k]) for k in left)
    if type(left) is list:
        return len(left) == len(right) and all(_exact_equal(a, b) for a, b in zip(left, right))
    return left == right


def _nonempty_text(value: Any, label: str, *, maximum_bytes: int = MAX_IDENTIFIER_BYTES) -> str:
    if type(value) is not str or not value or value != value.strip():
        raise _fail(f"{label} must be nonempty trimmed text")
    encoded = value.encode("utf-8")
    if len(encoded) > maximum_bytes:
        raise _fail(f"{label} exceeds {maximum_bytes} UTF-8 bytes")
    return value


def _iupac_bases(value: Any, label: str, *, allow_empty: bool = False) -> str:
    if type(value) is not str:
        raise _fail(f"{label} must be a string")
    if not value and not allow_empty:
        raise _fail(f"{label} must be nonempty")
    if len(value) > MAX_VARIANT_BASES:
        raise _fail(f"{label} exceeds {MAX_VARIANT_BASES} bases")
    if any(ch not in IUPAC_DNA for ch in value):
        raise _fail(f"{label} contains non-IUPAC or lowercase DNA characters")
    return value


def _nonneg_int(value: Any, label: str, *, maximum: int = SAFE_INTEGER) -> int:
    if type(value) is not int or not 0 <= value <= maximum:
        raise _fail(f"{label} must be an integer in [0, {maximum}]")
    return value


def vrs_location_digest(refget_id: str, start: int, end: int) -> str:
    """Compute the GA4GH VRS 2.x sha512t24u digest for a SequenceLocation."""
    if type(refget_id) is not str or _REFGET_RE.fullmatch(refget_id) is None:
        raise _fail("refget_id must match SQ.<32 base64url chars>")
    _nonneg_int(start, "start")
    _nonneg_int(end, "end")
    if end < start:
        raise _fail("SequenceLocation end must be >= start")
    digest_object = {
        "end": end,
        "sequenceReference": {
            "refgetAccession": refget_id,
            "type": "SequenceReference",
        },
        "start": start,
        "type": "SequenceLocation",
    }
    return _sha512t24u(canonical_bytes(digest_object))


def vrs_allele_digest(location_digest: str, alt: str) -> str:
    """Compute the GA4GH VRS 2.x sha512t24u digest for a LiteralSequenceExpression Allele."""
    if type(location_digest) is not str or _VRS_DIGEST_RE.fullmatch(location_digest) is None:
        raise _fail("location_digest must be a 32-character base64url sha512t24u digest")
    _iupac_bases(alt, "alt", allow_empty=True)
    digest_object = {
        "location": location_digest,
        "state": {
            "sequence": alt,
            "type": "LiteralSequenceExpression",
        },
        "type": "Allele",
    }
    return _sha512t24u(canonical_bytes(digest_object))


def vrs_sequence_location(refget_id: str, start: int, end: int) -> dict[str, Any]:
    """Build a closed GA4GH VRS 2.x SequenceLocation object."""
    loc_digest = vrs_location_digest(refget_id, start, end)
    return {
        "id": f"ga4gh:SL.{loc_digest}",
        "type": "SequenceLocation",
        "digest": loc_digest,
        "sequenceReference": {
            "type": "SequenceReference",
            "refgetAccession": refget_id,
        },
        "start": start,
        "end": end,
    }


def vrs_allele(refget_id: str, start: int, end: int, alt: str) -> dict[str, Any]:
    """Build a closed GA4GH VRS 2.x Allele object with computed SL and VA IDs."""
    location = vrs_sequence_location(refget_id, start, end)
    allele_digest = vrs_allele_digest(location["digest"], alt)
    return {
        "id": f"ga4gh:VA.{allele_digest}",
        "type": "Allele",
        "digest": allele_digest,
        "location": location,
        "state": {
            "type": "LiteralSequenceExpression",
            "sequence": alt,
        },
    }


def _bounded_gzip_decompress(raw: bytes) -> bytes:
    if not raw.startswith(_GZIP_HEADER):
        raise _fail("reference sequence does not have a gzip header")
    out = bytearray()
    view = memoryview(raw)
    offset = 0
    while offset < len(view):
        if len(view) - offset < 2 or view[offset : offset + 2].tobytes() != _GZIP_HEADER:
            raise _fail("trailing non-gzip data in reference sequence")
        decomp = zlib.decompressobj(16 + zlib.MAX_WBITS)
        try:
            while not decomp.eof:
                chunk = view[offset : offset + 65536]
                if not chunk:
                    raise _fail("truncated gzip reference sequence")
                remaining = MAX_DECOMPRESSED_BYTES - len(out)
                piece = decomp.decompress(chunk, remaining + 1)
                out.extend(piece)
                if len(out) > MAX_DECOMPRESSED_BYTES:
                    raise _fail("decompressed reference sequence exceeds limit")
                consumed = len(chunk) - len(decomp.unconsumed_tail)
                offset += consumed
                if consumed == 0 and not decomp.eof:
                    raise _fail("malformed gzip reference sequence")
            tail = decomp.flush()
            if len(out) + len(tail) > MAX_DECOMPRESSED_BYTES:
                raise _fail("decompressed reference sequence exceeds limit")
            out.extend(tail)
            if not decomp.eof:
                raise _fail("truncated gzip reference sequence")
        except zlib.error as failure:
            raise _fail(f"malformed gzip reference sequence: {failure}") from failure
    return bytes(out)


def _decode_fasta_bytes(raw: bytes, wrapper: str) -> dict[str, str]:
    if wrapper == "gzip":
        logical = _bounded_gzip_decompress(raw)
    elif wrapper == "identity":
        if raw.startswith(_GZIP_HEADER):
            raise _fail("identity reference sequence has unexpected gzip header")
        logical = raw
    else:
        raise _fail(f"unsupported reference wrapper: {wrapper}")
    try:
        text = logical.decode("ascii")
    except UnicodeDecodeError as failure:
        raise _fail(f"FASTA reference must be ASCII: {failure}") from failure
    records: dict[str, str] = {}
    current_id: str | None = None
    chunks: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if current_id is not None:
                records[current_id] = "".join(chunks)
            header = line[1:].strip()
            if not header:
                raise _fail("FASTA record has empty header")
            current_id = header.split()[0]
            if current_id in records:
                raise _fail(f"duplicate FASTA record id: {current_id}")
            chunks = []
        else:
            stripped = line.strip().upper()
            if not stripped:
                continue
            if current_id is None:
                raise _fail("FASTA bases appear before any header")
            if any(ch not in IUPAC_DNA for ch in stripped):
                raise _fail(f"FASTA record {current_id!r} contains non-IUPAC bases")
            chunks.append(stripped)
    if current_id is not None:
        records[current_id] = "".join(chunks)
    if not records:
        raise _fail("FASTA source contains no records")
    return records


def _extract_source_sequences(
    bundle: SourceBundle,
    source_inputs: Mapping[str, str | Path] | None = None,
) -> dict[str, str]:
    descriptor = bundle.to_dict()
    try:
        validate_source_descriptor(descriptor)
    except SourceError as failure:
        raise _fail(f"invalid source descriptor: {failure}") from failure

    source_ir = descriptor["source_ir"]
    profile = source_ir["profile"]
    artifacts = bundle.artifacts
    expected_records = {item["record_id"]: item for item in source_ir["records"]}
    sequences: dict[str, str] = {}

    if profile in {RAW_PROFILE, FASTA_PROFILE, GFF3_PROFILE}:
        seq_artifact = artifacts.get("sequence")
        if type(seq_artifact) is not dict:
            raise _fail("source bundle is missing native sequence artifact")
        for member in seq_artifact["members"]:
            seq_ir = member["sequence_artifact"]["sequence_ir"]
            rid = seq_ir["record_id"]
            if rid in sequences:
                raise _fail(f"duplicate record_id {rid!r} in sequence artifact")
            sequences[rid] = seq_ir["sequence"].upper()
    elif profile == GENBANK_PROFILE:
        gb_artifact = artifacts.get("genbank")
        if type(gb_artifact) is not dict:
            raise _fail("source bundle is missing native genbank artifact")
        for member in gb_artifact["sequence_collection"]["members"]:
            rid = member["record_id"]
            if rid in sequences:
                raise _fail(f"duplicate record_id {rid!r} in genbank artifact")
            storage = member["sequence"]["storage"]
            sequences[rid] = "".join(storage["chunks"]).upper()
    elif profile == REFERENCE_FASTA_PROFILE:
        if source_inputs is None or "sequence" not in source_inputs:
            raise _fail(
                "fasta-reference-dna/v1 source bundles omit DNA strings; pass source_inputs={'sequence': path}"
            )
        seq_path = Path(source_inputs["sequence"])
        try:
            raw = read_regular_file(seq_path, maximum_bytes=MAX_INPUT_BYTES, label="sequence input")
        except BoundedIOError as failure:
            raise _fail(str(failure)) from failure
        input_ref = source_ir["inputs"]["sequence"]
        if (
            hashlib.sha256(raw).hexdigest() != input_ref["sha256"]
            or len(raw) != input_ref["byte_length"]
        ):
            raise _fail("reference sequence input bytes do not match sealed source descriptor")
        wrapper = source_ir["parameters"]["wrapper"]
        sequences = _decode_fasta_bytes(raw, wrapper)
    elif profile == EXTERNAL_PROFILE:
        if source_inputs is not None and "sequence" in source_inputs:
            seq_path = Path(source_inputs["sequence"])
            try:
                raw = read_regular_file(seq_path, maximum_bytes=MAX_INPUT_BYTES, label="sequence input")
            except BoundedIOError as failure:
                raise _fail(str(failure)) from failure
            sequences = _decode_fasta_bytes(raw, "identity")
        else:
            raise _fail(
                "external-dna-source/v1 variant compilation requires explicit source_inputs={'sequence': fasta_path}"
            )
    else:
        raise _fail(f"unsupported source profile: {profile}")

    if set(sequences) != set(expected_records):
        raise _fail("extracted record identifiers do not match source descriptor")
    for record_id, expected in expected_records.items():
        seq = sequences[record_id]
        if any(ch not in IUPAC_DNA for ch in seq):
            raise _fail(f"record {record_id!r} contains non-IUPAC bases")
        if len(seq) != expected["bases"]:
            raise _fail(f"record {record_id!r} base count mismatch against source descriptor")
        seq_sha = hashlib.sha256(seq.encode("ascii")).hexdigest()
        if seq_sha != expected["sequence_sha256"]:
            raise _fail(f"record {record_id!r} sequence SHA-256 mismatch against source descriptor")
        if _refget_sequence_id(seq) != expected["refget_id"]:
            raise _fail(f"record {record_id!r} refget_id mismatch against source descriptor")
    return sequences


def _edit_kind(ref: str, alt: str) -> str:
    if len(ref) == 0 and len(alt) > 0:
        return "insertion"
    if len(ref) > 0 and len(alt) == 0:
        return "deletion"
    if len(ref) == 1 and len(alt) == 1:
        return "snv"
    if len(ref) == len(alt) and len(ref) > 1:
        return "mnv"
    return "delins"


def _variant_edit_parts(variant: Mapping[str, Any], max_coord: int) -> tuple[int, int, str, str]:
    if not isinstance(variant, Mapping):
        raise _fail("variant must be a mapping")
    loc = variant.get("location") if "location" in variant else variant
    if not isinstance(loc, Mapping) or "start" not in loc or "end" not in loc:
        raise _fail("variant location must provide 'start' and 'end'")
    if "ref" not in variant or "alt" not in variant:
        raise _fail("variant must provide 'ref' and 'alt'")
    start = _nonneg_int(loc["start"], "start", maximum=max_coord)
    end = _nonneg_int(loc["end"], "end", maximum=SAFE_INTEGER)
    if end < start:
        raise _fail("variant end must be >= start")
    ref = _iupac_bases(variant["ref"], "ref", allow_empty=(start == end))
    if start == end and ref != "":
        raise _fail("insertion requires empty ref")
    if start < end and len(ref) != end - start:
        raise _fail("variant ref length must equal end - start")
    alt = _iupac_bases(variant["alt"], "alt", allow_empty=(start < end))
    if start == end and not alt:
        raise _fail("insertion requires nonempty alt")
    if ref == alt:
        raise _fail("no-op edit (ref == alt) is not allowed")
    return start, end, ref, alt


def apply_variant_to_sequence(reference_sequence: str, variant: Mapping[str, Any]) -> str:
    """Apply a verified normalized variant record to a reference DNA string."""
    record_seq = _iupac_bases(reference_sequence, "reference_sequence", allow_empty=False)
    start, end, ref, alt = _variant_edit_parts(variant, len(record_seq))
    if end > len(record_seq):
        raise _fail(f"end must be an integer in [0, {len(record_seq)}]")
    actual_ref = record_seq[start:end]
    if actual_ref != ref:
        raise _fail(f"reference slice {actual_ref!r} at [{start}:{end}] does not match ref {ref!r}")
    return record_seq[:start] + alt + record_seq[end:]


def revert_variant(mutated_sequence: str, variant: Mapping[str, Any]) -> str:
    """Revert a single applied variant back to the exact reference DNA sequence."""
    mut_seq = _iupac_bases(mutated_sequence, "mutated_sequence", allow_empty=True)
    start, _, ref, alt = _variant_edit_parts(variant, len(mut_seq))
    alt_end = start + len(alt)
    if alt_end > len(mut_seq) or mut_seq[start:alt_end] != alt:
        raise _fail("mutated sequence does not contain expected alternate allele at edit coordinates")
    return mut_seq[:start] + ref + mut_seq[alt_end:]


def _source_summary(bundle: SourceBundle) -> dict[str, Any]:
    descriptor = bundle.to_dict()
    source_ir = descriptor["source_ir"]
    return {
        "format": descriptor["format"],
        "version": descriptor["version"],
        "profile": source_ir["profile"],
        "source_ir_sha256": descriptor["source_ir_sha256"],
        "artifact_sha256": descriptor["artifact_sha256"],
        "records": copy.deepcopy(source_ir["records"]),
    }


def _validate_variant_collection_payload(payload: dict[str, Any]) -> dict[str, Any]:
    _keys(
        payload,
        {
            "format",
            "version",
            "vrs_version",
            "producer",
            "source",
            "variants_ir",
            "variants_ir_sha256",
            "artifact_sha256",
        },
        "variant collection",
    )
    if payload["format"] != FORMAT:
        raise _fail(f"unexpected variant collection format: {payload['format']!r}")
    if type(payload["version"]) is not int or payload["version"] != VERSION:
        raise _fail(f"unexpected variant collection version: {payload['version']!r}")
    if payload["vrs_version"] != VRS_VERSION:
        raise _fail(f"unexpected VRS version pin: {payload['vrs_version']!r}")
    if not _exact_equal(payload["producer"], PRODUCER):
        raise _fail("unexpected variant producer metadata")

    source = _keys(
        payload["source"],
        {"format", "version", "profile", "source_ir_sha256", "artifact_sha256", "records"},
        "variant collection source",
    )
    if source["format"] != "brainc.source-descriptor" or type(source["version"]) is not int or source["version"] != 1:
        raise _fail("variant collection source descriptor metadata is invalid")
    if type(source["profile"]) is not str or not source["profile"]:
        raise _fail("variant collection source profile must be nonempty string")
    for key in ("source_ir_sha256", "artifact_sha256"):
        if type(source[key]) is not str or _SHA256_RE.fullmatch(source[key]) is None:
            raise _fail(f"source.{key} must be lowercase SHA-256")
    records = source["records"]
    if type(records) is not list or not records:
        raise _fail("source.records must be a nonempty array")
    record_by_id: dict[str, dict[str, Any]] = {}
    for idx, rec in enumerate(records):
        _keys(rec, {"record_id", "bases", "sequence_sha256", "refget_id"}, f"source.records[{idx}]")
        rid = _nonempty_text(rec["record_id"], f"source.records[{idx}].record_id")
        if rid in record_by_id:
            raise _fail(f"duplicate source record_id: {rid!r}")
        _nonneg_int(rec["bases"], f"source.records[{idx}].bases")
        if rec["bases"] < 1:
            raise _fail(f"source.records[{idx}].bases must be >= 1")
        if type(rec["sequence_sha256"]) is not str or _SHA256_RE.fullmatch(rec["sequence_sha256"]) is None:
            raise _fail(f"source.records[{idx}].sequence_sha256 must be lowercase SHA-256")
        if type(rec["refget_id"]) is not str or _REFGET_RE.fullmatch(rec["refget_id"]) is None:
            raise _fail(f"source.records[{idx}].refget_id must be valid refget SQ identifier")
        record_by_id[rid] = rec

    variants_ir = _keys(payload["variants_ir"], {"variants"}, "variants_ir")
    variants = variants_ir["variants"]
    if type(variants) is not list or not 1 <= len(variants) <= MAX_VARIANTS:
        raise _fail(f"variants_ir.variants must contain 1..{MAX_VARIANTS} entries")

    seen_ids: set[str] = set()
    for ordinal, entry in enumerate(variants):
        label = f"variants_ir.variants[{ordinal}]"
        _keys(
            entry,
            {
                "ordinal",
                "id",
                "record_id",
                "refget_id",
                "sequence_sha256",
                "edit_kind",
                "ref",
                "alt",
                "location",
                "state",
                "vrs_location_id",
                "vrs_allele_id",
                "vrs_allele",
                "mutated_bases",
                "mutated_sequence_sha256",
                "mutated_refget_id",
            },
            label,
        )
        if type(entry["ordinal"]) is not int or entry["ordinal"] != ordinal:
            raise _fail(f"{label}.ordinal must equal {ordinal}")
        vid = _nonempty_text(entry["id"], f"{label}.id")
        if vid in seen_ids:
            raise _fail(f"duplicate variant id: {vid!r}")
        seen_ids.add(vid)

        rid = _nonempty_text(entry["record_id"], f"{label}.record_id")
        if rid not in record_by_id:
            raise _fail(f"{label}.record_id {rid!r} is absent from bound source records")
        rec = record_by_id[rid]
        if entry["refget_id"] != rec["refget_id"]:
            raise _fail(f"{label}.refget_id does not match bound source record")
        if entry["sequence_sha256"] != rec["sequence_sha256"]:
            raise _fail(f"{label}.sequence_sha256 does not match bound source record")

        loc = _keys(
            entry["location"],
            {"type", "sequenceReference", "start", "end"},
            f"{label}.location",
        )
        if loc["type"] != "SequenceLocation":
            raise _fail(f"{label}.location.type must be SequenceLocation")
        seq_ref = _keys(
            loc["sequenceReference"],
            {"type", "refgetAccession"},
            f"{label}.location.sequenceReference",
        )
        if seq_ref["type"] != "SequenceReference" or seq_ref["refgetAccession"] != rec["refget_id"]:
            raise _fail(f"{label}.location.sequenceReference mismatch")
        start = _nonneg_int(loc["start"], f"{label}.location.start", maximum=rec["bases"])
        end = _nonneg_int(loc["end"], f"{label}.location.end", maximum=rec["bases"])
        if end < start:
            raise _fail(f"{label}.location.end must be >= start")

        ref = _iupac_bases(entry["ref"], f"{label}.ref", allow_empty=(start == end))
        if start == end and ref != "":
            raise _fail(f"{label}.ref must be empty when start == end (insertion)")
        if start < end and len(ref) != end - start:
            raise _fail(f"{label}.ref length must equal end - start")
        alt = _iupac_bases(entry["alt"], f"{label}.alt", allow_empty=(start < end))
        if start == end and not alt:
            raise _fail(f"{label}.alt must be nonempty for insertion")
        if ref == alt:
            raise _fail(f"{label} ref and alt cannot be identical")
        if entry["edit_kind"] != _edit_kind(ref, alt):
            raise _fail(f"{label}.edit_kind does not match ref/alt lengths")

        state = _keys(entry["state"], {"type", "sequence"}, f"{label}.state")
        if state["type"] != "LiteralSequenceExpression" or state["sequence"] != alt:
            raise _fail(f"{label}.state mismatch")

        expected_allele = vrs_allele(rec["refget_id"], start, end, alt)
        if entry["vrs_location_id"] != expected_allele["location"]["id"]:
            raise _fail(f"{label}.vrs_location_id mismatch")
        if entry["vrs_allele_id"] != expected_allele["id"]:
            raise _fail(f"{label}.vrs_allele_id mismatch")
        if not _exact_equal(entry["vrs_allele"], expected_allele):
            raise _fail(f"{label}.vrs_allele mismatch")

        expected_mut_bases = rec["bases"] - (end - start) + len(alt)
        if type(entry["mutated_bases"]) is not int or entry["mutated_bases"] != expected_mut_bases:
            raise _fail(f"{label}.mutated_bases mismatch")
        if type(entry["mutated_sequence_sha256"]) is not str or _SHA256_RE.fullmatch(entry["mutated_sequence_sha256"]) is None:
            raise _fail(f"{label}.mutated_sequence_sha256 must be lowercase SHA-256")
        if type(entry["mutated_refget_id"]) is not str or _REFGET_RE.fullmatch(entry["mutated_refget_id"]) is None:
            raise _fail(f"{label}.mutated_refget_id must be valid SQ identifier")

    if payload["variants_ir_sha256"] != digest(variants_ir):
        raise _fail("variants_ir_sha256 does not match RFC 8785 canonical digest")
    core = {key: value for key, value in payload.items() if key != "artifact_sha256"}
    if payload["artifact_sha256"] != digest(core):
        raise _fail("artifact_sha256 does not match RFC 8785 canonical digest")
    return payload


class VariantCollectionArtifact:
    """Immutable-by-interface sealed variant collection artifact."""

    __slots__ = ("_payload",)

    def __init__(self, *_: Any, **__: Any) -> None:
        raise TypeError("VariantCollectionArtifact is constructed by compile_variants or load_variant_collection")

    @classmethod
    def _create(cls, payload: dict[str, Any]) -> VariantCollectionArtifact:
        instance = object.__new__(cls)
        instance._payload = copy.deepcopy(payload)
        return instance

    @property
    def artifact_sha256(self) -> str:
        return self._payload["artifact_sha256"]

    @property
    def variants_ir_sha256(self) -> str:
        return self._payload["variants_ir_sha256"]

    def to_dict(self) -> dict[str, Any]:
        return copy.deepcopy(self._payload)

    def save(self, path: str | Path) -> Path:
        try:
            pretty_json_bytes(
                self._payload,
                ensure_ascii=False,
                maximum_bytes=MAX_JSON_BYTES,
                maximum_depth=MAX_JSON_DEPTH,
                maximum_members=MAX_JSON_MEMBERS,
                maximum_string_bytes=MAX_STRING_BYTES,
            )
            return save_artifact(self._payload, path)
        except (BoundedIOError, ContractError, OSError) as failure:
            raise _fail(f"cannot save variant collection artifact: {failure}") from failure


def compile_variants(
    source: SourceBundle | str | Path,
    variants: Sequence[Mapping[str, Any]],
    *,
    source_inputs: Mapping[str, str | Path] | None = None,
) -> VariantCollectionArtifact:
    """Compile and seal edits against a validated SourceBundle and exact reference bases."""
    if isinstance(source, (str, Path)):
        try:
            bundle = load_source_bundle(source)
        except SourceError as failure:
            raise _fail(f"cannot load source bundle: {failure}") from failure
    elif isinstance(source, SourceBundle):
        bundle = source
    else:
        raise _fail("source must be a SourceBundle or directory path")

    if not isinstance(variants, Sequence) or isinstance(variants, (str, bytes, bytearray)):
        raise _fail("variants must be a nonempty array of variant objects")
    if not 1 <= len(variants) <= MAX_VARIANTS:
        raise _fail(f"variants must contain 1..{MAX_VARIANTS} entries")

    sequences = _extract_source_sequences(bundle, source_inputs=source_inputs)
    source_meta = _source_summary(bundle)
    record_by_id = {item["record_id"]: item for item in source_meta["records"]}

    seen_ids: set[str] = set()
    compiled_entries: list[dict[str, Any]] = []
    for ordinal, raw_variant in enumerate(variants):
        label = f"variants[{ordinal}]"
        if not isinstance(raw_variant, Mapping):
            raise _fail(f"{label} must be a JSON object")
        item = dict(raw_variant)
        _keys(item, {"id", "record_id", "start", "end", "ref", "alt"}, label)
        vid = _nonempty_text(item["id"], f"{label}.id")
        if vid in seen_ids:
            raise _fail(f"duplicate variant id: {vid!r}")
        seen_ids.add(vid)

        rid = _nonempty_text(item["record_id"], f"{label}.record_id")
        if rid not in record_by_id:
            raise _fail(f"{label}.record_id {rid!r} is absent from bound source records")
        rec = record_by_id[rid]
        ref_seq = sequences[rid]

        start = _nonneg_int(item["start"], f"{label}.start", maximum=rec["bases"])
        end = _nonneg_int(item["end"], f"{label}.end", maximum=rec["bases"])
        if end < start:
            raise _fail(f"{label}.end must be >= start")

        ref = _iupac_bases(item["ref"], f"{label}.ref", allow_empty=(start == end))
        if start == end and ref != "":
            raise _fail(f"{label}.ref must be empty when start == end (insertion)")
        if start < end and len(ref) != end - start:
            raise _fail(f"{label}.ref length must equal end - start")
        actual_ref = ref_seq[start:end]
        if actual_ref != ref:
            raise _fail(
                f"{label} reference mismatch at {rid}:{start}-{end}: expected {actual_ref!r}, got {ref!r}"
            )

        alt = _iupac_bases(item["alt"], f"{label}.alt", allow_empty=(start < end))
        if start == end and not alt:
            raise _fail(f"{label}.alt must be nonempty for insertion")
        if ref == alt:
            raise _fail(f"{label} ref and alt cannot be identical")

        mutated_seq = ref_seq[:start] + alt + ref_seq[end:]
        mutated_bytes = mutated_seq.encode("ascii")
        mutated_sha256 = hashlib.sha256(mutated_bytes).hexdigest()
        mutated_refget = "SQ." + _sha512t24u(mutated_bytes)

        allele_obj = vrs_allele(rec["refget_id"], start, end, alt)
        location_obj = {
            "type": "SequenceLocation",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": rec["refget_id"],
            },
            "start": start,
            "end": end,
        }
        state_obj = {
            "type": "LiteralSequenceExpression",
            "sequence": alt,
        }
        compiled_entries.append(
            {
                "ordinal": ordinal,
                "id": vid,
                "record_id": rid,
                "refget_id": rec["refget_id"],
                "sequence_sha256": rec["sequence_sha256"],
                "edit_kind": _edit_kind(ref, alt),
                "ref": ref,
                "alt": alt,
                "location": location_obj,
                "state": state_obj,
                "vrs_location_id": allele_obj["location"]["id"],
                "vrs_allele_id": allele_obj["id"],
                "vrs_allele": allele_obj,
                "mutated_bases": len(mutated_seq),
                "mutated_sequence_sha256": mutated_sha256,
                "mutated_refget_id": mutated_refget,
            }
        )

    variants_ir = {"variants": compiled_entries}
    core = {
        "format": FORMAT,
        "version": VERSION,
        "vrs_version": VRS_VERSION,
        "producer": copy.deepcopy(PRODUCER),
        "source": source_meta,
        "variants_ir": variants_ir,
        "variants_ir_sha256": digest(variants_ir),
    }
    payload = {**core, "artifact_sha256": digest(core)}
    _validate_variant_collection_payload(payload)
    return VariantCollectionArtifact._create(payload)


def load_variant_collection(path: str | Path) -> VariantCollectionArtifact:
    """Load and structurally validate a sealed brainc.variant-collection/v1 JSON artifact."""
    try:
        raw = read_regular_file(path, maximum_bytes=MAX_JSON_BYTES, label="variant collection")
        payload = load_json_object(
            raw,
            "variant collection",
            maximum_bytes=MAX_JSON_BYTES,
            maximum_depth=MAX_JSON_DEPTH,
            maximum_members=MAX_JSON_MEMBERS,
            maximum_string_bytes=MAX_STRING_BYTES,
        )
    except (BoundedIOError, ContractError) as failure:
        raise _fail(f"cannot load variant collection: {failure}") from failure
    _validate_variant_collection_payload(payload)
    return VariantCollectionArtifact._create(payload)


def load_variant_inputs(path: str | Path) -> list[dict[str, Any]]:
    """Load an input variant list (JSON array or {'variants': [...]}) safely."""
    try:
        raw = read_regular_file(path, maximum_bytes=MAX_JSON_BYTES, label="variant input list")
        text = raw.decode("utf-8")
    except (BoundedIOError, UnicodeDecodeError) as failure:
        raise _fail(f"cannot read variant input list: {failure}") from failure

    import json

    def _no_dups(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        obj: dict[str, Any] = {}
        for k, v in pairs:
            if k in obj:
                raise _fail(f"duplicate JSON key in variant input list: {k!r}")
            obj[k] = v
        return obj

    def _no_const(val: str) -> None:
        raise _fail(f"non-finite JSON constant in variant input list: {val}")

    try:
        parsed = json.loads(text, object_pairs_hook=_no_dups, parse_constant=_no_const)
    except (ValueError, RecursionError) as failure:
        raise _fail(f"invalid variant input JSON: {failure}") from failure

    if type(parsed) is dict:
        _keys(parsed, {"variants"}, "variant input wrapper")
        parsed = parsed["variants"]
    if type(parsed) is not list:
        raise _fail("variant input file must contain a JSON array or {'variants': [...]}")
    return parsed


__all__ = [
    "FORMAT",
    "IUPAC_DNA",
    "PRODUCER",
    "VERSION",
    "VRS_VERSION",
    "VariantCollectionArtifact",
    "VariantError",
    "apply_variant_to_sequence",
    "compile_variants",
    "load_variant_collection",
    "load_variant_inputs",
    "revert_variant",
    "vrs_allele",
    "vrs_allele_digest",
    "vrs_location_digest",
    "vrs_sequence_location",
]

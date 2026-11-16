"""Producer-isolated replay validator for brainc.variant-collection/v1.

This module uses Python's standard library exclusively. It deliberately imports
no producer canonicalization, I/O, source, or variant compilation modules.
"""

from __future__ import annotations

import argparse
import base64
from collections.abc import Mapping
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
from typing import Any
import zlib


SAFE_INTEGER = 2**53 - 1
MAX_JSON_BYTES = 16 * 1024 * 1024
MAX_SOURCE_DESCRIPTOR_BYTES = 256 * 1024 * 1024
MAX_GENBANK_ARTIFACT_BYTES = 64 * 1024 * 1024
MAX_SEQUENCE_INPUT_BYTES = 16 * 1024 * 1024
MAX_DECOMPRESSED_BYTES = 64 * 1024 * 1024
MAX_JSON_DEPTH = 64
MAX_JSON_MEMBERS = 1_000_000
MAX_STRING_BYTES = 1 * 1024 * 1024
MAX_IDENTIFIER_BYTES = 256
MAX_VARIANTS = 100_000
MAX_VARIANT_BASES = 1_000_000

VARIANT_FORMAT = "brainc.variant-collection"
VARIANT_VERSION = 1
VRS_VERSION = "2.x"
SOURCE_FORMAT = "brainc.source-descriptor"
SOURCE_VERSION = 1
REPORT_FORMAT = "brainc.variant-validation-report"
REPORT_VERSION = 1

RAW_PROFILE = "raw-iupac-dna/v1"
FASTA_PROFILE = "fasta-dna/v1"
REFERENCE_FASTA_PROFILE = "fasta-reference-dna/v1"
GENBANK_PROFILE = "genbank-273-traditional-dna-physical-structural/v2"
GFF3_PROFILE = "gff3-external-sequence/v1"
EXTERNAL_PROFILE = "external-dna-source/v1"

EXPECTED_SOURCE_PRODUCER = {
    "name": "brainc-source",
    "version": "1.0.0",
    "passes": [
        "select-source-profile",
        "compile-native-source",
        "bind-source-evidence",
        "emit-source-descriptor",
    ],
}
EXPECTED_PRODUCER = {
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
VALIDATOR_INFO = {
    "name": "brainc-validator-variant",
    "version": "1.0.0",
}

IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")
_SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
_REFGET_RE = re.compile(r"SQ\.[A-Za-z0-9_-]{32}\Z")
_GZIP_HEADER = b"\x1f\x8b"

_HAS_DIR_FD = (
    os.name == "posix"
    and os.listdir in getattr(os, "supports_fd", set())
    and all(op in getattr(os, "supports_dir_fd", set()) for op in (os.link, os.open, os.stat, os.unlink))
    and getattr(os, "O_DIRECTORY", 0) != 0
    and getattr(os, "O_NOFOLLOW", 0) != 0
)


class VariantValidationError(ValueError):
    """Independent replay validation failure (`VARVAL001`)."""


def _fail(detail: str) -> VariantValidationError:
    return VariantValidationError(f"VARVAL001: {detail}")


def _sha512t24u(raw: bytes) -> str:
    return base64.urlsafe_b64encode(hashlib.sha512(raw).digest()[:24]).decode("ascii").rstrip("=")


def _jcs_string(value: str) -> str:
    if any(0xD800 <= ord(ch) <= 0xDFFF for ch in value):
        raise _fail("lone Unicode surrogate in JSON string")
    return json.dumps(value, ensure_ascii=False, allow_nan=False)


def _jcs_number(value: float) -> str:
    if not math.isfinite(value):
        raise _fail("non-finite number in JSON")
    if value == 0:
        return "0"
    sign = "-" if value < 0 else ""
    mantissa, sep, exp_text = repr(abs(value)).lower().partition("e")
    exponent = int(exp_text) if sep else 0
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
    norm_exp = point - 1
    rendered = digits[0] + (("." + digits[1:]) if len(digits) > 1 else "")
    return sign + rendered + ("e+" if norm_exp >= 0 else "e") + str(norm_exp)


def _canonical_bytes(value: Any, *, _depth: int = 0) -> bytes:
    if _depth > MAX_JSON_DEPTH:
        raise _fail("JSON exceeds maximum depth during canonicalization")
    if value is None:
        return b"null"
    if type(value) is bool:
        return b"true" if value else b"false"
    if type(value) is str:
        return _jcs_string(value).encode("utf-8")
    if type(value) is int:
        if not -SAFE_INTEGER <= value <= SAFE_INTEGER:
            raise _fail(f"integer outside safe range: {value}")
        return str(value).encode("ascii")
    if type(value) is float:
        return _jcs_number(value).encode("ascii")
    if type(value) is list:
        return b"[" + b",".join(_canonical_bytes(v, _depth=_depth + 1) for v in value) + b"]"
    if type(value) is dict:
        parts = []
        for k in sorted(value, key=lambda key: key.encode("utf-16be")):
            if type(k) is not str:
                raise _fail("JSON object keys must be strings")
            parts.append(_jcs_string(k).encode("utf-8") + b":" + _canonical_bytes(value[k], _depth=_depth + 1))
        return b"{" + b",".join(parts) + b"}"
    raise _fail(f"unsupported JSON type: {type(value).__name__}")


def _digest(value: Any) -> str:
    return hashlib.sha256(_canonical_bytes(value)).hexdigest()


def _utf8_bytes(value: str, label: str) -> bytes:
    try:
        return value.encode("utf-8")
    except UnicodeEncodeError as failure:
        raise _fail(f"{label} contains a lone Unicode surrogate") from failure


def _check_json_tree(value: Any, *, depth: int = 0, members: list[int] | None = None) -> int:
    if members is None:
        members = [0]
    if depth > MAX_JSON_DEPTH:
        raise _fail("JSON exceeds maximum nesting depth")
    members[0] += 1
    if members[0] > MAX_JSON_MEMBERS:
        raise _fail("JSON exceeds maximum member count")
    if value is None or type(value) is bool:
        return members[0]
    if type(value) is int:
        if not -SAFE_INTEGER <= value <= SAFE_INTEGER:
            raise _fail(f"integer outside safe range: {value}")
        return members[0]
    if type(value) is float:
        if not math.isfinite(value):
            raise _fail("non-finite float in JSON")
        return members[0]
    if type(value) is str:
        if len(_utf8_bytes(value, "JSON string")) > MAX_STRING_BYTES:
            raise _fail("JSON string exceeds maximum byte limit")
        return members[0]
    if type(value) is list:
        for item in value:
            _check_json_tree(item, depth=depth + 1, members=members)
        return members[0]
    if type(value) is dict:
        for k, v in value.items():
            if type(k) is not str:
                raise _fail("JSON object keys must be strings")
            if len(_utf8_bytes(k, "JSON key")) > MAX_STRING_BYTES:
                raise _fail("JSON key exceeds maximum byte limit")
            _check_json_tree(v, depth=depth + 1, members=members)
        return members[0]
    raise _fail(f"unsupported JSON value type: {type(value).__name__}")


def _parse_json_object(raw: bytes, label: str) -> dict[str, Any]:
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as failure:
        raise _fail(f"{label} is not valid UTF-8: {failure}") from failure

    def _no_dups(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        out: dict[str, Any] = {}
        for k, v in pairs:
            if k in out:
                raise _fail(f"duplicate JSON key in {label}: {k!r}")
            out[k] = v
        return out

    def _no_const(val: str) -> None:
        raise _fail(f"non-finite JSON constant in {label}: {val}")

    try:
        parsed = json.loads(text, object_pairs_hook=_no_dups, parse_constant=_no_const)
    except VariantValidationError:
        raise
    except (ValueError, RecursionError) as failure:
        raise _fail(f"malformed JSON in {label}: {failure}") from failure
    if type(parsed) is not dict:
        raise _fail(f"{label} must be a top-level JSON object")
    _check_json_tree(parsed)
    return parsed


def _stat_identity(info: os.stat_result) -> tuple[Any, ...]:
    return (
        info.st_dev,
        info.st_ino,
        info.st_size,
        info.st_nlink,
        getattr(info, "st_mtime_ns", int(info.st_mtime * 1e9)),
        getattr(info, "st_ctime_ns", int(info.st_ctime * 1e9)),
    )


def _read_fd_bounded(fd: int, label: str, *, maximum_bytes: int) -> bytes:
    info = os.fstat(fd)
    if not stat.S_ISREG(info.st_mode):
        raise _fail(f"{label} is not a regular file")
    if os.name == "posix" and info.st_nlink != 1:
        raise _fail(f"{label} must have link count 1 (found {info.st_nlink})")
    if info.st_size > maximum_bytes:
        raise _fail(f"{label} exceeds maximum byte size {maximum_bytes}")
    chunks: list[bytes] = []
    total = 0
    while True:
        block = os.read(fd, 65536)
        if not block:
            break
        total += len(block)
        if total > maximum_bytes:
            raise _fail(f"{label} exceeded maximum byte size during read")
        chunks.append(block)
    after = os.fstat(fd)
    if _stat_identity(after) != _stat_identity(info) or total != info.st_size:
        raise _fail(f"{label} changed during read")
    return b"".join(chunks)


def _read_regular_file(path: str | Path, label: str, *, maximum_bytes: int) -> bytes:
    target = Path(path)
    if not str(target):
        raise _fail(f"{label} path is empty")
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
    try:
        fd = os.open(target, flags)
    except OSError as failure:
        raise _fail(f"cannot open {label}: {failure}") from failure
    try:
        return _read_fd_bounded(fd, label, maximum_bytes=maximum_bytes)
    finally:
        os.close(fd)


def _read_source_bundle_files(directory: str | Path) -> dict[str, bytes]:
    bundle_dir = Path(directory)
    role_limits = {
        "source.json": MAX_SOURCE_DESCRIPTOR_BYTES,
        "sequence.json": MAX_JSON_BYTES,
        "reference.json": MAX_SOURCE_DESCRIPTOR_BYTES,
        "external.json": MAX_SOURCE_DESCRIPTOR_BYTES,
        "genbank.json": MAX_GENBANK_ARTIFACT_BYTES,
        "annotation.json": MAX_SOURCE_DESCRIPTOR_BYTES,
    }
    if not _HAS_DIR_FD:
        if bundle_dir.is_symlink() or not bundle_dir.is_dir():
            raise _fail("source bundle must be a non-symlink directory")
        names = [entry.name for entry in bundle_dir.iterdir()]
        if "source.json" not in names or not set(names) <= set(role_limits):
            raise _fail(f"unexpected children in source bundle directory: {sorted(names)}")
        return {name: _read_regular_file(bundle_dir / name, name, maximum_bytes=role_limits[name]) for name in names}

    dir_fd = -1
    try:
        dir_fd = os.open(bundle_dir, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
        names = os.listdir(dir_fd)
        if "source.json" not in names or not set(names) <= set(role_limits):
            raise _fail(f"unexpected children in source bundle directory: {sorted(names)}")
        out: dict[str, bytes] = {}
        child_flags = os.O_RDONLY | os.O_NOFOLLOW | getattr(os, "O_NONBLOCK", 0)
        for name in names:
            cfd = os.open(name, child_flags, dir_fd=dir_fd)
            try:
                out[name] = _read_fd_bounded(cfd, f"source bundle child {name!r}", maximum_bytes=role_limits[name])
            finally:
                os.close(cfd)
        return out
    except VariantValidationError:
        raise
    except OSError as failure:
        raise _fail(f"cannot read source bundle directory: {failure}") from failure
    finally:
        if dir_fd >= 0:
            os.close(dir_fd)


def _keys(obj: Any, expected: set[str], label: str) -> dict[str, Any]:
    if type(obj) is not dict:
        raise _fail(f"{label} must be an object")
    missing = sorted(expected - obj.keys())
    extra = sorted(obj.keys() - expected)
    if missing or extra:
        raise _fail(f"{label} keys invalid; missing={missing or 'none'}, unknown={extra or 'none'}")
    return obj


def _exact_equal(left: Any, right: Any) -> bool:
    if type(left) is not type(right):
        return False
    if type(left) is dict:
        return left.keys() == right.keys() and all(_exact_equal(left[k], right[k]) for k in left)
    if type(left) is list:
        return len(left) == len(right) and all(_exact_equal(a, b) for a, b in zip(left, right))
    return left == right


def _verify_bundle_child(
    bundle_files: dict[str, bytes],
    artifacts_meta: Mapping[str, Any],
    role: str,
    filename: str,
    *,
    expected_format: str,
    expected_version: int,
    ir_key: str,
) -> dict[str, Any]:
    if filename not in bundle_files:
        raise _fail(f"missing {filename} in source bundle")
    obj = _parse_json_object(bundle_files[filename], filename)
    meta = artifacts_meta.get(role)
    if type(meta) is not dict:
        raise _fail(f"source_ir.artifacts missing {role!r} entry")
    _keys(meta, {"format", "version", "artifact_sha256", "ir_sha256"}, f"source_ir.artifacts.{role}")
    if obj.get("format") != expected_format or type(obj.get("version")) is not int or obj.get("version") != expected_version:
        raise _fail(f"{filename} format or version invalid")
    if meta["format"] != expected_format or type(meta["version"]) is not int or meta["version"] != expected_version:
        raise _fail(f"source_ir.artifacts.{role} format or version mismatch")
    core = {k: v for k, v in obj.items() if k != "artifact_sha256"}
    if obj.get("artifact_sha256") != _digest(core) or obj.get("artifact_sha256") != meta["artifact_sha256"]:
        raise _fail(f"{filename} artifact_sha256 disagrees with source descriptor")
    if ir_key not in obj or obj[ir_key] != meta["ir_sha256"]:
        raise _fail(f"{filename} {ir_key} disagrees with source descriptor")
    return obj


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
                    raise _fail("decompressed reference sequence exceeds ceiling")
                consumed = len(chunk) - len(decomp.unconsumed_tail)
                offset += consumed
                if consumed == 0 and not decomp.eof:
                    raise _fail("malformed gzip reference sequence")
            tail = decomp.flush()
            if len(out) + len(tail) > MAX_DECOMPRESSED_BYTES:
                raise _fail("decompressed reference sequence exceeds ceiling")
            out.extend(tail)
            if not decomp.eof:
                raise _fail("truncated gzip reference sequence")
        except zlib.error as failure:
            raise _fail(f"malformed gzip reference sequence: {failure}") from failure
    return bytes(out)


def _parse_fasta_independent(raw: bytes, wrapper: str) -> dict[str, str]:
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
    buf: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if current_id is not None:
                records[current_id] = "".join(buf)
            hdr = line[1:].strip()
            if not hdr:
                raise _fail("empty FASTA header")
            current_id = hdr.split()[0]
            if current_id in records:
                raise _fail(f"duplicate FASTA record id: {current_id!r}")
            buf = []
        else:
            s = line.strip().upper()
            if not s:
                continue
            if current_id is None:
                raise _fail("FASTA bases before header")
            if any(c not in IUPAC_DNA for c in s):
                raise _fail(f"FASTA record {current_id!r} contains non-IUPAC bases")
            buf.append(s)
    if current_id is not None:
        records[current_id] = "".join(buf)
    if not records:
        raise _fail("FASTA contains no records")
    return records


def _extract_reference_sequences(
    bundle_files: dict[str, bytes],
    source_inputs: Mapping[str, str | Path] | None = None,
) -> tuple[dict[str, Any], dict[str, str]]:
    source_desc = _parse_json_object(bundle_files["source.json"], "source.json")
    _keys(
        source_desc,
        {"format", "version", "producer", "source_ir", "source_ir_sha256", "artifact_sha256"},
        "source.json",
    )
    if source_desc["format"] != SOURCE_FORMAT or type(source_desc["version"]) is not int or source_desc["version"] != SOURCE_VERSION:
        raise _fail("unsupported source descriptor format or version")
    if not _exact_equal(source_desc["producer"], EXPECTED_SOURCE_PRODUCER):
        raise _fail("unexpected source.json producer object")
    if source_desc["source_ir_sha256"] != _digest(source_desc["source_ir"]):
        raise _fail("source.json source_ir_sha256 mismatch")
    source_core = {k: v for k, v in source_desc.items() if k != "artifact_sha256"}
    if source_desc["artifact_sha256"] != _digest(source_core):
        raise _fail("source.json artifact_sha256 mismatch")

    source_ir = _keys(source_desc["source_ir"], {"profile", "parameters", "inputs", "artifacts", "records"}, "source_ir")
    profile = source_ir["profile"]
    if type(profile) is not str:
        raise _fail("source_ir.profile must be a string")
    if type(source_ir["parameters"]) is not dict or type(source_ir["inputs"]) is not dict:
        raise _fail("source_ir.parameters and source_ir.inputs must be objects")
    records_list = source_ir["records"]
    if type(records_list) is not list or not records_list:
        raise _fail("source_ir.records must be nonempty")
    expected_by_id: dict[str, dict[str, Any]] = {}
    for idx, rec in enumerate(records_list):
        _keys(rec, {"record_id", "bases", "sequence_sha256", "refget_id"}, f"source_ir.records[{idx}]")
        rid = rec["record_id"]
        if type(rid) is not str or not rid:
            raise _fail(f"source_ir.records[{idx}].record_id invalid")
        if rid in expected_by_id:
            raise _fail(f"duplicate record_id {rid!r} in source_ir.records")
        if type(rec["bases"]) is not int or rec["bases"] < 1:
            raise _fail(f"source_ir.records[{idx}].bases must be positive integer")
        if type(rec["sequence_sha256"]) is not str or type(rec["refget_id"]) is not str:
            raise _fail(f"source_ir.records[{idx}] digest fields must be strings")
        expected_by_id[rid] = rec

    artifacts_meta = source_ir["artifacts"]
    if type(artifacts_meta) is not dict:
        raise _fail("source_ir.artifacts must be an object")

    expected_file_closure: dict[str, set[str]] = {
        RAW_PROFILE: {"source.json", "sequence.json"},
        FASTA_PROFILE: {"source.json", "sequence.json"},
        GFF3_PROFILE: {"source.json", "sequence.json", "annotation.json"},
        GENBANK_PROFILE: {"source.json", "genbank.json"},
        REFERENCE_FASTA_PROFILE: {"source.json", "reference.json"},
        EXTERNAL_PROFILE: {"source.json", "external.json"},
    }
    if profile not in expected_file_closure:
        raise _fail(f"unsupported source profile: {profile!r}")
    if set(bundle_files) != expected_file_closure[profile]:
        raise _fail(f"source bundle files {sorted(bundle_files)} do not match profile {profile!r}")

    sequences: dict[str, str] = {}
    if profile in {RAW_PROFILE, FASTA_PROFILE, GFF3_PROFILE}:
        seq_obj = _verify_bundle_child(
            bundle_files,
            artifacts_meta,
            "sequence",
            "sequence.json",
            expected_format="brain01.sequence-collection-ir",
            expected_version=1,
            ir_key="collection_ir_sha256",
        )
        if profile == GFF3_PROFILE:
            _verify_bundle_child(
                bundle_files,
                artifacts_meta,
                "annotation",
                "annotation.json",
                expected_format="brainc.bio.gff3-ir",
                expected_version=1,
                ir_key="bio_ir_sha256",
            )
            if set(artifacts_meta) != {"sequence", "annotation"}:
                raise _fail("source_ir.artifacts roles mismatch for GFF3 profile")
        elif set(artifacts_meta) != {"sequence"}:
            raise _fail("source_ir.artifacts roles mismatch for sequence profile")
        members = seq_obj.get("members")
        if type(members) is not list:
            raise _fail("sequence.json.members must be an array")
        for idx, member in enumerate(members):
            if type(member) is not dict or type(member.get("sequence_artifact")) is not dict:
                raise _fail(f"sequence.json.members[{idx}] missing sequence_artifact")
            seq_ir = member["sequence_artifact"].get("sequence_ir")
            if type(seq_ir) is not dict or type(seq_ir.get("record_id")) is not str or type(seq_ir.get("sequence")) is not str:
                raise _fail(f"sequence.json.members[{idx}] missing sequence_ir record_id/sequence")
            rid = seq_ir["record_id"]
            if rid in sequences:
                raise _fail(f"duplicate record_id {rid!r} in sequence.json")
            sequences[rid] = seq_ir["sequence"].upper()
    elif profile == GENBANK_PROFILE:
        gb_obj = _verify_bundle_child(
            bundle_files,
            artifacts_meta,
            "genbank",
            "genbank.json",
            expected_format="brainc.bio.insdc-genbank-ir",
            expected_version=2,
            ir_key="bio_ir_sha256",
        )
        if set(artifacts_meta) != {"genbank"}:
            raise _fail("source_ir.artifacts roles mismatch for genbank profile")
        sc = gb_obj.get("sequence_collection")
        if type(sc) is not dict or type(sc.get("members")) is not list:
            raise _fail("genbank.json missing sequence_collection.members")
        for idx, member in enumerate(sc["members"]):
            if type(member) is not dict or type(member.get("record_id")) is not str:
                raise _fail(f"genbank.json.members[{idx}] missing record_id")
            seq_obj_m = member.get("sequence")
            storage = seq_obj_m.get("storage") if type(seq_obj_m) is dict else None
            chunks = storage.get("chunks") if type(storage) is dict else None
            if type(chunks) is not list or any(type(c) is not str for c in chunks):
                raise _fail(f"genbank.json.members[{idx}] invalid sequence chunks")
            rid = member["record_id"]
            if rid in sequences:
                raise _fail(f"duplicate record_id {rid!r} in genbank.json")
            sequences[rid] = "".join(chunks).upper()
    elif profile == REFERENCE_FASTA_PROFILE:
        _verify_bundle_child(
            bundle_files,
            artifacts_meta,
            "reference",
            "reference.json",
            expected_format="brainc.reference-sequence-catalog",
            expected_version=1,
            ir_key="catalog_sha256",
        )
        if set(artifacts_meta) != {"reference"}:
            raise _fail("source_ir.artifacts roles mismatch for reference profile")
        if source_inputs is None or "sequence" not in source_inputs:
            raise _fail("fasta-reference-dna/v1 validation requires --source-input sequence=PATH")
        raw = _read_regular_file(source_inputs["sequence"], "sequence input", maximum_bytes=MAX_SEQUENCE_INPUT_BYTES)
        inp_ref = source_ir["inputs"].get("sequence")
        if type(inp_ref) is not dict or type(inp_ref.get("sha256")) is not str or type(inp_ref.get("byte_length")) is not int:
            raise _fail("source_ir.inputs.sequence invalid")
        if hashlib.sha256(raw).hexdigest() != inp_ref["sha256"] or len(raw) != inp_ref["byte_length"]:
            raise _fail("sequence input SHA-256 or byte_length disagrees with source.json")
        wrapper = source_ir["parameters"].get("wrapper")
        if type(wrapper) is not str:
            raise _fail("source_ir.parameters.wrapper must be a string")
        sequences = _parse_fasta_independent(raw, wrapper)
    elif profile == EXTERNAL_PROFILE:
        _verify_bundle_child(
            bundle_files,
            artifacts_meta,
            "external",
            "external.json",
            expected_format="brainc.external-source-closure",
            expected_version=1,
            ir_key="closure_ir_sha256",
        )
        if set(artifacts_meta) != {"external"}:
            raise _fail("source_ir.artifacts roles mismatch for external profile")
        if source_inputs is None or "sequence" not in source_inputs:
            raise _fail("external-dna-source/v1 variant validation requires --source-input sequence=PATH")
        raw = _read_regular_file(source_inputs["sequence"], "sequence input", maximum_bytes=MAX_SEQUENCE_INPUT_BYTES)
        sequences = _parse_fasta_independent(raw, "identity")

    if set(sequences) != set(expected_by_id):
        raise _fail("reconstructed record ids do not match source_ir.records")
    for rid, rec in expected_by_id.items():
        seq = sequences[rid]
        if any(c not in IUPAC_DNA for c in seq):
            raise _fail(f"record {rid!r} contains non-IUPAC bases")
        try:
            seq_bytes = seq.encode("ascii")
        except UnicodeEncodeError as failure:
            raise _fail(f"record {rid!r} must be ASCII DNA") from failure
        if len(seq) != rec["bases"]:
            raise _fail(f"record {rid!r} length mismatch")
        if hashlib.sha256(seq_bytes).hexdigest() != rec["sequence_sha256"]:
            raise _fail(f"record {rid!r} sequence_sha256 mismatch")
        if "SQ." + _sha512t24u(seq_bytes) != rec["refget_id"]:
            raise _fail(f"record {rid!r} refget_id mismatch")
    return source_desc, sequences


def _independent_vrs_location_digest(refget_id: str, start: int, end: int) -> str:
    obj = {
        "end": end,
        "sequenceReference": {
            "refgetAccession": refget_id,
            "type": "SequenceReference",
        },
        "start": start,
        "type": "SequenceLocation",
    }
    return _sha512t24u(_canonical_bytes(obj))


def _independent_vrs_allele_digest(loc_digest: str, alt: str) -> str:
    obj = {
        "location": loc_digest,
        "state": {
            "sequence": alt,
            "type": "LiteralSequenceExpression",
        },
        "type": "Allele",
    }
    return _sha512t24u(_canonical_bytes(obj))


def _independent_edit_kind(ref: str, alt: str) -> str:
    if len(ref) == 0 and len(alt) > 0:
        return "insertion"
    if len(ref) > 0 and len(alt) == 0:
        return "deletion"
    if len(ref) == 1 and len(alt) == 1:
        return "snv"
    if len(ref) == len(alt) and len(ref) > 1:
        return "mnv"
    return "delins"


def _seal_report(core: dict[str, Any]) -> dict[str, Any]:
    return {**core, "report_sha256": _digest(core)}


def validate_variant_collection(
    bundle_files: dict[str, bytes],
    variant_payload: dict[str, Any],
    *,
    source_inputs: Mapping[str, str | Path] | None = None,
) -> dict[str, Any]:
    """Independently replay and verify a brainc.variant-collection/v1 artifact."""
    checks: list[dict[str, Any]] = []
    source_artifact_sha: str | None = None
    variant_artifact_sha: str | None = variant_payload.get("artifact_sha256") if isinstance(variant_payload, dict) else None
    try:
        source_desc, sequences = _extract_reference_sequences(bundle_files, source_inputs=source_inputs)
        source_artifact_sha = source_desc["artifact_sha256"]
        checks.append({"name": "source-bundle-integrity-and-sequences", "passed": True})

        _keys(
            variant_payload,
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
        if variant_payload["format"] != VARIANT_FORMAT:
            raise _fail(f"invalid variant collection format: {variant_payload['format']!r}")
        if type(variant_payload["version"]) is not int or variant_payload["version"] != VARIANT_VERSION:
            raise _fail(f"invalid variant collection version: {variant_payload['version']!r}")
        if variant_payload["vrs_version"] != VRS_VERSION:
            raise _fail(f"unsupported vrs_version pin: {variant_payload['vrs_version']!r}")
        if not _exact_equal(variant_payload["producer"], EXPECTED_PRODUCER):
            raise _fail("unexpected variant producer object")

        src_meta = _keys(
            variant_payload["source"],
            {"format", "version", "profile", "source_ir_sha256", "artifact_sha256", "records"},
            "variant.source",
        )
        if (
            src_meta["format"] != source_desc["format"]
            or type(src_meta["version"]) is not int
            or src_meta["version"] != source_desc["version"]
            or src_meta["profile"] != source_desc["source_ir"]["profile"]
            or src_meta["source_ir_sha256"] != source_desc["source_ir_sha256"]
            or src_meta["artifact_sha256"] != source_desc["artifact_sha256"]
            or not _exact_equal(src_meta["records"], source_desc["source_ir"]["records"])
        ):
            raise _fail("variant.source binding disagrees with source.json")
        checks.append({"name": "source-descriptor-binding", "passed": True})

        variants_ir = _keys(variant_payload["variants_ir"], {"variants"}, "variants_ir")
        variants = variants_ir["variants"]
        if type(variants) is not list or not 1 <= len(variants) <= MAX_VARIANTS:
            raise _fail(f"variants_ir.variants must contain 1..{MAX_VARIANTS} entries")

        rec_by_id = {rec["record_id"]: rec for rec in source_desc["source_ir"]["records"]}
        seen_ids: set[str] = set()
        allele_ids: list[str] = []

        for idx, entry in enumerate(variants):
            lbl = f"variants[{idx}]"
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
                lbl,
            )
            if type(entry["ordinal"]) is not int or entry["ordinal"] != idx:
                raise _fail(f"{lbl}.ordinal must equal {idx}")
            vid = entry["id"]
            if type(vid) is not str or not vid or vid != vid.strip() or len(_utf8_bytes(vid, f"{lbl}.id")) > MAX_IDENTIFIER_BYTES:
                raise _fail(f"{lbl}.id is invalid")
            if vid in seen_ids:
                raise _fail(f"duplicate variant id {vid!r}")
            seen_ids.add(vid)

            rid = entry["record_id"]
            if type(rid) is not str or rid not in rec_by_id:
                raise _fail(f"{lbl}.record_id {rid!r} absent from source")
            rec = rec_by_id[rid]
            ref_seq = sequences[rid]
            if entry["refget_id"] != rec["refget_id"] or entry["sequence_sha256"] != rec["sequence_sha256"]:
                raise _fail(f"{lbl} refget_id/sequence_sha256 mismatch")

            loc = _keys(entry["location"], {"type", "sequenceReference", "start", "end"}, f"{lbl}.location")
            if loc["type"] != "SequenceLocation":
                raise _fail(f"{lbl}.location.type invalid")
            seq_ref = _keys(loc["sequenceReference"], {"type", "refgetAccession"}, f"{lbl}.sequenceReference")
            if seq_ref["type"] != "SequenceReference" or seq_ref["refgetAccession"] != rec["refget_id"]:
                raise _fail(f"{lbl}.sequenceReference mismatch")
            start = loc["start"]
            end = loc["end"]
            if type(start) is not int or type(end) is not int or not 0 <= start <= end <= len(ref_seq):
                raise _fail(f"{lbl} coordinates [{start}, {end}) out of bounds [0, {len(ref_seq)}]")

            ref = entry["ref"]
            alt = entry["alt"]
            if type(ref) is not str or type(alt) is not str:
                raise _fail(f"{lbl} ref and alt must be strings")
            if len(ref) > MAX_VARIANT_BASES or len(alt) > MAX_VARIANT_BASES:
                raise _fail(f"{lbl} ref or alt exceeds {MAX_VARIANT_BASES} bases")
            if any(c not in IUPAC_DNA for c in ref) or any(c not in IUPAC_DNA for c in alt):
                raise _fail(f"{lbl} contains invalid DNA characters")
            if start == end and (ref != "" or not alt):
                raise _fail(f"{lbl} insertion requires empty ref and nonempty alt")
            if start < end and len(ref) != end - start:
                raise _fail(f"{lbl} ref length disagrees with end - start")
            if ref_seq[start:end] != ref:
                raise _fail(f"{lbl} reference sequence slice mismatch")
            if ref == alt:
                raise _fail(f"{lbl} no-op edit (ref == alt) is forbidden")
            if entry["edit_kind"] != _independent_edit_kind(ref, alt):
                raise _fail(f"{lbl}.edit_kind mismatch")

            state = _keys(entry["state"], {"type", "sequence"}, f"{lbl}.state")
            if state["type"] != "LiteralSequenceExpression" or state["sequence"] != alt:
                raise _fail(f"{lbl}.state mismatch")

            loc_digest = _independent_vrs_location_digest(rec["refget_id"], start, end)
            allele_digest = _independent_vrs_allele_digest(loc_digest, alt)
            expected_vrs_loc_id = f"ga4gh:SL.{loc_digest}"
            expected_vrs_allele_id = f"ga4gh:VA.{allele_digest}"
            expected_vrs_allele = {
                "id": expected_vrs_allele_id,
                "type": "Allele",
                "digest": allele_digest,
                "location": {
                    "id": expected_vrs_loc_id,
                    "type": "SequenceLocation",
                    "digest": loc_digest,
                    "sequenceReference": {
                        "type": "SequenceReference",
                        "refgetAccession": rec["refget_id"],
                    },
                    "start": start,
                    "end": end,
                },
                "state": {
                    "type": "LiteralSequenceExpression",
                    "sequence": alt,
                },
            }
            if entry["vrs_location_id"] != expected_vrs_loc_id:
                raise _fail(f"{lbl}.vrs_location_id mismatch")
            if entry["vrs_allele_id"] != expected_vrs_allele_id:
                raise _fail(f"{lbl}.vrs_allele_id mismatch")
            if not _exact_equal(entry["vrs_allele"], expected_vrs_allele):
                raise _fail(f"{lbl}.vrs_allele mismatch")

            mutated_seq = ref_seq[:start] + alt + ref_seq[end:]
            reverted_seq = mutated_seq[:start] + ref + mutated_seq[start + len(alt) :]
            if reverted_seq != ref_seq:
                raise _fail(f"{lbl} allele reversion failed exact reference identity")
            mutated_bytes = mutated_seq.encode("ascii")
            mut_sha = hashlib.sha256(mutated_bytes).hexdigest()
            mut_refget = "SQ." + _sha512t24u(mutated_bytes)
            if type(entry["mutated_bases"]) is not int or entry["mutated_bases"] != len(mutated_seq):
                raise _fail(f"{lbl}.mutated_bases mismatch")
            if entry["mutated_sequence_sha256"] != mut_sha:
                raise _fail(f"{lbl}.mutated_sequence_sha256 mismatch")
            if entry["mutated_refget_id"] != mut_refget:
                raise _fail(f"{lbl}.mutated_refget_id mismatch")
            allele_ids.append(expected_vrs_allele_id)

        checks.append({"name": "vrs-v2-identifiers-and-reversion", "passed": True})

        if variant_payload["variants_ir_sha256"] != _digest(variants_ir):
            raise _fail("variants_ir_sha256 mismatch")
        core = {k: v for k, v in variant_payload.items() if k != "artifact_sha256"}
        if variant_payload["artifact_sha256"] != _digest(core):
            raise _fail("variant artifact_sha256 mismatch")
        checks.append({"name": "rfc8785-canonical-seals", "passed": True})

        return _seal_report(
            {
                "format": REPORT_FORMAT,
                "version": REPORT_VERSION,
                "validator": VALIDATOR_INFO,
                "valid": True,
                "code": None,
                "detail": None,
                "source_artifact_sha256": source_artifact_sha,
                "variant_artifact_sha256": variant_payload["artifact_sha256"],
                "checks": checks,
                "summary": {
                    "variant_count": len(variants),
                    "vrs_version": VRS_VERSION,
                    "vrs_allele_ids": allele_ids,
                },
            }
        )
    except VariantValidationError as failure:
        return _seal_report(
            {
                "format": REPORT_FORMAT,
                "version": REPORT_VERSION,
                "validator": VALIDATOR_INFO,
                "valid": False,
                "code": "VARVAL001",
                "detail": str(failure),
                "source_artifact_sha256": source_artifact_sha,
                "variant_artifact_sha256": variant_artifact_sha,
                "checks": [*checks, {"name": "variant-replay", "passed": False}],
                "summary": None,
            }
        )


def _atomic_write(path: str | Path, content: bytes) -> None:
    target = Path(path)
    parent = target.parent
    if not parent.is_dir():
        raise _fail(f"output parent directory does not exist: {parent}")
    tmp_name = f".{target.name}.{secrets.token_hex(8)}.tmp"
    if _HAS_DIR_FD:
        try:
            pfd = os.open(parent, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
        except OSError as failure:
            raise _fail(f"cannot open output parent directory: {failure}") from failure
        try:
            fd = os.open(tmp_name, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o644, dir_fd=pfd)
            try:
                view = memoryview(content)
                written = 0
                while written < len(view):
                    written += os.write(fd, view[written:])
                os.fsync(fd)
            finally:
                os.close(fd)
            os.link(tmp_name, target.name, src_dir_fd=pfd, dst_dir_fd=pfd)
            os.unlink(tmp_name, dir_fd=pfd)
        except OSError as failure:
            try:
                os.unlink(tmp_name, dir_fd=pfd)
            except OSError:
                pass
            raise _fail(f"cannot write validation report: {failure}") from failure
        finally:
            os.close(pfd)
    else:
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_BINARY", 0)
        try:
            fd = os.open(target, flags, 0o644)
        except OSError as failure:
            raise _fail(f"cannot write validation report: {failure}") from failure
        try:
            view = memoryview(content)
            written = 0
            while written < len(view):
                written += os.write(fd, view[written:])
        finally:
            os.close(fd)


def validate_variant_paths(
    source_bundle_dir: str | Path,
    variant_collection_path: str | Path,
    *,
    source_inputs: Mapping[str, str | Path] | None = None,
    report_path: str | Path | None = None,
) -> dict[str, Any]:
    """Read sealed source bundle + variant collection from paths and independently validate."""
    try:
        bundle_files = _read_source_bundle_files(source_bundle_dir)
        variant_raw = _read_regular_file(variant_collection_path, "variant collection", maximum_bytes=MAX_JSON_BYTES)
        variant_payload = _parse_json_object(variant_raw, "variant collection")
        report = validate_variant_collection(bundle_files, variant_payload, source_inputs=source_inputs)
    except VariantValidationError as failure:
        report = _seal_report(
            {
                "format": REPORT_FORMAT,
                "version": REPORT_VERSION,
                "validator": VALIDATOR_INFO,
                "valid": False,
                "code": "VARVAL001",
                "detail": str(failure),
                "source_artifact_sha256": None,
                "variant_artifact_sha256": None,
                "checks": [{"name": "load-inputs", "passed": False}],
                "summary": None,
            }
        )
    if report_path is not None:
        encoded = (json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n").encode("utf-8")
        _atomic_write(report_path, encoded)
    return report


def _role_path(value: str) -> tuple[str, str]:
    role, sep, path = value.partition("=")
    if not sep or not role or not path:
        raise argparse.ArgumentTypeError("expected ROLE=PATH")
    return role, path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="brainc-validate-variants",
        description="Independently validate a sealed brainc.variant-collection/v1 artifact",
    )
    parser.add_argument("source_bundle", help="validated brainc.source-descriptor/v1 directory")
    parser.add_argument("variant_collection", help="sealed brainc.variant-collection/v1 JSON file")
    parser.add_argument(
        "--source-input",
        action="append",
        type=_role_path,
        default=[],
        metavar="ROLE=PATH",
        help="original sequence file for reference-catalog or external source profiles",
    )
    parser.add_argument(
        "-o",
        "--output",
        "--report",
        dest="output",
        help="atomically write the sealed validation report",
    )
    args = parser.parse_args(argv)
    source_inputs: dict[str, str] = {}
    for role, path in args.source_input:
        if role in source_inputs:
            parser.error(f"duplicate --source-input role: {role!r}")
        source_inputs[role] = path
    try:
        report = validate_variant_paths(
            args.source_bundle,
            args.variant_collection,
            source_inputs=source_inputs or None,
            report_path=args.output,
        )
    except VariantValidationError as failure:
        parser.exit(1, f"validation failed: {failure}\n")
    rendered = (json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n").encode("utf-8")
    if args.output is None:
        sys.stdout.buffer.write(rendered)
    return 0 if report["valid"] else 1


__all__ = [
    "REPORT_FORMAT",
    "REPORT_VERSION",
    "VARIANT_FORMAT",
    "VARIANT_VERSION",
    "VariantValidationError",
    "main",
    "validate_variant_collection",
    "validate_variant_paths",
]


if __name__ == "__main__":
    raise SystemExit(main())

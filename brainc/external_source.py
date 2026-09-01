"""Bind validated out-of-process frontend evidence into one compiler source."""

from __future__ import annotations

import copy
import hashlib
import os
from pathlib import Path
import re
import stat
from typing import Any, Mapping

from ._canonical import ContractError, artifact_digest, canonical_bytes, digest, loads
from .external_profile import (
    MAX_EVIDENCE_BYTES,
    MAX_PROFILE_BYTES,
    validate_external_evidence,
    validate_profile_manifest,
    validate_source_descriptor,
    validate_validation_report,
)


FORMAT = "brainc.external-source-closure"
VERSION = 1
PROFILE = "external-dna-source/v1"
MAX_CLOSURE_BYTES = MAX_PROFILE_BYTES + 2 * MAX_EVIDENCE_BYTES + 1024 * 1024
PRODUCER = {
    "name": "brainc-external-source",
    "version": "1.0.0",
    "passes": [
        "validate-profile-manifest",
        "validate-frontend-evidence",
        "bind-external-replay",
        "emit-external-source-closure",
    ],
}

_SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")


class ExternalSourceError(ContractError):
    """The external frontend closure is invalid or incomplete."""


def _fail(detail: str) -> ExternalSourceError:
    return ExternalSourceError(f"EXTSRC001: {detail}")


def _keys(value: Any, expected: set[str], label: str) -> dict[str, Any]:
    if type(value) is not dict:
        raise _fail(f"{label} must be an object")
    if any(type(key) is not str for key in value):
        raise _fail(f"{label} keys must be strings")
    missing = sorted(expected - value.keys())
    extra = sorted(value.keys() - expected)
    if missing or extra:
        raise _fail(
            f"{label} keys invalid; missing={missing or 'none'}, "
            f"unknown={extra or 'none'}"
        )
    return value


def _sha256(value: Any, label: str) -> str:
    if type(value) is not str or _SHA256_RE.fullmatch(value) is None:
        raise _fail(f"{label} must be a lowercase SHA-256 digest")
    return value


def _bounded(value: dict[str, Any]) -> None:
    try:
        size = len(canonical_bytes(value))
    except ContractError as failure:
        raise _fail(f"external source closure is not canonical I-JSON: {failure}") from failure
    if size > MAX_CLOSURE_BYTES:
        raise _fail(f"external source closure exceeds {MAX_CLOSURE_BYTES} bytes")


def _seal_closure(evidence: dict[str, dict[str, Any]]) -> dict[str, Any]:
    closure_ir = {
        "profile_manifest": evidence["profile_manifest"],
        "source_descriptor": evidence["source_descriptor"],
        "validation_report": evidence["validation_report"],
    }
    core = {
        "format": FORMAT,
        "version": VERSION,
        "producer": copy.deepcopy(PRODUCER),
        "closure_ir": closure_ir,
        "closure_ir_sha256": digest(closure_ir),
    }
    return validate_external_source_closure(
        {**core, "artifact_sha256": digest(core)}
    )


def build_external_source_closure(
    profile_manifest: dict[str, Any],
    source_descriptor: dict[str, Any],
    validation_report: dict[str, Any],
    *,
    original_payloads: Mapping[str, bytes],
    native_artifact_payloads: Mapping[str, bytes],
) -> dict[str, Any]:
    """Validate exact evidence bytes and seal the portable compiler-side closure."""

    try:
        evidence = validate_external_evidence(
            profile_manifest,
            source_descriptor,
            validation_report,
            original_payloads=original_payloads,
            native_artifact_payloads=native_artifact_payloads,
        )
        return _seal_closure(evidence)
    except ExternalSourceError:
        raise
    except ContractError as failure:
        raise _fail(f"external frontend evidence is invalid: {failure}") from failure
    except MemoryError as failure:
        raise _fail("external frontend evidence exceeds the memory ceiling") from failure


def _path_mapping(value: Any, expected: set[str], label: str) -> dict[str, Path]:
    if not isinstance(value, Mapping):
        raise _fail(f"{label} must be a mapping")
    try:
        observed_roles = len(value)
    except (TypeError, ValueError, OverflowError) as failure:
        raise _fail(f"{label} is not a valid role mapping") from failure
    if observed_roles != len(expected):
        raise _fail(f"{label} role closure does not match the descriptor")
    result: dict[str, Path] = {}
    for role in sorted(expected):
        try:
            raw_path = value[role]
        except (KeyError, TypeError, ValueError) as failure:
            raise _fail(f"{label} role closure does not match the descriptor") from failure
        if not isinstance(raw_path, (str, Path)) or not os.fspath(raw_path):
            raise _fail(f"{label}.{role} must be a filesystem path")
        result[role] = Path(raw_path)
    return result


def _read_path(
    path: Path,
    *,
    maximum_bytes: int,
    label: str,
    capture: bool,
) -> tuple[dict[str, Any], bytes | None]:
    descriptor = -1
    try:
        inspected = path.lstat()
        if (
            stat.S_ISLNK(inspected.st_mode)
            or not stat.S_ISREG(inspected.st_mode)
            or inspected.st_nlink != 1
        ):
            raise _fail(f"{label} must be one regular non-linked file")
        if inspected.st_size > maximum_bytes:
            raise _fail(f"{label} exceeds byte limit {maximum_bytes}")
        flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
        flags |= getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
        descriptor = os.open(path, flags)
        opened = os.fstat(descriptor)
        stable = (
            "st_dev",
            "st_ino",
            "st_nlink",
            "st_size",
            "st_mtime_ns",
            "st_ctime_ns",
        )
        if (
            not stat.S_ISREG(opened.st_mode)
            or opened.st_nlink != 1
            or any(
                getattr(inspected, field) != getattr(opened, field)
                for field in stable
            )
        ):
            raise _fail(f"{label} changed while opening")
        hasher = hashlib.sha256()
        chunks: list[bytes] | None = [] if capture else None
        byte_length = 0
        while True:
            raw = os.read(descriptor, min(64 * 1024, maximum_bytes + 1 - byte_length))
            if not raw:
                break
            byte_length += len(raw)
            if byte_length > maximum_bytes:
                raise _fail(f"{label} exceeds byte limit {maximum_bytes}")
            hasher.update(raw)
            if chunks is not None:
                chunks.append(raw)
        finished = os.fstat(descriptor)
        attached = path.lstat()
        if (
            byte_length != opened.st_size
            or stat.S_ISLNK(attached.st_mode)
            or any(getattr(opened, field) != getattr(finished, field) for field in stable)
            or any(getattr(opened, field) != getattr(attached, field) for field in stable)
        ):
            raise _fail(f"{label} changed while reading")
        return (
            {"sha256": hasher.hexdigest(), "byte_length": byte_length},
            None if chunks is None else b"".join(chunks),
        )
    except ExternalSourceError:
        raise
    except (OSError, TypeError, ValueError) as failure:
        raise _fail(f"cannot read {label}: {failure}") from failure
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def build_external_source_closure_from_paths(
    profile_manifest: dict[str, Any],
    source_descriptor: dict[str, Any],
    validation_report: dict[str, Any],
    *,
    original_paths: Mapping[str, str | Path],
    native_artifact_paths: Mapping[str, str | Path],
) -> dict[str, Any]:
    """Stream exact external inputs and seal their data-only frontend evidence."""

    try:
        manifest = validate_profile_manifest(profile_manifest)
        descriptor = validate_source_descriptor(manifest, source_descriptor)
        report = validate_validation_report(
            manifest,
            descriptor,
            validation_report,
        )
        input_references = descriptor["frontend_ir"]["inputs"]
        native_references = descriptor["frontend_ir"]["native_artifacts"]
        input_declarations = {
            declaration["role"]: declaration
            for declaration in manifest["profile_ir"]["inputs"]
        }
        native_declarations = {
            declaration["role"]: declaration
            for declaration in manifest["profile_ir"]["native_artifacts"]
        }
        input_paths = _path_mapping(
            original_paths,
            {reference["role"] for reference in input_references},
            "original_paths",
        )
        native_paths = _path_mapping(
            native_artifact_paths,
            {reference["role"] for reference in native_references},
            "native_artifact_paths",
        )
        for reference in input_references:
            role = reference["role"]
            observed, _ = _read_path(
                input_paths[role],
                maximum_bytes=input_declarations[role]["maximum_byte_length"],
                label=f"original source {role}",
                capture=False,
            )
            if observed != {
                "sha256": reference["sha256"],
                "byte_length": reference["byte_length"],
            }:
                raise _fail(f"original source {role} differs from its exact reference")
        for reference in native_references:
            role = reference["role"]
            observed, raw = _read_path(
                native_paths[role],
                maximum_bytes=native_declarations[role]["maximum_byte_length"],
                label=f"native artifact {role}",
                capture=True,
            )
            if observed != {
                "sha256": reference["sha256"],
                "byte_length": reference["byte_length"],
            }:
                raise _fail(f"native artifact {role} differs from its exact reference")
            assert raw is not None
            artifact = loads(raw, f"native artifact {role}")
            declaration = native_declarations[role]
            if artifact.get("format") != declaration["format"] or (
                type(artifact.get("version")) is not int
                or artifact["version"] != declaration["version"]
            ):
                raise _fail(f"native artifact {role} has the wrong format/version")
            ir_field = declaration["ir_digest_field"]
            if artifact.get(ir_field) != reference["ir_sha256"]:
                raise _fail(f"native artifact {role} exposes a different IR digest")
        return _seal_closure(
            {
                "profile_manifest": manifest,
                "source_descriptor": descriptor,
                "validation_report": report,
            }
        )
    except ExternalSourceError:
        raise
    except ContractError as failure:
        raise _fail(f"external frontend evidence is invalid: {failure}") from failure
    except MemoryError as failure:
        raise _fail("external frontend evidence exceeds the memory ceiling") from failure


def validate_external_source_closure(payload: dict[str, Any]) -> dict[str, Any]:
    """Validate the closed manifest/descriptor/replay relationship without code loading."""

    item = _keys(
        payload,
        {
            "format",
            "version",
            "producer",
            "closure_ir",
            "closure_ir_sha256",
            "artifact_sha256",
        },
        "external source closure",
    )
    if item["format"] != FORMAT or type(item["version"]) is not int or item["version"] != VERSION:
        raise _fail("external source closure format/version is unsupported")
    if type(item["producer"]) is not dict or item["producer"] != PRODUCER:
        raise _fail("external source closure producer is unsupported")
    closure = _keys(
        item["closure_ir"],
        {"profile_manifest", "source_descriptor", "validation_report"},
        "external source closure.closure_ir",
    )
    try:
        manifest = validate_profile_manifest(closure["profile_manifest"])
        descriptor = validate_source_descriptor(
            manifest,
            closure["source_descriptor"],
        )
        validate_validation_report(
            manifest,
            descriptor,
            closure["validation_report"],
        )
    except ContractError as failure:
        raise _fail(f"external source closure evidence is invalid: {failure}") from failure
    if _sha256(item["closure_ir_sha256"], "closure_ir_sha256") != digest(closure):
        raise _fail("closure_ir_sha256 does not match closure_ir")
    try:
        artifact_digest(item)
    except (ContractError, KeyError) as failure:
        raise _fail(f"external source closure digest is invalid: {failure}") from failure
    _bounded(item)
    return copy.deepcopy(item)


def profile_parameters(payload: dict[str, Any]) -> dict[str, Any]:
    closure = validate_external_source_closure(payload)["closure_ir"]
    manifest = closure["profile_manifest"]
    profile = manifest["profile_ir"]["profile"]
    return {
        "profile_id": profile["id"],
        "profile_version": profile["version"],
        "profile_manifest_sha256": manifest["artifact_sha256"],
    }


def input_references(payload: dict[str, Any]) -> dict[str, dict[str, Any]]:
    closure = validate_external_source_closure(payload)["closure_ir"]
    references = closure["source_descriptor"]["frontend_ir"]["inputs"]
    return {
        reference["role"]: {
            "sha256": reference["sha256"],
            "byte_length": reference["byte_length"],
        }
        for reference in references
    }


def source_records(payload: dict[str, Any]) -> list[dict[str, Any]]:
    closure = validate_external_source_closure(payload)["closure_ir"]
    records = closure["source_descriptor"]["frontend_ir"]["record_catalog"][
        "records"
    ]
    return [
        {
            "record_id": record["record_id"],
            "bases": record["bases"],
            "sequence_sha256": record["sequence_sha256"],
            "refget_id": record["refget_id"],
        }
        for record in records
    ]


__all__ = [
    "ExternalSourceError",
    "FORMAT",
    "MAX_CLOSURE_BYTES",
    "PRODUCER",
    "PROFILE",
    "VERSION",
    "build_external_source_closure",
    "build_external_source_closure_from_paths",
    "input_references",
    "profile_parameters",
    "source_records",
    "validate_external_source_closure",
]

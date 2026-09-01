"""Independent replay validator for compiler-only development bundles.

The trust boundary contains only Python's standard library and the existing
independent biological/v2 validators.  It deliberately imports no source,
bundle, provider, target, policy, tensor, or compiler producer module.
"""

from __future__ import annotations

import argparse
import functools
import hashlib
import json
import math
import os
from pathlib import Path
import re
import stat
from typing import Any

from . import (
    validator_bio,
    validator_external,
    validator_insdc,
    validator_reference,
    validator_v2,
)


SAFE_INTEGER = 2**53 - 1
MAX_INPUT_BYTES = 16 * 1024 * 1024
MAX_GFF3_BYTES = 64 * 1024 * 1024
MAX_JSON_BYTES = 16 * 1024 * 1024
MAX_GENBANK_ARTIFACT_BYTES = 64 * 1024 * 1024
MAX_GFF3_ARTIFACT_BYTES = 256 * 1024 * 1024
MAX_SOURCE_DESCRIPTOR_BYTES = MAX_GFF3_ARTIFACT_BYTES
MAX_JSON_DEPTH = 64
MAX_JSON_MEMBERS = 1_000_000
MAX_GFF3_JSON_MEMBERS = 5_000_000
MAX_STRING_BYTES = 1 * 1024 * 1024
MAX_IDENTIFIER_BYTES = 256
MAX_FAILURE_MESSAGE_BYTES = 16 * 1024
MAX_SOURCE_RECORDS = 100_000
MAX_EXTERNAL_BLOBS = 100_000
MAX_EXTERNAL_BLOB_BYTES = 256 * 1024 * 1024

SOURCE_FORMAT = "brainc.source-descriptor"
SOURCE_VERSION = 1
RAW_PROFILE = "raw-iupac-dna/v1"
FASTA_PROFILE = "fasta-dna/v1"
REFERENCE_FASTA_PROFILE = "fasta-reference-dna/v1"
EXTERNAL_PROFILE = "external-dna-source/v1"
GENBANK_PROFILE = "genbank-273-traditional-dna-physical-structural/v2"
GFF3_PROFILE = "gff3-external-sequence/v1"
PROFILES = frozenset(
    {
        RAW_PROFILE,
        FASTA_PROFILE,
        REFERENCE_FASTA_PROFILE,
        EXTERNAL_PROFILE,
        GENBANK_PROFILE,
        GFF3_PROFILE,
    }
)
SOURCE_PRODUCER = {
    "name": "brainc-source",
    "version": "1.0.0",
    "passes": [
        "select-source-profile",
        "compile-native-source",
        "bind-source-evidence",
        "emit-source-descriptor",
    ],
}
SOURCE_FILENAME = "source.json"
SOURCE_FILENAMES = {
    "sequence": "sequence.json",
    "reference": "reference.json",
    "external": "external.json",
    "genbank": "genbank.json",
    "annotation": "annotation.json",
}
SOURCE_FILE_LIMITS = {
    SOURCE_FILENAME: MAX_SOURCE_DESCRIPTOR_BYTES,
    "sequence.json": MAX_JSON_BYTES,
    "reference.json": validator_reference.MAX_ARTIFACT_BYTES,
    "external.json": validator_external.MAX_CLOSURE_BYTES,
    "genbank.json": MAX_GENBANK_ARTIFACT_BYTES,
    "annotation.json": MAX_GFF3_ARTIFACT_BYTES,
}
SOURCE_IDENTITIES = {
    "sequence": ("brain01.sequence-collection-ir", 1, "collection_ir_sha256"),
    "reference": (
        validator_reference.FORMAT,
        validator_reference.VERSION,
        "catalog_sha256",
    ),
    "external": (
        validator_external.CLOSURE_FORMAT,
        validator_external.CLOSURE_VERSION,
        "closure_ir_sha256",
    ),
    "genbank": ("brainc.bio.insdc-genbank-ir", 2, "bio_ir_sha256"),
    "annotation": ("brainc.bio.gff3-ir", 1, "bio_ir_sha256"),
}

BUNDLE_FORMAT = "brainc.development-bundle"
BUNDLE_VERSION = 1
COMPILATION_FORMAT = "brainc.development-compilation"
COMPILATION_VERSION = 1
PROFILE = "development-module/v1"
BUNDLE_PRODUCER = {
    "name": "brainc-development-bundle",
    "version": "1.0.0",
    "passes": [
        "validate-source-bundle",
        "bind-caller-interpretation",
        "type-check-development-module",
        "emit-reference-only-bundle",
    ],
}
CHILD_FILENAMES = {
    "compilation_record": "compilation_record.json",
    "development_module": "development_module.json",
    "lowering_policy": "lowering_policy.json",
    "prediction_request": "prediction_request.json",
    "prediction_response": "prediction_response.json",
    "provider_manifest": "provider_manifest.json",
    "target_contract": "target_contract.json",
}
CHILD_IDENTITIES = {
    "compilation_record": (COMPILATION_FORMAT, COMPILATION_VERSION),
    "development_module": ("brainc.development-module", 1),
    "lowering_policy": ("brainc.lowering-policy", 2),
    "prediction_request": ("brainc.prediction-request", 2),
    "prediction_response": ("brainc.prediction-response", 2),
    "provider_manifest": ("brainc.provider-manifest", 2),
    "target_contract": ("brainc.target-contract", 1),
}
BUNDLE_FILENAME = "bundle.json"
BUNDLE_FILENAMES = frozenset({BUNDLE_FILENAME, *CHILD_FILENAMES.values()})
REPORT_FORMAT = "brainc.development-validation-report"
REPORT_VERSION = 1
REPORT_CHECKS = (
    "original-source-and-native-artifact-replay",
    "closed-source-descriptor",
    "provider-request-response-binding",
    "target-policy-and-development-module-replay",
    "external-blob-identity-and-closure",
    "development-compilation-record-replay",
    "reference-only-development-bundle-replay",
)

SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
REFGET_RE = re.compile(r"SQ\.[A-Za-z0-9_-]{32}\Z")
EXTERNAL_ROLE_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9._:/+@-]{0,255}\Z")
_HAS_DIRECTORY_DESCRIPTOR = (
    os.name == "posix"
    and os.listdir in getattr(os, "supports_fd", set())
    and all(
        operation in getattr(os, "supports_dir_fd", set())
        for operation in (os.open, os.stat)
    )
    and getattr(os, "O_DIRECTORY", 0) != 0
    and getattr(os, "O_NOFOLLOW", 0) != 0
)


class DevelopmentValidationError(ValueError):
    """A development compilation differs from independent replay."""


def _fail(detail: str) -> DevelopmentValidationError:
    return DevelopmentValidationError(f"DEVVAL001: {detail}")


def _bounded(function: Any) -> Any:
    @functools.wraps(function)
    def guarded(*args: Any, **kwargs: Any) -> Any:
        try:
            return function(*args, **kwargs)
        except DevelopmentValidationError:
            raise
        except (MemoryError, RecursionError) as failure:
            raise _fail("validation resource ceiling exceeded") from failure

    return guarded


def _keys(value: Any, expected: set[str], label: str) -> dict[str, Any]:
    if type(value) is not dict:
        raise _fail(f"{label} must be an object")
    if len(value) != len(expected):
        raise _fail(
            f"{label} keys invalid; missing/unknown or non-string object key; "
            f"expected={len(expected)}, observed={len(value)}"
        )
    if any(type(key) is not str for key in value):
        raise _fail(f"{label} contains a non-string object key")
    missing = expected - value.keys()
    extra = value.keys() - expected
    if missing or extra:
        raise _fail(
            f"{label} keys invalid; missing={_key_summary(missing)}, "
            f"unknown={_key_summary(extra)}"
        )
    return value


def _key_summary(values: list[str] | set[str]) -> str:
    if not values:
        return "none"
    if len(values) > 8:
        return f"{len(values)} keys"
    try:
        oversized = any(len(value.encode("utf-8")) > 128 for value in values)
    except UnicodeEncodeError:
        return f"{len(values)} keys"
    if oversized:
        return f"{len(values)} keys"
    rendered = repr(sorted(values))
    if len(rendered.encode("utf-8")) > 1024:
        return f"{len(values)} keys"
    return rendered


def _external_roles(value: Any, label: str) -> set[str]:
    if type(value) is not dict:
        raise _fail(f"{label} must be an object")
    if not 1 <= len(value) <= validator_external.MAX_ROLES:
        raise _fail(
            f"{label} must contain 1..{validator_external.MAX_ROLES} roles"
        )
    roles = set(value)
    if any(
        type(role) is not str or EXTERNAL_ROLE_RE.fullmatch(role) is None
        for role in roles
    ):
        raise _fail(f"{label} contains an invalid role identifier")
    return roles


def _integer(
    value: Any,
    label: str,
    *,
    minimum: int = 0,
    maximum: int = SAFE_INTEGER,
) -> int:
    if type(value) is not int or not minimum <= value <= maximum:
        raise _fail(f"{label} must be an integer in [{minimum}, {maximum}]")
    return value


def _text(value: Any, label: str, *, identifier: bool = False) -> str:
    if type(value) is not str or not value or value != value.strip():
        raise _fail(f"{label} must be a nonempty trimmed string")
    if any(0xD800 <= ord(character) <= 0xDFFF for character in value):
        raise _fail(f"{label} contains a lone Unicode surrogate")
    maximum = MAX_IDENTIFIER_BYTES if identifier else MAX_STRING_BYTES
    if len(value.encode("utf-8")) > maximum:
        raise _fail(f"{label} exceeds {maximum} UTF-8 bytes")
    if identifier and any(character.isspace() for character in value):
        raise _fail(f"{label} must be whitespace-free")
    return value


def _sha(value: Any, label: str) -> str:
    if type(value) is not str or SHA256_RE.fullmatch(value) is None:
        raise _fail(f"{label} must be a lowercase SHA-256 digest")
    return value


def _exact_equal(left: Any, right: Any) -> bool:
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


def _first_difference(actual: Any, expected: Any, path: str = "$") -> str:
    if type(actual) is not type(expected):
        return f"{path} type {type(actual).__name__} != {type(expected).__name__}"
    if type(actual) is dict:
        actual_keys, expected_keys = set(actual), set(expected)
        if actual_keys != expected_keys:
            return (
                f"{path} keys differ; "
                f"missing={_key_summary(expected_keys - actual_keys)}, "
                f"unknown={_key_summary(actual_keys - expected_keys)}"
            )
        for key in sorted(actual, key=lambda item: item.encode("utf-8")):
            if not _exact_equal(actual[key], expected[key]):
                return _first_difference(actual[key], expected[key], f"{path}.{key}")
    elif type(actual) is list:
        if len(actual) != len(expected):
            return f"{path} length {len(actual)} != {len(expected)}"
        for index, (left, right) in enumerate(zip(actual, expected)):
            if not _exact_equal(left, right):
                return _first_difference(left, right, f"{path}[{index}]")
    elif actual != expected:
        return f"{path} value differs"
    return path


def _validate_tree(
    value: Any,
    label: str,
    *,
    maximum_members: int = MAX_JSON_MEMBERS,
) -> None:
    members = 0
    stack: list[tuple[Any, int]] = [(value, 0)]
    while stack:
        current, depth = stack.pop()
        if depth > MAX_JSON_DEPTH:
            raise _fail(f"{label} exceeds JSON depth {MAX_JSON_DEPTH}")
        if current is None or type(current) is bool:
            continue
        if type(current) is int:
            if not -SAFE_INTEGER <= current <= SAFE_INTEGER:
                raise _fail(f"{label} contains an unsafe JSON integer")
            continue
        if type(current) is float:
            if not math.isfinite(current):
                raise _fail(f"{label} contains a non-finite JSON number")
            continue
        if type(current) is str:
            if any(0xD800 <= ord(character) <= 0xDFFF for character in current):
                raise _fail(f"{label} contains a lone Unicode surrogate")
            if len(current.encode("utf-8")) > MAX_STRING_BYTES:
                raise _fail(f"{label} contains an oversized JSON string")
            continue
        if type(current) is list:
            members += len(current)
            stack.extend((item, depth + 1) for item in reversed(current))
        elif type(current) is dict:
            members += len(current)
            for key, item in reversed(list(current.items())):
                if type(key) is not str:
                    raise _fail(f"{label} contains a non-string object key")
                stack.append((key, depth + 1))
                stack.append((item, depth + 1))
        else:
            raise _fail(f"{label} contains unsupported type {type(current).__name__}")
        if members > maximum_members:
            raise _fail(f"{label} exceeds JSON member ceiling {maximum_members}")


def _parse_integer(token: str) -> int:
    value = int(token)
    if not -SAFE_INTEGER <= value <= SAFE_INTEGER:
        raise _fail("JSON contains an unsafe integer")
    return value


def _parse_float(token: str) -> float:
    value = float(token)
    if not math.isfinite(value):
        raise _fail("JSON contains a non-finite number")
    return value


def _json_object(
    raw: bytes,
    label: str,
    *,
    maximum_members: int = MAX_JSON_MEMBERS,
) -> dict[str, Any]:
    def pairs(values: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in values:
            if key in result:
                raise _fail(f"{label} contains duplicate JSON key {key!r}")
            result[key] = value
        return result

    def constant(token: str) -> None:
        raise _fail(f"{label} contains non-finite number {token}")

    try:
        value = json.loads(
            raw.decode("utf-8"),
            object_pairs_hook=pairs,
            parse_constant=constant,
            parse_int=_parse_integer,
            parse_float=_parse_float,
        )
    except DevelopmentValidationError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as failure:
        raise _fail(f"invalid {label} JSON: {failure}") from failure
    if type(value) is not dict:
        raise _fail(f"{label} must contain a JSON object")
    _validate_tree(value, label, maximum_members=maximum_members)
    return value


def _digest(value: Any) -> str:
    try:
        return validator_v2.digest(value)
    except ValueError as failure:
        raise _fail(f"value is not canonical I-JSON: {failure}") from failure


def _seal(core: dict[str, Any]) -> dict[str, Any]:
    return {**core, "artifact_sha256": _digest(core)}


def _artifact_seal(value: dict[str, Any], label: str) -> None:
    claimed = _sha(value.get("artifact_sha256"), f"{label}.artifact_sha256")
    expected = _digest(
        {key: item for key, item in value.items() if key != "artifact_sha256"}
    )
    if claimed != expected:
        raise _fail(f"{label}.artifact_sha256 does not match canonical content")


def _read_descriptor(
    descriptor: int,
    opened: os.stat_result,
    label: str,
    maximum: int,
) -> bytes:
    if not stat.S_ISREG(opened.st_mode) or opened.st_nlink != 1:
        raise _fail(f"{label} must be one regular non-linked file")
    if opened.st_size > maximum:
        raise _fail(f"{label} exceeds {maximum} bytes")
    chunks: list[bytes] = []
    total = 0
    while True:
        chunk = os.read(descriptor, min(64 * 1024, maximum + 1 - total))
        if not chunk:
            break
        chunks.append(chunk)
        total += len(chunk)
        if total > maximum:
            raise _fail(f"{label} exceeds {maximum} bytes")
    finished = os.fstat(descriptor)
    stable = ("st_dev", "st_ino", "st_nlink", "st_size", "st_mtime_ns", "st_ctime_ns")
    if any(getattr(opened, field) != getattr(finished, field) for field in stable):
        raise _fail(f"{label} changed while being read")
    if total != opened.st_size:
        raise _fail(f"{label} length changed while being read")
    return b"".join(chunks)


def _read_regular(path: str | Path, label: str, maximum: int) -> bytes:
    descriptor = -1
    try:
        source = Path(path)
        inspected = source.lstat()
        if stat.S_ISLNK(inspected.st_mode) or not stat.S_ISREG(inspected.st_mode):
            raise _fail(f"{label} must be a regular non-linked file")
        flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
        flags |= getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
        descriptor = os.open(source, flags)
        opened = os.fstat(descriptor)
        if (opened.st_dev, opened.st_ino) != (inspected.st_dev, inspected.st_ino):
            raise _fail(f"{label} changed while being opened")
        raw = _read_descriptor(descriptor, opened, label, maximum)
        attached = source.lstat()
        stable = ("st_dev", "st_ino", "st_nlink", "st_size", "st_mtime_ns", "st_ctime_ns")
        if stat.S_ISLNK(attached.st_mode) or any(
            getattr(opened, field) != getattr(attached, field) for field in stable
        ):
            raise _fail(f"{label} path changed while being read")
        return raw
    except DevelopmentValidationError:
        raise
    except (OSError, TypeError, ValueError) as failure:
        raise _fail(f"cannot read {label}: {failure}") from failure
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _read_directory(
    directory: str | Path,
    limits: dict[str, int],
    *,
    required: set[str],
    label: str,
) -> dict[str, bytes]:
    descriptor = -1
    try:
        if not _HAS_DIRECTORY_DESCRIPTOR:
            raise _fail(f"secure {label} input is unsupported on this platform")
        source = Path(directory)
        inspected = source.lstat()
        if stat.S_ISLNK(inspected.st_mode) or not stat.S_ISDIR(inspected.st_mode):
            raise _fail(f"{label} must be a regular non-linked directory")
        flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW
        flags |= getattr(os, "O_CLOEXEC", 0)
        descriptor = os.open(source, flags)
        opened = os.fstat(descriptor)
        identity = (opened.st_dev, opened.st_ino)
        if not stat.S_ISDIR(opened.st_mode) or identity != (
            inspected.st_dev,
            inspected.st_ino,
        ):
            raise _fail(f"{label} changed while being opened")
        names = set(os.listdir(descriptor))
        missing, unknown = sorted(required - names), sorted(names - limits.keys())
        if missing or unknown:
            raise _fail(
                f"{label} is not closed; missing={missing or 'none'}, "
                f"unknown={unknown or 'none'}"
            )

        snapshots: dict[str, bytes] = {}
        identities: dict[str, os.stat_result] = {}
        file_flags = os.O_RDONLY | os.O_NOFOLLOW
        file_flags |= getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NONBLOCK", 0)
        for filename in sorted(names):
            inspected_child = os.stat(
                filename,
                dir_fd=descriptor,
                follow_symlinks=False,
            )
            if (
                stat.S_ISLNK(inspected_child.st_mode)
                or not stat.S_ISREG(inspected_child.st_mode)
                or inspected_child.st_nlink != 1
            ):
                raise _fail(
                    f"{label} child {filename!r} must be one regular non-linked file"
                )
            child = os.open(filename, file_flags, dir_fd=descriptor)
            try:
                opened_child = os.fstat(child)
                if (opened_child.st_dev, opened_child.st_ino) != (
                    inspected_child.st_dev,
                    inspected_child.st_ino,
                ):
                    raise _fail(f"{label} child {filename!r} changed while opening")
                snapshots[filename] = _read_descriptor(
                    child,
                    opened_child,
                    f"{label} child {filename}",
                    limits[filename],
                )
                identities[filename] = inspected_child
            finally:
                os.close(child)

        if set(os.listdir(descriptor)) != names:
            raise _fail(f"{label} filename closure changed while reading")
        stable = ("st_dev", "st_ino", "st_nlink", "st_size", "st_mtime_ns", "st_ctime_ns")
        for filename, before in identities.items():
            after = os.stat(filename, dir_fd=descriptor, follow_symlinks=False)
            if stat.S_ISLNK(after.st_mode) or any(
                getattr(before, field) != getattr(after, field) for field in stable
            ):
                raise _fail(f"{label} child {filename!r} changed while reading")
        attached = source.lstat()
        finished = os.fstat(descriptor)
        directory_stable = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
        if (
            stat.S_ISLNK(attached.st_mode)
            or not stat.S_ISDIR(attached.st_mode)
            or (attached.st_dev, attached.st_ino) != identity
            or (finished.st_dev, finished.st_ino) != identity
            or any(
                getattr(opened, field) != getattr(finished, field)
                for field in directory_stable
            )
        ):
            raise _fail(f"{label} changed while being read")
        return snapshots
    except DevelopmentValidationError:
        raise
    except (OSError, TypeError, ValueError) as failure:
        raise _fail(f"cannot read {label}: {failure}") from failure
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _route(profile: Any, parameters: Any) -> tuple[tuple[str, ...], tuple[str, ...]]:
    if type(profile) is not str or profile not in PROFILES:
        raise _fail("source profile is unsupported")
    if type(parameters) is not dict:
        raise _fail("source parameters must be an object")
    if profile == RAW_PROFILE:
        parsed = _keys(parameters, {"record_id"}, "source parameters")
        _text(parsed["record_id"], "source parameters.record_id", identifier=True)
        return ("sequence",), ("sequence",)
    if profile in {FASTA_PROFILE, REFERENCE_FASTA_PROFILE}:
        parsed = _keys(parameters, {"wrapper"}, "source parameters")
        if type(parsed["wrapper"]) is not str or parsed["wrapper"] not in {
            "identity",
            "gzip",
        }:
            raise _fail("source parameters.wrapper must be identity or gzip")
        return (
            ("sequence",),
            ("reference",) if profile == REFERENCE_FASTA_PROFILE else ("sequence",),
        )
    if profile == GENBANK_PROFILE:
        _keys(parameters, set(), "source parameters")
        return ("genbank",), ("genbank",)
    if profile == EXTERNAL_PROFILE:
        parsed = _keys(
            parameters,
            {"profile_id", "profile_version", "profile_manifest_sha256"},
            "source parameters",
        )
        _text(parsed["profile_id"], "source parameters.profile_id", identifier=True)
        _integer(
            parsed["profile_version"],
            "source parameters.profile_version",
            minimum=1,
        )
        _sha(
            parsed["profile_manifest_sha256"],
            "source parameters.profile_manifest_sha256",
        )
        return (), ("external",)

    nested = _keys(
        parameters,
        {"sequence_profile", "sequence_parameters"},
        "source parameters",
    )
    if type(nested["sequence_profile"]) is not str or nested[
        "sequence_profile"
    ] not in {RAW_PROFILE, FASTA_PROFILE}:
        raise _fail("GFF3 source sequence profile must be raw IUPAC or FASTA")
    _route(nested["sequence_profile"], nested["sequence_parameters"])
    return ("sequence", "annotation"), ("sequence", "annotation")


def _source_reference(descriptor: dict[str, Any]) -> dict[str, Any]:
    return {
        "format": SOURCE_FORMAT,
        "version": SOURCE_VERSION,
        "artifact_sha256": descriptor["artifact_sha256"],
        "ir_sha256": descriptor["source_ir_sha256"],
    }


def _native_reference(role: str, artifact: dict[str, Any]) -> dict[str, Any]:
    expected_format, expected_version, ir_field = SOURCE_IDENTITIES[role]
    if (
        artifact.get("format") != expected_format
        or type(artifact.get("version")) is not int
        or artifact["version"] != expected_version
    ):
        raise _fail(f"native source artifact {role} has the wrong format/version")
    return {
        "format": expected_format,
        "version": expected_version,
        "artifact_sha256": _sha(
            artifact.get("artifact_sha256"),
            f"native source artifact {role}.artifact_sha256",
        ),
        "ir_sha256": _sha(
            artifact.get(ir_field),
            f"native source artifact {role}.{ir_field}",
        ),
    }


def _record_catalog(
    profile: str,
    artifacts: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    if profile == GENBANK_PROFILE:
        return [
            {
                "record_id": member["record_id"],
                "bases": member["sequence"]["bases"],
                "sequence_sha256": member["sequence"]["sha256"],
                "refget_id": member["sequence"]["refget_id"],
            }
            for member in artifacts["genbank"]["sequence_collection"]["members"]
        ]
    if profile == REFERENCE_FASTA_PROFILE:
        return artifacts["reference"]["sequence_catalog"]["records"]
    if profile == EXTERNAL_PROFILE:
        records = artifacts["external"]["closure_ir"]["source_descriptor"][
            "frontend_ir"
        ]["record_catalog"]["records"]
        return [
            {
                "record_id": record["record_id"],
                "bases": record["bases"],
                "sequence_sha256": record["sequence_sha256"],
                "refget_id": record["refget_id"],
            }
            for record in records
        ]
    return [
        {
            "record_id": member["sequence_artifact"]["sequence_ir"]["record_id"],
            "bases": len(member["sequence_artifact"]["sequence_ir"]["sequence"]),
            "sequence_sha256": member["sequence_artifact"]["sequence_ir"][
                "canonical_sha256"
            ],
            "refget_id": member["sequence_artifact"]["sequence_ir"]["refget_id"],
        }
        for member in artifacts["sequence"]["members"]
    ]


def _sequence_profile_matches(
    profile: str,
    parameters: dict[str, Any],
    sequence: dict[str, Any],
) -> None:
    if profile == RAW_PROFILE:
        expected = ("raw-iupac", None)
    else:
        expected = (
            "fasta",
            None if parameters["wrapper"] == "identity" else "gzip",
        )
    actual = (sequence["inputs"]["kind"], sequence["inputs"]["wrapper"])
    if actual != expected:
        raise _fail("declared sequence source profile does not match native artifact")


def _validate_native_sources(
    descriptor: dict[str, Any],
    artifacts: dict[str, dict[str, Any]],
    source_inputs: dict[str, Any],
    native_artifact_paths: Any,
) -> dict[str, Any] | None:
    source_ir = descriptor["source_ir"]
    profile = source_ir["profile"]
    parameters = source_ir["parameters"]
    external_report: dict[str, Any] | None = None
    try:
        if profile in {RAW_PROFILE, FASTA_PROFILE}:
            validator_v2.validate_source_artifact(
                artifacts["sequence"],
                raw_source=source_inputs["sequence"],
                record_id=parameters["record_id"] if profile == RAW_PROFILE else None,
            )
            _sequence_profile_matches(profile, parameters, artifacts["sequence"])
        elif profile == REFERENCE_FASTA_PROFILE:
            validator_reference.validate_reference_fasta(
                source_inputs["sequence"],
                artifacts["reference"],
            )
            reference_inputs = artifacts["reference"]["inputs"]
            if (
                reference_inputs["profile"] != REFERENCE_FASTA_PROFILE
                or reference_inputs["wrapper"] != parameters["wrapper"]
            ):
                raise _fail(
                    "declared reference FASTA profile does not match native artifact"
                )
        elif profile == GENBANK_PROFILE:
            validator_insdc.validate_genbank(
                artifacts["genbank"],
                genbank_source=source_inputs["genbank"],
            )
            if artifacts["genbank"].get("profile") != GENBANK_PROFILE:
                raise _fail("declared GenBank profile does not match native artifact")
        elif profile == EXTERNAL_PROFILE:
            if native_artifact_paths is None:
                raise _fail("external native artifact paths are required")
            external_report = validator_external.validate_external_source(
                artifacts["external"],
                original_paths=source_inputs,
                native_artifact_paths=native_artifact_paths,
            )
        else:
            nested_profile = parameters["sequence_profile"]
            nested_parameters = parameters["sequence_parameters"]
            validator_bio.validate_bio_chain(
                artifacts["annotation"],
                artifacts["sequence"],
                fasta_source=source_inputs["sequence"],
                gff3_source=source_inputs["annotation"],
                record_id=(
                    nested_parameters["record_id"]
                    if nested_profile == RAW_PROFILE
                    else None
                ),
            )
            _sequence_profile_matches(
                nested_profile,
                nested_parameters,
                artifacts["sequence"],
            )
    except DevelopmentValidationError:
        raise
    except (OSError, ValueError) as failure:
        raise _fail(f"independent native source replay failed: {failure}") from failure
    return external_report


def _validate_source_descriptor(
    descriptor_value: Any,
    artifacts_value: Any,
    source_inputs_value: Any,
    native_artifact_paths_value: Any = None,
) -> tuple[
    dict[str, Any],
    dict[str, dict[str, Any]],
    dict[str, dict[str, Any]],
    dict[str, int],
]:
    descriptor = _keys(
        descriptor_value,
        {
            "format",
            "version",
            "producer",
            "source_ir",
            "source_ir_sha256",
            "artifact_sha256",
        },
        "source descriptor",
    )
    _validate_tree(descriptor, "source descriptor")
    if (
        descriptor["format"] != SOURCE_FORMAT
        or type(descriptor["version"]) is not int
        or descriptor["version"] != SOURCE_VERSION
        or not _exact_equal(descriptor["producer"], SOURCE_PRODUCER)
    ):
        raise _fail("source descriptor identity or producer is unsupported")
    source_ir = _keys(
        descriptor["source_ir"],
        {"profile", "parameters", "inputs", "artifacts", "records"},
        "source descriptor.source_ir",
    )
    input_roles, artifact_roles = _route(source_ir["profile"], source_ir["parameters"])
    artifacts = _keys(
        artifacts_value,
        set(artifact_roles),
        "native source artifacts",
    )
    for role, artifact in artifacts.items():
        if type(artifact) is not dict:
            raise _fail(f"native source artifact {role} must be an object")
        _validate_tree(
            artifact,
            f"native source artifact {role}",
            maximum_members=(
                MAX_GFF3_JSON_MEMBERS if role == "annotation" else MAX_JSON_MEMBERS
            ),
        )
    if source_ir["profile"] == EXTERNAL_PROFILE:
        if type(source_inputs_value) is not dict:
            raise _fail("external original source inputs must be a path mapping")
        expected_input_roles = _external_roles(
            source_ir["inputs"],
            "external source descriptor inputs",
        )
        source_inputs = _keys(
            source_inputs_value,
            expected_input_roles,
            "external original source inputs",
        )
    else:
        source_inputs = _keys(
            source_inputs_value,
            set(input_roles),
            "original source inputs",
        )
    if source_ir["profile"] in {REFERENCE_FASTA_PROFILE, EXTERNAL_PROFILE}:
        if source_ir["profile"] == REFERENCE_FASTA_PROFILE:
            path_items = {"sequence": source_inputs["sequence"]}
        else:
            path_items = source_inputs
        for role, source_path in path_items.items():
            if not isinstance(source_path, (str, os.PathLike)):
                raise _fail(f"original source input {role} must be a path")
    else:
        for role in input_roles:
            raw = source_inputs[role]
            maximum = MAX_GFF3_BYTES if role == "annotation" else MAX_INPUT_BYTES
            if type(raw) is not bytes or len(raw) > maximum:
                raise _fail(f"original source input {role} must be bounded bytes")
    if source_ir["profile"] != EXTERNAL_PROFILE and native_artifact_paths_value is not None:
        _keys(
            native_artifact_paths_value,
            set(),
            "external native artifact paths",
        )
    external_report = _validate_native_sources(
        descriptor,
        artifacts,
        source_inputs,
        native_artifact_paths_value,
    )
    if source_ir["profile"] == EXTERNAL_PROFILE:
        assert external_report is not None
        profile_reference = external_report["profile"]
        expected_parameters = {
            "profile_id": profile_reference["id"],
            "profile_version": profile_reference["version"],
            "profile_manifest_sha256": profile_reference["manifest_sha256"],
        }
        if not _exact_equal(source_ir["parameters"], expected_parameters):
            raise _fail(
                "external source profile parameters differ from independent replay"
            )
        expected_inputs = {
            reference["role"]: {
                "sha256": reference["sha256"],
                "byte_length": reference["byte_length"],
            }
            for reference in external_report["result"]["input_references"]
        }
    elif source_ir["profile"] == REFERENCE_FASTA_PROFILE:
        expected_inputs = {
            "sequence": artifacts["reference"]["inputs"]["source"]
        }
    else:
        expected_inputs = {
            role: {
                "sha256": hashlib.sha256(source_inputs[role]).hexdigest(),
                "byte_length": len(source_inputs[role]),
            }
            for role in input_roles
        }
    if not _exact_equal(source_ir["inputs"], expected_inputs):
        detail = (
            "source descriptor input identities differ from original source"
            if source_ir["profile"]
            in {REFERENCE_FASTA_PROFILE, EXTERNAL_PROFILE}
            else "source descriptor input digests differ from original bytes"
        )
        raise _fail(detail)
    expected_artifacts = {
        role: _native_reference(role, artifacts[role]) for role in artifact_roles
    }
    if not _exact_equal(source_ir["artifacts"], expected_artifacts):
        raise _fail("source descriptor native references differ from replayed artifacts")
    expected_records = _record_catalog(source_ir["profile"], artifacts)
    if not expected_records or len(expected_records) > MAX_SOURCE_RECORDS:
        raise _fail("source descriptor record count is outside limits")
    seen: set[str] = set()
    for index, record in enumerate(expected_records):
        parsed = _keys(
            record,
            {"record_id", "bases", "sequence_sha256", "refget_id"},
            f"source record {index}",
        )
        record_id = _text(
            parsed["record_id"],
            f"source record {index}.record_id",
            identifier=True,
        )
        if record_id in seen:
            raise _fail("source descriptor records contain duplicate record identifiers")
        seen.add(record_id)
        _integer(parsed["bases"], f"source record {index}.bases", minimum=1)
        _sha(parsed["sequence_sha256"], f"source record {index}.sequence_sha256")
        if (
            type(parsed["refget_id"]) is not str
            or REFGET_RE.fullmatch(parsed["refget_id"]) is None
        ):
            raise _fail(f"source record {index}.refget_id is invalid")
    if not _exact_equal(source_ir["records"], expected_records):
        raise _fail("source descriptor record catalog differs from native artifacts")
    source_ir_sha = _sha(
        descriptor["source_ir_sha256"],
        "source descriptor.source_ir_sha256",
    )
    if source_ir_sha != _digest(source_ir):
        raise _fail("source descriptor.source_ir_sha256 does not match source_ir")
    _artifact_seal(descriptor, "source descriptor")
    records = {record["record_id"]: record["bases"] for record in expected_records}
    return descriptor, artifacts, expected_inputs, records


def _child_reference(role: str, artifact: dict[str, Any]) -> dict[str, Any]:
    expected_format, expected_version = CHILD_IDENTITIES[role]
    if (
        artifact.get("format") != expected_format
        or type(artifact.get("version")) is not int
        or artifact["version"] != expected_version
    ):
        raise _fail(f"development child {role} has the wrong format/version")
    return {
        "format": expected_format,
        "version": expected_version,
        "artifact_sha256": _sha(
            artifact.get("artifact_sha256"),
            f"development child {role}.artifact_sha256",
        ),
    }


def _blob_references(
    response: dict[str, Any],
    module: dict[str, Any],
) -> list[dict[str, Any]]:
    if type(response.get("outputs")) is not list:
        raise _fail("prediction response.outputs must be an array")
    module_body = module.get("module")
    if type(module_body) is not dict or type(module_body.get("tensors")) is not list:
        raise _fail("development module.module.tensors must be an array")
    locations = (
        (response["outputs"], "prediction response.outputs"),
        (module_body["tensors"], "development module.module.tensors"),
    )
    found: dict[str, int] = {}
    total = 0
    for tensors, label in locations:
        if len(tensors) > validator_v2.DEFAULT_LIMITS.tensors:
            raise _fail(
                f"{label} exceeds tensor reference ceiling "
                f"{validator_v2.DEFAULT_LIMITS.tensors}"
            )
        for index, tensor in enumerate(tensors):
            if type(tensor) is not dict or type(tensor.get("storage")) is not dict:
                raise _fail(f"{label}[{index}].storage must be an object")
            storage = tensor["storage"]
            kind = storage.get("kind")
            if kind == "inline-base64":
                continue
            if kind != "sha256-blob":
                raise _fail(f"{label}[{index}].storage.kind is unsupported")
            parsed = _keys(
                storage,
                {"kind", "byte_length", "sha256"},
                f"{label}[{index}].storage",
            )
            sha = _sha(parsed["sha256"], f"{label}[{index}].storage.sha256")
            length = _integer(
                parsed["byte_length"],
                f"{label}[{index}].storage.byte_length",
                maximum=validator_v2.DEFAULT_LIMITS.tensor_bytes,
            )
            previous = found.get(sha)
            if previous is not None and previous != length:
                raise _fail("one external blob digest has conflicting byte lengths")
            if previous is None:
                found[sha] = length
                total += length
                if len(found) > MAX_EXTERNAL_BLOBS:
                    raise _fail(
                        f"external blob set exceeds count ceiling {MAX_EXTERNAL_BLOBS}"
                    )
                if total > MAX_EXTERNAL_BLOB_BYTES:
                    raise _fail(
                        "external blob set exceeds cumulative declared byte ceiling "
                        f"{MAX_EXTERNAL_BLOB_BYTES}"
                    )
    return [
        {"sha256": sha, "byte_length": found[sha]}
        for sha in sorted(found)
    ]


def _read_external_blobs(
    blob_root: str | Path | None,
    references: list[dict[str, Any]],
) -> dict[str, bytes]:
    if not references:
        return {}
    if blob_root is None:
        raise _fail("external blobs require an explicit blob root")
    if not _HAS_DIRECTORY_DESCRIPTOR:
        raise _fail("secure external blob input is unsupported on this platform")

    root_descriptor = -1
    try:
        root = Path(blob_root)
        inspected_root = root.lstat()
        if stat.S_ISLNK(inspected_root.st_mode) or not stat.S_ISDIR(
            inspected_root.st_mode
        ):
            raise _fail("blob root must be a regular non-linked directory")
        flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW
        flags |= getattr(os, "O_CLOEXEC", 0)
        root_descriptor = os.open(root, flags)
        opened_root = os.fstat(root_descriptor)
        root_identity = (opened_root.st_dev, opened_root.st_ino)
        if not stat.S_ISDIR(opened_root.st_mode) or root_identity != (
            inspected_root.st_dev,
            inspected_root.st_ino,
        ):
            raise _fail("blob root changed while being opened")

        result: dict[str, bytes] = {}
        child_flags = os.O_RDONLY | os.O_NOFOLLOW
        child_flags |= getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NONBLOCK", 0)
        for index, reference in enumerate(references):
            parsed = _keys(
                reference,
                {"sha256", "byte_length"},
                f"external blob reference {index}",
            )
            sha = _sha(parsed["sha256"], f"external blob reference {index}.sha256")
            length = _integer(
                parsed["byte_length"],
                f"external blob reference {index}.byte_length",
                maximum=validator_v2.DEFAULT_LIMITS.tensor_bytes,
            )
            before = os.stat(sha, dir_fd=root_descriptor, follow_symlinks=False)
            if (
                stat.S_ISLNK(before.st_mode)
                or not stat.S_ISREG(before.st_mode)
                or before.st_nlink != 1
            ):
                raise _fail(f"external blob {sha} must be one regular non-linked file")
            if before.st_size != length:
                raise _fail(f"external blob {sha} length differs from its descriptor")
            descriptor = os.open(sha, child_flags, dir_fd=root_descriptor)
            try:
                opened = os.fstat(descriptor)
                if (opened.st_dev, opened.st_ino) != (before.st_dev, before.st_ino):
                    raise _fail(f"external blob {sha} changed while being opened")
                raw = _read_descriptor(
                    descriptor,
                    opened,
                    f"external blob {sha}",
                    length,
                )
            finally:
                os.close(descriptor)
            after = os.stat(sha, dir_fd=root_descriptor, follow_symlinks=False)
            stable = (
                "st_dev",
                "st_ino",
                "st_nlink",
                "st_size",
                "st_mtime_ns",
                "st_ctime_ns",
            )
            if stat.S_ISLNK(after.st_mode) or any(
                getattr(before, field) != getattr(after, field) for field in stable
            ):
                raise _fail(f"external blob {sha} changed while being read")
            if hashlib.sha256(raw).hexdigest() != sha:
                raise _fail(f"external blob {sha} digest mismatch")
            result[sha] = raw

        attached = root.lstat()
        finished = os.fstat(root_descriptor)
        if (
            stat.S_ISLNK(attached.st_mode)
            or not stat.S_ISDIR(attached.st_mode)
            or (attached.st_dev, attached.st_ino) != root_identity
            or (finished.st_dev, finished.st_ino) != root_identity
        ):
            raise _fail("blob root changed while being read")
        return result
    except DevelopmentValidationError:
        raise
    except (OSError, TypeError, ValueError) as failure:
        raise _fail(f"cannot read external blob set: {failure}") from failure
    finally:
        if root_descriptor >= 0:
            os.close(root_descriptor)


def _validate_v2_chain(
    descriptor: dict[str, Any],
    records: dict[str, int],
    artifacts: dict[str, dict[str, Any]],
    *,
    blob_root: str | Path | None,
    blob_cache: Any,
) -> None:
    try:
        manifest = artifacts["provider_manifest"]
        validator_v2.validate_manifest_artifact(manifest)
        request = artifacts["prediction_request"]
        validator_v2.validate_request_artifact(request)
        response, response_tensors = validator_v2.validate_response_artifact(
            artifacts["prediction_response"],
            sequence_sha256=descriptor["artifact_sha256"],
            records=records,
            blob_root=blob_root,
            _blob_cache=blob_cache,
        )
        policy = artifacts["lowering_policy"]
        validator_v2.validate_policy_artifact(policy)
        target = artifacts["target_contract"]
        validator_v2.validate_target_artifact(target)
        module = artifacts["development_module"]
        validator_v2.validate_module_artifact(
            module,
            target,
            blob_root=blob_root,
            _blob_cache=blob_cache,
        )
        validator_v2._chain_replay(
            descriptor,
            manifest,
            request,
            response,
            response_tensors,
            policy,
            target,
            module,
            blob_root=blob_root,
            blob_cache=blob_cache,
            limits=validator_v2.DEFAULT_LIMITS,
        )
    except (OSError, ValueError) as failure:
        raise _fail(f"independent v2 compilation replay failed: {failure}") from failure


def _validate_compilation_record(
    value: Any,
    descriptor: dict[str, Any],
    artifacts: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    record = _keys(
        value,
        {
            "format",
            "version",
            "profile",
            "producer",
            "source",
            "artifacts",
            "result",
            "artifact_sha256",
        },
        "development compilation record",
    )
    _validate_tree(record, "development compilation record")
    if (
        record["format"] != COMPILATION_FORMAT
        or type(record["version"]) is not int
        or record["version"] != COMPILATION_VERSION
        or record["profile"] != PROFILE
        or not _exact_equal(record["producer"], BUNDLE_PRODUCER)
    ):
        raise _fail("development compilation record identity is unsupported")
    _artifact_seal(record, "development compilation record")
    module = artifacts["development_module"]
    chain_roles = set(CHILD_FILENAMES) - {"compilation_record"}
    expected = _seal(
        {
            "format": COMPILATION_FORMAT,
            "version": COMPILATION_VERSION,
            "profile": PROFILE,
            "producer": BUNDLE_PRODUCER,
            "source": _source_reference(descriptor),
            "artifacts": {
                role: _child_reference(role, artifacts[role])
                for role in sorted(chain_roles)
            },
            "result": {
                "module_sha256": module["module_sha256"],
                "budgets": module["module"]["budgets"],
            },
        }
    )
    if not _exact_equal(record, expected):
        raise _fail(
            "development compilation record differs from independent replay at "
            + _first_difference(record, expected)
        )
    return record


@_bounded
def validate_development_bundle(
    bundle_value: Any,
    artifacts_value: Any,
    source_descriptor: Any,
    native_artifacts: Any,
    *,
    source_inputs: Any,
    native_artifact_paths: Any = None,
    blob_root: str | Path | None = None,
) -> dict[str, Any]:
    """Replay one compiler-only bundle from its original source evidence."""

    bundle = _keys(
        bundle_value,
        {
            "format",
            "version",
            "profile",
            "source",
            "artifacts",
            "blobs",
            "artifact_sha256",
        },
        "development bundle",
    )
    _validate_tree(bundle, "development bundle")
    if (
        bundle["format"] != BUNDLE_FORMAT
        or type(bundle["version"]) is not int
        or bundle["version"] != BUNDLE_VERSION
        or bundle["profile"] != PROFILE
    ):
        raise _fail("development bundle identity is unsupported")
    _artifact_seal(bundle, "development bundle")
    descriptor, _natives, input_references, records = _validate_source_descriptor(
        source_descriptor,
        native_artifacts,
        source_inputs,
        native_artifact_paths,
    )
    artifacts = _keys(
        artifacts_value,
        set(CHILD_FILENAMES),
        "development child artifacts",
    )
    for role, artifact in artifacts.items():
        if type(artifact) is not dict:
            raise _fail(f"development child {role} must be an object")
        _validate_tree(artifact, f"development child {role}")
        _child_reference(role, artifact)

    expected_blobs = _blob_references(
        artifacts["prediction_response"],
        artifacts["development_module"],
    )
    expected_bundle = _seal(
        {
            "format": BUNDLE_FORMAT,
            "version": BUNDLE_VERSION,
            "profile": PROFILE,
            "source": {
                "descriptor": _source_reference(descriptor),
                "native": descriptor["source_ir"]["artifacts"],
            },
            "artifacts": {
                role: _child_reference(role, artifacts[role])
                for role in sorted(CHILD_FILENAMES)
            },
            "blobs": expected_blobs,
        }
    )
    if not _exact_equal(bundle, expected_bundle):
        raise _fail(
            "development bundle differs from independent replay at "
            + _first_difference(bundle, expected_bundle)
        )

    blob_bytes = _read_external_blobs(blob_root, expected_blobs)
    blob_cache = validator_v2._validated_blob_cache(blob_bytes)
    _validate_v2_chain(
        descriptor,
        records,
        artifacts,
        blob_root=blob_root,
        blob_cache=blob_cache,
    )
    compilation = _validate_compilation_record(
        artifacts["compilation_record"],
        descriptor,
        artifacts,
    )

    report_core = {
        "format": REPORT_FORMAT,
        "version": REPORT_VERSION,
        "stage": "development-compilation",
        "valid": True,
        "inputs": {
            "source_profile": descriptor["source_ir"]["profile"],
            "source_descriptor_sha256": descriptor["artifact_sha256"],
            "source_ir_sha256": descriptor["source_ir_sha256"],
            "source_inputs": {
                role: dict(input_references[role])
                for role in sorted(input_references)
            },
            "bundle_artifact_sha256": bundle["artifact_sha256"],
        },
        "result": {
            "compilation_artifact_sha256": compilation["artifact_sha256"],
            "module_artifact_sha256": artifacts["development_module"][
                "artifact_sha256"
            ],
            "module_sha256": artifacts["development_module"]["module_sha256"],
            "budgets": artifacts["development_module"]["module"]["budgets"],
            "blobs": expected_blobs,
        },
        "checks": list(REPORT_CHECKS),
    }
    return {**report_core, "report_sha256": _digest(report_core)}


@_bounded
def validate_development_report(report_value: Any) -> dict[str, Any]:
    """Validate the exact sealed success-report contract."""

    if type(report_value) is not dict:
        raise _fail("development validation report must be an object")
    success_keys = {
        "format",
        "version",
        "stage",
        "valid",
        "inputs",
        "result",
        "checks",
        "report_sha256",
    }
    failure_keys = {
        "format",
        "version",
        "stage",
        "valid",
        "error",
        "report_sha256",
    }
    expected_keys = success_keys if report_value.get("valid") is True else failure_keys
    report = _keys(
        report_value,
        expected_keys,
        "development validation report",
    )
    _validate_tree(report, "development validation report")
    if (
        report["format"] != REPORT_FORMAT
        or type(report["version"]) is not int
        or report["version"] != REPORT_VERSION
        or report["stage"] != "development-compilation"
        or type(report["valid"]) is not bool
    ):
        raise _fail("development validation report identity is invalid")
    if report["valid"] is False:
        error = _keys(
            report["error"],
            {"code", "message"},
            "development validation report.error",
        )
        if error["code"] != "DEVVAL001":
            raise _fail("development validation report.error.code is invalid")
        _text(error["message"], "development validation report.error.message")
        _validate_report_seal(report)
        return report
    if report["checks"] != list(REPORT_CHECKS):
        raise _fail("development validation report checks are invalid")
    inputs = _keys(
        report["inputs"],
        {
            "source_profile",
            "source_descriptor_sha256",
            "source_ir_sha256",
            "source_inputs",
            "bundle_artifact_sha256",
        },
        "development validation report.inputs",
    )
    if inputs["source_profile"] not in PROFILES:
        raise _fail("development validation report source profile is unsupported")
    for name in (
        "source_descriptor_sha256",
        "source_ir_sha256",
        "bundle_artifact_sha256",
    ):
        _sha(inputs[name], f"development validation report.inputs.{name}")
    if inputs["source_profile"] == EXTERNAL_PROFILE:
        source_hashes = inputs["source_inputs"]
        if (
            type(source_hashes) is not dict
            or not source_hashes
            or len(source_hashes) > validator_external.MAX_ROLES
            or any(
                type(role) is not str or EXTERNAL_ROLE_RE.fullmatch(role) is None
                for role in source_hashes
            )
        ):
            raise _fail(
                "development validation report external source roles are invalid"
            )
    else:
        input_roles = _report_input_roles(inputs["source_profile"])
        source_hashes = _keys(
            inputs["source_inputs"],
            set(input_roles),
            "development validation report.inputs.source_inputs",
        )
    for role, value in source_hashes.items():
        reference = _keys(
            value,
            {"sha256", "byte_length"},
            f"development validation report source input {role}",
        )
        _sha(reference["sha256"], f"development validation report source input {role}")
        _integer(
            reference["byte_length"],
            f"development validation report source input {role}.byte_length",
        )
    result = _keys(
        report["result"],
        {
            "compilation_artifact_sha256",
            "module_artifact_sha256",
            "module_sha256",
            "budgets",
            "blobs",
        },
        "development validation report.result",
    )
    for name in (
        "compilation_artifact_sha256",
        "module_artifact_sha256",
        "module_sha256",
    ):
        _sha(result[name], f"development validation report.result.{name}")
    budgets = _keys(
        result["budgets"],
        {"operations", "tensor_bytes", "units", "edges", "attachments"},
        "development validation report.result.budgets",
    )
    for name, value in budgets.items():
        _integer(value, f"development validation report.result.budgets.{name}")
    if type(result["blobs"]) is not list:
        raise _fail("development validation report.result.blobs must be an array")
    if len(result["blobs"]) > MAX_EXTERNAL_BLOBS:
        raise _fail("development validation report exceeds external blob count ceiling")
    blob_digests: list[str] = []
    blob_bytes = 0
    for index, blob in enumerate(result["blobs"]):
        parsed = _keys(blob, {"sha256", "byte_length"}, f"report blob {index}")
        blob_digests.append(_sha(parsed["sha256"], f"report blob {index}.sha256"))
        blob_bytes += _integer(
            parsed["byte_length"],
            f"report blob {index}.byte_length",
            maximum=validator_v2.DEFAULT_LIMITS.tensor_bytes,
        )
    if blob_digests != sorted(blob_digests) or len(blob_digests) != len(
        set(blob_digests)
    ):
        raise _fail("development validation report blobs must be sorted and unique")
    if blob_bytes > MAX_EXTERNAL_BLOB_BYTES:
        raise _fail("development validation report exceeds external blob byte ceiling")
    _validate_report_seal(report)
    return report


def _validate_report_seal(report: dict[str, Any]) -> None:
    expected = _digest(
        {key: value for key, value in report.items() if key != "report_sha256"}
    )
    if _sha(
        report["report_sha256"],
        "development validation report.report_sha256",
    ) != expected:
        raise _fail("development validation report seal is invalid")


def _report_input_roles(profile: str) -> tuple[str, ...]:
    if profile == GENBANK_PROFILE:
        return ("genbank",)
    if profile == GFF3_PROFILE:
        return ("sequence", "annotation")
    return ("sequence",)


@_bounded
def load_source_bundle_directory(
    path: str | Path,
) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    """Load a source descriptor and its exact profile-selected native children."""

    raw = _read_directory(
        path,
        SOURCE_FILE_LIMITS,
        required={SOURCE_FILENAME},
        label="source bundle",
    )
    descriptor = _json_object(raw[SOURCE_FILENAME], "source descriptor")
    source_ir = _keys(
        descriptor.get("source_ir"),
        {"profile", "parameters", "inputs", "artifacts", "records"},
        "source descriptor.source_ir",
    )
    if source_ir["profile"] == EXTERNAL_PROFILE:
        _external_roles(
            source_ir["inputs"],
            "external source descriptor inputs",
        )
    _, artifact_roles = _route(source_ir["profile"], source_ir["parameters"])
    expected_names = {SOURCE_FILENAME, *(SOURCE_FILENAMES[role] for role in artifact_roles)}
    if set(raw) != expected_names:
        raise _fail("source bundle filename closure does not match its source profile")
    artifacts = {
        role: _json_object(
            raw[SOURCE_FILENAMES[role]],
            f"native source artifact {role}",
            maximum_members=(
                MAX_GFF3_JSON_MEMBERS if role == "annotation" else MAX_JSON_MEMBERS
            ),
        )
        for role in artifact_roles
    }
    return descriptor, artifacts


@_bounded
def load_development_bundle_directory(
    path: str | Path,
) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    """Load the exact flat index-and-seven-child development bundle."""

    limits = {filename: MAX_JSON_BYTES for filename in BUNDLE_FILENAMES}
    raw = _read_directory(
        path,
        limits,
        required=set(BUNDLE_FILENAMES),
        label="development bundle",
    )
    bundle = _json_object(raw[BUNDLE_FILENAME], "development bundle index")
    artifacts = {
        role: _json_object(raw[filename], f"development child {role}")
        for role, filename in CHILD_FILENAMES.items()
    }
    return bundle, artifacts


@_bounded
def validate_development_paths(
    source_bundle: str | Path,
    development_bundle: str | Path,
    *,
    source_inputs: dict[str, str | Path],
    native_artifact_paths: dict[str, str | Path] | None = None,
    blob_root: str | Path | None = None,
) -> dict[str, Any]:
    """Snapshot both closed bundles and replay compilation from source paths."""

    descriptor, natives = load_source_bundle_directory(source_bundle)
    source_ir = _keys(
        descriptor.get("source_ir"),
        {"profile", "parameters", "inputs", "artifacts", "records"},
        "source descriptor.source_ir",
    )
    input_roles, _ = _route(source_ir["profile"], source_ir["parameters"])
    if source_ir["profile"] == EXTERNAL_PROFILE:
        expected_input_roles = _external_roles(
            source_ir["inputs"],
            "external source descriptor inputs",
        )
    else:
        expected_input_roles = set(input_roles)
    paths = _keys(
        source_inputs,
        expected_input_roles,
        "original source input paths",
    )
    if source_ir["profile"] in {REFERENCE_FASTA_PROFILE, EXTERNAL_PROFILE}:
        original: dict[str, Any] = dict(paths)
    else:
        original = {
            role: _read_regular(
                paths[role],
                f"original source input {role}",
                MAX_GFF3_BYTES if role == "annotation" else MAX_INPUT_BYTES,
            )
            for role in input_roles
        }
    bundle, artifacts = load_development_bundle_directory(development_bundle)
    return validate_development_bundle(
        bundle,
        artifacts,
        descriptor,
        natives,
        source_inputs=original,
        native_artifact_paths=native_artifact_paths,
        blob_root=blob_root,
    )


def _failure_report(failure: DevelopmentValidationError) -> dict[str, Any]:
    message = str(failure)
    try:
        encoded = message.encode("utf-8")
    except UnicodeEncodeError:
        encoded = message.encode("utf-8", "backslashreplace")
    suffix = b"... [truncated]"
    if len(encoded) > MAX_FAILURE_MESSAGE_BYTES:
        encoded = encoded[: MAX_FAILURE_MESSAGE_BYTES - len(suffix)] + suffix
    message = encoded.decode("utf-8", "ignore")
    core = {
        "format": REPORT_FORMAT,
        "version": REPORT_VERSION,
        "stage": "development-compilation",
        "valid": False,
        "error": {
            "code": "DEVVAL001",
            "message": message,
        },
    }
    return {**core, "report_sha256": _digest(core)}


def _source_input_argument(value: str) -> tuple[str, Path]:
    role, separator, raw_path = value.partition("=")
    if (
        not separator
        or EXTERNAL_ROLE_RE.fullmatch(role) is None
        or not raw_path
    ):
        raise argparse.ArgumentTypeError(
            "input must be ROLE=PATH with a portable role identifier"
        )
    return role, Path(raw_path)


def main(argv: list[str] | None = None) -> int:
    """Validate one source/development pair and emit a sealed JSON report."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_bundle", type=Path)
    parser.add_argument("development_bundle", type=Path)
    parser.add_argument(
        "--source-input",
        action="append",
        type=_source_input_argument,
        required=True,
        metavar="ROLE=PATH",
        help="original source path; repeat for each descriptor input role",
    )
    parser.add_argument(
        "--native-input",
        action="append",
        type=_source_input_argument,
        default=[],
        metavar="ROLE=PATH",
        help="external native artifact path; repeat for each declared native role",
    )
    parser.add_argument("--blob-root", type=Path)
    parser.add_argument("-o", "--output", type=Path)
    arguments = parser.parse_args(argv)

    source_inputs: dict[str, Path] = {}
    for role, path in arguments.source_input:
        if role in source_inputs:
            parser.error(f"source input role {role!r} was supplied more than once")
        source_inputs[role] = path
    native_inputs: dict[str, Path] = {}
    for role, path in arguments.native_input:
        if role in native_inputs:
            parser.error(f"native input role {role!r} was supplied more than once")
        native_inputs[role] = path
    try:
        report = validate_development_paths(
            arguments.source_bundle,
            arguments.development_bundle,
            source_inputs=source_inputs,
            native_artifact_paths=native_inputs or None,
            blob_root=arguments.blob_root,
        )
    except DevelopmentValidationError as failure:
        report = _failure_report(failure)

    if arguments.output is not None:
        try:
            validator_v2.save_report(report, arguments.output)
        except (OSError, ValueError) as failure:
            report = _failure_report(
                _fail(f"cannot write validation report: {failure}")
            )
    print(
        json.dumps(
            report,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        )
    )
    return 0 if report["valid"] else 1


__all__ = [
    "BUNDLE_FILENAME",
    "CHILD_FILENAMES",
    "DevelopmentValidationError",
    "REPORT_FORMAT",
    "load_development_bundle_directory",
    "load_source_bundle_directory",
    "main",
    "validate_development_bundle",
    "validate_development_paths",
    "validate_development_report",
]


if __name__ == "__main__":
    raise SystemExit(main())

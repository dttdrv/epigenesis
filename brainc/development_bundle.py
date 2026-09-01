"""Reference-only packaging for the existing typed development compiler."""

from __future__ import annotations

import copy
from pathlib import Path
import tempfile
from typing import Any

from ._canonical import ContractError, artifact_digest, digest
from ._io import BoundedIOError, MAX_JSON_BYTES, atomic_write_file, read_regular_file
from ._publish import PublicationError, publish_directory
from .source import (
    FORMAT as SOURCE_FORMAT,
    VERSION as SOURCE_VERSION,
    SourceBundle,
    SourceError,
    load_source_bundle,
    validate_source_bundle,
)
from .v2._common import V2Error, loads, pretty_bytes
from .v2.compiler import _compile_module_from_validated_source_bundle
from .v2.policy import load_policy
from .v2.provider import load_manifest, load_request, load_response
from .v2.target import load_target


BUNDLE_FORMAT = "brainc.development-bundle"
BUNDLE_VERSION = 1
COMPILATION_FORMAT = "brainc.development-compilation"
COMPILATION_VERSION = 1
PROFILE = "development-module/v1"
PRODUCER = {
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
BUNDLE_FILENAME = "bundle.json"
_CHAIN_INPUTS = (
    ("provider_manifest", "manifest", "provider manifest"),
    ("prediction_request", "request", "prediction request"),
    ("prediction_response", "response", "prediction response"),
    ("lowering_policy", "policy", "lowering policy"),
    ("target_contract", "target", "target contract"),
)
_IDENTITIES = {
    "development_module": ("brainc.development-module", 1),
    "lowering_policy": ("brainc.lowering-policy", 2),
    "prediction_request": ("brainc.prediction-request", 2),
    "prediction_response": ("brainc.prediction-response", 2),
    "provider_manifest": ("brainc.provider-manifest", 2),
    "target_contract": ("brainc.target-contract", 1),
    "compilation_record": (COMPILATION_FORMAT, COMPILATION_VERSION),
}


class DevelopmentBundleError(ContractError):
    """The source, interpretation chain, or development bundle is invalid."""


def _fail(detail: str) -> DevelopmentBundleError:
    return DevelopmentBundleError(f"DEVB001: {detail}")


def _source(value: str | Path | SourceBundle) -> SourceBundle:
    if isinstance(value, SourceBundle):
        return validate_source_bundle(value)
    try:
        return load_source_bundle(value)
    except (SourceError, OSError, TypeError) as failure:
        raise _fail(f"invalid source bundle: {failure}") from failure


def _snapshot(path: str | Path, label: str) -> bytes:
    try:
        raw = read_regular_file(path, maximum_bytes=MAX_JSON_BYTES, label=label)
        loads(raw, label, maximum_bytes=MAX_JSON_BYTES)
        return raw
    except (BoundedIOError, V2Error, OSError, TypeError) as failure:
        raise _fail(f"cannot snapshot {label}: {failure}") from failure


def _write_snapshot(path: Path, raw: bytes, label: str) -> None:
    try:
        atomic_write_file(path, raw, maximum_bytes=MAX_JSON_BYTES, label=label)
    except BoundedIOError as failure:
        raise _fail(f"cannot stage {label}: {failure}") from failure


def _reference(role: str, payload: dict[str, Any]) -> dict[str, Any]:
    expected_format, expected_version = _IDENTITIES[role]
    if (
        payload.get("format") != expected_format
        or type(payload.get("version")) is not int
        or payload["version"] != expected_version
    ):
        raise _fail(f"{role} has the wrong format/version")
    try:
        artifact_sha = artifact_digest(payload)
    except (ContractError, KeyError) as failure:
        raise _fail(f"{role} artifact digest is invalid: {failure}") from failure
    return {
        "format": expected_format,
        "version": expected_version,
        "artifact_sha256": artifact_sha,
    }


def _source_reference(descriptor: dict[str, Any]) -> dict[str, Any]:
    return {
        "format": SOURCE_FORMAT,
        "version": SOURCE_VERSION,
        "artifact_sha256": descriptor["artifact_sha256"],
        "ir_sha256": descriptor["source_ir_sha256"],
    }


def _blob_references(*artifacts: dict[str, Any]) -> list[dict[str, Any]]:
    found: dict[str, int] = {}
    for artifact in artifacts:
        if artifact["format"] == "brainc.prediction-response":
            tensors = artifact["outputs"]
        elif artifact["format"] == "brainc.development-module":
            tensors = artifact["module"]["tensors"]
        else:
            continue
        for tensor in tensors:
            storage = tensor["storage"]
            if storage["kind"] != "sha256-blob":
                continue
            previous = found.setdefault(storage["sha256"], storage["byte_length"])
            if previous != storage["byte_length"]:
                raise _fail("one external blob digest has conflicting byte lengths")
    return [
        {"sha256": blob_sha, "byte_length": found[blob_sha]}
        for blob_sha in sorted(found)
    ]


def _seal(core: dict[str, Any]) -> dict[str, Any]:
    payload = {**core, "artifact_sha256": digest(core)}
    pretty_bytes(payload)
    return payload


class DevelopmentBundle:
    """Immutable-by-interface flat index and its seven JSON children."""

    __slots__ = ("_bundle", "_artifacts")

    def __init__(self, *_: Any, **__: Any) -> None:
        raise TypeError("DevelopmentBundle instances are created by compile_development")

    @classmethod
    def _create(
        cls,
        bundle: dict[str, Any],
        artifacts: dict[str, dict[str, Any]],
    ) -> DevelopmentBundle:
        instance = object.__new__(cls)
        instance._bundle = copy.deepcopy(bundle)
        instance._artifacts = copy.deepcopy(artifacts)
        return instance

    @property
    def artifacts(self) -> dict[str, dict[str, Any]]:
        return copy.deepcopy(self._artifacts)

    def to_dict(self) -> dict[str, Any]:
        return copy.deepcopy(self._bundle)

    def save(self, directory: str | Path) -> dict[str, Path]:
        entries = {BUNDLE_FILENAME: self._bundle}
        entries.update(
            {
                CHILD_FILENAMES[role]: artifact
                for role, artifact in self._artifacts.items()
            }
        )
        try:
            return publish_directory(directory, entries, pretty_bytes)
        except PublicationError as failure:
            state = " after commit" if failure.committed else ""
            raise _fail(
                f"cannot publish development bundle{state}: {failure.strerror}"
            ) from failure


def compile_development(
    source_bundle: str | Path | SourceBundle,
    manifest: str | Path,
    request: str | Path,
    response: str | Path,
    policy: str | Path,
    target: str | Path,
    *,
    blob_root: str | Path | None = None,
) -> DevelopmentBundle:
    """Compile one stable caller-supplied v2 chain and index its exact closure."""

    source = _source(source_bundle)
    descriptor = source.to_dict()
    supplied = {
        "manifest": manifest,
        "request": request,
        "response": response,
        "policy": policy,
        "target": target,
    }
    snapshots: dict[str, bytes] = {}
    for _, argument, label in _CHAIN_INPUTS:
        snapshots[argument] = _snapshot(supplied[argument], label)

    try:
        with tempfile.TemporaryDirectory(prefix="brainc-development-") as temporary:
            root = Path(temporary)
            paths: dict[str, Path] = {}
            for _, argument, label in _CHAIN_INPUTS:
                path = root / f"{argument}.json"
                _write_snapshot(path, snapshots[argument], f"{label} snapshot")
                paths[argument] = path

            manifest_payload = load_manifest(paths["manifest"])
            request_payload = load_request(paths["request"])
            response_payload = load_response(paths["response"])
            policy_payload = load_policy(paths["policy"])
            target_payload = load_target(paths["target"]).artifact
            module = _compile_module_from_validated_source_bundle(
                source,
                paths["manifest"],
                paths["request"],
                paths["response"],
                paths["policy"],
                paths["target"],
                blob_root=blob_root,
            )
    except (V2Error, OSError) as failure:
        raise _fail(f"development compilation failed: {failure}") from failure

    chain = {
        "development_module": module,
        "lowering_policy": policy_payload,
        "prediction_request": request_payload,
        "prediction_response": response_payload,
        "provider_manifest": manifest_payload,
        "target_contract": target_payload,
    }
    chain_references = {
        role: _reference(role, artifact) for role, artifact in chain.items()
    }
    source_reference = _source_reference(descriptor)
    compilation_record = _seal(
        {
            "format": COMPILATION_FORMAT,
            "version": COMPILATION_VERSION,
            "profile": PROFILE,
            "producer": copy.deepcopy(PRODUCER),
            "source": source_reference,
            "artifacts": copy.deepcopy(chain_references),
            "result": {
                "module_sha256": module["module_sha256"],
                "budgets": copy.deepcopy(module["module"]["budgets"]),
            },
        }
    )
    artifacts = {"compilation_record": compilation_record, **chain}
    bundle_references = {
        role: _reference(role, artifact) for role, artifact in artifacts.items()
    }
    bundle = _seal(
        {
            "format": BUNDLE_FORMAT,
            "version": BUNDLE_VERSION,
            "profile": PROFILE,
            "source": {
                "descriptor": source_reference,
                "native": copy.deepcopy(descriptor["source_ir"]["artifacts"]),
            },
            "artifacts": bundle_references,
            "blobs": _blob_references(response_payload, module),
        }
    )
    return DevelopmentBundle._create(bundle, artifacts)


__all__ = [
    "BUNDLE_FILENAME",
    "BUNDLE_FORMAT",
    "BUNDLE_VERSION",
    "CHILD_FILENAMES",
    "COMPILATION_FORMAT",
    "COMPILATION_VERSION",
    "DevelopmentBundle",
    "DevelopmentBundleError",
    "PROFILE",
    "compile_development",
]

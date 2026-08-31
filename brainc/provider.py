"""Offline file ABI between brainc and an independently operated predictor."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ._canonical import (
    ContractError, artifact_digest, digest, keys, load, number, save_artifact, sha256, text
)
from .sequence import _as_artifact as _as_sequence_artifact
from .sequence_collection import _as_collection as _as_sequence_collection


class ProviderError(ContractError):
    """A provider manifest, request, response, or binding is invalid."""


VALUE_TYPES = {"number", "integer", "boolean"}
SOURCE_FORMATS = {
    ("brain01.sequence-ir", 2): "ir_sha256",
    ("brain01.sequence-collection-ir", 1): "collection_ir_sha256",
}


def _checked_artifact(path: str | Path, label: str) -> tuple[dict[str, Any], bytes]:
    payload, raw = load(path, label)
    artifact_digest(payload)
    return payload, raw


def load_source(path: str | Path) -> dict[str, Any]:
    payload, _ = load(path, "sequence source")
    identity = (payload.get("format"), payload.get("version"))
    if identity == ("brain01.sequence-ir", 2):
        _as_sequence_artifact(payload)
    elif identity == ("brain01.sequence-collection-ir", 1):
        _as_sequence_collection(payload)
    else:
        raise ProviderError("unsupported sequence source format or version")
    return payload


def _source_binding(payload: dict[str, Any]) -> dict[str, Any]:
    identity = (payload["format"], payload["version"])
    ir_field = SOURCE_FORMATS.get(identity)
    if ir_field is None:
        raise ProviderError("unsupported sequence source format or version")
    return {
        "format": payload["format"],
        "version": payload["version"],
        "artifact_sha256": payload["artifact_sha256"],
        "ir_sha256": payload[ir_field],
    }


def _output_contract(value: Any, label: str) -> dict[str, Any]:
    item = keys(value, {"id", "type", "unit"}, label)
    output_id = text(item["id"], f"{label}.id")
    output_type = item["type"]
    if type(output_type) is not str or output_type not in VALUE_TYPES:
        raise ProviderError(f"{label}.type is unsupported")
    unit = item["unit"]
    if unit is not None:
        unit = text(unit, f"{label}.unit")
    return {"id": output_id, "type": output_type, "unit": unit}


def load_manifest(path: str | Path) -> dict[str, Any]:
    payload, _ = _checked_artifact(path, "provider manifest")
    keys(payload, {"format", "version", "provider", "model_identity", "accepts", "outputs", "artifact_sha256"}, "provider manifest")
    if payload["format"] != "brainc.provider-manifest" or type(payload["version"]) is not int or payload["version"] != 1:
        raise ProviderError("unsupported provider manifest format or version")
    provider = keys(payload["provider"], {"name", "version"}, "provider")
    text(provider["name"], "provider.name"); text(provider["version"], "provider.version")
    model = keys(payload["model_identity"], {"kind", "value"}, "model_identity")
    if model["kind"] == "content-sha256":
        sha256(model["value"], "model_identity.value")
    elif model["kind"] == "opaque":
        text(model["value"], "model_identity.value")
    else:
        raise ProviderError("model_identity.kind must be 'content-sha256' or 'opaque'")
    accepts = payload["accepts"]
    if type(accepts) is not list or not accepts:
        raise ProviderError("accepts must be a non-empty array")
    normalized_accepts = [text(item, f"accepts[{index}]") for index, item in enumerate(accepts)]
    if len(set(normalized_accepts)) != len(normalized_accepts):
        raise ProviderError("accepts must not contain duplicates")
    outputs = payload["outputs"]
    if type(outputs) is not list or not outputs:
        raise ProviderError("outputs must be a non-empty array")
    normalized = [_output_contract(item, f"outputs[{index}]") for index, item in enumerate(outputs)]
    if len({item["id"] for item in normalized}) != len(normalized):
        raise ProviderError("provider output ids must be unique")
    return payload


def make_request(source_path: str | Path, manifest_path: str | Path, output_ids: list[str]) -> dict[str, Any]:
    source = load_source(source_path)
    manifest = load_manifest(manifest_path)
    source_tag = f'{source["format"]}/v{source["version"]}'
    if source_tag not in manifest["accepts"]:
        raise ProviderError(f"provider does not declare support for {source_tag}")
    if not output_ids:
        raise ProviderError("at least one requested output id is required")
    normalized_ids = [text(output_id, "requested output id") for output_id in output_ids]
    if len(set(normalized_ids)) != len(normalized_ids):
        raise ProviderError("requested output ids must be unique")
    available = {item["id"]: item for item in manifest["outputs"]}
    requested: list[dict[str, Any]] = []
    for output_id in normalized_ids:
        if output_id not in available:
            raise ProviderError(f"provider does not declare output {output_id!r}")
        requested.append(dict(available[output_id]))
    core = {
        "format": "brainc.prediction-request", "version": 1,
        "source": _source_binding(source),
        "provider_manifest_sha256": manifest["artifact_sha256"],
        "requested_outputs": requested,
    }
    return {**core, "artifact_sha256": digest(core)}


def load_request(path: str | Path) -> dict[str, Any]:
    payload, _ = _checked_artifact(path, "prediction request")
    keys(payload, {"format", "version", "source", "provider_manifest_sha256", "requested_outputs", "artifact_sha256"}, "prediction request")
    if payload["format"] != "brainc.prediction-request" or type(payload["version"]) is not int or payload["version"] != 1:
        raise ProviderError("unsupported prediction request format or version")
    source = keys(payload["source"], {"format", "version", "artifact_sha256", "ir_sha256"}, "request source")
    identity = (source["format"], source["version"])
    if type(source["format"]) is not str or type(source["version"]) is not int or identity not in SOURCE_FORMATS:
        raise ProviderError("request source format/version is unsupported")
    sha256(source["artifact_sha256"], "request.source.artifact_sha256")
    sha256(source["ir_sha256"], "request.source.ir_sha256")
    sha256(payload["provider_manifest_sha256"], "provider_manifest_sha256")
    outputs = payload["requested_outputs"]
    if type(outputs) is not list or not outputs:
        raise ProviderError("requested_outputs must be a non-empty array")
    normalized = [_output_contract(item, f"requested_outputs[{index}]") for index, item in enumerate(outputs)]
    if len({item["id"] for item in normalized}) != len(normalized):
        raise ProviderError("requested output ids must be unique")
    return payload


def _typed_value(value: Any, kind: str, label: str) -> int | float | bool:
    if kind == "boolean":
        if type(value) is not bool:
            raise ProviderError(f"{label} must be a boolean")
    elif kind == "integer":
        if type(value) is not int:
            raise ProviderError(f"{label} must be an integer")
        number(value, label)
    else:
        number(value, label)
    return value


def load_response(path: str | Path) -> dict[str, Any]:
    payload, _ = _checked_artifact(path, "prediction response")
    keys(payload, {"format", "version", "request_artifact_sha256", "provider", "model_identity", "outputs", "artifact_sha256"}, "prediction response")
    if payload["format"] != "brainc.prediction-response" or type(payload["version"]) is not int or payload["version"] != 1:
        raise ProviderError("unsupported prediction response format or version")
    sha256(payload["request_artifact_sha256"], "request_artifact_sha256")
    provider = keys(payload["provider"], {"name", "version"}, "response provider")
    text(provider["name"], "response.provider.name"); text(provider["version"], "response.provider.version")
    model = keys(payload["model_identity"], {"kind", "value"}, "response model_identity")
    if model["kind"] == "content-sha256": sha256(model["value"], "response.model_identity.value")
    elif model["kind"] == "opaque": text(model["value"], "response.model_identity.value")
    else: raise ProviderError("response model_identity.kind is unsupported")
    outputs = payload["outputs"]
    if type(outputs) is not list or not outputs:
        raise ProviderError("response outputs must be a non-empty array")
    seen: set[str] = set()
    for index, value in enumerate(outputs):
        item = keys(value, {"id", "type", "unit", "value"}, f"response.outputs[{index}]")
        contract = _output_contract({key: item[key] for key in ("id", "type", "unit")}, f"response.outputs[{index}]")
        if contract["id"] in seen: raise ProviderError("response output ids must be unique")
        seen.add(contract["id"])
        _typed_value(item["value"], contract["type"], f"response.outputs[{index}].value")
    return payload


def validate_binding(source_path: str | Path, manifest_path: str | Path, request_path: str | Path, response_path: str | Path) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any]]:
    source = load_source(source_path); manifest = load_manifest(manifest_path)
    request = load_request(request_path); response = load_response(response_path)
    if request["source"] != _source_binding(source): raise ProviderError("request is not bound to the supplied sequence source")
    if request["provider_manifest_sha256"] != manifest["artifact_sha256"]: raise ProviderError("request is not bound to the supplied provider manifest")
    source_tag = f'{source["format"]}/v{source["version"]}'
    if source_tag not in manifest["accepts"]: raise ProviderError(f"provider does not declare support for {source_tag}")
    declared = {item["id"]: item for item in manifest["outputs"]}
    for item in request["requested_outputs"]:
        if item["id"] not in declared or item != declared[item["id"]]:
            raise ProviderError(f"request output {item['id']!r} is not declared by the provider manifest")
    if response["request_artifact_sha256"] != request["artifact_sha256"]: raise ProviderError("response is not bound to the supplied request")
    if response["provider"] != manifest["provider"] or response["model_identity"] != manifest["model_identity"]:
        raise ProviderError("response provider/model identity differs from the manifest")
    expected = request["requested_outputs"]
    observed = [{key: item[key] for key in ("id", "type", "unit")} for item in response["outputs"]]
    if observed != expected: raise ProviderError("response output contract/order differs from the request")
    return source, manifest, request, response


def save(payload: dict[str, Any], path: str | Path) -> None:
    save_artifact(payload, path)

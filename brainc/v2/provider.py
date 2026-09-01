"""Version-two tensor provider contract and request/response binding."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from brainc.insdc import (
    FORMAT as GENBANK_FORMAT,
    VERSION as GENBANK_VERSION,
    GenBankError,
    validate_genbank_artifact,
)
from brainc.sequence import SequenceCompilerError, _as_artifact
from brainc.sequence_collection import SequenceCollectionError, _as_collection
from brainc.source import (
    FORMAT as SOURCE_DESCRIPTOR_FORMAT,
    PROFILES as SOURCE_DESCRIPTOR_PROFILES,
    VERSION as SOURCE_DESCRIPTOR_VERSION,
    SourceBundle,
    SourceError,
    load_source_bundle,
    validate_source_bundle,
    validate_source_descriptor,
)

from ._common import (
    V2Error,
    _check_artifact,
    checked_artifact,
    exact_version,
    integer,
    keys,
    load,
    seal,
    sha256,
    sorted_unique,
    text,
)
from .limits import (
    MAX_JSON_BYTES,
    MAX_SOURCE_JSON_BYTES,
    MAX_SOURCE_RECORDS,
    MAX_TENSORS,
)
from .tensor import parse_axes, parse_storage, parse_type, parse_unit


SOURCE_FORMATS = {
    ("brain01.sequence-ir", 2): "ir_sha256",
    ("brain01.sequence-collection-ir", 1): "collection_ir_sha256",
    (GENBANK_FORMAT, GENBANK_VERSION): "bio_ir_sha256",
    (SOURCE_DESCRIPTOR_FORMAT, SOURCE_DESCRIPTOR_VERSION): "source_ir_sha256",
}


def _source_identity(source: dict[str, Any], label: str) -> tuple[str, int]:
    source_format = source.get("format")
    source_version = source.get("version")
    if type(source_format) is not str or type(source_version) is not int:
        raise V2Error(f"{label} format/version is unsupported")
    return source_format, source_version


def _reject_standalone_descriptor(
    source: dict[str, Any],
    *,
    use: str,
) -> None:
    if _source_identity(source, "sequence source") == (
        SOURCE_DESCRIPTOR_FORMAT,
        SOURCE_DESCRIPTOR_VERSION,
    ):
        raise V2Error(
            "standalone brainc.source-descriptor/v1 is not a complete source; "
            f"use {use} with a validated SourceBundle"
        )


def _source_binding(source: dict[str, Any]) -> dict[str, Any]:
    identity = _source_identity(source, "sequence source")
    ir_field = SOURCE_FORMATS.get(identity)
    if ir_field is None:
        raise V2Error("unsupported sequence source format/version")
    return {
        "format": source["format"],
        "version": source["version"],
        "artifact_sha256": source["artifact_sha256"],
        "ir_sha256": source[ir_field],
    }


def _source_acceptance_tag(source: dict[str, Any]) -> str:
    """Name the exact source language a provider promises to interpret."""

    identity = _source_identity(source, "sequence source")
    source_tag = f"{identity[0]}/v{identity[1]}"
    if identity != (SOURCE_DESCRIPTOR_FORMAT, SOURCE_DESCRIPTOR_VERSION):
        return source_tag

    source_ir = source.get("source_ir")
    if type(source_ir) is not dict:
        raise V2Error("source descriptor is missing its source IR")
    profile = source_ir.get("profile")
    if type(profile) is not str or profile not in SOURCE_DESCRIPTOR_PROFILES:
        raise V2Error("source descriptor profile is unsupported")
    return f"{source_tag};profile={profile}"


def load_source(path: str | Path) -> tuple[dict[str, Any], dict[str, int]]:
    source, raw = load(
        path,
        "sequence source",
        maximum_bytes=MAX_SOURCE_JSON_BYTES,
    )
    identity = _source_identity(source, "sequence source")
    if identity != (GENBANK_FORMAT, GENBANK_VERSION) and len(raw) > MAX_JSON_BYTES:
        raise V2Error(f"sequence source exceeds JSON byte limit {MAX_JSON_BYTES}")
    _check_artifact(source, "sequence source")
    try:
        if identity == ("brain01.sequence-ir", 2):
            loaded = _as_artifact(source)
            records = {loaded.record_id: len(loaded.sequence)}
        elif identity == ("brain01.sequence-collection-ir", 1):
            loaded = _as_collection(source)
            records = {
                member.artifact.record_id: len(member.artifact.sequence)
                for member in loaded.members
            }
        elif identity == (GENBANK_FORMAT, GENBANK_VERSION):
            validate_genbank_artifact(source)
            records = {
                member["record_id"]: member["sequence"]["bases"]
                for member in source["sequence_collection"]["members"]
            }
        elif identity == (SOURCE_DESCRIPTOR_FORMAT, SOURCE_DESCRIPTOR_VERSION):
            descriptor = validate_source_descriptor(source)
            records = {
                member["record_id"]: member["bases"]
                for member in descriptor["source_ir"]["records"]
            }
        else:
            raise V2Error("unsupported sequence source format/version")
    except (GenBankError, SequenceCompilerError, SequenceCollectionError, SourceError) as failure:
        raise V2Error(f"invalid sequence source: {failure}") from failure
    if not records or len(records) > MAX_SOURCE_RECORDS:
        raise V2Error("sequence source record count is outside compiler limits")
    if len(records) != len(set(records)):
        raise V2Error("sequence source contains duplicate record ids")
    return source, records


def _provider(value: Any, label: str) -> dict[str, str]:
    item = keys(value, {"name", "version"}, label)
    text(item["name"], f"{label}.name", identifier=True)
    text(item["version"], f"{label}.version", identifier=True)
    return dict(item)


def _model(value: Any, label: str) -> dict[str, str]:
    item = keys(value, {"kind", "value"}, label)
    if item["kind"] == "content-sha256":
        sha256(item["value"], f"{label}.value")
    elif item["kind"] == "opaque":
        text(item["value"], f"{label}.value")
    else:
        raise V2Error(f"{label}.kind is unsupported")
    return dict(item)


def parse_output_contract(
    value: Any,
    label: str,
    *,
    source_artifact_sha256: str | None = None,
    records: dict[str, int] | None = None,
) -> dict[str, Any]:
    item = keys(value, {"id", "type", "unit", "axes"}, label)
    output_id = text(item["id"], f"{label}.id", identifier=True)
    tensor_type = parse_type(item["type"], f"{label}.type")
    parse_unit(item["unit"], f"{label}.unit")
    parse_axes(
        item["axes"],
        tensor_type,
        f"{label}.axes",
        source_artifact_sha256=source_artifact_sha256,
        records=records,
    )
    return {"id": output_id, "type": item["type"], "unit": item["unit"], "axes": item["axes"]}


def load_manifest(path: str | Path) -> dict[str, Any]:
    artifact = checked_artifact(path, "provider manifest")
    item = keys(
        artifact,
        {"format", "version", "provider", "model_identity", "accepts", "outputs", "artifact_sha256"},
        "provider manifest",
    )
    if item["format"] != "brainc.provider-manifest":
        raise V2Error("unsupported provider manifest format")
    exact_version(item["version"], 2, "provider manifest.version")
    _provider(item["provider"], "provider manifest.provider")
    _model(item["model_identity"], "provider manifest.model_identity")
    accepts = item["accepts"]
    if type(accepts) is not list or not accepts:
        raise V2Error("provider manifest.accepts must be a nonempty array")
    normalized_accepts = [text(value, f"provider manifest.accepts[{index}]") for index, value in enumerate(accepts)]
    sorted_unique(normalized_accepts, "provider manifest.accepts")
    outputs = item["outputs"]
    if type(outputs) is not list or not outputs or len(outputs) > MAX_TENSORS:
        raise V2Error("provider manifest.outputs must be a finite nonempty array")
    parsed = [parse_output_contract(value, f"provider manifest.outputs[{index}]") for index, value in enumerate(outputs)]
    sorted_unique([value["id"] for value in parsed], "provider manifest.outputs")
    return artifact


def make_request(
    source_path: str | Path,
    manifest_path: str | Path,
    output_ids: list[str],
) -> dict[str, Any]:
    source, records = load_source(source_path)
    _reject_standalone_descriptor(source, use="make_development_request")
    manifest = load_manifest(manifest_path)
    return _make_request(source, records, manifest, output_ids)


def _validated_source_bundle_state(
    value: str | Path | SourceBundle,
) -> tuple[dict[str, Any], dict[str, int]]:
    try:
        bundle = (
            validate_source_bundle(value)
            if isinstance(value, SourceBundle)
            else load_source_bundle(value)
        )
    except (SourceError, OSError, TypeError) as failure:
        raise V2Error(f"invalid source bundle: {failure}") from failure
    source = bundle.to_dict()
    if _source_identity(source, "source bundle descriptor") != (
        SOURCE_DESCRIPTOR_FORMAT,
        SOURCE_DESCRIPTOR_VERSION,
    ):
        raise V2Error("source bundle descriptor format/version is unsupported")
    records = {
        record["record_id"]: record["bases"]
        for record in source["source_ir"]["records"]
    }
    return source, records


def make_development_request(
    source_bundle: str | Path | SourceBundle,
    manifest_path: str | Path,
    output_ids: list[str],
) -> dict[str, Any]:
    """Bind a request only after validating the complete native source closure."""

    source, records = _validated_source_bundle_state(source_bundle)
    return _make_request(source, records, load_manifest(manifest_path), output_ids)


def _require_validated_genbank_records(
    source: dict[str, Any], records: dict[str, int]
) -> None:
    if _source_identity(source, "validated GenBank source") != (
        GENBANK_FORMAT,
        GENBANK_VERSION,
    ):
        raise V2Error("internal validated-source path requires GenBank v2")
    expected = {
        member["record_id"]: member["sequence"]["bases"]
        for member in source["sequence_collection"]["members"]
    }
    if records != expected:
        raise V2Error("validated GenBank record lengths do not match the source")


def _make_request_from_validated_genbank_source(
    source: dict[str, Any],
    records: dict[str, int],
    manifest_path: str | Path,
    output_ids: list[str],
) -> dict[str, Any]:
    """Build a request for the same-stack validated GenBank snapshot."""

    _require_validated_genbank_records(source, records)
    return _make_request(source, records, load_manifest(manifest_path), output_ids)


def _make_request(
    source: dict[str, Any],
    records: dict[str, int],
    manifest: dict[str, Any],
    output_ids: list[str],
) -> dict[str, Any]:
    source_tag = _source_acceptance_tag(source)
    if source_tag not in manifest["accepts"]:
        raise V2Error(f"provider does not declare support for {source_tag}")
    if type(output_ids) is not list or not output_ids:
        raise V2Error("at least one requested output id is required")
    normalized_ids = [text(value, f"output_ids[{index}]", identifier=True) for index, value in enumerate(output_ids)]
    sorted_unique(normalized_ids, "output_ids")
    available = {value["id"]: value for value in manifest["outputs"]}
    try:
        requested = [dict(available[output_id]) for output_id in normalized_ids]
    except KeyError as failure:
        raise V2Error(f"provider does not declare output {failure.args[0]!r}") from failure
    for index, contract in enumerate(requested):
        parse_output_contract(
            contract,
            f"requested output {index}",
            source_artifact_sha256=source["artifact_sha256"],
            records=records,
        )
    return seal(
        {
            "format": "brainc.prediction-request",
            "version": 2,
            "source": _source_binding(source),
            "provider_manifest_sha256": manifest["artifact_sha256"],
            "requested_outputs": requested,
        }
    )


def load_request(path: str | Path) -> dict[str, Any]:
    artifact = checked_artifact(path, "prediction request")
    item = keys(
        artifact,
        {"format", "version", "source", "provider_manifest_sha256", "requested_outputs", "artifact_sha256"},
        "prediction request",
    )
    if item["format"] != "brainc.prediction-request":
        raise V2Error("unsupported prediction request format")
    exact_version(item["version"], 2, "prediction request.version")
    source = keys(item["source"], {"format", "version", "artifact_sha256", "ir_sha256"}, "prediction request.source")
    if _source_identity(source, "prediction request source") not in SOURCE_FORMATS:
        raise V2Error("prediction request source format/version is unsupported")
    sha256(source["artifact_sha256"], "prediction request.source.artifact_sha256")
    sha256(source["ir_sha256"], "prediction request.source.ir_sha256")
    sha256(item["provider_manifest_sha256"], "prediction request.provider_manifest_sha256")
    outputs = item["requested_outputs"]
    if type(outputs) is not list or not outputs or len(outputs) > MAX_TENSORS:
        raise V2Error("prediction request.requested_outputs must be a finite nonempty array")
    parsed = [
        parse_output_contract(
            value,
            f"prediction request.requested_outputs[{index}]",
            source_artifact_sha256=source["artifact_sha256"],
        )
        for index, value in enumerate(outputs)
    ]
    sorted_unique([value["id"] for value in parsed], "prediction request.requested_outputs")
    return artifact


def load_response(path: str | Path) -> dict[str, Any]:
    artifact = checked_artifact(path, "prediction response")
    item = keys(
        artifact,
        {"format", "version", "request_artifact_sha256", "provider", "model_identity", "outputs", "artifact_sha256"},
        "prediction response",
    )
    if item["format"] != "brainc.prediction-response":
        raise V2Error("unsupported prediction response format")
    exact_version(item["version"], 2, "prediction response.version")
    sha256(item["request_artifact_sha256"], "prediction response.request_artifact_sha256")
    _provider(item["provider"], "prediction response.provider")
    _model(item["model_identity"], "prediction response.model_identity")
    outputs = item["outputs"]
    if type(outputs) is not list or not outputs or len(outputs) > MAX_TENSORS:
        raise V2Error("prediction response.outputs must be a finite nonempty array")
    contracts: list[dict[str, Any]] = []
    for index, output in enumerate(outputs):
        label = f"prediction response.outputs[{index}]"
        parsed = keys(output, {"id", "type", "unit", "axes", "storage"}, label)
        contract = parse_output_contract(
            {key: parsed[key] for key in ("id", "type", "unit", "axes")}, label
        )
        contracts.append(contract)
        storage = parsed["storage"]
        storage_label = f"{label}.storage"
        if type(storage) is not dict or "kind" not in storage:
            raise V2Error(f"{storage_label} must be a descriptor")
        expected_length = parse_type(contract["type"], f"{label}.type").byte_length
        if storage["kind"] == "inline-base64":
            parse_storage(storage, expected_length, storage_label)
        elif storage["kind"] == "sha256-blob":
            item = keys(storage, {"kind", "byte_length", "sha256"}, storage_label)
            if integer(item["byte_length"], f"{storage_label}.byte_length") != expected_length:
                raise V2Error(f"{storage_label}.byte_length does not match tensor type")
            sha256(item["sha256"], f"{storage_label}.sha256")
        else:
            raise V2Error(f"{storage_label}.kind is unsupported")
    sorted_unique([value["id"] for value in contracts], "prediction response.outputs")
    return artifact


def validate_binding(
    source_path: str | Path,
    manifest_path: str | Path,
    request_path: str | Path,
    response_path: str | Path,
) -> tuple[dict[str, Any], dict[str, int], dict[str, Any], dict[str, Any], dict[str, Any]]:
    source, records = load_source(source_path)
    _reject_standalone_descriptor(source, use="compile_development")
    manifest = load_manifest(manifest_path)
    request = load_request(request_path)
    response = load_response(response_path)
    return _validate_binding(source, records, manifest, request, response)


def _validate_binding_from_validated_genbank_source(
    source: dict[str, Any],
    records: dict[str, int],
    manifest_path: str | Path,
    request_path: str | Path,
    response_path: str | Path,
) -> tuple[dict[str, Any], dict[str, int], dict[str, Any], dict[str, Any], dict[str, Any]]:
    """Validate a provider chain for the same-stack validated GenBank snapshot."""

    _require_validated_genbank_records(source, records)
    return _validate_binding(
        source,
        records,
        load_manifest(manifest_path),
        load_request(request_path),
        load_response(response_path),
    )


def _validate_binding_from_validated_source_bundle(
    source_bundle: str | Path | SourceBundle,
    manifest_path: str | Path,
    request_path: str | Path,
    response_path: str | Path,
) -> tuple[dict[str, Any], dict[str, int], dict[str, Any], dict[str, Any], dict[str, Any]]:
    """Validate a provider chain after replaying a complete SourceBundle closure."""

    source, records = _validated_source_bundle_state(source_bundle)
    return _validate_binding(
        source,
        records,
        load_manifest(manifest_path),
        load_request(request_path),
        load_response(response_path),
    )


def _validate_binding(
    source: dict[str, Any],
    records: dict[str, int],
    manifest: dict[str, Any],
    request: dict[str, Any],
    response: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, int], dict[str, Any], dict[str, Any], dict[str, Any]]:
    if request["source"] != _source_binding(source):
        raise V2Error("prediction request is not bound to the supplied sequence source")
    if request["provider_manifest_sha256"] != manifest["artifact_sha256"]:
        raise V2Error("prediction request is not bound to the supplied provider manifest")
    source_tag = _source_acceptance_tag(source)
    if source_tag not in manifest["accepts"]:
        raise V2Error(f"provider does not declare support for {source_tag}")
    declared = {value["id"]: value for value in manifest["outputs"]}
    for requested in request["requested_outputs"]:
        if requested["id"] not in declared or requested != declared[requested["id"]]:
            raise V2Error(f"request output {requested['id']!r} differs from provider declaration")
    if response["request_artifact_sha256"] != request["artifact_sha256"]:
        raise V2Error("prediction response is not bound to the supplied request")
    if response["provider"] != manifest["provider"] or response["model_identity"] != manifest["model_identity"]:
        raise V2Error("prediction response provider/model identity differs from manifest")
    expected = request["requested_outputs"]
    observed = [
        {key: value[key] for key in ("id", "type", "unit", "axes")}
        for value in response["outputs"]
    ]
    if observed != expected:
        raise V2Error("prediction response tensor contracts/order differ from request")
    return source, records, manifest, request, response

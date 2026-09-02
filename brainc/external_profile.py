"""Closed manifests and evidence for external DNA source frontends.

Version 1 is data-only. Version 2 binds explicitly selected frontend and
validator executables; execution is implemented by :mod:`external_runtime`.
"""

from __future__ import annotations

from collections.abc import Mapping
import copy
import hashlib
import re
from typing import Any

from ._canonical import ContractError, SAFE_INTEGER, canonical_bytes, digest, loads
from ._io import MAX_IDENTIFIER_BYTES, MAX_JSON_BYTES, MAX_STRING_BYTES


ABI = "brainc.external-frontend-data/v1"
EXECUTABLE_ABI = "brainc.external-frontend-executable/v2"
PROFILE_FORMAT = "brainc.external-frontend-profile"
PROFILE_VERSION = 1
EXECUTABLE_PROFILE_VERSION = 2
SOURCE_FORMAT = "brainc.external-source-descriptor"
SOURCE_VERSION = 1
VALIDATION_FORMAT = "brainc.external-frontend-validation"
VALIDATION_VERSION = 1
VALIDATOR_PROTOCOL = "brainc.external-frontend-validator/v1"
FRONTEND_PROTOCOL = "brainc.external-frontend-execution/v1"
RECORD_CATALOG_SCHEMA = "brainc.sequence-record-catalog/v1"

MAX_PROFILE_BYTES = 1 * 1024 * 1024
MAX_EVIDENCE_BYTES = MAX_JSON_BYTES
MAX_ROLES = 64
MAX_RECORDS = 100_000
MAX_RECORD_ID_BYTES = MAX_IDENTIFIER_BYTES

_SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
_REFGET_RE = re.compile(r"SQ\.[A-Za-z0-9_-]{32}\Z")
_IDENTIFIER_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9._:/+@-]{0,255}\Z")
_FIELD_RE = re.compile(r"[A-Za-z_][A-Za-z0-9_]{0,255}\Z")


class ExternalProfileError(ContractError):
    """An external-profile manifest or evidence closure is invalid."""


def _fail(detail: str) -> ExternalProfileError:
    return ExternalProfileError(f"EXT001: {detail}")


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


def _text(
    value: Any,
    label: str,
    *,
    maximum_bytes: int = MAX_STRING_BYTES,
    identifier: bool = False,
    field: bool = False,
) -> str:
    if type(value) is not str or not value or value != value.strip():
        raise _fail(f"{label} must be a nonempty trimmed string")
    if any(
        ord(character) < 0x20
        or 0x7F <= ord(character) <= 0x9F
        or 0xD800 <= ord(character) <= 0xDFFF
        for character in value
    ):
        raise _fail(f"{label} contains a control character or lone surrogate")
    if len(value.encode("utf-8")) > maximum_bytes:
        raise _fail(f"{label} exceeds {maximum_bytes} UTF-8 bytes")
    if identifier and _IDENTIFIER_RE.fullmatch(value) is None:
        raise _fail(f"{label} is not a portable identifier")
    if field and _FIELD_RE.fullmatch(value) is None:
        raise _fail(f"{label} is not a portable JSON field name")
    return value


def _sha256(value: Any, label: str) -> str:
    if type(value) is not str or _SHA256_RE.fullmatch(value) is None:
        raise _fail(f"{label} must be a lowercase SHA-256 digest")
    return value


def _bounded_artifact(value: dict[str, Any], maximum: int, label: str) -> None:
    try:
        size = len(canonical_bytes(value))
    except ContractError as failure:
        raise _fail(f"{label} is not canonical I-JSON: {failure}") from failure
    if size > maximum:
        raise _fail(f"{label} exceeds {maximum} canonical bytes")


def _artifact_digest(value: dict[str, Any], label: str) -> None:
    stored = _sha256(value["artifact_sha256"], f"{label}.artifact_sha256")
    expected = digest(
        {key: member for key, member in value.items() if key != "artifact_sha256"}
    )
    if stored != expected:
        raise _fail(f"{label}.artifact_sha256 does not match canonical content")


def _role_declarations(value: Any, label: str) -> list[dict[str, Any]]:
    if type(value) is not list or not value or len(value) > MAX_ROLES:
        raise _fail(f"{label} must be a nonempty array of at most {MAX_ROLES} roles")
    roles: set[str] = set()
    result: list[dict[str, Any]] = []
    for index, member in enumerate(value):
        item_label = f"{label}[{index}]"
        if label == "profile_ir.inputs":
            item = _keys(
                member,
                {"role", "media_type", "wrapper", "maximum_byte_length"},
                item_label,
            )
            _text(item["wrapper"], f"{item_label}.wrapper", identifier=True)
        else:
            item = _keys(
                member,
                {
                    "role",
                    "media_type",
                    "format",
                    "version",
                    "schema_sha256",
                    "ir_digest_field",
                    "maximum_byte_length",
                },
                item_label,
            )
            _text(item["format"], f"{item_label}.format", identifier=True)
            _integer(item["version"], f"{item_label}.version", minimum=1)
            _sha256(item["schema_sha256"], f"{item_label}.schema_sha256")
            _text(item["ir_digest_field"], f"{item_label}.ir_digest_field", field=True)
        role = _text(
            item["role"],
            f"{item_label}.role",
            maximum_bytes=MAX_IDENTIFIER_BYTES,
            identifier=True,
        )
        if role in roles:
            raise _fail(f"{label} contains duplicate role {role!r}")
        roles.add(role)
        _text(item["media_type"], f"{item_label}.media_type", identifier=True)
        maximum = _integer(
            item["maximum_byte_length"],
            f"{item_label}.maximum_byte_length",
            minimum=1,
        )
        if label == "profile_ir.native_artifacts":
            if item["media_type"] != "application/json":
                raise _fail(f"{item_label}.media_type must be application/json")
            if maximum > MAX_EVIDENCE_BYTES:
                raise _fail(
                    f"{item_label}.maximum_byte_length exceeds the JSON evidence ceiling"
                )
        result.append(item)
    return result


def _command(value: Any, label: str, *, executable: bool) -> dict[str, Any]:
    expected = {
        "id",
        "distribution",
        "version",
        "distribution_sha256",
        "executable_sha256",
    }
    if executable:
        expected.add("runtime")
    command = _keys(value, expected, label)
    _text(command["id"], f"{label}.id", identifier=True)
    _text(command["distribution"], f"{label}.distribution", identifier=True)
    _text(command["version"], f"{label}.version")
    _sha256(command["distribution_sha256"], f"{label}.distribution_sha256")
    _sha256(command["executable_sha256"], f"{label}.executable_sha256")
    if executable and command["runtime"] not in {"native", "python"}:
        raise _fail(f"{label}.runtime must be native or python")
    return command


def _profile_ir(value: Any, manifest_version: int = PROFILE_VERSION) -> dict[str, Any]:
    executable = manifest_version == EXECUTABLE_PROFILE_VERSION
    expected = {
        "abi",
        "profile",
        "grammar",
        "inputs",
        "native_artifacts",
        "catalog",
        "limits",
        "validator",
    }
    if executable:
        expected.add("frontend")
    item = _keys(
        value,
        expected,
        "profile_ir",
    )
    if item["abi"] != (EXECUTABLE_ABI if executable else ABI):
        raise _fail("profile_ir.abi is unsupported")

    profile = _keys(item["profile"], {"id", "version"}, "profile_ir.profile")
    _text(profile["id"], "profile_ir.profile.id", identifier=True)
    _integer(profile["version"], "profile_ir.profile.version", minimum=1)

    grammar = _keys(
        item["grammar"], {"id", "version", "authority"}, "profile_ir.grammar"
    )
    _text(grammar["id"], "profile_ir.grammar.id", identifier=True)
    _text(grammar["version"], "profile_ir.grammar.version")
    authority = _keys(
        grammar["authority"], {"uri", "sha256"}, "profile_ir.grammar.authority"
    )
    _text(authority["uri"], "profile_ir.grammar.authority.uri")
    _sha256(authority["sha256"], "profile_ir.grammar.authority.sha256")

    inputs = _role_declarations(item["inputs"], "profile_ir.inputs")
    native = _role_declarations(
        item["native_artifacts"], "profile_ir.native_artifacts"
    )

    catalog = _keys(item["catalog"], {"schema"}, "profile_ir.catalog")
    if catalog["schema"] != RECORD_CATALOG_SCHEMA:
        raise _fail("profile_ir.catalog.schema is unsupported")

    limits = _keys(
        item["limits"],
        {
            "maximum_total_input_bytes",
            "maximum_total_native_bytes",
            "maximum_records",
            "maximum_total_bases",
            "maximum_record_id_bytes",
        },
        "profile_ir.limits",
    )
    total_input = _integer(
        limits["maximum_total_input_bytes"],
        "profile_ir.limits.maximum_total_input_bytes",
        minimum=1,
    )
    total_native = _integer(
        limits["maximum_total_native_bytes"],
        "profile_ir.limits.maximum_total_native_bytes",
        minimum=1,
    )
    _integer(
        limits["maximum_records"],
        "profile_ir.limits.maximum_records",
        minimum=1,
        maximum=MAX_RECORDS,
    )
    _integer(
        limits["maximum_total_bases"],
        "profile_ir.limits.maximum_total_bases",
        minimum=1,
    )
    _integer(
        limits["maximum_record_id_bytes"],
        "profile_ir.limits.maximum_record_id_bytes",
        minimum=1,
        maximum=MAX_RECORD_ID_BYTES,
    )
    if any(member["maximum_byte_length"] > total_input for member in inputs):
        raise _fail("an input role byte ceiling exceeds the cumulative input ceiling")
    if any(member["maximum_byte_length"] > total_native for member in native):
        raise _fail("a native role byte ceiling exceeds the cumulative native ceiling")
    if executable and total_native > MAX_EVIDENCE_BYTES:
        raise _fail("executable profile native bytes exceed the closure ceiling")

    if executable:
        frontend = _keys(
            item["frontend"], {"protocol", "command"}, "profile_ir.frontend"
        )
        if frontend["protocol"] != FRONTEND_PROTOCOL:
            raise _fail("profile_ir.frontend.protocol is unsupported")
        frontend_command = _command(
            frontend["command"],
            "profile_ir.frontend.command",
            executable=True,
        )

    validator = _keys(
        item["validator"], {"protocol", "command"}, "profile_ir.validator"
    )
    if validator["protocol"] != VALIDATOR_PROTOCOL:
        raise _fail("profile_ir.validator.protocol is unsupported")
    validator_command = _command(
        validator["command"],
        "profile_ir.validator.command",
        executable=executable,
    )
    if executable and (
        frontend_command["executable_sha256"]
        == validator_command["executable_sha256"]
    ):
        raise _fail("frontend and validator executable digests must be distinct")
    _bounded_artifact(item, MAX_PROFILE_BYTES, "profile_ir")
    return item


def seal_profile_manifest(
    profile_ir: Mapping[str, Any], *, version: int = PROFILE_VERSION
) -> dict[str, Any]:
    """Seal a caller-authored profile manifest after exact contract validation."""

    if not isinstance(profile_ir, Mapping):
        raise _fail("profile_ir must be a mapping")
    if type(version) is not int or version not in {
        PROFILE_VERSION,
        EXECUTABLE_PROFILE_VERSION,
    }:
        raise _fail("unsupported profile manifest version")
    normalized = copy.deepcopy(dict(profile_ir))
    _profile_ir(normalized, version)
    core = {
        "format": PROFILE_FORMAT,
        "version": version,
        "profile_ir": normalized,
        "profile_ir_sha256": digest(normalized),
    }
    payload = {**core, "artifact_sha256": digest(core)}
    return validate_profile_manifest(payload)


def validate_profile_manifest(payload: dict[str, Any]) -> dict[str, Any]:
    """Validate a closed external-profile manifest without loading code."""

    item = _keys(
        payload,
        {
            "format",
            "version",
            "profile_ir",
            "profile_ir_sha256",
            "artifact_sha256",
        },
        "profile manifest",
    )
    if item["format"] != PROFILE_FORMAT or type(item["version"]) is not int or item[
        "version"
    ] not in {PROFILE_VERSION, EXECUTABLE_PROFILE_VERSION}:
        raise _fail("unsupported profile manifest format/version")
    profile_ir = _profile_ir(item["profile_ir"], item["version"])
    if _sha256(item["profile_ir_sha256"], "profile_ir_sha256") != digest(profile_ir):
        raise _fail("profile_ir_sha256 does not match profile_ir")
    _artifact_digest(item, "profile manifest")
    _bounded_artifact(item, MAX_PROFILE_BYTES, "profile manifest")
    return copy.deepcopy(item)


def _profile_reference(manifest: dict[str, Any]) -> dict[str, Any]:
    profile = manifest["profile_ir"]["profile"]
    return {
        "id": profile["id"],
        "version": profile["version"],
        "manifest_sha256": manifest["artifact_sha256"],
    }


def _input_references(
    value: Any,
    declarations: list[dict[str, Any]],
    maximum_total: int,
    label: str,
) -> list[dict[str, Any]]:
    if type(value) is not list or len(value) != len(declarations):
        raise _fail(f"{label} must contain exactly the declared input roles")
    total = 0
    result: list[dict[str, Any]] = []
    for index, (member, declaration) in enumerate(zip(value, declarations)):
        item_label = f"{label}[{index}]"
        item = _keys(member, {"role", "sha256", "byte_length"}, item_label)
        if item["role"] != declaration["role"]:
            raise _fail(f"{item_label}.role does not match the profile role order")
        _sha256(item["sha256"], f"{item_label}.sha256")
        length = _integer(item["byte_length"], f"{item_label}.byte_length")
        if length > declaration["maximum_byte_length"]:
            raise _fail(f"{item_label}.byte_length exceeds its role ceiling")
        total += length
        if total > maximum_total:
            raise _fail(f"{label} exceeds its cumulative byte ceiling")
        result.append(item)
    return result


def _native_references(
    value: Any,
    declarations: list[dict[str, Any]],
    maximum_total: int,
    label: str,
) -> list[dict[str, Any]]:
    if type(value) is not list or len(value) != len(declarations):
        raise _fail(f"{label} must contain exactly the declared native roles")
    total = 0
    result: list[dict[str, Any]] = []
    expected_keys = {
        "role",
        "format",
        "version",
        "schema_sha256",
        "sha256",
        "byte_length",
        "ir_sha256",
    }
    for index, (member, declaration) in enumerate(zip(value, declarations)):
        item_label = f"{label}[{index}]"
        item = _keys(member, expected_keys, item_label)
        expected_identity = {
            "role": declaration["role"],
            "format": declaration["format"],
            "version": declaration["version"],
            "schema_sha256": declaration["schema_sha256"],
        }
        if {key: item[key] for key in expected_identity} != expected_identity:
            raise _fail(f"{item_label} does not match its declared schema identity")
        if type(item["version"]) is not int:
            raise _fail(f"{item_label}.version must be an integer")
        _sha256(item["sha256"], f"{item_label}.sha256")
        _sha256(item["ir_sha256"], f"{item_label}.ir_sha256")
        length = _integer(item["byte_length"], f"{item_label}.byte_length")
        if length > declaration["maximum_byte_length"]:
            raise _fail(f"{item_label}.byte_length exceeds its role ceiling")
        total += length
        if total > maximum_total:
            raise _fail(f"{label} exceeds its cumulative byte ceiling")
        result.append(item)
    return result


def _record_catalog(
    value: Any,
    input_roles: set[str],
    limits: dict[str, Any],
) -> dict[str, Any]:
    item = _keys(
        value,
        {"schema", "records", "record_count", "total_bases", "catalog_sha256"},
        "frontend_ir.record_catalog",
    )
    if item["schema"] != RECORD_CATALOG_SCHEMA:
        raise _fail("frontend_ir.record_catalog.schema is unsupported")
    records = item["records"]
    if (
        type(records) is not list
        or not records
        or len(records) > limits["maximum_records"]
    ):
        raise _fail("frontend_ir.record_catalog.records is outside its record ceiling")
    total_bases = 0
    record_ids: set[str] = set()
    for index, member in enumerate(records):
        label = f"frontend_ir.record_catalog.records[{index}]"
        record = _keys(
            member,
            {
                "ordinal",
                "input_role",
                "record_id",
                "bases",
                "sequence_sha256",
                "refget_id",
            },
            label,
        )
        if type(record["ordinal"]) is not int or record["ordinal"] != index:
            raise _fail(f"{label}.ordinal must equal its zero-based array position")
        input_role = _text(
            record["input_role"],
            f"{label}.input_role",
            maximum_bytes=MAX_IDENTIFIER_BYTES,
            identifier=True,
        )
        if input_role not in input_roles:
            raise _fail(f"{label}.input_role is not declared by the profile")
        record_id = _text(
            record["record_id"],
            f"{label}.record_id",
            maximum_bytes=limits["maximum_record_id_bytes"],
            identifier=True,
        )
        if record_id in record_ids:
            raise _fail("frontend_ir.record_catalog contains duplicate record ids")
        record_ids.add(record_id)
        total_bases += _integer(record["bases"], f"{label}.bases", minimum=1)
        if total_bases > limits["maximum_total_bases"]:
            raise _fail("frontend_ir.record_catalog exceeds its total-base ceiling")
        _sha256(record["sequence_sha256"], f"{label}.sequence_sha256")
        if type(record["refget_id"]) is not str or _REFGET_RE.fullmatch(
            record["refget_id"]
        ) is None:
            raise _fail(f"{label}.refget_id is invalid")
    if type(item["record_count"]) is not int or item["record_count"] != len(records):
        raise _fail("frontend_ir.record_catalog.record_count does not match records")
    if type(item["total_bases"]) is not int or item["total_bases"] != total_bases:
        raise _fail("frontend_ir.record_catalog.total_bases does not match records")
    catalog_core = {
        key: member for key, member in item.items() if key != "catalog_sha256"
    }
    if _sha256(item["catalog_sha256"], "catalog_sha256") != digest(catalog_core):
        raise _fail("frontend_ir.record_catalog.catalog_sha256 does not match catalog")
    return item


def _frontend_ir(manifest: dict[str, Any], value: Any) -> dict[str, Any]:
    item = _keys(
        value,
        {"profile", "inputs", "native_artifacts", "record_catalog"},
        "frontend_ir",
    )
    expected_profile = _profile_reference(manifest)
    profile = _keys(
        item["profile"], {"id", "version", "manifest_sha256"}, "frontend_ir.profile"
    )
    if profile != expected_profile or type(profile["version"]) is not int:
        raise _fail("frontend_ir.profile does not exactly bind the profile manifest")

    profile_ir = manifest["profile_ir"]
    limits = profile_ir["limits"]
    inputs = _input_references(
        item["inputs"],
        profile_ir["inputs"],
        limits["maximum_total_input_bytes"],
        "frontend_ir.inputs",
    )
    _native_references(
        item["native_artifacts"],
        profile_ir["native_artifacts"],
        limits["maximum_total_native_bytes"],
        "frontend_ir.native_artifacts",
    )
    _record_catalog(
        item["record_catalog"],
        {reference["role"] for reference in inputs},
        limits,
    )
    _bounded_artifact(item, MAX_EVIDENCE_BYTES, "frontend_ir")
    return item


def _make_catalog(records: Any) -> dict[str, Any]:
    if type(records) is not list:
        raise _fail("records must be an array")
    copied = copy.deepcopy(records)
    core = {
        "schema": RECORD_CATALOG_SCHEMA,
        "records": copied,
        "record_count": len(copied),
        "total_bases": sum(
            member.get("bases", 0)
            if type(member) is dict and type(member.get("bases", 0)) is int
            else 0
            for member in copied
        ),
    }
    try:
        catalog_sha256 = digest(core)
    except ContractError as failure:
        raise _fail(f"records are not canonical I-JSON: {failure}") from failure
    return {**core, "catalog_sha256": catalog_sha256}


def seal_source_descriptor(
    profile_manifest: dict[str, Any],
    *,
    input_references: list[dict[str, Any]],
    native_artifact_references: list[dict[str, Any]],
    records: list[dict[str, Any]],
) -> dict[str, Any]:
    """Seal an external frontend's data-only source descriptor."""

    manifest = validate_profile_manifest(profile_manifest)
    frontend_ir = {
        "profile": _profile_reference(manifest),
        "inputs": copy.deepcopy(input_references),
        "native_artifacts": copy.deepcopy(native_artifact_references),
        "record_catalog": _make_catalog(records),
    }
    _frontend_ir(manifest, frontend_ir)
    core = {
        "format": SOURCE_FORMAT,
        "version": SOURCE_VERSION,
        "frontend_ir": frontend_ir,
        "frontend_ir_sha256": digest(frontend_ir),
    }
    payload = {**core, "artifact_sha256": digest(core)}
    return validate_source_descriptor(manifest, payload)


def validate_source_descriptor(
    profile_manifest: dict[str, Any], payload: dict[str, Any]
) -> dict[str, Any]:
    """Validate an external descriptor against one caller-selected manifest."""

    manifest = validate_profile_manifest(profile_manifest)
    item = _keys(
        payload,
        {
            "format",
            "version",
            "frontend_ir",
            "frontend_ir_sha256",
            "artifact_sha256",
        },
        "external source descriptor",
    )
    if item["format"] != SOURCE_FORMAT or (
        type(item["version"]) is not int or item["version"] != SOURCE_VERSION
    ):
        raise _fail("unsupported external source descriptor format/version")
    frontend_ir = _frontend_ir(manifest, item["frontend_ir"])
    if _sha256(item["frontend_ir_sha256"], "frontend_ir_sha256") != digest(
        frontend_ir
    ):
        raise _fail("frontend_ir_sha256 does not match frontend_ir")
    _artifact_digest(item, "external source descriptor")
    _bounded_artifact(item, MAX_EVIDENCE_BYTES, "external source descriptor")
    return copy.deepcopy(item)


def _validation_ir(
    manifest: dict[str, Any], descriptor: dict[str, Any], value: Any
) -> dict[str, Any]:
    item = _keys(
        value,
        {
            "protocol",
            "validator",
            "profile",
            "frontend_output_sha256",
            "replay",
            "result",
        },
        "validation_ir",
    )
    expected_replay = {
        "inputs": descriptor["frontend_ir"]["inputs"],
        "native_artifacts": descriptor["frontend_ir"]["native_artifacts"],
        "record_catalog_sha256": descriptor["frontend_ir"]["record_catalog"][
            "catalog_sha256"
        ],
        "frontend_ir_sha256": descriptor["frontend_ir_sha256"],
    }
    if item != {
        "protocol": VALIDATOR_PROTOCOL,
        "validator": manifest["profile_ir"]["validator"]["command"],
        "profile": _profile_reference(manifest),
        "frontend_output_sha256": descriptor["artifact_sha256"],
        "replay": expected_replay,
        "result": "valid",
    }:
        raise _fail("validation_ir does not exactly bind the replayed frontend closure")
    return item


def seal_validation_report(
    profile_manifest: dict[str, Any],
    source_descriptor: dict[str, Any],
    *,
    replayed_input_references: list[dict[str, Any]],
    replayed_native_artifact_references: list[dict[str, Any]],
    replayed_records: list[dict[str, Any]],
) -> dict[str, Any]:
    """Seal caller-attested references from an external validator run."""

    manifest = validate_profile_manifest(profile_manifest)
    descriptor = validate_source_descriptor(manifest, source_descriptor)
    observed = {
        "profile": _profile_reference(manifest),
        "inputs": copy.deepcopy(replayed_input_references),
        "native_artifacts": copy.deepcopy(replayed_native_artifact_references),
        "record_catalog": _make_catalog(replayed_records),
    }
    _frontend_ir(manifest, observed)
    if observed != descriptor["frontend_ir"]:
        raise _fail("validator replay does not reproduce the frontend output")
    replay = {
        "inputs": descriptor["frontend_ir"]["inputs"],
        "native_artifacts": descriptor["frontend_ir"]["native_artifacts"],
        "record_catalog_sha256": descriptor["frontend_ir"]["record_catalog"][
            "catalog_sha256"
        ],
        "frontend_ir_sha256": descriptor["frontend_ir_sha256"],
    }
    validation_ir = {
        "protocol": VALIDATOR_PROTOCOL,
        "validator": copy.deepcopy(manifest["profile_ir"]["validator"]["command"]),
        "profile": _profile_reference(manifest),
        "frontend_output_sha256": descriptor["artifact_sha256"],
        "replay": copy.deepcopy(replay),
        "result": "valid",
    }
    _validation_ir(manifest, descriptor, validation_ir)
    core = {
        "format": VALIDATION_FORMAT,
        "version": VALIDATION_VERSION,
        "validation_ir": validation_ir,
        "validation_ir_sha256": digest(validation_ir),
    }
    payload = {**core, "artifact_sha256": digest(core)}
    return validate_validation_report(manifest, descriptor, payload)


def validate_validation_report(
    profile_manifest: dict[str, Any],
    source_descriptor: dict[str, Any],
    payload: dict[str, Any],
) -> dict[str, Any]:
    """Validate a replay report and its exact manifest/output bindings."""

    manifest = validate_profile_manifest(profile_manifest)
    descriptor = validate_source_descriptor(manifest, source_descriptor)
    item = _keys(
        payload,
        {
            "format",
            "version",
            "validation_ir",
            "validation_ir_sha256",
            "artifact_sha256",
        },
        "external validation report",
    )
    if item["format"] != VALIDATION_FORMAT or (
        type(item["version"]) is not int or item["version"] != VALIDATION_VERSION
    ):
        raise _fail("unsupported external validation report format/version")
    validation_ir = _validation_ir(manifest, descriptor, item["validation_ir"])
    if _sha256(item["validation_ir_sha256"], "validation_ir_sha256") != digest(
        validation_ir
    ):
        raise _fail("validation_ir_sha256 does not match validation_ir")
    _artifact_digest(item, "external validation report")
    _bounded_artifact(item, MAX_EVIDENCE_BYTES, "external validation report")
    return copy.deepcopy(item)


def _payload_mapping(
    value: Any,
    expected_roles: set[str],
    label: str,
) -> dict[str, bytes]:
    if not isinstance(value, Mapping):
        raise _fail(f"{label} must be a mapping")
    try:
        observed_roles = len(value)
    except (TypeError, ValueError, OverflowError) as failure:
        raise _fail(f"{label} is not a valid role mapping") from failure
    if observed_roles != len(expected_roles):
        raise _fail(f"{label} role closure does not match the descriptor")
    result: dict[str, bytes] = {}
    for role in sorted(expected_roles):
        try:
            raw = value[role]
        except (KeyError, TypeError, ValueError) as failure:
            raise _fail(f"{label} role closure does not match the descriptor") from failure
        if type(raw) is not bytes:
            raise _fail(f"{label} values must be immutable bytes")
        result[role] = raw
    return result


def _verify_payloads(
    value: Any,
    references: list[dict[str, Any]],
    label: str,
    declarations: list[dict[str, Any]] | None = None,
) -> None:
    expected_roles = [reference["role"] for reference in references]
    payloads = _payload_mapping(value, set(expected_roles), label)
    for index, reference in enumerate(references):
        raw = payloads[reference["role"]]
        if len(raw) != reference["byte_length"] or hashlib.sha256(raw).hexdigest() != reference["sha256"]:
            raise _fail(f"{label}.{reference['role']} bytes do not match their exact reference")
        if declarations is not None:
            declaration = declarations[index]
            try:
                artifact = loads(raw, f"{label}.{reference['role']}")
            except ContractError as failure:
                raise _fail(
                    f"{label}.{reference['role']} is not bounded duplicate-free JSON: {failure}"
                ) from failure
            if artifact.get("format") != declaration["format"] or (
                type(artifact.get("version")) is not int
                or artifact["version"] != declaration["version"]
            ):
                raise _fail(
                    f"{label}.{reference['role']} has the wrong format/version identity"
                )
            ir_digest_field = declaration["ir_digest_field"]
            if _sha256(
                artifact.get(ir_digest_field),
                f"{label}.{reference['role']}.{ir_digest_field}",
            ) != reference["ir_sha256"]:
                raise _fail(
                    f"{label}.{reference['role']} does not expose its referenced IR digest"
                )


def validate_external_evidence(
    profile_manifest: dict[str, Any],
    source_descriptor: dict[str, Any],
    validation_report: dict[str, Any],
    *,
    original_payloads: Mapping[str, bytes],
    native_artifact_payloads: Mapping[str, bytes],
) -> dict[str, dict[str, Any]]:
    """Validate the complete data closure, including exact original/native bytes."""

    manifest = validate_profile_manifest(profile_manifest)
    descriptor = validate_source_descriptor(manifest, source_descriptor)
    report = validate_validation_report(manifest, descriptor, validation_report)
    _verify_payloads(
        original_payloads, descriptor["frontend_ir"]["inputs"], "original_payloads"
    )
    _verify_payloads(
        native_artifact_payloads,
        descriptor["frontend_ir"]["native_artifacts"],
        "native_artifact_payloads",
        manifest["profile_ir"]["native_artifacts"],
    )
    return {
        "profile_manifest": manifest,
        "source_descriptor": descriptor,
        "validation_report": report,
    }


__all__ = [
    "ABI",
    "EXECUTABLE_ABI",
    "EXECUTABLE_PROFILE_VERSION",
    "ExternalProfileError",
    "MAX_EVIDENCE_BYTES",
    "MAX_PROFILE_BYTES",
    "MAX_RECORD_ID_BYTES",
    "MAX_RECORDS",
    "MAX_ROLES",
    "PROFILE_FORMAT",
    "PROFILE_VERSION",
    "FRONTEND_PROTOCOL",
    "RECORD_CATALOG_SCHEMA",
    "SOURCE_FORMAT",
    "SOURCE_VERSION",
    "VALIDATION_FORMAT",
    "VALIDATION_VERSION",
    "VALIDATOR_PROTOCOL",
    "seal_profile_manifest",
    "seal_source_descriptor",
    "seal_validation_report",
    "validate_external_evidence",
    "validate_profile_manifest",
    "validate_source_descriptor",
    "validate_validation_report",
]

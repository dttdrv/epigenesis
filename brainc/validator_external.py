"""Independent validator for data-only external DNA frontend closures.

Only the Python standard library is inside this trust boundary.  The module
does not import or execute the compiler, a frontend, a profile implementation,
or shared Epigenesis helpers.
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping
import copy
import functools
import hashlib
import json
import math
import os
from pathlib import Path
import re
import stat
from typing import Any


SAFE_INTEGER = 2**53 - 1
MAX_IDENTIFIER_BYTES = 256
MAX_STRING_BYTES = 1024 * 1024
MAX_JSON_DEPTH = 64
MAX_JSON_MEMBERS = 1_000_000
MAX_PROFILE_BYTES = 1024 * 1024
MAX_EVIDENCE_BYTES = 16 * 1024 * 1024
MAX_CLOSURE_BYTES = MAX_PROFILE_BYTES + 2 * MAX_EVIDENCE_BYTES + 1024 * 1024
MAX_REPORT_BYTES = 1024 * 1024
MAX_ROLES = 64
MAX_RECORDS = 100_000
MAX_RECORD_ID_BYTES = MAX_IDENTIFIER_BYTES
CHUNK_BYTES = 64 * 1024

ABI = "brainc.external-frontend-data/v1"
PROFILE_FORMAT = "brainc.external-frontend-profile"
PROFILE_VERSION = 1
SOURCE_FORMAT = "brainc.external-source-descriptor"
SOURCE_VERSION = 1
FRONTEND_VALIDATION_FORMAT = "brainc.external-frontend-validation"
FRONTEND_VALIDATION_VERSION = 1
VALIDATOR_PROTOCOL = "brainc.external-frontend-validator/v1"
RECORD_CATALOG_SCHEMA = "brainc.sequence-record-catalog/v1"

CLOSURE_FORMAT = "brainc.external-source-closure"
CLOSURE_VERSION = 1
CLOSURE_PRODUCER = {
    "name": "brainc-external-source",
    "version": "1.0.0",
    "passes": [
        "validate-profile-manifest",
        "validate-frontend-evidence",
        "bind-external-replay",
        "emit-external-source-closure",
    ],
}

REPORT_FORMAT = "brainc.external-source-validation-report"
REPORT_VERSION = 1
REPORT_STAGE = "external-frontend-admission"
REPORT_VALIDATOR = {
    "id": "brainc-independent-external-source-validator",
    "version": "1.0.0",
    "protocol": "brainc.external-source-independent-validator/v1",
}
REPORT_CHECKS = (
    "closed-profile-manifest-schema",
    "closed-source-descriptor-schema",
    "closed-frontend-validation-schema",
    "manifest-qualified-profile-identity",
    "canonical-evidence-and-closure-seals",
    "stable-single-link-original-paths",
    "exact-original-byte-identities",
    "stable-single-link-native-artifact-paths",
    "bounded-native-json-identities-and-ir-digests",
)

_SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
_REFGET_RE = re.compile(r"SQ\.[A-Za-z0-9_-]{32}\Z")
_IDENTIFIER_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9._:/+@-]{0,255}\Z")
_FIELD_RE = re.compile(r"[A-Za-z_][A-Za-z0-9_]{0,255}\Z")


class ExternalValidationError(ValueError):
    """The external closure or its independently replayed files are invalid."""


def _fail(detail: str) -> ExternalValidationError:
    return ExternalValidationError(f"EXTVAL001: {detail}")


def _bounded(function: Any) -> Any:
    @functools.wraps(function)
    def guarded(*args: Any, **kwargs: Any) -> Any:
        try:
            return function(*args, **kwargs)
        except ExternalValidationError:
            raise
        except MemoryError as failure:
            raise _fail("validation memory ceiling exceeded") from failure
        except RecursionError as failure:
            raise _fail("validation recursion ceiling exceeded") from failure

    return guarded


def _jcs_string(value: str) -> str:
    if any(0xD800 <= ord(character) <= 0xDFFF for character in value):
        raise ValueError("lone Unicode surrogate")
    return json.dumps(value, ensure_ascii=False, allow_nan=False)


def _jcs_number(value: float) -> str:
    if not math.isfinite(value):
        raise ValueError("non-finite number")
    if value == 0:
        return "0"
    sign = "-" if value < 0 else ""
    mantissa, separator, exponent_text = repr(abs(value)).lower().partition("e")
    exponent = int(exponent_text) if separator else 0
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
    normalized_exponent = point - 1
    rendered = digits[0] + ("." + digits[1:] if len(digits) > 1 else "")
    return sign + rendered + ("e+" if normalized_exponent >= 0 else "e") + str(
        normalized_exponent
    )


def _jcs_chunks(value: Any) -> Iterator[str]:
    if value is None:
        yield "null"
    elif type(value) is bool:
        yield "true" if value else "false"
    elif type(value) is str:
        yield _jcs_string(value)
    elif type(value) is int:
        if not -SAFE_INTEGER <= value <= SAFE_INTEGER:
            raise ValueError("unsafe integer")
        yield str(value)
    elif type(value) is float:
        yield _jcs_number(value)
    elif type(value) is list:
        yield "["
        for index, member in enumerate(value):
            if index:
                yield ","
            yield from _jcs_chunks(member)
        yield "]"
    elif type(value) is dict:
        if any(type(key) is not str for key in value):
            raise TypeError("non-string object key")
        yield "{"
        keys = sorted(value, key=lambda key: key.encode("utf-16be"))
        for index, key in enumerate(keys):
            if index:
                yield ","
            yield _jcs_string(key)
            yield ":"
            yield from _jcs_chunks(value[key])
        yield "}"
    else:
        raise TypeError(f"unsupported JSON value type {type(value).__name__}")


def digest(value: Any) -> str:
    """Return an independently computed RFC 8785 SHA-256 digest."""

    hasher = hashlib.sha256()
    try:
        for chunk in _jcs_chunks(value):
            hasher.update(chunk.encode("utf-8"))
    except ExternalValidationError:
        raise
    except MemoryError as failure:
        raise _fail("canonical JSON exceeds the validation memory ceiling") from failure
    except (TypeError, ValueError, UnicodeError, RecursionError) as failure:
        raise _fail(f"value is not RFC 8785 canonical JSON: {failure}") from failure
    return hasher.hexdigest()


def _same_json(left: Any, right: Any) -> bool:
    if type(left) is not type(right):
        return False
    if type(left) is dict:
        return left.keys() == right.keys() and all(
            _same_json(left[key], right[key]) for key in left
        )
    if type(left) is list:
        return len(left) == len(right) and all(
            _same_json(first, second) for first, second in zip(left, right)
        )
    return bool(left == right)


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
    try:
        encoded = value.encode("utf-8")
    except UnicodeEncodeError as failure:
        raise _fail(f"{label} is not valid Unicode") from failure
    if len(encoded) > maximum_bytes:
        raise _fail(f"{label} exceeds {maximum_bytes} UTF-8 bytes")
    if identifier and _IDENTIFIER_RE.fullmatch(value) is None:
        raise _fail(f"{label} is not a portable identifier")
    if field and _FIELD_RE.fullmatch(value) is None:
        raise _fail(f"{label} is not a portable JSON field name")
    return value


def _sha(value: Any, label: str) -> str:
    if type(value) is not str or _SHA256_RE.fullmatch(value) is None:
        raise _fail(f"{label} must be a lowercase SHA-256 digest")
    return value


def _validate_json_tree(value: Any, label: str) -> None:
    members = 0
    pending: list[tuple[Any, int]] = [(value, 0)]
    while pending:
        current, depth = pending.pop()
        if depth > MAX_JSON_DEPTH:
            raise _fail(f"{label} exceeds JSON depth ceiling {MAX_JSON_DEPTH}")
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
            try:
                encoded = current.encode("utf-8")
            except UnicodeEncodeError as failure:
                raise _fail(f"{label} contains a lone Unicode surrogate") from failure
            if len(encoded) > MAX_STRING_BYTES:
                raise _fail(
                    f"{label} exceeds JSON string byte ceiling {MAX_STRING_BYTES}"
                )
            continue
        if type(current) is list:
            members += len(current)
            pending.extend((member, depth + 1) for member in current)
        elif type(current) is dict:
            members += len(current)
            for key, member in current.items():
                if type(key) is not str:
                    raise _fail(f"{label} contains a non-string JSON key")
                pending.append((key, depth + 1))
                pending.append((member, depth + 1))
        else:
            raise _fail(f"{label} contains unsupported JSON value")
        if members > MAX_JSON_MEMBERS:
            raise _fail(f"{label} exceeds JSON member ceiling {MAX_JSON_MEMBERS}")


def _bounded_artifact(value: Any, maximum: int, label: str) -> None:
    _validate_json_tree(value, label)
    size = 0
    try:
        for chunk in _jcs_chunks(value):
            size += len(chunk.encode("utf-8"))
            if size > maximum:
                raise _fail(f"{label} exceeds {maximum} canonical bytes")
    except ExternalValidationError:
        raise
    except MemoryError as failure:
        raise _fail(f"{label} exceeds the validation memory ceiling") from failure
    except (TypeError, ValueError, UnicodeError, RecursionError) as failure:
        raise _fail(f"{label} is not canonical I-JSON: {failure}") from failure


def _artifact_seal(value: dict[str, Any], label: str) -> None:
    claimed = _sha(value["artifact_sha256"], f"{label}.artifact_sha256")
    expected = digest(
        {key: member for key, member in value.items() if key != "artifact_sha256"}
    )
    if claimed != expected:
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
            _sha(item["schema_sha256"], f"{item_label}.schema_sha256")
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


def _profile_ir(value: Any) -> dict[str, Any]:
    item = _keys(
        value,
        {
            "abi",
            "profile",
            "grammar",
            "inputs",
            "native_artifacts",
            "catalog",
            "limits",
            "validator",
        },
        "profile_ir",
    )
    if type(item["abi"]) is not str or item["abi"] != ABI:
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
    _sha(authority["sha256"], "profile_ir.grammar.authority.sha256")

    inputs = _role_declarations(item["inputs"], "profile_ir.inputs")
    native = _role_declarations(
        item["native_artifacts"], "profile_ir.native_artifacts"
    )

    catalog = _keys(item["catalog"], {"schema"}, "profile_ir.catalog")
    if type(catalog["schema"]) is not str or catalog["schema"] != RECORD_CATALOG_SCHEMA:
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
    if any(declaration["maximum_byte_length"] > total_input for declaration in inputs):
        raise _fail("an input role byte ceiling exceeds the cumulative input ceiling")
    if any(declaration["maximum_byte_length"] > total_native for declaration in native):
        raise _fail("a native role byte ceiling exceeds the cumulative native ceiling")

    validator = _keys(
        item["validator"], {"protocol", "command"}, "profile_ir.validator"
    )
    if type(validator["protocol"]) is not str or validator["protocol"] != VALIDATOR_PROTOCOL:
        raise _fail("profile_ir.validator.protocol is unsupported")
    command = _keys(
        validator["command"],
        {
            "id",
            "distribution",
            "version",
            "distribution_sha256",
            "executable_sha256",
        },
        "profile_ir.validator.command",
    )
    _text(command["id"], "profile_ir.validator.command.id", identifier=True)
    _text(
        command["distribution"],
        "profile_ir.validator.command.distribution",
        identifier=True,
    )
    _text(command["version"], "profile_ir.validator.command.version")
    _sha(
        command["distribution_sha256"],
        "profile_ir.validator.command.distribution_sha256",
    )
    _sha(
        command["executable_sha256"],
        "profile_ir.validator.command.executable_sha256",
    )
    _bounded_artifact(item, MAX_PROFILE_BYTES, "profile_ir")
    return item


def _validate_manifest(value: Any) -> dict[str, Any]:
    item = _keys(
        value,
        {
            "format",
            "version",
            "profile_ir",
            "profile_ir_sha256",
            "artifact_sha256",
        },
        "profile manifest",
    )
    if (
        type(item["format"]) is not str
        or item["format"] != PROFILE_FORMAT
        or type(item["version"]) is not int
        or item["version"] != PROFILE_VERSION
    ):
        raise _fail("unsupported profile manifest format/version")
    profile_ir = _profile_ir(item["profile_ir"])
    if _sha(item["profile_ir_sha256"], "profile_ir_sha256") != digest(profile_ir):
        raise _fail("profile_ir_sha256 does not match profile_ir")
    _artifact_seal(item, "profile manifest")
    _bounded_artifact(item, MAX_PROFILE_BYTES, "profile manifest")
    return item


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
        if type(item["role"]) is not str or item["role"] != declaration["role"]:
            raise _fail(f"{item_label}.role does not match the profile role order")
        _sha(item["sha256"], f"{item_label}.sha256")
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
    expected_keys = {
        "role",
        "format",
        "version",
        "schema_sha256",
        "sha256",
        "byte_length",
        "ir_sha256",
    }
    result: list[dict[str, Any]] = []
    for index, (member, declaration) in enumerate(zip(value, declarations)):
        item_label = f"{label}[{index}]"
        item = _keys(member, expected_keys, item_label)
        expected_identity = {
            "role": declaration["role"],
            "format": declaration["format"],
            "version": declaration["version"],
            "schema_sha256": declaration["schema_sha256"],
        }
        if not _same_json(
            {key: item[key] for key in expected_identity}, expected_identity
        ):
            raise _fail(f"{item_label} does not match its declared schema identity")
        _integer(item["version"], f"{item_label}.version", minimum=1)
        _sha(item["schema_sha256"], f"{item_label}.schema_sha256")
        _sha(item["sha256"], f"{item_label}.sha256")
        _sha(item["ir_sha256"], f"{item_label}.ir_sha256")
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
    if type(item["schema"]) is not str or item["schema"] != RECORD_CATALOG_SCHEMA:
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
        _sha(record["sequence_sha256"], f"{label}.sequence_sha256")
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
    if _sha(item["catalog_sha256"], "catalog_sha256") != digest(catalog_core):
        raise _fail("frontend_ir.record_catalog.catalog_sha256 does not match catalog")
    return item


def _frontend_ir(manifest: dict[str, Any], value: Any) -> dict[str, Any]:
    item = _keys(
        value,
        {"profile", "inputs", "native_artifacts", "record_catalog"},
        "frontend_ir",
    )
    profile = _keys(
        item["profile"], {"id", "version", "manifest_sha256"}, "frontend_ir.profile"
    )
    if not _same_json(profile, _profile_reference(manifest)):
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


def _validate_descriptor(
    manifest: dict[str, Any], value: Any
) -> dict[str, Any]:
    item = _keys(
        value,
        {
            "format",
            "version",
            "frontend_ir",
            "frontend_ir_sha256",
            "artifact_sha256",
        },
        "external source descriptor",
    )
    if (
        type(item["format"]) is not str
        or item["format"] != SOURCE_FORMAT
        or type(item["version"]) is not int
        or item["version"] != SOURCE_VERSION
    ):
        raise _fail("unsupported external source descriptor format/version")
    frontend_ir = _frontend_ir(manifest, item["frontend_ir"])
    if _sha(item["frontend_ir_sha256"], "frontend_ir_sha256") != digest(frontend_ir):
        raise _fail("frontend_ir_sha256 does not match frontend_ir")
    _artifact_seal(item, "external source descriptor")
    _bounded_artifact(item, MAX_EVIDENCE_BYTES, "external source descriptor")
    return item


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
    replay = _keys(
        item["replay"],
        {
            "inputs",
            "native_artifacts",
            "record_catalog_sha256",
            "frontend_ir_sha256",
        },
        "validation_ir.replay",
    )
    expected_replay = {
        "inputs": descriptor["frontend_ir"]["inputs"],
        "native_artifacts": descriptor["frontend_ir"]["native_artifacts"],
        "record_catalog_sha256": descriptor["frontend_ir"]["record_catalog"][
            "catalog_sha256"
        ],
        "frontend_ir_sha256": descriptor["frontend_ir_sha256"],
    }
    expected = {
        "protocol": VALIDATOR_PROTOCOL,
        "validator": manifest["profile_ir"]["validator"]["command"],
        "profile": _profile_reference(manifest),
        "frontend_output_sha256": descriptor["artifact_sha256"],
        "replay": expected_replay,
        "result": "valid",
    }
    if not _same_json(item, expected) or not _same_json(replay, expected_replay):
        raise _fail("validation_ir does not exactly bind the replayed frontend closure")
    return item


def _validate_frontend_report(
    manifest: dict[str, Any], descriptor: dict[str, Any], value: Any
) -> dict[str, Any]:
    item = _keys(
        value,
        {
            "format",
            "version",
            "validation_ir",
            "validation_ir_sha256",
            "artifact_sha256",
        },
        "external frontend validation report",
    )
    if (
        type(item["format"]) is not str
        or item["format"] != FRONTEND_VALIDATION_FORMAT
        or type(item["version"]) is not int
        or item["version"] != FRONTEND_VALIDATION_VERSION
    ):
        raise _fail("unsupported external frontend validation report format/version")
    validation_ir = _validation_ir(manifest, descriptor, item["validation_ir"])
    if _sha(item["validation_ir_sha256"], "validation_ir_sha256") != digest(
        validation_ir
    ):
        raise _fail("validation_ir_sha256 does not match validation_ir")
    _artifact_seal(item, "external frontend validation report")
    _bounded_artifact(item, MAX_EVIDENCE_BYTES, "external frontend validation report")
    return item


def _validate_closure(value: Any) -> dict[str, Any]:
    item = _keys(
        value,
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
    if (
        type(item["format"]) is not str
        or item["format"] != CLOSURE_FORMAT
        or type(item["version"]) is not int
        or item["version"] != CLOSURE_VERSION
    ):
        raise _fail("external source closure format/version is unsupported")
    if not _same_json(item["producer"], CLOSURE_PRODUCER):
        raise _fail("external source closure producer is unsupported")
    closure_ir = _keys(
        item["closure_ir"],
        {"profile_manifest", "source_descriptor", "validation_report"},
        "external source closure.closure_ir",
    )
    manifest = _validate_manifest(closure_ir["profile_manifest"])
    descriptor = _validate_descriptor(manifest, closure_ir["source_descriptor"])
    _validate_frontend_report(manifest, descriptor, closure_ir["validation_report"])
    if _sha(item["closure_ir_sha256"], "closure_ir_sha256") != digest(closure_ir):
        raise _fail("closure_ir_sha256 does not match closure_ir")
    _artifact_seal(item, "external source closure")
    _bounded_artifact(item, MAX_CLOSURE_BYTES, "external source closure")
    return item


@_bounded
def validate_external_source_closure(payload: dict[str, Any]) -> dict[str, Any]:
    """Validate every closed nested schema and seal without producer imports."""

    return copy.deepcopy(_validate_closure(payload))


def _json_preflight(raw: bytes, label: str) -> None:
    stack: list[list[Any]] = []
    members = 0
    in_string = False
    escaped = False
    string_start = 0

    def add_member() -> None:
        nonlocal members
        members += 1
        if members > MAX_JSON_MEMBERS:
            raise _fail(f"{label} exceeds JSON member ceiling {MAX_JSON_MEMBERS}")

    for index, byte in enumerate(raw):
        if in_string:
            if escaped:
                escaped = False
            elif byte == 0x5C:
                escaped = True
            elif byte == 0x22:
                in_string = False
            if index - string_start > MAX_STRING_BYTES:
                raise _fail(
                    f"{label} exceeds JSON string byte ceiling {MAX_STRING_BYTES}"
                )
            continue
        if byte in b" \t\r\n":
            continue
        if byte == 0x22:
            if stack and stack[-1] == ["array", True]:
                add_member()
                stack[-1][1] = False
            in_string = True
            string_start = index + 1
        elif byte in (0x7B, 0x5B):
            if stack and stack[-1] == ["array", True]:
                add_member()
                stack[-1][1] = False
            stack.append(["object" if byte == 0x7B else "array", byte == 0x5B])
            if len(stack) - 1 > MAX_JSON_DEPTH:
                raise _fail(f"{label} exceeds JSON depth ceiling {MAX_JSON_DEPTH}")
        elif byte == 0x3A and stack and stack[-1][0] == "object":
            add_member()
        elif byte == 0x2C and stack and stack[-1][0] == "array":
            stack[-1][1] = True
        elif byte in (0x7D, 0x5D):
            if stack:
                stack.pop()
        elif stack and stack[-1] == ["array", True]:
            add_member()
            stack[-1][1] = False


def _load_json(raw: bytes, label: str) -> dict[str, Any]:
    def pairs(items: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, member in items:
            if key in result:
                raise _fail(f"{label} contains duplicate JSON key {key!r}")
            result[key] = member
        return result

    def constant(token: str) -> None:
        raise _fail(f"{label} contains non-finite number {token}")

    def integer(token: str) -> int:
        digits = token[1:] if token.startswith("-") else token
        if len(digits) > 16:
            raise _fail(f"{label} contains an unsafe JSON integer")
        value = int(token)
        if not -SAFE_INTEGER <= value <= SAFE_INTEGER:
            raise _fail(f"{label} contains an unsafe JSON integer")
        return value

    def number(token: str) -> float:
        value = float(token)
        if not math.isfinite(value):
            raise _fail(f"{label} contains a non-finite JSON number")
        return value

    try:
        _json_preflight(raw, label)
        value = json.loads(
            raw.decode("utf-8"),
            object_pairs_hook=pairs,
            parse_constant=constant,
            parse_int=integer,
            parse_float=number,
        )
    except ExternalValidationError:
        raise
    except MemoryError as failure:
        raise _fail(f"invalid {label} JSON: memory ceiling exceeded") from failure
    except (
        UnicodeDecodeError,
        json.JSONDecodeError,
        TypeError,
        ValueError,
        RecursionError,
    ) as failure:
        raise _fail(f"invalid {label} JSON: {failure}") from failure
    if type(value) is not dict:
        raise _fail(f"{label} must contain a JSON object")
    _validate_json_tree(value, label)
    return value


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
            raw_value = value[role]
        except KeyError as failure:
            raise _fail(f"{label} role closure does not match the descriptor") from failure
        try:
            raw = os.fspath(raw_value)
        except (TypeError, ValueError, OSError) as failure:
            raise _fail(f"{label}.{role} must be a filesystem path") from failure
        if type(raw) is not str or not raw or "\0" in raw:
            raise _fail(f"{label}.{role} must be a filesystem path")
        result[role] = Path(raw)
    return result


def _read_regular(
    path: Path,
    *,
    label: str,
    maximum_bytes: int,
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
            raise _fail(f"{label} must be one regular single-link file")
        if inspected.st_size > maximum_bytes:
            raise _fail(f"{label} exceeds byte ceiling {maximum_bytes}")
        flags = os.O_RDONLY
        flags |= getattr(os, "O_CLOEXEC", 0)
        flags |= getattr(os, "O_NOFOLLOW", 0)
        flags |= getattr(os, "O_NONBLOCK", 0)
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
            or any(getattr(inspected, field) != getattr(opened, field) for field in stable)
        ):
            raise _fail(f"{label} changed while opening")
        hasher = hashlib.sha256()
        chunks: list[bytes] | None = [] if capture else None
        byte_length = 0
        while True:
            raw = os.read(
                descriptor,
                min(CHUNK_BYTES, maximum_bytes + 1 - byte_length),
            )
            if not raw:
                break
            byte_length += len(raw)
            if byte_length > maximum_bytes:
                raise _fail(f"{label} exceeds byte ceiling {maximum_bytes}")
            hasher.update(raw)
            if chunks is not None:
                chunks.append(raw)
        finished = os.fstat(descriptor)
        attached = path.lstat()
        if (
            byte_length != opened.st_size
            or stat.S_ISLNK(attached.st_mode)
            or not stat.S_ISREG(attached.st_mode)
            or attached.st_nlink != 1
            or any(getattr(opened, field) != getattr(finished, field) for field in stable)
            or any(getattr(opened, field) != getattr(attached, field) for field in stable)
        ):
            raise _fail(f"{label} changed while reading")
        return (
            {"sha256": hasher.hexdigest(), "byte_length": byte_length},
            None if chunks is None else b"".join(chunks),
        )
    except ExternalValidationError:
        raise
    except MemoryError:
        raise
    except (OSError, TypeError, ValueError) as failure:
        raise _fail(f"cannot read {label}: {failure}") from failure
    finally:
        if descriptor >= 0:
            try:
                os.close(descriptor)
            except OSError:
                pass


def _validate_native_json(
    raw: bytes,
    declaration: dict[str, Any],
    reference: dict[str, Any],
    label: str,
) -> None:
    artifact = _load_json(raw, label)
    if (
        type(artifact.get("format")) is not str
        or artifact.get("format") != declaration["format"]
        or type(artifact.get("version")) is not int
        or artifact.get("version") != declaration["version"]
    ):
        raise _fail(f"{label} has the wrong format/version")
    digest_field = declaration["ir_digest_field"]
    if _sha(artifact.get(digest_field), f"{label}.{digest_field}") != reference[
        "ir_sha256"
    ]:
        raise _fail(f"{label} exposes a different IR digest")


def _validation_report(closure: dict[str, Any]) -> dict[str, Any]:
    closure_ir = closure["closure_ir"]
    manifest = closure_ir["profile_manifest"]
    descriptor = closure_ir["source_descriptor"]
    frontend_report = closure_ir["validation_report"]
    frontend_ir = descriptor["frontend_ir"]
    catalog = frontend_ir["record_catalog"]
    inputs = frontend_ir["inputs"]
    native = frontend_ir["native_artifacts"]
    core = {
        "format": REPORT_FORMAT,
        "version": REPORT_VERSION,
        "stage": REPORT_STAGE,
        "valid": True,
        "validator": copy.deepcopy(REPORT_VALIDATOR),
        "profile": _profile_reference(manifest),
        "evidence": {
            "closure_sha256": closure["artifact_sha256"],
            "profile_manifest_sha256": manifest["artifact_sha256"],
            "source_descriptor_sha256": descriptor["artifact_sha256"],
            "frontend_validation_sha256": frontend_report["artifact_sha256"],
        },
        "result": {
            "frontend_ir_sha256": descriptor["frontend_ir_sha256"],
            "record_catalog_sha256": catalog["catalog_sha256"],
            "input_references": copy.deepcopy(inputs),
            "native_artifact_references": copy.deepcopy(native),
            "records": catalog["record_count"],
            "total_bases": catalog["total_bases"],
            "original_bytes": sum(reference["byte_length"] for reference in inputs),
            "native_bytes": sum(reference["byte_length"] for reference in native),
        },
        "checks": list(REPORT_CHECKS),
    }
    report = {**core, "report_sha256": digest(core)}
    _validate_report_schema(report)
    return report


def _report_input_references(value: Any, label: str) -> list[dict[str, Any]]:
    if type(value) is not list or not value or len(value) > MAX_ROLES:
        raise _fail(f"{label} must be a nonempty role array")
    roles: set[str] = set()
    for index, member in enumerate(value):
        item_label = f"{label}[{index}]"
        item = _keys(member, {"role", "sha256", "byte_length"}, item_label)
        role = _text(item["role"], f"{item_label}.role", identifier=True)
        if role in roles:
            raise _fail(f"{label} contains a duplicate role")
        roles.add(role)
        _sha(item["sha256"], f"{item_label}.sha256")
        _integer(item["byte_length"], f"{item_label}.byte_length")
    return value


def _report_native_references(value: Any, label: str) -> list[dict[str, Any]]:
    if type(value) is not list or not value or len(value) > MAX_ROLES:
        raise _fail(f"{label} must be a nonempty role array")
    roles: set[str] = set()
    expected = {
        "role",
        "format",
        "version",
        "schema_sha256",
        "sha256",
        "byte_length",
        "ir_sha256",
    }
    for index, member in enumerate(value):
        item_label = f"{label}[{index}]"
        item = _keys(member, expected, item_label)
        role = _text(item["role"], f"{item_label}.role", identifier=True)
        if role in roles:
            raise _fail(f"{label} contains a duplicate role")
        roles.add(role)
        _text(item["format"], f"{item_label}.format", identifier=True)
        _integer(item["version"], f"{item_label}.version", minimum=1)
        for field in ("schema_sha256", "sha256", "ir_sha256"):
            _sha(item[field], f"{item_label}.{field}")
        _integer(item["byte_length"], f"{item_label}.byte_length")
    return value


def _validate_report_schema(value: Any) -> dict[str, Any]:
    item = _keys(
        value,
        {
            "format",
            "version",
            "stage",
            "valid",
            "validator",
            "profile",
            "evidence",
            "result",
            "checks",
            "report_sha256",
        },
        "external validation report",
    )
    identity = (
        type(item["format"]) is str
        and item["format"] == REPORT_FORMAT
        and type(item["version"]) is int
        and item["version"] == REPORT_VERSION
        and type(item["stage"]) is str
        and item["stage"] == REPORT_STAGE
        and type(item["valid"]) is bool
        and item["valid"] is True
        and _same_json(item["validator"], REPORT_VALIDATOR)
    )
    if not identity:
        raise _fail("external validation report identity is invalid")
    profile = _keys(
        item["profile"], {"id", "version", "manifest_sha256"}, "report.profile"
    )
    _text(profile["id"], "report.profile.id", identifier=True)
    _integer(profile["version"], "report.profile.version", minimum=1)
    _sha(profile["manifest_sha256"], "report.profile.manifest_sha256")
    evidence = _keys(
        item["evidence"],
        {
            "closure_sha256",
            "profile_manifest_sha256",
            "source_descriptor_sha256",
            "frontend_validation_sha256",
        },
        "report.evidence",
    )
    for name, member in evidence.items():
        _sha(member, f"report.evidence.{name}")
    if evidence["profile_manifest_sha256"] != profile["manifest_sha256"]:
        raise _fail("report profile does not bind its manifest evidence")
    result = _keys(
        item["result"],
        {
            "frontend_ir_sha256",
            "record_catalog_sha256",
            "input_references",
            "native_artifact_references",
            "records",
            "total_bases",
            "original_bytes",
            "native_bytes",
        },
        "report.result",
    )
    _sha(result["frontend_ir_sha256"], "report.result.frontend_ir_sha256")
    _sha(result["record_catalog_sha256"], "report.result.record_catalog_sha256")
    inputs = _report_input_references(
        result["input_references"], "report.result.input_references"
    )
    native = _report_native_references(
        result["native_artifact_references"],
        "report.result.native_artifact_references",
    )
    _integer(result["records"], "report.result.records", minimum=1, maximum=MAX_RECORDS)
    _integer(result["total_bases"], "report.result.total_bases", minimum=1)
    original_bytes = _integer(result["original_bytes"], "report.result.original_bytes")
    native_bytes = _integer(result["native_bytes"], "report.result.native_bytes")
    if original_bytes != sum(member["byte_length"] for member in inputs):
        raise _fail("report.result.original_bytes does not match input references")
    if native_bytes != sum(member["byte_length"] for member in native):
        raise _fail("report.result.native_bytes does not match native references")
    if not _same_json(item["checks"], list(REPORT_CHECKS)):
        raise _fail("external validation report checks differ from the closed contract")
    claimed = _sha(item["report_sha256"], "report.report_sha256")
    expected = digest(
        {key: member for key, member in item.items() if key != "report_sha256"}
    )
    if claimed != expected:
        raise _fail("external validation report seal is invalid")
    _bounded_artifact(item, MAX_REPORT_BYTES, "external validation report")
    return item


@_bounded
def validate_external_report(
    report: dict[str, Any], closure: dict[str, Any] | None = None
) -> dict[str, Any]:
    """Validate a deterministic report, optionally binding it to one closure."""

    observed = _validate_report_schema(report)
    if closure is not None:
        expected = _validation_report(_validate_closure(closure))
        if not _same_json(observed, expected):
            raise _fail("external validation report does not bind the supplied closure")
    return copy.deepcopy(observed)


@_bounded
def validate_external_source(
    payload: dict[str, Any],
    *,
    original_paths: Mapping[str, str | os.PathLike[str]],
    native_artifact_paths: Mapping[str, str | os.PathLike[str]],
) -> dict[str, Any]:
    """Replay exact original/native files and return a sealed success report."""

    closure = _validate_closure(payload)
    closure_ir = closure["closure_ir"]
    manifest = closure_ir["profile_manifest"]
    descriptor = closure_ir["source_descriptor"]
    frontend_ir = descriptor["frontend_ir"]
    input_references = frontend_ir["inputs"]
    native_references = frontend_ir["native_artifacts"]
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
        observed, _ = _read_regular(
            input_paths[role],
            label=f"original source {role}",
            maximum_bytes=min(
                input_declarations[role]["maximum_byte_length"],
                reference["byte_length"],
            ),
            capture=False,
        )
        expected = {
            "sha256": reference["sha256"],
            "byte_length": reference["byte_length"],
        }
        if not _same_json(observed, expected):
            raise _fail(f"original source {role} differs from its exact reference")
    for reference in native_references:
        role = reference["role"]
        declaration = native_declarations[role]
        observed, raw = _read_regular(
            native_paths[role],
            label=f"native artifact {role}",
            maximum_bytes=min(
                declaration["maximum_byte_length"], reference["byte_length"]
            ),
            capture=True,
        )
        expected = {
            "sha256": reference["sha256"],
            "byte_length": reference["byte_length"],
        }
        if not _same_json(observed, expected):
            raise _fail(f"native artifact {role} differs from its exact reference")
        assert raw is not None
        _validate_native_json(raw, declaration, reference, f"native artifact {role}")
    report = _validation_report(closure)
    _validate_report_schema(report)
    return report


@_bounded
def validate_external_paths(
    closure_path: str | os.PathLike[str],
    *,
    original_paths: Mapping[str, str | os.PathLike[str]],
    native_artifact_paths: Mapping[str, str | os.PathLike[str]],
    maximum_closure_bytes: int = MAX_CLOSURE_BYTES,
) -> dict[str, Any]:
    """Load a safe closure path, replay its files, and return a sealed report."""

    limit = _integer(
        maximum_closure_bytes,
        "maximum_closure_bytes",
        minimum=1,
        maximum=MAX_CLOSURE_BYTES,
    )
    paths = _path_mapping({"closure": closure_path}, {"closure"}, "closure_path")
    _, raw = _read_regular(
        paths["closure"],
        label="external source closure",
        maximum_bytes=limit,
        capture=True,
    )
    assert raw is not None
    payload = _load_json(raw, "external source closure")
    return validate_external_source(
        payload,
        original_paths=original_paths,
        native_artifact_paths=native_artifact_paths,
    )


__all__ = [
    "CLOSURE_FORMAT",
    "CLOSURE_PRODUCER",
    "CLOSURE_VERSION",
    "ExternalValidationError",
    "MAX_CLOSURE_BYTES",
    "REPORT_CHECKS",
    "REPORT_FORMAT",
    "REPORT_STAGE",
    "REPORT_VALIDATOR",
    "REPORT_VERSION",
    "digest",
    "validate_external_paths",
    "validate_external_report",
    "validate_external_source",
    "validate_external_source_closure",
]

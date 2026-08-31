"""Closed target-contract v1 producer and parser."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Any

from brainc._canonical import artifact_digest, digest, keys, sha256

from ._common import V2Error, checked_artifact, exact_version, integer, seal, sorted_unique, text
from .limits import MAX_FIELDS_PER_SCHEMA, MAX_OPERATIONS, MAX_PORTS, MAX_RULES, MAX_SCHEMAS
from .tensor import TensorType, parse_type


DEV_DOMAIN = "io.github.dttdrv.epigenesis.dev"
DEV_VERSION = 1
OP_UNIT_CREATE = f"{DEV_DOMAIN}.unit.create"
OP_EDGE_CREATE = f"{DEV_DOMAIN}.edge.create"
OP_RULE_ATTACH = f"{DEV_DOMAIN}.rule.attach"
OP_PORT_BIND = f"{DEV_DOMAIN}.port.bind"
SUPPORTED_OPERATIONS = {OP_UNIT_CREATE, OP_EDGE_CREATE, OP_RULE_ATTACH, OP_PORT_BIND}

_DOMAIN_SEGMENT = r"[a-z0-9](?:[a-z0-9-]*[a-z0-9])?"
_DOMAIN_RE = re.compile(rf"(?:{_DOMAIN_SEGMENT}\.)+{_DOMAIN_SEGMENT}\Z")
_PHASES = {"INPUT", "DELIVER", "UPDATE", "EMIT", "LEARN", "OUTPUT", "SEAL"}
_SCOPES = {"subject", "source", "target", "event"}


@dataclass(frozen=True)
class FieldInfo:
    id: str
    type: TensorType
    unit: str | None
    mutability: str


@dataclass(frozen=True)
class SchemaInfo:
    id: str
    fields: tuple[FieldInfo, ...]

    def field(self, field_id: str) -> FieldInfo:
        for field in self.fields:
            if field.id == field_id:
                return field
        raise V2Error(f"schema {self.id!r} has no field {field_id!r}")


@dataclass(frozen=True)
class PortInfo:
    id: str
    direction: str
    type: TensorType
    unit: str | None


@dataclass(frozen=True)
class ParameterInfo:
    id: str
    type: TensorType
    unit: str | None


@dataclass(frozen=True)
class RuleInfo:
    id: str
    version: int
    subject: str
    schema: str
    phase: str
    parameters: tuple[ParameterInfo, ...]
    reads: tuple[tuple[str, str], ...]
    writes: tuple[tuple[str, str], ...]


@dataclass(frozen=True)
class TargetInfo:
    artifact: dict[str, Any]
    id: str
    abi_major: int
    contract_sha256: str
    artifact_sha256: str
    opsets: tuple[tuple[str, int], ...]
    unit_schemas: tuple[SchemaInfo, ...]
    edge_schemas: tuple[SchemaInfo, ...]
    ports: tuple[PortInfo, ...]
    rules: tuple[RuleInfo, ...]

    def unit_schema(self, schema_id: str) -> SchemaInfo:
        return _find(self.unit_schemas, schema_id, "unit schema")

    def edge_schema(self, schema_id: str) -> SchemaInfo:
        return _find(self.edge_schemas, schema_id, "edge schema")

    def port(self, port_id: str) -> PortInfo:
        return _find(self.ports, port_id, "port")

    def rule(self, rule_id: str) -> RuleInfo:
        return _find(self.rules, rule_id, "rule")


def _find(values: tuple[Any, ...], wanted: str, label: str) -> Any:
    for value in values:
        if value.id == wanted:
            return value
    raise V2Error(f"unknown {label}: {wanted!r}")


def _domain(value: Any, label: str) -> str:
    result = text(value, label, identifier=True)
    if _DOMAIN_RE.fullmatch(result) is None:
        raise V2Error(f"{label} must be a lowercase dotted namespaced identifier")
    return result


def _unit(value: Any, label: str) -> str | None:
    if value is None:
        return None
    return _domain(value, label)


def _numeric(value: Any, dtype: str, label: str) -> None:
    if type(value) is not dict or "kind" not in value:
        raise V2Error(f"{label} must be a numeric contract")
    kind = value["kind"]
    if kind in {"boolean", "integer", "float"}:
        keys(value, {"kind"}, label)
        if kind == "boolean" and dtype != "bool":
            raise V2Error(f"{label} boolean requires bool")
        if kind == "integer" and dtype not in {"i64", "u64"}:
            raise V2Error(f"{label} integer requires i64 or u64")
        if kind == "float" and dtype != "f64":
            raise V2Error(f"{label} float requires f64")
        return
    if kind == "binary-fixed":
        item = keys(value, {"kind", "fraction_bits", "rounding", "overflow"}, label)
        if dtype != "i64":
            raise V2Error(f"{label} binary-fixed requires i64")
        integer(item["fraction_bits"], f"{label}.fraction_bits", maximum=62)
        if item["rounding"] != "nearest-even" or item["overflow"] != "saturate":
            raise V2Error(f"{label} has unsupported fixed-point semantics")
        return
    raise V2Error(f"{label}.kind is unsupported: {kind!r}")


def _field(value: Any, label: str) -> FieldInfo:
    item = keys(value, {"id", "type", "unit", "mutability", "numeric"}, label)
    field_type = parse_type(item["type"], f"{label}.type")
    mutability = item["mutability"]
    if type(mutability) is not str or mutability not in {"constant", "state"}:
        raise V2Error(f"{label}.mutability must be constant or state")
    _numeric(item["numeric"], field_type.dtype, f"{label}.numeric")
    return FieldInfo(
        text(item["id"], f"{label}.id", identifier=True),
        field_type,
        _unit(item["unit"], f"{label}.unit"),
        mutability,
    )


def _schema(value: Any, label: str) -> SchemaInfo:
    item = keys(value, {"id", "fields"}, label)
    raw_fields = item["fields"]
    if type(raw_fields) is not list or not raw_fields or len(raw_fields) > MAX_FIELDS_PER_SCHEMA:
        raise V2Error(f"{label}.fields must be a finite nonempty array")
    fields = tuple(_field(value, f"{label}.fields[{index}]") for index, value in enumerate(raw_fields))
    sorted_unique([value.id for value in fields], f"{label}.fields")
    return SchemaInfo(text(item["id"], f"{label}.id", identifier=True), fields)


def _port(value: Any, label: str) -> PortInfo:
    item = keys(value, {"id", "direction", "type", "unit", "codec"}, label)
    if type(item["direction"]) is not str or item["direction"] not in {"input", "output"}:
        raise V2Error(f"{label}.direction must be input or output")
    codec = keys(item["codec"], {"id", "version", "semantics_sha256"}, f"{label}.codec")
    _domain(codec["id"], f"{label}.codec.id")
    integer(codec["version"], f"{label}.codec.version", minimum=1)
    sha256(codec["semantics_sha256"], f"{label}.codec.semantics_sha256")
    return PortInfo(
        text(item["id"], f"{label}.id", identifier=True),
        item["direction"],
        parse_type(item["type"], f"{label}.type"),
        _unit(item["unit"], f"{label}.unit"),
    )


def _access(value: Any, label: str) -> tuple[str, str]:
    item = keys(value, {"scope", "field"}, label)
    if type(item["scope"]) is not str or item["scope"] not in _SCOPES:
        raise V2Error(f"{label}.scope is unsupported")
    return item["scope"], text(item["field"], f"{label}.field", identifier=True)


def _parameter(value: Any, label: str) -> ParameterInfo:
    item = keys(value, {"id", "type", "unit"}, label)
    return ParameterInfo(
        text(item["id"], f"{label}.id", identifier=True),
        parse_type(item["type"], f"{label}.type"),
        _unit(item["unit"], f"{label}.unit"),
    )


def _rule(value: Any, label: str) -> RuleInfo:
    item = keys(
        value,
        {"id", "version", "semantics_sha256", "subject", "schema", "phase", "triggers", "parameters", "reads", "writes"},
        label,
    )
    subject = item["subject"]
    if type(subject) is not str or subject not in {"unit", "edge"}:
        raise V2Error(f"{label}.subject must be unit or edge")
    phase = item["phase"]
    if type(phase) is not str or phase not in _PHASES:
        raise V2Error(f"{label}.phase is unsupported")
    triggers_raw = item["triggers"]
    if type(triggers_raw) is not list or not triggers_raw:
        raise V2Error(f"{label}.triggers must be a nonempty array")
    triggers = [text(value, f"{label}.triggers[{index}]", identifier=True) for index, value in enumerate(triggers_raw)]
    sorted_unique(triggers, f"{label}.triggers")
    for collection in ("parameters", "reads", "writes"):
        if type(item[collection]) is not list:
            raise V2Error(f"{label}.{collection} must be an array")
    parameters = tuple(_parameter(value, f"{label}.parameters[{index}]") for index, value in enumerate(item["parameters"]))
    sorted_unique([value.id for value in parameters], f"{label}.parameters")
    reads = tuple(_access(value, f"{label}.reads[{index}]") for index, value in enumerate(item["reads"]))
    writes = tuple(_access(value, f"{label}.writes[{index}]") for index, value in enumerate(item["writes"]))
    if list(reads) != sorted(reads) or len(reads) != len(set(reads)):
        raise V2Error(f"{label}.reads must be sorted and unique")
    if list(writes) != sorted(writes) or len(writes) != len(set(writes)):
        raise V2Error(f"{label}.writes must be sorted and unique")
    if any(scope != "subject" for scope, _ in writes):
        raise V2Error(f"{label}.writes must be subject-scoped")
    if subject == "unit" and any(scope in {"source", "target"} for scope, _ in reads):
        raise V2Error(f"{label} unit rule cannot read source/target scope")
    _domain(item["id"], f"{label}.id")
    sha256(item["semantics_sha256"], f"{label}.semantics_sha256")
    return RuleInfo(
        item["id"],
        integer(item["version"], f"{label}.version", minimum=1),
        subject,
        text(item["schema"], f"{label}.schema", identifier=True),
        phase,
        parameters,
        reads,
        writes,
    )


def parse_target_artifact(artifact: dict[str, Any]) -> TargetInfo:
    item = keys(
        artifact,
        {"format", "version", "contract", "contract_sha256", "artifact_sha256"},
        "target contract",
    )
    if item["format"] != "brainc.target-contract":
        raise V2Error("unsupported target contract format")
    exact_version(item["version"], 1, "target contract.version")
    body = keys(
        item["contract"],
        {"id", "abi_major", "opsets", "unit_schemas", "edge_schemas", "ports", "rules"},
        "target contract body",
    )
    contract_sha = sha256(item["contract_sha256"], "target contract.contract_sha256")
    if contract_sha != digest(body):
        raise V2Error("target contract.contract_sha256 does not match contract")
    try:
        artifact_sha = artifact_digest(item)
    except Exception as failure:
        raise V2Error(str(failure)) from failure
    if integer(body["abi_major"], "target contract.abi_major", minimum=1) != 1:
        raise V2Error("compiler supports target ABI major 1 only")
    opsets_raw = body["opsets"]
    if type(opsets_raw) is not list or not opsets_raw or len(opsets_raw) > MAX_OPERATIONS:
        raise V2Error("target contract.opsets must be a nonempty array")
    opsets: list[tuple[str, int]] = []
    for index, value in enumerate(opsets_raw):
        operation_set = keys(value, {"domain", "version"}, f"target contract.opsets[{index}]")
        opsets.append((
            _domain(operation_set["domain"], f"target contract.opsets[{index}].domain"),
            integer(operation_set["version"], f"target contract.opsets[{index}].version", minimum=1),
        ))
    if opsets != sorted(opsets) or len({domain for domain, _ in opsets}) != len(opsets):
        raise V2Error("target contract.opsets must be sorted with unique domains")
    parsed_collections: dict[str, tuple[Any, ...]] = {}
    for name, parser, limit, allow_empty in (
        ("unit_schemas", _schema, MAX_SCHEMAS, False),
        ("edge_schemas", _schema, MAX_SCHEMAS, True),
        ("ports", _port, MAX_PORTS, True),
        ("rules", _rule, MAX_RULES, True),
    ):
        values = body[name]
        if type(values) is not list or (not allow_empty and not values) or len(values) > limit:
            raise V2Error(f"target contract.{name} is outside compiler limits")
        parsed = tuple(parser(value, f"target contract.{name}[{index}]") for index, value in enumerate(values))
        sorted_unique([value.id for value in parsed], f"target contract.{name}")
        parsed_collections[name] = parsed
    result = TargetInfo(
        artifact,
        _domain(body["id"], "target contract.id"),
        body["abi_major"],
        contract_sha,
        artifact_sha,
        tuple(opsets),
        parsed_collections["unit_schemas"],
        parsed_collections["edge_schemas"],
        parsed_collections["ports"],
        parsed_collections["rules"],
    )
    all_unit_fields = {field.id for schema in result.unit_schemas for field in schema.fields}
    for rule in result.rules:
        schema = result.unit_schema(rule.schema) if rule.subject == "unit" else result.edge_schema(rule.schema)
        for scope, field_id in rule.reads + rule.writes:
            if scope == "event":
                continue
            if scope == "subject":
                field = schema.field(field_id)
                if (scope, field_id) in rule.writes and field.mutability != "state":
                    raise V2Error(f"rule {rule.id!r} writes constant field {field_id!r}")
            elif field_id not in all_unit_fields:
                raise V2Error(f"rule {rule.id!r} {scope} field does not resolve: {field_id}")
    return result


def target_artifact(contract: dict[str, Any]) -> dict[str, Any]:
    if type(contract) is not dict:
        raise V2Error("target contract body must be an object")
    contract_copy = deepcopy(contract)
    core = {
        "format": "brainc.target-contract",
        "version": 1,
        "contract": contract_copy,
        "contract_sha256": digest(contract_copy),
    }
    artifact = seal(core)
    parse_target_artifact(artifact)
    return artifact


def load_target(path: str | Path) -> TargetInfo:
    return parse_target_artifact(checked_artifact(path, "target contract"))


__all__ = [
    "DEV_DOMAIN",
    "DEV_VERSION",
    "OP_EDGE_CREATE",
    "OP_PORT_BIND",
    "OP_RULE_ATTACH",
    "OP_UNIT_CREATE",
    "SUPPORTED_OPERATIONS",
    "TargetInfo",
    "load_target",
    "parse_target_artifact",
    "target_artifact",
]

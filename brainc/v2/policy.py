"""Target-bound tensor lowering policy v2."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any

from ._common import V2Error, checked_artifact, exact_version, integer, keys, seal, sha256, sorted_unique, text
from .limits import MAX_OPERATIONS, MAX_TENSORS
from .target import (
    OP_EDGE_CREATE,
    OP_PORT_BIND,
    OP_RULE_ATTACH,
    OP_UNIT_CREATE,
    SUPPORTED_OPERATIONS,
)


def _target_binding(value: Any, label: str) -> dict[str, Any]:
    item = keys(value, {"id", "abi_major", "contract_sha256"}, label)
    text(item["id"], f"{label}.id", identifier=True)
    integer(item["abi_major"], f"{label}.abi_major", minimum=1)
    sha256(item["contract_sha256"], f"{label}.contract_sha256")
    return dict(item)


def _value(value: Any, label: str) -> dict[str, str]:
    item = keys(value, {"id", "from_output"}, label)
    return {
        "id": text(item["id"], f"{label}.id", identifier=True),
        "from_output": text(item["from_output"], f"{label}.from_output", identifier=True),
    }


def _operation(value: Any, label: str) -> dict[str, Any]:
    if type(value) is not dict or "op" not in value:
        raise V2Error(f"{label} must be an operation object")
    op_name = value["op"]
    if type(op_name) is not str or op_name not in SUPPORTED_OPERATIONS:
        raise V2Error(f"{label}.op is unsupported")
    common = {"id", "op", "version"}
    if op_name == OP_UNIT_CREATE:
        expected = common | {"schema", "count", "initializers"}
    elif op_name == OP_EDGE_CREATE:
        expected = common | {"schema", "sources", "targets", "pairs", "initializers"}
    elif op_name == OP_RULE_ATTACH:
        expected = common | {"rule", "rule_version", "subject", "parameters"}
    else:
        expected = common | {"port", "subject", "field", "indices"}
    item = keys(value, expected, label)
    text(item["id"], f"{label}.id", identifier=True)
    exact_version(item["version"], 1, f"{label}.version")
    for name in expected - common - {"initializers", "parameters", "rule_version"}:
        text(item[name], f"{label}.{name}", identifier=True)
    if "rule_version" in item:
        integer(item["rule_version"], f"{label}.rule_version", minimum=1)
    if "initializers" in item:
        bindings = item["initializers"]
        if type(bindings) is not list:
            raise V2Error(f"{label}.initializers must be an array")
        fields: list[str] = []
        for index, binding in enumerate(bindings):
            parsed = keys(binding, {"field", "tensor"}, f"{label}.initializers[{index}]")
            fields.append(text(parsed["field"], f"{label}.initializers[{index}].field", identifier=True))
            text(parsed["tensor"], f"{label}.initializers[{index}].tensor", identifier=True)
        sorted_unique(fields, f"{label}.initializers")
    if "parameters" in item:
        bindings = item["parameters"]
        if type(bindings) is not list:
            raise V2Error(f"{label}.parameters must be an array")
        ids: list[str] = []
        for index, binding in enumerate(bindings):
            parsed = keys(binding, {"id", "tensor"}, f"{label}.parameters[{index}]")
            ids.append(text(parsed["id"], f"{label}.parameters[{index}].id", identifier=True))
            text(parsed["tensor"], f"{label}.parameters[{index}].tensor", identifier=True)
        sorted_unique(ids, f"{label}.parameters")
    return dict(item)


def parse_policy(artifact: dict[str, Any]) -> dict[str, Any]:
    item = keys(
        artifact,
        {"format", "version", "id", "target", "values", "operations", "artifact_sha256"},
        "lowering policy",
    )
    if item["format"] != "brainc.lowering-policy":
        raise V2Error("unsupported lowering policy format")
    exact_version(item["version"], 2, "lowering policy.version")
    text(item["id"], "lowering policy.id", identifier=True)
    _target_binding(item["target"], "lowering policy.target")
    values = item["values"]
    if type(values) is not list or not values or len(values) > MAX_TENSORS:
        raise V2Error("lowering policy.values must be a finite nonempty array")
    parsed_values = [_value(value, f"lowering policy.values[{index}]") for index, value in enumerate(values)]
    sorted_unique([value["id"] for value in parsed_values], "lowering policy.values")
    outputs = [value["from_output"] for value in parsed_values]
    if len(outputs) != len(set(outputs)):
        raise V2Error("lowering policy.values must map provider outputs one-to-one")
    operations = item["operations"]
    if type(operations) is not list or len(operations) > MAX_OPERATIONS:
        raise V2Error("lowering policy.operations must be a finite array")
    parsed_operations = [_operation(value, f"lowering policy.operations[{index}]") for index, value in enumerate(operations)]
    ids = [value["id"] for value in parsed_operations]
    if len(ids) != len(set(ids)):
        raise V2Error("lowering policy operation ids must be unique")
    return artifact


def policy_artifact(
    policy_id: str,
    target: dict[str, Any],
    values: list[dict[str, Any]],
    operations: list[dict[str, Any]],
) -> dict[str, Any]:
    core = {
        "format": "brainc.lowering-policy",
        "version": 2,
        "id": policy_id,
        "target": deepcopy(target),
        "values": deepcopy(values),
        "operations": deepcopy(operations),
    }
    artifact = seal(core)
    parse_policy(artifact)
    return artifact


def load_policy(path: str | Path) -> dict[str, Any]:
    return parse_policy(checked_artifact(path, "lowering policy"))


__all__ = [
    "OP_EDGE_CREATE",
    "OP_PORT_BIND",
    "OP_RULE_ATTACH",
    "OP_UNIT_CREATE",
    "load_policy",
    "parse_policy",
    "policy_artifact",
]

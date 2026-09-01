"""Typed lowering from tensor-provider artifacts to a developmental module."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from brainc._canonical import SAFE_INTEGER, digest, keys
from brainc.source import SourceBundle

from ._common import V2Error, pretty_bytes, seal, sorted_unique, text
from .limits import (
    MAX_ATTACHMENTS,
    MAX_EDGES,
    MAX_OPERATIONS,
    MAX_PORTS,
    MAX_TENSORS,
    MAX_TOTAL_TENSOR_BYTES,
    MAX_UNITS,
)
from .policy import load_policy
from .provider import (
    _source_binding,
    _validate_binding_from_validated_genbank_source,
    _validate_binding_from_validated_source_bundle,
    validate_binding,
)
from .target import (
    DEV_DOMAIN,
    DEV_VERSION,
    OP_EDGE_CREATE,
    OP_PORT_BIND,
    OP_RULE_ATTACH,
    OP_UNIT_CREATE,
    FieldInfo,
    RuleInfo,
    SchemaInfo,
    TargetInfo,
    load_target,
)
from .tensor import (
    ParsedTensor,
    TensorType,
    parse_tensor_output,
    unpack_u64_pairs,
    unpack_u64_scalar,
    unpack_u64_vector,
)


class CompilerError(V2Error):
    """A v2 compilation request cannot be lowered to the declared target."""


COMPILER = {
    "name": "brainc",
    "version": "0.5.0",
    "passes": [
        "validate-sequence-source",
        "bind-provider-tensor-contract",
        "validate-target-contract",
        "validate-lowering-policy",
        "type-check-tensors",
        "validate-development-operations",
        "compute-exact-budgets",
        "emit-development-module",
    ],
}


@dataclass(frozen=True)
class UnitSet:
    id: str
    schema: str
    count: int
    first_id: int


@dataclass(frozen=True)
class EdgeSet:
    id: str
    schema: str
    sources: str
    targets: str
    pairs: tuple[tuple[int, int], ...]
    first_id: int


OperationSet = UnitSet | EdgeSet


def _require_tensor(tensors: dict[str, ParsedTensor], tensor_id: str, label: str) -> ParsedTensor:
    try:
        return tensors[tensor_id]
    except KeyError as failure:
        raise CompilerError(f"{label} references unknown tensor {tensor_id!r}") from failure


def _initializers(value: Any, label: str) -> tuple[tuple[str, str], ...]:
    if type(value) is not list:
        raise CompilerError(f"{label} must be an array")
    result: list[tuple[str, str]] = []
    for index, binding in enumerate(value):
        item = keys(binding, {"field", "tensor"}, f"{label}[{index}]")
        result.append(
            (
                text(item["field"], f"{label}[{index}].field", identifier=True),
                text(item["tensor"], f"{label}[{index}].tensor", identifier=True),
            )
        )
    sorted_unique([field for field, _ in result], label)
    return tuple(result)


def _parameters(value: Any, label: str) -> tuple[tuple[str, str], ...]:
    if type(value) is not list:
        raise CompilerError(f"{label} must be an array")
    result: list[tuple[str, str]] = []
    for index, binding in enumerate(value):
        item = keys(binding, {"id", "tensor"}, f"{label}[{index}]")
        result.append(
            (
                text(item["id"], f"{label}[{index}].id", identifier=True),
                text(item["tensor"], f"{label}[{index}].tensor", identifier=True),
            )
        )
    sorted_unique([identifier for identifier, _ in result], label)
    return tuple(result)


def _validate_initializers(
    initializers: tuple[tuple[str, str], ...],
    schema: SchemaInfo,
    count: int,
    tensors: dict[str, ParsedTensor],
    label: str,
) -> None:
    if tuple(field for field, _ in initializers) != tuple(field.id for field in schema.fields):
        raise CompilerError(f"{label} must initialize every schema field exactly once")
    for index, (field_id, tensor_id) in enumerate(initializers):
        field = schema.field(field_id)
        tensor = _require_tensor(tensors, tensor_id, f"{label}[{index}]")
        if tensor.type.dtype != field.type.dtype or tensor.unit != field.unit:
            raise CompilerError(f"{label}[{index}] dtype/unit does not match field")
        suffix = field.type.shape
        if tensor.type.shape not in (suffix, (count,) + suffix):
            raise CompilerError(f"{label}[{index}] uses an unsupported initializer shape")


def _validate_rule_accesses(
    rule: RuleInfo,
    subject: OperationSet,
    operation_sets: dict[str, OperationSet],
    target: TargetInfo,
) -> None:
    if isinstance(subject, UnitSet):
        return
    source_set = operation_sets[subject.sources]
    target_set = operation_sets[subject.targets]
    assert isinstance(source_set, UnitSet) and isinstance(target_set, UnitSet)
    source_schema = target.unit_schema(source_set.schema)
    target_schema = target.unit_schema(target_set.schema)
    for scope, field_id in rule.reads:
        if scope == "source":
            source_schema.field(field_id)
        elif scope == "target":
            target_schema.field(field_id)


def _compile_operations(
    operations: list[dict[str, Any]],
    tensors: dict[str, ParsedTensor],
    target: TargetInfo,
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    if len(operations) > MAX_OPERATIONS:
        raise CompilerError("operation count exceeds compiler ceiling")
    emitted: list[dict[str, Any]] = []
    operation_sets: dict[str, OperationSet] = {}
    all_ids: set[str] = set()
    next_unit_id = 0
    next_edge_id = 0
    attachments = 0
    seen_edges: set[tuple[str, int, int]] = set()
    attached: set[tuple[str, str]] = set()
    writes: set[tuple[str, str, str]] = set()
    bound_ports: set[str] = set()
    input_writes: set[tuple[int, str]] = set()

    for index, raw in enumerate(operations):
        label = f"lowering policy.operations[{index}]"
        operation_id = raw["id"]
        if operation_id in all_ids:
            raise CompilerError(f"duplicate operation id: {operation_id!r}")
        all_ids.add(operation_id)
        op_name = raw["op"]
        if op_name == OP_UNIT_CREATE:
            schema = target.unit_schema(raw["schema"])
            count = unpack_u64_scalar(
                _require_tensor(tensors, raw["count"], label), f"{label}.count"
            )
            initializers = _initializers(raw["initializers"], f"{label}.initializers")
            _validate_initializers(initializers, schema, count, tensors, f"{label}.initializers")
            if next_unit_id + count > MAX_UNITS or next_unit_id + count > SAFE_INTEGER:
                raise CompilerError(f"{label} exceeds unit ceiling")
            operation_sets[operation_id] = UnitSet(operation_id, schema.id, count, next_unit_id)
            next_unit_id += count
        elif op_name == OP_EDGE_CREATE:
            schema = target.edge_schema(raw["schema"])
            source_set = operation_sets.get(raw["sources"])
            target_set = operation_sets.get(raw["targets"])
            if not isinstance(source_set, UnitSet) or not isinstance(target_set, UnitSet):
                raise CompilerError(f"{label} source/target must reference earlier unit.create operations")
            pairs = unpack_u64_pairs(
                _require_tensor(tensors, raw["pairs"], label), f"{label}.pairs"
            )
            for row, (source_local, target_local) in enumerate(pairs):
                if source_local >= source_set.count or target_local >= target_set.count:
                    raise CompilerError(f"{label}.pairs[{row}] endpoint is out of range")
                edge_key = (
                    schema.id,
                    source_set.first_id + source_local,
                    target_set.first_id + target_local,
                )
                if edge_key in seen_edges:
                    raise CompilerError(f"{label}.pairs[{row}] duplicates a typed edge")
                seen_edges.add(edge_key)
            initializers = _initializers(raw["initializers"], f"{label}.initializers")
            _validate_initializers(initializers, schema, len(pairs), tensors, f"{label}.initializers")
            if next_edge_id + len(pairs) > MAX_EDGES or next_edge_id + len(pairs) > SAFE_INTEGER:
                raise CompilerError(f"{label} exceeds edge ceiling")
            operation_sets[operation_id] = EdgeSet(
                operation_id,
                schema.id,
                source_set.id,
                target_set.id,
                pairs,
                next_edge_id,
            )
            next_edge_id += len(pairs)
        elif op_name == OP_RULE_ATTACH:
            rule = target.rule(raw["rule"])
            if raw["rule_version"] != rule.version:
                raise CompilerError(f"{label}.rule_version does not match target rule")
            subject = operation_sets.get(raw["subject"])
            expected_type = UnitSet if rule.subject == "unit" else EdgeSet
            if not isinstance(subject, expected_type) or subject.schema != rule.schema:
                raise CompilerError(f"{label}.subject does not match rule subject/schema")
            parameters = _parameters(raw["parameters"], f"{label}.parameters")
            if tuple(identifier for identifier, _ in parameters) != tuple(value.id for value in rule.parameters):
                raise CompilerError(f"{label}.parameters must bind every target parameter exactly once")
            for (identifier, tensor_id), contract in zip(parameters, rule.parameters):
                tensor = _require_tensor(tensors, tensor_id, label)
                if tensor.type != contract.type or tensor.unit != contract.unit:
                    raise CompilerError(f"{label} parameter {identifier!r} type/unit mismatch")
            attach_key = (rule.id, subject.id)
            if attach_key in attached:
                raise CompilerError(f"{label} duplicates a rule attachment")
            attached.add(attach_key)
            for _, field_id in rule.writes:
                write_key = (subject.id, rule.phase, field_id)
                if write_key in writes:
                    raise CompilerError(f"{label} creates a same-phase rule write conflict")
                writes.add(write_key)
            _validate_rule_accesses(rule, subject, operation_sets, target)
            attachments += 1
            if attachments > MAX_ATTACHMENTS:
                raise CompilerError(f"{label} exceeds attachment ceiling")
        elif op_name == OP_PORT_BIND:
            port = target.port(raw["port"])
            if port.id in bound_ports:
                raise CompilerError(f"{label} binds a target port more than once")
            subject = operation_sets.get(raw["subject"])
            if not isinstance(subject, UnitSet):
                raise CompilerError(f"{label}.subject must reference an earlier unit.create")
            schema = target.unit_schema(subject.schema)
            field = schema.field(raw["field"])
            if port.direction == "input" and field.mutability != "state":
                raise CompilerError(f"{label} input port cannot bind a constant field")
            indices = unpack_u64_vector(
                _require_tensor(tensors, raw["indices"], label), f"{label}.indices"
            )
            if len(indices) != len(set(indices)):
                raise CompilerError(f"{label}.indices must be unique")
            if any(value >= subject.count for value in indices):
                raise CompilerError(f"{label}.indices contains an out-of-range index")
            expected_type = TensorType(field.type.dtype, (len(indices),) + field.type.shape)
            if port.type != expected_type or port.unit != field.unit:
                raise CompilerError(f"{label} port type/unit does not match bound field")
            if port.direction == "input":
                for local_id in indices:
                    write_key = (subject.first_id + local_id, field.id)
                    if write_key in input_writes:
                        raise CompilerError(f"{label} overlaps another input binding")
                    input_writes.add(write_key)
            bound_ports.add(port.id)
            if len(bound_ports) > MAX_PORTS:
                raise CompilerError(f"{label} exceeds port ceiling")
        else:  # policy parsing makes this unreachable, but retain a closed compiler boundary.
            raise CompilerError(f"{label}.op is unsupported")
        emitted.append(deepcopy(raw))

    if bound_ports != {port.id for port in target.ports}:
        raise CompilerError("development module must bind every target port exactly once")
    return emitted, {
        "operations": len(emitted),
        "units": next_unit_id,
        "edges": next_edge_id,
        "attachments": attachments,
    }


def compile_module(
    source_path: str | Path,
    manifest_path: str | Path,
    request_path: str | Path,
    response_path: str | Path,
    policy_path: str | Path,
    target_path: str | Path,
    *,
    blob_root: str | Path | None = None,
) -> dict[str, Any]:
    source, records, manifest, request, response = validate_binding(
        source_path, manifest_path, request_path, response_path
    )
    return _compile_module(
        source,
        records,
        manifest,
        request,
        response,
        policy_path,
        target_path,
        blob_root=blob_root,
    )


def _compile_module_from_validated_genbank_source(
    source: dict[str, Any],
    records: dict[str, int],
    manifest_path: str | Path,
    request_path: str | Path,
    response_path: str | Path,
    policy_path: str | Path,
    target_path: str | Path,
    *,
    blob_root: str | Path | None = None,
) -> dict[str, Any]:
    """Lower a chain for the same-stack validated GenBank snapshot."""

    source, records, manifest, request, response = (
        _validate_binding_from_validated_genbank_source(
            source,
            records,
            manifest_path,
            request_path,
            response_path,
        )
    )
    return _compile_module(
        source,
        records,
        manifest,
        request,
        response,
        policy_path,
        target_path,
        blob_root=blob_root,
    )


def _compile_module_from_validated_source_bundle(
    source_bundle: str | Path | SourceBundle,
    manifest_path: str | Path,
    request_path: str | Path,
    response_path: str | Path,
    policy_path: str | Path,
    target_path: str | Path,
    *,
    blob_root: str | Path | None = None,
) -> dict[str, Any]:
    """Lower a chain only after replaying its complete native source closure."""

    source, records, manifest, request, response = (
        _validate_binding_from_validated_source_bundle(
            source_bundle,
            manifest_path,
            request_path,
            response_path,
        )
    )
    return _compile_module(
        source,
        records,
        manifest,
        request,
        response,
        policy_path,
        target_path,
        blob_root=blob_root,
    )


def _compile_module(
    source: dict[str, Any],
    records: dict[str, int],
    manifest: dict[str, Any],
    request: dict[str, Any],
    response: dict[str, Any],
    policy_path: str | Path,
    target_path: str | Path,
    *,
    blob_root: str | Path | None = None,
) -> dict[str, Any]:
    target = load_target(target_path)
    policy = load_policy(policy_path)
    expected_target = {
        "id": target.id,
        "abi_major": target.abi_major,
        "contract_sha256": target.contract_sha256,
    }
    if policy["target"] != expected_target:
        raise CompilerError("lowering policy is not bound to the supplied target contract")
    values = policy["values"]
    response_by_id = {
        value["id"]: (index, value)
        for index, value in enumerate(response["outputs"])
    }
    if {value["from_output"] for value in values} != set(response_by_id):
        raise CompilerError("lowering policy must consume every response output exactly once")
    parsed_tensors: list[ParsedTensor] = []
    emitted_tensors: list[dict[str, Any]] = []
    total_tensor_bytes = 0
    for index, value in enumerate(values):
        response_index, response_output = response_by_id[value["from_output"]]
        parsed = parse_tensor_output(
            response_output,
            f"prediction response.outputs[{response_index}]",
            blob_root=blob_root,
            source_artifact_sha256=source["artifact_sha256"],
            records=records,
        )
        emitted = {
            "id": value["id"],
            "type": deepcopy(response_output["type"]),
            "unit": response_output["unit"],
            "axes": deepcopy(response_output["axes"]),
            "storage": deepcopy(response_output["storage"]),
            "lineage": {
                "response_output": value["from_output"],
                "lowering_value": value["id"],
            },
        }
        lowered = ParsedTensor(
            emitted,
            value["id"],
            parsed.type,
            parsed.unit,
            parsed.axes,
            parsed.raw,
        )
        parsed_tensors.append(lowered)
        emitted_tensors.append(emitted)
        total_tensor_bytes += parsed.type.byte_length
        if total_tensor_bytes > SAFE_INTEGER or total_tensor_bytes > MAX_TOTAL_TENSOR_BYTES:
            raise CompilerError("total tensor bytes exceed compiler ceiling")
    if len(parsed_tensors) > MAX_TENSORS:
        raise CompilerError("tensor count exceeds compiler ceiling")
    sorted_unique([tensor.id for tensor in parsed_tensors], "development module.tensors")
    tensor_by_id = {tensor.id: tensor for tensor in parsed_tensors}
    operations, budget_counts = _compile_operations(policy["operations"], tensor_by_id, target)
    requirements = [{"domain": DEV_DOMAIN, "version": DEV_VERSION}] if operations else []
    if operations and (DEV_DOMAIN, DEV_VERSION) not in target.opsets:
        raise CompilerError("target contract does not admit the developmental operation set")
    budgets = {
        "operations": budget_counts["operations"],
        "tensor_bytes": total_tensor_bytes,
        "units": budget_counts["units"],
        "edges": budget_counts["edges"],
        "attachments": budget_counts["attachments"],
    }
    body = {
        "target": expected_target,
        "requirements": requirements,
        "budgets": budgets,
        "tensors": emitted_tensors,
        "entrypoint": {"name": "develop", "operations": operations},
    }
    core = {
        "format": "brainc.development-module",
        "version": 1,
        "producer": deepcopy(COMPILER),
        "sources": {
            "sequence": _source_binding(source),
            "provider_manifest": {"artifact_sha256": manifest["artifact_sha256"]},
            "prediction_request": {"artifact_sha256": request["artifact_sha256"]},
            "prediction_response": {"artifact_sha256": response["artifact_sha256"]},
            "lowering_policy": {"artifact_sha256": policy["artifact_sha256"]},
            "target_contract": {"artifact_sha256": target.artifact_sha256},
        },
        "module": body,
        "module_sha256": digest(body),
    }
    artifact = seal(core)
    # The artifact returned by compile_module must be writable and reloadable by
    # every conforming consumer; reject oversize inline expansion before return.
    pretty_bytes(artifact)
    return artifact


__all__ = ["COMPILER", "CompilerError", "compile_module"]

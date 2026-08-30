"""Typed lowering from validated DNA/provider artifacts to a state-program IR."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

from ._canonical import SAFE_INTEGER, ContractError, artifact_digest, digest, keys, load, number, save_artifact, sha256, text
from .provider import validate_binding


class CompilerError(ContractError):
    """Compilation or policy validation failed."""


COMPILER = {
    "name": "brainc",
    "version": "0.4.0",
    "passes": [
        "validate-sequence-source",
        "bind-provider-contract",
        "type-check-response",
        "validate-lowering-policy",
        "lower-state-expressions",
        "emit-state-program",
    ],
}
STATE_TYPES = {"number", "integer", "boolean"}


def _string_array(value: Any, label: str) -> list[str]:
    if type(value) is not list:
        raise CompilerError(f"{label} must be an array")
    result = [text(item, f"{label}[{index}]") for index, item in enumerate(value)]
    if len(set(result)) != len(result):
        raise CompilerError(f"{label} must not contain duplicates")
    return result


def _parse_transform(value: Any, label: str, state_type: str) -> dict[str, Any] | None:
    if state_type == "boolean":
        if value is not None:
            raise CompilerError(f"{label} must be null for boolean states")
        return None
    item = keys(value, {"scale", "offset", "clamp", "rounding"}, label)
    scale = float(number(item["scale"], f"{label}.scale"))
    offset = float(number(item["offset"], f"{label}.offset"))
    clamp = item["clamp"]
    normalized_clamp = None
    if clamp is not None:
        clamp = keys(clamp, {"minimum", "maximum"}, f"{label}.clamp")
        minimum = float(number(clamp["minimum"], f"{label}.clamp.minimum"))
        maximum = float(number(clamp["maximum"], f"{label}.clamp.maximum"))
        if minimum > maximum:
            raise CompilerError(f"{label}.clamp minimum exceeds maximum")
        normalized_clamp = {"minimum": minimum, "maximum": maximum}
    expected_rounding = "nearest" if state_type == "integer" else None
    if item["rounding"] != expected_rounding:
        raise CompilerError(f"{label}.rounding must be {expected_rounding!r}")
    return {"scale": scale, "offset": offset, "clamp": normalized_clamp, "rounding": expected_rounding}


def load_policy(path: str | Path) -> dict[str, Any]:
    payload, _ = load(path, "lowering policy")
    keys(payload, {"format", "version", "id", "target", "states", "links", "ports", "artifact_sha256"}, "lowering policy")
    if payload["format"] != "brainc.lowering-policy" or type(payload["version"]) is not int or payload["version"] != 1:
        raise CompilerError("unsupported lowering policy format or version")
    text(payload["id"], "policy.id")
    target = keys(payload["target"], {"name", "version"}, "policy.target")
    text(target["name"], "policy.target.name"); text(target["version"], "policy.target.version")
    states = payload["states"]
    if type(states) is not list or not states:
        raise CompilerError("policy.states must be a non-empty array")
    state_ids: list[str] = []
    for index, value in enumerate(states):
        item = keys(value, {"id", "type", "from_output", "transform"}, f"policy.states[{index}]")
        state_id = text(item["id"], f"policy.states[{index}].id")
        state_type = item["type"]
        if type(state_type) is not str or state_type not in STATE_TYPES:
            raise CompilerError(f"policy.states[{index}].type is unsupported")
        text(item["from_output"], f"policy.states[{index}].from_output")
        _parse_transform(item["transform"], f"policy.states[{index}].transform", state_type)
        state_ids.append(state_id)
    if len(set(state_ids)) != len(state_ids):
        raise CompilerError("policy state ids must be unique")
    links = payload["links"]
    if type(links) is not list:
        raise CompilerError("policy.links must be an array")
    seen_links: set[tuple[str, str, str]] = set()
    for index, value in enumerate(links):
        item = keys(value, {"source", "target", "kind"}, f"policy.links[{index}]")
        source = text(item["source"], f"policy.links[{index}].source")
        target_id = text(item["target"], f"policy.links[{index}].target")
        kind = text(item["kind"], f"policy.links[{index}].kind")
        if source not in state_ids or target_id not in state_ids or source == target_id:
            raise CompilerError(f"policy.links[{index}] has invalid state endpoints")
        edge = (source, target_id, kind)
        if edge in seen_links:
            raise CompilerError(f"policy.links[{index}] is duplicated")
        seen_links.add(edge)
    ports = keys(payload["ports"], {"inputs", "outputs"}, "policy.ports")
    for direction in ("inputs", "outputs"):
        ids = _string_array(ports[direction], f"policy.ports.{direction}")
        if any(item not in state_ids for item in ids):
            raise CompilerError(f"policy.ports.{direction} references an unknown state")
    artifact_digest(payload)
    return payload


def _resolve(state: dict[str, Any], output: dict[str, Any]) -> dict[str, Any]:
    state_type = state["type"]
    input_type = output["type"]
    input_value = output["value"]
    if state_type == "boolean":
        if input_type != "boolean":
            raise CompilerError(f"state {state['id']!r} requires a boolean provider output")
        return {
            "id": state["id"], "type": "boolean", "value": input_value,
            "derivation": {"output_id": output["id"], "input_value": input_value},
        }
    if input_type not in {"number", "integer"}:
        raise CompilerError(f"state {state['id']!r} requires a numeric provider output")
    transform = _parse_transform(state["transform"], f"state {state['id']} transform", state_type)
    assert transform is not None
    input_number = float(input_value)
    contribution = transform["scale"] * input_number
    unclamped = transform["offset"] + contribution
    if not math.isfinite(contribution) or not math.isfinite(unclamped):
        raise CompilerError(f"state {state['id']!r} arithmetic overflowed")
    resolved = unclamped
    if transform["clamp"] is not None:
        resolved = min(max(resolved, transform["clamp"]["minimum"]), transform["clamp"]["maximum"])
    if state_type == "integer":
        if resolved < 0:
            raise CompilerError(f"integer state {state['id']!r} resolves below zero")
        lower = math.floor(resolved)
        final: int | float = lower + (1 if resolved - lower >= 0.5 else 0)
        if final > SAFE_INTEGER:
            raise CompilerError(f"integer state {state['id']!r} exceeds the I-JSON safe integer range")
    else:
        final = 0.0 if resolved == 0.0 else resolved
    return {
        "id": state["id"], "type": state_type, "value": final,
        "derivation": {
            "output_id": output["id"], "input_value": input_value,
            "scale": transform["scale"], "offset": transform["offset"],
            "contribution": contribution, "unclamped_value": unclamped,
            "clamp": transform["clamp"], "rounding": transform["rounding"],
        },
    }


def compile_program(source_path: str | Path, manifest_path: str | Path, request_path: str | Path,
                    response_path: str | Path, policy_path: str | Path) -> dict[str, Any]:
    source, manifest, request, response = validate_binding(source_path, manifest_path, request_path, response_path)
    policy = load_policy(policy_path)
    outputs = {item["id"]: item for item in response["outputs"]}
    states: list[dict[str, Any]] = []
    for state in policy["states"]:
        output_id = state["from_output"]
        if output_id not in outputs:
            raise CompilerError(f"policy references unavailable provider output {output_id!r}")
        states.append(_resolve(state, outputs[output_id]))
    source_ir_field = "ir_sha256" if source["format"] == "brain01.sequence-ir" else "collection_ir_sha256"
    program_ir = {
        "target": dict(policy["target"]),
        "states": states,
        "links": [dict(item) for item in policy["links"]],
        "ports": {"inputs": list(policy["ports"]["inputs"]), "outputs": list(policy["ports"]["outputs"])},
    }
    sources = {
        "sequence": {"format": source["format"], "version": source["version"],
                     "artifact_sha256": source["artifact_sha256"], "ir_sha256": source[source_ir_field]},
        "provider_manifest": {"artifact_sha256": manifest["artifact_sha256"]},
        "prediction_request": {"artifact_sha256": request["artifact_sha256"]},
        "prediction_response": {"artifact_sha256": response["artifact_sha256"]},
        "lowering_policy": {"artifact_sha256": policy["artifact_sha256"]},
    }
    core = {
        "format": "brainc.state-program", "version": 1,
        "compiler": COMPILER, "sources": sources,
        "program_ir": program_ir, "ir_sha256": digest(program_ir),
    }
    return {**core, "artifact_sha256": digest(core)}


def _validate_program_ir(value: Any) -> None:
    program = keys(value, {"target", "states", "links", "ports"}, "program IR")
    target = keys(program["target"], {"name", "version"}, "program target")
    text(target["name"], "program target.name"); text(target["version"], "program target.version")
    states = program["states"]
    if type(states) is not list or not states:
        raise CompilerError("program states must be a non-empty array")
    state_ids: list[str] = []
    for index, raw in enumerate(states):
        item = keys(raw, {"id", "type", "value", "derivation"}, f"program states[{index}]")
        state_id = text(item["id"], f"program states[{index}].id")
        state_type = item["type"]
        if type(state_type) is not str or state_type not in STATE_TYPES:
            raise CompilerError(f"program states[{index}].type is unsupported")
        derivation = item["derivation"]
        if state_type == "boolean":
            if type(item["value"]) is not bool:
                raise CompilerError(f"program states[{index}].value must be boolean")
            derivation = keys(derivation, {"output_id", "input_value"}, f"program states[{index}].derivation")
            if type(derivation["input_value"]) is not bool:
                raise CompilerError(f"program states[{index}] boolean derivation is invalid")
        else:
            if state_type == "integer":
                if type(item["value"]) is not int or item["value"] < 0:
                    raise CompilerError(f"program states[{index}].value must be a nonnegative integer")
            else:
                number(item["value"], f"program states[{index}].value")
            derivation = keys(
                derivation,
                {"output_id", "input_value", "scale", "offset", "contribution", "unclamped_value", "clamp", "rounding"},
                f"program states[{index}].derivation",
            )
            for field in ("input_value", "scale", "offset", "contribution", "unclamped_value"):
                number(derivation[field], f"program states[{index}].derivation.{field}")
            _parse_transform(
                {key: derivation[key] for key in ("scale", "offset", "clamp", "rounding")},
                f"program states[{index}].derivation",
                state_type,
            )
        text(derivation["output_id"], f"program states[{index}].derivation.output_id")
        state_ids.append(state_id)
    if len(set(state_ids)) != len(state_ids):
        raise CompilerError("program state ids must be unique")
    links = program["links"]
    if type(links) is not list:
        raise CompilerError("program links must be an array")
    seen_links: set[tuple[str, str, str]] = set()
    for index, raw in enumerate(links):
        item = keys(raw, {"source", "target", "kind"}, f"program links[{index}]")
        edge = (
            text(item["source"], f"program links[{index}].source"),
            text(item["target"], f"program links[{index}].target"),
            text(item["kind"], f"program links[{index}].kind"),
        )
        if edge[0] not in state_ids or edge[1] not in state_ids or edge[0] == edge[1]:
            raise CompilerError(f"program links[{index}] has invalid state endpoints")
        if edge in seen_links:
            raise CompilerError(f"program links[{index}] is duplicated")
        seen_links.add(edge)
    ports = keys(program["ports"], {"inputs", "outputs"}, "program ports")
    for direction in ("inputs", "outputs"):
        ids = _string_array(ports[direction], f"program ports.{direction}")
        if any(item not in state_ids for item in ids):
            raise CompilerError(f"program ports.{direction} references an unknown state")


def load_program(path: str | Path) -> dict[str, Any]:
    payload, _ = load(path, "state program")
    keys(payload, {"format", "version", "compiler", "sources", "program_ir", "ir_sha256", "artifact_sha256"}, "state program")
    if payload["format"] != "brainc.state-program" or type(payload["version"]) is not int or payload["version"] != 1:
        raise CompilerError("unsupported state program format or version")
    if payload["compiler"] != COMPILER:
        raise CompilerError("unsupported state-program compiler identity")
    sources = keys(payload["sources"], {"sequence", "provider_manifest", "prediction_request", "prediction_response", "lowering_policy"}, "program sources")
    sequence = keys(sources["sequence"], {"format", "version", "artifact_sha256", "ir_sha256"}, "program sequence source")
    if (type(sequence["format"]) is not str or type(sequence["version"]) is not int
            or (sequence["format"], sequence["version"]) not in {("brain01.sequence-ir", 2), ("brain01.sequence-collection-ir", 1)}):
        raise CompilerError("program sequence source is unsupported")
    sha256(sequence["artifact_sha256"], "program sequence artifact"); sha256(sequence["ir_sha256"], "program sequence IR")
    for name in ("provider_manifest", "prediction_request", "prediction_response", "lowering_policy"):
        item = keys(sources[name], {"artifact_sha256"}, f"program source {name}")
        sha256(item["artifact_sha256"], f"program source {name}")
    _validate_program_ir(payload["program_ir"])
    if payload["ir_sha256"] != digest(payload["program_ir"]):
        raise CompilerError("state-program IR digest mismatch")
    artifact_digest(payload)
    return payload


def save(payload: dict[str, Any], path: str | Path) -> None:
    save_artifact(payload, path)

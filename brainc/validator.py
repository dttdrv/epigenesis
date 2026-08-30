"""Independent standard-library replay validator for the compiler chain.

This module deliberately imports no brainc compiler modules.  It validates
serialization, identities, bindings, types, lowering arithmetic, and emitted
program equality.  It cannot establish the biological truth of an external
provider's values.
"""

from __future__ import annotations

import base64
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any


SHA_RE = re.compile(r"[0-9a-f]{64}\Z")
IUPAC = frozenset("ACGTRYSWKMBDHVN")
SEQUENCE_COMPILER = {
    "name": "brainc-dna", "version": "0.2.0",
    "passes": ["parse-fasta", "parse-context", "validate-iupac", "resolve-reference",
               "canonicalize-sequence", "emit-sequence-ir"],
}
PROGRAM_COMPILER = {
    "name": "brainc", "version": "0.3.0",
    "passes": ["validate-sequence-source", "bind-provider-contract", "type-check-response",
               "validate-lowering-policy", "lower-state-expressions", "emit-state-program"],
}


class ValidationError(ValueError):
    """An independently evaluated compiler invariant failed."""


def _reject(value: str) -> None:
    raise ValidationError(f"non-finite number {value}")


def _pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result: raise ValidationError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def _decode(raw: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(raw.decode("utf-8"), object_pairs_hook=_pairs, parse_constant=_reject)
    except ValidationError: raise
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as failure:
        raise ValidationError(f"invalid {label}: {failure}") from failure
    if type(value) is not dict: raise ValidationError(f"{label} must be an object")
    return value


def _read(path: str | Path, label: str) -> tuple[dict[str, Any], bytes]:
    try: raw = Path(path).read_bytes()
    except OSError as failure: raise ValidationError(f"cannot read {label}: {failure}") from failure
    return _decode(raw, label), raw


def _canonical(value: Any) -> bytes:
    try:
        return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
                          allow_nan=False).encode("utf-8")
    except (TypeError, ValueError) as failure: raise ValidationError(str(failure)) from failure


def _digest(value: Any) -> str: return hashlib.sha256(_canonical(value)).hexdigest()
def _raw_digest(value: bytes) -> str: return hashlib.sha256(value).hexdigest()


def _keys(value: Any, expected: set[str], label: str) -> dict[str, Any]:
    if type(value) is not dict: raise ValidationError(f"{label} must be an object")
    missing, extra = sorted(expected - value.keys()), sorted(value.keys() - expected)
    if missing or extra: raise ValidationError(f"{label} keys: missing={missing or 'none'}, unknown={extra or 'none'}")
    return value


def _sha(value: Any, label: str) -> str:
    if type(value) is not str or SHA_RE.fullmatch(value) is None: raise ValidationError(f"{label} is not SHA-256")
    return value


def _text(value: Any, label: str, empty: bool = False) -> str:
    if type(value) is not str or (not empty and not value) or value != value.strip():
        raise ValidationError(f"{label} is not a valid trimmed string")
    return value


def _number(value: Any, label: str) -> float:
    if type(value) not in (int, float) or not math.isfinite(value): raise ValidationError(f"{label} is not finite numeric")
    return float(value)


def _artifact(value: dict[str, Any], label: str) -> None:
    stored = _sha(value.get("artifact_sha256"), f"{label}.artifact_sha256")
    expected = _digest({key: item for key, item in value.items() if key != "artifact_sha256"})
    if stored != expected: raise ValidationError(f"{label} artifact digest mismatch")


def _parse_fasta(raw: bytes) -> tuple[str, str, str, list[dict[str, int]]]:
    try: lines = raw.decode("utf-8").splitlines()
    except UnicodeDecodeError as failure: raise ValidationError("FASTA must be UTF-8") from failure
    if not lines or not lines[0].startswith(">"): raise ValidationError("FASTA must begin with a defline")
    head = lines[0][1:].split(maxsplit=1)
    if not head: raise ValidationError("FASTA record id is empty")
    record_id, description = head[0], head[1] if len(head) == 2 else ""
    parts: list[str] = []; segments: list[dict[str, int]] = []; offset = 0
    for line_number, line in enumerate(lines[1:], 2):
        if line.startswith(">"): raise ValidationError("single-record FASTA contains another defline")
        if not line: continue
        if any(character.isspace() or character.upper() not in IUPAC for character in line):
            raise ValidationError(f"invalid IUPAC sequence at FASTA line {line_number}")
        parts.append(line); segments.append({"line": line_number, "normalized_start": offset, "length": len(line)})
        offset += len(line)
    sequence = "".join(parts)
    if not sequence: raise ValidationError("FASTA sequence is empty")
    return record_id, description, sequence, segments


def _context(raw: bytes, record_id: str, sequence_length: int) -> tuple[Any, Any]:
    value = _decode(raw, "sequence context")
    _keys(value, {"format", "version", "record_id", "reference", "provenance"}, "sequence context")
    if value["format"] != "brain01.sequence-context" or value["version"] != 1 or value["record_id"] != record_id:
        raise ValidationError("sequence context identity mismatch")
    reference = value["reference"]
    if reference is not None:
        reference = _keys(reference, {"assembly", "contig", "start", "end", "coordinate_system", "orientation", "aliases"}, "context reference")
        _text(reference["assembly"], "reference assembly"); _text(reference["contig"], "reference contig")
        if type(reference["start"]) is not int or type(reference["end"]) is not int or reference["start"] < 0 or reference["end"] - reference["start"] != sequence_length:
            raise ValidationError("context reference span mismatch")
        if reference["coordinate_system"] != "0-based-half-open" or reference["orientation"] not in ("forward", "reverse"):
            raise ValidationError("context coordinate convention mismatch")
        if type(reference["aliases"]) is not list or any(type(item) is not str or not item for item in reference["aliases"]):
            raise ValidationError("context aliases invalid")
    provenance = value["provenance"]
    if type(provenance) is not list: raise ValidationError("context provenance must be an array")
    return reference, provenance


def _validate_sequence(value: dict[str, Any], raw: bytes, context_raw: bytes | None) -> None:
    _keys(value, {"format", "version", "compiler", "inputs", "sequence_ir", "ir_sha256", "source_map", "artifact_sha256"}, "Sequence IR")
    if value["format"] != "brain01.sequence-ir" or value["version"] != 2 or value["compiler"] != SEQUENCE_COMPILER:
        raise ValidationError("Sequence IR format/compiler mismatch")
    record_id, description, sequence, segments = _parse_fasta(raw)
    inputs = _keys(value["inputs"], {"fasta_sha256", "context_sha256"}, "sequence inputs")
    if inputs["fasta_sha256"] != _raw_digest(raw): raise ValidationError("FASTA digest mismatch")
    reference: Any = None; provenance: Any = []
    if context_raw is None:
        if inputs["context_sha256"] is not None: raise ValidationError("unexpected context binding")
    else:
        if inputs["context_sha256"] != _raw_digest(context_raw): raise ValidationError("context digest mismatch")
        reference, provenance = _context(context_raw, record_id, len(sequence))
    canonical = sequence.upper()
    ir = {
        "record_id": record_id, "description": description, "sequence": sequence,
        "sequence_sha256": _raw_digest(sequence.encode("ascii")),
        "canonical_sha256": _raw_digest(canonical.encode("ascii")),
        "refget_id": "SQ." + base64.urlsafe_b64encode(hashlib.sha512(canonical.encode("ascii")).digest()[:24]).decode("ascii").rstrip("="),
        "reference": reference, "provenance": provenance,
    }
    if value["sequence_ir"] != ir: raise ValidationError("Sequence IR does not replay from FASTA/context")
    if value["source_map"] != {"sequence_segments": segments}: raise ValidationError("Sequence IR source map mismatch")
    if value["ir_sha256"] != _digest(ir): raise ValidationError("Sequence IR body digest mismatch")
    _artifact(value, "Sequence IR")


def _validate_manifest(value: dict[str, Any]) -> None:
    _keys(value, {"format", "version", "provider", "model_identity", "accepts", "outputs", "artifact_sha256"}, "manifest")
    if value["format"] != "brainc.provider-manifest" or value["version"] != 1: raise ValidationError("manifest version mismatch")
    provider = _keys(value["provider"], {"name", "version"}, "provider"); _text(provider["name"], "provider name"); _text(provider["version"], "provider version")
    model = _keys(value["model_identity"], {"kind", "value"}, "model identity")
    if model["kind"] == "content-sha256": _sha(model["value"], "model value")
    elif model["kind"] == "opaque": _text(model["value"], "model value")
    else: raise ValidationError("model identity kind unsupported")
    if type(value["accepts"]) is not list or not value["accepts"]: raise ValidationError("manifest accepts invalid")
    _output_contracts(value["outputs"], "manifest outputs")
    _artifact(value, "manifest")


def _output_contracts(values: Any, label: str, response: bool = False) -> list[dict[str, Any]]:
    if type(values) is not list or not values: raise ValidationError(f"{label} must be a non-empty array")
    result: list[dict[str, Any]] = []; seen: set[str] = set()
    expected = {"id", "type", "unit", "value"} if response else {"id", "type", "unit"}
    for index, raw in enumerate(values):
        item = _keys(raw, expected, f"{label}[{index}]")
        output_id = _text(item["id"], f"{label}[{index}].id")
        if output_id in seen: raise ValidationError(f"{label} ids duplicated")
        seen.add(output_id)
        kind = item["type"]
        if kind not in ("number", "integer", "boolean"): raise ValidationError(f"{label} type unsupported")
        if item["unit"] is not None: _text(item["unit"], f"{label} unit")
        if response:
            val = item["value"]
            if kind == "boolean" and type(val) is not bool: raise ValidationError(f"{label} boolean value mismatch")
            if kind == "integer" and type(val) is not int: raise ValidationError(f"{label} integer value mismatch")
            if kind == "number": _number(val, f"{label} numeric value")
        result.append(item)
    return result


def _validate_request(value: dict[str, Any]) -> None:
    _keys(value, {"format", "version", "source", "provider_manifest_sha256", "requested_outputs", "artifact_sha256"}, "request")
    if value["format"] != "brainc.prediction-request" or value["version"] != 1: raise ValidationError("request version mismatch")
    source = _keys(value["source"], {"format", "version", "artifact_sha256", "ir_sha256"}, "request source")
    if (source["format"], source["version"]) not in (("brain01.sequence-ir", 2), ("brain01.sequence-collection-ir", 1)):
        raise ValidationError("request source unsupported")
    _sha(source["artifact_sha256"], "request source artifact"); _sha(source["ir_sha256"], "request source IR")
    _sha(value["provider_manifest_sha256"], "request manifest"); _output_contracts(value["requested_outputs"], "request outputs")
    _artifact(value, "request")


def _validate_response(value: dict[str, Any]) -> None:
    _keys(value, {"format", "version", "request_artifact_sha256", "provider", "model_identity", "outputs", "artifact_sha256"}, "response")
    if value["format"] != "brainc.prediction-response" or value["version"] != 1: raise ValidationError("response version mismatch")
    _sha(value["request_artifact_sha256"], "response request")
    provider = _keys(value["provider"], {"name", "version"}, "response provider"); _text(provider["name"], "response provider name"); _text(provider["version"], "response provider version")
    model = _keys(value["model_identity"], {"kind", "value"}, "response model")
    if model["kind"] == "content-sha256": _sha(model["value"], "response model value")
    elif model["kind"] == "opaque": _text(model["value"], "response model value")
    else: raise ValidationError("response model identity unsupported")
    _output_contracts(value["outputs"], "response outputs", True); _artifact(value, "response")


def _validate_policy(value: dict[str, Any]) -> None:
    _keys(value, {"format", "version", "id", "target", "states", "links", "ports", "artifact_sha256"}, "policy")
    if value["format"] != "brainc.lowering-policy" or value["version"] != 1: raise ValidationError("policy version mismatch")
    _text(value["id"], "policy id"); target = _keys(value["target"], {"name", "version"}, "policy target"); _text(target["name"], "target name"); _text(target["version"], "target version")
    if type(value["states"]) is not list or not value["states"]: raise ValidationError("policy states invalid")
    state_ids: list[str] = []
    for index, raw in enumerate(value["states"]):
        item = _keys(raw, {"id", "type", "from_output", "transform"}, f"policy state {index}")
        state_ids.append(_text(item["id"], f"policy state {index} id")); _text(item["from_output"], f"policy state {index} output")
        if item["type"] == "boolean":
            if item["transform"] is not None: raise ValidationError("boolean transform must be null")
        elif item["type"] in ("number", "integer"):
            transform = _keys(item["transform"], {"scale", "offset", "clamp", "rounding"}, f"state {index} transform")
            _number(transform["scale"], "transform scale"); _number(transform["offset"], "transform offset")
            if transform["clamp"] is not None:
                clamp = _keys(transform["clamp"], {"minimum", "maximum"}, "transform clamp")
                if _number(clamp["minimum"], "clamp min") > _number(clamp["maximum"], "clamp max"): raise ValidationError("clamp inverted")
            if transform["rounding"] != ("nearest" if item["type"] == "integer" else None): raise ValidationError("rounding/type mismatch")
        else: raise ValidationError("state type unsupported")
    if len(set(state_ids)) != len(state_ids): raise ValidationError("state ids duplicated")
    if type(value["links"]) is not list: raise ValidationError("links must be an array")
    for index, raw in enumerate(value["links"]):
        item = _keys(raw, {"source", "target", "kind"}, f"link {index}")
        if item["source"] not in state_ids or item["target"] not in state_ids or item["source"] == item["target"]: raise ValidationError("link endpoint invalid")
        _text(item["kind"], f"link {index} kind")
    ports = _keys(value["ports"], {"inputs", "outputs"}, "ports")
    for direction in ("inputs", "outputs"):
        if type(ports[direction]) is not list or any(item not in state_ids for item in ports[direction]): raise ValidationError(f"{direction} ports invalid")
    _artifact(value, "policy")


def _expected_program(source: dict[str, Any], manifest: dict[str, Any], request: dict[str, Any], response: dict[str, Any], policy: dict[str, Any]) -> dict[str, Any]:
    source_ir = source["ir_sha256"] if source["format"] == "brain01.sequence-ir" else source["collection_ir_sha256"]
    expected_binding = {"format": source["format"], "version": source["version"], "artifact_sha256": source["artifact_sha256"], "ir_sha256": source_ir}
    if request["source"] != expected_binding or request["provider_manifest_sha256"] != manifest["artifact_sha256"]: raise ValidationError("request binding mismatch")
    if response["request_artifact_sha256"] != request["artifact_sha256"] or response["provider"] != manifest["provider"] or response["model_identity"] != manifest["model_identity"]: raise ValidationError("response binding mismatch")
    contracts = [{key: item[key] for key in ("id", "type", "unit")} for item in response["outputs"]]
    if contracts != request["requested_outputs"]: raise ValidationError("response contract mismatch")
    outputs = {item["id"]: item for item in response["outputs"]}; states: list[dict[str, Any]] = []
    for state in policy["states"]:
        if state["from_output"] not in outputs: raise ValidationError("policy output unavailable")
        output = outputs[state["from_output"]]
        if state["type"] == "boolean":
            if output["type"] != "boolean": raise ValidationError("boolean lowering type mismatch")
            states.append({"id": state["id"], "type": "boolean", "value": output["value"], "derivation": {"output_id": output["id"], "input_value": output["value"]}})
            continue
        if output["type"] not in ("number", "integer"): raise ValidationError("numeric lowering type mismatch")
        transform = state["transform"]; source_value = float(output["value"])
        contribution = float(transform["scale"]) * source_value; unclamped = float(transform["offset"]) + contribution; resolved = unclamped
        if transform["clamp"] is not None: resolved = min(max(resolved, float(transform["clamp"]["minimum"])), float(transform["clamp"]["maximum"]))
        final: int | float = math.floor(resolved + 0.5) if state["type"] == "integer" else (0.0 if resolved == 0.0 else resolved)
        if state["type"] == "integer" and resolved < 0: raise ValidationError("integer lowering below zero")
        states.append({"id": state["id"], "type": state["type"], "value": final,
                       "derivation": {"output_id": output["id"], "input_value": output["value"], "scale": float(transform["scale"]), "offset": float(transform["offset"]), "contribution": contribution, "unclamped_value": unclamped, "clamp": transform["clamp"], "rounding": transform["rounding"]}})
    ir = {"target": policy["target"], "states": states, "links": policy["links"], "ports": policy["ports"]}
    sources = {"sequence": expected_binding, "provider_manifest": {"artifact_sha256": manifest["artifact_sha256"]},
               "prediction_request": {"artifact_sha256": request["artifact_sha256"]}, "prediction_response": {"artifact_sha256": response["artifact_sha256"]},
               "lowering_policy": {"artifact_sha256": policy["artifact_sha256"]}}
    core = {"format": "brainc.state-program", "version": 1, "compiler": PROGRAM_COMPILER, "sources": sources, "program_ir": ir, "ir_sha256": _digest(ir)}
    return {**core, "artifact_sha256": _digest(core)}


def validate_chain(*, fasta: str | Path, sequence: str | Path, manifest: str | Path,
                   request: str | Path, response: str | Path, policy: str | Path,
                   program: str | Path, context: str | Path | None = None) -> dict[str, Any]:
    named = {"fasta": fasta, "sequence": sequence, "manifest": manifest, "request": request,
             "response": response, "policy": policy, "program": program}
    if context is not None: named["context"] = context
    checks: list[dict[str, Any]] = []
    hashes: dict[str, str] = {}
    try:
        fasta_raw = Path(fasta).read_bytes(); hashes["fasta"] = _raw_digest(fasta_raw)
        context_raw = None if context is None else Path(context).read_bytes()
        if context_raw is not None: hashes["context"] = _raw_digest(context_raw)
        source, source_raw = _read(sequence, "sequence artifact"); hashes["sequence"] = _raw_digest(source_raw)
        if source.get("format") != "brain01.sequence-ir": raise ValidationError("full replay currently requires a single-record Sequence IR")
        _validate_sequence(source, fasta_raw, context_raw); checks.append({"name": "sequence-replay", "passed": True})
        manifest_value, raw = _read(manifest, "manifest"); hashes["manifest"] = _raw_digest(raw); _validate_manifest(manifest_value)
        request_value, raw = _read(request, "request"); hashes["request"] = _raw_digest(raw); _validate_request(request_value)
        response_value, raw = _read(response, "response"); hashes["response"] = _raw_digest(raw); _validate_response(response_value)
        checks.append({"name": "provider-contract-and-binding", "passed": True, "scope": "identity/types/bindings only; not biological truth"})
        policy_value, raw = _read(policy, "policy"); hashes["policy"] = _raw_digest(raw); _validate_policy(policy_value)
        program_value, raw = _read(program, "program"); hashes["program"] = _raw_digest(raw); _artifact(program_value, "program")
        expected = _expected_program(source, manifest_value, request_value, response_value, policy_value)
        if _canonical(program_value) != _canonical(expected): raise ValidationError("program differs from independent lowering replay")
        checks.append({"name": "lowering-and-emission-replay", "passed": True})
    except (OSError, ValidationError) as failure:
        checks.append({"name": "compiler-chain", "passed": False, "detail": str(failure)})
    valid = all(item["passed"] for item in checks) and len(checks) == 3
    core = {"format": "brainc.validation-report", "version": 1, "stage": "compiler-chain",
            "inputs": {key: str(value) for key, value in sorted(named.items())},
            "input_sha256": dict(sorted(hashes.items())), "valid": valid,
            "passed_checks": sum(item["passed"] for item in checks), "total_checks": len(checks),
            "checks": checks,
            "trust_boundary": "Provider schema, identity, request binding, output types, and artifact integrity are checked. Provider execution and biological correctness are not replayed."}
    return {**core, "report_sha256": _digest(core)}


def save_report(report: dict[str, Any], path: str | Path) -> None:
    destination = Path(path); destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")

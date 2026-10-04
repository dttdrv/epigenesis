"""Published neural-fate ODEs driven by an explicitly artificial DNA recipe.

This example owns its parameter interpretation and target execution. Epigenesis
only compiles and validates. Natural DNA inference is not implemented here.
Equations: Balaskas et al. 2012, doi:10.1016/j.cell.2011.10.047, equations 1-3;
Panovska-Griffiths et al. 2013, doi:10.1098/rsif.2012.0826, Appendix C.
"""

from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import json
import math
from pathlib import Path
import struct
import sys
from typing import NamedTuple

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from brainc._io import load_json_object, read_regular_file
from brainc.development_bundle import compile_development
from brainc.source import RAW_PROFILE, compile_source
from brainc.v2 import inline_storage, make_development_request, pack, policy_artifact, save, target_artifact
from brainc.v2._common import seal
from brainc.v2.target import DEV_DOMAIN, OP_RULE_ATTACH, OP_UNIT_CREATE
from brainc.validator_development import load_development_bundle_directory, validate_development_paths


class Parameters(NamedTuple):
    alpha: float
    beta: float
    gamma: float
    h1: float
    h2: float
    h3: float
    h4: float
    h5: float
    k1: float
    k2: float
    k3: float
    ncrit_p: float
    ocrit_p: float
    ncrit_o: float
    ocrit_n: float
    pcrit_n: float


PARAMETER_STORAGE = struct.Struct("<" + "d" * len(Parameters._fields))
RECIPE_BASES = 4 * PARAMETER_STORAGE.size
TARGET_ID = "org.epigenesis.example.neural-fate"
RULE_ID = TARGET_ID + ".hill"
MAX_STEPS = 1_000_000  # explicit per-trajectory resource ceiling for this example.


def _validate_parameters(parameters: Parameters) -> None:
    for name, value in parameters._asdict().items():
        if type(value) not in (int, float) or not math.isfinite(value):
            raise ValueError(f"{name} must be finite")
        if value < 0 or (value == 0 and name not in ("alpha", "beta", "gamma")):
            raise ValueError(f"{name} must be positive (production rates may be zero)")


def load_models() -> dict[str, Parameters]:
    document = load_json_object(read_regular_file(Path(__file__).with_name("models.json")), "published models")
    result = {name: Parameters(**values) for name, values in document["models"].items()}
    for parameters in result.values():
        _validate_parameters(parameters)
    return result


def encode_recipe(parameters: Parameters) -> bytes:
    _validate_parameters(parameters)
    return bytes(b"ACGT"[(byte >> shift) & 3]
                 for byte in PARAMETER_STORAGE.pack(*parameters) for shift in (6, 4, 2, 0))


def decode_recipe(sequence: bytes) -> Parameters:
    if len(sequence) != RECIPE_BASES or set(sequence) - set(b"ACGT"):
        raise ValueError(f"artificial recipe requires exactly {RECIPE_BASES} uppercase A/C/G/T bases")
    digits = sequence.translate(bytes.maketrans(b"ACGT", bytes(range(4))))
    raw = bytes((a << 6) | (b << 4) | (c << 2) | d for a, b, c, d in zip(*[iter(digits)] * 4))
    parameters = Parameters(*PARAMETER_STORAGE.unpack(raw))
    _validate_parameters(parameters)
    return parameters


def derivative(parameters: Parameters, state: tuple[float, ...], signal: float) -> tuple[float, float, float]:
    p, o, n = state
    c = parameters
    activation = signal / (1 + signal)
    return (
        c.alpha / (1 + (n / c.ncrit_p)**c.h1 + (o / c.ocrit_p)**c.h2) - c.k1 * p,
        c.beta * activation / (1 + (n / c.ncrit_o)**c.h3) - c.k2 * o,
        c.gamma * activation / (1 + (o / c.ocrit_n)**c.h4 + (p / c.pcrit_n)**c.h5) - c.k3 * n,
    )


def _state(values) -> tuple[float, float, float]:
    values = tuple(values)
    if len(values) != 3 or any(type(x) not in (int, float) or not math.isfinite(x) or x < 0 for x in values):
        raise ValueError("state must contain three finite nonnegative concentrations; reduce dt if integration left the domain")
    return values


def simulate(parameters: Parameters, schedule, *, initial=None, dt: float = 0.02) -> list[tuple[float, ...]]:
    _validate_parameters(parameters)
    if type(dt) not in (int, float) or not math.isfinite(dt) or dt <= 0:
        raise ValueError("dt must be finite and positive")
    segments = []
    total_steps = 0
    end = 0.0
    for duration, signal in schedule:
        if any(type(x) not in (int, float) or not math.isfinite(x) for x in (duration, signal)) or duration <= 0 or signal < 0:
            raise ValueError("schedule requires finite positive durations and nonnegative signals")
        count = duration / dt
        if not math.isfinite(count) or count > MAX_STEPS:
            raise ValueError("trajectory exceeds step budget")
        count = max(1, math.ceil(count))
        total_steps += count
        if total_steps > MAX_STEPS:
            raise ValueError("trajectory exceeds step budget")
        next_end = end + duration
        if not math.isfinite(next_end) or next_end <= end or duration / count < math.ulp(next_end):
            raise ValueError("schedule time or time step cannot be represented finitely and distinctly")
        segments.append((duration, signal, count))
        end = next_end
    if not segments:
        raise ValueError("signal schedule is empty")
    state = _state((parameters.alpha / parameters.k1, 0.0, 0.0) if initial is None else initial)
    trajectory = [(0.0, *state)]
    start = 0.0
    for duration, signal, count in segments:
        step = duration / count
        for index in range(count):
            a = derivative(parameters, state, signal)
            b = derivative(parameters, _state(x + step * dx / 2 for x, dx in zip(state, a)), signal)
            c = derivative(parameters, _state(x + step * dx / 2 for x, dx in zip(state, b)), signal)
            d = derivative(parameters, _state(x + step * dx for x, dx in zip(state, c)), signal)
            state = _state(x + step * (u + 2*v + 2*w + z) / 6 for x, u, v, w, z in zip(state, a, b, c, d))
            trajectory.append((start + duration * ((index + 1) / count), *state))
        start += duration
    return trajectory


def _target() -> dict:
    fields = [{"id": name, "type": {"dtype": "f64", "shape": [size]}, "unit": None,
               "mutability": mutability, "numeric": {"kind": "float"}}
              for name, size, mutability in (("parameters", len(Parameters._fields), "constant"), ("state", 3, "state"))]
    return target_artifact({
        "id": TARGET_ID, "abi_major": 1, "opsets": [{"domain": DEV_DOMAIN, "version": 1}],
        "unit_schemas": [{"id": "progenitor", "fields": fields}], "edge_schemas": [], "ports": [],
        "rules": [{
            "id": RULE_ID, "version": 1,
            "semantics_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "subject": "unit", "schema": "progenitor", "phase": "UPDATE", "triggers": ["step"],
            "parameters": [],
            "reads": [{"scope": "event", "field": "signal"}, {"scope": "subject", "field": "parameters"}, {"scope": "subject", "field": "state"}],
            "writes": [{"scope": "subject", "field": "state"}],
        }],
    })


def _operations() -> list[dict]:
    return [
        {"id": "create.cell", "op": OP_UNIT_CREATE, "version": 1, "schema": "progenitor", "count": "count",
         "initializers": [{"field": name, "tensor": name} for name in ("parameters", "state")]},
        {"id": "attach.dynamics", "op": OP_RULE_ATTACH, "version": 1, "rule": RULE_ID,
         "rule_version": 1, "subject": "create.cell", "parameters": []},
    ]


def compile_recipe(recipe: Path, output: Path) -> None:
    sequence = read_regular_file(recipe, maximum_bytes=RECIPE_BASES, label="artificial recipe")
    parameters = decode_recipe(sequence)
    output.mkdir()
    snapshot = output / "recipe.dna"
    snapshot.write_bytes(sequence)
    source = compile_source(RAW_PROFILE, {"sequence": snapshot}, parameters={"record_id": "synthetic-neural-recipe"})
    source.save(output / "source")
    values = (("count", "u64", [], [1]),
              ("parameters", "f64", [len(parameters)], list(parameters)),
              ("state", "f64", [3], [parameters.alpha / parameters.k1, 0.0, 0.0]))
    outputs = [{"id": name, "type": {"dtype": dtype, "shape": shape}, "unit": None, "axes": [None] * len(shape),
                "storage": inline_storage(pack(dtype, numbers))} for name, dtype, shape, numbers in values]
    target = _target()
    save(target, output / "target.json")
    manifest = seal({
        "format": "brainc.provider-manifest", "version": 2, "provider": {"name": TARGET_ID, "version": "1"},
        "model_identity": {"kind": "content-sha256", "value": target["contract"]["rules"][0]["semantics_sha256"]},
        "accepts": [f"brainc.source-descriptor/v1;profile={RAW_PROFILE}"],
        "outputs": [{key: tensor[key] for key in ("id", "type", "unit", "axes")} for tensor in outputs],
    })
    save(manifest, output / "manifest.json")
    request = make_development_request(source, output / "manifest.json", [tensor["id"] for tensor in outputs])
    save(request, output / "request.json")
    save(seal({"format": "brainc.prediction-response", "version": 2,
               "request_artifact_sha256": request["artifact_sha256"], "provider": manifest["provider"],
               "model_identity": manifest["model_identity"], "outputs": outputs}), output / "response.json")
    policy = policy_artifact(TARGET_ID,
        {"id": TARGET_ID, "abi_major": 1, "contract_sha256": target["contract_sha256"]},
        [{"id": tensor["id"], "from_output": tensor["id"]} for tensor in outputs], _operations())
    save(policy, output / "policy.json")
    compile_development(source, output / "manifest.json", output / "request.json", output / "response.json",
                        output / "policy.json", output / "target.json").save(output / "development")


def load_compiled(directory: Path) -> tuple[Parameters, tuple[float, ...], dict]:
    validation = validate_development_paths(directory / "source", directory / "development",
                                             source_inputs={"sequence": directory / "recipe.dna"})
    if not validation["valid"]:
        raise ValueError(f"compiled recipe failed validation: {validation}")
    _, artifacts = load_development_bundle_directory(directory / "development")
    module = artifacts["development_module"]["module"]
    if artifacts["target_contract"] != _target() or module["entrypoint"]["operations"] != _operations():
        raise ValueError("unsupported target semantics or operations")
    tensors = {tensor["id"]: base64.b64decode(tensor["storage"]["data"], validate=True) for tensor in module["tensors"]}
    if struct.unpack("<Q", tensors["count"]) != (1,):
        raise ValueError("this example executes exactly one progenitor")
    parameters = Parameters(*PARAMETER_STORAGE.unpack(tensors["parameters"]))
    expected = decode_recipe(read_regular_file(directory / "recipe.dna", maximum_bytes=RECIPE_BASES))
    initial = struct.unpack("<3d", tensors["state"])
    if parameters != expected or initial != (expected.alpha / expected.k1, 0.0, 0.0):
        raise ValueError("interpreted parameters or initial state disagree with the recipe")
    return parameters, initial, validation


def _oscillation(trajectory: list[tuple[float, ...]]) -> dict:
    early = [row[2] for row in trajectory if 200 <= row[0] <= 300]
    late = [row[2] for row in trajectory if 400 <= row[0] <= 500]
    peaks = [trajectory[i][0] for i in range(1, len(trajectory) - 1)
             if trajectory[i][0] >= 400 and trajectory[i-1][2] < trajectory[i][2] >= trajectory[i+1][2]]
    return {"early_amplitude": max(early) - min(early), "late_amplitude": max(late) - min(late),
            "late_peaks": len(peaks), "mean_period": sum(b-a for a, b in zip(peaks, peaks[1:])) / (len(peaks)-1) if len(peaks) > 1 else None}


def run_experiment(output: Path) -> dict:
    output.mkdir()
    models = load_models()
    models["balaskas2012-double-knockout"] = models["balaskas2012"]._replace(alpha=0.0, beta=0.0)
    compiled = {}
    provenance = {}
    for name, parameters in models.items():
        recipe = output / (name + ".dna")
        recipe.write_bytes(encode_recipe(parameters))
        compile_recipe(recipe, output / name)
        loaded, initial, validation = load_compiled(output / name)
        compiled[name] = loaded, initial
        (output / name / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
        provenance[name] = load_json_object(read_regular_file(output / name / "development/bundle.json"), "compiled bundle")["artifact_sha256"]
    protocols = {
        "naive": ("balaskas2012", [(40.0, 1.0)]),
        "primed": ("balaskas2012", [(20.0, 5.0), (20.0, 1.0)]),
        "washout": ("balaskas2012", [(20.0, 5.0), (20.0, 0.0)]),
        "wild-type": ("balaskas2012", [(20.0, 0.5)]),
        "double-knockout": ("balaskas2012-double-knockout", [(20.0, 0.5)]),
        "steady": ("panovska2013-steady", [(500.0, 5.0)]),
        "oscillation": ("panovska2013-oscillating", [(500.0, 5.0)]),
    }
    trajectories = {name: simulate(compiled[model][0], schedule, initial=compiled[model][1])
                    for name, (model, schedule) in protocols.items()}
    with (output / "trajectories.csv").open("x", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("condition", "time_model_units", "Pax6", "Olig2", "Nkx2.2"))
        for name, trajectory in trajectories.items():
            writer.writerows((name, *row) for row in trajectory)
    steady = _oscillation(trajectories["steady"])
    oscillation = _oscillation(trajectories["oscillation"])
    fine = simulate(compiled["panovska2013-oscillating"][0], [(500.0, 5.0)], dt=0.01)
    refined_oscillation = _oscillation(fine)
    error = max(abs(x-y) for a, b in zip(trajectories["oscillation"], fine[::2]) for x, y in zip(a[1:], b[1:]))
    period_error = (abs(oscillation["mean_period"] - refined_oscillation["mean_period"])
                    if oscillation["mean_period"] is not None and refined_oscillation["mean_period"] is not None else None)
    checks = [
        ("history changes fate at identical final signal", trajectories["naive"][-1][3] < 1 < trajectories["primed"][-1][3]),
        ("withdrawal extinguishes Nkx2.2", trajectories["washout"][-1][3] < 1e-7),
        ("double gene loss increases low-signal Nkx2.2", trajectories["double-knockout"][-1][3] > 1 > trajectories["wild-type"][-1][3]),
        ("2013 steady case settles", steady["late_amplitude"] < 1e-4),
        ("2013 oscillator persists in two late windows", min(oscillation["early_amplitude"], oscillation["late_amplitude"]) > 0.1 and oscillation["late_peaks"] >= 3),
        ("oscillations survive timestep refinement", error < 1e-4),
        ("oscillation period agrees within sampling resolution", period_error is not None and period_error <= 2 * (0.02 + 0.01)),
    ]
    with (output / "trajectories.csv").open("rb") as stream:
        trajectory_digest = hashlib.file_digest(stream, "sha256").hexdigest()
    report = {
        "status": "published-model-reproduced" if all(passed for _, passed in checks) else "reproduction-failed",
        "natural_dna_inference": False, "independent_biological_validation": False,
        "evidence_scope": "qualitative published-model reproduction; time and concentration are model units",
        "implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "model_document_sha256": hashlib.sha256(Path(__file__).with_name("models.json").read_bytes()).hexdigest(),
        "compiled_bundles": provenance, "protocols": protocols, "dt": 0.02, "refined_dt": 0.01,
        "trajectories_sha256": trajectory_digest,
        "sample_counts": {name: len(rows) for name, rows in trajectories.items()},
        "maximum_refinement_error": error, "period_refinement_error": period_error,
        "steady": steady, "oscillation": oscillation,
        "endpoints": {name: trajectory[-1] for name, trajectory in trajectories.items()},
        "checks": [{"claim": claim, "passed": passed} for claim, passed in checks],
    }
    (output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="new directory for compiled recipes, trajectories, and evidence")
    arguments = parser.parse_args()
    try:
        report = run_experiment(arguments.output)
    except (ValueError, OSError, OverflowError) as failure:
        print(f"experiment failed: {failure}", file=sys.stderr)
        return 1
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0 if report["status"] == "published-model-reproduced" else 1


if __name__ == "__main__":
    raise SystemExit(main())

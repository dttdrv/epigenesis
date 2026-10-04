"""Compile a shared artificial developmental recipe and execute it with CX3D.

The recipe controls a new local construction experiment. It is not a natural
genome interpreter or a reproduction of the complete 2013 cortical model.
"""

from __future__ import annotations

import argparse
import base64
import hashlib
import io
import json
import math
from pathlib import Path, PurePosixPath
import shutil
import struct
import subprocess
import sys
from typing import NamedTuple
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from brainc._io import load_json_object, read_regular_file
from brainc.development_bundle import compile_development
from brainc.source import RAW_PROFILE, compile_source
from brainc.v2 import inline_storage, make_development_request, pack, policy_artifact, save, target_artifact
from brainc.v2._common import seal
from brainc.v2.target import DEV_DOMAIN, OP_RULE_ATTACH, OP_UNIT_CREATE
from brainc.validator_development import load_development_bundle_directory, validate_development_paths


HERE = Path(__file__).resolve().parent
TARGET = "org.epigenesis.example.neural-construction"
RULE = TARGET + ".cx3d-local-program"
MAX_STEPS = 20_000
MAX_CELLS = 256
MAX_SEGMENTS = 10_000
MAX_PROTRUSIONS = 100_000
MAX_TRACE_BYTES = 64*1024*1024


class Program(NamedTuple):
    division_signal: float
    division_threshold: float
    division_volume: float
    volume_rate: float
    axon_speed: float
    dendrite_speed: float
    guidance_gain: float
    branch_rate: float
    diameter: float
    taper: float
    stop_diameter: float
    contact_spacing: float
    contact_reach: float
    recognition: float
    noise: float


STORAGE = struct.Struct("<" + "d"*len(Program._fields))
RECIPE_BASES = STORAGE.size * 4


class Backend(NamedTuple):
    classes: Path
    manifest_sha256: str


def default_program() -> Program:
    return Program(4.0, 1.0, math.pi*20**3/6, 3000.0, 80.0, 55.0, 2.0,
                   0.3, 1.0, 0.05, 0.3, 2.0, 5.0, 1.0, 0.2)


def validate_program(program: Program) -> None:
    positive = {"division_threshold", "division_volume", "diameter", "stop_diameter", "contact_spacing", "contact_reach"}
    for name, value in program._asdict().items():
        if type(value) not in (int, float) or not math.isfinite(value) or value < 0 or (name in positive and value == 0):
            raise ValueError(f"invalid program gene: {name}")
    if program.stop_diameter > program.diameter or program.recognition not in (0.0, 1.0):
        raise ValueError("unsupported diameter or recognition gene")


def encode_recipe(program: Program) -> bytes:
    validate_program(program)
    return bytes(b"ACGT"[(byte >> shift) & 3] for byte in STORAGE.pack(*program) for shift in (6, 4, 2, 0))


def decode_recipe(sequence: bytes) -> Program:
    if len(sequence) != RECIPE_BASES or set(sequence) - set(b"ACGT"):
        raise ValueError(f"artificial construction recipe requires {RECIPE_BASES} uppercase ACGT bases")
    digits = sequence.translate(bytes.maketrans(b"ACGT", bytes(range(4))))
    raw = bytes((a<<6) | (b<<4) | (c<<2) | d for a,b,c,d in zip(*[iter(digits)]*4))
    program = Program(*STORAGE.unpack(raw))
    validate_program(program)
    return program


def semantics() -> str:
    identity = hashlib.sha256()
    for name in ("construction.py", "Construction.java", "vendor/provenance.json", "vendor/compatibility.patch"):
        identity.update(name.encode()+b"\0"+read_regular_file(HERE/name, maximum_bytes=2*1024*1024))
    return identity.hexdigest()


def _target() -> dict:
    return target_artifact({"id": TARGET, "abi_major": 1, "opsets": [{"domain": DEV_DOMAIN, "version": 1}],
        "unit_schemas": [{"id": "founder", "fields": [{"id": "program", "type": {"dtype": "f64", "shape": [len(Program._fields)]},
          "unit": None, "mutability": "constant", "numeric": {"kind": "float"}}]}],
        "edge_schemas": [], "ports": [], "rules": [{"id": RULE, "version": 1, "semantics_sha256": semantics(),
          "subject": "unit", "schema": "founder", "phase": "UPDATE", "triggers": ["step"], "parameters": [],
          "reads": [{"scope": "subject", "field": "program"}], "writes": []}]})


def _operations() -> list[dict]:
    return [{"id": "create.founder", "op": OP_UNIT_CREATE, "version": 1, "schema": "founder", "count": "count",
             "initializers": [{"field": "program", "tensor": "program"}]},
            {"id": "attach.program", "op": OP_RULE_ATTACH, "version": 1, "rule": RULE,
             "rule_version": 1, "subject": "create.founder", "parameters": []}]


def compile_recipe(recipe: Path, output: Path) -> None:
    sequence = read_regular_file(recipe, maximum_bytes=RECIPE_BASES)
    program = decode_recipe(sequence)
    output.mkdir()
    snapshot = output/"recipe.dna"
    snapshot.write_bytes(sequence)
    source = compile_source(RAW_PROFILE, {"sequence": snapshot}, parameters={"record_id": "synthetic-construction-recipe"})
    source.save(output/"source")
    tensors = [{"id": name, "type": {"dtype": dtype, "shape": shape}, "unit": None,
                "axes": [None]*len(shape), "storage": inline_storage(pack(dtype, values))}
               for name,dtype,shape,values in (("count","u64",[],[1]),("program","f64",[len(program)],list(program)))]
    manifest = seal({"format": "brainc.provider-manifest", "version": 2, "provider": {"name": TARGET, "version": "1"},
                     "model_identity": {"kind": "content-sha256", "value": semantics()},
                     "accepts": [f"brainc.source-descriptor/v1;profile={RAW_PROFILE}"],
                     "outputs": [{key:tensor[key] for key in ("id","type","unit","axes")} for tensor in tensors]})
    save(manifest, output/"manifest.json")
    request = make_development_request(source, output/"manifest.json", [t["id"] for t in tensors])
    save(request, output/"request.json")
    save(seal({"format": "brainc.prediction-response", "version": 2, "request_artifact_sha256": request["artifact_sha256"],
               "provider": manifest["provider"], "model_identity": manifest["model_identity"], "outputs": tensors}), output/"response.json")
    target = _target()
    save(target, output/"target.json")
    save(policy_artifact(TARGET, {"id": TARGET, "abi_major": 1, "contract_sha256": target["contract_sha256"]},
          [{"id": t["id"], "from_output": t["id"]} for t in tensors], _operations()), output/"policy.json")
    compile_development(source, output/"manifest.json", output/"request.json", output/"response.json",
                        output/"policy.json", output/"target.json").save(output/"development")
    _, validation = load_compiled(output)
    with (output/"validation.json").open("x") as stream:
        stream.write(json.dumps(validation, indent=2)+"\n")


def load_compiled(directory: Path) -> tuple[Program, dict]:
    validation = validate_development_paths(directory/"source", directory/"development", source_inputs={"sequence": directory/"recipe.dna"})
    if not validation["valid"]:
        raise ValueError(f"construction artifact validation failed: {validation}")
    _, artifacts = load_development_bundle_directory(directory/"development")
    if artifacts["target_contract"] != _target():
        raise ValueError("unsupported construction target semantics")
    module = artifacts["development_module"]["module"]
    if module["entrypoint"]["operations"] != _operations():
        raise ValueError("unsupported construction operations")
    tensors = {tensor["id"]:base64.b64decode(tensor["storage"]["data"], validate=True) for tensor in module["tensors"]}
    if set(tensors) != {"count", "program"} or struct.unpack("<Q", tensors["count"]) != (1,):
        raise ValueError("construction must start from exactly one founder and a shared program")
    program = Program(*STORAGE.unpack(tensors["program"]))
    if program != decode_recipe(read_regular_file(directory/"recipe.dna", maximum_bytes=RECIPE_BASES)):
        raise ValueError("linked construction program disagrees with the recipe")
    validate_program(program)
    return program, validation


def build_backend(output: Path) -> Backend:
    output = output.resolve()
    provenance = load_json_object(read_regular_file(HERE/"vendor/provenance.json"), "CX3D provenance")
    archive = read_regular_file(HERE/"vendor/cx3d-0.03.zip", maximum_bytes=2*1024*1024)
    if len(archive) != provenance["source_bytes"] or hashlib.sha256(archive).hexdigest() != provenance["source_sha256"]:
        raise ValueError("CX3D source identity mismatch")
    patch = HERE/"vendor/compatibility.patch"
    patch_bytes = read_regular_file(patch)
    if hashlib.sha256(patch_bytes).hexdigest() != provenance["patch_sha256"]:
        raise ValueError("CX3D compatibility patch identity mismatch")
    compiler = shutil.which("javac")
    if compiler is None:
        raise ValueError("construction requires an installed JDK with javac and java")
    output.mkdir()
    patch = output/"compatibility.patch"
    patch.write_bytes(patch_bytes)
    source = output/"source"
    source.mkdir()
    with zipfile.ZipFile(io.BytesIO(archive)) as bundle:
        for item in bundle.infolist():
            path = PurePosixPath(item.filename)
            if path.is_absolute() or ".." in path.parts or item.file_size > 2*1024*1024:
                raise ValueError("invalid CX3D archive member")
            if path.suffix != ".java":
                continue
            destination = source.joinpath(*path.parts)
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(bundle.read(item))
    patcher = shutil.which("patch")
    if patcher is None:
        raise ValueError("construction requires the standard patch utility")
    with (output/"patch.log").open("x") as log:
        result = subprocess.run([patcher,"--batch","--fuzz=0","-p1","-i",str(patch)],cwd=source,stdout=log,stderr=subprocess.STDOUT,timeout=30)
    if result.returncode:
        raise ValueError(f"CX3D compatibility patch failed; evidence: {output/'patch.log'}")
    shutil.copyfile(HERE/"Construction.java", source/"Construction.java")
    classes = output/"classes"
    classes.mkdir()
    sources = sorted(str(path) for path in source.rglob("*.java"))
    command = [compiler, "--release", "11", "-d", str(classes), *sources]
    with (output/"build.log").open("x") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=120)
    if result.returncode:
        raise ValueError(f"CX3D build failed; evidence: {output/'build.log'}")
    (output/"build.json").write_text(json.dumps({"semantics_sha256": semantics(), "command": command,
        "source_provenance": provenance, "classes": {str(p.relative_to(classes)):hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(classes.rglob("*.class"))}}, indent=2)+"\n")
    return Backend(classes, hashlib.sha256((output/"build.json").read_bytes()).hexdigest())


def develop(program: Program, backend: Backend, output: Path, *, duration: float = 3.0, dt: float = 0.01,
            seed: int = 7, cue: float = 40.0, origin=(0.0,0.0,0.0), max_cells: int = MAX_CELLS,
            max_segments: int = MAX_SEGMENTS, max_protrusions: int = MAX_PROTRUSIONS,
            cue_switch: tuple[int,float] | None = None) -> dict:
    validate_program(program)
    if any(type(x) not in (int,float) or not math.isfinite(x) for x in (duration,dt,cue,*origin)) or duration <= 0 or dt <= 0 or len(origin) != 3:
        raise ValueError("invalid development environment or duration")
    count = duration/dt
    if not math.isfinite(count) or count > MAX_STEPS or count < 1 or type(seed) is not int or not 0 <= seed < 2**63:
        raise ValueError("development exceeds supported step or seed range")
    if cue_switch is not None and (len(cue_switch) != 2 or type(cue_switch[0]) is not int or not 1 <= cue_switch[0] <= math.ceil(count)
            or type(cue_switch[1]) not in (int,float) or not math.isfinite(cue_switch[1])):
        raise ValueError("invalid timed guidance intervention")
    for value, maximum in ((max_cells,MAX_CELLS),(max_segments,MAX_SEGMENTS),(max_protrusions,MAX_PROTRUSIONS)):
        if type(value) is not int or not 1 <= value <= maximum:
            raise ValueError("invalid construction resource budget")
    classes = backend.classes
    manifest = read_regular_file(classes.parent/"build.json", maximum_bytes=1024*1024)
    if hashlib.sha256(manifest).hexdigest() != backend.manifest_sha256:
        raise ValueError("backend manifest changed since the verified build")
    build = load_json_object(manifest, "CX3D build")
    if build["semantics_sha256"] != semantics():
        raise ValueError("construction backend semantics have changed; rebuild")
    if {str(p.relative_to(classes)) for p in classes.rglob("*") if not p.is_dir()} != set(build["classes"]):
        raise ValueError("compiled backend inventory mismatch")
    for name,digest in build["classes"].items():
        if hashlib.sha256(read_regular_file(classes/name, maximum_bytes=1024*1024)).hexdigest() != digest:
            raise ValueError("compiled backend identity mismatch")
    java = shutil.which("java")
    if java is None:
        raise ValueError("construction requires an installed Java runtime")
    output = output.resolve()
    output.mkdir()
    inputs = {**program._asdict(), "steps": math.ceil(count), "dt": duration/math.ceil(count), "seed": seed,
              "cue": cue, "origin_x": origin[0], "origin_y": origin[1], "origin_z": origin[2],
              "cue_switch_step": cue_switch[0] if cue_switch else 0, "cue_after": cue_switch[1] if cue_switch else cue,
              "max_cells": max_cells, "max_segments": max_segments, "max_protrusions": max_protrusions}
    (output/"input.json").write_text(json.dumps(inputs, indent=2, allow_nan=False)+"\n")
    (output/"input.properties").write_text("".join(f"{key}={value}\n" for key,value in inputs.items()))
    command = [java,"-ea","-Djava.awt.headless=true","-Xmx512m","-cp",str(classes),"Construction",
               str(output/"input.properties"),str(output/"trace.jsonl"),str(output/"summary.json")]
    with (output/"runtime.log").open("x") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=120)
    if result.returncode:
        raise ValueError(f"construction failed; partial evidence: {output}")
    report = verify_trace(output)
    report.update({"biological_acceptance": False, "semantics_sha256": semantics(), "command": command,
                   "backend_manifest_sha256": backend.manifest_sha256})
    (output/"report.json").write_text(json.dumps(report, indent=2)+"\n")
    return report


def verify_trace(directory: Path) -> dict:
    inputs = load_json_object(read_regular_file(directory/"input.json", maximum_bytes=65536), "construction input")
    program = Program(*(inputs[name] for name in Program._fields))
    validate_program(program)
    dt, steps = inputs["dt"], inputs["steps"]
    if type(steps) is not int or not 1 <= steps <= MAX_STEPS or not math.isfinite(dt) or dt <= 0:
        raise ValueError("invalid trace clock")
    schemas = {
        "founder": "cell position factor", "cell": "cell soma position volume",
        "segment": "segment cell parent parent_is_soma position proximal diameter axon",
        "division": "parent children factor_before factor_after volume_before volumes",
        "differentiate": "cell soma factor axon dendrite",
        "growth": "cell segment position direction prior_direction noise axon gradient speed dt",
        "branch": "cell parent children diameter_before diameter_after", "stop": "segment diameter",
        "synapse": "contact pre post pre_segment post_segment pre_position post_position reach",
        "connection": "contact pre post pre_segment post_segment pre_axial post_axial pre_length post_length", "environment": "cue",
        "subdivide": "distal proximal_segment parent", "merge": "distal proximal_segment parent",
    }
    data = read_regular_file(directory/"trace.jsonl", maximum_bytes=MAX_TRACE_BYTES)
    events = [load_json_object(line, "construction event") for line in data.splitlines()]
    if not events or events[0].get("event") != "founder":
        raise ValueError("construction trace must begin with a single founder")
    live, born, differentiated, segments = {}, set(), set(), {}
    cells, frame_segments = {}, {}
    current_step, contacts = 0, {}
    retained, pending_contacts, interventions = {}, [], 0
    active, expected_actions, actions = set(), set(), set()

    def close_frame():
        if set(cells) != set(live):
            raise ValueError("snapshot cells disagree with lineage")
        if set(frame_segments) != set(segments) or actions != expected_actions:
            raise ValueError("snapshot lacks a complete construction history")
        for segment in frame_segments.values():
            visited = set()
            current = segment
            while not current["parent_is_soma"]:
                parent = current["parent"]
                if parent in visited or parent not in frame_segments:
                    raise ValueError("segment tree has a cycle or missing parent")
                visited.add(parent)
                current = frame_segments[parent]
                if current["cell"] != segment["cell"] or current["axon"] != segment["axon"]:
                    raise ValueError("segment tree crosses cells or process types")
            if current["parent"] != cells[segment["cell"]]["soma"]:
                raise ValueError("segment is not rooted at its own soma")
        for contact in pending_contacts:
            for side, axon in (("pre",True),("post",False)):
                segment = frame_segments.get(contact[side+"_segment"])
                if segment is None or segment["cell"] != contact[side] or segment["axon"] != axon:
                    raise ValueError("synapse is not attached to the developed morphology")
                if side+"_position" in contact:
                    point = contact[side+"_position"]
                    axis = [b-a for a,b in zip(segment["proximal"],segment["position"])]
                    offset = [b-a for a,b in zip(segment["proximal"],point)]
                    length_squared = math.fsum(x*x for x in axis)
                    if length_squared == 0:
                        raise ValueError("contact belongs to a degenerate segment")
                    fraction = math.fsum(a*b for a,b in zip(axis,offset))/length_squared
                    radius = math.hypot(*(x-fraction*a for x,a in zip(offset,axis)))
                    if not -1e-10 <= fraction <= 1+1e-10 or not math.isclose(radius,segment["diameter"]/2,rel_tol=1e-9,abs_tol=1e-9):
                        raise ValueError("contact is not on its attached cylinder surface")
        pending_contacts.clear()

    for event in events:
        kind = event.get("event")
        if kind not in schemas or set(event) != {"event","step","time",*schemas[kind].split()}:
            raise ValueError("unknown or malformed construction event")
        step = event["step"]
        if type(step) is not int or step not in (current_step,current_step+1) or step > steps:
            raise ValueError("invalid construction event order")
        if step != current_step:
            close_frame()
            segments = frame_segments.copy()
            cells, frame_segments = {}, {}
            current_step = step
            expected_actions, actions = active.copy(), set()
        snapshot = kind in ("cell","segment","synapse","connection")
        expected_time = max(0, step if snapshot else step-1)*dt
        if not math.isclose(event["time"],expected_time,rel_tol=1e-10,abs_tol=1e-12):
            raise ValueError("construction event time disagrees with clock")
        for field in ("position","proximal","pre_position","post_position","direction","prior_direction","noise","gradient"):
            if field in event and (not isinstance(event[field],list) or len(event[field]) != 3 or
                    any(type(x) not in (int,float) or not math.isfinite(x) for x in event[field])):
                raise ValueError("invalid physical vector")
        if kind == "founder":
            if born or step != 0 or event["factor"] != program.division_signal:
                raise ValueError("invalid founder state")
            if event["position"] != [inputs["origin_"+axis] for axis in "xyz"]:
                raise ValueError("founder position disagrees with initial conditions")
            live[event["cell"]] = event["factor"]
            born.add(event["cell"])
        elif kind == "division":
            parent, children = event["parent"], event["children"]
            if parent not in live or parent in differentiated or len(children) != 2 or len(set(children)) != 2 or born.intersection(children):
                raise ValueError("invalid developmental lineage")
            if event["factor_before"] != live[parent] or live[parent] <= program.division_threshold:
                raise ValueError("division disagrees with the inherited program")
            if event["volume_before"] < program.division_volume or len(event["volumes"]) != 2 or any(x <= 0 for x in event["volumes"]):
                raise ValueError("division lacks the required cell volume")
            if not math.isclose(sum(event["volumes"]),event["volume_before"],rel_tol=1e-12):
                raise ValueError("division does not conserve volume")
            if event["factor_after"]*2 != event["factor_before"]:
                raise ValueError("division does not conserve the inherited factor")
            del live[parent]
            live.update(dict.fromkeys(children,event["factor_after"]))
            born.update(children)
        elif kind == "differentiate":
            cell = event["cell"]
            if cell not in live or cell in differentiated or event["factor"] != live[cell] or live[cell] > program.division_threshold:
                raise ValueError("differentiation disagrees with the inherited program")
            differentiated.add(cell)
            for channel in ("axon","dendrite"):
                identity = event[channel]
                if identity in segments:
                    raise ValueError("duplicate new neurite")
                segments[identity] = {"cell":cell,"axon":channel == "axon","parent":event["soma"],"parent_is_soma":True}
                active.add(identity)
        elif kind == "cell":
            cell = event["cell"]
            if cell not in live or cell in cells or event["volume"] <= 0:
                raise ValueError("invalid or duplicate cell snapshot")
            cells[cell] = event
        elif kind == "segment":
            if event["cell"] not in differentiated or event["segment"] in frame_segments or event["diameter"] <= 0:
                raise ValueError("invalid segment snapshot")
            if type(event["axon"]) is not bool or type(event["parent_is_soma"]) is not bool:
                raise ValueError("invalid segment type")
            known = segments.get(event["segment"])
            if known is None or any(known[key] != event[key] for key in ("cell","axon","parent","parent_is_soma")):
                raise ValueError("segment snapshot disagrees with construction history")
            frame_segments[event["segment"]] = event
        elif kind == "growth":
            known = segments.get(event["segment"])
            if known is None or event["cell"] != known["cell"] or event["axon"] != known["axon"]:
                raise ValueError("growth is not attached to a known process")
            if event["segment"] not in expected_actions or event["segment"] in actions:
                raise ValueError("unexpected growth action")
            actions.add(event["segment"])
            if event["speed"] != (program.axon_speed if event["axon"] else program.dendrite_speed) or event["dt"] != dt:
                raise ValueError("growth disagrees with the declared speed or clock")
            prior, gradient, noise = (event[key] for key in ("prior_direction","gradient","noise"))
            magnitude = math.hypot(*gradient)
            combined = [p + (program.guidance_gain*g/magnitude if magnitude else 0) + n for p,g,n in zip(prior,gradient,noise)]
            length = math.hypot(*combined)
            expected = [x/length for x in combined] if length else prior
            if any(abs(x) > program.noise for x in noise) or math.dist(event["direction"],expected) > 1e-12:
                raise ValueError("growth direction disagrees with local inputs")
        elif kind == "branch":
            parent = segments.get(event["parent"])
            if parent is None or parent["cell"] != event["cell"] or program.branch_rate == 0 or event["parent"] not in active or event["parent"] not in actions:
                raise ValueError("branch is not attached to a growing process")
            if not math.isclose(2*event["diameter_after"]**2,event["diameter_before"]**2,rel_tol=1e-12):
                raise ValueError("branch does not conserve cross-sectional area")
            children = event["children"]
            if len(children) != 2 or len(set(children)) != 2 or any(child in segments for child in children):
                raise ValueError("invalid branch children")
            segments.update({child:{**parent,"parent":event["parent"],"parent_is_soma":False} for child in children})
            active.remove(event["parent"])
            active.update(children)
        elif kind == "stop":
            if event["segment"] not in expected_actions or event["segment"] in actions or event["diameter"] > program.stop_diameter:
                raise ValueError("invalid growth stop")
            actions.add(event["segment"])
            active.remove(event["segment"])
        elif kind in ("subdivide","merge"):
            distal, proximal = event["distal"], event["proximal_segment"]
            if distal not in segments:
                raise ValueError("physical remeshing has no existing distal segment")
            if kind == "subdivide":
                if proximal in segments or segments[distal]["parent"] != event["parent"]:
                    raise ValueError("invalid physical subdivision")
                segments[proximal] = segments[distal].copy()
                segments[distal] = {**segments[distal],"parent":proximal,"parent_is_soma":False}
            else:
                if proximal not in segments or proximal in active or segments[distal]["parent"] != proximal or segments[proximal]["parent"] != event["parent"]:
                    raise ValueError("invalid physical segment merger")
                segments[distal] = {**segments[distal],"parent":event["parent"],"parent_is_soma":segments[proximal]["parent_is_soma"]}
                del segments[proximal]
        elif kind == "synapse":
            if program.recognition != 1 or event["reach"] != program.contact_reach:
                raise ValueError("contact disagrees with the recognition program")
            if event["pre"] not in live or event["post"] not in live or event["pre"] == event["post"]:
                raise ValueError("synapse has missing or identical cells")
            if math.dist(event["pre_position"],event["post_position"]) > event["reach"] + 1e-12:
                raise ValueError("synapse formed outside physical reach")
            if event["contact"] in contacts:
                raise ValueError("duplicate contact identity")
            contacts[event["contact"]] = (event["pre"],event["post"])
            pending_contacts.append(event)
        elif kind == "connection":
            if step != steps or event["contact"] in retained or contacts.get(event["contact"]) != (event["pre"],event["post"]):
                raise ValueError("final connection disagrees with formation history")
            retained[event["contact"]] = event
            pending_contacts.append(event)
        elif kind == "environment":
            if step != inputs["cue_switch_step"] or event["cue"] != inputs["cue_after"] or interventions:
                raise ValueError("undeclared guidance intervention")
            interventions += 1
    close_frame()
    if current_step != steps or interventions != bool(inputs["cue_switch_step"]) or set(retained) != set(contacts):
        raise ValueError("construction trace ends before the declared duration")
    outside = 0
    for connection in retained.values():
        for side in ("pre","post"):
            segment = frame_segments[connection[side+"_segment"]]
            length = math.dist(segment["position"],segment["proximal"])
            if not math.isclose(connection[side+"_length"],length,rel_tol=1e-12,abs_tol=1e-12):
                raise ValueError("retained attachment length disagrees with morphology")
            outside += not 0 <= connection[side+"_axial"] <= connection[side+"_length"]
    summary = load_json_object(read_regular_file(directory/"summary.json", maximum_bytes=1024*1024), "construction summary")
    if summary["initial_cells"] != 1 or summary["cells"] != len(live) or summary["segments"] != len(frame_segments) or summary["physics"] is not True or summary["diffusion"] is not True:
        raise ValueError("construction summary disagrees with lineage, morphology or physics")
    if summary["retained_sites_outside_segments"] != outside:
        raise ValueError("retained attachment diagnostic disagrees with geometry")
    length = math.fsum(math.dist(segment["position"],segment["proximal"]) for segment in frame_segments.values())
    if summary["synapses"] != len(contacts) or not math.isclose(summary["length"],length,rel_tol=1e-12,abs_tol=1e-12) or not math.isclose(summary["time"],steps*dt,rel_tol=1e-10):
        raise ValueError("construction summary disagrees with physical events")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--recipe", type=Path)
    parser.add_argument("--duration", type=float, default=3.0)
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args()
    try:
        args.output.mkdir()
        recipe = args.recipe
        if recipe is None:
            recipe = args.output/"recipe.dna"
            recipe.write_bytes(encode_recipe(default_program()))
        compile_recipe(recipe, args.output/"compiled")
        program, validation = load_compiled(args.output/"compiled")
        classes = build_backend(args.output/"backend")
        report = develop(program,classes,args.output/"development-run",duration=args.duration,seed=args.seed)
        from activity import assay, write_viewer
        activity = assay(args.output/"development-run",verify_trace)
        (args.output/"activity.json").write_text(json.dumps(activity,indent=2)+"\n")
        write_viewer(args.output/"development-run",activity,HERE/"viewer.html",args.output/"viewer.html")
        names = ("compiled/recipe.dna","compiled/development/bundle.json","compiled/validation.json",
                 "backend/build.json","development-run/input.json","development-run/trace.jsonl",
                 "development-run/summary.json","development-run/report.json","activity.json","viewer.html")
        experiment = {"construction":report,"activity_controls_pass":activity["transmission_controls_pass"],
            "compiled_valid":validation["valid"],"biological_acceptance":False,
            "evidence_sha256":{name:hashlib.sha256(read_regular_file(args.output/name,maximum_bytes=MAX_TRACE_BYTES)).hexdigest() for name in names}}
        (args.output/"experiment.json").write_text(json.dumps(experiment,indent=2)+"\n")
        if not activity["transmission_controls_pass"]:
            raise ValueError("activity transmission controls failed; evidence retained")
        print(json.dumps(experiment,indent=2))
    except (ValueError,OSError,subprocess.TimeoutExpired) as failure:
        print(f"construction experiment failed: {failure}",file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

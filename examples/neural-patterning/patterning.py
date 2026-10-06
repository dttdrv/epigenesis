"""Measured regulatory DNA profiles coupled to a published neural-patterning model."""

from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import struct
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from brainc._canonical import canonical_bytes
from brainc._io import load_json_object, read_regular_file
from brainc.development_bundle import compile_development
from brainc.sequence_collection import SequenceCollectionCompiler
from brainc.source import FASTA_PROFILE, compile_source
from brainc.v2 import inline_storage, make_development_request, pack, policy_artifact, save, target_artifact
from brainc.v2._common import seal
from brainc.v2.target import DEV_DOMAIN, OP_RULE_ATTACH, OP_UNIT_CREATE
from brainc.validator_development import (
    load_development_bundle_directory, load_source_bundle_directory, validate_development_bundle,
)

HERE = Path(__file__).resolve().parent
DATA_SHA256 = "b9b8407e8489d8f9c1207cc9b5001c031d536f801a2090ed025fe4ee3f1a578a"
IMPLEMENTATION = read_regular_file(HERE/"patterning.py", maximum_bytes=128*1024)
TARGET = "org.epigenesis.example.neural-patterning"
RULE = TARGET + ".gli-site"
MAX_SOURCE_BYTES = 1024*1024
MAX_STEPS = 1_000_000  # per-trajectory resource ceiling, not a biological limit.
DEFAULT_DT = .1
LOCAL_ERROR_TOLERANCE = 1e-10


def digest(value: dict) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def load_data() -> dict:
    raw = read_regular_file(HERE/"data.json", maximum_bytes=128*1024)
    if hashlib.sha256(raw).hexdigest() != DATA_SHA256:
        raise ValueError("scientific data identity mismatch")
    return load_json_object(raw, "patterning data")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def profile(sequence: str, protein: str) -> dict:
    binding = load_data()["binding"]
    reference = binding["consensus"]
    if (type(sequence) is not str or len(sequence) != len(reference) or set(sequence)-set("ACGT")
            or type(protein) is not str or protein not in binding["profiles"]):
        raise ValueError("expected a nine-base unambiguous site and a measured GLI protein")
    changes = [i for i, (a,b) in enumerate(zip(sequence,reference)) if a != b]
    if len(changes) > 1:
        raise ValueError("multiple substitutions require an unvalidated independence assumption")
    result = {"protein":protein, "oriented_sequence":sequence, "reference":reference,
              "changed_position_1based":None, "profile_ratio":1.0}
    if changes:
        i = changes[0]
        weights = binding["profiles"][protein]
        ratio = weights[sequence[i]][i]/weights[reference[i]][i]
        if ratio == 0:
            raise ValueError("a zero assay weight does not establish zero physical affinity")
        result.update(changed_position_1based=i+1, profile_ratio=ratio)
    return result


def selection(fasta: Path, record: str, start: int, strand: int, protein: str) -> dict:
    return {"format":"epigenesis.gli-site-selection", "version":1,
            "source_sha256":hashlib.sha256(read_regular_file(fasta, maximum_bytes=MAX_SOURCE_BYTES)).hexdigest(),
            "record_id":record, "start":start, "strand":strand, "protein":protein, "context":"Nkx2.2"}


def _derive(raw: bytes, contract: dict) -> str:
    fields = {"format","version","source_sha256","record_id","start","strand","protein","context"}
    if (type(contract) is not dict or set(contract) != fields
            or contract["format"] != "epigenesis.gli-site-selection"
            or type(contract["version"]) is not int or contract["version"] != 1
            or type(contract["start"]) is not int or contract["start"] < 0
            or type(contract["strand"]) is not int or contract["strand"] not in (-1,1)
            or type(contract["record_id"]) is not str or not contract["record_id"]
            or contract["context"] != "Nkx2.2"
            or type(contract["source_sha256"]) is not str
            or re.fullmatch("[0-9a-f]{64}", contract["source_sha256"]) is None):
        raise ValueError("invalid site selection contract")
    if len(raw) > MAX_SOURCE_BYTES or hashlib.sha256(raw).hexdigest() != contract["source_sha256"]:
        raise ValueError("selected source identity mismatch")
    native = SequenceCollectionCompiler().compile_fasta_bytes(raw).to_dict()
    if native["inputs"]["wrapper"] is not None:
        raise ValueError("the site consumer requires identity-wrapped FASTA")
    records = [m["sequence_artifact"]["sequence_ir"] for m in native["members"]]
    chosen = [r for r in records if r["record_id"] == contract["record_id"]]
    if len(chosen) != 1:
        raise ValueError("selected FASTA record must be present exactly once")
    count = len(load_data()["binding"]["consensus"])
    sequence = chosen[0]["sequence"][contract["start"]:contract["start"]+count]
    if contract["strand"] == -1:
        sequence = reverse_complement(sequence)
    profile(sequence,contract["protein"])
    return sequence


def semantics(contract: dict) -> str:
    return hashlib.sha256(IMPLEMENTATION+b"\0"+DATA_SHA256.encode()+b"\0"+canonical_bytes(contract)).hexdigest()


def _target(contract: dict, count: int) -> dict:
    return target_artifact({"id":TARGET,"abi_major":1,"opsets":[{"domain":DEV_DOMAIN,"version":1}],
        "unit_schemas":[{"id":"gli-site","fields":[{"id":"nucleotides","type":{"dtype":"u64","shape":[count]},
            "unit":None,"mutability":"constant","numeric":{"kind":"integer"}}]}],
        "edge_schemas":[],"ports":[],"rules":[{"id":RULE,"version":1,"semantics_sha256":semantics(contract),
            "subject":"unit","schema":"gli-site","phase":"UPDATE","triggers":["step"],"parameters":[],
            "reads":[{"scope":"subject","field":"nucleotides"}],"writes":[]}]})


def _operations() -> list[dict]:
    return [{"id":"create.site","op":OP_UNIT_CREATE,"version":1,"schema":"gli-site","count":"count",
             "initializers":[{"field":"nucleotides","tensor":"nucleotides"}]},
            {"id":"attach.profile","op":OP_RULE_ATTACH,"version":1,"rule":RULE,"rule_version":1,
             "subject":"create.site","parameters":[]}]


def compile_site(fasta: Path, contract: dict, output: Path) -> str:
    raw = read_regular_file(fasta, maximum_bytes=MAX_SOURCE_BYTES)
    sequence = _derive(raw,contract)
    contract = load_json_object(canonical_bytes(contract), "site selection")
    receipt = digest(contract)
    output.mkdir()
    (output/"sequence.fasta").write_bytes(raw)
    (output/"selection.json").write_bytes(canonical_bytes(contract)+b"\n")
    source = compile_source(FASTA_PROFILE,{"sequence":output/"sequence.fasta"},parameters={"wrapper":"identity"})
    source.save(output/"source")
    tensors = [{"id":name,"type":{"dtype":"u64","shape":shape},"unit":None,"axes":[None]*len(shape),
                "storage":inline_storage(pack("u64",values))}
               for name,shape,values in (("count",[],[1]),("nucleotides",[len(sequence)],list(sequence.encode())))]
    manifest = seal({"format":"brainc.provider-manifest","version":2,"provider":{"name":TARGET,"version":"1"},
        "model_identity":{"kind":"content-sha256","value":semantics(contract)},
        "accepts":[f"brainc.source-descriptor/v1;profile={FASTA_PROFILE}"],
        "outputs":[{k:t[k] for k in ("id","type","unit","axes")} for t in tensors]})
    save(manifest,output/"manifest.json")
    request = make_development_request(source,output/"manifest.json",[t["id"] for t in tensors])
    save(request,output/"request.json")
    save(seal({"format":"brainc.prediction-response","version":2,"request_artifact_sha256":request["artifact_sha256"],
               "provider":manifest["provider"],"model_identity":manifest["model_identity"],"outputs":tensors}),output/"response.json")
    target = _target(contract,len(sequence))
    save(target,output/"target.json")
    save(policy_artifact(TARGET,{"id":TARGET,"abi_major":1,"contract_sha256":target["contract_sha256"]},
         [{"id":t["id"],"from_output":t["id"]} for t in tensors],_operations()),output/"policy.json")
    compile_development(source,output/"manifest.json",output/"request.json",output/"response.json",
                        output/"policy.json",output/"target.json").save(output/"development")
    consume(output,receipt)
    return receipt


def consume(directory: Path, expected_contract_sha256: str) -> dict:
    contract = load_json_object(read_regular_file(directory/"selection.json",maximum_bytes=128*1024),"site selection")
    if digest(contract) != expected_contract_sha256:
        raise ValueError("selection differs from the caller's expected receipt")
    sequence = read_regular_file(directory/"sequence.fasta",maximum_bytes=MAX_SOURCE_BYTES)
    expected = _derive(sequence,contract)
    descriptor,native = load_source_bundle_directory(directory/"source")
    bundle,artifacts = load_development_bundle_directory(directory/"development")
    validation = validate_development_bundle(bundle,artifacts,descriptor,native,source_inputs={"sequence":sequence})
    if not validation["valid"]:
        raise ValueError(f"site artifact validation failed: {validation}")
    manifest = artifacts["provider_manifest"]
    if (manifest["model_identity"] != {"kind":"content-sha256","value":semantics(contract)}
            or manifest["provider"] != {"name":TARGET,"version":"1"}):
        raise ValueError("unsupported site interpreter identity")
    if artifacts["target_contract"] != _target(contract,len(expected)):
        raise ValueError("unsupported site target or rule semantics")
    module = artifacts["development_module"]["module"]
    if module["entrypoint"]["operations"] != _operations():
        raise ValueError("unsupported site operations")
    tensors = {t["id"]:t for t in module["tensors"]}
    types = {"count":{"dtype":"u64","shape":[]},"nucleotides":{"dtype":"u64","shape":[len(expected)]}}
    if set(tensors) != set(types) or any(tensors[k]["type"] != v for k,v in types.items()):
        raise ValueError("unsupported site tensors")
    if any(t["storage"]["kind"] != "inline-base64" for t in tensors.values()):
        raise ValueError("site consumer requires inline tensors")
    raw = {k:base64.b64decode(t["storage"]["data"],validate=True) for k,t in tensors.items()}
    if raw["count"] != pack("u64",[1]) or raw["nucleotides"] != pack("u64",list(expected.encode())):
        raise ValueError("linked nucleotide tensor disagrees with selected source bases")
    loaded = "".join(chr(value[0]) for value in struct.iter_unpack("<Q",raw["nucleotides"]))
    return {"profile":profile(loaded,contract["protein"]),"contract_sha256":expected_contract_sha256,
            "source_artifact_sha256":descriptor["artifact_sha256"],"development_bundle_sha256":bundle["artifact_sha256"],
            "semantics_sha256":semantics(contract),"data_sha256":DATA_SHA256}


def occupancy(model: dict, state, signal: float, ratio: float = 1, mechanism: str = "competitive") -> tuple:
    target = model["genes"].index("Nkx2.2")
    result = []
    for i,basal in enumerate(model["basal"]):
        q = ratio if i == target else 1
        k = model["gli_binding"][i]
        on = basal*(1+k)*(1+q*model["input"][i]*signal)
        # effective basal weights already absorb reference GliA/GliR competition.
        off = (1+(q if mechanism == "competitive" else 1)*k)*math.prod(
            (1+w*x)**model["sites_per_repressor"] for w,x in zip(model["repression"][i],state) if w)
        result.append(on/(on+off))
    return tuple(result)


def simulate(protocol, *, ratio: float = 1, dt: float = DEFAULT_DT, initial=(0.,0.,0.,0.),
             mechanism: str = "competitive") -> list[tuple]:
    def finite(value):
        return type(value) in (int,float) and math.isfinite(value)

    model = load_data()["model"]
    if (not finite(dt) or dt <= 0 or not finite(ratio) or not 0 < ratio <= 1
            or mechanism not in ("competitive","activation-only")
            or len(initial) != len(model["genes"])
            or any(type(v) not in (int,float) or not 0 <= v <= 1 for v in initial)):
        raise ValueError("invalid timestep, affinity ratio, mechanism or initial state")
    phases = []
    total,steps = 0.,0
    for phase in protocol:
        if (not isinstance(phase,(list,tuple)) or len(phase) != 2
                or any(not finite(v) for v in phase) or phase[0] <= 0 or not 0 <= phase[1] <= 1):
            raise ValueError("protocol needs positive durations and Gli activator fractions in [0,1]")
        duration,signal = phase
        count = duration/dt
        end = total+duration
        if not math.isfinite(count) or not math.isfinite(end) or end <= total or count > MAX_STEPS-steps:
            raise ValueError("protocol exceeds the step budget or representable time")
        count = max(1,math.ceil(count))
        steps += count
        phases.append((total,end,signal,count,duration/count))
        total = end
    if not phases:
        raise ValueError("protocol must not be empty")
    state = tuple(initial)
    trajectory = [(0.,*state)]

    def rhs(y,s):
        if any(not math.isfinite(v) or not 0 <= v <= 1 for v in y):
            raise ArithmeticError("integration stage outside the physical state range")
        return tuple(model["production"]*p-model["decay"]*x
                     for p,x in zip(occupancy(model,y,s,ratio,mechanism),y))

    def step(y,h,signal):
        a = rhs(y,signal)
        b = rhs(tuple(x+h*v/2 for x,v in zip(y,a)),signal)
        c = rhs(tuple(x+h*v/2 for x,v in zip(y,b)),signal)
        d = rhs(tuple(x+h*v for x,v in zip(y,c)),signal)
        result = tuple(x+h*(v+2*w+2*z+u)/6 for x,v,w,z,u in zip(y,a,b,c,d))
        rhs(result,signal)
        return result

    attempts = 0
    for start,end,signal,count,interval in phases:
        elapsed,h = 0.,interval
        for index in range(count):
            stop = (index+1)*interval
            while elapsed < stop:
                h = min(h,stop-elapsed)
                if elapsed+h <= elapsed or attempts >= MAX_STEPS:
                    raise ValueError("adaptive integration exceeds the step budget or representable time")
                attempts += 1
                try:
                    coarse = step(state,h,signal)
                    fine = step(step(state,h/2,signal),h/2,signal)
                    error = max(abs(a-b) for a,b in zip(fine,coarse))/15
                except ArithmeticError:
                    h /= 2
                    continue
                if error <= LOCAL_ERROR_TOLERANCE:
                    state = fine
                    elapsed = stop if h == stop-elapsed else elapsed+h
                # RK4 step doubling estimates the local fifth-order error.
                factor = 2 if error == 0 else min(2,max(.2,.9*(LOCAL_ERROR_TOLERANCE/error)**.2))
                h = min(interval,h*factor)
            trajectory.append((end if index == count-1 else start+stop,*state))
    return trajectory


def named_state(initial: dict) -> tuple:
    genes = load_data()["model"]["genes"]
    if type(initial) is not dict or set(initial) != set(genes):
        raise ValueError("initial_state must specify exactly " + ", ".join(genes))
    return tuple(initial[gene] for gene in genes)


def run_assay(protocol: list, *, ratio: float, mode: str, dt: float = DEFAULT_DT,
              mechanism: str = "competitive", initial: dict | None = None) -> dict:
    if mode not in ("reporter","endogenous"):
        raise ValueError("mode must be reporter or endogenous")
    if type(ratio) not in (int,float) or not math.isfinite(ratio) or not 0 < ratio <= 1:
        raise ValueError("invalid affinity ratio")
    model = load_data()["model"]
    if initial is None:
        initial = dict.fromkeys(model["genes"],0.)
    state = named_state(initial)
    trajectory = simulate(protocol,ratio=1 if mode == "reporter" else ratio,dt=dt,
                          mechanism=mechanism,initial=state)
    phase,boundary = 0,protocol[0][0]
    reporter,signals = [],[]
    for time,*state in trajectory:
        while time > boundary and phase+1 < len(protocol):
            phase += 1
            boundary += protocol[phase][0]
        signal = protocol[phase][1]
        signals.append(signal)
        reporter.append(occupancy(model,state,signal,ratio,mechanism)[model["genes"].index("Nkx2.2")])
    return {"mode":mode,"mechanism":mechanism,"protocol":protocol,"ratio":ratio,
            "initial_state":dict(zip(model["genes"],trajectory[0][1:])),
            "trajectory":trajectory,"signal":signals,"reporter_occupancy":reporter}


def write_assay(path: Path, assay: dict) -> None:
    with path.open("w",newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["time",*load_data()["model"]["genes"],"GliA_fraction","Nkx2.2_promoter_occupancy"])
        for row,signal,reporter in zip(assay["trajectory"],assay["signal"],assay["reporter_occupancy"]):
            writer.writerow([*row,signal,reporter])


def demo(output: Path) -> dict:
    data = load_data()
    output.mkdir()
    compiled = {}
    for name,source in (("reference","reference.fasta"),("variant","variant.fasta"),("allele-rescue","reference.fasta")):
        fasta = HERE/source
        receipt = compile_site(fasta,selection(fasta,"nkx2_2_site",0,-1,"GLI2"),output/name)
        compiled[name] = consume(output/name,receipt)
    for name,row in zip(("reference","variant"),data["observations"]["nkx2_2_reporter"]):
        if compiled[name]["profile"]["oriented_sequence"] != reverse_complement(row["site"]):
            raise ValueError("demo core sequence differs from the published reporter intervention")
    reference_ratio = compiled["reference"]["profile"]["profile_ratio"]
    q = compiled["variant"]["profile"]["profile_ratio"]
    q3 = profile(compiled["variant"]["profile"]["oriented_sequence"],"GLI3")["profile_ratio"]
    protocol = [(20.,.2)]
    cases = {"reference":(reference_ratio,"endogenous","competitive"),
             "variant":(q,"endogenous","competitive"),
             "allele-rescue":(compiled["allele-rescue"]["profile"]["profile_ratio"],"endogenous","competitive"),
             "mediator-rescue":(reference_ratio,"endogenous","competitive"),
             "phenocopy":(q,"endogenous","competitive"),
             "activation-only":(q,"endogenous","activation-only"),
             "GLI3-sensitivity":(q3,"endogenous","competitive"),
             "reporter-reference":(reference_ratio,"reporter","competitive"),
             "reporter-variant":(q,"reporter","competitive")}
    assays = {name:run_assay(protocol,ratio=ratio,mode=mode,mechanism=mechanism)
              for name,(ratio,mode,mechanism) in cases.items()}
    assays["history-naive"] = run_assay([(20,.05)],ratio=reference_ratio,mode="endogenous")
    assays["history-primed"] = run_assay([(5,1),(15,.05)],ratio=reference_ratio,mode="endogenous")
    for name,assay in assays.items():
        write_assay(output/(name+".csv"),assay)
    checks = []

    def check(name,passed,evidence):
        checks.append({"name":name,"passed":bool(passed),"evidence":evidence})

    def difference(a,b):
        return max(abs(x-y) for row,other in zip(a,b) for x,y in zip(row[1:],other[1:]))

    for rescue,reference in (("allele-rescue","reference"),("mediator-rescue","reference"),("phenocopy","variant")):
        error = difference(assays[rescue]["trajectory"],assays[reference]["trajectory"])
        check(rescue,error == 0,{"max_state_difference":error})
    effect = difference(assays["reference"]["trajectory"],assays["variant"]["trajectory"])
    check("sequence changes coupled dynamics",effect > .1,{"max_state_difference":effect})
    check("reporter leaves endogenous trajectory unchanged",
          assays["reporter-reference"]["trajectory"] == assays["reporter-variant"]["trajectory"],
          {"readout":"passive promoter occupancy, not reporter abundance"})
    naive = assays["history-naive"]["trajectory"][-1][3]
    primed = assays["history-primed"]["trajectory"][-1][3]
    check("signal history changes the final Nkx2.2 state",naive < .5 < primed,{"naive":naive,"primed":primed})
    low,high = .01,.6
    model = data["model"]
    state = (.2,.5,.1,.7)
    effects = [occupancy(model,state,s,q)[2]-occupancy(model,state,s,1)[2] for s in (low,high)]
    ablated = occupancy(model,state,low,q,"activation-only")[2]-occupancy(model,state,low,1)[2]
    check("competitive binding reverses the effect across environments",effects[0] > 0 > effects[1],
          {"signals":[low,high],"occupancy_changes":effects})
    check("input-only ablation loses low-signal derepression",ablated < 0,{"occupancy_change":ablated})
    check("neutral environment follows the analytic threshold",
          abs(occupancy(model,state,.1,q)[2]-occupancy(model,state,.1,1)[2]) < 1e-14,
          {"GliA_fraction":model["gli_binding"][2]/model["input"][2]})

    # refine every spatial trajectory, including its transient, before accepting the result.
    positions = [index/40 for index in range(41)]
    times = list(range(101))
    spatial = {"positions":positions,"times":times,"initial_state":dict.fromkeys(model["genes"],0.),
               "reference":[],"variant":[],
               "reporter-reference":[],"reporter-variant":[]}
    errors = {"reference":0.,"variant":0.}
    states = {"reference":[],"variant":[]}
    for name,ratio in (("reference",reference_ratio),("variant",q)):
        for position in positions:
            signal = math.exp(-position/model["gradient_length"])
            trajectory = simulate([(times[-1],signal)],ratio=ratio)
            refined = simulate([(times[-1],signal)],ratio=ratio,dt=DEFAULT_DT/2)
            errors[name] = max(errors[name],difference(trajectory,refined[::2]))
            states[name].append([list(trajectory[round(t/DEFAULT_DT)][1:]) for t in times])
        spatial[name] = [list(frame) for frame in zip(*states[name])]
    for t in range(len(times)):
        reference = spatial["reference"][t]
        for name,ratio in (("reporter-reference",reference_ratio),("reporter-variant",q)):
            spatial[name].append([occupancy(model,y,math.exp(-x/model["gradient_length"]),ratio)[2]
                                  for x,y in zip(positions,reference)])
    check("all spatial transients converge under timestep refinement",max(errors.values()) < 1e-4,
          {"dt":DEFAULT_DT,"refined_dt":DEFAULT_DT/2,"maximum_errors":errors,"tolerance":1e-4})
    final = spatial["reference"][-1]
    domains = {gene:[x for x,y in zip(positions,final) if y[i] > .5]
               for i,gene in enumerate(model["genes"]) if gene != "Pax6"}
    ordered = (all(domains.values()) and max(domains["Nkx2.2"]) < min(domains["Olig2"])
               and max(domains["Olig2"]) < min(domains["Irx3"]))
    check("published baseline forms ordered neural domains",ordered,domains)
    assay_refinement = {}
    for name,assay in assays.items():
        fine = run_assay(assay["protocol"],ratio=assay["ratio"],mode=assay["mode"],
                         mechanism=assay["mechanism"],initial=assay["initial_state"],dt=DEFAULT_DT/2)
        assay_refinement[name] = difference(assay["trajectory"],fine["trajectory"][::2])
    check("intervention trajectories converge under timestep refinement",max(assay_refinement.values()) < 1e-4,
          {"maximum_errors":assay_refinement,"tolerance":1e-4})
    report = {"format":"epigenesis.neural-patterning-experiment","version":1,
              "implementation_sha256":hashlib.sha256(IMPLEMENTATION).hexdigest(),"data_sha256":DATA_SHA256,
              "model_source":model,"binding_source":data["binding"]["source"],
              "assumptions":data["assumptions"],"compiled":compiled,"checks":checks,
              "mechanism_checks_passed":all(c["passed"] for c in checks),
              "biological_validation":"not established; no embryo-level prediction is scored",
              "biological_comparison":data["observations"],"spatial":spatial,
              "numerics":{"method":"adaptive RK4 step doubling","output_interval":DEFAULT_DT,
                  "local_absolute_tolerance":LOCAL_ERROR_TOLERANCE,"switch_readout":"left-limit"},
              "interventions":{name:{"source":"variant" if name in ("variant","mediator-rescue","activation-only","GLI3-sensitivity","reporter-variant") else "reference",
                  "mediator_override":name in ("mediator-rescue","phenocopy"),
                  **{k:v for k,v in assay.items() if k != "trajectory"},
                  "trajectory":assay["trajectory"]} for name,assay in assays.items()}}
    write_report(output,report)
    return report


def write_report(output: Path, report: dict) -> None:
    (output/"report.json").write_text(json.dumps(report,separators=(",",":"),allow_nan=False)+"\n")
    if "spatial" in report:
        template = (HERE/"viewer.html").read_text()
        (output/"index.html").write_text(template.replace("__EXPERIMENT_DATA__",json.dumps(report,allow_nan=False).replace("<","\\u003c")))
    files = {str(path.relative_to(output)):hashlib.sha256(path.read_bytes()).hexdigest()
             for path in sorted(output.rglob("*")) if path.is_file() and path.name != "files.sha256.json"}
    (output/"files.sha256.json").write_text(json.dumps(files,indent=2)+"\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command",required=True)
    demonstration = commands.add_parser("demo",help="run the complete measured-site experiment and controls")
    demonstration.add_argument("--output",type=Path,required=True)
    runner = commands.add_parser("run",help="compile a selected site and run a custom signal protocol")
    runner.add_argument("--fasta",type=Path,required=True)
    runner.add_argument("--record",required=True)
    runner.add_argument("--start",type=int,default=0,help="zero-based start of the nine-base interval")
    runner.add_argument("--strand",type=int,choices=(-1,1),required=True,help="orientation into Hallikas GACCACCCA")
    runner.add_argument("--protein",choices=("GLI2","GLI3"),default="GLI2")
    runner.add_argument("--mode",choices=("reporter","endogenous"),required=True)
    runner.add_argument("--protocol",type=Path,required=True,
                        help='JSON object with "phases" and optional named "initial_state" concentrations')
    runner.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    try:
        if args.command == "demo":
            report = demo(args.output)
            print(json.dumps({"output":str(args.output),"mechanism_checks_passed":report["mechanism_checks_passed"],
                              "checks":report["checks"]},indent=2))
            if not report["mechanism_checks_passed"]:
                parser.exit(1,"patterning: one or more mechanism checks failed; results preserved\n")
        else:
            protocol = load_json_object(read_regular_file(args.protocol,maximum_bytes=128*1024),"protocol")
            if ("phases" not in protocol or type(protocol["phases"]) is not list
                    or set(protocol)-{"phases","initial_state"}):
                raise ValueError('protocol requires a "phases" list and permits an "initial_state" object')
            initial = protocol.get("initial_state",dict.fromkeys(load_data()["model"]["genes"],0.))
            state = named_state(initial)
            contract = selection(args.fasta,args.record,args.start,args.strand,args.protein)
            ratio = profile(_derive(read_regular_file(args.fasta,maximum_bytes=MAX_SOURCE_BYTES),contract),args.protein)["profile_ratio"]
            simulate(protocol["phases"],ratio=ratio,initial=state)
            args.output.mkdir()
            receipt = compile_site(args.fasta,contract,args.output/"compiled")
            compiled = consume(args.output/"compiled",receipt)
            assay = run_assay(protocol["phases"],ratio=compiled["profile"]["profile_ratio"],mode=args.mode,
                              initial=initial)
            write_assay(args.output/"trajectory.csv",assay)
            report = {"compiled":compiled,"assay":assay,"assumptions":load_data()["assumptions"],
                      "biological_validation":"not established"}
            write_report(args.output,report)
            print(json.dumps({"output":str(args.output),"compiled":compiled},indent=2))
    except (ValueError,OSError) as failure:
        parser.exit(2,f"patterning: {failure}\n")


if __name__ == "__main__":
    main()

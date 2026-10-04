"""Source-bound coding predictions under an explicit CDS selection contract."""

from __future__ import annotations

import argparse
import base64
import hashlib
import json
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
IMPLEMENTATION_BYTES = read_regular_file(HERE/"coding.py",maximum_bytes=256*1024)
TARGET = "org.epigenesis.example.coding-consequences"
RULE = TARGET + ".standard-translation"
CODE_SHA256 = "2aecbdd09ecc5e0d233c105090a593737298ef93d88f6bc0c1bef26df6051a04"
MAX_SOURCE_BYTES = 1024*1024
MAX_TEMPLATE_BASES = 10_000


def digest(value: dict) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def genetic_code() -> dict[str, str]:
    raw = read_regular_file(HERE/"data/ncbi-gc.prt", maximum_bytes=32*1024)
    if hashlib.sha256(raw).hexdigest() != CODE_SHA256:
        raise ValueError("NCBI genetic-code authority identity mismatch")
    standard = re.search(rb'\bid 1\s*,(.*?)\n \}', raw, re.S).group(1)
    amino = re.search(rb'ncbieaa\s+"([A-Z*]{64})"', standard).group(1).decode()
    bases = [re.search(rb'-- Base'+str(i).encode()+rb'\s+([TCAG]{64})', standard).group(1).decode() for i in (1,2,3)]
    return dict(zip(map("".join,zip(*bases)),amino))


def translate(template: str) -> dict:
    if type(template) is not str or not 3 <= len(template) <= MAX_TEMPLATE_BASES or set(template)-set("ACGT"):
        raise ValueError("coding template must be bounded unambiguous uppercase DNA")
    if not template.startswith("ATG"):
        raise ValueError("only canonical ATG initiation is supported")
    code = genetic_code()
    peptide = []
    stop_offset = None
    for offset in range(0,len(template)-2,3):
        aa = code[template[offset:offset+3]]
        if aa == "*":
            stop_offset = offset
            break
        peptide.append(aa)
    protein = "".join(peptide)
    return {"peptide":protein, "peptide_sha256":hashlib.sha256(protein.encode()).hexdigest(),
            "template_sha256":hashlib.sha256(template.encode()).hexdigest(), "template_bases":len(template),
            "stop_offset":stop_offset, "stop_codon":None if stop_offset is None else template[stop_offset:stop_offset+3],
            "status":"no-stop-in-template" if stop_offset is None else "terminated",
            "trailing_bases":template[len(template)//3*3:] if stop_offset is None else "",
            "untranslated_bases":0 if stop_offset is None else len(template)-stop_offset-3}


def _check_contract(contract: dict) -> None:
    fields = {"format","version","source_sha256","record_id","start","end","strand","genetic_code","initiation","edit"}
    if type(contract) is dict and type(contract.get("version")) is int and contract["version"] == 2:
        fields = fields-{"start","end"}|{"segments"}
    if type(contract) is not dict or set(contract) != fields:
        raise ValueError("coding contract has missing or unknown fields")
    if contract["format"] != "epigenesis.coding-contract" or type(contract["version"]) is not int or contract["version"] not in (1,2):
        raise ValueError("unsupported coding contract")
    if type(contract["source_sha256"]) is not str or re.fullmatch("[0-9a-f]{64}",contract["source_sha256"]) is None:
        raise ValueError("invalid reference source identity")
    if type(contract["record_id"]) is not str or not contract["record_id"] or len(contract["record_id"]) > 1024:
        raise ValueError("invalid selected record")
    segments = _segments(contract)
    if (type(segments) is not list or not 1 <= len(segments) <= MAX_TEMPLATE_BASES
            or any(type(s) is not dict or set(s) != {"start","end"}
                   or any(type(s[k]) is not int for k in ("start","end"))
                   or not 0 <= s["start"] < s["end"] for s in segments)
            or any(a["end"] > b["start"] for a,b in zip(segments,segments[1:]))
            or not 3 <= sum(s["end"]-s["start"] for s in segments) <= MAX_TEMPLATE_BASES
            or any(type(contract[k]) is not int for k in ("strand","genetic_code"))
            or contract["strand"] not in (-1,1) or contract["genetic_code"] != 1 or contract["initiation"] != "ATG"):
        raise ValueError("unsupported coding coordinates, strand, code or initiation")
    edit = contract["edit"]
    if edit is not None:
        if type(edit) is not dict or set(edit) != {"start","deleted","inserted"}:
            raise ValueError("invalid reference-bound edit fields")
        if (type(edit["start"]) is not int or edit["start"] < 0
                or any(type(edit[k]) is not str or len(edit[k]) > MAX_TEMPLATE_BASES or set(edit[k])-set("ACGT") for k in ("deleted","inserted"))
                or edit["deleted"] == edit["inserted"]):
            raise ValueError("invalid reference-bound edit")


def _segments(contract: dict) -> list[dict]:
    return ([{"start":contract["start"],"end":contract["end"]}]
            if contract["version"] == 1 else contract["segments"])


def _reverse(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT","TGCA"))[::-1]


def _derive(raw: bytes, contract: dict) -> tuple[bytes, str]:
    _check_contract(contract)
    if len(raw) > MAX_SOURCE_BYTES or hashlib.sha256(raw).hexdigest() != contract["source_sha256"]:
        raise ValueError("reference source identity mismatch")
    native = SequenceCollectionCompiler().compile_fasta_bytes(raw).to_dict()
    if native["inputs"]["wrapper"] is not None:
        raise ValueError("only identity-wrapped FASTA is supported by this consumer")
    records = [m["sequence_artifact"]["sequence_ir"] for m in native["members"]]
    selected = next((r for r in records if r["record_id"] == contract["record_id"]),None)
    segments = _segments(contract)
    if selected is None or segments[-1]["end"] > len(selected["sequence"]):
        raise ValueError("selected record or coding interval is absent")
    sequence = selected["sequence"]
    strand = contract["strand"]
    reference = "".join(sequence[s["start"]:s["end"]] for s in segments)
    if set(reference)-set("ACGT"):
        raise ValueError("ambiguous coding template is unsupported")
    if strand == -1:
        reference = _reverse(reference)
    if translate(reference)["stop_offset"] != len(reference)-3:
        raise ValueError("complete reference CDS must end at its first in-frame stop")
    template = reference
    admitted = raw
    edit = contract["edit"]
    if edit is not None:
        lo,hi = edit["start"],edit["start"]+len(edit["deleted"])
        if hi > len(reference) or reference[lo:hi] != edit["deleted"]:
            raise ValueError("edit does not match reference coding bases")
        template = reference[:lo]+edit["inserted"]+reference[hi:]
        oriented = segments if strand == 1 else segments[::-1]
        offset = 0
        owner = None
        for i,segment in enumerate(oriented):
            stop = offset+segment["end"]-segment["start"]
            if not edit["deleted"] and lo == stop and i != len(oriented)-1:
                raise ValueError("insertion at an internal exon junction is ambiguous")
            if offset <= lo and hi <= stop and (lo < stop or i == len(oriented)-1):
                owner = segment
                break
            offset = stop
        if owner is None:
            raise ValueError("a coding edit cannot delete across an exon junction")
        if strand == 1:
            start,end = owner["start"]+lo-offset,owner["start"]+hi-offset
            inserted = edit["inserted"]
        else:
            start,end = owner["end"]-(hi-offset),owner["end"]-(lo-offset)
            inserted = _reverse(edit["inserted"])
        changed = sequence[:start]+inserted+sequence[end:]
        # only the selected record is reconstructed; its remaining bases are assumed reference.
        admitted = (">"+contract["record_id"]+" reconstructed local allele; reference flanks assumed\n"
                    +"\n".join(changed[i:i+70] for i in range(0,len(changed),70))+"\n").encode()
    if not 3 <= len(template) <= MAX_TEMPLATE_BASES or not template.startswith("ATG"):
        raise ValueError("edited template exceeds the size or canonical-initiation contract")
    if len(admitted) > MAX_SOURCE_BYTES:
        raise ValueError("reconstructed source exceeds the input budget")
    return admitted,template


def semantics(contract: dict) -> str:
    return hashlib.sha256(IMPLEMENTATION_BYTES+b"\0"+CODE_SHA256.encode()+b"\0"+canonical_bytes(contract)).hexdigest()


def _target(contract: dict, count: int) -> dict:
    return target_artifact({"id":TARGET, "abi_major":1, "opsets":[{"domain":DEV_DOMAIN,"version":1}],
        "unit_schemas":[{"id":"coding-template","fields":[{"id":"nucleotides","type":{"dtype":"u64","shape":[count]},
            "unit":None,"mutability":"constant","numeric":{"kind":"integer"}}]}],
        "edge_schemas":[],"ports":[],"rules":[{"id":RULE,"version":1,"semantics_sha256":semantics(contract),
            "subject":"unit","schema":"coding-template","phase":"UPDATE","triggers":["step"],"parameters":[],
            "reads":[{"scope":"subject","field":"nucleotides"}],"writes":[]}]})


def _operations() -> list[dict]:
    return [{"id":"create.template","op":OP_UNIT_CREATE,"version":1,"schema":"coding-template","count":"count",
             "initializers":[{"field":"nucleotides","tensor":"nucleotides"}]},
            {"id":"attach.translation","op":OP_RULE_ATTACH,"version":1,"rule":RULE,"rule_version":1,
             "subject":"create.template","parameters":[]}]


def compile_coding(fasta: Path, contract: dict, output: Path) -> str:
    raw = read_regular_file(fasta,maximum_bytes=MAX_SOURCE_BYTES)
    _check_contract(contract)
    contract = load_json_object(canonical_bytes(contract),"coding contract")
    admitted,template = _derive(raw,contract)
    receipt = digest(contract)
    output.mkdir()
    (output/"reference.fasta").write_bytes(raw)
    (output/"sequence.fasta").write_bytes(admitted)
    (output/"contract.json").write_bytes(canonical_bytes(contract)+b"\n")
    source = compile_source(FASTA_PROFILE,{"sequence":output/"sequence.fasta"},parameters={"wrapper":"identity"})
    source.save(output/"source")
    tensors = [{"id":name,"type":{"dtype":"u64","shape":shape},"unit":None,"axes":[None]*len(shape),
                "storage":inline_storage(pack("u64",values))}
               for name,shape,values in (("count",[],[1]),("nucleotides",[len(template)],list(template.encode())))]
    manifest = seal({"format":"brainc.provider-manifest","version":2,"provider":{"name":TARGET,"version":"1"},
        "model_identity":{"kind":"content-sha256","value":semantics(contract)},
        "accepts":[f"brainc.source-descriptor/v1;profile={FASTA_PROFILE}"],
        "outputs":[{k:t[k] for k in ("id","type","unit","axes")} for t in tensors]})
    save(manifest,output/"manifest.json")
    request = make_development_request(source,output/"manifest.json",[t["id"] for t in tensors])
    save(request,output/"request.json")
    save(seal({"format":"brainc.prediction-response","version":2,"request_artifact_sha256":request["artifact_sha256"],
               "provider":manifest["provider"],"model_identity":manifest["model_identity"],"outputs":tensors}),output/"response.json")
    target = _target(contract,len(template))
    save(target,output/"target.json")
    save(policy_artifact(TARGET,{"id":TARGET,"abi_major":1,"contract_sha256":target["contract_sha256"]},
         [{"id":t["id"],"from_output":t["id"]} for t in tensors],_operations()),output/"policy.json")
    compile_development(source,output/"manifest.json",output/"request.json",output/"response.json",
                        output/"policy.json",output/"target.json").save(output/"development")
    report = consume(output,receipt)
    (output/"prediction.json").write_text(json.dumps(report,indent=2)+"\n")
    return receipt


def consume(directory: Path, expected_contract_sha256: str) -> dict:
    contract = load_json_object(read_regular_file(directory/"contract.json",maximum_bytes=128*1024),"coding contract")
    if digest(contract) != expected_contract_sha256:
        raise ValueError("coding contract identity differs from the caller's expected receipt")
    reference = read_regular_file(directory/"reference.fasta",maximum_bytes=MAX_SOURCE_BYTES)
    expected,template = _derive(reference,contract)
    sequence = read_regular_file(directory/"sequence.fasta",maximum_bytes=MAX_SOURCE_BYTES)
    if sequence != expected:
        raise ValueError("admitted sequence does not reconstruct the selected reference and edit")
    descriptor,native = load_source_bundle_directory(directory/"source")
    bundle,artifacts = load_development_bundle_directory(directory/"development")
    validation = validate_development_bundle(bundle,artifacts,descriptor,native,source_inputs={"sequence":sequence})
    if not validation["valid"]:
        raise ValueError(f"coding artifact validation failed: {validation}")
    if artifacts["target_contract"] != _target(contract,len(template)):
        raise ValueError("unsupported coding target, contract or rule semantics")
    module = artifacts["development_module"]["module"]
    if module["entrypoint"]["operations"] != _operations():
        raise ValueError("unsupported coding operations")
    tensors = {t["id"]:t for t in module["tensors"]}
    expected_types = {"count":{"dtype":"u64","shape":[]},"nucleotides":{"dtype":"u64","shape":[len(template)]}}
    if set(tensors) != set(expected_types) or any(tensors[k]["type"] != v for k,v in expected_types.items()):
        raise ValueError("unsupported coding tensor inventory or types")
    if any(t["storage"]["kind"] != "inline-base64" for t in tensors.values()):
        raise ValueError("coding consumer requires inline nucleotide tensors")
    raw = {k:base64.b64decode(t["storage"]["data"],validate=True) for k,t in tensors.items()}
    if raw["count"] != pack("u64",[1]) or raw["nucleotides"] != pack("u64",list(template.encode())):
        raise ValueError("linked nucleotide tensor disagrees with source extraction")
    loaded = "".join(chr(value[0]) for value in struct.iter_unpack("<Q",raw["nucleotides"]))
    return {"format":"epigenesis.coding-prediction","version":1,**translate(loaded),
            "contract_sha256":expected_contract_sha256,"reference_source_sha256":contract["source_sha256"],
            "admitted_source_sha256":hashlib.sha256(sequence).hexdigest(),"source_artifact_sha256":descriptor["artifact_sha256"],
            "development_bundle_sha256":bundle["artifact_sha256"],"semantics_sha256":semantics(contract),
            "genetic_code_sha256":CODE_SHA256,"reconstructed_allele":contract["edit"] is not None,
            "biological_acceptance":False}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command",required=True)
    compiler = commands.add_parser("compile")
    compiler.add_argument("--fasta",type=Path,required=True)
    compiler.add_argument("--contract",type=Path,required=True)
    compiler.add_argument("--output",type=Path,required=True)
    consumer = commands.add_parser("consume")
    consumer.add_argument("directory",type=Path)
    consumer.add_argument("--contract-sha256",required=True)
    args = parser.parse_args()
    try:
        if args.command == "compile":
            contract = load_json_object(read_regular_file(args.contract,maximum_bytes=128*1024),"coding contract")
            receipt = compile_coding(args.fasta,contract,args.output)
            report = consume(args.output,receipt)
        else:
            report = consume(args.directory,args.contract_sha256)
        print(json.dumps(report,indent=2))
    except (ValueError,OSError) as failure:
        parser.exit(2,f"coding: {failure}\n")


if __name__ == "__main__":
    main()

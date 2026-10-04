"""Source-bound FamilyCode single-residue relative DNA preferences."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from brainc._io import load_json_object, read_regular_file


HERE = Path(__file__).resolve().parent
IMPLEMENTATION_SHA256 = hashlib.sha256(read_regular_file(HERE/"binding.py",maximum_bytes=128*1024)).hexdigest()
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
BASES = "ACGT"
FLOOR = 0.01
PINS = {
    "homeodomain.json":"510cd52a527ee1f9b0f871dd0822a1a65a3f26f4daaeaae14baa085282b441ba",
    "BLOSUM62.json":"2d09a35756264ef08df133694c3fd17fd5909a0d64f5d36393868d5079b12cef",
}
spec = importlib.util.spec_from_file_location("binding_coding",HERE.parent/"coding-consequences/coding.py")
coding = importlib.util.module_from_spec(spec)
spec.loader.exec_module(coding)


def _read_data(name: str, directory: Path | None = None) -> dict:
    raw = read_regular_file((HERE/"data" if directory is None else directory)/name,maximum_bytes=1024*1024)
    if hashlib.sha256(raw).hexdigest() != PINS[name]:
        raise ValueError(f"binding model data identity mismatch: {name}")
    return load_json_object(raw,name)


def load_model(directory: Path | None = None) -> list[dict]:
    _read_data("BLOSUM62.json",directory)
    return _read_data("homeodomain.json",directory)["motifs"]


def _predict(rows: list[dict], position: int, residue: str) -> dict:
    if not rows or type(position) is not int or not 1 <= position <= len(rows[0]["protein_alignment"]):
        raise ValueError("conditioning position must be a one-based alignment column")
    if type(residue) is not str or len(residue) != 1 or residue not in AMINO_ACIDS:
        raise ValueError("conditioning residue must be one standard amino acid")
    groups = {}
    for row in rows:
        aa = row["protein_alignment"][position-1]
        groups.setdefault(aa,[]).append(row)
    matched = len(groups.get(residue,[]))
    if matched:
        weights = {residue:1.0}
    else:
        table = _read_data("BLOSUM62.json")
        scores = table["values"][table["axes"][0].index(residue)]
        diagonal = scores[table["axes"][1].index(residue)]
        weights = {aa:1/(diagonal-scores[table["axes"][1].index(aa)]) for aa in groups}
    columns = len(rows[0]["relative_affinity_by_base"]["A"])
    psam = {base:[] for base in BASES}
    centered = {base:[] for base in BASES}
    for column in range(columns):
        contributions = {base:[] for base in BASES}
        for aa,weight in weights.items():
            for row in groups[aa]:
                values = row["relative_affinity_by_base"]
                total = math.fsum(values[base][column] for base in BASES)
                for base in BASES:
                    contributions[base].append(weight*values[base][column]/total/len(groups[aa]))
        # all three PCs with the same single residue reduce to conditional means.
        predicted = {base:math.fsum(contributions[base])/math.fsum(weights.values()) for base in BASES}
        maximum = max(predicted.values())
        values = {base:max(FLOOR,predicted[base]/maximum) for base in BASES}
        log_mean = math.fsum(math.log(values[base]) for base in BASES)/4
        for base in BASES:
            psam[base].append(values[base])
            centered[base].append(math.log(values[base])-log_mean)
    return {"psam":psam,"centered_log_preference":centered,"conditioning_position":position,
            "conditioning_residue":residue,"matched_training_rows":matched,"interpolated":not matched}


def consume(directory: Path, selection: dict, expected_selection_sha256: str) -> dict:
    fields = {"format","version","model_sha256","coding_contract_sha256","reference_alignment",
              "peptide_positions","conditioning_position"}
    if type(selection) is not dict or set(selection) != fields:
        raise ValueError("binding selection has missing or unknown fields")
    selection = dict(selection)
    if type(selection["peptide_positions"]) is list:
        selection["peptide_positions"] = selection["peptide_positions"].copy()
    if (selection["format"] != "epigenesis.binding-selection" or type(selection["version"]) is not int
            or selection["version"] != 1 or selection["model_sha256"] != PINS["homeodomain.json"]):
        raise ValueError("unsupported binding selection or model identity")
    reference = selection["reference_alignment"]
    positions = selection["peptide_positions"]
    condition = selection["conditioning_position"]
    model = load_model()
    if (type(reference) is not str or len(reference) != len(model[0]["protein_alignment"]) or set(reference)-set(AMINO_ACIDS+"-")
            or type(positions) is not list or len(positions) != len(reference)
            or type(condition) is not int or not 1 <= condition <= len(reference)
            or any(p is not None and (type(p) is not int or p < 0) for p in positions)
            or any((p is None) != (aa == "-") for p,aa in zip(positions,reference))):
        raise ValueError("invalid aligned domain or peptide coordinate map")
    mapped = [p for p in positions if p is not None]
    if mapped != sorted(set(mapped)) or positions[condition-1] is None:
        raise ValueError("conditioning column must be mapped and peptide coordinates strictly increase")
    receipt = selection["coding_contract_sha256"]
    if type(receipt) is not str or re.fullmatch("[0-9a-f]{64}",receipt) is None:
        raise ValueError("invalid coding receipt")
    if coding.digest(selection) != expected_selection_sha256:
        raise ValueError("binding selection differs from the caller's expected receipt")
    translated = coding.consume(directory,receipt)
    peptide = translated["peptide"]
    if mapped[-1] >= len(peptide):
        raise ValueError("domain mapping exceeds the translated peptide")
    alignment = "".join("-" if p is None else peptide[p] for p in positions)
    if any(actual != expected and i != condition-1 for i,(actual,expected) in enumerate(zip(alignment,reference))):
        raise ValueError("unmodelled domain residues differ from the selected reference")
    prediction = _predict(model,condition,alignment[condition-1])
    return {"format":"epigenesis.binding-prediction","version":1,"coding":translated,
            "selection_sha256":expected_selection_sha256,"model_sha256":PINS["homeodomain.json"],
            "substitution_table_sha256":PINS["BLOSUM62.json"],"alignment":alignment,
            "implementation_sha256":IMPLEMENTATION_SHA256,
            "dna_positions":model[0]["DNA_positions"],"prediction":prediction,"biological_acceptance":False}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("coding_directory",type=Path)
    parser.add_argument("--selection",type=Path,required=True)
    parser.add_argument("--selection-sha256",required=True)
    args = parser.parse_args()
    try:
        selection = load_json_object(read_regular_file(args.selection,maximum_bytes=128*1024),"binding selection")
        print(json.dumps(consume(args.coding_directory,selection,args.selection_sha256),indent=2))
    except (ValueError,OSError) as failure:
        parser.exit(2,f"binding: {failure}\n")


if __name__ == "__main__":
    main()

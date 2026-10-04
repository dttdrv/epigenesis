"""Learn and audit sequence effects in the GSE188264 neural-induction reporter assay.

Outputs are differences in asinh-transformed reporter activity at seven assayed
times. They are not endogenous expression, biochemical rates, or cell fates.
"""

from __future__ import annotations

import argparse
import base64
from collections import Counter, defaultdict
import csv
from functools import lru_cache
import gzip
import hashlib
import io
import itertools
import json
import math
from pathlib import Path
import random
import re
import statistics
import struct
import sys
from typing import NamedTuple

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from brainc._io import load_json_object, read_regular_file
from brainc.development_bundle import compile_development
from brainc.source import FASTA_PROFILE, compile_source
from brainc.v2 import inline_storage, make_development_request, pack, policy_artifact, save, target_artifact
from brainc.v2._common import seal
from brainc.v2.target import DEV_DOMAIN, OP_UNIT_CREATE
from brainc.validator_development import load_development_bundle_directory, validate_development_paths


TIMES = [0, 3, 6, 12, 24, 48, 72]
PENALTIES = (0.01, 0.1, 1.0, 10.0)
BASES = 171
TARGET = "org.epigenesis.example.neural-regulation"
MODEL_FORMAT = "epigenesis.neural-regulation.kmer-ridge/v1"
COMPLEMENT = str.maketrans("ACGT", "TGCA")
LOCUS = re.compile(r"chr[^_:]+:[0-9]+-[0-9]+")


class Row(NamedTuple):
    id: str
    locus: str
    group: str
    split: str
    kind: str
    reference: str
    sequence: str
    alpha: tuple[float, ...]
    wt_alphas: tuple[tuple[float, ...], ...]


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def words(sequence: str, k: int) -> dict[str, int]:
    if not sequence or set(sequence) - set("ACGT") or type(k) is not int or not 1 <= k <= len(sequence):
        raise ValueError("word counts require uppercase ACGT and a positive supported word length")
    return dict(Counter(min(word, reverse_complement(word))
                        for i in range(len(sequence)-k+1) for word in (sequence[i:i+k],)))


@lru_cache(maxsize=4)
def vocabulary(k: int) -> tuple[str, ...]:
    if type(k) is not int or k not in (2, 3, 4, 5):
        raise ValueError("supported word lengths are 2, 3, 4, 5")
    return tuple(sorted({min(word, reverse_complement(word))
                         for letters in itertools.product("ACGT", repeat=k) for word in ("".join(letters),)}))


def family_groups(parents: dict[str, str], sequences: dict[str, list[str]]) -> dict[str, tuple[str, str]]:
    roots = {locus: locus for locus in parents}

    def root(locus):
        while roots[locus] != locus:
            roots[locus] = roots[roots[locus]]
            locus = roots[locus]
        return locus

    def merge(a, b):
        a, b = sorted((root(a), root(b)))
        roots[b] = a

    intervals = []
    for locus in parents:
        chromosome, span = locus.split(":")
        start, end = map(int, span.split("-"))
        if start < 0 or end-start != BASES:
            raise ValueError("invalid parent interval")
        intervals.append((chromosome, start, end, locus))
    active = []
    for chromosome, start, end, locus in sorted(intervals):
        active = [entry for entry in active if entry[0] == chromosome and entry[2] > start]
        for entry in active:
            merge(locus, entry[3])
        active.append((chromosome, start, end, locus))
    seen = {}
    for locus in sorted(sequences):
        for sequence in sequences[locus]:
            key = min(sequence, reverse_complement(sequence))
            if key in seen:
                merge(locus, seen[key])
            seen[key] = locus
    members = defaultdict(list)
    for locus in sorted(parents):
        members[root(locus)].append(locus)
    result = {}
    for loci in members.values():
        identity = hashlib.sha256(("epigenesis-neural-regulation-v1\n" + "\n".join(loci)).encode()).hexdigest()
        bucket = int(identity, 16) % 10
        split = "train" if bucket < 6 else "validation" if bucket < 8 else "test"
        for locus in loci:
            result[locus] = identity, split
    return result


def load_data(directory: Path | None = None) -> tuple[list[Row], dict]:
    directory = Path(directory) if directory is not None else Path(__file__).with_name("data")
    provenance = load_json_object(read_regular_file(directory / "provenance.json", maximum_bytes=65536), "provenance")
    sources = {}
    for name in ("design.fasta", "alphas.csv.gz"):
        data = read_regular_file(directory / name, maximum_bytes=8*1024*1024)
        pin = provenance["files"][name]
        if len(data) != pin["bytes"] or hashlib.sha256(data).hexdigest() != pin["sha256"]:
            raise ValueError(f"source identity mismatch: {name}")
        sources[name] = data
    if provenance["insert_bases"] != BASES or provenance["times_hours"] != TIMES:
        raise ValueError("unsupported assay coordinates or timepoints")
    designs = defaultdict(set)
    physical = 0
    exact = 0
    header = None
    pieces = []

    def finish():
        nonlocal physical, exact
        if header is None:
            return
        sequence = "".join(pieces)
        prefix, suffix = provenance["prefix"], provenance["suffix"]
        if not prefix or not suffix or not sequence.startswith(prefix) or not sequence.endswith(suffix):
            raise ValueError(f"assay flanks disagree: {header}")
        sequence = sequence[len(prefix):-len(suffix)]
        if len(sequence) != BASES or set(sequence) - set("ACGT"):
            raise ValueError(f"invalid insert: {header}")
        physical += 1
        exact += sequence in designs[header]
        designs[header].add(sequence)

    for line in sources["design.fasta"].decode("ascii").splitlines():
        if line.startswith(">"):
            finish()
            header, pieces = line[1:], []
            if not header or len(header) > 512:
                raise ValueError("invalid design identifier")
        elif header is None:
            raise ValueError("sequence before FASTA header")
        else:
            pieces.append(line)
    finish()
    ambiguous = sorted(name for name, sequences in designs.items() if len(sequences) != 1)
    unique = {name: next(iter(sequences)) for name, sequences in designs.items() if len(sequences) == 1}
    with gzip.GzipFile(fileobj=io.BytesIO(sources["alphas.csv.gz"])) as stream:
        rates = stream.read(8*1024*1024 + 1)
    if len(rates) > 8*1024*1024:
        raise ValueError("activity table exceeds decoded byte budget")
    reader = csv.reader(io.StringIO(rates.decode("ascii")))
    if next(reader) != [""] + [f"{time}h" for time in TIMES]:
        raise ValueError("activity columns disagree with assay times")
    alphas = {}
    for record in reader:
        if len(record) != 8 or record[0] not in designs or record[0] in alphas:
            raise ValueError("unmatched, duplicate, or malformed activity row")
        values = tuple(float(value) for value in record[1:])
        if any(not math.isfinite(value) or value < 0 for value in values):
            raise ValueError("activity must be finite and nonnegative")
        alphas[record[0]] = values
    regions, parents, wt_values = {}, {}, defaultdict(list)
    for name, sequence in unique.items():
        if name.startswith("WT_region"):
            match = re.fullmatch(r"WT_region(\d+)_(chr[^_:]+:\d+-\d+)", name)
            if match is None:
                raise ValueError("malformed parent name")
            number, locus = match.groups()
            if locus in parents and parents[locus] != sequence:
                raise ValueError("same locus has different oriented parent sequences")
            regions[number], parents[locus] = locus, sequence
            if name in alphas:
                wt_values[locus].append(alphas[name])

    def location(name):
        match = LOCUS.search(name)
        if match:
            locus = match.group()
        else:
            control = re.match(r"(?:scram|nonmotif)_seq(\d+)(?:_|$)", name)
            if control is None or control.group(1) not in regions:
                raise ValueError(f"unmatched control parent: {name}")
            locus = regions[control.group(1)]
        if locus not in parents:
            raise ValueError(f"unmatched parent locus: {name}")
        return locus

    by_locus = defaultdict(list)
    for name, sequence in unique.items():
        by_locus[location(name)].append(sequence)
    groups = family_groups(parents, by_locus)
    rows, missing = [], []
    for name, alpha in sorted(alphas.items()):
        if name in ambiguous:
            continue
        locus = location(name)
        if not wt_values[locus]:
            missing.append(name)
            continue
        if name.startswith("WT_region"):
            continue
        kind = "random" if name.startswith("nonmotif_") else "scrambled" if name.startswith("scram_") else "targeted"
        group, split = groups[locus]
        rows.append(Row(name, locus, group, split, kind, parents[locus], unique[name], alpha, tuple(wt_values[locus])))
    qc = {"physical_design_records": physical, "exact_duplicate_records": exact,
          "ambiguous_ids": ambiguous, "activity_rows": len(alphas), "missing_parent_rows": len(missing),
          "missing_parent_ids": missing, "zero_estimates": sum(x == 0 for a in alphas.values() for x in a),
          "eligible_variant_rows": len(rows), "source_provenance": provenance,
          "split_groups": {split: len({row.group for row in rows if row.split == split}) for split in ("train", "validation", "test")}}
    return rows, qc


def transform(value: float, scale: float) -> float:
    if not math.isfinite(value) or value < 0 or not math.isfinite(scale) or scale <= 0:
        raise ValueError("transformation needs a finite nonnegative estimate and positive scale")
    ratio = value / scale
    if math.isinf(ratio):
        # asinh(x) = log(2x) to floating-point precision at this magnitude.
        return (math.log(value) - math.log(scale)) / math.log(2) + 1
    return math.asinh(ratio) / math.log(2)


def training_scale(rows: list[Row]) -> float:
    parents = {row.locus: row.wt_alphas for row in rows}
    values = [value for rates in parents.values() for curve in rates for value in curve if value > 0]
    if not values:
        raise ValueError("no positive training parent estimates")
    return statistics.median(values)


def target(row: Row, scale: float) -> list[float]:
    return [transform(value, scale) - statistics.mean(transform(parent[i], scale) for parent in row.wt_alphas)
            for i, value in enumerate(row.alpha)]


def solve_ridge(gram: list[list[float]], cross: list[list[float]], penalty: float) -> list[list[float]]:
    if not math.isfinite(penalty) or penalty <= 0:
        raise ValueError("ridge penalty must be finite and positive")
    size, outputs = len(gram), len(cross[0])
    lower = [[0.0] * size for _ in gram]
    for i in range(size):
        for j in range(i+1):
            value = gram[i][j] + (penalty if i == j else 0) - sum(lower[i][q]*lower[j][q] for q in range(j))
            if i == j:
                if not math.isfinite(value) or value <= 0:
                    raise ValueError("ridge matrix is not finite positive definite")
                lower[i][j] = math.sqrt(value)
            else:
                lower[i][j] = value / lower[j][j]
    solved = [[0.0] * outputs for _ in gram]
    for h in range(outputs):
        intermediate = [0.0] * size
        for i in range(size):
            intermediate[i] = (cross[i][h] - sum(lower[i][j]*intermediate[j] for j in range(i))) / lower[i][i]
        for i in range(size-1, -1, -1):
            solved[i][h] = (intermediate[i] - sum(lower[j][i]*solved[j][h] for j in range(i+1, size))) / lower[i][i]
    return solved


def moments(rows: list[Row], k: int, scale: float, shuffle_seed: int | None = None):
    index = {word: i for i, word in enumerate(vocabulary(k))}
    counts = Counter(row.group for row in rows)
    if not counts:
        raise ValueError("training rows are empty")
    responses = [target(row, scale) for row in rows]
    if shuffle_seed is not None:
        random.Random(shuffle_seed).shuffle(responses)
    gram = [[0.0] * len(index) for _ in index]
    cross = [[0.0] * len(TIMES) for _ in index]
    for row, response in zip(rows, responses):
        difference = Counter(words(row.sequence, k))
        difference.subtract(words(row.reference, k))
        sparse = sorted((index[word], value) for word, value in difference.items() if value)
        weight = 1 / (len(counts) * counts[row.group])
        for offset, (i, x) in enumerate(sparse):
            weighted = weight * x
            for h, y in enumerate(response):
                cross[i][h] += weighted * y
            for j, y in sparse[:offset+1]:
                gram[i][j] += weighted * y
    for i in range(len(gram)):
        for j in range(i):
            gram[j][i] = gram[i][j]
    return gram, cross


def fit(rows: list[Row], k: int, penalty: float, scale: float, *, sufficient=None, shuffle_seed=None) -> dict:
    gram, cross = moments(rows, k, scale, shuffle_seed) if sufficient is None else sufficient
    return {"format": MODEL_FORMAT, "k": k, "penalty": penalty, "scale": scale, "times_hours": TIMES,
            "coefficients": solve_ridge(gram, cross, len(TIMES) * penalty), "training_rows": len(rows),
            "training_groups": sorted({row.group for row in rows}), "provenance": {}}


def validate_model(model: dict) -> None:
    keys = {"format", "k", "penalty", "scale", "times_hours", "coefficients", "training_rows", "training_groups", "provenance"}
    if set(model) != keys or model["format"] != MODEL_FORMAT or model["times_hours"] != TIMES:
        raise ValueError("unsupported model contract")
    for key in ("scale", "penalty"):
        if type(model[key]) not in (int, float) or not math.isfinite(model[key]) or model[key] <= 0:
            raise ValueError(f"invalid model {key}")
    coefficients = model["coefficients"]
    if not isinstance(coefficients, list) or len(coefficients) != len(vocabulary(model["k"])):
        raise ValueError("coefficient shape mismatch")
    if any(not isinstance(curve, list) or len(curve) != len(TIMES) or
           any(type(x) not in (float, int) or not math.isfinite(x) for x in curve) for curve in coefficients):
        raise ValueError("coefficients must be finite seven-timepoint vectors")


def predict(model: dict, reference: str, sequence: str) -> list[float]:
    validate_model(model)
    if len(reference) != BASES or len(sequence) != BASES:
        raise ValueError("prediction requires two 171-base inserts")
    difference = Counter(words(sequence, model["k"]))
    difference.subtract(words(reference, model["k"]))
    prediction = [math.fsum(difference.get(word, 0) * curve[h] for word, curve in zip(vocabulary(model["k"]), model["coefficients"]))
                  for h in range(len(TIMES))]
    if any(not math.isfinite(value) for value in prediction):
        raise ValueError("sequence-model prediction must be finite")
    return prediction


def family_losses(rows: list[Row], predictions: dict[str, list[float]], scale: float, *, temporal=False) -> dict[str, float]:
    grouped = defaultdict(list)
    for row in rows:
        observed = target(row, scale)
        predicted = predictions[row.id]
        if len(predicted) != len(TIMES) or any(not math.isfinite(x) for x in predicted):
            raise ValueError("prediction must contain seven finite values")
        if temporal:
            observed = [x-statistics.mean(observed) for x in observed]
            predicted = [x-statistics.mean(predicted) for x in predicted]
        grouped[row.group].append(statistics.mean((x-y)**2 for x,y in zip(observed, predicted)))
    if not grouped:
        raise ValueError("evaluation set is empty")
    return {group: statistics.mean(values) for group, values in sorted(grouped.items())}


def improvement_interval(baseline: dict[str, float], candidate: dict[str, float]) -> list[float]:
    if not baseline or baseline.keys() != candidate.keys():
        raise ValueError("paired family scores must have the same nonempty keys")
    differences = [baseline[key]-candidate[key] for key in sorted(baseline)]
    rng = random.Random(20261004)
    bootstraps = sorted(statistics.mean(rng.choices(differences, k=len(differences))) for _ in range(2000))
    return [bootstraps[49], bootstraps[1949]]


def benchmark(output: Path, data: Path | None = None) -> dict:
    output.mkdir()
    rows, qc = load_data(data)
    (output / "data-quality.json").write_text(json.dumps(qc, indent=2)+"\n")
    with (output / "split.csv").open("x", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("id", "locus", "group", "split", "kind"))
        writer.writerows((row.id, row.locus, row.group, row.split, row.kind) for row in rows)
    training = [row for row in rows if row.split == "train"]
    validation = [row for row in rows if row.split == "validation" and row.kind == "targeted"]
    testing = [row for row in rows if row.split == "test"]
    selection_scale = training_scale(training)
    choices = []
    for k in (2, 3, 4, 5):
        sufficient = moments(training, k, selection_scale)
        for penalty in PENALTIES:
            model = fit(training, k, penalty, selection_scale, sufficient=sufficient)
            predictions = {row.id: predict(model, row.reference, row.sequence) for row in validation}
            loss = statistics.mean(family_losses(validation, predictions, selection_scale).values())
            choices.append({"k": k, "penalty": penalty, "validation_mse": loss})
            print(f"validation k={k} penalty={penalty:g} mse={loss:.6g}", flush=True)
    order = lambda choice: (choice["validation_mse"], choice["k"], -choice["penalty"])
    chosen = min((choice for choice in choices if choice["k"] > 2), key=order)
    baseline_choice = min((choice for choice in choices if choice["k"] == 2), key=order)
    fitting = [row for row in rows if row.split != "test"]
    scale = training_scale(fitting)
    model = fit(fitting, chosen["k"], chosen["penalty"], scale)
    model["provenance"] = {"sources": qc["source_provenance"], "selection_scale": selection_scale,
                            "selection": choices, "implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    baseline = fit(fitting, 2, baseline_choice["penalty"], scale)
    (output / "model.json").write_text(json.dumps(model, indent=2, allow_nan=False)+"\n")
    (output / "dinucleotide-model.json").write_text(json.dumps(baseline, indent=2, allow_nan=False)+"\n")
    predicted = {row.id: predict(model, row.reference, row.sequence) for row in testing}
    simple = {row.id: predict(baseline, row.reference, row.sequence) for row in testing}
    zero = {row.id: [0.0]*len(TIMES) for row in testing}
    averaged = {key: [statistics.mean(values)]*len(TIMES) for key,values in predicted.items()}
    primary = [row for row in testing if row.kind == "targeted"]
    scores = family_losses(primary, predicted, scale)
    intervals = {name: improvement_interval(family_losses(primary, other, scale), scores)
                 for name, other in (("zero_effect", zero), ("dinucleotide", simple))}
    temporal = improvement_interval(family_losses(primary, averaged, scale), scores)
    permutations = []
    for seed in (19, 73, 211):
        shuffled = fit(fitting, chosen["k"], chosen["penalty"], scale, shuffle_seed=seed)
        values = {row.id: predict(shuffled, row.reference, row.sequence) for row in primary}
        permutations.append({"seed": seed, "mse": statistics.mean(family_losses(primary, values, scale).values())})
    with (output / "test-predictions.csv").open("x", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("id", "locus", "group", "kind", "hour", "alpha", "wt_alphas", "observed_effect", "predicted_effect", "dinucleotide_effect"))
        for row in testing:
            observed = target(row, scale)
            for i,hour in enumerate(TIMES):
                writer.writerow((row.id,row.locus,row.group,row.kind,hour,row.alpha[i],
                                 json.dumps([curve[i] for curve in row.wt_alphas]),observed[i],predicted[row.id][i],simple[row.id][i]))
    report = {"status": "evaluation-complete", "predictive_acceptance": all(interval[0] > 0 for interval in intervals.values()),
              "biological_development_accepted": False, "scope": "retrospective held-out-family reporter effects under published preprocessing",
              "selection_scale": selection_scale, "final_scale": scale, "selected": chosen, "baseline_selected": baseline_choice,
              "test_groups": len(scores), "test_targeted_rows": len(primary), "test_targeted_mse": statistics.mean(scores.values()),
              "improvement_intervals_95": intervals, "temporal_improvement_interval_95": temporal,
              "temporal_contrast_mse": {name: statistics.mean(family_losses(primary, values, scale, temporal=True).values())
                                        for name,values in (("model",predicted),("dinucleotide",simple),("time_averaged",averaged))},
              "label_permutation_diagnostics": permutations,
              "categories": {kind: {name: statistics.mean(family_losses([row for row in testing if row.kind == kind], values, scale).values())
                                     for name,values in (("model",predicted),("dinucleotide",simple),("zero_effect",zero))}
                             for kind in ("targeted","random","scrambled")},
              "artifacts": {name: hashlib.sha256((output/name).read_bytes()).hexdigest()
                            for name in ("model.json","dinucleotide-model.json","data-quality.json","split.csv","test-predictions.csv")}}
    (output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    return report


def _target() -> dict:
    return target_artifact({"id": TARGET, "abi_major": 1, "opsets": [{"domain": DEV_DOMAIN, "version": 1}],
        "unit_schemas": [{"id": "reporter-effect", "fields": [{"id": "effect", "type": {"dtype": "f64", "shape": [len(TIMES)]},
          "unit": None, "mutability": "constant", "numeric": {"kind": "float"}}]}], "edge_schemas": [], "rules": [], "ports": []})


def compile_prediction(model_path: Path, reference_path: Path, sequence_path: Path, output: Path) -> None:
    model_bytes = read_regular_file(model_path, maximum_bytes=2*1024*1024)
    model = load_json_object(model_bytes, "sequence model")
    reference = read_regular_file(reference_path, maximum_bytes=BASES).decode("ascii")
    sequence = read_regular_file(sequence_path, maximum_bytes=BASES).decode("ascii")
    effect = predict(model, reference, sequence)
    output.mkdir()
    (output / "model.json").write_bytes(model_bytes)
    pair = output / "pair.fasta"
    pair.write_text(f">reference\n{reference}\n>variant\n{sequence}\n")
    source = compile_source(FASTA_PROFILE, {"sequence": pair}, parameters={"wrapper": "identity"})
    source.save(output / "source")
    tensors = [{"id": name, "type": {"dtype": dtype, "shape": shape}, "unit": None,
                "axes": [None]*len(shape), "storage": inline_storage(pack(dtype, values))}
               for name,dtype,shape,values in (("count","u64",[],[1]),("effect","f64",[len(TIMES)],effect))]
    identity = {"kind": "content-sha256", "value": hashlib.sha256(model_bytes).hexdigest()}
    manifest = seal({"format": "brainc.provider-manifest", "version": 2,
                     "provider": {"name": TARGET, "version": "1"}, "model_identity": identity,
                     "accepts": [f"brainc.source-descriptor/v1;profile={FASTA_PROFILE}"],
                     "outputs": [{key:tensor[key] for key in ("id","type","unit","axes")} for tensor in tensors]})
    save(manifest, output / "manifest.json")
    request = make_development_request(source, output / "manifest.json", [tensor["id"] for tensor in tensors])
    save(request, output / "request.json")
    save(seal({"format": "brainc.prediction-response", "version": 2, "request_artifact_sha256": request["artifact_sha256"],
               "provider": manifest["provider"], "model_identity": identity, "outputs": tensors}), output / "response.json")
    target_contract = _target()
    save(target_contract, output / "target.json")
    policy = policy_artifact(TARGET, {"id": TARGET, "abi_major": 1, "contract_sha256": target_contract["contract_sha256"]},
        [{"id": tensor["id"], "from_output": tensor["id"]} for tensor in tensors],
        [{"id": "create.effect", "op": OP_UNIT_CREATE, "version": 1, "schema": "reporter-effect", "count": "count",
          "initializers": [{"field": "effect", "tensor": "effect"}]}])
    save(policy, output / "policy.json")
    compile_development(source, output / "manifest.json", output / "request.json", output / "response.json",
                        output / "policy.json", output / "target.json").save(output / "development")
    load_prediction(output)
    validation = validate_development_paths(output / "source", output / "development", source_inputs={"sequence": pair})
    with (output / "validation.json").open("x") as stream:
        stream.write(json.dumps(validation, indent=2)+"\n")


def load_prediction(directory: Path) -> list[float]:
    validation = validate_development_paths(directory / "source", directory / "development", source_inputs={"sequence": directory / "pair.fasta"})
    if not validation["valid"]:
        raise ValueError(f"prediction artifact validation failed: {validation}")
    _, artifacts = load_development_bundle_directory(directory / "development")
    model_bytes = read_regular_file(directory / "model.json", maximum_bytes=2*1024*1024)
    if artifacts["provider_manifest"]["model_identity"]["value"] != hashlib.sha256(model_bytes).hexdigest():
        raise ValueError("prediction model identity mismatch")
    if artifacts["target_contract"] != _target():
        raise ValueError("unsupported prediction target")
    lines = read_regular_file(directory / "pair.fasta", maximum_bytes=2*BASES+64).decode("ascii").splitlines()
    if len(lines) != 4 or lines[0] != ">reference" or lines[2] != ">variant":
        raise ValueError("unsupported prediction source layout")
    module = artifacts["development_module"]["module"]
    operations = module["entrypoint"]["operations"]
    expected = [{"id": "create.effect", "op": OP_UNIT_CREATE, "version": 1, "schema": "reporter-effect", "count": "count",
                 "initializers": [{"field": "effect", "tensor": "effect"}]}]
    if operations != expected:
        raise ValueError("unsupported prediction operations")
    tensors = {tensor["id"]: base64.b64decode(tensor["storage"]["data"], validate=True) for tensor in module["tensors"]}
    if struct.unpack("<Q", tensors["count"]) != (1,):
        raise ValueError("expected one perturbation effect")
    effect = list(struct.unpack("<7d", tensors["effect"]))
    calculated = predict(load_json_object(model_bytes, "sequence model"), lines[1], lines[3])
    if effect != calculated:
        raise ValueError("compiled effect disagrees with sequence-model prediction")
    return effect


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    trial = commands.add_parser("benchmark")
    trial.add_argument("--output", type=Path, required=True)
    trial.add_argument("--data", type=Path)
    prediction = commands.add_parser("predict")
    for name in ("model", "reference", "variant", "output"):
        prediction.add_argument("--"+name, type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == "benchmark":
            print(json.dumps(benchmark(args.output, args.data), indent=2))
        else:
            compile_prediction(args.model, args.reference, args.variant, args.output)
            print(json.dumps({"times_hours": TIMES, "predicted_effect": load_prediction(args.output)}))
    except (ValueError, OSError, OverflowError) as failure:
        print(f"regulation experiment failed: {failure}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

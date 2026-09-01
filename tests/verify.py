"""Release gates for the compiler-only Epigenesis distribution."""

from __future__ import annotations

import argparse
import ast
import base64
from collections.abc import Iterator
import copy
import gc
import hashlib
import json
import math
import os
from pathlib import Path
from pathlib import PurePosixPath
import shutil
import stat
import subprocess
import sys
import tarfile
import tempfile
import zipfile

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT))

from brainc._canonical import (
    ContractError,
    artifact_digest,
    canonical_bytes,
    digest,
    save_artifact,
)
from brainc._io import MAX_INPUT_BYTES, MAX_JSON_BYTES
from brainc.bio import GFF3Compiler
from brainc.bio_graph import compile_feature_graph
from brainc.compiler import CompilerError, compile_program, load_policy, load_program, save as save_program
from brainc.development_bundle import (
    BUNDLE_FILENAME as DEVELOPMENT_BUNDLE_FILENAME,
    CHILD_FILENAMES as DEVELOPMENT_CHILD_FILENAMES,
    DevelopmentBundleError,
    compile_development,
)
from brainc.insdc import GenBankCompiler, load_genbank_artifact
from brainc.provider import ProviderError, load_response, make_request, save as save_provider
from brainc.sequence_collection import (
    SequenceCollectionCompiler,
    SequenceCollectionError,
    load_sequence_collection,
)
from brainc.source import (
    FASTA_PROFILE,
    RAW_PROFILE,
    REFERENCE_FASTA_PROFILE,
    SourceError,
    compile_source,
    load_source_bundle,
)
from brainc.validator import validate_chain
from brainc.validator_bio import validate_bio_chain
from brainc.validator_bio_graph import validate_feature_graph_bundle
from brainc.v2 import (
    inline_storage,
    make_development_request,
    pack,
    policy_artifact,
    save as save_v2,
    target_artifact,
)
from brainc.v2._common import seal as seal_v2
from brainc.v2.target import DEV_DOMAIN, OP_UNIT_CREATE


U49845_EXPECTED = {
    "raw_bytes": 10541,
    "raw_sha256": "1b41f0096dece0626236d49e0bede43b07ab4fd1aa576c1284cc4fb50bea2fd4",
    "artifact_bytes": 31047,
    "artifact_sha256": "88644e5230b7b568ea60e1c85dd4f0cebf6a020b7c87c849a50c05fb918afd11",
    "bio_ir_sha256": "16d430edfd507bc00ec3c1a6b077a7430713b3ecc1cb0484d9250d97c3f60d8f",
    "sequence_collection_artifact_sha256": "d2c28a3b9fc441c4049875c560aa129a4bc923e45548e7edbc9618b35c3ce5fa",
    "refget_seqcol_digest": "EnzSK21xMemjRN-WyXOuo6mb5DylpTpE",
    "report_sha256": "ec1d09607820fbfa71d47a874b7051bca7d31d5c2381e8694d3ef8acbb1d7f6a",
    "record_ids": ["U49845.1"],
    "sequence_bases": 5028,
    "sequence_sha256": "36203848f0560d3bc23561205438cbb3dc22b068ecfc92f9e3d77334b9ed9e9d",
    "sequence_refget": "SQ.lDmukOs0TZpjopBrnhAwI775FEAhj379",
    "feature_key_counts": {"CDS": 3, "gene": 2, "source": 1},
    "summary": {
        "records": 1,
        "source_bytes": 10541,
        "sequence_bases": 5028,
        "sequence_chunks": 1,
        "features": 6,
        "segments": 6,
        "unresolved_references": [],
    },
}


def _artifact(core: dict) -> dict:
    return {**core, "artifact_sha256": digest(core)}


def _write_json(value: dict, path: Path) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _development_case(
    root: Path,
    source_bundle: object,
    *,
    node_count: int = 1,
) -> dict[str, object]:
    """Build the smallest caller-supplied interpretation over a source bundle."""

    root.mkdir(parents=True)
    source_profile = source_bundle.to_dict()["source_ir"]["profile"]
    source_directory = root / "source"
    source_bundle.save(source_directory)
    source_path = source_directory / "source.json"
    target = target_artifact(
        {
            "id": "org.epigenesis.verify",
            "abi_major": 1,
            "opsets": [{"domain": DEV_DOMAIN, "version": 1}],
            "unit_schemas": [
                {
                    "id": "node",
                    "fields": [
                        {
                            "id": "value",
                            "type": {"dtype": "u64", "shape": []},
                            "unit": None,
                            "mutability": "state",
                            "numeric": {"kind": "integer"},
                        }
                    ],
                }
            ],
            "edge_schemas": [],
            "ports": [],
            "rules": [],
        }
    )
    target_path = root / "target.json"
    save_v2(target, target_path)
    output = {
        "id": "development.node-count",
        "type": {"dtype": "u64", "shape": []},
        "unit": None,
        "axes": [],
        "storage": inline_storage(pack("u64", [node_count])),
    }
    manifest = seal_v2(
        {
            "format": "brainc.provider-manifest",
            "version": 2,
            "provider": {"name": "org.epigenesis.verify", "version": "1"},
            "model_identity": {
                "kind": "content-sha256",
                "value": digest({"algorithm": "caller-supplied-count", "version": 1}),
            },
            "accepts": [
                f"brainc.source-descriptor/v1;profile={source_profile}"
            ],
            "outputs": [
                {key: copy.deepcopy(output[key]) for key in ("id", "type", "unit", "axes")}
            ],
        }
    )
    manifest_path = root / "manifest.json"
    save_v2(manifest, manifest_path)
    request = make_development_request(
        source_bundle,
        manifest_path,
        [output["id"]],
    )
    request_path = root / "request.json"
    save_v2(request, request_path)
    response = seal_v2(
        {
            "format": "brainc.prediction-response",
            "version": 2,
            "request_artifact_sha256": request["artifact_sha256"],
            "provider": copy.deepcopy(manifest["provider"]),
            "model_identity": copy.deepcopy(manifest["model_identity"]),
            "outputs": [output],
        }
    )
    response_path = root / "response.json"
    save_v2(response, response_path)
    policy = policy_artifact(
        "org.epigenesis.verify",
        {
            "id": target["contract"]["id"],
            "abi_major": target["contract"]["abi_major"],
            "contract_sha256": target["contract_sha256"],
        },
        [{"id": output["id"], "from_output": output["id"]}],
        [
            {
                "id": "create.nodes",
                "op": OP_UNIT_CREATE,
                "version": 1,
                "schema": "node",
                "count": output["id"],
                "initializers": [{"field": "value", "tensor": output["id"]}],
            }
        ],
    )
    policy_path = root / "policy.json"
    save_v2(policy, policy_path)
    compiled = compile_development(
        source_bundle,
        manifest_path,
        request_path,
        response_path,
        policy_path,
        target_path,
    )
    return {
        "compiled": compiled,
        "source_directory": source_directory,
        "source_path": source_path,
        "manifest_path": manifest_path,
        "request_path": request_path,
        "response_path": response_path,
        "policy_path": policy_path,
        "target_path": target_path,
    }


def _build_real_chain(root: Path) -> dict[str, Path | dict]:
    source_input = ROOT / "tests" / "data" / "J02482.1.fasta"
    source = root / "source.json"
    SequenceCollectionCompiler().compile_file(source_input).save(source)
    manifest = root / "manifest.json"
    save_artifact(_artifact({
        "format": "brainc.provider-manifest", "version": 1,
        "provider": {"name": "fixture-provider", "version": "1"},
        "model_identity": {"kind": "opaque", "value": "external-test-fixture"},
        "accepts": ["brain01.sequence-collection-ir/v1"],
        "outputs": [{"id": "external.scalar", "type": "number", "unit": None}],
    }), manifest)
    request = root / "request.json"
    request_value = make_request(source, manifest, ["external.scalar"])
    save_provider(request_value, request)
    response = root / "response.json"
    save_artifact(_artifact({
        "format": "brainc.prediction-response", "version": 1,
        "request_artifact_sha256": request_value["artifact_sha256"],
        "provider": {"name": "fixture-provider", "version": "1"},
        "model_identity": {"kind": "opaque", "value": "external-test-fixture"},
        "outputs": [{"id": "external.scalar", "type": "number", "unit": None, "value": 0.375}],
    }), response)
    policy = root / "policy.json"
    save_artifact(_artifact({
        "format": "brainc.lowering-policy", "version": 1, "id": "fixture-lowering",
        "target": {"name": "generic-state-vm", "version": "1"},
        "states": [{
            "id": "fixture.level", "type": "number", "from_output": "external.scalar",
            "transform": {"scale": 2.0, "offset": 0.25,
                          "clamp": {"minimum": 0.0, "maximum": 2.0}, "rounding": None},
        }],
        "links": [], "ports": {"inputs": ["fixture.level"], "outputs": ["fixture.level"]},
    }), policy)
    program = root / "program.json"
    save_program(compile_program(source, manifest, request, response, policy), program)
    report = validate_chain(fasta=source_input, sequence=source, manifest=manifest, request=request,
                            response=response, policy=policy, program=program)
    return {"source_input": source_input, "source": source, "manifest": manifest, "request": request,
            "response": response, "policy": policy, "program": program, "report": report}


def real_dna() -> None:
    source_input = ROOT / "tests" / "data" / "J02482.1.fasta"
    raw = source_input.read_bytes()
    metadata = json.loads((ROOT / "tests" / "data" / "J02482.1.source.json").read_text(encoding="utf-8"))
    assert metadata["accession"] == "J02482.1"
    assert hashlib.sha256(raw).hexdigest() == metadata["fasta_sha256"] == "2826ee08e3506154cdeb7ab734ec6f9ebbea6e5d16d0bec488030192faf492a7"
    bases = b"".join(raw.splitlines()[1:])
    assert len(bases) == metadata["bases"] == 5386
    assert hashlib.sha256(bases).hexdigest() == metadata["normalized_sequence_sha256"] == "97038c7e1edea2297667d7f0426ba942b322c74cb30e072ec66ba47f9c0448d0"
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        chain = _build_real_chain(root / "scalar")
        assert chain["report"]["valid"], chain["report"]
        source = json.loads(Path(chain["source"]).read_text(encoding="utf-8"))
        program = json.loads(Path(chain["program"]).read_text(encoding="utf-8"))
        assert source["members"][0]["sequence_artifact"]["sequence_ir"]["refget_id"] == "SQ.IIXILYBQCpHdC4qpI3sOQ_HAeAm9bmeF"
        assert program["program_ir"]["states"][0]["value"] == 1.0

        gff3_path = ROOT / "tests" / "data" / "J02482.1.gff3"
        gff3 = gff3_path.read_bytes()
        collection = SequenceCollectionCompiler().compile_file(source_input)
        biological_source = root / "biological-source.json"
        collection.save(biological_source)
        bio = GFF3Compiler().compile_file(gff3_path, collection)
        bio_report = validate_bio_chain(
            bio.to_dict(),
            collection.to_dict(),
            fasta_source=raw,
            gff3_source=gff3,
            source_metadata=(ROOT / "tests" / "data" / "J02482.1.source.json").read_bytes(),
        )
        assert bio_report["valid"], bio_report
        assert bio_report["summary"]["sequence_bases"] == 5386
        assert bio_report["summary"]["features"] == 23
        assert bio_report["summary"]["relationships"] == 4

        bundle = compile_feature_graph(
            biological_source, bio, gff3_source=gff3
        )
        graph_report = validate_feature_graph_bundle(
            bundle.to_dict(), collection.to_dict(), bio.to_dict()
        )
        assert graph_report["valid"], graph_report
        assert graph_report["result"]["units"] == 23
        assert graph_report["result"]["edges"] == 4
        bundle.save(root / "feature-graph")
    print("REAL-DNA-E2E-PASSED")


def _walk_key_paths(
    value: object,
    path: tuple[str, ...] = (),
) -> Iterator[tuple[tuple[str, ...], object]]:
    pending = [(path, value)]
    while pending:
        current_path, current = pending.pop()
        if type(current) is dict:
            for key, child in current.items():
                child_path = (*current_path, key)
                yield child_path, child
                pending.append((child_path, child))
        elif type(current) is list:
            pending.extend(((*current_path, str(index)), child) for index, child in enumerate(current))


def _assert_genbank_wire(payload: dict, raw: bytes, expected: dict) -> None:
    collection = payload["sequence_collection"]
    bio_ir = payload["bio_ir"]
    records = bio_ir["records"]
    members = collection["members"]

    assert payload["version"] == 2
    assert payload["artifact_sha256"] == expected["artifact_sha256"]
    assert payload["bio_ir_sha256"] == expected["bio_ir_sha256"]
    assert collection["artifact_sha256"] == expected["sequence_collection_artifact_sha256"]
    assert bio_ir["sequence_collection_artifact_sha256"] == collection["artifact_sha256"]
    assert collection["inputs"]["source"]["sha256"] == expected["raw_sha256"]
    assert collection["inputs"]["source"]["byte_length"] == len(raw)
    assert collection["refget_seqcol"]["digest"] == expected["refget_seqcol_digest"]
    assert [record["record_id"] for record in records] == expected["record_ids"]
    assert [member["record_id"] for member in members] == expected["record_ids"]
    assert sum(member["sequence"]["bases"] for member in members) == expected["sequence_bases"]
    assert [member["sequence"]["refget_id"] for member in members] == [expected["sequence_refget"]]
    assert [member["sequence"]["sha256"] for member in members] == [expected["sequence_sha256"]]
    assert collection["refget_seqcol"]["level_2"] == {
        "lengths": [member["sequence"]["bases"] for member in members],
        "names": expected["record_ids"],
        "sequences": [expected["sequence_refget"]],
    }

    feature_key_counts: dict[str, int] = {}
    for record in records:
        for feature in record["features"]:
            feature_key_counts[feature["key"]] = feature_key_counts.get(feature["key"], 0) + 1
    assert feature_key_counts == expected["feature_key_counts"]

    sequence_paths: list[tuple[str, ...]] = []
    storage_paths: list[tuple[str, ...]] = []
    chunk_paths: list[tuple[str, ...]] = []
    source_chunks_found = False
    raw_sha_count = 0
    sequence_sha_count = 0
    for path, value in _walk_key_paths(payload):
        if path[-1] == "sequence":
            sequence_paths.append(path)
        elif path[-1] == "storage":
            storage_paths.append(path)
        elif path[-1] == "chunks":
            chunk_paths.append(path)
        elif path[-1] == "source_chunks":
            source_chunks_found = True
        if value == expected["raw_sha256"]:
            raw_sha_count += 1
        if value == expected["sequence_sha256"]:
            sequence_sha_count += 1
        if path[-1] == "source_map":
            for nested_path, _ in _walk_key_paths(value, path):
                assert nested_path[-1] not in {"text", "eol", "chunks", "source_chunks"}, nested_path
    expected_sequence_paths = [
        ("sequence_collection", "members", str(index), "sequence")
        for index in range(len(members))
    ]
    expected_storage_paths = [
        ("sequence_collection", "inputs", "source", "storage"),
        *(
            ("sequence_collection", "members", str(index), "sequence", "storage")
            for index in range(len(members))
        ),
    ]
    expected_chunk_paths = [(*path, "chunks") for path in expected_storage_paths]
    assert sorted(sequence_paths) == sorted(expected_sequence_paths)
    assert sorted(storage_paths) == sorted(expected_storage_paths)
    assert sorted(chunk_paths) == sorted(expected_chunk_paths)
    assert not source_chunks_found
    assert "inputs" not in payload
    assert all("sequence" not in record and "source_chunks" not in record for record in records)
    assert raw_sha_count == 1
    assert sequence_sha_count == 1


def _independent_genbank_report(source_path: Path, artifact_path: Path) -> dict:
    validator = ROOT / "brainc" / "validator_insdc.py"
    with tempfile.TemporaryDirectory() as directory:
        result = subprocess.run(
            [
                sys.executable,
                "-I",
                str(validator),
                str(source_path),
                str(artifact_path),
            ],
            check=False,
            capture_output=True,
            text=True,
            cwd=directory,
            timeout=180,
            env=dict(os.environ),
        )
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stderr == "", result.stderr
    return json.loads(result.stdout)


def _run_genbank_gate(source_path: Path, expected: dict) -> dict:
    inspected = source_path.lstat()
    assert stat.S_ISREG(inspected.st_mode), f"GenBank source must be a regular non-linked file: {source_path}"
    assert inspected.st_size == expected["raw_bytes"]
    with source_path.open("rb") as source:
        raw = source.read(expected["raw_bytes"] + 1)
    assert len(raw) == expected["raw_bytes"]
    assert hashlib.sha256(raw).hexdigest() == expected["raw_sha256"]
    if "raw_line_count" in expected:
        assert len(raw.splitlines()) == expected["raw_line_count"]
    if expected.get("final_blank_line") is True:
        assert raw.endswith(b"//\n\n")

    artifact = GenBankCompiler().compile_file(source_path)
    payload = artifact.to_dict()
    _assert_genbank_wire(payload, raw, expected)

    with tempfile.TemporaryDirectory() as directory:
        artifact_path = Path(directory) / "compiled-genbank.json"
        artifact.save(artifact_path)
        assert artifact_path.stat().st_size == expected["artifact_bytes"]

        del payload, artifact
        gc.collect()
        loaded = load_genbank_artifact(artifact_path, genbank_source=raw)
        assert loaded.digest == expected["artifact_sha256"]
        del loaded
        gc.collect()

        report = _independent_genbank_report(source_path, artifact_path)

    assert report["valid"] is True
    assert report["report_sha256"] == expected["report_sha256"]
    assert report["inputs"] == {
        "genbank_source_sha256": expected["raw_sha256"],
        "genbank_artifact_sha256": expected["artifact_sha256"],
        "sequence_collection_artifact_sha256": expected["sequence_collection_artifact_sha256"],
    }
    assert report["replay"] == {
        "source_sha256": expected["raw_sha256"],
        "sequence_collection_artifact_sha256": expected["sequence_collection_artifact_sha256"],
        "refget_seqcol_digest": expected["refget_seqcol_digest"],
        "bio_ir_sha256": expected["bio_ir_sha256"],
        "genbank_artifact_sha256": expected["artifact_sha256"],
    }
    assert report["summary"] == expected["summary"]
    return report


def genbank_real_dna() -> None:
    source_path = ROOT / "tests" / "data" / "U49845.1.gb"
    provenance = json.loads(
        (ROOT / "tests" / "data" / "U49845.1.gb.source.json").read_text(encoding="utf-8")
    )
    assert provenance["record"] == "U49845.1"
    assert provenance["bytes"] == U49845_EXPECTED["raw_bytes"]
    assert provenance["sha256"] == U49845_EXPECTED["raw_sha256"]
    assert provenance["normalized_sequence_sha256"] == U49845_EXPECTED["sequence_sha256"]
    assert provenance["refget_id"] == U49845_EXPECTED["sequence_refget"]
    _run_genbank_gate(source_path, U49845_EXPECTED)
    print("GENBANK-REAL-DNA-PASSED")


def genbank_scale(source_path: Path) -> None:
    metadata = json.loads(
        (ROOT / "tests" / "data" / "CP032762.1.source.json").read_text(encoding="utf-8")
    )
    source = metadata["source"]
    sequence = metadata["sequence"]
    compiled = metadata["compiler_release_gate"]
    expected = {
        **compiled,
        "raw_bytes": source["raw_byte_length"],
        "raw_sha256": source["raw_sha256"],
        "raw_line_count": source["raw_line_count"],
        "final_blank_line": source["final_blank_line"],
        "record_ids": [metadata["version"]],
        "sequence_bases": sequence["length"],
        "sequence_sha256": sequence["sha256"],
        "sequence_refget": sequence["refget"],
        "feature_key_counts": metadata["structure"]["feature_key_counts"],
    }
    _run_genbank_gate(source_path, expected)
    print("GENBANK-SCALE-PASSED")


def runtime_integration() -> None:
    """Materialize the real-DNA module in the separately released runtime."""

    runtime_root = ROOT.parent / "epigenesis-runtime"
    assert runtime_root.is_dir(), "sibling Epigenesis runtime checkout is required"
    source_input = ROOT / "tests" / "data" / "J02482.1.fasta"
    gff3_path = ROOT / "tests" / "data" / "J02482.1.gff3"
    gff3 = gff3_path.read_bytes()
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        collection = SequenceCollectionCompiler().compile_file(source_input)
        source = root / "source.json"
        collection.save(source)
        bio = GFF3Compiler().compile_file(gff3_path, collection)
        bundle = compile_feature_graph(source, bio, gff3_source=gff3)
        paths = bundle.save(root / "feature-graph")

        sys.path.insert(0, str(runtime_root))
        try:
            from epirun.contracts import load_development_module, load_target_contract
            from epirun.development import develop

            target = load_target_contract(paths["target_contract"])
            module = load_development_module(paths["development_module"], target)
            developed = develop(target, module)
        finally:
            if sys.path[0] == str(runtime_root):
                sys.path.pop(0)
        assert developed.final_state_sha256 == (
            "95cd313a51e6c43c29582c6834db0452430468a12975e0429434c329c313c2aa"
        )
        parent_edges = developed.state.edge_set("create.relationship.parent")
        assert parent_edges.source_ids == (6, 9, 12, 13)
        assert parent_edges.target_ids == (22, 21, 21, 22)
    print("RUNTIME-INTEGRATION-PASSED")


def tamper() -> None:
    with tempfile.TemporaryDirectory() as directory:
        chain = _build_real_chain(Path(directory))
        program_path = Path(chain["program"])
        forged = json.loads(program_path.read_text(encoding="utf-8"))
        forged["program_ir"]["states"][0]["value"] = 1.5
        forged["ir_sha256"] = digest(forged["program_ir"])
        forged["artifact_sha256"] = digest({key: value for key, value in forged.items() if key != "artifact_sha256"})
        _write_json(forged, program_path)
        report = validate_chain(fasta=chain["source_input"], sequence=chain["source"], manifest=chain["manifest"],
                                request=chain["request"], response=chain["response"], policy=chain["policy"],
                                program=program_path)
        assert not report["valid"]
        assert "independent lowering replay" in report["checks"][-1]["detail"]

        fresh = _build_real_chain(Path(directory) / "source-map")
        source_path = Path(fresh["source"])
        forged_source = json.loads(source_path.read_text(encoding="utf-8"))
        for segment in forged_source["members"][0]["input_source_map"]["sequence_segments"]:
            segment["line"] += 1
        forged_source["artifact_sha256"] = digest(
            {key: value for key, value in forged_source.items() if key != "artifact_sha256"}
        )
        _write_json(forged_source, source_path)
        try:
            load_sequence_collection(source_path)
        except SequenceCollectionError:
            pass
        else:
            raise AssertionError("producer accepted a shifted collection source map")
        source_report = validate_chain(
            fasta=fresh["source_input"], sequence=source_path, manifest=fresh["manifest"],
            request=fresh["request"], response=fresh["response"], policy=fresh["policy"],
            program=fresh["program"],
        )
        assert not source_report["valid"]
        assert "does not replay from source input" in source_report["checks"][-1]["detail"]
    print("TAMPER-E2E-PASSED")


def source_causality() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        original_path = root / "original.dna"
        copied_path = root / "copied-at-an-unrelated-path.dna"
        changed_path = root / "changed.dna"
        original_path.write_bytes(b"ACGTACGT")
        copied_path.write_bytes(original_path.read_bytes())
        changed_path.write_bytes(b"ACGTACGA")

        original = compile_source(
            RAW_PROFILE,
            {"sequence": original_path},
            parameters={"record_id": "arbitrary-one"},
        )
        copied = compile_source(
            RAW_PROFILE,
            {"sequence": copied_path},
            parameters={"record_id": "arbitrary-one"},
        )
        changed = compile_source(
            RAW_PROFILE,
            {"sequence": changed_path},
            parameters={"record_id": "arbitrary-one"},
        )
        renamed = compile_source(
            RAW_PROFILE,
            {"sequence": copied_path},
            parameters={"record_id": "unrecognized-accession"},
        )

        original_descriptor = original.to_dict()
        changed_descriptor = changed.to_dict()
        renamed_descriptor = renamed.to_dict()
        assert original_descriptor == copied.to_dict(), "filesystem paths affected source output"
        assert original_descriptor["source_ir"]["profile"] == RAW_PROFILE
        assert changed_descriptor["source_ir"]["profile"] == RAW_PROFILE
        assert renamed_descriptor["source_ir"]["profile"] == RAW_PROFILE
        assert (
            original_descriptor["source_ir"]["artifacts"]["sequence"]["format"]
            == changed_descriptor["source_ir"]["artifacts"]["sequence"]["format"]
            == renamed_descriptor["source_ir"]["artifacts"]["sequence"]["format"]
        )
        assert (
            original_descriptor["source_ir"]["inputs"]["sequence"]["sha256"]
            != changed_descriptor["source_ir"]["inputs"]["sequence"]["sha256"]
        )
        assert (
            original_descriptor["source_ir"]["records"][0]["sequence_sha256"]
            != changed_descriptor["source_ir"]["records"][0]["sequence_sha256"]
        )
        assert original_descriptor["source_ir_sha256"] != changed_descriptor["source_ir_sha256"]
        assert original_descriptor["artifact_sha256"] != changed_descriptor["artifact_sha256"]
        assert renamed_descriptor["artifact_sha256"] != original_descriptor["artifact_sha256"]
        assert original.profile == changed.profile == renamed.profile == RAW_PROFILE

        original_case = _development_case(root / "original-case", original)
        copied_case = _development_case(root / "copied-case", copied)
        changed_case = _development_case(root / "changed-case", changed)
        interpreted_case = _development_case(
            root / "interpreted-case",
            original,
            node_count=2,
        )
        original_bundle = original_case["compiled"]
        copied_bundle = copied_case["compiled"]
        changed_bundle = changed_case["compiled"]
        interpreted_bundle = interpreted_case["compiled"]
        assert original_bundle.to_dict() == copied_bundle.to_dict()

        original_module = original_bundle.artifacts["development_module"]
        changed_module = changed_bundle.artifacts["development_module"]
        interpreted_module = interpreted_bundle.artifacts["development_module"]
        assert original_module["module"] == changed_module["module"], (
            "source bytes selected compiler behavior despite identical caller interpretation"
        )
        assert original_module["module_sha256"] == changed_module["module_sha256"]
        assert original_module["sources"]["sequence"] != changed_module["sources"]["sequence"]
        assert original_module["artifact_sha256"] != changed_module["artifact_sha256"]
        assert original_bundle.to_dict()["artifact_sha256"] != changed_bundle.to_dict()["artifact_sha256"]

        assert original_module["sources"]["sequence"] == interpreted_module["sources"]["sequence"]
        assert original_module["module_sha256"] != interpreted_module["module_sha256"]
        assert original_module["module"]["budgets"]["units"] == 1
        assert interpreted_module["module"]["budgets"]["units"] == 2
        assert (
            original_module["module"]["tensors"][0]["storage"]["sha256"]
            != interpreted_module["module"]["tensors"][0]["storage"]["sha256"]
        )
    print("SOURCE-CAUSALITY-PASSED")


def boundary() -> None:
    denied = {"socket", "urllib", "requests", "httpx", "aiohttp", "subprocess", "importlib",
              "torch", "tensorflow", "jax", "sklearn"}
    biological_fixtures = {
        "j02482",
        "nc_001422",
        "phix",
        "u49845",
        "cp032762",
        "5386",
        "5028",
        "5868661",
        "2826ee08e3506154cdeb7ab734ec6f9ebbea6e5d16d0bec488030192faf492a7",
        "97038c7e1edea2297667d7f0426ba942b322c74cb30e072ec66ba47f9c0448d0",
        "sq.iixilybqcphdc4qpi3soq_haeam9bmef",
        "95cd313a51e6c43c29582c6834db0452430468a12975e0429434c329c313c2aa",
        "1b41f0096dece0626236d49e0bede43b07ab4fd1aa576c1284cc4fb50bea2fd4",
        "36203848f0560d3bc23561205438cbb3dc22b068ecfc92f9e3d77334b9ed9e9d",
        "sq.ldmukos0tzpjopbrnhawi775feahj379",
        "c54efd1ee2a811c8efbe12e3f8ea9757013877ce7e5faf8c06d39ac326a91582",
        "b2c317b26275d822aae0063819b0b4012f247404edf87982d1eed303416bde62",
        "sq.sfu4yxxy_mskaqhe06_smfwklbobmb-n",
        "97b13016b899da15a76b66c5d2a94bf839f225e42a3f92fdc20bb63e873768b1",
        "5367765fa8310e0491ec39e8e496750a18b2a51b5f52daaab5ad681b7f929035",
        "447f3c4e32714ae9efc637657ed777f80a8faaceec2302c38c7ba4b9ecee5c0c",
    }
    compiler_root = ROOT / "brainc"
    forbidden_modules = {
        "environment",
        "learner",
        "optimizer",
        "runtime",
        "simulation",
        "trainer",
        "world",
    }
    assert not {path.stem for path in compiler_root.rglob("*.py")} & forbidden_modules
    for path in compiler_root.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        imports: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import): imports.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module: imports.add(node.module.split(".")[0])
        assert not imports & denied, (path, imports & denied)
        text = path.read_text(encoding="utf-8").lower()
        for marker in ("train_model(", "optimizer.step(", "requests.get(", "subprocess.run("):
            assert marker not in text, (path, marker)
        for marker in biological_fixtures:
            assert marker not in text, (path, marker)
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, (str, bytes)):
                literal = node.value.decode("ascii", "ignore") if isinstance(node.value, bytes) else node.value
                normalized = "".join(literal.split()).upper()
                assert not (
                    len(normalized) >= 80 and set(normalized) <= set("ACGTUNRYKMSWBDHV-.?")
                ), (path, "embedded biological sequence literal")

    import tomllib

    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    assert project["dependencies"] == [], "compiler distribution acquired runtime/model dependencies"
    public = set(__import__("brainc").__all__)
    assert not public & {
        "DevelopmentRuntime",
        "Environment",
        "Learner",
        "Simulation",
        "Trainer",
        "World",
    }
    print("COMPILER-BOUNDARY-PASSED")


def resource_path_safety() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        source_path = root / "input.dna"
        source_path.write_bytes(b"ACGTACGT")
        source = compile_source(
            RAW_PROFILE,
            {"sequence": source_path},
            parameters={"record_id": "resource-case"},
        )

        oversized = root / "oversized.dna"
        with oversized.open("xb") as stream:
            stream.truncate(MAX_INPUT_BYTES + 1)
        try:
            compile_source(
                RAW_PROFILE,
                {"sequence": oversized},
                parameters={"record_id": "oversized"},
            )
        except SourceError:
            pass
        else:
            raise AssertionError("source reader accepted an input above its byte ceiling")

        occupied_source = root / "occupied-source"
        occupied_source.mkdir()
        source_sentinel = occupied_source / "keep"
        source_sentinel.write_bytes(b"unchanged")
        try:
            source.save(occupied_source)
        except SourceError:
            pass
        else:
            raise AssertionError("source publisher overwrote an existing directory")
        assert source_sentinel.read_bytes() == b"unchanged"
        assert {item.name for item in occupied_source.iterdir()} == {"keep"}

        unknown_directory = root / "source-with-unknown-child"
        source.save(unknown_directory)
        (unknown_directory / "unknown.json").write_text("{}", encoding="utf-8")
        try:
            load_source_bundle(unknown_directory)
        except SourceError:
            pass
        else:
            raise AssertionError("source loader accepted an unknown bundle child")

        if os.name == "posix":
            linked_input = root / "linked-input.dna"
            linked_input.symlink_to(source_path)
            try:
                compile_source(
                    RAW_PROFILE,
                    {"sequence": linked_input},
                    parameters={"record_id": "linked"},
                )
            except SourceError:
                pass
            else:
                raise AssertionError("source compiler followed an input symlink")

            linked_bundle = root / "linked-source-bundle"
            clean_bundle = root / "clean-source-bundle"
            source.save(clean_bundle)
            linked_bundle.symlink_to(clean_bundle, target_is_directory=True)
            try:
                load_source_bundle(linked_bundle)
            except SourceError:
                pass
            else:
                raise AssertionError("source loader followed a bundle-directory symlink")

            linked_child_bundle = root / "linked-child-bundle"
            source.save(linked_child_bundle)
            child = linked_child_bundle / "sequence.json"
            child.unlink()
            child.symlink_to(clean_bundle / "sequence.json")
            try:
                load_source_bundle(linked_child_bundle)
            except SourceError:
                pass
            else:
                raise AssertionError("source loader followed a bundle-child symlink")

            hardlinked_child_bundle = root / "hardlinked-child-bundle"
            source.save(hardlinked_child_bundle)
            hardlinked_child = hardlinked_child_bundle / "sequence.json"
            external_child = root / "external-sequence.json"
            shutil.copyfile(hardlinked_child, external_child)
            hardlinked_child.unlink()
            os.link(external_child, hardlinked_child)
            try:
                load_source_bundle(hardlinked_child_bundle)
            except SourceError:
                pass
            else:
                raise AssertionError("source loader accepted a multiply linked bundle child")

            source_destination = root / "source-destination"
            source_destination.mkdir()
            destination_alias = root / "source-destination-alias"
            destination_alias.symlink_to(source_destination, target_is_directory=True)
            try:
                source.save(destination_alias)
            except SourceError:
                pass
            else:
                raise AssertionError("source publisher accepted a linked destination")

        case = _development_case(root / "development-case", source)
        compiled = case["compiled"]
        safe_output = root / "development-output"
        compiled.save(safe_output)
        assert {item.name for item in safe_output.iterdir()} == {
            DEVELOPMENT_BUNDLE_FILENAME,
            *DEVELOPMENT_CHILD_FILENAMES.values(),
        }
        published_text = "".join(
            path.read_text(encoding="utf-8") for path in sorted(safe_output.iterdir())
        )
        assert "ACGTACGT" not in published_text
        assert str(root) not in published_text

        occupied_output = root / "occupied-development"
        occupied_output.mkdir()
        development_sentinel = occupied_output / "keep"
        development_sentinel.write_bytes(b"unchanged")
        try:
            compiled.save(occupied_output)
        except DevelopmentBundleError:
            pass
        else:
            raise AssertionError("development publisher overwrote an existing directory")
        assert development_sentinel.read_bytes() == b"unchanged"
        assert {item.name for item in occupied_output.iterdir()} == {"keep"}

        oversized_json = root / "oversized.json"
        with oversized_json.open("xb") as stream:
            stream.truncate(MAX_JSON_BYTES + 1)
        try:
            compile_development(
                source,
                oversized_json,
                case["request_path"],
                case["response_path"],
                case["policy_path"],
                case["target_path"],
            )
        except DevelopmentBundleError:
            pass
        else:
            raise AssertionError("development compiler accepted an oversized chain artifact")

        if os.name == "posix":
            linked_response = root / "response-link.json"
            linked_response.symlink_to(case["response_path"])
            try:
                compile_development(
                    source,
                    case["manifest_path"],
                    case["request_path"],
                    linked_response,
                    case["policy_path"],
                    case["target_path"],
                )
            except DevelopmentBundleError:
                pass
            else:
                raise AssertionError("development compiler followed a chain-artifact symlink")
    print("RESOURCE-PATH-SAFETY-PASSED")


def development_real_dna() -> None:
    """Compile pinned DNA and independently replay it from an isolated wheel install."""

    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        wheels = root / "wheels"
        installed = root / "installed"
        wheels.mkdir()
        installed.mkdir()
        subprocess.run(
            [
                sys.executable,
                "-m",
                "pip",
                "wheel",
                "--no-index",
                "--no-deps",
                "--no-build-isolation",
                "--wheel-dir",
                str(wheels),
                str(ROOT),
            ],
            check=True,
            capture_output=True,
            text=True,
            cwd=root,
        )
        wheel_path, = wheels.glob("*.whl")
        with zipfile.ZipFile(wheel_path) as archive:
            names = set(archive.namelist())
            metadata_name, = [name for name in names if name.endswith(".dist-info/METADATA")]
            entry_points_name, = [
                name for name in names if name.endswith(".dist-info/entry_points.txt")
            ]
            metadata = archive.read(metadata_name).decode("utf-8")
            entry_points = archive.read(entry_points_name).decode("utf-8")
            assert "Version: 1.0.0\n" in metadata
            assert "brainc/source.py" in names
            assert "brainc/development_bundle.py" in names
            assert "brainc/validator_development.py" in names
            assert "brainc = brainc.cli:main" in entry_points
            assert (
                "brainc-validate-development = brainc.validator_development:main"
                in entry_points
            )
            assert not any("tests/data" in name for name in names)

        subprocess.run(
            [
                sys.executable,
                "-m",
                "pip",
                "install",
                "--no-index",
                "--no-deps",
                "--target",
                str(installed),
                str(wheel_path),
            ],
            check=True,
            capture_output=True,
            text=True,
            cwd=root,
        )
        environment = dict(os.environ)
        environment["PYTHONPATH"] = str(installed)
        dna = root / "J02482.1.fasta"
        shutil.copyfile(ROOT / "tests" / "data" / "J02482.1.fasta", dna)
        raw = dna.read_bytes()
        assert hashlib.sha256(raw).hexdigest() == (
            "2826ee08e3506154cdeb7ab734ec6f9ebbea6e5d16d0bec488030192faf492a7"
        )

        source_directory = root / "source"
        source_result = subprocess.run(
            [
                sys.executable,
                "-m",
                "brainc",
                "compile-source",
                "--profile",
                REFERENCE_FASTA_PROFILE,
                "--sequence",
                str(dna),
                "--wrapper",
                "identity",
                "--output",
                str(source_directory),
            ],
            check=True,
            capture_output=True,
            text=True,
            cwd=root,
            env=environment,
        )
        assert json.loads(source_result.stdout)["source"] == str(
            source_directory / "source.json"
        )
        assert {item.name for item in source_directory.iterdir()} == {
            "source.json",
            "reference.json",
        }
        descriptor = json.loads(
            (source_directory / "source.json").read_text(encoding="utf-8")
        )
        assert artifact_digest(descriptor) == descriptor["artifact_sha256"]
        assert descriptor["source_ir"]["profile"] == REFERENCE_FASTA_PROFILE
        assert descriptor["source_ir"]["inputs"]["sequence"] == {
            "sha256": "2826ee08e3506154cdeb7ab734ec6f9ebbea6e5d16d0bec488030192faf492a7",
            "byte_length": len(raw),
        }
        assert descriptor["source_ir"]["records"] == [
            {
                "record_id": "J02482.1",
                "bases": 5386,
                "sequence_sha256": "97038c7e1edea2297667d7f0426ba942b322c74cb30e072ec66ba47f9c0448d0",
                "refget_id": "SQ.IIXILYBQCpHdC4qpI3sOQ_HAeAm9bmeF",
            }
        ]

        chain_root = root / "caller-chain"
        chain_root.mkdir()
        chain_script = r'''
from copy import deepcopy
import json
from pathlib import Path
import sys

from brainc._canonical import digest
from brainc.v2 import inline_storage, make_development_request, pack, policy_artifact, save, target_artifact
from brainc.v2._common import seal
from brainc.v2.target import DEV_DOMAIN, OP_UNIT_CREATE

root = Path(sys.argv[1])
source = Path(sys.argv[2])
source_descriptor = json.loads((source / "source.json").read_text(encoding="utf-8"))
source_profile = source_descriptor["source_ir"]["profile"]
source_bases = sum(record["bases"] for record in source_descriptor["source_ir"]["records"])
target = target_artifact({
    "id": "org.epigenesis.installed-real-dna",
    "abi_major": 1,
    "opsets": [{"domain": DEV_DOMAIN, "version": 1}],
    "unit_schemas": [{
        "id": "dna-unit",
        "fields": [{
            "id": "source-bases",
            "type": {"dtype": "u64", "shape": []},
            "unit": "org.ga4gh.base",
            "mutability": "constant",
            "numeric": {"kind": "integer"},
        }],
    }],
    "edge_schemas": [],
    "ports": [],
    "rules": [],
})
target_path = root / "target.json"
save(target, target_path)
source_output = {
    "id": "dna.source-bases",
    "type": {"dtype": "u64", "shape": []},
    "unit": "org.ga4gh.base",
    "axes": [],
    "storage": inline_storage(pack("u64", [source_bases])),
}
count_output = {
    "id": "dna.unit-count",
    "type": {"dtype": "u64", "shape": []},
    "unit": None,
    "axes": [],
    "storage": inline_storage(pack("u64", [source_bases])),
}
outputs = [source_output, count_output]
manifest = seal({
    "format": "brainc.provider-manifest",
    "version": 2,
    "provider": {"name": "org.epigenesis.installed-gate", "version": "1"},
    "model_identity": {
        "kind": "content-sha256",
        "value": digest({"algorithm": "caller-supplied-base-count", "version": 1}),
    },
    "accepts": [f"brainc.source-descriptor/v1;profile={source_profile}"],
    "outputs": [
        {key: deepcopy(output[key]) for key in ("id", "type", "unit", "axes")}
        for output in outputs
    ],
})
manifest_path = root / "manifest.json"
save(manifest, manifest_path)
request = make_development_request(
    source,
    manifest_path,
    [output["id"] for output in outputs],
)
request_path = root / "request.json"
save(request, request_path)
response = seal({
    "format": "brainc.prediction-response",
    "version": 2,
    "request_artifact_sha256": request["artifact_sha256"],
    "provider": deepcopy(manifest["provider"]),
    "model_identity": deepcopy(manifest["model_identity"]),
    "outputs": outputs,
})
response_path = root / "response.json"
save(response, response_path)
policy = policy_artifact(
    "org.epigenesis.installed-real-dna",
    {
        "id": target["contract"]["id"],
        "abi_major": target["contract"]["abi_major"],
        "contract_sha256": target["contract_sha256"],
    },
    [
        {"id": output["id"], "from_output": output["id"]}
        for output in outputs
    ],
    [{
        "id": "create.dna-units",
        "op": OP_UNIT_CREATE,
        "version": 1,
        "schema": "dna-unit",
        "count": count_output["id"],
        "initializers": [
            {"field": "source-bases", "tensor": source_output["id"]}
        ],
    }],
)
save(policy, root / "policy.json")
'''
        chain_result = subprocess.run(
            [
                sys.executable,
                "-c",
                chain_script,
                str(chain_root),
                str(source_directory),
            ],
            check=False,
            capture_output=True,
            text=True,
            cwd=root,
            env=environment,
        )
        assert chain_result.returncode == 0, chain_result.stdout + chain_result.stderr

        compile_arguments = [
            sys.executable,
            "-m",
            "brainc",
            "compile-development",
            str(source_directory),
            "--manifest",
            str(chain_root / "manifest.json"),
            "--request",
            str(chain_root / "request.json"),
            "--response",
            str(chain_root / "response.json"),
            "--policy",
            str(chain_root / "policy.json"),
            "--target",
            str(chain_root / "target.json"),
        ]
        development_one = root / "development-one"
        development_two = root / "development-two"
        for output in (development_one, development_two):
            result = subprocess.run(
                [*compile_arguments, "--output", str(output)],
                check=False,
                capture_output=True,
                text=True,
                cwd=root,
                env=environment,
            )
            assert result.returncode == 0, result.stdout + result.stderr
            assert json.loads(result.stdout)["bundle"] == str(
                output / DEVELOPMENT_BUNDLE_FILENAME
            )
        first_bytes = {
            path.name: path.read_bytes() for path in sorted(development_one.iterdir())
        }
        second_bytes = {
            path.name: path.read_bytes() for path in sorted(development_two.iterdir())
        }
        assert first_bytes == second_bytes
        assert set(first_bytes) == {
            DEVELOPMENT_BUNDLE_FILENAME,
            *DEVELOPMENT_CHILD_FILENAMES.values(),
        }
        assert b"J02482" not in b"".join(first_bytes.values())
        normalized_dna = b"".join(raw.splitlines()[1:])
        assert normalized_dna[:128] not in b"".join(first_bytes.values())

        bundle = json.loads(first_bytes[DEVELOPMENT_BUNDLE_FILENAME])
        artifacts = {
            role: json.loads(first_bytes[filename])
            for role, filename in DEVELOPMENT_CHILD_FILENAMES.items()
        }
        assert artifact_digest(bundle) == bundle["artifact_sha256"]
        assert bundle["source"]["descriptor"]["artifact_sha256"] == descriptor[
            "artifact_sha256"
        ]
        for role, artifact in artifacts.items():
            assert artifact_digest(artifact) == artifact["artifact_sha256"], role
            assert bundle["artifacts"][role]["artifact_sha256"] == artifact[
                "artifact_sha256"
            ]
        module = artifacts["development_module"]
        record = artifacts["compilation_record"]
        assert module["module"]["budgets"] == {
            "operations": 1,
            "tensor_bytes": 16,
            "units": 5386,
            "edges": 0,
            "attachments": 0,
        }
        assert module["sources"]["sequence"] == bundle["source"]["descriptor"]
        assert record["result"]["module_sha256"] == module["module_sha256"]
        assert record["result"]["budgets"] == module["module"]["budgets"]

        scripts = installed / ("Scripts" if os.name == "nt" else "bin")
        validator = scripts / (
            "brainc-validate-development.exe"
            if os.name == "nt"
            else "brainc-validate-development"
        )
        assert validator.is_file(), validator
        report_one = root / "validation-one.json"
        report_two = root / "validation-two.json"
        for report_path in (report_one, report_two):
            subprocess.run(
                [
                    str(validator),
                    str(source_directory),
                    str(development_one),
                    "--source-input",
                    f"sequence={dna}",
                    "-o",
                    str(report_path),
                ],
                check=True,
                capture_output=True,
                text=True,
                cwd=root,
                env=environment,
            )
        assert report_one.read_bytes() == report_two.read_bytes()
        report = json.loads(report_one.read_text(encoding="utf-8"))
        assert report["format"] == "brainc.development-validation-report"
        assert type(report["version"]) is int and report["version"] == 1
        assert report["valid"] is True
        assert digest(
            {key: value for key, value in report.items() if key != "report_sha256"}
        ) == report["report_sha256"]
    print("DEVELOPMENT-REAL-DNA-PASSED")


def wheel() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory); wheels = root / "wheels"; installed = root / "installed"
        wheels.mkdir(); installed.mkdir()
        subprocess.run([sys.executable, "-m", "pip", "wheel", "--no-index", "--no-deps", "--no-build-isolation",
                        "--wheel-dir", str(wheels), str(ROOT)], check=True, capture_output=True, text=True)
        wheel_path, = wheels.glob("*.whl")
        with zipfile.ZipFile(wheel_path) as archive:
            names = archive.namelist()
            metadata_name, = [name for name in names if name.endswith(".dist-info/METADATA")]
            metadata_text = archive.read(metadata_name).decode("utf-8")
            entry_points_name, = [
                name for name in names if name.endswith(".dist-info/entry_points.txt")
            ]
            entry_points_text = archive.read(entry_points_name).decode("utf-8")
            assert "Version: 1.0.0\n" in metadata_text
            assert "License-Expression: Apache-2.0\n" in metadata_text
            assert "brainc/validator_bio_graph.py" in names
            for required_module in (
                "brainc/insdc.py",
                "brainc/insdc_graph.py",
                "brainc/sequence_collection_v2.py",
                "brainc/validator_insdc.py",
                "brainc/source.py",
                "brainc/source_scale.py",
                "brainc/external_profile.py",
                "brainc/external_source.py",
                "brainc/development_bundle.py",
                "brainc/validator_reference.py",
                "brainc/validator_external.py",
                "brainc/validator_development.py",
                "brainc/validator_insdc_graph.py",
                "brainc/standards/genbank-273-insdc-ft-11.4.authority.json",
            ):
                assert required_module in names, required_module
            assert "brainc = brainc.cli:main" in entry_points_text
            assert (
                "brainc-validate-genbank = brainc.validator_insdc:main"
                in entry_points_text
            )
            assert (
                "brainc-validate-insdc-graph = brainc.validator_insdc_graph:main"
                in entry_points_text
            )
            assert (
                "brainc-validate-development = brainc.validator_development:main"
                in entry_points_text
            )
            assert not any(
                PurePosixPath(name).name.casefold() == "cp032762.1.gb"
                for name in names
            ), "large CP032762.1 source bytes must not be bundled in the wheel"
            assert any(name.endswith(".dist-info/licenses/LICENSE") for name in names)
        subprocess.run([sys.executable, "-m", "pip", "install", "--no-index", "--no-deps", "--target", str(installed),
                        str(wheel_path)], check=True, capture_output=True, text=True)
        environment = dict(os.environ); environment["PYTHONPATH"] = str(installed)
        result = subprocess.run([sys.executable, "-m", "brainc", "--help"], check=True,
                                capture_output=True, text=True, cwd=root, env=environment)
        for command in (
            "compile-collection",
            "validate",
            "compile-v2",
            "check-v2",
            "validate-v2",
            "compile-gff3",
            "validate-gff3",
            "compile-feature-graph",
            "validate-feature-graph",
            "compile-genbank",
            "validate-genbank",
            "compile-insdc-graph",
            "compile-source",
            "compile-external-source",
            "compile-development",
        ):
            assert command in result.stdout, (command, result.stdout)
        subprocess.run(
            [
                sys.executable,
                "-c",
                (
                    "import brainc, brainc.v2, brainc.validator_v2, brainc.bio, "
                    "brainc.bio_graph, brainc.validator_bio, "
                    "brainc.validator_bio_graph, brainc.insdc, "
                    "brainc.insdc_graph, brainc.sequence_collection_v2, "
                    "brainc.validator_insdc, brainc.validator_insdc_graph, "
                    "brainc.source_scale, brainc.external_profile, "
                    "brainc.external_source, brainc.development_bundle, "
                    "brainc.validator_reference, brainc.validator_external, "
                    "brainc.validator_development; "
                    "assert brainc.__version__ == '1.0.0'"
                ),
            ],
            check=True,
            capture_output=True,
            text=True,
            cwd=root,
            env=environment,
        )

        genbank_source = ROOT / "tests" / "data" / "U49845.1.gb"
        genbank_artifact = root / "U49845.1.genbank.json"
        subprocess.run(
            [
                sys.executable,
                "-m",
                "brainc",
                "compile-genbank",
                str(genbank_source),
                "-o",
                str(genbank_artifact),
            ],
            check=True,
            capture_output=True,
            text=True,
            cwd=root,
            env=environment,
        )
        compiled_genbank = json.loads(genbank_artifact.read_text(encoding="utf-8"))
        assert compiled_genbank["format"] == "brainc.bio.insdc-genbank-ir"
        assert compiled_genbank["version"] == 2
        assert compiled_genbank["bio_ir"]["records"][0]["record_id"] == "U49845.1"
        assert (
            compiled_genbank["sequence_collection"]["members"][0]["sequence"]["bases"]
            == 5028
        )

        genbank_report_path = root / "U49845.1.validation.json"
        genbank_validation = subprocess.run(
            [
                sys.executable,
                "-m",
                "brainc",
                "validate-genbank",
                str(genbank_source),
                str(genbank_artifact),
                "--report",
                str(genbank_report_path),
            ],
            check=True,
            capture_output=True,
            text=True,
            cwd=root,
            env=environment,
        )
        genbank_report = json.loads(genbank_validation.stdout)
        assert genbank_report["valid"] is True
        assert genbank_report["summary"]["records"] == 1
        assert genbank_report["summary"]["sequence_bases"] == 5028
        assert json.loads(genbank_report_path.read_text(encoding="utf-8")) == genbank_report

        scripts = installed / ("Scripts" if os.name == "nt" else "bin")
        validator_script = scripts / (
            "brainc-validate-genbank.exe" if os.name == "nt" else "brainc-validate-genbank"
        )
        assert validator_script.is_file(), validator_script
        standalone_report_path = root / "U49845.1.standalone-validation.json"
        subprocess.run(
            [
                str(validator_script),
                str(genbank_source),
                str(genbank_artifact),
                "-o",
                str(standalone_report_path),
            ],
            check=True,
            capture_output=True,
            text=True,
            cwd=root,
            env=environment,
        )
        assert (
            json.loads(standalone_report_path.read_text(encoding="utf-8"))
            == genbank_report
        )

        insdc_graph_directory = root / "U49845.1.insdc-graph"
        insdc_graph_compile = subprocess.run(
            [
                sys.executable,
                "-m",
                "brainc",
                "compile-insdc-graph",
                str(genbank_artifact),
                "-o",
                str(insdc_graph_directory),
            ],
            check=True,
            capture_output=True,
            text=True,
            cwd=root,
            env=environment,
        )
        insdc_graph_summary = json.loads(insdc_graph_compile.stdout)
        assert Path(insdc_graph_summary["bundle"]) == (
            insdc_graph_directory / "bundle.json"
        )
        insdc_graph_validator = scripts / (
            "brainc-validate-insdc-graph.exe"
            if os.name == "nt"
            else "brainc-validate-insdc-graph"
        )
        assert insdc_graph_validator.is_file(), insdc_graph_validator
        insdc_graph_report_path = root / "U49845.1.insdc-graph-validation.json"
        subprocess.run(
            [
                str(insdc_graph_validator),
                str(genbank_source),
                str(genbank_artifact),
                str(insdc_graph_directory),
                "-o",
                str(insdc_graph_report_path),
            ],
            check=True,
            capture_output=True,
            text=True,
            cwd=root,
            env=environment,
        )
        insdc_graph_report = json.loads(
            insdc_graph_report_path.read_text(encoding="utf-8")
        )
        assert insdc_graph_report["valid"] is True
        assert insdc_graph_report["result"]["units"] == 6
        assert insdc_graph_report["result"]["edges"] == 0
        assert insdc_graph_report["result"]["tensor_bytes"] == 536

        source_input = ROOT / "tests" / "data" / "J02482.1.fasta"
        source = root / "source.json"
        subprocess.run(
            [sys.executable, "-m", "brainc", "compile-collection", str(source_input), "-o", str(source)],
            check=True, capture_output=True, text=True, cwd=root, env=environment,
        )
        target = target_artifact(
            {
                "id": "org.example.wheel-gate",
                "abi_major": 1,
                "opsets": [{"domain": DEV_DOMAIN, "version": 1}],
                "unit_schemas": [
                    {
                        "id": "node",
                        "fields": [
                            {
                                "id": "value",
                                "type": {"dtype": "u64", "shape": []},
                                "unit": None,
                                "mutability": "constant",
                                "numeric": {"kind": "integer"},
                            }
                        ],
                    }
                ],
                "edge_schemas": [],
                "ports": [],
                "rules": [],
            }
        )
        target_path = root / "target.json"
        save_v2(target, target_path)
        outputs = [
            {
                "id": "t.count",
                "type": {"dtype": "u64", "shape": []},
                "unit": None,
                "axes": [],
                "storage": inline_storage(pack("u64", [1])),
            },
            {
                "id": "t.value",
                "type": {"dtype": "u64", "shape": []},
                "unit": None,
                "axes": [],
                "storage": inline_storage(pack("u64", [7])),
            },
        ]
        manifest = seal_v2(
            {
                "format": "brainc.provider-manifest",
                "version": 2,
                "provider": {"name": "wheel-gate", "version": "1"},
                "model_identity": {"kind": "content-sha256", "value": digest("wheel-gate")},
                "accepts": ["brain01.sequence-collection-ir/v1"],
                "outputs": [
                    {key: output[key] for key in ("id", "type", "unit", "axes")}
                    for output in outputs
                ],
            }
        )
        manifest_path = root / "manifest.json"
        save_v2(manifest, manifest_path)
        request = root / "request.json"
        subprocess.run(
            [
                sys.executable, "-m", "brainc", "make-request-v2", str(source),
                "--manifest", str(manifest_path), "--output-id", "t.count",
                "--output-id", "t.value", "-o", str(request),
            ],
            check=True, capture_output=True, text=True, cwd=root, env=environment,
        )
        request_value = json.loads(request.read_text(encoding="utf-8"))
        response = seal_v2(
            {
                "format": "brainc.prediction-response",
                "version": 2,
                "request_artifact_sha256": request_value["artifact_sha256"],
                "provider": manifest["provider"],
                "model_identity": manifest["model_identity"],
                "outputs": outputs,
            }
        )
        response_path = root / "response.json"
        save_v2(response, response_path)
        target_binding = {
            "id": target["contract"]["id"],
            "abi_major": target["contract"]["abi_major"],
            "contract_sha256": target["contract_sha256"],
        }
        policy = policy_artifact(
            "wheel-gate",
            target_binding,
            [
                {"id": "t.count", "from_output": "t.count"},
                {"id": "t.value", "from_output": "t.value"},
            ],
            [
                {
                    "id": "create.nodes",
                    "op": OP_UNIT_CREATE,
                    "version": 1,
                    "schema": "node",
                    "count": "t.count",
                    "initializers": [{"field": "value", "tensor": "t.value"}],
                }
            ],
        )
        policy_path = root / "policy.json"
        save_v2(policy, policy_path)
        module = root / "module.json"
        subprocess.run(
            [
                sys.executable, "-m", "brainc", "compile-v2", str(source),
                "--manifest", str(manifest_path), "--request", str(request),
                "--response", str(response_path), "--policy", str(policy_path),
                "--target", str(target_path), "-o", str(module),
            ],
            check=True, capture_output=True, text=True, cwd=root, env=environment,
        )
        validation = subprocess.run(
            [
                sys.executable, "-m", "brainc", "validate-v2",
                "--source-input", str(source_input), "--source-artifact", str(source),
                "--manifest", str(manifest_path), "--request", str(request),
                "--response", str(response_path), "--policy", str(policy_path),
                "--target", str(target_path), "--module", str(module),
            ],
            check=True, capture_output=True, text=True, cwd=root, env=environment,
        )
        assert json.loads(validation.stdout)["valid"]

        gff3_source = ROOT / "tests" / "data" / "J02482.1.gff3"
        bio_ir = root / "bio-ir.json"
        subprocess.run(
            [
                sys.executable, "-m", "brainc", "compile-gff3", str(gff3_source),
                "--sequence-collection", str(source), "-o", str(bio_ir),
            ],
            check=True, capture_output=True, text=True, cwd=root, env=environment,
        )
        bio_report_path = root / "bio-validation.json"
        bio_validation = subprocess.run(
            [
                sys.executable, "-m", "brainc", "validate-gff3", str(bio_ir),
                "--sequence-collection", str(source),
                "--sequence-input", str(source_input),
                "--gff3-source", str(gff3_source),
                "--report", str(bio_report_path),
            ],
            check=True, capture_output=True, text=True, cwd=root, env=environment,
        )
        bio_report = json.loads(bio_validation.stdout)
        assert bio_report["valid"]
        assert bio_report["summary"]["features"] == 23
        assert bio_report["summary"]["relationships"] == 4
        assert json.loads(bio_report_path.read_text(encoding="utf-8")) == bio_report

        graph_directory = root / "feature-graph"
        graph_compile = subprocess.run(
            [
                sys.executable, "-m", "brainc", "compile-feature-graph", str(bio_ir),
                "--sequence-collection", str(source),
                "--gff3-source", str(gff3_source),
                "-o", str(graph_directory),
            ],
            check=True, capture_output=True, text=True, cwd=root, env=environment,
        )
        graph_summary = json.loads(graph_compile.stdout)
        bundle_path = Path(graph_summary["bundle"])
        assert bundle_path == graph_directory / "bundle.json"
        graph_report_path = root / "graph-validation.json"
        graph_validation = subprocess.run(
            [
                sys.executable, "-m", "brainc", "validate-feature-graph", str(bundle_path),
                "--sequence-collection", str(source), "--bio-ir", str(bio_ir),
                "--report", str(graph_report_path),
            ],
            check=True, capture_output=True, text=True, cwd=root, env=environment,
        )
        graph_report = json.loads(graph_validation.stdout)
        assert graph_report["valid"]
        assert graph_report["result"]["units"] == 23
        assert graph_report["result"]["edges"] == 4
        assert json.loads(graph_report_path.read_text(encoding="utf-8")) == graph_report
    print("WHEEL-SMOKE-PASSED")


def sdist() -> None:
    """Build the published source archive and test it outside the checkout."""

    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        distributions = root / "dist"
        distributions.mkdir()
        script = (
            "import pathlib, setuptools.build_meta as backend; "
            f"print(backend.build_sdist({str(distributions)!r}))"
        )
        subprocess.run(
            [sys.executable, "-c", script],
            check=True,
            capture_output=True,
            text=True,
            cwd=ROOT,
            env=dict(os.environ),
        )
        archives = list(distributions.glob("*.tar.gz"))
        assert len(archives) == 1, archives
        archive_path = archives[0]
        extracted = root / "extracted"
        extracted.mkdir()
        with tarfile.open(archive_path, "r:gz") as archive:
            names = archive.getnames()
            top_levels = {name.split("/", 1)[0] for name in names if name}
            assert len(top_levels) == 1, top_levels
            package_root = next(iter(top_levels))
            required = {
                "LICENSE",
                "README.md",
                "SPEC.md",
                "docs/PRIOR_ART.md",
                "docs/STANDARDS.md",
                "pyproject.toml",
                "tests/verify.py",
                "tests/data/J02482.1.fasta",
                "tests/data/J02482.1.gff3",
                "tests/data/J02482.1.source.json",
                "tests/data/J02482.1.gff3.source.json",
                "tests/data/U49845.1.gb",
                "brainc/insdc.py",
                "brainc/sequence_collection_v2.py",
                "brainc/validator_insdc.py",
                "brainc/source.py",
                "brainc/source_scale.py",
                "brainc/external_profile.py",
                "brainc/external_source.py",
                "brainc/development_bundle.py",
                "brainc/validator_reference.py",
                "brainc/validator_external.py",
                "brainc/validator_development.py",
                "tests/test_source_scale.py",
                "tests/test_validator_reference.py",
                "tests/test_external_profile.py",
                "tests/test_external_source.py",
                "tests/test_validator_external.py",
                "tests/test_validator_development.py",
                "brainc/standards/genbank-273-insdc-ft-11.4.authority.json",
            }
            missing = sorted(
                relative
                for relative in required
                if f"{package_root}/{relative}" not in names
            )
            assert not missing, f"sdist is missing publication inputs: {missing}"
            assert not any(
                PurePosixPath(name).name.casefold() == "cp032762.1.gb"
                for name in names
            ), "large CP032762.1 source bytes must not be bundled in the sdist"
            for member in archive.getmembers():
                relative = PurePosixPath(member.name)
                assert not relative.is_absolute(), member.name
                assert relative.parts and ".." not in relative.parts, member.name
                assert "\\" not in member.name, member.name
                assert all(":" not in part for part in relative.parts), member.name
                target = extracted.joinpath(*relative.parts)
                assert target.resolve(strict=False).is_relative_to(
                    extracted.resolve()
                ), member.name
                if member.isdir():
                    target.mkdir(parents=True, exist_ok=True)
                    continue
                assert member.isfile(), (
                    "sdist contains a non-regular archive member",
                    member.name,
                    member.type,
                )
                target.parent.mkdir(parents=True, exist_ok=True)
                source = archive.extractfile(member)
                assert source is not None, member.name
                with source, target.open("xb") as output:
                    shutil.copyfileobj(source, output, length=1024 * 1024)
        source_root = extracted / package_root
        environment = dict(os.environ)
        environment.pop("PYTHONPATH", None)
        result = subprocess.run(
            [sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"],
            check=False,
            capture_output=True,
            text=True,
            cwd=source_root,
            env=environment,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        assert "OK" in result.stderr, result.stderr

        installed = root / "installed"
        subprocess.run(
            [
                sys.executable,
                "-m",
                "pip",
                "install",
                "--no-index",
                "--no-deps",
                "--no-build-isolation",
                "--target",
                str(installed),
                str(source_root),
            ],
            check=True,
            capture_output=True,
            text=True,
            cwd=root,
            env=environment,
        )
        installed_environment = dict(environment)
        installed_environment["PYTHONPATH"] = str(installed)
        help_result = subprocess.run(
            [sys.executable, "-m", "brainc", "--help"],
            check=True,
            capture_output=True,
            text=True,
            cwd=root,
            env=installed_environment,
        )
        for command in (
            "compile-genbank",
            "validate-genbank",
            "compile-source",
            "compile-external-source",
            "compile-development",
        ):
            assert command in help_result.stdout, (command, help_result.stdout)
        subprocess.run(
            [
                sys.executable,
                "-c",
                (
                    "import brainc, brainc.insdc, brainc.sequence_collection_v2, "
                    "brainc.validator_insdc, brainc.source_scale, "
                    "brainc.external_profile, brainc.external_source, "
                    "brainc.development_bundle, brainc.validator_reference, "
                    "brainc.validator_external, brainc.validator_development; "
                    "assert brainc.__version__ == '1.0.0'"
                ),
            ],
            check=True,
            capture_output=True,
            text=True,
            cwd=root,
            env=installed_environment,
        )

        genbank_source = source_root / "tests" / "data" / "U49845.1.gb"
        genbank_artifact = root / "sdist-U49845.1.genbank.json"
        subprocess.run(
            [
                sys.executable,
                "-m",
                "brainc",
                "compile-genbank",
                str(genbank_source),
                "-o",
                str(genbank_artifact),
            ],
            check=True,
            capture_output=True,
            text=True,
            cwd=root,
            env=installed_environment,
        )
        validation = subprocess.run(
            [
                sys.executable,
                "-m",
                "brainc",
                "validate-genbank",
                str(genbank_source),
                str(genbank_artifact),
            ],
            check=True,
            capture_output=True,
            text=True,
            cwd=root,
            env=installed_environment,
        )
        report = json.loads(validation.stdout)
        assert report["valid"] is True
        assert report["summary"]["records"] == 1
        assert report["summary"]["sequence_bases"] == 5028

        scripts = installed / ("Scripts" if os.name == "nt" else "bin")
        validator_script = scripts / (
            "brainc-validate-genbank.exe" if os.name == "nt" else "brainc-validate-genbank"
        )
        assert validator_script.is_file(), validator_script
        standalone = subprocess.run(
            [str(validator_script), str(genbank_source), str(genbank_artifact)],
            check=True,
            capture_output=True,
            text=True,
            cwd=root,
            env=installed_environment,
        )
        assert json.loads(standalone.stdout) == report
    print("SDIST-SELF-TEST-PASSED")


def canonical() -> None:
    sample = {
        "numbers": [333333333.33333329, 1e30, 4.50, 2e-3, 1e-27],
        "string": "€$\u000f\nA'B\"\\\"/",
        "literals": [None, True, False],
    }
    expected = '{"literals":[null,true,false],"numbers":[333333333.3333333,1e+30,4.5,0.002,1e-27],"string":"€$\\u000f\\nA\'B\\\"\\\\\\\"/"}'
    assert canonical_bytes(sample).decode("utf-8") == expected
    vectors = {
        0.0: "0", -0.0: "0", 1e-7: "1e-7", 1e-6: "0.000001", 1e20: "100000000000000000000",
        1e21: "1e+21", 4.9406564584124654e-324: "5e-324", 1.7976931348623157e308: "1.7976931348623157e+308",
    }
    for value, rendered in vectors.items(): assert canonical_bytes(value).decode("ascii") == rendered
    for value in (2**53, math.nan, math.inf, "\ud800"):
        try: canonical_bytes(value)
        except ContractError: pass
        else: raise AssertionError(f"non-I-JSON value accepted: {value!r}")
    seqcol_level_2 = {
        "lengths": [248956422, 133797422, 135086622],
        "names": ["chr1", "chr2", "chr3"],
        "sequences": ["SQ.2648ae1bacce4ec4b6cf337dcae37816", "SQ.907112d17fcb73bcab1ed1c72b97ce68",
                      "SQ.1511375dc2dd1b633af8cf439ae90cec"],
    }
    def sha512t24u(value: bytes) -> str:
        return base64.urlsafe_b64encode(hashlib.sha512(value).digest()[:24]).decode("ascii")
    seqcol_level_1 = {name: sha512t24u(canonical_bytes(value)) for name, value in seqcol_level_2.items()}
    assert seqcol_level_1 == {
        "lengths": "IOlarejnLTmdv3-CqehLpcxAR9yNeR1i",
        "names": "g04lKdxiYtG3dOGeUC5AdKEifw65G0Wp",
        "sequences": "ixJdEJlNBgz5U49vfIUqmq3kD4oOtLpd",
    }
    assert sha512t24u(canonical_bytes({name: seqcol_level_1[name] for name in ("names", "sequences")})) == "KxZO6qIbVNCIKtQj0WR3fwzg2rsJLlC3"
    print("CANONICAL-RFC8785-PASSED")


def contract_attacks() -> None:
    with tempfile.TemporaryDirectory() as directory:
        chain = _build_real_chain(Path(directory))
        manifest_path = Path(chain["manifest"]); request_path = Path(chain["request"])
        response_path = Path(chain["response"]); program_path = Path(chain["program"])
        manifest = json.loads(manifest_path.read_text()); request = json.loads(request_path.read_text())
        response = json.loads(response_path.read_text()); program = json.loads(program_path.read_text())
        manifest["accepts"] = ["brain01.sequence-ir/v2"]
        manifest["artifact_sha256"] = digest({key: value for key, value in manifest.items() if key != "artifact_sha256"})
        request["provider_manifest_sha256"] = manifest["artifact_sha256"]
        request["artifact_sha256"] = digest({key: value for key, value in request.items() if key != "artifact_sha256"})
        response["request_artifact_sha256"] = request["artifact_sha256"]
        response["artifact_sha256"] = digest({key: value for key, value in response.items() if key != "artifact_sha256"})
        program["sources"]["provider_manifest"]["artifact_sha256"] = manifest["artifact_sha256"]
        program["sources"]["prediction_request"]["artifact_sha256"] = request["artifact_sha256"]
        program["sources"]["prediction_response"]["artifact_sha256"] = response["artifact_sha256"]
        program["artifact_sha256"] = digest({key: value for key, value in program.items() if key != "artifact_sha256"})
        for value, path in ((manifest, manifest_path), (request, request_path), (response, response_path), (program, program_path)):
            _write_json(value, path)
        try: compile_program(chain["source"], manifest_path, request_path, response_path, chain["policy"])
        except ProviderError: pass
        else: raise AssertionError("compiler accepted undeclared source format")
        report = validate_chain(fasta=chain["source_input"], sequence=chain["source"], manifest=manifest_path,
                                request=request_path, response=response_path, policy=chain["policy"], program=program_path)
        assert not report["valid"]

        output_chain = _build_real_chain(Path(directory) / "output-contract")
        output_manifest_path = Path(output_chain["manifest"]); output_request_path = Path(output_chain["request"])
        output_response_path = Path(output_chain["response"]); output_program_path = Path(output_chain["program"])
        output_manifest = json.loads(output_manifest_path.read_text()); output_request = json.loads(output_request_path.read_text())
        output_response = json.loads(output_response_path.read_text()); output_program = json.loads(output_program_path.read_text())
        output_manifest["outputs"] = [{"id": "other", "type": "number", "unit": None}]
        output_manifest["artifact_sha256"] = digest({key: value for key, value in output_manifest.items() if key != "artifact_sha256"})
        output_request["provider_manifest_sha256"] = output_manifest["artifact_sha256"]
        output_request["artifact_sha256"] = digest({key: value for key, value in output_request.items() if key != "artifact_sha256"})
        output_response["request_artifact_sha256"] = output_request["artifact_sha256"]
        output_response["artifact_sha256"] = digest({key: value for key, value in output_response.items() if key != "artifact_sha256"})
        output_program["sources"]["provider_manifest"]["artifact_sha256"] = output_manifest["artifact_sha256"]
        output_program["sources"]["prediction_request"]["artifact_sha256"] = output_request["artifact_sha256"]
        output_program["sources"]["prediction_response"]["artifact_sha256"] = output_response["artifact_sha256"]
        output_program["artifact_sha256"] = digest({key: value for key, value in output_program.items() if key != "artifact_sha256"})
        for value, path in ((output_manifest, output_manifest_path), (output_request, output_request_path),
                            (output_response, output_response_path), (output_program, output_program_path)):
            _write_json(value, path)
        try:
            compile_program(output_chain["source"], output_manifest_path, output_request_path,
                            output_response_path, output_chain["policy"])
        except ProviderError: pass
        else: raise AssertionError("compiler accepted request output absent from manifest")
        output_report = validate_chain(
            fasta=output_chain["source_input"], sequence=output_chain["source"], manifest=output_manifest_path,
            request=output_request_path, response=output_response_path, policy=output_chain["policy"],
            program=output_program_path,
        )
        assert not output_report["valid"]

        unsafe = {"format": "brainc.prediction-response", "version": 1, "request_artifact_sha256": "0" * 64,
                  "provider": {"name": "p", "version": "1"}, "model_identity": {"kind": "opaque", "value": "m"},
                  "outputs": [{"id": "n", "type": "integer", "unit": None, "value": 2**53}],
                  "artifact_sha256": "0" * 64}
        _write_json(unsafe, response_path)
        try: load_response(response_path)
        except (ContractError, ProviderError): pass
        else: raise AssertionError("unsafe integer accepted")

        integer_root = Path(directory) / "integer-boundary"; integer_root.mkdir()
        integer_source = integer_root / "source.json"
        SequenceCollectionCompiler().compile_file(chain["source_input"]).save(integer_source)
        integer_manifest = integer_root / "manifest.json"
        save_artifact(_artifact({
            "format": "brainc.provider-manifest", "version": 1,
            "provider": {"name": "integer-provider", "version": "1"},
            "model_identity": {"kind": "opaque", "value": "integer-fixture"},
            "accepts": ["brain01.sequence-collection-ir/v1"],
            "outputs": [{"id": "count", "type": "integer", "unit": None}],
        }), integer_manifest)
        integer_request = integer_root / "request.json"
        integer_request_value = make_request(integer_source, integer_manifest, ["count"])
        save_provider(integer_request_value, integer_request)
        integer_response = integer_root / "response.json"
        save_artifact(_artifact({
            "format": "brainc.prediction-response", "version": 1,
            "request_artifact_sha256": integer_request_value["artifact_sha256"],
            "provider": {"name": "integer-provider", "version": "1"},
            "model_identity": {"kind": "opaque", "value": "integer-fixture"},
            "outputs": [{"id": "count", "type": "integer", "unit": None, "value": 2**53 - 1}],
        }), integer_response)
        integer_policy = integer_root / "policy.json"
        integer_policy_value = _artifact({
            "format": "brainc.lowering-policy", "version": 1, "id": "integer-boundary",
            "target": {"name": "generic-state-vm", "version": "1"},
            "states": [{"id": "count", "type": "integer", "from_output": "count",
                        "transform": {"scale": 1, "offset": 0, "clamp": None, "rounding": "nearest"}}],
            "links": [], "ports": {"inputs": ["count"], "outputs": ["count"]},
        })
        save_artifact(integer_policy_value, integer_policy)
        integer_program = integer_root / "program.json"
        save_program(compile_program(integer_source, integer_manifest, integer_request, integer_response, integer_policy), integer_program)
        assert json.loads(integer_program.read_text())["program_ir"]["states"][0]["value"] == 2**53 - 1
        integer_report = validate_chain(
            fasta=chain["source_input"], sequence=integer_source, manifest=integer_manifest,
            request=integer_request, response=integer_response, policy=integer_policy, program=integer_program,
        )
        assert integer_report["valid"], integer_report
        integer_policy_value["states"][0]["transform"]["scale"] = 2
        integer_policy_value["artifact_sha256"] = digest(
            {key: value for key, value in integer_policy_value.items() if key != "artifact_sha256"}
        )
        save_artifact(integer_policy_value, integer_policy)
        try: compile_program(integer_source, integer_manifest, integer_request, integer_response, integer_policy)
        except CompilerError: pass
        else: raise AssertionError("derived unsafe integer accepted")

        original = _build_real_chain(Path(directory) / "fresh")
        program_path = Path(original["program"]); boolean_version = json.loads(program_path.read_text())
        boolean_version["sources"]["sequence"]["version"] = True
        boolean_version["artifact_sha256"] = digest({key: value for key, value in boolean_version.items() if key != "artifact_sha256"})
        _write_json(boolean_version, program_path)
        try: load_program(program_path)
        except CompilerError: pass
        else: raise AssertionError("boolean source version accepted")

        corrupt_gzip = Path(directory) / "corrupt.gz"
        corrupt_gzip.write_bytes(b"\x1f\x8b\x08\x00broken")
        malformed_report = validate_chain(
            fasta=corrupt_gzip, sequence=original["source"], manifest=original["manifest"],
            request=original["request"], response=original["response"], policy=original["policy"],
            program=original["program"],
        )
        assert not malformed_report["valid"]
        assert "malformed gzip" in malformed_report["checks"][-1]["detail"]

        policy_path = Path(original["policy"]); duplicate = json.loads(policy_path.read_text())
        duplicate["links"] = [{"source": "fixture.level", "target": "fixture.level", "kind": "bad"}]
        duplicate["artifact_sha256"] = digest({key: value for key, value in duplicate.items() if key != "artifact_sha256"})
        _write_json(duplicate, policy_path)
        try: load_policy(policy_path)
        except CompilerError: pass
        else: raise AssertionError("invalid policy link accepted")
    print("CONTRACT-ATTACKS-PASSED")


COMMANDS = {
    "boundary": boundary,
    "canonical": canonical,
    "contract-attacks": contract_attacks,
    "development-real-dna": development_real_dna,
    "genbank-real-dna": genbank_real_dna,
    "real-dna": real_dna,
    "resource-path-safety": resource_path_safety,
    "runtime-integration": runtime_integration,
    "sdist": sdist,
    "source-causality": source_causality,
    "tamper": tamper,
    "wheel": wheel,
}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gate", choices=sorted([*COMMANDS, "genbank-scale"]))
    parser.add_argument("--source", type=Path, help="local CP032762.1 GenBank source for genbank-scale")
    arguments = parser.parse_args()
    if arguments.gate == "genbank-scale":
        if arguments.source is None:
            parser.error("genbank-scale requires --source PATH; it never downloads source data")
        genbank_scale(arguments.source)
    else:
        if arguments.source is not None:
            parser.error("--source is only valid with genbank-scale")
        COMMANDS[arguments.gate]()

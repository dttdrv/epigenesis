"""Release gates for the compiler-only Epigenesis distribution."""

from __future__ import annotations

import argparse
import ast
import base64
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT))

from brainc._canonical import ContractError, canonical_bytes, digest, save_artifact
from brainc.compiler import CompilerError, compile_program, load_policy, load_program, save as save_program
from brainc.provider import ProviderError, load_response, make_request, save as save_provider
from brainc.sequence_collection import SequenceCollectionCompiler, load_sequence_collection
from brainc.validator import validate_chain


def _artifact(core: dict) -> dict:
    return {**core, "artifact_sha256": digest(core)}


def _write_json(value: dict, path: Path) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


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
    raw = (ROOT / "tests" / "data" / "J02482.1.fasta").read_bytes()
    metadata = json.loads((ROOT / "tests" / "data" / "J02482.1.source.json").read_text(encoding="utf-8"))
    assert metadata["accession"] == "J02482.1"
    assert hashlib.sha256(raw).hexdigest() == metadata["fasta_sha256"] == "2826ee08e3506154cdeb7ab734ec6f9ebbea6e5d16d0bec488030192faf492a7"
    bases = b"".join(raw.splitlines()[1:])
    assert len(bases) == metadata["bases"] == 5386
    assert hashlib.sha256(bases).hexdigest() == metadata["normalized_sequence_sha256"] == "97038c7e1edea2297667d7f0426ba942b322c74cb30e072ec66ba47f9c0448d0"
    with tempfile.TemporaryDirectory() as directory:
        chain = _build_real_chain(Path(directory))
        assert chain["report"]["valid"], chain["report"]
        source = json.loads(Path(chain["source"]).read_text(encoding="utf-8"))
        program = json.loads(Path(chain["program"]).read_text(encoding="utf-8"))
        assert source["members"][0]["sequence_artifact"]["sequence_ir"]["refget_id"] == "SQ.IIXILYBQCpHdC4qpI3sOQ_HAeAm9bmeF"
        assert program["program_ir"]["states"][0]["value"] == 1.0
    print("REAL-DNA-E2E-PASSED")


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
        load_sequence_collection(source_path)
        source_report = validate_chain(
            fasta=fresh["source_input"], sequence=source_path, manifest=fresh["manifest"],
            request=fresh["request"], response=fresh["response"], policy=fresh["policy"],
            program=fresh["program"],
        )
        assert not source_report["valid"]
        assert "does not replay from source input" in source_report["checks"][-1]["detail"]
    print("TAMPER-E2E-PASSED")


def boundary() -> None:
    denied = {"socket", "urllib", "requests", "httpx", "aiohttp", "subprocess", "importlib",
              "torch", "tensorflow", "jax", "sklearn"}
    for path in (ROOT / "brainc").glob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        imports: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import): imports.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module: imports.add(node.module.split(".")[0])
        assert not imports & denied, (path, imports & denied)
        text = path.read_text(encoding="utf-8").lower()
        for marker in ("train_model(", "optimizer.step(", "requests.get(", "subprocess.run("):
            assert marker not in text, (path, marker)
    print("COMPILER-BOUNDARY-PASSED")


def wheel() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory); wheels = root / "wheels"; installed = root / "installed"
        wheels.mkdir(); installed.mkdir()
        subprocess.run([sys.executable, "-m", "pip", "wheel", "--no-index", "--no-deps", "--no-build-isolation",
                        "--wheel-dir", str(wheels), str(ROOT)], check=True, capture_output=True, text=True)
        wheel_path, = wheels.glob("*.whl")
        subprocess.run([sys.executable, "-m", "pip", "install", "--no-index", "--no-deps", "--target", str(installed),
                        str(wheel_path)], check=True, capture_output=True, text=True)
        environment = dict(os.environ); environment["PYTHONPATH"] = str(installed)
        result = subprocess.run([sys.executable, "-m", "brainc", "--help"], check=True,
                                capture_output=True, text=True, cwd=root, env=environment)
        assert "compile-collection" in result.stdout and "validate" in result.stdout
    print("WHEEL-SMOKE-PASSED")


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


COMMANDS = {"real-dna": real_dna, "tamper": tamper, "boundary": boundary, "wheel": wheel,
            "canonical": canonical, "contract-attacks": contract_attacks}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("gate", choices=COMMANDS)
    COMMANDS[parser.parse_args().gate]()

from __future__ import annotations

import ast
import copy
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from brainc._canonical import digest, save_artifact
from brainc.compiler import CompilerError, compile_program, load_program, save as save_program
from brainc.provider import ProviderError, load_response, make_request, save as save_provider
from brainc.sequence import SequenceCompiler, SequenceCompilerError
from brainc.sequence_collection import SequenceCollectionCompiler, load_sequence_collection
from brainc.validator import validate_chain


def artifact(core: dict) -> dict:
    return {**core, "artifact_sha256": digest(core)}


def manifest(name: str, model_value: str, output_id: str) -> dict:
    return artifact({
        "format": "brainc.provider-manifest",
        "version": 1,
        "provider": {"name": name, "version": "1"},
        "model_identity": {"kind": "opaque", "value": model_value},
        "accepts": ["brain01.sequence-ir/v2", "brain01.sequence-collection-ir/v1"],
        "outputs": [{"id": output_id, "type": "number", "unit": None}],
    })


def policy(output_id: str, state_prefix: str) -> dict:
    return artifact({
        "format": "brainc.lowering-policy",
        "version": 1,
        "id": f"{state_prefix}-policy",
        "target": {"name": "generic-state-vm", "version": "1"},
        "states": [
            {
                "id": f"{state_prefix}.level",
                "type": "number",
                "from_output": output_id,
                "transform": {"scale": 2.0, "offset": 0.25, "clamp": {"minimum": 0.0, "maximum": 3.0}, "rounding": None},
            },
            {
                "id": f"{state_prefix}.slots",
                "type": "integer",
                "from_output": output_id,
                "transform": {"scale": 10.0, "offset": 1.0, "clamp": {"minimum": 1.0, "maximum": 64.0}, "rounding": "nearest"},
            },
        ],
        "links": [{"source": f"{state_prefix}.level", "target": f"{state_prefix}.slots", "kind": "declared-dependency"}],
        "ports": {"inputs": [f"{state_prefix}.level"], "outputs": [f"{state_prefix}.slots"]},
    })


class CompilerTests(unittest.TestCase):
    def run_chain(self, root: Path, provider_name: str, output_id: str, value: float) -> dict:
        root.mkdir(parents=True, exist_ok=True)
        fasta = root / "input.fa"
        fasta.write_text(">record-x independent sample\nACGTryswkmbdhvnACGT\n", encoding="utf-8")
        sequence_path = root / "sequence.json"
        SequenceCompiler().compile_file(fasta).save(sequence_path)
        manifest_path = root / "manifest.json"
        save_artifact(manifest(provider_name, f"{provider_name}-model", output_id), manifest_path)
        request_path = root / "request.json"
        request_value = make_request(sequence_path, manifest_path, [output_id])
        save_provider(request_value, request_path)
        response_path = root / "response.json"
        response_value = artifact({
            "format": "brainc.prediction-response",
            "version": 1,
            "request_artifact_sha256": request_value["artifact_sha256"],
            "provider": {"name": provider_name, "version": "1"},
            "model_identity": {"kind": "opaque", "value": f"{provider_name}-model"},
            "outputs": [{"id": output_id, "type": "number", "unit": None, "value": value}],
        })
        save_artifact(response_value, response_path)
        policy_path = root / "policy.json"
        save_artifact(policy(output_id, provider_name), policy_path)
        program_path = root / "program.json"
        save_program(compile_program(sequence_path, manifest_path, request_path, response_path, policy_path), program_path)
        report = validate_chain(fasta=fasta, sequence=sequence_path, manifest=manifest_path,
                                request=request_path, response=response_path, policy=policy_path,
                                program=program_path)
        self.assertTrue(report["valid"], report)
        return {"paths": (fasta, sequence_path, manifest_path, request_path, response_path, policy_path, program_path),
                "program": json.loads(program_path.read_text()), "report": report}

    def test_two_unrelated_external_providers(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            one = self.run_chain(root / "one", "provider-alpha", "alpha.signal", 0.2)
            two = self.run_chain(root / "two", "provider-zeta", "zeta.score", 0.8)
            self.assertNotEqual(one["program"]["artifact_sha256"], two["program"]["artifact_sha256"])

    def test_rehashed_program_forgery_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            result = self.run_chain(Path(directory), "provider-forge", "score", 0.4)
            fasta, sequence, manifest_path, request, response, policy_path, program_path = result["paths"]
            forged = copy.deepcopy(result["program"])
            forged["program_ir"]["states"][0]["value"] = 999.0
            forged["ir_sha256"] = digest(forged["program_ir"])
            forged["artifact_sha256"] = digest({key: value for key, value in forged.items() if key != "artifact_sha256"})
            program_path.write_text(json.dumps(forged), encoding="utf-8")
            report = validate_chain(fasta=fasta, sequence=sequence, manifest=manifest_path,
                                    request=request, response=response, policy=policy_path,
                                    program=program_path)
            self.assertFalse(report["valid"])
            self.assertIn("independent lowering replay", report["checks"][-1]["detail"])

    def test_sequence_sensitivity_and_determinism(self) -> None:
        one = SequenceCompiler().compile_text(">x\nACGT\n").to_dict()
        again = SequenceCompiler().compile_text(">x\nACGT\n").to_dict()
        changed = SequenceCompiler().compile_text(">x\nACGA\n").to_dict()
        self.assertEqual(one, again)
        self.assertNotEqual(one["artifact_sha256"], changed["artifact_sha256"])

    def test_context_rejects_ambiguous_text_and_null_digest(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); fasta = root / "input.fa"; context = root / "context.json"
            fasta.write_text(">x\nACGT\n", encoding="utf-8")
            base = {
                "format": "brain01.sequence-context", "version": 1, "record_id": "x",
                "reference": {"assembly": "a", "contig": "c", "start": 0, "end": 4,
                              "coordinate_system": "0-based-half-open", "orientation": "forward",
                              "aliases": [" x "]},
                "provenance": [],
            }
            context.write_text(json.dumps(base), encoding="utf-8")
            with self.assertRaises(SequenceCompilerError): SequenceCompiler().compile_file(fasta, context)
            base["reference"]["aliases"] = []
            base["provenance"] = [{"id": "p", "kind": "database", "uri": "https://example.test",
                                   "version": "1", "sha256": None}]
            context.write_text(json.dumps(base), encoding="utf-8")
            with self.assertRaises(SequenceCompilerError): SequenceCompiler().compile_file(fasta, context)

    def test_multi_fasta_and_gzip_ingress(self) -> None:
        source = b">one\nAC\nGT\n>two desc\nNN\nry\n>three\nA\n"
        plain = SequenceCollectionCompiler().compile_fasta_bytes(source)
        wrapped = SequenceCollectionCompiler().compile_fasta_bytes(gzip.compress(source, mtime=0))
        self.assertEqual([item.artifact.record_id for item in plain.members], ["one", "two", "three"])
        self.assertEqual([[segment.line for segment in item.artifact.source_map] for item in plain.members], [[2, 3], [2, 3], [2]])
        self.assertEqual([[segment.line for segment in item.input_source_map] for item in plain.members], [[2, 3], [5, 6], [8]])
        self.assertEqual(plain.to_dict()["collection_ir"], wrapped.to_dict()["collection_ir"])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "collection.json"
            wrapped.save(path)
            self.assertEqual(load_sequence_collection(path).to_dict(), wrapped.to_dict())

    def test_collection_reaches_target_without_bundled_predictor(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = root / "input.fa"
            input_path.write_bytes(b">chrA\nACGT\n>chrB\nNNNN\n")
            source_path = root / "collection.json"
            SequenceCollectionCompiler().compile_file(input_path).save(source_path)
            manifest_path = root / "manifest.json"
            save_artifact(manifest("collection-provider", "whole-input-v1", "collection.score"), manifest_path)
            request_path = root / "request.json"
            request_value = make_request(source_path, manifest_path, ["collection.score"])
            save_provider(request_value, request_path)
            response_path = root / "response.json"
            save_artifact(artifact({
                "format": "brainc.prediction-response",
                "version": 1,
                "request_artifact_sha256": request_value["artifact_sha256"],
                "provider": {"name": "collection-provider", "version": "1"},
                "model_identity": {"kind": "opaque", "value": "whole-input-v1"},
                "outputs": [{"id": "collection.score", "type": "number", "unit": None, "value": 0.6}],
            }), response_path)
            policy_path = root / "policy.json"
            save_artifact(policy("collection.score", "collection"), policy_path)
            program_path = root / "program.json"
            program = compile_program(source_path, manifest_path, request_path, response_path, policy_path)
            save_program(program, program_path)
            self.assertEqual(program["sources"]["sequence"]["format"], "brain01.sequence-collection-ir")
            self.assertEqual(program["program_ir"]["states"][0]["value"], 1.45)
            report = validate_chain(fasta=input_path, sequence=source_path, manifest=manifest_path,
                                    request=request_path, response=response_path, policy=policy_path,
                                    program=program_path)
            self.assertTrue(report["valid"], report)

    def test_raw_gzip_collection_replays_independently(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = root / "raw.dna.gz"
            input_path.write_bytes(gzip.compress(b"ACGTRYSWKMBDHVN", mtime=0))
            source_path = root / "collection.json"
            SequenceCollectionCompiler().compile_raw(input_path.read_bytes(), "raw-1").save(source_path)
            manifest_path = root / "manifest.json"
            save_artifact(manifest("raw-provider", "raw-model", "score"), manifest_path)
            request_path = root / "request.json"
            request_value = make_request(source_path, manifest_path, ["score"])
            save_provider(request_value, request_path)
            response_path = root / "response.json"
            save_artifact(artifact({
                "format": "brainc.prediction-response", "version": 1,
                "request_artifact_sha256": request_value["artifact_sha256"],
                "provider": {"name": "raw-provider", "version": "1"},
                "model_identity": {"kind": "opaque", "value": "raw-model"},
                "outputs": [{"id": "score", "type": "number", "unit": None, "value": 0.25}],
            }), response_path)
            policy_path = root / "policy.json"; save_artifact(policy("score", "raw"), policy_path)
            program_path = root / "program.json"
            save_program(compile_program(source_path, manifest_path, request_path, response_path, policy_path), program_path)
            report = validate_chain(fasta=input_path, record_id="raw-1", sequence=source_path,
                                    manifest=manifest_path, request=request_path, response=response_path,
                                    policy=policy_path, program=program_path)
            self.assertTrue(report["valid"], report)

    def test_safe_integer_and_manifest_contract_boundaries(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            result = self.run_chain(Path(directory), "provider-boundary", "score", 0.4)
            fasta, sequence, manifest_path, request_path, response_path, policy_path, _ = result["paths"]
            invalid_manifest = manifest("provider-boundary", "provider-boundary-model", "other")
            save_artifact(invalid_manifest, manifest_path)
            with self.assertRaises(ProviderError):
                compile_program(sequence, manifest_path, request_path, response_path, policy_path)
            unsafe = {
                "format": "brainc.prediction-response", "version": 1,
                "request_artifact_sha256": "0" * 64,
                "provider": {"name": "p", "version": "1"},
                "model_identity": {"kind": "opaque", "value": "m"},
                "outputs": [{"id": "n", "type": "integer", "unit": None, "value": 2**53}],
                "artifact_sha256": "0" * 64,
            }
            response_path.write_text(json.dumps(unsafe), encoding="utf-8")
            with self.assertRaises((ProviderError, ValueError)):
                load_response(response_path)

    def test_closed_program_schema_rejects_rehashed_nonsense(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            result = self.run_chain(Path(directory), "provider-schema", "score", 0.4)
            program_path = result["paths"][-1]
            forged = copy.deepcopy(result["program"])
            forged["program_ir"]["ports"]["outputs"].append("missing")
            forged["ir_sha256"] = digest(forged["program_ir"])
            forged["artifact_sha256"] = digest({key: value for key, value in forged.items() if key != "artifact_sha256"})
            program_path.write_text(json.dumps(forged), encoding="utf-8")
            with self.assertRaises(CompilerError): load_program(program_path)

    def test_validator_has_no_compiler_imports(self) -> None:
        tree = ast.parse((Path(__file__).parents[1] / "brainc" / "validator.py").read_text())
        imports = {
            node.module for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom) and node.module is not None
        }
        self.assertFalse(any(name == "brainc" or name.startswith("brainc.") for name in imports))

    def test_core_has_no_network_dynamic_loading_or_unscoped_process_execution(self) -> None:
        denied_imports = {"socket", "urllib", "requests", "httpx", "aiohttp", "importlib"}
        process_owners = {"external_runtime.py", "validator_external.py"}
        for path in (Path(__file__).parents[1] / "brainc").glob("*.py"):
            tree = ast.parse(path.read_text())
            imports = set()
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    imports.update(alias.name.split(".")[0] for alias in node.names)
                elif isinstance(node, ast.ImportFrom) and node.module:
                    imports.add(node.module.split(".")[0])
            self.assertFalse(imports & denied_imports, (path, imports & denied_imports))
            if "subprocess" in imports:
                self.assertIn(path.name, process_owners)

    def test_cli_exposes_no_training_or_runtime_commands(self) -> None:
        result = subprocess.run([sys.executable, "-m", "brainc", "--help"], check=True,
                                text=True, capture_output=True)
        for command in ("train-model", "run-world", "simulate", "finalize", "checkpoint"):
            self.assertNotIn(command, result.stdout)


if __name__ == "__main__":
    unittest.main()

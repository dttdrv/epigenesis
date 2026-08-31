from __future__ import annotations

import ast
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

from brainc import validator as legacy_validator
from brainc._canonical import (
    ContractError,
    digest,
    load as load_artifact,
    loads as loads_artifact,
    save_artifact,
)
from brainc.compiler import compile_program, load_policy, load_program, save as save_program
from brainc.provider import (
    load_manifest,
    load_request,
    load_response,
    load_source,
    make_request,
    save as save_provider,
)
from brainc.sequence import SequenceCompiler
from brainc.sequence_collection import SequenceCollectionCompiler
from brainc.validator import ValidationError, save_report, validate_chain
import brainc._canonical as canonical_io


def _seal(core: dict) -> dict:
    return {**core, "artifact_sha256": digest(core)}


def _build_chain(root: Path, *, collection: bool = False, compressed: bool = False) -> dict:
    root.mkdir(parents=True, exist_ok=True)
    fasta = b">legacy sample\nACGTN\n"
    source_input = gzip.compress(fasta, mtime=0) if compressed else fasta
    input_path = root / ("input.fa.gz" if compressed else "input.fa")
    input_path.write_bytes(source_input)
    source_path = root / "source.json"
    if collection:
        SequenceCollectionCompiler().compile_fasta_bytes(source_input).save(source_path)
        source_tag = "brain01.sequence-collection-ir/v1"
    else:
        if compressed:
            raise ValueError("single-sequence ingress does not accept gzip")
        SequenceCompiler().compile_file(input_path).save(source_path)
        source_tag = "brain01.sequence-ir/v2"

    manifest = _seal(
        {
            "format": "brainc.provider-manifest",
            "version": 1,
            "provider": {"name": "legacy-provider", "version": "1"},
            "model_identity": {"kind": "opaque", "value": "legacy-model"},
            "accepts": [source_tag],
            "outputs": [{"id": "score", "type": "number", "unit": None}],
        }
    )
    manifest_path = root / "manifest.json"
    save_artifact(manifest, manifest_path)
    request = make_request(source_path, manifest_path, ["score"])
    request_path = root / "request.json"
    save_provider(request, request_path)
    response = _seal(
        {
            "format": "brainc.prediction-response",
            "version": 1,
            "request_artifact_sha256": request["artifact_sha256"],
            "provider": manifest["provider"],
            "model_identity": manifest["model_identity"],
            "outputs": [
                {
                    "id": "score",
                    "type": "number",
                    "unit": None,
                    "value": 0.25,
                }
            ],
        }
    )
    response_path = root / "response.json"
    save_artifact(response, response_path)
    policy = _seal(
        {
            "format": "brainc.lowering-policy",
            "version": 1,
            "id": "legacy-policy",
            "target": {"name": "generic-state-vm", "version": "1"},
            "states": [
                {
                    "id": "level",
                    "type": "number",
                    "from_output": "score",
                    "transform": {
                        "scale": 2.0,
                        "offset": 0.25,
                        "clamp": {"minimum": 0.0, "maximum": 3.0},
                        "rounding": None,
                    },
                }
            ],
            "links": [],
            "ports": {"inputs": ["level"], "outputs": ["level"]},
        }
    )
    policy_path = root / "policy.json"
    save_artifact(policy, policy_path)
    program = compile_program(
        source_path, manifest_path, request_path, response_path, policy_path
    )
    program_path = root / "program.json"
    save_program(program, program_path)
    return {
        "source_input": input_path,
        "source": source_path,
        "manifest": manifest_path,
        "request": request_path,
        "response": response_path,
        "policy": policy_path,
        "program": program_path,
        "values": {
            "manifest": manifest,
            "request": request,
            "response": response,
            "policy": policy,
            "program": program,
        },
        "logical_source": fasta,
    }


class LegacyIOLimitTests(unittest.TestCase):
    def test_legacy_limit_constants_match_the_public_compiler_contract(self) -> None:
        self.assertEqual(canonical_io.MAX_JSON_BYTES, 16 * 1024 * 1024)
        self.assertEqual(legacy_validator.MAX_INPUT_BYTES, 16 * 1024 * 1024)
        self.assertEqual(legacy_validator.MAX_JSON_BYTES, 16 * 1024 * 1024)
        self.assertEqual(legacy_validator.MAX_DECOMPRESSED_BYTES, 64 * 1024 * 1024)
        self.assertEqual(legacy_validator.MAX_STRING_BYTES, 1024 * 1024)
        self.assertEqual(legacy_validator.MAX_JSON_DEPTH, 64)
        self.assertEqual(legacy_validator.MAX_JSON_MEMBERS, 1_000_000)

    def test_canonical_json_exact_boundaries_and_finite_i_json(self) -> None:
        with mock.patch.object(canonical_io, "MAX_JSON_BYTES", 2):
            self.assertEqual(loads_artifact(b"{}", "artifact"), {})
            with self.assertRaisesRegex(ContractError, "JSON byte limit 2"):
                loads_artifact(b"{} ", "artifact")
        with mock.patch.object(canonical_io, "MAX_STRING_BYTES", 2):
            self.assertEqual(loads_artifact(b'{"x":"ab"}', "artifact"), {"x": "ab"})
            with self.assertRaisesRegex(ContractError, "string byte limit 2"):
                loads_artifact(b'{"x":"abc"}', "artifact")
        with mock.patch.object(canonical_io, "MAX_JSON_MEMBERS", 1):
            self.assertEqual(loads_artifact(b'{"x":[]}', "artifact"), {"x": []})
            with self.assertRaisesRegex(ContractError, "member limit 1"):
                loads_artifact(b'{"x":[],"y":[]}', "artifact")
        with mock.patch.object(canonical_io, "MAX_JSON_DEPTH", 1):
            self.assertEqual(loads_artifact(b'{"x":[]}', "artifact"), {"x": []})
            with self.assertRaisesRegex(ContractError, "depth limit 1"):
                loads_artifact(b'{"x":[[]]}', "artifact")
        with self.assertRaisesRegex(ContractError, "duplicate JSON key"):
            loads_artifact(b'{"x":1,"x":2}', "artifact")
        with self.assertRaisesRegex(ContractError, "safe range"):
            loads_artifact(b'{"x":9007199254740992}', "artifact")
        with self.assertRaisesRegex(ContractError, "non-finite"):
            loads_artifact(b'{"x":1e400}', "artifact")

    def test_independent_json_decoder_has_the_same_exact_boundaries(self) -> None:
        with mock.patch.object(legacy_validator, "MAX_JSON_BYTES", 2):
            self.assertEqual(legacy_validator._decode(b"{}", "artifact"), {})
            with self.assertRaisesRegex(ValidationError, "at most 2 JSON bytes"):
                legacy_validator._decode(b"{} ", "artifact")
        with mock.patch.object(legacy_validator, "MAX_STRING_BYTES", 2):
            self.assertEqual(
                legacy_validator._decode(b'{"x":"ab"}', "artifact"),
                {"x": "ab"},
            )
            with self.assertRaisesRegex(ValidationError, "string byte limit 2"):
                legacy_validator._decode(b'{"x":"abc"}', "artifact")
        with self.assertRaisesRegex(ValidationError, "duplicate JSON key"):
            legacy_validator._decode(b'{"x":1,"x":2}', "artifact")
        with self.assertRaisesRegex(ValidationError, "safe range"):
            legacy_validator._decode(b'{"x":9007199254740992}', "artifact")

    @unittest.skipUnless(hasattr(os, "symlink"), "symbolic links are unavailable")
    def test_every_public_legacy_artifact_loader_rejects_symlinks(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            chain = _build_chain(root / "chain")
            loaders = {
                "source": load_source,
                "manifest": load_manifest,
                "request": load_request,
                "response": load_response,
                "policy": load_policy,
                "program": load_program,
            }
            for name, loader in loaders.items():
                with self.subTest(name=name):
                    linked = root / f"{name}-linked.json"
                    linked.symlink_to(chain[name])
                    with self.assertRaisesRegex(ContractError, "regular non-linked file"):
                        loader(linked)

    def test_every_public_legacy_artifact_loader_rejects_oversized_input(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            chain = _build_chain(root / "chain")
            loaders = {
                "source": load_source,
                "manifest": load_manifest,
                "request": load_request,
                "response": load_response,
                "policy": load_policy,
                "program": load_program,
            }
            for name, loader in loaders.items():
                with self.subTest(name=name), mock.patch.object(
                    canonical_io, "MAX_JSON_BYTES", 16
                ):
                    oversized = root / f"{name}-oversized.json"
                    oversized.write_bytes(b"{" + b" " * 16 + b"}")
                    with self.assertRaisesRegex(ContractError, "byte limit 16"):
                        loader(oversized)

    def test_canonical_regular_file_detects_path_replacement(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "artifact.json"
            replacement = root / "replacement.json"
            source.write_bytes(b"{}")
            replacement.write_bytes(b"{}")
            real_open = os.open

            def replace_then_open(path: str | bytes | os.PathLike[str], flags: int) -> int:
                os.replace(replacement, source)
                return real_open(path, flags)

            with mock.patch("brainc._io.os.open", side_effect=replace_then_open):
                with self.assertRaisesRegex(ContractError, "changed while being opened"):
                    load_artifact(source, "artifact")

    def test_legacy_provider_validates_one_descriptor_snapshot_in_memory(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            single = _build_chain(root / "single")["source"]
            collection_input = b">one\nACGT\n>two\nNN\n"
            collection = root / "collection.json"
            SequenceCollectionCompiler().compile_fasta_bytes(collection_input).save(
                collection
            )
            real_open = os.open
            for source in (single, collection):
                with self.subTest(source=source.name):
                    opens = 0

                    def counted_open(
                        path: str | bytes | os.PathLike[str], flags: int
                    ) -> int:
                        nonlocal opens
                        opens += 1
                        return real_open(path, flags)

                    with mock.patch("brainc._io.os.open", side_effect=counted_open):
                        loaded = load_source(source)
                    self.assertEqual(opens, 1)
                    stored = json.loads(source.read_text(encoding="utf-8"))
                    self.assertEqual(
                        loaded["artifact_sha256"], stored["artifact_sha256"]
                    )

    def test_legacy_artifact_output_is_bounded_atomic_and_never_follows_symlink(self) -> None:
        payload = {"format": "example", "version": 1}
        expected = (
            json.dumps(
                payload,
                indent=2,
                sort_keys=True,
                ensure_ascii=False,
                allow_nan=False,
            )
            + "\n"
        ).encode("utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "artifact.json"
            with mock.patch.object(canonical_io, "MAX_JSON_BYTES", len(expected)):
                save_artifact(payload, output)
            self.assertEqual(output.read_bytes(), expected)
            with mock.patch.object(canonical_io, "MAX_JSON_BYTES", len(expected) - 1):
                with self.assertRaisesRegex(ContractError, "serialized output exceeds"):
                    save_artifact(payload, output)
            self.assertEqual(output.read_bytes(), expected)

            if hasattr(os, "symlink"):
                target = root / "target.json"
                target.write_bytes(b"sentinel")
                linked = root / "linked.json"
                linked.symlink_to(target)
                with self.assertRaisesRegex(ContractError, "regular non-linked file"):
                    save_artifact(payload, linked)
                self.assertEqual(target.read_bytes(), b"sentinel")

    def test_independent_report_output_is_atomic_and_never_follows_symlink(self) -> None:
        report = {"format": "brainc.validation-report", "valid": False}
        expected = (
            json.dumps(
                report,
                indent=2,
                sort_keys=True,
                ensure_ascii=False,
                allow_nan=False,
            )
            + "\n"
        ).encode("utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "report.json"
            output.write_bytes(b"sentinel")
            with mock.patch(
                "brainc.validator.os.replace", side_effect=OSError("injected")
            ):
                with self.assertRaisesRegex(ValidationError, "cannot write"):
                    save_report(report, output)
            self.assertEqual(output.read_bytes(), b"sentinel")
            self.assertEqual(list(root.glob(".report.json.*")), [])

            with mock.patch.object(legacy_validator, "MAX_JSON_BYTES", len(expected)):
                save_report(report, output)
            self.assertEqual(output.read_bytes(), expected)
            with mock.patch.object(
                legacy_validator, "MAX_JSON_BYTES", len(expected) - 1
            ):
                with self.assertRaisesRegex(ValidationError, "byte limit"):
                    save_report(report, output)
            self.assertEqual(output.read_bytes(), expected)

            if hasattr(os, "symlink"):
                linked = root / "linked-report.json"
                linked.symlink_to(output)
                with self.assertRaisesRegex(ValidationError, "regular non-linked file"):
                    save_report(report, linked)
                self.assertEqual(output.read_bytes(), expected)

    def test_independent_regular_file_detects_path_replacement(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "artifact.json"
            replacement = root / "replacement.json"
            source.write_bytes(b"{}")
            replacement.write_bytes(b"{}")
            real_open = os.open

            def replace_then_open(path: str | bytes | os.PathLike[str], flags: int) -> int:
                os.replace(replacement, source)
                return real_open(path, flags)

            with mock.patch("brainc.validator.os.open", side_effect=replace_then_open):
                with self.assertRaisesRegex(ValidationError, "changed or is not"):
                    legacy_validator._read(source, "artifact")

    def test_independent_gzip_replay_boundaries_multistream_and_trailing_policy(self) -> None:
        logical = b">legacy\nACGT\n"
        compressed = gzip.compress(logical, mtime=0)
        with mock.patch.object(
            legacy_validator, "MAX_INPUT_BYTES", len(compressed)
        ), mock.patch.object(
            legacy_validator, "MAX_DECOMPRESSED_BYTES", len(logical)
        ):
            self.assertEqual(legacy_validator._unwrap(compressed)[0], logical)
        with mock.patch.object(
            legacy_validator, "MAX_INPUT_BYTES", len(compressed) - 1
        ):
            with self.assertRaisesRegex(ValidationError, "input exceeds byte limit"):
                legacy_validator._unwrap(compressed)
        with mock.patch.object(
            legacy_validator, "MAX_DECOMPRESSED_BYTES", len(logical) - 1
        ):
            with self.assertRaisesRegex(ValidationError, "decompressed byte limit"):
                legacy_validator._unwrap(compressed)

        first = gzip.compress(b">legacy\n", mtime=0)
        second = gzip.compress(b"ACGT\n", mtime=0)
        self.assertEqual(legacy_validator._unwrap(first + second)[0], logical)
        self.assertEqual(legacy_validator._unwrap(compressed + b"\0\0")[0], logical)
        with self.assertRaisesRegex(ValidationError, "malformed gzip source"):
            legacy_validator._unwrap(compressed + b"junk")

    def test_public_validation_replay_enforces_source_limits_and_safe_files(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            chain = _build_chain(root / "chain", collection=True, compressed=True)
            arguments = {
                "fasta": chain["source_input"],
                "sequence": chain["source"],
                "manifest": chain["manifest"],
                "request": chain["request"],
                "response": chain["response"],
                "policy": chain["policy"],
                "program": chain["program"],
            }
            report = validate_chain(**arguments)
            self.assertTrue(report["valid"], report)

            compressed_size = chain["source_input"].stat().st_size
            with mock.patch.object(
                legacy_validator, "MAX_INPUT_BYTES", compressed_size - 1
            ):
                rejected = validate_chain(**arguments)
            self.assertFalse(rejected["valid"])
            self.assertIn("byte limit", rejected["checks"][-1]["detail"])

            if hasattr(os, "symlink"):
                targets = {
                    "fasta": chain["source_input"],
                    "sequence": chain["source"],
                    "manifest": chain["manifest"],
                    "request": chain["request"],
                    "response": chain["response"],
                    "policy": chain["policy"],
                    "program": chain["program"],
                }
                for argument, target in targets.items():
                    with self.subTest(argument=argument):
                        linked = root / f"linked-{argument}"
                        linked.symlink_to(target)
                        linked_arguments = dict(arguments)
                        linked_arguments[argument] = linked
                        rejected = validate_chain(**linked_arguments)
                        self.assertFalse(rejected["valid"])
                        self.assertIn(
                            "regular non-linked file",
                            rejected["checks"][-1]["detail"],
                        )

                context_target = root / "context.json"
                context_target.write_text("{}", encoding="utf-8")
                linked_context = root / "linked-context.json"
                linked_context.symlink_to(context_target)
                rejected = validate_chain(**arguments, context=linked_context)
                self.assertFalse(rejected["valid"])
                self.assertIn(
                    "regular non-linked file", rejected["checks"][-1]["detail"]
                )

    def test_v04_hashes_serialized_bytes_and_independent_validation_are_unchanged(self) -> None:
        expected_hashes = {
            "manifest": "333ed3f26df65af7100a79da76d35290d7f3a38447a105943a006b501a68a071",
            "request": "4bb2079e1e593ddd20d2573067535ab77b9ecb21e218650b2e1e123683cee85a",
            "response": "e5973737a214a7628a082e6d9cd93e224a0f60310c825136426b1a6b6297d9c5",
            "policy": "0576feefa04efd977355a0ab678bb9ac245a8932e8411c5e3e8907d19a51c404",
            "program": "487d1703916a30f17e2b8d76db1b92a9eda334448453bef094acd3edde4d485a",
        }
        with tempfile.TemporaryDirectory() as directory:
            chain = _build_chain(Path(directory))
            for name, expected_hash in expected_hashes.items():
                with self.subTest(name=name):
                    self.assertEqual(
                        chain["values"][name]["artifact_sha256"], expected_hash
                    )
                    expected_bytes = (
                        json.dumps(
                            chain["values"][name],
                            indent=2,
                            sort_keys=True,
                            ensure_ascii=False,
                            allow_nan=False,
                        )
                        + "\n"
                    ).encode("utf-8")
                    self.assertEqual(chain[name].read_bytes(), expected_bytes)
            report = validate_chain(
                fasta=chain["source_input"],
                sequence=chain["source"],
                manifest=chain["manifest"],
                request=chain["request"],
                response=chain["response"],
                policy=chain["policy"],
                program=chain["program"],
            )
            self.assertTrue(report["valid"], report)

    def test_legacy_cli_check_compile_and_validate_paths_remain_compatible(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            chain = _build_chain(root)
            for kind in (
                "sequence",
                "manifest",
                "request",
                "response",
                "policy",
                "program",
            ):
                path = chain["source"] if kind == "sequence" else chain[kind]
                result = subprocess.run(
                    [sys.executable, "-m", "brainc", "check", kind, str(path)],
                    text=True,
                    capture_output=True,
                    check=False,
                )
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "brainc",
                    "validate",
                    "--source-input",
                    str(chain["source_input"]),
                    "--source-artifact",
                    str(chain["source"]),
                    "--manifest",
                    str(chain["manifest"]),
                    "--request",
                    str(chain["request"]),
                    "--response",
                    str(chain["response"]),
                    "--policy",
                    str(chain["policy"]),
                    "--program",
                    str(chain["program"]),
                    "--report",
                    str(root / "validation.json"),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertTrue(json.loads(result.stdout)["valid"])

    def test_independent_validator_still_has_no_brainc_imports(self) -> None:
        source = Path(legacy_validator.__file__).read_text(encoding="utf-8")
        tree = ast.parse(source)
        imports = {
            node.module
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom) and node.module is not None
        }
        self.assertFalse(any(name == "brainc" or name.startswith("brainc.") for name in imports))


if __name__ == "__main__":
    unittest.main()

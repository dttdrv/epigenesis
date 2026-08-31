"""Adversarial tests for the independent v2 sequence-source boundary."""

from __future__ import annotations

import base64
from copy import deepcopy
from dataclasses import replace
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

from brainc.sequence import SequenceCompiler
from brainc.sequence_collection import SequenceCollectionCompiler
from brainc.insdc import GenBankCompiler
from brainc.v2 import policy_artifact, target_artifact
from brainc.validator_v2 import (
    DEFAULT_LIMITS,
    MODULE_PRODUCER,
    ValidationError,
    digest,
    save_report,
    source_record_lengths,
    validate_chain,
    validate_module_artifact,
    validate_request_artifact,
    validate_response_artifact,
    validate_source_artifact,
    validate_target_artifact,
)
import brainc.validator_v2 as independent_validator


ROOT = Path(__file__).parents[1]
U49845 = ROOT / "tests/data/U49845.1.gb"


def _seal(value: dict) -> dict:
    result = deepcopy(value)
    result.pop("artifact_sha256", None)
    result["artifact_sha256"] = digest(result)
    return result


def _reseal_source(value: dict) -> dict:
    result = deepcopy(value)
    if result.get("format") == "brain01.sequence-ir":
        result["ir_sha256"] = digest(result["sequence_ir"])
    else:
        result["collection_ir_sha256"] = digest(result["collection_ir"])
    return _seal(result)


def _write_json(path: Path, value: dict, *, padding: int = 0) -> Path:
    path.write_bytes(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
        + b" " * padding
    )
    return path


class IndependentSourceValidatorTests(unittest.TestCase):
    def setUp(self) -> None:
        self.fasta = b">alpha description\nACGT\nTG\n>beta\nNNNN\n"
        self.collection = SequenceCollectionCompiler().compile_fasta_bytes(
            self.fasta
        ).to_dict()
        self.sequence_fasta = b">single arbitrary description\nacgTRYSW\nKMBDHVN\n"
        self.sequence = SequenceCompiler().compile_bytes(self.sequence_fasta).to_dict()

    def _genbank(self) -> dict:
        return GenBankCompiler().compile_bytes(U49845.read_bytes()).to_dict()

    def _genbank_chain(self, root: Path) -> tuple[dict, dict[str, Path]]:
        source = self._genbank()
        source_binding = {
            "format": source["format"],
            "version": source["version"],
            "artifact_sha256": source["artifact_sha256"],
            "ir_sha256": source["bio_ir_sha256"],
        }
        axis = {
            "kind": "source-coordinate",
            "source_artifact_sha256": source["artifact_sha256"],
            "record_id": "U49845.1",
            "start": 0,
            "step": 1,
            "span": 1,
            "coordinate_system": "0-based-half-open",
        }
        contract = {
            "id": "coordinate.values",
            "type": {"dtype": "bool", "shape": [2]},
            "unit": None,
            "axes": [axis],
        }
        manifest = _seal(
            {
                "format": "brainc.provider-manifest",
                "version": 2,
                "provider": {"name": "source-test", "version": "1"},
                "model_identity": {
                    "kind": "content-sha256",
                    "value": "2" * 64,
                },
                "accepts": ["brainc.bio.insdc-genbank-ir/v2"],
                "outputs": [contract],
            }
        )
        request = _seal(
            {
                "format": "brainc.prediction-request",
                "version": 2,
                "source": source_binding,
                "provider_manifest_sha256": manifest["artifact_sha256"],
                "requested_outputs": [contract],
            }
        )
        raw_tensor = b"\x01\x00"
        output = {
            **contract,
            "storage": {
                "kind": "inline-base64",
                "data": base64.b64encode(raw_tensor).decode("ascii"),
                "byte_length": len(raw_tensor),
                "sha256": hashlib.sha256(raw_tensor).hexdigest(),
            },
        }
        response = _seal(
            {
                "format": "brainc.prediction-response",
                "version": 2,
                "request_artifact_sha256": request["artifact_sha256"],
                "provider": manifest["provider"],
                "model_identity": manifest["model_identity"],
                "outputs": [output],
            }
        )
        target = target_artifact(
            {
                "id": "org.example.source-test",
                "abi_major": 1,
                "opsets": [{"domain": "org.example.noop", "version": 1}],
                "unit_schemas": [
                    {
                        "id": "node",
                        "fields": [
                            {
                                "id": "state",
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
        target_binding = {
            "id": target["contract"]["id"],
            "abi_major": target["contract"]["abi_major"],
            "contract_sha256": target["contract_sha256"],
        }
        policy = policy_artifact(
            "org.example.source-test",
            target_binding,
            [{"id": contract["id"], "from_output": contract["id"]}],
            [],
        )
        body = {
            "target": target_binding,
            "requirements": [],
            "budgets": {
                "operations": 0,
                "tensor_bytes": len(raw_tensor),
                "units": 0,
                "edges": 0,
                "attachments": 0,
            },
            "tensors": [
                {
                    **output,
                    "lineage": {
                        "response_output": contract["id"],
                        "lowering_value": contract["id"],
                    },
                }
            ],
            "entrypoint": {"name": "develop", "operations": []},
        }
        module = _seal(
            {
                "format": "brainc.development-module",
                "version": 1,
                "producer": deepcopy(MODULE_PRODUCER),
                "sources": {
                    "sequence": source_binding,
                    "provider_manifest": {
                        "artifact_sha256": manifest["artifact_sha256"]
                    },
                    "prediction_request": {
                        "artifact_sha256": request["artifact_sha256"]
                    },
                    "prediction_response": {
                        "artifact_sha256": response["artifact_sha256"]
                    },
                    "lowering_policy": {
                        "artifact_sha256": policy["artifact_sha256"]
                    },
                    "target_contract": {
                        "artifact_sha256": target["artifact_sha256"]
                    },
                },
                "module": body,
                "module_sha256": digest(body),
            }
        )
        paths = {
            "fasta": U49845,
            "sequence": _write_json(root / "source.json", source),
            "manifest": _write_json(root / "manifest.json", manifest),
            "request": _write_json(root / "request.json", request),
            "response": _write_json(root / "response.json", response),
            "policy": _write_json(root / "policy.json", policy),
            "target": _write_json(root / "target.json", target),
            "module": _write_json(root / "module.json", module),
        }
        return {
            "source": source,
            "request": request,
            "target": target,
            "module": module,
        }, paths

    def assertSourceRejected(self, value: dict) -> None:
        with self.assertRaises(ValidationError):
            validate_source_artifact(value)

    def test_real_sequence_and_collection_validate_and_expose_exact_lengths(self) -> None:
        self.assertEqual(validate_source_artifact(self.sequence), self.sequence)
        self.assertEqual(source_record_lengths(self.sequence), {"single": 15})
        self.assertEqual(validate_source_artifact(self.collection), self.collection)
        lengths = source_record_lengths(self.collection)
        self.assertEqual(lengths, {"alpha": 6, "beta": 4})
        lengths["alpha"] = 99
        self.assertEqual(source_record_lengths(self.collection)["alpha"], 6)

    def test_u49845_genbank_validates_across_every_source_binding(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            artifacts, paths = self._genbank_chain(Path(temporary))
            source = artifacts["source"]
            self.assertEqual(validate_source_artifact(source), source)
            self.assertEqual(
                validate_source_artifact(source, raw_source=U49845.read_bytes()),
                source,
            )
            self.assertEqual(source_record_lengths(source), {"U49845.1": 5028})
            self.assertEqual(validate_request_artifact(artifacts["request"]), artifacts["request"])
            self.assertEqual(
                validate_module_artifact(artifacts["module"], artifacts["target"]),
                artifacts["module"],
            )
            report = validate_chain(**paths)
            self.assertTrue(report["valid"], report)
            self.assertEqual(report["passed_checks"], report["total_checks"])

    def test_coherently_resealed_genbank_semantic_tamper_is_rejected(self) -> None:
        tampered = self._genbank()
        record = tampered["bio_ir"]["records"][0]
        record["definition"] = "coherently resealed forgery"
        record["record_ir_sha256"] = digest(
            {key: value for key, value in record.items() if key != "record_ir_sha256"}
        )
        tampered["bio_ir_sha256"] = digest(tampered["bio_ir"])
        tampered = _seal(tampered)
        with self.assertRaisesRegex(ValidationError, "independent raw GenBank replay"):
            validate_source_artifact(tampered, raw_source=U49845.read_bytes())

    def test_source_format_specific_admission_limits(self) -> None:
        self.assertEqual(DEFAULT_LIMITS.input_bytes, 16 * 1024 * 1024)
        self.assertEqual(DEFAULT_LIMITS.source_artifact_bytes, 64 * 1024 * 1024)
        genbank = self._genbank()
        genbank_raw = json.dumps(genbank, separators=(",", ":")).encode("utf-8")
        generic_raw = json.dumps(self.collection, separators=(",", ":")).encode("utf-8")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genbank_path = root / "genbank.json"
            generic_path = root / "generic.json"
            genbank_path.write_bytes(genbank_raw)
            generic_path.write_bytes(generic_raw)
            genbank_limits = replace(
                DEFAULT_LIMITS,
                input_bytes=len(genbank_raw) - 1,
                source_artifact_bytes=len(genbank_raw),
            )
            loaded, _ = independent_validator._load_source_artifact(
                genbank_path, "source artifact", genbank_limits
            )
            self.assertEqual(loaded["artifact_sha256"], genbank["artifact_sha256"])

            generic_limits = replace(
                DEFAULT_LIMITS,
                input_bytes=len(generic_raw) - 1,
                source_artifact_bytes=len(generic_raw),
            )
            with self.assertRaisesRegex(ValidationError, "JSON byte limit"):
                independent_validator._load_source_artifact(
                    generic_path, "source artifact", generic_limits
                )
            with self.assertRaisesRegex(ValidationError, "byte limit"):
                independent_validator.load(
                    generic_path, "generic chain artifact", limits=generic_limits
                )
            with self.assertRaisesRegex(ValidationError, "byte limit"):
                independent_validator._load_source_artifact(
                    genbank_path,
                    "source artifact",
                    replace(genbank_limits, source_artifact_bytes=len(genbank_raw) - 1),
                )

        with self.assertRaisesRegex(ValidationError, "raw sequence source"):
            validate_source_artifact(
                genbank,
                raw_source=U49845.read_bytes(),
                limits=replace(DEFAULT_LIMITS, input_bytes=U49845.stat().st_size - 1),
            )

    def test_sparse_coherently_rehashed_auditor_forgery_is_rejected(self) -> None:
        forged = {
            "format": "brain01.sequence-ir",
            "version": 2,
            "compiler": {},
            "inputs": {},
            "sequence_ir": {
                "record_id": "x",
                "sequence": "A",
                "surprise": "accepted",
            },
            "source_map": {},
        }
        forged["ir_sha256"] = digest(forged["sequence_ir"])
        forged = _seal(forged)
        self.assertSourceRejected(forged)

    def test_fully_shaped_resealed_semantic_forgeries_are_rejected(self) -> None:
        attacks = []

        wrong_compiler = deepcopy(self.sequence)
        wrong_compiler["compiler"]["version"] = "999.0.0"
        attacks.append(_reseal_source(wrong_compiler))

        wrong_derived_digest = deepcopy(self.sequence)
        wrong_derived_digest["sequence_ir"]["sequence_sha256"] = "0" * 64
        attacks.append(_reseal_source(wrong_derived_digest))

        wrong_refget = deepcopy(self.sequence)
        wrong_refget["sequence_ir"]["refget_id"] = "SQ." + "A" * 32
        attacks.append(_reseal_source(wrong_refget))

        incomplete_map = deepcopy(self.sequence)
        incomplete_map["source_map"]["sequence_segments"][-1]["length"] -= 1
        attacks.append(_reseal_source(incomplete_map))

        unbound_context = deepcopy(self.sequence)
        unbound_context["sequence_ir"]["provenance"] = [
            {
                "id": "source-1",
                "kind": "database",
                "uri": "urn:example:source",
                "version": "1",
            }
        ]
        attacks.append(_reseal_source(unbound_context))

        for index, attack in enumerate(attacks):
            with self.subTest(attack=index):
                self.assertSourceRejected(attack)

    def test_collection_rejects_resealed_nested_global_and_refget_tampering(self) -> None:
        nested_digest = deepcopy(self.collection)
        nested = nested_digest["members"][0]["sequence_artifact"]
        nested["sequence_ir"]["canonical_sha256"] = "0" * 64
        nested["ir_sha256"] = digest(nested["sequence_ir"])
        nested = _seal(nested)
        nested_digest["members"][0]["sequence_artifact"] = nested
        nested_digest["members"][0]["ir_sha256"] = nested["ir_sha256"]
        nested_digest["members"][0]["artifact_sha256"] = nested["artifact_sha256"]
        nested_digest["collection_ir"]["members"][0]["ir_sha256"] = nested[
            "ir_sha256"
        ]
        nested_digest = _reseal_source(nested_digest)

        overlapping_map = deepcopy(self.collection)
        overlapping_map["members"][1]["input_source_map"]["sequence_segments"][0][
            "line"
        ] = 3
        overlapping_map = _reseal_source(overlapping_map)

        wrong_refget = deepcopy(self.collection)
        wrong_refget["refget_seqcol"]["level_2"]["lengths"][0] += 1
        wrong_refget = _reseal_source(wrong_refget)

        wrong_member_binding = deepcopy(self.collection)
        wrong_member_binding["members"][0]["record_id"] = "other"
        wrong_member_binding = _reseal_source(wrong_member_binding)

        for name, attack in (
            ("nested-derived-digest", nested_digest),
            ("global-map-overlap", overlapping_map),
            ("refget", wrong_refget),
            ("member-binding", wrong_member_binding),
        ):
            with self.subTest(attack=name):
                self.assertSourceRejected(attack)

    def test_raw_iupac_structure_and_full_raw_replay(self) -> None:
        raw = b"acgTRYSWKMBDHVN"
        for supplied in (raw, gzip.compress(raw, mtime=0)):
            with self.subTest(gzip=supplied is not raw):
                source = SequenceCollectionCompiler().compile_raw(
                    supplied, "raw-1"
                ).to_dict()
                validate_source_artifact(source)
                validate_source_artifact(
                    source, raw_source=supplied, record_id="raw-1"
                )
        source = SequenceCollectionCompiler().compile_raw(raw, "raw-1").to_dict()
        with self.assertRaises(ValidationError):
            validate_source_artifact(source, raw_source=raw + b"A", record_id="raw-1")

        wrong_logical = deepcopy(source)
        wrong_logical["inputs"]["logical_sha256"] = hashlib.sha256(b"other").hexdigest()
        self.assertSourceRejected(_reseal_source(wrong_logical))

    def test_single_sequence_full_replay_with_context(self) -> None:
        fasta = b">contextual\nACGTAC\n"
        context = json.dumps(
            {
                "format": "brain01.sequence-context",
                "version": 1,
                "record_id": "contextual",
                "reference": {
                    "assembly": "assembly-1",
                    "contig": "contig-1",
                    "start": 10,
                    "end": 16,
                    "coordinate_system": "0-based-half-open",
                    "orientation": "forward",
                    "aliases": ["alias-1"],
                },
                "provenance": [
                    {
                        "id": "reference-source",
                        "kind": "database",
                        "uri": "urn:example:reference",
                        "version": "1",
                    }
                ],
            },
            separators=(",", ":"),
        ).encode("utf-8")
        source = SequenceCompiler().compile_bytes(fasta, context).to_dict()
        validate_source_artifact(source, raw_source=fasta, context_raw=context)
        with self.assertRaises(ValidationError):
            validate_source_artifact(
                source, raw_source=fasta.replace(b"ACGTAC", b"ACGTAA"), context_raw=context
            )

    def test_source_record_ceiling_is_enforced_before_axis_use(self) -> None:
        limits = replace(DEFAULT_LIMITS, source_records=1)
        with self.assertRaises(ValidationError):
            source_record_lengths(self.collection, limits=limits)
        with self.assertRaises(ValidationError):
            source_record_lengths(
                self.sequence, limits=replace(DEFAULT_LIMITS, source_records=0)
            )
        with self.assertRaises(ValidationError):
            source_record_lengths(
                self.sequence, limits=replace(DEFAULT_LIMITS, decompressed_bytes=14)
            )

    def _response(self, record_id: str, start: int) -> dict:
        raw = b"\x01\x00"
        return _seal(
            {
                "format": "brainc.prediction-response",
                "version": 2,
                "request_artifact_sha256": "1" * 64,
                "provider": {"name": "source-test", "version": "1"},
                "model_identity": {"kind": "content-sha256", "value": "2" * 64},
                "outputs": [
                    {
                        "id": "coordinate.values",
                        "type": {"dtype": "bool", "shape": [2]},
                        "unit": None,
                        "axes": [
                            {
                                "kind": "source-coordinate",
                                "source_artifact_sha256": self.collection[
                                    "artifact_sha256"
                                ],
                                "record_id": record_id,
                                "start": start,
                                "step": 1,
                                "span": 1,
                                "coordinate_system": "0-based-half-open",
                            }
                        ],
                        "storage": {
                            "kind": "inline-base64",
                            "data": base64.b64encode(raw).decode("ascii"),
                            "byte_length": 2,
                            "sha256": hashlib.sha256(raw).hexdigest(),
                        },
                    }
                ],
            }
        )

    def test_response_axes_use_validated_record_ids_and_exact_lengths(self) -> None:
        records = source_record_lengths(self.collection)
        validate_response_artifact(
            self._response("alpha", 4),
            sequence_sha256=self.collection["artifact_sha256"],
            records=records,
        )
        for record_id, start in (("unknown", 0), ("alpha", 5)):
            with self.subTest(record_id=record_id, start=start):
                with self.assertRaises(ValidationError):
                    validate_response_artifact(
                        self._response(record_id, start),
                        sequence_sha256=self.collection["artifact_sha256"],
                        records=records,
                    )

    def test_external_blob_fails_closed_without_posix_descriptor_guarantees(self) -> None:
        response = self._response("alpha", 4)
        raw = b"\x01\x00"
        response["outputs"][0]["storage"] = {
            "kind": "sha256-blob",
            "byte_length": len(raw),
            "sha256": hashlib.sha256(raw).hexdigest(),
        }
        response = _seal(response)
        with mock.patch.object(independent_validator.os, "name", "nt"):
            with self.assertRaisesRegex(ValidationError, "requires a POSIX host"):
                validate_response_artifact(
                    response,
                    sequence_sha256=self.collection["artifact_sha256"],
                    records=source_record_lengths(self.collection),
                    blob_root="unused",
                )

    def test_cli_check_response_passes_validated_record_lengths(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source_path = root / "source.json"
            response_path = root / "response.json"
            source_path.write_text(json.dumps(self.collection), encoding="utf-8")
            response_path.write_text(
                json.dumps(self._response("alpha", 5)), encoding="utf-8"
            )
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "brainc",
                    "check-v2",
                    "response",
                    str(response_path),
                    "--source",
                    str(source_path),
                ],
                cwd=ROOT,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
            self.assertIn("coordinate exceeds", result.stderr)

    def test_public_module_validator_rejects_spoofed_internal_target_metadata(self) -> None:
        target = target_artifact(
            {
                "id": "org.example.source-test",
                "abi_major": 1,
                "opsets": [{"domain": "org.example.noop", "version": 1}],
                "unit_schemas": [
                    {
                        "id": "node",
                        "fields": [
                            {
                                "id": "state",
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
        target_binding = {
            "id": target["contract"]["id"],
            "abi_major": target["contract"]["abi_major"],
            "contract_sha256": target["contract_sha256"],
        }
        body = {
            "target": target_binding,
            "requirements": [],
            "budgets": {
                "operations": 0,
                "tensor_bytes": 0,
                "units": 0,
                "edges": 0,
                "attachments": 0,
            },
            "tensors": [],
            "entrypoint": {"name": "develop", "operations": []},
        }
        module = _seal(
            {
                "format": "brainc.development-module",
                "version": 1,
                "producer": {"name": "source-test", "version": "1", "passes": []},
                "sources": {
                    "sequence": {
                        "format": self.collection["format"],
                        "version": self.collection["version"],
                        "artifact_sha256": self.collection["artifact_sha256"],
                        "ir_sha256": self.collection["collection_ir_sha256"],
                    },
                    "provider_manifest": {"artifact_sha256": "1" * 64},
                    "prediction_request": {"artifact_sha256": "2" * 64},
                    "prediction_response": {"artifact_sha256": "3" * 64},
                    "lowering_policy": {"artifact_sha256": "4" * 64},
                    "target_contract": {
                        "artifact_sha256": target["artifact_sha256"]
                    },
                },
                "module": body,
                "module_sha256": digest(body),
            }
        )
        validate_module_artifact(module, target)
        internal_target_info = validate_target_artifact(target)
        with self.assertRaises(ValidationError):
            validate_module_artifact(module, internal_target_info)

    def test_report_save_is_atomic_and_never_follows_existing_symlink(self) -> None:
        report = {"format": "brainc.validation-report", "valid": True}
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            victim = root / "victim.json"
            victim.write_bytes(b"do-not-touch")
            linked_report = root / "linked-report.json"
            linked_report.symlink_to(victim)
            with self.assertRaises(ValidationError):
                save_report(report, linked_report)
            self.assertTrue(linked_report.is_symlink())
            self.assertEqual(victim.read_bytes(), b"do-not-touch")

            report_path = root / "report.json"
            report_path.write_bytes(b"previous-complete-report")
            with mock.patch.object(
                independent_validator.os,
                "fsync",
                side_effect=OSError("injected incomplete write"),
            ):
                with self.assertRaises(ValidationError):
                    save_report(report, report_path)
            self.assertEqual(report_path.read_bytes(), b"previous-complete-report")
            self.assertEqual(list(root.glob(f".{report_path.name}.*")), [])

            save_report(report, report_path)
            self.assertEqual(json.loads(report_path.read_text(encoding="utf-8")), report)

    def test_independent_file_reader_rejects_same_size_in_place_change(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "source.json"
            path.write_text(json.dumps(self.collection), encoding="utf-8")
            real_fstat = independent_validator.os.fstat
            calls = 0

            def changing_fstat(descriptor: int):
                nonlocal calls
                calls += 1
                value = real_fstat(descriptor)
                if calls != 2:
                    return value

                class Changed:
                    st_mode = value.st_mode
                    st_dev = value.st_dev
                    st_ino = value.st_ino
                    st_size = value.st_size
                    st_mtime_ns = value.st_mtime_ns + 1
                    st_ctime_ns = value.st_ctime_ns

                return Changed()

            with mock.patch.object(
                independent_validator.os, "fstat", side_effect=changing_fstat
            ):
                with self.assertRaisesRegex(ValidationError, "changed while being read"):
                    independent_validator.load(path, "sequence source")

    def test_genbank_validation_import_isolation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            artifact_path = _write_json(root / "source.json", self._genbank())
            command = (
                "import json,sys; "
                f"sys.path.insert(0,{str(ROOT)!r}); "
                "from brainc.validator_v2 import validate_source_artifact; "
                f"artifact=json.loads(open({str(artifact_path)!r},encoding='utf-8').read()); "
                f"raw=open({str(U49845)!r},'rb').read(); "
                "validate_source_artifact(artifact,raw_source=raw); "
                "forbidden={'brainc.insdc','brainc.insdc_location','brainc.compiler',"
                "'brainc.provider','brainc.sequence','brainc.sequence_collection',"
                "'brainc.sequence_collection_v2','brainc.v2','brainc._canonical',"
                "'brainc._io'}; "
                "loaded=sorted(forbidden.intersection(sys.modules)); "
                "assert not loaded,loaded"
            )
            subprocess.run(
                [sys.executable, "-I", "-c", command],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
            )


if __name__ == "__main__":
    unittest.main()

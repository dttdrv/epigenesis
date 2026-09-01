from __future__ import annotations

import ast
from contextlib import redirect_stdout
from copy import deepcopy
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

from brainc._canonical import digest
from brainc.development_bundle import compile_development
from brainc.source import (
    FASTA_PROFILE,
    GFF3_PROFILE,
    GENBANK_PROFILE,
    RAW_PROFILE,
    REFERENCE_FASTA_PROFILE,
    compile_source,
)
from brainc.v2 import (
    inline_storage,
    make_development_request,
    pack,
    policy_artifact,
    save,
    target_artifact,
)
from brainc.v2.target import DEV_DOMAIN, OP_UNIT_CREATE
import brainc.validator_development as independent
from brainc.validator_development import (
    DevelopmentValidationError,
    load_development_bundle_directory,
    load_source_bundle_directory,
    main,
    validate_development_bundle,
    validate_development_paths,
    validate_development_report,
)


ROOT = Path(__file__).resolve().parents[1]
FASTA = ROOT / "tests/data/J02482.1.fasta"
GENBANK = ROOT / "tests/data/U49845.1.gb"
GFF3 = ROOT / "tests/data/J02482.1.gff3"


def _seal(core: dict) -> dict:
    return {**core, "artifact_sha256": digest(core)}


def _reseal(value: dict) -> dict:
    return _seal({key: item for key, item in value.items() if key != "artifact_sha256"})


def _source_reference(descriptor: dict) -> dict:
    return {
        "format": descriptor["format"],
        "version": descriptor["version"],
        "artifact_sha256": descriptor["artifact_sha256"],
        "ir_sha256": descriptor["source_ir_sha256"],
    }


def _child_reference(value: dict) -> dict:
    return {
        "format": value["format"],
        "version": value["version"],
        "artifact_sha256": value["artifact_sha256"],
    }


class _Fixture:
    def __init__(
        self,
        root: Path,
        profile: str = FASTA_PROFILE,
        inputs: dict[str, Path] | None = None,
        parameters: dict | None = None,
        *,
        external_blob: bool = False,
    ) -> None:
        self.root = root
        self.root.mkdir(parents=True, exist_ok=True)
        self.source_inputs = {"sequence": FASTA} if inputs is None else inputs
        self.source = compile_source(
            profile,
            self.source_inputs,
            parameters=(
                {"wrapper": "identity"} if parameters is None else parameters
            ),
        )
        self.source_directory = root / "source"
        self.source.save(self.source_directory)
        self.source_path = self.source_directory / "source.json"

        self.target = target_artifact(
            {
                "id": "org.example.development",
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
        self.target_path = self._write("target", self.target)
        raw = pack("u64", [1])
        self.blob_root: Path | None = None
        if external_blob:
            blob_sha = hashlib.sha256(raw).hexdigest()
            self.blob_root = root / "blobs"
            self.blob_root.mkdir()
            (self.blob_root / blob_sha).write_bytes(raw)
            storage = {
                "kind": "sha256-blob",
                "byte_length": len(raw),
                "sha256": blob_sha,
            }
        else:
            storage = inline_storage(raw)
        output = {
            "id": "t.value",
            "type": {"dtype": "u64", "shape": []},
            "unit": None,
            "axes": [],
            "storage": storage,
        }
        self.manifest = _seal(
            {
                "format": "brainc.provider-manifest",
                "version": 2,
                "provider": {"name": "org.example.provider", "version": "1"},
                "model_identity": {
                    "kind": "content-sha256",
                    "value": digest({"algorithm": "caller-supplied", "version": 1}),
                },
                "accepts": [f"brainc.source-descriptor/v1;profile={profile}"],
                "outputs": [
                    {
                        key: deepcopy(output[key])
                        for key in ("id", "type", "unit", "axes")
                    }
                ],
            }
        )
        self.manifest_path = self._write("manifest", self.manifest)
        self.request = make_development_request(
            self.source_directory,
            self.manifest_path,
            ["t.value"],
        )
        self.request_path = self._write("request", self.request)
        self.response = _seal(
            {
                "format": "brainc.prediction-response",
                "version": 2,
                "request_artifact_sha256": self.request["artifact_sha256"],
                "provider": deepcopy(self.manifest["provider"]),
                "model_identity": deepcopy(self.manifest["model_identity"]),
                "outputs": [output],
            }
        )
        self.response_path = self._write("response", self.response)
        target_binding = {
            "id": self.target["contract"]["id"],
            "abi_major": self.target["contract"]["abi_major"],
            "contract_sha256": self.target["contract_sha256"],
        }
        self.policy = policy_artifact(
            "org.example.lowering",
            target_binding,
            [{"id": "t.value", "from_output": "t.value"}],
            [
                {
                    "id": "create.node",
                    "op": OP_UNIT_CREATE,
                    "version": 1,
                    "schema": "node",
                    "count": "t.value",
                    "initializers": [{"field": "value", "tensor": "t.value"}],
                }
            ],
        )
        self.policy_path = self._write("policy", self.policy)
        self.compiled = compile_development(
            self.source_directory,
            self.manifest_path,
            self.request_path,
            self.response_path,
            self.policy_path,
            self.target_path,
            blob_root=self.blob_root,
        )
        self.development_directory = root / "development"
        self.compiled.save(self.development_directory)

    def _write(self, name: str, value: dict) -> Path:
        path = self.root / f"{name}.json"
        save(value, path)
        return path

    @property
    def source_bytes(self) -> dict[str, bytes]:
        return {role: path.read_bytes() for role, path in self.source_inputs.items()}

    def validate_paths(self) -> dict:
        return validate_development_paths(
            self.source_directory,
            self.development_directory,
            source_inputs=self.source_inputs,
            blob_root=self.blob_root,
        )


class DevelopmentValidatorTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.fixture = _Fixture(self.root)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_full_path_and_in_memory_replay_are_deterministic_and_sealed(self) -> None:
        first = self.fixture.validate_paths()
        second = self.fixture.validate_paths()
        self.assertEqual(first, second)
        self.assertIs(validate_development_report(first), first)
        self.assertTrue(first["valid"])
        self.assertEqual(first["result"]["budgets"]["units"], 1)
        self.assertEqual(first["result"]["blobs"], [])

        bundle, artifacts = load_development_bundle_directory(
            self.fixture.development_directory
        )
        descriptor, natives = load_source_bundle_directory(
            self.fixture.source_directory
        )
        direct = validate_development_bundle(
            bundle,
            artifacts,
            descriptor,
            natives,
            source_inputs=self.fixture.source_bytes,
        )
        self.assertEqual(direct, first)

    def test_every_source_profile_replays_from_original_sources(self) -> None:
        multi = self.root / "multi.fasta"
        multi.write_bytes(b">a\nACGT\n>b\nNNRY\n")
        compressed = self.root / "multi.fasta.gz"
        compressed.write_bytes(gzip.compress(multi.read_bytes(), mtime=0))
        phix_compressed = self.root / "phix.fasta.gz"
        phix_compressed.write_bytes(gzip.compress(FASTA.read_bytes(), mtime=0))
        raw = self.root / "raw.dna"
        raw.write_bytes(b"ACGTRYSWKMBDHVN")
        phix_raw = self.root / "phix.dna"
        phix_raw.write_bytes(
            b"".join(
                line.strip()
                for line in FASTA.read_bytes().splitlines()
                if not line.startswith(b">")
            )
        )
        cases = (
            (RAW_PROFILE, {"sequence": raw}, {"record_id": "raw-1"}),
            (FASTA_PROFILE, {"sequence": multi}, {"wrapper": "identity"}),
            (FASTA_PROFILE, {"sequence": compressed}, {"wrapper": "gzip"}),
            (
                REFERENCE_FASTA_PROFILE,
                {"sequence": multi},
                {"wrapper": "identity"},
            ),
            (
                REFERENCE_FASTA_PROFILE,
                {"sequence": compressed},
                {"wrapper": "gzip"},
            ),
            (GENBANK_PROFILE, {"genbank": GENBANK}, {}),
            (
                GFF3_PROFILE,
                {"sequence": FASTA, "annotation": GFF3},
                {
                    "sequence_profile": FASTA_PROFILE,
                    "sequence_parameters": {"wrapper": "identity"},
                },
            ),
            (
                GFF3_PROFILE,
                {"sequence": phix_compressed, "annotation": GFF3},
                {
                    "sequence_profile": FASTA_PROFILE,
                    "sequence_parameters": {"wrapper": "gzip"},
                },
            ),
            (
                GFF3_PROFILE,
                {"sequence": phix_raw, "annotation": GFF3},
                {
                    "sequence_profile": RAW_PROFILE,
                    "sequence_parameters": {"record_id": "J02482.1"},
                },
            ),
        )
        for index, (profile, inputs, parameters) in enumerate(cases):
            with self.subTest(profile=profile, parameters=parameters):
                fixture = _Fixture(
                    self.root / f"case-{index}",
                    profile,
                    inputs,
                    parameters,
                )
                self.assertTrue(fixture.validate_paths()["valid"])

    def test_reference_fasta_replay_streams_original_path_and_fails_closed(self) -> None:
        reference = self.root / "reference.fasta"
        reference.write_bytes(b">chr1\nACGTACGT\n>chr2\nNNRY\n")
        fixture = _Fixture(
            self.root / "reference-route",
            REFERENCE_FASTA_PROFILE,
            {"sequence": reference},
            {"wrapper": "identity"},
        )
        with mock.patch.object(
            independent,
            "_read_regular",
            side_effect=AssertionError("reference FASTA must stream by path"),
        ):
            report = fixture.validate_paths()
        self.assertEqual(
            report["inputs"]["source_inputs"]["sequence"],
            fixture.source.artifacts["reference"]["inputs"]["source"],
        )

        stdout = io.StringIO()
        with redirect_stdout(stdout):
            status = main(
                [
                    str(fixture.source_directory),
                    str(fixture.development_directory),
                    "--source-input",
                    f"sequence={reference}",
                ]
            )
        self.assertEqual(status, 0)
        cli_report = json.loads(stdout.getvalue())
        self.assertIs(validate_development_report(cli_report), cli_report)

        with self.assertRaisesRegex(DevelopmentValidationError, "must be a path"):
            validate_development_bundle(
                fixture.compiled.to_dict(),
                fixture.compiled.artifacts,
                fixture.source.to_dict(),
                fixture.source.artifacts,
                source_inputs={"sequence": reference.read_bytes()},
            )

        wrong = self.root / "wrong-reference.fasta"
        wrong.write_bytes(b">chr1\nACGT\n")
        with self.assertRaisesRegex(DevelopmentValidationError, "reference FASTA replay"):
            validate_development_paths(
                fixture.source_directory,
                fixture.development_directory,
                source_inputs={"sequence": wrong},
            )

        alias = self.root / "reference-link"
        alias.symlink_to(reference)
        with self.assertRaisesRegex(DevelopmentValidationError, "non-linked"):
            validate_development_paths(
                fixture.source_directory,
                fixture.development_directory,
                source_inputs={"sequence": alias},
            )

        hardlink = self.root / "reference-hardlink"
        hardlink.hardlink_to(reference)
        with self.assertRaisesRegex(DevelopmentValidationError, "non-linked"):
            validate_development_paths(
                fixture.source_directory,
                fixture.development_directory,
                source_inputs={"sequence": hardlink},
            )
        hardlink.unlink()

        descriptor = deepcopy(fixture.source.to_dict())
        natives = deepcopy(fixture.source.artifacts)
        native = natives["reference"]
        native["inputs"]["source"]["sha256"] = "0" * 64
        native["inputs"]["logical"]["sha256"] = "0" * 64
        natives["reference"] = _reseal(native)
        descriptor["source_ir"]["inputs"]["sequence"] = deepcopy(
            natives["reference"]["inputs"]["source"]
        )
        descriptor["source_ir"]["artifacts"]["reference"][
            "artifact_sha256"
        ] = natives["reference"]["artifact_sha256"]
        descriptor["source_ir_sha256"] = digest(descriptor["source_ir"])
        descriptor = _reseal(descriptor)
        with self.assertRaisesRegex(DevelopmentValidationError, "reference FASTA replay"):
            validate_development_bundle(
                fixture.compiled.to_dict(),
                fixture.compiled.artifacts,
                descriptor,
                natives,
                source_inputs={"sequence": reference},
            )

    def test_path_loaders_require_exact_flat_nonlinked_closure_and_bounds(self) -> None:
        extra = self.fixture.development_directory / "extra.json"
        extra.write_text("{}", encoding="utf-8")
        with self.assertRaisesRegex(DevelopmentValidationError, "not closed"):
            self.fixture.validate_paths()
        extra.unlink()

        child = self.fixture.development_directory / "target_contract.json"
        raw = child.read_bytes()
        child.unlink()
        with self.assertRaisesRegex(DevelopmentValidationError, "not closed"):
            self.fixture.validate_paths()
        child.write_bytes(raw)

        outside = self.root / "outside-target.json"
        outside.write_bytes(raw)
        child.unlink()
        child.symlink_to(outside)
        with self.assertRaisesRegex(DevelopmentValidationError, "non-linked"):
            self.fixture.validate_paths()
        child.unlink()
        child.write_bytes(raw)

        source_extra = self.fixture.source_directory / "extra.json"
        source_extra.write_text("{}", encoding="utf-8")
        with self.assertRaisesRegex(DevelopmentValidationError, "not closed"):
            self.fixture.validate_paths()
        source_extra.unlink()

        with mock.patch.object(independent, "MAX_JSON_BYTES", 1):
            with self.assertRaisesRegex(DevelopmentValidationError, "exceeds 1 bytes"):
                load_development_bundle_directory(self.fixture.development_directory)

        source_alias = self.root / "source-alias"
        source_alias.symlink_to(self.fixture.source_directory, target_is_directory=True)
        with self.assertRaisesRegex(DevelopmentValidationError, "non-linked directory"):
            validate_development_paths(
                source_alias,
                self.fixture.development_directory,
                source_inputs=self.fixture.source_inputs,
            )

    def test_closed_schemas_and_exact_json_types_reject_coherent_seals(self) -> None:
        bundle = deepcopy(self.fixture.compiled.to_dict())
        artifacts = deepcopy(self.fixture.compiled.artifacts)
        descriptor = self.fixture.source.to_dict()
        natives = self.fixture.source.artifacts

        wrong_type = deepcopy(bundle)
        wrong_type["version"] = 1.0
        wrong_type = _reseal(wrong_type)
        with self.assertRaisesRegex(DevelopmentValidationError, "identity"):
            validate_development_bundle(
                wrong_type,
                artifacts,
                descriptor,
                natives,
                source_inputs=self.fixture.source_bytes,
            )

        extra = deepcopy(bundle)
        extra["surprise"] = True
        extra = _reseal(extra)
        with self.assertRaisesRegex(DevelopmentValidationError, "unknown"):
            validate_development_bundle(
                extra,
                artifacts,
                descriptor,
                natives,
                source_inputs=self.fixture.source_bytes,
            )

        non_string_key = deepcopy(bundle)
        non_string_key[1] = "not-json"
        with self.assertRaisesRegex(DevelopmentValidationError, "non-string object key"):
            validate_development_bundle(
                non_string_key,
                artifacts,
                descriptor,
                natives,
                source_inputs=self.fixture.source_bytes,
            )

        descriptor_extra = deepcopy(descriptor)
        descriptor_extra["source_ir"]["inputs"]["sequence"]["byte_length"] = True
        descriptor_extra["source_ir_sha256"] = digest(descriptor_extra["source_ir"])
        descriptor_extra = _reseal(descriptor_extra)
        with self.assertRaisesRegex(DevelopmentValidationError, "input digests"):
            validate_development_bundle(
                bundle,
                artifacts,
                descriptor_extra,
                natives,
                source_inputs=self.fixture.source_bytes,
            )

        unhashable_route = deepcopy(descriptor)
        unhashable_route["source_ir"]["profile"] = GFF3_PROFILE
        unhashable_route["source_ir"]["parameters"] = {
            "sequence_profile": [],
            "sequence_parameters": {},
        }
        with self.assertRaisesRegex(DevelopmentValidationError, "sequence profile"):
            validate_development_bundle(
                bundle,
                artifacts,
                unhashable_route,
                natives,
                source_inputs=self.fixture.source_bytes,
            )

        with self.assertRaisesRegex(DevelopmentValidationError, "must be an object"):
            validate_development_bundle(
                bundle,
                artifacts,
                descriptor,
                {"sequence": True},
                source_inputs=self.fixture.source_bytes,
            )

    def test_coherently_resealed_compilation_and_module_forgeries_fail(self) -> None:
        descriptor = self.fixture.source.to_dict()
        natives = self.fixture.source.artifacts
        original = self.fixture.source_bytes

        artifacts = deepcopy(self.fixture.compiled.artifacts)
        bundle = deepcopy(self.fixture.compiled.to_dict())
        record = artifacts["compilation_record"]
        record["result"]["budgets"]["units"] += 1
        artifacts["compilation_record"] = _reseal(record)
        bundle["artifacts"]["compilation_record"] = _child_reference(
            artifacts["compilation_record"]
        )
        bundle = _reseal(bundle)
        with self.assertRaisesRegex(DevelopmentValidationError, "compilation record differs"):
            validate_development_bundle(
                bundle,
                artifacts,
                descriptor,
                natives,
                source_inputs=original,
            )

        artifacts = deepcopy(self.fixture.compiled.artifacts)
        bundle = deepcopy(self.fixture.compiled.to_dict())
        module = artifacts["development_module"]
        module["module"]["budgets"]["units"] += 1
        module["module_sha256"] = digest(module["module"])
        artifacts["development_module"] = _reseal(module)
        record = artifacts["compilation_record"]
        record["artifacts"]["development_module"] = _child_reference(
            artifacts["development_module"]
        )
        record["result"]["module_sha256"] = module["module_sha256"]
        record["result"]["budgets"] = deepcopy(module["module"]["budgets"])
        artifacts["compilation_record"] = _reseal(record)
        for role in ("development_module", "compilation_record"):
            bundle["artifacts"][role] = _child_reference(artifacts[role])
        bundle = _reseal(bundle)
        with self.assertRaisesRegex(DevelopmentValidationError, "v2 compilation replay"):
            validate_development_bundle(
                bundle,
                artifacts,
                descriptor,
                natives,
                source_inputs=original,
            )

    def test_descriptor_provider_acceptance_is_exactly_profile_qualified(self) -> None:
        descriptor = self.fixture.source.to_dict()
        for acceptance in (
            "brainc.source-descriptor/v1",
            f"brainc.source-descriptor/v1;profile={RAW_PROFILE}",
        ):
            with self.subTest(acceptance=acceptance):
                artifacts = deepcopy(self.fixture.compiled.artifacts)
                manifest = artifacts["provider_manifest"]
                manifest["accepts"] = [acceptance]
                artifacts["provider_manifest"] = _reseal(manifest)
                request = artifacts["prediction_request"]
                request["provider_manifest_sha256"] = artifacts[
                    "provider_manifest"
                ]["artifact_sha256"]
                artifacts["prediction_request"] = _reseal(request)
                response = artifacts["prediction_response"]
                response["request_artifact_sha256"] = artifacts[
                    "prediction_request"
                ]["artifact_sha256"]
                artifacts["prediction_response"] = _reseal(response)
                module = artifacts["development_module"]
                for role in (
                    "provider_manifest",
                    "prediction_request",
                    "prediction_response",
                ):
                    module["sources"][role]["artifact_sha256"] = artifacts[role][
                        "artifact_sha256"
                    ]
                artifacts["development_module"] = _reseal(module)
                record = artifacts["compilation_record"]
                for role in (
                    "provider_manifest",
                    "prediction_request",
                    "prediction_response",
                    "development_module",
                ):
                    record["artifacts"][role] = _child_reference(artifacts[role])
                artifacts["compilation_record"] = _reseal(record)
                bundle = deepcopy(self.fixture.compiled.to_dict())
                for role in (
                    "provider_manifest",
                    "prediction_request",
                    "prediction_response",
                    "development_module",
                    "compilation_record",
                ):
                    bundle["artifacts"][role] = _child_reference(artifacts[role])
                bundle = _reseal(bundle)

                with self.assertRaisesRegex(
                    DevelopmentValidationError,
                    "does not declare source format support",
                ):
                    validate_development_bundle(
                        bundle,
                        artifacts,
                        descriptor,
                        self.fixture.source.artifacts,
                        source_inputs=self.fixture.source_bytes,
                    )

    def test_fully_coherent_source_digest_forgery_is_bound_to_original_bytes(self) -> None:
        length_forgery = deepcopy(self.fixture.source.to_dict())
        length_forgery["source_ir"]["inputs"]["sequence"]["byte_length"] += 1
        length_forgery["source_ir_sha256"] = digest(length_forgery["source_ir"])
        length_forgery = _reseal(length_forgery)
        with self.assertRaisesRegex(DevelopmentValidationError, "original bytes"):
            validate_development_bundle(
                self.fixture.compiled.to_dict(),
                self.fixture.compiled.artifacts,
                length_forgery,
                self.fixture.source.artifacts,
                source_inputs=self.fixture.source_bytes,
            )

        descriptor = deepcopy(self.fixture.source.to_dict())
        descriptor["source_ir"]["inputs"]["sequence"]["sha256"] = "0" * 64
        descriptor["source_ir_sha256"] = digest(descriptor["source_ir"])
        descriptor = _reseal(descriptor)
        source_reference = _source_reference(descriptor)

        artifacts = deepcopy(self.fixture.compiled.artifacts)
        request = artifacts["prediction_request"]
        request["source"] = deepcopy(source_reference)
        artifacts["prediction_request"] = _reseal(request)
        response = artifacts["prediction_response"]
        response["request_artifact_sha256"] = artifacts["prediction_request"][
            "artifact_sha256"
        ]
        artifacts["prediction_response"] = _reseal(response)
        module = artifacts["development_module"]
        module["sources"]["sequence"] = deepcopy(source_reference)
        module["sources"]["prediction_request"]["artifact_sha256"] = artifacts[
            "prediction_request"
        ]["artifact_sha256"]
        module["sources"]["prediction_response"]["artifact_sha256"] = artifacts[
            "prediction_response"
        ]["artifact_sha256"]
        artifacts["development_module"] = _reseal(module)

        record = artifacts["compilation_record"]
        record["source"] = deepcopy(source_reference)
        for role in (
            "prediction_request",
            "prediction_response",
            "development_module",
        ):
            record["artifacts"][role] = _child_reference(artifacts[role])
        artifacts["compilation_record"] = _reseal(record)
        bundle = deepcopy(self.fixture.compiled.to_dict())
        bundle["source"]["descriptor"] = deepcopy(source_reference)
        for role in (
            "prediction_request",
            "prediction_response",
            "development_module",
            "compilation_record",
        ):
            bundle["artifacts"][role] = _child_reference(artifacts[role])
        bundle = _reseal(bundle)

        with self.assertRaisesRegex(DevelopmentValidationError, "original bytes"):
            validate_development_bundle(
                bundle,
                artifacts,
                descriptor,
                self.fixture.source.artifacts,
                source_inputs=self.fixture.source_bytes,
            )

    def test_external_blob_bytes_and_bundle_list_are_both_validated(self) -> None:
        fixture = _Fixture(self.root / "external", external_blob=True)
        with mock.patch.object(
            independent.validator_v2,
            "_read_blob",
            side_effect=AssertionError("v2 must consume the prevalidated cache"),
        ):
            report = fixture.validate_paths()
        self.assertEqual(len(report["result"]["blobs"]), 1)
        blob = report["result"]["blobs"][0]
        assert fixture.blob_root is not None
        blob_path = fixture.blob_root / blob["sha256"]
        raw = blob_path.read_bytes()

        blob_path.unlink()
        with self.assertRaisesRegex(DevelopmentValidationError, "blob"):
            fixture.validate_paths()
        blob_path.write_bytes(raw)

        target = self.root / "blob-outside"
        target.write_bytes(raw)
        blob_path.unlink()
        blob_path.symlink_to(target)
        with self.assertRaisesRegex(DevelopmentValidationError, "blob"):
            fixture.validate_paths()
        blob_path.unlink()
        blob_path.write_bytes(raw)

        hardlink = fixture.blob_root / "hardlink"
        hardlink.hardlink_to(blob_path)
        with self.assertRaisesRegex(DevelopmentValidationError, "regular non-linked"):
            fixture.validate_paths()
        hardlink.unlink()

        bundle = deepcopy(fixture.compiled.to_dict())
        bundle["blobs"] = []
        bundle = _reseal(bundle)
        with self.assertRaisesRegex(DevelopmentValidationError, "bundle differs"):
            validate_development_bundle(
                bundle,
                fixture.compiled.artifacts,
                fixture.source.to_dict(),
                fixture.source.artifacts,
                source_inputs=fixture.source_bytes,
                blob_root=fixture.blob_root,
            )

    def test_external_blob_unique_count_and_cumulative_byte_boundaries(self) -> None:
        fixture = _Fixture(self.root / "blob-boundaries", external_blob=True)
        with (
            mock.patch.object(independent, "MAX_EXTERNAL_BLOBS", 1),
            mock.patch.object(independent, "MAX_EXTERNAL_BLOB_BYTES", 8),
        ):
            self.assertTrue(fixture.validate_paths()["valid"])
        with mock.patch.object(independent, "MAX_EXTERNAL_BLOBS", 0):
            with self.assertRaisesRegex(DevelopmentValidationError, "count ceiling 0"):
                fixture.validate_paths()
        with mock.patch.object(independent, "MAX_EXTERNAL_BLOB_BYTES", 7):
            with self.assertRaisesRegex(
                DevelopmentValidationError,
                "cumulative declared byte ceiling 7",
            ):
                fixture.validate_paths()

    def test_standalone_cli_emits_and_optionally_saves_sealed_reports(self) -> None:
        output = self.root / "validation-report.json"
        arguments = [
            str(self.fixture.source_directory),
            str(self.fixture.development_directory),
            "--source-input",
            f"sequence={self.fixture.source_inputs['sequence']}",
            "--output",
            str(output),
        ]
        stdout = io.StringIO()
        with redirect_stdout(stdout):
            status = main(arguments)
        self.assertEqual(status, 0)
        report = json.loads(stdout.getvalue())
        self.assertEqual(report, json.loads(output.read_text(encoding="utf-8")))
        self.assertIs(validate_development_report(report), report)

        wrong = self.root / "wrong.fasta"
        wrong.write_bytes(b">wrong\nACGT\n")
        failed_output = self.root / "failed-report.json"
        failed_stdout = io.StringIO()
        with redirect_stdout(failed_stdout):
            status = main(
                [
                    str(self.fixture.source_directory),
                    str(self.fixture.development_directory),
                    "--source-input",
                    f"sequence={wrong}",
                    "-o",
                    str(failed_output),
                ]
            )
        self.assertEqual(status, 1)
        failure = json.loads(failed_stdout.getvalue())
        self.assertFalse(failure["valid"])
        self.assertEqual(failure["error"]["code"], "DEVVAL001")
        self.assertEqual(
            failure,
            json.loads(failed_output.read_text(encoding="utf-8")),
        )
        self.assertIs(validate_development_report(failure), failure)

    def test_original_input_roles_and_paths_are_exact_and_nonlinked(self) -> None:
        with self.assertRaisesRegex(DevelopmentValidationError, "keys invalid"):
            validate_development_paths(
                self.fixture.source_directory,
                self.fixture.development_directory,
                source_inputs={**self.fixture.source_inputs, "annotation": GFF3},
            )
        alias = self.root / "fasta-link"
        alias.symlink_to(FASTA)
        with self.assertRaisesRegex(DevelopmentValidationError, "non-linked"):
            validate_development_paths(
                self.fixture.source_directory,
                self.fixture.development_directory,
                source_inputs={"sequence": alias},
            )


class DevelopmentValidatorIsolationTests(unittest.TestCase):
    def test_import_loads_no_source_bundle_or_v2_producer_modules(self) -> None:
        command = (
            "import sys; "
            f"sys.path.insert(0,{str(ROOT)!r}); "
            "import brainc.validator_development as validator; "
            "assert callable(validator.main); "
            "forbidden={'brainc.source','brainc.development_bundle','brainc.bio',"
            "'brainc.insdc','brainc.sequence','brainc.sequence_collection',"
            "'brainc.source_scale',"
            "'brainc.sequence_collection_v2','brainc._canonical','brainc._io',"
            "'brainc.v2','brainc.v2.compiler','brainc.v2.provider',"
            "'brainc.v2.policy','brainc.v2.target','brainc.v2.tensor'}; "
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

    def test_source_contains_only_allowed_independent_relative_imports(self) -> None:
        path = ROOT / "brainc/validator_development.py"
        tree = ast.parse(path.read_text(encoding="utf-8"))
        relative: set[str] = set()
        forbidden_calls: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.level:
                relative.update(alias.name for alias in node.names)
            elif isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                if node.func.id in {"eval", "exec", "compile", "__import__"}:
                    forbidden_calls.append(node.func.id)
        self.assertEqual(
            relative,
            {
                "validator_bio",
                "validator_insdc",
                "validator_reference",
                "validator_v2",
            },
        )
        self.assertEqual(forbidden_calls, [])


if __name__ == "__main__":
    unittest.main()

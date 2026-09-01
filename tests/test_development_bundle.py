from __future__ import annotations

from copy import deepcopy
import hashlib
import os
from pathlib import Path
import tempfile
import unittest

from brainc._canonical import artifact_digest, digest
from brainc.development_bundle import (
    BUNDLE_FILENAME,
    CHILD_FILENAMES,
    DevelopmentBundleError,
    compile_development,
)
from brainc.source import (
    EXTERNAL_PROFILE,
    FASTA_PROFILE,
    RAW_PROFILE,
    REFERENCE_FASTA_PROFILE,
    SourceBundle,
    compile_external_source,
    compile_source,
)
from brainc.v2 import (
    V2Error,
    inline_storage,
    make_development_request,
    make_request,
    pack,
    policy_artifact,
    save,
    target_artifact,
)
from brainc.v2.target import DEV_DOMAIN, OP_UNIT_CREATE
from tests.test_external_profile import FASTQ, NATIVE_INDEX, _closure


ROOT = Path(__file__).resolve().parents[1]
FASTA = ROOT / "tests/data/J02482.1.fasta"


def _seal(core: dict) -> dict:
    return {**core, "artifact_sha256": digest(core)}


class _Chain:
    def __init__(
        self,
        root: Path,
        profile: str = FASTA_PROFILE,
        source: SourceBundle | None = None,
    ) -> None:
        self.root = root
        self.root.mkdir(parents=True, exist_ok=True)
        self.source = source or compile_source(
            profile,
            {"sequence": FASTA},
            parameters={"wrapper": "identity"},
        )
        profile = self.source.profile
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
        output = {
            "id": "t.value",
            "type": {"dtype": "u64", "shape": []},
            "unit": None,
            "axes": [],
            "storage": inline_storage(pack("u64", [1])),
        }
        acceptance = f"brainc.source-descriptor/v1;profile={profile}"
        if profile == EXTERNAL_PROFILE:
            acceptance += ";manifest=" + self.source.to_dict()["source_ir"][
                "parameters"
            ]["profile_manifest_sha256"]
        self.manifest = _seal(
            {
                "format": "brainc.provider-manifest",
                "version": 2,
                "provider": {"name": "org.example.provider", "version": "1"},
                "model_identity": {
                    "kind": "content-sha256",
                    "value": digest({"algorithm": "caller-supplied", "version": 1}),
                },
                "accepts": [acceptance],
                "outputs": [
                    {key: deepcopy(output[key]) for key in ("id", "type", "unit", "axes")}
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

    def _write(self, name: str, value: dict) -> Path:
        path = self.root / f"{name}.json"
        save(value, path)
        return path

    def compile(self, **overrides):
        values = {
            "source_bundle": self.source_directory,
            "manifest": self.manifest_path,
            "request": self.request_path,
            "response": self.response_path,
            "policy": self.policy_path,
            "target": self.target_path,
            "blob_root": None,
        }
        values.update(overrides)
        return compile_development(**values)


class DevelopmentBundleTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.chain = _Chain(self.root)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_bundle_is_deterministic_flat_and_exactly_bound(self) -> None:
        first = self.chain.compile()
        second = self.chain.compile()
        self.assertEqual(first.to_dict(), second.to_dict())
        bundle = first.to_dict()
        artifacts = first.artifacts
        self.assertEqual(artifact_digest(bundle), bundle["artifact_sha256"])
        self.assertEqual(bundle["format"], "brainc.development-bundle")
        self.assertEqual(type(bundle["version"]), int)
        self.assertEqual(bundle["version"], 1)
        self.assertEqual(set(artifacts), set(CHILD_FILENAMES))
        self.assertEqual(set(bundle["artifacts"]), set(CHILD_FILENAMES))
        self.assertEqual(bundle["source"]["descriptor"]["artifact_sha256"], self.chain.source.to_dict()["artifact_sha256"])
        self.assertEqual(bundle["source"]["native"], self.chain.source.to_dict()["source_ir"]["artifacts"])
        for role, reference in bundle["artifacts"].items():
            self.assertEqual(
                set(reference),
                {"format", "version", "artifact_sha256"},
            )
            self.assertEqual(reference["artifact_sha256"], artifacts[role]["artifact_sha256"])
            self.assertNotIn("module", reference)
            self.assertNotIn("outputs", reference)

        record = artifacts["compilation_record"]
        module = artifacts["development_module"]
        self.assertEqual(record["result"]["module_sha256"], module["module_sha256"])
        self.assertEqual(record["result"]["budgets"], module["module"]["budgets"])
        self.assertEqual(module["module"]["budgets"], {"operations": 1, "tensor_bytes": 8, "units": 1, "edges": 0, "attachments": 0})
        self.assertEqual(module["sources"]["sequence"], bundle["source"]["descriptor"])
        self.assertEqual(bundle["blobs"], [])
        self.assertNotIn(str(self.root), repr(bundle))

    def test_reference_fasta_catalog_reaches_the_typed_development_target(self) -> None:
        chain = _Chain(self.root / "reference-chain", REFERENCE_FASTA_PROFILE)
        compiled = chain.compile()
        bundle = compiled.to_dict()
        module = compiled.artifacts["development_module"]
        self.assertEqual(chain.source.profile, REFERENCE_FASTA_PROFILE)
        self.assertEqual(set(chain.source.artifacts), {"reference"})
        self.assertEqual(
            bundle["source"]["native"],
            chain.source.to_dict()["source_ir"]["artifacts"],
        )
        self.assertEqual(
            module["sources"]["sequence"],
            bundle["source"]["descriptor"],
        )
        self.assertEqual(module["module"]["budgets"]["units"], 1)

    def test_external_profile_reaches_target_with_exact_manifest_acceptance(self) -> None:
        manifest, descriptor, report, _, _, _ = _closure()
        source = compile_external_source(
            manifest,
            descriptor,
            report,
            original_payloads={"reads": FASTQ},
            native_artifact_payloads={"read-index": NATIVE_INDEX},
        )
        chain = _Chain(self.root / "external-chain", source=source)
        compiled = chain.compile()
        expected_tag = (
            "brainc.source-descriptor/v1;profile=external-dna-source/v1;manifest="
            + manifest["artifact_sha256"]
        )
        self.assertEqual(chain.manifest["accepts"], [expected_tag])
        self.assertEqual(compiled.artifacts["development_module"]["sources"]["sequence"], compiled.to_dict()["source"]["descriptor"])

        for rejected in (
            "brainc.source-descriptor/v1;profile=external-dna-source/v1",
            expected_tag[:-64] + "0" * 64,
        ):
            attacked = deepcopy(chain.manifest)
            attacked["accepts"] = [rejected]
            attacked = _seal(
                {
                    key: value
                    for key, value in attacked.items()
                    if key != "artifact_sha256"
                }
            )
            path = chain._write("attacked-external-manifest", attacked)
            with self.assertRaises(V2Error):
                make_development_request(source, path, ["t.value"])

    def test_request_builder_requires_and_replays_the_native_source_closure(self) -> None:
        self.assertEqual(
            make_development_request(
                self.chain.source,
                self.chain.manifest_path,
                ["t.value"],
            ),
            self.chain.request,
        )
        with self.assertRaisesRegex(V2Error, "not a complete source"):
            make_request(
                self.chain.source_path,
                self.chain.manifest_path,
                ["t.value"],
            )
        (self.chain.source_directory / "sequence.json").unlink()
        with self.assertRaisesRegex(V2Error, "source bundle"):
            make_development_request(
                self.chain.source_directory,
                self.chain.manifest_path,
                ["t.value"],
            )

    def test_provider_must_accept_the_exact_descriptor_profile(self) -> None:
        exact_tag = f"brainc.source-descriptor/v1;profile={FASTA_PROFILE}"
        self.assertEqual(self.chain.manifest["accepts"], [exact_tag])

        rejected_tags = (
            "brainc.source-descriptor/v1",
            f"brainc.source-descriptor/v1;profile={RAW_PROFILE}",
        )
        for index, rejected_tag in enumerate(rejected_tags):
            with self.subTest(accepts=rejected_tag):
                manifest = deepcopy(self.chain.manifest)
                manifest["accepts"] = [rejected_tag]
                manifest = _seal(
                    {
                        key: value
                        for key, value in manifest.items()
                        if key != "artifact_sha256"
                    }
                )
                manifest_path = self.chain._write(f"wrong-profile-manifest-{index}", manifest)

                with self.assertRaisesRegex(V2Error, exact_tag):
                    make_development_request(
                        self.chain.source_directory,
                        manifest_path,
                        ["t.value"],
                    )

                request = deepcopy(self.chain.request)
                request["provider_manifest_sha256"] = manifest["artifact_sha256"]
                request = _seal(
                    {
                        key: value
                        for key, value in request.items()
                        if key != "artifact_sha256"
                    }
                )
                request_path = self.chain._write(f"wrong-profile-request-{index}", request)
                response = deepcopy(self.chain.response)
                response["request_artifact_sha256"] = request["artifact_sha256"]
                response = _seal(
                    {
                        key: value
                        for key, value in response.items()
                        if key != "artifact_sha256"
                    }
                )
                response_path = self.chain._write(
                    f"wrong-profile-response-{index}", response
                )
                with self.assertRaisesRegex(DevelopmentBundleError, exact_tag):
                    self.chain.compile(
                        manifest=manifest_path,
                        request=request_path,
                        response=response_path,
                    )

    def test_bundle_publication_contains_only_index_and_seven_children(self) -> None:
        compiled = self.chain.compile(source_bundle=self.chain.source)
        output = self.root / "development"
        paths = compiled.save(output)
        self.assertEqual(
            {item.name for item in output.iterdir()},
            {BUNDLE_FILENAME, *CHILD_FILENAMES.values()},
        )
        self.assertEqual(set(paths), {BUNDLE_FILENAME, *CHILD_FILENAMES.values()})
        self.assertNotIn("sequence.json", {item.name for item in output.iterdir()})
        self.assertNotIn("source.json", {item.name for item in output.iterdir()})

    def test_external_blob_set_is_derived_sorted_and_deduplicated(self) -> None:
        raw = pack("u64", [1])
        blob_sha = hashlib.sha256(raw).hexdigest()
        blob_root = self.root / "blobs"
        blob_root.mkdir()
        (blob_root / blob_sha).write_bytes(raw)
        response = deepcopy(self.chain.response)
        response["outputs"][0]["storage"] = {
            "kind": "sha256-blob",
            "byte_length": len(raw),
            "sha256": blob_sha,
        }
        response = _seal({key: value for key, value in response.items() if key != "artifact_sha256"})
        response_path = self.chain._write("external-response", response)
        bundle = self.chain.compile(response=response_path, blob_root=blob_root).to_dict()
        self.assertEqual(bundle["blobs"], [{"sha256": blob_sha, "byte_length": len(raw)}])

    def test_cross_artifact_splice_and_linked_inputs_fail_closed(self) -> None:
        request = deepcopy(self.chain.request)
        request["source"]["artifact_sha256"] = "0" * 64
        request_path = self.chain._write(
            "spliced-request",
            _seal({key: value for key, value in request.items() if key != "artifact_sha256"}),
        )
        with self.assertRaises(DevelopmentBundleError):
            self.chain.compile(request=request_path)

        if os.name == "posix":
            alias = self.root / "response-link.json"
            alias.symlink_to(self.chain.response_path)
            with self.assertRaises(DevelopmentBundleError):
                self.chain.compile(response=alias)


if __name__ == "__main__":
    unittest.main()

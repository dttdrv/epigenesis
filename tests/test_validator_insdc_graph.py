from __future__ import annotations

import ast
import base64
import copy
import json
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import unittest

import brainc.validator_insdc_graph as validator
from brainc.insdc import GenBankCompiler
from brainc.insdc_graph import compile_insdc_graph


ROOT = Path(__file__).resolve().parents[1]


def _reseal(artifact: dict) -> None:
    artifact["artifact_sha256"] = validator.digest(
        {key: value for key, value in artifact.items() if key != "artifact_sha256"}
    )


def _feature(
    feature_id: str,
    ordinal: int,
    key: str,
    segments: list[dict],
    qualifier_count: int,
) -> dict:
    return {
        "feature_id": feature_id,
        "ordinal": ordinal,
        "key": key,
        "location": {"segments": segments},
        "qualifiers": [{} for _ in range(qualifier_count)],
    }


class IndependentINSDCGraphValidatorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.raw = (ROOT / "tests/data/U49845.1.gb").read_bytes()
        cls.source = GenBankCompiler().compile_bytes(cls.raw).to_dict()
        compiled = compile_insdc_graph(cls.source)
        cls.index = compiled.to_dict()
        cls.artifacts = compiled.artifacts

    def test_actual_genbank_dna_reaches_independent_graph_and_v2_validators(self) -> None:
        report = validator.validate_insdc_graph_bundle(
            self.index,
            self.artifacts,
            self.source,
            genbank_source=self.raw,
        )
        self.assertTrue(report["valid"])
        self.assertEqual(
            report["result"],
            {
                "operations": 1,
                "tensor_bytes": 536,
                "units": 6,
                "edges": 0,
                "attachments": 0,
                "semantics_sha256": (
                    self.artifacts["backend_semantics"]["semantics_sha256"]
                ),
                "module_sha256": self.artifacts["development_module"][
                    "module_sha256"
                ],
            },
        )
        self.assertEqual(len(report["checks"]), 9)

    def test_embedded_exact_source_is_sufficient_for_independent_replay(self) -> None:
        report = validator.validate_insdc_graph_bundle(
            self.index,
            self.artifacts,
            self.source,
        )
        self.assertTrue(report["valid"])

    def test_bundle_index_is_reference_only_and_exactly_bound(self) -> None:
        self.assertEqual(
            set(self.index),
            {"format", "version", "backend", "source", "artifacts", "artifact_sha256"},
        )
        self.assertEqual(set(self.index["artifacts"]), set(validator.CHILD_ROLES))
        self.assertTrue(
            all(set(reference) == {"artifact_sha256"} for reference in self.index["artifacts"].values())
        )

    def test_coherently_resealed_tensor_substitution_is_rejected(self) -> None:
        artifacts = copy.deepcopy(self.artifacts)
        index = copy.deepcopy(self.index)
        response = artifacts["prediction_response"]
        outputs = {output["id"]: output for output in response["outputs"]}
        left, right = "t.feature.ordinal", "t.feature.orientation_mask"
        outputs[left]["storage"], outputs[right]["storage"] = (
            outputs[right]["storage"],
            outputs[left]["storage"],
        )
        _reseal(response)

        module = artifacts["development_module"]
        tensors = {tensor["id"]: tensor for tensor in module["module"]["tensors"]}
        tensors[left]["storage"], tensors[right]["storage"] = (
            tensors[right]["storage"],
            tensors[left]["storage"],
        )
        module["sources"]["prediction_response"]["artifact_sha256"] = response[
            "artifact_sha256"
        ]
        module["module_sha256"] = validator.digest(module["module"])
        _reseal(module)

        record = artifacts["compilation_record"]
        for role in ("prediction_response", "development_module"):
            record["artifacts"][role]["artifact_sha256"] = artifacts[role][
                "artifact_sha256"
            ]
        record["result"]["module_sha256"] = module["module_sha256"]
        _reseal(record)
        for role in ("prediction_response", "development_module", "compilation_record"):
            index["artifacts"][role]["artifact_sha256"] = artifacts[role][
                "artifact_sha256"
            ]
        _reseal(index)

        with self.assertRaisesRegex(
            validator.INSDCGraphValidationError,
            "differs from independent adapter replay",
        ):
            validator.validate_insdc_graph_bundle(
                index, artifacts, self.source, genbank_source=self.raw
            )

    def test_coherently_resealed_backend_spec_substitution_is_rejected(self) -> None:
        artifacts = copy.deepcopy(self.artifacts)
        index = copy.deepcopy(self.index)
        specification = artifacts["backend_spec"]
        specification["capabilities"]["learning"] = True
        _reseal(specification)
        changed_sha = specification["artifact_sha256"]

        semantics = artifacts["backend_semantics"]
        semantics["backend"]["spec_sha256"] = changed_sha
        _reseal(semantics)
        manifest = artifacts["provider_manifest"]
        manifest["model_identity"]["value"] = changed_sha
        _reseal(manifest)
        request = artifacts["prediction_request"]
        request["provider_manifest_sha256"] = manifest["artifact_sha256"]
        _reseal(request)
        response = artifacts["prediction_response"]
        response["request_artifact_sha256"] = request["artifact_sha256"]
        response["model_identity"]["value"] = changed_sha
        _reseal(response)
        module = artifacts["development_module"]
        for role in ("provider_manifest", "prediction_request", "prediction_response"):
            module["sources"][role]["artifact_sha256"] = artifacts[role][
                "artifact_sha256"
            ]
        _reseal(module)
        record = artifacts["compilation_record"]
        record["backend"]["spec_sha256"] = changed_sha
        for role, artifact in artifacts.items():
            if role != "compilation_record" and role in record["artifacts"]:
                record["artifacts"][role]["artifact_sha256"] = artifact[
                    "artifact_sha256"
                ]
        _reseal(record)
        index["backend"]["spec_sha256"] = changed_sha
        for role, artifact in artifacts.items():
            index["artifacts"][role]["artifact_sha256"] = artifact[
                "artifact_sha256"
            ]
        _reseal(index)

        with self.assertRaisesRegex(
            validator.INSDCGraphValidationError, "normative commitment"
        ):
            validator.validate_insdc_graph_bundle(
                index, artifacts, self.source, genbank_source=self.raw
            )

    def test_resealed_boolean_as_integer_is_rejected_without_python_coercion(self) -> None:
        artifacts = copy.deepcopy(self.artifacts)
        index = copy.deepcopy(self.index)
        record = artifacts["compilation_record"]
        record["result"]["inference"] = 0
        _reseal(record)
        index["artifacts"]["compilation_record"]["artifact_sha256"] = record[
            "artifact_sha256"
        ]
        _reseal(index)
        with self.assertRaisesRegex(
            validator.INSDCGraphValidationError,
            r"compilation_record differs.*\.inference type int != bool",
        ):
            validator.validate_insdc_graph_bundle(
                index, artifacts, self.source, genbank_source=self.raw
            )

    def test_graph_closure_rejects_jcs_equivalent_floats(self) -> None:
        cases = []

        artifacts = copy.deepcopy(self.artifacts)
        artifacts["backend_spec"]["target_contract"]["abi_major"] = 1.0
        cases.append(("backend specification", self.index, artifacts))

        artifacts = copy.deepcopy(self.artifacts)
        artifacts["compilation_record"]["result"]["ports"] = -0.0
        cases.append(("compilation record", self.index, artifacts))

        index = copy.deepcopy(self.index)
        index["version"] = 1.0
        cases.append(("bundle index", index, self.artifacts))

        for label, index, artifacts in cases:
            with self.subTest(label=label):
                with self.assertRaisesRegex(
                    validator.INSDCGraphValidationError,
                    "floating-point JSON number",
                ):
                    validator.validate_insdc_graph_bundle(
                        index,
                        artifacts,
                        self.source,
                        genbank_source=self.raw,
                    )

    def test_mismatched_raw_genbank_evidence_is_rejected(self) -> None:
        with self.assertRaisesRegex(
            validator.INSDCGraphValidationError,
            "independent GenBank source replay failed",
        ):
            validator.validate_insdc_graph_bundle(
                self.index,
                self.artifacts,
                self.source,
                genbank_source=self.raw.replace(b"Saccharomyces", b"Schizosacchar", 1),
            )

    def test_coherently_resealed_source_semantics_and_index_tampering_fail(self) -> None:
        source = copy.deepcopy(self.source)
        source["bio_ir"]["records"][0]["definition"] += " forged"
        source["bio_ir_sha256"] = validator.digest(source["bio_ir"])
        _reseal(source)
        with self.assertRaisesRegex(
            validator.INSDCGraphValidationError,
            "independent GenBank source replay failed",
        ):
            validator.validate_insdc_graph_bundle(
                self.index, self.artifacts, source, genbank_source=self.raw
            )

        artifacts = copy.deepcopy(self.artifacts)
        index = copy.deepcopy(self.index)
        semantics = artifacts["backend_semantics"]
        semantics["dictionaries"]["keys"].reverse()
        semantic_core = {
            key: semantics[key]
            for key in ("source", "dictionaries", "unit_order", "edge_order", "counts")
        }
        semantics["semantics_sha256"] = validator.digest(semantic_core)
        _reseal(semantics)
        record = artifacts["compilation_record"]
        record["backend"]["semantics_sha256"] = semantics["semantics_sha256"]
        record["artifacts"]["backend_semantics"]["artifact_sha256"] = semantics[
            "artifact_sha256"
        ]
        _reseal(record)
        for role in ("backend_semantics", "compilation_record"):
            index["artifacts"][role]["artifact_sha256"] = artifacts[role][
                "artifact_sha256"
            ]
        _reseal(index)
        with self.assertRaisesRegex(
            validator.INSDCGraphValidationError,
            "backend_semantics differs from independent adapter replay",
        ):
            validator.validate_insdc_graph_bundle(
                index, artifacts, self.source, genbank_source=self.raw
            )

        index = copy.deepcopy(self.index)
        index["source"]["bio_ir_sha256"] = "0" * 64
        _reseal(index)
        with self.assertRaisesRegex(
            validator.INSDCGraphValidationError,
            "bundle index differs from independent reference replay",
        ):
            validator.validate_insdc_graph_bundle(
                index, self.artifacts, self.source, genbank_source=self.raw
            )

    def test_runtime_validation_loads_no_producer_modules(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            payload_path = root / "payload.json"
            source_path = root / "source.gb"
            payload_path.write_text(
                json.dumps(
                    {
                        "index": self.index,
                        "artifacts": self.artifacts,
                        "source": self.source,
                    }
                ),
                encoding="utf-8",
            )
            source_path.write_bytes(self.raw)
            command = (
                "import json,sys; "
                f"sys.path.insert(0,{str(ROOT)!r}); "
                "from pathlib import Path; "
                "from brainc.validator_insdc_graph import validate_insdc_graph_bundle; "
                f"p=json.loads(Path({str(payload_path)!r}).read_text()); "
                f"raw=Path({str(source_path)!r}).read_bytes(); "
                "assert validate_insdc_graph_bundle(p['index'],p['artifacts'],"
                "p['source'],genbank_source=raw)['valid']; "
                "forbidden={'brainc.insdc','brainc.insdc_graph','brainc.insdc_location',"
                "'brainc._canonical','brainc._io','brainc.v2','brainc.v2.compiler',"
                "'brainc.v2.provider','brainc.v2.target'}; "
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

    def test_feature_state_replays_dictionaries_order_and_all_fields(self) -> None:
        bio_ir = {
            "records": [
                {
                    "record_id": "Z.1",
                    "features": [
                        _feature(
                            "insdc:z",
                            1,
                            "gene",
                            [
                                {
                                    "reference": None,
                                    "start": 10,
                                    "end": 20,
                                    "orientation": 1,
                                    "kind": "interval",
                                    "start_fuzz": "<",
                                    "end_fuzz": None,
                                    "bounds_status": "verified",
                                },
                                {
                                    "reference": "REMOTE.1",
                                    "start": 20,
                                    "end": 20,
                                    "orientation": -1,
                                    "kind": "between",
                                    "start_fuzz": None,
                                    "end_fuzz": ">",
                                    "bounds_status": "unresolved",
                                },
                            ],
                            2,
                        )
                    ],
                },
                {
                    "record_id": "A.1",
                    "features": [
                        _feature(
                            "insdc:a",
                            1,
                            "CDS",
                            [
                                {
                                    "reference": None,
                                    "start": 5,
                                    "end": 8,
                                    "orientation": 1,
                                    "kind": "uncertain-point",
                                    "start_fuzz": None,
                                    "end_fuzz": None,
                                    "bounds_status": "verified",
                                }
                            ],
                            0,
                        )
                    ],
                },
            ]
        }

        state = validator._feature_state(bio_ir)

        self.assertEqual(
            state["dictionaries"],
            {
                "keys": [{"code": 0, "value": "CDS"}, {"code": 1, "value": "gene"}],
                "records": [
                    {"code": 0, "value": "A.1"},
                    {"code": 1, "value": "Z.1"},
                ],
            },
        )
        self.assertEqual(state["feature_ids"], ["insdc:z", "insdc:a"])
        self.assertEqual(
            state["values"],
            {
                "fuzzy_boundary_count": [2, 0],
                "key_code": [1, 0],
                "kind_mask": [3, 4],
                "orientation_mask": [3, 1],
                "ordinal": [1, 1],
                "qualifier_count": [2, 0],
                "record_code": [1, 0],
                "remote_segment_count": [1, 0],
                "segment_extent_sum": [10, 3],
                "segment_count": [2, 1],
                "unresolved_segment_count": [1, 0],
            },
        )

    def test_u64_tensor_is_exact_little_endian_inline_storage(self) -> None:
        tensor = validator._tensor("t.feature.ordinal", [1, 2**32 + 3])
        raw = struct.pack("<QQ", 1, 2**32 + 3)
        self.assertEqual(tensor["type"], {"dtype": "u64", "shape": [2]})
        self.assertEqual(tensor["axes"], [None])
        self.assertEqual(
            base64.b64decode(tensor["storage"]["data"], validate=True), raw
        )
        self.assertEqual(tensor["storage"]["byte_length"], len(raw))
        self.assertEqual(
            validator._pack_u64([2**64 - 1]), struct.pack("<Q", 2**64 - 1)
        )
        with self.assertRaises(validator.INSDCGraphValidationError):
            validator._pack_u64([2**64])

    def test_nonordinal_feature_array_is_rejected(self) -> None:
        bio_ir = {
            "records": [
                {
                    "record_id": "X.1",
                    "features": [
                        _feature(
                            "insdc:x",
                            2,
                            "source",
                            [
                                {
                                    "reference": None,
                                    "start": 0,
                                    "end": 1,
                                    "orientation": 1,
                                    "kind": "interval",
                                    "start_fuzz": None,
                                    "end_fuzz": None,
                                    "bounds_status": "verified",
                                }
                            ],
                            0,
                        )
                    ],
                }
            ]
        }
        with self.assertRaisesRegex(
            validator.INSDCGraphValidationError, "record-array order then ordinal"
        ):
            validator._feature_state(bio_ir)

    def test_validator_source_has_no_producer_or_shared_canonical_imports(self) -> None:
        path = ROOT / "brainc/validator_insdc_graph.py"
        tree = ast.parse(path.read_text(encoding="utf-8"))
        imported: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module is not None:
                imported.append(node.module)
        forbidden = {
            "brainc.insdc",
            "brainc.insdc_graph",
            "brainc._canonical",
            "brainc._io",
            "brainc.v2",
        }
        self.assertFalse(forbidden.intersection(imported), imported)


if __name__ == "__main__":
    unittest.main()

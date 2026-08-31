from __future__ import annotations

import ast
import copy
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import brainc.validator_bio_graph as bridge_validator
from brainc._canonical import digest
from brainc.bio import GFF3Compiler
from brainc.bio_graph import compile_feature_graph
from brainc.sequence_collection import SequenceCollectionCompiler
from brainc.validator_bio_graph import (
    BioGraphValidationError,
    validate_feature_graph_bundle,
)
from brainc.validator_v2 import validate_chain


ROOT = Path(__file__).resolve().parents[1]


def _reseal(artifact: dict) -> None:
    artifact["artifact_sha256"] = digest(
        {key: value for key, value in artifact.items() if key != "artifact_sha256"}
    )


def _reseal_collection(collection: dict) -> None:
    for member in collection["members"]:
        sequence = member["sequence_artifact"]
        sequence["ir_sha256"] = digest(sequence["sequence_ir"])
        _reseal(sequence)
        member["record_id"] = sequence["sequence_ir"]["record_id"]
        member["ir_sha256"] = sequence["ir_sha256"]
        member["artifact_sha256"] = sequence["artifact_sha256"]
    collection["collection_ir"] = {
        "members": [
            {
                "record_id": member["record_id"],
                "ir_sha256": member["ir_sha256"],
            }
            for member in collection["members"]
        ]
    }
    collection["collection_ir_sha256"] = digest(collection["collection_ir"])
    _reseal(collection)


def _rewire_semantics(bundle: dict) -> None:
    semantics = bundle["artifacts"]["backend_semantics"]
    semantics["semantics_sha256"] = digest(semantics["semantics"])
    _reseal(semantics)
    record = bundle["artifacts"]["compilation_record"]
    record["backend"]["semantics_sha256"] = semantics["semantics_sha256"]
    record["artifacts"]["backend_semantics"]["artifact_sha256"] = semantics[
        "artifact_sha256"
    ]
    _reseal(record)
    _reseal(bundle)


def _rewire_response_and_module(bundle: dict) -> None:
    artifacts = bundle["artifacts"]
    response = artifacts["prediction_response"]
    _reseal(response)
    module = artifacts["development_module"]
    module["sources"]["prediction_response"]["artifact_sha256"] = response[
        "artifact_sha256"
    ]
    module["module_sha256"] = digest(module["module"])
    _reseal(module)
    record = artifacts["compilation_record"]
    record["artifacts"]["prediction_response"]["artifact_sha256"] = response[
        "artifact_sha256"
    ]
    record["artifacts"]["development_module"]["artifact_sha256"] = module[
        "artifact_sha256"
    ]
    record["result"]["module_sha256"] = module["module_sha256"]
    _reseal(record)
    _reseal(bundle)


class IndependentBioGraphValidatorTests(unittest.TestCase):
    def _chain(self, root: Path, fasta: bytes, gff3: bytes):
        collection = SequenceCollectionCompiler().compile_fasta_bytes(fasta)
        source_path = root / "sequence.json"
        source_path.parent.mkdir(parents=True, exist_ok=True)
        collection.save(source_path)
        bio = GFF3Compiler().compile_bytes(gff3, collection)
        bundle = compile_feature_graph(source_path, bio, gff3_source=gff3)
        return collection, source_path, bio, bundle

    def test_actual_ncbi_dna_reaches_independent_bridge_and_v2_validators(self) -> None:
        fasta_path = ROOT / "tests/data/J02482.1.fasta"
        gff3_path = ROOT / "tests/data/J02482.1.gff3"
        gff3 = gff3_path.read_bytes()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            collection = SequenceCollectionCompiler().compile_file(fasta_path)
            source_path = root / "sequence.json"
            collection.save(source_path)
            bio = GFF3Compiler().compile_file(gff3_path, collection)
            bundle = compile_feature_graph(source_path, bio, gff3_source=gff3)

            report = validate_feature_graph_bundle(
                bundle.to_dict(), collection.to_dict(), bio.to_dict()
            )
            self.assertTrue(report["valid"])
            self.assertEqual(
                report["result"],
                {
                    "operations": 3,
                    "tensor_bytes": 2158,
                    "units": 23,
                    "edges": 4,
                    "attachments": 0,
                    "relationship_edges": 4,
                    "semantics_sha256": bundle.semantics["semantics_sha256"],
                    "module_sha256": bundle.development_module["module_sha256"],
                },
            )

            paths = bundle.save(root / "bundle")
            v2_report = validate_chain(
                fasta=fasta_path,
                sequence=source_path,
                manifest=paths["provider_manifest"],
                request=paths["prediction_request"],
                response=paths["prediction_response"],
                policy=paths["lowering_policy"],
                target=paths["target_contract"],
                module=paths["development_module"],
            )
            self.assertTrue(v2_report["valid"], v2_report)
            self.assertEqual(v2_report["passed_checks"], 4)

    def test_discontinuous_features_do_not_create_quadratic_edges(self) -> None:
        fasta = b">arbitrary\n" + b"A" * 120 + b"\n"
        gff3 = (
            b"##gff-version 3\n"
            b"arbitrary\t.\tgene\t1\t40\t.\t+\t.\tID=parent\n"
            b"arbitrary\t.\tgene\t60\t100\t.\t+\t.\tID=parent\n"
            b"arbitrary\t.\texon\t2\t10\t.\t+\t.\tID=child;Parent=parent\n"
            b"arbitrary\t.\texon\t20\t30\t.\t+\t.\tID=child;Parent=parent\n"
        )
        with tempfile.TemporaryDirectory() as temporary:
            collection, _, bio, bundle = self._chain(
                Path(temporary), fasta, gff3
            )
            report = validate_feature_graph_bundle(
                bundle.to_dict(), collection.to_dict(), bio.to_dict()
            )
            self.assertEqual(len(bio.to_dict()["bio_ir"]["features"]), 2)
            self.assertEqual(len(bio.to_dict()["bio_ir"]["relationships"]), 1)
            self.assertEqual(report["result"]["units"], 2)
            self.assertEqual(report["result"]["edges"], 1)
            pairs = next(
                tensor
                for tensor in bundle.development_module["module"]["tensors"]
                if tensor["id"] == "t.relationship.parent.pairs"
            )
            self.assertEqual(pairs["type"]["shape"], [1, 2])

    def test_empty_annotation_is_a_valid_zero_state_graph(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            collection, _, bio, bundle = self._chain(
                Path(temporary),
                b">unannotated\nACGTACGT\n",
                b"##gff-version 3\n",
            )
            report = validate_feature_graph_bundle(
                bundle.to_dict(), collection.to_dict(), bio.to_dict()
            )
            self.assertEqual(report["result"]["units"], 0)
            self.assertEqual(report["result"]["edges"], 0)
            self.assertEqual(report["result"]["tensor_bytes"], 24)

    def test_descriptor_safe_paths_reject_symlink_and_oversized_source(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            collection, source_path, bio, bundle = self._chain(
                root,
                b">x\nACGTACGT\n",
                b"##gff-version 3\nx\t.\tgene\t1\t8\t.\t+\t.\tID=g\n",
            )
            bio_path = root / "bio.json"
            bio.save(bio_path)
            paths = bundle.save(root / "bundle")
            report = validate_feature_graph_bundle(
                paths["bundle"], source_path, bio_path
            )
            self.assertTrue(report["valid"])

            replacement = root / "replacement-bundle.json"
            replacement.write_bytes(paths["bundle"].read_bytes())
            real_open = bridge_validator.os.open
            swapped = False

            def replace_between_inspection_and_open(path, flags):
                nonlocal swapped
                if Path(path) == paths["bundle"] and not swapped:
                    swapped = True
                    replacement.replace(paths["bundle"])
                return real_open(path, flags)

            with mock.patch.object(
                bridge_validator.os,
                "open",
                side_effect=replace_between_inspection_and_open,
            ):
                with self.assertRaisesRegex(
                    BioGraphValidationError, "changed while being opened"
                ):
                    validate_feature_graph_bundle(
                        paths["bundle"], source_path, bio_path
                    )

            linked = root / "linked-bundle.json"
            linked.symlink_to(paths["bundle"])
            with self.assertRaisesRegex(
                BioGraphValidationError, "regular non-linked file"
            ):
                validate_feature_graph_bundle(linked, source_path, bio_path)

            oversized = root / "oversized-source.json"
            with oversized.open("wb") as stream:
                stream.truncate(bridge_validator.MAX_SOURCE_JSON_BYTES + 1)
            with self.assertRaisesRegex(
                BioGraphValidationError,
                f"exceeds {bridge_validator.MAX_SOURCE_JSON_BYTES} bytes",
            ):
                validate_feature_graph_bundle(
                    bundle.to_dict(), oversized, bio.to_dict()
                )

    def test_collection_rejects_coherently_resealed_nested_context(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            collection, _, bio, bundle = self._chain(
                Path(temporary),
                b">x\nACGTACGT\n",
                b"##gff-version 3\nx\t.\tgene\t1\t8\t.\t+\t.\tID=g\n",
            )
            base = collection.to_dict()
            with_context = copy.deepcopy(base)
            nested = with_context["members"][0]["sequence_artifact"]
            nested["inputs"]["context_sha256"] = "1" * 64
            nested["sequence_ir"]["reference"] = {
                "assembly": "explicit",
                "contig": "x",
                "start": 0,
                "end": 8,
                "coordinate_system": "0-based-half-open",
                "orientation": "forward",
                "aliases": [],
            }
            nested["sequence_ir"]["provenance"] = [
                {
                    "id": "source",
                    "kind": "test",
                    "uri": "urn:test",
                    "version": "1",
                }
            ]
            _reseal_collection(with_context)

            assertions_without_context = copy.deepcopy(base)
            nested = assertions_without_context["members"][0]["sequence_artifact"]
            nested["sequence_ir"]["provenance"] = [
                {
                    "id": "source",
                    "kind": "test",
                    "uri": "urn:test",
                    "version": "1",
                }
            ]
            _reseal_collection(assertions_without_context)

            cases = (
                (with_context, "context not emitted by collection ingress"),
                (assertions_without_context, "context assertions without a context input"),
            )
            for forged, message in cases:
                with self.subTest(message=message):
                    with self.assertRaisesRegex(BioGraphValidationError, message):
                        validate_feature_graph_bundle(
                            bundle.to_dict(), forged, bio.to_dict()
                        )

    def test_raw_collection_rejects_all_coherently_resealed_ingress_mismatches(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            collection = SequenceCollectionCompiler().compile_raw("ACGTACGT", "raw")
            source_path = root / "source.json"
            collection.save(source_path)
            gff3 = b"##gff-version 3\nraw\t.\tgene\t1\t8\t.\t+\t.\tID=g\n"
            bio = GFF3Compiler().compile_bytes(gff3, collection)
            bundle = compile_feature_graph(
                source_path, bio, gff3_source=gff3
            ).to_dict()

            description = collection.to_dict()
            description["members"][0]["sequence_artifact"]["sequence_ir"][
                "description"
            ] = "forged"
            _reseal_collection(description)

            synthetic_fasta = collection.to_dict()
            synthetic_fasta["members"][0]["sequence_artifact"]["inputs"][
                "fasta_sha256"
            ] = "2" * 64
            _reseal_collection(synthetic_fasta)

            logical = collection.to_dict()
            logical["inputs"]["raw_sha256"] = "3" * 64
            logical["inputs"]["logical_sha256"] = "3" * 64
            _reseal_collection(logical)

            cases = (
                (description, "cannot have a description"),
                (synthetic_fasta, "synthetic FASTA digest is inconsistent"),
                (logical, "logical digest does not match sequence"),
            )
            for forged, message in cases:
                with self.subTest(message=message):
                    with self.assertRaisesRegex(BioGraphValidationError, message):
                        validate_feature_graph_bundle(bundle, forged, bio.to_dict())

    def test_collection_rejects_overlapping_global_maps_and_total_bases(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            collection, _, bio, bundle = self._chain(
                Path(temporary),
                b">one\nAAAA\n>two\nCCCC\n",
                b"##gff-version 3\none\t.\tgene\t1\t4\t.\t+\t.\tID=g\n",
            )
            overlapping = collection.to_dict()
            overlapping["members"][1]["input_source_map"]["sequence_segments"][0][
                "line"
            ] = overlapping["members"][0]["input_source_map"][
                "sequence_segments"
            ][-1]["line"]
            _reseal_collection(overlapping)
            with self.assertRaisesRegex(
                BioGraphValidationError, "overlaps the preceding record"
            ):
                validate_feature_graph_bundle(
                    bundle.to_dict(), overlapping, bio.to_dict()
                )

            with mock.patch.object(
                bridge_validator, "MAX_LOGICAL_DNA_BYTES", 7
            ):
                with self.assertRaisesRegex(
                    BioGraphValidationError, "logical DNA bytes"
                ):
                    validate_feature_graph_bundle(
                        bundle.to_dict(), collection.to_dict(), bio.to_dict()
                    )

    def test_collection_rejects_impossible_source_map_line(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            collection, _, bio, bundle = self._chain(
                Path(temporary),
                b">x\nACGTACGT\n",
                b"##gff-version 3\nx\t.\tgene\t1\t8\t.\t+\t.\tID=g\n",
            )
            forged = collection.to_dict()
            forged["members"][0]["input_source_map"]["sequence_segments"][0][
                "line"
            ] = 2**53 - 1
            forged["members"][0]["sequence_artifact"]["source_map"][
                "sequence_segments"
            ][0]["line"] = 2**53 - 1
            _reseal_collection(forged)
            with self.assertRaisesRegex(BioGraphValidationError, "1000001"):
                validate_feature_graph_bundle(
                    bundle.to_dict(), forged, bio.to_dict()
                )

    def test_total_feature_segment_limit_is_not_entity_count(self) -> None:
        gff3 = (
            b"##gff-version 3\n"
            b"x\t.\tgene\t1\t2\t.\t+\t.\tID=one\n"
            b"x\t.\tgene\t4\t5\t.\t+\t.\tID=one\n"
            b"x\t.\tgene\t7\t8\t.\t+\t.\tID=one\n"
        )
        with tempfile.TemporaryDirectory() as temporary:
            collection, _, bio, bundle = self._chain(
                Path(temporary), b">x\nACGTACGT\n", gff3
            )
            self.assertEqual(len(bio.to_dict()["bio_ir"]["features"]), 1)
            self.assertEqual(
                len(bio.to_dict()["bio_ir"]["features"][0]["segments"]), 3
            )
            with mock.patch.object(bridge_validator, "MAX_FEATURE_ROWS", 2):
                with self.assertRaisesRegex(
                    BioGraphValidationError, "total feature segments"
                ):
                    validate_feature_graph_bundle(
                        bundle.to_dict(), collection.to_dict(), bio.to_dict()
                    )

    def test_opaque_attribute_value_order_and_duplicates_are_preserved(self) -> None:
        gff3 = (
            b"##gff-version 3\n"
            b"x\t.\tgene\t1\t8\t.\t+\t.\tID=g;custom=z,a,z\n"
        )
        with tempfile.TemporaryDirectory() as temporary:
            collection, _, bio, bundle = self._chain(
                Path(temporary), b">x\nACGTACGT\n", gff3
            )
            attributes = bio.to_dict()["bio_ir"]["features"][0]["segments"][0][
                "attributes"
            ]
            custom = next(item for item in attributes if item["tag"] == "custom")
            self.assertEqual(custom["values"], ["z", "a", "z"])
            report = validate_feature_graph_bundle(
                bundle.to_dict(), collection.to_dict(), bio.to_dict()
            )
            self.assertTrue(report["valid"])

    def test_coherently_resealed_dictionary_and_semantics_substitutions_fail(self) -> None:
        fasta = b">x\nACGTACGTACGT\n"
        gff3 = (
            b"##gff-version 3\n"
            b"x\ttool-a\tgene\t1\t12\t.\t+\t.\tID=g\n"
            b"x\ttool-b\texon\t2\t5\t.\t+\t.\tID=e;Parent=g\n"
        )
        with tempfile.TemporaryDirectory() as temporary:
            collection, _, bio, compiled = self._chain(Path(temporary), fasta, gff3)

            swapped = copy.deepcopy(compiled.to_dict())
            dictionaries = swapped["artifacts"]["backend_semantics"]["semantics"][
                "dictionaries"
            ]
            dictionaries["types"][0]["value"], dictionaries["types"][1]["value"] = (
                dictionaries["types"][1]["value"],
                dictionaries["types"][0]["value"],
            )
            _rewire_semantics(swapped)

            substituted = copy.deepcopy(compiled.to_dict())
            substituted["artifacts"]["backend_semantics"]["semantics"]["bio_ir"][
                "features"
            ][0]["type"] = "substituted-type"
            _rewire_semantics(substituted)

            for forged in (swapped, substituted):
                with self.subTest(forged=forged["artifact_sha256"]):
                    with self.assertRaisesRegex(
                        BioGraphValidationError, "independent BioIR-to-module replay"
                    ):
                        validate_feature_graph_bundle(
                            forged, collection.to_dict(), bio.to_dict()
                        )

    def test_resealed_backend_specification_substitution_fails_commitment(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            collection, _, bio, compiled = self._chain(
                Path(temporary),
                b">x\nACGTACGT\n",
                b"##gff-version 3\nx\t.\tgene\t1\t8\t.\t+\t.\tID=g\n",
            )
            forged = copy.deepcopy(compiled.to_dict())
            specification = forged["artifacts"]["backend_spec"]
            specification["capabilities"]["learning"] = True
            _reseal(specification)
            _reseal(forged)
            with self.assertRaisesRegex(
                BioGraphValidationError, "normative commitment"
            ):
                validate_feature_graph_bundle(
                    forged, collection.to_dict(), bio.to_dict()
                )

    def test_coherent_tensor_and_compilation_record_substitutions_fail(self) -> None:
        fasta = b">x\nACGTACGTACGT\n"
        gff3 = (
            b"##gff-version 3\n"
            b"x\t.\tgene\t1\t12\t.\t+\t.\tID=g\n"
            b"x\t.\texon\t2\t5\t.\t+\t.\tID=e;Parent=g\n"
        )
        with tempfile.TemporaryDirectory() as temporary:
            collection, _, bio, compiled = self._chain(Path(temporary), fasta, gff3)

            tensor_forgery = copy.deepcopy(compiled.to_dict())
            response_outputs = {
                output["id"]: output
                for output in tensor_forgery["artifacts"]["prediction_response"][
                    "outputs"
                ]
            }
            module_tensors = {
                tensor["id"]: tensor
                for tensor in tensor_forgery["artifacts"]["development_module"][
                    "module"
                ]["tensors"]
            }
            left, right = "t.feature.attribute_count", "t.feature.interval_count"
            response_outputs[left]["storage"], response_outputs[right]["storage"] = (
                response_outputs[right]["storage"],
                response_outputs[left]["storage"],
            )
            module_tensors[left]["storage"], module_tensors[right]["storage"] = (
                module_tensors[right]["storage"],
                module_tensors[left]["storage"],
            )
            _rewire_response_and_module(tensor_forgery)

            record_forgery = copy.deepcopy(compiled.to_dict())
            record = record_forgery["artifacts"]["compilation_record"]
            record["result"]["units"] += 1
            _reseal(record)
            _reseal(record_forgery)

            for forged in (tensor_forgery, record_forgery):
                with self.subTest(forged=forged["artifact_sha256"]):
                    with self.assertRaisesRegex(
                        BioGraphValidationError, "independent BioIR-to-module replay"
                    ):
                        validate_feature_graph_bundle(
                            forged, collection.to_dict(), bio.to_dict()
                        )

    def test_production_validator_ast_has_no_compiler_imports(self) -> None:
        path = ROOT / "brainc/validator_bio_graph.py"
        tree = ast.parse(path.read_text(encoding="utf-8"))
        imported: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module is not None:
                imported.append(node.module)
        self.assertFalse(
            [name for name in imported if name == "brainc" or name.startswith("brainc.")],
            imported,
        )


if __name__ == "__main__":
    unittest.main()

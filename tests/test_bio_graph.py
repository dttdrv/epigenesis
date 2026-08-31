from __future__ import annotations

import copy
import errno
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import brainc.bio_graph as bio_graph_module
from brainc._canonical import digest
from brainc.bio import GFF3Compiler
from brainc.bio_graph import (
    BACKEND_ID,
    BACKEND_SPEC,
    BACKEND_SPEC_SHA256,
    FEATURE_FIELD_SPECS,
    FeatureGraphError,
    compile_feature_graph,
    validate_feature_graph_bundle,
)
from brainc.sequence_collection import SequenceCollectionCompiler


ROOT = Path(__file__).resolve().parents[1]
RUNTIME_ROOT = ROOT.parent / "epigenesis-runtime"


def _reseal(payload: dict) -> dict:
    payload["artifact_sha256"] = digest(
        {key: value for key, value in payload.items() if key != "artifact_sha256"}
    )
    return payload


class FeatureGraphBackendTests(unittest.TestCase):
    def _chain(
        self,
        directory: Path,
        fasta: bytes,
        gff3: bytes,
    ):
        collection = SequenceCollectionCompiler().compile_fasta_bytes(fasta)
        source_path = directory / "source.json"
        collection.save(source_path)
        bio = GFF3Compiler().compile_bytes(gff3, collection)
        bundle = compile_feature_graph(
            source_path, bio, gff3_source=gff3
        )
        return collection, source_path, bio, bundle

    def test_complete_bio_ir_sidecar_mapping_has_no_biological_inference(self) -> None:
        fasta = b">arbitrary.record\nACGTACGTACGTACGT\n"
        gff3 = (
            b"##gff-version 3\n"
            b"##pipeline opaque-value\n"
            b"arbitrary.record\ttool alpha\tregion\t1\t16\t.\t.\t.\t"
            b"ID=reference;custom_tag=opaque value\n"
            b"arbitrary.record\ttool alpha\topaque_type\t2\t8\t1.25\t?\t.\t"
            b"ID=child;Parent=reference;Derives_from=reference;Note=kept\n"
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, source_path, bio, first = self._chain(root, fasta, gff3)
            second = compile_feature_graph(
                source_path, bio, gff3_source=gff3
            )

            self.assertEqual(first.to_dict(), second.to_dict())
            ir = bio.to_dict()["bio_ir"]
            semantics = first.semantics
            self.assertEqual(semantics["semantics"]["bio_ir"], ir)
            self.assertEqual(semantics["semantics_sha256"], digest(semantics["semantics"]))
            self.assertEqual(
                semantics["backend"]["id"], BACKEND_ID
            )
            self.assertEqual(
                [entry["value"] for entry in semantics["semantics"]["dictionaries"]["types"]],
                ["opaque_type", "region"],
            )

            module = first.development_module
            self.assertEqual(module["module"]["budgets"]["units"], len(ir["features"]))
            self.assertEqual(module["module"]["budgets"]["edges"], len(ir["relationships"]))
            self.assertEqual(module["module"]["budgets"]["attachments"], 0)
            self.assertEqual(first.target_contract["contract"]["ports"], [])
            self.assertEqual(first.target_contract["contract"]["rules"], [])

            # Both explicit relationship kinds may share endpoints.  They remain
            # distinct typed edges rather than being collapsed by a graph set.
            edge_order = semantics["semantics"]["edge_order"]
            self.assertEqual(len(edge_order), 2)
            self.assertEqual({edge["kind"] for edge in edge_order}, {"parent", "derives-from"})
            self.assertEqual(
                {(edge["source_unit_index"], edge["target_unit_index"]) for edge in edge_order},
                {(0, 1)},
            )

            record = first.compilation_record
            artifact_map = {
                "backend_spec": first.backend_spec,
                "backend_semantics": first.semantics,
                "development_module": first.development_module,
                "lowering_policy": first.lowering_policy,
                "prediction_request": first.prediction_request,
                "prediction_response": first.prediction_response,
                "provider_manifest": first.provider_manifest,
                "target_contract": first.target_contract,
            }
            for role, artifact in artifact_map.items():
                self.assertEqual(
                    record["artifacts"][role]["artifact_sha256"],
                    artifact["artifact_sha256"],
                )
            self.assertFalse(record["result"]["inference"])
            self.assertFalse(record["result"]["learning"])
            validate_feature_graph_bundle(
                first, source_path, bio, gff3_source=gff3
            )

    def test_arbitrary_content_changes_outputs_without_fixture_branches(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, _, first_bio, first = self._chain(
                root / "first",
                b">one\nAAAAAAAAAA\n",
                b"##gff-version 3\n"
                b"one\t.\talpha\t2\t4\t.\t+\t.\tID=a\n",
            )
            _, _, second_bio, second = self._chain(
                root / "second",
                b">totally-different\nCCCCCCCCCCCC\n",
                b"##gff-version 3\n"
                b"totally-different\t.\tbeta\t1\t3\t.\t-\t.\tID=b\n"
                b"totally-different\t.\tgamma\t5\t9\t.\t?\t.\tID=c;Parent=b\n",
            )
            self.assertNotEqual(
                first_bio.to_dict()["artifact_sha256"],
                second_bio.to_dict()["artifact_sha256"],
            )
            self.assertNotEqual(
                first.development_module["artifact_sha256"],
                second.development_module["artifact_sha256"],
            )
            self.assertEqual(first.development_module["module"]["budgets"]["units"], 1)
            self.assertEqual(second.development_module["module"]["budgets"]["units"], 2)
            self.assertEqual(second.development_module["module"]["budgets"]["edges"], 1)

    def test_backend_identity_commits_to_the_complete_normative_mapping(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, _, _, bundle = self._chain(
                root,
                b">x\nACGTACGT\n",
                b"##gff-version 3\nx\t.\tgene\t1\t8\t.\t+\t.\tID=g\n",
            )
            self.assertEqual(BACKEND_SPEC_SHA256, digest(BACKEND_SPEC))
            self.assertEqual(
                bundle.backend_spec,
                {**BACKEND_SPEC, "artifact_sha256": BACKEND_SPEC_SHA256},
            )
            self.assertEqual(
                [field["id"] for field in FEATURE_FIELD_SPECS],
                [
                    "attribute_count",
                    "declared_id",
                    "interval_count",
                    "location_max_end",
                    "location_min_start",
                    "phase_mask",
                    "segment_bases",
                    "segment_count",
                    "seqid_code",
                    "source_code",
                    "strand_code",
                    "type_code",
                    "wraps_origin",
                ],
            )
            self.assertTrue(
                all(type(field["formula"]) is dict for field in FEATURE_FIELD_SPECS)
            )
            self.assertEqual(
                BACKEND_SPEC["target_contract"],
                bundle.target_contract["contract"],
            )
            self.assertEqual(
                BACKEND_SPEC["operations"], bundle.lowering_policy["operations"]
            )
            self.assertEqual(
                BACKEND_SPEC["semantics_sidecar"]["dictionary_wire"]["container"],
                "object-keyed-by-dictionary-id",
            )
            self.assertEqual(
                BACKEND_SPEC["semantics_sidecar"]["edge_order"]["indices"],
                {
                    "edge_index": "zero-based-position-in-emitted-edge-order",
                    "bio_ir_relationship_index": (
                        "zero-based-position-in-bio-ir-relationships"
                    ),
                },
            )
            self.assertTrue(
                all(
                    relationship["relationship_count"]
                    == {"op": "length", "of": "selection"}
                    for relationship in BACKEND_SPEC["tensors"]["relationships"]
                )
            )
            tensor_ids = {"t.feature.count"}
            tensor_ids.update(field["tensor_id"] for field in FEATURE_FIELD_SPECS)
            for relationship in BACKEND_SPEC["tensors"]["relationships"]:
                tensor_ids.add(relationship["kind_tensor"]["id"])
                tensor_ids.add(relationship["pairs_tensor"]["id"])
            self.assertEqual(
                tensor_ids,
                {output["id"] for output in bundle.prediction_response["outputs"]},
            )
            for artifact in (
                bundle.semantics,
                bundle.compilation_record,
            ):
                self.assertEqual(
                    artifact["backend"]["spec_sha256"], BACKEND_SPEC_SHA256
                )
            self.assertEqual(
                bundle.provider_manifest["model_identity"]["value"],
                BACKEND_SPEC_SHA256,
            )

            changed = copy.deepcopy(BACKEND_SPEC)
            changed["tensors"]["feature_fields"][0]["formula"]["op"] = "product"
            self.assertNotEqual(digest(changed), BACKEND_SPEC_SHA256)
            import brainc.validator_bio_graph as independent

            self.assertEqual(independent.BACKEND_SPEC_SHA256, BACKEND_SPEC_SHA256)

    def test_public_spec_mutation_cannot_taint_compilation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            collection = SequenceCollectionCompiler().compile_raw("ACGT", "x")
            source_path = root / "source.json"
            collection.save(source_path)
            source = b"##gff-version 3\nx\t.\tgene\t1\t4\t.\t+\t.\tID=g\n"
            bio = GFF3Compiler().compile_bytes(source, collection)

            original = BACKEND_SPEC["capabilities"]["learning"]
            BACKEND_SPEC["capabilities"]["learning"] = True
            try:
                bundle = compile_feature_graph(
                    source_path, bio, gff3_source=source
                )
            finally:
                BACKEND_SPEC["capabilities"]["learning"] = original
            self.assertFalse(bundle.backend_spec["capabilities"]["learning"])
            self.assertEqual(
                bundle.backend_spec["artifact_sha256"], BACKEND_SPEC_SHA256
            )
            forged = copy.deepcopy(bundle.to_dict())
            forged["artifacts"]["backend_spec"]["capabilities"][
                "learning"
            ] = True
            _reseal(forged)
            with self.assertRaisesRegex(
                FeatureGraphError, "invalid backend_spec artifact digest"
            ):
                validate_feature_graph_bundle(
                    forged, source_path, bio, gff3_source=source
                )

            previous_strands = bio_graph_module.STRANDS
            bio_graph_module.STRANDS = ("-", "+", ".", "?")
            try:
                with self.assertRaisesRegex(
                    FeatureGraphError, "implementation constants drifted"
                ):
                    compile_feature_graph(
                        source_path, bio, gff3_source=source
                    )
            finally:
                bio_graph_module.STRANDS = previous_strands

    def test_compile_snapshots_source_for_the_entire_v2_transaction(self) -> None:
        gff3 = b"##gff-version 3\nx\t.\tgene\t1\t4\t.\t+\t.\tID=g\n"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            original = SequenceCollectionCompiler().compile_raw("ACGT", "x")
            replacement = SequenceCollectionCompiler().compile_raw("TGCA", "x")
            source_path = root / "source.json"
            original.save(source_path)
            bio = GFF3Compiler().compile_bytes(gff3, original)
            from brainc.bio_graph import make_request as graph_make_request

            observed_snapshot: Path | None = None

            def swap_before_request(*args, **kwargs):
                nonlocal observed_snapshot
                observed_snapshot = Path(args[0])
                replacement.save(source_path)
                return graph_make_request(*args, **kwargs)

            with mock.patch(
                "brainc.bio_graph.make_request", side_effect=swap_before_request
            ):
                bundle = compile_feature_graph(
                    source_path, bio, gff3_source=gff3
                )
            self.assertIsNotNone(observed_snapshot)
            self.assertNotEqual(observed_snapshot, source_path)
            wanted = original.to_dict()["artifact_sha256"]
            self.assertEqual(
                bundle.prediction_request["source"]["artifact_sha256"], wanted
            )
            self.assertEqual(
                bundle.compilation_record["inputs"]["sequence"][
                    "artifact_sha256"
                ],
                wanted,
            )

    def test_coherently_resealed_semantics_tampering_fails_replay(self) -> None:
        fasta = b">x\nACGTACGT\n"
        gff3 = b"##gff-version 3\nx\t.\tgene\t1\t8\t.\t+\t.\tID=g\n"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, source_path, bio, bundle = self._chain(root, fasta, gff3)
            forged = copy.deepcopy(bundle.to_dict())
            semantics = forged["artifacts"]["backend_semantics"]
            semantics["semantics"]["bio_ir"]["features"][0]["type"] = "forged"
            semantics["semantics_sha256"] = digest(semantics["semantics"])
            _reseal(semantics)
            record = forged["artifacts"]["compilation_record"]
            record["backend"]["semantics_sha256"] = semantics["semantics_sha256"]
            record["artifacts"]["backend_semantics"]["artifact_sha256"] = semantics[
                "artifact_sha256"
            ]
            _reseal(record)
            _reseal(forged)
            with self.assertRaisesRegex(FeatureGraphError, "deterministic lowering replay"):
                validate_feature_graph_bundle(
                    forged, source_path, bio, gff3_source=gff3
                )

    def test_bio_ir_must_match_supplied_sequence_collection(self) -> None:
        gff3 = b"##gff-version 3\nx\t.\tgene\t1\t4\t.\t+\t.\tID=g\n"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            original = SequenceCollectionCompiler().compile_raw("ACGT", "x")
            bio = GFF3Compiler().compile_bytes(gff3, original)
            changed = SequenceCollectionCompiler().compile_raw("ACGA", "x")
            changed_path = root / "changed.json"
            changed.save(changed_path)
            with self.assertRaisesRegex(FeatureGraphError, "binding"):
                compile_feature_graph(changed_path, bio, gff3_source=gff3)

    def test_compile_rejects_before_return_when_graph_wire_is_unsavable(self) -> None:
        fasta = b">x\nACGTACGT\n"
        gff3 = b"##gff-version 3\nx\t.\tgene\t1\t8\t.\t+\t.\tID=g\n"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            collection = SequenceCollectionCompiler().compile_fasta_bytes(fasta)
            source_path = root / "source.json"
            collection.save(source_path)
            bio = GFF3Compiler().compile_bytes(gff3, collection)
            from brainc.v2._common import V2Error, pretty_bytes

            def reject_semantics(payload: dict) -> bytes:
                if payload.get("format") == "brainc.bio.feature-graph-semantics":
                    raise V2Error("serialized output exceeds JSON byte limit")
                return pretty_bytes(payload)

            with mock.patch(
                "brainc.bio_graph.pretty_bytes", side_effect=reject_semantics
            ):
                with self.assertRaisesRegex(FeatureGraphError, "graph wire contract"):
                    compile_feature_graph(source_path, bio, gff3_source=gff3)

    def test_bundle_directory_publish_rejects_symlinks_and_is_transactional(self) -> None:
        fasta = b">x\nACGTACGT\n"
        gff3 = b"##gff-version 3\nx\t.\tgene\t1\t8\t.\t+\t.\tID=g\n"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, _, _, bundle = self._chain(root, fasta, gff3)
            victim = root / "victim"
            victim.mkdir()
            linked = root / "linked"
            linked.symlink_to(victim, target_is_directory=True)
            with self.assertRaisesRegex(FeatureGraphError, "must be absent"):
                bundle.save(linked)
            self.assertEqual(list(victim.iterdir()), [])

            from brainc.bio_graph import _save_staged_artifact as graph_stage_save
            from brainc.v2._common import V2Error

            calls = 0

            def fail_after_one(
                payload: dict, role: str, staging: Path, descriptor: int
            ) -> None:
                nonlocal calls
                calls += 1
                if calls == 2:
                    raise V2Error("injected artifact write failure")
                graph_stage_save(payload, role, staging, descriptor)

            destination = root / "transactional"
            with mock.patch(
                "brainc.bio_graph._save_staged_artifact",
                side_effect=fail_after_one,
            ):
                with self.assertRaisesRegex(FeatureGraphError, "cannot save"):
                    bundle.save(destination)
            self.assertFalse(destination.exists())
            self.assertEqual(list(root.glob(f".{destination.name}.*")), [])

            calls = 0
            portable = root / "portable-transactional"
            with mock.patch("brainc.bio_graph._HAS_DIRECTORY_DESCRIPTOR", False):
                with mock.patch(
                    "brainc.bio_graph._save_staged_artifact",
                    side_effect=fail_after_one,
                ):
                    with self.assertRaisesRegex(FeatureGraphError, "cannot save"):
                        bundle.save(portable)
            self.assertFalse(portable.exists())
            self.assertEqual(list(root.glob(f".{portable.name}.*")), [])

            durable = root / "unsupported-directory-fsync"
            unsupported = OSError(errno.EINVAL, "directory fsync unsupported")
            with mock.patch(
                "brainc.bio_graph._directory_fsync", side_effect=unsupported
            ):
                paths = bundle.save(durable)
            self.assertTrue(paths["backend_spec"].is_file())
            self.assertTrue(paths["bundle"].is_file())

    def test_bundle_cleanup_cannot_follow_a_substituted_staging_path(self) -> None:
        if not bio_graph_module._HAS_DIRECTORY_DESCRIPTOR:
            self.skipTest("staging substitution defense requires directory descriptors")
        fasta = b">x\nACGTACGT\n"
        gff3 = b"##gff-version 3\nx\t.\tgene\t1\t8\t.\t+\t.\tID=g\n"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, _, _, bundle = self._chain(root, fasta, gff3)
            victim = root / "victim"
            victim.mkdir()
            sentinel = victim / "backend_spec.json"
            sentinel.write_text("do not remove", encoding="utf-8")

            from brainc.bio_graph import _save_staged_artifact as graph_stage_save

            displaced: Path | None = None
            staging_link: Path | None = None
            calls = 0

            def substitute_before_first_write(
                payload: dict, role: str, staging: Path, descriptor: int
            ) -> None:
                nonlocal calls, displaced, staging_link
                calls += 1
                if calls == 1:
                    staging_link = staging
                    displaced = root / "displaced-staging"
                    staging_link.rename(displaced)
                    staging_link.symlink_to(victim, target_is_directory=True)
                graph_stage_save(payload, role, staging, descriptor)

            destination = root / "substitution-target"
            with mock.patch(
                "brainc.bio_graph._save_staged_artifact",
                side_effect=substitute_before_first_write,
            ):
                with self.assertRaisesRegex(FeatureGraphError, "changed identity"):
                    bundle.save(destination)
            self.assertFalse(destination.exists())
            self.assertEqual(sentinel.read_text(encoding="utf-8"), "do not remove")
            self.assertIsNotNone(displaced)
            self.assertEqual(list(displaced.iterdir()), [])
            self.assertIsNotNone(staging_link)
            self.assertTrue(staging_link.is_symlink())

    def test_bundle_final_publish_detects_staging_path_substitution(self) -> None:
        if not bio_graph_module._HAS_DIRECTORY_DESCRIPTOR:
            self.skipTest("staging substitution defense requires directory descriptors")
        fasta = b">x\nACGTACGT\n"
        gff3 = b"##gff-version 3\nx\t.\tgene\t1\t8\t.\t+\t.\tID=g\n"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, _, _, bundle = self._chain(root, fasta, gff3)
            destination = root / "final-substitution"
            real_fsync = os.fsync
            displaced = root / "displaced-final-staging"
            substituted: Path | None = None
            calls = 0

            def substitute_at_first_directory_sync(descriptor: int) -> None:
                nonlocal calls, substituted
                calls += 1
                if calls == 1:
                    candidates = list(root.glob(f".{destination.name}.*"))
                    self.assertEqual(len(candidates), 1)
                    substituted = candidates[0]
                    substituted.rename(displaced)
                    substituted.mkdir()
                    (substituted / "ATTACKER").write_text(
                        "unowned", encoding="utf-8"
                    )
                real_fsync(descriptor)

            with mock.patch(
                "brainc.bio_graph._directory_fsync",
                side_effect=substitute_at_first_directory_sync,
            ):
                with self.assertRaisesRegex(
                    FeatureGraphError,
                    "staging identity changed during publication",
                ):
                    bundle.save(destination)
            self.assertEqual(
                (destination / "ATTACKER").read_text(encoding="utf-8"),
                "unowned",
            )
            self.assertEqual(list(displaced.iterdir()), [])

    def test_bundle_destination_race_is_no_clobber(self) -> None:
        if not bio_graph_module._HAS_DIRECTORY_DESCRIPTOR:
            self.skipTest("destination reservation requires directory descriptors")
        fasta = b">x\nACGTACGT\n"
        gff3 = b"##gff-version 3\nx\t.\tgene\t1\t8\t.\t+\t.\tID=g\n"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, _, _, bundle = self._chain(root, fasta, gff3)
            destination = root / "appearing-destination"
            real_fsync = os.fsync
            calls = 0

            def create_destination_at_first_directory_sync(
                descriptor: int,
            ) -> None:
                nonlocal calls
                calls += 1
                if calls == 1:
                    destination.mkdir()
                    (destination / "SENTINEL").write_text(
                        "preserve", encoding="utf-8"
                    )
                real_fsync(descriptor)

            with mock.patch(
                "brainc.bio_graph._directory_fsync",
                side_effect=create_destination_at_first_directory_sync,
            ):
                with self.assertRaisesRegex(
                    FeatureGraphError,
                    "output path appeared during staging",
                ):
                    bundle.save(destination)
            self.assertEqual(
                (destination / "SENTINEL").read_text(encoding="utf-8"),
                "preserve",
            )
            self.assertEqual(list(root.glob(f".{destination.name}.*")), [])

    def test_portable_bundle_publish_uses_absent_destination_rename(self) -> None:
        fasta = b">x\nACGTACGT\n"
        gff3 = b"##gff-version 3\nx\t.\tgene\t1\t8\t.\t+\t.\tID=g\n"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, _, _, bundle = self._chain(root, fasta, gff3)
            destination = root / "portable-publish"
            with mock.patch("brainc.bio_graph._HAS_DIRECTORY_DESCRIPTOR", False):
                paths = bundle.save(destination)
            self.assertTrue(paths["backend_spec"].is_file())
            self.assertTrue(paths["bundle"].is_file())

    def test_bundle_children_are_defensive_immutable_snapshots(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, _, _, bundle = self._chain(
                root,
                b">x\nACGT\n",
                b"##gff-version 3\nx\t.\tgene\t1\t4\t.\t+\t.\tID=g\n",
            )
            expected_format = bundle.semantics["format"]
            exposed = bundle.semantics
            exposed["format"] = "corrupted"
            exposed["semantics"]["bio_ir"]["features"][0]["type"] = "forged"
            self.assertEqual(bundle.semantics["format"], expected_format)
            self.assertNotEqual(
                bundle.semantics["semantics"]["bio_ir"]["features"][0]["type"],
                "forged",
            )
            output = root / "immutable-bundle"
            paths = bundle.save(output)
            saved = json.loads(paths["backend_semantics"].read_bytes())
            self.assertEqual(saved["format"], expected_format)

    def test_actual_ncbi_annotation_materializes_exact_structural_runtime_state(self) -> None:
        if not RUNTIME_ROOT.is_dir():
            self.skipTest("sibling Epigenesis runtime checkout is unavailable")
        fasta_path = ROOT / "tests/data/J02482.1.fasta"
        gff3_path = ROOT / "tests/data/J02482.1.gff3"
        gff3 = gff3_path.read_bytes()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            collection = SequenceCollectionCompiler().compile_file(fasta_path)
            source_path = root / "source.json"
            collection.save(source_path)
            bio = GFF3Compiler().compile_file(gff3_path, collection)
            bundle = compile_feature_graph(
                source_path, bio, gff3_source=gff3
            )
            paths = bundle.save(root / "bundle")

            sys.path.insert(0, str(RUNTIME_ROOT))
            try:
                from epirun.contracts import (
                    load_development_module,
                    load_target_contract,
                )
                from epirun.development import develop

                target = load_target_contract(paths["target_contract"])
                module = load_development_module(
                    paths["development_module"], target
                )
                result = develop(target, module)
            finally:
                if sys.path[0] == str(RUNTIME_ROOT):
                    sys.path.pop(0)

            ir = bio.to_dict()["bio_ir"]
            self.assertEqual(len(ir["relationships"]), 4)
            self.assertEqual(result.state.unit_count, len(ir["features"]))
            self.assertEqual(result.state.edge_count, len(ir["relationships"]))
            self.assertEqual(len(result.state.attachments), 0)
            self.assertEqual(len(result.state.ports), 0)
            self.assertEqual(
                result.final_state_sha256,
                "95cd313a51e6c43c29582c6834db0452430468a12975e0429434c329c313c2aa",
            )
            feature_state = result.state.unit_set("create.features")
            self.assertEqual(
                {field.id: field.scalar(0) for field in feature_state.fields},
                {
                    "attribute_count": 6,
                    "declared_id": True,
                    "interval_count": 1,
                    "location_max_end": 5386,
                    "location_min_start": 0,
                    "phase_mask": 0,
                    "segment_bases": 5386,
                    "segment_count": 1,
                    "seqid_code": 0,
                    "source_code": 1,
                    "strand_code": 0,
                    "type_code": 4,
                    "wraps_origin": False,
                },
            )
            parent_edges = result.state.edge_set("create.relationship.parent")
            self.assertEqual(parent_edges.ids, (0, 1, 2, 3))
            self.assertEqual(parent_edges.source_ids, (6, 9, 12, 13))
            self.assertEqual(parent_edges.target_ids, (22, 21, 21, 22))
            self.assertEqual(
                [
                    parent_edges.field("kind_code").scalar(index)
                    for index in range(len(parent_edges.ids))
                ],
                [1, 1, 1, 1],
            )
            derives_edges = result.state.edge_set(
                "create.relationship.derives-from"
            )
            self.assertEqual(
                (derives_edges.ids, derives_edges.source_ids, derives_edges.target_ids),
                ((), (), ()),
            )
            self.assertEqual(
                bundle.compilation_record["inputs"]["bio_ir"]["artifact_sha256"],
                bio.to_dict()["artifact_sha256"],
            )
            validate_feature_graph_bundle(
                bundle, source_path, bio, gff3_source=gff3
            )


if __name__ == "__main__":
    unittest.main()

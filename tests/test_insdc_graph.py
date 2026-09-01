from __future__ import annotations

import base64
import copy
import errno
import json
from pathlib import Path
import stat
import struct
import tempfile
import unittest
from unittest import mock

from brainc._canonical import digest
import brainc.insdc as insdc_module
from brainc.insdc import GenBankCompiler
import brainc.insdc_graph as graph_module
from brainc.insdc_graph import (
    BACKEND_ID,
    BACKEND_SPEC,
    BACKEND_SPEC_SHA256,
    CHILD_ROLES,
    FEATURE_FIELD_SPECS,
    compile_insdc_graph,
)
from brainc.v2._common import pretty_bytes


ROOT = Path(__file__).resolve().parents[1]

U49845_CLOSURE_SHA256 = {
    "bundle": "fcb748c89e4a9b81f06c08cea2e43cc28cb9ab5394a6cdfc08cd5a1f1702b4de",
    "backend_semantics": "ee8dbe0b98db3bd3866318717d3bf80cdb5ac0a026167b9f832ea4f9caa673b6",
    "backend_spec": "a8f2494dddd259f7f41604001f5dcae9e19b0e491eac7c13aec298f2f99a5da4",
    "compilation_record": "7b63de123eed7e349ce65114919418c210a8b0aece350df6c775d7e0eba83527",
    "development_module": "46764eae6feb0f018ff0611ef074b37181d442de2a8db04d86f1a57d9ea0ac71",
    "lowering_policy": "0eee90dea3d62eb44ce9ab1b1704ff3178d944048df2ad7714e57e4094a8c070",
    "prediction_request": "7a8d5361b5a3b5635232bfc8b32f98a6184ffcbaa7dbbeadbfbeb58cfd678939",
    "prediction_response": "2026a9d667bbe6b2bf0577a0f11a4b1d2a0994de2d9b2feec765beffd464d4bc",
    "provider_manifest": "d4d6541d1e917addeaaefb0566cd554ea28eef048c9953864a3c4f450329739b",
    "target_contract": "2c7d7cc926d3e0125a867e1d3ab777105311e8dec67b95cf282c6d95e27c6748",
}


def _origin(sequence: str) -> list[str]:
    lines: list[str] = []
    for offset in range(0, len(sequence), 60):
        part = sequence[offset : offset + 60].lower()
        groups = " ".join(
            part[index : index + 10] for index in range(0, len(part), 10)
        )
        lines.append(f"{offset + 1:>9} {groups}")
    return lines


def _record(
    accession: str,
    sequence: str,
    *,
    features: list[str] | None = None,
) -> bytes:
    record_id = f"{accession}.1"
    feature_lines = [
        f"     {'source':<15} 1..{len(sequence)}",
        '                     /organism="Arbitrary organism"',
        '                     /mol_type="genomic DNA"',
        *(features or []),
    ]
    lines = [
        f"LOCUS       GENERATED_{accession} {len(sequence)} bp DNA linear SYN 31-AUG-2026",
        "DEFINITION  Generated adapter input.",
        f"ACCESSION   {accession}",
        f"VERSION     {record_id}",
        "KEYWORDS    .",
        "SOURCE      synthetic construct",
        "  ORGANISM  Arbitrary organism",
        "            artificial sequences.",
        f"REFERENCE   1  (bases 1 to {len(sequence)})",
        "  AUTHORS   Example,A.",
        "  TITLE     Direct Submission",
        "  JOURNAL   Submitted (31-AUG-2026)",
        "FEATURES             Location/Qualifiers",
        *feature_lines,
        "ORIGIN",
        *_origin(sequence),
        "//",
    ]
    return ("\n".join(lines) + "\n").encode("ascii")


def _synthetic_source():
    raw = _record(
        "AX000001",
        "ACGT" * 5,
        features=[
            "     misc_feature    join(<2..>4,complement(6..8),",
            "                     9^10,12.15,AX000002.1:3..5,",
            "                     ZZ000001.1:7..9)",
            '                     /note="one"',
            '                     /note="two"',
            "                     /pseudo",
        ],
    ) + _record(
        "AX000002",
        "TGCA" * 5,
        features=["     misc_feature    complement(1..2)"],
    )
    return GenBankCompiler().compile_bytes(raw)


def _u64_values(tensor: dict) -> list[int]:
    storage = tensor["storage"]
    raw = base64.b64decode(storage["data"], validate=True)
    if not raw:
        return []
    return list(struct.unpack("<" + "Q" * (len(raw) // 8), raw))


def _tensors(bundle) -> dict[str, dict]:
    return {
        tensor["id"]: tensor
        for tensor in bundle.development_module["module"]["tensors"]
    }


class GenBankFeatureStateTests(unittest.TestCase):
    def test_u49845_replays_once_and_preserves_the_exact_artifact_closure(self) -> None:
        source = GenBankCompiler().compile_bytes(
            (ROOT / "tests" / "data" / "U49845.1.gb").read_bytes()
        )
        with tempfile.TemporaryDirectory() as temporary:
            source_path = Path(temporary) / "source.json"
            source.save(source_path)
            for value in (source, source.to_dict(), source_path):
                with self.subTest(input_type=type(value).__name__), mock.patch.object(
                    insdc_module,
                    "_compile_parts",
                    wraps=insdc_module._compile_parts,
                ) as replay:
                    bundle = compile_insdc_graph(value)
                self.assertEqual(replay.call_count, 1)
                closure = {
                    "bundle": bundle.to_dict()["artifact_sha256"],
                    **{
                        role: artifact["artifact_sha256"]
                        for role, artifact in bundle.artifacts.items()
                    },
                }
                self.assertEqual(closure, U49845_CLOSURE_SHA256)

    def test_u49845_compiles_six_units_without_edges(self) -> None:
        source = GenBankCompiler().compile_bytes(
            (ROOT / "tests" / "data" / "U49845.1.gb").read_bytes()
        )
        bundle = compile_insdc_graph(source)

        budgets = bundle.development_module["module"]["budgets"]
        self.assertEqual(
            budgets,
            {
                "operations": 1,
                "tensor_bytes": 8 * (1 + 11 * 6),
                "units": 6,
                "edges": 0,
                "attachments": 0,
            },
        )
        operations = bundle.development_module["module"]["entrypoint"][
            "operations"
        ]
        self.assertEqual(len(operations), 1)
        self.assertEqual(operations[0]["id"], "create.features")
        self.assertEqual(bundle.semantics["counts"], {"units": 6, "edges": 0})
        self.assertEqual(bundle.semantics["edge_order"], [])
        self.assertEqual(
            [entry["value"] for entry in bundle.semantics["dictionaries"]["keys"]],
            ["CDS", "gene", "source"],
        )
        self.assertEqual(
            [
                entry["value"]
                for entry in bundle.semantics["dictionaries"]["records"]
            ],
            ["U49845.1"],
        )
        expected_order = [
            feature["feature_id"]
            for record in source.to_dict()["bio_ir"]["records"]
            for feature in record["features"]
        ]
        self.assertEqual(bundle.semantics["unit_order"], expected_order)
        self.assertEqual(_u64_values(_tensors(bundle)["t.feature.count"]), [6])

    def test_all_field_formulas_on_arbitrary_structural_features(self) -> None:
        bundle = compile_insdc_graph(_synthetic_source())
        tensors = _tensors(bundle)
        expected = {
            "fuzzy_boundary_count": [0, 2, 0, 0],
            "key_code": [1, 0, 1, 0],
            "kind_mask": [1, 7, 1, 1],
            "ordinal": [1, 2, 1, 2],
            "orientation_mask": [1, 3, 1, 2],
            "qualifier_count": [2, 3, 2, 0],
            "record_code": [0, 0, 1, 1],
            "remote_segment_count": [0, 2, 0, 0],
            "segment_count": [1, 6, 1, 1],
            "segment_extent_sum": [20, 16, 20, 2],
            "unresolved_segment_count": [0, 1, 0, 0],
        }
        self.assertEqual(
            [field["id"] for field in FEATURE_FIELD_SPECS], list(expected)
        )
        for field, values in expected.items():
            tensor = tensors[f"t.feature.{field}"]
            self.assertEqual(tensor["type"], {"dtype": "u64", "shape": [4]})
            self.assertIsNone(tensor["unit"])
            self.assertEqual(tensor["axes"], [None])
            self.assertEqual(_u64_values(tensor), values)

        self.assertEqual(
            [entry["value"] for entry in bundle.semantics["dictionaries"]["keys"]],
            ["misc_feature", "source"],
        )
        self.assertEqual(
            [
                entry["value"]
                for entry in bundle.semantics["dictionaries"]["records"]
            ],
            ["AX000001.1", "AX000002.1"],
        )

    def test_bundle_save_uses_the_shared_publication_owner(self) -> None:
        bundle = compile_insdc_graph(_synthetic_source())
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / "state"
            shared = graph_module.publish_directory
            with mock.patch.object(
                graph_module,
                "publish_directory",
                wraps=shared,
            ) as publisher:
                paths = bundle.save(destination)
            publisher.assert_called_once()
            entries = publisher.call_args.args[1]
            self.assertEqual(len(entries), 10)
            self.assertEqual(set(entries), {path.name for path in paths.values()})

    def test_bundle_is_deterministic_reference_only_and_bounded(self) -> None:
        source = _synthetic_source()
        first = compile_insdc_graph(source)
        second = compile_insdc_graph(source.to_dict())
        self.assertEqual(first.to_dict(), second.to_dict())
        self.assertEqual(first.artifacts, second.artifacts)

        index = first.to_dict()
        self.assertEqual(set(index["artifacts"]), set(CHILD_ROLES))
        self.assertTrue(
            all(
                set(reference) == {"artifact_sha256"}
                for reference in index["artifacts"].values()
            )
        )
        self.assertNotIn("bio_ir", index)
        self.assertNotIn("sequence_collection", index)
        self.assertNotIn("bio_ir", first.semantics)
        self.assertNotIn("sequence_collection", first.semantics)
        for role, artifact in first.artifacts.items():
            self.assertLessEqual(len(pretty_bytes(artifact)), graph_module.MAX_CHILD_BYTES)
            self.assertEqual(
                index["artifacts"][role]["artifact_sha256"],
                artifact["artifact_sha256"],
            )
        self.assertLessEqual(len(pretty_bytes(index)), graph_module.MAX_CHILD_BYTES)
        self.assertEqual(index["backend"]["id"], BACKEND_ID)
        self.assertEqual(index["backend"]["spec_sha256"], BACKEND_SPEC_SHA256)
        self.assertEqual(BACKEND_SPEC_SHA256, digest(BACKEND_SPEC))
        self.assertEqual(
            BACKEND_SPEC["input"],
            {
                "format": "brainc.bio.insdc-genbank-ir",
                "version": 2,
                "profile": "genbank-273-traditional-dna-physical-structural/v2",
                "authority_manifest_sha256": (
                    "6dc3658a5f7d6a774ce9cd08ff6dc5327cbef8a9951c7679fd9a34650f912bf9"
                ),
                "coordinate_system": "0-based-half-open",
                "binding": "whole-artifact-sha256-and-bio-ir-sha256",
            },
        )
        projection = BACKEND_SPEC["projection"]
        self.assertEqual(
            projection,
            {
                "unit": "one-per-ordered-insdc-feature-table-entry",
                "columns": {"count": 11, "kind": "structural-summary"},
                "semantic_authority": {
                    "artifact": "bound-genbank-bio-ir",
                    "role": "source-dependency",
                },
                "dictionary_codes": "bundle-local-and-may-renumber-with-membership",
                "masks": "discard-multiplicity-order-and-segment-association",
                "fuzzy_boundary_count": "excludes-within-uncertainty",
                "unresolved_segment_count": (
                    "reference-resolution-and-bounds-status"
                ),
                "edges": "none-inferred-or-emitted",
            },
        )
        self.assertTrue(
            all(
                dictionary["stability"]
                == "bundle-local; may-renumber-when-membership-changes"
                for dictionary in BACKEND_SPEC["dictionaries"]
            )
        )
        field_specs = {field["id"]: field for field in FEATURE_FIELD_SPECS}
        self.assertEqual(
            field_specs["segment_extent_sum"]["formula"],
            {
                "op": "sum",
                "over": "feature.location.segments",
                "value": {
                    "op": "subtract",
                    "left": {"path": "segment.end"},
                    "right": {"path": "segment.start"},
                },
                "coordinate_system": "normalized-0-based-half-open",
                "between_contribution": 0,
                "uncertain_point_contribution": "normalized-envelope-width",
                "overlap_policy": "count-each-segment-extent-with-multiplicity",
            },
        )

        encoded = json.dumps(index, sort_keys=True)
        self.assertNotIn("AX000001", encoded)
        self.assertNotIn("misc_feature", encoded)

        original = BACKEND_SPEC["capabilities"]["learning"]
        BACKEND_SPEC["capabilities"]["learning"] = True
        original_tensor_id = FEATURE_FIELD_SPECS[0]["tensor_id"]
        FEATURE_FIELD_SPECS[0]["tensor_id"] = "t.feature.tainted"
        original_kind_bits = graph_module.KIND_BITS
        original_orientation_bits = graph_module.ORIENTATION_BITS
        graph_module.KIND_BITS = (("between", 8), ("interval", 16), ("uncertain-point", 32))
        graph_module.ORIENTATION_BITS = ((-1, 8), (1, 16))
        try:
            unaffected = compile_insdc_graph(source)
        finally:
            BACKEND_SPEC["capabilities"]["learning"] = original
            FEATURE_FIELD_SPECS[0]["tensor_id"] = original_tensor_id
            graph_module.KIND_BITS = original_kind_bits
            graph_module.ORIENTATION_BITS = original_orientation_bits
        self.assertFalse(unaffected.backend_spec["capabilities"]["learning"])
        self.assertIn(original_tensor_id, _tensors(unaffected))
        self.assertNotIn("t.feature.tainted", _tensors(unaffected))
        self.assertEqual(
            _u64_values(_tensors(unaffected)["t.feature.kind_mask"]),
            [1, 7, 1, 1],
        )
        self.assertEqual(
            _u64_values(_tensors(unaffected)["t.feature.orientation_mask"]),
            [1, 3, 1, 2],
        )
        copy_out = unaffected.artifacts
        copy_out["backend_semantics"]["counts"]["units"] = 99
        self.assertEqual(unaffected.semantics["counts"]["units"], 4)

    def test_bundle_directory_contains_only_index_and_sealed_children(self) -> None:
        bundle = compile_insdc_graph(_synthetic_source())
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / "state"
            paths = bundle.save(destination)
            self.assertEqual(
                {path.name for path in destination.iterdir()},
                {
                    graph_module.BUNDLE_FILENAME,
                    *(filename for _, filename in graph_module.CHILD_FILENAMES),
                },
            )
            self.assertEqual(
                json.loads(paths["bundle"].read_text(encoding="utf-8")),
                bundle.to_dict(),
            )
            for role in CHILD_ROLES:
                self.assertEqual(
                    json.loads(paths[role].read_text(encoding="utf-8")),
                    bundle.artifact(role),
                )
            with self.assertRaisesRegex(
                graph_module.INSDCGraphError, "output path"
            ):
                bundle.save(destination)

            empty = Path(temporary) / "empty"
            empty.mkdir()
            with self.assertRaisesRegex(
                graph_module.INSDCGraphError, "output path"
            ):
                bundle.save(empty)
            self.assertEqual(list(empty.iterdir()), [])

    @unittest.skipUnless(
        graph_module._HAS_DIRECTORY_DESCRIPTOR,
        "directory descriptors are unavailable",
    )
    def _legacy_bundle_publication_is_staged_and_retryable_after_write_failure(self) -> None:
        bundle = compile_insdc_graph(_synthetic_source())
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            destination = Path(temporary) / "state"
            original_save = graph_module._save_staged_artifact
            calls = 0

            def fail_after_first_child(payload, filename, staging, descriptor):
                nonlocal calls
                calls += 1
                self.assertFalse(destination.exists())
                if calls == 2:
                    raise graph_module.V2Error("injected staged write failure")
                return original_save(payload, filename, staging, descriptor)

            with mock.patch.object(
                graph_module,
                "_save_staged_artifact",
                side_effect=fail_after_first_child,
            ):
                with self.assertRaisesRegex(
                    graph_module.INSDCGraphError, "injected staged write failure"
                ):
                    bundle.save(destination)
            self.assertFalse(destination.exists())
            self.assertEqual(list(root.glob(f".{destination.name}.*")), [])
            self.assertTrue(bundle.save(destination)["bundle"].is_file())

    @unittest.skipUnless(
        graph_module._HAS_DIRECTORY_DESCRIPTOR,
        "directory descriptors are unavailable",
    )
    def _legacy_bundle_publication_cleans_up_after_child_fstat_failure(self) -> None:
        bundle = compile_insdc_graph(_synthetic_source())
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            destination = root / "state"
            original_fstat = graph_module.os.fstat
            calls = 0

            def fail_first_child_fstat(descriptor):
                nonlocal calls
                metadata = original_fstat(descriptor)
                if stat.S_ISREG(metadata.st_mode):
                    calls += 1
                    if calls == 1:
                        raise OSError(errno.EIO, "injected child fstat failure")
                return metadata

            with mock.patch.object(
                graph_module.os,
                "fstat",
                side_effect=fail_first_child_fstat,
            ):
                with self.assertRaisesRegex(
                    graph_module.INSDCGraphError,
                    "injected child fstat failure",
                ):
                    bundle.save(destination)
            self.assertFalse(destination.exists())
            self.assertEqual(list(root.glob(f".{destination.name}.*")), [])
            self.assertTrue(bundle.save(destination)["bundle"].is_file())

    @unittest.skipUnless(
        graph_module._HAS_DIRECTORY_DESCRIPTOR,
        "directory descriptors are unavailable",
    )
    def _legacy_bundle_publication_detects_child_swap_at_directory_sync(self) -> None:
        bundle = compile_insdc_graph(_synthetic_source())
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            destination = root / "state"
            displaced = root / "displaced-child.json"
            original_fsync = graph_module.os.fsync
            raced = False

            def swap_at_directory_sync(descriptor):
                nonlocal raced
                metadata = graph_module.os.fstat(descriptor)
                if stat.S_ISDIR(metadata.st_mode) and not raced:
                    candidates = list(root.glob(f".{destination.name}.*"))
                    self.assertEqual(len(candidates), 1)
                    child = candidates[0] / graph_module.CHILD_FILENAMES[0][1]
                    child.rename(displaced)
                    child.write_text('{"attacker":true}\n', encoding="utf-8")
                    raced = True
                return original_fsync(descriptor)

            with mock.patch.object(
                graph_module.os,
                "fsync",
                side_effect=swap_at_directory_sync,
            ):
                with self.assertRaisesRegex(
                    graph_module.INSDCGraphError,
                    "changed during publication",
                ):
                    bundle.save(destination)
            self.assertTrue(raced)
            self.assertFalse(destination.exists())
            self.assertTrue(displaced.is_file())
            retained = list(root.glob(f".{destination.name}.*"))
            self.assertEqual(len(retained), 1)
            self.assertEqual(
                (retained[0] / graph_module.CHILD_FILENAMES[0][1]).read_text(
                    encoding="utf-8"
                ),
                '{"attacker":true}\n',
            )

    @unittest.skipUnless(
        graph_module._HAS_DIRECTORY_DESCRIPTOR,
        "directory descriptors are unavailable",
    )
    def _legacy_bundle_destination_race_preserves_appearing_path(self) -> None:
        bundle = compile_insdc_graph(_synthetic_source())
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            destination = root / "state"
            original_fsync = graph_module.os.fsync
            raced = False

            def create_destination_at_directory_sync(descriptor):
                nonlocal raced
                metadata = graph_module.os.fstat(descriptor)
                if stat.S_ISDIR(metadata.st_mode) and not raced:
                    destination.mkdir()
                    (destination / "SENTINEL").write_text(
                        "preserve",
                        encoding="utf-8",
                    )
                    raced = True
                return original_fsync(descriptor)

            with mock.patch.object(
                graph_module.os,
                "fsync",
                side_effect=create_destination_at_directory_sync,
            ):
                with self.assertRaisesRegex(
                    graph_module.INSDCGraphError,
                    "appeared during staging",
                ):
                    bundle.save(destination)
            self.assertTrue(raced)
            self.assertEqual(
                (destination / "SENTINEL").read_text(encoding="utf-8"),
                "preserve",
            )
            self.assertEqual(list(root.glob(f".{destination.name}.*")), [])


if __name__ == "__main__":
    unittest.main()

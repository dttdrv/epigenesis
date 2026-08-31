from __future__ import annotations

from copy import deepcopy
import gzip
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from brainc._io import (
    BoundedIOError,
    atomic_write_file,
    load_json_object,
    pretty_json_bytes,
    read_regular_file,
    unwrap_gzip,
)
from brainc.sequence import SequenceCompiler, SequenceCompilerError, load_sequence_artifact
from brainc.sequence_collection import (
    SequenceCollectionCompiler,
    SequenceCollectionError,
)
from brainc._canonical import digest
from brainc.v2._common import V2Error
from brainc.v2.provider import load_source
from brainc.validator_v2 import ValidationError, validate_source_artifact


def _reseal(value: dict) -> dict:
    result = deepcopy(value)
    result.pop("artifact_sha256", None)
    result["artifact_sha256"] = digest(result)
    return result


class SourceLimitTests(unittest.TestCase):
    def test_regular_file_exact_boundary_and_plus_one(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "source.bin"
            source.write_bytes(b"12345678")
            self.assertEqual(
                read_regular_file(source, maximum_bytes=8, label="test source"),
                b"12345678",
            )
            source.write_bytes(b"123456789")
            with self.assertRaisesRegex(BoundedIOError, "exceeds byte limit 8"):
                read_regular_file(source, maximum_bytes=8, label="test source")

    @unittest.skipUnless(hasattr(os, "symlink"), "symbolic links are unavailable")
    def test_sequence_file_ingress_rejects_symlinks(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / "target.fa"
            target.write_bytes(b">x\nACGT\n")
            linked = root / "linked.fa"
            linked.symlink_to(target)
            with self.assertRaisesRegex(SequenceCompilerError, "regular non-linked file"):
                SequenceCompiler().compile_file(linked)
            with self.assertRaisesRegex(SequenceCollectionError, "regular non-linked file"):
                SequenceCollectionCompiler().compile_file(linked)

            with self.assertRaisesRegex(SequenceCompilerError, "regular non-linked file"):
                SequenceCompiler().compile_file(root)

    def test_path_replacement_between_inspection_and_open_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.bin"
            replacement = root / "replacement.bin"
            source.write_bytes(b"original")
            replacement.write_bytes(b"replaced")
            real_open = os.open

            def replace_then_open(path: str | bytes | os.PathLike[str], flags: int) -> int:
                os.replace(replacement, source)
                return real_open(path, flags)

            with mock.patch("brainc._io.os.open", side_effect=replace_then_open):
                with self.assertRaisesRegex(BoundedIOError, "changed while being opened"):
                    read_regular_file(source, maximum_bytes=8, label="test source")

    def test_gzip_streaming_exact_decompressed_boundary_and_plus_one(self) -> None:
        exact = gzip.compress(b"A" * 32, mtime=0)
        logical, wrapped = unwrap_gzip(
            exact,
            maximum_input_bytes=len(exact),
            maximum_decompressed_bytes=32,
        )
        self.assertTrue(wrapped)
        self.assertEqual(logical, b"A" * 32)

        over = gzip.compress(b"A" * 33, mtime=0)
        with self.assertRaisesRegex(BoundedIOError, "decompressed byte limit 32"):
            unwrap_gzip(
                over,
                maximum_input_bytes=len(over),
                maximum_decompressed_bytes=32,
            )

    def test_high_ratio_gzip_is_stopped_at_logical_ceiling(self) -> None:
        compressed = gzip.compress(b"A" * (1024 * 1024), mtime=0)
        self.assertLess(len(compressed), 2048)
        with self.assertRaisesRegex(BoundedIOError, "decompressed byte limit 65536"):
            unwrap_gzip(
                compressed,
                maximum_input_bytes=2048,
                maximum_decompressed_bytes=65536,
            )

    def test_gzip_multistream_and_trailing_data_policy_matches_prior_ingress(self) -> None:
        first = gzip.compress(b"one", mtime=0)
        second = gzip.compress(b"two", mtime=0)
        logical, wrapped = unwrap_gzip(
            first + second,
            maximum_input_bytes=len(first + second),
            maximum_decompressed_bytes=6,
        )
        self.assertTrue(wrapped)
        self.assertEqual(logical, b"onetwo")

        padded, _ = unwrap_gzip(
            first + b"\0\0",
            maximum_input_bytes=len(first) + 2,
            maximum_decompressed_bytes=3,
        )
        self.assertEqual(padded, b"one")
        with self.assertRaisesRegex(BoundedIOError, "malformed gzip source"):
            unwrap_gzip(
                first + b"junk",
                maximum_input_bytes=len(first) + 4,
                maximum_decompressed_bytes=3,
            )
        with self.assertRaisesRegex(BoundedIOError, "byte limit"):
            unwrap_gzip(
                first,
                maximum_input_bytes=len(first) - 1,
                maximum_decompressed_bytes=3,
            )

    def test_bounded_json_exact_boundaries_duplicates_and_i_json(self) -> None:
        self.assertEqual(load_json_object(b"{}", "test", maximum_bytes=2), {})
        with self.assertRaisesRegex(BoundedIOError, "JSON byte limit 1"):
            load_json_object(b"{}", "test", maximum_bytes=1)

        self.assertEqual(
            load_json_object(b'{"a":"xy"}', "test", maximum_string_bytes=2),
            {"a": "xy"},
        )
        with self.assertRaisesRegex(BoundedIOError, "JSON string byte limit 1"):
            load_json_object(b'{"a":"xy"}', "test", maximum_string_bytes=1)
        with self.assertRaisesRegex(BoundedIOError, "duplicate JSON key"):
            load_json_object(b'{"a":1,"a":2}', "test")
        with self.assertRaisesRegex(BoundedIOError, "safe range"):
            load_json_object(b'{"a":9007199254740992}', "test")

    def test_bounded_json_member_and_depth_exact_boundaries(self) -> None:
        self.assertEqual(
            load_json_object(b'{"a":[]}', "test", maximum_members=1, maximum_depth=1),
            {"a": []},
        )
        with self.assertRaisesRegex(BoundedIOError, "member limit 0"):
            load_json_object(b'{"a":[]}', "test", maximum_members=0)
        with self.assertRaisesRegex(BoundedIOError, "depth limit 0"):
            load_json_object(b'{"a":[]}', "test", maximum_depth=0)

    def test_pretty_output_exact_byte_boundary_and_plus_one(self) -> None:
        raw = pretty_json_bytes({}, ensure_ascii=False)
        self.assertEqual(
            pretty_json_bytes({}, ensure_ascii=False, maximum_bytes=len(raw)),
            raw,
        )
        with self.assertRaisesRegex(BoundedIOError, "serialized output exceeds"):
            pretty_json_bytes({}, ensure_ascii=False, maximum_bytes=len(raw) - 1)

    def test_atomic_output_exact_boundary_and_failed_replace_preserves_target(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            destination = root / "artifact.json"
            atomic_write_file(destination, b"1234", maximum_bytes=4, label="artifact")
            self.assertEqual(destination.read_bytes(), b"1234")
            with self.assertRaisesRegex(BoundedIOError, "byte limit 4"):
                atomic_write_file(destination, b"12345", maximum_bytes=4, label="artifact")
            self.assertEqual(destination.read_bytes(), b"1234")

            with mock.patch("brainc._io.os.replace", side_effect=OSError("injected")):
                with self.assertRaisesRegex(BoundedIOError, "cannot write artifact"):
                    atomic_write_file(destination, b"next", maximum_bytes=4, label="artifact")
            self.assertEqual(destination.read_bytes(), b"1234")
            self.assertEqual(list(root.glob(".artifact.json.*")), [])

    def test_atomic_output_parent_substitution_never_reports_false_success(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            parent = root / "output-parent"
            parent.mkdir()
            destination = parent / "artifact.json"
            displaced = root / "original-parent"
            original_replace = os.replace
            substituted = False

            def substitute_parent(source, target, *args, **kwargs):
                nonlocal substituted
                if not substituted:
                    substituted = True
                    original_replace(parent, displaced)
                    parent.mkdir()
                    destination.write_bytes(b"attacker")
                return original_replace(source, target, *args, **kwargs)

            with mock.patch(
                "brainc._io.os.replace", side_effect=substitute_parent
            ):
                with self.assertRaisesRegex(
                    BoundedIOError, "parent changed while publishing"
                ):
                    atomic_write_file(
                        destination,
                        b"intended",
                        maximum_bytes=8,
                        label="artifact",
                    )
            self.assertEqual(destination.read_bytes(), b"attacker")
            self.assertEqual(
                (displaced / "artifact.json").read_bytes(), b"intended"
            )

    @unittest.skipUnless(hasattr(os, "symlink"), "symbolic links are unavailable")
    def test_sequence_saves_reject_symlink_and_nonregular_outputs(self) -> None:
        sequence = SequenceCompiler().compile_bytes(b">x\nACGT\n")
        collection = SequenceCollectionCompiler().compile_fasta_bytes(b">x\nACGT\n")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / "target.json"
            target.write_bytes(b"sentinel")
            sequence_link = root / "sequence-link.json"
            sequence_link.symlink_to(target)
            with self.assertRaisesRegex(SequenceCompilerError, "regular non-linked file"):
                sequence.save(sequence_link)
            self.assertEqual(target.read_bytes(), b"sentinel")

            collection_link = root / "collection-link.json"
            collection_link.symlink_to(target)
            with self.assertRaisesRegex(SequenceCollectionError, "regular non-linked file"):
                collection.save(collection_link)
            self.assertEqual(target.read_bytes(), b"sentinel")

            nonregular = root / "directory-output"
            nonregular.mkdir()
            with self.assertRaisesRegex(SequenceCompilerError, "regular non-linked file"):
                sequence.save(nonregular)

    def test_sequence_producers_preflight_v2_output_limits(self) -> None:
        with mock.patch("brainc.sequence.MAX_STRING_BYTES", 4):
            with self.assertRaisesRegex(SequenceCompilerError, "emitted JSON string limit 4"):
                SequenceCompiler().compile_bytes(b">x\nAAAAA\n")
        with mock.patch("brainc.sequence.MAX_JSON_BYTES", 128):
            with self.assertRaisesRegex(SequenceCompilerError, "serialized output exceeds"):
                SequenceCompiler().compile_bytes(b">x\nA\n")
        with mock.patch("brainc.sequence_collection.MAX_JSON_BYTES", 128):
            with self.assertRaisesRegex(SequenceCollectionError, "serialized output exceeds"):
                SequenceCollectionCompiler().compile_fasta_bytes(b">x\nA\n")

    def test_raw_ingress_exact_boundary_and_plus_one(self) -> None:
        with mock.patch("brainc.sequence_collection.MAX_INPUT_BYTES", 4):
            artifact = SequenceCollectionCompiler().compile_raw(b"ACGT", "raw")
            self.assertEqual(artifact.members[0].artifact.sequence, "ACGT")
            with self.assertRaisesRegex(SequenceCollectionError, "byte limit 4"):
                SequenceCollectionCompiler().compile_raw(b"ACGTA", "raw")

    def test_oversized_and_linked_artifact_loads_are_rejected_before_json(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            artifact = root / "sequence.json"
            artifact.write_bytes(b"{" + b" " * 16 + b"}")
            with mock.patch("brainc.sequence.MAX_JSON_BYTES", 16):
                with self.assertRaisesRegex(SequenceCompilerError, "exceeds byte limit 16"):
                    load_sequence_artifact(artifact)

            if hasattr(os, "symlink"):
                linked = root / "linked.json"
                linked.symlink_to(artifact)
                with self.assertRaisesRegex(SequenceCompilerError, "regular non-linked file"):
                    load_sequence_artifact(linked)

    def test_compiled_collection_passes_public_loader_and_independent_replay(self) -> None:
        fasta = b">alpha description\nACGT\nTG\n>beta\nNNNN\n"
        collection = SequenceCollectionCompiler().compile_fasta_bytes(fasta)
        with tempfile.TemporaryDirectory() as directory:
            source_path = Path(directory) / "collection.json"
            collection.save(source_path)
            loaded, lengths = load_source(source_path)
        self.assertEqual(lengths, {"alpha": 6, "beta": 4})
        self.assertEqual(validate_source_artifact(loaded, raw_source=fasta), loaded)

    def test_public_loader_rejects_collection_member_context(self) -> None:
        fasta = b">alpha description\nACGT\nTG\n>beta\nNNNN\n"
        collection = SequenceCollectionCompiler().compile_fasta_bytes(fasta).to_dict()
        context = json.dumps(
            {
                "format": "brain01.sequence-context",
                "version": 1,
                "record_id": "alpha",
                "reference": None,
                "provenance": [
                    {
                        "id": "source-1",
                        "kind": "database",
                        "uri": "urn:example:source",
                        "version": "1",
                    }
                ],
            }
        ).encode("utf-8")
        nested = SequenceCompiler().compile_bytes(
            b">alpha description\nACGT\nTG\n", context
        ).to_dict()
        member = collection["members"][0]
        member["sequence_artifact"] = nested
        member["ir_sha256"] = nested["ir_sha256"]
        member["artifact_sha256"] = nested["artifact_sha256"]
        collection["collection_ir"]["members"][0]["ir_sha256"] = nested["ir_sha256"]
        collection["collection_ir_sha256"] = digest(collection["collection_ir"])
        forged = _reseal(collection)

        with self.assertRaisesRegex(ValidationError, "context not emitted"):
            validate_source_artifact(forged)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "forged.json"
            path.write_text(json.dumps(forged), encoding="utf-8")
            with self.assertRaisesRegex(V2Error, "context not emitted"):
                load_source(path)

    def test_public_loader_rejects_unbound_single_sequence_assertions(self) -> None:
        forged = SequenceCompiler().compile_bytes(b">single\nACGT\n").to_dict()
        forged["sequence_ir"]["provenance"] = [
            {
                "id": "source-1",
                "kind": "database",
                "uri": "urn:example:source",
                "version": "1",
            }
        ]
        forged["ir_sha256"] = digest(forged["sequence_ir"])
        forged = _reseal(forged)

        with self.assertRaisesRegex(ValidationError, "context assertions"):
            validate_source_artifact(forged)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "forged.json"
            path.write_text(json.dumps(forged), encoding="utf-8")
            with self.assertRaisesRegex(V2Error, "context assertions require"):
                load_source(path)

    def test_public_loader_rejects_overlapping_collection_source_maps(self) -> None:
        collection = SequenceCollectionCompiler().compile_fasta_bytes(
            b">one\nAC\n>two\nGT\n"
        ).to_dict()
        collection["members"][1]["input_source_map"]["sequence_segments"][0]["line"] = 2
        forged = _reseal(collection)

        with self.assertRaisesRegex(ValidationError, "overlaps"):
            validate_source_artifact(forged)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "forged.json"
            path.write_text(json.dumps(forged), encoding="utf-8")
            with self.assertRaisesRegex(V2Error, "overlaps"):
                load_source(path)

    def test_phi_x_and_multi_fasta_artifact_identities_are_unchanged(self) -> None:
        data = Path(__file__).parent / "data"
        single = SequenceCompiler().compile_file(data / "J02482.1.fasta")
        collection = SequenceCollectionCompiler().compile_file(data / "J02482.1.fasta")
        self.assertEqual(
            single.digest,
            "04d59484935bd73cec5499f44f59be9db30fe4d5d381ca28e567cb0e50ae9579",
        )
        self.assertEqual(
            collection.digest,
            "117cb7f85fa9ab6fa27d0618d2666764015136394b179a524c4850ee47eacd91",
        )

        multi = SequenceCollectionCompiler().compile_fasta_bytes(
            b">one\nAC\nGT\n>two desc\nNN\nry\n>three\nA\n"
        )
        self.assertEqual(
            multi.digest,
            "b334e0104a9a46bffb5ff7a3572a5e05146e27035781637909b095a0b08b5755",
        )

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            single_path = root / "single.json"
            collection_path = root / "collection.json"
            single.save(single_path)
            collection.save(collection_path)
            self.assertEqual(
                single_path.read_bytes(),
                (json.dumps(single.to_dict(), indent=2, sort_keys=True) + "\n").encode("utf-8"),
            )
            self.assertEqual(
                collection_path.read_bytes(),
                (
                    json.dumps(
                        collection.to_dict(),
                        indent=2,
                        sort_keys=True,
                        ensure_ascii=False,
                        allow_nan=False,
                    )
                    + "\n"
                ).encode("utf-8"),
            )


if __name__ == "__main__":
    unittest.main()

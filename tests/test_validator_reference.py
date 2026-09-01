from __future__ import annotations

import copy
import gzip
import hashlib
import inspect
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import tracemalloc
import unittest
from unittest import mock

from brainc.source_scale import compile_reference_fasta
import brainc.validator_reference as validator


ROOT = Path(__file__).resolve().parents[1]


def _reseal(artifact: dict) -> dict:
    artifact["catalog_sha256"] = validator.digest(artifact["sequence_catalog"])
    core = {
        key: value
        for key, value in artifact.items()
        if key != "artifact_sha256"
    }
    artifact["artifact_sha256"] = validator.digest(core)
    return artifact


class IndependentReferenceValidatorTests(unittest.TestCase):
    def test_identity_gzip_path_replay_and_reports_are_exact(self) -> None:
        raw = b">chr1 description\r\nacgt\n>chr2\nNNry\n"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            identity_path = root / "genome.fasta"
            identity_path.write_bytes(raw)
            gzip_path = root / "genome.fasta.gz"
            gzip_path.write_bytes(gzip.compress(raw, mtime=0))

            identity = compile_reference_fasta(identity_path, wrapper="identity")
            compressed = compile_reference_fasta(gzip_path, wrapper="gzip")
            self.assertEqual(
                identity["sequence_catalog"],
                compressed["sequence_catalog"],
            )
            self.assertEqual(
                validator.validate_reference_catalog(identity),
                identity,
            )

            identity_report = validator.validate_reference_fasta(
                identity_path,
                identity,
            )
            gzip_report = validator.validate_reference_fasta(gzip_path, compressed)
            self.assertTrue(identity_report["valid"])
            self.assertTrue(gzip_report["valid"])
            self.assertEqual(identity_report["result"]["records"], 2)
            self.assertEqual(identity_report["result"]["total_bases"], 8)
            self.assertEqual(
                validator.validate_reference_report(identity_report),
                identity_report,
            )

            artifact_path = root / "catalog.json"
            artifact_path.write_text(
                json.dumps(identity, ensure_ascii=False),
                encoding="utf-8",
            )
            first = validator.validate_reference_paths(identity_path, artifact_path)
            second = validator.validate_reference_paths(identity_path, artifact_path)
            self.assertEqual(first, second)
            self.assertEqual(first, identity_report)

    def test_import_isolation_loads_no_producer_or_shared_helpers(self) -> None:
        command = (
            "import sys; "
            f"sys.path.insert(0, {str(ROOT)!r}); "
            "import brainc.validator_reference; "
            "forbidden={'brainc.source_scale','brainc.source','brainc.compiler',"
            "'brainc.provider','brainc._canonical','brainc._io','brainc.v2',"
            "'brainc.v2.compiler'}; "
            "loaded=sorted(forbidden.intersection(sys.modules)); "
            "assert not loaded, loaded"
        )
        subprocess.run(
            [sys.executable, "-I", "-c", command],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )

    def test_direct_resealed_and_closed_schema_tampering_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "source.fasta"
            source.write_bytes(b">a\nACGT\n>b\nNNRY\n")
            artifact = compile_reference_fasta(source, wrapper="identity")

            direct = copy.deepcopy(artifact)
            direct["sequence_catalog"]["records"][0]["sequence_sha256"] = "0" * 64
            with self.assertRaisesRegex(
                validator.ReferenceValidationError,
                "catalog_sha256",
            ):
                validator.validate_reference_catalog(direct)

            resealed = _reseal(copy.deepcopy(direct))
            self.assertEqual(
                validator.validate_reference_catalog(resealed),
                resealed,
            )
            with self.assertRaisesRegex(
                validator.ReferenceValidationError,
                "replay does not match",
            ):
                validator.validate_reference_fasta(source, resealed)

            attacks = []
            version_bool = copy.deepcopy(artifact)
            version_bool["version"] = True
            attacks.append(_reseal(version_bool))
            length_float = copy.deepcopy(artifact)
            length_float["inputs"]["source"]["byte_length"] = float(
                length_float["inputs"]["source"]["byte_length"]
            )
            attacks.append(_reseal(length_float))
            unknown = copy.deepcopy(artifact)
            unknown["training_data"] = "forbidden"
            attacks.append(_reseal(unknown))
            seqcol_float = copy.deepcopy(artifact)
            seqcol_float["sequence_catalog"]["refget_seqcol"]["level_2"][
                "lengths"
            ][0] = 4.0
            attacks.append(_reseal(seqcol_float))
            identity_divergence = copy.deepcopy(artifact)
            identity_divergence["inputs"]["logical"]["sha256"] = "0" * 64
            attacks.append(_reseal(identity_divergence))
            for attack in attacks:
                with self.subTest(attack=attack), self.assertRaises(
                    validator.ReferenceValidationError
                ):
                    validator.validate_reference_catalog(attack)

            report = validator.validate_reference_fasta(source, artifact)
            report["result"]["records"] += 1
            with self.assertRaisesRegex(
                validator.ReferenceValidationError,
                "report seal",
            ):
                validator.validate_reference_report(report)

    def test_wrapper_and_gzip_corruption_fail_closed(self) -> None:
        raw = b">chr\nACGTNNRY\n"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            identity_path = root / "source.fasta"
            identity_path.write_bytes(raw)
            gzip_path = root / "source.fasta.gz"
            gzip_path.write_bytes(gzip.compress(raw, mtime=0))
            identity = compile_reference_fasta(identity_path, wrapper="identity")
            compressed = compile_reference_fasta(gzip_path, wrapper="gzip")

            with self.assertRaises(validator.ReferenceValidationError):
                validator.validate_reference_fasta(identity_path, compressed)
            with self.assertRaises(validator.ReferenceValidationError):
                validator.validate_reference_fasta(gzip_path, identity)

            truncated_path = root / "truncated.fasta.gz"
            truncated_path.write_bytes(gzip_path.read_bytes()[:-4])
            truncated = copy.deepcopy(compressed)
            truncated["inputs"]["source"]["sha256"] = "0" * 64
            truncated["inputs"]["source"]["byte_length"] = truncated_path.stat().st_size
            truncated = _reseal(truncated)
            with self.assertRaises(validator.ReferenceValidationError):
                validator.validate_reference_fasta(truncated_path, truncated)

    def test_strict_concatenated_gzip_members_and_member_ceiling(self) -> None:
        members = (
            gzip.compress(b"", mtime=0),
            gzip.compress(b">chr\nACGTNNRY\n", mtime=0),
            gzip.compress(b"", mtime=0),
        )
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "concatenated.fasta.gz"
            source.write_bytes(b"".join(members))
            artifact = compile_reference_fasta(source, wrapper="gzip")
            report = validator.validate_reference_fasta(
                source,
                artifact,
                maximum_gzip_members=3,
            )
            self.assertEqual(report["result"]["total_bases"], 8)
            with self.assertRaisesRegex(
                validator.ReferenceValidationError,
                "member ceiling 2",
            ):
                validator.validate_reference_fasta(
                    source,
                    artifact,
                    maximum_gzip_members=2,
                )
            with self.assertRaisesRegex(
                validator.ReferenceValidationError,
                "maximum_gzip_members",
            ):
                validator.validate_reference_fasta(
                    source,
                    artifact,
                    maximum_gzip_members=validator.MAXIMUM_GZIP_MEMBERS + 1,
                )

    def test_strict_gzip_rejects_padding_and_trailing_junk(self) -> None:
        raw = b">chr\nACGT\n"
        member = gzip.compress(raw, mtime=0)
        empty = gzip.compress(b"", mtime=0)
        attacks = (
            (member, b"J" + member),
            (member, member + b"\0"),
            (member, member + b"J"),
            (member, member + b"\x1f"),
            (member, member + b"\0" * 32),
            (member, member + b"junk"),
            (member + empty, member + b"\0" + empty),
            (member + empty, (member + empty)[:-1]),
        )
        for valid_raw, attacked_raw in attacks:
            with self.subTest(attacked_raw=attacked_raw[-32:]), tempfile.TemporaryDirectory() as temporary:
                source = Path(temporary) / "invalid.fasta.gz"
                source.write_bytes(valid_raw)
                artifact = compile_reference_fasta(source, wrapper="gzip")
                source.write_bytes(attacked_raw)
                artifact["inputs"]["source"] = {
                    "sha256": hashlib.sha256(attacked_raw).hexdigest(),
                    "byte_length": len(attacked_raw),
                }
                artifact = _reseal(artifact)
                with self.assertRaisesRegex(
                    validator.ReferenceValidationError,
                    "outside complete members|invalid gzip FASTA member|incomplete member",
                ):
                    validator.validate_reference_fasta(source, artifact)

    def test_crlf_header_at_chunk_boundary_replays_exactly(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "boundary.fasta"
            prefix = b">a\n"
            first_bases = validator.CHUNK_BYTES - len(prefix) - 1
            source.write_bytes(prefix + b"A" * first_bases + b"\r\n>b\nC\n")
            artifact = compile_reference_fasta(source, wrapper="identity")
            report = validator.validate_reference_fasta(source, artifact)
            self.assertEqual(report["result"]["records"], 2)
            self.assertEqual(report["result"]["total_bases"], first_bases + 1)

    @unittest.skipUnless(os.name == "posix", "link attacks require POSIX")
    def test_source_and_catalog_paths_reject_symbolic_and_hard_links(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.fasta"
            source.write_bytes(b">a\nACGT\n")
            artifact = compile_reference_fasta(source, wrapper="identity")
            artifact_path = root / "catalog.json"
            artifact_path.write_text(json.dumps(artifact), encoding="utf-8")

            source_symlink = root / "source-symlink.fasta"
            source_symlink.symlink_to(source)
            artifact_symlink = root / "catalog-symlink.json"
            artifact_symlink.symlink_to(artifact_path)
            with self.assertRaisesRegex(
                validator.ReferenceValidationError,
                "single-link",
            ):
                validator.validate_reference_paths(source_symlink, artifact_path)
            with self.assertRaisesRegex(
                validator.ReferenceValidationError,
                "single-link",
            ):
                validator.validate_reference_paths(source, artifact_symlink)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.fasta"
            source.write_bytes(b">a\nACGT\n")
            artifact = compile_reference_fasta(source, wrapper="identity")
            artifact_path = root / "catalog.json"
            artifact_path.write_text(json.dumps(artifact), encoding="utf-8")
            source_hardlink = root / "source-hardlink.fasta"
            os.link(source, source_hardlink)
            with self.assertRaisesRegex(
                validator.ReferenceValidationError,
                "single-link",
            ):
                validator.validate_reference_paths(source_hardlink, artifact_path)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.fasta"
            source.write_bytes(b">a\nACGT\n")
            artifact = compile_reference_fasta(source, wrapper="identity")
            artifact_path = root / "catalog.json"
            artifact_path.write_text(json.dumps(artifact), encoding="utf-8")
            artifact_hardlink = root / "catalog-hardlink.json"
            os.link(artifact_path, artifact_hardlink)
            with self.assertRaisesRegex(
                validator.ReferenceValidationError,
                "single-link",
            ):
                validator.validate_reference_paths(source, artifact_hardlink)

    def test_duplicate_json_and_resource_ceilings_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.fasta"
            source.write_bytes(b">a\nACGT\n")
            artifact = compile_reference_fasta(source, wrapper="identity")
            artifact_path = root / "catalog.json"
            artifact_path.write_text(json.dumps(artifact), encoding="utf-8")
            duplicate_path = root / "duplicate.json"
            duplicate_path.write_text('{"format":"a","format":"b"}', encoding="utf-8")

            with self.assertRaisesRegex(
                validator.ReferenceValidationError,
                "duplicate JSON key",
            ):
                validator.validate_reference_paths(source, duplicate_path)
            with self.assertRaisesRegex(
                validator.ReferenceValidationError,
                "physical source exceeds|reference FASTA exceeds",
            ):
                validator.validate_reference_fasta(
                    source,
                    artifact,
                    maximum_source_bytes=source.stat().st_size - 1,
                )
            with self.assertRaisesRegex(
                validator.ReferenceValidationError,
                "logical FASTA exceeds",
            ):
                validator.validate_reference_fasta(
                    source,
                    artifact,
                    maximum_logical_bytes=source.stat().st_size - 1,
                )
            with self.assertRaisesRegex(
                validator.ReferenceValidationError,
                "reference catalog exceeds",
            ):
                validator.validate_reference_paths(
                    source,
                    artifact_path,
                    maximum_artifact_bytes=artifact_path.stat().st_size - 1,
                )

    def test_default_and_hard_override_source_budgets_are_exact(self) -> None:
        hard_maximum = 2**53 - 1
        self.assertEqual(validator.DEFAULT_MAX_INPUT_BYTES, 64 * 1024**3)
        self.assertEqual(validator.DEFAULT_MAX_LOGICAL_BYTES, 64 * 1024**3)
        for api in (
            validator.validate_reference_fasta,
            validator.validate_reference_paths,
        ):
            parameters = inspect.signature(api).parameters
            self.assertEqual(
                parameters["maximum_source_bytes"].default,
                validator.DEFAULT_MAX_INPUT_BYTES,
            )
            self.assertEqual(
                parameters["maximum_logical_bytes"].default,
                validator.DEFAULT_MAX_LOGICAL_BYTES,
            )
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "source.fasta"
            source.write_bytes(b">a\nACGT\n")
            artifact = compile_reference_fasta(source, wrapper="identity")

            with mock.patch.object(
                validator,
                "_observed_artifact",
                wraps=validator._observed_artifact,
            ) as replay:
                validator.validate_reference_fasta(source, artifact)
            self.assertEqual(
                replay.call_args.kwargs["maximum_source_bytes"],
                validator.DEFAULT_MAX_INPUT_BYTES,
            )
            self.assertEqual(
                replay.call_args.kwargs["maximum_logical_bytes"],
                validator.DEFAULT_MAX_LOGICAL_BYTES,
            )

            with mock.patch.object(
                validator,
                "_observed_artifact",
                wraps=validator._observed_artifact,
            ) as replay:
                validator.validate_reference_fasta(
                    source,
                    artifact,
                    maximum_source_bytes=hard_maximum,
                    maximum_logical_bytes=hard_maximum,
                )
            self.assertEqual(
                replay.call_args.kwargs["maximum_source_bytes"],
                hard_maximum,
            )
            self.assertEqual(
                replay.call_args.kwargs["maximum_logical_bytes"],
                hard_maximum,
            )

            for name in ("maximum_source_bytes", "maximum_logical_bytes"):
                with self.subTest(name=name), self.assertRaisesRegex(
                    validator.ReferenceValidationError,
                    name,
                ):
                    validator.validate_reference_fasta(
                        source,
                        artifact,
                        **{name: hard_maximum + 1},
                    )

    def test_invalid_and_embedded_nul_paths_are_normalized(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.fasta"
            source.write_bytes(b">a\nACGT\n")
            artifact = compile_reference_fasta(source, wrapper="identity")
            artifact_path = root / "catalog.json"
            artifact_path.write_text(json.dumps(artifact), encoding="utf-8")
            for invalid in (None, 1, b"bytes", "embedded\0nul"):
                with self.subTest(invalid=invalid):
                    with self.assertRaises(validator.ReferenceValidationError) as caught:
                        validator.validate_reference_fasta(invalid, artifact)
                    self.assertTrue(str(caught.exception).startswith("REFVAL001:"))
                    with self.assertRaises(validator.ReferenceValidationError) as caught:
                        validator.validate_reference_paths(source, invalid)
                    self.assertTrue(str(caught.exception).startswith("REFVAL001:"))

    def test_catalog_json_is_bounded_before_decode(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.fasta"
            source.write_bytes(b">a\nACGT\n")

            oversized = root / "oversized.json"
            with oversized.open("wb") as stream:
                stream.truncate(validator.MAX_ARTIFACT_BYTES + 1)
            with self.assertRaisesRegex(
                validator.ReferenceValidationError,
                "reference catalog exceeds byte ceiling",
            ):
                validator.validate_reference_paths(source, oversized)

            deep = root / "deep.json"
            deep.write_bytes(
                b'{"x":' + b"[" * 65 + b"0" + b"]" * 65 + b"}"
            )
            with self.assertRaisesRegex(
                validator.ReferenceValidationError,
                "JSON depth ceiling",
            ):
                validator.validate_reference_paths(source, deep)

            with mock.patch.object(validator, "MAX_JSON_MEMBERS", 3):
                with self.assertRaisesRegex(
                    validator.ReferenceValidationError,
                    "JSON member ceiling",
                ):
                    validator._load_json(b'{"x":[1,2,3]}', "catalog")
            with mock.patch.object(validator, "MAX_STRING_BYTES", 3):
                with self.assertRaisesRegex(
                    validator.ReferenceValidationError,
                    "JSON string byte ceiling",
                ):
                    validator._load_json(b'{"x":"1234"}', "catalog")

    def test_memory_errors_are_normalized_at_public_boundaries(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "source.fasta"
            source.write_bytes(b">a\nACGT\n")
            artifact = compile_reference_fasta(source, wrapper="identity")
            with mock.patch.object(
                validator,
                "_observed_artifact",
                side_effect=MemoryError("injected"),
            ), self.assertRaisesRegex(
                validator.ReferenceValidationError,
                "memory ceiling exceeded",
            ):
                validator.validate_reference_fasta(source, artifact)

    def test_streaming_replay_exceeds_legacy_physical_and_logical_limits(self) -> None:
        block = b"A" * validator.CHUNK_BYTES
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            identity_path = root / "large.fasta"
            identity_bases = 17 * 1024 * 1024 + 1
            with identity_path.open("wb") as stream:
                stream.write(b">chromosome-identity\n")
                remaining = identity_bases
                while remaining:
                    piece = block[: min(len(block), remaining)]
                    stream.write(piece)
                    remaining -= len(piece)
                stream.write(b"\n")
            identity = compile_reference_fasta(identity_path, wrapper="identity")

            tracemalloc.start()
            identity_report = validator.validate_reference_fasta(
                identity_path,
                identity,
            )
            _, identity_peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            self.assertEqual(identity_report["result"]["total_bases"], identity_bases)
            self.assertGreater(identity_report["result"]["source_bytes"], 16 * 1024 * 1024)
            self.assertLess(identity_peak, 8 * 1024 * 1024)

            gzip_path = root / "large.fasta.gz"
            gzip_bases = 64 * 1024 * 1024 + 1
            with gzip_path.open("wb") as physical:
                with gzip.GzipFile(fileobj=physical, mode="wb", mtime=0) as stream:
                    stream.write(b">chromosome-gzip\n")
                    remaining = gzip_bases
                    while remaining:
                        piece = block[: min(len(block), remaining)]
                        stream.write(piece)
                        remaining -= len(piece)
                    stream.write(b"\n")
            compressed = compile_reference_fasta(gzip_path, wrapper="gzip")

            tracemalloc.start()
            gzip_report = validator.validate_reference_fasta(gzip_path, compressed)
            _, gzip_peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            self.assertEqual(gzip_report["result"]["total_bases"], gzip_bases)
            self.assertGreater(gzip_report["result"]["logical_bytes"], 64 * 1024 * 1024)
            self.assertLess(gzip_peak, 8 * 1024 * 1024)


if __name__ == "__main__":
    unittest.main()

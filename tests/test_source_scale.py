from __future__ import annotations

from copy import deepcopy
import gzip
import hashlib
import json
import os
from pathlib import Path
import tempfile
import tracemalloc
import unittest
from unittest import mock

from brainc._canonical import digest
import brainc.source_scale as source_scale
from brainc.source_scale import (
    DEFAULT_LIMITS,
    DEFAULT_MAX_INPUT_BYTES,
    DEFAULT_MAX_LOGICAL_BYTES,
    PROFILE,
    ScaleLimits,
    SourceScaleError,
    compile_reference_fasta,
    replay_reference_fasta,
    validate_reference_catalog,
)


def _reseal(value: dict) -> dict:
    value["catalog_sha256"] = digest(value["sequence_catalog"])
    core = {key: item for key, item in value.items() if key != "artifact_sha256"}
    value["artifact_sha256"] = digest(core)
    return value


class ReferenceFastaTests(unittest.TestCase):
    def test_default_byte_budgets_are_finite_and_explicitly_overridable(self) -> None:
        self.assertEqual(DEFAULT_MAX_INPUT_BYTES, 64 * 1024**3)
        self.assertEqual(DEFAULT_MAX_LOGICAL_BYTES, 64 * 1024**3)
        self.assertEqual(DEFAULT_LIMITS.maximum_input_bytes, DEFAULT_MAX_INPUT_BYTES)
        self.assertEqual(
            DEFAULT_LIMITS.maximum_logical_bytes,
            DEFAULT_MAX_LOGICAL_BYTES,
        )

        hard_maximum = 2**53 - 1
        exceptional = ScaleLimits(
            maximum_input_bytes=hard_maximum,
            maximum_logical_bytes=hard_maximum,
        )
        self.assertEqual(exceptional.maximum_input_bytes, hard_maximum)
        self.assertEqual(exceptional.maximum_logical_bytes, hard_maximum)
        with self.assertRaises(SourceScaleError):
            ScaleLimits(maximum_logical_bytes=hard_maximum + 1)

    def test_line_ending_grammar_and_all_byte_splits_are_deterministic(self) -> None:
        variants = []
        for separator in (b"\n", b"\r\n", b"\r"):
            variants.append(
                separator.join(
                    (b">a description", b"ac", b"GT", b"", b">b", b"nR", b"")
                )
            )
        variants.append(b">a description\r\nac\nGT\r\r>b\r\n\r\nnR\n")

        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "grammar.fasta"
            catalogs = []
            for index, raw in enumerate(variants):
                with self.subTest(variant=index):
                    source.write_bytes(raw)
                    expected = compile_reference_fasta(source, wrapper="identity")
                    catalogs.append(expected["sequence_catalog"])
                    with mock.patch.object(source_scale, "_CHUNK_BYTES", 1):
                        observed = compile_reference_fasta(source, wrapper="identity")
                    self.assertEqual(observed, expected)

            self.assertTrue(all(catalog == catalogs[0] for catalog in catalogs[1:]))
            self.assertEqual(
                [(record["record_id"], record["bases"]) for record in catalogs[0]["records"]],
                [("a", 4), ("b", 2)],
            )

            interior = b">a\nAC>GT\n"
            source.write_bytes(interior)
            for chunk_bytes in (1, 2, 3, 64 * 1024):
                with self.subTest(interior_chunk_bytes=chunk_bytes):
                    with mock.patch.object(
                        source_scale,
                        "_CHUNK_BYTES",
                        chunk_bytes,
                    ), self.assertRaises(SourceScaleError):
                        compile_reference_fasta(source, wrapper="identity")

    def test_crlf_header_marker_split_at_stream_chunk_boundary(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "boundary.fasta"
            prefix = b">a\n"
            first_bases = 64 * 1024 - len(prefix) - 1
            source.write_bytes(
                prefix + b"A" * first_bases + b"\r\n>b\nC\n"
            )
            records = compile_reference_fasta(
                source,
                wrapper="identity",
            )["sequence_catalog"]["records"]
            self.assertEqual(
                [(record["record_id"], record["bases"]) for record in records],
                [("a", first_bases), ("b", 1)],
            )

    def test_large_unwrapped_record_is_streamed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "chromosome.fasta"
            bases = 17 * 1024 * 1024 + 1
            block = b"A" * (64 * 1024)
            with source.open("wb") as stream:
                stream.write(b">chromosome-1\n")
                remaining = bases
                while remaining:
                    piece = block[: min(len(block), remaining)]
                    stream.write(piece)
                    remaining -= len(piece)
                stream.write(b"\n")

            tracemalloc.start()
            artifact = compile_reference_fasta(source, wrapper="identity")
            _, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()

            record = artifact["sequence_catalog"]["records"][0]
            self.assertEqual(record["bases"], bases)
            self.assertEqual(
                record["sequence_sha256"],
                hashlib.sha256(b"A" * bases).hexdigest(),
            )
            self.assertLess(peak, 8 * 1024 * 1024)
            wire = json.dumps(artifact, separators=(",", ":"))
            self.assertLess(len(wire), 4096)
            self.assertNotIn("sequence\"", wire)

    def test_large_gzip_record_is_streamed_past_legacy_logical_limit(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "chromosome.fasta.gz"
            bases = 64 * 1024 * 1024 + 1
            block = b"A" * (64 * 1024)
            with source.open("wb") as physical:
                with gzip.GzipFile(fileobj=physical, mode="wb", mtime=0) as stream:
                    stream.write(b">chromosome-1\n")
                    remaining = bases
                    while remaining:
                        piece = block[: min(len(block), remaining)]
                        stream.write(piece)
                        remaining -= len(piece)
                    stream.write(b"\n")

            tracemalloc.start()
            artifact = compile_reference_fasta(source, wrapper="gzip")
            _, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()

            self.assertEqual(artifact["sequence_catalog"]["total_bases"], bases)
            self.assertLess(source.stat().st_size, 1024 * 1024)
            self.assertLess(peak, 1024 * 1024)
            self.assertLess(len(json.dumps(artifact, separators=(",", ":"))), 4096)

    def test_identity_and_gzip_are_exact_and_deterministic(self) -> None:
        raw = b">chr1 description\r\nacgt\n>chr2\nNNry\n"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "genome.fasta"
            source.write_bytes(raw)
            compressed = root / "genome.fasta.gz"
            compressed.write_bytes(gzip.compress(raw, mtime=0))

            identity = compile_reference_fasta(source, wrapper="identity")
            self.assertEqual(
                identity,
                compile_reference_fasta(source, wrapper="identity"),
            )
            zipped = compile_reference_fasta(compressed, wrapper="gzip")
            self.assertEqual(identity["sequence_catalog"], zipped["sequence_catalog"])
            self.assertEqual(identity["inputs"]["profile"], PROFILE)
            self.assertEqual(identity["inputs"]["source"]["sha256"], hashlib.sha256(raw).hexdigest())
            self.assertEqual(identity["inputs"]["source"]["byte_length"], len(raw))
            self.assertEqual(zipped["inputs"]["logical"]["sha256"], hashlib.sha256(raw).hexdigest())
            self.assertEqual(zipped["inputs"]["logical"]["byte_length"], len(raw))
            self.assertEqual(
                identity["sequence_catalog"]["records"],
                [
                    {
                        "record_id": "chr1",
                        "bases": 4,
                        "sequence_sha256": hashlib.sha256(b"ACGT").hexdigest(),
                        "refget_id": "SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2",
                    },
                    {
                        "record_id": "chr2",
                        "bases": 4,
                        "sequence_sha256": hashlib.sha256(b"NNRY").hexdigest(),
                        "refget_id": "SQ.pdpAuSJhwmprg7AdEwy26bZrWim9uP4M",
                    },
                ],
            )
            self.assertEqual(replay_reference_fasta(source, identity), identity)
            self.assertEqual(replay_reference_fasta(compressed, zipped), zipped)

            concatenated = root / "concatenated.fasta.gz"
            split = len(raw) // 2
            concatenated.write_bytes(
                gzip.compress(raw[:split], mtime=0)
                + gzip.compress(raw[split:], mtime=0)
            )
            concatenated_artifact = compile_reference_fasta(
                concatenated,
                wrapper="gzip",
            )
            self.assertEqual(
                concatenated_artifact["sequence_catalog"],
                identity["sequence_catalog"],
            )
            with mock.patch.object(source_scale, "_CHUNK_BYTES", 1):
                self.assertEqual(
                    compile_reference_fasta(concatenated, wrapper="gzip"),
                    concatenated_artifact,
                )

            many_members = root / "many-members.fasta.gz"
            many_members.write_bytes(
                b"".join(gzip.compress(b"", mtime=0) for _ in range(1025))
                + gzip.compress(raw, mtime=0)
            )
            self.assertEqual(
                compile_reference_fasta(many_members, wrapper="gzip")[
                    "sequence_catalog"
                ],
                identity["sequence_catalog"],
            )

    def test_attacks_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.fasta"
            source.write_bytes(b">aa\nACGT\n>b\nNNNN\n")
            compressed = root / "source.fasta.gz"
            compressed.write_bytes(gzip.compress(source.read_bytes(), mtime=0))

            for path, wrapper in (
                (source, "auto"),
                (source, "gzip"),
                (compressed, "identity"),
            ):
                with self.subTest(path=path.name, wrapper=wrapper):
                    with self.assertRaises(SourceScaleError):
                        compile_reference_fasta(path, wrapper=wrapper)

            duplicate = root / "duplicate.fasta"
            duplicate.write_bytes(b">x\nA\n>x\nC\n")
            invalid = root / "invalid.fasta"
            invalid.write_bytes(b">x\nACGT!\n")
            truncated = root / "truncated.fasta.gz"
            truncated.write_bytes(compressed.read_bytes()[:-4])
            compressed_bytes = compressed.read_bytes()
            corrupt_crc = root / "corrupt-crc.fasta.gz"
            corrupt_crc_bytes = bytearray(compressed_bytes)
            corrupt_crc_bytes[-8] ^= 1
            corrupt_crc.write_bytes(corrupt_crc_bytes)
            corrupt_size = root / "corrupt-size.fasta.gz"
            corrupt_size_bytes = bytearray(compressed_bytes)
            corrupt_size_bytes[-4] ^= 1
            corrupt_size.write_bytes(corrupt_size_bytes)
            reserved_flag = root / "reserved-flag.fasta.gz"
            reserved_flag_bytes = bytearray(compressed_bytes)
            reserved_flag_bytes[3] |= 0x20
            reserved_flag.write_bytes(reserved_flag_bytes)
            leading = root / "leading-junk.fasta.gz"
            leading.write_bytes(b"J" + compressed_bytes)
            between = root / "between-junk.fasta.gz"
            between.write_bytes(
                compressed_bytes + b"J" + gzip.compress(b"", mtime=0)
            )
            too_many_members = root / "too-many-members.fasta.gz"
            too_many_members.write_bytes(
                b"".join(
                    gzip.compress(b"", mtime=0)
                    for _ in range(3)
                )
            )
            trailing = []
            for index, suffix in enumerate((b"\0", b"J", b"\x1f", b"\0" * 16)):
                path = root / f"trailing-{index}.fasta.gz"
                path.write_bytes(compressed.read_bytes() + suffix)
                trailing.append(path)
            for path, wrapper in (
                (duplicate, "identity"),
                (invalid, "identity"),
                (truncated, "gzip"),
                (corrupt_crc, "gzip"),
                (corrupt_size, "gzip"),
                (reserved_flag, "gzip"),
                (leading, "gzip"),
                (between, "gzip"),
                *((path, "gzip") for path in trailing),
            ):
                with self.subTest(path=path.name):
                    with self.assertRaises(SourceScaleError):
                        compile_reference_fasta(path, wrapper=wrapper)

            with self.assertRaisesRegex(SourceScaleError, "member ceiling 2"):
                compile_reference_fasta(
                    too_many_members,
                    wrapper="gzip",
                    limits=ScaleLimits(maximum_gzip_members=2),
                )

            with self.assertRaises(SourceScaleError):
                compile_reference_fasta(
                    source,
                    wrapper="identity",
                    limits=ScaleLimits(maximum_logical_bytes=8),
                )
            resource_limits = (
                ScaleLimits(maximum_input_bytes=8),
                ScaleLimits(maximum_records=1),
                ScaleLimits(maximum_record_bases=3),
                ScaleLimits(maximum_header_bytes=1),
            )
            for limits in resource_limits:
                with self.subTest(limits=limits):
                    with self.assertRaises(SourceScaleError):
                        compile_reference_fasta(source, wrapper="identity", limits=limits)

            if os.name == "posix":
                if getattr(os, "O_NONBLOCK", 0):
                    flags = []
                    real_open = os.open

                    def checked_open(path: object, value: int) -> int:
                        flags.append(value)
                        return real_open(path, value)  # type: ignore[arg-type]

                    with mock.patch.object(
                        source_scale.os,
                        "open",
                        side_effect=checked_open,
                    ):
                        compile_reference_fasta(source, wrapper="identity")
                    self.assertTrue(flags[0] & os.O_NONBLOCK)

                alias = root / "alias.fasta"
                alias.symlink_to(source)
                with self.assertRaises(SourceScaleError):
                    compile_reference_fasta(alias, wrapper="identity")

                hardlink = root / "hardlink.fasta"
                hardlink.hardlink_to(source)
                with self.assertRaises(SourceScaleError):
                    compile_reference_fasta(source, wrapper="identity")
                hardlink.unlink()

            with self.assertRaises(SourceScaleError):
                compile_reference_fasta(source, wrapper=[])

            for invalid_path in (None, 1, b"source.fasta", "bad\0path"):
                with self.subTest(invalid_path=invalid_path):
                    with self.assertRaisesRegex(SourceScaleError, "SCALE001"):
                        compile_reference_fasta(
                            invalid_path,  # type: ignore[arg-type]
                            wrapper="identity",
                        )

            artifact = compile_reference_fasta(source, wrapper="identity")
            typed_attacks = []
            version_float = deepcopy(artifact)
            version_float["version"] = 1.0
            typed_attacks.append(_reseal(version_float))
            length_float = deepcopy(artifact)
            length_float["inputs"]["source"]["byte_length"] = float(
                length_float["inputs"]["source"]["byte_length"]
            )
            typed_attacks.append(_reseal(length_float))
            total_float = deepcopy(artifact)
            total_float["sequence_catalog"]["total_bases"] = float(
                total_float["sequence_catalog"]["total_bases"]
            )
            typed_attacks.append(_reseal(total_float))
            seqcol_float = deepcopy(artifact)
            seqcol_float["sequence_catalog"]["refget_seqcol"]["level_2"]["lengths"][0] = 4.0
            typed_attacks.append(_reseal(seqcol_float))
            for attack in typed_attacks:
                with self.subTest(attack=attack):
                    with self.assertRaises(SourceScaleError):
                        validate_reference_catalog(attack)

            detached_logical = deepcopy(artifact)
            detached_logical["inputs"]["logical"]["sha256"] = "0" * 64
            detached_logical = _reseal(detached_logical)
            with self.assertRaisesRegex(SourceScaleError, "must be identical"):
                validate_reference_catalog(detached_logical)

            forged = deepcopy(artifact)
            forged["sequence_catalog"]["records"][0]["sequence_sha256"] = "0" * 64
            forged = _reseal(forged)
            self.assertEqual(validate_reference_catalog(forged), forged)
            with self.assertRaisesRegex(SourceScaleError, "replay does not match"):
                replay_reference_fasta(source, forged)

            changed = root / "changed.fasta"
            changed.write_bytes(b">aa\r\nACGT\r\n>b\r\nNNNN\r\n")
            with self.assertRaisesRegex(SourceScaleError, "replay does not match"):
                replay_reference_fasta(changed, artifact)

    def test_memory_errors_are_normalized_at_public_boundaries(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "source.fasta"
            source.write_bytes(b">a\nACGT\n")
            artifact = compile_reference_fasta(source, wrapper="identity")

            with mock.patch.object(
                source_scale,
                "_open_regular",
                side_effect=MemoryError("probe"),
            ), self.assertRaisesRegex(SourceScaleError, "memory ceiling"):
                compile_reference_fasta(source, wrapper="identity")

            with mock.patch.object(
                source_scale,
                "_keys",
                side_effect=MemoryError("probe"),
            ), self.assertRaisesRegex(SourceScaleError, "memory ceiling"):
                validate_reference_catalog(artifact)

            with mock.patch.object(
                source_scale,
                "validate_reference_catalog",
                return_value=artifact,
            ), mock.patch.object(
                source_scale,
                "compile_reference_fasta",
                side_effect=MemoryError("probe"),
            ), self.assertRaisesRegex(SourceScaleError, "memory ceiling"):
                replay_reference_fasta(source, artifact)


if __name__ == "__main__":
    unittest.main()

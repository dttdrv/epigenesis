from __future__ import annotations

import ast
import copy
import gzip
import hashlib
import json
from pathlib import Path
import unittest
from unittest import mock

from brainc._canonical import digest
from brainc.bio import GFF3Compiler
import brainc.sequence as sequence_frontend
import brainc.sequence_collection as collection_frontend
from brainc.sequence_collection import (
    SequenceCollectionCompiler,
    SequenceCollectionError,
)
from brainc.validator_bio import (
    BioValidationError,
    validate_bio_chain,
    validate_bio_report,
)
import brainc.validator_bio as independent_validator


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data"


def rehash_bio(payload: dict) -> dict:
    payload["bio_ir_sha256"] = digest(payload["bio_ir"])
    payload["artifact_sha256"] = digest(
        {key: value for key, value in payload.items() if key != "artifact_sha256"}
    )
    return payload


class IndependentBioValidatorTests(unittest.TestCase):
    def compile(self, fasta: bytes, gff3: bytes) -> tuple[dict, dict]:
        collection = SequenceCollectionCompiler().compile_fasta_bytes(fasta).to_dict()
        bio = GFF3Compiler().compile_bytes(gff3, collection).to_dict()
        return collection, bio

    def minimal_chain(self) -> tuple[bytes, bytes, dict, dict]:
        fasta = b">circle\nACGTACGTAC\n"
        gff3 = (
            b"##gff-version 3\n"
            b"##sequence-region circle 1 10\n"
            b"circle\t.\tregion\t1\t10\t.\t.\t.\tID=circle;Is_circular=true\n"
            b"circle\t.\tgene\t8\t13\t.\t+\t.\tID=gene-1;old-name=opaque\n"
            b"circle\t.\tCDS\t8\t13\t.\t-\t2\tParent=gene-1\n"
        )
        collection, bio = self.compile(fasta, gff3)
        return fasta, gff3, collection, bio

    def test_actual_phix_raw_sources_and_metadata_replay(self) -> None:
        fasta = (DATA / "J02482.1.fasta").read_bytes()
        gff3 = (DATA / "J02482.1.gff3").read_bytes()
        metadata = (DATA / "J02482.1.source.json").read_bytes()
        gff3_metadata = json.loads(
            (DATA / "J02482.1.gff3.source.json").read_text(encoding="utf-8")
        )
        self.assertEqual(gff3_metadata["record"], "J02482.1")
        self.assertEqual(gff3_metadata["sha256"], hashlib.sha256(gff3).hexdigest())
        collection, bio = self.compile(fasta, gff3)

        first = validate_bio_chain(
            bio,
            collection,
            fasta_source=fasta,
            gff3_source=gff3,
            source_metadata=metadata,
        )
        second = validate_bio_chain(
            bio,
            collection,
            fasta_source=fasta,
            gff3_source=gff3,
            source_metadata=metadata,
        )

        self.assertEqual(first, second)
        self.assertIs(validate_bio_report(first), first)
        self.assertTrue(first["valid"])
        self.assertEqual(first["summary"]["records"], 1)
        self.assertEqual(first["summary"]["sequence_bases"], 5386)
        self.assertEqual(first["summary"]["features"], 23)
        self.assertEqual(first["summary"]["segments"], 23)
        self.assertEqual(first["summary"]["relationships"], 4)
        self.assertEqual(first["summary"]["circular_records"], ["J02482.1"])
        self.assertTrue(first["summary"]["source_metadata_verified"])
        self.assertEqual(first["replay"]["bio_artifact_sha256"], bio["artifact_sha256"])
        self.assertEqual(
            first["replay"]["sequence_collection_sha256"],
            collection["artifact_sha256"],
        )

    def test_validator_has_no_compiler_import_or_dynamic_import_escape(self) -> None:
        path = ROOT / "brainc" / "validator_bio.py"
        source = path.read_text(encoding="utf-8")
        tree = ast.parse(source)
        allowed = {
            "__future__",
            "base64",
            "decimal",
            "gzip",
            "hashlib",
            "io",
            "json",
            "math",
            "re",
            "typing",
        }
        imported: set[str] = set()
        forbidden_calls: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.update(alias.name.split(".", 1)[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                self.assertIsNotNone(node.module)
                imported.add(node.module.split(".", 1)[0])
            elif isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                if node.func.id in {"eval", "exec", "compile", "__import__"}:
                    forbidden_calls.append(node.func.id)
        self.assertLessEqual(imported, allowed)
        self.assertNotIn("brainc", imported)
        self.assertEqual(forbidden_calls, [])

    def test_source_limits_match_the_published_collection_ingress(self) -> None:
        self.assertEqual(independent_validator.MAX_SOURCE_INPUT_BYTES, 16 * 1024 * 1024)
        self.assertEqual(independent_validator.MAX_LOGICAL_DNA_BYTES, 64 * 1024 * 1024)
        self.assertEqual(independent_validator.MAX_SOURCE_ARTIFACT_BYTES, 16 * 1024 * 1024)
        self.assertEqual(independent_validator.MAX_TEXT_BYTES, 1024 * 1024)
        with mock.patch.object(independent_validator, "MAX_TEXT_BYTES", 3):
            with self.assertRaisesRegex(BioValidationError, "sequence exceeds 3 bytes"):
                independent_validator._replay_collection(b">x\nAAAA\n")
            with self.assertRaisesRegex(BioValidationError, "sequence exceeds 3 bytes"):
                independent_validator._replay_collection(b"AAAA", record_id="x")

    def test_sequence_identifier_and_line_boundaries_match_compiler(self) -> None:
        accepted_id = "x" * 256
        fasta = f">{accepted_id}\nACGT\n".encode("ascii")
        gff3 = (
            f"##gff-version 3\n{accepted_id}\t.\tgene\t1\t4\t.\t+\t.\tID=g\n"
        ).encode("ascii")
        collection, bio = self.compile(fasta, gff3)
        self.assertTrue(
            validate_bio_chain(
                bio,
                collection,
                fasta_source=fasta,
                gff3_source=gff3,
            )["valid"]
        )

        too_long = f">{'x' * 257}\nACGT\n".encode("ascii")
        with self.assertRaisesRegex(SequenceCollectionError, "256 UTF-8 bytes"):
            SequenceCollectionCompiler().compile_fasta_bytes(too_long)
        with self.assertRaisesRegex(BioValidationError, "256 UTF-8 bytes"):
            independent_validator._replay_collection(too_long)

        for source in (b">x\x00y\nACGT\n", b">x\x7fy\nACGT\n"):
            with self.subTest(source=source):
                with self.assertRaisesRegex(
                    SequenceCollectionError, "control character"
                ):
                    SequenceCollectionCompiler().compile_fasta_bytes(source)
                with self.assertRaisesRegex(BioValidationError, "control character"):
                    independent_validator._replay_collection(source)

        collection = SequenceCollectionCompiler().compile_fasta_bytes(b">x\nACGT\n")
        for escaped in (b"%00", b"%7F", b"%C2%85"):
            source = (
                b"##gff-version 3\n"
                + escaped
                + b"\t.\tgene\t1\t4\t.\t+\t.\tID=g\n"
            )
            with self.subTest(seqid=escaped):
                with self.assertRaisesRegex(Exception, "control character"):
                    GFF3Compiler().compile_bytes(source, collection)

        at_boundary = b">x\n\nACGT"
        above_boundary = b">x\n\n\nACGT"
        patches = (
            mock.patch.object(sequence_frontend, "MAX_SOURCE_LINE_BREAKS", 2),
            mock.patch.object(collection_frontend, "MAX_SOURCE_LINE_BREAKS", 2),
            mock.patch.object(independent_validator, "MAX_SEQUENCE_LINE_BREAKS", 2),
        )
        with patches[0], patches[1], patches[2]:
            compiled = SequenceCollectionCompiler().compile_fasta_bytes(at_boundary)
            replayed, _ = independent_validator._replay_collection(at_boundary)
            self.assertEqual(compiled.to_dict(), replayed)
            with self.assertRaisesRegex(SequenceCollectionError, "line breaks"):
                SequenceCollectionCompiler().compile_fasta_bytes(above_boundary)
            with self.assertRaisesRegex(BioValidationError, "line breaks"):
                independent_validator._replay_collection(above_boundary)

    def test_deep_object_metadata_fails_closed_before_canonicalization(self) -> None:
        metadata: dict = {}
        for _ in range(independent_validator.MAX_JSON_DEPTH + 2):
            metadata = {"nested": metadata}
        with self.assertRaisesRegex(BioValidationError, "JSON depth"):
            independent_validator._source_metadata(metadata, b"", [])

    def test_coherently_rehashed_bio_tampering_is_rejected(self) -> None:
        fasta, gff3, collection, original = self.minimal_chain()
        attacks = []

        interval = copy.deepcopy(original)
        crossing = next(
            feature
            for feature in interval["bio_ir"]["features"]
            if feature["declared_id"] == "gene-1"
        )
        crossing["segments"][0]["intervals"][1]["end"] -= 1
        attacks.append(rehash_bio(interval))

        relationship = copy.deepcopy(original)
        relationship["bio_ir"]["relationships"] = []
        attacks.append(rehash_bio(relationship))

        source_line = copy.deepcopy(original)
        feature = next(
            feature
            for feature in source_line["bio_ir"]["features"]
            if feature["declared_id"] == "gene-1"
        )
        feature["segments"][0]["line"] = 99
        attacks.append(rehash_bio(source_line))

        identity = copy.deepcopy(original)
        anonymous = next(
            feature
            for feature in identity["bio_ir"]["features"]
            if feature["declared_id"] is None
        )
        old_id = anonymous["entity_id"]
        anonymous["entity_id"] = "anon:" + "0" * 64
        identity["bio_ir"]["features"].sort(key=lambda item: item["entity_id"])
        for edge in identity["bio_ir"]["relationships"]:
            if edge["source"] == old_id:
                edge["source"] = anonymous["entity_id"]
        identity["bio_ir"]["relationships"].sort(
            key=lambda item: (item["source"], item["kind"], item["target"])
        )
        attacks.append(rehash_bio(identity))

        unknown = copy.deepcopy(original)
        unknown["bio_ir"]["features"][0]["invented"] = True
        attacks.append(rehash_bio(unknown))

        for attack in attacks:
            with self.subTest(artifact=attack["artifact_sha256"]):
                with self.assertRaisesRegex(BioValidationError, "source replay"):
                    validate_bio_chain(
                        attack,
                        collection,
                        fasta_source=fasta,
                        gff3_source=gff3,
                    )

    def test_external_source_and_collection_substitution_is_rejected(self) -> None:
        fasta, gff3, collection, bio = self.minimal_chain()
        other_fasta = fasta.replace(b"ACGTACGTAC", b"ACGTACGTAA")
        other_collection = SequenceCollectionCompiler().compile_fasta_bytes(other_fasta).to_dict()

        with self.assertRaisesRegex(BioValidationError, "source replay"):
            validate_bio_chain(
                bio,
                other_collection,
                fasta_source=fasta,
                gff3_source=gff3,
            )
        with self.assertRaises(BioValidationError):
            validate_bio_chain(
                bio,
                collection,
                fasta_source=other_fasta,
                gff3_source=gff3,
            )
        with self.assertRaises(BioValidationError):
            validate_bio_chain(
                bio,
                collection,
                fasta_source=fasta,
                gff3_source=gff3.replace(b"old-name=opaque", b"old-name=changed"),
            )

    def test_strict_escape_rules_are_independently_enforced(self) -> None:
        fasta, gff3, collection, bio = self.minimal_chain()
        invalid = [
            gff3.replace(b"gene-1", b"gene-%31", 1),
            gff3.replace(b"old-name=opaque", b"old-name=bad%20space"),
            gff3.replace(b"old-name=opaque", b"old-name=bad&value"),
            gff3.replace(b"old-name=opaque", b"old-name=bad\x01value"),
        ]
        for source in invalid:
            with self.subTest(source=source):
                with self.assertRaises(BioValidationError):
                    validate_bio_chain(
                        bio,
                        collection,
                        fasta_source=fasta,
                        gff3_source=source,
                    )

    def test_ampersand_is_reserved_only_in_attribute_column(self) -> None:
        fasta, gff3, _, _ = self.minimal_chain()
        source = gff3.replace(
            b"circle\t.\tgene\t8\t13",
            b"circle\ttool&x\topaque&type\t8\t13",
        )
        collection, bio = self.compile(fasta, source)
        report = validate_bio_chain(
            bio,
            collection,
            fasta_source=fasta,
            gff3_source=source,
        )
        self.assertTrue(report["valid"])

    def test_hash_fragments_and_ordered_application_values_replay(self) -> None:
        fasta = b">x\nAAAAAAAAAAAAAAAAAAAA\n"
        gff3 = (
            b"##gff-version 3\n"
            b"##species https://example.test/taxon#fragment\n"
            b"# whole-line comment\n"
            b"x\t.\tgene\t1\t10\t.\t+\t.\tID=a\n"
            b"x\t.\tgene\t1\t10\t.\t+\t.\tID=b\n"
            b"x\t.\texon\t2\t5\t.\t+\t.\tID=child;Parent=b,a;Note=C#;old-name=former#1;app-values=z,a,z\n"
        )
        collection, bio = self.compile(fasta, gff3)
        report = validate_bio_chain(
            bio,
            collection,
            fasta_source=fasta,
            gff3_source=gff3,
        )
        self.assertTrue(report["valid"])
        child = next(
            feature for feature in bio["bio_ir"]["features"]
            if feature["declared_id"] == "child"
        )
        attributes = {
            item["tag"]: item["values"] for item in child["segments"][0]["attributes"]
        }
        self.assertEqual(attributes["app-values"], ["z", "a", "z"])
        self.assertEqual(attributes["Parent"], ["a", "b"])
        self.assertEqual(attributes["Note"], ["C#"])

    def test_coordinate_lexemes_are_rejected_by_independent_replay(self) -> None:
        fasta = b">x\nAAAAAAAAAAAAAAAAAAAA\n"
        valid_gff3 = b"##gff-version 3\nx\t.\tgene\t1\t2\t.\t+\t.\tID=g\n"
        collection, bio = self.compile(fasta, valid_gff3)
        for token in ("+1", " 1", "1 ", "١", "１", "0", "-1", "1.0"):
            invalid = f"##gff-version 3\nx\t.\tgene\t{token}\t2\t.\t+\t.\tID=g\n".encode("utf-8")
            with self.subTest(feature_coordinate=token):
                with self.assertRaisesRegex(
                    BioValidationError, "ASCII positive-decimal|integer in"
                ):
                    validate_bio_chain(
                        bio,
                        collection,
                        fasta_source=fasta,
                        gff3_source=invalid,
                    )

        region_gff3 = (
            b"##gff-version 3\n"
            b"##sequence-region x 1 20\n"
            b"x\t.\tgene\t1\t2\t.\t+\t.\tID=g\n"
        )
        region_collection, region_bio = self.compile(fasta, region_gff3)
        for token in ("+1", "١", "１", "0", "-1", "1.0"):
            invalid = region_gff3.replace(b"sequence-region x 1 20", f"sequence-region x {token} 20".encode("utf-8"))
            with self.subTest(sequence_region_coordinate=token):
                with self.assertRaisesRegex(
                    BioValidationError, "ASCII positive-decimal|integer in"
                ):
                    validate_bio_chain(
                        region_bio,
                        region_collection,
                        fasta_source=fasta,
                        gff3_source=invalid,
                    )

    def test_tiny_circular_feature_cannot_authorize_origin_crossing(self) -> None:
        fasta, gff3, collection, bio = self.minimal_chain()
        forged_source = gff3.replace(
            b"region\t1\t10\t.\t.\t.\tID=circle;Is_circular=true",
            b"region\t2\t3\t.\t.\t.\tID=circle;Is_circular=true",
        )
        with self.assertRaisesRegex(BioValidationError, "non-circular sequence origin"):
            validate_bio_chain(
                bio,
                collection,
                fasta_source=fasta,
                gff3_source=forged_source,
            )

    def test_gzip_collection_replay_and_report_seal(self) -> None:
        fasta, gff3, _, _ = self.minimal_chain()
        wrapped = gzip.compress(fasta, mtime=0)
        collection, bio = self.compile(wrapped, gff3)
        report = validate_bio_chain(
            bio,
            collection,
            fasta_source=wrapped,
            gff3_source=gff3,
        )
        self.assertFalse(report["summary"]["source_metadata_verified"])

        altered = copy.deepcopy(report)
        altered["summary"]["features"] += 1
        with self.assertRaisesRegex(BioValidationError, "report_sha256"):
            validate_bio_report(altered)

    def test_raw_iupac_collection_requires_exact_record_id(self) -> None:
        raw = b"Acgtryswkmbdhvn"
        record_id = "raw-record"
        gff3 = (
            b"##gff-version 3\n"
            b"raw-record\t.\tregion\t1\t15\t.\t.\t.\tID=raw-record\n"
        )
        collection = SequenceCollectionCompiler().compile_raw(raw, record_id).to_dict()
        bio = GFF3Compiler().compile_bytes(gff3, collection).to_dict()
        report = validate_bio_chain(
            bio,
            collection,
            fasta_source=raw,
            gff3_source=gff3,
            record_id=record_id,
        )
        self.assertTrue(report["valid"])
        self.assertEqual(report["summary"]["sequence_bases"], 15)
        with self.assertRaisesRegex(BioValidationError, "sequence collection"):
            validate_bio_chain(
                bio,
                collection,
                fasta_source=raw,
                gff3_source=gff3,
                record_id="wrong-record",
            )
        with self.assertRaisesRegex(BioValidationError, "requires.*record_id"):
            validate_bio_chain(
                bio,
                collection,
                fasta_source=raw,
                gff3_source=gff3,
            )

    def test_source_metadata_is_checked_as_evidence_not_compiler_input(self) -> None:
        fasta = (DATA / "J02482.1.fasta").read_bytes()
        gff3 = (DATA / "J02482.1.gff3").read_bytes()
        metadata = (DATA / "J02482.1.source.json").read_bytes()
        collection, bio = self.compile(fasta, gff3)
        forged = metadata.replace(b'"bases": 5386', b'"bases": 5385')
        with self.assertRaisesRegex(BioValidationError, "base count"):
            validate_bio_chain(
                bio,
                collection,
                fasta_source=fasta,
                gff3_source=gff3,
                source_metadata=forged,
            )


if __name__ == "__main__":
    unittest.main()

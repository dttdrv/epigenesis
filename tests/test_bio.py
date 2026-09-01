from __future__ import annotations

import copy
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from brainc._canonical import digest
from brainc.bio import (
    GFF3Compiler,
    GFF3Error,
    MAX_ATTRIBUTES_PER_ROW,
    MAX_OUTPUT_BYTES,
    load_gff3_artifact,
    validate_gff3_artifact,
)
from brainc.sequence_collection import SequenceCollectionCompiler


def rehash(payload: dict) -> dict:
    payload["bio_ir_sha256"] = digest(payload["bio_ir"])
    payload["artifact_sha256"] = digest(
        {key: value for key, value in payload.items() if key != "artifact_sha256"}
    )
    return payload


class GFF3CompilerTests(unittest.TestCase):
    def collection(self, *records: tuple[str, str]):
        fasta = "".join(f">{record_id}\n{sequence}\n" for record_id, sequence in records)
        return SequenceCollectionCompiler().compile_fasta_bytes(fasta.encode("ascii"))

    def test_compile_and_validation_preflight_the_exact_saved_wire(self) -> None:
        collection = self.collection(("x", "ACGT"))
        source = b"##gff-version 3\nx\t.\tgene\t1\t4\t.\t+\t.\tID=g\n"
        artifact = GFF3Compiler().compile_bytes(source, collection).to_dict()
        with mock.patch(
            "brainc.bio._pretty_wire_size", return_value=MAX_OUTPUT_BYTES + 1
        ):
            with self.assertRaisesRegex(GFF3Error, "serialized artifact exceeds"):
                GFF3Compiler().compile_bytes(source, collection)
            with self.assertRaisesRegex(GFF3Error, "serialized artifact exceeds"):
                validate_gff3_artifact(
                    artifact, collection, gff3_source=source
                )

    def test_structural_profile_preserves_opaque_feature_types(self) -> None:
        collection = self.collection(("x", "ACGT"))
        source = (
            b"##gff-version 3\n"
            b"x\t.\tnot_an_SO_term\t1\t4\t.\t+\t.\tID=opaque\n"
        )
        artifact = GFF3Compiler().compile_bytes(source, collection).to_dict()
        self.assertEqual(artifact["bio_ir"]["features"][0]["type"], "not_an_SO_term")

    def test_forward_links_discontinuous_ids_and_opaque_extensions(self) -> None:
        collection = self.collection(("chr,1", "ACGT" * 8), ("unused", "NNNN"))
        source = (
            "##gff-version 3.1.26\n"
            "##sequence-region chr%2C1 1 32\n"
            "##pipeline opaque α value\n"
            "chr%2C1\tsource one\texon\t2\t5\t.\t+\t.\tID=exon%3B1;Parent=tx1;custom=a%2Cb,c\n"
            "chr%2C1\tsource one\tCDS\t2\t5\t-2.50e+1\t+\t0\tID=cds1;Parent=tx1\n"
            "chr%2C1\tsource one\texon\t9\t12\t.\t+\t.\tID=exon%3B1;Parent=tx1\n"
            "chr%2C1\tsource one\tmRNA\t1\t20\t.\t+\t.\tID=tx1;Parent=gene1;Derives_from=template1\n"
            "chr%2C1\t.\tgene\t1\t20\t.\t+\t.\tID=gene1;Name=opaque name\n"
            "chr%2C1\t.\tregion\t1\t32\t.\t.\t.\tID=template1\n"
        ).encode("utf-8")

        first = GFF3Compiler().compile_bytes(source, collection).to_dict()
        second = GFF3Compiler().compile_bytes(source, collection).to_dict()
        self.assertEqual(first, second)
        ir = first["bio_ir"]
        self.assertEqual(ir["coordinate_system"], {
            "source": "1-based-closed", "normalized": "0-based-half-open"
        })
        self.assertEqual(ir["directives"], [
            {"line": 3, "name": "pipeline", "value": "opaque α value"}
        ])
        self.assertEqual([record["seqid"] for record in ir["records"]], ["chr,1", "unused"])
        exon = next(feature for feature in ir["features"] if feature["entity_id"] == "id:exon;1")
        self.assertEqual(len(exon["segments"]), 2)
        self.assertEqual(exon["source"], "source one")
        custom = next(
            attribute
            for attribute in exon["segments"][0]["attributes"]
            if attribute["tag"] == "custom"
        )
        self.assertEqual(custom["values"], ["a,b", "c"])
        cds = next(feature for feature in ir["features"] if feature["entity_id"] == "id:cds1")
        self.assertEqual(cds["segments"][0]["score"], "-2.50e+1")
        self.assertEqual(cds["segments"][0]["phase"], 0)
        self.assertIn(
            {"source": "id:tx1", "kind": "derives-from", "target": "id:template1"},
            ir["relationships"],
        )

    def test_circular_origin_crossing_and_all_strand_states(self) -> None:
        collection = self.collection(("circle", "ACGTACGTAC"))
        source = (
            "##gff-version 3\n"
            "##sequence-region circle 1 10\n"
            "circle\t.\tregion\t1\t10\t.\t.\t.\tID=circle;Is_circular=true\n"
            "circle\t.\tgene\t8\t13\t.\t+\t.\tID=wrap\n"
            "circle\t.\tCDS\t8\t13\t.\t-\t2\tParent=wrap\n"
            "circle\t.\tmarker\t2\t2\t.\t.\t.\tlowercase=opaque\n"
            "circle\t.\tmarker\t2\t2\t.\t?\t.\tlowercase=opaque\n"
        ).encode()
        payload = GFF3Compiler().compile_bytes(source, collection).to_dict()
        wrap = next(feature for feature in payload["bio_ir"]["features"] if feature["entity_id"] == "id:wrap")
        self.assertEqual(wrap["segments"][0]["raw_start"], 8)
        self.assertEqual(wrap["segments"][0]["raw_end"], 13)
        self.assertEqual(wrap["segments"][0]["intervals"], [
            {"start": 7, "end": 10}, {"start": 0, "end": 3}
        ])
        anonymous = [feature for feature in payload["bio_ir"]["features"] if feature["declared_id"] is None]
        self.assertEqual({feature["strand"] for feature in anonymous}, {"-", ".", "?"})
        self.assertEqual(len({feature["entity_id"] for feature in anonymous}), 3)

    def test_crossing_requires_explicit_circularity_and_one_lap(self) -> None:
        collection = self.collection(("x", "ACGTACGTAC"))
        without_marker = (
            b"##gff-version 3\n"
            b"x\t.\tgene\t8\t13\t.\t+\t.\tID=wrap\n"
        )
        with self.assertRaisesRegex(GFF3Error, "without Is_circular=true"):
            GFF3Compiler().compile_bytes(without_marker, collection)
        over_one_lap = (
            b"##gff-version 3\n"
            b"x\t.\tregion\t1\t10\t.\t.\t.\tID=x;Is_circular=true\n"
            b"x\t.\tgene\t2\t12\t.\t+\t.\tID=too-long\n"
        )
        with self.assertRaisesRegex(GFF3Error, "more than one circular"):
            GFF3Compiler().compile_bytes(over_one_lap, collection)
        non_landmark_marker = (
            b"##gff-version 3\n"
            b"x\t.\tregion\t2\t3\t.\t.\t.\tID=not-landmark;Is_circular=true\n"
            b"x\t.\tgene\t8\t13\t.\t+\t.\tID=wrap\n"
        )
        with self.assertRaisesRegex(GFF3Error, "without Is_circular=true"):
            GFF3Compiler().compile_bytes(non_landmark_marker, collection)
        backwards = b"##gff-version 3\nx\t.\tgene\t8\t3\t.\t+\t.\tID=no\n"
        with self.assertRaisesRegex(GFF3Error, "end precedes start"):
            GFF3Compiler().compile_bytes(backwards, collection)

    def test_sequence_region_and_record_boundaries(self) -> None:
        collection = self.collection(("a", "ACGTACGT"), ("b", "NNNN"))
        cases = [
            b"##gff-version 3\nmissing\t.\tgene\t1\t2\t.\t+\t.\tID=g\n",
            b"##gff-version 3\n##sequence-region a 1 9\na\t.\tgene\t1\t2\t.\t+\t.\tID=g\n",
            b"##gff-version 3\n##sequence-region a 2 7\na\t.\tgene\t1\t2\t.\t+\t.\tID=g\n",
            b"##gff-version 3\n##sequence-region a 1 8\n##sequence-region a 1 8\n",
        ]
        for source in cases:
            with self.subTest(source=source):
                with self.assertRaises(GFF3Error):
                    GFF3Compiler().compile_bytes(source, collection)

    def test_relationship_forward_references_dangling_and_cycles(self) -> None:
        collection = self.collection(("x", "A" * 20))
        forward = (
            b"##gff-version 3\n"
            b"x\t.\texon\t2\t3\t.\t+\t.\tParent=parent\n"
            b"x\t.\tmRNA\t1\t10\t.\t+\t.\tID=parent\n"
        )
        self.assertEqual(
            len(GFF3Compiler().compile_bytes(forward, collection).to_dict()["bio_ir"]["relationships"]),
            1,
        )
        dangling = b"##gff-version 3\nx\t.\texon\t2\t3\t.\t+\t.\tParent=absent\n"
        with self.assertRaisesRegex(GFF3Error, "dangling Parent"):
            GFF3Compiler().compile_bytes(dangling, collection)
        cycle = (
            b"##gff-version 3\n"
            b"x\t.\tgene\t1\t5\t.\t+\t.\tID=a;Parent=b\n"
            b"x\t.\tmRNA\t1\t5\t.\t+\t.\tID=b;Derives_from=a\n"
        )
        with self.assertRaisesRegex(GFF3Error, "relationship cycle"):
            GFF3Compiler().compile_bytes(cycle, collection)

    def test_bounded_profile_syntax_phase_score_and_extensions(self) -> None:
        collection = self.collection(("x", "A" * 20))
        invalid_sources = [
            b"# comment\n##gff-version 3\n",
            b"##gff-version 3\nx . gene 1 2 . + . ID=g\n",
            b"##gff-version 3\nx\t.\tgene\t1\t2\tNaN\t+\t.\tID=g\n",
            b"##gff-version 3\nx\t.\tCDS\t1\t2\t.\t+\t.\tID=g\n",
            b"##gff-version 3\nx\t.\tgene\t1\t2\t.\t+\t0\tID=g\n",
            b"##gff-version 3\nx\t.\tgene\t1\t2\t.\tx\t.\tID=g\n",
            b"##gff-version 3\nx\t.\tgene\t1\t2\t.\t+\t.\tID=%ZZ\n",
            b"##gff-version 3\n%78\t.\tgene\t1\t2\t.\t+\t.\tID=g\n",
            b"##gff-version 3\nx\tsource%20one\tgene\t1\t2\t.\t+\t.\tID=g\n",
            b"##gff-version 3\nx\t.\tgene\t1\t2\t.\t+\t.\t1bad=x\n",
            b"##gff-version 3\nx\t.\tgene\t1\t2\t.\t+\t.\tID=a&b\n",
            b"##gff-version 3\nx\tsource\x01bad\tgene\t1\t2\t.\t+\t.\tID=g\n",
            b"##gff-version 3\n##FASTA\n>x\nAAAA\n",
            b"##gff-version 3\nx\t.\tgene\t1\t2\t.\t+\t.\tID=g;ID=h\n",
            b"##gff-version 3\nx\t.\tgene\t1\t2\t.\t+\t.\tFuture=x\n",
            b"##gff-version 3\n##Future opaque\n",
            b"##gff-version 3\nx\t.\tgene\t1\t2\t.\t+\t.\tIs_circular=maybe\n",
        ]
        for source in invalid_sources:
            with self.subTest(source=source):
                with self.assertRaises(GFF3Error):
                    GFF3Compiler().compile_bytes(source, collection)
        with self.assertRaisesRegex(GFF3Error, "strict UTF-8"):
            GFF3Compiler().compile_bytes(b"##gff-version 3\n\xff", collection)
        excessive_attributes = ";".join(
            f"tag{index}=v" for index in range(MAX_ATTRIBUTES_PER_ROW + 1)
        )
        with self.assertRaisesRegex(GFF3Error, "exceeds .* attributes"):
            GFF3Compiler().compile_bytes(
                f"##gff-version 3\nx\t.\tgene\t1\t2\t.\t+\t.\t{excessive_attributes}\n".encode(),
                collection,
            )

    def test_hash_data_and_application_value_order_are_preserved(self) -> None:
        collection = self.collection(("x", "A" * 20))
        source = (
            b"##gff-version 3\n"
            b"##species https://example.test/taxon#fragment\n"
            b"# a whole-line comment containing # is ignored\n"
            b"x\t.\tgene\t1\t10\t.\t+\t.\tID=a\n"
            b"x\t.\tgene\t1\t10\t.\t+\t.\tID=b\n"
            b"x\t.\texon\t2\t5\t.\t+\t.\tID=child;Parent=b,a;Note=C#;old-name=former#1;app-values=z,a,z\n"
        )
        payload = GFF3Compiler().compile_bytes(source, collection).to_dict()
        self.assertEqual(
            payload["bio_ir"]["directives"],
            [{"line": 2, "name": "species", "value": "https://example.test/taxon#fragment"}],
        )
        child = next(
            feature for feature in payload["bio_ir"]["features"]
            if feature["declared_id"] == "child"
        )
        attributes = {
            item["tag"]: item["values"] for item in child["segments"][0]["attributes"]
        }
        self.assertEqual(attributes["Note"], ["C#"])
        self.assertEqual(attributes["old-name"], ["former#1"])
        self.assertEqual(attributes["app-values"], ["z", "a", "z"])
        self.assertEqual(attributes["Parent"], ["a", "b"])
        self.assertIs(
            validate_gff3_artifact(payload, collection, gff3_source=source),
            payload,
        )

    def test_coordinates_require_ascii_positive_decimal_lexemes(self) -> None:
        collection = self.collection(("x", "A" * 20))
        feature_template = "##gff-version 3\nx\t.\tgene\t{}\t2\t.\t+\t.\tID=g\n"
        for token in ("+1", " 1", "1 ", "١", "１", "0", "-1", "1.0"):
            with self.subTest(feature_coordinate=token):
                with self.assertRaisesRegex(GFF3Error, "ASCII positive-decimal|integer in"):
                    GFF3Compiler().compile_bytes(feature_template.format(token).encode("utf-8"), collection)
        for token in ("+1", "١", "１", "0", "-1", "1.0"):
            source = f"##gff-version 3\n##sequence-region x {token} 20\n".encode("utf-8")
            with self.subTest(sequence_region_coordinate=token):
                with self.assertRaisesRegex(GFF3Error, "ASCII positive-decimal|integer in"):
                    GFF3Compiler().compile_bytes(source, collection)

    def test_rehashed_semantic_tampering_is_rejected(self) -> None:
        collection = self.collection(("x", "ACGTACGT"))
        source = (
            b"##gff-version 3\n"
            b"x\t.\tgene\t1\t8\t.\t+\t.\tID=parent\n"
            b"x\t.\tCDS\t2\t7\t.\t+\t0\tParent=parent\n"
        )
        original = GFF3Compiler().compile_bytes(source, collection).to_dict()
        mutations = []

        forged_interval = copy.deepcopy(original)
        forged_interval["bio_ir"]["features"][0]["segments"][0]["intervals"][0]["end"] -= 1
        mutations.append(rehash(forged_interval))

        forged_relationship = copy.deepcopy(original)
        forged_relationship["bio_ir"]["relationships"] = []
        mutations.append(rehash(forged_relationship))

        forged_anonymous_id = copy.deepcopy(original)
        anonymous = next(
            feature for feature in forged_anonymous_id["bio_ir"]["features"]
            if feature["declared_id"] is None
        )
        anonymous["entity_id"] = "anon:" + "0" * 64
        forged_anonymous_id["bio_ir"]["features"].sort(key=lambda feature: feature["entity_id"])
        forged_anonymous_id["bio_ir"]["relationships"][0]["source"] = anonymous["entity_id"]
        mutations.append(rehash(forged_anonymous_id))

        forged_binding = copy.deepcopy(original)
        forged_binding["inputs"]["sequence_collection"]["artifact_sha256"] = "0" * 64
        forged_binding["artifact_sha256"] = digest({
            key: value for key, value in forged_binding.items() if key != "artifact_sha256"
        })
        mutations.append(forged_binding)

        forged_unknown = copy.deepcopy(original)
        forged_unknown["bio_ir"]["features"][0]["nonsense"] = True
        mutations.append(rehash(forged_unknown))

        for forged in mutations:
            with self.subTest(forged=forged):
                with self.assertRaises(GFF3Error):
                    validate_gff3_artifact(forged, collection)

    def test_save_load_and_exact_collection_binding(self) -> None:
        collection = self.collection(("x", "ACGT"))
        other = self.collection(("x", "ACGA"))
        artifact = GFF3Compiler().compile_bytes(
            b"##gff-version 3\nx\t.\tgene\t1\t4\t.\t+\t.\tID=g\n",
            collection,
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bio.json"
            artifact.save(path)
            self.assertEqual(load_gff3_artifact(path, collection).to_dict(), artifact.to_dict())
            self.assertEqual(
                load_gff3_artifact(
                    path,
                    collection,
                    gff3_source=b"##gff-version 3\nx\t.\tgene\t1\t4\t.\t+\t.\tID=g\n",
                ).to_dict(),
                artifact.to_dict(),
            )
            with self.assertRaisesRegex(GFF3Error, "does not match supplied GFF3"):
                load_gff3_artifact(path, collection, gff3_source=b"##gff-version 3\n")
            with self.assertRaisesRegex(GFF3Error, "binding"):
                load_gff3_artifact(path, other)

            target = Path(directory) / "target.json"
            target.write_text("preserve", encoding="utf-8")
            linked = Path(directory) / "linked.json"
            linked.symlink_to(target)
            with self.assertRaisesRegex(GFF3Error, "regular non-linked"):
                artifact.save(linked)
            self.assertEqual(target.read_text(encoding="utf-8"), "preserve")

            duplicate = path.read_text().replace(
                '"format": "brainc.bio.gff3-ir",',
                '"format": "brainc.bio.gff3-ir", "format": "brainc.bio.gff3-ir",',
                1,
            )
            path.write_text(duplicate)
            with self.assertRaisesRegex(GFF3Error, "duplicate JSON key"):
                load_gff3_artifact(path, collection)


if __name__ == "__main__":
    unittest.main()

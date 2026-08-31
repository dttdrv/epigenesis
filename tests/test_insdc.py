from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import brainc.insdc as insdc_module
from brainc._canonical import digest
from brainc.insdc import (
    AUTHORITY,
    AUTHORITY_MANIFEST_SHA256,
    FEATURE_ID_PREFIX,
    FORMAT,
    PROFILE,
    GenBankCompiler,
    GenBankError,
    load_genbank_artifact,
    validate_genbank_artifact,
)


def _origin(sequence: str) -> list[str]:
    lines: list[str] = []
    for offset in range(0, len(sequence), 60):
        part = sequence[offset : offset + 60].lower()
        groups = " ".join(part[index : index + 10] for index in range(0, len(part), 10))
        lines.append(f"{offset + 1:>9} {groups}")
    return lines


def _record(
    accession: str,
    version: int,
    sequence: str,
    *,
    topology: str = "linear",
    features: list[str] | None = None,
    newline: str = "\n",
) -> bytes:
    record_id = f"{accession}.{version}"
    feature_lines = [
        f"     {'source':<15} {1}..{len(sequence)}",
        '                     /organism="Arbitrary organism"',
        '                     /mol_type="genomic DNA"',
    ]
    feature_lines.extend(features or [])
    lines = [
        f"LOCUS       GENERATED_{accession} {len(sequence)} bp DNA {topology} SYN 31-AUG-2026",
        "DEFINITION  Generated structural compiler input.",
        f"ACCESSION   {accession}",
        f"VERSION     {record_id}",
        "DBLINK      BioProject: PRJNA1",
        "KEYWORDS    .",
        "SOURCE      synthetic construct",
        "  ORGANISM  Arbitrary organism",
        "            artificial sequences.",
        "REFERENCE   1  (bases 1 to %d)" % len(sequence),
        "  AUTHORS   Example,A.",
        "  TITLE     Direct Submission",
        "  JOURNAL   Submitted (31-AUG-2026)",
        "FEATURES             Location/Qualifiers",
        *feature_lines,
        "ORIGIN",
        *_origin(sequence),
        "//",
    ]
    return (newline.join(lines) + newline).encode("ascii")


def _records() -> bytes:
    first_sequence = ("ACGTRYSWKMBDHVN" * 8)[:120]
    second_sequence = ("TTGCAAN" * 13)[:84]
    first = _record(
        "AA000001",
        2,
        first_sequence,
        topology="circular",
        newline="\r\n",
        features=[
            "     misc_feature    complement(join(<2..10,40..>50))",
            '                     /note="alpha ""quoted',
            '                     value"" continues"',
            '                     /db_xref="arbitrary:1"',
            '                     /db_xref="arbitrary:2"',
            "                     /pseudo",
            "     variation       120^1",
            '                     /note="circular boundary"',
            "     misc_feature    BB000002.3:3..8",
            '                     /note="same-input remote"',
        ],
    )
    second = _record(
        "BB000002",
        3,
        second_sequence,
        features=[
            "     misc_feature    order(2,4)",
            '                     /note="ordered sites"',
        ],
    )
    return first + second


class GenBankFrontendTests(unittest.TestCase):
    def test_feature_and_qualifier_limits_are_aggregate_across_records(self) -> None:
        raw = _record("XA000001", 1, "ACGT") + _record(
            "XB000002", 1, "TGCA"
        )
        with mock.patch.object(insdc_module, "MAX_FEATURES", 1):
            with self.assertRaisesRegex(GenBankError, "aggregate features"):
                GenBankCompiler().compile_bytes(raw)
        with mock.patch.object(insdc_module, "MAX_QUALIFIERS", 2):
            with self.assertRaisesRegex(GenBankError, "aggregate qualifiers"):
                GenBankCompiler().compile_bytes(raw)

    def test_combined_json_member_budget_preflights_before_serialization(self) -> None:
        raw = _record(
            "XC000003",
            1,
            "ACGT" * 4,
            features=[
                "     misc_feature    1..4",
                *['                     /note="x"'] * 5,
                "     misc_feature    5..8",
                *['                     /note="y"'] * 5,
            ],
        )
        payload = GenBankCompiler().compile_bytes(raw).to_dict()
        exact_members = insdc_module._json_member_count(payload)

        with (
            mock.patch.object(insdc_module, "MAX_JSON_MEMBERS", exact_members - 1),
            mock.patch.object(
                insdc_module,
                "compact_json_bytes",
                side_effect=AssertionError("serialization must not run"),
            ),
        ):
            with self.assertRaisesRegex(GenBankError, "aggregate JSON members"):
                GenBankCompiler().compile_bytes(raw)

    def test_profile_is_bound_to_the_packaged_authority_manifest(self) -> None:
        manifest_path = (
            Path(insdc_module.__file__).with_name("standards")
            / "genbank-273-insdc-ft-11.4.authority.json"
        )
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        claimed = manifest.pop("artifact_sha256")
        self.assertEqual(claimed, AUTHORITY_MANIFEST_SHA256)
        self.assertEqual(digest(manifest), AUTHORITY_MANIFEST_SHA256)
        self.assertEqual(
            AUTHORITY["authority_manifest_sha256"],
            AUTHORITY_MANIFEST_SHA256,
        )

        with tempfile.TemporaryDirectory() as directory:
            tampered_path = Path(directory) / "authority.json"
            tampered = copy.deepcopy(manifest)
            tampered["sources"]["genbank_flatfile"]["declared_version"] = "273.1"
            tampered["artifact_sha256"] = claimed
            tampered_path.write_text(json.dumps(tampered), encoding="utf-8")
            with mock.patch.object(
                insdc_module,
                "_AUTHORITY_MANIFEST_PATH",
                tampered_path,
            ):
                with self.assertRaisesRegex(GenBankError, "content pin"):
                    GenBankCompiler().compile_bytes(_record("AZ000001", 1, "ACGT"))

    def test_compiles_official_ncbi_sample_record_without_identity_rules(self) -> None:
        fixture = Path(__file__).parent / "data" / "U49845.1.gb"
        raw = fixture.read_bytes()
        artifact = GenBankCompiler().compile_bytes(raw)
        payload = artifact.to_dict()
        record = payload["bio_ir"]["records"][0]
        sequence_ir = payload["sequence_collection"]["members"][0]["sequence"]

        self.assertEqual(hashlib.sha256(raw).hexdigest(), "1b41f0096dece0626236d49e0bede43b07ab4fd1aa576c1284cc4fb50bea2fd4")
        self.assertEqual(record["record_id"], "U49845.1")
        self.assertEqual(record["legacy_gi"], 1293613)
        self.assertEqual(sequence_ir["bases"], 5028)
        self.assertEqual(len(record["features"]), 6)
        self.assertTrue(
            all(
                feature["feature_id"].startswith(FEATURE_ID_PREFIX)
                and len(feature["feature_id"]) == len(FEATURE_ID_PREFIX) + 64
                for feature in record["features"]
            )
        )
        self.assertNotIn(":", FEATURE_ID_PREFIX)
        translations = [
            qualifier["value"]
            for feature in record["features"]
            for qualifier in feature["qualifiers"]
            if qualifier["name"] == "translation"
        ]
        self.assertEqual([len(value) for value in translations], [67, 823, 245])
        self.assertTrue(all(" " not in value for value in translations))
        self.assertEqual(sequence_ir["sha256"], "36203848f0560d3bc23561205438cbb3dc22b068ecfc92f9e3d77334b9ed9e9d")
        self.assertEqual(sequence_ir["refget_id"], "SQ.lDmukOs0TZpjopBrnhAwI775FEAhj379")
        validate_genbank_artifact(payload, genbank_source=raw)

    def test_compiles_two_records_and_preserves_closed_structural_semantics(self) -> None:
        raw = _records()
        artifact = GenBankCompiler().compile_bytes(raw)
        payload = artifact.to_dict()

        self.assertEqual(payload["format"], FORMAT)
        self.assertEqual(payload["profile"], PROFILE)
        self.assertEqual(payload["authority"], AUTHORITY)
        self.assertNotIn("inputs", payload)
        self.assertEqual(
            payload["sequence_collection"]["inputs"]["source"]["sha256"],
            hashlib.sha256(raw).hexdigest(),
        )
        self.assertEqual(
            payload["bio_ir"]["sequence_collection_artifact_sha256"],
            payload["sequence_collection"]["artifact_sha256"],
        )
        records = payload["bio_ir"]["records"]
        self.assertEqual([record["record_id"] for record in records], ["AA000001.2", "BB000002.3"])
        self.assertEqual(
            [member["record_id"] for member in payload["sequence_collection"]["members"]],
            ["AA000001.2", "BB000002.3"],
        )
        self.assertEqual(
            "".join(
                payload["sequence_collection"]["members"][0]["sequence"]["storage"]["chunks"]
            ),
            ("ACGTRYSWKMBDHVN" * 8)[:120],
        )

        joined = records[0]["features"][1]
        self.assertEqual(joined["location"]["ast"]["kind"], "complement")
        self.assertEqual(joined["location"]["ast"]["location"]["kind"], "join")
        self.assertEqual(joined["location"]["segments"][0]["orientation"], -1)
        self.assertEqual(joined["location"]["segments"][0]["end_fuzz"], ">")
        self.assertEqual(joined["qualifiers"][0]["value"], 'alpha "quoted value" continues')
        self.assertEqual(
            [qualifier["value"] for qualifier in joined["qualifiers"] if qualifier["name"] == "db_xref"],
            ["arbitrary:1", "arbitrary:2"],
        )
        self.assertIsNone(joined["qualifiers"][-1]["value"])

        between = records[0]["features"][2]["location"]["segments"][0]
        self.assertEqual((between["kind"], between["start"], between["end"]), ("between", 0, 0))
        remote = records[0]["features"][3]["location"]["segments"][0]
        self.assertEqual((remote["reference"], remote["start"], remote["end"]), ("BB000002.3", 2, 8))
        ordered = records[1]["features"][1]["location"]["ast"]
        self.assertEqual(ordered["kind"], "order")
        self.assertEqual(records[0]["relationships"], [])
        self.assertNotIn("sequence", records[0])
        self.assertNotIn("source_chunks", records[0])
        self.assertNotIn("chunks", records[0]["headers"][0])
        self.assertNotIn("text", records[0]["headers"][0]["source_map"][0])
        self.assertNotIn("chunks", joined["location"])
        self.assertNotIn("chunks", joined["qualifiers"][0])
        self.assertEqual(artifact.sequence_collection.source_bytes(), raw)
        validate_genbank_artifact(payload, genbank_source=raw)

    def test_compilation_is_deterministic_and_accessors_are_defensive(self) -> None:
        raw = _records()
        with mock.patch.object(
            insdc_module,
            "validate_genbank_artifact",
            side_effect=AssertionError("producer must not run consumer replay"),
        ):
            first = GenBankCompiler().compile_bytes(raw)
        second = GenBankCompiler().compile_bytes(raw)
        self.assertEqual(first.to_dict(), second.to_dict())
        self.assertEqual(first.digest, second.digest)
        copy_out = first.to_dict()
        copy_out["bio_ir"]["records"][0]["definition"] = "mutated"
        self.assertNotEqual(copy_out, first.to_dict())
        collection = first.sequence_collection
        self.assertEqual(collection.to_dict(), first.to_dict()["sequence_collection"])

    def test_save_load_and_full_source_replay(self) -> None:
        raw = _records()
        artifact = GenBankCompiler().compile_bytes(raw)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "compiled.json"
            artifact.save(path)
            loaded = load_genbank_artifact(path, genbank_source=raw)
            self.assertEqual(loaded.to_dict(), artifact.to_dict())
            with self.assertRaisesRegex(GenBankError, "does not match embedded"):
                load_genbank_artifact(path, genbank_source=raw.replace(b"Generated", b"Altered__", 1))

    def test_tamper_fails_at_digest_and_at_semantic_replay(self) -> None:
        payload = GenBankCompiler().compile_bytes(_records()).to_dict()
        direct = copy.deepcopy(payload)
        direct["bio_ir"]["records"][0]["definition"] = "tampered"
        with self.assertRaisesRegex(GenBankError, "artifact_sha256"):
            validate_genbank_artifact(direct)

        resealed = copy.deepcopy(payload)
        resealed["bio_ir"]["records"][0]["definition"] = "tampered"
        resealed["bio_ir_sha256"] = digest(resealed["bio_ir"])
        core = {key: value for key, value in resealed.items() if key != "artifact_sha256"}
        resealed["artifact_sha256"] = digest(core)
        with self.assertRaisesRegex(GenBankError, "frontend replay"):
            validate_genbank_artifact(resealed)

        for replacement in (True, 1.0):
            with self.subTest(type_confusion=type(replacement).__name__):
                confused = copy.deepcopy(payload)
                confused["bio_ir"]["records"][0]["features"][0]["ordinal"] = replacement
                confused["bio_ir_sha256"] = digest(confused["bio_ir"])
                confused["artifact_sha256"] = digest(
                    {
                        key: value
                        for key, value in confused.items()
                        if key != "artifact_sha256"
                    }
                )
                with self.assertRaisesRegex(GenBankError, "frontend replay"):
                    validate_genbank_artifact(confused)

        invalid_keys = copy.deepcopy(payload)
        invalid_keys[0] = "not JSON"
        invalid_keys[()] = "not JSON"
        with self.assertRaisesRegex(GenBankError, "non-string object key"):
            validate_genbank_artifact(invalid_keys)

    def test_rejects_malformed_records_and_profile_exits(self) -> None:
        valid = _record("CC000003", 1, "ACGT" * 20)
        cases = {
            "bare carriage return": valid.replace(b"\n", b"\r", 1),
            "tab": valid.replace(b"DEFINITION  ", b"DEFINITION\t ", 1),
            "unknown top-level": valid.replace(b"FEATURES", b"UNSUPPORTED\nFEATURES", 1),
            "misordered required": valid.replace(b"DEFINITION  Generated structural compiler input.\n", b"", 1),
            "duplicate required": valid.replace(b"ACCESSION   CC000003\n", b"ACCESSION   CC000003\nACCESSION   CC000003\n", 1),
            "reference numbering": valid.replace(b"REFERENCE   1", b"REFERENCE   2", 1),
            "wrong terminator": valid.replace(b"//\n", b"// \n"),
            "protein": valid.replace(b" bp DNA linear ", b" aa PRT linear ", 1),
            "RNA": valid.replace(b" bp DNA linear ", b" bp RNA linear ", 1),
            "CONTIG": valid.replace(b"ORIGIN\n", b"CONTIG     join(1..80)\n", 1),
            "uppercase origin": valid.replace(b"acgtacgtac", b"ACGTACGTAC", 1),
            "bad origin index": valid.replace(b"        1 ", b"        2 ", 1),
            "length mismatch": valid.replace(b" 80 bp ", b" 81 bp ", 1),
            "numeric feature key": valid.replace(b"     source          ", b"     123             ", 1),
            "numeric qualifier name": valid.replace(b"/organism=", b"/123=====", 1),
            "duplicate source mol_type": valid.replace(
                b'                     /mol_type="genomic DNA"\n',
                b'                     /mol_type="genomic DNA"\n'
                b'                     /mol_type="genomic DNA"\n',
                1,
            ),
            "blank after optional header": valid.replace(
                b"DBLINK      BioProject: PRJNA1\n",
                b"DBLINK      BioProject: PRJNA1\n\n",
                1,
            ),
            "blank after organism": valid.replace(
                b"            artificial sequences.\n",
                b"            artificial sequences.\n\n",
                1,
            ),
        }
        for label, raw in cases.items():
            with self.subTest(label=label):
                with self.assertRaises(GenBankError):
                    GenBankCompiler().compile_bytes(raw)

    def test_segment_uses_the_documented_pre_source_order_and_shape(self) -> None:
        valid = _record("CS000018", 1, "ACGT" * 20)
        segmented = valid.replace(
            b"KEYWORDS    .\n",
            b"KEYWORDS    .\nSEGMENT     1 of 2\n",
            1,
        )
        record = GenBankCompiler().compile_bytes(segmented).to_dict()["bio_ir"][
            "records"
        ][0]
        segment_headers = [
            header for header in record["headers"] if header["name"] == "SEGMENT"
        ]
        self.assertEqual([header["value"] for header in segment_headers], ["1 of 2"])

        malformed = (
            segmented.replace(b"SEGMENT     1 of 2", b"SEGMENT     0 of 2", 1),
            segmented.replace(b"SEGMENT     1 of 2", b"SEGMENT     2 of 1", 1),
            segmented.replace(b"SEGMENT     1 of 2", b"SEGMENT     1 of 1", 1),
            segmented.replace(b"SEGMENT     1 of 2", b"SEGMENT     1 OF 2", 1),
            segmented.replace(
                b"SEGMENT     1 of 2",
                b"SEGMENT     1 of\n            2",
                1,
            ),
            valid.replace(
                b"            artificial sequences.\n",
                b"            artificial sequences.\nSEGMENT     1 of 2\n",
                1,
            ),
        )
        for raw in malformed:
            with self.subTest(raw=raw[:180]):
                with self.assertRaises(GenBankError):
                    GenBankCompiler().compile_bytes(raw)

    def test_accession_version_identity_has_an_explicit_remote_reference_ceiling(self) -> None:
        raw = _record("AL000001", 1, "ACGT" * 20)
        accepted_accession = "A" + "B" * 60 + "1"
        accepted_version = accepted_accession + ".1"
        accepted = raw.replace(b"ACCESSION   AL000001", f"ACCESSION   {accepted_accession}".encode("ascii"), 1)
        accepted = accepted.replace(b"VERSION     AL000001.1", f"VERSION     {accepted_version}".encode("ascii"), 1)
        artifact = GenBankCompiler().compile_bytes(accepted)
        self.assertEqual(artifact.to_dict()["bio_ir"]["records"][0]["record_id"], accepted_version)

        rejected_accession = "A" + "B" * 61 + "1"
        rejected_version = rejected_accession + ".1"
        rejected = raw.replace(b"ACCESSION   AL000001", f"ACCESSION   {rejected_accession}".encode("ascii"), 1)
        rejected = rejected.replace(b"VERSION     AL000001.1", f"VERSION     {rejected_version}".encode("ascii"), 1)
        with self.assertRaisesRegex(GenBankError, "64 ASCII bytes"):
            GenBankCompiler().compile_bytes(rejected)

    def test_fixed_columns_qualifiers_source_coverage_and_limits_fail_closed(self) -> None:
        valid = _record("DD000004", 4, "ACGT" * 20)
        malformed = [
            valid.replace(b"     source          1..80", b"    source           1..80", 1),
            valid.replace(b"     source          1..80", b"     source          2..80", 1),
            valid.replace(
                b'                     /organism="Arbitrary organism"\n',
                b'                     /organism="Arbitrary organism"\n'
                b'                     /organism="Different organism"\n',
                1,
            ),
            valid.replace(b'/mol_type="genomic DNA"', b'/mol_type="unterminated', 1),
            valid.replace(b"FEATURES             Location/Qualifiers", b"FEATURES             Wrong", 1),
            valid.replace(
                b'/organism="Arbitrary organism"',
                b'/organism="' + b"x" * 60 + b'"',
                1,
            ),
        ]
        for raw in malformed:
            with self.subTest(raw=raw[:100]):
                with self.assertRaises(GenBankError):
                    GenBankCompiler().compile_bytes(raw)

    def test_multiline_qualifier_is_scanned_linearly(self) -> None:
        continuations = 600
        feature_lines = [
            "     misc_feature    1..80",
            '                     /note="first',
            *("                     continued" for _ in range(continuations)),
            '                     last"',
        ]
        raw = _record("EE000005", 5, "ACGT" * 20, features=feature_lines)
        with mock.patch.object(
            insdc_module,
            "_quoted_value",
            wraps=insdc_module._quoted_value,
        ) as observed:
            artifact = GenBankCompiler().compile_bytes(raw)

        scanned_parts = sum(len(call.args[0]) for call in observed.call_args_list)
        self.assertLess(scanned_parts, continuations * 3)
        note = artifact.to_dict()["bio_ir"]["records"][0]["features"][1]["qualifiers"][0]
        self.assertEqual(note["value"].count("continued"), continuations)

        structured = _record(
            "EF000005",
            5,
            "ACGT" * 20,
            features=[
                "     misc_feature    2..4",
                "                     /transl_except=(pos: 2..4,",
                "                     aa: Trp)",
            ],
        )
        qualifier = GenBankCompiler().compile_bytes(structured).to_dict()["bio_ir"][
            "records"
        ][0]["features"][1]["qualifiers"][0]
        self.assertEqual(qualifier["value"], "(pos: 2..4, aa: Trp)")

    def test_multiple_source_features_cover_by_union_and_may_overlap(self) -> None:
        valid = _record("FF000006", 6, "ACGT" * 20)
        original = (
            b"     source          1..80\n"
            b'                     /organism="Arbitrary organism"\n'
            b'                     /mol_type="genomic DNA"\n'
        )
        split_sources = (
            b"     source          1..40\n"
            b'                     /organism="First organism"\n'
            b'                     /mol_type="genomic DNA"\n'
            b"     source          41..80\n"
            b'                     /organism="Second organism"\n'
            b'                     /mol_type="genomic DNA"\n'
        )
        artifact = GenBankCompiler().compile_bytes(valid.replace(original, split_sources, 1))
        sources = [
            feature
            for feature in artifact.to_dict()["bio_ir"]["records"][0]["features"]
            if feature["key"] == "source"
        ]
        self.assertEqual(len(sources), 2)

        overlapping = split_sources + (
            b"     source          10..20\n"
            b'                     /organism="Third organism"\n'
            b'                     /mol_type="genomic DNA"\n'
        )
        GenBankCompiler().compile_bytes(valid.replace(original, overlapping, 1))

        gapped = split_sources.replace(b"41..80", b"42..80", 1)
        with self.assertRaisesRegex(GenBankError, "leave bases uncovered"):
            GenBankCompiler().compile_bytes(valid.replace(original, gapped, 1))

        external = _record(
            "FH000006",
            6,
            "ACGT" * 20,
            features=[
                "     misc_feature    ZZ999999.1:2..9",
                '                     /note="external bounds supplied later"',
            ],
        )
        external_location = GenBankCompiler().compile_bytes(external).to_dict()[
            "bio_ir"
        ]["records"][0]["features"][1]["location"]
        self.assertEqual(external_location["unresolved_references"], ["ZZ999999.1"])
        self.assertEqual(external_location["segments"][0]["bounds_status"], "unresolved")

        unresolved_circular = external.replace(b"ZZ999999.1:2..9", b"ZZ999999.1:10^1", 1)
        with self.assertRaisesRegex(GenBankError, "adjacent or explicit circular"):
            GenBankCompiler().compile_bytes(unresolved_circular)

    def test_wrapped_organism_reference_structure_and_origin_values(self) -> None:
        valid = _record("GG000007", 7, "ACGT" * 20)
        first_name_chunk = b"A" * 68
        wrapped = valid.replace(
            b"  ORGANISM  Arbitrary organism\n            artificial sequences.\n",
            b"  ORGANISM  "
            + first_name_chunk
            + b"\n            organism-name-continuation\n"
            + b"            Bacteria; Example lineage.\n",
            1,
        )
        record = GenBankCompiler().compile_bytes(wrapped).to_dict()["bio_ir"]["records"][0]
        self.assertEqual(
            record["organism"]["name"],
            first_name_chunk.decode("ascii") + " organism-name-continuation",
        )
        self.assertEqual(record["organism"]["lineage"], "Bacteria; Example lineage.")

        for value in (b"", b"Unreported.", b"EcoRI restriction site."):
            with self.subTest(value=value):
                artifact = GenBankCompiler().compile_bytes(
                    valid.replace(b"ORIGIN\n", b"ORIGIN      " + value + b"\n", 1)
                )
                origin = artifact.to_dict()["bio_ir"]["records"][0]["headers"][-1]
                self.assertEqual(origin["name"], "ORIGIN")
                self.assertEqual(origin["value"], value.decode("ascii"))

        without_reference = valid.replace(
            b"REFERENCE   1  (bases 1 to 80)\n"
            b"  AUTHORS   Example,A.\n"
            b"  TITLE     Direct Submission\n"
            b"  JOURNAL   Submitted (31-AUG-2026)\n",
            b"",
            1,
        )
        with self.assertRaisesRegex(GenBankError, "must contain a REFERENCE"):
            GenBankCompiler().compile_bytes(without_reference)
        with self.assertRaisesRegex(GenBankError, "exactly one JOURNAL"):
            GenBankCompiler().compile_bytes(
                valid.replace(b"  JOURNAL   Submitted (31-AUG-2026)\n", b"", 1)
            )
        for indentation in (b" ", b"            "):
            with self.subTest(journal_indentation=indentation):
                malformed_journal = valid.replace(
                    b"  JOURNAL   Submitted (31-AUG-2026)\n",
                    indentation + b"JOURNAL   Submitted (31-AUG-2026)\n",
                    1,
                )
                with self.assertRaisesRegex(GenBankError, "exactly one JOURNAL"):
                    GenBankCompiler().compile_bytes(malformed_journal)
        with self.assertRaisesRegex(GenBankError, "ending in a period"):
            GenBankCompiler().compile_bytes(
                valid.replace(b"ORIGIN\n", b"ORIGIN      Unreported\n", 1)
            )

        with_blank_comment = valid.replace(
            b"FEATURES             Location/Qualifiers\n",
            b"COMMENT     first paragraph\n"
            b"            \n"
            b"            second paragraph\n"
            b"FEATURES             Location/Qualifiers\n",
            1,
        )
        comment_record = GenBankCompiler().compile_bytes(with_blank_comment).to_dict()[
            "bio_ir"
        ]["records"][0]
        comment = next(
            header for header in comment_record["headers"] if header["name"] == "COMMENT"
        )
        self.assertEqual(comment["value"], "first paragraph second paragraph")
        self.assertEqual(
            comment["source_map"][1],
            {"line": 15, "column_start": 13, "column_end_exclusive": 13},
        )

    def test_single_blank_record_separator_is_skipped_without_renumbering(self) -> None:
        first = _record("AH000008", 1, "ACGT" * 20)
        second = _record("AI000009", 1, "TGCA" * 20)
        raw = first + b"\n" + second
        payload = GenBankCompiler().compile_bytes(raw).to_dict()
        records = payload["bio_ir"]["records"]
        self.assertEqual(records[0]["source_map"]["line_start"], 1)
        self.assertEqual(
            records[1]["source_map"]["line_start"],
            first.count(b"\n") + 2,
        )
        self.assertEqual(
            payload["sequence_collection"]["members"][1]["source_map"]["line_start"],
            records[1]["headers"][-1]["source_map"][0]["line"] + 1,
        )
        with self.assertRaises(GenBankError):
            GenBankCompiler().compile_bytes(first + b"\n\n" + second)
        with self.assertRaises(GenBankError):
            GenBankCompiler().compile_bytes(first + b" \n" + second)

    def test_record_length_can_exceed_the_legacy_one_megabyte_string_ceiling(self) -> None:
        length = insdc_module.MAX_STRING_BYTES + 1
        sequence = ("ACGT" * ((length + 3) // 4))[:length]
        artifact = GenBankCompiler().compile_bytes(
            _record("AJ000010", 1, sequence)
        )
        member = artifact.to_dict()["sequence_collection"]["members"][0]
        self.assertEqual(member["sequence"]["bases"], length)
        self.assertEqual(len(member["sequence"]["storage"]["chunks"]), 2)

    def test_production_frontend_contains_no_fixture_identity_tokens(self) -> None:
        source = (Path(__file__).parents[1] / "brainc" / "insdc.py").read_text(encoding="utf-8")
        for token in ("AA000001", "BB000002", "J02482", "5386"):
            self.assertNotIn(token, source)


if __name__ == "__main__":
    unittest.main()

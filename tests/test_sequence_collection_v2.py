from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from brainc._canonical import digest
from brainc.sequence_collection import SequenceCollectionError, load_sequence_collection
import brainc.sequence_collection_v2 as v2
from brainc.sequence_collection_v2 import (
    CHUNK_BYTES,
    SequenceCollectionV2Compiler,
    SequenceCollectionV2Error,
    SequenceMemberInput,
    load_sequence_collection_v2,
    validate_sequence_collection_v2,
)


def _origin(sequence: str) -> list[str]:
    lines: list[str] = []
    for offset in range(0, len(sequence), 60):
        part = sequence[offset : offset + 60].lower()
        groups = " ".join(part[index : index + 10] for index in range(0, len(part), 10))
        lines.append(f"{offset + 1:>9} {groups}")
    return lines


def _record(record_id: str, sequence: str, newline: str = "\n") -> tuple[bytes, SequenceMemberInput]:
    lines = [
        f"LOCUS       GENERATED {len(sequence)} bp DNA linear SYN 31-AUG-2026",
        f"VERSION     {record_id}",
        "ORIGIN",
        *_origin(sequence),
        "//",
    ]
    first_sequence_line = 4
    last_sequence_line = first_sequence_line + len(_origin(sequence)) - 1
    raw = (newline.join(lines) + newline).encode("ascii")
    return raw, SequenceMemberInput(record_id, sequence, first_sequence_line, last_sequence_line)


def _reseal(payload: dict[str, object]) -> None:
    payload["artifact_sha256"] = digest(
        {key: value for key, value in payload.items() if key != "artifact_sha256"}
    )


class SequenceCollectionV2Tests(unittest.TestCase):
    def test_small_exact_source_collection_round_trip_and_legacy_rejection(self) -> None:
        raw, member = _record("AR000001.1", "ACGT")
        artifact = SequenceCollectionV2Compiler().compile_genbank(raw, [member])
        payload = artifact.to_dict()

        self.assertEqual(payload["version"], 2)
        self.assertEqual(payload["inputs"]["source"]["byte_length"], len(raw))
        self.assertEqual(payload["inputs"]["source"]["sha256"], hashlib.sha256(raw).hexdigest())
        sequence = payload["members"][0]["sequence"]
        self.assertEqual(sequence["sha256"], hashlib.sha256(b"ACGT").hexdigest())
        self.assertEqual(sequence["refget_id"], "SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2")
        self.assertEqual(artifact.source_bytes(), raw)
        validate_sequence_collection_v2(payload)

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "collection.json"
            artifact.save(path)
            self.assertEqual(load_sequence_collection_v2(path).to_dict(), payload)
            with self.assertRaises(SequenceCollectionError):
                load_sequence_collection(path)

    def test_canonical_chunk_boundaries(self) -> None:
        for length, wanted_chunks in (
            (1, 1),
            (CHUNK_BYTES - 1, 1),
            (CHUNK_BYTES, 1),
            (CHUNK_BYTES + 1, 2),
        ):
            with self.subTest(length=length):
                storage = v2._storage(b"A" * length, "test")
                self.assertEqual(len(storage["chunks"]), wanted_chunks)
                self.assertEqual(v2._storage_bytes(storage, "test", maximum=CHUNK_BYTES * 2), b"A" * length)

        malformed = v2._storage(b"A" * (CHUNK_BYTES + 1), "test")
        malformed["chunks"][0] = malformed["chunks"][0][:-1]
        with self.assertRaisesRegex(SequenceCollectionV2Error, "short non-final"):
            v2._storage_bytes(malformed, "test", maximum=CHUNK_BYTES * 2)
        malformed = v2._storage(b"A", "test")
        malformed["chunks"][-1] = ""
        with self.assertRaisesRegex(SequenceCollectionV2Error, "empty or oversized final"):
            v2._storage_bytes(malformed, "test", maximum=CHUNK_BYTES * 2)

    def test_lf_and_crlf_bind_different_sources_but_identical_sequences(self) -> None:
        lf, lf_member = _record("AR000002.1", "ACGTRYSWKMBDHVN")
        crlf, crlf_member = _record("AR000002.1", "ACGTRYSWKMBDHVN", "\r\n")
        left = SequenceCollectionV2Compiler().compile_genbank(lf, [lf_member]).to_dict()
        right = SequenceCollectionV2Compiler().compile_genbank(crlf, [crlf_member]).to_dict()

        self.assertNotEqual(left["artifact_sha256"], right["artifact_sha256"])
        self.assertNotEqual(left["inputs"]["source"]["sha256"], right["inputs"]["source"]["sha256"])
        self.assertEqual(left["members"][0]["sequence"]["sha256"], right["members"][0]["sequence"]["sha256"])
        self.assertEqual(left["members"][0]["sequence"]["refget_id"], right["members"][0]["sequence"]["refget_id"])

    def test_coherently_resealed_sequence_map_and_chunk_forgery_fail(self) -> None:
        raw, member = _record("AR000003.1", "ACGT" * 20)
        payload = SequenceCollectionV2Compiler().compile_genbank(raw, [member]).to_dict()

        sequence_forgery = copy.deepcopy(payload)
        sequence = sequence_forgery["members"][0]["sequence"]
        sequence["storage"]["chunks"][0] = "T" + sequence["storage"]["chunks"][0][1:]
        forged_bytes = sequence["storage"]["chunks"][0].encode("ascii")
        sequence["sha256"] = hashlib.sha256(forged_bytes).hexdigest()
        sequence["refget_id"] = v2._refget(forged_bytes)
        sequence_forgery["refget_seqcol"] = v2._refget_seqcol(sequence_forgery["members"])
        _reseal(sequence_forgery)
        with self.assertRaisesRegex(SequenceCollectionV2Error, "source lines"):
            validate_sequence_collection_v2(sequence_forgery)

        map_forgery = copy.deepcopy(payload)
        map_forgery["members"][0]["source_map"]["line_start"] -= 1
        _reseal(map_forgery)
        with self.assertRaises(SequenceCollectionV2Error):
            validate_sequence_collection_v2(map_forgery)

        chunk_forgery = copy.deepcopy(payload)
        chunk_forgery["members"][0]["sequence"]["storage"]["chunk_bytes"] = CHUNK_BYTES - 1
        _reseal(chunk_forgery)
        with self.assertRaisesRegex(SequenceCollectionV2Error, "chunk contract"):
            validate_sequence_collection_v2(chunk_forgery)

    def test_source_record_identity_length_order_and_completeness_are_bound(self) -> None:
        raw, member = _record("AR000004.1", "ACGT" * 20)
        wrong_id = SequenceMemberInput("AR000005.1", member.sequence, member.line_start, member.line_end)
        with self.assertRaisesRegex(SequenceCollectionV2Error, "VERSION"):
            SequenceCollectionV2Compiler().compile_genbank(raw, [wrong_id])

        wrong_length = raw.replace(b" 80 bp ", b" 81 bp ", 1)
        with self.assertRaisesRegex(SequenceCollectionV2Error, "LOCUS length"):
            SequenceCollectionV2Compiler().compile_genbank(wrong_length, [member])

        second_raw, second_member = _record("AR000006.1", "TGCA" * 20)
        combined = raw + second_raw
        offset = raw.count(b"\n")
        shifted_second = SequenceMemberInput(
            second_member.record_id,
            second_member.sequence,
            second_member.line_start + offset,
            second_member.line_end + offset,
        )
        with self.assertRaisesRegex(SequenceCollectionV2Error, "every source LOCUS"):
            SequenceCollectionV2Compiler().compile_genbank(combined, [member])
        collection = SequenceCollectionV2Compiler().compile_genbank(combined, [member, shifted_second])
        self.assertEqual([item["record_id"] for item in collection.to_dict()["members"]], ["AR000004.1", "AR000006.1"])

    def test_source_envelope_and_molecule_profile_fail_closed(self) -> None:
        raw, member = _record("AR000008.1", "ACGT" * 20)
        malformed = {
            "amino acids": raw.replace(b" bp DNA ", b" aa DNA ", 1),
            "RNA": raw.replace(b" bp DNA ", b" bp RNA ", 1),
            "ORIGIN text": raw.replace(b"ORIGIN\n", b"ORIGIN      malicious header\n", 1),
            "trailing junk": raw + b"junk\n",
            "extra terminal blanks": raw + b"\n\n",
            "bare final carriage return": raw[:-1] + b"\r",
            "tab in physical line": raw.replace(b"GENERATED", b"GENERATED\t", 1),
        }
        for label, source in malformed.items():
            with self.subTest(label=label):
                with self.assertRaises(SequenceCollectionV2Error):
                    SequenceCollectionV2Compiler().compile_genbank(source, [member])

        leading = b"junk\n" + raw
        shifted = SequenceMemberInput(
            member.record_id,
            member.sequence,
            member.line_start + 1,
            member.line_end + 1,
        )
        with self.assertRaisesRegex(SequenceCollectionV2Error, "begin"):
            SequenceCollectionV2Compiler().compile_genbank(leading, [shifted])

        official_terminal_blank = raw + b"\n"
        SequenceCollectionV2Compiler().compile_genbank(official_terminal_blank, [member])

    def test_refget_types_and_deep_input_are_rejected_before_digest_recursion(self) -> None:
        raw, member = _record("AR000009.1", "A")
        payload = SequenceCollectionV2Compiler().compile_genbank(raw, [member]).to_dict()
        bool_length = copy.deepcopy(payload)
        bool_length["refget_seqcol"]["level_2"]["lengths"][0] = True
        _reseal(bool_length)
        with self.assertRaisesRegex(SequenceCollectionV2Error, "level_2"):
            validate_sequence_collection_v2(bool_length)

        deep = copy.deepcopy(payload)
        nested: dict[str, object] = {}
        deep["unexpected"] = nested
        for _ in range(v2.MAX_JSON_DEPTH + 2):
            child: dict[str, object] = {}
            nested["child"] = child
            nested = child
        with self.assertRaisesRegex(SequenceCollectionV2Error, "depth"):
            validate_sequence_collection_v2(deep)

    def test_large_sequence_uses_a_constant_number_of_canonical_chunks(self) -> None:
        length = CHUNK_BYTES + 1
        sequence = ("ACGT" * ((length + 3) // 4))[:length]
        raw, member = _record("AR000007.1", sequence)
        artifact = SequenceCollectionV2Compiler().compile_genbank(raw, [member])
        payload = artifact.to_dict()
        self.assertEqual(len(payload["members"][0]["sequence"]["storage"]["chunks"]), 2)
        self.assertEqual(payload["members"][0]["sequence"]["bases"], length)
        self.assertLess(len(json.dumps(payload).encode("utf-8")), v2.MAX_ARTIFACT_BYTES)


if __name__ == "__main__":
    unittest.main()

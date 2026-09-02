from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from brainc.source import (
    EXTERNAL_PROFILE,
    FASTA_PROFILE,
    GFF3_PROFILE,
    GENBANK_PROFILE,
    PROFILES,
    RAW_PROFILE,
    REFERENCE_FASTA_PROFILE,
    compile_external_source_executable,
    compile_source,
)
from tests.test_external_execution import FRONTEND, VALIDATOR, _manifest
from tests.test_external_profile import _refget, _sha
from tests.test_insdc import _record
from tests.test_validator_development import _Fixture


RECORD_ID = "SAMP0001.1"
SEQUENCE = b"ACGTRYSWKMBDHVN"


class UniversalTranslationTests(unittest.TestCase):
    def test_every_builtin_and_an_executable_external_grammar_reach_development(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sources = self._sources(root)
            records = [source.to_dict()["source_ir"]["records"] for source, _, _ in sources]
            self.assertEqual({source.profile for source, _, _ in sources}, PROFILES)
            self.assertTrue(all(value == records[0] for value in records[1:]))
            self.assertEqual(
                records[0],
                [
                    {
                        "record_id": RECORD_ID,
                        "bases": len(SEQUENCE),
                        "sequence_sha256": _sha(SEQUENCE),
                        "refget_id": _refget(SEQUENCE),
                    }
                ],
            )
            for index, (source, inputs, validator) in enumerate(sources):
                with self.subTest(profile=source.profile):
                    fixture = _Fixture(
                        root / f"development-{index}",
                        inputs=inputs,
                        source=source,
                        external_validator=validator,
                    )
                    report = fixture.validate_paths()
                    self.assertTrue(report["valid"])
                    self.assertEqual(
                        fixture.compiled.artifacts["development_module"]["format"],
                        "brainc.development-module",
                    )

    def _sources(self, root: Path):
        raw = root / "sample.dna"
        raw.write_bytes(SEQUENCE.lower())
        fasta = root / "sample.fasta"
        fasta.write_bytes(b">" + RECORD_ID.encode() + b"\n" + SEQUENCE.lower() + b"\n")
        genbank = root / "sample.gb"
        genbank.write_bytes(_record("SAMP0001", 1, SEQUENCE.decode()))
        annotation = root / "sample.gff3"
        annotation.write_bytes(
            b"##gff-version 3\n"
            + RECORD_ID.encode()
            + b"\t.\tregion\t1\t15\t.\t+\t.\tID=whole\n"
        )
        fastq = root / "sample.fastq"
        fastq.write_bytes(
            b"@" + RECORD_ID.encode() + b"\n" + SEQUENCE.lower() + b"\n+\n" + b"I" * len(SEQUENCE) + b"\n"
        )
        external = compile_external_source_executable(
            _manifest(),
            original_paths={"reads": fastq},
            frontend_executable=FRONTEND,
            validator_executable=VALIDATOR,
        )
        return [
            (
                compile_source(RAW_PROFILE, {"sequence": raw}, parameters={"record_id": RECORD_ID}),
                {"sequence": raw},
                None,
            ),
            (
                compile_source(FASTA_PROFILE, {"sequence": fasta}, parameters={"wrapper": "identity"}),
                {"sequence": fasta},
                None,
            ),
            (
                compile_source(REFERENCE_FASTA_PROFILE, {"sequence": fasta}, parameters={"wrapper": "identity"}),
                {"sequence": fasta},
                None,
            ),
            (
                compile_source(GENBANK_PROFILE, {"genbank": genbank}, parameters={}),
                {"genbank": genbank},
                None,
            ),
            (
                compile_source(
                    GFF3_PROFILE,
                    {"sequence": fasta, "annotation": annotation},
                    parameters={
                        "sequence_profile": FASTA_PROFILE,
                        "sequence_parameters": {"wrapper": "identity"},
                    },
                ),
                {"sequence": fasta, "annotation": annotation},
                None,
            ),
            (external, {"reads": fastq}, VALIDATOR),
        ]


if __name__ == "__main__":
    unittest.main()

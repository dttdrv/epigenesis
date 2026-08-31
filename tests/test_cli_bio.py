"""End-to-end tests for the public v0.6 biological compiler CLI."""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class BioCompilerCliTests(unittest.TestCase):
    def _run(
        self, *arguments: str | Path, expected: int = 0
    ) -> subprocess.CompletedProcess[str]:
        result = subprocess.run(
            [sys.executable, "-m", "brainc", *(str(value) for value in arguments)],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(
            result.returncode,
            expected,
            msg=f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}",
        )
        return result

    def _compile_bio(self, root: Path) -> dict[str, Path]:
        fasta = root / "input.fasta"
        fasta.write_bytes(b">arbitrary\nACGTACGTACGT\n")
        gff3 = root / "annotation.gff3"
        gff3.write_bytes(
            b"##gff-version 3\n"
            b"arbitrary\t.\tgene\t1\t12\t.\t+\t.\tID=parent\n"
            b"arbitrary\t.\texon\t2\t5\t.\t+\t.\tID=child;Parent=parent\n"
        )
        source = root / "sequence-collection.json"
        bio_ir = root / "bio-ir.json"
        self._run("compile-collection", fasta, "-o", source)
        self._run(
            "compile-gff3",
            gff3,
            "--sequence-collection",
            source,
            "-o",
            bio_ir,
        )
        return {"fasta": fasta, "gff3": gff3, "source": source, "bio_ir": bio_ir}

    def test_complete_gff3_and_feature_graph_round_trip_with_reports(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = self._compile_bio(root)

            bio_report_path = root / "gff3-validation.json"
            result = self._run(
                "validate-gff3",
                paths["bio_ir"],
                "--sequence-collection",
                paths["source"],
                "--sequence-input",
                paths["fasta"],
                "--gff3-source",
                paths["gff3"],
                "--report",
                bio_report_path,
            )
            bio_report = json.loads(result.stdout)
            self.assertTrue(bio_report["valid"])
            self.assertEqual(bio_report["summary"]["features"], 2)
            self.assertEqual(bio_report["summary"]["relationships"], 1)
            self.assertEqual(
                json.loads(bio_report_path.read_text(encoding="utf-8")), bio_report
            )

            graph_directory = root / "feature-graph"
            compile_result = self._run(
                "compile-feature-graph",
                paths["bio_ir"],
                "--sequence-collection",
                paths["source"],
                "--gff3-source",
                paths["gff3"],
                "-o",
                graph_directory,
            )
            compile_summary = json.loads(compile_result.stdout)
            bundle_path = Path(compile_summary["bundle"])
            self.assertEqual(bundle_path, graph_directory / "bundle.json")
            self.assertTrue(bundle_path.is_file())
            self.assertEqual(len(list(graph_directory.glob("*.json"))), 10)

            graph_report_path = root / "feature-graph-validation.json"
            result = self._run(
                "validate-feature-graph",
                bundle_path,
                "--sequence-collection",
                paths["source"],
                "--bio-ir",
                paths["bio_ir"],
                "--report",
                graph_report_path,
            )
            graph_report = json.loads(result.stdout)
            self.assertTrue(graph_report["valid"])
            self.assertEqual(graph_report["result"]["units"], 2)
            self.assertEqual(graph_report["result"]["edges"], 1)
            self.assertEqual(graph_report["result"]["relationship_edges"], 1)
            self.assertEqual(
                json.loads(graph_report_path.read_text(encoding="utf-8")),
                graph_report,
            )

    def test_source_replay_and_bundle_failures_do_not_write_success_reports(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = self._compile_bio(root)
            changed_gff3 = root / "changed.gff3"
            changed_gff3.write_bytes(
                paths["gff3"].read_bytes().replace(b"\texon\t", b"\tCDS\t")
            )
            failed_bio_report = root / "must-not-exist-bio.json"
            failure = self._run(
                "validate-gff3",
                paths["bio_ir"],
                "--sequence-collection",
                paths["source"],
                "--sequence-input",
                paths["fasta"],
                "--gff3-source",
                changed_gff3,
                "--report",
                failed_bio_report,
                expected=2,
            )
            self.assertIn("brainc:", failure.stderr)
            self.assertFalse(failed_bio_report.exists())

            graph_directory = root / "graph"
            self._run(
                "compile-feature-graph",
                paths["bio_ir"],
                "--sequence-collection",
                paths["source"],
                "-o",
                graph_directory,
            )
            bundle_path = graph_directory / "bundle.json"
            bundle = json.loads(bundle_path.read_text(encoding="utf-8"))
            bundle["artifacts"]["compilation_record"]["result"]["units"] += 1
            forged = root / "forged-bundle.json"
            forged.write_text(json.dumps(bundle), encoding="utf-8")
            failed_graph_report = root / "must-not-exist-graph.json"
            failure = self._run(
                "validate-feature-graph",
                forged,
                "--sequence-collection",
                paths["source"],
                "--bio-ir",
                paths["bio_ir"],
                "--report",
                failed_graph_report,
                expected=2,
            )
            self.assertIn("brainc:", failure.stderr)
            self.assertFalse(failed_graph_report.exists())

    def test_embedded_fasta_is_rejected_by_external_sequence_profile(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = self._compile_bio(root)
            embedded = root / "embedded.gff3"
            embedded.write_bytes(
                paths["gff3"].read_bytes() + b"##FASTA\n>arbitrary\nACGT\n"
            )
            output = root / "embedded.json"
            failure = self._run(
                "compile-gff3",
                embedded,
                "--sequence-collection",
                paths["source"],
                "-o",
                output,
                expected=2,
            )
            self.assertIn("embedded FASTA", failure.stderr)
            self.assertFalse(output.exists())

    def test_compile_raw_uses_safe_file_ingress(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            raw = root / "raw.txt"
            raw.write_bytes(b"ACGTRYSWKMBDHVN")
            output = root / "raw.json"
            self._run(
                "compile-raw", raw, "--record-id", "raw-record", "-o", output
            )
            payload = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(payload["members"][0]["record_id"], "raw-record")

            if hasattr(os, "symlink"):
                linked = root / "linked-raw.txt"
                linked.symlink_to(raw)
                failure = self._run(
                    "compile-raw",
                    linked,
                    "--record-id",
                    "raw-record",
                    "-o",
                    root / "linked-output.json",
                    expected=2,
                )
                self.assertIn("regular non-linked file", failure.stderr)


if __name__ == "__main__":
    unittest.main()

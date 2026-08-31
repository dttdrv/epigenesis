from __future__ import annotations

import ast
from copy import deepcopy
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import tomllib
import unittest
from unittest import mock

from brainc.bio import GFF3Compiler
from brainc.bio_graph import compile_feature_graph
from brainc.sequence_collection import SequenceCollectionCompiler
import brainc.validator_bio_chain as chain_validator
from brainc.validator_bio_chain import (
    BioChainValidationError,
    validate_bio_feature_graph_paths,
    validate_chain_report,
)


ROOT = Path(__file__).resolve().parents[1]


class StandaloneBioChainValidatorTests(unittest.TestCase):
    def _chain(self, root: Path) -> dict[str, Path]:
        root.mkdir(parents=True, exist_ok=True)
        fasta = root / "input.fasta"
        fasta.write_bytes(b">arbitrary\nACGTACGTACGT\n")
        gff3 = root / "input.gff3"
        gff3.write_bytes(
            b"##gff-version 3\n"
            b"arbitrary\t.\tgene\t1\t12\t.\t+\t.\tID=parent\n"
            b"arbitrary\t.\texon\t2\t5\t.\t+\t.\tID=child;Parent=parent\n"
        )
        collection = SequenceCollectionCompiler().compile_file(fasta)
        collection_path = root / "sequence-collection.json"
        collection.save(collection_path)
        bio = GFF3Compiler().compile_file(gff3, collection)
        bio_path = root / "bio-ir.json"
        bio.save(bio_path)
        bundle = compile_feature_graph(
            collection_path,
            bio,
            gff3_source=gff3.read_bytes(),
        )
        bundle_path = bundle.save(root / "feature-graph")["bundle"]
        return {
            "fasta": fasta,
            "gff3": gff3,
            "collection": collection_path,
            "bio": bio_path,
            "bundle": bundle_path,
        }

    def _validate(self, paths: dict[str, Path]) -> dict:
        return validate_bio_feature_graph_paths(
            sequence_input=paths["fasta"],
            gff3_source=paths["gff3"],
            sequence_collection=paths["collection"],
            bio_artifact=paths["bio"],
            feature_graph_bundle=paths["bundle"],
        )

    def test_one_snapshot_per_input_and_deterministic_closed_report(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = self._chain(Path(temporary))
            original = chain_validator._read_regular
            reads: list[Path] = []

            def observe(path, label, maximum):
                reads.append(Path(path))
                return original(path, label, maximum)

            with mock.patch.object(
                chain_validator, "_read_regular", side_effect=observe
            ):
                first = self._validate(paths)
            second = self._validate(paths)

            self.assertEqual(reads, list(paths.values()))
            self.assertEqual(first, second)
            self.assertIs(validate_chain_report(first), first)
            self.assertTrue(first["valid"])
            self.assertEqual(first["result"]["features"], 2)
            self.assertEqual(first["result"]["relationships"], 1)
            self.assertEqual(first["result"]["units"], 2)
            self.assertEqual(first["result"]["edges"], 1)

            mutable = deepcopy(first)
            mutable["validator"]["name"] = "attacker-chain"
            mutable["stages"]["source_to_bio_ir"]["validator"]["name"] = "attacker-bio"
            self.assertEqual(
                chain_validator.VALIDATOR["name"],
                "brainc-independent-bio-chain-validator",
            )
            self.assertEqual(
                chain_validator.validator_bio.VALIDATOR["name"],
                "brainc-independent-bio-validator",
            )
            with self.assertRaises(BioChainValidationError):
                validate_chain_report(mutable)

            forged = deepcopy(first)
            forged["result"]["edges"] = 2
            with self.assertRaisesRegex(
                BioChainValidationError,
                "combined result does not match",
            ):
                validate_chain_report(forged)

    def test_returned_report_identity_cannot_mutate_validator_contract(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = self._chain(Path(temporary))
            expected_chain_identity = dict(chain_validator.VALIDATOR)
            expected_bio_identity = dict(chain_validator.validator_bio.VALIDATOR)

            forged = self._validate(paths)
            forged["validator"]["name"] = "attacker-chain-validator"
            forged["report_sha256"] = chain_validator.validator_bio_graph.digest(
                {key: value for key, value in forged.items() if key != "report_sha256"}
            )
            self.assertEqual(dict(chain_validator.VALIDATOR), expected_chain_identity)
            with self.assertRaisesRegex(
                BioChainValidationError,
                "combined validation report identity",
            ):
                validate_chain_report(forged)

            forged = self._validate(paths)
            biological = forged["stages"]["source_to_bio_ir"]
            biological["validator"]["name"] = "attacker-bio-validator"
            biological["report_sha256"] = chain_validator.validator_bio_graph.digest(
                {
                    key: value
                    for key, value in biological.items()
                    if key != "report_sha256"
                }
            )
            forged["report_sha256"] = chain_validator.validator_bio_graph.digest(
                {key: value for key, value in forged.items() if key != "report_sha256"}
            )
            self.assertEqual(
                dict(chain_validator.validator_bio.VALIDATOR),
                expected_bio_identity,
            )
            with self.assertRaisesRegex(ValueError, "validation report identity"):
                validate_chain_report(forged)

    def test_actual_ncbi_record_passes_the_complete_public_chain(self) -> None:
        fasta = ROOT / "tests" / "data" / "J02482.1.fasta"
        gff3 = ROOT / "tests" / "data" / "J02482.1.gff3"
        metadata = ROOT / "tests" / "data" / "J02482.1.source.json"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            collection = SequenceCollectionCompiler().compile_file(fasta)
            collection_path = root / "sequence-collection.json"
            collection.save(collection_path)
            bio = GFF3Compiler().compile_file(gff3, collection)
            bio_path = root / "bio-ir.json"
            bio.save(bio_path)
            bundle_path = compile_feature_graph(
                collection_path,
                bio,
                gff3_source=gff3.read_bytes(),
            ).save(root / "feature-graph")["bundle"]

            report = validate_bio_feature_graph_paths(
                sequence_input=fasta,
                gff3_source=gff3,
                sequence_collection=collection_path,
                bio_artifact=bio_path,
                feature_graph_bundle=bundle_path,
                source_metadata=metadata,
            )
            self.assertEqual(report["result"]["sequence_bases"], 5386)
            self.assertEqual(report["result"]["features"], 23)
            self.assertEqual(report["result"]["relationships"], 4)
            self.assertEqual(report["result"]["units"], 23)
            self.assertEqual(report["result"]["edges"], 4)
            self.assertTrue(
                report["stages"]["source_to_bio_ir"]["summary"][
                    "source_metadata_verified"
                ]
            )

    def test_isolated_public_executable_loads_no_compiler_modules(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = self._chain(root / "inputs")
            installed = root / "installed"
            shutil.copytree(ROOT / "brainc", installed / "brainc")
            report_path = root / "combined-report.json"
            allowed = {
                "brainc",
                "brainc.validator_bio",
                "brainc.validator_bio_chain",
                "brainc.validator_bio_graph",
            }
            command = (
                "import sys; "
                f"sys.path.insert(0, {str(installed)!r}); "
                "from brainc.validator_bio_chain import main; "
                f"code=main({[str(paths[name]) for name in ('fasta', 'gff3', 'collection', 'bio', 'bundle')] + ['--report', str(report_path)]!r}); "
                f"allowed={allowed!r}; "
                "loaded={name for name in sys.modules if name.startswith('brainc')}; "
                "assert loaded == allowed, sorted(loaded); "
                "raise SystemExit(code)"
            )
            result = subprocess.run(
                [sys.executable, "-I", "-c", command],
                cwd=root,
                check=False,
                capture_output=True,
            )
            self.assertEqual(
                result.returncode,
                0,
                msg=(result.stdout + result.stderr).decode("utf-8", errors="replace"),
            )
            self.assertEqual(report_path.read_bytes(), result.stdout)
            self.assertTrue(json.loads(result.stdout)["valid"])

    def test_cli_tamper_fails_closed_without_writing_report(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = self._chain(root)
            paths["gff3"].write_bytes(
                paths["gff3"].read_bytes().replace(b"\texon\t", b"\tCDS\t")
            )
            report = root / "must-not-exist.json"
            result = subprocess.run(
                [
                    sys.executable,
                    "-I",
                    "-c",
                    (
                        "import sys; "
                        f"sys.path.insert(0, {str(ROOT)!r}); "
                        "from brainc.validator_bio_chain import main; "
                        f"raise SystemExit(main({[str(paths[name]) for name in ('fasta', 'gff3', 'collection', 'bio', 'bundle')] + ['--report', str(report)]!r}))"
                    ),
                ],
                cwd=root,
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertIn("validation failed:", result.stderr)
            self.assertFalse(report.exists())

    def test_imports_and_console_entry_are_public_and_isolated(self) -> None:
        source = (ROOT / "brainc" / "validator_bio_chain.py").read_text(
            encoding="utf-8"
        )
        tree = ast.parse(source)
        imported_brainc_modules = {
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, ast.Import)
            for alias in node.names
            if alias.name.startswith("brainc")
        }
        imported_from_brainc = {
            node.module
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom)
            and node.module is not None
            and node.module.startswith("brainc")
        }
        self.assertEqual(
            imported_brainc_modules,
            {"brainc.validator_bio", "brainc.validator_bio_graph"},
        )
        self.assertEqual(imported_from_brainc, set())
        project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
        self.assertEqual(
            project["project"]["scripts"]["brainc-validate-bio-chain"],
            "brainc.validator_bio_chain:main",
        )

    def test_input_symlink_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = self._chain(root)
            linked = root / "linked.fasta"
            try:
                linked.symlink_to(paths["fasta"])
            except (OSError, NotImplementedError):
                self.skipTest("symbolic links unavailable")
            paths["fasta"] = linked
            with self.assertRaisesRegex(
                BioChainValidationError,
                "regular non-linked file",
            ):
                self._validate(paths)


if __name__ == "__main__":
    unittest.main()

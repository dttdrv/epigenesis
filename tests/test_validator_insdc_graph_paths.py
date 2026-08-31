from __future__ import annotations

import copy
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
import tempfile
import tomllib
import unittest
from unittest import mock

import brainc.validator_insdc_graph as validator
from brainc.insdc import GenBankCompiler
from brainc.insdc_graph import compile_insdc_graph


ROOT = Path(__file__).resolve().parents[1]


class IndependentINSDCGraphPathTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.raw = (ROOT / "tests/data/U49845.1.gb").read_bytes()
        cls.source = GenBankCompiler().compile_bytes(cls.raw).to_dict()
        cls.bundle = compile_insdc_graph(cls.source)

    def _paths(self, root: Path) -> tuple[Path, Path, Path]:
        genbank = root / "U49845.1.gb"
        artifact = root / "U49845.1.genbank.json"
        bundle = root / "feature-state"
        genbank.write_bytes(self.raw)
        artifact.write_text(
            json.dumps(self.source, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        self.bundle.save(bundle)
        return genbank, artifact, bundle

    @unittest.skipUnless(
        validator._HAS_DIRECTORY_DESCRIPTOR,
        "directory descriptors are unavailable",
    )
    def test_complete_loader_path_api_cli_and_atomic_report(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genbank, artifact, bundle = self._paths(root)
            index, children = validator.load_insdc_graph_bundle_directory(bundle)
            self.assertEqual(index, self.bundle.to_dict())
            self.assertEqual(children, self.bundle.artifacts)

            report_path = root / "report.json"
            report = validator.validate_insdc_graph_paths(
                genbank,
                artifact,
                bundle,
                report_path=report_path,
            )
            self.assertTrue(report["valid"])
            self.assertEqual(
                json.loads(report_path.read_text(encoding="utf-8")),
                report,
            )
            self.assertEqual(validator.validate_insdc_graph_report(report), report)

            completed = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "brainc.validator_insdc_graph",
                    str(genbank),
                    str(artifact),
                    str(bundle),
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertEqual(json.loads(completed.stdout), report)
            for option in ("-o", "--output", "--report"):
                output = root / f"{option.lstrip('-') or 'short'}.json"
                completed = subprocess.run(
                    [
                        sys.executable,
                        "-m",
                        "brainc.validator_insdc_graph",
                        str(genbank),
                        str(artifact),
                        str(bundle),
                        option,
                        str(output),
                    ],
                    cwd=ROOT,
                    check=True,
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(completed.stdout, "")
                self.assertEqual(json.loads(output.read_text(encoding="utf-8")), report)
            scripts = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"][
                "scripts"
            ]
            self.assertEqual(
                scripts["brainc-validate-insdc-graph"],
                "brainc.validator_insdc_graph:main",
            )

    @unittest.skipUnless(
        validator._HAS_DIRECTORY_DESCRIPTOR,
        "directory descriptors are unavailable",
    )
    def test_raw_source_and_child_tampering_publish_no_success_report(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genbank, artifact, bundle = self._paths(root)
            genbank.write_bytes(
                self.raw.replace(b"Saccharomyces", b"Schizosacchar", 1)
            )
            report = root / "report.json"
            with self.assertRaisesRegex(
                validator.INSDCGraphValidationError,
                "independent GenBank source replay failed",
            ):
                validator.validate_insdc_graph_paths(
                    genbank,
                    artifact,
                    bundle,
                    report_path=report,
                )
            self.assertFalse(report.exists())

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genbank, artifact, bundle = self._paths(root)
            semantic_path = bundle / "backend_semantics.json"
            semantic = json.loads(semantic_path.read_text(encoding="utf-8"))
            semantic["dictionaries"]["keys"].reverse()
            semantic["artifact_sha256"] = validator.digest(
                {
                    key: value
                    for key, value in semantic.items()
                    if key != "artifact_sha256"
                }
            )
            semantic_path.write_text(json.dumps(semantic), encoding="utf-8")
            report = root / "report.json"
            report.write_bytes(b"preserve-existing-report")
            with self.assertRaises(validator.INSDCGraphValidationError):
                validator.validate_insdc_graph_paths(
                    genbank,
                    artifact,
                    bundle,
                    report_path=report,
                )
            self.assertEqual(report.read_bytes(), b"preserve-existing-report")

    @unittest.skipUnless(
        validator._HAS_DIRECTORY_DESCRIPTOR,
        "directory descriptors are unavailable",
    )
    def test_closed_directory_rejects_unknown_missing_and_links(self) -> None:
        for case in ("unknown", "missing", "linked-child", "linked-directory"):
            with self.subTest(case=case), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                _, _, bundle = self._paths(root)
                if case == "unknown":
                    (bundle / "extra.json").write_text("{}", encoding="utf-8")
                    target = bundle
                    pattern = "unknown"
                elif case == "missing":
                    (bundle / "target_contract.json").unlink()
                    target = bundle
                    pattern = "missing"
                elif case == "linked-child":
                    child = bundle / "target_contract.json"
                    victim = root / "victim.json"
                    child.rename(victim)
                    child.symlink_to(victim)
                    target = bundle
                    pattern = "regular non-linked file"
                else:
                    target = root / "linked-bundle"
                    target.symlink_to(bundle, target_is_directory=True)
                    pattern = "non-linked directory"
                with self.assertRaisesRegex(
                    validator.INSDCGraphValidationError,
                    pattern,
                ):
                    validator.load_insdc_graph_bundle_directory(target)

    @unittest.skipUnless(
        validator._HAS_DIRECTORY_DESCRIPTOR,
        "directory descriptors are unavailable",
    )
    def test_duplicate_json_keys_fail_at_every_path_boundary(self) -> None:
        cases = (
            ("bundle.json", "bundle index"),
            ("target_contract.json", "child artifact target_contract"),
        )
        for filename, label in cases:
            with self.subTest(filename=filename), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                _, _, bundle = self._paths(root)
                path = bundle / filename
                raw = path.read_bytes().replace(
                    b"{",
                    b'{"format":"forged",',
                    1,
                )
                path.write_bytes(raw)
                with self.assertRaisesRegex(
                    validator.INSDCGraphValidationError,
                    rf"{label} contains duplicate JSON key 'format'",
                ):
                    validator.load_insdc_graph_bundle_directory(bundle)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genbank, artifact, bundle = self._paths(root)
            artifact.write_bytes(
                artifact.read_bytes().replace(
                    b"{",
                    b'{"format":"forged",',
                    1,
                )
            )
            with self.assertRaisesRegex(
                validator.INSDCGraphValidationError,
                "compiled GenBank artifact contains duplicate JSON key 'format'",
            ):
                validator.validate_insdc_graph_paths(genbank, artifact, bundle)

    @unittest.skipUnless(
        validator._HAS_DIRECTORY_DESCRIPTOR,
        "directory descriptors are unavailable",
    )
    def test_exact_wire_ceilings_and_plus_one(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            exact_source = root / "exact-source.json"
            with exact_source.open("wb") as stream:
                stream.truncate(validator.MAX_SOURCE_ARTIFACT_BYTES)
            self.assertEqual(
                len(
                    validator._read_regular(
                        exact_source,
                        "compiled GenBank artifact",
                        validator.MAX_SOURCE_ARTIFACT_BYTES,
                    )
                ),
                validator.MAX_SOURCE_ARTIFACT_BYTES,
            )
            with exact_source.open("r+b") as stream:
                stream.truncate(validator.MAX_SOURCE_ARTIFACT_BYTES + 1)
            with self.assertRaisesRegex(
                validator.INSDCGraphValidationError,
                "exceeds 67108864 bytes",
            ):
                validator._read_regular(
                    exact_source,
                    "compiled GenBank artifact",
                    validator.MAX_SOURCE_ARTIFACT_BYTES,
                )

            _, _, bundle = self._paths(root)
            index_path = bundle / "bundle.json"
            index_raw = index_path.read_bytes()
            index_path.write_bytes(
                index_raw + b" " * (validator.MAX_CHILD_BYTES - len(index_raw))
            )
            validator.load_insdc_graph_bundle_directory(bundle)
            with index_path.open("ab") as stream:
                stream.write(b" ")
            with self.assertRaisesRegex(
                validator.INSDCGraphValidationError,
                "exceeds 16777216 bytes",
            ):
                validator.load_insdc_graph_bundle_directory(bundle)

    @unittest.skipUnless(
        validator._HAS_DIRECTORY_DESCRIPTOR and hasattr(os, "symlink"),
        "secure symbolic-link checks are unavailable",
    )
    def test_report_output_never_follows_symlink(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genbank, artifact, bundle = self._paths(root)
            victim = root / "victim"
            victim.write_bytes(b"preserved")
            linked = root / "report.json"
            linked.symlink_to(victim)
            with self.assertRaisesRegex(
                validator.INSDCGraphValidationError,
                "regular non-linked file",
            ):
                validator.validate_insdc_graph_paths(
                    genbank,
                    artifact,
                    bundle,
                    report_path=linked,
                )
            self.assertTrue(linked.is_symlink())
            self.assertEqual(victim.read_bytes(), b"preserved")

    @unittest.skipUnless(
        validator._HAS_DIRECTORY_DESCRIPTOR,
        "directory descriptors are unavailable",
    )
    def test_report_temp_substitution_rolls_back_and_cleans_up(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genbank, artifact, bundle = self._paths(root)
            report = root / "report.json"
            report.write_bytes(b"PRESERVE")
            original_replace = validator.os.replace
            attacked = False

            def substitute(source, target, *, src_dir_fd=None, dst_dir_fd=None):
                nonlocal attacked
                if target == report.name and not attacked:
                    attacked = True
                    validator.os.rename(
                        source,
                        ".stolen-valid-report",
                        src_dir_fd=src_dir_fd,
                        dst_dir_fd=src_dir_fd,
                    )
                    descriptor = validator.os.open(
                        source,
                        validator.os.O_WRONLY
                        | validator.os.O_CREAT
                        | validator.os.O_EXCL,
                        0o600,
                        dir_fd=src_dir_fd,
                    )
                    try:
                        validator.os.write(descriptor, b"ATTACKER")
                    finally:
                        validator.os.close(descriptor)
                return original_replace(
                    source,
                    target,
                    src_dir_fd=src_dir_fd,
                    dst_dir_fd=dst_dir_fd,
                )

            with mock.patch.object(validator.os, "replace", substitute):
                with self.assertRaisesRegex(
                    validator.INSDCGraphValidationError,
                    "changed while publishing",
                ):
                    validator.validate_insdc_graph_paths(
                        genbank,
                        artifact,
                        bundle,
                        report_path=report,
                    )
            self.assertTrue(attacked)
            self.assertEqual(report.read_bytes(), b"PRESERVE")
            self.assertEqual(
                sorted(path.name for path in root.iterdir() if path.name.startswith(".")),
                [],
            )

    @unittest.skipUnless(
        validator._HAS_DIRECTORY_DESCRIPTOR,
        "directory descriptors are unavailable",
    )
    def test_report_final_fsync_substitution_never_returns_success(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genbank, artifact, bundle = self._paths(root)
            report = root / "report.json"
            report.write_bytes(b"PRESERVE")
            original_fsync = validator.os.fsync
            attacked = False

            def substitute(descriptor):
                nonlocal attacked
                metadata = validator.os.fstat(descriptor)
                if stat.S_ISDIR(metadata.st_mode) and report.exists() and not attacked:
                    attacked = True
                    report.rename(root / ".stolen-valid-report")
                    report.write_bytes(b"ATTACKER")
                return original_fsync(descriptor)

            with mock.patch.object(validator.os, "fsync", substitute):
                with self.assertRaisesRegex(
                    validator.INSDCGraphValidationError,
                    "changed before completion",
                ):
                    validator.validate_insdc_graph_paths(
                        genbank,
                        artifact,
                        bundle,
                        report_path=report,
                    )
            self.assertTrue(attacked)
            self.assertEqual(report.read_bytes(), b"PRESERVE")
            self.assertEqual(
                sorted(path.name for path in root.iterdir() if path.name.startswith(".")),
                [],
            )

    @unittest.skipUnless(
        validator._HAS_DIRECTORY_DESCRIPTOR,
        "directory descriptors are unavailable",
    )
    def test_report_in_place_fsync_tamper_is_replayed_and_rolled_back(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genbank, artifact, bundle = self._paths(root)
            report = root / "report.json"
            report.write_bytes(b"PRESERVE")
            original_fsync = validator.os.fsync
            attacked = False

            def substitute(descriptor):
                nonlocal attacked
                metadata = validator.os.fstat(descriptor)
                if stat.S_ISDIR(metadata.st_mode) and report.exists() and not attacked:
                    attacked = True
                    report.write_bytes(b"ATTACKER")
                return original_fsync(descriptor)

            with mock.patch.object(validator.os, "fsync", substitute):
                with self.assertRaisesRegex(
                    validator.INSDCGraphValidationError,
                    "content changed during publication",
                ):
                    validator.validate_insdc_graph_paths(
                        genbank,
                        artifact,
                        bundle,
                        report_path=report,
                    )
            self.assertTrue(attacked)
            self.assertEqual(report.read_bytes(), b"PRESERVE")

    @unittest.skipUnless(
        validator._HAS_DIRECTORY_DESCRIPTOR,
        "directory descriptors are unavailable",
    )
    def test_bundle_directory_substitution_is_detected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, _, bundle = self._paths(root)
            displaced = root / "displaced"
            original_open = validator.os.open
            raced = False

            def substitute(path, flags, mode=0o777, *, dir_fd=None):
                nonlocal raced
                if Path(path) == bundle and dir_fd is None and not raced:
                    raced = True
                    bundle.rename(displaced)
                    bundle.mkdir()
                return original_open(path, flags, mode, dir_fd=dir_fd)

            with mock.patch.object(validator.os, "open", substitute):
                with self.assertRaisesRegex(
                    validator.INSDCGraphValidationError,
                    "changed while being opened",
                ):
                    validator.load_insdc_graph_bundle_directory(bundle)
            self.assertTrue(raced)

    @unittest.skipUnless(
        validator._HAS_DIRECTORY_DESCRIPTOR,
        "directory descriptors are unavailable",
    )
    def test_source_and_child_post_snapshot_substitution_is_detected(self) -> None:
        for label in ("original GenBank source", "bundle artifact target_contract.json"):
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                genbank, artifact, bundle = self._paths(root)
                original_read = validator._read_descriptor
                attacked = False

                def substitute(descriptor, opened, current_label, maximum):
                    nonlocal attacked
                    raw = original_read(descriptor, opened, current_label, maximum)
                    if current_label == label and not attacked:
                        attacked = True
                        if current_label == "original GenBank source":
                            genbank.rename(root / "displaced.gb")
                            genbank.write_bytes(b"FORGED")
                        else:
                            child = bundle / "target_contract.json"
                            child.rename(root / "displaced-child.json")
                            child.write_bytes(b"{}")
                    return raw

                with mock.patch.object(validator, "_read_descriptor", substitute):
                    with self.assertRaisesRegex(
                        validator.INSDCGraphValidationError,
                        "path changed while",
                    ):
                        validator.validate_insdc_graph_paths(
                            genbank,
                            artifact,
                            bundle,
                        )
                self.assertTrue(attacked)

    @unittest.skipUnless(
        validator._HAS_DIRECTORY_DESCRIPTOR,
        "directory descriptors are unavailable",
    )
    def test_path_runtime_loads_no_producer_modules(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genbank, artifact, bundle = self._paths(root)
            command = (
                "import sys; "
                f"sys.path.insert(0,{str(ROOT)!r}); "
                "from brainc.validator_insdc_graph import validate_insdc_graph_paths; "
                f"assert validate_insdc_graph_paths({str(genbank)!r},"
                f"{str(artifact)!r},{str(bundle)!r})['valid']; "
                "forbidden={'brainc.insdc','brainc.insdc_graph','brainc.compiler',"
                "'brainc.provider','brainc._canonical','brainc._io','brainc.v2',"
                "'brainc.v2.compiler','brainc.v2.provider','brainc.v2.target'}; "
                "loaded=sorted(forbidden.intersection(sys.modules)); "
                "assert not loaded,loaded"
            )
            subprocess.run(
                [sys.executable, "-I", "-c", command],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
            )

            public_command = (
                "import sys; "
                f"sys.path.insert(0,{str(ROOT)!r}); "
                "from brainc import (load_insdc_graph_bundle_directory,"
                "validate_insdc_graph_paths,validate_insdc_graph_report); "
                "forbidden={'brainc.insdc','brainc.insdc_graph','brainc.compiler',"
                "'brainc.provider','brainc._canonical','brainc._io','brainc.v2'}; "
                "loaded=sorted(forbidden.intersection(sys.modules)); "
                "assert not loaded,loaded"
            )
            subprocess.run(
                [sys.executable, "-I", "-c", public_command],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
            )

    def test_report_seal_and_closed_shape_are_enforced(self) -> None:
        report = validator.validate_insdc_graph_bundle(
            self.bundle.to_dict(),
            self.bundle.artifacts,
            self.source,
            genbank_source=self.raw,
        )
        changed = copy.deepcopy(report)
        changed["result"]["units"] += 1
        with self.assertRaisesRegex(
            validator.INSDCGraphValidationError,
            "seal is invalid",
        ):
            validator.validate_insdc_graph_report(changed)
        changed = copy.deepcopy(report)
        changed["unexpected"] = True
        with self.assertRaisesRegex(
            validator.INSDCGraphValidationError,
            r"unknown=\['unexpected'\]",
        ):
            validator.validate_insdc_graph_report(changed)


if __name__ == "__main__":
    unittest.main()

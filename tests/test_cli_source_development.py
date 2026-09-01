from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
import io
import json
from pathlib import Path
import tempfile
import unittest

import brainc
from brainc.cli import main
from brainc.development_bundle import CHILD_FILENAMES
from brainc.source import (
    FASTA_PROFILE,
    GFF3_PROFILE,
    GENBANK_PROFILE,
    PROFILES,
    RAW_PROFILE,
    load_source_bundle,
)
from brainc.v2 import V2Error, compile_module
from tests.test_development_bundle import _Chain


ROOT = Path(__file__).resolve().parents[1]
FASTA = ROOT / "tests/data/J02482.1.fasta"
GENBANK = ROOT / "tests/data/U49845.1.gb"
GFF3 = ROOT / "tests/data/J02482.1.gff3"


class CompilerOneZeroCliTests(unittest.TestCase):
    def test_root_api_exports_source_and_development_surfaces(self) -> None:
        source = __import__("brainc.source", fromlist=["compile_source"])
        development = __import__(
            "brainc.development_bundle",
            fromlist=["compile_development"],
        )
        self.assertIs(brainc.compile_source, source.compile_source)
        self.assertIs(
            brainc.compile_development,
            development.compile_development,
        )
        self.assertIs(brainc.SourceBundle, source.SourceBundle)
        self.assertIs(brainc.SourceError, source.SourceError)
        self.assertIs(brainc.DevelopmentBundle, development.DevelopmentBundle)
        self.assertIs(brainc.DevelopmentBundleError, development.DevelopmentBundleError)
        self.assertIs(brainc.validate_source_bundle, source.validate_source_bundle)
        self.assertEqual(brainc.PROFILES, PROFILES)
        self.assertIs(
            brainc.make_development_request,
            __import__(
                "brainc.v2.provider",
                fromlist=["make_development_request"],
            ).make_development_request,
        )

    def test_compile_source_requires_an_explicit_exact_profile(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sequence = root / "input.dna"
            sequence.write_bytes(b"ACGTRYSWKMBDHVN")
            output = root / "source"
            stdout = io.StringIO()
            with redirect_stdout(stdout):
                result = main(
                    [
                        "compile-source",
                        "--profile",
                        RAW_PROFILE,
                        "--sequence",
                        str(sequence),
                        "--record-id",
                        "raw-1",
                        "--output",
                        str(output),
                    ]
                )
            self.assertEqual(result, 0)
            self.assertEqual(load_source_bundle(output).profile, RAW_PROFILE)
            self.assertEqual(
                {path.name for path in output.iterdir()},
                {"source.json", "sequence.json"},
            )
            self.assertEqual(json.loads(stdout.getvalue())["source"], str(output / "source.json"))

            rejected = root / "rejected"
            stderr = io.StringIO()
            with redirect_stderr(stderr):
                result = main(
                    [
                        "compile-source",
                        "--profile",
                        RAW_PROFILE,
                        "--sequence",
                        str(sequence),
                        "--record-id",
                        "raw-1",
                        "--wrapper",
                        "identity",
                        "--output",
                        str(rejected),
                    ]
                )
            self.assertEqual(result, 2)
            self.assertIn("not-applicable=['wrapper']", stderr.getvalue())
            self.assertFalse(rejected.exists())

            with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as caught:
                main(
                    [
                        "compile-source",
                        "--profile",
                        "auto",
                        "--sequence",
                        str(sequence),
                        "--output",
                        str(root / "auto"),
                    ]
                )
            self.assertEqual(caught.exception.code, 2)

    def test_compile_source_routes_every_declared_profile_without_sniffing(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cases = (
                (
                    FASTA_PROFILE,
                    ["--sequence", str(FASTA), "--wrapper", "identity"],
                    {"source.json", "sequence.json"},
                ),
                (
                    GENBANK_PROFILE,
                    ["--genbank", str(GENBANK)],
                    {"source.json", "genbank.json"},
                ),
                (
                    GFF3_PROFILE,
                    [
                        "--sequence",
                        str(FASTA),
                        "--annotation",
                        str(GFF3),
                        "--sequence-profile",
                        FASTA_PROFILE,
                        "--wrapper",
                        "identity",
                    ],
                    {"source.json", "sequence.json", "annotation.json"},
                ),
            )
            for index, (profile, profile_args, filenames) in enumerate(cases):
                with self.subTest(profile=profile):
                    output = root / f"source-{index}"
                    with redirect_stdout(io.StringIO()):
                        result = main(
                            [
                                "compile-source",
                                "--profile",
                                profile,
                                *profile_args,
                                "--output",
                                str(output),
                            ]
                        )
                    self.assertEqual(result, 0)
                    self.assertEqual(load_source_bundle(output).profile, profile)
                    self.assertEqual(
                        {path.name for path in output.iterdir()},
                        filenames,
                    )

            rejected = root / "sniffed"
            raw = root / "raw.dna"
            raw.write_bytes(b"ACGT")
            stderr = io.StringIO()
            with redirect_stderr(stderr):
                result = main(
                    [
                        "compile-source",
                        "--profile",
                        FASTA_PROFILE,
                        "--sequence",
                        str(raw),
                        "--wrapper",
                        "identity",
                        "--output",
                        str(rejected),
                    ]
                )
            self.assertEqual(result, 2)
            self.assertIn("brainc: SOURCE001:", stderr.getvalue())
            self.assertNotIn("Traceback", stderr.getvalue())
            self.assertFalse(rejected.exists())

    def test_compile_development_publishes_the_exact_flat_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            chain = _Chain(root)
            output = root / "development-cli"
            stdout = io.StringIO()
            with redirect_stdout(stdout):
                result = main(
                    [
                        "compile-development",
                        str(chain.source_directory),
                        "--manifest",
                        str(chain.manifest_path),
                        "--request",
                        str(chain.request_path),
                        "--response",
                        str(chain.response_path),
                        "--policy",
                        str(chain.policy_path),
                        "--target",
                        str(chain.target_path),
                        "--output",
                        str(output),
                    ]
                )
            self.assertEqual(result, 0)
            self.assertEqual(
                {path.name for path in output.iterdir()},
                {"bundle.json", *CHILD_FILENAMES.values()},
            )
            bundle = json.loads((output / "bundle.json").read_text(encoding="utf-8"))
            self.assertEqual((bundle["format"], bundle["version"]), ("brainc.development-bundle", 1))
            self.assertEqual(json.loads(stdout.getvalue())["bundle"], str(output / "bundle.json"))

            stderr = io.StringIO()
            with redirect_stderr(stderr):
                result = main(
                    [
                        "compile-development",
                        str(chain.source_directory),
                        "--manifest",
                        str(chain.manifest_path),
                        "--request",
                        str(chain.request_path),
                        "--response",
                        str(root / "missing-response.json"),
                        "--policy",
                        str(chain.policy_path),
                        "--target",
                        str(chain.target_path),
                        "--output",
                        str(root / "failed-development"),
                    ]
                )
            self.assertEqual(result, 2)
            self.assertIn("brainc: DEVB001:", stderr.getvalue())
            self.assertNotIn("Traceback", stderr.getvalue())

    def test_detached_descriptor_is_rejected_by_low_level_api_and_cli(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            chain = _Chain(root)
            arguments = (
                chain.source_path,
                chain.manifest_path,
                chain.request_path,
                chain.response_path,
                chain.policy_path,
                chain.target_path,
            )
            with self.assertRaisesRegex(V2Error, "not a complete source"):
                compile_module(*arguments)

            output = root / "detached-module.json"
            stderr = io.StringIO()
            with redirect_stderr(stderr):
                result = main(
                    [
                        "compile-v2",
                        str(chain.source_path),
                        "--manifest",
                        str(chain.manifest_path),
                        "--request",
                        str(chain.request_path),
                        "--response",
                        str(chain.response_path),
                        "--policy",
                        str(chain.policy_path),
                        "--target",
                        str(chain.target_path),
                        "--output",
                        str(output),
                    ]
                )
            self.assertEqual(result, 2)
            self.assertIn("not a complete source", stderr.getvalue())
            self.assertNotIn("Traceback", stderr.getvalue())
            self.assertFalse(output.exists())

            generic_request = root / "detached-request.json"
            stderr = io.StringIO()
            with redirect_stderr(stderr):
                result = main(
                    [
                        "make-request-v2",
                        str(chain.source_path),
                        "--manifest",
                        str(chain.manifest_path),
                        "--output-id",
                        "t.value",
                        "--output",
                        str(generic_request),
                    ]
                )
            self.assertEqual(result, 2)
            self.assertIn("make_development_request", stderr.getvalue())
            self.assertNotIn("Traceback", stderr.getvalue())
            self.assertFalse(generic_request.exists())

            request_output = root / "closed-request.json"
            result = main(
                [
                    "make-development-request",
                    str(chain.source_directory),
                    "--manifest",
                    str(chain.manifest_path),
                    "--output-id",
                    "t.value",
                    "--output",
                    str(request_output),
                ]
            )
            self.assertEqual(result, 0)
            self.assertEqual(
                json.loads(request_output.read_text(encoding="utf-8")),
                chain.request,
            )


if __name__ == "__main__":
    unittest.main()

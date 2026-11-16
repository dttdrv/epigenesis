from __future__ import annotations

import gzip
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from brainc.development_bundle import compile_development
from brainc.source import (
    FASTA_PROFILE,
    GENBANK_PROFILE,
    GFF3_PROFILE,
    REFERENCE_FASTA_PROFILE,
    compile_source,
)
from brainc.variant import compile_variants


ROOT = Path(__file__).resolve().parents[1]
RUST_DIR = ROOT / "validators" / "rust"


def _build_rust_validator(cache_dir: Path) -> Path:
    binary = RUST_DIR / "target" / "debug" / "epigenesis-validate-rs"
    if binary.is_file():
        return binary
    subprocess.run(
        ["cargo", "build", "--offline", "--manifest-path", str(RUST_DIR / "Cargo.toml")],
        check=True,
        capture_output=True,
        text=True,
        cwd=ROOT,
    )
    return binary


@unittest.skipUnless(
    shutil.which("cargo") is not None and (RUST_DIR / "Cargo.toml").is_file(),
    "Rust cargo toolchain or validators/rust is unavailable",
)
class RustValidatorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._tmp = tempfile.TemporaryDirectory()
        cls.binary = _build_rust_validator(Path(cls._tmp.name))

    @classmethod
    def tearDownClass(cls) -> None:
        cls._tmp.cleanup()

    def test_rust_unit_tests_pass(self) -> None:
        res = subprocess.run(
            ["cargo", "test", "--offline", "--manifest-path", str(RUST_DIR / "Cargo.toml")],
            check=False,
            capture_output=True,
            text=True,
            cwd=ROOT,
        )
        self.assertEqual(res.returncode, 0, res.stdout + res.stderr)

    def test_replays_minimal_development_to_identical_hashes(self) -> None:
        example_dir = ROOT / "examples" / "minimal-development"
        expected = json.loads((example_dir / "expected.json").read_text(encoding="utf-8"))
        seq_path = example_dir / "sequence.fasta"
        interp = example_dir / "interpretation"
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            src_dir = root / "source"
            dev_dir = root / "development"
            compile_source(
                REFERENCE_FASTA_PROFILE,
                {"sequence": seq_path},
                parameters={"wrapper": "identity"},
            ).save(src_dir)
            compile_development(
                src_dir,
                interp / "manifest.json",
                interp / "request.json",
                interp / "response.json",
                interp / "lowering-policy.json",
                interp / "target-contract.json",
            ).save(dev_dir)
            res = subprocess.run(
                [
                    str(self.binary),
                    "validate-development",
                    str(src_dir),
                    str(dev_dir),
                    "--source-input",
                    f"sequence={seq_path}",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(res.returncode, 0, res.stderr + res.stdout)
            report = json.loads(res.stdout)
            self.assertTrue(report["valid"])
            self.assertEqual(report["source_artifact_sha256"], expected["source_artifact_sha256"])
            self.assertEqual(
                report["development_bundle_artifact_sha256"],
                expected["development_bundle_artifact_sha256"],
            )
            self.assertEqual(
                report["development_module_artifact_sha256"],
                expected["development_module_artifact_sha256"],
            )

    def test_replays_variant_collections_and_rejects_tampering(self) -> None:
        vectors = json.loads((ROOT / "tests" / "conformance" / "vrs-v2-vectors.json").read_text(encoding="utf-8"))
        jv = vectors["j02482_variant_collection"]
        fasta_path = ROOT / jv["source_path"]
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            # 1) FASTA profile + J02482 vectors
            src_fasta = root / "src_fasta"
            vc_fasta = root / "vc_fasta.json"
            compile_source(FASTA_PROFILE, {"sequence": fasta_path}, parameters={"wrapper": "identity"}).save(src_fasta)
            compile_variants(src_fasta, jv["inputs"]).save(vc_fasta)
            res = subprocess.run(
                [str(self.binary), "validate-variants", str(src_fasta), str(vc_fasta)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(res.returncode, 0, res.stderr + res.stdout)
            rep = json.loads(res.stdout)
            self.assertTrue(rep["valid"])
            self.assertEqual(rep["variant_artifact_sha256"], jv["expected_artifact_sha256"])
            self.assertEqual(
                rep["summary"]["vrs_allele_ids"],
                [e["vrs_allele_id"] for e in jv["expected_entries"]],
            )

            # 2) Reference gzip FASTA profile (verifies pure-Rust DEFLATE/gzip decoder)
            gz_path = root / "j02482.fasta.gz"
            gz_path.write_bytes(gzip.compress(fasta_path.read_bytes(), mtime=0))
            src_gz = root / "src_gz"
            vc_gz = root / "vc_gz.json"
            compile_source(REFERENCE_FASTA_PROFILE, {"sequence": gz_path}, parameters={"wrapper": "gzip"}).save(src_gz)
            compile_variants(src_gz, jv["inputs"], source_inputs={"sequence": gz_path}).save(vc_gz)
            res_gz = subprocess.run(
                [
                    str(self.binary),
                    "validate-variants",
                    str(src_gz),
                    str(vc_gz),
                    "--source-input",
                    f"sequence={gz_path}",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(res_gz.returncode, 0, res_gz.stderr + res_gz.stdout)
            self.assertTrue(json.loads(res_gz.stdout)["valid"])

            # 3) GenBank + GFF3 profiles
            src_gb = root / "src_gb"
            vc_gb = root / "vc_gb.json"
            compile_source(
                GENBANK_PROFILE,
                {"genbank": ROOT / "tests" / "data" / "U49845.1.gb"},
                parameters={},
            ).save(src_gb)
            compile_variants(
                src_gb,
                [{"id": "u49845-0", "record_id": "U49845.1", "start": 0, "end": 1, "ref": "G", "alt": "A"}],
            ).save(vc_gb)
            res_gb = subprocess.run(
                [str(self.binary), "validate-variants", str(src_gb), str(vc_gb)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(res_gb.returncode, 0, res_gb.stderr + res_gb.stdout)
            self.assertTrue(json.loads(res_gb.stdout)["valid"])

            src_gff = root / "src_gff"
            vc_gff = root / "vc_gff.json"
            compile_source(
                GFF3_PROFILE,
                {"sequence": fasta_path, "annotation": ROOT / "tests" / "data" / "J02482.1.gff3"},
                parameters={"sequence_profile": FASTA_PROFILE, "sequence_parameters": {"wrapper": "identity"}},
            ).save(src_gff)
            compile_variants(src_gff, jv["inputs"]).save(vc_gff)
            res_gff = subprocess.run(
                [str(self.binary), "validate-variants", str(src_gff), str(vc_gff)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(res_gff.returncode, 0, res_gff.stderr + res_gff.stdout)
            self.assertTrue(json.loads(res_gff.stdout)["valid"])

            # 4) Tampered variant collection must fail with exit code 3 & VARVAL001
            tampered = json.loads(vc_fasta.read_text(encoding="utf-8"))
            tampered["variants_ir"]["variants"][0]["alt"] = "C"
            tampered_path = root / "tampered.json"
            tampered_path.write_text(json.dumps(tampered), encoding="utf-8")
            bad = subprocess.run(
                [str(self.binary), "validate-variants", str(src_fasta), str(tampered_path)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(bad.returncode, 3, bad.stderr + bad.stdout)
            bad_rep = json.loads(bad.stdout)
            self.assertFalse(bad_rep["valid"])
            self.assertEqual(bad_rep["code"], "VARVAL001")


if __name__ == "__main__":
    unittest.main()

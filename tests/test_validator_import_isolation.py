from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]


class ValidatorImportIsolationTests(unittest.TestCase):
    def test_validator_import_loads_no_compiler_modules(self) -> None:
        command = (
            "import sys; "
            f"sys.path.insert(0, {str(ROOT)!r}); "
            "import brainc.validator_bio, brainc.validator_bio_graph; "
            "forbidden={'brainc.compiler','brainc.provider','brainc.sequence',"
            "'brainc.sequence_collection','brainc._canonical','brainc._io'}; "
            "loaded=sorted(forbidden.intersection(sys.modules)); "
            "assert not loaded, loaded"
        )
        subprocess.run(
            [sys.executable, "-I", "-c", command],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )

    def test_lazy_public_api_remains_compatible(self) -> None:
        import brainc

        self.assertEqual(brainc.SequenceCompiler.__name__, "SequenceCompiler")
        self.assertEqual(
            brainc.SequenceCollectionCompiler.__name__,
            "SequenceCollectionCompiler",
        )

    def test_lazy_public_genbank_validator_loads_no_compiler_modules(self) -> None:
        command = (
            "import sys; "
            f"sys.path.insert(0, {str(ROOT)!r}); "
            "from brainc import (INSDCValidationError, validate_genbank, "
            "validate_genbank_paths, validate_genbank_report); "
            "assert INSDCValidationError.__module__ == 'brainc.validator_insdc'; "
            "assert validate_genbank.__module__ == 'brainc.validator_insdc'; "
            "forbidden={'brainc.insdc','brainc.insdc_location','brainc.compiler',"
            "'brainc.provider','brainc.sequence','brainc.sequence_collection',"
            "'brainc.sequence_collection_v2','brainc._canonical','brainc._io'}; "
            "loaded=sorted(forbidden.intersection(sys.modules)); "
            "assert not loaded, loaded"
        )
        subprocess.run(
            [sys.executable, "-I", "-c", command],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )

    def test_insdc_graph_validator_import_loads_no_producer_modules(self) -> None:
        command = (
            "import sys; "
            f"sys.path.insert(0, {str(ROOT)!r}); "
            "import brainc.validator_insdc_graph; "
            "forbidden={'brainc.insdc','brainc.insdc_graph','brainc.insdc_location',"
            "'brainc._canonical','brainc._io','brainc.v2','brainc.v2.compiler',"
            "'brainc.v2.provider','brainc.v2.target'}; "
            "loaded=sorted(forbidden.intersection(sys.modules)); "
            "assert not loaded, loaded"
        )
        subprocess.run(
            [sys.executable, "-I", "-c", command],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )

    def test_lazy_public_genbank_compiler_api(self) -> None:
        command = (
            "import sys; "
            f"sys.path.insert(0, {str(ROOT)!r}); "
            "from brainc import (GenBankArtifact, GenBankCompiler, GenBankError, "
            "load_genbank_artifact, validate_genbank_artifact); "
            "assert GenBankCompiler.__module__ == 'brainc.insdc'; "
            "assert GenBankArtifact.__module__ == 'brainc.insdc'; "
            "assert GenBankError.__module__ == 'brainc.insdc'; "
            "assert load_genbank_artifact.__module__ == 'brainc.insdc'; "
            "assert validate_genbank_artifact.__module__ == 'brainc.insdc'"
        )
        subprocess.run(
            [sys.executable, "-I", "-c", command],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )


if __name__ == "__main__":
    unittest.main()

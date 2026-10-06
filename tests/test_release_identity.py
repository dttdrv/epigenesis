from __future__ import annotations

import json
from pathlib import Path
import unittest

import brainc
from brainc import compiler, development_bundle, external_source, source, source_scale
from brainc import validator_development, validator_external, validator_reference
from brainc.v2 import compiler as tensor_compiler


ROOT = Path(__file__).parents[1]


class ReleaseIdentityTests(unittest.TestCase):
    def test_package_release_preserves_compiler_artifact_versions(self) -> None:
        pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
        self.assertIn('version = "1.6.0"', pyproject)
        self.assertIn(
            'description = "Universal DNA translation compiler targeting Development Module IR"',
            pyproject,
        )
        self.assertEqual(brainc.__version__, "1.6.0")

        final_producers = (
            source.PRODUCER,
            source_scale.PRODUCER,
            external_source.PRODUCER,
            external_source.EXECUTABLE_PRODUCER,
            development_bundle.PRODUCER,
            validator_reference.PRODUCER,
            validator_external.CLOSURE_PRODUCER,
            validator_external.EXECUTABLE_CLOSURE_PRODUCER,
            validator_development.SOURCE_PRODUCER,
            validator_development.BUNDLE_PRODUCER,
        )
        self.assertEqual(
            {producer["version"] for producer in final_producers},
            {"1.0.0"},
        )
        self.assertEqual(compiler.COMPILER["version"], "0.4.0")
        self.assertEqual(tensor_compiler.COMPILER["version"], "0.5.0")

        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        introduction = readme.split("## Complete source-to-module example", 1)[0]
        normalized_introduction = " ".join(introduction.split())
        self.assertIn("universal DNA translation compiler", normalized_introduction)
        self.assertIn("caller-supplied interpretation", normalized_introduction)
        self.assertIn("executable external frontend", normalized_introduction)
        self.assertIn("separate executable validator", normalized_introduction)

        public_contracts = [
            ROOT / "README.md",
            ROOT / "SPEC.md",
            ROOT / "pyproject.toml",
            ROOT / "docs" / "RELEASE-1.6.md",
            ROOT / "docs" / "STANDARDS.md",
            ROOT / "docs" / "PRIOR_ART.md",
        ]
        combined = "\n".join(
            path.read_text(encoding="utf-8")
            for path in public_contracts
            if path.is_file()
        )
        self.assertIn("Epigenesis 1.6", combined)
        self.assertIn("universal dna translation compiler", combined.lower())
        self.assertIn("common typed source boundary", combined)
        self.assertIn("profile-native", combined)
        self.assertIn("caller-supplied interpretation", combined)
        self.assertIn("version-2 executable", combined.lower())
        self.assertIn("version-1", combined.lower())
        self.assertIn("digest-pinned", combined.lower())

        example = ROOT / "examples" / "minimal-development"
        for relative in (
            "sequence.fasta",
            "expected.json",
            "interpretation/manifest.json",
            "interpretation/request.json",
            "interpretation/response.json",
            "interpretation/lowering-policy.json",
            "interpretation/target-contract.json",
        ):
            self.assertTrue((example / relative).is_file(), relative)
        expected = json.loads((example / "expected.json").read_text(encoding="utf-8"))
        for identity in expected.values():
            self.assertIn(identity, readme)

        external = ROOT / "examples" / "external-fastq"
        for relative in ("frontend.py", "validator.py", "profile.py", "reads.fastq", "GRAMMAR.md"):
            self.assertTrue((external / relative).is_file(), relative)
        frontend = (external / "frontend.py").read_text(encoding="utf-8")
        validator = (external / "validator.py").read_text(encoding="utf-8")
        for semantic_field in ("header", "sequence", "quality", "quality_encoding"):
            self.assertIn(semantic_field, frontend)
            self.assertIn(semantic_field, validator)


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

from pathlib import Path
import unittest

import brainc
from brainc import compiler, development_bundle, external_source, source, source_scale
from brainc import validator_development, validator_external, validator_reference
from brainc.v2 import compiler as tensor_compiler


ROOT = Path(__file__).parents[1]


class ReleaseIdentityTests(unittest.TestCase):
    def test_public_release_is_universal_dna_compiler_one_zero(self) -> None:
        pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
        self.assertIn('version = "1.0.0"', pyproject)
        self.assertIn(
            'description = "Universal DNA translation compiler to Development Module IR"',
            pyproject,
        )
        self.assertEqual(brainc.__version__, "1.0.0")

        final_producers = (
            source.PRODUCER,
            source_scale.PRODUCER,
            external_source.PRODUCER,
            development_bundle.PRODUCER,
            validator_reference.PRODUCER,
            validator_external.CLOSURE_PRODUCER,
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
        introduction = readme.split("## Compile a DNA source", 1)[0]
        self.assertIn("one deterministic, verifiable pipeline", introduction)
        self.assertNotIn("does not", introduction.lower())
        self.assertNotIn("not mean", introduction.lower())

        public_contracts = [
            ROOT / "README.md",
            ROOT / "SPEC.md",
            ROOT / ".github" / "workflows" / "ci.yml",
            ROOT / "docs" / "RELEASE-1.0.md",
        ]
        combined = "\n".join(
            path.read_text(encoding="utf-8")
            for path in public_contracts
            if path.is_file()
        )
        self.assertIn("Epigenesis 1.0", combined)
        self.assertIn("universal DNA translation", combined)
        self.assertIn("common typed source boundary", combined)
        self.assertIn("profile-native", combined)
        self.assertIn("caller-supplied interpretation", combined)
        self.assertIn("data-only", combined)


if __name__ == "__main__":
    unittest.main()

"""Publication copies must not pass as byte-identical frozen originals."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


class PublicationChecks(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory()
        self.addCleanup(self.scratch.cleanup)
        self.root = Path(self.scratch.name)
        shutil.copy(Path(__file__).with_name("check_freeze.py"), self.root)
        self.trial = self.root / "frozen-tests" / "trial"
        self.trial.mkdir(parents=True)
        self.original = hashlib.sha256(b"original plan\n").hexdigest()
        self.published = hashlib.sha256(b"public plan\n").hexdigest()
        (self.trial / "PLAN.md").write_bytes(b"original plan\n")
        self.write(self.trial / "FREEZE.json", {"files": {"PLAN.md": self.original}})

    def write(self, path, value):
        path.write_text(json.dumps(value))

    def publication(self, **changes):
        self.write(self.root / "PUBLICATION.json", {
            "format": "epigenesis-publication-copies-v1",
            "files": {"trial/PLAN.md": {
                "original_sha256": self.original,
                "published_sha256": self.published,
                **changes,
            }},
        })
        (self.trial / "PLAN.md").write_bytes(b"public plan\n")

    def run_check(self):
        return subprocess.run([sys.executable, str(self.root / "check_freeze.py")],
                              capture_output=True, text=True)

    def test_original_bytes_pass_without_a_publication_manifest(self):
        self.assertEqual(self.run_check().returncode, 0)

    def test_declared_copy_is_reported_separately(self):
        self.publication()
        result = self.run_check()
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("0 original digests verified", result.stdout)
        self.assertIn("1 declared publication copies", result.stdout)

    def test_tampered_copy_fails(self):
        self.publication()
        (self.trial / "PLAN.md").write_bytes(b"tampered plan\n")
        self.assertNotEqual(self.run_check().returncode, 0)

    def test_wrong_original_identity_fails(self):
        self.publication(original_sha256="0" * 64)
        self.assertNotEqual(self.run_check().returncode, 0)

    def test_missing_declared_copy_fails(self):
        self.publication()
        (self.trial / "PLAN.md").unlink()
        self.assertNotEqual(self.run_check().returncode, 0)

    def test_absent_external_input_is_reported_not_verified(self):
        (self.trial / "PLAN.md").unlink()
        result = self.run_check()
        self.assertEqual(result.returncode, 0)
        self.assertIn("1 not present locally", result.stdout)

    def test_unknown_publication_format_fails(self):
        self.publication()
        path = self.root / "PUBLICATION.json"
        value = json.loads(path.read_text())
        value["format"] = "unsupported"
        self.write(path, value)
        self.assertNotEqual(self.run_check().returncode, 0)


if __name__ == "__main__":
    unittest.main()

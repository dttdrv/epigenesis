import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class ReleaseLauncherTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="epigenesis release ")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        shutil.copyfile(ROOT / "tools/release/epigenesis.py", self.root / "epigenesis.py")
        # dependencies are irrelevant to dispatch; the real model is checked separately.
        for name in ("numpy", "torch"):
            (self.root / f"{name}.py").write_text("")
        for name, path in (
            ("huber", "studio/epigenesis_studio/model/predict.py"),
            ("activity", "examples/neural-predictor/activity.py"),
            ("attention", "examples/neural-predictor/predictor.py"),
        ):
            script = self.root / path
            script.parent.mkdir(parents=True, exist_ok=True)
            script.write_text(f"import json, sys\nprint(json.dumps([{name!r}, sys.argv[1:]]))\n")

    def invoke(self, *options):
        return subprocess.run(
            [sys.executable, str(self.root / "epigenesis.py"), "predict",
             "--input", "edits with spaces.json", "--output", "scores.json", *options],
            cwd=self.root, capture_output=True, text=True,
        )

    def test_default_dispatches_to_delivered_model(self):
        result = self.invoke()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout),
                         ["huber", ["--input", "edits with spaces.json", "--output", "scores.json"]])

    def test_panel_budget_reaches_delivered_model(self):
        result = self.invoke("--budget", "3")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout),
                         ["huber", ["--input", "edits with spaces.json", "--output", "scores.json",
                                    "--budget", "3"]])

    def test_explicit_legacy_models_keep_their_arguments(self):
        for flag, model, extra in (
            ("--activity", "activity", []),
            ("--attention", "attention", ["--device", "cpu"]),
        ):
            with self.subTest(model=model):
                result = self.invoke(flag, "--families", "families.json")
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(json.loads(result.stdout),
                                 [model, ["--input", "edits with spaces.json", "--output", "scores.json",
                                          "--families", "families.json", *extra]])

    def test_incompatible_options_are_refused_before_dispatch(self):
        for options in (
            ("--families", "families.json"), ("--device", "mps"),
            ("--activity", "--budget", "3"), ("--attention", "--budget", "3"),
            ("--activity", "--attention"), ("--activity", "--device", "cpu"),
        ):
            with self.subTest(options=options):
                result = self.invoke(*options)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertEqual(result.stdout, "")

    def test_backend_failure_is_not_reported_as_success(self):
        (self.root / "studio/epigenesis_studio/model/predict.py").write_text(
            "import sys\nsys.exit(7)\n")
        self.assertEqual(self.invoke().returncode, 7)

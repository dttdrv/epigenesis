from __future__ import annotations

import importlib.util
import contextlib
import csv
import io
import json
import math
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples/neural-regulation/regulation.py"


class NeuralRegulationTests(unittest.TestCase):
    def setUp(self):
        self.assertTrue(EXAMPLE.is_file(), "the measured sequence-effect interpreter is missing")
        spec = importlib.util.spec_from_file_location("neural_regulation_example", EXAMPLE)
        self.model = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = self.model
        spec.loader.exec_module(self.model)

    def test_word_counts_preserve_reverse_complement_and_sequence_order(self):
        self.assertEqual(self.model.words("AACG", 2), {"AA": 1, "AC": 1, "CG": 1})
        self.assertEqual(self.model.words("CGTT", 2), self.model.words("AACG", 2))
        self.assertNotEqual(self.model.words("AGAC", 2), self.model.words("AACG", 2))
        with self.assertRaises(ValueError):
            self.model.words("ACNT", 2)

    def test_ridge_solution_matches_orthogonal_closed_form(self):
        # (X'WX + lambda I) B = X'WY, with independent orthogonal columns.
        result = self.model.solve_ridge([[1.0, 0.0], [0.0, 4.0]], [[2.0, -1.0], [10.0, 5.0]], 1.0)
        for observed, expected in zip(result, [[1.0, -0.5], [2.0, 1.0]]):
            for a, b in zip(observed, expected):
                self.assertAlmostEqual(a, b, places=14)
        coupled = self.model.solve_ridge([[2.0, 1.0], [1.0, 3.0]], [[5.0], [6.0]], 1.0)
        for observed, expected in zip(coupled, (14/11, 13/11)):
            self.assertAlmostEqual(observed[0], expected, places=14)
        with self.assertRaises(ValueError):
            self.model.solve_ridge([[1.0]], [[2.0]], 0.0)

    def test_transform_keeps_zero_and_has_training_defined_units(self):
        self.assertEqual(self.model.transform(0.0, 2.0), 0.0)
        expected = math.log(1 + math.sqrt(2)) / math.log(2)
        self.assertAlmostEqual(self.model.transform(2.0, 2.0), expected)
        for invalid in (-1.0, math.nan, math.inf):
            with self.assertRaises(ValueError):
                self.model.transform(invalid, 2.0)
        with self.assertRaises(ValueError):
            self.model.transform(1.0, 0.0)
        expected = (math.log(1e308) - math.log(1e-308)) / math.log(2) + 1
        self.assertAlmostEqual(self.model.transform(1e308, 1e-308), expected, places=12)

    def test_finite_model_coefficients_cannot_emit_infinite_predictions(self):
        trained = self.model.fit(self._rows()[:4], 2, 0.01, 1.0)
        trained["coefficients"] = [[1e308]*7] + [[0.0]*7 for _ in trained["coefficients"][1:]]
        with self.assertRaisesRegex(ValueError, "finite"):
            self.model.predict(trained, "A"*171, "C"*171)

    def test_families_merge_overlap_and_identical_sequences_transitively(self):
        parents = {"chr1:0-171": "A"*171, "chr1:170-341": "C"*171,
                   "chr2:500-671": "G"*171, "chr3:0-171": "AC"*85 + "A"}
        groups = self.model.family_groups(parents, {key: [value] for key, value in parents.items()})
        self.assertEqual(groups["chr1:0-171"], groups["chr1:170-341"])
        self.assertEqual(groups["chr1:0-171"], groups["chr2:500-671"])
        self.assertNotEqual(groups["chr1:0-171"][0], groups["chr3:0-171"][0])
        self.assertEqual(groups, self.model.family_groups(dict(reversed(list(parents.items()))), {key: [value] for key, value in parents.items()}))

    def test_published_data_conflicts_and_missing_parents_are_visible(self):
        rows, qc = self.model.load_data()
        self.assertEqual(qc["physical_design_records"], 10041)
        self.assertEqual(qc["exact_duplicate_records"], 8)
        self.assertEqual(len(qc["ambiguous_ids"]), 4)
        self.assertEqual(qc["missing_parent_rows"], 116)
        self.assertEqual(len(rows), 9241)
        self.assertEqual(len({row.locus for row in rows}), 585)
        self.assertTrue(any(0.0 in row.alpha for row in rows))
        self.assertFalse(set(qc["ambiguous_ids"]) & {row.id for row in rows})
        for group in {row.group for row in rows}:
            self.assertEqual(len({row.split for row in rows if row.group == group}), 1)

    def test_source_identity_mismatch_fails_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            data = Path(temporary) / "data"
            shutil.copytree(EXAMPLE.parent / "data", data)
            with (data / "design.fasta").open("ab") as stream:
                stream.write(b"A")
            with self.assertRaisesRegex(ValueError, "identity"):
                self.model.load_data(data)

    def _rows(self):
        rows = []
        reference = "A"*171
        for family, split in (("train-a", "train"), ("train-b", "train"), ("test-c", "test")):
            for base, alpha in (("C", 3.0), ("G", 0.0)):
                rows.append(self.model.Row(family+base, family, family, split, "targeted",
                    reference, reference[:-1]+base, (alpha,)*7, ((1.0,)*7,)))
        return rows

    def test_fitting_predicts_sequence_effects_without_held_out_labels(self):
        rows = self._rows()
        trained = self.model.fit([row for row in rows if row.split == "train"], 3, 1e-8, 1.0)
        for row in rows[-2:]:
            prediction = self.model.predict(trained, row.reference, row.sequence)
            expected = self.model.transform(row.alpha[0], 1.0) - self.model.transform(1.0, 1.0)
            self.assertLess(max(abs(value-expected) for value in prediction), 1e-6)
            self.assertEqual(self.model.predict(trained, row.reference, row.reference), [0.0]*7)
            self.assertEqual(self.model.predict(trained, row.sequence, row.reference), [-x for x in prediction])
        with self.assertRaises(ValueError):
            self.model.predict(trained, "A", "C")

    def test_declared_penalty_uses_mean_over_seven_outputs(self):
        row = self._rows()[0]
        trained = self.model.fit([row], 3, 1.0, 1.0)
        # one AAA is lost and one AAC gained: squared feature norm is 2.
        effect = self.model.transform(3.0, 1.0) - self.model.transform(1.0, 1.0)
        expected = 2 / (2 + 7) * effect
        for value in self.model.predict(trained, row.reference, row.sequence):
            self.assertAlmostEqual(value, expected, places=14)

    def test_family_scores_and_bootstrap_do_not_count_rows_as_independent(self):
        rows = self._rows()
        predictions = {row.id: [0.0]*7 for row in rows}
        scores = self.model.family_losses(rows, predictions, 1.0)
        self.assertEqual(len(scores), 3)
        interval = self.model.improvement_interval({key: value+1 for key,value in scores.items()}, scores)
        for endpoint in interval:
            self.assertAlmostEqual(endpoint, 1.0, places=14)

    def test_repeating_an_entire_family_cannot_increase_its_training_weight(self):
        rows = self._rows()[:4]
        original = self.model.fit(rows, 3, 0.1, self.model.training_scale(rows))
        repeated = rows + [row._replace(id=row.id+"copy") for row in rows[:2]]
        duplicate = self.model.fit(repeated, 3, 0.1, self.model.training_scale(repeated))
        for row in rows:
            for a, b in zip(self.model.predict(original, row.reference, row.sequence),
                            self.model.predict(duplicate, row.reference, row.sequence)):
                self.assertAlmostEqual(a, b, places=14)

    def test_benchmark_test_labels_change_scores_but_never_the_fitted_model(self):
        rows = self._rows()
        rows += [row._replace(id="validation-"+row.id, group="validation", split="validation") for row in rows[:2]]
        rows += [row._replace(id=kind+row.id, kind=kind) for row in rows[-4:-2] for kind in ("random", "scrambled")]
        altered = [row._replace(alpha=(100.0,)*7, wt_alphas=((50.0,)*7,)) if row.split == "test" else row for row in rows]
        qc = {"source_provenance": {"fixture": "analytic linear sequence effects"}}
        with tempfile.TemporaryDirectory() as temporary, contextlib.redirect_stdout(io.StringIO()):
            outputs = [Path(temporary) / name for name in ("original", "changed-test-labels")]
            reports = []
            for data, output in zip((rows, altered), outputs):
                with mock.patch.object(self.model, "load_data", return_value=(data, qc)), mock.patch.object(self.model, "PENALTIES", (0.01,)):
                    reports.append(self.model.benchmark(output))
            self.assertEqual((outputs[0]/"model.json").read_bytes(), (outputs[1]/"model.json").read_bytes())
            self.assertEqual((outputs[0]/"dinucleotide-model.json").read_bytes(), (outputs[1]/"dinucleotide-model.json").read_bytes())
            self.assertNotEqual(reports[0]["test_targeted_mse"], reports[1]["test_targeted_mse"])
            for report, output in zip(reports, outputs):
                self.assertEqual(report["status"], "evaluation-complete")
                self.assertFalse(report["biological_development_accepted"])
                self.assertEqual(report["predictive_acceptance"], all(interval[0] > 0 for interval in report["improvement_intervals_95"].values()))
                with (output/"test-predictions.csv").open() as stream:
                    values = [row for row in csv.DictReader(stream) if row["kind"] == "targeted"]
                # this fixture has one equally weighted test family.
                expected = sum((float(row["observed_effect"])-float(row["predicted_effect"]))**2 for row in values)/len(values)
                self.assertAlmostEqual(report["test_targeted_mse"], expected, places=12)
            with self.assertRaises(FileExistsError):
                self.model.benchmark(outputs[0])

    def test_compiled_prediction_uses_actual_source_and_rejects_substitution(self):
        trained = self.model.fit(self._rows()[:4], 3, 0.01, 1.0)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            model_path = root / "model.json"
            model_path.write_text(json.dumps(trained))
            reference, variant = root / "reference.dna", root / "variant.dna"
            reference.write_text("A"*171)
            variant.write_text("A"*170+"C")
            output = root / "compiled"
            self.model.compile_prediction(model_path, reference, variant, output)
            observed = self.model.load_prediction(output)
            self.assertEqual(observed, self.model.predict(trained, reference.read_text(), variant.read_text()))
            with self.assertRaises(FileExistsError):
                self.model.compile_prediction(model_path, reference, variant, output)
            sentinel = root/"unrelated.txt"
            sentinel.write_text("preserve me")
            (output/"validation.json").unlink()
            (output/"validation.json").symlink_to(sentinel)
            self.assertEqual(self.model.load_prediction(output), observed)
            self.assertEqual(sentinel.read_text(), "preserve me")
            (output / "pair.fasta").write_text(">reference\n"+"A"*171+"\n>variant\n"+"A"*170+"G"+"\n")
            with self.assertRaises(ValueError):
                self.model.load_prediction(output)

    def test_structurally_valid_wrong_effect_fails_sequence_replay(self):
        trained = self.model.fit(self._rows()[:4], 3, 0.01, 1.0)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            model_path, reference, variant = (root/name for name in ("model.json", "reference.dna", "variant.dna"))
            model_path.write_text(json.dumps(trained))
            reference.write_text("A"*171)
            variant.write_text("A"*170+"C")
            output = root/"compiled"
            self.model.compile_prediction(model_path, reference, variant, output)
            response = json.loads((output/"response.json").read_text())
            response.pop("artifact_sha256")
            response["outputs"][1]["storage"] = self.model.inline_storage(self.model.pack("f64", [10.0]*7))
            self.model.save(self.model.seal(response), output/"response.json")
            (output/"development").rename(output/"original-development")
            source = self.model.compile_source(self.model.FASTA_PROFILE, {"sequence": output/"pair.fasta"}, parameters={"wrapper": "identity"})
            self.model.compile_development(source, output/"manifest.json", output/"request.json", output/"response.json",
                                           output/"policy.json", output/"target.json").save(output/"development")
            validation = self.model.validate_development_paths(output/"source", output/"development", source_inputs={"sequence": output/"pair.fasta"})
            self.assertTrue(validation["valid"], validation)
            with self.assertRaisesRegex(ValueError, "disagrees"):
                self.model.load_prediction(output)


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples/neural-fate/experiment.py"


class NeuralFateExperimentTests(unittest.TestCase):
    def setUp(self) -> None:
        self.assertTrue(EXAMPLE.is_file(), "the executable neural-fate mechanism is missing")
        spec = importlib.util.spec_from_file_location("neural_fate_example", EXAMPLE)
        self.model = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = self.model
        spec.loader.exec_module(self.model)
        self.parameters = self.model.load_models()["balaskas2012"]

    def test_recipe_roundtrip_and_invalid_recipes(self) -> None:
        recipe = self.model.encode_recipe(self.parameters)
        self.assertEqual(set(recipe), set(b"ACGT"))
        self.assertEqual(self.model.decode_recipe(recipe), self.parameters)
        asymmetric = self.model.Parameters(*(index + 0.25 for index in range(1, len(self.parameters) + 1)))
        self.assertEqual(self.model.decode_recipe(self.model.encode_recipe(asymmetric)), asymmetric)
        for invalid in (recipe[:-1], recipe + b"A", b"N" + recipe[1:]):
            with self.subTest(length=len(invalid)), self.assertRaises(ValueError):
                self.model.decode_recipe(invalid)
        for replacement in (math.nan, math.inf, -1.0):
            with self.subTest(value=replacement), self.assertRaises(ValueError):
                self.model.encode_recipe(self.parameters._replace(k1=replacement))

    def test_exact_derivative_from_published_equations(self) -> None:
        # Balaskas 2012 equations 1-3, Table S2, evaluated independently at P=O=N=G=1.
        observed = self.model.derivative(self.parameters, (1.0, 1.0, 1.0), 1.0)
        for actual, expected in zip(observed, (0.0, 0.25, -1.0 / 6.0)):
            self.assertAlmostEqual(actual, expected, places=14)
        observed = self.model.derivative(self.parameters, (2.0, 3.0, 0.5), 2.0)
        for actual, expected in zip(observed, (-1090 / 641, 23 / 99, 1 / 18)):
            self.assertAlmostEqual(actual, expected, places=14)

    def test_double_knockout_matches_analytic_solution_with_nonunit_decay(self) -> None:
        parameters = self.parameters._replace(alpha=0.0, beta=0.0, k3=0.7)
        signal, initial_n = 1.3, 0.4
        trajectory = self.model.simulate(parameters, [(2.137, signal)], initial=(0.0, 0.0, initial_n), dt=0.02)
        equilibrium = parameters.gamma * signal / ((1 + signal) * parameters.k3)
        self.assertEqual(trajectory[-1][0], 2.137)
        for time, pax, olig, nkx in trajectory:
            self.assertEqual((pax, olig), (0.0, 0.0))
            expected = equilibrium + (initial_n - equilibrium) * math.exp(-parameters.k3 * time)
            self.assertLess(abs(nkx - expected), 1e-8)

    def test_zero_signal_matches_independent_exponential_limits(self) -> None:
        parameters = self.parameters._replace(k1=0.7, k2=1.3, k3=0.4)
        trajectory = self.model.simulate(parameters, [(3.137, 0.0)], initial=(0.2, 0.7, 1.1), dt=0.02)
        for time, _, olig, nkx in trajectory:
            self.assertLess(abs(olig - 0.7 * math.exp(-1.3 * time)), 1e-8)
            self.assertLess(abs(nkx - 1.1 * math.exp(-0.4 * time)), 1e-8)
        trajectory = self.model.simulate(parameters, [(3.137, 0.0)], initial=(0.2, 0.0, 0.0), dt=0.02)
        equilibrium = parameters.alpha / parameters.k1
        self.assertAlmostEqual(trajectory[-1][1], equilibrium + (0.2 - equilibrium) * math.exp(-0.7 * 3.137), places=8)

    def test_piecewise_signal_uses_exact_boundaries_and_restarts_from_state(self) -> None:
        first = self.model.simulate(self.parameters, [(1.137, 5.0)], dt=0.02)
        second = self.model.simulate(self.parameters, [(2.173, 0.5)], initial=first[-1][1:], dt=0.02)
        together = self.model.simulate(self.parameters, [(1.137, 5.0), (2.173, 0.5)], dt=0.02)
        self.assertEqual(together[len(first) - 1][0], 1.137)
        self.assertEqual(together[-1][1:], second[-1][1:])
        self.assertAlmostEqual(together[-1][0], 3.31)

    def test_invalid_protocols_fail_instead_of_clipping_or_falling_back(self) -> None:
        for schedule in ([], [(0.0, 1.0)], [(-1.0, 1.0)], [(1.0, -1.0)], [(math.inf, 1.0)], [(1.0, math.nan)]):
            with self.subTest(schedule=schedule), self.assertRaises(ValueError):
                self.model.simulate(self.parameters, schedule)
        for step in (0.0, -1.0, math.nan, math.inf):
            with self.subTest(step=step), self.assertRaises(ValueError):
                self.model.simulate(self.parameters, [(1.0, 1.0)], dt=step)
        with self.assertRaises(ValueError):
            self.model.simulate(self.parameters, [(1.0, 1.0)], initial=(-1.0, 0.0, 0.0))

    def test_finite_time_extremes_do_not_produce_division_by_zero_or_infinity(self) -> None:
        inactive = self.parameters._replace(alpha=0.0, beta=0.0, gamma=0.0)
        tiny = self.model.simulate(inactive, [(5e-324, 0.0)], dt=1e308)
        self.assertEqual(tiny[-1], (5e-324, 0.0, 0.0, 0.0))
        huge = self.model.simulate(inactive, [(1e308, 0.0)], dt=5e307)
        self.assertEqual(huge[-1], (1e308, 0.0, 0.0, 0.0))
        with self.assertRaisesRegex(ValueError, "time"):
            self.model.simulate(inactive, [(1e308, 0.0), (1e308, 0.0)], dt=1e308)
        with self.assertRaisesRegex(ValueError, "time"):
            self.model.simulate(inactive, [(1e30, 0.0), (1.0, 0.0)], dt=1e30)

    def test_step_budget_bounds_schedule_consumption_before_integration(self) -> None:
        consumed = []

        def schedule():
            for index in range(100):
                consumed.append(index)
                yield (1.0, 0.0)

        with mock.patch.object(self.model, "MAX_STEPS", 3):
            with self.assertRaisesRegex(ValueError, "budget"):
                self.model.simulate(self.parameters, schedule(), dt=1.0)
        self.assertLessEqual(len(consumed), 4)

    def test_history_dependent_fate_and_signal_withdrawal(self) -> None:
        # Balaskas 2012 Figure 6: maintenance needs less signal than induction; HIGH is >1.
        naive = self.model.simulate(self.parameters, [(40.0, 1.0)])
        primed = self.model.simulate(self.parameters, [(20.0, 5.0), (20.0, 1.0)])
        washout = self.model.simulate(self.parameters, [(20.0, 5.0), (20.0, 0.0)])
        self.assertLess(naive[-1][3], 1.0)
        self.assertGreater(primed[-1][3], 1.0)
        self.assertLess(washout[-1][3], 1e-7)
        self.assertGreater(primed[-1][3] - naive[-1][3], 1.0)

    def test_coupled_dynamics_converge_towards_independent_midpoint_reference(self) -> None:
        # separately transcribed Table S2 equations; never call the production RHS.
        def rhs(state):
            p, o, n = state
            return (3 / (1 + n**6 + o**2) - p,
                    (25 / 6) / (1 + n**5) - o,
                    (25 / 6) / (1 + o + p) - n)

        def midpoint(step):
            state = (3.0, 0.0, 0.0)
            snapshots = {}
            wanted = {round(time / step): time for time in (1, 3, 10)}
            for index in range(round(10 / step)):
                start = rhs(state)
                middle = tuple(x + step * dx / 2 for x, dx in zip(state, start))
                state = tuple(x + step * dx for x, dx in zip(state, rhs(middle)))
                if index + 1 in wanted:
                    snapshots[wanted[index + 1]] = state
            return snapshots

        reference, refined = midpoint(0.001), midpoint(0.0005)
        coarse = self.model.simulate(self.parameters, [(10.0, 5.0)], dt=0.04)
        fine = self.model.simulate(self.parameters, [(10.0, 5.0)], dt=0.02)
        for time in (1, 3, 10):
            for a, b, c, d in zip(reference[time], refined[time], coarse[round(time / 0.04)][1:], fine[round(time / 0.02)][1:]):
                self.assertLess(abs(a - b), 2e-6)
                self.assertLess(abs(d - b), 2e-6)
                self.assertLess(abs(c - d), 2e-6)

    def test_valid_compiler_artifacts_cannot_substitute_different_recipe_parameters(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            recipe = root / "recipe.dna"
            recipe.write_bytes(self.model.encode_recipe(self.parameters))
            output = root / "compiled"
            self.model.compile_recipe(recipe, output)
            response = json.loads((output / "response.json").read_text())
            response.pop("artifact_sha256")
            wrong = self.parameters._replace(beta=0.0)
            for tensor in response["outputs"]:
                if tensor["id"] == "parameters":
                    tensor["storage"] = self.model.inline_storage(self.model.pack("f64", [float(value) for value in wrong]))
            self.model.save(self.model.seal(response), output / "response.json")
            (output / "development").rename(output / "honest-development")
            self.model.compile_development(
                output / "source", output / "manifest.json", output / "request.json",
                output / "response.json", output / "policy.json", output / "target.json",
            ).save(output / "development")
            validation = self.model.validate_development_paths(
                output / "source", output / "development", source_inputs={"sequence": output / "recipe.dna"})
            self.assertTrue(validation["valid"])
            with self.assertRaisesRegex(ValueError, "disagree with the recipe"):
                self.model.load_compiled(output)

    def test_compiled_parameters_drive_gene_loss_and_exact_recipe_rescue(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            observed = []
            identities = []
            for index, parameters in enumerate((self.parameters, self.parameters._replace(alpha=0.0, beta=0.0), self.parameters)):
                recipe = root / f"recipe-{index}.dna"
                recipe.write_bytes(self.model.encode_recipe(parameters))
                output = root / f"compiled-{index}"
                self.model.compile_recipe(recipe, output)
                loaded, initial, validation = self.model.load_compiled(output)
                self.assertEqual(loaded, parameters)
                self.assertTrue(validation["valid"])
                observed.append(self.model.simulate(loaded, [(20.0, 0.5)], initial=initial)[-1][3])
                identities.append(json.loads((output / "source/source.json").read_text())["artifact_sha256"])
            self.assertGreater(observed[1], observed[0])
            self.assertEqual(observed[0], observed[2])
            self.assertNotEqual(identities[0], identities[1])
            self.assertEqual(identities[0], identities[2])
            with self.assertRaises(FileExistsError):
                self.model.compile_recipe(recipe, output)
            (output / "recipe.dna").write_bytes(b"A" * len(recipe.read_bytes()))
            with self.assertRaises(ValueError):
                self.model.load_compiled(output)

    def test_cli_reproduces_published_steady_and_oscillatory_dynamics(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "experiment"
            result = subprocess.run([sys.executable, str(EXAMPLE), "--output", str(output)], capture_output=True, text=True, cwd=ROOT)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            report = json.loads((output / "report.json").read_text())
            self.assertEqual(report["status"], "published-model-reproduced")
            self.assertTrue(all(check["passed"] for check in report["checks"]))
            self.assertFalse(report["natural_dna_inference"])
            self.assertFalse(report["independent_biological_validation"])
            self.assertGreater(report["oscillation"]["late_amplitude"], 0.1)
            self.assertLess(report["steady"]["late_amplitude"], 1e-4)
            self.assertTrue((output / "trajectories.csv").is_file())
            before = (output / "report.json").read_bytes()
            repeated = subprocess.run([sys.executable, str(EXAMPLE), "--output", str(output)], capture_output=True, text=True, cwd=ROOT)
            self.assertNotEqual(repeated.returncode, 0)
            self.assertEqual((output / "report.json").read_bytes(), before)


if __name__ == "__main__":
    unittest.main()

"""Artificial assay controls, not biological sequence interpretation or development.

The specified equilibrium is y = A*S/(K+S). Two base-4 nucleotide pairs encode
A and K (1..16); the remaining bases are neutral by construction. The oracle
uses the experimental parameters, never the candidate's sequence decoder.
"""

from __future__ import annotations

import base64
from collections.abc import Callable
from fractions import Fraction
import hashlib
import math
from pathlib import Path
import random
import struct
import tempfile
from typing import NamedTuple
import unittest

from brainc._canonical import digest
from brainc.development_bundle import compile_development
from brainc.source import FASTA_PROFILE, compile_source
from brainc.v2 import (
    inline_storage,
    make_development_request,
    pack,
    policy_artifact,
    save,
    target_artifact,
)
from brainc.v2._common import seal
from brainc.v2.target import DEV_DOMAIN, OP_UNIT_CREATE
from brainc.validator_development import load_development_bundle_directory, validate_development_paths


class Case(NamedTuple):
    name: str
    sequence: bytes
    signal: Fraction
    amplitude: int
    affinity: int
    clamp: int | None = None

    @property
    def expected(self) -> Fraction:
        return self.amplitude * self.signal / (self.affinity + self.signal)


def _encode(value: int) -> bytes:
    high, low = divmod(value - 1, 4)
    return bytes((b"ACGT"[high], b"ACGT"[low]))


def _cases(seed: int) -> list[Case]:
    rng = random.Random(seed)
    cases = []
    for family in range(4):
        amplitude = rng.randint(1, 16)
        affinity = rng.randint(1, 15)
        neutral = bytes(rng.choices(b"ACGT", k=8))
        original = _encode(amplitude) + _encode(affinity) + neutral
        mutant = _encode(amplitude) + _encode(affinity + 1) + neutral
        for signal in (Fraction(0), Fraction(1, 2), Fraction(7)):
            prefix = f"family={family},S={signal}"
            cases.extend([
                Case(prefix + ":baseline", original, signal, amplitude, affinity),
                Case(prefix + ":mutation", mutant, signal, amplitude, affinity + 1),
                Case(prefix + ":rescue", original, signal, amplitude, affinity),
                Case(prefix + ":mediation", mutant, signal, amplitude, affinity, affinity),
                Case(prefix + ":phenocopy", original, signal, amplitude, affinity + 1, affinity + 1),
            ])
            for position in (4, 7, 11):
                replacement = bytes([b"ACGT"[(b"ACGT".index(original[position]) + 1) % 4]])
                changed = original[:position] + replacement + original[position + 1:]
                cases.append(Case(prefix + f":neutral={position}", changed, signal, amplitude, affinity))
    return cases


def _synthetic_response(sequence: bytes, signal: float, clamp: int | None) -> float:
    digits = sequence[:4].translate(bytes.maketrans(b"ACGT", b"0123"))
    amplitude = int(digits[:2], 4) + 1
    affinity = int(digits[2:], 4) + 1 if clamp is None else clamp
    return 0.0 if signal == 0 else amplitude / (1 + affinity / signal)


def _compiled_response(root: Path, sequence: bytes, value: float) -> tuple[float, dict, dict]:
    """Return the validated linked initializer, module body, and source descriptor."""
    root.mkdir()
    sequence_path = root / "input.fasta"
    sequence_path.write_bytes(b">" + root.name.encode("ascii") + b"\n" + sequence + b"\n")
    source = compile_source(FASTA_PROFILE, {"sequence": sequence_path}, parameters={"wrapper": "identity"})
    source.save(root / "source")
    outputs = [
        {"id": name, "type": {"dtype": dtype, "shape": []}, "unit": None,
         "axes": [], "storage": inline_storage(pack(dtype, [number]))}
        for name, dtype, number in (("count", "u64", 1), ("response", "f64", value))
    ]
    manifest = seal({
        "format": "brainc.provider-manifest", "version": 2,
        "provider": {"name": "org.epigenesis.synthetic-assay", "version": "1"},
        "model_identity": {"kind": "content-sha256", "value": digest({"kind": "caller-supplied-assay-response"})},
        "accepts": [f"brainc.source-descriptor/v1;profile={FASTA_PROFILE}"],
        "outputs": [{key: output[key] for key in ("id", "type", "unit", "axes")} for output in outputs],
    })
    save(manifest, root / "manifest.json")
    request = make_development_request(source, root / "manifest.json", [output["id"] for output in outputs])
    save(request, root / "request.json")
    save(seal({
        "format": "brainc.prediction-response", "version": 2,
        "request_artifact_sha256": request["artifact_sha256"],
        "provider": manifest["provider"], "model_identity": manifest["model_identity"],
        "outputs": outputs,
    }), root / "response.json")
    target = target_artifact({
        "id": "org.epigenesis.synthetic-assay", "abi_major": 1,
        "opsets": [{"domain": DEV_DOMAIN, "version": 1}],
        "unit_schemas": [{"id": "assay", "fields": [{
            "id": "response", "type": {"dtype": "f64", "shape": []}, "unit": None,
            "mutability": "constant", "numeric": {"kind": "float"},
        }]}],
        "edge_schemas": [], "ports": [], "rules": [],
    })
    save(target, root / "target.json")
    policy = policy_artifact(
        "org.epigenesis.synthetic-assay",
        {"id": target["contract"]["id"], "abi_major": 1, "contract_sha256": target["contract_sha256"]},
        [{"id": output["id"], "from_output": output["id"]} for output in outputs],
        [{"id": "create.assay", "op": OP_UNIT_CREATE, "version": 1, "schema": "assay",
          "count": "count", "initializers": [{"field": "response", "tensor": "response"}]}],
    )
    save(policy, root / "policy.json")
    development = compile_development(
        source, root / "manifest.json", root / "request.json", root / "response.json",
        root / "policy.json", root / "target.json",
    )
    development.save(root / "development")
    report = validate_development_paths(root / "source", root / "development", source_inputs={"sequence": sequence_path})
    if report["valid"] is not True:
        raise AssertionError(f"artifact validation failed: {report}")
    _, artifacts = load_development_bundle_directory(root / "development")
    module = artifacts["development_module"]["module"]
    initializer = module["entrypoint"]["operations"][0]["initializers"][0]["tensor"]
    tensor, = [tensor for tensor in module["tensors"] if tensor["id"] == initializer]
    observed, = struct.unpack("<d", base64.b64decode(tensor["storage"]["data"], validate=True))
    return observed, module, source.to_dict()


def _check_prediction(case: Case, observed: float) -> None:
    # these are arithmetic tolerances, not biological uncertainty or effect sizes.
    if not math.isfinite(observed) or not math.isclose(observed, float(case.expected), rel_tol=1e-12, abs_tol=1e-12):
        raise AssertionError(
            f"prediction mismatch: {case.name}, DNA={case.sequence.decode()}, "
            f"clamp={case.clamp}, expected={case.expected}, observed={observed!r}"
        )


class NeuralMechanismAcceptanceTests(unittest.TestCase):
    def _challenge(self, candidate: Callable, *, seed: int = 20261004, reverse: bool = False) -> dict[str, float]:
        cases = _cases(seed)
        random.Random(seed ^ 0xFFFF).shuffle(cases)
        if reverse:
            cases.reverse()
        observations = {}
        with tempfile.TemporaryDirectory() as directory:
            rng = random.Random(seed + 1)
            for case in cases:
                supplied = candidate(case.sequence, float(case.signal), case.clamp)
                if not math.isfinite(supplied):
                    _check_prediction(case, supplied)
                path = Path(directory) / f"sample-{rng.getrandbits(128):032x}"
                observed, _, _ = _compiled_response(path, case.sequence, supplied)
                _check_prediction(case, observed)
                observations[case.name] = observed
        return observations

    def test_positive_control_survives_new_recipes_order_paths_and_record_labels(self) -> None:
        first = self._challenge(_synthetic_response)
        self.assertEqual(len(first), 96)
        self.assertEqual(first, self._challenge(_synthetic_response, reverse=True))
        self.assertEqual(len(self._challenge(_synthetic_response, seed=39107)), 96)

    def test_case_design_separates_mutation_neutrality_rescue_and_mediation(self) -> None:
        cases = {case.name: case for case in _cases(20261004)}
        for family in range(4):
            prefix = f"family={family},S=1/2:"
            baseline = cases[prefix + "baseline"]
            mutant = cases[prefix + "mutation"]
            self.assertNotEqual(baseline.sequence, mutant.sequence)
            self.assertGreater(baseline.expected, mutant.expected)
            self.assertEqual(cases[prefix + "rescue"].sequence, baseline.sequence)
            self.assertEqual(cases[prefix + "mediation"].sequence, mutant.sequence)
            self.assertEqual(cases[prefix + "mediation"].expected, baseline.expected)
            self.assertEqual(cases[prefix + "phenocopy"].sequence, baseline.sequence)
            self.assertEqual(cases[prefix + "phenocopy"].expected, mutant.expected)
            for position in (4, 7, 11):
                neutral = cases[prefix + f"neutral={position}"]
                self.assertNotEqual(neutral.sequence, baseline.sequence)
                self.assertEqual(neutral.expected, baseline.expected)

    def test_valid_artifacts_and_changed_source_do_not_establish_causality(self) -> None:
        cases = _cases(20261004)
        baseline, mutant = [case for case in cases if case.signal == Fraction(1, 2)][:2]
        with tempfile.TemporaryDirectory() as directory:
            first = _compiled_response(Path(directory) / "first", baseline.sequence, float(baseline.expected))
            second = _compiled_response(Path(directory) / "second", mutant.sequence, float(baseline.expected))
        self.assertNotEqual(first[2], second[2])
        self.assertNotEqual(
            first[2]["source_ir"]["records"][0]["sequence_sha256"],
            second[2]["source_ir"]["records"][0]["sequence_sha256"],
        )
        self.assertEqual(first[1], second[1])
        _check_prediction(baseline, first[0])
        with self.assertRaisesRegex(AssertionError, "prediction mismatch"):
            _check_prediction(mutant, second[0])

    def test_rejects_frozen_interpretation(self) -> None:
        with self.assertRaisesRegex(AssertionError, "prediction mismatch"):
            self._challenge(lambda sequence, signal, clamp: 1.0)

    def test_rejects_hash_as_phenotype(self) -> None:
        with self.assertRaisesRegex(AssertionError, "prediction mismatch"):
            self._challenge(lambda sequence, signal, clamp: int.from_bytes(hashlib.sha256(sequence).digest()[:4], "big") / 2**32)

    def test_rejects_correct_direction_with_wrong_effect_size(self) -> None:
        with self.assertRaisesRegex(AssertionError, "prediction mismatch"):
            self._challenge(lambda *args: 1.1 * _synthetic_response(*args))

    def test_rejects_environment_blind_candidate(self) -> None:
        with self.assertRaisesRegex(AssertionError, "prediction mismatch"):
            self._challenge(lambda sequence, signal, clamp: _synthetic_response(sequence, 1.0, clamp))

    def test_rejects_neutral_sequence_effect(self) -> None:
        neutral_region = {}
        for case in _cases(20261004):
            if case.name.endswith((":baseline", ":mutation")):
                key = case.sequence[:4]
                self.assertEqual(neutral_region.setdefault(key, case.sequence[4:]), case.sequence[4:])

        def neutral_sensitive(sequence: bytes, signal: float, clamp: int | None) -> float:
            response = _synthetic_response(sequence, signal, clamp)
            return response + (0.1 if sequence[4:] != neutral_region[sequence[:4]] else 0.0)

        with self.assertRaisesRegex(AssertionError, "prediction mismatch: .*:neutral="):
            self._challenge(neutral_sensitive)

    def test_rejects_candidate_that_ignores_mechanism_intervention(self) -> None:
        with self.assertRaisesRegex(AssertionError, "prediction mismatch: .*:(mediation|phenocopy)"):
            self._challenge(lambda sequence, signal, clamp: _synthetic_response(sequence, signal, None))

    def test_rejects_lookup_learned_from_previous_recipes(self) -> None:
        table = {(case.sequence, float(case.signal), case.clamp): float(case.expected) for case in _cases(20261004)}
        with self.assertRaisesRegex(AssertionError, "prediction mismatch"):
            self._challenge(lambda *args: table.get(args, 0.0), seed=39107)

    def test_rejects_nonfinite_predictions(self) -> None:
        case = _cases(20261004)[0]
        for value in (math.nan, math.inf, -math.inf):
            with self.subTest(value=value), self.assertRaisesRegex(AssertionError, "prediction mismatch"):
                _check_prediction(case, value)


if __name__ == "__main__":
    unittest.main()

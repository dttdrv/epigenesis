"""Independent arithmetic and input contracts, using synthetic DNA and weights."""

import copy
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import torch

import predictor


BASE = 'A'*270
ROOT = Path(__file__).resolve().parent


def rc(sequence):
    return ''.join({'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}[base]
                   for base in sequence[::-1])


def edit(position, base='G'):
    return BASE[:position]+base+BASE[position+1:]


def pair(name, reference=BASE, sequence=BASE):
    return dict(id=name, reference=reference, sequence=sequence)


def scalar(sequence, factor=1):
    value = sum({'A': 0, 'G': 2, 'C': 4, 'T': 0}[base]*(i+1)
                for i, base in enumerate(sequence))
    return float(np.float32(2**24+factor*value))


def expected(reference, sequence, factor=1):
    score = lambda s: (scalar(s, factor)+scalar(rc(s), factor))/2
    return score(sequence)-score(reference)


class SequenceSpy(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.calls = []

    def forward(self, value):
        if value.dtype != torch.float32 or value.shape[1:] != (4, 270):
            raise ValueError('wrong encoded shape or dtype')
        if not torch.equal(value.sum(1), torch.ones((len(value), 270))):
            raise ValueError('invalid one-hot encoding')
        sequences = [''.join('AGCT'[int(c)] for c in row)
                     for row in value.argmax(1).cpu().numpy()]
        self.calls.append(sequences)
        return torch.tensor([scalar(s) for s in sequences], dtype=torch.float32)


class WeightedSpy(predictor.RegulatoryAttention):
    def forward(self, value):
        factor = float(self.readout.weight[0, 0])
        sequences = [''.join('AGCT'[int(c)] for c in row)
                     for row in value.argmax(1).cpu().numpy()]
        return torch.tensor([scalar(s, factor) for s in sequences], dtype=torch.float32)


class PredictorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        torch.set_num_threads(1)

    def test_kmer_raw_overlapping_counts(self):
        words = [''.join(chars) for chars in itertools.product('ACGT', repeat=5)]
        vocabulary = sorted(word for word in words if word <= rc(word))
        self.assertEqual(predictor.vocabulary(5), vocabulary)
        coefficients = {'AAAAA': 2, 'AAAAG': 3, 'AAAGA': 5,
                        'AAGAA': 7, 'AGAAA': 11, 'GAAAA': 13}
        beta = np.array([coefficients.get(word, 0) for word in vocabulary], dtype=np.float64)
        rows = [pair(str(i), sequence=edit(i)) for i in (0, 100, 269)]
        np.testing.assert_array_equal(predictor.features(rows, 5) @ beta, [11, 29, 1])

    def test_strands_precision_and_padding_boundaries(self):
        for count in (63, 64, 65):
            with self.subTest(unique_sequences=count):
                rows = [pair(str(i), sequence=edit(i)) for i in range(count-1)]
                spy = SequenceSpy()
                actual = predictor.pair_corrections(spy, rows, 'cpu', 64, 270)
                canonical = sorted({min(s, rc(s)) for r in rows
                                    for s in (r['reference'], r['sequence'])})
                self.assertEqual(len(canonical), count)
                calls = []
                for start in range(0, count, 64):
                    batch = canonical[start:start+64]
                    batch += batch[-1:]*(64-len(batch))
                    calls.extend([batch, [rc(s) for s in batch]])
                self.assertEqual(spy.calls, calls)
                np.testing.assert_array_equal(actual, [expected(r['reference'], r['sequence']) for r in rows])
                reordered = predictor.pair_corrections(SequenceSpy(), rows[::-1], 'cpu', 64, 270)
                np.testing.assert_array_equal(actual, reordered[::-1])
        # averaging in float32 instead of float64 changes this difference to 540.
        self.assertEqual(expected(BASE, edit(0)), 541)
        a, b, c = BASE, edit(0), edit(0, 'C')
        rows = [pair('zero'), pair('ab', a, b), pair('ba', b, a),
                pair('bc', b, c), pair('ac', a, c), pair('rc', rc(a), rc(b))]
        value = predictor.pair_corrections(SequenceSpy(), rows, 'cpu', 64, 270)
        self.assertEqual(value[0], 0)
        self.assertEqual(value[1], -value[2])
        self.assertEqual(value[1], value[5])
        self.assertEqual(value[1]+value[3], value[4])

    def fixture(self, directory):
        directory.mkdir()
        manifest = copy.deepcopy(json.loads((ROOT/'data/manifest.json').read_bytes()))
        for name, factor in [('long', 1), ('short', 2)]:
            model = predictor.RegulatoryAttention(**manifest['architecture'])
            arrays = {key: np.zeros(tuple(value.shape), dtype=np.float32)
                      for key, value in model.state_dict().items()}
            arrays['readout.weight'][0, 0] = factor
            np.savez(directory/f'{name}.npz', **arrays)
            manifest['weights'][name]['sha256'] = hashlib.sha256((directory/f'{name}.npz').read_bytes()).hexdigest()
        beta = np.zeros(manifest['baseline']['parameters'], dtype=np.float64)
        np.savez(directory/'baseline.npz', coefficients=beta)
        manifest['baseline']['sha256'] = hashlib.sha256((directory/'baseline.npz').read_bytes()).hexdigest()
        manifest['mixture']['coefficients'] = [.25, .5]
        (directory/'manifest.json').write_text(json.dumps(manifest))
        return manifest

    def test_loader_mixture_input_and_weight_failures(self):
        with tempfile.TemporaryDirectory() as temporary:
            data = Path(temporary)/'data'
            manifest = self.fixture(data)
            rows = [pair('a', sequence=edit(0)), pair('b', sequence=edit(100))]
            rng = torch.get_rng_state().clone()
            with patch.object(predictor, 'RegulatoryAttention', WeightedSpy):
                result = predictor.predict(rows, model_dir=data)
            self.assertTrue(torch.equal(rng, torch.get_rng_state()))
            for row, prediction in zip(rows, result):
                long = expected(row['reference'], row['sequence'])
                short = expected(row['reference'], row['sequence'], 2)
                self.assertEqual(prediction['long_correction'], long)
                self.assertEqual(prediction['short_correction'], short)
                self.assertEqual(prediction['effect'], .25*long+.5*short)
            self.assertTrue(all(r['effect'] == 0 for r in predictor.predict(rows, model_dir=data)))
            invalid = [[], {}, [pair('')], [pair('x'), pair('x')], [pair(5)],
                       [pair('x', sequence='A'*269)], [pair('x', sequence='A'*271)],
                       [pair('x', sequence='a'*270)], [pair('x', sequence='N'+'A'*269)],
                       [pair('x', sequence='GG'+'A'*268)], [dict(pair('x'), effect=1)]]
            for rows_invalid in invalid:
                with self.subTest(input=repr(rows_invalid)[:60]), self.assertRaises(ValueError):
                    predictor.predict(rows_invalid, model_dir=data)
            original = (data/'long.npz').read_bytes()
            (data/'long.npz').write_bytes(original+b'tamper')
            with self.assertRaisesRegex(ValueError, 'identity mismatch'):
                predictor.predict(rows, model_dir=data)
            (data/'long.npz').write_bytes(original)
            with np.load(data/'long.npz', allow_pickle=False) as archive:
                arrays = {key: archive[key] for key in archive.files}
            malformed = [dict(arrays, **{'readout.weight': arrays['readout.weight'].astype(np.float64)}),
                         dict(arrays, **{'readout.weight': np.zeros((1, 1), np.float32)}),
                         dict(arrays, **{'readout.weight': np.full_like(arrays['readout.weight'], np.nan)}),
                         {k: v for k, v in arrays.items() if k != 'readout.weight'}]
            for bad in malformed:
                np.savez(data/'long.npz', **bad)
                manifest['weights']['long']['sha256'] = hashlib.sha256((data/'long.npz').read_bytes()).hexdigest()
                (data/'manifest.json').write_text(json.dumps(manifest))
                with self.assertRaisesRegex(ValueError, 'invalid model tensors'):
                    predictor.predict(rows, model_dir=data)
            (data/'short.npz').unlink()
            (data/'long.npz').write_bytes(original)
            manifest['weights']['long']['sha256'] = hashlib.sha256(original).hexdigest()
            (data/'manifest.json').write_text(json.dumps(manifest))
            with self.assertRaises(FileNotFoundError):
                predictor.predict(rows, model_dir=data)

    def test_bundled_weights_and_command_outside_checkout(self):
        rows = json.loads((ROOT/'edits.json').read_bytes())
        results = predictor.predict(rows)
        self.assertNotEqual(results[0]['effect'], 0)
        for key in ('effect', 'baseline', 'long_correction', 'short_correction'):
            self.assertEqual(results[1][key], 0)
            self.assertEqual(results[0][key], -results[2][key])
            self.assertEqual(results[0][key], results[3][key])
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            for name in ('predictor.py', 'model.py', 'edits.json'):
                shutil.copyfile(ROOT/name, directory/name)
            shutil.copytree(ROOT/'data', directory/'data')
            output = directory/'predictions.json'
            command = [sys.executable, str(directory/'predictor.py'), '--input',
                       str(directory/'edits.json'), '--output', str(output)]
            process = subprocess.run(command, cwd=directory, capture_output=True, text=True)
            self.assertEqual(process.returncode, 0, process.stderr)
            raw = output.read_bytes()
            self.assertEqual(json.loads(raw)['predictions'], results)
            process = subprocess.run(command, cwd=directory, capture_output=True, text=True)
            self.assertNotEqual(process.returncode, 0)
            self.assertEqual(output.read_bytes(), raw)
            self.assertIn('already exists', process.stderr)

    def test_command_records_the_manifest_used_for_scoring(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            for name in ('predictor.py', 'model.py', 'edits.json'):
                shutil.copyfile(ROOT/name, directory/name)
            shutil.copytree(ROOT/'data', directory/'data')
            manifest = directory/'data/manifest.json'
            original = manifest.read_bytes()
            digest = hashlib.sha256(original).hexdigest()
            rows = json.loads((directory/'edits.json').read_bytes())
            with self.assertRaisesRegex(ValueError, 'manifest identity mismatch'):
                predictor.predict(rows, model_dir=directory/'data', expected_manifest_sha256='0'*64)
            output = directory/'predictions.json'
            original_predict = predictor.predict

            def concurrent_edit(rows, device, **kwargs):
                result = original_predict(rows, device, model_dir=directory/'data', **kwargs)
                changed = json.loads(original)
                changed['mixture']['coefficients'] = [0, 0]
                manifest.write_text(json.dumps(changed))
                return result

            arguments = ['predictor.py', '--input', str(directory/'edits.json'), '--output', str(output)]
            with patch.object(predictor, 'ROOT', directory), patch.object(sys, 'argv', arguments), \
                    patch.object(predictor, 'predict', side_effect=concurrent_edit):
                predictor.main()
            result = json.loads(output.read_bytes())
            self.assertNotEqual(manifest.read_bytes(), original)
            self.assertEqual(result['manifest_sha256'], digest)
            self.assertNotEqual(result['predictions'][0]['effect'], result['predictions'][0]['baseline'])


if __name__ == '__main__':
    unittest.main()

"""Portable activity predictor checks against literal scalar sequence potentials."""
import hashlib
import itertools
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import shutil
import unittest
from unittest.mock import patch

import numpy as np

SOURCE = Path(os.environ.get('ACTIVITY_SOURCE_DIR', Path(__file__).resolve().parent))
sys.path.insert(0, str(SOURCE))
import activity


def reverse(sequence):
    complement = dict(A='T', C='G', G='C', T='A')
    return ''.join(complement[base] for base in reversed(sequence))


def reference_effect(row, arrays, k):
    words = [''.join(word) for word in itertools.product('ACGT', repeat=k)]
    canonical = sorted({min(word, reverse(word)) for word in words})
    activity_weights = dict(zip(words, arrays['activity_coefficients']))
    baseline_weights = dict(zip(canonical, arrays['baseline_coefficients']))
    def potential(sequence):
        baseline = sum(float(baseline_weights[min(sequence[i:i+k], reverse(sequence[i:i+k]))])
                       for i in range(len(sequence)-k+1))
        direct = sum(float(activity_weights[sequence[i:i+k]]) for i in range(len(sequence)-k+1))
        flipped = reverse(sequence)
        opposite = sum(float(activity_weights[flipped[i:i+k]]) for i in range(len(sequence)-k+1))
        return baseline, float(arrays['head'][0])*(direct+opposite)/2
    alternate, reference = potential(row['sequence']), potential(row['reference'])
    return tuple(a-b for a,b in zip(alternate, reference))


class ActivityTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='epigenesis-activity-test-')
        self.addCleanup(temporary.cleanup)
        self.path = Path(temporary.name)
        self.k = 5
        self.arrays = dict(activity_coefficients=((np.arange(4**self.k)*37)%101-50).astype(np.float64)/8,
            baseline_coefficients=((np.arange(512)*13)%47-23).astype(np.float64)/8,
            head=np.array([.375], dtype=np.float64))
        self.manifest = dict(format='epigenesis-activity-transfer-v1', selected_arm='linear',
            source_sha256={name:hashlib.sha256((SOURCE/name).read_bytes()).hexdigest()
                           for name in ('activity.py','predictor.py','model.py')},
            insert_bases=270, k=self.k, prediction_batch=64, ranking_budget=10,
            status='synthetic', endpoint='invented', context='invented')
        self.save()

    def save(self):
        np.savez(self.path/'activity.npz', **self.arrays)
        self.manifest['weights_sha256'] = hashlib.sha256((self.path/'activity.npz').read_bytes()).hexdigest()
        (self.path/'activity-manifest.json').write_text(json.dumps(self.manifest))
        return hashlib.sha256((self.path/'activity-manifest.json').read_bytes()).hexdigest()

    @staticmethod
    def rows(number):
        template = ('AACGATTCGCGTACCAAGGTCCTATGCTAC'*10)[:270]
        result = []
        for i in range(number):
            rotated = template[i%270:]+template[:i%270]
            position = i%270
            alternate = rotated[:position]+{'A':'C','C':'G','G':'T','T':'A'}[rotated[position]]+rotated[position+1:]
            result.append(dict(id=f'edit-{i:03}', reference=rotated, sequence=alternate))
        return result

    def assert_oracle(self, rows):
        output = activity.predict(rows, self.path)
        self.assertEqual([r['id'] for r in output], [r['id'] for r in rows])
        expected = np.array([reference_effect(row,self.arrays,self.k) for row in rows])
        np.testing.assert_allclose([r['baseline'] for r in output], expected[:,0], rtol=0, atol=1e-12)
        np.testing.assert_allclose([r['activity_correction'] for r in output], expected[:,1], rtol=0, atol=1e-12)
        np.testing.assert_allclose([r['effect'] for r in output], expected.sum(axis=1), rtol=0, atol=1e-12)
        return np.array([r['effect'] for r in output])

    def test_literal_potential_and_batch_edges(self):
        for count in (63,64,65):
            with self.subTest(count=count):
                self.assert_oracle(self.rows(count))
        rows = []
        for n,position in enumerate((0,1,4,5,264,265,268,269)):
            reference = 'A'*270
            rows.append(dict(id=f'boundary-{n}', reference=reference,
                             sequence=reference[:position]+'C'+reference[position+1:]))
        self.assert_oracle(rows)

    def test_zero_swap_reverse_and_cycle(self):
        row = self.rows(1)[0]
        reverse_row = dict(id='reverse', reference=reverse(row['reference']), sequence=reverse(row['sequence']))
        swap = dict(id='swap', reference=row['sequence'], sequence=row['reference'])
        same = dict(id='self', reference=row['reference'], sequence=row['reference'])
        values = self.assert_oracle([row,reverse_row,swap,same])
        np.testing.assert_array_equal(values, [values[0],values[0],-values[0],0.])
        a = 'A'*135+'A'+'G'*134
        b = 'A'*135+'C'+'G'*134
        c = 'A'*135+'T'+'G'*134
        cycle = [dict(id='ab',reference=a,sequence=b),dict(id='bc',reference=b,sequence=c),
                 dict(id='ac',reference=a,sequence=c)]
        scores = self.assert_oracle(cycle)
        self.assertEqual(scores[0]+scores[1],scores[2])

    def test_artifact_and_manifest_rejections(self):
        rows = self.rows(1)
        with self.assertRaisesRegex(ValueError,'manifest identity'):
            activity.predict(rows,self.path,'0'*64)
        for field,value in [('selected_arm','quadratic'),('prediction_batch',0),('k',True)]:
            old = self.manifest[field]
            self.manifest[field] = value
            self.save()
            with self.subTest(field=field), self.assertRaises(ValueError):
                activity.predict(rows,self.path)
            self.manifest[field] = old
        self.manifest['source_sha256']['activity.py'] = '0'*64
        self.save()
        with self.assertRaisesRegex(ValueError,'source identity'):
            activity.predict(rows,self.path)
        self.manifest['source_sha256']['activity.py'] = hashlib.sha256((SOURCE/'activity.py').read_bytes()).hexdigest()
        self.save()
        with (self.path/'activity.npz').open('ab') as stream:
            stream.write(b'changed')
        with self.assertRaisesRegex(ValueError,'weight identity'):
            activity.predict(rows,self.path)
        for name,value in [('head',np.array([np.nan])),('head',np.zeros(2)),
                           ('activity_coefficients',np.zeros(1024,dtype=np.float32)),('extra',np.zeros(1))]:
            original = self.arrays.get(name)
            self.arrays[name] = value
            self.save()
            with self.subTest(array=name,shape=value.shape), self.assertRaisesRegex(ValueError,'weights'):
                activity.predict(rows,self.path)
            if original is None:
                del self.arrays[name]
            else:
                self.arrays[name] = original

    @staticmethod
    def selection_fixture():
        literal = [('a0',-20,'f0'),('z0',20,'f0'),('a1',19,'f1'),('z1',18,'f1'),
                   ('b2',-17,'f2'),('a3',17,'f3')]+[(f'a{i}',20-i,f'f{i}') for i in range(4,11)]
        return [dict(id=i,effect=float(p)) for i,p,f in reversed(literal)], {i:f for i,p,f in literal}

    def assert_selection(self):
        predictions,families = self.selection_fixture()
        selected = activity.select(predictions,families,10)
        expected = ['a0','a1','a3','b2']+[f'a{i}' for i in range(4,10)]
        self.assertEqual([row['id'] for row in selected],expected)
        self.assertEqual([row['direction'] for row in selected],[-1,1,1,-1]+[1]*6)
        self.assertEqual(len({row['family'] for row in selected}),10)
        self.assertEqual([row['family'] for row in selected],[families[i] for i in expected])

    def test_family_selection_contract(self):
        self.assert_selection()
        predictions,families = self.selection_fixture()
        wrong = [None,{},dict(families,extra='other'),dict(families,a0=' '),dict(families,a0=7)]
        for mapping in wrong:
            with self.subTest(mapping=mapping), self.assertRaises(ValueError):
                activity.select(predictions,mapping,10)
        for budget in (True,0,12):
            with self.subTest(budget=budget), self.assertRaises(ValueError):
                activity.select(predictions,families,budget)
        limited = [row for row in predictions if families[row['id']] in {f'f{i}' for i in range(9)}]
        with self.assertRaises(ValueError):
            activity.select(limited,{r['id']:families[r['id']] for r in limited},10)
        with self.assertRaises(ValueError):
            activity.select(predictions+[predictions[0]],families,10)
        nonfinite = [dict(row,effect=float('inf')) if i == 0 else row for i,row in enumerate(predictions)]
        with self.assertRaises(ValueError):
            activity.select(nonfinite,families,10)

    def test_bundled_weights_and_cli_outside_checkout(self):
        deployed = self.path/'portable'
        (deployed/'data').mkdir(parents=True)
        for name in ('activity.py','predictor.py','model.py'):
            shutil.copyfile(SOURCE/name,deployed/name)
        for name in ('activity.npz','activity-manifest.json'):
            shutil.copyfile(SOURCE/'data'/name,deployed/'data'/name)
        raw = (deployed/'data/activity-manifest.json').read_bytes()
        manifest = json.loads(raw)
        with np.load(deployed/'data/activity.npz',allow_pickle=False) as archive:
            weights = {key:archive[key].copy() for key in archive.files}
        rows = self.rows(65)
        (self.path/'input.json').write_text(json.dumps(rows))
        mapping = {row['id']:f'family-{i%11}' for i,row in enumerate(rows)}
        (self.path/'families.json').write_text(json.dumps(mapping))
        output = self.path/'output.json'
        command = [sys.executable,str(deployed/'activity.py'),'--input',str(self.path/'input.json'),
                   '--output',str(output),'--families',str(self.path/'families.json')]
        env = dict(os.environ)
        env.pop('PYTHONPATH',None)
        def run(arguments):
            return subprocess.run(arguments,cwd=self.path,env=env,capture_output=True,text=True,timeout=30)
        completed = run(command)
        self.assertEqual(completed.returncode,0,completed.stderr)
        result = json.loads(output.read_bytes())
        self.assertEqual(result['manifest_sha256'],hashlib.sha256(raw).hexdigest())
        for key in ('status','endpoint','context'):
            self.assertEqual(result[key],manifest[key])
        expected = np.array([reference_effect(row,weights,manifest['k']) for row in rows])
        for column,key in enumerate(('baseline','activity_correction')):
            np.testing.assert_allclose([p[key] for p in result['predictions']],expected[:,column],rtol=0,atol=1e-12)
        np.testing.assert_allclose([p['effect'] for p in result['predictions']],expected.sum(axis=1),rtol=0,atol=1e-12)
        winners = []
        for family in sorted(set(mapping.values())):
            indices = [i for i,row in enumerate(rows) if mapping[row['id']] == family]
            winners.append(min(indices,key=lambda i:(-abs(result['predictions'][i]['effect']),rows[i]['id'])))
        chosen = sorted(winners,key=lambda i:(-abs(result['predictions'][i]['effect']),rows[i]['id']))[:manifest['ranking_budget']]
        self.assertEqual([p['id'] for p in result['selection']],[rows[i]['id'] for i in chosen])
        retained = output.read_bytes()
        repeated = run(command)
        self.assertNotEqual(repeated.returncode,0)
        self.assertIn('output already exists',repeated.stderr)
        self.assertEqual(output.read_bytes(),retained)
        (self.path/'families.json').write_text(json.dumps({row['id']:'one-family' for row in rows}))
        command[command.index('--output')+1] = str(self.path/'invalid.json')
        rejected = run(command)
        self.assertNotEqual(rejected.returncode,0)
        self.assertFalse((self.path/'invalid.json').exists())

    def test_regression_controls(self):
        original = activity.effects
        def stripped(*args):
            baseline,correction = original(*args)
            return baseline,np.zeros_like(correction)
        with patch.object(activity,'effects',stripped), self.assertRaises(AssertionError):
            self.assert_oracle(self.rows(3))
        namespace = dict(activity.__dict__)
        mutant = (SOURCE/'activity.py').read_text().replace(
            'strands = np.column_stack((beta, beta[reverse]))', 'strands = np.column_stack((beta, beta))')
        self.assertNotEqual(mutant,(SOURCE/'activity.py').read_text())
        exec(compile(mutant,str(SOURCE/'activity.py'),'exec'),namespace)
        with patch.object(activity,'effects',namespace['effects']), self.assertRaises(AssertionError):
            self.assert_oracle(self.rows(3))
        namespace = dict(activity.__dict__)
        mutant = (SOURCE/'activity.py').read_text().replace('if family not in seen:', 'if True:')
        self.assertNotEqual(mutant,(SOURCE/'activity.py').read_text())
        exec(compile(mutant,str(SOURCE/'activity.py'),'exec'),namespace)
        with patch.object(activity,'select',namespace['select']), self.assertRaises(AssertionError):
            self.assert_selection()


if __name__ == '__main__':
    unittest.main()

from __future__ import annotations

import hashlib
import importlib.util
import math
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT/"examples/neural-binding/benchmark.py"


class NeuralBindingBenchmarkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location("binding_benchmark",EXAMPLE)
        cls.benchmark = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.benchmark)

    def test_log_weight_preparation_has_fixed_gauge_floor_and_orientation(self):
        # Frozen protocol observations: exp(W - column max), then floor 0.01.
        result = self.benchmark.prepare([[1000,0],[1000-math.log(2),-1000],
                                         [1000-math.log(4),0],[0,0]],log_weights=True)
        for actual,expected in zip(result,[[1,1],[.5,.01],[.25,1],[.01,1]]):
            for a,e in zip(actual,expected):
                self.assertAlmostEqual(a,e,places=12)
        for invalid in ([[1],[2],[3]],[[0],[0],[0],[0]],[[1],[2],[3],[math.nan]],
                        [[1],[2],[3],[-1]],[[1,2],[2],[3],[4]]):
            with self.subTest(invalid=invalid),self.assertRaises(ValueError):
                self.benchmark.prepare(invalid,log_weights=False)

    def test_wt_replicates_are_all_used_and_centered_before_averaging(self):
        wt = [[[1],[.1],[.1],[.1]],[[.1],[1],[.1],[.1]]]
        mutant = [[1],[.1],[.1],[.1]]
        result = self.benchmark.score(mutant,mutant,[mutant],wt,mutant)
        self.assertEqual(result['wt_mean'],[[.55],[.55],[.1],[.1]])
        expected = [[math.log(10)/2],[-math.log(10)/2],[0],[0]]
        for row,wanted in zip(result['replicates'][0]['observed_contrast'],expected):
            self.assertAlmostEqual(row[0],wanted[0],places=14)
        self.assertAlmostEqual(result['losses']['raw_wt'],2*.45**2/4)
        self.assertAlmostEqual(result['losses']['contrast_zero'],math.log(10)**2/8)
        self.assertEqual(result['losses']['contrast_model'],result['losses']['contrast_zero'])
        self.assertEqual(result['predicted_contrast'],[[0],[0],[0],[0]])

    def test_mutation_contrast_subtracts_reference_prediction(self):
        ref = [[1],[.1],[.1],[.1]]
        mutant = [[.1],[1],[.1],[.1]]
        scored = self.benchmark.score(mutant,ref,[mutant],[ref],ref)
        self.assertEqual(scored['losses']['raw_model'],0)
        self.assertEqual(scored['losses']['contrast_model'],0)
        self.assertAlmostEqual(scored['losses']['contrast_zero'],math.log(10)**2/2)
        reverse = self.benchmark.score(ref,mutant,[ref],[mutant],mutant)
        for a,b in zip(scored['predicted_contrast'],reverse['predicted_contrast']):
            self.assertEqual(a[0],-b[0])

    def test_all_mutant_replicates_and_tied_motifs_have_equal_weight(self):
        a = [[1],[.1],[.1],[.1]]
        b = [[.1],[1],[.1],[.1]]
        nearest = self.benchmark.matrix_mean([a,b])
        self.assertEqual(nearest,[[.55],[.55],[.1],[.1]])
        result = self.benchmark.score(a,a,[a,b,b],[a],nearest)
        self.assertAlmostEqual(result['losses']['raw_model'],(.9**2*2/4)*2/3)
        self.assertAlmostEqual(result['losses']['raw_nearest'],.45**2*2/4)

    def test_nearest_baseline_excludes_all_self_copies_and_keeps_all_ties(self):
        rows = [{'name':name,'protein_alignment':alignment} for name,alignment in
                [('Target','AAA'),('target','CAA'),('Else','AAA'),
                 ('First','CCA'),('Second','CAC'),('Distant','CCC')]]
        self.assertEqual(self.benchmark.nearest_rows(rows,'TARGET','AAA','CAA'),([4,5],1))
        with self.assertRaisesRegex(ValueError,'non-self'):
            self.benchmark.nearest_rows(rows[:3],'Target','AAA','CAA')

    def test_factor_weighting_is_not_variant_or_replicate_weighting(self):
        keys = self.benchmark.LOSSES
        records = [{'allele':str(i),'factor':'many','losses':dict.fromkeys(keys,0.)} for i in range(10)]
        records.append({'allele':'last','factor':'one','losses':dict.fromkeys(keys,2.)})
        result = self.benchmark.summarize(records,['many','one'],20,73)
        self.assertEqual(result['losses'],dict.fromkeys(keys,1.))
        self.assertFalse(result['all_comparisons_pass'])
        records[-1]['losses']['raw_wt'] = math.inf
        with self.assertRaisesRegex(ValueError,'finite'):
            self.benchmark.summarize(records,['many','one'],20,73)

    def test_bootstrap_pairing_quantiles_and_joint_acceptance(self):
        # Two constant-loss factors give a degenerate independently known interval.
        losses = dict(raw_model=1.,raw_wt=2.,raw_nearest=3.,contrast_model=4.,contrast_zero=7.)
        records = [{'allele':factor,'factor':factor,'losses':losses.copy()} for factor in ['a','b']]
        result = self.benchmark.summarize(records,['a','b'],101,73)
        self.assertTrue(result['all_comparisons_pass'])
        for name,expected in zip(self.benchmark.COMPARISONS,[1.,2.,3.]):
            self.assertEqual(result['comparisons'][name]['interval95'],[expected,expected])
            self.assertEqual(result['comparisons'][name]['leave_one_factor_out_range'],[expected,expected])
        for record in records:
            record['losses']['contrast_zero'] = record['losses']['contrast_model']
        self.assertFalse(self.benchmark.summarize(records,['a','b'],101,73)['all_comparisons_pass'])
        self.assertEqual(self.benchmark.quantile([0,10,20,30],.25),7.5)
        records[0]['losses']['raw_wt'] = 0
        result = self.benchmark.summarize(records,['a','b'],101,73)
        self.assertEqual(result['comparisons']['raw_vs_wt']['interval95'],[-1,1])
        self.assertFalse(result['all_comparisons_pass'])

    def test_frozen_inputs_and_predictor_fail_before_metrics_if_changed(self):
        self.benchmark.load_inputs()
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            for name in self.benchmark.PINS:
                shutil.copyfile(EXAMPLE.parent/'benchmark-data'/name,directory/name)
            for name in self.benchmark.PINS:
                path = directory/name
                original = path.read_bytes()
                path.write_bytes(original+b' ')
                with self.subTest(name=name),self.assertRaisesRegex(ValueError,'identity'):
                    self.benchmark.load_inputs(directory)
                path.write_bytes(original)
        with mock.patch.object(self.benchmark.binding,'IMPLEMENTATION_SHA256','0'*64):
            with self.assertRaisesRegex(ValueError,'implementation'):
                self.benchmark.load_inputs()

    def test_full_frozen_panel_has_no_dropped_targets_or_reference_matrices(self):
        report = self.benchmark.evaluate()
        self.assertEqual(len(report['variants']),92)
        self.assertEqual(len(report['observations']),316)
        self.assertEqual(sum(len(v['replicates']) for v in report['variants']),222)
        for name,count,factors in [('primary',66,27),('secondary_metadata',26,8)]:
            self.assertEqual(report['cohorts'][name]['variants'],count)
            self.assertEqual(len(report['cohorts'][name]['factors']),factors)
        self.assertEqual(report['predictive_acceptance'],report['cohorts']['primary']['all_comparisons_pass'])
        self.assertNotIn('all_comparisons_pass',report['cohorts']['secondary_metadata'])
        self.assertFalse(report['neural_development_acceptance'])
        self.assertEqual(report['implementation_sha256'],hashlib.sha256(EXAMPLE.read_bytes()).hexdigest())
        samples = set()
        for variant in report['variants']:
            samples.update(variant['measured_reference_samples'])
            for replicate in variant['replicates']:
                samples.add(replicate['sample'])
                for key in self.benchmark.LOSSES:
                    residual = replicate['residuals'][key]
                    self.assertEqual([len(row) for row in residual],[8]*4)
                    self.assertAlmostEqual(replicate['losses'][key],
                                           sum(v*v for row in residual for v in row)/32,places=14)
        self.assertEqual(samples,{r['sample'] for r in report['observations']})


if __name__ == '__main__':
    unittest.main()

import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]/'examples/neural-predictor'


class PredictorBundleTests(unittest.TestCase):
    def test_activity_assets_and_evidence(self):
        raw = (ROOT/'data/activity-manifest.json').read_bytes()
        manifest = json.loads(raw)
        for name, digest in dict(manifest['source_sha256'], **{'data/activity.npz':manifest['weights_sha256']}).items():
            with self.subTest(path=name):
                self.assertEqual(hashlib.sha256((ROOT/name).read_bytes()).hexdigest(),digest)
        evidence = json.loads((ROOT/'results/activity.json').read_bytes())
        export = json.loads((ROOT/'results/activity-export.json').read_bytes())
        self.assertEqual(export['model_manifest_sha256'],hashlib.sha256(raw).hexdigest())
        self.assertEqual(evidence['identities']['candidate_sha256'],manifest['source_candidate_sha256'])
        confirmation = evidence['confirmation']
        original = (json.dumps(confirmation,indent=2,allow_nan=False)+'\n').encode()
        self.assertEqual(hashlib.sha256(original).hexdigest(),evidence['identities']['confirmation_result_sha256'])
        self.assertEqual(confirmation['input_identity'],manifest['confirmation_freeze_sha256'])
        self.assertEqual(set(confirmation['comparisons']),
                         {'mse_zero','mse_dinucleotide','mse_kmer','ranking_zero','ranking_dinucleotide','ranking_kmer','ranking_random'})
        self.assertEqual(confirmation['accepted'],all(v['interval'][0] > 0 for v in confirmation['comparisons'].values()))
        for record in export['populations'].values():
            self.assertEqual(record['rows'],record['exact_prediction_values'])
            self.assertTrue(record['selected_ids_exact'])

    def test_frozen_inference_assets(self):
        manifest = json.loads((ROOT/'data/manifest.json').read_bytes())
        paths = {'model.py': manifest['model_sha256'],
                 'data/baseline.npz': manifest['baseline']['sha256']}
        paths.update({f'data/{name}.npz': spec['sha256']
                      for name, spec in manifest['weights'].items()})
        for path, digest in paths.items():
            with self.subTest(path=path):
                self.assertEqual(hashlib.sha256((ROOT/path).read_bytes()).hexdigest(), digest)
        for path in ('requirements.txt', 'LICENSE.Salomon', 'NOTICE.md', 'edits.json'):
            self.assertTrue((ROOT/path).is_file(), path)

    def test_evidence_links_to_packaged_predictor(self):
        read = lambda name: json.loads((ROOT/'results'/name).read_bytes())
        validation = read('validation.json')
        review = read('mixture-review.json')
        result = validation['models']['mixture']
        self.assertEqual(result['source_state_sha256'], review['experiment_state_sha256'])
        self.assertEqual(result['validation_mse'], review['results']['interlace']['validation_mse'])
        self.assertEqual(result['ranking_utility'], review['results']['interlace']['ranking_utility'])
        manifest = json.loads((ROOT/'data/manifest.json').read_bytes())
        self.assertEqual(manifest['mixture'], review['mixture'])
        for path, digest in read('inference-export.json')['files'].items():
            self.assertEqual(hashlib.sha256((ROOT/path).read_bytes()).hexdigest(), digest, path)
        uncertainty = read('uncertainty-review.json')
        self.assertEqual(hashlib.sha256((ROOT/'results/uncertainty.json').read_bytes()).hexdigest(),
                         uncertainty['report_sha256'])
        self.assertFalse(result['confirmation_evaluated'])
        self.assertIsNone(result['nomination'])


if __name__ == '__main__':
    unittest.main()

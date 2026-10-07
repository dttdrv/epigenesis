import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]/'examples/neural-predictor'


class PredictorBundleTests(unittest.TestCase):
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

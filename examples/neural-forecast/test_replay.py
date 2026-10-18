"""Exercise relocation, archive binding and a deliberately broken calculation."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
import hashlib

ROOT=Path(__file__).resolve().parent


class Replay(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='epigenesis-forecast-replay-')
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)/'archive'
        shutil.copytree(ROOT,self.root,ignore=shutil.ignore_patterns('__pycache__'))

    def run_replay(self):
        env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
        return subprocess.run([sys.executable,str(self.root/'replay.py')],text=True,capture_output=True,env=env)

    def test_relocated_archive_reproduces_failed_decision(self):
        result=self.run_replay()
        self.assertEqual(result.returncode,0,result.stderr)
        result=json.loads(result.stdout)
        self.assertEqual(result['decision'],'not_passed')
        self.assertEqual(result['strata'],15)
        self.assertEqual(result['verified_prediction_files'],30)
        self.assertLess(result['maximum_absolute_numeric_difference'],2e-12)

    def test_altered_counts_are_rejected(self):
        path=self.root/'result/COUNTS.npz'
        with path.open('r+b') as stream:stream.write(b'wrong')
        result=self.run_replay()
        self.assertNotEqual(result.returncode,0)
        self.assertIn('archive identity differs: result/COUNTS.npz',result.stderr)

    def test_wrong_feature_mean_is_detected_after_resealing_code(self):
        path=self.root/'replay.py'
        original=path.read_text()
        self.assertEqual(original.count('mean(prefix[:3])'),1)
        path.write_text(original.replace('mean(prefix[:3])','mean(prefix[:2])'))
        manifest=json.loads((self.root/'MANIFEST.json').read_text())
        manifest['files']['replay.py']=hashlib.sha256(path.read_bytes()).hexdigest()
        (self.root/'MANIFEST.json').write_text(json.dumps(manifest))
        result=self.run_replay()
        self.assertNotEqual(result.returncode,0)
        self.assertIn('AssertionError',result.stderr)
        self.assertIn('center',result.stderr)


if __name__=='__main__':unittest.main()

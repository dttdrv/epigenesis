"""add S2 and Ceq6 = z(N1)+z(S2) to the CNN predictions; seal-ready PREDICTIONS.npz."""
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parent.parent; D6 = ROOT.parent / 'hnp-loso6-20261019/derived'
raw, p6 = np.load(ROOT / 'derived/PREDICTIONS_RAW.npz'), np.load(D6 / 'PREDICTIONS.npz'); out = {}
z = lambda a: a / a.std()
for s in ['lee', 'mcafee', 'matoba', 'gse225815', 'gse221384', 'weiss']:
    assert list(raw[f'{s}/ids']) == list(p6[f'{s}/ids'])
    out[f'{s}/ids'], out[f'{s}/family'] = raw[f'{s}/ids'], raw[f'{s}/family']
    out[f'{s}/N1'], out[f'{s}/S2'] = raw[f'{s}/N1'], p6[f'{s}/S2']
    out[f'{s}/Ceq6'] = z(out[f'{s}/N1']) + z(out[f'{s}/S2'])
np.savez(ROOT / 'derived/PREDICTIONS.npz', **out)

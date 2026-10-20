"""frozen Salomon-trained activity-transfer model applied unchanged to the six studies' held-out variants. nothing is fitted here.
A4 = baseline + correction (the model); Abase = k5 baseline part; Aact = activity-correction part; Ceq5 = z(A4)+z(S2), S2 = siraj-trained head."""
import json, sys
from pathlib import Path
import numpy as np
BUILD = Path(__file__).resolve().parent.parent.parent
EX = BUILD.parent.parent / 'examples/neural-predictor'
sys.path.insert(0, str(EX)); import predictor
D6 = BUILD / 'hnp-loso6-20261019/derived'
W = np.load(EX / 'data/activity.npz'); K = 5
words = [''.join(b) for b in __import__('itertools').product('ACGT', repeat=K)]
col = {w: i for i, w in enumerate(words)}
rc = {w: predictor.reverse_complement(w) for w in words}
rev = np.array([col[rc[w]] for w in words])
beta = W['activity_coefficients']; strands = np.column_stack((beta, beta[rev]))
def counts(s):
    v = np.zeros(len(words))
    for i in range(len(s) - K + 1): v[col[s[i:i + K]]] += 1
    return v
P6 = np.load(D6 / 'PREDICTIONS.npz'); out = {}
z = lambda a: a / a.std()
for s in ['lee', 'mcafee', 'matoba', 'gse225815', 'gse221384', 'weiss']:
    rows = json.load(open(D6 / f'inputs_{s}.json')); idx = {r['id']: i for i, r in enumerate(rows)}
    sel = [rows[idx[i]] for i in P6[f'{s}/ids']]
    delta = np.array([counts(r['sequence']) - counts(r['reference']) for r in sel])
    act = (delta @ strands).mean(axis=1) * W['head'][0]
    base = predictor.features(sel, K) @ W['baseline_coefficients']
    assert np.isfinite(act).all() and np.isfinite(base).all()
    S2 = P6[f'{s}/S2']
    for k, v in {'ids': P6[f'{s}/ids'], 'family': P6[f'{s}/family'], 'A4': base + act, 'Aact': act, 'Abase': base, 'S2': S2, 'Ceq5': z(base + act) + z(S2)}.items():
        out[f'{s}/{k}'] = v
    print(s, len(sel), 'std base %.4g act %.4g' % (base.std(), act.std()), flush=True)
np.savez(Path(__file__).resolve().parent.parent / 'derived/PREDICTIONS.npz', **out)

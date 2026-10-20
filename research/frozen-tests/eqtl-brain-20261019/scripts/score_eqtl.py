"""score sealed eQTL-direction predictions once. --synthetic [--null] replaces the labels to qualify the scorer."""
import hashlib, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parent.parent; DER = ROOT / 'derived'
SYNTH, NULL = '--synthetic' in sys.argv, '--null' in sys.argv
frz = json.load(open(ROOT / 'FREEZE.json'))['files']
for f in ('PLAN.md', 'derived/PREDICTIONS.npz'):
    assert hashlib.sha256((ROOT / f).read_bytes()).hexdigest() == frz[f], 'sealed input changed: ' + f
P = np.load(DER / 'PREDICTIONS.npz'); MODELS = ('M4', 'M5', 'B1', 'B2', 'B3')
tab = {(r['dataset'], r['chrom'], r['pos'], r['alt'], r['gene']): r['beta'] for r in json.load(open(DER / 'eqtl_table.json')) if r['role'] == 'test'}
y = np.array([1.0 if tab[(d, c, int(p), a, g)] > 0 else 0.0 for d, c, p, a, g in zip(P['dataset'], P['chrom'], P['pos'], P['alt'], P['gene'])])
rng = np.random.Generator(np.random.PCG64(20261019)); srng = np.random.default_rng(1)
if SYNTH:
    z = P['M4'] - 0.5; z = z / z.std()
    y = (srng.random(len(y)) < (0.5 if NULL else 1 / (1 + np.exp(-2.0 * z)))).astype(float)
def acc(s, yy):  # ties (score exactly 0.5) get half credit
    return float(np.mean(np.where(s == 0.5, 0.5, ((s > 0.5) == (yy == 1)).astype(float))))
def auc(s, yy):
    from scipy.stats import rankdata
    r = rankdata(s); n1, n0 = yy.sum(), (1 - yy).sum(); return float((r[yy == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))
res = {'synthetic': SYNTH, 'null': NULL, 'sealed_predictions_sha256': frz['derived/PREDICTIONS.npz']}
for ds in sorted(set(P['dataset'])):
    m = P['dataset'] == ds; genes = P['gene'][m]; ug, inv = np.unique(genes, return_inverse=True)
    idx = [np.where(inv == k)[0] for k in range(len(ug))]
    S = {k: P[k][m] for k in MODELS}; yy = y[m]
    draws = [np.concatenate([idx[j] for j in rng.integers(0, len(ug), len(ug))]) for _ in range(2000)]
    out = {'n': int(m.sum()), 'genes': int(len(ug)), 'positive_rate': float(yy.mean())}
    for k in MODELS:
        pt = acc(S[k], yy); bs = np.array([acc(S[k][d], yy[d]) for d in draws])
        out['acc_' + k] = {'point': pt, 'q05_50_95': [float(x) for x in np.quantile(bs, [.05, .5, .95])], 'auc': auc(S[k], yy) if len(set(S[k])) > 1 else None}
    base = max(('B1', 'B2', 'B3'), key=lambda k: out['acc_' + k]['point']); out['best_baseline'] = base
    for k in ('M4', 'M5'):
        bs = np.array([acc(S[k][d], yy[d]) - acc(S[base][d], yy[d]) for d in draws])
        out[f'gain_{k}_over_{base}'] = {'point': out['acc_' + k]['point'] - out['acc_' + base]['point'], 'q05_50_95': [float(x) for x in np.quantile(bs, [.05, .5, .95])]}
    res[ds] = out
w = res['walker_fetal_neocortex']; base = w['best_baseline']
g1 = w['acc_M4']['q05_50_95'][0] > 0.5; g2 = w[f'gain_M4_over_{base}']['q05_50_95'][0] > 0; g3 = w[f'gain_M5_over_{base}']['q05_50_95'][0] > 0
res['gates'] = {'G1_M4_acc_bound_gt_0.5': g1, 'G2_M4_gain_over_best_baseline_bound_gt_0': g2, 'G3_M5_gain_over_best_baseline_bound_gt_0': g3}
ay = [res[k]['acc_M4']['point'] > 0.5 for k in res if k.startswith('aygun')]
res['label_fetal_brain_eqtl_direction'] = bool(g1 and g2 and all(ay))
name = 'RESULT' + ('_SYNTHETIC_NULL' if NULL else '_SYNTHETIC' if SYNTH else '') + '.json'
json.dump(res, open(DER / name, 'w'), indent=1)
for ds in sorted(set(P['dataset'])):
    o = res[ds]; print(ds, o['n'], {k: round(o['acc_' + k]['point'], 3) for k in MODELS}, 'auc M4', o['acc_M4']['auc'])
print(res['gates'], res['label_fetal_brain_eqtl_direction'])

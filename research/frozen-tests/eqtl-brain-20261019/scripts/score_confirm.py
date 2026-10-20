"""score sealed Jerber-neuron eQTL-direction predictions once. --synthetic [--null] replaces the labels to qualify the scorer."""
import hashlib, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parent.parent; DER = ROOT / 'derived'
SYNTH, NULL = '--synthetic' in sys.argv, '--null' in sys.argv
frz = json.load(open(ROOT / 'FREEZE_CONFIRM.json'))['files']
for f in ('PLAN_CONFIRM.md', 'derived/CONFIRM_PREDICTIONS.npz'):
    assert hashlib.sha256((ROOT / f).read_bytes()).hexdigest() == frz[f], 'sealed input changed: ' + f
P = np.load(DER / 'CONFIRM_PREDICTIONS.npz'); MODELS = ('M4', 'M5', 'B1', 'B2', 'B3')
tab = {(r['dataset'], r['chrom'], r['pos'], r['alt'], r['gene']): r['beta'] for r in json.load(open(DER / 'confirm_table.json'))}
y = np.array([1.0 if tab[(d, c, int(p), a, g)] > 0 else 0.0 for d, c, p, a, g in zip(P['dataset'], P['chrom'], P['pos'], P['alt'], P['gene'])])
rng = np.random.Generator(np.random.PCG64(20261019)); srng = np.random.default_rng(1)
if SYNTH:
    z = P['M5'] - 0.5; z = z / z.std()
    y = (srng.random(len(y)) < (0.5 if NULL else 1 / (1 + np.exp(-2.0 * z)))).astype(float)
def acc(s, yy):  # ties (score exactly 0.5) get half credit
    return float(np.mean(np.where(s == 0.5, 0.5, ((s > 0.5) == (yy == 1)).astype(float))))
def auc(s, yy):
    from scipy.stats import rankdata
    r = rankdata(s); n1, n0 = yy.sum(), (1 - yy).sum(); return float((r[yy == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))
q = lambda a: [float(x) for x in np.quantile(a, [.05, .5, .95])]
genes = P['gene']; ug, inv = np.unique(genes, return_inverse=True); idx = [np.where(inv == k)[0] for k in range(len(ug))]
S = {k: P[k] for k in MODELS}
draws = [np.concatenate([idx[j] for j in rng.integers(0, len(ug), len(ug))]) for _ in range(2000)]
res = {'synthetic': SYNTH, 'null': NULL, 'sealed_predictions_sha256': frz['derived/CONFIRM_PREDICTIONS.npz'], 'n': int(len(y)), 'genes': int(len(ug)), 'positive_rate': float(y.mean())}
for k in MODELS:
    bs = np.array([acc(S[k][d], y[d]) for d in draws])
    res['acc_' + k] = {'point': acc(S[k], y), 'q05_50_95': q(bs), 'auc': auc(S[k], y) if len(set(S[k])) > 1 else None}
base = max(('B1', 'B2', 'B3'), key=lambda k: res['acc_' + k]['point']); res['best_baseline'] = base
for k in ('M4', 'M5'):
    bs = np.array([acc(S[k][d], y[d]) - acc(S[base][d], y[d]) for d in draws])
    res[f'gain_{k}_over_{base}'] = {'point': res['acc_' + k]['point'] - res['acc_' + base]['point'], 'q05_50_95': q(bs)}
res['gates'] = {'C1_M5_acc_bound_gt_0.5': res['acc_M5']['q05_50_95'][0] > 0.5, 'C2_M4_acc_bound_gt_0.5': res['acc_M4']['q05_50_95'][0] > 0.5,
                'C3_M5_gain_over_best_baseline_bound_gt_0': res[f'gain_M5_over_{base}']['q05_50_95'][0] > 0, 'C4_M4_gain_over_best_baseline_bound_gt_0': res[f'gain_M4_over_{base}']['q05_50_95'][0] > 0}
res['label_eqtl_direction_confirmed_in_differentiating_neurons'] = bool(res['gates']['C1_M5_acc_bound_gt_0.5'] and res['gates']['C3_M5_gain_over_best_baseline_bound_gt_0'])
name = 'CONFIRM_RESULT' + ('_SYNTHETIC_NULL' if NULL else '_SYNTHETIC' if SYNTH else '') + '.json'
json.dump(res, open(DER / name, 'w'), indent=1)
print(res['n'], res['genes'], {k: [round(res['acc_' + k]['point'], 3), res['acc_' + k]['auc'] and round(res['acc_' + k]['auc'], 3)] for k in MODELS}, 'best baseline', base)
print({k: res[k] for k in res if k.startswith('gain_')}); print(res['gates'], res['label_eqtl_direction_confirmed_in_differentiating_neurons'])

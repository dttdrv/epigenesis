"""score sealed credible-set ranking predictions once. --synthetic [--null] replaces the labels to qualify the scorer."""
import hashlib, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parent.parent; DER = ROOT / 'derived'
SYNTH, NULL = '--synthetic' in sys.argv, '--null' in sys.argv
frz = json.load(open(ROOT / 'FREEZE2.json'))['files']
for f in ('PLAN2.md', 'derived/PREDICTIONS2.npz'):
    assert hashlib.sha256((ROOT / f).read_bytes()).hexdigest() == frz[f], 'sealed input changed: ' + f
P = np.load(DER / 'PREDICTIONS2.npz'); T = json.load(open(DER / 'cs_table.json'))
MODELS = ('B1', 'B2', 'M4', 'M6', 'M7')
rng = np.random.Generator(np.random.PCG64(20261019)); srng = np.random.default_rng(1)

def labels(name):  # index of the top-PIP member of every set, in the order run_cs.py stored them
    ptr = P[name + '/ptr']; chrom, pos, alt = P[name + '/chrom'], P[name + '/pos'], P[name + '/alt']
    sets = {tuple(sorted((m['chrom'], m['pos'], m['alt']) for m in s['members'])): max(s['members'], key=lambda m: m['pip']) for s in T[name if name == 'walker' else 'secondary']}
    y = []
    for a, b in zip(ptr[:-1], ptr[1:]):
        key = tuple(sorted((str(chrom[i]), int(pos[i]), str(alt[i])) for i in range(a, b))); t = sets[key]
        y.append(next(i - a for i in range(a, b) if (str(chrom[i]), int(pos[i]), str(alt[i])) == (t['chrom'], t['pos'], t['alt'])))
    return np.array(y)
def sid(ptr): return np.repeat(np.arange(len(ptr) - 1), np.diff(ptr))
def per_set(s, ptr, y):  # top-1 credit (ties share) and pairwise accuracy against the other members, per set
    n = len(ptr) - 1; si = sid(ptr); mx = np.maximum.reduceat(s, ptr[:-1]); win = s == mx[si]; k = np.add.reduceat(win.astype(float), ptr[:-1])
    lab = np.zeros(len(s), bool); lab[ptr[:-1] + y] = True; t1 = np.add.reduceat((win & lab).astype(float), ptr[:-1]) / k
    sl = s[ptr[:-1] + y]; gt = np.add.reduceat((s < sl[si]).astype(float), ptr[:-1]); eq = np.add.reduceat((s == sl[si]).astype(float), ptr[:-1]) - 1
    return t1, (gt + 0.5 * eq) / (np.diff(ptr) - 1)
q = lambda a: [float(x) for x in np.quantile(a, [.05, .5, .95])]
res = {'synthetic': SYNTH, 'null': NULL, 'sealed_predictions_sha256': frz['derived/PREDICTIONS2.npz']}
cache = {}
for name in ('walker', 'secondary'):
    ptr = P[name + '/ptr']; y = labels(name); S = {k: P[f'{name}/{k}'].astype(float) for k in MODELS}
    if SYNTH:
        z = S['M6'] - np.repeat(np.add.reduceat(S['M6'], ptr[:-1]) / np.diff(ptr), np.diff(ptr)); z = z / z.std(); g = srng.gumbel(size=len(z))
        sc = (0.0 if NULL else 2.0) * z + g; y = np.array([int(np.argmax(sc[a:b])) for a, b in zip(ptr[:-1], ptr[1:])])
    genes = P[name + '/gene'][ptr[:-1]]; ug, inv = np.unique(genes, return_inverse=True); idx = [np.where(inv == k)[0] for k in range(len(ug))]
    draws = [np.concatenate([idx[j] for j in rng.integers(0, len(ug), len(ug))]) for _ in range(2000)]
    base = {k: per_set(S[k], ptr, y) for k in MODELS}; rand = float((1 / np.diff(ptr)).mean())
    out = {'sets': int(len(ptr) - 1), 'members': int(len(y) and ptr[-1]), 'genes': int(len(ug)), 'random_top1': rand}
    for k in MODELS:
        out['top1_' + k] = {'point': float(base[k][0].mean()), 'q05_50_95': q([base[k][0][d].mean() for d in draws]), 'pairwise_auc': float(base[k][1].mean())}
    for a, b in (('M6', 'B1'), ('M7', 'B1'), ('M6', 'M4'), ('M6', 'B2')):
        out[f'gain_top1_{a}_minus_{b}'] = {'point': float(base[a][0].mean() - base[b][0].mean()), 'q05_50_95': q([base[a][0][d].mean() - base[b][0][d].mean() for d in draws])}
    res[name] = out; print(name, out['sets'], {k: round(out['top1_' + k]['point'], 3) for k in MODELS}, 'random', round(rand, 3))
w, s = res['walker'], res['secondary']
res['gates'] = {'G1_M6_gain_over_nearest_TSS_bound_gt_0': w['gain_top1_M6_minus_B1']['q05_50_95'][0] > 0, 'G2_M7_gain_over_nearest_TSS_bound_gt_0': w['gain_top1_M7_minus_B1']['q05_50_95'][0] > 0,
                'G3_M6_gain_over_M4_bound_gt_0': w['gain_top1_M6_minus_M4']['q05_50_95'][0] > 0, 'secondary_M6_gain_over_nearest_TSS_point_gt_0': s['gain_top1_M6_minus_B1']['point'] > 0}
LAB = 'label_fetal_accessibility_plus_sequence_model_prioritises_causal_variants'
res[LAB] = bool(res['gates']['G1_M6_gain_over_nearest_TSS_bound_gt_0'] and res['gates']['secondary_M6_gain_over_nearest_TSS_point_gt_0'])
json.dump(res, open(DER / ('RESULT2' + ('_SYNTHETIC_NULL' if NULL else '_SYNTHETIC' if SYNTH else '') + '.json'), 'w'), indent=1)
print({k: res[k] for k in ('gates', LAB)})
for n_ in ('walker', 'secondary'): print(n_, {k: res[n_][k] for k in res[n_] if k.startswith('gain')})

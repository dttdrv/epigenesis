"""score sealed sensory-neuron / microglia eQTL-direction predictions once. --synthetic [--null] replaces the labels to qualify the scorer."""
import hashlib, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parent.parent; DER = ROOT / 'derived'
SYNTH, NULL = '--synthetic' in sys.argv, '--null' in sys.argv
frz = json.load(open(ROOT / 'FREEZE_REP.json'))['files']
for f in ('PLAN_REP.md', 'derived/REP_PREDICTIONS.npz'):
    assert hashlib.sha256((ROOT / f).read_bytes()).hexdigest() == frz[f], 'sealed input changed: ' + f
MODELS = ('M4', 'M5', 'B1', 'B2', 'B3')
def load(pred, table, roles=None):
    P = np.load(DER / pred); tab = {(r['dataset'], r['chrom'], r['pos'], r['alt'], r['gene']): r['beta'] for r in json.load(open(DER / table)) if roles is None or r['role'] in roles}
    y = np.array([1.0 if tab[(d, c, int(p), a, g)] > 0 else 0.0 for d, c, p, a, g in zip(P['dataset'], P['chrom'], P['pos'], P['alt'], P['gene'])])
    return P, y
P, y = load('REP_PREDICTIONS.npz', 'rep_table.json')
rng = np.random.Generator(np.random.PCG64(20261019)); srng = np.random.default_rng(1)
if SYNTH:
    z = P['M5'] - 0.5; z = z / z.std(); y = (srng.random(len(y)) < (0.5 if NULL else 1 / (1 + np.exp(-2.0 * z)))).astype(float)
def acc(s, yy): return float(np.mean(np.where(s == 0.5, 0.5, ((s > 0.5) == (yy == 1)).astype(float))))
def auc(s, yy):
    from scipy.stats import rankdata
    r = rankdata(s); n1, n0 = yy.sum(), (1 - yy).sum(); return float((r[yy == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))
q = lambda a: [float(x) for x in np.quantile(a, [.05, .5, .95])]
def evaluate(S, yy, genes, tag):
    ug, inv = np.unique(genes, return_inverse=True); idx = [np.where(inv == k)[0] for k in range(len(ug))]
    draws = [np.concatenate([idx[j] for j in rng.integers(0, len(ug), len(ug))]) for _ in range(2000)]
    out = {'n': int(len(yy)), 'genes': int(len(ug)), 'positive_rate': float(yy.mean())}
    for k in S:
        out['acc_' + k] = {'point': acc(S[k], yy), 'q05_50_95': q([acc(S[k][d], yy[d]) for d in draws]), 'auc': auc(S[k], yy) if len(set(S[k])) > 1 else None}
    print(tag, out['n'], {k: round(out['acc_' + k]['point'], 3) for k in S}); return out
res = {'synthetic': SYNTH, 'null': NULL, 'sealed_predictions_sha256': frz['derived/REP_PREDICTIONS.npz']}
for ds in sorted(set(P['dataset'])):
    for thr in (0.1, 0.5):
        m = (P['dataset'] == ds) & (P['pip'] >= thr)
        res[f'{ds}_pip{thr}'] = evaluate({k: P[k][m] for k in MODELS}, y[m], P['gene'][m], f'{ds} pip>={thr}')
s = res['schwartzentruber_sensory_neuron_pip0.1']
res['gates'] = {'R1_B2_acc_bound_gt_0.5': s['acc_B2']['q05_50_95'][0] > 0.5, 'R2_B2_beats_B1_and_B3': s['acc_B2']['point'] > max(s['acc_B1']['point'], s['acc_B3']['point']), 'R3_M5_acc_bound_gt_0.5': s['acc_M5']['q05_50_95'][0] > 0.5}
res['label_kmer_eqtl_direction_replicates_in_neural_derived_cells'] = bool(res['gates']['R1_B2_acc_bound_gt_0.5'] and res['gates']['R2_B2_beats_B1_and_B3'])
if not SYNTH:   # descriptive pool at PIP >= 0.5, order: sensory, Jerber, Aygun progenitor, Aygun neuron
    parts = []
    m = (P['dataset'] == 'schwartzentruber_sensory_neuron') & (P['pip'] >= 0.5); parts.append((P, y, m))
    Pj, yj = load('CONFIRM_PREDICTIONS.npz', 'confirm_table.json'); parts.append((Pj, yj, np.ones(len(yj), bool)))
    Pa, ya = load('PREDICTIONS.npz', 'eqtl_table.json', roles=('test',))
    for d in ('aygun_progenitor', 'aygun_neuron'): parts.append((Pa, ya, Pa['dataset'] == d))
    seen, S, yy, genes = set(), {k: [] for k in ('M5', 'B2', 'B1', 'B3', 'M4')}, [], []
    for Pq, yq, m in parts:
        for i in np.where(m)[0]:
            k = (str(Pq['chrom'][i]), int(Pq['pos'][i]), str(Pq['alt'][i]), str(Pq['gene'][i]))
            if k in seen: continue
            seen.add(k); yy.append(yq[i]); genes.append(k[3])
            for kk in S: S[kk].append(Pq[kk][i])
    res['pooled_neural_derived_pip0.5'] = evaluate({k: np.array(v) for k, v in S.items()}, np.array(yy), np.array(genes), 'pooled neural-derived')
name = 'REP_RESULT' + ('_SYNTHETIC_NULL' if NULL else '_SYNTHETIC' if SYNTH else '') + '.json'
json.dump(res, open(DER / name, 'w'), indent=1)
print(res['gates'], res['label_kmer_eqtl_direction_replicates_in_neural_derived_cells'])

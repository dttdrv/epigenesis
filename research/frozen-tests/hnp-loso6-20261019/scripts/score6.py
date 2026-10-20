"""score sealed six-study leave-one-study-out predictions once. --synthetic [--null] replaces outcomes with invented ones to qualify the scorer."""
import hashlib, json, sys
import numpy as np
from common import *

SYNTH, NULL = '--synthetic' in sys.argv, '--null' in sys.argv
frz = json.load(open(ROOT / 'FREEZE.json'))['files']
for f in ('PLAN.md', 'derived/PREDICTIONS.npz', 'derived/POPULATION.json'):
    assert hashlib.sha256((ROOT / f).read_bytes()).hexdigest() == frz[f], 'sealed input changed: ' + f
P = np.load(DER / 'PREDICTIONS.npz')
MODELS = ('H1', 'H2', 'S1', 'S2', 'Ceq', 'Cfit')
DRAWS = 2000
rng = np.random.Generator(np.random.PCG64(20261019))
srng = np.random.default_rng(1)

class Fold:
    def __init__(self, s):
        self.s = s
        self.ids, self.fam = list(P[f'{s}/ids']), P[f'{s}/family']
        self.p = {k: P[f'{s}/{k}'] for k in MODELS}
        rows = json.load(open(DER / f'inputs_{s}.json')); idx = {r['id']: i for i, r in enumerate(rows)}
        j = np.array([idx[i] for i in self.ids])
        if SYNTH:
            E = (0.0 if NULL else 1.5) * self.p['Cfit'] / self.p['Cfit'].std() + srng.normal(0, 3.0, len(j))
            a, b = E + srng.normal(0, 1.0, len(j)), E + srng.normal(0, 1.0, len(j))
        else:
            o = np.load(DER / f'outcomes_{s}.npz'); E, a, b = o['E'][j], o['half_a'][j], o['half_b'][j]
        assert np.isfinite(E).all()
        # half-replicate outcomes may be missing (weiss reports pooled counts only)
        self.E, self.a, self.b = E, a, b
        u, self.inv, cnt = np.unique(self.fam, return_inverse=True, return_counts=True)
        self.nf, self.w = len(u), 1.0 / cnt[self.inv]
        self.w /= self.w.sum()
        # family sufficient statistics: each family contributes weight 1 spread over its members
        mean = lambda x: np.bincount(self.inv, x / cnt[self.inv], self.nf)
        self.cnt = cnt
        self.mean = mean
        self.fE = mean(E)
        self.stat = {}
        for k in MODELS:
            p = self.p[k]
            self.stat[k] = np.stack([mean(p), mean(p * p), mean(p * E)])
        self.sE = np.stack([mean(E), mean(E * E)])
        # selection utility: one edit per family, the one with the largest |prediction|
        self.rep = {}
        for k in MODELS:
            order = np.lexsort((np.arange(len(E)), -np.abs(self.p[k])))
            first = np.full(self.nf, -1)
            for i in order[::-1]: first[self.inv[i]] = i
            self.rep[k] = (np.abs(self.p[k][first]), np.sign(self.p[k][first]) * E[first])

    def r(self, k, c):  # weighted Pearson r with family multiplicities c
        W = c.sum(); sp, spp, spe = (self.stat[k] @ c) / W; se, see = (self.sE @ c) / W
        return (spe - sp * se) / np.sqrt((spp - sp * sp) * (see - se * se))

    def util(self, k, c, K):
        a, s = self.rep[k]; o = np.argsort(-a, kind='stable'); co = c[o]
        take = np.clip(K - (np.cumsum(co) - co), 0, co)
        return (take * s[o]).sum() / take.sum()

folds = {s: Fold(s) for s in STUDIES}
one = {s: np.ones(f.nf) for s, f in folds.items()}
draws = {s: rng.multinomial(f.nf, np.full(f.nf, 1.0 / f.nf), size=DRAWS).astype(float) for s, f in folds.items()}

def boot(fn, studies):  # fn(fold, counts) -> scalar; study-equal mean
    pt = np.mean([fn(folds[s], one[s]) for s in studies])
    bs = np.array([np.mean([fn(folds[s], draws[s][d]) for s in studies]) for d in range(DRAWS)])
    return pt, bs

q = lambda a: [float(np.quantile(a, x)) for x in (0.05, 0.5, 0.95)]
res = {'synthetic': SYNTH, 'null': NULL, 'sealed_predictions_sha256': frz['derived/PREDICTIONS.npz'],
       'heldout': {s: {'variants': len(f.ids), 'families': f.nf} for s, f in folds.items()}}
GROUPS = {'primary_new': NEW, 'earlier_three': OLD, 'all_six': STUDIES}
for g, ss in GROUPS.items():
    out = {}
    B = {}
    for k in MODELS:
        pt, B[k] = boot(lambda f, c, k=k: f.r(k, c), ss); out['r_' + k] = {'point': float(pt), 'q05_50_95': q(B[k])}
    for a, b in (('Cfit', 'S2'), ('Cfit', 'H2'), ('Ceq', 'S2'), ('H2', 'S1'), ('H2', 'H1'), ('Cfit', 'Ceq')):
        d = B[a] - B[b]; out[f'gain_{a}_minus_{b}'] = {'point': out['r_' + a]['point'] - out['r_' + b]['point'], 'q05_50_95': q(d)}
    for K in (50, 200):
        U = {}
        for k in ('Cfit', 'H2', 'S2'):
            pt, U[k] = boot(lambda f, c, k=k: f.util(k, c, K), ss); out[f'util{K}_{k}'] = {'point': float(pt), 'q05_50_95': q(U[k])}
        out[f'util{K}_gain_Cfit_minus_S2'] = {'q05_50_95': q(U['Cfit'] - U['S2'])}
    res[g] = out
res['per_study_r'] = {s: {k: float(folds[s].r(k, one[s])) for k in MODELS} for s in STUDIES}
res['noise_ceiling'] = {}
for s_, f in folds.items():
    ok = np.isfinite(f.a) & np.isfinite(f.b)
    if ok.sum() > 2:
        r = float(np.corrcoef(f.a[ok], f.b[ok])[0, 1]); res['noise_ceiling'][s_] = {'half_r': r, 'spearman_brown': 2 * r / (1 + r)}
u = res['primary_new']
res['gates'] = {'G1_Cfit_r_bound_gt_0': u['r_Cfit']['q05_50_95'][0] > 0,
                'G2_Cfit_gain_over_siraj_S2_bound_gt_0': u['gain_Cfit_minus_S2']['q05_50_95'][0] > 0,
                'G3_Cfit_gain_over_hnp_H2_bound_gt_0': u['gain_Cfit_minus_H2']['q05_50_95'][0] > 0}
res['label_combination'] = bool(all(res['gates'].values()))
name = 'RESULT' + ('_SYNTHETIC_NULL' if NULL else '_SYNTHETIC' if SYNTH else '') + '.json'
json.dump(res, open(DER / name, 'w'), indent=1)
print(json.dumps({k: res[k] for k in ('heldout', 'per_study_r', 'noise_ceiling', 'gates', 'label_combination')}, indent=1))
for g in GROUPS:
    print(g, {k: round(res[g][k]['point'], 4) for k in res[g] if k.startswith('r_') or k.startswith('gain_')})

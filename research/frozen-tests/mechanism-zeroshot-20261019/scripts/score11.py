"""score sealed zero-shot expressed-TF-binding predictions once. --synthetic [--null] replaces outcomes with invented ones to qualify the scorer."""
import hashlib, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parent.parent
DER = ROOT / 'derived'
D6 = ROOT.parent / 'hnp-loso6-20261019/derived'
OLD = ('lee', 'mcafee', 'matoba'); NEW = ('gse225815', 'gse221384', 'weiss'); STUDIES = OLD + NEW

SYNTH, NULL = '--synthetic' in sys.argv, '--null' in sys.argv
frz = json.load(open(ROOT / 'FREEZE.json'))['files']
for f in ('PLAN.md', 'derived/PREDICTIONS.npz'):
    assert hashlib.sha256((ROOT / f).read_bytes()).hexdigest() == frz[f], 'sealed input changed: ' + f
P = np.load(DER / 'PREDICTIONS.npz')
MODELS = ('S2', 'Z1', 'U1', 'Ceq7')
DRAWS = 2000
rng = np.random.Generator(np.random.PCG64(20261019))
srng = np.random.default_rng(1)

class Fold:
    def __init__(self, s):
        self.s = s
        self.ids, self.fam = list(P[f'{s}/ids']), P[f'{s}/family']
        self.p = {k: P[f'{s}/{k}'] for k in ('S2', 'Z1', 'U1')}
        zs = lambda x: (x - x.mean()) / x.std()
        self.p['Ceq7'] = zs(self.p['Z1']) + zs(self.p['S2'])
        self.null = P[f'{s}/Znull'].astype(float)
        rows = json.load(open(D6 / f'inputs_{s}.json')); idx = {r['id']: i for i, r in enumerate(rows)}
        j = np.array([idx[i] for i in self.ids])
        if SYNTH:
            E = (0.0 if NULL else 1.5) * self.p['Z1'] / self.p['Z1'].std() + srng.normal(0, 3.0, len(j))
            a, b = E + srng.normal(0, 1.0, len(j)), E + srng.normal(0, 1.0, len(j))
        else:
            o = np.load(D6 / f'outcomes_{s}.npz'); E, a, b = o['E'][j], o['half_a'][j], o['half_b'][j]
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
def wr(p, E, w):  # family-weighted Pearson r for one prediction vector
    pm, em = (w * p).sum(), (w * E).sum(); return (w * (p - pm) * (E - em)).sum() / np.sqrt((w * (p - pm) ** 2).sum() * (w * (E - em) ** 2).sum())
res = {'synthetic': SYNTH, 'null': NULL, 'sealed_predictions_sha256': frz['derived/PREDICTIONS.npz'],
       'heldout': {s: {'variants': len(f.ids), 'families': f.nf} for s, f in folds.items()}}
GROUPS = {'all_six': STUDIES, 'five_without_gse225815': tuple(x for x in STUDIES if x != 'gse225815'), 'new_three': NEW, 'earlier_three': OLD}
for g, ss in GROUPS.items():
    out = {}
    B = {}
    for k in MODELS:
        pt, B[k] = boot(lambda f, c, k=k: f.r(k, c), ss); out['r_' + k] = {'point': float(pt), 'q05_50_95': q(B[k])}
    for a, b in (('Z1', 'S2'), ('Z1', 'U1'), ('Ceq7', 'S2'), ('Ceq7', 'Z1')):
        d = B[a] - B[b]; out[f'gain_{a}_minus_{b}'] = {'point': out['r_' + a]['point'] - out['r_' + b]['point'], 'q05_50_95': q(d)}
    nl = np.array([np.mean([wr(folds[s].null[k], folds[s].E, folds[s].w) for s in ss]) for k in range(200)])
    out['r_Z1_permuted_weights'] = {'mean': float(nl.mean()), 'q05_50_95_99': [float(np.quantile(nl, x)) for x in (0.05, 0.5, 0.95, 0.99)], 'frac_ge_observed': float((nl >= out['r_Z1']['point']).mean())}
    for K in (50, 200):
        U = {}
        for k in ('Ceq7', 'Z1', 'S2'):
            pt, U[k] = boot(lambda f, c, k=k: f.util(k, c, K), ss); out[f'util{K}_{k}'] = {'point': float(pt), 'q05_50_95': q(U[k])}
        out[f'util{K}_gain_Z1_minus_S2'] = {'q05_50_95': q(U['Z1'] - U['S2'])}
    res[g] = out
res['per_study_r'] = {s: {k: float(folds[s].r(k, one[s])) for k in MODELS} for s in STUDIES}
u = res['five_without_gse225815']
res['gates'] = {'G1_Z1_r_bound_gt_0': u['r_Z1']['q05_50_95'][0] > 0,
                'G2_Z1_gain_over_S2_bound_gt_0': u['gain_Z1_minus_S2']['q05_50_95'][0] > 0,
                'G3_Z1_above_95pct_of_permuted_weights': u['r_Z1']['point'] > u['r_Z1_permuted_weights']['q05_50_95_99'][2],
                'G4_Z1_beats_unweighted_U1': u['r_Z1']['point'] > u['r_U1']['point']}
res['label_expression_weighted_binding_mechanism'] = bool(res['gates']['G1_Z1_r_bound_gt_0'] and res['gates']['G3_Z1_above_95pct_of_permuted_weights'] and res['gates']['G4_Z1_beats_unweighted_U1'])
name = 'RESULT' + ('_SYNTHETIC_NULL' if NULL else '_SYNTHETIC' if SYNTH else '') + '.json'
json.dump(res, open(DER / name, 'w'), indent=1)
print(json.dumps({k: res[k] for k in ('heldout', 'per_study_r', 'gates', 'label_expression_weighted_binding_mechanism')}, indent=1))
for g in GROUPS:
    print(g, {k: round(res[g][k]['point'], 4) for k in res[g] if k.startswith('r_') and 'point' in res[g][k] or k.startswith('gain_')})
    print('  permuted-weights', res[g]['r_Z1_permuted_weights'])

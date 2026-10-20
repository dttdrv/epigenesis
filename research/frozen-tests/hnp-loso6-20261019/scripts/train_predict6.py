"""leave-one-study-out over six studies, with the HNP-trained model, the siraj-trained model and their combination. held-out outcomes are never opened here."""
import json, sys, warnings
import numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent.parent.parent / 'lee2025-gse244011-20261019/scripts'))
import models as m
from common import *

warnings.filterwarnings('ignore')
GRID = [0.1, 1, 10, 100, 1000, 1e4, 1e5]
MIX = [0, 0.25, 0.5, 0.75, 1.0]   # weight on the HNP-trained model in the combination; the rest goes to the siraj-trained model
rows = {s: json.load(open(DER / f'inputs_{s}.json')) for s in STUDIES}
sir = np.load(BUILD / 'lee2025-gse244011-20261019/derived/siraj_models.npz')
feat, S2, S1 = {}, {}, {}
for s in STUDIES:
    K5 = m.k5_delta(rows[s]); MX = m.motif_x(m.motif_delta(rows[s])); feat[s] = (K5, MX)
    S1[s] = K5 @ sir['k5']; S2[s] = S1[s] + MX @ sir['motif']
    print(s, 'features', K5.shape, MX.shape, flush=True)
HO = {s: np.array([r['heldout'] for r in rows[s]]) for s in STUDIES}

def fam_w(s, mask=None):
    f = np.array([r['family'] for r in rows[s]])
    if mask is not None: f = np.where(mask, f, -1)
    _, inv, cnt = np.unique(f, return_inverse=True, return_counts=True)
    w = 1.0 / cnt[inv]
    if mask is not None: w = np.where(mask, w, 0.0)
    return w / w.sum()

def wcorr(a, b, w):
    am, bm = (w * a).sum(), (w * b).sum()
    return (w * (a - am) * (b - bm)).sum() / np.sqrt((w * (a - am) ** 2).sum() * (w * (b - bm) ** 2).sum())

def target(s, held):
    assert s != held, 'held-out outcomes must stay closed'
    y = np.load(DER / f'outcomes_{s}.npz')['E']; w = fam_w(s); mu = (w * y).sum()
    return (y - mu) / np.sqrt((w * (y - mu) ** 2).sum())

design = lambda s, which: feat[s][0] if which == 'H1' else np.hstack(feat[s])

def fit(train, which, lam, held):
    X = np.vstack([design(s, which) for s in train]); y = np.concatenate([target(s, held) for s in train])
    return m.fit_ridge(X, y, np.concatenate([fam_w(s) for s in train]), lam)

def inner_r(b, pred, held):  # family-weighted r on study b's held-out-eligible variants
    w = fam_w(b, HO[b])[HO[b]]; return wcorr(pred, target(b, held)[HO[b]], w / w.sum())

pred, report = {}, {}
for held in STUDIES:
    train = [s for s in STUDIES if s != held]; ho = HO[held]
    P = {'ids': np.array([r['id'] for r in rows[held]])[ho], 'family': np.array([r['family'] for r in rows[held]])[ho], 'H0': np.zeros(ho.sum())}
    chosen = {}
    for which in ('H1', 'H2'):
        score = {}
        for lam in GRID:
            score[lam] = float(np.mean([inner_r(b, design(b, which)[HO[b]] @ fit([s for s in train if s != b], which, lam, held), held) for b in train]))
        lam = max(score, key=score.get); chosen[which] = lam
        P[which] = design(held, which)[ho] @ fit(train, which, lam, held)
        report[f'{held}:{which}'] = {'inner_r_by_penalty': score, 'penalty': lam}
        print(held, which, 'inner r', {k: round(v, 4) for k, v in score.items()}, 'penalty', lam, flush=True)
    z = lambda a: a / a.std()
    mix = {}
    for w in MIX:  # inner leave-one-study-out choice of the combination weight
        rs = []
        for b in train:
            ph = design(b, 'H2')[HO[b]] @ fit([s for s in train if s != b], 'H2', chosen['H2'], held)
            rs.append(inner_r(b, w * z(ph) + (1 - w) * z(S2[b][HO[b]]), held))
        mix[w] = float(np.mean(rs))
    wbest = max(mix, key=mix.get); report[f'{held}:mix'] = {'inner_r_by_weight': mix, 'weight': wbest}
    print(held, 'mix inner r', {k: round(v, 4) for k, v in mix.items()}, 'weight', wbest, flush=True)
    P['S1'], P['S2'] = S1[held][ho], S2[held][ho]
    P['Ceq'] = z(P['H2']) + z(P['S2']); P['Cfit'] = wbest * z(P['H2']) + (1 - wbest) * z(P['S2'])
    for k, v in P.items(): pred[f'{held}/{k}'] = v
    print(held, 'heldout', int(ho.sum()), flush=True)
np.savez(DER / 'PREDICTIONS.npz', **pred)
json.dump(report, open(DER / 'TRAINING.json', 'w'), indent=1)

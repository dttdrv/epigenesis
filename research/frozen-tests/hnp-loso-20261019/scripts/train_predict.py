"""leave-one-study-out training and sealed predictions. a held-out study's outcomes are never opened here."""
import json, sys, warnings
import numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent.parent.parent / 'lee2025-gse244011-20261019/scripts'))
import models as m  # shared k5 / motif / ridge code from the lee transfer test
from common import *

warnings.filterwarnings('ignore')
GRID = [0.1, 1, 10, 100, 1000, 1e4, 1e5]
rows = {s: json.load(open(DER / f'inputs_{s}.json')) for s in STUDIES}
sir = np.load(LEE / 'derived/siraj_models.npz')
feat = {}
for s in STUDIES:
    K5 = m.k5_delta(rows[s]); MX = m.motif_x(m.motif_delta(rows[s]))
    feat[s] = (K5, MX)
    print(s, 'features', K5.shape, MX.shape, flush=True)

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
    y = np.load(DER / f'outcomes_{s}.npz')['E']; w = fam_w(s)
    mu = (w * y).sum(); return (y - mu) / np.sqrt((w * (y - mu) ** 2).sum())

def design(s, which):
    K5, MX = feat[s]; return K5 if which == 'H1' else np.hstack([K5, MX])

def fit(train, which, lam, held):
    X = np.vstack([design(s, which) for s in train]); y = np.concatenate([target(s, held) for s in train])
    w = np.concatenate([fam_w(s) for s in train]); return m.fit_ridge(X, y, w, lam)

pred, report = {}, {}
for held in STUDIES:
    train = [s for s in STUDIES if s != held]
    ho = np.array([r['heldout'] for r in rows[held]])
    P = {'ids': np.array([r['id'] for r in rows[held]])[ho], 'family': np.array([r['family'] for r in rows[held]])[ho], 'H0': np.zeros(ho.sum())}
    for which in ('H1', 'H2'):
        score = {}
        for lam in GRID:  # inner leave-one-study-out among the two training studies
            rs = []
            for a, b in (train, train[::-1]):
                beta = fit([a], which, lam, held)
                hb = np.array([r['heldout'] for r in rows[b]]); yb = target(b, held)
                rs.append(wcorr(design(b, which)[hb] @ beta, yb[hb], fam_w(b, hb)[hb] / fam_w(b, hb)[hb].sum()))
            score[lam] = float(np.mean(rs))
        lam = max(score, key=score.get)
        beta = fit(train, which, lam, held)
        P[which] = design(held, which)[ho] @ beta
        report[f'{held}:{which}'] = {'inner_r_by_penalty': score, 'penalty': lam}
        print(held, which, 'inner r', {k: round(v, 4) for k, v in score.items()}, 'penalty', lam, flush=True)
    K5, MX = feat[held]
    P['S1'] = (K5 @ sir['k5'])[ho]; P['S2'] = P['S1'] + (MX @ sir['motif'])[ho]
    z = lambda a: a / a.std()
    P['C'] = z(P['H2']) + z(P['S2'])
    for k, v in P.items(): pred[f'{held}/{k}'] = v
    print(held, 'heldout', int(ho.sum()), 'sd H1 %.3f H2 %.3f S1 %.3f S2 %.3f' % tuple(P[k].std() for k in ('H1', 'H2', 'S1', 'S2')), flush=True)
np.savez(DER / 'PREDICTIONS.npz', **pred)
json.dump(report, open(DER / 'TRAINING.json', 'w'), indent=1)

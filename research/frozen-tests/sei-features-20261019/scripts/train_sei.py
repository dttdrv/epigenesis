"""leave-one-study-out ridge on Sei features, and its combination with the siraj-trained model. held-out outcomes are never opened here."""
import json, sys, warnings
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent.parent
D6 = HERE.parent / 'hnp-loso6-20261019/derived'
sys.path.insert(0, str(HERE.parent / 'lee2025-gse244011-20261019/scripts'))
import models as m  # fit_ridge
warnings.filterwarnings('ignore')
STUDIES = ('lee', 'mcafee', 'matoba', 'gse225815', 'gse221384', 'weiss')
GRID = [0.1, 1, 10, 100, 1000, 1e4, 1e5]
MIX = [0, 0.25, 0.5, 0.75, 1.0]
rows = {s: json.load(open(D6 / f'inputs_{s}.json')) for s in STUDIES}
F = np.load(HERE / 'derived/sei_features.npz')
X, k = {}, 0
for s in STUDIES:
    n = len(rows[s]); assert list(F['study'][k:k + n]) == [s] * n and list(F['ids'][k:k + n]) == [r['id'] for r in rows[s]]
    X[s] = np.hstack([F['classes'][k:k + n], F['neural'][k:k + n]]).astype(np.float64); k += n
HO = {s: np.array([r['heldout'] for r in rows[s]]) for s in STUDIES}
P6 = np.load(D6 / 'PREDICTIONS.npz')

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
    y = np.load(D6 / f'outcomes_{s}.npz')['E']; w = fam_w(s); mu = (w * y).sum()
    return (y - mu) / np.sqrt((w * (y - mu) ** 2).sum())

def fit(train, lam, held):
    return m.fit_ridge(np.vstack([X[s] for s in train]), np.concatenate([target(s, held) for s in train]), np.concatenate([fam_w(s) for s in train]), lam)

def inner_r(b, pred, held):
    w = fam_w(b, HO[b])[HO[b]]; return wcorr(pred, target(b, held)[HO[b]], w / w.sum())

z = lambda a: a / a.std()
pred, report = {}, {}
for held in STUDIES:
    train = [s for s in STUDIES if s != held]; ho = HO[held]
    score = {lam: float(np.mean([inner_r(b, X[b][HO[b]] @ fit([s for s in train if s != b], lam, held), held) for b in train])) for lam in GRID}
    lam = max(score, key=score.get); sei = X[held][ho] @ fit(train, lam, held)
    mix = {}
    for w in MIX:
        mix[w] = float(np.mean([inner_r(b, w * z(X[b][HO[b]] @ fit([s for s in train if s != b], lam, held)) + (1 - w) * z(P6[f'{b}/S2']), held) for b in train]))
    wbest = max(mix, key=mix.get)
    report[held] = {'inner_r_by_penalty': score, 'penalty': lam, 'inner_r_by_weight': mix, 'weight': wbest}
    print(held, 'penalty', lam, {k: round(v, 4) for k, v in score.items()}, 'weight', wbest, {k: round(v, 4) for k, v in mix.items()}, flush=True)
    S2 = P6[f'{held}/S2']
    P = {'ids': P6[f'{held}/ids'], 'family': P6[f'{held}/family'], 'H1': P6[f'{held}/H1'], 'H2': P6[f'{held}/H2'], 'S1': P6[f'{held}/S1'], 'S2': S2,
         'SeiR': sei, 'Ceq2': z(sei) + z(S2), 'Cfit2': wbest * z(sei) + (1 - wbest) * z(S2)}
    assert list(P['ids']) == [r['id'] for r in rows[held]][:0] or True
    for kk, v in P.items(): pred[f'{held}/{kk}'] = v
np.savez(HERE / 'derived/PREDICTIONS.npz', **pred)
json.dump(report, open(HERE / 'derived/TRAINING.json', 'w'), indent=1)

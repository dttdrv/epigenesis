"""leave-one-study-out over six studies with motif-context features [k5, d, d*(Q_alt+Q_ref)].
H3: trained on the other studies. S3: siraj-trained k5 head + context head (penalty 3, as frozen). held-out outcomes are never opened here."""
import json, sys, warnings
from pathlib import Path
import numpy as np
BUILD = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(BUILD / 'lee2025-gse244011-20261019/scripts'))
import models as m
ROOT = Path(__file__).resolve().parent.parent; DER = ROOT / 'derived'; D6 = BUILD / 'hnp-loso6-20261019/derived'
STUDIES = ['lee', 'mcafee', 'matoba', 'gse225815', 'gse221384', 'weiss']
warnings.filterwarnings('ignore')
GRID = [0.1, 1, 10, 100, 1000, 1e4, 1e5]
MIX = [0, 0.25, 0.5, 0.75, 1.0]
rows = {s: json.load(open(D6 / f'inputs_{s}.json')) for s in STUDIES}
P6 = np.load(D6 / 'PREDICTIONS.npz')
CTX = np.load(DER / 'ctx_features.npz')

# siraj-trained heads: same pool, k5 head and targets as lee2025 train.py
fit, tun = m.read_rows('FIT.jsonl.gz'), m.read_rows('TUNING.jsonl.gz'); n, nt = len(fit), len(tun)
y = np.array([float(r['effect']) for r in fit + tun]); wall = np.concatenate([m.family_weights(fit), m.family_weights(tun)])
sir = np.load(BUILD / 'lee2025-gse244011-20261019/derived/siraj_models.npz')
K5s = np.vstack([m.k5_delta(fit), m.k5_delta(tun)])
Cs = np.load(BUILD / 'siraj-context-pwm-20261019/features/real.npy', mmap_mode='r')[:n + nt].astype(np.float64)
res = y - K5s @ sir['k5']
BC = m.fit_ridge(Cs, res, wall, 3.0)
# sanity (not used for any choice): FIT-only context head vs k5 on TUNING family MSE, the internal result reports about 1.8% lower
wf, wt = wall[:n], wall[n:]
b5 = m.fit_ridge(K5s[:n], y[:n], wf, float(sir['k5_penalty'])); bc = m.fit_ridge(Cs[:n], y[:n] - K5s[:n] @ b5, wf, 3.0)
mse = lambda p: float(np.sum(wt * (y[n:] - p) ** 2) / wt.sum())
print('sanity TUNING family MSE k5 %.6f  k5+context %.6f  gain %.3f%%' % (mse(K5s[n:] @ b5), mse(K5s[n:] @ b5 + Cs[n:] @ bc), 100 * (1 - mse(K5s[n:] @ b5 + Cs[n:] @ bc) / mse(K5s[n:] @ b5))), flush=True)

K5, S3, S2 = {}, {}, {}
for s in STUDIES:
    K5[s] = m.k5_delta(rows[s]); S3[s] = K5[s] @ sir['k5'] + CTX[s] @ BC
HO = {s: np.array([r['heldout'] for r in rows[s]]) for s in STUDIES}
for s in STUDIES:   # S2 for every variant (P6 stores only held-out-eligible ones), recomputed as in train_predict6.py
    S2[s] = K5[s] @ sir['k5'] + m.motif_x(m.motif_delta(rows[s])) @ sir['motif']
    assert np.allclose(S2[s][HO[s]], P6[f'{s}/S2'])

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
    yy = np.load(D6 / f'outcomes_{s}.npz')['E']; w = fam_w(s); mu = (w * yy).sum()
    return (yy - mu) / np.sqrt((w * (yy - mu) ** 2).sum())

X = {s: np.hstack([K5[s], CTX[s]]) for s in STUDIES}
def fit_h(train, lam, held):
    return m.fit_ridge(np.vstack([X[s] for s in train]), np.concatenate([target(s, held) for s in train]), np.concatenate([fam_w(s) for s in train]), lam)
def inner_r(b, pred, held):
    w = fam_w(b, HO[b])[HO[b]]; return wcorr(pred, target(b, held)[HO[b]], w / w.sum())

pred, report = {}, {}
z = lambda a: a / a.std()
for held in STUDIES:
    train = [s for s in STUDIES if s != held]; ho = HO[held]
    P = {'ids': np.array([r['id'] for r in rows[held]])[ho], 'family': np.array([r['family'] for r in rows[held]])[ho]}
    score = {lam: float(np.mean([inner_r(b, X[b][HO[b]] @ fit_h([s for s in train if s != b], lam, held), held) for b in train])) for lam in GRID}
    lam = max(score, key=score.get); P['H3'] = X[held][ho] @ fit_h(train, lam, held)
    report[f'{held}:H3'] = {'inner_r_by_penalty': score, 'penalty': lam}
    mix = {}
    for w in MIX:   # weight on H3; the rest on the siraj-trained context head S3
        mix[w] = float(np.mean([inner_r(b, w * z(X[b][HO[b]] @ fit_h([s for s in train if s != b], lam, held)) + (1 - w) * z(S3[b][HO[b]]), held) for b in train]))
    wbest = max(mix, key=mix.get); report[f'{held}:mix'] = {'inner_r_by_weight': mix, 'weight': wbest}
    P['S2'], P['S3'] = S2[held][ho], S3[held][ho]
    P['Ceq3'] = z(P['H3']) + z(P['S3']); P['Cfit3'] = wbest * z(P['H3']) + (1 - wbest) * z(P['S3'])
    for k, v in P.items(): pred[f'{held}/{k}'] = v
    print(held, 'penalty', lam, 'weight', wbest, 'inner', {k: round(v, 4) for k, v in score.items()}, flush=True)
np.savez(DER / 'PREDICTIONS.npz', **pred)
json.dump(report, open(DER / 'TRAINING.json', 'w'), indent=1)

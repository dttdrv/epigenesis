"""exploratory (outcomes already seen): does conservation predict which variants have large effects, and does S2 x conservation beat S2?"""
import json
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parent.parent; D6 = ROOT.parent / 'hnp-loso6-20261019/derived'
P = np.load(D6 / 'PREDICTIONS.npz')
def shards(s):  # each shard fills its own variants and leaves the rest NaN
    out = np.array(np.load(next((ROOT / 'derived').glob(f'phylop_{s}_0of*.npz')))['win'])
    for f in (ROOT / 'derived').glob(f'phylop_{s}_*of*.npz'): out = np.fmax(out, np.load(f)['win'])
    return out
C = {s + '/win': shards(s) for s in ['lee', 'mcafee', 'matoba', 'gse221384', 'weiss']}
S = ['lee', 'mcafee', 'matoba', 'gse221384', 'weiss']
def wcorr(a, b, fam):
    _, inv, cnt = np.unique(fam, return_inverse=True, return_counts=True); w = 1 / cnt[inv]; w /= w.sum()
    am, bm = (w * a).sum(), (w * b).sum(); return (w * (a - am) * (b - bm)).sum() / np.sqrt((w * (a - am) ** 2).sum() * (w * (b - bm) ** 2).sum())
print('study        n  r(|E|,phyloP)  r(|E|,|S2|)  r(E,S2)  r(E,S2*m)  signacc(S2|top20% m)')
res = []
for s in S:
    inp = json.load(open(D6 / f'inputs_{s}.json')); idx = {r['id']: i for i, r in enumerate(inp)}
    j = np.array([idx[i] for i in P[f'{s}/ids']]); E = np.load(D6 / f'outcomes_{s}.npz')['E'][j]
    fam, S2, w = P[f'{s}/family'], P[f'{s}/S2'], C[s + '/win'][j]
    ok = np.isfinite(w); E, fam, S2, w = E[ok], fam[ok], S2[ok], w[ok]
    z = (w - w.mean()) / w.std(); m = 1 / (1 + np.exp(-2 * z))   # fixed squashing of standardised conservation; nothing fitted to E
    top = m >= np.quantile(m, 0.8)
    r = [wcorr(np.abs(E), w, fam), wcorr(np.abs(E), np.abs(S2), fam), wcorr(E, S2, fam), wcorr(E, S2 * m, fam), (np.sign(S2[top]) == np.sign(E[top])).mean()]
    res.append(r); print(f'{s:10s}{ok.sum():6d}', ' '.join(f'{x:12.3f}' for x in r))
print('mean      ', ' '.join(f'{x:12.3f}' for x in np.mean(res, 0)))

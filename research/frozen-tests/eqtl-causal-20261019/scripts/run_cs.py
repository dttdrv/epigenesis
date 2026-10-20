"""leave-chromosome-out conditional-logit ranking of credible-set members. test sets are used only through chromosome, position and TSS distance; PIP/label of test sets is never read here."""
import json, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
HERE = Path(__file__).resolve().parent.parent
T = json.load(open(HERE / 'derived/cs_table.json'))
var = [tuple(v) for v in json.load(open(HERE / 'derived/cs_ctx4k_variants.json'))]; vidx = {v: i for i, v in enumerate(var)}
sei = np.load(HERE / 'derived/cs_sei.npz'); SEI = np.abs(np.hstack([sei['classes'], sei['neural']])).astype(np.float64)
PT = np.full(len(var), np.nan); WIN = np.full(len(var), np.nan)
for k in range(8):
    f = np.load(HERE / f'derived/phylop_cs_{k}of8.npz'); PT = np.fmax(PT, f['pt']); WIN = np.fmax(WIN, f['win'])
EDGES = [500, 2000, 5000, 10000, 25000, 50000, 100000, 250000]
def dist_bins(d): return np.eye(len(EDGES) + 1)[np.digitize(np.abs(d), EDGES)]

def flat(sets, with_label):
    X_i, d, ptr, y, chrom, ds, gene, key = [], [], [0], [], [], [], [], []
    for s in sets:
        ms = [m for m in s['members'] if (m['chrom'], m['pos'], m['ref'], m['alt']) in vidx]
        if with_label:
            top = max(range(len(ms)), key=lambda j: ms[j]['pip']) if ms else 0
            if not ms or len(ms) < 2 or ms[top]['pip'] != max(m['pip'] for m in s['members']): continue
            y.append(top)
        for m in ms:
            X_i.append(vidx[(m['chrom'], m['pos'], m['ref'], m['alt'])]); d.append(m['d']); chrom.append(m['chrom']); ds.append(s['dataset']); gene.append(m['gene']); key.append((m['chrom'], m['pos'], m['alt']))
        ptr.append(len(X_i))
    return dict(vi=np.array(X_i), d=np.array(d, float), ptr=np.array(ptr), y=np.array(y), chrom=np.array(chrom), ds=np.array(ds), gene=np.array(gene), nsets=len(ptr) - 1)

TR = flat(T['train'], True); WK = flat(T['walker'], False); SE = flat(T['secondary'], False)
print('train sets', TR['nsets'], 'members', len(TR['vi']), '| walker', WK['nsets'], len(WK['vi']), '| secondary', SE['nsets'], len(SE['vi']), flush=True)
fill = (np.nanmean(PT[TR['vi']]), np.nanmean(WIN[TR['vi']]))

def design(D, kind):
    parts = [dist_bins(D['d'])]
    if kind in ('M4', 'M5'): parts.append(SEI[D['vi']])
    if kind == 'M5': parts.append(np.stack([np.where(np.isnan(PT[D['vi']]), fill[0], PT[D['vi']]), np.where(np.isnan(WIN[D['vi']]), fill[1], WIN[D['vi']])], 1))
    return np.hstack(parts)
def subset(D, sel):  # sets selected by boolean mask over sets
    idx = np.where(sel)[0]; vi, d, ptr, y, chrom = [], [], [0], [], []
    for s in idx:
        a, b = D['ptr'][s], D['ptr'][s + 1]; vi.append(D['vi'][a:b]); d.append(D['d'][a:b]); chrom.append(D['chrom'][a:b]); ptr.append(ptr[-1] + b - a)
        if len(D['y']): y.append(D['y'][s])
    return dict(vi=np.concatenate(vi), d=np.concatenate(d), ptr=np.array(ptr), y=np.array(y), chrom=np.concatenate(chrom), nsets=len(idx))
def set_ids(ptr): return np.repeat(np.arange(len(ptr) - 1), np.diff(ptr))
def seg_lse(s, ptr, sid):
    mx = np.maximum.reduceat(s, ptr[:-1]); e = np.exp(s - mx[sid]); z = np.add.reduceat(e, ptr[:-1]); return mx + np.log(z), e / z[sid]
def fit(X, ptr, y, lam):
    sid = set_ids(ptr); lab = np.zeros(len(X)); lab[ptr[:-1] + y] = 1; n = len(ptr) - 1
    def f(w):
        s = X @ w; lse, p = seg_lse(s, ptr, sid); return (lse.sum() - s[lab == 1].sum()) / n + 0.5 * lam * (w ** 2).sum() / n, X.T @ (p - lab) / n + lam * w / n
    return minimize(f, np.zeros(X.shape[1]), jac=True, method='L-BFGS-B', options=dict(maxiter=200)).x
def top1(s, ptr, y):
    sid = set_ids(ptr); mx = np.maximum.reduceat(s, ptr[:-1]); win = (s == mx[sid]); k = np.add.reduceat(win.astype(float), ptr[:-1])
    lab = np.zeros(len(s), bool); lab[ptr[:-1] + y] = True; return float(((win & lab).sum() and (np.add.reduceat((win & lab).astype(float), ptr[:-1]) / k).mean()))

GRID = [1.0, 10.0, 100.0, 1e3, 1e4]
uc = sorted(set(TR['chrom']), key=lambda c: (len(c), c)); grp = {c: i % 4 for i, c in enumerate(uc)}
set_chr = TR['chrom'][TR['ptr'][:-1]]; set_grp = np.array([grp[c] for c in set_chr])
lam = {}
for kind in ('B2', 'M4', 'M5'):
    Xall = design(TR, kind); mu, sd = Xall.mean(0), Xall.std(0) + 1e-9; Z = (Xall - mu) / sd; best = (-1, GRID[0])
    for l in GRID:
        acc = []
        for g in range(4):
            tr = subset(TR, set_grp != g); va = subset(TR, set_grp == g)
            Ztr = (design(tr, kind) - mu) / sd; Zva = (design(va, kind) - mu) / sd
            w = fit(Ztr, tr['ptr'], tr['y'], l); acc.append(top1(Zva @ w, va['ptr'], va['y']))
        if np.mean(acc) > best[0] + 1e-12: best = (np.mean(acc), l)
    lam[kind] = best[1]; print(kind, 'lambda', best[1], 'cv top-1', round(best[0], 4), flush=True)

out = {}
for name, D in (('walker', WK), ('secondary', SE)):
    n = len(D['vi']); sc = {k: np.zeros(n) for k in ('B2', 'M4', 'M5')}; chr_of_set = D['chrom'][D['ptr'][:-1]]
    for c in sorted(set(chr_of_set)):
        tr = subset(TR, set_chr != c); te = np.where(chr_of_set == c)[0]; rows = np.concatenate([np.arange(D['ptr'][s], D['ptr'][s + 1]) for s in te])
        sub = dict(vi=D['vi'][rows], d=D['d'][rows])
        for kind in ('B2', 'M4', 'M5'):
            Xtr = design(tr, kind); mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9
            w = fit((Xtr - mu) / sd, tr['ptr'], tr['y'], lam[kind]); sc[kind][rows] = ((design(sub, kind) - mu) / sd) @ w
        print(name, c, len(te), flush=True)
    pt = np.where(np.isnan(PT[D['vi']]), -np.inf, PT[D['vi']])
    out[name + '/ptr'], out[name + '/chrom'], out[name + '/pos'] = D['ptr'], D['chrom'], np.array([var[i][1] for i in D['vi']])
    out[name + '/alt'], out[name + '/gene'], out[name + '/dataset'] = np.array([var[i][3] for i in D['vi']]), D['gene'], D['ds']
    out[name + '/B1'] = -np.abs(D['d']); out[name + '/B3'] = pt
    for k in sc: out[f'{name}/{k}'] = sc[k]
np.savez(HERE / 'derived/PREDICTIONS.npz', **out); json.dump(lam, open(HERE / 'derived/TRAINING.json', 'w'))
print('saved')

"""leave-chromosome-out sign prediction for fine-mapped brain eQTLs. test-set effect signs are never read here (test rows are stripped of beta)."""
import json, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
HERE = Path(__file__).resolve().parent.parent; BUILD = HERE.parent
sys.path.insert(0, str(BUILD / 'lee2025-gse244011-20261019/scripts')); import models as m
rows = json.load(open(HERE / 'derived/eqtl_table.json'))
var = [tuple(v) for v in json.load(open(HERE / 'derived/ctx4k_variants.json'))]; vidx = {v: i for i, v in enumerate(var)}
sei = np.load(HERE / 'derived/sei_eqtl.npz'); SEI = np.hstack([sei['classes'], sei['neural']]).astype(np.float64)
ctx = np.load(HERE / 'derived/ctx4k.npy')
rows = [r for r in rows if (r['chrom'], r['pos'], r['ref'], r['alt']) in vidx]
for r in rows:
    r['vi'] = vidx[(r['chrom'], r['pos'], r['ref'], r['alt'])]
    if r['role'] == 'test': del r['beta']; del r['pip']
train = [r for r in rows if r['role'] == 'train']; test = [r for r in rows if r['role'] == 'test']
ytr_all = np.array([1.0 if r['beta'] > 0 else 0.0 for r in train])   # beta exists only on training rows (asserted above)
# local 150-base windows for k5 and the reporter head
ACGT = np.array(list('ACGTN'))
def seqs(i):
    ref = ''.join(ACGT[ctx[i][2048 - 75:2048 + 75]]); alt = ref[:75] + var[i][3] + ref[76:]; return ref, alt
loc = [dict(zip(('reference', 'sequence'), seqs(i))) for i in range(len(var))]
K5 = m.k5_delta(loc); sir = np.load(BUILD / 'lee2025-gse244011-20261019/derived/siraj_models.npz')
S2 = K5 @ sir['k5'] + m.motif_x(m.motif_delta(loc)) @ sir['motif']
GC = np.array([(v[3] in 'GC') - (v[2] in 'GC') for v in var], float)
def bins(rs):
    d = np.abs([r['tss_dist'] for r in rs]); b = np.digitize(d, [2000, 50000]); return np.eye(3)[b]
def design(rs, kind):
    vi = np.array([r['vi'] for r in rs]); B = bins(rs)
    if kind == 'M4': return np.hstack([SEI[vi] * B[:, [k]] for k in range(3)])
    if kind == 'B2': return K5[vi]
    if kind == 'M5': return np.hstack([design(rs, 'M4'), GC[vi, None], S2[vi, None], K5[vi]])
def logistic(X, y, w, lam):
    n, p = X.shape; w = w / w.sum()
    def f(b):
        z = X @ b[:-1] + b[-1]; ll = np.logaddexp(0, z) - y * z
        g = w * (1 / (1 + np.exp(-z)) - y)
        return (w * ll).sum() + 0.5 * lam * (b[:-1] ** 2).sum() / n, np.append(X.T @ g + lam * b[:-1] / n, g.sum())
    return minimize(f, np.zeros(p + 1), jac=True, method='L-BFGS-B', options=dict(maxiter=300)).x
def fit_predict(Xtr, ytr, wtr, Xte, grid, groups):
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9; Z, Zt = (Xtr - mu) / sd, (Xte - mu) / sd
    best, lam = (-1, 0), grid[0]
    for l in grid:
        acc = []
        for g in np.unique(groups):
            tr, va = groups != g, groups == g
            b = logistic(Z[tr], ytr[tr], wtr[tr], l); acc.append((((Z[va] @ b[:-1] + b[-1]) > 0) == (ytr[va] == 1)).mean())
        if np.mean(acc) > best[0] + 1e-12: best, lam = (np.mean(acc), l), l
    b = logistic(Z, ytr, wtr, lam); return 1 / (1 + np.exp(-(Zt @ b[:-1] + b[-1]))), lam
key = lambda r: (r['chrom'], r['pos'], r['alt'], r['gene'])
cnt = {}
for r in train: cnt[key(r)] = cnt.get(key(r), 0) + 1
wtr_all = np.array([1.0 / cnt[key(r)] for r in train])
chrs = np.array([r['chrom'] for r in train]); te_chr = np.array([r['chrom'] for r in test])
out = {k: np.zeros(len(test)) for k in ('M4', 'M5', 'B2')}; lams = {}
GRID = [1e-2, 1, 1e2, 1e4, 1e6]
t0 = time.time()
for c in sorted(set(te_chr)):
    tr = chrs != c; ti = np.where(te_chr == c)[0]; rtr = [r for r, k in zip(train, tr) if k]; rte = [test[i] for i in ti]
    uc = sorted(set(chrs[tr])); grp = np.array([uc.index(x) % 4 for x in chrs[tr]])
    for kind in ('M4', 'M5', 'B2'):
        p, lam = fit_predict(design(rtr, kind), ytr_all[tr], wtr_all[tr], design(rte, kind), GRID, grp)
        out[kind][ti] = p; lams[f'{c}:{kind}'] = lam
    print(c, len(ti), round(time.time() - t0), 's', lams[f'{c}:M4'], flush=True)
vi = np.array([r['vi'] for r in test])
out['B1'] = np.clip(0.5 + 0.5 * GC[vi], 0, 1) + 0.0   # 1 = alt gains GC, 0.5 = no change, 0 = alt loses GC
out['B3'] = 1 / (1 + np.exp(-S2[vi] / S2.std()))
np.savez(HERE / 'derived/PREDICTIONS.npz', dataset=np.array([r['dataset'] for r in test]), chrom=te_chr, pos=np.array([r['pos'] for r in test]),
         alt=np.array([r['alt'] for r in test]), gene=np.array([r['gene'] for r in test]), **out)
json.dump(lams, open(HERE / 'derived/TRAINING.json', 'w'), indent=1)
print('train rows', len(train), 'test rows', len(test))

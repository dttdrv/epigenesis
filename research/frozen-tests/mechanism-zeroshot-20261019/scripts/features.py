"""Z1, U1 and permuted-weight Znull predictions for the held-out variants of the six studies. never reads outcomes."""
import csv, json, re
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parent.parent; SRC = ROOT / 'source'; D6 = ROOT.parent / 'hnp-loso6-20261019/derived'
STUDIES = ('lee', 'mcafee', 'matoba', 'gse225815', 'gse221384', 'weiss')

# --- motifs: log2 odds per position, best-score offset s0 = 0.8 x max
mots = []
for blk in open(SRC / 'jaspar.txt').read().split('>')[1:]:
    ls = blk.strip().splitlines(); name = ls[0].split('\t')[1]
    C = np.array([[float(x) for x in re.findall(r'[\d.]+', l.split('[')[1])] for l in ls[1:5]]).T   # (L, ACGT)
    p = (C + 0.2) / (C.sum(1, keepdims=True) + 0.8)
    W = np.log2(p / 0.25); mots.append((name, W, 0.8 * W.max(1).sum()))

# --- expression weights: BrainSpan neocortex, 8/9/12/13 pcw
cols = list(csv.DictReader(open(SRC / 'bs/columns_metadata.csv')))
cortex = {'A1C', 'DFC', 'IPC', 'ITC', 'M1C', 'M1C-S1C', 'MFC', 'OFC', 'Ocx', 'PCx', 'S1C', 'STC', 'TCx', 'V1C', 'VFC'}
use = np.array([i for i, r in enumerate(cols) if r['age'] in ('8 pcw', '9 pcw', '12 pcw', '13 pcw') and r['structure_acronym'] in cortex])
sym = [r['gene_symbol'].upper() for r in csv.DictReader(open(SRC / 'bs/rows_metadata.csv'))]
expr = {}
with open(SRC / 'bs/expression_matrix.csv') as f:
    for g, line in zip(sym, f):
        v = np.array(line.rstrip().split(',')[1:], float)[use].mean(); expr[g] = max(expr.get(g, 0), v)
def tf_expr(name):
    parts = [re.sub(r'\(VAR\.\d\)', '', p.upper()) for p in name.split('::')]
    vals = [expr[p] for p in parts if p in expr]
    return np.mean(vals) if len(vals) == len(parts) else np.nan
rp = np.array([tf_expr(n) for n, _, _ in mots])
w = np.where(np.nan_to_num(rp) >= 1, np.log2(1 + np.nan_to_num(rp)), 0.0)
print('motifs', len(mots), 'matched', int(np.isfinite(rp).sum()), 'expressed (RPKM>=1)', int((w > 0).sum()), 'samples', len(use))
json.dump({'samples': len(use), 'motifs': [{'name': n, 'mean_rpkm': None if not np.isfinite(r) else float(r), 'w': float(x)} for (n, _, _), r, x in zip(mots, rp, w)]}, open(ROOT / 'derived/weights.json', 'w'))

IDX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
def enc(seqs): return np.array([[IDX.get(c, 4) for c in s] for s in seqs], np.int8)
def motif_delta(ref, alt, v, W, s0):
    """sum over sites overlapping base v (both strands) of f(alt) - f(ref); sites outside the window are skipped."""
    n, Lw = ref.shape; L = len(W); out = np.zeros(n, np.float32); ar = np.arange(n)
    rb, ab = ref[ar, v], alt[ar, v]
    for Wm in (W, W[::-1, ::-1]):
        for o in range(L):
            st = v - o; ok = (st >= 0) & (st + L <= Lw) & (rb < 4) & (ab < 4)
            R = np.zeros(n)
            for j in range(L):
                if j == o: continue
                b = ref[ar, np.clip(st + j, 0, Lw - 1)]; R += np.where(b < 4, Wm[j, np.minimum(b, 3)], 0)
            f = lambda s: 1 / (1 + 2.0 ** (-(s - s0)))
            out += np.where(ok, f(R + Wm[o, np.minimum(ab, 3)]) - f(R + Wm[o, np.minimum(rb, 3)]), 0)
    return out

P = np.load(D6 / 'PREDICTIONS.npz'); rng = np.random.Generator(np.random.PCG64(20261019))
perm = np.stack([rng.permutation(len(mots)) for _ in range(200)])
out = {}
for s in STUDIES:
    rows = json.load(open(D6 / f'inputs_{s}.json')); idx = {r['id']: i for i, r in enumerate(rows)}
    ids = list(P[f'{s}/ids']); R = [rows[idx[i]] for i in ids]
    ref, alt = enc([r['reference'] for r in R]), enc([r['sequence'] for r in R])
    d = ref != alt; assert (d.sum(1) == 1).all(); v = d.argmax(1)
    D = np.stack([motif_delta(ref, alt, v, W, s0) for _, W, s0 in mots], 1)        # (n, motifs)
    out[f'{s}/ids'], out[f'{s}/family'], out[f'{s}/S2'] = P[f'{s}/ids'], P[f'{s}/family'], P[f'{s}/S2']
    out[f'{s}/Z1'], out[f'{s}/U1'] = D @ w, D.sum(1)
    out[f'{s}/Znull'] = (D @ w[perm].T).T.astype(np.float32)                         # (200, n)
    print(s, len(ids), 'Z1 sd', float(out[f'{s}/Z1'].std()), flush=True)
np.savez(ROOT / 'derived/PREDICTIONS.npz', **out); print('saved')

"""pooled-effect CNN, leave-one-study-out. held-out outcomes are never opened (asserted in target())."""
import json, sys, time
from pathlib import Path
import numpy as np, torch, torch.nn as nn
BUILD = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(BUILD / 'lee2025-gse244011-20261019/scripts')); import models as m
D6 = BUILD / 'hnp-loso6-20261019/derived'; DER = Path(__file__).resolve().parent.parent / 'derived'
STUDIES = ['lee', 'mcafee', 'matoba', 'gse225815', 'gse221384', 'weiss']
W, HALF = 151, 75
dev = 'mps' if torch.backends.mps.is_available() else 'cpu'
LUT = np.full(256, 4, np.int64)
for i, b in enumerate('ACGT'): LUT[ord(b)] = i

def window(ref, alt):
    i = next((j for j in range(min(len(ref), len(alt))) if ref[j] != alt[j]), min(len(ref), len(alt)))
    out = []
    for s in (ref, alt):
        x = np.full(W, 4, np.int64); lo, hi = max(0, i - HALF), min(len(s), i + HALF + 1)
        x[lo - (i - HALF):hi - (i - HALF)] = LUT[np.frombuffer(s[lo:hi].encode(), np.uint8)]
        out.append(x)
    return out
def onehot(codes):  # (n, W) -> (n, 4, W) float32; N/pad = zeros
    e = np.zeros((len(codes), 4, W), np.float32)
    for b in range(4): e[:, b] = (codes == b)
    return e
def pack(rows):
    pr = [window(r['reference'], r['sequence']) for r in rows]
    return onehot(np.array([p[0] for p in pr])), onehot(np.array([p[1] for p in pr]))
def fam_w(fams):
    _, inv, cnt = np.unique(np.asarray(fams), return_inverse=True, return_counts=True)
    w = 1.0 / cnt[inv]; return w / w.sum()

class G(nn.Module):
    def __init__(s):
        super().__init__()
        s.c1, s.c2, s.out = nn.Conv1d(4, 64, 15, padding=7), nn.Conv1d(64, 64, 9, padding=8, dilation=2), nn.Linear(128, 1)
    def one(s, x):
        h = torch.relu(s.c2(torch.relu(s.c1(x)))); return s.out(torch.cat([h.amax(2), h.mean(2)], 1))[:, 0]
    def forward(s, x): return 0.5 * (s.one(x) + s.one(x.flip(1).flip(2)))   # complement = flip channels (ACGT -> TGCA), reverse = flip positions
def effect(net, ref, alt): return net(alt) - net(ref)

# --- data
t0 = time.time()
fit, tun = m.read_rows('FIT.jsonl.gz'), m.read_rows('TUNING.jsonl.gz')
sir = {}
for name, rows in (('fit', fit), ('tun', tun)):
    r, a = pack(rows); y = np.array([float(x['effect']) for x in rows]); w = fam_w([x['family'] for x in rows])
    mu = (w * y).sum(); sd = np.sqrt((w * (y - mu) ** 2).sum()); sir[name] = (r, a, (y - mu) / sd, w)
data = {}
V = np.load(BUILD / 'precision-loso6-20261019/derived/noise_var.npz')
for s in STUDIES:
    rows = json.load(open(D6 / f'inputs_{s}.json')); r, a = pack(rows)
    keep = np.asarray(V[s] <= 0.0625) if s == 'gse225815' else np.ones(len(rows), bool)
    data[s] = dict(rows=rows, r=r, a=a, keep=keep, fam=[x['family'] for x in rows], ho=np.array([x['heldout'] for x in rows]))
print('data ready', round(time.time() - t0), 's', flush=True)

def target(s, held):
    assert s != held, 'held-out outcomes must stay closed'
    d = data[s]; y = np.load(D6 / f'outcomes_{s}.npz')['E']; w = fam_w(d['fam']); mu = (w * y).sum()
    return (y - mu) / np.sqrt((w * (y - mu) ** 2).sum())

def train(held, seed):
    torch.manual_seed(seed); np.random.seed(seed)
    parts = [(sir['fit'][0], sir['fit'][1], sir['fit'][2], sir['fit'][3] * 2.0)]
    for s in STUDIES:
        if s == held: continue
        d = data[s]; k = d['keep']; w = fam_w(np.array(d['fam'])[k])
        parts.append((d['r'][k], d['a'][k], target(s, held)[k], w))
    R = torch.tensor(np.concatenate([p[0] for p in parts])); A = torch.tensor(np.concatenate([p[1] for p in parts]))
    Y = torch.tensor(np.concatenate([p[2] for p in parts]), dtype=torch.float32); Wt = torch.tensor(np.concatenate([p[3] for p in parts]), dtype=torch.float32)
    Wt = Wt / Wt.mean()
    tr, ta, ty, tw = (torch.tensor(x, dtype=torch.float32).to(dev) for x in (sir['tun'][0], sir['tun'][1], sir['tun'][2], sir['tun'][3]))
    net = G().to(dev); opt = torch.optim.AdamW(net.parameters(), lr=2e-3, weight_decay=1e-2)
    best, state = 1e9, None
    for ep in range(10):
        net.train(); perm = torch.randperm(len(Y))
        for i in range(0, len(Y), 256):
            b = perm[i:i + 256]; r_, a_ = R[b].to(dev), A[b].to(dev)
            loss = (Wt[b].to(dev) * (effect(net, r_, a_) - Y[b].to(dev)) ** 2).sum() / Wt[b].sum().to(dev)
            opt.zero_grad(); loss.backward(); opt.step()
        net.eval()
        with torch.no_grad(): v = float(sum(((effect(net, tr[i:i + 2048], ta[i:i + 2048]) - ty[i:i + 2048]) ** 2 * tw[i:i + 2048]).sum() for i in range(0, len(ty), 2048)) / tw.sum())
        if v < best: best, state = v, {k: x.detach().clone() for k, x in net.state_dict().items()}
        print(held, seed, 'epoch', ep, 'siraj-tuning wMSE %.4f' % v, flush=True)
    net.load_state_dict(state); net.eval(); return net, best

P = {}; log = {}
for held in STUDIES:
    d = data[held]; ho = d['ho']; preds = []
    for seed in (1, 2):
        net, best = train(held, seed)
        with torch.no_grad():
            preds.append(np.concatenate([effect(net, torch.tensor(d['r'][ho][i:i + 2048]).to(dev), torch.tensor(d['a'][ho][i:i + 2048]).to(dev)).cpu().numpy() for i in range(0, ho.sum(), 2048)]))
        log[f'{held}:{seed}'] = best
    P[f'{held}/ids'] = np.array([r['id'] for r in d['rows']])[ho]; P[f'{held}/family'] = np.array(d['fam'])[ho]
    P[f'{held}/N1'] = np.mean(preds, 0); P[f'{held}/seed_corr'] = np.corrcoef(*preds)[0, 1]
    print(held, 'done', 'seed-seed corr %.3f' % P[f'{held}/seed_corr'], flush=True)
    np.savez(DER / 'PREDICTIONS_PARTIAL.npz', **P)
np.savez(DER / 'PREDICTIONS_RAW.npz', **P); json.dump(log, open(DER / 'TRAINING.json', 'w'), indent=1)

"""Sei ref/alt predictions for every eligible variant (forward strand, no reverse-complement averaging); keeps the 40 histone-normalised sequence-class
variant scores and the alt-minus-ref probabilities of a fixed neural/brain profile subset. outcomes are never read here."""
import importlib.util, json, re, sys, time
from pathlib import Path
import numpy as np, torch

HERE = Path(__file__).resolve().parent.parent
D6 = HERE.parent / 'hnp-loso6-20261019/derived'
W = HERE / 'source/weights/model/'
OUT = HERE / 'derived'
STUDIES = ('lee', 'mcafee', 'matoba', 'gse225815', 'gse221384', 'weiss')
NEURAL = re.compile(r'Neural_Progenitor|Neuronal_Progenitor|Neural_Stem|Neuroepithelial|Neuroectoderm|Neurosphere|Fetal_Brain|Germinal_Matrix|Cortical_Neuron|^Neuron|_Neuron|Astrocyte|Cortex_Prefrontal|Hippocampus|Midbrain|SK-N-SH|SH-SY5Y|ReN-VM|Brain_')

targets = [l.strip() for l in open(W / 'target.names')]
sel = np.array([i for i, t in enumerate(targets) if NEURAL.search(t.split(' | ')[0]) and 'Kidney' not in t and 'Bone_Marrow' not in t])
proj, hist = np.load(W / 'projvec_targets.npy'), np.load(W / 'histone_inds.npy')
print('neural profiles', len(sel), flush=True)

spec = importlib.util.spec_from_file_location('sei', W / 'sei.py'); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
net = mod.Sei(); sd = torch.load(W / 'sei.pth', map_location='cpu')
net.load_state_dict({k.replace('module.model.', '').replace('module.', ''): v for k, v in sd.items()}); net.eval()
with torch.no_grad(): net(torch.zeros(1, 4, 4096))  # builds the lazy spline matrix
dev = 'mps'; net = net.to(dev)
for m in net.modules():
    if getattr(m, '_spline_tr', None) is not None: m._spline_tr = m._spline_tr.to(dev)

rows = [(s, r) for s in STUDIES for r in json.load(open(D6 / f'inputs_{s}.json'))]
ctx = np.load(HERE / 'derived/context4k.npy', mmap_mode='r'); assert len(ctx) == len(rows)
LUT = np.full(256, 4, np.uint8)
for i, b in enumerate(b'ACGT'): LUT[b] = i
EYE = np.vstack([np.eye(4, dtype=np.float32), np.zeros((1, 4), np.float32)])

def alleles(i):
    r = rows[i][1]; ref, alt = r['reference'], r['sequence']
    v = next(j for j in range(len(ref)) if ref[j] != alt[j])
    a, b = np.array(ctx[i]), np.array(ctx[i])
    s = 2048 - v; a[s:s + len(ref)] = LUT[np.frombuffer(ref.encode(), np.uint8)]; b[s:s + len(alt)] = LUT[np.frombuffer(alt.encode(), np.uint8)]
    return a, b

def classes(ref, alt):  # sei-framework utils.sc_hnorm_varianteffect
    mid = 0.5 * ref[:, hist].sum(1) + 0.5 * alt[:, hist].sum(1)
    r, a = ref.copy(), alt.copy()
    r[:, hist] *= (mid / ref[:, hist].sum(1))[:, None]; a[:, hist] *= (mid / alt[:, hist].sum(1))[:, None]
    p = lambda x: x @ proj.T / np.linalg.norm(proj, axis=1)
    return (p(a) - p(r))[:, :40]

N, B = len(rows), 32
cls = np.zeros((N, 40), np.float32); neu = np.zeros((N, len(sel)), np.float32)
t0 = time.time()
for s in range(0, N, B):
    idx = range(s, min(s + B, N)); pair = [alleles(i) for i in idx]
    x = torch.from_numpy(EYE[np.stack([p[k] for p in pair for k in (0, 1)])]).permute(0, 2, 1).contiguous().to(dev)
    with torch.no_grad(): y = net(x).float().cpu().numpy().astype(np.float64)
    ref, alt = y[0::2], y[1::2]
    cls[s:s + len(idx)] = classes(ref, alt); neu[s:s + len(idx)] = (alt - ref)[:, sel]
    if (s // B) % 50 == 0: print(s, round(time.time() - t0), 's', flush=True)
np.savez(OUT / 'sei_features.npz', classes=cls, neural=neu, ids=np.array([r['id'] for _, r in rows]), study=np.array([s for s, _ in rows]), profiles=np.array(targets)[sel])
print('done', N, round(time.time() - t0), 's')

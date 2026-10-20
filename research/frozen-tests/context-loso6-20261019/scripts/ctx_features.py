"""motif-context features [d, d*(Q_alt+Q_ref)] (532+532, scale .003) for the six studies' variants.
run with --check to reproduce the cached siraj rows first."""
import json, sys, time
from pathlib import Path
import numpy as np
BUILD = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(BUILD / 'siraj-context-pwm-20261019'))
import context_features as cf
sys.path.insert(0, str(BUILD / 'lee2025-gse244011-20261019/scripts'))
import models as m
SCALE = 0.003
real = [m.motif_features.log_odds(p['counts']) for p in m.pwm.profiles()]

def features(rows):
    d = cf.linear.deltas(rows, real, SCALE)
    return np.hstack([d, cf.quadratic(rows, real, SCALE, d)])

if sys.argv[1:] == ['--check']:
    rows = (m.read_rows('FIT.jsonl.gz') + m.read_rows('TUNING.jsonl.gz'))
    idx = np.r_[0:64, 70000:70064]
    cache = np.load(BUILD / 'siraj-context-pwm-20261019/features/real.npy', mmap_mode='r')[idx]
    mine = features([rows[i] for i in idx])
    print('max abs diff vs reference cache', np.abs(mine - cache).max(), 'cache scale', np.abs(cache).max())
    sys.exit()
DER = Path(__file__).resolve().parent.parent / 'derived'
D6 = BUILD / 'hnp-loso6-20261019/derived'
out = {}
for s in ['lee', 'mcafee', 'matoba', 'gse225815', 'gse221384', 'weiss']:
    rows = json.load(open(D6 / f'inputs_{s}.json')); t = time.time()
    out[s] = np.vstack([features(rows[i:i + 64]) for i in range(0, len(rows), 64)])
    print(s, out[s].shape, round(time.time() - t), 's', flush=True)
np.savez(DER / 'ctx_features.npz', **out)

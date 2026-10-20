"""phyloP100way (hg38) at each credible-set member and its mean over +-10 bases; one shard per process. neighbouring variants share one remote read. no outcome is read."""
import json, sys, time
from pathlib import Path
import numpy as np, pyBigWig
HERE = Path(__file__).resolve().parent.parent
k, n = int(sys.argv[1]), int(sys.argv[2])
var = json.load(open(HERE / 'derived/cs_ctx4k_variants.json'))
order = sorted(range(len(var)), key=lambda i: (var[i][0], var[i][1]))
groups, cur = [], []
for i in order:
    if cur and (var[i][0] != var[cur[0]][0] or var[i][1] - var[cur[0]][1] > 5000): groups.append(cur); cur = []
    cur.append(i)
groups.append(cur)
for attempt in range(8):
    try: bw = pyBigWig.open("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw"); break
    except RuntimeError: time.sleep(5)
pt = np.full(len(var), np.nan); win = np.full(len(var), np.nan); t0 = time.time()
for gi in range(k, len(groups), n):
    g = groups[gi]; c = var[g[0]][0]; lo = var[g[0]][1] - 1 - 10; hi = var[g[-1]][1] - 1 + 11
    try: v = np.array(bw.values(c, lo, hi), float)
    except RuntimeError: continue
    for i in g:
        j = var[i][1] - 1 - 10 - lo; w = v[j:j + 21]; pt[i] = w[10]; win[i] = np.nanmean(w) if np.isfinite(w).any() else np.nan
    if gi % (100 * n) == k: print(gi, len(groups), round(time.time() - t0), 's', flush=True)
np.savez(HERE / f'derived/phylop_cs_{k}of{n}.npz', pt=pt, win=win); print('done', k, flush=True)

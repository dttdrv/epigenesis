"""phyloP100way (hg19) at each variant and the mean over +-10 bp; one study shard per process. neighbouring variants share one remote read. outcomes are never read here."""
import json, sys, time
from pathlib import Path
import numpy as np, pyBigWig
ROOT = Path(__file__).resolve().parent.parent; D6 = ROOT.parent / 'hnp-loso6-20261019/derived'
s, k, n = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
rows = json.load(open(D6 / f'inputs_{s}.json'))
order = sorted(range(len(rows)), key=lambda i: (str(rows[i]['chrom']), int(rows[i]['pos'])))
groups, cur = [], []
for i in order:
    if cur and (str(rows[i]['chrom']) != str(rows[cur[0]]['chrom']) or int(rows[i]['pos']) - int(rows[cur[0]]['pos']) > 5000): groups.append(cur); cur = []
    cur.append(i)
groups.append(cur)
bw = pyBigWig.open("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way.bw")
pt = np.full(len(rows), np.nan); win = np.full(len(rows), np.nan); t0 = time.time()
for gi in range(k, len(groups), n):
    g = groups[gi]; c = 'chr' + str(rows[g[0]]['chrom']); lo = int(rows[g[0]]['pos']) - 11; hi = int(rows[g[-1]]['pos']) + 10
    try: v = np.array(bw.values(c, lo, hi), float)
    except RuntimeError: continue
    for i in g:
        j = int(rows[i]['pos']) - 11 - lo; w = v[j:j + 21]; pt[i] = w[10]; win[i] = np.nanmean(w) if np.isfinite(w).any() else np.nan
    if gi % (50 * n) == k: print(s, gi, len(groups), round(time.time() - t0), 's', flush=True)
np.savez(ROOT / f'derived/phylop_{s}_{k}of{n}.npz', pt=pt, win=win); print('done', s, k, flush=True)

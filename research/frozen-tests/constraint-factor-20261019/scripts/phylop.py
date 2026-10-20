"""phyloP100way (hg19) at each variant and the mean over +-10 bp, for the six studies' variants. outcomes are never read here."""
import json, time
from pathlib import Path
import numpy as np, pyBigWig
D6 = Path(__file__).resolve().parent.parent.parent / 'hnp-loso6-20261019/derived'
bw = pyBigWig.open("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way.bw")
out = {}; t0 = time.time()
for s in ['lee', 'mcafee', 'matoba', 'gse225815', 'gse221384', 'weiss']:
    rows = json.load(open(D6 / f'inputs_{s}.json')); pt = np.full(len(rows), np.nan); win = np.full(len(rows), np.nan)
    for i, r in enumerate(rows):
        c = 'chr' + str(r['chrom']); p = int(r['pos'])
        try:
            v = bw.values(c, p - 11, p + 10)
        except RuntimeError: continue
        v = np.array(v, float); pt[i] = v[10]; win[i] = np.nanmean(v) if np.isfinite(v).any() else np.nan
        if i % 2000 == 0: print(s, i, round(time.time() - t0), 's', flush=True)
    out[s + '/pt'], out[s + '/win'] = pt, win
np.savez(Path(__file__).resolve().parent.parent / 'derived/phylop.npz', **out); print('done')

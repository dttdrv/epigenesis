"""Roadmap DNase peak features for every credible-set member (variant lifted GRCh38 -> GRCh37). no outcome is read."""
import gzip, json
from pathlib import Path
import numpy as np
from pyliftover import LiftOver
HERE = Path(__file__).resolve().parent.parent
var = [tuple(v) for v in json.load(open(HERE / 'derived/cs_ctx4k_variants.json'))]
lo = LiftOver('hg38', 'hg19')
pos19 = []
for c, p, _, _ in var:
    r = lo.convert_coordinate(c, p - 1)
    pos19.append((r[0][0], r[0][1]) if r and r[0][0] == c else (None, None))
print('lifted', sum(x[0] is not None for x in pos19), 'of', len(var), flush=True)
TR = ('E081', 'E082', 'E007')
inpeak = np.zeros((len(var), 3)); sig = np.zeros((len(var), 3)); near = np.full(len(var), 1e5)
for t, e in enumerate(TR):
    peaks = {}
    with gzip.open(HERE / f'source/roadmap/{e}-DNase.narrowPeak.gz', 'rt') as fh:
        for l in fh:
            x = l.split('\t'); peaks.setdefault(x[0], []).append((int(x[1]), int(x[2]), float(x[6])))
    arr = {c: np.array(sorted(v)) for c, v in peaks.items()}
    for i, (c, p) in enumerate(pos19):
        if c is None or c not in arr: continue
        a = arr[c]; j = np.searchsorted(a[:, 0], p, side='right')            # peaks starting at or before p are a[:j]
        best = None
        for k in range(max(0, j - 30), j):                                    # peaks are non-overlapping per chromosome up to rare nesting; scan a short window
            if a[k, 0] <= p < a[k, 1] and (best is None or a[k, 2] > best): best = a[k, 2]
        if best is not None: inpeak[i, t] = 1; sig[i, t] = np.log1p(best)
        d = min([abs(p - a[j - 1, 1]) if j > 0 else 1e9, abs(a[j, 0] - p) if j < len(a) else 1e9]); near[i] = min(near[i], d)
    print(e, 'variants in peak', int(inpeak[:, t].sum()), flush=True)
np.savez(HERE / 'derived/peaks_cs.npz', inpeak=inpeak, sig=sig, near=np.log10(1 + near))

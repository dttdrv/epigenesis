"""4,096-base hg19 contexts (variant at 0-based index 2048) for every eligible variant of the six studies, from Ensembl GRCh37; the study's own allele sequences are written over the centre in build_sei_inputs.py."""
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent.parent
D6 = HERE.parent / 'hnp-loso6-20261019'
sys.path.insert(0, str(D6 / 'scripts'))
from ensembl import fetch
rows = [(s, r) for s in ('lee', 'mcafee', 'matoba', 'gse225815', 'gse221384', 'weiss') for r in json.load(open(D6 / f'derived/inputs_{s}.json'))]
print(len(rows), 'variants', flush=True)
out = HERE / 'derived/context4k.npy'
LUT = np.full(256, 4, np.uint8)
for i, b in enumerate(b'ACGT'): LUT[b] = i
arr = np.full((len(rows), 4096), 4, np.uint8)
for i in range(0, len(rows), 500):
    regs = [f"{r['chrom']}:{r['pos'] - 2048}..{r['pos'] + 2047}:1" for _, r in rows[i:i + 500]]
    for j, s in enumerate(fetch(regs)):
        assert len(s) == 4096, (rows[i + j][1]['id'], len(s))
        arr[i + j] = LUT[np.frombuffer(s.encode(), np.uint8)]
    print(i, flush=True)
np.save(out, arr)
json.dump([[s, r['id']] for s, r in rows], open(HERE / 'derived/context4k_ids.json', 'w'))

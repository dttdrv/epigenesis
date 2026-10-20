"""4,096-base GRCh38 contexts (variant at 0-based index 2048) for the unique eQTL variants, sliced from local UCSC hg38 chromosome files; the table's ref allele must match."""
import gzip, json
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent.parent
rows = json.load(open(HERE / 'derived/eqtl_table.json'))
var = sorted({(r['chrom'], r['pos'], r['ref'], r['alt']) for r in rows})
LUT = np.full(256, 4, np.uint8)
for i, b in enumerate(b'ACGT'): LUT[b] = i; LUT[b + 32] = i
arr = np.full((len(var), 4096), 4, np.uint8); bad = 0
for chrom in sorted({v[0] for v in var}):
    with gzip.open(HERE / f'source/hg38/{chrom}.fa.gz', 'rb') as fh:
        fh.readline(); seq = LUT[np.frombuffer(fh.read().replace(b'\n', b''), np.uint8)]
    for i, (c, p, ref, alt) in enumerate(var):
        if c != chrom: continue
        lo = p - 1 - 2048                      # p is 1-based; variant at index 2048
        if lo < 0 or lo + 4096 > len(seq): continue
        arr[i] = seq[lo:lo + 4096]; bad += arr[i, 2048] != LUT[ord(ref)]
    print(chrom, 'ref mismatches so far', bad, flush=True)
keep = [i for i in range(len(var)) if (arr[i] != 4).all() and arr[i, 2048] == LUT[ord(var[i][2])]]
print('variants', len(var), 'usable (ref matches, full window)', len(keep))
np.save(HERE / 'derived/ctx4k.npy', arr[keep]); json.dump([var[i] for i in keep], open(HERE / 'derived/ctx4k_variants.json', 'w'))

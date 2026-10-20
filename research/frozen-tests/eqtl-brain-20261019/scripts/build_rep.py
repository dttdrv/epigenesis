"""fine-mapped eQTL variants (SNV, PIP >= 0.5, one best row per variant-gene-dataset) with effect sign and strand-aware TSS distance."""
import csv, gzip, json, collections
from pathlib import Path
import numpy as np
SRC = Path(__file__).resolve().parent.parent / 'source'; DER = Path(__file__).resolve().parent.parent / 'derived'
TEST = {'schwartzentruber_sensory_neuron': 'QTD000514', 'young_microglia': 'QTD000559'}
TRAIN = {}
tss = {}
with gzip.open(SRC / 'gencode.v39.gtf.gz', 'rt') as fh:
    for l in fh:
        if l[0] == '#': continue
        x = l.split('\t')
        if x[2] != 'gene': continue
        gid = x[8].split('gene_id "')[1].split('"')[0].split('.')[0]
        tss[gid] = (x[0], int(x[3]) if x[6] == '+' else int(x[4]), x[6])
rows = []
for role, group in (('test', TEST),):
    for name, q in group.items():
        with gzip.open(SRC / f'{q}.credible_sets.tsv.gz', 'rt') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                if float(r['pip']) < 0.1: continue
                c, p, ref, alt = r['variant'].split('_')
                if len(ref) != 1 or len(alt) != 1 or r['gene_id'].split('.')[0] not in tss: continue
                ch, t, strand = tss[r['gene_id'].split('.')[0]]
                if ch != c: continue
                d = (int(p) - t) * (1 if strand == '+' else -1)
                rows.append(dict(role=role, dataset=name, chrom=c, pos=int(p), ref=ref, alt=alt, gene=r['gene_id'], beta=float(r['beta']), pip=float(r['pip']), tss_dist=d, cs=r['cs_id']))
json.dump(rows, open(DER / 'rep_table.json', 'w'))
for role in ('test',):
    c = collections.Counter(r['dataset'] for r in rows if r['role'] == role); print(role, sum(c.values()), dict(c))
print('unique test variants', len({(r['chrom'], r['pos'], r['alt']) for r in rows}))

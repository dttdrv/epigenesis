"""fine-mapped eQTL variants (SNV, PIP >= 0.5, one best row per variant-gene-dataset) with effect sign and strand-aware TSS distance."""
import csv, gzip, json, collections
from pathlib import Path
import numpy as np
SRC = Path(__file__).resolve().parent.parent / 'source'; DER = Path(__file__).resolve().parent.parent / 'derived'
TEST = {'jerber_'+n: q for n, q in {'NB_Naive_D11': 'QTD000630', 'P_FPP_Naive_D30': 'QTD000635', 'U_Neur1_Naive_D30': 'QTD000636', 'FPP_Naive_D11': 'QTD000637', 'FPP_Naive_D30': 'QTD000638', 'P_Sert_Naive_D52': 'QTD000639', 'Astro_Naive_D52': 'QTD000640', 'Epen2_Naive_D52': 'QTD000641', 'FPP_Naive_D52': 'QTD000642', 'Sert_Naive_D52': 'QTD000643', 'P_FPP_Naive_D11': 'QTD000644', 'Epen1_Naive_D30': 'QTD000649', 'P_FPP_Naive_D52': 'QTD000650', 'Sert_Naive_D30': 'QTD000651', 'U_Neur1_Naive_D52': 'QTD000652', 'U_Neur2_Naive_D30': 'QTD000654', 'U_Neur3_Naive_D52': 'QTD000656', 'DA_Naive_D30': 'QTD000657', 'Epen1_Naive_D52': 'QTD000658', 'DA_Naive_D52': 'QTD000659'}.items()}
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
                if float(r['pip']) < 0.5: continue
                c, p, ref, alt = r['variant'].split('_')
                if len(ref) != 1 or len(alt) != 1 or r['gene_id'].split('.')[0] not in tss: continue
                ch, t, strand = tss[r['gene_id'].split('.')[0]]
                if ch != c: continue
                d = (int(p) - t) * (1 if strand == '+' else -1)
                rows.append(dict(role=role, dataset=name, chrom=c, pos=int(p), ref=ref, alt=alt, gene=r['gene_id'], beta=float(r['beta']), pip=float(r['pip']), tss_dist=d, cs=r['cs_id']))
json.dump(rows, open(DER / 'confirm_table.json', 'w'))
for role in ('test',):
    c = collections.Counter(r['dataset'] for r in rows if r['role'] == role); print(role, sum(c.values()), dict(c))
print('unique test variants', len({(r['chrom'], r['pos'], r['alt']) for r in rows}))

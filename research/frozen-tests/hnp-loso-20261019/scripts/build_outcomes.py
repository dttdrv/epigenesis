"""allelic effects (alt minus ref) with odd/even replicate halves per study, aligned to inputs_<study>.json. prints counts only."""
import json
import numpy as np
from common import *

def save(s, E, a, b):
    np.savez(DER / f'outcomes_{s}.npz', E=E, half_a=a, half_b=b)
    print(s, len(E), 'outcomes written', flush=True)

rows = json.load(open(DER / 'inputs_lee.json'))
_, D = table(LEE / 'source/GSE244011_dnaMatrix.csv.gz'); _, R = table(LEE / 'source/GSE244011_rnaMatrix.csv.gz')
d, r = np.stack([D[x['id']] for x in rows]), np.stack([R[x['id']] for x in rows])
cols = lambda sl: (np.arange(22)[sl])
def lee(sl):
    c = cols(sl); return ratio(r[:, c + 22], d[:, c + 22]) - ratio(r[:, c], d[:, c])
save('lee', lee(slice(0, 22)), lee(slice(0, 22, 2)), lee(slice(1, 22, 2)))

rows = json.load(open(DER / 'inputs_mcafee.json'))
_, D = table(MCA / 'source/GSE211045_dna_count_matrix.tsv.gz', '\t'); _, R = table(MCA / 'source/GSE211045_rna_count_matrix.tsv.gz', '\t')
d, r = np.stack([D[x['id']] for x in rows]), np.stack([R[x['id']] for x in rows])
def mca(sl):
    c = np.arange(10)[sl]; return ratio(r[:, c + 10], d[:, c + 10]) - ratio(r[:, c], d[:, c])
save('mcafee', mca(slice(0, 10)), mca(slice(0, 10, 2)), mca(slice(1, 10, 2)))

cont = {}
for l in open(MAT / 'CONTRASTS.jsonl'):
    t = json.loads(l); cont[t['design_key']] = t['conditions']['Vehicle']
v = np.array([cont[x['id']] for x in json.load(open(DER / 'inputs_matoba.json'))])
save('matoba', v.mean(1), v[:, 0], v[:, 1])

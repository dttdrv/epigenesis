"""eligible variants per study with REF/ALT windows, positions, families and held-out overlap flags. outcomes are not read."""
import bisect, collections, json
import numpy as np
from common import *

out = {}

# lee: author exclude list, DNA support, hg19 windows with the labelled alleles placed at index 75
_, dna = table(LEE / 'source/GSE244011_dnaMatrix.csv.gz')
excl = {r[0].replace('_', ':', 1) for r in list(csv.reader(gzip.open(LEE / 'source/GSE244011_excludelist.csv.gz', 'rt')))[1:]}
win = {l.split('\t')[0]: l.split('\t')[1] for l in open(LEE / 'derived/windows_hg19.tsv').read().splitlines()[1:]}
rows = []
for v, d in dna.items():
    if v in excl or d[:22].sum() < 1000 or d[22:].sum() < 1000: continue
    c, p, ref, alt = lee_pos(v); w = win[v]
    assert w[75] in (ref, alt)
    rows.append({'id': v, 'chrom': c, 'pos': p, 'reference': w[:75] + ref + w[76:], 'sequence': w[:75] + alt + w[76:]})
out['lee'] = rows

# mcafee: single-base differences only, DNA support, GRCh37 position of the rsID
seq = dict(csv.reader(gzip.open(MCA / 'source/GSE211045_cre_sequences.tsv.gz', 'rt'), delimiter='\t')); seq.pop('names')
_, dna = table(MCA / 'source/GSE211045_dna_count_matrix.tsv.gz', '\t')
coord = json.load(open(MCA / 'rsid_coordinates_grch37.json'))
rows = []
for rs, d in dna.items():
    ref, alt = seq.get(rs + '_scz_ref'), seq.get(rs + '_scz_alt')
    m = [x for x in coord.get(rs, []) if x['chr'] in STD]
    if not (ref and alt and len(m) == 1 and d[:10].sum() >= 1000 and d[10:].sum() >= 1000): continue
    diff = [i for i in range(len(ref)) if ref[i] != alt[i]]
    if diff == [75]: rows.append({'id': rs, 'chrom': m[0]['chr'], 'pos': m[0]['start'], 'reference': ref, 'sequence': alt})
out['mcafee'] = rows

# matoba: TRAIN pairs only (validation and test stay closed); allele1 is reference, allele2 alternate
pair = {}
for l in open(BUILD / 'matoba-source-admission/PAIR_DESIGNS.jsonl'):
    r = json.loads(l); pair[r['allele1_id'].split('_A1_')[0]] = r
coord = json.load(open(ROOT.parent / 'matoba-hnp-baseline-20261019/rsid_coordinates_grch37.json'))
rows = []
for l in open(MAT / 'CONTRASTS.jsonl'):
    t = json.loads(l); p = pair[t['design_key']]
    m = [x for x in coord.get(p['snp'], []) if x['chr'] in STD]
    assert p['change_position_zero_based'] == 75
    if len(m) != 1: continue  # rsID retired or unmapped in GRCh37
    rows.append({'id': t['design_key'], 'chrom': m[0]['chr'], 'pos': m[0]['start'], 'reference': p['allele1'], 'sequence': p['allele2'], 'family_hash': t['family']})
out['matoba'] = rows

# siraj design positions (all roles)
siraj = collections.defaultdict(list)
with gzip.open(BUILD / 'siraj-design-inventory-20261019/inventory/PHYSICAL_DESIGN.jsonl.gz', 'rt') as fh:
    for line in fh:
        c, p, _, _ = json.loads(line)['variant'].split(':'); siraj[c.removeprefix('chr')].append(int(p))

def near(arr, p):
    j = bisect.bisect_left(arr, p - WINDOW); return j < len(arr) and arr[j] <= p + WINDOW

pos = {s: collections.defaultdict(list) for s in STUDIES}
for s in STUDIES:
    for r in out[s]: pos[s][r['chrom']].append(r['pos'])
for d in list(pos.values()) + [siraj]:
    for a in d.values(): a.sort()
pos['siraj'] = siraj

summary = {}
for s in STUDIES:
    rows = out[s]
    # families: 2,000-base connected components; matoba keeps its genomic-locus hash
    last, fid = {}, 0
    if s == 'matoba':
        ids = {h: i + 1 for i, h in enumerate(sorted({r['family_hash'] for r in rows}))}
        for r in rows: r['family'] = ids[r['family_hash']]; del r['family_hash']
    else:
        for r in sorted(rows, key=lambda r: (r['chrom'], r['pos'])):
            if r['chrom'] not in last or r['pos'] - last[r['chrom']] > WINDOW: fid += 1
            last[r['chrom']] = r['pos']; r['family'] = fid
    others = [pos[o] for o in list(STUDIES) + ['siraj'] if o != s]
    for r in rows: r['heldout'] = not any(near(o.get(r['chrom'], []), r['pos']) for o in others)
    summary[s] = {'eligible': len(rows), 'families': len({r['family'] for r in rows}), 'heldout_eligible': sum(r['heldout'] for r in rows),
                  'heldout_families': len({r['family'] for r in rows if r['heldout']})}
    json.dump(rows, open(DER / f'inputs_{s}.json', 'w'))
print(json.dumps(summary, indent=1))
json.dump(summary, open(DER / 'POPULATION.json', 'w'), indent=1)

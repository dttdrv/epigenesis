"""inputs and outcomes for the three new studies, then unified families and held-out flags for all six."""
import collections, glob, json, shutil
from common import *
from ensembl import fetch

AB = 'ACGT'
def window(chrom, pos):
    return f'{chrom}:{pos - 75}..{pos + 74}:1'

# --- gse225815: schizophrenia variants, iPSC-derived neural progenitors, barcode counts, GRCh37, position in the eid
agg = collections.defaultdict(lambda: collections.defaultdict(float)); meta = {}
for f in sorted(glob.glob(str(BUILD / 'gse225815-scz-20261019/source/*.gz'))):
    with gzip.open(f, 'rt') as fh:
        hdr = next(fh).rstrip('\n').split('\t'); npc = [i for i, h in enumerate(hdr) if h.startswith('NPC_')]
        for l in fh:
            x = l.rstrip('\n').split('\t')
            if x[12] != 'SNP' or x[2] != '0': continue  # single-base variants, exact barcode match
            e = x[1].split('_')
            if len(e) != 3 or not e[1].isdigit() or e[2] not in ('Ref', 'Alt'): continue  # eids of other designs
            c, p, a = e; k = (c[3:], int(p)); meta[k] = (x[10], x[11])
            g = agg[(k, a)]; g['dna'] += float(x[13])
            for i in npc: g['rna:' + hdr[i]] += float(x[i] or 0)
sites = [k for k in meta if (k, 'Ref') in agg and (k, 'Alt') in agg and min(agg[(k, 'Ref')]['dna'], agg[(k, 'Alt')]['dna']) >= 100]
print('gse225815 sites with DNA >= 100 per allele', len(sites), flush=True)
win = fetch([window(c, p) for c, p in sites])
rows, E, A, B = [], [], [], []
for (c, p), w in zip(sites, win):
    ref, alt = meta[(c, p)]
    if len(w) != 150 or w[75] != ref: continue
    reps = sorted({k[4:] for k in agg[((c, p), 'Ref')] if k.startswith('rna:')})
    def eff(sub): return float(ratio(sum(agg[((c, p), 'Alt')]['rna:' + r] for r in sub), agg[((c, p), 'Alt')]['dna']) - ratio(sum(agg[((c, p), 'Ref')]['rna:' + r] for r in sub), agg[((c, p), 'Ref')]['dna']))
    rows.append({'id': f'chr{c}:{p}_{ref}:{alt}', 'chrom': c, 'pos': p, 'reference': w[:75] + ref + w[76:], 'sequence': w[:75] + alt + w[76:]})
    E.append(eff(reps)); A.append(eff(reps[0::2])); B.append(eff(reps[1::2]))
json.dump(rows, open(DER / 'inputs_gse225815.json', 'w')); np.savez(DER / 'outcomes_gse225815.npz', E=E, half_a=A, half_b=B)
print('gse225815 written', len(rows), flush=True)

# --- gse221384: neuropsychiatric variants, anterior and posterior neural progenitors, plasmid DNA, GRCh37
with gzip.open(BUILD / 'gse221384-neuropsych-20261019/source/GSE221384_count_matrix_combined_wdna.csv.gz', 'rt') as fh:
    r = csv.reader(fh); hdr = next(r)
    npc = [i for i, h in enumerate(hdr) if h.startswith(('P-NPC_', 'A-NPC_'))]; pl = hdr.index('plasmid')
    agg = collections.defaultdict(lambda: np.zeros(len(npc) + 1))
    for x in r:
        c, p, al, b, _ = x[1].split('.')
        agg[(c[3:], int(p), 'Ref' if al == 'Ref' else 'alt' + b)] += [float(x[i]) for i in npc] + [float(x[pl])]
refs = {(c, p): a for (c, p, a) in agg if a == 'Ref'}
pairs = [(c, p, a) for (c, p, a) in agg if a != 'Ref' and (c, p, 'Ref') in agg and min(agg[(c, p, a)][-1], agg[(c, p, 'Ref')][-1]) >= 500]
# the ref base at a site is the allele with the Ref tag in the barcode name
refbase = {}
with gzip.open(BUILD / 'gse221384-neuropsych-20261019/source/GSE221384_count_matrix_combined_wdna.csv.gz', 'rt') as fh:
    next(fh)
    for l in fh:
        n = l.split(',', 2)[1].split('.')
        if n[2] == 'Ref': refbase[(n[0][3:], int(n[1]))] = n[3]
sites = sorted({(c, p) for c, p, a in pairs})
print('gse221384 pairs with plasmid DNA >= 500 per allele', len(pairs), 'sites', len(sites), flush=True)
sw = dict(zip(sites, fetch([window(c, p) for c, p in sites])))
rows, E, A, B = [], [], [], []
for c, p, a in pairs:
    w, ref, alt = sw[(c, p)], refbase[(c, p)], a[3:]
    if len(w) != 150 or w[75] != ref: continue
    R, L = agg[(c, p, 'Ref')], agg[(c, p, a)]
    eff = lambda idx: float(ratio(L[idx].sum(), L[-1]) - ratio(R[idx].sum(), R[-1]))
    n = len(npc)
    rows.append({'id': f'chr{c}:{p}_{ref}:{alt}', 'chrom': c, 'pos': p, 'reference': w[:75] + ref + w[76:], 'sequence': w[:75] + alt + w[76:]})
    E.append(eff(list(range(n)))); A.append(eff(list(range(0, n, 2)))); B.append(eff(list(range(1, n, 2))))
json.dump(rows, open(DER / 'inputs_gse221384.json', 'w')); np.savez(DER / 'outcomes_gse221384.npz', E=E, half_a=A, half_b=B)
print('gse221384 written', len(rows), flush=True)

# --- weiss: modern (derived) versus archaic (ancestral) alleles in neural progenitors, GRCh37, 200-base oligos
fa, name = {}, None
with gzip.open(BUILD / 'weiss2021-gse152404-20261019/source/GSE152404_Oligos_library_joint_noDups.fasta.gz', 'rt') as fh:
    for l in fh:
        l = l.strip()
        if l.startswith('>'): name = l[1:]
        else: fa[name] = l.upper()
tab = [l.split('\t') for l in open(BUILD / 'weiss2021-gse152404-20261019/source/ProcessedData.tsv').read().splitlines()]
h = tab[0]; ix = {k: h.index(k) for k in ('NPC_raw_anc_DNA_cpm', 'NPC_raw_der_DNA_cpm', 'NPC_raw_anc_RNA_cpm', 'NPC_raw_der_RNA_cpm', 'chromosome', 'oligo_start_hg19')}
def num(s):
    try: return float(s)
    except ValueError: return float('nan')
rows, E = [], []
for t in tab[1:]:
    oid = t[0]; head, tail = oid.split('_[', 1)
    anc, der = fa.get(f'{head}_Ancestral_[{tail}'), fa.get(f'{head}_Derived_[{tail}')
    v = [num(t[ix[k]]) for k in ('NPC_raw_anc_DNA_cpm', 'NPC_raw_der_DNA_cpm', 'NPC_raw_anc_RNA_cpm', 'NPC_raw_der_RNA_cpm')]
    if not (anc and der and len(anc) == len(der)) or np.isnan(v).any() or min(v[:2]) < 30: continue
    d = [i for i in range(len(anc)) if anc[i] != der[i]]
    if len(d) != 1: continue
    rows.append({'id': oid, 'chrom': t[ix['chromosome']].removeprefix('chr'), 'pos': int(t[ix['oligo_start_hg19']]) + d[0], 'reference': anc, 'sequence': der})
    E.append(float(np.log2((v[3] + 1) / (v[1] + 1)) - np.log2((v[2] + 1) / (v[0] + 1))))
json.dump(rows, open(DER / 'inputs_weiss.json', 'w')); np.savez(DER / 'outcomes_weiss.npz', E=E, half_a=np.full(len(E), np.nan), half_b=np.full(len(E), np.nan))
print('weiss written', len(rows), flush=True)

# --- the three earlier studies: inputs and outcomes copied from the three-study test
for s in OLD:
    shutil.copy(V1 / f'inputs_{s}.json', DER / f'inputs_{s}.json'); shutil.copy(V1 / f'outcomes_{s}.npz', DER / f'outcomes_{s}.npz')

# --- families and held-out flags over all six studies plus siraj designs
siraj = collections.defaultdict(list)
with gzip.open(BUILD / 'siraj-design-inventory-20261019/inventory/PHYSICAL_DESIGN.jsonl.gz', 'rt') as fh:
    for line in fh:
        c, p, _, _ = json.loads(line)['variant'].split(':'); siraj[c.removeprefix('chr')].append(int(p))
data = {s: json.load(open(DER / f'inputs_{s}.json')) for s in STUDIES}
pos = {s: collections.defaultdict(list) for s in STUDIES}
for s in STUDIES:
    for r in data[s]: pos[s][r['chrom']].append(r['pos'])
for d in list(pos.values()) + [siraj]:
    for a in d.values(): a.sort()
pos['siraj'] = siraj
near = lambda arr, p: (lambda j: j < len(arr) and arr[j] <= p + WINDOW)(bisect.bisect_left(arr, p - WINDOW))
summary = {}
for s in STUDIES:
    rows = data[s]
    if s != 'matoba':
        last, fid = {}, 0
        for r in sorted(rows, key=lambda r: (r['chrom'], r['pos'])):
            if r['chrom'] not in last or r['pos'] - last[r['chrom']] > WINDOW: fid += 1
            last[r['chrom']] = r['pos']; r['family'] = fid
    others = [pos[o] for o in list(STUDIES) + ['siraj'] if o != s]
    for r in rows: r['heldout'] = not any(near(o.get(r['chrom'], []), r['pos']) for o in others)
    summary[s] = {'eligible': len(rows), 'families': len({r['family'] for r in rows}), 'heldout_eligible': sum(r['heldout'] for r in rows), 'heldout_families': len({r['family'] for r in rows if r['heldout']})}
    json.dump(rows, open(DER / f'inputs_{s}.json', 'w'))
print(json.dumps(summary, indent=1)); json.dump(summary, open(DER / 'POPULATION.json', 'w'), indent=1)

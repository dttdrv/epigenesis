"""qualifying SuSiE credible sets (2-30 single-base members, unique top PIP >= 0.5) with strand-aware TSS distance. writes derived/cs_table.json (role, dataset, set id, members)."""
import csv, gzip, json, collections
from pathlib import Path
HERE = Path(__file__).resolve().parent.parent; SRC = HERE.parent / 'eqtl-brain-20261019/source'
TRAIN = {'Braineac2_putamen': 'QTD000041', 'Braineac2_sn': 'QTD000046', 'BrainSeq': 'QTD000051', 'CommonMind': 'QTD000075', 'ROSMAP': 'QTD000434',
         **{f'GTEx_{i}': f'QTD000{i}' for i in (146, 151, 156, 161, 166, 171, 176, 181, 186, 191, 196, 201, 206)}}
JERBER = {'NB_Naive_D11': 'QTD000630', 'P_FPP_Naive_D30': 'QTD000635', 'U_Neur1_Naive_D30': 'QTD000636', 'FPP_Naive_D11': 'QTD000637', 'FPP_Naive_D30': 'QTD000638', 'P_Sert_Naive_D52': 'QTD000639', 'Astro_Naive_D52': 'QTD000640', 'Epen2_Naive_D52': 'QTD000641', 'FPP_Naive_D52': 'QTD000642', 'Sert_Naive_D52': 'QTD000643', 'P_FPP_Naive_D11': 'QTD000644', 'Epen1_Naive_D30': 'QTD000649', 'P_FPP_Naive_D52': 'QTD000650', 'Sert_Naive_D30': 'QTD000651', 'U_Neur1_Naive_D52': 'QTD000652', 'U_Neur2_Naive_D30': 'QTD000654', 'U_Neur3_Naive_D52': 'QTD000656', 'DA_Naive_D30': 'QTD000657', 'Epen1_Naive_D52': 'QTD000658', 'DA_Naive_D52': 'QTD000659'}
GROUPS = {'train': TRAIN, 'walker': {'walker': 'QTD000579'}, 'secondary': {'aygun_progenitor': 'QTD000564', 'aygun_neuron': 'QTD000569', 'sensory': 'QTD000514', 'microglia': 'QTD000559', **{'jerber_' + k: v for k, v in JERBER.items()}}}
tss = {}
with gzip.open(SRC / 'gencode.v39.gtf.gz', 'rt') as fh:
    for l in fh:
        if l[0] == '#': continue
        x = l.split('\t')
        if x[2] != 'gene': continue
        tss[x[8].split('gene_id "')[1].split('"')[0].split('.')[0]] = (x[0], int(x[3]) if x[6] == '+' else int(x[4]), x[6])
out = {}
for role, group in GROUPS.items():
    best = {}
    for name, q in group.items():
        sets = collections.defaultdict(list)
        with gzip.open(SRC / f'{q}.credible_sets.tsv.gz', 'rt') as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                g = r['gene_id'].split('.')[0]
                if g not in tss: continue
                c, p, ref, alt = r['variant'].split('_'); ch, t, strand = tss[g]
                if ch != c: continue
                snv = len(ref) == 1 and len(alt) == 1
                sets[(r['molecular_trait_id'], r['cs_id'])].append(dict(chrom=c, pos=int(p), ref=ref, alt=alt, pip=float(r['pip']), d=(int(p) - t) * (1 if strand == '+' else -1), snv=snv, gene=g))
        for key, m in sets.items():
            top = sorted(m, key=lambda x: -x['pip']); mm = [x for x in m if x['snv']]
            if not (top[0]['snv'] and top[0]['pip'] >= 0.5 and (len(top) == 1 or top[1]['pip'] < top[0]['pip']) and 2 <= len(mm) <= 30): continue
            k = (top[0]['gene'], top[0]['chrom'], top[0]['pos'], top[0]['alt'])
            if k not in best or top[0]['pip'] > best[k][0]: best[k] = (top[0]['pip'], name, key, mm)
    out[role] = [dict(dataset=n, cs=list(key), members=mm) for _, n, key, mm in best.values()]
    print(role, 'sets', len(out[role]), 'members', sum(len(s['members']) for s in out[role]))
# a variant set shared between roles never crosses the chromosome rule, but report overlap for the record
json.dump(out, open(HERE / 'derived/cs_table.json', 'w'))

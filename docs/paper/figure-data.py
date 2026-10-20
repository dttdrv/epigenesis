"""Rebuilds the predictor figure data of the technical paper and prints the numbers quoted in its section 5.

Run with the Studio interpreter, from any folder:
    $(uv tool dir)/epigenesis-studio/bin/python docs/paper/figure-data.py

figure-scan.csv     Studio Scan of the 200-base example sequence (600 edits); `nan` marks the reference base.
figure-studies.csv  per-study results copied from research/frozen-tests/huber-transfer-20261020/RESULTS.md.
"""
import copy
import csv
import json
import re
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT/'studio')]
from epigenesis_studio import server
from brainc import _canonical as canonical

HERE = Path(__file__).resolve().parent
EDITS = json.loads((ROOT/'studio/epigenesis_studio/example-edits.json').read_text())
SEQUENCE = EDITS[0]['reference']

# Scan
predictions = server.run(server.normalize(dict(mode='scan', sequence=SEQUENCE)))['predictions']
effect = {p['id']: p['effect'] for p in predictions}
with open(HERE/'figure-scan.csv', 'w', newline='') as f:
    out = csv.writer(f)
    out.writerow(['pos', 'altidx', 'alt', 'effect'])
    for j, alt in enumerate('ACGT'):
        for i, base in enumerate(SEQUENCE):
            value = effect.get(f'{i+1}:{base}>{alt}')
            out.writerow([i+1, j, alt, 'nan' if value is None else f'{value:.6f}'])
best, worst = max(predictions, key=lambda p: p['effect']), min(predictions, key=lambda p: p['effect'])
print(f'scan: {len(predictions)} edits, {sum(p["effect"] > 0 for p in predictions)} up, '
      f'{sum(p["effect"] < 0 for p in predictions)} down, median |effect| {statistics.median(abs(p["effect"]) for p in predictions):.4f}')
print(f'  largest increase {best["id"]} {best["effect"]:+.4f}; largest decrease {worst["id"]} {worst["effect"]:+.4f}')
first = EDITS[0]
changed = next(i for i, (a, b) in enumerate(zip(first['reference'], first['sequence'])) if a != b)
pair = server.run(server.normalize(dict(mode='pairs', rows=[{k: first[k] for k in ('id', 'reference', 'sequence')}], budget=None)))
print(f'  example edit {first["id"]} at position {changed+1}: Edits {pair["predictions"][0]["effect"]:+.6f}, '
      f'Scan {effect[f"{changed+1}:{first["reference"][changed]}>{first["sequence"][changed]}"]:+.6f}')

# Edits panel and receipt
BUDGET = 5
receipt = server.score(dict(mode='pairs', rows=EDITS, budget=BUDGET))
ranked = sorted(receipt['result']['predictions'], key=lambda p: (-abs(p['effect']), p['id']))
print(f'edits: {len(EDITS)} edits, {len({r["family"] for r in EDITS})} families; receipt {receipt["receipt_sha256"]}')
print('  strongest five by |effect|:', [(p['id'].split(':', 1)[1], round(p['effect'], 4)) for p in ranked[:5]])
for s in receipt['result']['selection']:
    print(f'  panel {s["id"]} {s["family"]} {s["effect"]:+.4f} (words {s["word_effect"]:+.4f}, motifs {s["motif_effect"]:+.4f})')


def outcome(r):
    verdict = server.verify(r)
    return [c['ok'] for c in verdict['checks']], verdict['max_abs_difference']


def reseal(r):
    r = copy.deepcopy(r)
    r['receipt_sha256'] = canonical.digest({k: v for k, v in r.items() if k != 'receipt_sha256'})
    return r


edited = copy.deepcopy(receipt)
edited['result']['selection'][0]['effect'] += 1e-3
renamed = copy.deepcopy(receipt)
renamed['model'] = 'other-model'
print('verify (digest, identity, replay):')
for name, r in [('original', receipt), ('result edited', edited), ('result edited, digest recomputed', reseal(edited)),
                ('model renamed, digest recomputed', reseal(renamed))]:
    print(f'  {name}: {outcome(r)}')

# per-study table
text = (ROOT/'research/frozen-tests/huber-transfer-20261020/RESULTS.md').read_text()
number = r'(-?[\d.]+)'
rows = {m[0]: m[1:] for m in re.findall(rf'^\| (\w+) \| ([\d,]+) \| ([\d,]+) \| {number} \| {number} \| {number} \| {number} \| {number} \|$', text, re.M)}
zero = {m[0]: m[1:] for m in re.findall(rf'^\| (\w+) \| {number} \| {number} \| {number}% \|$', text, re.M)}
with open(HERE/'figure-studies.csv', 'w', newline='') as f:
    out = csv.writer(f)
    out.writerow(['study', 'primary', 'edits', 'families', 'huber', 'jl', 'k5', 's2', 'zero_mse', 'huber_mse', 'excess'])
    for study in sorted(rows, key=lambda s: s == 'gse225815'):  # primary five first
        edits, families, huber, jl, k5, s2, _s3 = rows[study]
        out.writerow([study, int(study != 'gse225815'), edits.replace(',', ''), families.replace(',', ''), huber, jl, k5, s2, *zero[study]])
primary = [s for s in rows if s != 'gse225815']
for name, column in (('Huber', 2), ('K5', 4), ('S2', 5)):
    print(f'mean r, primary five, {name}: {statistics.mean(float(rows[s][column]) for s in primary):.6f}')
print('mean MSE, primary five, zero / Huber:', *(f'{statistics.mean(float(zero[s][i]) for s in primary):.8f}' for i in (0, 1)))

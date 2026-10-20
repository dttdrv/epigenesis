"""Apply the frozen Siraj Huber-linear model to native DNA pairs; no training."""
import argparse
import hashlib
import json
import os
from pathlib import Path

# the lecture command must not occupy every CPU core by default.
for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[name] = '1'

import numpy as np
import features

HERE = Path(__file__).resolve().parent


def load():
    manifest = json.loads((HERE/'model.json').read_text())
    if manifest['format'] != 'epigenesis-siraj-huber-linear-v1':
        raise ValueError('unsupported model format')
    for name, expected in manifest['sha256'].items():
        if Path(name).name != name or hashlib.sha256((HERE/name).read_bytes()).hexdigest() != expected:
            raise ValueError('artifact identity mismatch: '+name)
    with np.load(HERE/'weights.npz', allow_pickle=False) as saved:
        arrays = dict(saved)
    n = len(manifest['motif_ids'])
    shapes = {'word': (len(features.vocabulary(manifest['k'])),), 'motif': (n,)}
    shapes.update({f'motif_{i}': (4, width) for i, width in enumerate(manifest['motif_widths'])})
    if (len(manifest['motif_widths']) != n or len(set(manifest['motif_ids'])) != n
            or set(arrays) != set(shapes) or any(a.shape != shapes[key]
            or a.dtype != np.float64 or not np.isfinite(a).all() for key, a in arrays.items())):
        raise ValueError('invalid model arrays')
    if type(manifest['batch']) is not int or manifest['batch'] <= 0:
        raise ValueError('invalid model batch')
    return manifest, arrays


def validate(rows, k):
    if not isinstance(rows, list) or not rows:
        raise ValueError('input must be a nonempty JSON array')
    ids = set()
    for row in rows:
        if (not isinstance(row, dict) or set(row) not in ({'id', 'reference', 'sequence'}, {'id', 'reference', 'sequence', 'family'})
                or not isinstance(row.get('id'), str) or not row['id'].strip() or row['id'] in ids):
            raise ValueError('each edit needs a unique ID, reference and sequence; optional family only')
        ids.add(row['id'])
        for key in ('reference', 'sequence'):
            sequence = row[key]
            if not isinstance(sequence, str) or len(sequence) < k or set(sequence)-set('ACGT'):
                raise ValueError(f'{row["id"]}: {key} must contain at least {k} uppercase ACGT bases')
        if 'family' in row and (not isinstance(row['family'], str) or not row['family'].strip()):
            raise ValueError('family must be a nonempty string')


def predict(rows):
    manifest, arrays = load()
    validate(rows, manifest['k'])
    matrices = [arrays[f'motif_{i}'] for i in range(len(manifest['motif_ids']))]
    result = []
    for start in range(0, len(rows), manifest['batch']):
        batch = rows[start:start+manifest['batch']]
        delta = features.deltas(batch, matrices, manifest['scale'])
        words = features.whole_counts([r['sequence'] for r in batch], manifest['k']) - features.whole_counts([r['reference'] for r in batch], manifest['k'])
        word_effect = words @ arrays['word']
        motif_effect = delta @ arrays['motif']
        predicted = word_effect + motif_effect
        if not np.isfinite(predicted).all():
            raise ValueError('nonfinite predictions')
        for row, b, d, y in zip(batch, word_effect, motif_effect, predicted, strict=True):
            result.append(dict(id=row['id'], effect=float(y), word_effect=float(b), motif_effect=float(d)))
    return result


def select(predictions, rows, budget):
    families = {row['id']: row.get('family') for row in rows}
    ids = [row['id'] for row in predictions]
    if (not predictions or len(set(ids)) != len(ids) or len(families) != len(rows)
            or set(ids) != set(families) or any(not np.isfinite(row['effect']) for row in predictions)):
        raise ValueError('finite predictions must match every unique input ID exactly')
    if (type(budget) is not int or budget <= 0 or any(not value for value in families.values())
            or len(set(families.values())) < budget):
        raise ValueError('selection needs explicit families and enough distinct families for the positive budget')
    result, seen = [], set()
    for row in sorted(predictions, key=lambda r: (-abs(r['effect']), r['id'])):
        family = families[row['id']]
        if family not in seen:
            seen.add(family)
            result.append(dict(row, family=family, direction=int(np.sign(row['effect']))))
            if len(result) == budget:
                break
    assert len(result) == budget
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--budget', type=int, help='optional panel size; requires a family for every input edit')
    args = parser.parse_args()
    try:
        if args.output.exists():
            raise ValueError('output already exists')
        rows = json.loads(args.input.read_text())
        predictions = predict(rows)
        manifest, _ = load()
        result = dict(model=manifest['format'], status=manifest['status'], endpoint=manifest['endpoint'],
                      manifest_sha256=hashlib.sha256((HERE/'model.json').read_bytes()).hexdigest(), predictions=predictions)
        if args.budget is not None:
            result['selection'] = select(predictions, rows, args.budget)
        with args.output.open('x', encoding='utf-8') as stream:
            stream.write(json.dumps(result, indent=2, allow_nan=False)+'\n')
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        parser.exit(1, str(error)+'\n')


if __name__ == '__main__':
    main()

"""Score neuronal reporter edits with the frozen linear activity-transfer model."""

import argparse
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

import predictor


ROOT = Path(__file__).resolve().parent


def effects(rows, arrays, k, batch):
    words = [''.join(bases) for bases in itertools.product('ACGT', repeat=k)]
    columns = {word: i for i, word in enumerate(words)}
    reverse = np.array([columns[predictor.reverse_complement(word)] for word in words])
    canonical = lambda s: min(s, predictor.reverse_complement(s))
    sequences = sorted({canonical(row[key]) for row in rows for key in ('reference', 'sequence')})
    indices = {sequence: i for i, sequence in enumerate(sequences)}
    counts = np.zeros((len(sequences), len(words)), dtype=np.min_scalar_type(max(map(len,sequences))-k+1))
    for i, sequence in enumerate(sequences):
        for word, count in Counter(sequence[j:j+k] for j in range(len(sequence)-k+1)).items():
            counts[i, columns[word]] = count
    reference = np.array([indices[canonical(row['reference'])] for row in rows])
    alternate = np.array([indices[canonical(row['sequence'])] for row in rows])
    beta = arrays['activity_coefficients']
    strands = np.column_stack((beta, beta[reverse]))
    feature = np.empty((len(rows), 1), dtype=np.float64)
    for start in range(0, len(rows), batch):
        end = start+batch
        delta = counts[alternate[start:end]].astype(np.float64)-counts[reference[start:end]].astype(np.float64)
        feature[start:end, 0] = (delta@strands).mean(axis=1)
    baseline = predictor.features(rows, k)@arrays['baseline_coefficients']
    correction = feature@arrays['head']
    return baseline, correction


def predict(rows, model_dir=ROOT/'data', expected_manifest_sha256=None):
    model_dir = Path(model_dir)
    raw = (model_dir/'activity-manifest.json').read_bytes()
    if expected_manifest_sha256 is not None and hashlib.sha256(raw).hexdigest() != expected_manifest_sha256:
        raise ValueError('model manifest identity mismatch')
    manifest = json.loads(raw)
    if manifest['format'] != 'epigenesis-activity-transfer-v1' or manifest['selected_arm'] != 'linear':
        raise ValueError('unsupported activity-transfer model')
    for name in ('activity.py', 'predictor.py', 'model.py'):
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != manifest['source_sha256'][name]:
            raise ValueError('inference source identity mismatch: '+name)
    length, k, batch = (manifest[key] for key in ('insert_bases', 'k', 'prediction_batch'))
    if any(type(value) is not int for value in (length,k,batch)) or not 0 < k <= length or batch <= 0:
        raise ValueError('invalid model dimensions')
    predictor.validate_rows(rows, length)
    arrays = predictor.checked_arrays(model_dir/'activity.npz', manifest['weights_sha256'])
    shapes = dict(activity_coefficients=(4**k,), baseline_coefficients=(len(predictor.vocabulary(k)),), head=(1,))
    if set(arrays) != set(shapes) or any(value.dtype != np.float64 or value.shape != shapes[key]
            or not np.isfinite(value).all() for key,value in arrays.items()):
        raise ValueError('invalid activity-transfer weights')
    baseline, correction = effects(rows, arrays, k, batch)
    predicted = baseline+correction
    if not np.isfinite(predicted).all():
        raise ValueError('nonfinite activity-transfer predictions')
    return [dict(id=row['id'], effect=float(predicted[i]), baseline=float(baseline[i]),
                 activity_correction=float(correction[i])) for i,row in enumerate(rows)]


def select(predictions, families, budget):
    ids = [row['id'] for row in predictions]
    if (not predictions or len(set(ids)) != len(ids) or not isinstance(families,dict)
            or set(families) != set(ids) or any(not isinstance(v,str) or not v.strip() for v in families.values())):
        raise ValueError('families must map every unique prediction ID to a nonempty component ID')
    if type(budget) is not int or budget <= 0 or len(set(families.values())) < budget:
        raise ValueError('selection requires at least the requested number of distinct families')
    if any(not np.isfinite(row['effect']) for row in predictions):
        raise ValueError('finite effects required for selection')
    selected, seen = [], set()
    for row in sorted(predictions, key=lambda row:(-abs(row['effect']),row['id'])):
        family = families[row['id']]
        if family not in seen:
            seen.add(family)
            selected.append(dict(row, family=family, direction=int(np.sign(row['effect']))))
            if len(selected) == budget:
                break
    return selected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--families', type=Path, help='optional JSON mapping of edit IDs to independent sequence components')
    args = parser.parse_args()
    try:
        if args.output.exists():
            raise ValueError('output already exists')
        raw = (ROOT/'data/activity-manifest.json').read_bytes()
        manifest = json.loads(raw)
        identity = hashlib.sha256(raw).hexdigest()
        rows = json.loads(args.input.read_bytes())
        predictions = predict(rows, expected_manifest_sha256=identity)
        result = dict(status=manifest['status'], endpoint=manifest['endpoint'], context=manifest['context'],
                      manifest_sha256=identity, predictions=predictions)
        if args.families:
            result['selection'] = select(predictions,json.loads(args.families.read_bytes()),manifest['ranking_budget'])
        with args.output.open('x',encoding='utf-8') as stream:
            stream.write(json.dumps(result,indent=2,allow_nan=False)+'\n')
    except (OSError,ValueError,KeyError,RuntimeError) as error:
        parser.exit(1,f'{error}\n')


if __name__ == '__main__':
    main()

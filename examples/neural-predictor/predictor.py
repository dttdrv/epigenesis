"""Predict reporter effects of 270-base DNA pairs with the frozen experimental model."""

import argparse
from collections import Counter
import hashlib
import io
import itertools
import json
from pathlib import Path

import numpy as np
import torch
from torch.nn import functional as F

from model import RegulatoryAttention


ROOT = Path(__file__).resolve().parent


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def vocabulary(k):
    return sorted({min(word, reverse_complement(word))
                   for letters in itertools.product('ACGT', repeat=k)
                   for word in [''.join(letters)]})


def features(rows, k):
    columns = {word: index for index, word in enumerate(vocabulary(k))}
    result = np.zeros((len(rows), len(columns)), dtype=np.float64)
    for index, row in enumerate(rows):
        for key, sign in (('sequence', 1), ('reference', -1)):
            sequence = row[key]
            counts = Counter(min(sequence[i:i+k], reverse_complement(sequence[i:i+k]))
                             for i in range(len(sequence)-k+1))
            for word, count in counts.items():
                result[index, columns[word]] += sign*count
    return result


def validate_rows(rows, length):
    if not isinstance(rows, list) or not rows:
        raise ValueError('input must be a nonempty JSON array of sequence pairs')
    ids = set()
    for row in rows:
        if not isinstance(row, dict) or set(row) != {'id', 'reference', 'sequence'}:
            raise ValueError('each pair must contain only id, reference and sequence')
        if not isinstance(row['id'], str) or not row['id'].strip() or row['id'] in ids:
            raise ValueError('pair identifiers must be unique nonempty strings')
        ids.add(row['id'])
        for key in ('reference', 'sequence'):
            sequence = row[key]
            if not isinstance(sequence, str) or len(sequence) != length or set(sequence)-set('ACGT'):
                raise ValueError(f'{row["id"]}: expected {length} uppercase ACGT bases')
        if sum(a != b for a, b in zip(row['reference'], row['sequence'])) > 1:
            raise ValueError(f'{row["id"]}: only single substitutions or identical controls are supported')


def checked_arrays(path, expected):
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != expected:
        raise ValueError(f'{path.name}: weight identity mismatch')
    with np.load(io.BytesIO(raw), allow_pickle=False) as archive:
        return {key: archive[key].copy() for key in archive.files}


def pair_corrections(model, rows, device, batch_size, length):
    canonical = {s: min(s, reverse_complement(s))
                 for row in rows for s in (row['reference'], row['sequence'])}
    sequences = sorted(set(canonical.values()))
    codes = {base: index for index, base in enumerate('AGCT')}
    encoded = torch.tensor([[codes[base] for base in sequence] for sequence in sequences],
                           dtype=torch.int64)
    if encoded.shape[1] != length:
        raise ValueError('incorrect sequence length')
    model.eval()
    values = []
    with torch.inference_mode():
        for batch in encoded.split(batch_size):
            size = len(batch)
            if size < batch_size:
                batch = torch.cat((batch, batch[-1:].expand(batch_size-size, -1)))
            forward = model(F.one_hot(batch, 4).permute(0, 2, 1).to(device=device, dtype=torch.float32)).cpu().double()
            reverse = model(F.one_hot(3-batch.flip(-1), 4).permute(0, 2, 1).to(device=device, dtype=torch.float32)).cpu().double()
            if forward.shape != (batch_size,) or reverse.shape != (batch_size,):
                raise ValueError('invalid sequence score shape')
            averaged = (forward+reverse)/2
            if not torch.isfinite(averaged).all():
                raise ValueError('nonfinite sequence score')
            values.extend(averaged[:size].tolist())
    scores = dict(zip(sequences, values))
    return np.array([scores[canonical[row['sequence']]]-scores[canonical[row['reference']]]
                     for row in rows], dtype=np.float64)


def predict(rows, device='cpu', model_dir=ROOT/'data', expected_manifest_sha256=None):
    model_dir = Path(model_dir)
    manifest_raw = (model_dir/'manifest.json').read_bytes()
    if expected_manifest_sha256 is not None and hashlib.sha256(manifest_raw).hexdigest() != expected_manifest_sha256:
        raise ValueError('model manifest identity mismatch')
    manifest = json.loads(manifest_raw)
    if manifest['format'] != 'epigenesis-neuronal-predictor-v1' or manifest['alphabet'] != 'AGCT':
        raise ValueError('unsupported predictor format')
    if hashlib.sha256((ROOT/'model.py').read_bytes()).hexdigest() != manifest['model_sha256']:
        raise ValueError('model source identity mismatch')
    architecture = manifest['architecture']
    validate_rows(rows, architecture['length'])
    if device not in ('cpu', 'mps'):
        raise ValueError('supported devices are cpu and mps')
    if device == 'mps' and not torch.backends.mps.is_available():
        raise ValueError('MPS is not available; use --device cpu')
    batch_size = manifest['prediction_batch']
    if type(batch_size) is not int or batch_size <= 0:
        raise ValueError('invalid prediction batch')
    coefficients = np.asarray(manifest['mixture']['coefficients'], dtype=np.float64)
    if (coefficients.shape != (2,) or not np.isfinite(coefficients).all()
            or np.any(coefficients < 0) or coefficients.sum() > 1
            or manifest['mixture']['selection_data'] != 'TRAIN excluded-family predictions only'):
        raise ValueError('invalid frozen mixture')
    baseline_spec = manifest['baseline']
    baseline_arrays = checked_arrays(model_dir/'baseline.npz', baseline_spec['sha256'])
    if set(baseline_arrays) != {'coefficients'}:
        raise ValueError('unexpected baseline arrays')
    prior = baseline_arrays['coefficients']
    if (prior.dtype != np.float64 or prior.shape != (len(vocabulary(baseline_spec['k'])),)
            or not np.isfinite(prior).all()):
        raise ValueError('invalid baseline coefficients')
    baseline = features(rows, baseline_spec['k']) @ prior
    corrections = []
    for name in ('long', 'short'):
        with torch.random.fork_rng(devices=[]):
            model = RegulatoryAttention(**architecture)
        arrays = checked_arrays(model_dir/f'{name}.npz', manifest['weights'][name]['sha256'])
        expected = model.state_dict()
        if set(arrays) != set(expected) or any(
                array.dtype != np.float32 or array.shape != tuple(expected[key].shape)
                or not np.isfinite(array).all() for key, array in arrays.items()):
            raise ValueError(f'{name}: invalid model tensors')
        model.load_state_dict({key: torch.from_numpy(array) for key, array in arrays.items()})
        model.to(device)
        corrections.append(pair_corrections(model, rows, device, batch_size, architecture['length']))
    long, short = corrections
    effect = baseline+coefficients[0]*long+coefficients[1]*short
    if not np.isfinite(effect).all():
        raise ValueError('nonfinite predictions')
    return [dict(id=row['id'], effect=float(effect[i]), baseline=float(baseline[i]),
                 long_correction=float(long[i]), short_correction=float(short[i]))
            for i, row in enumerate(rows)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--device', choices=('cpu', 'mps'), default='cpu')
    args = parser.parse_args()
    try:
        if args.output.exists():
            raise ValueError('output already exists')
        torch.set_num_threads(1)
        rows = json.loads(args.input.read_bytes())
        manifest_sha256 = hashlib.sha256((ROOT/'data/manifest.json').read_bytes()).hexdigest()
        predictions = predict(rows, args.device, expected_manifest_sha256=manifest_sha256)
        result = dict(status='experimental; joint biological validation gate unmet',
                      endpoint='ALT minus REF log2 reporter activity coefficient',
                      manifest_sha256=manifest_sha256,
                      device=args.device, torch=str(torch.__version__), numpy=str(np.__version__),
                      predictions=predictions)
        with args.output.open('x', encoding='utf-8') as stream:
            stream.write(json.dumps(result, indent=2, allow_nan=False)+'\n')
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        parser.exit(1, f'{error}\n')


if __name__ == '__main__':
    main()

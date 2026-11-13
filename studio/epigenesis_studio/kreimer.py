"""Bundled Kreimer H1 neural-induction 0-72h predictor for Epigenesis Studio 1.7.

Executes the deterministic motif-temporal ridge checkpoint fitted on Kreimer et al.
(2022) GSE188264 train+validation families using 5-mer counts + 532 JASPAR motifs.
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if (REPO / 'brainc').is_dir() and str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from brainc import _canonical as canonical

MODEL_DIR = Path(__file__).resolve().parent / 'model'
if str(MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(MODEL_DIR))

import features  # noqa: E402
import predict as siraj_predict  # noqa: E402

MANIFEST_PATH = MODEL_DIR / 'kreimer_model.json'
MANIFEST_BYTES = MANIFEST_PATH.read_bytes()
MANIFEST = json.loads(MANIFEST_BYTES)
MANIFEST_SHA256 = canonical.bytes_digest(MANIFEST_BYTES)
IDENTITY = dict(
    model=MANIFEST['format'],
    manifest_sha256=MANIFEST_SHA256,
    endpoint=MANIFEST['endpoint'],
)

_CACHE = None


def load():
    global _CACHE
    if _CACHE is not None:
        return _CACHE
    if MANIFEST.get('format') != 'epigenesis-kreimer-motif-temporal-v1':
        raise ValueError('unsupported Kreimer model format')
    for name, expected in MANIFEST['sha256'].items():
        if Path(name).name != name or hashlib.sha256((MODEL_DIR / name).read_bytes()).hexdigest() != expected:
            raise ValueError('Kreimer artifact identity mismatch: ' + name)
    siraj_bytes = (MODEL_DIR / 'model.json').read_bytes()
    weights_bytes = (MODEL_DIR / 'weights.npz').read_bytes()
    if hashlib.sha256(siraj_bytes).hexdigest() != MANIFEST['motif_manifest_sha256']:
        raise ValueError('shared motif manifest identity mismatch')
    if hashlib.sha256(weights_bytes).hexdigest() != MANIFEST['motif_weights_sha256']:
        raise ValueError('shared motif weights identity mismatch')
    with np.load(MODEL_DIR / 'kreimer_weights.npz', allow_pickle=False) as ks:
        if 'coefficients' not in ks:
            raise ValueError('missing coefficients array')
        coefs = np.array(ks['coefficients'], dtype=np.float64)
    n_k5 = len(features.vocabulary(MANIFEST['k']))
    n_motifs = len(MANIFEST['motif_ids'])
    n_times = len(MANIFEST['times_hours'])
    if coefs.shape != (n_k5 + n_motifs, n_times) or not np.isfinite(coefs).all():
        raise ValueError('invalid Kreimer coefficient shape or values')
    with np.load(MODEL_DIR / 'weights.npz', allow_pickle=False) as ms:
        matrices = [np.array(ms[f'motif_{i}'], dtype=np.float64) for i in range(n_motifs)]
    _CACHE = (MANIFEST, coefs, matrices, n_k5)
    return _CACHE


def predict(rows):
    manifest, coefs, matrices, n_k5 = load()
    siraj_predict.validate(rows, manifest['k'])
    w_coefs = coefs[:n_k5, :]
    m_coefs = coefs[n_k5:, :]
    motif_ids = manifest['motif_ids']
    motif_mean_coefs = m_coefs.mean(axis=1)
    motif_h72_coefs = m_coefs[:, -1]
    out = []
    for start in range(0, len(rows), manifest['batch']):
        batch = rows[start : start + manifest['batch']]
        delta = features.deltas(batch, matrices, manifest['scale'])
        words = (
            features.whole_counts([r['sequence'] for r in batch], manifest['k'])
            - features.whole_counts([r['reference'] for r in batch], manifest['k'])
        )
        word_mat = words @ w_coefs
        motif_mat = delta @ m_coefs
        total_mat = word_mat + motif_mat
        if not np.isfinite(total_mat).all():
            raise ValueError('nonfinite Kreimer predictions')
        for idx, (row, wm, mm, tm) in enumerate(zip(batch, word_mat, motif_mat, total_mat, strict=True)):
            contribs = delta[idx] * motif_mean_coefs
            nonzero = [j for j in range(len(motif_ids)) if contribs[j] != 0.0]
            nonzero.sort(key=lambda j: (-abs(float(contribs[j])), motif_ids[j]))
            top_motifs = [
                dict(
                    id=motif_ids[j],
                    delta=float(delta[idx, j]),
                    effect=float(contribs[j]),
                    h72_effect=float(delta[idx, j] * motif_h72_coefs[j]),
                )
                for j in nonzero[:5]
            ]
            out.append(dict(
                id=row['id'],
                effect=float(tm.mean()),
                word_effect=float(wm.mean()),
                motif_effect=float(mm.mean()),
                temporal=[float(v) for v in tm],
                top_motifs=top_motifs,
            ))
    return out


def select(predictions, rows, budget):
    return siraj_predict.select(predictions, rows, budget)


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
        manifest, *_ = load()
        result = dict(
            model=manifest['format'],
            status=manifest['status'],
            endpoint=manifest['endpoint'],
            manifest_sha256=MANIFEST_SHA256,
            times_hours=manifest['times_hours'],
            predictions=predictions,
        )
        if args.budget is not None:
            result['selection'] = select(predictions, rows, args.budget)
        with args.output.open('x', encoding='utf-8') as stream:
            stream.write(json.dumps(result, indent=2, allow_nan=False) + '\n')
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        parser.exit(1, str(error) + '\n')


if __name__ == '__main__':
    main()

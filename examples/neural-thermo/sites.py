"""Per-site PWM occupancy sums for the thermodynamic enhancer model.

For a sequence s and motif m with log-odds matrix M (4 x width), each window i
has strand-averaged odds o_i = (exp(f_i) + exp(r_i)) / 2, where f and r are the
forward and reverse-complement scores (the convention of studio features.py).
The saturating site occupancy is theta_i = kappa * o_i / (1 + kappa * o_i), and
the model input is Theta_m(s) = sum_i theta_i.
"""

from __future__ import annotations

import numpy as np

LOOKUP = np.full(256, -1, dtype=np.int16)
LOOKUP[np.frombuffer(b"ACGT", dtype=np.uint8)] = np.arange(4)


def encode(sequences: list[str]) -> np.ndarray:
    if not sequences or any(len(s) != len(sequences[0]) for s in sequences):
        raise ValueError("sequences must be nonempty and equally long")
    result = LOOKUP[np.frombuffer("".join(sequences).encode("ascii"), dtype=np.uint8)]
    result = result.reshape(len(sequences), -1)
    if (result < 0).any():
        raise ValueError("sequence alphabet must be ACGT")
    return result


def site_sums(sequences: list[str], matrices: list[np.ndarray], kappas, *, batch: int = 512) -> np.ndarray:
    """Return Theta with shape (len(kappas), len(sequences), len(matrices))."""
    kappas = np.asarray(kappas, dtype=np.float64)
    if kappas.ndim != 1 or (kappas <= 0).any():
        raise ValueError("kappas must be a positive vector")
    out = np.zeros((len(kappas), len(sequences), len(matrices)))
    widths: dict[int, list[int]] = {}
    for j, m in enumerate(matrices):
        if m.ndim != 2 or m.shape[0] != 4 or not np.isfinite(m).all():
            raise ValueError("invalid motif matrix")
        widths.setdefault(m.shape[1], []).append(j)
    encoded = encode(sequences)
    length = encoded.shape[1]
    for width, columns in sorted(widths.items()):
        if width > length:
            continue
        group = np.stack([matrices[j] for j in columns])
        kernels = np.concatenate((group, group[:, ::-1, ::-1])).transpose(2, 1, 0).reshape(4 * width, -1)
        for start in range(0, len(sequences), batch):
            letters = np.eye(4, dtype=np.float64)[encoded[start:start + batch]]
            windows = np.lib.stride_tricks.sliding_window_view(letters, width, axis=1).transpose(0, 1, 3, 2)
            scores = windows.reshape(-1, 4 * width) @ kernels
            forward, reverse = np.split(scores, 2, axis=1)
            with np.errstate(over="raise"):
                odds = (np.exp(forward) + np.exp(reverse)) / 2
            odds = odds.reshape(windows.shape[0], windows.shape[1], len(columns))
            for k, kappa in enumerate(kappas):
                x = kappa * odds
                out[k, start:start + batch][:, columns] = (x / (1.0 + x)).sum(axis=1)
    if not np.isfinite(out).all():
        raise ValueError("nonfinite site sums")
    return out


def shuffle_columns(matrices: list[np.ndarray], seed: int) -> list[np.ndarray]:
    """PWM-shuffle null: permute the position columns of every motif, keeping its column composition."""
    rng = np.random.default_rng(seed)
    result = []
    for m in matrices:
        result.append(m[:, rng.permutation(m.shape[1])].copy())
    return result

"""Qualified frozen inference functions; only private module bindings are flattened."""

from collections import Counter, defaultdict
import itertools
import numpy as np
COMPLEMENT = str.maketrans("ACGT", "TGCA")

def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]

def vocabulary(k):
    return sorted({min(word, reverse_complement(word))
                   for letters in itertools.product('ACGT', repeat=k)
                   for word in [''.join(letters)]})

def whole_counts(sequences, k):
    columns = {word: index for index, word in enumerate(vocabulary(k))}
    matrix = np.zeros((len(sequences), len(columns)))
    for index, sequence in enumerate(sequences):
        if len(sequence) < k or set(sequence)-set('ACGT'):
            raise ValueError('invalid sequence')
        counts = Counter(min(sequence[i:i+k], reverse_complement(sequence[i:i+k]))
                         for i in range(len(sequence)-k+1))
        for word, count in counts.items():
            matrix[index, columns[word]] = count
    return matrix

def encode(sequences):
    if not sequences or not sequences[0] or any(len(sequence) != len(sequences[0]) for sequence in sequences):
        raise ValueError('sequences must be nonempty and equally long')
    lookup = np.full(256, -1, dtype=np.int16)
    lookup[np.frombuffer(b'ACGT', dtype=np.uint8)] = np.arange(4)
    result = lookup[np.frombuffer(''.join(sequences).encode('ascii'), dtype=np.uint8)].reshape(len(sequences), -1)
    if np.any(result < 0):
        raise ValueError('sequence alphabet must be ACGT')
    return result

def changed_words(reference, alternate, width):
    if type(width) is not int or width <= 0 or not reference or not alternate or set(reference+alternate)-set('ACGT'):
        raise ValueError('invalid native pair or word width')
    def counts(sequence):
        return Counter(min(word, reverse_complement(word))
                       for start in range(len(sequence)-width+1)
                       for word in [sequence[start:start+width]])
    delta = counts(alternate)
    delta.subtract(counts(reference))
    return [(word, count) for word, count in sorted(delta.items()) if count]

def word_odds(words, matrices):
    encoded = encode(words)
    group = np.stack(matrices)
    if group.ndim != 3 or group.shape[1:] != (4, encoded.shape[1]) or not np.isfinite(group).all():
        raise ValueError('invalid same-width motif group')
    forward = np.zeros((len(words), len(group)))
    reverse = np.zeros_like(forward)
    for offset in range(encoded.shape[1]):
        letters = encoded[:, offset]
        forward += group[:, letters, offset].T
        reverse += group[:, 3-letters, encoded.shape[1]-1-offset].T
    result = (np.exp(forward)+np.exp(reverse))/2
    if not np.isfinite(result).all():
        raise ValueError('nonfinite motif likelihood ratio')
    return result

def deltas(rows, matrices, scale):
    if not rows or not matrices or not np.isfinite(scale) or scale <= 0:
        raise ValueError('invalid motif input or compression scale')
    grouped = defaultdict(list)
    for column, matrix in enumerate(matrices):
        if matrix.ndim != 2 or matrix.shape[0] != 4 or matrix.shape[1] <= 0 or not np.isfinite(matrix).all():
            raise ValueError('invalid motif matrix')
        grouped[matrix.shape[1]].append(column)
    result = np.zeros((len(rows), len(matrices)))
    for width, columns in sorted(grouped.items()):
        changes = [changed_words(row['reference'], row['sequence'], width) for row in rows]
        words = sorted({word for change in changes for word, count in change})
        if not words:
            continue
        index = {word: i for i, word in enumerate(words)}
        odds = word_odds(words, [matrices[c] for c in columns])
        baseline = np.log1p(scale)
        uncentered = np.log1p(scale*odds)
        centered = np.log1p((scale/(1+scale))*(odds-1))
        # use the smaller representation and cancel weak-site offsets as integers.
        weak = np.abs(uncentered) < np.abs(centered)
        values = np.where(weak, uncentered, centered)
        for row, change in enumerate(changes):
            if change:
                ids = [index[word] for word, count in change]
                counts = np.array([count for word, count in change], dtype=np.int64)
                result[row, columns] = counts @ values[ids]-(counts @ weak[ids])*baseline
    if not np.isfinite(result).all():
        raise ValueError('nonfinite logarithmic motif delta')
    return result

def endpoints(sequences, matrices, scale):
    if not sequences or not matrices or not np.isfinite(scale) or scale <= 0:
        raise ValueError('invalid motif input or compression scale')
    if any(not s or set(s)-set('ACGT') for s in sequences):
        raise ValueError('invalid native sequence')
    lengths = defaultdict(list); widths = defaultdict(list)
    sequences = [min(s, reverse_complement(s)) for s in sequences]
    for i, sequence in enumerate(sequences):
        lengths[len(sequence)].append(i)
    for j, matrix in enumerate(matrices):
        if matrix.ndim != 2 or matrix.shape[0] != 4 or matrix.shape[1] <= 0 or not np.isfinite(matrix).all():
            raise ValueError('invalid motif matrix')
        widths[matrix.shape[1]].append(j)
    result = np.zeros((len(sequences), len(matrices)))
    with np.errstate(over='raise', invalid='raise'):
        for length, indices in sorted(lengths.items()):
            encoded = encode([sequences[i] for i in indices])
            letters = np.eye(4)[encoded]
            for width, columns in sorted(widths.items()):
                if width > length:
                    continue
                group = np.stack([matrices[j] for j in columns])
                kernels = np.concatenate((group, group[:, ::-1, ::-1])).transpose(2, 1, 0).reshape(4*width, -1)
                windows = np.lib.stride_tricks.sliding_window_view(letters, width, axis=1).transpose(0, 1, 3, 2)
                scores = windows.reshape(-1, 4*width) @ kernels
                forward, reverse = np.split(scores, 2, axis=1)
                odds = (np.exp(forward)+np.exp(reverse))/2
                uncentered = np.log1p(scale*odds)
                centered = np.log1p((scale/(1+scale))*((np.expm1(forward)+np.expm1(reverse))/2))
                weak = np.abs(uncentered) < np.abs(centered)
                values = np.where(weak, uncentered, centered).reshape(len(indices), length-width+1, -1)
                totals = values.sum(axis=1)-weak.reshape(values.shape).sum(axis=1)*np.log1p(scale)
                result[np.ix_(indices, columns)] = totals
    if not np.isfinite(result).all():
        raise ValueError('nonfinite motif potential')
    return result

def quadratic(rows, matrices, scale, delta):
    if delta.shape != (len(rows), len(matrices)) or not np.isfinite(delta).all():
        raise ValueError('invalid linear motif deltas')
    totals = endpoints([r[key] for key in ('reference', 'sequence') for r in rows], matrices, scale)
    reference, alternate = np.split(totals, 2)
    result = delta*(alternate+reference)
    if not np.isfinite(result).all():
        raise ValueError('nonfinite quadratic motif delta')
    return result

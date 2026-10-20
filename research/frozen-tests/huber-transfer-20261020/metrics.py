"""Paired whole-family-copy summaries for a fixed transfer population."""
import numpy as np


def evaluate(ids, families, y, predictions, counts, budgets):
    ids, families, y = np.asarray(ids), np.asarray(families), np.asarray(y)
    groups, inverse, sizes = np.unique(families, return_inverse=True, return_counts=True)
    n, f = len(ids), len(groups)
    if (len(set(ids.tolist())) != n or families.shape != (n,) or y.shape != (n,)
            or counts.ndim != 2 or counts.shape[1] != f or f < 2
            or not np.isfinite(y).all() or not np.isfinite(counts).all()
            or np.any(counts < 0) or np.any(counts != np.floor(counts))
            or not np.all(counts.sum(axis=1) == f)):
        raise ValueError('invalid whole-family-copy population')
    mean = lambda value: np.bincount(inverse, weights=value/sizes[inverse], minlength=f)
    result = {}
    for name, p in predictions.items():
        p = np.asarray(p)
        if p.shape != (n,) or not np.isfinite(p).all():
            raise ValueError('invalid prediction vector')
        yc, pc = y-y.mean(), p-p.mean()
        moments = counts@np.stack([mean(yc), mean(yc*yc), mean(pc), mean(pc*pc), mean(yc*pc)], axis=1)/f
        ym, yy, pm, pp, yp = moments.T
        vy, vp = yy-ym*ym, pp-pm*pm
        if (not np.isfinite(moments).all() or not np.isfinite(vy).all()
                or not np.isfinite(vp).all() or np.any(vy <= 0) or np.any(vp <= 0)):
            raise ValueError('undefined correlation in a family draw')
        order = np.lexsort((ids, -np.abs(p)))
        first = np.full(f, -1, dtype=int)
        for i in order[::-1]:
            first[inverse[i]] = i
        chosen = first[np.lexsort((ids[first], -np.abs(p[first])))]
        copies = counts[:, inverse[chosen]]
        before = np.cumsum(copies, axis=1)-copies
        utility = {}
        for budget in budgets:
            if not isinstance(budget, int) or budget <= 0:
                raise ValueError('invalid panel budget')
            utility[budget] = (None if f < budget else
                np.minimum(np.maximum(budget-before, 0), copies)@(np.sign(p[chosen])*y[chosen])/budget)
        result[name] = dict(r=(yp-ym*pm)/(np.sqrt(vy)*np.sqrt(vp)),
                            mse=counts@mean((y-p)**2)/f, utility=utility)
        if not all(np.isfinite(value).all() for value in
                   [result[name]['r'], result[name]['mse'], *[v for v in utility.values() if v is not None]]):
            raise ValueError('nonfinite transfer metric')
    return result


def summary(values):
    values = np.asarray(values)
    if values.ndim != 1 or len(values) < 2 or not np.isfinite(values).all():
        raise ValueError('invalid point/bootstrap vector')
    return dict(point=float(values[0]), q05_50_95=np.quantile(values[1:], [.05,.5,.95]).tolist())

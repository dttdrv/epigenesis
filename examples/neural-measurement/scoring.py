"""Unchanged scoring functions extracted from the archived independent checker."""
import math
import sys
import numpy as np

N, JOBS, DELTA, TOLERANCE = 2000000, 12, .01, 1e-8

def same(actual, expected, label, tolerance=2e-11):
    left, right = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
    if left.shape != right.shape or not np.isfinite(left).all() or not np.isfinite(right).all():
        raise AssertionError(label + ': shape/finiteness differs')
    error = float(np.max(np.abs(left-right))) if left.size else 0.
    if error > tolerance:
        raise AssertionError(f'{label}: maximum absolute difference {error}')
    return error


def pooled_mass(training, total):
    alpha = 1+sum(pair[1] for pair in training)
    beta = 1+sum(pair[0] for pair in training)
    log_mass, running = [0.], 0.
    for k in range(total):
        running += math.fsum((math.log(total-k), math.log(alpha+k),
                              -math.log(k+1), -math.log(beta+total-k-1)))
        log_mass.append(running)
    peak = max(log_mass)
    scaled = [math.exp(value-peak) for value in log_mass]
    normalizer = math.fsum(scaled)
    return [value/normalizer for value in scaled]


def prefix(values):
    running, result = 0., []
    for value in values:
        running += float(value)
        result.append(running)
    return result


def pair_loss(masses, observed):
    total = len(masses)-1
    first = math.fsum(p*abs(k-observed) for k,p in enumerate(masses))
    preceding_mass, preceding_count, terms = 0., 0., []
    for k,p in enumerate(masses):
        terms.append(p*(k*preceding_mass-preceding_count))
        preceding_mass += p
        preceding_count += k*p
    return (first-math.fsum(terms))/total


def enclosure_loss(lower, upper, observed):
    low, high = [], []
    for k,(left,right) in enumerate(zip(lower, upper, strict=True)):
        low.append(left*left if k < observed else (1-right)**2)
        high.append(right*right if k < observed else (1-left)**2)
    return [math.fsum(low)/len(low), math.fsum(high)/len(high)]


def summarize(sums, squares, size, envelope, observed, *, jobs=JOBS, delta=DELTA, tolerance=TOLERANCE):
    sums, squares = list(map(float, sums)), list(map(float, squares))
    if (size < 2 or len(sums) < 2 or len(sums) != len(squares)
            or any(not math.isfinite(x) or x < 0 for x in sums+squares)
            or not 0 < envelope <= 1 or not 0 <= observed < len(sums)):
        raise ValueError('invalid moment contract')
    total = len(sums)-1
    cap = envelope*math.exp(tolerance)
    rounding = 4*size*sys.float_info.epsilon
    first = [math.fsum(sums), *prefix(sums[:-1])]
    second = [math.fsum(squares), *prefix(squares[:-1])]
    if first[0] > size*cap*(1+rounding):
        raise ValueError('sum exceeds raw population bound')
    for a,b in zip(first+sums, second+squares, strict=True):
        if b > cap*a*(1+rounding) or b+rounding*size*cap**2 < a*a/size:
            raise ValueError('impossible bounded moments')
    level = math.log(4*jobs*(total+1)/delta)
    intervals, variances = [], []
    for a,b in zip(first, second, strict=True):
        average = a/size
        variance = max(0., math.fsum((b, -a*a/size))/(size-1))+rounding*cap**2
        radius = math.sqrt(2*variance*level/size)+7*cap*level/(3*(size-1))+rounding*cap
        intervals.append([max(0., average-radius), min(cap, average+radius)])
        variances.append(variance)
    dlow, dhigh = intervals[0]
    distortion = math.expm1(2*tolerance)
    lower, upper = [], []
    for a,b in intervals[1:]:
        lower.append(max(0., (a/dhigh if dhigh else 0.)-distortion))
        upper.append(min(1., (min(1., b/dlow) if dlow else 1.)+distortion))
    mass = [value/first[0] for value in sums] if first[0] else None
    cdf = [value/first[0] for value in first[1:]] if first[0] else None
    return dict(raw_particles=size, mean_weight=first[0]/size,
                effective_particles=first[0]**2/second[0] if second[0] else 0.,
                cdf=cdf, cdf_lower=lower, cdf_upper=upper,
                loss=pair_loss(mass, observed) if mass is not None else None,
                loss_interval=enclosure_loss(lower, upper, observed),
                denominator_interval=intervals[0], variance=variances)

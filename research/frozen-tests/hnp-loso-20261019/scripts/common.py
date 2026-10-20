"""shared paths and loaders for the three human-neural-progenitor reporter studies."""
import csv, gzip, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
BUILD = ROOT.parent
LEE, MCA = BUILD / 'lee2025-gse244011-20261019', BUILD / 'mcafee2023-gse211045-20261019'
MAT = BUILD / 'matoba-measurement/train-measurements'
DER = ROOT / 'derived'
STUDIES = ('lee', 'mcafee', 'matoba')
STD = {*map(str, range(1, 23)), 'X'}
WINDOW = 2000  # bases; same rule as the Lee transfer test


def table(path, delim=','):
    with gzip.open(path, 'rt') as fh:
        t = list(csv.reader(fh, delimiter=delim))
    return t[0], {r[0]: np.array([float(x) for x in r[1:]]) for r in t[1:]}


def lee_pos(v):  # chr1:104789431_A:G
    c, rest = v.split(':', 1); p, al = rest.split('_'); r, a = al.split(':')
    return c.removeprefix('chr'), int(p), r, a


def ratio(rna, dna):
    return np.log2((rna.sum(1) + 0.5) / (dna.sum(1) + 0.5))

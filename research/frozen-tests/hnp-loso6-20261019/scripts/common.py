"""shared paths for the six neural-cell variant studies."""
import bisect, csv, gzip, json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
BUILD = ROOT.parent
V1 = BUILD / 'hnp-loso-20261019/derived'
DER = ROOT / 'derived'
OLD = ('lee', 'mcafee', 'matoba')            # scored in the three-study test; outcomes already seen
NEW = ('gse225815', 'gse221384', 'weiss')    # first scored here
STUDIES = OLD + NEW
WINDOW = 2000

def ratio(rna, dna):
    return np.log2((np.asarray(rna, float) + 0.5) / (np.asarray(dna, float) + 0.5))

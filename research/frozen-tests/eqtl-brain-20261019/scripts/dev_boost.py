"""development only (Walker + Aygun are spent as test sets): compare model families by leave-chromosome-out accuracy/AUC on the pooled development pairs."""
import json, sys, time
from pathlib import Path
import numpy as np
exec(open(Path(__file__).resolve().parent / 'run_eqtl.py').read().split("key = lambda r:")[0])   # reuse data + feature code (no fitting)
from sklearn.ensemble import HistGradientBoostingClassifier as HGB
from sklearn.metrics import roc_auc_score
tab = {(r['dataset'], r['chrom'], r['pos'], r['alt'], r['gene']): r['beta'] for r in json.load(open(HERE / 'derived/eqtl_table.json')) if r['role'] == 'test'}
# run_eqtl strips beta from `test`; the lookup above supplies labels for development scoring only
yte = np.array([1.0 if tab[(r['dataset'], r['chrom'], r['pos'], r['alt'], r['gene'])] > 0 else 0.0 for r in test])
def dense(rs):
    vi = np.array([r['vi'] for r in rs]); d = np.array([r['tss_dist'] for r in rs], float)
    return np.hstack([SEI[vi], np.c_[np.log1p(np.abs(d)), np.sign(d), GC[vi], S2[vi]]])
cnt = {}
for r in train: k = (r['chrom'], r['pos'], r['alt'], r['gene']); cnt[k] = cnt.get(k, 0) + 1
w = np.array([1.0 / cnt[(r['chrom'], r['pos'], r['alt'], r['gene'])] for r in train])
chrs = np.array([r['chrom'] for r in train]); te_chr = np.array([r['chrom'] for r in test])
Xtr_all, Xte_all = dense(train), dense(test)
cfgs = {'gbm_small': dict(max_iter=200, learning_rate=0.03, max_leaf_nodes=6, min_samples_leaf=40, l2_regularization=5.0, max_features=0.3),
        'gbm_mid': dict(max_iter=400, learning_rate=0.03, max_leaf_nodes=8, min_samples_leaf=30, l2_regularization=5.0, max_features=0.3)}
P = {k: np.zeros(len(test)) for k in cfgs}
for c in sorted(set(te_chr)):
    tr = chrs != c; ti = te_chr == c
    for k, cfg in cfgs.items():
        mdl = HGB(**cfg, random_state=1).fit(Xtr_all[tr], ytr_all[tr], sample_weight=w[tr])
        P[k][ti] = mdl.predict_proba(Xte_all[ti])[:, 1]
base = np.load(HERE / 'derived/PREDICTIONS.npz')
ds = np.array([r['dataset'] for r in test])
for name, sc in {**P, 'M4_logistic': base['M4'], 'M5_logistic': base['M5'], 'B3_reporter': base['B3']}.items():
    print(name, 'pooled acc %.3f auc %.3f' % (((sc > 0.5) == (yte == 1)).mean(), roc_auc_score(yte, sc)), ' walker auc %.3f' % roc_auc_score(yte[ds == 'walker_fetal_neocortex'], sc[ds == 'walker_fetal_neocortex']))

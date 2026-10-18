"""Replay saved aggregate predictions only; no source imports or model fits."""
import argparse
from collections import defaultdict
import hashlib
import gzip
import json
import math
from pathlib import Path
import sys

import numpy as np

METHODS = ('level','flat','latest','zero','population')
CONTROLS = METHODS[1:]
HOURS = (0,3,6,12,24,48,72)
PENALTIES = (.001,.01,.1,1.,10.)
DEFAULT_ROOT = Path(__file__).resolve().parent


def mean(values):
    values = list(values)
    assert values, 'cannot average an unobserved cell'
    return math.fsum(values)/len(values)


def fold(group):
    return int.from_bytes(hashlib.sha256(('epigenesis-temporal-v1\n'+group).encode()).digest(),'big')%5


def archived_bytes(path):
    if path.exists():
        return path.read_bytes()
    with gzip.open(path.with_suffix(path.suffix+'.gz'),'rb') as stream:
        return stream.read()


def digest(path):
    return hashlib.sha256(archived_bytes(path)).hexdigest()


def features(prefix, mode):
    if mode == 'level':
        return [mean(prefix[:3]),float(prefix[3])]
    if mode == 'flat':
        return [mean(prefix)]
    assert mode == 'latest'
    return [float(prefix[3])]


def main(root, output):
    out = root/'result'
    assert (out/'RESULT.json.gz').exists(), 'complete result is not present; do not score a running attempt'
    assert not list(out.glob('FAILURE.json*')), 'attempt recorded a failure'
    if not __debug__:
        raise RuntimeError('replay requires assertions; omit Python -O')
    pins = json.loads((root/'MANIFEST.json').read_text())
    paths = {path.relative_to(root).as_posix() for path in root.rglob('*') if path.is_file() and '__pycache__' not in path.parts and path.name!='MANIFEST.json'}
    assert paths == set(pins['files']), 'archive is incomplete or contains extra files'
    for name, expected in pins['files'].items():
        path = (root/name).resolve()
        assert path.is_relative_to(root.resolve()), 'archive path leaves its root'
        assert digest(path) == expected, 'archive identity differs: '+name
    result = json.loads(archived_bytes(out/'RESULT.json'))
    summary = json.loads((root/'SUMMARY.json').read_text())
    assert summary == {key:result[key] for key in summary}, 'public summary differs'
    attempt = json.loads(archived_bytes(out/'ATTEMPT.json'))
    assert attempt['freeze_sha256'] == pins['original_freeze_sha256']
    rows = json.loads(archived_bytes(out/'ROWS.json'))
    assert rows == json.loads(archived_bytes(root/'COHORT.json')), 'fixed cohort differs'
    assert len({r['id'] for r in rows}) == len(rows), 'duplicate cohort identity'
    folds = np.array([fold(r['group']) for r in rows])
    assert folds.tolist() == [r['fold'] for r in rows]
    assert set(folds) == set(range(5))
    admission = json.loads(archived_bytes(out/'ADMISSION.json'))
    assert admission['times'] == list(HOURS) and admission['cultures'] == [1,2,3]
    assert admission['rows'] == len(rows)
    assert admission['families'] == len({r['group'] for r in rows})
    with np.load(out/'COUNTS.npz',allow_pickle=False) as saved:
        assert set(saved.files) == {'effects','variant','parent'}
        effects,variant,parent = (saved[k] for k in ('effects','variant','parent'))
    assert effects.shape == (len(rows),3,7)
    assert variant.shape == parent.shape == (len(rows),3,7,3)
    for totals in (variant,parent):
        assert np.isfinite(totals).all() and np.all(totals>=0)
        assert np.all(totals<2**53) and np.all(totals==np.floor(totals))
        assert np.array_equal(totals[:,:, :,2]>0,(totals[:,:,:,0]>0)&(totals[:,:,:,1]>0))
    valid = (variant[:,:,:,:2]>0).all(axis=3)&(parent[:,:,:,:2]>0).all(axis=3)
    assert np.array_equal(np.isfinite(effects),valid)
    assert admission['rows_with_support_by_culture_time'] == valid.sum(axis=0).tolist()
    differences = []
    def close(actual, expected, label):
        a,b = np.asarray(actual,dtype=float),np.asarray(expected,dtype=float)
        assert a.shape == b.shape, (label,a.shape,b.shape)
        assert np.isfinite(a).all() and np.isfinite(b).all(), label
        if a.size:
            differences.append(float(np.max(np.abs(a-b))))
        np.testing.assert_allclose(a,b,rtol=1e-11,atol=2e-12,err_msg=label)
    reconstructed = []
    for i,c,t in zip(*np.nonzero(valid)):
        d,r,_ = variant[i,c,t]; pd,pr,_ = parent[i,c,t]
        reconstructed.append(math.fsum((math.log2(r),-math.log2(d),-math.log2(pr),math.log2(pd))))
    close(effects[valid],reconstructed,'saved log-count observable')
    absent = set(admission['absent_ids'])
    def missing(indices,cultures,times):
        ledger = []
        for i in indices:
            for c in cultures:
                for t in times:
                    reasons = {}
                    if variant[i,c,t,2] == 0:
                        reasons['variant'] = 'absent_row' if rows[i]['id'] in absent else 'no_matched_positive_slots'
                    if parent[i,c,t,2] == 0:
                        ids = admission['parent_ids_by_locus'][rows[i]['locus']]
                        reasons['parent'] = 'all_parent_rows_absent' if all(k in absent for k in ids) else 'no_matched_positive_slots'
                    if reasons:
                        ledger.append({'row':int(i),'culture':c+1,'hour':HOURS[t],**reasons})
        return ledger
    seals = json.loads(archived_bytes(out/'PREDICTION_SEALS.json'))
    expected_seals = {f'culture-{c+1}-fold-{k}{ext}' for c in range(3) for k in range(5) for ext in ('.npz','.json')}
    assert set(seals) == expected_seals
    for name, expected in seals.items():
        assert digest(out/name) == expected, 'prediction seal changed: '+name
    cells = defaultdict(dict)
    strata,availability = [],[]
    for c in range(3):
        other = [k for k in range(3) if k != c]
        for k in range(5):
            stem = f'culture-{c+1}-fold-{k}'
            receipt = json.loads(archived_bytes(out/(stem+'.json')))
            assert receipt['omitted_culture'] == c+1 and receipt['held_fold'] == k
            ix = np.flatnonzero(folds==k)
            eligible = np.array([all(valid[i,a,:4]) and all(valid[i,b,:4]) for i in ix for a,b in [other]],dtype=bool)
            with np.load(out/(stem+'.npz'),allow_pickle=False) as saved:
                assert set(saved.files) == {'rows','eligible',*METHODS}
                assert np.array_equal(saved['rows'],ix)
                assert np.array_equal(saved['eligible'],eligible)
                predictions = {m:saved[m] for m in METHODS}
            for p in predictions.values():
                assert p.shape == (int(eligible.sum()),3) and np.isfinite(p).all()
            assert receipt['input_missing'] == missing(ix,other,range(4))
            training = []
            missing_training = []
            for a,b in (other,other[::-1]):
                rejected = []
                for i in np.flatnonzero(folds!=k):
                    if all(valid[i,a,:4]) and all(valid[i,b,3:]):
                        training.append((int(i),a,b,rows[i]['group']))
                    else:
                        rejected.append(int(i))
                missing_training.append({'source':a+1,'target':b+1,'rows':rejected})
            assert receipt['training_rows'] == [v[0] for v in training]
            assert receipt['training_source_culture'] == [v[1]+1 for v in training]
            assert receipt['training_target_culture'] == [v[2]+1 for v in training]
            assert receipt['training_missing'] == missing_training
            assert {fold(v[3]) for v in training} == set(range(5))-{k}
            tcells = defaultdict(lambda:defaultdict(list))
            for n,(_,_,target,g) in enumerate(training):
                tcells[g][target].append(n)
            weights = [0.]*len(training)
            for family in tcells.values():
                for obs in family.values():
                    for n in obs:
                        weights[n] = 1/len(tcells)/len(family)/len(obs)
            targets = [[float(effects[i,b,t]-effects[i,b,3]) for t in (4,5,6)] for i,a,b,g in training]
            target_mean = [math.fsum(w*y[h] for w,y in zip(weights,targets)) for h in range(3)]
            close(predictions['zero'],np.zeros_like(predictions['zero']),stem+' zero')
            close(predictions['population'],np.tile(target_mean,(int(eligible.sum()),1)),stem+' population')
            for mode in ('level','flat','latest'):
                model = receipt['models'][mode]
                assert model['penalty'] in PENALTIES
                inner = model['inner_mse']
                assert len(inner) == 5 and all(math.isfinite(x) and x>=0 for x in inner)
                assert model['penalty'] == PENALTIES[min(range(5),key=lambda n:inner[n])]
                x = [features(effects[i,a,:4],mode) for i,a,b,g in training]
                center = [math.fsum(w*v[j] for w,v in zip(weights,x)) for j in range(len(x[0]))]
                for j in range(len(center)):
                    if min(v[j] for v in x) == max(v[j] for v in x):
                        center[j] = x[0][j]
                scale = [math.sqrt(math.fsum(w*(v[j]-center[j])**2 for w,v in zip(weights,x))) for j in range(len(center))]
                scale = [v if v else 1. for v in scale]
                close(model['center'],center,stem+' '+mode+' center')
                close(model['scale'],scale,stem+' '+mode+' scale')
                close(model['mean'],target_mean,stem+' '+mode+' intercept')
                predicted = []
                for i in ix[eligible]:
                    directions = []
                    for a in other:
                        xx = features(effects[i,a,:4],mode)
                        directions.append([math.fsum([model['mean'][h]]+[(xx[j]-model['center'][j])/model['scale'][j]*model['beta'][j][h] for j in range(len(xx))]) for h in range(3)])
                    predicted.append([mean(v[h] for v in directions) for h in range(3)])
                close(predictions[mode],np.array(predicted).reshape((-1,3)),stem+' '+mode+' sealed coefficients')
            candidates = ix[eligible]
            observed = np.array([all(valid[i,c,3:]) for i in candidates],dtype=bool)
            observed_families = sorted({rows[i]['group'] for i in candidates[observed]})
            assert observed_families, 'unobserved outer stratum'
            strata.append({'culture':c+1,'fold':k,'population':len(ix),'input_eligible':int(eligible.sum()),'complete_targets':int(observed.sum()),'unscored_predictions':int((~observed).sum()),'families':len(observed_families)})
            availability.append({'culture':c+1,'fold':k,'rows':ix.tolist(),'source_eligible':eligible.tolist(),'target_support_by_time':valid[ix,c,3:].tolist()})
            row_losses = defaultdict(lambda:defaultdict(list))
            for n,i in enumerate(candidates):
                if not observed[n]:
                    continue
                truth = [float(effects[i,c,t]-effects[i,c,3]) for t in (4,5,6)]
                for m in METHODS:
                    row_losses[rows[i]['group']][m].append(mean((float(predictions[m][n,h])-truth[h])**2 for h in range(3)))
            for g in observed_families:
                cells[g][c] = {m:mean(row_losses[g][m]) for m in METHODS}
    assert result['status'] == 'complete'
    assert result['strata'] == strata and result['availability'] == availability
    assert result['target_missing'] == missing(range(len(rows)),range(3),range(3,7))
    families = sorted(cells)
    assert result['families'] == len(families)
    losses = [[mean(cell[m] for cell in cells[g].values()) for m in METHODS] for g in families]
    for j,m in enumerate(METHODS):
        close(result['mse'][m],mean(v[j] for v in losses),'overall '+m)
    assert set(result['family_cells']) == set(cells)
    for g in families:
        assert set(result['family_cells'][g]) == {str(c+1) for c in cells[g]}
        for c,cell in cells[g].items():
            for m in METHODS:
                close(result['family_cells'][g][str(c+1)][m],cell[m],'family/culture '+m)
    culture_gains = []
    assert len(result['cultures']) == 3
    for c,record in enumerate(result['cultures']):
        selected = [cells[g][c] for g in families if c in cells[g]]
        assert record['culture'] == c+1 and record['families'] == len(selected)
        cm = {m:mean(v[m] for v in selected) for m in METHODS}
        cg = {m:cm[m]-cm['level'] for m in CONTROLS}
        for m in METHODS:
            close(record['mse'][m],cm[m],'culture '+str(c+1)+' '+m)
        for m in CONTROLS:
            close(record['gain'][m],cg[m],'culture gain '+m)
        culture_gains.extend(cg.values())
    gains = np.array([[v[j]-v[0] for j in range(1,5)] for v in losses])
    rng = np.random.Generator(np.random.PCG64(20261018))
    boot = []
    for _ in range(20):
        samples = rng.integers(len(families),size=(1000,len(families)))
        boot.extend([[mean(gains[i,j] for i in sample) for j in range(4)] for sample in samples])
    bounds = np.quantile(np.asarray(boot),[.01/8,1-.01/8],axis=0,method='linear')
    assert len(result['contrasts']) == 4
    for j,(m,record) in enumerate(zip(CONTROLS,result['contrasts'])):
        assert record['control'] == m
        close(record['gain'],mean(gains[:,j]),'gain '+m)
        close(record['interval'],bounds[:,j],'bootstrap '+m)
    passed = bool(np.all(bounds[0]>0)) and all(g>0 for g in culture_gains)
    expected = 'conditional_on_supported_held_culture_measurements' if passed else 'not_passed'
    assert result['decision'] == expected
    receipt = {'scope':'saved aggregates and sealed predictions only; no count-matrix conversion or model fitting','root':str(root),'runtime':{'python':sys.version.split()[0],'numpy':np.__version__},'files':{str(p):digest(p) for p in (root/'MANIFEST.json',out/'COUNTS.npz',out/'ROWS.json',out/'PREDICTION_SEALS.json',out/'RESULT.json',Path(__file__))},'rows':len(rows),'families':len(families),'strata':len(strata),'verified_prediction_files':len(seals),'maximum_absolute_numeric_difference':max(differences),'decision':expected,'mse':result['mse'],'contrasts':result['contrasts'],'cultures':result['cultures']}
    if output is not None:
        with output.open('x') as stream:
            stream.write(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ('rows','families','strata','verified_prediction_files','maximum_absolute_numeric_difference','decision','mse')},indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=DEFAULT_ROOT)
    parser.add_argument('--output',type=Path,default=None)
    args = parser.parse_args()
    main(args.root,args.output)

"""Whole-family and whole-culture exclusion for early-to-late changes."""
from collections import Counter,defaultdict
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

import counts

HERE=Path(__file__).resolve().parent
PENALTIES=(.001,.01,.1,1.,10.)
METHODS=('level','flat','latest','zero','population')
CONTROLS=METHODS[1:]


def fold(group):
    return int(hashlib.sha256(('epigenesis-temporal-v1\n'+group).encode()).hexdigest(),16)%5


def weights(groups,cultures):
    sizes=Counter(zip(groups,cultures))
    cells=Counter(g for g,c in sizes)
    return np.asarray([1/(len(cells)*cells[g]*sizes[g,c]) for g,c in zip(groups,cultures)])


def features(x,mode):
    if mode=='level':return np.column_stack((x[:,:3].mean(axis=1),x[:,-1]))
    if mode=='flat':return x.mean(axis=1,keepdims=True)
    if mode=='latest':return x[:,-1:]
    raise ValueError('unknown model')


def fit(x,y,w,penalty):
    xm,ym=w@x,w@y
    constant=np.ptp(x,axis=0)==0
    xm[constant]=x[0,constant]
    sd=np.sqrt(w@((x-xm)**2));sd[sd==0]=1.
    z=(x-xm)/sd
    beta=np.linalg.solve(z.T@(w[:,None]*z)+penalty*np.eye(x.shape[1]),z.T@(w[:,None]*(y-ym)))
    return dict(center=xm,scale=sd,mean=ym,beta=beta)


def predict(model,x):return (x-model['center'])/model['scale']@model['beta']+model['mean']


def outer(effects,groups,culture,heldfold):
    families=np.asarray(groups)
    folds=np.asarray([fold(g) for g in groups])
    other=[c for c in range(3) if c!=culture]
    train_rows=np.flatnonzero(folds!=heldfold)
    index=[];sources=[];targets=[];missing=[]
    xs=[];ys=[]
    for a,b in (other,other[::-1]):
        x=effects[train_rows,a,:4]
        target=effects[train_rows,b,3:]
        valid=np.isfinite(x).all(axis=1)&np.isfinite(target).all(axis=1)
        used=train_rows[valid]
        xs.append(x[valid]);ys.append(target[valid,1:]-target[valid,:1])
        index.extend(used.tolist());sources.extend([a]*len(used));targets.extend([b]*len(used))
        missing.append({'source':a+1,'target':b+1,'rows':train_rows[~valid].tolist()})
    index=np.asarray(index,dtype=int);targets=np.asarray(targets,dtype=int)
    x,y=np.concatenate(xs),np.concatenate(ys)
    if not len(x):raise ValueError('no training observations')
    g=families[index];f=folds[index];w=weights(g,targets)
    if set(f)!=set(range(5))-{heldfold}:raise ValueError('missing required inner family fold')
    evaluation=np.flatnonzero(folds==heldfold)
    source_prefixes=[effects[evaluation,c,:4] for c in other]
    eligible=np.logical_and.reduce([np.isfinite(a).all(axis=1) for a in source_prefixes])
    contexts=[a[eligible] for a in source_prefixes]
    predictions={'zero':np.zeros((sum(eligible),3)),'population':np.tile(w@y,(sum(eligible),1))}
    receipt={'omitted_culture':culture+1,'held_fold':heldfold,'models':{},'training_rows':index.tolist(),
             'training_source_culture':(np.asarray(sources)+1).tolist(),'training_target_culture':(targets+1).tolist(),
             'training_missing':missing}
    for mode in ('level','flat','latest'):
        scores=np.zeros(len(PENALTIES));nfamily=0
        for held in sorted(set(f)):
            train=f!=held;check=f==held
            if not train.any() or not check.any():raise ValueError('missing inner family fold')
            tw=weights(g[train],targets[train]);cw=weights(g[check],targets[check]);n=len(set(g[check]))
            for j,penalty in enumerate(PENALTIES):
                model=fit(features(x[train],mode),y[train],tw,penalty)
                error=np.mean((predict(model,features(x[check],mode))-y[check])**2,axis=1)
                scores[j]+=n*(cw@error)
            nfamily+=n
        penalty=PENALTIES[int(np.argmin(scores))]
        model=fit(features(x,mode),y,w,penalty)
        predictions[mode]=np.mean([predict(model,features(v,mode)) for v in contexts],axis=0)
        receipt['models'][mode]={'penalty':penalty,'inner_mse':(scores/nfamily).tolist(),**{k:v.tolist() for k,v in model.items()}}
    if any(not np.isfinite(v).all() for v in predictions.values()):raise ValueError('nonfinite forecast')
    return evaluation,eligible,predictions,receipt


def score(effects,groups,stages):
    if set(stages)!={(c,f) for c in range(3) for f in range(5)}:
        raise ValueError('all 15 outer strata required')
    cells=defaultdict(dict);strata=[];availability=[]
    for (culture,heldfold),entry in stages.items():
        evaluation,eligible,predictions=entry
        expected=np.flatnonzero(np.asarray([fold(g) for g in groups])==heldfold)
        if not np.array_equal(evaluation,expected) or eligible.shape!=evaluation.shape or eligible.dtype!=bool:
            raise ValueError('outer population differs')
        if set(predictions)!=set(METHODS) or any(v.shape!=(int(eligible.sum()),3) or not np.isfinite(v).all() for v in predictions.values()):
            raise ValueError('incomplete or invalid forecasts')
        target=effects[evaluation,culture,3:]
        complete=np.isfinite(target).all(axis=1)
        scored=eligible&complete
        candidates=evaluation[eligible];observed=complete[eligible]
        y=target[eligible,1:]-target[eligible,:1]
        fam=np.asarray(groups)[candidates]
        observed_families=sorted(set(fam[observed]))
        if not observed_families:raise ValueError('unobserved outer stratum')
        strata.append({'culture':culture+1,'fold':heldfold,'population':len(evaluation),'input_eligible':int(eligible.sum()),
                       'complete_targets':int(scored.sum()),'unscored_predictions':int((eligible&~complete).sum()),'families':len(observed_families)})
        availability.append({'culture':culture+1,'fold':heldfold,'rows':evaluation.tolist(),'source_eligible':eligible.tolist(),
                             'target_support_by_time':np.isfinite(target).tolist()})
        for family in observed_families:
            mask=(fam==family)&observed
            cells[family][culture]={m:float(np.mean((predictions[m][mask]-y[mask])**2)) for m in METHODS}
    families=sorted(cells)
    loss=np.asarray([[np.mean([cell[m] for cell in cells[g].values()]) for m in METHODS] for g in families])
    gains=loss[:,1:]-loss[:,:1]
    rng=np.random.Generator(np.random.PCG64(20261018))
    boot=np.concatenate([gains[rng.integers(len(families),size=(1000,len(families)))].mean(axis=1) for _ in range(20)])
    interval=np.quantile(boot,[.01/8,1-.01/8],axis=0,method='linear')
    cultures=[]
    for c in range(3):
        rows=[cells[g][c] for g in families if c in cells[g]]
        means={m:float(np.mean([v[m] for v in rows])) for m in METHODS}
        cultures.append({'culture':c+1,'families':len(rows),'mse':means,'gain':{m:means[m]-means['level'] for m in CONTROLS}})
    contrasts=[{'control':m,'gain':float(gains[:,i].mean()),'interval':interval[:,i].tolist()} for i,m in enumerate(CONTROLS)]
    passed=all(c['interval'][0]>0 for c in contrasts) and all(all(v>0 for v in c['gain'].values()) for c in cultures)
    return {'status':'complete','decision':'conditional_on_supported_held_culture_measurements' if passed else 'not_passed',
            'families':len(families),'mse':dict(zip(METHODS,loss.mean(axis=0).tolist())),
            'contrasts':contrasts,'cultures':cultures,'strata':strata,'availability':availability,
            'family_cells':{g:{str(c+1):v for c,v in cells[g].items()} for g in families}}


def run():
    frozen=(HERE/'FREEZE.json').read_bytes();pins=json.loads(frozen)
    def verify():
        if (HERE/'FREEZE.json').read_bytes()!=frozen:raise ValueError('freeze changed')
        if pins['runtime']!={'python':sys.version.split()[0],'numpy':np.__version__}:raise ValueError('runtime differs')
        required=[HERE/name for name in ('PLAN.md','counts.py','forecast.py','COHORT.json','test_forecast.py','qualify.py','QUALIFICATION.json')]
        required.extend(counts.ADMISSION/name for name in ('GSE188264_dnaCounts.tsv.gz','GSE188264_rnaCounts.tsv.gz','GSE188264_dnaAnnot.tsv.gz','GSE188264_rnaAnnot.tsv.gz'))
        if any(str(path) not in pins['files'] for path in required):raise ValueError('freeze lacks required source')
        for path,digest in pins['files'].items():
            with Path(path).open('rb') as stream:
                if hashlib.file_digest(stream,'sha256').hexdigest()!=digest:raise ValueError('pinned source changed: '+path)
    verify();out=HERE/'result';out.mkdir();start=time.perf_counter()
    (out/'ATTEMPT.json').write_text(json.dumps({'freeze_sha256':hashlib.sha256(frozen).hexdigest()})+'\n')
    try:
        rows,effects,variant,parent,qc=counts.aggregate();groups=[r.group for r in rows]
        identities=[{'id':r.id,'group':r.group,'locus':r.locus,'fold':fold(r.group)} for r in rows]
        if identities!=json.loads((HERE/'COHORT.json').read_text()):raise ValueError('fixed cohort differs')
        np.savez(out/'COUNTS.npz',effects=effects,variant=variant,parent=parent)
        (out/'ROWS.json').write_text(json.dumps(identities,indent=2)+'\n')
        (out/'ADMISSION.json').write_text(json.dumps(qc,indent=2)+'\n')
        stages={};seals={}
        for c in range(3):
            for f in range(5):
                evaluation,eligible,predictions,receipt=outer(effects,groups,c,f)
                receipt['input_missing']=counts.missing_samples(rows,variant,parent,evaluation,[a for a in range(3) if a!=c],range(4),qc)
                stem=f'culture-{c+1}-fold-{f}'
                np.savez(out/(stem+'.npz'),rows=evaluation,eligible=eligible,**predictions)
                (out/(stem+'.json')).write_text(json.dumps(receipt,indent=2)+'\n')
                for suffix in ('.npz','.json'):
                    file=out/(stem+suffix);seals[file.name]=hashlib.sha256(file.read_bytes()).hexdigest()
                stages[c,f]=evaluation,eligible,predictions
        (out/'PREDICTION_SEALS.json').write_text(json.dumps(seals,indent=2,sort_keys=True)+'\n')
        result=score(effects,groups,stages)
        result['target_missing']=counts.missing_samples(rows,variant,parent,range(len(rows)),range(3),range(3,7),qc)
        for name,digest in seals.items():
            if hashlib.sha256((out/name).read_bytes()).hexdigest()!=digest:raise ValueError('prediction seal changed')
        verify();result['seconds']=time.perf_counter()-start
        (out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps({k:v for k,v in result.items() if k not in ('family_cells','availability')},indent=2))
    except Exception as exc:
        (out/'FAILURE.json').write_text(json.dumps({'error':repr(exc)})+'\n');raise


if __name__=='__main__':run()

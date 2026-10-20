"""No-refit transfer of frozen Siraj heads to the existing neural assay exports."""
import importlib.util
import json
from pathlib import Path
import signal
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
HUBER = BASE/'siraj-huber-joint-20261019'
SIX = BASE/'hnp-loso6-20261019'
CONTEXT = BASE/'context-loso6-20261019'
STUDIES = ('lee','mcafee','matoba','gse225815','gse221384','weiss')
METHODS = ('H','JL','K5','S2','S3')
BUDGETS = (10,50,200)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name,path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


parent = load('completed_huber_transfer_parent',HUBER/'run.py')
trial = parent.trial
linear = parent.previous.parent.context_features.linear
metrics = load('transfer_metrics',HERE/'metrics.py')


def legacy_pins(root):
    frozen = trial.read(root/'FREEZE.json')
    pins = {name:trial.pin(root/name) for name in frozen['files']}
    assert all(item['sha256'] == frozen['files'][name] for name,item in pins.items())
    return pins


def freeze():
    inherited = trial.read(HUBER/'FREEZE.json')
    trial.verify(inherited['pins'])
    pins = {name:item for name,item in inherited['pins'].items()
            if item['path'].endswith('.py') or name.startswith('motif/')}
    pins['huber_freeze'] = trial.pin(HUBER/'FREEZE.json')
    for name in ('STATE.json','RESULT.json','PREDICTION_FREEZE.json'):
        pins['huber/'+name] = trial.pin(HUBER/'experiment'/name)
    state = trial.read(pins['huber/STATE.json']['path'])
    seal = trial.read(pins['huber/PREDICTION_FREEZE.json']['path'])
    result = trial.read(pins['huber/RESULT.json']['path'])
    assert state['status'] == 'completed' and state['freeze'] == pins['huber_freeze']
    assert state['result'] == pins['huber/RESULT.json'] and seal['fit_ledger'] == state['fits']
    assert result['freeze'] == seal['experiment_freeze'] == pins['huber_freeze']
    assert result['prediction_freeze'] == pins['huber/PREDICTION_FREEZE.json'] and result['fits'] == len(state['fits'])
    head = next(r for r in state['fits'] if r['arm']=='huber' and r['fold'] is None)
    pins['H'] = head['artifact']
    pins['JL'] = inherited['pins']['initial_head/huber/None']
    assert head['initial'] == pins['JL'] and head['optimization']['converged'] is True
    pins['K5'] = inherited['pins']['k5_inner']
    pins['huber/head_review'] = trial.pin(HUBER/'reviews/heads/HEAD_REPLAY_RECEIPT.json')
    reviewed = trial.read(pins['huber/head_review']['path'])
    assert reviewed['status'] == 'passed' and reviewed['freeze'] == pins['huber_freeze']
    assert reviewed['heads']['huber-None'] == pins['H'] and reviewed['all_prediction_vectors_and_full_orders_verified'] is True
    for name in ('STATE.json','RESULT.json','PREDICTION_FREEZE.json'):
        assert reviewed['inputs'][name] == pins['huber/'+name]
    pins['huber/score_review'] = trial.pin(HUBER/'reviews/metrics/REVIEW.json')
    reviewed = trial.read(pins['huber/score_review']['path'])
    assert reviewed['status'] == 'passed' and reviewed['inputs']['freeze'] == pins['huber_freeze']
    for key,name in (('state','STATE.json'),('result','RESULT.json'),('seal','PREDICTION_FREEZE.json')):
        assert reviewed['inputs'][key] == pins['huber/'+name]
    for name,root in (('six',SIX),('context',CONTEXT)):
        pins[name+'/freeze'] = trial.pin(root/'FREEZE.json')
        pins.update({name+'/'+key:item for key,item in legacy_pins(root).items()})
        pins[name+'/result'] = trial.pin(root/'derived/RESULT.json')
        result = trial.read(pins[name+'/result']['path'])
        assert result['synthetic'] is False and result['null'] is False
        assert result['sealed_predictions_sha256'] == pins[name+'/derived/PREDICTIONS.npz']['sha256']
    for study in STUDIES:
        pins['outcome/'+study] = trial.pin(SIX/f'derived/outcomes_{study}.npz')
    for name in ('run.py','metrics.py','PLAN.md'):
        pins['transfer/'+name] = trial.pin(HERE/name)
    for kind in ('metrics','runner'):
        receipt = trial.read(HERE/f'reviews/{kind}/RECEIPT.json')
        assert receipt['status'] == 'passed'
        for path,sha in receipt['source_pins'].items():
            assert trial.pin(path)['sha256'] == sha,path
        for name in ('PLAN.md', 'metrics.py' if kind=='metrics' else 'run.py'):
            assert receipt['source_pins'][str(HERE/name)] == pins['transfer/'+name]['sha256']
        pins['review/'+kind] = trial.pin(HERE/f'reviews/{kind}/RECEIPT.json')
    plan = trial.read(inherited['pins']['plan_json']['path'])
    trial.verify(pins)
    trial.save(HERE/'FREEZE.json',dict(pins=pins,motif_ids=plan['motif_ids'],scale=plan['scale'],
        studies=list(STUDIES),methods=list(METHODS),budgets=list(BUDGETS),resamples=2000,seed=20261019,
        primary=[s for s in STUDIES if s!='gse225815'],biological_fits=0))


def inputs(frozen, study, old, context):
    rows = trial.read(frozen['pins'][f'six/derived/inputs_{study}.json']['path'])
    if len({r['id'] for r in rows}) != len(rows) or any(type(r['heldout']) is not bool for r in rows):
        raise ValueError('invalid legacy input identities')
    take = np.array([i for i,r in enumerate(rows) if r['heldout']],dtype=int)
    selected = [rows[i] for i in take]
    ids = np.array([r['id'] for r in selected]); families = np.array([r['family'] for r in selected])
    for archive in (old,context):
        assert np.array_equal(ids,archive[study+'/ids']) and np.array_equal(families,archive[study+'/family'])
    return rows,take,selected,ids,families


def predict():
    frozen = trial.read(HERE/'FREEZE.json');trial.verify(frozen['pins'])
    output = HERE/'prediction';output.mkdir(exist_ok=False)
    state = dict(status='started',freeze=trial.pin(HERE/'FREEZE.json'),studies=[],biological_fits=0)
    trial.save(output/'STATE.json',state);started=time.monotonic()
    try:
        profiles = trial.features.profiles()
        assert [p['id'] for p in profiles] == frozen['motif_ids']
        matrices = [trial.features.motif_features.log_odds(p['counts']) for p in profiles]
        heads = {name:dict(np.load(frozen['pins'][name]['path'],allow_pickle=False)) for name in ('H','JL','K5')}
        with np.load(frozen['pins']['six/derived/PREDICTIONS.npz']['path'],allow_pickle=False) as old, \
             np.load(frozen['pins']['context/derived/PREDICTIONS.npz']['path'],allow_pickle=False) as context, \
             np.load(frozen['pins']['context/derived/ctx_features.npz']['path'],allow_pickle=False) as cache:
            vectors = {};rankings = {}
            for study in STUDIES:
                record=dict(study=study,status='started');state['studies'].append(record)
                trial.save(output/'STATE.json',state,replace=True)
                rows,take,selected,ids,families = inputs(frozen,study,old,context)
                stored = cache[study]
                assert stored.shape == (len(rows),2*len(profiles))
                z = linear.deltas(selected,matrices,frozen['scale'])
                if not np.array_equal(z,stored[take,:len(profiles)]):
                    raise ValueError('native linear feature/cache mismatch: '+study)
                x = (trial.context_model.whole_counts([r['sequence'] for r in selected],5)
                     -trial.context_model.whole_counts([r['reference'] for r in selected],5))
                pred = {name:x@heads[name]['word']+z@heads[name]['motif'] for name in ('H','JL')}
                pred.update(K5=x@heads['K5']['coefficients'],S2=old[study+'/S2'].copy(),S3=context[study+'/S3'].copy())
                assert np.array_equal(pred['S2'],context[study+'/S2'])
                assert all(p.shape==(len(selected),) and np.isfinite(p).all() for p in pred.values())
                vectors.update({study+'/'+name:p for name,p in dict(ids=ids,family=families,**pred).items()})
                rankings[study] = {name:ids[np.lexsort((ids,-abs(p)))].tolist() for name,p in pred.items()}
                record.update(status='completed',rows=len(ids),families=len(set(families.tolist())),
                              native_motif_cache_bit_identical=True,seconds=time.monotonic()-started)
                trial.save(output/'STATE.json',state,replace=True)
                print(json.dumps(record),flush=True)
        trial.verify(frozen['pins'])
        trial.arrays(output/'PREDICTIONS.npz',**vectors)
        trial.save(output/'SEAL.json',dict(freeze=state['freeze'],predictions=trial.pin(output/'PREDICTIONS.npz'),
                                         studies=state['studies'],ranking_orders=rankings))
        state.update(status='completed',seal=trial.pin(output/'SEAL.json'),seconds=time.monotonic()-started)
        trial.save(output/'STATE.json',state,replace=True)
    except BaseException as error:
        state.update(status='failed',error=type(error).__name__+': '+str(error));trial.save(output/'STATE.json',state,replace=True)
        raise


def score():
    frozen = trial.read(HERE/'FREEZE.json');trial.verify(frozen['pins'])
    state = trial.read(HERE/'prediction/STATE.json');seal = trial.read(HERE/'prediction/SEAL.json')
    assert state['status']=='completed' and state['seal']==trial.pin(HERE/'prediction/SEAL.json')
    assert state['freeze']==seal['freeze']==trial.pin(HERE/'FREEZE.json') and state['studies']==seal['studies']
    trial.verify({'predictions':seal['predictions']})
    output=HERE/'scoring';output.mkdir(exist_ok=False)
    status=dict(status='started',seal=state['seal'],biological_fits=0);trial.save(output/'STATE.json',status)
    started=time.monotonic()
    try:
        historical = {name:trial.read(frozen['pins'][name+'/result']['path']) for name in ('six','context')}
        rng=np.random.Generator(np.random.PCG64(frozen['seed']));values={};populations={}
        with np.load(seal['predictions']['path'],allow_pickle=False) as saved:
            for study in STUDIES:
                rows=trial.read(frozen['pins'][f'six/derived/inputs_{study}.json']['path'])
                take=np.array([i for i,r in enumerate(rows) if r['heldout']],dtype=int)
                ids=saved[study+'/ids'];families=saved[study+'/family']
                assert ids.tolist()==[rows[i]['id'] for i in take] and families.tolist()==[rows[i]['family'] for i in take]
                with np.load(frozen['pins']['outcome/'+study]['path'],allow_pickle=False) as target:
                    assert target['E'].shape==(len(rows),)
                    y=target['E'][take]
                f=len(set(families.tolist()))
                counts=np.vstack((np.ones(f,dtype=int),rng.multinomial(f,np.full(f,1/f),size=frozen['resamples'])))
                pred={name:saved[study+'/'+name] for name in METHODS}
                for name,p in pred.items():assert seal['ranking_orders'][study][name]==ids[np.lexsort((ids,-abs(p)))].tolist()
                scored=metrics.evaluate(ids,families,y,pred,counts,BUDGETS)
                for name,source in (('S2','six'),('S3','context')):
                    np.testing.assert_allclose(scored[name]['r'][0],historical[source]['per_study_r'][study][name],rtol=1e-10,atol=1e-12)
                values[study]={}
                for name,record in scored.items():
                    for metric in ('r','mse'):values[study][metric+'/'+name]=record[metric]
                    for budget,vector in record['utility'].items():values[study][f'u{budget}/'+name]=vector
                for metric in ('r','mse',*[f'u{k}' for k in BUDGETS]):
                    for control in METHODS[1:]:
                        h,c=values[study][metric+'/H'],values[study][metric+'/'+control]
                        values[study][metric+'/gain_'+control]=None if h is None else ((c-h) if metric=='mse' else (h-c))
                populations[study]=dict(rows=len(ids),families=f,historical_target_order_reproduced=True)
        groups={};sets={'primary_five':frozen['primary'],'all_six':list(STUDIES)}
        for group,studies in sets.items():
            groups[group]={}
            for key in values[studies[0]]:
                eligible=[s for s in studies if values[s][key] is not None]
                averaged=np.mean([values[s][key] for s in eligible],axis=0)
                groups[group][key]=dict(studies=eligible,**metrics.summary(averaged))
        per_study={s:{key:None if vector is None else metrics.summary(vector) for key,vector in result.items()}
                   for s,result in values.items()}
        leave_one_out={s:{key:float(np.mean([values[t][key][0] for t in frozen['primary'] if t!=s]))
                          for key in ('r/H','r/gain_JL','r/gain_K5','r/gain_S2')} for s in frozen['primary']}
        result=dict(status='completed_exposed_transfer_diagnostic',groups=groups,per_study=per_study,
                    populations=populations,leave_one_study_out=leave_one_out,biological_fits=0,
                    original_siraj_acceptance=False,confirmation_evaluated=False,predictor_goal_complete=False,
                    freeze=trial.pin(HERE/'FREEZE.json'),prediction_seal=state['seal'],seconds=time.monotonic()-started)
        trial.verify(frozen['pins']);trial.verify({'predictions':seal['predictions']})
        trial.save(output/'RESULT.json',result);status.update(status='completed',result=trial.pin(output/'RESULT.json'))
        trial.save(output/'STATE.json',status,replace=True)
        print(json.dumps({k:v for k,v in result.items() if k not in ('groups','per_study')}),flush=True)
    except BaseException as error:
        status.update(status='failed',error=type(error).__name__+': '+str(error));trial.save(output/'STATE.json',status,replace=True)
        raise


if __name__=='__main__':
    action=sys.argv[1:]
    if action not in (['freeze'],['predict'],['score']):raise SystemExit('usage: run.py freeze|predict|score')
    if action!=['freeze']:
        def expired(signum,frame):raise TimeoutError('fixed 180-second transfer stage cap')
        signal.signal(signal.SIGALRM,expired);signal.setitimer(signal.ITIMER_REAL,180)
    {'freeze':freeze,'predict':predict,'score':score}[action[0]]()

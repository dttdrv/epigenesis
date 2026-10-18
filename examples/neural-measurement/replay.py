"""Replay the archived twelve-job RNA-repeat pilot without draws or fitting."""
import hashlib
import json
import math
from pathlib import Path

import numpy as np

import scoring


def read(path):
    return json.loads(path.read_text())


def replay(root):
    manifest = read(root/'MANIFEST.json')
    for name, expected in manifest['files'].items():
        path = Path(name)
        if path.is_absolute() or '..' in path.parts:
            raise ValueError('manifest path leaves the archive')
        if hashlib.sha256((root/path).read_bytes()).hexdigest() != expected:
            raise ValueError('archive identity differs: '+name)
    settings = manifest['calculation']
    selected = manifest['selected']
    omissions = range(settings['rna_repeats'])
    expected_jobs = {f'dev-{row["source_row"]}-{i}': (row, i) for row in selected for i in omissions}
    expected_contents = set(expected_jobs)|{'ATTEMPT.json', 'RESULT.json'}
    if {p.name for p in (root/'result').iterdir()} != expected_contents:
        raise ValueError('archive has incomplete or extra jobs')
    if len({row['component'] for row in selected}) != len(selected):
        raise ValueError('selected genomic components are not distinct')
    jobs, errors = [], []
    for name, (row, omission) in expected_jobs.items():
        folder = root/'result'/name
        inputs, saved = read(folder/'INPUT.json'), read(folder/'RESULT.json')
        if inputs['source'] != row or inputs['omission'] != omission:
            raise ValueError('source or omission differs')
        if any(saved[k] != v for k, v in inputs.items()):
            raise ValueError('saved input binding differs')
        counts = inputs['dna']+sum(inputs['training_rna'], [])+inputs['omitted_rna']
        if any(type(value) is not int or value <= 0 for value in counts):
            raise ValueError('counts violate the archived admission rule')
        total, observed = sum(inputs['omitted_rna']), inputs['omitted_rna'][1]
        if inputs['total'] != total:
            raise ValueError('omitted total differs')
        with np.load(folder/'moments.npz', allow_pickle=False) as arrays:
            if set(arrays.files) != {'sums', 'squares', 'baseline_cdf', 'cdf', 'cdf_lower', 'cdf_upper'}:
                raise ValueError('moment fields differ')
            for key in arrays.files:
                expected_shape = (total+1,) if key in ('sums', 'squares') else (total,)
                if arrays[key].shape != expected_shape:
                    raise ValueError('moment dimensions differ')
            bird = scoring.summarize(arrays['sums'], arrays['squares'], settings['raw_particles_per_job'],
                                     saved['envelope'], observed, jobs=len(expected_jobs),
                                     delta=settings['mc_delta'], tolerance=settings['log_tolerance'])
            for key in ('cdf', 'cdf_lower', 'cdf_upper'):
                errors.append(scoring.same(arrays[key], bird[key], key))
            pooled_mass = scoring.pooled_mass(inputs['training_rna'], total)
            errors.append(scoring.same(arrays['baseline_cdf'], scoring.prefix(pooled_mass)[:-1], 'baseline_cdf'))
        for key in ('raw_particles', 'mean_weight', 'effective_particles', 'loss', 'loss_interval', 'denominator_interval'):
            errors.append(scoring.same(saved['bird'][key], bird[key], key))
        pooled = scoring.pair_loss(pooled_mass, observed)
        gain = pooled-bird['loss']
        interval = [pooled-bird['loss_interval'][1], pooled-bird['loss_interval'][0]]
        for key, value in (('pooled_loss', pooled), ('gain', gain), ('gain_interval', interval)):
            errors.append(scoring.same(saved[key], value, key))
        jobs.append(dict(source=row, omission=omission, bird=bird['loss'], pooled=pooled,
                         gain=gain, gain_interval=interval))
    for row in selected:
        inputs = [read(root/'result'/f'dev-{row["source_row"]}-{i}'/'INPUT.json') for i in omissions]
        rna = [item['omitted_rna'] for item in inputs]
        for i, item in enumerate(inputs):
            if item['dna'] != inputs[0]['dna'] or item['training_rna'] != rna[:i]+rna[i+1:]:
                raise ValueError('omission count or training order differs')
    gain = math.fsum(job['gain'] for job in jobs)/len(jobs)
    interval = [math.fsum(job['gain_interval'][i] for job in jobs)/len(jobs) for i in (0, 1)]
    decision = 'bird_better' if interval[0] > 0 else 'pooled_better' if interval[1] < 0 else 'unresolved'
    saved = read(root/'result/RESULT.json')
    if (saved['status'], saved['jobs'], saved['raw_particles'], saved['decision']) != (
        'complete', len(jobs), settings['raw_particles_per_job']*len(jobs), decision
    ):
        raise ValueError('aggregate contract differs')
    errors.append(scoring.same(saved['gain'], gain, 'aggregate_gain'))
    errors.append(scoring.same(saved['gain_interval'], interval, 'aggregate_interval'))
    for component in saved['components']:
        selected_jobs = [job for job in jobs if job['source'] == component['source']]
        if len(selected_jobs) != settings['rna_repeats']:
            raise ValueError('component population differs')
        errors.append(scoring.same(component['gain'], math.fsum(job['gain'] for job in selected_jobs)/len(selected_jobs), 'component_gain'))
    bird = math.fsum(job['bird'] for job in jobs)/len(jobs)
    pooled = math.fsum(job['pooled'] for job in jobs)/len(jobs)
    return dict(status='saved_calculation_reproduced', jobs=len(jobs), bird_crps=bird,
                pooled_crps=pooled, gain=gain, relative_loss_reduction=100*gain/pooled,
                conditional_99_percent_simulation_interval=interval, archived_decision=decision,
                maximum_absolute_difference=max(errors), new_random_draws=0, new_model_fits=0,
                scope='fixed admitted TRAIN cohort; simulation uncertainty only; no biological generalization claim')


if __name__ == '__main__':
    print(json.dumps(replay(Path(__file__).resolve().parent), indent=2, allow_nan=False))

"""Frozen FamilyCode mutant-panel comparison, without fitting or tuning."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import random
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]))
from brainc._io import load_json_object, read_regular_file

IMPLEMENTATION_SHA256 = hashlib.sha256(read_regular_file(HERE/'benchmark.py',maximum_bytes=128*1024)).hexdigest()
spec = importlib.util.spec_from_file_location('benchmark_binding',HERE/'binding.py')
binding = importlib.util.module_from_spec(spec)
spec.loader.exec_module(binding)

PINS = {
    'protocol.json':'cb1747852f7e4a164d05d283389314f65dfdec822be7838160f0856965a4a6bc',
    'variant-mapping.json':'9ca227398575f3d2f11cc5fa2d95d1087955e8a2cfe93eae87351af02a27ff3b',
    'input-manifest.json':'1cf4a8c1f95344d4126b6814ac85374e161fd66cea3ea43bbd787e440dc7214b',
    'kock-panel.json':'c96d2f5f84634780ef17acc186632f643bf16ec799296a5c86b58cf7ae1ab500',
}
LOSSES = ('raw_model','raw_wt','raw_nearest','contrast_model','contrast_zero')
COMPARISONS = {'raw_vs_wt':('raw_wt','raw_model'),
               'raw_vs_nearest':('raw_nearest','raw_model'),
               'contrast_vs_zero':('contrast_zero','contrast_model')}


def prepare(matrix: list[list[float]], *, log_weights: bool) -> list[list[float]]:
    if (len(matrix) != 4 or not matrix[0] or any(len(row) != len(matrix[0]) for row in matrix)
            or any(not math.isfinite(v) for row in matrix for v in row)
            or (not log_weights and any(v < 0 for row in matrix for v in row))):
        raise ValueError('matrix must have four finite rectangular base rows')
    maxima = [max(column) for column in zip(*matrix)]
    if not log_weights and min(maxima) <= 0:
        raise ValueError('affinity column maximum must be positive')
    return [[max(binding.FLOOR,math.exp(v-top) if log_weights else v/top)
             for v,top in zip(row,maxima)] for row in matrix]


def matrix_mean(matrices: list[list[list[float]]]) -> list[list[float]]:
    return [[math.fsum(m[b][c] for m in matrices)/len(matrices)
             for c in range(len(matrices[0][0]))] for b in range(4)]


def centered(matrix: list[list[float]]) -> list[list[float]]:
    means = [math.fsum(math.log(v) for v in column)/4 for column in zip(*matrix)]
    return [[math.log(v)-mean for v,mean in zip(row,means)] for row in matrix]


def difference(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    return [[x-y for x,y in zip(ar,br)] for ar,br in zip(a,b)]


def mean_losses(records: list[dict]) -> dict:
    losses = {key:math.fsum(record['losses'][key] for record in records)/len(records) for key in LOSSES}
    if any(not math.isfinite(v) for v in losses.values()):
        raise ValueError('losses must be finite')
    return losses


def improvements(losses: dict) -> dict:
    return {name:losses[baseline]-losses[model] for name,(baseline,model) in COMPARISONS.items()}


def nearest_rows(rows: list[dict], factor: str, reference: str, mutant: str) -> tuple[list[int],int]:
    distances = [(i+1,sum(a != b for a,b in zip(row['protein_alignment'],mutant)))
                 for i,row in enumerate(rows) if row['name'].casefold() != factor.casefold()
                 and row['protein_alignment'] != reference]
    if not distances:
        raise ValueError('no non-self homeodomain available')
    minimum = min(d for _,d in distances)
    return [i for i,d in distances if d == minimum],minimum


def score(mutant: list, reference: list, measured: list, measured_wt: list, nearest: list) -> dict:
    wt_mean = matrix_mean(measured_wt)
    wt_log_mean = matrix_mean([centered(m) for m in measured_wt])
    predicted_contrast = difference(centered(mutant),centered(reference))
    zero = [[0.]*len(mutant[0]) for _ in range(4)]
    replicates = []
    for target in measured:
        observed = difference(centered(target),wt_log_mean)
        residuals = {
            'raw_model':difference(mutant,target), 'raw_wt':difference(wt_mean,target),
            'raw_nearest':difference(nearest,target),
            'contrast_model':difference(predicted_contrast,observed),
            'contrast_zero':difference(zero,observed),
        }
        losses = {key:math.fsum(v*v for row in matrix for v in row)/(4*len(matrix[0]))
                  for key,matrix in residuals.items()}
        replicates.append({'observed_contrast':observed,'residuals':residuals,'losses':losses})
    losses = mean_losses(replicates)
    return {'wt_mean':wt_mean,'wt_centered_log_mean':wt_log_mean,'nearest_mean':nearest,
            'predicted_contrast':predicted_contrast,'replicates':replicates,
            'losses':losses,'improvements':improvements(losses)}


def quantile(ordered: list[float], q: float) -> float:
    position = q*(len(ordered)-1)
    lo,hi = math.floor(position),math.ceil(position)
    return ordered[lo]+(position-lo)*(ordered[hi]-ordered[lo])


def summarize(records: list[dict], factor_order: list[str], resamples: int, seed: int) -> dict:
    if len(set(factor_order)) != len(factor_order) or set(factor_order) != {v['factor'] for v in records}:
        raise ValueError('factor cohort differs from frozen protocol')
    factors = []
    for factor in factor_order:
        variants = [v for v in records if v['factor'] == factor]
        losses = mean_losses(variants)
        factors.append({'factor':factor,'variants':[v['allele'] for v in variants],
                        'losses':losses,'improvements':improvements(losses)})
    losses = mean_losses(factors)
    count = len(factors)
    rng = random.Random(seed)
    bootstrap = {name:[] for name in COMPARISONS}
    for _ in range(resamples):
        indices = [rng.randrange(count) for _ in range(count)]
        for name in COMPARISONS:
            bootstrap[name].append(math.fsum(factors[i]['improvements'][name] for i in indices)/count)
    comparisons = {}
    for name,values in bootstrap.items():
        values.sort()
        interval = [quantile(values,q) for q in (.025,.975)]
        gains = [f['improvements'][name] for f in factors]
        leave_one_out = [math.fsum(gains[:i]+gains[i+1:])/(count-1) for i in range(count)]
        comparisons[name] = {'improvement':math.fsum(gains)/count,'interval95':interval,
                             'lower95_positive':interval[0] > 0,
                             'factor_signs':{'positive':sum(g > 0 for g in gains),
                                             'zero':sum(g == 0 for g in gains),
                                             'negative':sum(g < 0 for g in gains)},
                             'leave_one_factor_out_range':[min(leave_one_out),max(leave_one_out)]}
    return {'variants':len(records),'factors':factors,'losses':losses,'comparisons':comparisons,
            'all_comparisons_pass':all(c['lower95_positive'] for c in comparisons.values()),
            'bootstrap':{'resamples':resamples,'seed':seed,'unit':'factor','paired':True}}


def load_inputs(directory: Path | None = None) -> tuple[dict,dict,dict,list]:
    directory = HERE/'benchmark-data' if directory is None else directory
    data = {}
    for name,expected in PINS.items():
        raw = read_regular_file(directory/name,maximum_bytes=2*1024*1024)
        if hashlib.sha256(raw).hexdigest() != expected:
            raise ValueError(f'frozen benchmark input identity mismatch: {name}')
        data[name] = load_json_object(raw,name)
    protocol = data['protocol.json']
    bound = protocol['bound_inputs']
    if binding.IMPLEMENTATION_SHA256 != bound['prediction_implementation_sha256']:
        raise ValueError('prediction implementation differs from frozen protocol')
    if (binding.PINS['homeodomain.json'] != bound['executed_model_JSON_sha256']
            or binding.PINS['BLOSUM62.json'] != bound['portable_BLOSUM_JSON_sha256']):
        raise ValueError('model identity differs from frozen protocol')
    rows = binding.load_model()
    return protocol,data['variant-mapping.json'],data['kock-panel.json'],rows


def evaluate(directory: Path | None = None) -> dict:
    protocol,mapping,panel,rows = load_inputs(directory)
    observations = {r['sample']:dict(r,psam=prepare(r['original_log_weight_matrix'],log_weights=True))
                    for r in panel['replicates']}
    if len(observations) != len(panel['replicates']):
        raise ValueError('duplicate measured sample')
    if (len(rows) != protocol['model']['training_rows'] or panel['axes']['bases'] != list(binding.BASES)
            or any(r['DNA_positions'] != panel['axes']['dna_positions'] for r in rows)):
        raise ValueError('training or observation axes differ from protocol')
    variants,used = [],set()
    for entry in mapping['variants']:
        primary = entry['primary_clone_cohort']['status'] == 'eligible'
        cohort = 'primary' if primary else 'secondary_metadata'
        source = entry['primary_clone_cohort' if primary else 'secondary_metadata_cohort']
        if source['status'] not in ('eligible','eligible_descriptive_only'):
            raise ValueError(f"unsupported frozen case: {entry['allele']}")
        reference,mutant = source['reference_alignment'],source['mutant_alignment']
        position = entry['mutation']['alignment_column_one_based']
        if (len(reference) != protocol['model']['query_domain_length'] or len(mutant) != len(reference)
                or [i+1 for i,(a,b) in enumerate(zip(reference,mutant)) if a != b] != [position]
                or reference[position-1] != entry['mutation']['reference_residue']
                or mutant[position-1] != entry['mutation']['mutant_residue']):
            raise ValueError(f"invalid frozen mutation: {entry['allele']}")
        nearest,distance = nearest_rows(rows,entry['factor'],reference,mutant)
        expected = source['nearest_nonself_domain_baseline']
        if nearest != expected['selected_source_rows_one_based'] or distance != expected['minimum_distance']:
            raise ValueError(f"nearest baseline differs from frozen mapping: {entry['allele']}")
        nearest_mean = matrix_mean([prepare([rows[i-1]['relative_affinity_by_base'][b] for b in binding.BASES],
                                            log_weights=False) for i in nearest])
        measured = entry['measured_mutant_samples']
        wt = entry['measured_reference_samples']
        for sample_names,allele in ((measured,entry['allele']),(wt,entry['factor']+'-REF')):
            if not sample_names or set(sample_names) != {k for k,v in observations.items() if v['allele'] == allele}:
                raise ValueError(f'measured samples differ from frozen mapping: {allele}')
            used.update(sample_names)
        ref_prediction = binding._predict(rows,position,reference[position-1])
        mutant_prediction = binding._predict(rows,position,mutant[position-1])
        result = score([mutant_prediction['psam'][b] for b in binding.BASES],
                       [ref_prediction['psam'][b] for b in binding.BASES],
                       [observations[s]['psam'] for s in measured],[observations[s]['psam'] for s in wt],nearest_mean)
        for sample,replicate in zip(measured,result['replicates']):
            replicate['sample'] = sample
        variants.append(dict(result,allele=entry['allele'],factor=entry['factor'],cohort=cohort,status='evaluated',
                             source_provenance_status=entry['source_provenance_status'],
                             primary_exclusion_reasons=entry['primary_clone_cohort']['exclusion_reasons'],
                             source_row_correction=source.get('source_row_differs_from_original_script',False),
                             domain_mapping=source,mutation=entry['mutation'],
                             measured_reference_samples=wt,reference_prediction=ref_prediction,
                             mutant_prediction=mutant_prediction))
    if used != set(observations) or sorted(v['allele'] for v in variants) != protocol['cohorts']['accounting']['variant_ids']:
        raise ValueError('benchmark omitted frozen cases or observations')
    summaries = {}
    uncertainty = protocol['uncertainty']
    for name in ('primary','secondary_metadata'):
        selected = [v for v in variants if v['cohort'] == name]
        contract = protocol['cohorts'][name]
        if sorted(v['allele'] for v in selected) != contract['variant_ids']:
            raise ValueError('cohort differs from frozen protocol')
        summaries[name] = summarize(selected,contract['factor_order'],uncertainty['resamples'],uncertainty['seed'])
        summaries[name]['contributes_to_acceptance'] = name == 'primary'
        if name == 'secondary_metadata':
            del summaries[name]['all_comparisons_pass']
            for comparison in summaries[name]['comparisons'].values():
                del comparison['lower95_positive']
    return {'format':'epigenesis.binding-benchmark','version':1,'input_sha256':PINS,
            'implementation_sha256':IMPLEMENTATION_SHA256,'prediction_implementation_sha256':binding.IMPLEMENTATION_SHA256,
            'model_sha256':binding.PINS,'axes':panel['axes'],'observations':list(observations.values()),
            'variants':variants,'cohorts':summaries,
            'predictive_acceptance':summaries['primary']['all_comparisons_pass'],
            'acceptance_boundary':protocol['acceptance']['boundary'],
            'neural_development_acceptance':False,'limitations':protocol['limitations']}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True,help='absent JSON output path')
    args = parser.parse_args()
    try:
        report = evaluate()
        encoded = json.dumps(report,indent=2,allow_nan=False)+'\n'
        with args.output.open('x',encoding='utf-8') as output:
            output.write(encoded)
        print(json.dumps({'report':str(args.output),'predictive_acceptance':report['predictive_acceptance'],
                          'primary':report['cohorts']['primary']['comparisons']},indent=2))
    except (ValueError,OSError) as failure:
        parser.exit(2,f'benchmark: {failure}\n')


if __name__ == '__main__':
    main()

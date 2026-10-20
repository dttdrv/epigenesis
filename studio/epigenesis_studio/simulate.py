"""Live DNA-site to neural-dynamics runs through the repository's own backend, examples/neural-patterning.

Each run compiles the chosen site into a validated Development Module, reads the nucleotides
back, derives their binding-profile ratio, and drives the backend's four-gene model. Nothing here
models development: the equations, data and numerics are the backend's.
"""

import hashlib
import importlib.util
import math
import sys
import tempfile
from pathlib import Path

from brainc import _canonical as canonical

from .refusal import Refusal

EXAMPLE = Path(__file__).resolve().parents[2]/'examples/neural-patterning'
MODEL = 'epigenesis-neural-patterning-v1'
ENDPOINT = 'Nkx2.2 regulatory context: four-gene dynamics in model units'
MAX_PHASES = 10
MAX_TIME = 100.0  # model units; a resource limit so one request cannot run for hours
SITE_ID = 'nkx2_2_site'


def _load():
    path = EXAMPLE/'patterning.py'
    if not path.is_file():
        return None
    spec = importlib.util.spec_from_file_location('patterning', path)
    module = importlib.util.module_from_spec(spec)
    sys.modules['patterning'] = module
    spec.loader.exec_module(module)
    return module


patterning = _load()
IDENTITY = patterning and dict(
    model=MODEL, endpoint=ENDPOINT,
    manifest_sha256=canonical.digest(dict(implementation_sha256=hashlib.sha256(patterning.IMPLEMENTATION).hexdigest(),
                                          data_sha256=patterning.DATA_SHA256)))


def available():
    return patterning is not None


def _backend():
    if patterning is None:
        raise Refusal('sim_unavailable', 'the simulation backend examples/neural-patterning was not found')
    return patterning


def site_table():
    backend = _backend()
    data = backend.load_data()
    consensus = data['binding']['consensus']
    proteins = {}
    for protein in data['binding']['profiles']:
        rows = []
        for i, base in enumerate(consensus):
            row = {}
            for alt in 'ACGT':
                if alt != base:
                    try:
                        row[alt] = backend.profile(consensus[:i]+alt+consensus[i+1:], protein)['profile_ratio']
                    except ValueError:
                        row[alt] = None  # a published zero is not a measurement of zero binding
            rows.append(row)
        proteins[protein] = rows
    wild, engineered = data['observations']['nkx2_2_reporter']
    return dict(consensus=consensus, genome_strand=backend.reverse_complement(consensus), proteins=proteins,
                genes=data['model']['genes'], max_time=MAX_TIME, max_phases=MAX_PHASES,
                published=dict(site=backend.reverse_complement(engineered['site']), protein='GLI2',
                               reference_embryos=[wild['expressing'], wild['transgenic']],
                               variant_embryos=[engineered['expressing'], engineered['transgenic']]))


def normalize(request):
    backend = _backend()
    if not {'mode', 'site', 'protein', 'phases'} <= set(request) <= {'mode', 'site', 'protein', 'phases', 'initial_state'}:
        raise Refusal('unsupported', 'unsupported request')
    data = backend.load_data()
    consensus, genes = data['binding']['consensus'], data['model']['genes']
    site, protein, phases, initial = request['site'], request['protein'], request['phases'], request.get('initial_state')
    if type(protein) is not str or protein not in data['binding']['profiles']:
        raise Refusal('sim_protein', 'choose GLI2 or GLI3')
    if type(site) is not str or len(site) != len(consensus) or set(site)-set('ACGT'):
        raise Refusal('sim_site_shape', f'the site must be {len(consensus)} bases of A, C, G and T', length=len(consensus))
    changed = [i for i, (a, b) in enumerate(zip(site, consensus)) if a != b]
    if not changed:
        raise Refusal('sim_same', 'the site equals the reference; choose a substitution')
    if len(changed) > 1:
        raise Refusal('sim_site_multi', 'only one substitution from the reference site is supported')
    try:
        backend.profile(site, protein)
    except ValueError:
        raise Refusal('sim_site_zero', 'that substitution has a published zero weight, which the model cannot use') from None
    number = lambda v: type(v) in (int, float) and math.isfinite(v)
    if (not isinstance(phases, list) or not 1 <= len(phases) <= MAX_PHASES
            or any(not (isinstance(p, list) and len(p) == 2 and all(map(number, p)) and p[0] > 0 and 0 <= p[1] <= 1)
                   for p in phases)):
        raise Refusal('sim_phases', f'give 1 to {MAX_PHASES} phases, each with a positive duration and a signal from 0 to 1',
                      max=MAX_PHASES)
    if sum(p[0] for p in phases) > MAX_TIME:
        raise Refusal('sim_time', f'the total duration is at most {MAX_TIME:g} model units', max=MAX_TIME)
    if initial is None:
        initial = dict.fromkeys(genes, 0.0)
    if type(initial) is not dict or set(initial) != set(genes) or not all(number(v) and 0 <= v <= 1 for v in initial.values()):
        raise Refusal('sim_initial', 'the starting state needs all four genes, each from 0 to 1')
    return dict(mode='simulate', site=site, protein=protein, phases=[[float(d), float(s)] for d, s in phases],
                initial_state={gene: float(initial[gene]) for gene in genes})


def run(request):
    backend = _backend()
    data = backend.load_data()
    consensus, genes = data['binding']['consensus'], data['model']['genes']
    phases, initial, protein = request['phases'], request['initial_state'], request['protein']
    compiled = {}
    with tempfile.TemporaryDirectory() as scratch:
        for name, oriented in (('reference', consensus), ('variant', request['site']), ('reversion', consensus)):
            fasta = Path(scratch)/f'{name}.fasta'
            fasta.write_text(f'>{SITE_ID}\n{backend.reverse_complement(oriented)}\n')
            receipt = backend.compile_site(fasta, backend.selection(fasta, SITE_ID, 0, -1, protein), Path(scratch)/name)
            found = backend.consume(Path(scratch)/name, receipt)
            compiled[name] = dict(oriented_site=found['profile']['oriented_sequence'], ratio=found['profile']['profile_ratio'],
                                  changed_position=found['profile']['changed_position_1based'],
                                  selection_sha256=found['contract_sha256'], source_sha256=found['source_artifact_sha256'],
                                  development_sha256=found['development_bundle_sha256'])
    ratio = compiled['variant']['ratio']
    assay = lambda name, mode: backend.run_assay(phases, ratio=compiled[name]['ratio'], mode=mode, initial=initial)
    reference, variant, reporter, reverted = (assay('reference', 'endogenous'), assay('variant', 'endogenous'),
                                              assay('variant', 'reporter'), assay('reversion', 'endogenous'))
    gap = lambda a, b: max(abs(x-y) for ra, rb in zip(a['trajectory'], b['trajectory']) for x, y in zip(ra, rb))
    series = lambda a: dict(states=[list(row[1:]) for row in a['trajectory']], occupancy=list(a['reporter_occupancy']))
    return dict(mode='simulate', genes=genes, time=[row[0] for row in reference['trajectory']], signal=list(reference['signal']),
                site=dict(protein=protein, reference=consensus, variant=request['site'],
                          position=compiled['variant']['changed_position'], ratio=ratio),
                compiled=compiled, runs=dict(reference=series(reference), variant=series(variant), reporter=series(reporter)),
                controls=dict(reversion_max_difference=gap(reference, reverted), reporter_max_difference=gap(reference, reporter),
                              reversion_same_module=compiled['reversion']['development_sha256'] == compiled['reference']['development_sha256']),
                assumptions=data['assumptions'])

"""Live full-enhancer Variant-to-Circuit backend for Epigenesis Studio 1.7.

Reuses examples/neural-circuit/circuit.py (Boltzmann GLI partition sums over
Hallikas et al. 2006 weights + Herrera-Delgado et al. 2020 four-gene neural-tube
ODE + brainc Development Module compilation) and attaches passive reporter
comparisons (promoters, Kreimer 0-72h motif_temporal, Siraj Huber-linear).
"""

import hashlib
import importlib.util
import math
import sys
import tempfile
from pathlib import Path

from brainc import _canonical as canonical

from . import kreimer
from .model import predict as siraj_predict
from .refusal import Refusal

REPO = Path(__file__).resolve().parents[2]
CIRCUIT_DIR = REPO / 'examples/neural-circuit'
REG_DIR = REPO / 'examples/neural-regulation'
PAT_DIR = REPO / 'examples/neural-patterning'
MODEL = 'epigenesis-neural-circuit-v1'
ENDPOINT = 'Boltzmann GLI enhancer scan -> four-gene neural-tube bifurcation (model units)'
MIN_BASES = 9
MAX_BASES = 500
MAX_PHASES = 10
MAX_TIME = 100.0
DEFAULT_EPSILON = 1e-4
RADIUS = max(siraj_predict.load()[0]['k'], *siraj_predict.load()[0]['motif_widths']) - 1 if (Path(__file__).resolve().parent / 'model/model.json').is_file() else 18


def _load_backend():
    path = CIRCUIT_DIR / 'circuit.py'
    if not path.is_file() or not (PAT_DIR / 'patterning.py').is_file():
        return None
    if str(CIRCUIT_DIR) not in sys.path:
        sys.path.insert(0, str(CIRCUIT_DIR))
    spec = importlib.util.spec_from_file_location('neural_circuit_backend', path)
    module = importlib.util.module_from_spec(spec)
    sys.modules['neural_circuit_backend'] = module
    spec.loader.exec_module(module)
    return module


backend = _load_backend()
IDENTITY = backend and dict(
    model=MODEL,
    endpoint=ENDPOINT,
    manifest_sha256=canonical.digest(dict(
        implementation_sha256=hashlib.sha256(Path(backend.__file__).read_bytes()).hexdigest(),
        patterning_implementation_sha256=hashlib.sha256(backend.pat.IMPLEMENTATION).hexdigest(),
        data_sha256=backend.pat.DATA_SHA256,
    )),
)


def available():
    return backend is not None


def _backend():
    if backend is None:
        raise Refusal('circ_unavailable', 'the circuit backend (examples/neural-circuit) was not found')
    return backend


def example_enhancer():
    _backend()
    if str(REG_DIR) not in sys.path:
        sys.path.insert(0, str(REG_DIR))
    import regulation as reg
    rows, _ = reg.load_data()
    locus = 'chr15:90936956-90937127'
    reference = next(r.reference for r in rows if r.locus == locus)
    return dict(locus=locus, sequence=reference, length=len(reference), protein='GLI2',
                epsilon=DEFAULT_EPSILON, phases=[[40.0, 0.2]])


def _clean_sequence(raw):
    if not isinstance(raw, str):
        raise Refusal('not_text', 'sequence must be text', field='sequence', id=None)
    lines = raw.splitlines()
    records = sum(1 for line in lines if line.startswith('>'))
    if records > 1:
        raise Refusal('circ_multi_record', f'the FASTA has {records} records; scan one enhancer at a time', n=records)
    kept = [c for c in '\n'.join(l for l in lines if not l.startswith('>')) if not c.isspace() and not c.isdigit()]
    seq = ''.join(kept).upper()
    bad = sorted(set(seq) - set('ACGT'))
    if bad:
        raise Refusal('bad_chars', f'sequence contains characters other than A, C, G and T: {" ".join(bad)}',
                      chars=''.join(bad), field='sequence', id=None)
    if not MIN_BASES <= len(seq) <= MAX_BASES:
        raise Refusal('circ_len', f'enhancer sequence needs {MIN_BASES} to {MAX_BASES} bases',
                      min=MIN_BASES, max=MAX_BASES)
    return seq


def normalize(request):
    cb = _backend()
    allowed = {'mode', 'sequence', 'protein', 'phases', 'epsilon', 'initial_state'}
    needed = {'mode', 'sequence', 'protein', 'phases'}
    if not needed <= set(request) <= allowed:
        raise Refusal('unsupported', 'unsupported request')
    seq = _clean_sequence(request['sequence'])
    data = cb.pat.load_data()
    genes = data['model']['genes']
    protein = request['protein']
    if type(protein) is not str or protein not in data['binding']['profiles']:
        raise Refusal('circ_protein', 'choose GLI2 or GLI3')
    epsilon = request.get('epsilon', DEFAULT_EPSILON)
    number = lambda v: type(v) in (int, float) and math.isfinite(v)
    if not number(epsilon) or not 0.0 <= epsilon <= 1.0:
        raise Refusal('circ_epsilon', 'epsilon must be between 0 and 1')
    phases = request['phases']
    if (not isinstance(phases, list) or not 1 <= len(phases) <= MAX_PHASES
            or any(not (isinstance(p, list) and len(p) == 2 and all(map(number, p)) and p[0] > 0 and 0 <= p[1] <= 1)
                   for p in phases)):
        raise Refusal('circ_phases', f'give 1 to {MAX_PHASES} phases, each with a positive duration and a signal from 0 to 1',
                      max=MAX_PHASES)
    if sum(p[0] for p in phases) > MAX_TIME:
        raise Refusal('circ_time', f'the total duration is at most {MAX_TIME:g} model units', max=MAX_TIME)
    initial = request.get('initial_state')
    if initial is None:
        initial = dict.fromkeys(genes, 0.0)
    if type(initial) is not dict or set(initial) != set(genes) or not all(number(v) and 0 <= v <= 1 for v in initial.values()):
        raise Refusal('circ_initial', 'the starting state needs all four genes, each from 0 to 1')
    return dict(
        mode='circuit',
        sequence=seq,
        protein=protein,
        epsilon=float(epsilon),
        phases=[[float(d), float(s)] for d, s in phases],
        initial_state={gene: float(initial[gene]) for gene in genes},
    )


def _crop_rows(reference):
    rows = []
    for i, base in enumerate(reference):
        lo, hi = max(0, i - RADIUS), i + RADIUS + 1
        window = reference[lo:hi]
        for alt in 'ACGT':
            if alt != base:
                rows.append(dict(id=f'{i + 1}:{base}>{alt}', reference=window,
                                 sequence=window[: i - lo] + alt + window[i - lo + 1 :]))
    return rows


def run(request):
    cb = _backend()
    reference = request['sequence']
    protein = request['protein']
    epsilon = request['epsilon']
    phases = request['phases']
    initial_dict = request['initial_state']
    model_data = cb.pat.load_data()['model']
    genes = model_data['genes']
    initial_tuple = tuple(initial_dict[g] for g in genes)
    nkx_idx = genes.index('Nkx2.2')
    olig_idx = genes.index('Olig2')
    end_signal = phases[-1][1]

    consensus, ratios = cb._binding_tables(protein, epsilon)
    width = len(consensus)
    wt_info = cb.boltzmann_gli_affinity(reference, protein, epsilon)
    z_wt = wt_info['partition_sum']

    wt_assay = cb.run_circuit_assay(phases, ratio=1.0, mode='endogenous', initial=initial_dict)
    wt_end = tuple(wt_assay['trajectory'][-1][1:])
    wt_fate = cb.fate_label(wt_end, genes)
    wt_rep_end = wt_assay['reporter_occupancy'][-1]

    # Cache unique q values at 1e-8 resolution so a full 171-bp scan runs in <0.5s
    cache = {}
    edits = []
    for pos, base in enumerate(reference):
        for alt in 'ACGT':
            if alt == base:
                continue
            eid = f'{pos + 1}:{base}>{alt}'
            z_mut = z_wt + cb._delta_partition_single(reference, pos, alt, ratios, width)
            q = z_mut / z_wt if z_wt > 0 else 1.0
            key = round(q, 8)
            if key not in cache:
                end_state = cb.simulate_circuit(phases, ratio=q, initial=initial_tuple)[-1][1:]
                rep_end = cb.pat.occupancy(model_data, wt_end, end_signal, q)[nkx_idx]
                cache[key] = (end_state, rep_end)
            end_state, rep_end = cache[key]
            d_rep = rep_end - wt_rep_end
            d_nkx = end_state[nkx_idx] - wt_end[nkx_idx]
            d_olig = end_state[olig_idx] - wt_end[olig_idx]
            mut_fate = cb.fate_label(end_state, genes)
            fate_flipped = mut_fate != wt_fate
            if abs(d_rep) < 0.02 and fate_flipped:
                category = 'cryptic_bifurcation'
            elif abs(d_rep) >= 0.02 and abs(d_nkx) < 0.05 and not fate_flipped:
                category = 'canalized_buffered'
            elif fate_flipped:
                category = 'concordant_switch'
            else:
                category = 'neutral'
            edits.append(dict(
                id=eid, pos=pos, ref=base, alt=alt,
                affinity_ratio=float(q),
                reporter_delta_end=float(d_rep),
                nkx2_2_delta=float(d_nkx),
                olig2_delta=float(d_olig),
                wt_fate=wt_fate, mut_fate=mut_fate, fate_flipped=fate_flipped,
                amplification_ratio=float(abs(d_nkx) / max(abs(d_rep), 1e-12)),
                category=category,
            ))

    # Attach sequence-based passive reporter predictions using cropped local windows
    cropped = _crop_rows(reference)
    siraj_map = {r['id']: r['effect'] for r in siraj_predict.predict(cropped)}
    kreimer_map = {r['id']: r for r in kreimer.predict(cropped)}
    for e in edits:
        e['siraj_score'] = siraj_map[e['id']]
        e['kreimer_mean'] = kreimer_map[e['id']]['effect']
        e['kreimer_72h'] = kreimer_map[e['id']]['temporal'][-1]

    cryptic = sorted([e for e in edits if e['category'] == 'cryptic_bifurcation'],
                     key=lambda e: (-abs(e['nkx2_2_delta']), abs(e['reporter_delta_end']), e['id']))
    featured = cryptic[0] if cryptic else max(edits, key=lambda e: (abs(e['nkx2_2_delta']), -abs(e['reporter_delta_end'])))
    pos_f = featured['pos']
    variant_seq = reference[:pos_f] + featured['alt'] + reference[pos_f + 1 :]
    q_feat = cb.affinity_ratio(reference, variant_seq, protein, epsilon)

    # Compile and verify real brainc Development Module for reference -> featured switch, plus exact reversion
    with tempfile.TemporaryDirectory() as scratch:
        comp_dir = Path(scratch) / 'variant_bundle'
        receipt_var = cb.compile_circuit_bundle(reference, variant_seq, protein=protein, epsilon=epsilon, output=comp_dir)
        verified_var = cb.verify_circuit_bundle(comp_dir, receipt_var)
        rev_dir = Path(scratch) / 'reversion_bundle'
        receipt_rev = cb.compile_circuit_bundle(reference, reference, protein=protein, epsilon=epsilon, output=rev_dir)
        verified_rev = cb.verify_circuit_bundle(rev_dir, receipt_rev)

    mut_assay = cb.run_circuit_assay(phases, ratio=q_feat, mode='endogenous', initial=initial_dict)
    rep_assay = cb.run_circuit_assay(phases, ratio=q_feat, mode='reporter', initial=initial_dict)
    rev_assay = cb.run_circuit_assay(phases, ratio=verified_rev['affinity_ratio'], mode='endogenous', initial=initial_dict)
    gap = lambda a, b: max(abs(x - y) for ra, rb in zip(a['trajectory'], b['trajectory']) for x, y in zip(ra[1:], rb[1:]))
    series = lambda a: dict(states=[list(row[1:]) for row in a['trajectory']], occupancy=list(a['reporter_occupancy']))

    # 21-point spatial Shh gradient comparison
    positions = [i / 20 for i in range(21)]
    sp_ref = [list(cb.simulate_circuit([(100.0, math.exp(-x / model_data['gradient_length']))], ratio=1.0)[-1][1:]) for x in positions]
    sp_mut = [list(cb.simulate_circuit([(100.0, math.exp(-x / model_data['gradient_length']))], ratio=q_feat)[-1][1:]) for x in positions]
    p3_wt = sum(1 for y in sp_ref if y[nkx_idx] > 0.5)
    p3_mut = sum(1 for y in sp_mut if y[nkx_idx] > 0.5)

    # Pulse-chase GliR competitive derepression control
    hist_wt = cb.simulate_circuit([(5.0, 1.0), (15.0, 0.05)], ratio=1.0, mechanism='competitive')[-1][nkx_idx + 1]
    hist_mut = cb.simulate_circuit([(5.0, 1.0), (15.0, 0.05)], ratio=q_feat, mechanism='competitive')[-1][nkx_idx + 1]
    hist_act = cb.simulate_circuit([(5.0, 1.0), (15.0, 0.05)], ratio=q_feat, mechanism='activation-only')[-1][nkx_idx + 1]

    counts = {cat: sum(e['category'] == cat for e in edits)
              for cat in ('cryptic_bifurcation', 'concordant_switch', 'canalized_buffered', 'neutral')}
    ranked = sorted(edits, key=lambda e: (-abs(e['nkx2_2_delta']), -e['amplification_ratio'], e['id']))

    return dict(
        mode='circuit',
        genes=genes,
        length=len(reference),
        protein=protein,
        epsilon=epsilon,
        unique_affinity_solves=len(cache),
        wt_partition=wt_info,
        wt_end_state=dict(zip(genes, wt_end)),
        wt_fate=wt_fate,
        category_counts=counts,
        featured=featured,
        compiled=dict(variant_contract_sha256=verified_var['contract_sha256'],
                      reversion_contract_sha256=verified_rev['contract_sha256'],
                      variant_bundle_sha256=verified_var['bundle_sha256'],
                      reversion_bundle_sha256=verified_rev['bundle_sha256']),
        time=[row[0] for row in wt_assay['trajectory']],
        signal=list(wt_assay['signal']),
        runs=dict(reference=series(wt_assay), variant=series(mut_assay), reporter=series(rep_assay)),
        spatial=dict(positions=positions, reference=sp_ref, variant=sp_mut, p3_wt_points=p3_wt, p3_mut_points=p3_mut),
        pulse_chase=dict(wt_nkx=float(hist_wt), mut_competitive_nkx=float(hist_mut), mut_activation_only_nkx=float(hist_act)),
        controls=dict(
            reversion_max_difference=float(gap(wt_assay, rev_assay)),
            reporter_max_difference=float(gap(wt_assay, rep_assay)),
            reversion_ratio_exact=verified_rev['affinity_ratio'] == 1.0,
        ),
        edits=ranked,
    )

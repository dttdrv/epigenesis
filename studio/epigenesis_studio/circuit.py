"""Live full-enhancer Variant-to-Circuit backend for Epigenesis Studio 1.7.

Reuses examples/neural-circuit/circuit.py (Boltzmann TF partition sums over
Hallikas et al. 2006 GLI2/GLI3, 532 JASPAR 2024 motifs, or custom PWMs + arbitrary
N-gene GRN ODE presets + brainc Development Module compilation) and attaches passive
reporter comparisons (promoters, Kreimer 0-72h motif_temporal, Siraj Huber-linear).
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
ENDPOINT = 'Boltzmann TF enhancer scan -> N-gene regulatory bifurcation (model units)'
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
    cb = _backend()
    if str(REG_DIR) not in sys.path:
        sys.path.insert(0, str(REG_DIR))
    import regulation as reg
    rows, _ = reg.load_data()
    by_locus = {}
    for row in rows:
        if row.locus not in by_locus:
            by_locus[row.locus] = row.reference
    presets = [
        dict(locus=loc, sequence=by_locus[loc], length=len(by_locus[loc]))
        for loc in cb.CURATED_LOCI
        if loc in by_locus
    ]
    first = presets[0]
    circuits = [
        dict(
            name=name,
            genes=spec['genes'],
            default_target=spec['default_target'],
            default_protein=spec['default_protein'],
            default_signal=spec['default_signal'],
            edge_bindings=spec['edge_bindings'],
        )
        for name, spec in cb.list_circuit_presets().items()
    ]
    motifs = cb.list_available_motifs()
    return dict(
        locus=first['locus'], sequence=first['sequence'], length=first['length'],
        circuit=cb.DEFAULT_CIRCUIT, protein='GLI2', target=cb.DEFAULT_TARGET,
        epsilon=DEFAULT_EPSILON, phases=[[40.0, 0.2]], presets=presets,
        circuits=circuits, motifs=motifs,
    )


def _clean_sequence(raw, *, field='sequence', min_bases=MIN_BASES):
    if not isinstance(raw, str):
        raise Refusal('not_text', f'{field} must be text', field=field, id=None)
    lines = raw.splitlines()
    records = sum(1 for line in lines if line.startswith('>'))
    if records > 1:
        raise Refusal('circ_multi_record', f'the FASTA has {records} records; scan one enhancer at a time', n=records)
    kept = [c for c in '\n'.join(l for l in lines if not l.startswith('>')) if not c.isspace() and not c.isdigit()]
    seq = ''.join(kept).upper()
    bad = sorted(set(seq) - set('ACGT'))
    if bad:
        raise Refusal('bad_chars', f'{field} contains characters other than A, C, G and T: {" ".join(bad)}',
                      chars=''.join(bad), field=field, id=None)
    if not min_bases <= len(seq) <= MAX_BASES:
        raise Refusal('circ_len', f'enhancer sequence needs {min_bases} to {MAX_BASES} bases',
                      min=min_bases, max=MAX_BASES)
    return seq


def _resolve_circuit_choice(cb, raw_circuit):
    if raw_circuit is None:
        return cb.DEFAULT_CIRCUIT, cb.resolve_circuit(cb.DEFAULT_CIRCUIT)
    if isinstance(raw_circuit, str):
        if raw_circuit not in cb.SUPPORTED_PRESETS:
            raise Refusal('circ_circuit', f'choose one of {", ".join(cb.SUPPORTED_PRESETS)}')
        return raw_circuit, cb.resolve_circuit(raw_circuit)
    if isinstance(raw_circuit, dict):
        try:
            validated = cb.validate_circuit_spec(raw_circuit)
        except ValueError as err:
            raise Refusal('circ_circuit', str(err)) from err
        return validated, validated
    raise Refusal('circ_circuit', f'choose one of {", ".join(cb.SUPPORTED_PRESETS)} or a custom circuit object')


def normalize(request):
    cb = _backend()
    allowed = {'mode', 'sequence', 'variant_sequence', 'protein', 'target', 'circuit', 'phases', 'epsilon', 'initial_state'}
    needed = {'mode', 'sequence', 'phases'} if 'circuit' in request else {'mode', 'sequence', 'protein', 'phases'}
    if not needed <= set(request) <= allowed:
        raise Refusal('unsupported', 'unsupported request')

    circuit_val, model = _resolve_circuit_choice(cb, request.get('circuit'))
    genes = model['genes']

    epsilon = request.get('epsilon', DEFAULT_EPSILON)
    number = lambda v: type(v) in (int, float) and math.isfinite(v)
    if not number(epsilon) or not 0.0 <= epsilon <= 1.0:
        raise Refusal('circ_epsilon', 'epsilon must be between 0 and 1')

    protein = request.get('protein', model.get('default_protein', 'GLI2'))
    if type(protein) not in (str, dict):
        raise Refusal('circ_protein', 'choose GLI2, GLI3, or a JASPAR motif ID')
    try:
        prof = cb.load_motif_profile(protein, float(epsilon))
    except ValueError as err:
        raise Refusal('circ_protein', 'choose GLI2, GLI3, or a JASPAR motif ID') from err

    min_len = max(MIN_BASES, prof['width'])
    is_default_nt = (circuit_val == cb.DEFAULT_CIRCUIT)
    if not is_default_nt and model.get('edge_bindings'):
        for eb in model['edge_bindings']:
            try:
                ep = cb.load_motif_profile(eb['protein'], float(epsilon))
                min_len = max(min_len, ep['width'])
            except ValueError as err:
                raise Refusal('circ_circuit', str(err)) from err

    seq = _clean_sequence(request['sequence'], field='sequence', min_bases=min_len)
    variant_seq = None
    if request.get('variant_sequence') is not None:
        raw_var = request['variant_sequence']
        if type(raw_var) is not str:
            raise Refusal('not_text', 'variant_sequence must be text', field='variant_sequence', id=None)
        if raw_var.strip():
            variant_seq = _clean_sequence(raw_var, field='variant_sequence', min_bases=min_len)

    default_target = cb.DEFAULT_TARGET if is_default_nt else model['default_target']
    target = request.get('target', default_target)
    if is_default_nt:
        if type(target) is not str or target not in cb.SUPPORTED_TARGETS:
            raise Refusal('circ_target', 'choose Nkx2.2, Olig2 or Both')
    else:
        valid_targets = set(genes) | {'All'}
        if {'Nkx2.2', 'Olig2'} <= set(genes):
            valid_targets.add('Both')
        if type(target) is not str or target not in valid_targets:
            raise Refusal('circ_target', f'choose one of {", ".join(sorted(valid_targets))}')

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
        raise Refusal('circ_initial', f'the starting state needs all genes ({", ".join(genes)}), each from 0 to 1')

    out = dict(
        mode='circuit',
        sequence=seq,
        protein=protein,
        target=target,
        epsilon=float(epsilon),
        phases=[[float(d), float(s)] for d, s in phases],
        initial_state={gene: float(initial[gene]) for gene in genes},
    )
    if 'circuit' in request:
        out['circuit'] = circuit_val
    if variant_seq is not None:
        out['variant_sequence'] = variant_seq
    return out


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
    variant_override = request.get('variant_sequence')
    circuit_val = request.get('circuit', cb.DEFAULT_CIRCUIT)
    c_arg = None if circuit_val == cb.DEFAULT_CIRCUIT else circuit_val
    protein = request['protein']
    target = request.get('target', cb.DEFAULT_TARGET)
    epsilon = request['epsilon']
    phases = request['phases']
    initial_dict = request['initial_state']

    model_peek = cb._pat_data()['model'] if c_arg is None else cb.resolve_circuit(c_arg)
    multi_edge = (
        c_arg is not None
        and bool(model_peek.get('edge_bindings'))
        and protein == model_peek.get('default_protein')
    )
    ctx = cb._resolve_context(
        c_arg,
        None if multi_edge else protein,
        target,
        phases[0][1],
        None,
        None,
    )
    model_data = ctx['model']
    genes = ctx['genes']
    bindings = ctx['bindings']
    bound_genes = sorted(bindings)
    initial_tuple = tuple(initial_dict[g] for g in genes)
    rep_idx = cb._target_gene_index(genes, ctx['target'])
    has_nkx_olig = 'Nkx2.2' in genes and 'Olig2' in genes
    nkx_idx = genes.index('Nkx2.2') if has_nkx_olig else rep_idx
    olig_idx = genes.index('Olig2') if has_nkx_olig else (1 if len(genes) > 1 and rep_idx != 1 else (0 if len(genes) == 1 else (rep_idx + 1) % len(genes)))
    end_signal = phases[-1][1]
    fates_cfg = model_data.get('fates')
    kw = dict(circuit=ctx['circuit_arg'], roles=ctx['roles'], mechanisms=ctx['mechanisms'])

    tf_tables = {}
    for lbl, spec in ctx['unique_specs']:
        prof = cb.load_motif_profile(spec, epsilon)
        wt_part = cb.boltzmann_affinity(reference, spec, epsilon)
        tf_tables[lbl] = dict(width=prof['width'], ratios=prof['ratios'], wt_part=wt_part, z_wt=wt_part['partition_sum'])
    primary_label = ctx['primary_label']
    wt_info = tf_tables[primary_label]['wt_part']

    wt_ratios = dict.fromkeys(bound_genes, 1.0)
    wt_assay = cb.run_circuit_assay(phases, ratio=wt_ratios, mode='endogenous', target=ctx['target'], initial=initial_dict, **kw)
    wt_end = tuple(wt_assay['trajectory'][-1][1:])
    wt_fate = cb.fate_label(wt_end, genes, fates=fates_cfg)
    wt_rep_end = wt_assay['reporter_occupancy'][-1]

    # Cache unique ratio vectors at 1e-8 resolution so a full 171-bp scan runs in <0.5s
    cache = {}
    edits = []
    for pos, base in enumerate(reference):
        for alt in 'ACGT':
            if alt == base:
                continue
            eid = f'{pos + 1}:{base}>{alt}'
            q_by_tf = {}
            for lbl, tbl in tf_tables.items():
                z_mut = tbl['z_wt'] + cb._delta_partition_single(reference, pos, alt, tbl['ratios'], tbl['width'])
                q_by_tf[lbl] = z_mut / tbl['z_wt'] if tbl['z_wt'] > 0 else 1.0
            ratio_map = {g: q_by_tf[cb._motif_spec_label(bindings[g]['protein'])] for g in bound_genes}
            q_primary = q_by_tf[primary_label]
            key = tuple(round(ratio_map[g], 8) for g in bound_genes)
            if key not in cache:
                end_state = cb.simulate_circuit(phases, ratios=ratio_map, target=ctx['target'], initial=initial_tuple, **kw)[-1][1:]
                rep_end = cb.occupancy_ratios(model_data, wt_end, end_signal, ratio_map, roles=ctx['roles'], mechanisms=ctx['mechanisms'])[rep_idx]
                cache[key] = (end_state, rep_end)
            end_state, rep_end = cache[key]
            d_rep = rep_end - wt_rep_end
            gene_deltas = {g: float(end_state[i] - wt_end[i]) for i, g in enumerate(genes)}
            target_delta = float(end_state[rep_idx] - wt_end[rep_idx])
            mut_fate = cb.fate_label(end_state, genes, fates=fates_cfg)
            fate_flipped = mut_fate != wt_fate
            category = cb._classify_edit(d_rep, target_delta, fate_flipped)
            entry = dict(
                id=eid, pos=pos, ref=base, alt=alt,
                affinity_ratio=float(q_primary),
                affinity_ratios={k: float(v) for k, v in q_by_tf.items()},
                reporter_delta_end=float(d_rep),
                gene_deltas=gene_deltas,
                target_delta=float(target_delta),
                wt_fate=wt_fate, mut_fate=mut_fate, fate_flipped=fate_flipped,
                amplification_ratio=float(abs(target_delta) / max(abs(d_rep), 1e-12)),
                category=category,
            )
            if has_nkx_olig:
                entry['nkx2_2_delta'] = gene_deltas['Nkx2.2']
                entry['olig2_delta'] = gene_deltas['Olig2']
            edits.append(entry)

    # Attach sequence-based passive reporter predictions using cropped local windows
    cropped = _crop_rows(reference)
    siraj_map = {r['id']: r['effect'] for r in siraj_predict.predict(cropped)}
    kreimer_map = {r['id']: r for r in kreimer.predict(cropped)}
    for e in edits:
        e['siraj_score'] = siraj_map[e['id']]
        e['kreimer_mean'] = kreimer_map[e['id']]['effect']
        e['kreimer_72h'] = kreimer_map[e['id']]['temporal'][-1]

    cryptic = sorted([e for e in edits if e['category'] == 'cryptic_bifurcation'],
                     key=lambda e: (-abs(e['target_delta']), abs(e['reporter_delta_end']), e['id']))
    featured = cryptic[0] if cryptic else max(edits, key=lambda e: (abs(e['target_delta']), -abs(e['reporter_delta_end']), e['id']))
    pos_f = featured['pos']
    auto_variant_seq = reference[:pos_f] + featured['alt'] + reference[pos_f + 1 :]
    second_site = cb._find_second_site_rescue(
        reference, dict(featured, variant=auto_variant_seq),
        None if multi_edge else protein, epsilon,
        target=ctx['target'], phases=phases, initial=initial_dict,
        circuit=ctx['circuit_arg'], roles=ctx['roles'], mechanisms=ctx['mechanisms'],
    )
    if second_site is not None:
        second_site = {k: v for k, v in second_site.items() if k != 'double_variant'}

    variant_seq = variant_override if variant_override is not None else auto_variant_seq
    multi_var = cb.multi_tf_affinity(reference, variant_seq, dict(ctx['unique_specs']), epsilon)
    var_q_by_tf = multi_var['ratios']
    q_feat = var_q_by_tf[primary_label]
    var_ratio_map = {g: var_q_by_tf[cb._motif_spec_label(b['protein'])] for g, b in bindings.items()}

    # Compile and verify real brainc Development Module (v2 for GLI2/GLI3 neural-tube, v3 for JASPAR/multi-TF/custom circuit)
    motif_list = [spec for _, spec in ctx['unique_specs']]
    with tempfile.TemporaryDirectory() as scratch:
        comp_dir = Path(scratch) / 'variant_bundle'
        receipt_var = cb.compile_circuit_bundle(
            reference, variant_seq, protein=ctx['primary_protein'], epsilon=epsilon,
            output=comp_dir, motif_specs=motif_list, circuit=ctx['circuit_arg'],
        )
        verified_var = cb.verify_circuit_bundle(comp_dir, receipt_var)
        rev_dir = Path(scratch) / 'reversion_bundle'
        receipt_rev = cb.compile_circuit_bundle(
            reference, reference, protein=ctx['primary_protein'], epsilon=epsilon,
            output=rev_dir, motif_specs=motif_list, circuit=ctx['circuit_arg'],
        )
        verified_rev = cb.verify_circuit_bundle(rev_dir, receipt_rev)

    mut_assay = cb.run_circuit_assay(phases, ratio=var_ratio_map, mode='endogenous', target=ctx['target'], initial=initial_dict, **kw)
    rep_assay = cb.run_circuit_assay(phases, ratio=var_ratio_map, mode='reporter', target=ctx['target'], initial=initial_dict, **kw)
    rev_assay = cb.run_circuit_assay(phases, ratio=wt_ratios, mode='endogenous', target=ctx['target'], initial=initial_dict, **kw)
    gap = lambda a, b: max(abs(x - y) for ra, rb in zip(a['trajectory'], b['trajectory']) for x, y in zip(ra[1:], rb[1:]))
    series = lambda a: dict(states=[list(row[1:]) for row in a['trajectory']], occupancy=list(a['reporter_occupancy']))

    # 21-point spatial morphogen gradient comparison
    grad_len = float(model_data['gradient_length'])
    positions = [i / 20 for i in range(21)]
    sp_ref = [list(cb.simulate_circuit([(100.0, math.exp(-x / grad_len))], ratios=wt_ratios, target=ctx['target'], **kw)[-1][1:]) for x in positions]
    sp_mut = [list(cb.simulate_circuit([(100.0, math.exp(-x / grad_len))], ratios=var_ratio_map, target=ctx['target'], **kw)[-1][1:]) for x in positions]
    p3_wt = sum(1 for y in sp_ref if y[rep_idx] > 0.5)
    p3_mut = sum(1 for y in sp_mut if y[rep_idx] > 0.5)

    # Pulse-chase competitive derepression control
    hist_wt = cb.simulate_circuit([(5.0, 1.0), (15.0, 0.05)], ratios=wt_ratios, target=ctx['target'], mechanism='competitive', **kw)[-1][rep_idx + 1]
    hist_mut = cb.simulate_circuit([(5.0, 1.0), (15.0, 0.05)], ratios=var_ratio_map, target=ctx['target'], mechanism='competitive', **kw)[-1][rep_idx + 1]
    hist_act = cb.simulate_circuit([(5.0, 1.0), (15.0, 0.05)], ratios=var_ratio_map, target=ctx['target'], mechanism='activation-only', **kw)[-1][rep_idx + 1]

    counts = {cat: sum(e['category'] == cat for e in edits)
              for cat in ('cryptic_bifurcation', 'concordant_switch', 'canalized_buffered', 'neutral')}
    ranked = sorted(edits, key=lambda e: (-abs(e['target_delta']), -e['amplification_ratio'], e['id']))

    custom_summary = None
    if variant_override is not None:
        var_end = tuple(mut_assay['trajectory'][-1][1:])
        var_rep_end = cb.occupancy_ratios(model_data, wt_end, end_signal, var_ratio_map, roles=ctx['roles'], mechanisms=ctx['mechanisms'])[rep_idx]
        d_rep_c = var_rep_end - wt_rep_end
        c_gene_deltas = {g: float(var_end[i] - wt_end[i]) for i, g in enumerate(genes)}
        target_delta_c = float(var_end[rep_idx] - wt_end[rep_idx])
        mut_fate_c = cb.fate_label(var_end, genes, fates=fates_cfg)
        pair_row = [dict(id='custom', reference=reference, sequence=variant_override)]
        siraj_c = siraj_predict.predict(pair_row)[0]['effect']
        kreimer_c = kreimer.predict(pair_row)[0]
        custom_summary = dict(
            edit_kind=cb._edit_kind(reference, variant_override),
            reference_length=len(reference),
            variant_length=len(variant_override),
            affinity_ratio=float(q_feat),
            affinity_ratios={k: float(v) for k, v in var_q_by_tf.items()},
            reporter_delta_end=float(d_rep_c),
            gene_deltas=c_gene_deltas,
            target_delta=float(target_delta_c),
            wt_fate=wt_fate,
            mut_fate=mut_fate_c,
            fate_flipped=mut_fate_c != wt_fate,
            amplification_ratio=float(abs(target_delta_c) / max(abs(d_rep_c), 1e-12)),
            category=cb._classify_edit(d_rep_c, target_delta_c, mut_fate_c != wt_fate),
            siraj_score=siraj_c,
            kreimer_mean=kreimer_c['effect'],
            kreimer_72h=kreimer_c['temporal'][-1],
        )
        if has_nkx_olig:
            custom_summary['nkx2_2_delta'] = c_gene_deltas['Nkx2.2']
            custom_summary['olig2_delta'] = c_gene_deltas['Olig2']

    out = dict(
        mode='circuit',
        circuit=ctx['circuit_name'],
        genes=genes,
        target=ctx['target'],
        target_gene=genes[rep_idx],
        alt_gene=genes[olig_idx],
        length=len(reference),
        protein=primary_label,
        proteins=[lbl for lbl, _ in ctx['unique_specs']],
        epsilon=epsilon,
        unique_affinity_solves=len(cache),
        wt_partition=wt_info,
        wt_partitions={lbl: tbl['wt_part'] for lbl, tbl in tf_tables.items()},
        wt_end_state=dict(zip(genes, wt_end)),
        wt_fate=wt_fate,
        category_counts=counts,
        featured=featured,
        second_site_rescue=second_site,
        compiled=dict(variant_contract_sha256=verified_var['contract_sha256'],
                      reversion_contract_sha256=verified_rev['contract_sha256'],
                      variant_bundle_sha256=verified_var['bundle_sha256'],
                      reversion_bundle_sha256=verified_rev['bundle_sha256'],
                      reference_length=verified_var['reference_length'],
                      variant_length=verified_var['variant_length']),
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
    if custom_summary is not None:
        out['custom_variant'] = custom_summary
    return out

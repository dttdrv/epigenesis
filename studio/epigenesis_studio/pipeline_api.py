"""Studio backend for unified developmental pipeline (`develop`), inverse DNA compiler (`synthesize`),
and 532-TF human biophysical discovery (`discover`).

All responses are deterministic and free of host-specific paths or wall-clock timings
so Epigenesis Studio receipts replay identically (`server.verify(receipt)["ok"] == True`).
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import tempfile
from pathlib import Path

from brainc import _canonical as canonical
from .refusal import Refusal

REPO = Path(__file__).resolve().parents[2]
PIPE_DIR = REPO / 'examples/developmental-pipeline'
REG_DIR = REPO / 'examples/neural-regulation'

DEVELOP_MODEL = 'epigenesis-developmental-pipeline-v1'
DEVELOP_ENDPOINT = 'Multi-locus DNA + VRS 2.x -> chromatin/coding -> GRN uncertainty -> 3D CX3D+LIF'
SYNTH_MODEL = 'epigenesis-inverse-compiler-v1'
SYNTH_ENDPOINT = 'Inverse developmental DNA synthesis (fate reprogramming & epistatic 3D rescue)'
DISCOVER_MODEL = 'epigenesis-biophysical-discovery-v1'
DISCOVER_ENDPOINT = '532-TF Boltzmann thermodynamics discovery on GSE188264 (0-72h)'


def _load(name: str, filename: str):
    path = PIPE_DIR / filename
    if not path.is_file():
        return None
    if str(PIPE_DIR) not in sys.path:
        sys.path.insert(0, str(PIPE_DIR))
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


pipeline_mod = _load('epigenesis_dev_pipeline_backend', 'pipeline.py')
synth_mod = _load('epigenesis_inverse_synth_backend', 'synthesize.py')
discover_mod = _load('epigenesis_biophys_discover_backend', 'discover.py')


def _file_sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


DEVELOP_IDENTITY = pipeline_mod and dict(
    model=DEVELOP_MODEL,
    endpoint=DEVELOP_ENDPOINT,
    manifest_sha256=canonical.digest(dict(pipeline_sha256=_file_sha(PIPE_DIR / 'pipeline.py'))),
)
SYNTH_IDENTITY = synth_mod and dict(
    model=SYNTH_MODEL,
    endpoint=SYNTH_ENDPOINT,
    manifest_sha256=canonical.digest(dict(
        pipeline_sha256=_file_sha(PIPE_DIR / 'pipeline.py'),
        synthesize_sha256=_file_sha(PIPE_DIR / 'synthesize.py'),
    )),
)
DISCOVER_IDENTITY = discover_mod and dict(
    model=DISCOVER_MODEL,
    endpoint=DISCOVER_ENDPOINT,
    manifest_sha256=canonical.digest(dict(
        pipeline_sha256=_file_sha(PIPE_DIR / 'pipeline.py'),
        discover_sha256=_file_sha(PIPE_DIR / 'discover.py'),
    )),
)


def available() -> bool:
    return pipeline_mod is not None and synth_mod is not None


def _strip_host_fields(prov: dict) -> dict:
    clean = json.loads(json.dumps(prov))
    rc = clean.get('rust_cross_check')
    if isinstance(rc, dict):
        rc.pop('binary', None)
    return clean


def _strip_cx3d_paths(morph: dict) -> dict:
    clean = json.loads(json.dumps(morph))
    for b in (clean.get('branches') or {}).values():
        cx = b.get('construction')
        if isinstance(cx, dict):
            cx.pop('command', None)
    return clean


def presets_info() -> dict:
    if not available():
        raise Refusal('pipe_unavailable', 'developmental pipeline backend not found')
    return dict(
        develop_presets=['neural-tube-multi-locus', 'tri-lineage-chromatin'],
        rescue_presets=[
            dict(locus='chr18:56475717-56475888', breaking_edit='165:G>T', circuit='neural-tube', protein='GLI2', target_gene='Nkx2.2', signal=0.20),
            dict(locus='chr1:175099416-175099587', breaking_edit='149:G>T', circuit='neural-tube', protein='GLI2', target_gene='Nkx2.2', signal=0.20),
        ],
        reprogram_presets=[
            dict(locus='chr15:90936956-90937127', target_fate='Olig2', circuit='neural-tube', protein='GLI2', target_gene='Nkx2.2', signal=0.20),
            dict(locus='chr18:56475717-56475888', target_fate='Olig2', circuit='neural-tube', protein='GLI2', target_gene='Nkx2.2', signal=0.20),
        ],
        discover_available=discover_mod is not None,
    )


def normalize_develop(request: dict) -> dict:
    if pipeline_mod is None:
        raise Refusal('pipe_unavailable', 'developmental pipeline backend not found')
    allowed = {'mode', 'preset', 'accessibility', 'include_coding_truncation', 'samples', 'run_cx3d'}
    if not set(request) <= allowed:
        raise Refusal('unsupported', 'unsupported develop request fields')
    preset = request.get('preset', 'neural-tube-multi-locus')
    if preset not in ('neural-tube-multi-locus', 'tri-lineage-chromatin'):
        raise Refusal('bad_preset', 'choose neural-tube-multi-locus or tri-lineage-chromatin')
    acc = request.get('accessibility', 1.0)
    if type(acc) not in (int, float) or not 0.0 <= float(acc) <= 1.0:
        raise Refusal('bad_accessibility', 'chromatin accessibility must be between 0 and 1')
    samples = request.get('samples', 16)
    if type(samples) is not int or not 0 <= samples <= 100:
        raise Refusal('bad_samples', 'uncertainty samples must be an integer from 0 to 100')
    inc_coding = bool(request.get('include_coding_truncation', False))
    run_cx3d = bool(request.get('run_cx3d', False))
    return dict(
        mode='develop',
        preset=preset,
        accessibility=float(acc),
        include_coding_truncation=inc_coding,
        samples=samples,
        run_cx3d=run_cx3d,
    )


def run_develop(request: dict) -> dict:
    if pipeline_mod is None:
        raise Refusal('pipe_unavailable', 'developmental pipeline backend not found')
    preset = pipeline_mod.builtin_preset(request['preset'])
    spec = json.loads(json.dumps(preset['spec']))
    spec['loci'][0]['accessibility'] = request['accessibility']
    spec['uncertainty']['samples'] = request['samples']
    variants = list(preset['variants'])
    if request['include_coding_truncation'] and request['preset'] == 'neural-tube-multi-locus':
        variants.append(dict(id='pax6_early_stop_C4T', record_id='orf_pax6', start=3, end=4, ref='C', alt='T'))

    with tempfile.TemporaryDirectory() as tmp:
        troot = Path(tmp)
        fasta = troot / 'loci.fasta'
        pipeline_mod._write_fasta(preset['records'], fasta)
        rep = pipeline_mod.run_pipeline(
            fasta,
            spec,
            variants,
            troot / 'out',
            run_cx3d=request['run_cx3d'],
            cx3d_duration=1.2,
            cx3d_seed=7,
            export_bundles=False,
        )
    return dict(
        mode='develop',
        program_name=rep['program_name'],
        circuit=rep['circuit'],
        genes=rep['genes'],
        biological_acceptance=rep['biological_acceptance'],
        provenance=_strip_host_fields(rep['provenance']),
        variants=rep['variants'],
        loci=rep['loci'],
        nominal=rep['nominal'],
        uncertainty_ensemble=rep['uncertainty_ensemble'],
        controls=rep['controls'],
        morphogenesis=_strip_cx3d_paths(rep['morphogenesis']),
    )


def normalize_synthesize(request: dict) -> dict:
    if synth_mod is None:
        raise Refusal('pipe_unavailable', 'inverse compiler backend not found')
    allowed = {'mode', 'goal', 'locus', 'breaking_edit', 'target_fate', 'samples', 'run_cx3d'}
    if not set(request) <= allowed:
        raise Refusal('unsupported', 'unsupported synthesize request fields')
    goal = request.get('goal', 'rescue')
    if goal not in ('rescue', 'reprogram'):
        raise Refusal('bad_goal', 'choose rescue or reprogram')
    locus = request.get('locus', 'chr18:56475717-56475888' if goal == 'rescue' else 'chr15:90936956-90937127')
    if not isinstance(locus, str) or not locus:
        raise Refusal('bad_locus', 'locus must be a nonempty string')
    samples = request.get('samples', 12)
    if type(samples) is not int or not 2 <= samples <= 64:
        raise Refusal('bad_samples', 'uncertainty samples must be an integer from 2 to 64')
    run_cx3d = bool(request.get('run_cx3d', False))
    norm = dict(mode='synthesize', goal=goal, locus=locus, samples=samples, run_cx3d=run_cx3d)
    if goal == 'rescue':
        edit = request.get('breaking_edit', '165:G>T')
        if not isinstance(edit, str) or ':' not in edit or '>' not in edit:
            raise Refusal('bad_edit', 'breaking_edit must use POS:REF>ALT syntax (e.g. 165:G>T)')
        norm['breaking_edit'] = edit.strip()
    else:
        tfate = request.get('target_fate', 'Olig2')
        if not isinstance(tfate, str) or not tfate:
            raise Refusal('bad_fate', 'target_fate must be a nonempty string')
        norm['target_fate'] = tfate.strip()
    return norm


def run_synthesize(request: dict) -> dict:
    if synth_mod is None:
        raise Refusal('pipe_unavailable', 'inverse compiler backend not found')
    try:
        seq = synth_mod._load_gse_locus(request['locus'])
    except ValueError as err:
        raise Refusal('bad_locus', str(err)) from err
    rec_id = request['locus'].replace(':', '_').replace('-', '_')
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / 'syn'
        try:
            if request['goal'] == 'rescue':
                rep = synth_mod.synthesize_epistatic_rescue(
                    seq,
                    request['breaking_edit'],
                    out,
                    record_id=rec_id,
                    ensemble_samples=request['samples'],
                    run_cx3d=request['run_cx3d'],
                    cx3d_duration=1.2,
                    export_bundles=False,
                )
                return dict(
                    mode='synthesize',
                    goal='rescue',
                    locus=request['locus'],
                    biological_acceptance=rep['biological_acceptance'],
                    wt_fate=rep['wt_fate'],
                    primary_breaking_edit=rep['primary_breaking_edit'],
                    synthesized_epistatic_rescue=rep['synthesized_epistatic_rescue'],
                    provenance=dict(
                        broken=_strip_host_fields(rep['provenance']['broken']),
                        epistatic_rescue=_strip_host_fields(rep['provenance']['epistatic_rescue']),
                    ),
                    controls=rep['controls'],
                    cx3d_3d_restoration=rep['cx3d_3d_restoration'],
                )
            rep = synth_mod.synthesize_reprogram(
                seq,
                request['target_fate'],
                out,
                record_id=rec_id,
                ensemble_samples=request['samples'],
                run_cx3d=request['run_cx3d'],
                cx3d_duration=1.2,
                export_bundles=False,
            )
            return dict(
                mode='synthesize',
                goal='reprogram',
                locus=request['locus'],
                biological_acceptance=rep['biological_acceptance'],
                wt_fate=rep['wt_fate'],
                target_fate=rep['target_fate'],
                candidates_scanned=rep['candidates_scanned'],
                single_base_hits=rep['single_base_hits'],
                synthesized_edit=rep['synthesized_edit'],
                provenance=_strip_host_fields(rep['provenance']),
                controls=rep['controls'],
                morphogenesis=_strip_cx3d_paths(rep['morphogenesis']),
            )
        except ValueError as err:
            raise Refusal('synth_failed', str(err)) from err


def normalize_discover(request: dict) -> dict:
    if discover_mod is None:
        raise Refusal('pipe_unavailable', 'biophysical discovery backend not found')
    if not set(request) <= {'mode', 'bootstraps'}:
        raise Refusal('unsupported', 'unsupported discover request fields')
    boot = request.get('bootstraps', 200)
    if type(boot) is not int or not 20 <= boot <= 1000:
        raise Refusal('bad_bootstraps', 'bootstraps must be an integer from 20 to 1000')
    return dict(mode='discover', bootstraps=boot)


def run_discover(request: dict) -> dict:
    if discover_mod is None:
        raise Refusal('pipe_unavailable', 'biophysical discovery backend not found')
    with tempfile.TemporaryDirectory() as tmp:
        rep = discover_mod.discover_gse188264(
            Path(tmp) / 'disc',
            bootstraps=request['bootstraps'],
            run_cx3d=False,
            export_bundles=False,
        )
    top = json.loads(json.dumps(rep['top_corroborated_unseen_edit']))
    top['pipeline_provenance'] = _strip_host_fields(top['pipeline_provenance'])
    return dict(
        mode='discover',
        dataset=rep['dataset'],
        biological_acceptance=rep['biological_acceptance'],
        cohort=rep['cohort'],
        validation_selection=rep['validation_selection'],
        heldout_test_metrics=rep['heldout_test_metrics'],
        discovered_top_tfs=rep['discovered_top_tfs'][:10],
        top_corroborated_unseen_edit=top,
    )

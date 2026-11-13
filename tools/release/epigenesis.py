#!/usr/bin/env python3
"""One entry point for the Epigenesis demos. Run `python3 epigenesis.py` for the list."""
import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
EX = ROOT / 'examples'
OUT = Path.cwd() / 'epigenesis-out'


def run(*cmd):
    return subprocess.call([sys.executable, *map(str, cmd)])


def report(script, name, check=False):
    out = OUT / name
    if out.exists():
        sys.exit(f'{out} already exists; delete it or run from another folder.')
    OUT.mkdir(exist_ok=True)
    if run(script, 'demo', '--output', out) or (check and run(script, 'verify', '--output', out)):
        sys.exit('demo failed')
    print(f'\nOpen {out / "index.html"} in your browser.')


def predict(args):
    legacy = args.activity or args.attention
    if args.budget is not None and legacy:
        p.error('--budget belongs to the default Huber-linear and Kreimer models')
    if args.families and not legacy:
        p.error('the default models use family fields in the input and --budget; --families requires --activity or --attention')
    if args.device is not None and not args.attention:
        p.error('--device requires --attention')
    if args.kreimer:
        script = ROOT / 'studio/epigenesis_studio/kreimer.py'
    elif legacy:
        script = EX / 'neural-predictor' / ('predictor.py' if args.attention else 'activity.py')
    else:
        script = ROOT / 'studio/epigenesis_studio/model/predict.py'
    cmd = [script, '--input', args.input, '--output', args.output]
    if args.budget is not None:
        cmd += ['--budget', args.budget]
    if args.families:
        cmd += ['--families', args.families]
    if args.attention:
        cmd += ['--device', args.device or 'cpu']
    return run(*cmd)


def _circuit_opts(args):
    extra = []
    if getattr(args, 'circuit', None) is not None:
        extra += ['--circuit', args.circuit]
    if getattr(args, 'protein', None) is not None:
        extra += ['--protein', args.protein]
    if getattr(args, 'target', None) is not None:
        extra += ['--target', args.target]
    if getattr(args, 'signal', None) is not None:
        extra += ['--signal', args.signal]
    return extra


def scan(args):
    out = args.output or (OUT / 'enhancer-scan')
    if out.exists() and args.output is not None:
        sys.exit(f'{out} already exists; delete it or choose another --output.')
    while out.exists():
        out = out.parent / f'{out.name}-next'
    OUT.mkdir(exist_ok=True)
    cmd = [EX / 'neural-circuit' / 'circuit.py', 'scan', '--fasta', args.fasta, *_circuit_opts(args), '--output', out]
    rc = run(*cmd)
    if rc == 0:
        print(f'\nScan output: {out / "scan.json"}')
    return rc


def radar(args):
    script = EX / 'neural-v2f' / 'radar.py'
    modes = sum([bool(args.fasta is not None), bool(args.variants is not None), bool(args.benchmark)])
    if modes > 1:
        p.error('radar accepts at most one of --fasta, --variants, or --benchmark')
    out = args.output or (OUT / ('v2f-benchmark' if args.benchmark else 'v2f-radar'))
    if out.exists() and args.output is not None:
        sys.exit(f'{out} already exists; delete it or choose another --output.')
    while out.exists():
        out = out.parent / f'{out.name}-next'
    OUT.mkdir(exist_ok=True)
    opts = []
    if args.circuit is not None:
        opts += ['--circuit', args.circuit]
    if args.protein is not None:
        opts += ['--protein', args.protein]
    if args.target is not None:
        opts += ['--target', args.target]
    if args.signal_min is not None:
        opts += ['--signal-min', args.signal_min]
    if args.signal_max is not None:
        opts += ['--signal-max', args.signal_max]
    if args.steps is not None:
        opts += ['--steps', args.steps]
    if args.no_sequence_models:
        opts.append('--no-sequence-models')
    if args.benchmark:
        if args.shuffles is not None:
            opts += ['--shuffles', args.shuffles]
        return run(script, 'benchmark', *opts, '--output', out)
    if args.variants is not None:
        return run(script, 'variants', '--variants', args.variants, *opts, '--output', out)
    if args.fasta is not None:
        if args.shuffles is not None:
            opts += ['--shuffles', args.shuffles]
        return run(script, 'scan', '--fasta', args.fasta, *opts, '--output', out)
    if args.shuffles is not None:
        opts += ['--shuffles', args.shuffles]
    rc = run(script, 'demo', '--locus', args.locus, *opts, '--output', out)
    if rc == 0:
        print(f'\nOpen {out / "index.html"} in your browser.')
    return rc


def circuit(args):
    script = EX / 'neural-circuit' / 'circuit.py'
    if args.list_circuits:
        return run(script, 'list-circuits')
    if args.list_motifs:
        return run(script, 'list-motifs')
    out = args.output or (OUT / 'neural-circuit')
    if out.exists() and args.output is not None:
        sys.exit(f'{out} already exists; delete it or choose another --output.')
    while out.exists():
        out = out.parent / f'{out.name}-next'
    OUT.mkdir(exist_ok=True)
    opts = _circuit_opts(args)
    if args.benchmark:
        return run(script, 'benchmark', *opts, '--output', out)
    if bool(args.reference) != bool(args.variant):
        p.error('circuit requires both --reference and --variant')
    if args.reference and args.variant:
        subcmd = 'bridge' if args.with_construction else 'compare'
        cmd = [script, subcmd, '--reference', args.reference, '--variant', args.variant, *opts, '--output', out]
        if args.with_construction and args.mapping is not None:
            cmd += ['--mapping', args.mapping]
        return run(*cmd)
    cmd = [script, 'demo', '--locus', args.locus, *opts, '--output', out]
    if args.with_construction:
        cmd.append('--with-construction')
    rc = run(*cmd)
    if rc == 0:
        print(f'\nOpen {out / "index.html"} in your browser.')
    return rc


def verify(_):
    return run('-m', 'unittest', 'discover', '-s', ROOT / 'tests') or run(ROOT / 'tests' / 'verify.py', 'minimal-example')


p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
sub = p.add_subparsers(dest='cmd')
sub.add_parser('evidence', help='open a report of 109 published CHE-1 measurements traced to DNA edits')
sub.add_parser('patterning', help='explore an Nkx2.2 enhancer edit in a developmental model')
cc = sub.add_parser('circuit', help='scan an enhancer through Boltzmann TF affinity and an N-gene regulatory ODE')
cc.add_argument('--circuit', help='circuit preset (neural-tube, mutual-toggle, tri-lineage-switch) or JSON spec path')
cc.add_argument('--locus', default='chr15:90936956-90937127')
cc.add_argument('--target')
cc.add_argument('--protein')
cc.add_argument('--signal', type=float)
cc.add_argument('--reference')
cc.add_argument('--variant')
cc.add_argument('--mapping', type=Path)
cc.add_argument('--with-construction', action='store_true')
cc.add_argument('--benchmark', action='store_true')
cc.add_argument('--list-circuits', action='store_true')
cc.add_argument('--list-motifs', action='store_true')
cc.add_argument('--output', type=Path)
sc = sub.add_parser('scan', help='scan a FASTA enhancer across single-base substitutions')
sc.add_argument('--fasta', type=Path, required=True)
sc.add_argument('--circuit')
sc.add_argument('--protein')
sc.add_argument('--target')
sc.add_argument('--signal', type=float)
sc.add_argument('--output', type=Path)
rd = sub.add_parser('radar', help='run V2F tipping-point, dual rescue, and null-control radar')
rd.add_argument('--locus', default='chr15:90936956-90937127')
rd.add_argument('--fasta', type=Path)
rd.add_argument('--variants', type=Path)
rd.add_argument('--benchmark', action='store_true')
rd.add_argument('--circuit')
rd.add_argument('--protein')
rd.add_argument('--target')
rd.add_argument('--signal-min', type=float)
rd.add_argument('--signal-max', type=float)
rd.add_argument('--steps', type=int)
rd.add_argument('--shuffles', type=int)
rd.add_argument('--no-sequence-models', action='store_true')
rd.add_argument('--output', type=Path)
sub.add_parser('construction', help='compile a 480-bp developmental recipe and run the CX3D + LIF assay')
q = sub.add_parser('predict', help='score reporter DNA edits with Huber-linear (default) or Kreimer 0-72h (needs numpy)')
q.add_argument('--input', required=True, help='JSON list of {id, reference, sequence, optional family}; uppercase ACGT, at least 5 bases per sequence (see studio/epigenesis_studio/example-edits.json)')
q.add_argument('--output', required=True, help='new JSON file to write')
q.add_argument('--budget', type=int, help='panel size for Huber-linear or Kreimer; requires a family field on every edit')
legacy = q.add_mutually_exclusive_group()
legacy.add_argument('--kreimer', action='store_true', help='use the bundled Kreimer H1 neural-induction 0-72h motif-temporal model')
legacy.add_argument('--activity', action='store_true', help='use the archived Salomon activity-transfer model (270-base inputs; needs torch)')
legacy.add_argument('--attention', action='store_true', help='use the archived Salomon attention mixture (270-base inputs; needs torch)')
q.add_argument('--families', help='legacy models only: JSON {edit id: family id}; selects ten experiments')
q.add_argument('--device', choices=('cpu', 'mps'), help='attention mixture only; defaults to cpu')
sub.add_parser('verify', help='run the release self-checks')
args = p.parse_args()
if sys.version_info < (3, 11):
    sys.exit('Python 3.11 or newer is required.')
if args.cmd == 'evidence':
    report(EX / 'neural-evidence' / 'evidence.py', 'neural-evidence', check=True)
elif args.cmd == 'patterning':
    report(EX / 'neural-patterning' / 'patterning.py', 'neural-patterning')
elif args.cmd == 'circuit':
    sys.exit(circuit(args))
elif args.cmd == 'scan':
    sys.exit(scan(args))
elif args.cmd == 'radar':
    sys.exit(radar(args))
elif args.cmd == 'construction':
    out = OUT / 'neural-construction'
    if out.exists():
        sys.exit(f'{out} already exists; delete it or run from another folder.')
    OUT.mkdir(exist_ok=True)
    if run(EX / 'neural-construction' / 'construction.py', '--output', out):
        sys.exit('demo failed')
    print(f'\nOpen {out / "viewer.html"} in your browser.')
elif args.cmd == 'predict':
    sys.exit(predict(args))
elif args.cmd == 'verify':
    sys.exit(verify(args))
else:
    p.print_help()

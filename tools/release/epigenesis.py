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


def scan(args):
    out = args.output or (OUT / 'enhancer-scan')
    if out.exists() and args.output is not None:
        sys.exit(f'{out} already exists; delete it or choose another --output.')
    while out.exists():
        out = out.parent / f'{out.name}-next'
    OUT.mkdir(exist_ok=True)
    cmd = [EX / 'neural-circuit' / 'circuit.py', 'scan', '--fasta', args.fasta,
           '--protein', args.protein, '--target', args.target, '--signal', args.signal, '--output', out]
    rc = run(*cmd)
    if rc == 0:
        print(f'\nScan output: {out / "scan.json"}')
    return rc


def circuit(args):
    out = args.output or (OUT / 'neural-circuit')
    if out.exists() and args.output is not None:
        sys.exit(f'{out} already exists; delete it or choose another --output.')
    while out.exists():
        out = out.parent / f'{out.name}-next'
    OUT.mkdir(exist_ok=True)
    script = EX / 'neural-circuit' / 'circuit.py'
    if args.benchmark:
        return run(script, 'benchmark', '--protein', args.protein, '--target', args.target,
                   '--signal', args.signal, '--output', out)
    if bool(args.reference) != bool(args.variant):
        p.error('circuit requires both --reference and --variant')
    if args.reference and args.variant:
        subcmd = 'bridge' if args.with_construction else 'compare'
        return run(script, subcmd, '--reference', args.reference, '--variant', args.variant,
                   '--protein', args.protein, '--target', args.target, '--signal', args.signal, '--output', out)
    cmd = [script, 'demo', '--locus', args.locus, '--target', args.target,
           '--signal', args.signal, '--protein', args.protein, '--output', out]
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
cc = sub.add_parser('circuit', help='scan a human enhancer through Boltzmann GLI affinity and the neural-tube ODE')
cc.add_argument('--locus', default='chr15:90936956-90937127')
cc.add_argument('--target', choices=('Nkx2.2', 'Olig2', 'Both'), default='Nkx2.2')
cc.add_argument('--protein', choices=('GLI2', 'GLI3'), default='GLI2')
cc.add_argument('--signal', type=float, default=0.20)
cc.add_argument('--reference')
cc.add_argument('--variant')
cc.add_argument('--with-construction', action='store_true')
cc.add_argument('--benchmark', action='store_true')
cc.add_argument('--output', type=Path)
sc = sub.add_parser('scan', help='scan a FASTA enhancer across single-base GLI substitutions')
sc.add_argument('--fasta', type=Path, required=True)
sc.add_argument('--protein', choices=('GLI2', 'GLI3'), default='GLI2')
sc.add_argument('--target', choices=('Nkx2.2', 'Olig2', 'Both'), default='Nkx2.2')
sc.add_argument('--signal', type=float, default=0.20)
sc.add_argument('--output', type=Path)
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

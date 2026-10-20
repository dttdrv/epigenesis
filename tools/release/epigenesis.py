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
        p.error('--budget belongs to the default Huber-linear model')
    if args.families and not legacy:
        p.error('the default model uses family fields in the input and --budget; --families requires --activity or --attention')
    if args.device is not None and not args.attention:
        p.error('--device requires --attention')
    script = ROOT / 'studio/epigenesis_studio/model/predict.py'
    if legacy:
        script = EX / 'neural-predictor' / ('predictor.py' if args.attention else 'activity.py')
    cmd = [script, '--input', args.input, '--output', args.output]
    if args.budget is not None:
        cmd += ['--budget', args.budget]
    if args.families:
        cmd += ['--families', args.families]
    if args.attention:
        cmd += ['--device', args.device or 'cpu']
    return run(*cmd)


def verify(_):
    return run('-m', 'unittest', 'discover', '-s', ROOT / 'tests') or run(ROOT / 'tests' / 'verify.py', 'minimal-example')


p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
sub = p.add_subparsers(dest='cmd')
sub.add_parser('evidence', help='open a report of 109 published CHE-1 measurements traced to DNA edits')
sub.add_parser('patterning', help='explore an Nkx2.2 enhancer edit in a developmental model')
q = sub.add_parser('predict', help='score reporter DNA edits with the delivered Huber-linear model (needs numpy)')
q.add_argument('--input', required=True, help='JSON list of {id, reference, sequence, optional family}; uppercase ACGT, at least 5 bases per sequence (see studio/epigenesis_studio/example-edits.json)')
q.add_argument('--output', required=True, help='new JSON file to write')
q.add_argument('--budget', type=int, help='panel size for Huber-linear; requires a family field on every edit')
legacy = q.add_mutually_exclusive_group()
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
elif args.cmd == 'predict':
    sys.exit(predict(args))
elif args.cmd == 'verify':
    sys.exit(verify(args))
else:
    p.print_help()

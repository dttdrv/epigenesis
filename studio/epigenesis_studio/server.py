"""Epigenesis Studio: a local app for the frozen Siraj Huber-linear edit model.

`epigenesis-studio [FILE]` opens the app, optionally with a file preloaded.
`epigenesis-studio verify RECEIPT` replays a saved receipt without the app.
The app scans every single-base change in a sequence, scores reference/edited
pairs, selects one edit per family for an experimental panel, and writes
receipts that anyone can replay. The server listens on the loopback interface only.
"""

import argparse
import json
import re
import sys
import threading
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE/'model'))

import predict  # frozen inference code; model/model.json pins its digest
from brainc import _canonical as canonical  # the compiler's RFC 8785 digests

from . import simulate
from .refusal import Refusal

MANIFEST_BYTES = (HERE/'model/model.json').read_bytes()
MANIFEST = json.loads(MANIFEST_BYTES)
MANIFEST_SHA256 = canonical.bytes_digest(MANIFEST_BYTES)
RECEIPT_FORMAT = 'epigenesis-studio-receipt-v1'
MAX_ROWS = 5000
MAX_BASES = 2000
PRELOAD = None  # (name, text, kind) of the file named on the command line
FOLDER = Path.cwd()  # where the app looks for files to open
MAX_SCAN = 1000
TOLERANCE = 1e-9
# a single-base change only alters windows that overlap it, so a scan scores each
# change inside a crop of this radius; every other window cancels in the difference.
RADIUS = max(MANIFEST['k'], *MANIFEST['motif_widths']) - 1


def bases(text, field, id=None):
    label = field if id is None else f'edit {id!r} {field}'
    where = dict(field=field, id=id)
    if not isinstance(text, str):
        raise Refusal('not_text', f'{label} must be text', **where)
    kept = [c for c in text if not c.isspace() and not c.isdigit()]
    sequence = ''.join(kept).upper()
    bad = sorted(set(sequence)-set('ACGT'))
    if bad:
        raise Refusal('bad_chars', f'{label} contains characters other than A, C, G and T: {" ".join(bad)}',
                      chars=''.join(bad), **where)
    if not MANIFEST['k'] <= len(sequence) <= MAX_BASES:
        raise Refusal('bad_length', f'{label} needs {MANIFEST["k"]} to {MAX_BASES} bases',
                      min=MANIFEST['k'], max=MAX_BASES, **where)
    return sequence


def fasta_sequence(text):
    if not isinstance(text, str):
        raise Refusal('not_text', 'sequence must be text', field='sequence', id=None)
    lines = text.splitlines()
    records = sum(l.startswith('>') for l in lines)
    if records > 1:
        raise Refusal('multi_record', f'the FASTA has {records} records; scan one record at a time', n=records)
    sequence = bases('\n'.join(l for l in lines if not l.startswith('>')), 'sequence')
    if len(sequence) > MAX_SCAN:
        raise Refusal('scan_max', f'a scan takes at most {MAX_SCAN} bases', max=MAX_SCAN)
    return sequence


def scan_rows(reference):
    rows = []
    for i, base in enumerate(reference):
        lo, hi = max(0, i-RADIUS), i+RADIUS+1
        window = reference[lo:hi]
        for alt in 'ACGT':
            if alt != base:
                rows.append(dict(id=f'{i+1}:{base}>{alt}', reference=window,
                                 sequence=window[:i-lo]+alt+window[i-lo+1:]))
    return rows


def normalize(request):
    mode = request.get('mode')
    if mode == 'simulate':
        return simulate.normalize(request)
    if mode == 'scan' and set(request) == {'mode', 'sequence'}:
        return dict(mode='scan', sequence=fasta_sequence(request['sequence']))
    if mode == 'pairs' and set(request) <= {'mode', 'rows', 'budget'} and 'rows' in request:
        rows = request['rows']
        if not isinstance(rows, list) or not 0 < len(rows) <= MAX_ROWS:
            raise Refusal('row_count', f'provide 1 to {MAX_ROWS} edits', max=MAX_ROWS)
        clean = []
        for row in rows:
            if not isinstance(row, dict) or not isinstance(row.get('id'), str) \
                    or not {'id', 'reference', 'sequence'} <= set(row) <= {'id', 'reference', 'sequence', 'family'}:
                raise Refusal('row_shape', 'each edit needs an id, a reference and a sequence, and may add a family')
            if 'family' in row and not (isinstance(row['family'], str) and row['family']):
                raise Refusal('bad_family', f'edit {row["id"]!r}: family must be a nonempty string', id=row['id'])
            clean.append({key: bases(value, key, row['id']) if key in ('reference', 'sequence') else value
                          for key, value in row.items()})
        ids = [row['id'] for row in clean]
        if len(set(ids)) != len(ids):
            raise Refusal('dup_id', 'edit IDs must be unique', id=next(i for i in ids if ids.count(i) > 1))
        budget = request.get('budget')
        if budget is not None and (type(budget) is not int or budget <= 0):
            raise Refusal('bad_budget', 'panel size must be a positive whole number')
        if budget is not None and (any('family' not in row for row in clean)
                                   or len({row['family'] for row in clean}) < budget):
            raise Refusal('need_families', 'a panel needs a family on every edit and at least as many families as its size',
                          budget=budget)
        return dict(mode='pairs', rows=clean, budget=budget)
    raise Refusal('unsupported', 'unsupported request')


def run(request):
    if request['mode'] == 'simulate':
        return simulate.run(request)
    if request['mode'] == 'scan':
        reference = request['sequence']
        return dict(mode='scan', length=len(reference), predictions=predict.predict(scan_rows(reference)))
    predictions = predict.predict(request['rows'])
    result = dict(mode='pairs', predictions=predictions)
    if request['budget'] is not None:
        result['selection'] = predict.select(predictions, request['rows'], request['budget'])
    return result


def identity(mode):
    if mode == 'simulate':
        return simulate.IDENTITY
    return dict(model=MANIFEST['format'], manifest_sha256=MANIFEST_SHA256, endpoint=MANIFEST['endpoint'])


def score(raw_request):
    request = normalize(raw_request)
    receipt = dict(format=RECEIPT_FORMAT, **identity(request['mode']), request=request, result=run(request))
    receipt['receipt_sha256'] = canonical.digest(receipt)
    return receipt


def agree(stored, again, accord):
    """Deep equality of a stored result and its replay; numbers match within TOLERANCE."""
    if isinstance(stored, dict) and isinstance(again, dict):
        same = stored.keys() == again.keys()
        accord['shape'] &= same
        return all([agree(stored[k], again[k], accord) for k in stored.keys() & again.keys()]) and same
    if isinstance(stored, list) and isinstance(again, list):
        same = len(stored) == len(again)
        accord['shape'] &= same
        return all([agree(a, b, accord) for a, b in zip(stored, again)]) and same
    if type(stored) in (int, float) and type(again) in (int, float):
        accord['gap'] = max(accord['gap'], abs(stored-again))
        return abs(stored-again) <= TOLERANCE
    accord['shape'] &= type(stored) is type(again)
    return stored == again


def verify(receipt):
    try:
        if set(receipt) != {'format', 'model', 'manifest_sha256', 'endpoint', 'request', 'result', 'receipt_sha256'} \
                or receipt['format'] != RECEIPT_FORMAT:
            raise Refusal('not_receipt', 'not an Epigenesis Studio receipt')
        request = normalize(receipt['request'])
        expected = identity(request['mode'])
        content = {key: value for key, value in receipt.items() if key != 'receipt_sha256'}
        accord = dict(gap=0.0, shape=True)
        replays = agree(receipt['result'], run(request), accord)
        checks = [('receipt digest matches its content', canonical.digest(content) == receipt['receipt_sha256']),
                  ('model identity matches this installation', all(receipt[key] == value for key, value in expected.items())),
                  (f'results replay within {TOLERANCE:g}', replays)]
    except (KeyError, TypeError, AttributeError) as failure:
        raise Refusal('malformed_receipt', 'malformed receipt') from failure
    return dict(ok=all(ok for _, ok in checks), max_abs_difference=accord['gap'] if accord['shape'] else None,
                tolerance=TOLERANCE, checks=[dict(check=name, ok=ok) for name, ok in checks])


SEQUENCE_SUFFIXES = ('.fa', '.fasta', '.fna', '.fas', '.seq', '.txt')
EDIT_SUFFIXES = ('.json', '.csv', '.tsv')
MAX_FILES = 200


def classify(name, text):
    """What the app can do with a file: 'receipt', 'edits', 'sequence', or None. Only the head is read."""
    head, suffix = text[:2048], Path(name).suffix.lower()
    if suffix == '.json' and RECEIPT_FORMAT in head:
        return 'receipt'
    if suffix in EDIT_SUFFIXES and 'reference' in head.lower() and 'sequence' in head.lower():
        return 'edits'
    if suffix in SEQUENCE_SUFFIXES and (head.lstrip().startswith('>') or re.fullmatch(r'\s*[ACGTacgt][ACGTacgt\s]*', head)):
        return 'sequence'


def kind(request):
    name, text = (request.get(key) for key in ('name', 'text')) if type(request) is dict else (None, None)
    found = type(name) is str and type(text) is str and classify(name, text)
    if not found:
        raise Refusal('unknown_file', 'not a sequence, an edit table or a receipt')
    return dict(kind=found)


def folder_files():
    found = []
    for path in sorted(FOLDER.iterdir(), key=lambda p: p.name.lower()):
        if path.name.startswith('.') or path.suffix.lower() not in SEQUENCE_SUFFIXES+EDIT_SUFFIXES:
            continue
        try:
            size = path.stat().st_size
            if path.is_file() and 0 < size <= canonical.MAX_JSON_BYTES:
                with path.open('rb') as handle:
                    found.append(dict(name=path.name, size=size, kind=classify(path.name, handle.read(2048).decode('utf-8', 'ignore'))))
        except OSError:
            continue
    return [f for f in found if f['kind']][:MAX_FILES]


def open_file(name):
    """Only a file the listing offers can be opened, so a name can never reach outside the folder."""
    listed = next((f for f in folder_files() if f['name'] == name), None)
    if listed is None:
        raise Refusal('no_file', 'no such file in this folder', name=name)
    try:
        return dict(listed, text=(FOLDER/name).read_text('utf-8'))
    except UnicodeDecodeError as failure:
        raise Refusal('not_text', 'the file is not text') from failure


def model_info():
    return dict(format=MANIFEST['format'], endpoint=MANIFEST['endpoint'], scope=MANIFEST['biological_scope'],
                status=MANIFEST['status'], manifest_sha256=MANIFEST_SHA256, simulation=simulate.available(), motifs=len(MANIFEST['motif_ids']),
                motif_source=MANIFEST['motif_source']['license'], radius=RADIUS, max_scan=MAX_SCAN,
                max_rows=MAX_ROWS, max_bases=MAX_BASES)


class Handler(BaseHTTPRequestHandler):
    server_version = 'EpigenesisStudio'

    def log_message(self, *_):
        pass

    def send(self, status, body, kind='application/json'):
        if not isinstance(body, bytes):
            body = json.dumps(body, allow_nan=False).encode()
        self.send_response(status)
        self.send_header('Content-Type', kind+'; charset=utf-8')
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
        self.send_header('X-Content-Type-Options', 'nosniff')
        self.end_headers()
        self.wfile.write(body)

    def allowed(self):
        port = self.server.server_address[1]
        return self.headers.get('Host') in (f'127.0.0.1:{port}', f'localhost:{port}')

    def refuse(self, failure):
        self.send(400, dict(error=str(failure), code=getattr(failure, 'code', None), params=getattr(failure, 'params', {})))

    def do_GET(self):
        if not self.allowed():
            return self.send(403, dict(error='unexpected host'))
        path, _, query = self.path.partition('?')
        if path in ('/', '/index.html'):
            return self.send(200, (HERE/'index.html').read_bytes(), 'text/html')
        if path == '/api/model':
            return self.send(200, model_info())
        if path == '/api/example':
            return self.send(200, (HERE/'example-edits.json').read_bytes())
        try:
            if path == '/api/site':
                return self.send(200, simulate.site_table())
            if path == '/api/files':
                return self.send(200, dict(folder=str(FOLDER), files=folder_files()))
            if path == '/api/file':
                return self.send(200, open_file(parse_qs(query).get('name', [''])[0]))
        except ValueError as failure:
            return self.refuse(failure)
        if path == '/api/preload':
            return self.send(200, PRELOAD and dict(zip(('name', 'text', 'kind'), PRELOAD)))
        self.send(404, dict(error='not found'))

    def do_POST(self):
        handler = dict(zip(('/api/score', '/api/verify', '/api/kind'), (score, verify, kind))).get(self.path)
        try:
            if not self.allowed():
                return self.send(403, dict(error='unexpected host'))
            if handler is None:
                return self.send(404, dict(error='not found'))
            if self.headers.get('Content-Type', '').split(';')[0] != 'application/json':
                return self.send(415, dict(error='send application/json'))
            size = int(self.headers.get('Content-Length', ''))
            if not 0 < size <= canonical.MAX_JSON_BYTES:
                return self.send(413, dict(error='request body too large'))
            self.send(200, handler(canonical.loads(self.rfile.read(size), 'request')))
        except ValueError as failure:
            self.refuse(failure)
        except Exception as failure:
            print(f'{type(failure).__name__}: {failure}', file=sys.stderr)
            self.send(500, dict(error='internal error'))


def verify_command(path):
    outcome = verify(canonical.loads(Path(path).read_bytes(), 'receipt'))
    for check in outcome['checks']:
        print(('ok    ' if check['ok'] else 'FAIL  ')+check['check'])
    return 0 if outcome['ok'] else 1


def main():
    global PRELOAD, FOLDER
    if sys.argv[1:2] == ['verify']:
        parser = argparse.ArgumentParser(prog='epigenesis-studio verify')
        parser.add_argument('receipt', type=Path)
        try:
            return verify_command(parser.parse_args(sys.argv[2:]).receipt)
        except (OSError, ValueError) as failure:
            return f'epigenesis-studio: {failure}'
    parser = argparse.ArgumentParser(prog='epigenesis-studio', description=__doc__.split('\n\n')[0])
    parser.add_argument('file', nargs='?', type=Path, help='a file to open (edits, a sequence or a receipt), or a folder to look in instead of the current one')
    parser.add_argument('--port', type=int, help='default 8765, or any free port when that one is taken')
    parser.add_argument('--no-browser', action='store_true')
    args = parser.parse_args()
    if args.file and args.file.is_dir():
        FOLDER = args.file.resolve()
    elif args.file:
        try:
            raw = args.file.read_bytes()
            if len(raw) > canonical.MAX_JSON_BYTES:
                raise ValueError('file is too large')
            text = raw.decode('utf-8')
            PRELOAD = (args.file.name, text, classify(args.file.name, text) or 'sequence')
        except (OSError, ValueError) as failure:
            return f'epigenesis-studio: {args.file}: {failure}'
    predict.load()
    for port in [args.port] if args.port is not None else [8765, 0]:
        try:
            server = ThreadingHTTPServer(('127.0.0.1', port), Handler)
            break
        except OSError as failure:
            error = failure
    else:
        return f'epigenesis-studio: {error}'
    url = f'http://127.0.0.1:{server.server_address[1]}/'
    print(f'Epigenesis Studio at {url}  (Ctrl+C to stop)', flush=True)
    print(f'{len(folder_files())} files it can open in {FOLDER}', flush=True)
    if not args.no_browser:
        threading.Timer(0.3, webbrowser.open, (url,)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == '__main__':
    sys.exit(main())

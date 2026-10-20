import csv
import json
import re
import subprocess
import sys
import tempfile
import threading
import unittest
import unittest.mock
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer
from pathlib import Path

try:
    import numpy  # noqa: F401
    from epigenesis_studio import server, simulate
except ImportError as failure:
    raise unittest.SkipTest(f'Epigenesis Studio and its dependencies are required: {failure}')

EXAMPLE = json.loads((Path(server.HERE)/'example-edits.json').read_text())
REFERENCE = EXAMPLE[0]['reference']
# the three-edit panel documented for the frozen model's own run on these 64 edits
PANEL = [('chr10:101119242:G:A', -0.213959), ('chr10:119062236:C:T', -0.203347), ('chr10:102798346:T:C', 0.168665)]


class Engine(unittest.TestCase):
    def test_cropped_scan_equals_scoring_whole_sequences(self):
        scanned = server.score(dict(mode='scan', sequence=REFERENCE))['result']['predictions']
        whole = [dict(id=f'{i+1}:{b}>{a}', reference=REFERENCE, sequence=REFERENCE[:i]+a+REFERENCE[i+1:])
                 for i, b in enumerate(REFERENCE) for a in 'ACGT' if a != b]
        expected = server.predict.predict(whole)
        self.assertEqual([p['id'] for p in scanned], [p['id'] for p in expected])
        self.assertLess(max(abs(a['effect']-b['effect']) for a, b in zip(scanned, expected)), 1e-12)

    def test_scan_scores_every_change_and_strips_fasta(self):
        result = server.score(dict(mode='scan', sequence='>x\n'+REFERENCE[:60]+'\n'+REFERENCE[60:]))['result']
        self.assertEqual(result['length'], len(REFERENCE))
        self.assertEqual(len(result['predictions']), 3*len(REFERENCE))

    def test_panel_matches_the_frozen_model_run(self):
        selection = server.score(dict(mode='pairs', rows=EXAMPLE, budget=3))['result']['selection']
        self.assertEqual([s['id'] for s in selection], [i for i, _ in PANEL])
        for chosen, (_, effect) in zip(selection, PANEL):
            self.assertAlmostEqual(chosen['effect'], effect, places=6)
        self.assertEqual(len({s['family'] for s in selection}), 3)

    def test_receipt_replays_and_each_check_catches_its_own_tampering(self):
        receipt = server.score(dict(mode='pairs', rows=EXAMPLE[:8], budget=2))
        copy = lambda: json.loads(json.dumps(receipt))
        def sealed(edit):
            changed = copy()
            edit(changed)
            changed['receipt_sha256'] = server.canonical.digest({k: v for k, v in changed.items() if k != 'receipt_sha256'})
            return changed
        def verdict(r):
            return [c['ok'] for c in server.verify(r)['checks']]
        self.assertEqual(verdict(receipt), [True, True, True])
        self.assertTrue(server.verify(receipt)['ok'])
        stale = copy()
        stale['receipt_sha256'] = '0'*64
        self.assertEqual(verdict(stale), [False, True, True])
        for field, value in (('manifest_sha256', '0'*64), ('endpoint', 'wrong endpoint'), ('model', 'another model')):
            with self.subTest(field=field):
                self.assertEqual(verdict(sealed(lambda r: r.update({field: value}))), [True, False, True])
        wordless = sealed(lambda r: r['result']['predictions'][0].update(word_effect=r['result']['predictions'][0]['word_effect']+.5))
        self.assertEqual(verdict(wordless), [True, True, False])
        dropped = sealed(lambda r: r['result']['predictions'].pop())
        self.assertEqual(verdict(dropped), [True, True, False])
        rescored = sealed(lambda r: r['result']['predictions'][0].update(effect=r['result']['predictions'][0]['effect']+1e-6))
        self.assertEqual(verdict(rescored), [True, True, False])
        reseeded = sealed(lambda r: r['request']['rows'][0].update(sequence=r['request']['rows'][0]['reference']))
        self.assertEqual(verdict(reseeded), [True, True, False])
        with self.assertRaises(ValueError):
            server.verify({'format': 'something else'})

    def test_a_fasta_with_several_records_is_refused_not_joined(self):
        with self.assertRaises(server.Refusal) as caught:
            server.score(dict(mode='scan', sequence='>one\nACGTA\n>two\nTGCAA'))
        self.assertEqual((caught.exception.code, caught.exception.params), ('multi_record', {'n': 2}))
        single = server.score(dict(mode='scan', sequence='>only\nACGTACGTAC\n'))
        self.assertEqual(single['request']['sequence'], 'ACGTACGTAC')

    def test_bad_requests_are_refused(self):
        row = EXAMPLE[0]
        for request in (dict(mode='scan', sequence='ACGTN'*10), dict(mode='scan', sequence='ACGT'*300),
                        dict(mode='scan', sequence='ACG'), dict(mode='nope'), dict(mode='scan'),
                        dict(mode='pairs', rows=[]), dict(mode='pairs', rows=[row, row]),
                        dict(mode='pairs', rows=[row], budget=2), dict(mode='pairs', rows=[row], budget=0),
                        dict(mode='pairs', rows=[dict(row, extra=1)])):
            with self.subTest(request=str(request)[:60]), self.assertRaises(ValueError):
                server.score(request)


class Files(unittest.TestCase):
    """The app offers what the folder holds and opens only what it offers."""
    def test_files_are_told_apart_by_content_not_just_name(self):
        receipt = json.dumps(server.score(dict(mode='pairs', rows=EXAMPLE[:2])))
        csv_text = 'id,reference,sequence\na,'+REFERENCE+','+REFERENCE+'\n'
        for name, text, kind in [('x.json', receipt, 'receipt'), ('x.json', json.dumps(EXAMPLE[:2]), 'edits'), ('x.csv', csv_text, 'edits'),
                                 ('x.tsv', csv_text.replace(',', '\t'), 'edits'), ('x.fasta', '>a\nACGT', 'sequence'), ('x.txt', ' acgt\nacgt\n', 'sequence'),
                                 ('package.json', '{"name": "app", "version": "1"}', None), ('notes.txt', 'hello world', None),
                                 ('x.csv', 'a,b\n1,2\n', None), ('x.py', '>a\nACGT', None), ('empty.txt', '', None), ('blank.txt', ' \n', None)]:
            with self.subTest(name=name, text=text[:20]):
                self.assertEqual(server.classify(name, text), kind)

    def test_folder_listing_and_opening(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)/'work'
            folder.mkdir()
            (folder/'b.fasta').write_text('>b\n'+REFERENCE)
            (folder/'a-edits.json').write_text(json.dumps(EXAMPLE[:3]))
            (folder/'run.receipt.json').write_text(json.dumps(server.score(dict(mode='scan', sequence=REFERENCE[:40]))))
            (folder/'.hidden.fasta').write_text('>h\nACGTACGT')
            (folder/'notes.txt').write_text('not a sequence')
            (folder/'binary.txt').write_bytes(b'\xff\xfe\x00')
            (folder/'nested').mkdir()
            (folder/'nested'/'deep.fasta').write_text('>d\nACGTACGT')
            (Path(tmp)/'outside.fasta').write_text('>o\nACGTACGT')
            with unittest.mock.patch.object(server, 'FOLDER', folder):
                self.assertEqual([(f['name'], f['kind']) for f in server.folder_files()],
                                 [('a-edits.json', 'edits'), ('b.fasta', 'sequence'), ('run.receipt.json', 'receipt')])
                self.assertEqual(server.open_file('b.fasta')['text'], '>b\n'+REFERENCE)
                for name in ('notes.txt', '.hidden.fasta', 'nested/deep.fasta', '../outside.fasta', str(Path(tmp)/'outside.fasta'), 'missing.fasta', ''):
                    with self.subTest(name=name), self.assertRaises(server.Refusal) as refused:
                        server.open_file(name)
                    self.assertEqual(refused.exception.code, 'no_file')
                (folder/'b.fasta').write_bytes(b'>b\n\xff\xfeACGT')
                with self.assertRaises(server.Refusal) as refused:
                    server.open_file('b.fasta')
                self.assertEqual(refused.exception.code, 'not_text')

    def test_kind_endpoint_refuses_what_the_app_cannot_open(self):
        self.assertEqual(server.kind(dict(name='r.json', text='{"format": "%s"}' % server.RECEIPT_FORMAT)), dict(kind='receipt'))
        for request in (dict(name='n.txt', text='hello'), dict(name='n.txt'), dict(name=1, text='ACGT'), [], None):
            with self.subTest(request=request), self.assertRaises(server.Refusal) as refused:
                server.kind(request)
            self.assertEqual(refused.exception.code, 'unknown_file')


class Http(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.httpd = ThreadingHTTPServer(('127.0.0.1', 0), server.Handler)
        cls.base = f'http://127.0.0.1:{cls.httpd.server_address[1]}'
        threading.Thread(target=cls.httpd.serve_forever, daemon=True).start()

    @classmethod
    def tearDownClass(cls):
        cls.httpd.shutdown()
        cls.httpd.server_close()

    def call(self, path, body=None, headers=None, raw=None):
        data = raw if raw is not None else None if body is None else json.dumps(body).encode()
        request = urllib.request.Request(self.base+path, data=data, headers=headers or {})
        if data is not None:
            request.add_header('Content-Type', (headers or {}).get('Content-Type', 'application/json'))
        try:
            with urllib.request.urlopen(request) as response:
                return response.status, response.read()
        except urllib.error.HTTPError as failure:
            return failure.code, failure.read()

    def test_page_model_and_score_round_trip(self):
        status, page = self.call('/')
        self.assertEqual(status, 200)
        self.assertIn(b'Epigenesis Studio', page)
        self.assertEqual(json.loads(self.call('/api/model')[1])['manifest_sha256'], server.MANIFEST_SHA256)
        status, body = self.call('/api/score', dict(mode='pairs', rows=EXAMPLE[:4]))
        self.assertEqual(status, 200)
        status, body = self.call('/api/verify', json.loads(body))
        self.assertTrue(status == 200 and json.loads(body)['ok'])

    def test_guards(self):
        self.assertEqual(self.call('/api/model', headers={'Host': 'evil.example'})[0], 403)
        self.assertEqual(self.call('/api/score', headers={'Host': 'evil.example', 'Content-Type': 'application/json'}, raw=b'{}')[0], 403)
        self.assertEqual(self.call('/api/score', raw=b'{}', headers={'Content-Type': 'text/plain'})[0], 415)
        self.assertEqual(self.call('/api/score', raw=b'not json')[0], 400)
        self.assertEqual(self.call('/api/score', raw=b'{"mode":"scan","mode":"pairs"}')[0], 400)
        status, body = self.call('/api/score', dict(mode='scan', sequence='XYZ'))
        self.assertEqual((status, json.loads(body)['code']), (400, 'bad_chars'))
        self.assertEqual(self.call('/nowhere')[0], 404)

    def test_json_that_is_not_an_object_is_refused_not_a_server_error(self):
        for body in (b'null', b'[]', b'true', b'1', b'"scan"'):
            for path in ('/api/score', '/api/verify'):
                with self.subTest(path=path, body=body):
                    status, reply = self.call(path, raw=body)
                    self.assertEqual((status, json.loads(reply)['error']), (400, 'invalid request JSON: request must be a JSON object'))
        receipt = server.score(dict(mode='pairs', rows=EXAMPLE[:2]))
        for request in (None, [], 'scan', 1):
            with self.subTest(request=request):
                status, reply = self.call('/api/verify', dict(receipt, request=request))
                self.assertEqual((status, json.loads(reply)['code']), (400, 'malformed_receipt'))

    @unittest.skipUnless(simulate.available(), 'the repository checkout with examples/neural-patterning is required')
    def test_a_protein_that_is_not_a_name_is_refused(self):
        for protein in ([], {}, ['GLI2'], 1, None):
            with self.subTest(protein=protein):
                status, reply = self.call('/api/score', dict(mode='simulate', site='GAACACCCA', protein=protein, phases=[[20, 0.2]]))
                self.assertEqual((status, json.loads(reply)['code']), (400, 'sim_protein'))

    def test_folder_and_import_endpoints(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp)/'one.fasta').write_text('>one\nACGTACGTAC')
            with unittest.mock.patch.object(server, 'FOLDER', Path(tmp)):
                status, body = self.call('/api/files')
                self.assertEqual((status, json.loads(body)), (200, dict(folder=tmp, files=[dict(name='one.fasta', size=len('>one\nACGTACGTAC'), kind='sequence')])))
                status, body = self.call('/api/file?name=one.fasta')
                self.assertEqual((status, json.loads(body)['text']), (200, '>one\nACGTACGTAC'))
                status, body = self.call('/api/file?name=..%2Fetc%2Fhosts')
                self.assertEqual((status, json.loads(body)['code']), (400, 'no_file'))
        status, body = self.call('/api/kind', dict(name='x.fa', text='>x\nACGT'))
        self.assertEqual((status, json.loads(body)), (200, dict(kind='sequence')))
        status, body = self.call('/api/kind', dict(name='x.py', text='print(1)'))
        self.assertEqual((status, json.loads(body)['code']), (400, 'unknown_file'))
        self.assertEqual(self.call('/api/files', headers={'Host': 'evil.example'})[0], 403)


@unittest.skipUnless(simulate.available(), 'the repository checkout with examples/neural-patterning is required')
class Simulation(unittest.TestCase):
    request = dict(mode='simulate', site='GAACACCCA', protein='GLI2', phases=[[20, .2]])
    # Hallikas et al. (2006) Table S1 via examples/neural-patterning/README.md: GLI2 mutant/reference ratio
    PUBLISHED_RATIO = 0.04235549749669733

    @classmethod
    def setUpClass(cls):
        cls.receipt = server.score(cls.request)
        cls.result = cls.receipt['result']

    def test_the_page_translates_every_backend_assumption(self):
        self.assertEqual(len(self.result['assumptions']), 6)

    def test_published_substitution_uses_the_published_ratio_and_reverts_exactly(self):
        self.assertEqual(self.result['site']['ratio'], self.PUBLISHED_RATIO)
        self.assertEqual(self.result['site']['position'], 3)
        controls = self.result['controls']
        self.assertEqual(controls['reversion_max_difference'], 0.0)
        self.assertTrue(controls['reversion_same_module'])
        self.assertEqual(controls['reporter_max_difference'], 0.0)
        compiled = self.result['compiled']
        self.assertNotEqual(compiled['variant']['development_sha256'], compiled['reference']['development_sha256'])
        self.assertEqual(compiled['reversion'], compiled['reference'])

    def test_variant_changes_the_dynamics_and_reference_forms_the_nkx22_domain(self):
        genes = self.result['genes']
        end = lambda run: dict(zip(genes, self.result['runs'][run]['states'][-1]))
        self.assertGreater(end('reference')['Nkx2.2'], .5)
        self.assertLess(end('variant')['Nkx2.2'], .5)

    def test_endogenous_run_matches_the_backends_own_command_line(self):
        backend = Path(simulate.patterning.__file__)
        with tempfile.TemporaryDirectory() as scratch:
            protocol = Path(scratch)/'protocol.json'
            protocol.write_text(json.dumps(dict(phases=[[5, 1], [15, .05]], initial_state=dict(zip(self.result['genes'], [.1, 0, 0, .1])))))
            variant_fasta = backend.parent/'variant.fasta'
            command = [sys.executable, str(backend), 'run', '--fasta', str(variant_fasta), '--record', 'nkx2_2_site', '--start', '0',
                       '--strand', '-1', '--protein', 'GLI2', '--mode', 'endogenous', '--protocol', str(protocol), '--output', str(Path(scratch)/'out')]
            subprocess.run(command, check=True, capture_output=True)
            rows = list(csv.reader((Path(scratch)/'out/trajectory.csv').open()))[1:]
        mine = server.score(dict(self.request, phases=[[5, 1], [15, .05]], initial_state=dict(zip(self.result['genes'], [.1, 0, 0, .1]))))['result']
        self.assertEqual(len(rows), len(mine['time']))
        for row, time, state, occupancy in zip(rows, mine['time'], mine['runs']['variant']['states'], mine['runs']['variant']['occupancy']):
            self.assertEqual([float(v) for v in row[:5]], [time, *state])
            self.assertEqual(float(row[6]), occupancy)

    def test_the_receipt_replays_and_catches_resealed_edits(self):
        copy = lambda: json.loads(json.dumps(self.receipt))
        def sealed(edit):
            changed = copy()
            edit(changed)
            changed['receipt_sha256'] = server.canonical.digest({k: v for k, v in changed.items() if k != 'receipt_sha256'})
            return [c['ok'] for c in server.verify(changed)['checks']]
        self.assertTrue(server.verify(self.receipt)['ok'])
        self.assertEqual(sealed(lambda r: r['result']['runs']['variant']['states'][-1].__setitem__(0, .5)), [True, True, False])
        self.assertEqual(sealed(lambda r: r['result']['compiled']['variant'].update(development_sha256='0'*64)), [True, True, False])
        self.assertEqual(sealed(lambda r: r['request'].update(site='GACCACCAA')), [True, True, False])
        self.assertEqual(sealed(lambda r: r.update(manifest_sha256='0'*64)), [True, False, True])
        self.assertEqual(sealed(lambda r: r.update(endpoint='wrong endpoint')), [True, False, True])

    def test_unsupported_inputs_are_refused_with_codes(self):
        refused = {'sim_same': dict(site='GACCACCCA'), 'sim_site_shape': dict(site='GACC'), 'sim_site_multi': dict(site='GAACACCCC'),
                   'sim_site_zero': dict(site='GACAACCCA'), 'sim_protein': dict(protein='GLI1'),
                   'sim_phases': dict(phases=[[0, .2]]), 'sim_time': dict(phases=[[simulate.MAX_TIME+1, .2]]),
                   'sim_initial': dict(initial_state=dict(Pax6=.1))}
        for code, change in refused.items():
            with self.subTest(code=code), self.assertRaises(server.Refusal) as caught:
                server.score(dict(self.request, **change))
            self.assertEqual(caught.exception.code, code)
        for phases in ([[5, 1.5]], [[5]], [], [[5, .1]]*(simulate.MAX_PHASES+1), [['5', .1]]):
            with self.subTest(phases=phases), self.assertRaises(server.Refusal):
                server.score(dict(self.request, phases=phases))

    def test_every_supported_substitution_in_the_table_runs(self):
        table = simulate.site_table()
        self.assertEqual(table['genome_strand'], 'TGGGTGGTC')
        for protein, rows in table['proteins'].items():
            for i, row in enumerate(rows):
                for alt, ratio in row.items():
                    site = table['consensus'][:i]+alt+table['consensus'][i+1:]
                    if ratio is None:
                        with self.assertRaises(server.Refusal):
                            simulate.normalize(dict(self.request, site=site, protein=protein))
                    else:
                        self.assertTrue(0 < ratio <= 1)
                        simulate.normalize(dict(self.request, site=site, protein=protein))
        self.assertEqual(table['published']['site'], 'GAACACCCA')


class Page(unittest.TestCase):
    """The page has no JavaScript engine here, so its translation tables are checked as text."""
    PAGE = (Path(server.HERE)/'index.html').read_text()

    def table(self, name):
        block = re.search(r'/\*%s\*/(.*?)/\*end\*/' % name, self.PAGE, re.S).group(1)
        return set(re.findall(r"'([\w.]+)':", block))

    def test_english_and_bulgarian_have_the_same_keys(self):
        en, bg = self.table('en'), self.table('bg')
        self.assertEqual(en - bg, set())
        self.assertEqual(bg - en, set())

    def test_every_key_the_page_uses_is_translated(self):
        keys = self.table('en')
        used = set(re.findall(r'data-(?:k|ph|label)="([\w.]+)"', self.PAGE)) | set(re.findall(r"\bt\('([\w.]+)'", self.PAGE))
        self.assertEqual({k for k in used if not k.endswith('.')} - keys, set())
        for code in set(re.findall(r"Refusal\('(\w+)'", (Path(server.HERE)/'server.py').read_text()+(Path(server.HERE)/'simulate.py').read_text())):
            self.assertIn('err.'+code, keys, code)
        for family in ('kind.sequence', 'kind.edits', 'kind.receipt', 'dir.up', 'dir.down', 'dir.none', 'verify.c.0', 'verify.c.1', 'verify.c.2', 'field.reference', 'field.sequence', *('sim.a.%d' % i for i in range(6))):
            self.assertIn(family, keys)


if __name__ == '__main__':
    unittest.main()

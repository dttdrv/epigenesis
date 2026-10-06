"""Independent Figure 6F arithmetic and donor-geometry regression checks."""
import copy
import csv
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
import zipfile

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / 'examples/neural-evidence'
spec = importlib.util.spec_from_file_location('neural_evidence', HERE / 'evidence.py')
evidence = importlib.util.module_from_spec(spec)
spec.loader.exec_module(evidence)


class NeuralEvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = evidence.analyze()

    def test_original_workbook_arithmetic(self):
        # Independent selection from the original Figure 6 notebook, cells 7–11.
        expected = {('core','control'):(22,602), ('core','auxin'):(21,50),
                    ('fragment','control'):(12,316), ('fragment','auxin'):(30,524),
                    ('native','control'):(14,358), ('native','auxin'):(10,30)}
        for group in self.report['groups']:
            count, total = expected[group['construct'],group['condition']]
            self.assertEqual((group['n'],group['sum']), (count,total))
            self.assertEqual(group['mean'], total/count)
        rows = self.report['observations']
        self.assertEqual(len(rows),109)
        self.assertEqual(sum(r['count']==0 for r in rows),6)
        self.assertEqual(len({(r['sheet'],r['cell']) for r in rows}),109)
        self.assertTrue(any(r['batch']=='bacth 1' for r in rows))
        self.assertTrue(any(r['batch'] is None for r in rows))
        native = next(r for r in rows if r['sheet']=='degron_l3l4_con' and r['cell']=='E2')
        self.assertEqual((native['count'],native['specimen']), (20,'worm 2'))

    def test_original_donor_geometry(self):
        core, fragment, native = self.report['constructs']
        self.assertEqual((core['reference_bp'],core['insert_bp'],core['core_relative_to_gene']),(12,12,'opposite'))
        self.assertEqual((fragment['reference_bp'],fragment['insert_bp'],fragment['core_relative_to_gene']),(160,130,'same'))
        self.assertEqual(core['core'], fragment['core'])
        self.assertEqual(core['reference_interval'],{'accession':'NC_003283.11','start':20825221,'end':20825232,'strand':1})
        self.assertEqual(fragment['source_interval'],{'accession':'NC_003279.8','start':6521063,'end':6521192,'strand':-1})
        self.assertEqual(fragment['insert'][59:71], 'GAAGCCACAATT')
        self.assertEqual(native['core'], 'GAAGCCCTTCAA')

    def test_equal_input_constraint_is_descriptive(self):
        comparison = self.report['comparisons']['core:fragment']
        self.assertAlmostEqual(comparison['auxin_mean_difference'],528/35)
        self.assertAlmostEqual(comparison['constraints']['normalized_core']['minimum_mean_rmse'],264/35)
        self.assertTrue(comparison['constraints']['normalized_core']['equal_inputs'])
        self.assertFalse(comparison['constraints']['oriented_core']['equal_inputs'])
        self.assertEqual(comparison['constraints']['oriented_core']['minimum_mean_rmse'],0)
        self.assertAlmostEqual(self.report['constructs'][1]['retention'], (524/30)/(316/12))

    def test_representation_uses_physical_inputs_not_names(self):
        original = self.report['constructs'][0]
        renamed = dict(original, donor='renamed', id='another-name', label='Another label')
        for kind in ('normalized_core','oriented_core','full_design'):
            self.assertEqual(evidence.representation(original,kind), evidence.representation(renamed,kind))
        altered = dict(original, donor_sequence='C'+original['donor_sequence'][1:])
        self.assertNotEqual(evidence.representation(original,'full_design'), evidence.representation(altered,'full_design'))
        with self.assertRaisesRegex(ValueError,'unknown sequence'):
            evidence.representation(original,'unknown')

    def test_imported_implementation_identity_survives_disk_change(self):
        with tempfile.TemporaryDirectory() as td:
            source = Path(td)/'evidence.py'
            source.write_bytes((HERE/'evidence.py').read_bytes())
            local_spec = importlib.util.spec_from_file_location('identity_test',source)
            loaded = importlib.util.module_from_spec(local_spec)
            local_spec.loader.exec_module(loaded)
            original_hash = loaded.IMPLEMENTATION_SHA256
            source.write_bytes(source.read_bytes()+b'\n# changed on disk\n')
            self.assertEqual(loaded.analyze(HERE)['implementation_sha256'],original_hash)
            self.assertNotEqual(evidence.sha(source.read_bytes()),original_hash)

    def test_missing_and_invalid_counts_fail(self):
        case, sources = evidence.load_case()
        sheet = evidence.workbook(sources['figure6f.xlsx'])
        selector = case['groups'][0]
        for invalid in ('text', -1, 1.5, True, float('nan')):
            broken = copy.deepcopy(sheet)
            broken[selector['sheet']]['C2'] = invalid
            with self.subTest(invalid=invalid), self.assertRaisesRegex(ValueError,'molecule count'):
                evidence.observations(broken, case['groups'])
        missing = copy.deepcopy(sheet)
        del missing[selector['sheet']]
        with self.assertRaisesRegex(ValueError,'sheet'):
            evidence.observations(missing, case['groups'])
        broken = copy.deepcopy(case['groups'])
        del broken[0]['batch_column']
        with self.assertRaisesRegex(ValueError,'selection metadata'):
            evidence.observations(sheet, broken)

    def test_source_and_selection_substitution_fail_before_output(self):
        with tempfile.TemporaryDirectory() as td:
            clone = Path(td)/'case'
            import shutil
            shutil.copytree(HERE,clone)
            source = clone/'data/figure6f.xlsx'
            source.write_bytes(source.read_bytes()+b'changed')
            with self.assertRaisesRegex(ValueError,'source identity'):
                evidence.analyze(clone)
            shutil.copyfile(HERE/'data/figure6f.xlsx',source)
            case = json.loads((clone/'case.json').read_text())
            case['groups'][3]['column']='D'
            (clone/'case.json').write_text(json.dumps(case))
            with self.assertRaisesRegex(ValueError,'case identity'):
                evidence.analyze(clone)

    def test_office_formulas_rejected(self):
        raw = (HERE/'data/figure6f.xlsx').read_bytes()
        original = zipfile.ZipFile(io.BytesIO(raw))
        target = next(n for n in original.namelist() if n.startswith('xl/worksheets/sheet'))
        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer,'w') as z:
            for name in original.namelist():
                content=original.read(name)
                if name==target:
                    content=content.replace(b'<sheetData>',b'<sheetData><row r="999"><c r="A999"><f>1+1</f><v>2</v></c></row>')
                z.writestr(name,content)
        with self.assertRaisesRegex(ValueError,'formula'):
            evidence.workbook(buffer.getvalue())

    def test_reconstruction_requires_unique_arms(self):
        case,sources=evidence.load_case()
        reference=evidence.fasta(sources['gcy-22.fasta'])
        donors=evidence.donors(sources['donors.docx'])
        broken=copy.deepcopy(reference)
        broken['sequence'] *= 2
        with self.assertRaisesRegex(ValueError,'uniquely'):
            evidence.reconstruct(donors['3418'],broken,case['sequence'])

    def test_output_is_portable_and_verifiable(self):
        with tempfile.TemporaryDirectory() as td:
            output=Path(td)/'report'
            evidence.write_report(output)
            evidence.verify_report(output)
            rows=list(csv.DictReader(io.StringIO((output/'observations.csv').read_text())))
            self.assertEqual(len(rows),109)
            self.assertEqual(sum(int(r['count']) for r in rows),1880)
            html=(output/'index.html').read_text()
            self.assertNotIn('__EVIDENCE_DATA__',html)
            self.assertIn('application/json',html)
            self.assertTrue((output/'sources/figure6f.xlsx').is_file())
            with self.assertRaises(FileExistsError):
                evidence.write_report(output)
            (output/'observations.csv').write_text('wrong')
            with self.assertRaisesRegex(ValueError,'output differs'):
                evidence.verify_report(output)


if __name__=='__main__':
    unittest.main()

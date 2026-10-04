from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT/"examples/neural-binding/binding.py"


class NeuralBindingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location("neural_binding", EXAMPLE)
        cls.binding = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.binding)

    def test_conditional_mean_normalizes_each_experiment_before_averaging(self):
        # FCpackage matrixSVD + trainSVD, one feature in each of all three PCs.
        rows = [{"protein_alignment":aa, "relative_affinity_by_base":dict(zip("ACGT",([v] for v in values)))}
                for aa,values in [("A",[8,1,.5,.5]),("A",[.2,.2,.2,.4]),("C",[.1,.1,.7,.1])]]
        result = self.binding._predict(rows,1,"A")
        for base,expected in zip("ACGT",[1,.3,.25,.45]):
            self.assertAlmostEqual(result["psam"][base][0],expected,places=14)
        self.assertEqual(result["matched_training_rows"],2)
        self.assertFalse(result["interpolated"])
        self.assertAlmostEqual(sum(result["centered_log_preference"][b][0] for b in "ACGT"),0,places=14)
        self.assertGreater(result["centered_log_preference"]["A"][0],0)
        self.assertNotAlmostEqual(sum(result["psam"][b][0] for b in "ACGT"),1)

    def test_unseen_amino_acid_weights_groups_not_sample_counts(self):
        rows = [{"protein_alignment":aa,"relative_affinity_by_base":dict(zip("ACGT",([v] for v in values)))}
                for aa,values in [("A",[1,0,0,0]),("A",[1,0,0,0]),("C",[0,1,0,0])]]
        # Published BLOSUM62: W/W=11, W/A=-3, W/C=-2; weights 1/14 and 1/13.
        result = self.binding._predict(rows,1,"W")
        expected = {"A":13/14,"C":1,"G":.01,"T":.01}
        for base,value in expected.items():
            self.assertAlmostEqual(result["psam"][base][0],value,places=14)
        self.assertTrue(result["interpolated"])
        self.assertEqual(result["matched_training_rows"],0)

    def test_original_model_rows_axes_and_duplicate_names_are_preserved(self):
        model = self.binding.load_model()
        self.assertEqual(len(model),414)
        self.assertEqual({len(r["protein_alignment"]) for r in model},{57})
        for name in ("CRX","SIX6"):
            rows = [r for r in model if r["name"] == name]
            self.assertEqual(len(rows),2)
            self.assertEqual(len({r["protein_alignment"] for r in rows}),2)
        for position in (1,28,57):
            for aa in "ACDEFGHIKLMNPQRSTVWY":
                report = self.binding._predict(model,position,aa)
                for col in range(8):
                    values = [report["psam"][b][col] for b in "ACGT"]
                    self.assertEqual(max(values),1)
                    self.assertTrue(all(.01 <= v <= 1 for v in values))

    def test_model_and_substitution_table_cannot_be_replaced(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            for name in ("homeodomain.json","BLOSUM62.json"):
                (directory/name).write_bytes((EXAMPLE.parent/"data"/name).read_bytes())
            self.binding.load_model(directory)
            (directory/"homeodomain.json").write_bytes((directory/"homeodomain.json").read_bytes()+b" ")
            with self.assertRaisesRegex(ValueError,"identity"):
                self.binding.load_model(directory)
            (directory/"homeodomain.json").write_bytes((EXAMPLE.parent/"data/homeodomain.json").read_bytes())
            table = json.loads((directory/"BLOSUM62.json").read_text())
            table["values"][0][0] += 1
            (directory/"BLOSUM62.json").write_text(json.dumps(table))
            with self.assertRaisesRegex(ValueError,"identity"):
                self.binding.load_model(directory)

    def test_shipped_source_files_match_provenance(self):
        directory = EXAMPLE.parent/"data"
        provenance = json.loads((directory/"provenance.json").read_text())
        for name,pin in provenance["files"].items():
            raw = (directory/name).read_bytes()
            self.assertEqual(len(raw),pin["bytes"],name)
            self.assertEqual(hashlib.sha256(raw).hexdigest(),pin["sha256"],name)

    def test_position_and_residue_contract_rejects_unsupported_inputs(self):
        model = self.binding.load_model()
        for position in (0,58,1.0,True):
            with self.subTest(position=position),self.assertRaises(ValueError):
                self.binding._predict(model,position,"A")
        for aa in ("","AA","X","-","*","a",None):
            with self.subTest(aa=aa),self.assertRaises(ValueError):
                self.binding._predict(model,1,aa)

    def compile_control(self, directory, name, residue_codon="GCT", label="same-label"):
        # an explicitly synthetic peptide tests dataflow, with one modelled residue.
        sequence = "ATG"+"GCT"*27+residue_codon+"GCT"*29+"TAA"
        source = directory/(name+".fasta")
        source.write_text(f">{label}\n{sequence}\n")
        contract = {"format":"epigenesis.coding-contract","version":1,
                    "source_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),"record_id":label,
                    "start":0,"end":len(sequence),"strand":1,"genetic_code":1,"initiation":"ATG","edit":None}
        coding = self.binding.coding
        receipt = coding.compile_coding(source,contract,directory/name)
        selection = {"format":"epigenesis.binding-selection","version":1,
                     "model_sha256":self.binding.PINS["homeodomain.json"],"coding_contract_sha256":receipt,
                     "reference_alignment":"A"*57,"peptide_positions":list(range(1,58)),
                     "conditioning_position":28}
        return directory/name,selection

    def test_actual_coding_tensor_drives_binding_and_synonymous_recoding_is_invariant(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            reports = []
            for name,codon,label in [("wt","GCT","same-label"),("syn","GCC","same-label"),
                                      ("mut","TGG","same-label"),("rename","GCT","renamed")]:
                module,selection = self.compile_control(directory,name,codon,label)
                report = self.binding.consume(module,selection,self.binding.coding.digest(selection))
                self.assertEqual(report["alignment"][27],"W" if name == "mut" else "A")
                self.assertFalse(report["biological_acceptance"])
                self.assertEqual(report["coding"]["contract_sha256"],selection["coding_contract_sha256"])
                reports.append(report)
            self.assertEqual(reports[0]["prediction"],reports[1]["prediction"])
            self.assertEqual(reports[0]["prediction"],reports[3]["prediction"])
            self.assertNotEqual(reports[0]["prediction"]["psam"],reports[2]["prediction"]["psam"])

    def test_receipts_domain_coordinates_and_unmodelled_residues_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            module,selection = self.compile_control(Path(tmp),"input")
            receipt = self.binding.coding.digest(selection)
            self.binding.consume(module,selection,receipt)
            changed = {**selection,"conditioning_position":29}
            with self.assertRaisesRegex(ValueError,"receipt"):
                self.binding.consume(module,changed,receipt)
            for change in ({"version":1.0},{"conditioning_position":28.0},{"extra":True},
                           {"model_sha256":"0"*64},{"coding_contract_sha256":"0"*64},
                           {"reference_alignment":"C"+"A"*56},
                           {"peptide_positions":[0]+list(range(2,58))},
                           {"peptide_positions":[1,1]+list(range(3,58))},
                           {"peptide_positions":list(range(1,57))+[59]},
                           {"peptide_positions":[None]+list(range(2,58))}):
                candidate = {**selection,**change}
                with self.subTest(change=change),self.assertRaises(ValueError):
                    self.binding.consume(module,candidate,self.binding.coding.digest(candidate))
            path = module/"prediction.json"
            report = json.loads(path.read_text())
            report["peptide"] = "wrong stored prediction"
            path.write_text(json.dumps(report))
            loaded = self.binding.consume(module,selection,receipt)
            self.assertEqual(loaded["coding"]["peptide"],"M"+"A"*57)

    def test_discontinuous_domain_map_checks_gaps_and_native_insertions(self):
        with tempfile.TemporaryDirectory() as tmp:
            module,selection = self.compile_control(Path(tmp),"gapped")
            selection["reference_alignment"] = "-"+"A"*56
            selection["peptide_positions"] = [None]+list(range(2,58))
            report = self.binding.consume(module,selection,self.binding.coding.digest(selection))
            self.assertEqual(report["alignment"],selection["reference_alignment"])
            selection["conditioning_position"] = 1
            with self.assertRaisesRegex(ValueError,"mapped"):
                self.binding.consume(module,selection,self.binding.coding.digest(selection))

    def test_natural_joined_nkx2_2_cds_matches_external_protein_and_drives_variant_prediction(self):
        data = EXAMPLE.parent/"data"
        fasta = data/"NC_000068.8_147025815_147028038.fasta"
        oracle_raw = (data/"Nkx2-2_NP_035049.1.fasta").read_bytes()
        self.assertEqual(hashlib.sha256(oracle_raw).hexdigest(),
                         "f2f6f8a1d6cfa68b3106817998ceaf4712be9788b3ed4dbb514e79882bd6f44c")
        oracle = "".join(oracle_raw.decode().splitlines()[1:])
        reference = json.loads((data/"nkx2-2-wt.json").read_text())
        with tempfile.TemporaryDirectory() as tmp:
            reports = []
            for name,edit in [("wt",None),("p156l",{"start":466,"deleted":"C","inserted":"T"}),
                              ("synonymous",{"start":467,"deleted":"C","inserted":"T"})]:
                directory = Path(tmp)/name
                receipt = self.binding.coding.compile_coding(fasta,{**reference,"edit":edit},directory)
                selection = {"format":"epigenesis.binding-selection","version":1,
                    "model_sha256":self.binding.PINS["homeodomain.json"],"coding_contract_sha256":receipt,
                    "reference_alignment":oracle[128:185],"peptide_positions":list(range(128,185)),
                    "conditioning_position":28}
                report = self.binding.consume(directory,selection,self.binding.coding.digest(selection))
                expected = oracle[:155]+"L"+oracle[156:] if name == "p156l" else oracle
                self.assertEqual(report["coding"]["peptide"],expected)
                self.assertEqual(report["coding"]["stop_offset"],819)
                self.assertEqual(report["dna_positions"],["N1","N2","T3","D4","A5","Y6","N7","N8"])
                reports.append(report)
            self.assertEqual(reports[0]["prediction"],reports[2]["prediction"])
            self.assertNotEqual(reports[0]["prediction"]["psam"],reports[1]["prediction"]["psam"])

    def test_selection_is_snapshotted_before_consuming_coding_module(self):
        with tempfile.TemporaryDirectory() as tmp:
            module,selection = self.compile_control(Path(tmp),"snapshot")
            receipt = self.binding.coding.digest(selection)
            expected = self.binding.consume(module,selection,receipt)
            consume = self.binding.coding.consume
            def mutate_caller(*args):
                selection["peptide_positions"][27] = 0
                return consume(*args)
            with mock.patch.object(self.binding.coding,"consume",side_effect=mutate_caller):
                result = self.binding.consume(module,selection,receipt)
            self.assertEqual(result,expected)
            self.assertNotEqual(self.binding.coding.digest(selection),receipt)
            with self.assertRaises(ValueError):
                self.binding.consume(module,selection,receipt)

    def test_loaded_source_identity_survives_later_file_replacement(self):
        with tempfile.TemporaryDirectory() as tmp:
            for example in ("coding-consequences","neural-binding"):
                shutil.copytree(ROOT/"examples"/example,Path(tmp)/"examples"/example,
                                ignore=shutil.ignore_patterns("__pycache__"))
            source = Path(tmp)/"examples/neural-binding/binding.py"
            spec = importlib.util.spec_from_file_location("binding_identity_probe",source)
            candidate = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(candidate)
            self.binding = candidate
            module,selection = self.compile_control(Path(tmp),"original")
            receipt = candidate.coding.digest(selection)
            expected = candidate.consume(module,selection,receipt)
            for path in (source,Path(tmp)/"examples/coding-consequences/coding.py"):
                with self.subTest(source=path.name):
                    path.write_bytes(path.read_bytes()+b"\n# deliberately replaced after import\n")
                    actual = candidate.consume(module,selection,receipt)
                    self.assertEqual(actual,expected)


if __name__ == "__main__":
    unittest.main()

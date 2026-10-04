from __future__ import annotations

import base64
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

from brainc.development_bundle import compile_development
from brainc.source import FASTA_PROFILE, GENBANK_PROFILE, compile_source
from brainc.v2 import inline_storage, pack, save
from brainc.v2._common import seal
from brainc.validator_development import validate_development_paths


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT/"examples/coding-consequences/coding.py"
DATA = EXAMPLE.parent/"data"


class CodingConsequencesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location("coding_consequences", EXAMPLE)
        cls.model = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.model)
        cls.workspace = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.workspace.cleanup)
        cls.root = Path(cls.workspace.name)

    def contract(self, raw, record="same-label", start=0, end=None, strand=1, edit=None):
        if end is None:
            end = len(b"".join(raw.splitlines()[1:]))
        return {"format":"epigenesis.coding-contract", "version":1,
                "source_sha256":hashlib.sha256(raw).hexdigest(), "record_id":record,
                "start":start, "end":end, "strand":strand, "genetic_code":1,
                "initiation":"ATG", "edit":edit}

    def compile(self, raw, contract, name):
        fasta = self.root/(name+".fasta")
        fasta.write_bytes(raw)
        output = self.root/name
        receipt = self.model.compile_coding(fasta, contract, output)
        loaded = self.model.consume(output, receipt)
        return output, receipt, loaded

    def genomic(self, newer=False, edit=None):
        name, record, start, end = (
            ("NC_141042.1_40458511_40460773.fasta", "NC_141042.1:40458511-40460773", 992, 1940)
            if newer else
            ("NC_007133.7_37347793_37350078.fasta", "NC_007133.7:37347793-37350078", 1015, 1963))
        raw = (DATA/name).read_bytes()
        return raw, self.contract(raw,record,start,end,-1,edit)

    def test_all_standard_codons_first_stop_and_bounded_partial_template(self):
        # NCBI code 1, read by amino-acid groups independently of gc.prt ordering.
        groups = {"F":"TTT TTC", "L":"TTA TTG CTT CTC CTA CTG", "I":"ATT ATC ATA", "M":"ATG",
                  "V":"GTT GTC GTA GTG", "S":"TCT TCC TCA TCG AGT AGC", "P":"CCT CCC CCA CCG",
                  "T":"ACT ACC ACA ACG", "A":"GCT GCC GCA GCG", "Y":"TAT TAC", "*":"TAA TAG TGA",
                  "H":"CAT CAC", "Q":"CAA CAG", "N":"AAT AAC", "K":"AAA AAG", "D":"GAT GAC",
                  "E":"GAA GAG", "C":"TGT TGC", "W":"TGG", "R":"CGT CGC CGA CGG AGA AGG",
                  "G":"GGT GGC GGA GGG"}
        covered = set()
        for aa,codons in groups.items():
            for codon in codons.split():
                covered.add(codon)
                result = self.model.translate("ATG"+codon+"TAA"+"TGG")
                self.assertEqual(result["peptide"], "M" if aa == "*" else "M"+aa)
                self.assertEqual(result["stop_offset"], 3 if aa == "*" else 6)
                self.assertEqual(result["stop_codon"], codon if aa == "*" else "TAA")
                self.assertEqual(result["status"], "terminated")
        self.assertEqual(covered,{"".join(c) for c in itertools.product("ACGT",repeat=3)})
        for suffix in ("", "A", "AC"):
            result = self.model.translate("ATGGCT"+suffix)
            self.assertEqual(result["peptide"],"MA")
            self.assertEqual(result["status"],"no-stop-in-template")
            self.assertEqual(result["trailing_bases"],suffix)
            self.assertIsNone(result["stop_offset"])
        for bad in ("", "AT", "TTGAAA", "ATGNAA", "atgtaa", "ATGTAA\n"):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                self.model.translate(bad)
        for size in (self.model.MAX_TEMPLATE_BASES-1, self.model.MAX_TEMPLATE_BASES):
            self.model.translate("ATG"+"A"*(size-3))
        with self.assertRaises(ValueError):
            self.model.translate("ATG"+"A"*(self.model.MAX_TEMPLATE_BASES-2))

    def test_compiled_template_budget_and_authority_identity(self):
        length = self.model.MAX_TEMPLATE_BASES-1
        sequence = "ATG"+"A"*(length-6)+"TAA"
        raw = f">same-label\n{sequence}\n".encode()
        for extra in (0,1,2):
            contract = self.contract(raw,edit=None if extra == 0 else {"start":length,"deleted":"","inserted":"A"*extra})
            if extra < 2:
                _,_,report = self.compile(raw,contract,f"budget-{extra}")
                self.assertEqual(report["template_bases"],length+extra)
            else:
                with self.assertRaisesRegex(ValueError,"size"):
                    self.compile(raw,contract,f"budget-{extra}")
        original = self.model.read_regular_file
        def altered(path, **kwargs):
            raw = original(path,**kwargs)
            return raw+b"\n" if Path(path).name == "ncbi-gc.prt" else raw
        with mock.patch.object(self.model,"read_regular_file",side_effect=altered):
            with self.assertRaisesRegex(ValueError,"authority identity"):
                self.model.translate("ATGTAA")

    def test_compiled_tensor_is_nucleotides_and_labels_cannot_select_peptides(self):
        reports = []
        for index,(record,sequence) in enumerate((("same-label","ATGGCTTAA"),("same-label","ATGGCCTAA"),
                                                  ("same-label","ATGGAATAA"),("other-label","ATGGCTTAA"))):
            raw = f">{record}\n{sequence}\n".encode()
            output,receipt,report = self.compile(raw,self.contract(raw,record),f"labels-{index}")
            reports.append(report)
            module = json.loads((output/"development/development_module.json").read_text())["module"]
            tensor = next(t for t in module["tensors"] if t["id"] == "nucleotides")
            self.assertEqual(base64.b64decode(tensor["storage"]["data"]),pack("u64",list(sequence.encode())))
            self.assertEqual({t["id"] for t in module["tensors"]},{"count","nucleotides"})
            self.assertFalse(report["biological_acceptance"])
            self.assertEqual(report["contract_sha256"],receipt)
            with self.assertRaises(FileExistsError):
                self.model.compile_coding(output/"reference.fasta",self.contract(raw,record),output)
        self.assertEqual([r["peptide"] for r in reports],["MA","MA","ME","MA"])

    def test_orientation_flanks_and_record_selection(self):
        for index,(raw,record,start,end,strand) in enumerate((
            (b">r\nATGGAATAA\n","r",0,9,1),
            (b">r\nCCCTTATTCCATGG\n","r",3,12,-1),
            (b">decoy\nATGGCTTAA\n>r\nNNATGGAATAAN\n","r",2,11,1))):
            _,_,report = self.compile(raw,self.contract(raw,record,start,end,strand),f"orientation-{index}")
            self.assertEqual(report["peptide"],"ME")

    def joined(self, raw, segments, strand=1, edit=None):
        contract = self.contract(raw,strand=strand,edit=edit)
        del contract["start"],contract["end"]
        return {**contract,"version":2,"segments":[{"start":lo,"end":hi} for lo,hi in segments]}

    def test_joined_cds_preserves_introns_and_translates_split_codons_on_both_strands(self):
        # CDS ATGGA + ATAA = ME*, with a phase-two split codon across the intron.
        sequence = "CCATGGACCCTAGATAAGG"
        for strand in (1,-1):
            genomic = sequence if strand == 1 else self.model._reverse(sequence)
            spans = [(2,7),(13,17)] if strand == 1 else [(2,6),(12,17)]
            raw = f">same-label\n{genomic}\n".encode()
            contract = self.joined(raw,spans,strand)
            output,_,report = self.compile(raw,contract,f"joined-{strand}")
            self.assertEqual(report["peptide"],"ME")
            self.assertEqual(report["stop_offset"],6)
            self.assertEqual((output/"sequence.fasta").read_bytes(),raw)
            # c.6A>T changes the split GAA codon to GAT, preserving the whole intron.
            edited = self.joined(raw,spans,strand,{"start":5,"deleted":"A","inserted":"T"})
            output,_,report = self.compile(raw,edited,f"joined-edit-{strand}")
            self.assertEqual(report["peptide"],"MD")
            actual = b"".join((output/"sequence.fasta").read_bytes().splitlines()[1:]).decode()
            oriented = actual if strand == 1 else self.model._reverse(actual)
            self.assertEqual(oriented,"CCATGGACCCTAGTTAAGG")

    def test_joined_cds_rejects_invalid_segments_and_ambiguous_junction_edits(self):
        raw = b">same-label\nCCATGGACCCTAGATAAGG\n"
        valid = self.joined(raw,[(2,7),(13,17)])
        for index,segments in enumerate(([],[{"start":2,"end":7},{"start":6,"end":17}],
                [{"start":13,"end":17},{"start":2,"end":7}], [{"start":2.0,"end":7}],
                [{"start":2,"end":99}], [{"start":2,"end":2}], [{"start":2,"end":7,"extra":0}])):
            with self.subTest(segments=segments),self.assertRaises(ValueError):
                self.compile(raw,{**valid,"segments":segments},f"joined-bad-segments-{index}")
        for index,edit in enumerate(({"start":4,"deleted":"AA","inserted":"T"},
                                    {"start":5,"deleted":"","inserted":"G"})):
            with self.subTest(edit=edit),self.assertRaisesRegex(ValueError,"junction"):
                self.compile(raw,{**valid,"edit":edit},f"joined-junction-{index}")

    def test_edits_check_reference_bases_and_reverse_strand_coordinates(self):
        edit = {"start":3,"deleted":"GCT","inserted":"GAA"}
        for index,(raw,strand,start,end) in enumerate(((b">same-label\nCCATGGCTTAAGG\n",1,2,11),
                                                      (b">same-label\nCCTTAAGCCATGG\n",-1,2,11))):
            out,_,report = self.compile(raw,self.contract(raw,start=start,end=end,strand=strand,edit=edit),f"edit-{index}")
            self.assertEqual(report["peptide"],"ME")
            actual = b"".join((out/"sequence.fasta").read_bytes().splitlines()[1:])
            self.assertEqual(actual,b"CCATGGAATAAGG" if strand == 1 else b"CCTTATTCCATGG")
        raw = b">same-label\nATGGCTTAA\n"
        for index,edit in enumerate(({"start":3,"deleted":"CCC","inserted":""},
                                     {"start":9,"deleted":"T","inserted":""},
                                     {"start":-1,"deleted":"","inserted":"A"},
                                     {"start":3,"deleted":"","inserted":"N"},
                                     {"start":True,"deleted":"","inserted":"A"})):
            with self.subTest(edit=edit), self.assertRaises(ValueError):
                self.compile(raw,self.contract(raw,edit=edit),f"invalid-edit-{index}")

    def test_closed_contract_rejects_unsupported_or_mismatched_inputs(self):
        raw = b">same-label\nATGGCTTAA\n"
        changes = ({"strand":0},{"strand":True},{"start":True},{"start":-1},{"end":10},{"start":4,"end":3},
                   {"genetic_code":2},{"genetic_code":True},{"initiation":"TTG"},{"record_id":"absent"},
                   {"source_sha256":"0"*64},{"extra":"ignored?"},{"version":True})
        for index,change in enumerate(changes):
            with self.subTest(change=change), self.assertRaises(ValueError):
                self.compile(raw,{**self.contract(raw),**change},f"invalid-contract-{index}")

    def test_reference_is_a_complete_cds_before_editing(self):
        for index,sequence in enumerate(("ATGAAA","ATGTAAGCT","ATGTAAC","ATGTAACC")):
            raw = f">same-label\n{sequence}\n".encode()
            with self.subTest(sequence=sequence), self.assertRaisesRegex(ValueError,"complete reference CDS"):
                self.compile(raw,self.contract(raw),f"incomplete-reference-{index}")

    def test_compile_rejects_integral_floats_before_canonicalization(self):
        raw = b">same-label\nATGGCTTAA\n"
        for index,change in enumerate(({"version":1.0},{"end":9.0},{"strand":1.0},{"genetic_code":1.0},
            {"edit":{"start":3.0,"deleted":"GCT","inserted":"GAA"}})):
            with self.subTest(change=change), self.assertRaises(ValueError):
                self.compile(raw,{**self.contract(raw),**change},f"float-contract-{index}")

    def test_pinned_primary_evidence_and_genomic_orientation(self):
        provenance = json.loads((DATA/"provenance.json").read_text())
        for item in provenance["files"]:
            raw = (DATA/item["file"]).read_bytes()
            self.assertEqual(len(raw),item["bytes"])
            self.assertEqual(hashlib.sha256(raw).hexdigest(),item["sha256"])
        for newer in (False,True):
            raw,contract = self.genomic(newer)
            stem = "NC_141042.1_40458511_40460773" if newer else "NC_007133.7_37347793_37350078"
            genbank = (DATA/(stem+".gb")).read_text()
            self.assertIn(" DNA ",genbank.splitlines()[0])
            self.assertIn('/mol_type="genomic DNA"',genbank)
            self.assertIn(f'CDS             complement({contract["start"]+1}..{contract["end"]})',genbank)
            origin = re.sub('[^acgt]', '',genbank.split("ORIGIN",1)[1].split("//",1)[0]).upper()
            self.assertEqual(b"".join(raw.splitlines()[1:]).decode(),origin)
            with self.assertRaisesRegex(ValueError,"INSDC026"):
                compile_source(GENBANK_PROFILE,{"genbank":DATA/(stem+".gb")},parameters={})
        with self.assertRaises(ValueError):
            compile_source(GENBANK_PROFILE,{"genbank":DATA/"NM_213118.1.gb"},parameters={})

    def test_same_size_indels_change_actual_tail_and_report_missing_stop(self):
        raw = b">same-label\nATGAAACCCTAA\n"
        for offset,deleted,peptide in ((3,"A","MNP"),(6,"C","MKP")):
            _,_,report = self.compile(raw,self.contract(raw,edit={"start":offset,"deleted":deleted,"inserted":""}),f"frameshift-{offset}")
            self.assertEqual(report["peptide"],peptide)
            self.assertEqual(report["template_bases"],11)
            self.assertEqual(report["status"],"no-stop-in-template")
            self.assertEqual(report["trailing_bases"],"AA")

    def test_external_receipt_binds_identical_orfs_and_noncoding_flanks(self):
        raw = b">same-label\nATGGCTTAAATGGCTTAA\n"
        _,receipt,report = self.compile(raw,self.contract(raw,end=9),"first-orf")
        other,_,result = self.compile(raw,self.contract(raw,start=9,end=18),"second-orf")
        self.assertEqual(report["peptide"],result["peptide"])
        with self.assertRaisesRegex(ValueError,"contract identity"):
            self.model.consume(other,receipt)
        changed = raw[:-2]+b"C\n"
        other,_,result = self.compile(changed,self.contract(changed,end=9),"noncoding-change")
        self.assertEqual(report["peptide"],result["peptide"])
        with self.assertRaisesRegex(ValueError,"contract identity"):
            self.model.consume(other,receipt)

    def test_published_alleles_and_reference_haplotypes(self):
        # Gong et al. 2020 Fig S1C/D: local sequenced edits and predicted -4 protein.
        nm = (DATA/"NM_213118.1.gb").read_text()
        wt = "".join(re.search(r'/translation="([^"]+)"',nm).group(1).split())
        edits = {"WT":None,"minus4":{"start":208,"deleted":"AGCG","inserted":""},
                 "equivalent":{"start":209,"deleted":"GCGA","inserted":""},
                 "minus11":{"start":213,"deleted":"CTCGGGGCCGA","inserted":""},
                 "minus3":{"start":211,"deleted":"GAC","inserted":""},
                 "plus3":{"start":213,"deleted":"","inserted":"ATA"}}
        results = {}
        for newer in (False,True):
            for name,edit in edits.items():
                raw,contract = self.genomic(newer,edit)
                out,_,report = self.compile(raw,contract,f"genomic-{newer}-{name}")
                results[newer,name] = report["peptide"]
                if name == "WT":
                    self.assertEqual(report["peptide"],wt)
                elif name in ("minus4","equivalent"):
                    self.assertEqual(report["peptide"],wt[:69]+"NSGPSGNFCPRARSDHSSTKPNAFGLCT")
                    self.assertEqual(report["stop_offset"],291)
                elif name == "minus3":
                    self.assertEqual(report["peptide"],wt[:70]+wt[71:])
                elif name == "plus3":
                    self.assertEqual(report["peptide"],wt[:71]+"I"+wt[71:])
            a = self.root/f"genomic-{newer}-minus4/sequence.fasta"
            b = self.root/f"genomic-{newer}-equivalent/sequence.fasta"
            self.assertEqual(a.read_bytes(),b.read_bytes())
        older,newer = results[False,"minus11"],results[True,"minus11"]
        self.assertEqual((len(older),len(newer)),(241,241))
        self.assertEqual([(i+1,a,b) for i,(a,b) in enumerate(zip(older,newer)) if a != b],[(116,"G","R"),(230,"A","T")])

    def test_source_contract_and_coherent_tensor_substitution_fail(self):
        raw = b">same-label\nATGGCTTGGGCTTAA\n"
        edit = {"start":6,"deleted":"TGG","inserted":"TGA"}
        original,receipt,_ = self.compile(raw,self.contract(raw,edit=edit),"trusted")
        for attack in ("source","contract","tensor","after-stop","rule"):
            with self.subTest(attack=attack):
                out = self.root/("attack-"+attack)
                shutil.copytree(original,out)
                if attack == "source":
                    (out/"reference.fasta").write_bytes(raw.replace(b"GCT",b"GAA"))
                elif attack == "contract":
                    contract = json.loads((out/"contract.json").read_text())
                    contract["strand"] = -1
                    (out/"contract.json").write_text(json.dumps(contract))
                else:
                    response = json.loads((out/"response.json").read_text())
                    target = json.loads((out/"target.json").read_text())
                    policy = json.loads((out/"policy.json").read_text())
                    if attack in ("tensor","after-stop"):
                        sequence = b"ATGGAATGAGCTTAA" if attack == "tensor" else b"ATGGCTTGAGAATAA"
                        next(t for t in response["outputs"] if t["id"] == "nucleotides")["storage"] = inline_storage(pack("u64",list(sequence)))
                        save(seal({k:v for k,v in response.items() if k != "artifact_sha256"}),out/"changed-response.json")
                    else:
                        # change a real rule while preserving compiler-valid closure.
                        target["contract"]["rules"][0]["semantics_sha256"] = "0"*64
                        target.pop("artifact_sha256")
                        target["contract_sha256"] = self.model.digest(target["contract"])
                        save(seal(target),out/"target.json.new")
                        (out/"target.json.new").replace(out/"target.json")
                        policy["target"]["contract_sha256"] = target["contract_sha256"]
                        save(seal({k:v for k,v in policy.items() if k != "artifact_sha256"}),out/"policy.json.new")
                        (out/"policy.json.new").replace(out/"policy.json")
                    source = compile_source(FASTA_PROFILE,{"sequence":out/"sequence.fasta"},parameters={"wrapper":"identity"})
                    shutil.rmtree(out/"development")
                    compile_development(source,out/"manifest.json",out/"request.json",
                        out/("changed-response.json" if attack in ("tensor","after-stop") else "response.json"),
                        out/"policy.json",out/"target.json").save(out/"development")
                    self.assertTrue(validate_development_paths(out/"source",out/"development",source_inputs={"sequence":out/"sequence.fasta"})["valid"])
                with self.assertRaises(ValueError):
                    self.model.consume(out,receipt)
        alternate = raw.replace(b"GCT",b"GAA")
        other,other_receipt,_ = self.compile(alternate,self.contract(alternate),"coherent-alternative")
        self.assertNotEqual(other_receipt,receipt)
        with self.assertRaisesRegex(ValueError,"contract identity"):
            self.model.consume(other,receipt)

    def test_public_cli_and_external_receipt(self):
        raw,contract = self.genomic(edit={"start":208,"deleted":"AGCG","inserted":""})
        fasta,selection = self.root/"cli.fasta",self.root/"cli-contract.json"
        fasta.write_bytes(raw)
        selection.write_text(json.dumps(contract))
        output = self.root/"cli"
        command = [sys.executable,str(EXAMPLE),"compile","--fasta",str(fasta),"--contract",str(selection),"--output",str(output)]
        result = subprocess.run(command,capture_output=True,text=True,timeout=30)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        report = json.loads(result.stdout)
        replay = subprocess.run([sys.executable,str(EXAMPLE),"consume",str(output),"--contract-sha256",report["contract_sha256"]],capture_output=True,text=True,timeout=30)
        self.assertEqual(replay.returncode,0,replay.stdout+replay.stderr)
        self.assertEqual(json.loads(replay.stdout),report)
        self.assertEqual(len(report["peptide"]),97)


if __name__ == "__main__":
    unittest.main()

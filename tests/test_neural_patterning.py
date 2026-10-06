from __future__ import annotations

import importlib.util
import itertools
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples/neural-patterning/patterning.py"


class NeuralPatterningTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location("neural_patterning", EXAMPLE)
        cls.m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.m)
        cls.data = cls.m.load_data()

    def test_published_single_substitution_and_strand(self):
        # Hallikas 2006 Table S1: D10/D11 and D16/D17, not Peterson PBM E-scores.
        for protein, expected in (("GLI2", .04063440697383719/.9593655930261629),
                                  ("GLI3", .03409141625199168/.9659085837480084)):
            self.assertEqual(self.m.profile("GACCACCCA", protein)["profile_ratio"], 1)
            prediction = self.m.profile("GAACACCCA", protein)
            self.assertEqual(prediction["profile_ratio"], expected)
            self.assertEqual(prediction["changed_position_1based"], 3)
        self.assertEqual(self.m.reverse_complement("TGGGTGTTC"), "GAACACCCA")

    def test_unsupported_sequence_inference_fails(self):
        for sequence in ("GAATACCCA", "NACCACCCA", "GACCACCC", "gaccaccca", "GTCCACCCA"):
            with self.subTest(sequence=sequence), self.assertRaises(ValueError):
                self.m.profile(sequence, "GLI2")
        with self.assertRaises(ValueError):
            self.m.profile("GACCACCCA", "unmeasured")

    def test_all_measured_single_base_profiles(self):
        binding=self.data["binding"]
        for protein,weights in binding["profiles"].items():
            reference=binding["consensus"]
            for i,base in enumerate(reference):
                for alternate in "ACGT":
                    site=reference[:i]+alternate+reference[i+1:]
                    with self.subTest(protein=protein,position=i,base=alternate):
                        if weights[alternate][i] == 0:
                            with self.assertRaises(ValueError):
                                self.m.profile(site,protein)
                        else:
                            self.assertEqual(self.m.profile(site,protein)["profile_ratio"],
                                             weights[alternate][i]/weights[base][i])

    def test_occupancy_matches_independent_microstate_enumeration(self):
        # Cohen 2014 Eq9: two independent sites for each repressor, one Gli site.
        y, signal, q = (.21, .43, .17, .31), .37, .19
        weights = [4.8*y[0], 27.1*y[1], 47.1*y[3]]
        off = 0.0
        for gli_weight in (1.0, q*37.3*signal, q*37.3*(1-signal)):
            for sites in itertools.product((0, 1), repeat=2*len(weights)):
                off += gli_weight * math.prod(weights[j//2]**occupied for j, occupied in enumerate(sites))
        on = .572324*(1+37.3)*(1 + 10*q*37.3*signal)
        observed = self.m.occupancy(self.data["model"], y, signal, q)[2]
        self.assertAlmostEqual(observed, on/(on+off), places=15)

    def test_affinity_neutral_point_and_sign_reversal(self):
        model = self.data["model"]
        for state in ((0,0,0,0), (.2,.5,.1,.7), (1,1,1,1)):
            for q in (.01,.3,.8):
                for signal, sign in ((.01,1),(.1,0),(.9,-1)):
                    reference = self.m.occupancy(model,state,signal,1)[2]
                    changed = self.m.occupancy(model,state,signal,q)[2]
                    if sign:
                        self.assertGreater(sign*(changed-reference),0)
                    else:
                        self.assertAlmostEqual(changed,reference,places=15)
                ablated = self.m.occupancy(model,state,.01,q,"activation-only")[2]
                self.assertLess(ablated,self.m.occupancy(model,state,.01,1)[2])

    def test_published_baseline_against_independent_replay(self):
        # Untuned Herrera-Delgado 2020 Appendix J; separate transcription and solver.
        reference = self.data["reproduction"]["states"]
        for row in reference:
            actual = self.m.simulate([(100,math.exp(-row["position"]/.15))])[-1][1:]
            for a,b in zip(actual,row["state"]):
                self.assertLess(abs(a-b),2e-8)

    def test_numerical_refinement_and_history(self):
        coarse = self.m.simulate([(2.137,.6),(3.173,.025)],dt=.01)
        fine = self.m.simulate([(2.137,.6),(3.173,.025)],dt=.005)
        first = self.m.simulate([(2.137,.6)],dt=.005)
        second = self.m.simulate([(3.173,.025)],initial=first[-1][1:],dt=.005)
        self.assertEqual(fine[-1][1:],second[-1][1:])
        self.assertAlmostEqual(fine[-1][0],5.31)
        self.assertLess(max(abs(a-b) for a,b in zip(coarse[-1],fine[-1])),2e-5)
        reverse = self.m.simulate([(3.173,.025),(2.137,.6)],dt=.02)
        self.assertGreater(max(abs(a-b) for a,b in zip(fine[-1][1:],reverse[-1][1:])),.1)

    def test_invalid_protocols_and_integration_budget(self):
        for protocol in ([],[(0,1)],[(-1,.3)],[(1,-1)],[(1,1.01)],[(math.inf,.3)],[(1,math.nan)]):
            with self.subTest(protocol=protocol),self.assertRaises(ValueError):
                self.m.simulate(protocol)
        for dt in (0,-1,math.inf,math.nan,1e-12):
            with self.subTest(dt=dt),self.assertRaises(ValueError):
                self.m.simulate([(1,.5)],dt=dt)
        for ratio in (0,-1,math.inf,math.nan,1.01):
            with self.subTest(ratio=ratio),self.assertRaises(ValueError):
                self.m.simulate([(1,.5)],ratio=ratio)
            for mode in ("reporter","endogenous"):
                with self.subTest(mode=mode,ratio=ratio),self.assertRaises(ValueError):
                    self.m.run_assay([(1,.5)],ratio=ratio,mode=mode)

    def test_adaptive_rejection_and_budget(self):
        reference=self.m.simulate([(2,.2)],dt=.01)
        large_step=self.m.simulate([(2,.2)],dt=1)
        self.assertLess(max(abs(a-b) for a,b in zip(reference[-1],large_step[-1])),1e-7)
        with mock.patch.object(self.m,"MAX_STEPS",1):
            with self.assertRaisesRegex(ValueError,"budget"):
                self.m.simulate([(1,.2)],dt=1)

    def compile(self,root,sequence,name):
        fasta=root/(name+".fasta")
        fasta.write_text(">site\n"+sequence+"\n")
        contract=self.m.selection(fasta,"site",0,-1,"GLI2")
        output=root/name
        receipt=self.m.compile_site(fasta,contract,output)
        return output,receipt

    def test_actual_compiled_bases_drive_effect_and_exact_rescue(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            results=[]
            for name,sequence in (("wt","TGGGTGGTC"),("mutant","TGGGTGTTC"),("rescue","TGGGTGGTC")):
                output,receipt=self.compile(root,sequence,name)
                results.append(self.m.consume(output,receipt))
            self.assertEqual(results[0],results[2])
            self.assertNotEqual(results[0]["source_artifact_sha256"],results[1]["source_artifact_sha256"])
            trajectories=[self.m.simulate([(12,.2)],ratio=r["profile"]["profile_ratio"]) for r in results]
            self.assertEqual(trajectories[0],trajectories[2])
            self.assertGreater(max(abs(a-b) for a,b in zip(trajectories[0][-1],trajectories[1][-1])),.1)
            (root/"mutant/sequence.fasta").write_text(">site\nTGGGTGGTC\n")
            with self.assertRaises(ValueError):
                self.m.consume(root/"mutant",results[1]["contract_sha256"])

    def test_coherently_resealed_wrong_nucleotides_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            output,receipt=self.compile(Path(tmp),"TGGGTGGTC","compiled")
            response=json.loads((output/"response.json").read_text())
            response.pop("artifact_sha256")
            for tensor in response["outputs"]:
                if tensor["id"]=="nucleotides":
                    tensor["storage"]=self.m.inline_storage(self.m.pack("u64",list(b"GAACACCCA")))
            self.m.save(self.m.seal(response),output/"response.json")
            (output/"development").rename(output/"honest-development")
            self.m.compile_development(output/"source",output/"manifest.json",output/"request.json",
                output/"response.json",output/"policy.json",output/"target.json").save(output/"development")
            with self.assertRaisesRegex(ValueError,"nucleotide tensor"):
                self.m.consume(output,receipt)

    def test_reporter_does_not_edit_endogenous_network(self):
        wt=self.m.run_assay([(10,.6)],ratio=1,mode="reporter")
        mutant=self.m.run_assay([(10,.6)],ratio=.04,mode="reporter")
        self.assertEqual(wt["trajectory"],mutant["trajectory"])
        self.assertNotEqual(wt["reporter_occupancy"],mutant["reporter_occupancy"])
        edited=self.m.run_assay([(10,.6)],ratio=.04,mode="endogenous")
        self.assertNotEqual(wt["trajectory"],edited["trajectory"])

    def test_named_initial_state_drives_dynamics_and_reporter(self):
        initial = {"Irx3":0., "Nkx2.2":.0001, "Olig2":.1, "Pax6":.9}
        expected = tuple(initial[g] for g in self.data["model"]["genes"])
        wt = self.m.run_assay([(.001,0)],ratio=1,mode="reporter",initial=initial)
        mutant = self.m.run_assay([(.001,0)],ratio=.04,mode="reporter",initial=initial)
        self.assertEqual(wt["initial_state"],initial)
        self.assertEqual(wt["trajectory"][0],(0.,*expected))
        self.assertEqual(wt["trajectory"],mutant["trajectory"])
        self.assertNotEqual(wt["reporter_occupancy"],mutant["reporter_occupancy"])
        # Appendix J equations permit this direction from a nonzero state at S=0.
        change = [a-b for a,b in zip(wt["trajectory"][-1][1:],expected)]
        self.assertLess(change[0],0)
        self.assertGreater(change[1],0)
        self.assertGreater(change[2],0)
        naive = self.m.run_assay([(.001,0)],ratio=1,mode="reporter")
        self.assertGreater(naive["trajectory"][-1][1],0)
        fine = self.m.run_assay([(.001,0)],ratio=1,mode="reporter",initial=initial,dt=.0005)
        self.assertLess(max(abs(a-b) for a,b in zip(fine["trajectory"][-1],wt["trajectory"][-1])),1e-10)
        for invalid in ({}, {**initial,"unknown":0}, [0,0,0,0],
                        {**initial,"Pax6":True}, {**initial,"Pax6":-1},
                        {**initial,"Pax6":math.nan}, {**initial,"Irx3":1.01}, {**initial,"Pax6":10**400}):
            with self.subTest(initial=invalid),self.assertRaises(ValueError):
                self.m.run_assay([(.001,0)],ratio=1,mode="reporter",initial=invalid)

    def test_cli_initial_state_validation_precedes_publication(self):
        initial = {"Irx3":.1,"Nkx2.2":0,"Olig2":0,"Pax6":.1}
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            protocol = root/"protocol.json"
            for index,state in enumerate((initial,{},None,{**initial,"Pax6":True})):
                output = root/str(index)
                protocol.write_text(json.dumps({"phases":[[.01,.2]],"initial_state":state}))
                run = subprocess.run([sys.executable,str(EXAMPLE),"run",
                    "--fasta",str(EXAMPLE.parent/"reference.fasta"),"--record","nkx2_2_site",
                    "--strand","-1","--mode","endogenous","--protocol",str(protocol),
                    "--output",str(output)],text=True,capture_output=True)
                if index:
                    self.assertEqual(run.returncode,2,run.stderr)
                    self.assertFalse(output.exists())
                else:
                    self.assertEqual(run.returncode,0,run.stderr)
                    assay = json.loads((output/"report.json").read_text())["assay"]
                    self.assertEqual(assay["initial_state"],initial)
                    self.assertEqual(assay["trajectory"][0],[0,*[initial[g] for g in self.data["model"]["genes"]]])

    def test_coherently_resealed_wrong_interpreter_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            output,receipt=self.compile(Path(tmp),"TGGGTGGTC","compiled")
            manifest=json.loads((output/"manifest.json").read_text())
            manifest.pop("artifact_sha256")
            manifest["model_identity"]["value"]="f"*64
            self.m.save(self.m.seal(manifest),output/"manifest.json")
            request=self.m.make_development_request(output/"source",output/"manifest.json",["count","nucleotides"])
            self.m.save(request,output/"request.json")
            response=json.loads((output/"response.json").read_text())
            response.pop("artifact_sha256")
            response["request_artifact_sha256"]=request["artifact_sha256"]
            response["model_identity"]=manifest["model_identity"]
            self.m.save(self.m.seal(response),output/"response.json")
            (output/"development").rename(output/"original-development")
            self.m.compile_development(output/"source",output/"manifest.json",output/"request.json",
                output/"response.json",output/"policy.json",output/"target.json").save(output/"development")
            with self.assertRaisesRegex(ValueError,"interpreter identity"):
                self.m.consume(output,receipt)

    def test_cli_selected_interval_protocol_and_export(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            fasta=root/"sequence.fasta"
            fasta.write_text(">other\nAAAA\n>selected\nCCTGGGTGTTCAA\n")
            protocol=root/"protocol.json"
            protocol.write_text(json.dumps({"phases":[[.137,.6],[.173,.025]]}))
            for mode in ("reporter","endogenous"):
                output=root/mode
                run=subprocess.run([sys.executable,str(EXAMPLE),"run","--fasta",str(fasta),
                    "--record","selected","--start","2","--strand","-1","--mode",mode,
                    "--protocol",str(protocol),"--output",str(output)],cwd=root,text=True,capture_output=True)
                self.assertEqual(run.returncode,0,run.stderr)
                receipt=json.loads(run.stdout)["compiled"]["contract_sha256"]
                compiled=self.m.consume(output/"compiled",receipt)
                self.assertEqual(compiled["profile"]["oriented_sequence"],"GAACACCCA")
                report=json.loads((output/"report.json").read_text())
                self.assertAlmostEqual(report["assay"]["trajectory"][-1][0],.31)
                for name,expected in json.loads((output/"files.sha256.json").read_text()).items():
                    self.assertEqual(self.m.hashlib.sha256((output/name).read_bytes()).hexdigest(),expected)
                with self.assertRaises(ValueError):
                    self.m.consume(output/"compiled","0"*64)


if __name__ == "__main__":
    unittest.main()

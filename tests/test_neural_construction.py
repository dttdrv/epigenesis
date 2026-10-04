from __future__ import annotations

import importlib.util
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import sys
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples/neural-construction/construction.py"


class NeuralConstructionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location("neural_construction_example", EXAMPLE)
        cls.model = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = cls.model
        spec.loader.exec_module(cls.model)
        cls.workspace = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.workspace.cleanup)
        cls.root = Path(cls.workspace.name)
        cls.classes = cls.model.build_backend(cls.root/"backend")
        cls.sample = cls.root/"sample"
        cls.model.develop(cls.model.default_program(),cls.classes,cls.sample,duration=1.5,seed=7)

    def activity(self):
        spec = importlib.util.spec_from_file_location("construction_activity",EXAMPLE.with_name("activity.py"))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def test_activity_has_exact_temporal_integration_and_refractory_boundary(self):
        activity = self.activity()
        # NEST iaf_psc_delta: 8 + 8*exp(-interval/10), threshold gap 15 mV.
        for interval, expected in ((1.0,[[1.0,1]]),(2.0,[])):
            result = activity.simulate([1],[],[(0.0,1,8.0),(interval,1,8.0)],duration=5)
            self.assertEqual(result["spikes"],expected)
        self.assertAlmostEqual(activity.simulate([1],[],[(0,1,8)],duration=10)["voltages"]["1"],-70+8/math.e)
        result = activity.simulate([1],[],[(0,1,15),(1.999,1,100),(2,1,15)],duration=3)
        self.assertEqual(result["spikes"],[[0,1],[2,1]])
        events = [(0,1,20),(0,1,-10)]
        self.assertEqual(activity.simulate([1],[],events,duration=1),activity.simulate([1],[],list(reversed(events)),duration=1))
        self.assertEqual(activity.simulate([1],[],events,duration=1)["spikes"],[])

    def test_activity_uses_grown_contacts_and_lesions_block_only_transmission(self):
        activity = self.activity()
        report = activity.assay(self.sample,self.model.verify_trace)
        self.assertGreater(report["contacts"],0)
        self.assertTrue(report["transmission_controls_pass"])
        self.assertFalse(report["biological_acceptance"])
        self.assertTrue(any(trial["downstream_cells"] for trial in report["trials"]))
        for trial in report["trials"]:
            self.assertEqual(trial["lesioned_spikes"],[[1.0,trial["stimulated_cell"]]])
        with self.assertRaisesRegex(ValueError,"event budget"):
            activity.simulate([1,2],[(1,2),(2,1)],[(0,1,15)],duration=100,max_events=3)

    def test_recipe_roundtrip_rejects_missing_extra_and_nonfinite_genes(self):
        program = self.model.default_program()
        sequence = self.model.encode_recipe(program)
        self.assertEqual(set(sequence), set(b"ACGT"))
        self.assertEqual(self.model.decode_recipe(sequence), program)
        for wrong in (sequence[:-1], sequence+b"A", b"N"+sequence[1:]):
            with self.assertRaises(ValueError):
                self.model.decode_recipe(wrong)
        for value in (math.nan, math.inf, -1.0):
            with self.assertRaises(ValueError):
                self.model.encode_recipe(program._replace(volume_rate=value))

    def test_module_contains_one_founder_program_and_no_imported_topology(self):
        program = self.model.default_program()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            recipe = root/"recipe.dna"
            recipe.write_bytes(self.model.encode_recipe(program))
            self.model.compile_recipe(recipe, root/"compiled")
            loaded, report = self.model.load_compiled(root/"compiled")
            self.assertTrue(report["valid"])
            self.assertEqual(loaded, program)
            module = json.loads((root/"compiled/development/development_module.json").read_text())["module"]
            self.assertEqual({t["id"] for t in module["tensors"]}, {"count", "program"})
            with self.assertRaises(FileExistsError):
                self.model.compile_recipe(recipe, root/"compiled")
            recipe_in_bundle = root/"compiled/recipe.dna"
            sequence = recipe_in_bundle.read_bytes()
            recipe_in_bundle.write_bytes(b"C"+sequence[1:])
            with self.assertRaises(ValueError):
                self.model.load_compiled(root/"compiled")

    def test_backend_uses_pinned_source_and_real_physics(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            program = self.model.default_program()._replace(division_signal=2.0, branch_rate=0.0)
            report = self.model.develop(program, self.classes, root/"run", duration=0.3, seed=7)
            self.assertEqual(report["initial_cells"], 1)
            self.assertGreater(report["cells"], 1)
            self.assertGreater(report["segments"], 0)
            self.assertTrue(report["physics"])
            self.assertTrue(report["diffusion"])
            self.assertFalse(report["biological_acceptance"])
            self.model.verify_trace(root/"run")
            # after growth, two further divisions must still conserve volume.
            program = program._replace(division_signal=4.0)
            report = self.model.develop(program, self.classes, root/"grown", duration=1.5, seed=7)
            self.assertEqual(report["cells"], 4)

    def test_trace_rejects_injected_cells_wrong_counts_and_undeclared_actions(self):
        events = [json.loads(line) for line in (self.sample/"trace.jsonl").read_text().splitlines()]
        for mutation in ("invented-cell", "wrong-segments", "wrong-speed", "undeclared-action", "invented-reach", "missing-growth", "distant-contact", "retained-length", "retained-diagnostic"):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)/"run"
                shutil.copytree(self.sample,root)
                changed = json.loads(json.dumps(events))
                if mutation == "invented-cell":
                    cell = next(e for e in reversed(changed) if e["event"] == "cell").copy()
                    cell["cell"] = 99999
                    changed.append(cell)
                elif mutation == "wrong-segments":
                    summary = json.loads((root/"summary.json").read_text())
                    summary["segments"] += 1
                    (root/"summary.json").write_text(json.dumps(summary))
                elif mutation == "wrong-speed":
                    next(e for e in changed if e["event"] == "growth")["speed"] = 1e12
                elif mutation == "undeclared-action":
                    changed.append({"event":"import-connectome", "step":150, "time":1.5})
                elif mutation == "invented-reach":
                    next(e for e in changed if e["event"] == "synapse")["reach"] = 1e9
                elif mutation == "missing-growth":
                    changed = [e for e in changed if e["event"] not in ("growth","branch","stop")]
                elif mutation == "distant-contact":
                    contact = next(e for e in changed if e["event"] == "synapse")
                    contact["pre_position"] = contact["post_position"] = [1e12,1e12,1e12]
                elif mutation == "retained-length":
                    next(e for e in changed if e["event"] == "connection")["pre_length"] += 1
                else:
                    summary = json.loads((root/"summary.json").read_text())
                    summary["retained_sites_outside_segments"] += 1
                    (root/"summary.json").write_text(json.dumps(summary))
                (root/"trace.jsonl").write_text("".join(json.dumps(e)+"\n" for e in changed))
                with self.assertRaises(ValueError):
                    self.model.verify_trace(root)

    def test_recognition_changes_growth_only_after_physical_contact(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.model.develop(self.model.default_program()._replace(recognition=0),self.classes,root/"off",duration=1.5)
            on = [json.loads(line) for line in (self.sample/"trace.jsonl").read_text().splitlines()]
            off = [json.loads(line) for line in (root/"off/trace.jsonl").read_text().splitlines()]
            first_contact = min(event["step"] for event in on if event["event"] == "synapse")
            before_contact = lambda events: [event for event in events if event["step"] <= first_contact and event["event"] != "synapse"]
            self.assertEqual(before_contact(on),before_contact(off))
            self.assertFalse(any(event["event"] == "synapse" for event in off))

    def test_backend_manifest_and_inventory_are_bound_to_build_result(self):
        for mutation in ("omitted", "coherent", "extra"):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                shutil.copytree(self.classes.classes.parent,root/"backend")
                backend = self.classes._replace(classes=root/"backend/classes")
                manifest = root/"backend/build.json"
                data = json.loads(manifest.read_text())
                if mutation == "omitted":
                    del data["classes"]["Construction.class"]
                    manifest.write_text(json.dumps(data))
                elif mutation == "coherent":
                    (backend.classes/"Construction.class").write_bytes(b"changed bytecode")
                    data["classes"]["Construction.class"] = hashlib.sha256(b"changed bytecode").hexdigest()
                    manifest.write_text(json.dumps(data))
                else:
                    (backend.classes/"Extra.class").write_bytes(b"extra")
                with self.assertRaisesRegex(ValueError,"backend"):
                    self.model.develop(self.model.default_program(),backend,root/"run")
                self.assertFalse((root/"run").exists())

    def test_contact_resource_budget_checks_counts_before_allocation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            probe = root/"BudgetProbe.java"
            probe.write_text('''public class BudgetProbe {
                public static void main(String[] args) {
                    assert Construction.protrusionCount(9,1,10) == 9;
                    assert Construction.protrusionCount(10,1,10) == 10;
                    for (double spacing : new double[]{10.0/11,1e-300,Double.MIN_VALUE}) {
                        try { Construction.protrusionCount(10,spacing,10); throw new AssertionError("overflow accepted"); }
                        catch (IllegalStateException expected) { }
                    }
                }
            }''')
            subprocess.run([shutil.which("javac"),"--release","11","-cp",str(self.classes.classes),str(probe)],check=True,capture_output=True)
            subprocess.run([shutil.which("java"),"-ea","-Djava.awt.headless=true","-cp",str(root)+":"+str(self.classes.classes),"BudgetProbe"],check=True,capture_output=True)
            with self.assertRaisesRegex(ValueError,"partial evidence"):
                self.model.develop(self.model.default_program()._replace(contact_spacing=1e-300),self.classes,root/"run",duration=1.5)
            self.assertIn("protrusion budget exhausted",(root/"run/runtime.log").read_text())
            self.assertTrue((root/"run/trace.jsonl").exists())
            self.assertFalse((root/"run/report.json").exists())

    def test_division_depends_on_accumulated_volume_not_future_growth_rate(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            program = self.model.default_program()._replace(volume_rate=0.0)
            result = self.model.develop(program,self.classes,root/"run",duration=1.5)
            self.assertEqual(result["cells"],2)
            self.assertEqual(result["segments"],0)

    def test_timed_guidance_has_identical_prefix_and_local_effect(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            program = self.model.default_program()._replace(division_signal=2.0,recognition=0,branch_rate=0)
            self.model.develop(program,self.classes,root/"control",duration=.3)
            self.model.develop(program,self.classes,root/"change",duration=.3,cue_switch=(15,-40.0))
            control = [json.loads(line) for line in (root/"control/trace.jsonl").read_text().splitlines()]
            changed = [json.loads(line) for line in (root/"change/trace.jsonl").read_text().splitlines()]
            self.assertEqual([e for e in control if e["step"] < 15],[e for e in changed if e["step"] < 15])
            a = next(e for e in control if e["step"] == 15 and e["event"] == "growth")
            b = next(e for e in changed if e["step"] == 15 and e["event"] == "growth")
            self.assertEqual(a["position"],b["position"])
            self.assertEqual(a["noise"],b["noise"])
            self.assertGreater(a["gradient"][2],0)
            self.assertLess(b["gradient"][2],0)
            self.assertNotEqual(a["direction"],b["direction"])

    def test_repeatability_unseen_programs_and_translation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.model.develop(self.model.default_program(),self.classes,root/"repeat",duration=1.5)
            self.assertEqual((self.sample/"trace.jsonl").read_bytes(),(root/"repeat/trace.jsonl").read_bytes())
            for signal,seed in ((2.0,11),(3.0,23)):
                program = self.model.default_program()._replace(division_signal=signal,recognition=0,branch_rate=0,axon_speed=65,dendrite_speed=35)
                result = self.model.develop(program,self.classes,root/f"unseen-{seed}",duration=1.5,seed=seed)
                self.assertEqual(result["cells"],2**math.ceil(math.log2(signal)))
                self.assertGreater(result["length"],0)
            program = self.model.default_program()._replace(division_signal=2,recognition=0,branch_rate=0)
            for name,origin in (("zero",(0,0,0)),("moved",(123,-67,49))):
                self.model.develop(program,self.classes,root/name,duration=.2,origin=origin)
            def final_positions(name):
                events = [json.loads(line) for line in (root/name/"trace.jsonl").read_text().splitlines()]
                return [e["position"] for e in events if e["step"] == 20 and e["event"] in ("cell","segment")]
            a,b = final_positions("zero"),final_positions("moved")
            self.assertEqual(len(a),len(b))
            for first,second in zip(a,b):
                self.assertLess(math.dist(first,[x-offset for x,offset in zip(second,(123,-67,49))]),1e-8)
            with self.assertRaisesRegex(ValueError,"partial evidence"):
                self.model.develop(program,self.classes,root/"limited",duration=.2,max_cells=1)
            self.assertIn("resource budget exhausted",(root/"limited/runtime.log").read_text())

    def test_public_cli_binds_recipe_runtime_activity_and_viewer(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)/"experiment"
            result = subprocess.run([sys.executable,str(EXAMPLE),"--output",str(root),"--duration","1.5"],capture_output=True,text=True,timeout=120)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            report = json.loads((root/"experiment.json").read_text())
            self.assertTrue(report["compiled_valid"])
            self.assertTrue(report["activity_controls_pass"])
            self.assertGreater(report["construction"]["synapses"],0)
            activity = json.loads((root/"activity.json").read_text())
            self.assertEqual(activity["implementation_sha256"],hashlib.sha256(EXAMPLE.with_name("activity.py").read_bytes()).hexdigest())
            for name,digest in report["evidence_sha256"].items():
                self.assertEqual(hashlib.sha256((root/name).read_bytes()).hexdigest(),digest)
            viewer = (root/"viewer.html").read_text()
            self.assertNotIn("/* CONSTRUCTION_DATA */ null",viewer)
            self.assertIn("biological_acceptance",viewer)

    def test_activity_oracle_respects_arrival_horizon(self):
        # graph fixture tests the assay oracle; the CLI test covers physical admission.
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            events = [{"event":"cell","cell":cell,"step":1} for cell in range(22)]
            events += [{"event":"connection","pre":cell,"post":cell+1,"step":1} for cell in range(21)]
            (root/"trace.jsonl").write_text("".join(json.dumps(event)+"\n" for event in events))
            report = self.activity().assay(root,lambda _: {"cells":22,"retained_sites_outside_segments":0})
            self.assertTrue(report["transmission_controls_pass"])
            self.assertEqual(report["trials"][0]["downstream_cells"],list(range(1,20)))

    def test_contact_formation_requires_both_sites_within_their_cylinders(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            probe = root/"ContactProbe.java"
            probe.write_text('''import ini.cx3d.cells.*;
                import ini.cx3d.localBiology.*;
                import ini.cx3d.physics.*;
                import ini.cx3d.simulations.*;
                import ini.cx3d.synapses.*;
                import static ini.cx3d.utilities.Matrix.*;
                public class ContactProbe {
                    static double[] site(PhysicalCylinder c, double fraction, double side) {
                        double[] proximal = subtract(c.getMassLocation(),c.getSpringAxis());
                        double[] point = add(proximal,scalarMult(fraction*c.getActualLength(),c.getUnitaryAxisDirectionVector()),
                            new double[]{0,side*c.getDiameter()/2,0});
                        return c.transformCoordinatesGlobalToPolar(point);
                    }
                    public static void main(String[] args) {
                        ECM world = ECM.getInstance(); ECM.setRandomSeed(7);
                        for (int i=0;i<18;i++) world.getPhysicalNodeInstance(randomNoise(1000,3));
                        NeuriteElement pre = CellFactory.getCellInstance(new double[]{0,0,0}).getSomaElement().extendNewNeurite(new double[]{1,0,0});
                        NeuriteElement post = CellFactory.getCellInstance(new double[]{0,3,0}).getSomaElement().extendNewNeurite(new double[]{1,0,0});
                        PhysicalCylinder a = pre.getPhysicalCylinder(), b = post.getPhysicalCylinder();
                        PhysicalBouton bouton = new PhysicalBouton(a,site(a,Double.parseDouble(args[0]),1),2);
                        PhysicalSpine spine = new PhysicalSpine(b,site(b,Double.parseDouble(args[1]),-1),2);
                        a.addExcrescence(bouton); b.addExcrescence(spine);
                        pre.synapseBetweenExistingBS(1);
                        boolean linked = bouton.getEx()==spine && spine.getEx()==bouton;
                        assert linked == Boolean.parseBoolean(args[2]) : "incorrect physical contact admission";
                    }
                }''')
            subprocess.run([shutil.which("javac"),"--release","11","-cp",str(self.classes.classes),str(probe)],check=True,capture_output=True)
            for side in (0,1):
                for fraction in (-.01,0,.999999999,1,1.04389785):
                    with self.subTest(side=side,fraction=fraction):
                        sites = [.5,.5]
                        sites[side] = fraction
                        result = subprocess.run([shutil.which("java"),"-ea","-Djava.awt.headless=true","-cp",
                            os.pathsep.join((str(root),str(self.classes.classes))),"ContactProbe",
                            *map(str,sites),str(0 <= fraction <= 1).lower()],capture_output=True,text=True,timeout=30)
                        self.assertEqual(result.returncode,0,result.stdout+result.stderr)

    def test_larger_populations_audit_formation_and_report_retained_geometry(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for signal,seed in ((8,7),(8,19),(16,7),(16,19)):
                with self.subTest(signal=signal,seed=seed):
                    program = self.model.default_program()._replace(division_signal=signal)
                    result = self.model.develop(program,self.classes,root/f"{signal}-{seed}",duration=3,seed=seed)
                    self.assertEqual(result["cells"],signal)
                    self.assertGreater(result["synapses"],0)
                    self.assertEqual(result["retained_sites_outside_segments"],0)
                    activity = self.activity().assay(root/f"{signal}-{seed}",self.model.verify_trace)
                    self.assertEqual(activity["retained_sites_outside_segments"],result["retained_sites_outside_segments"])
                    self.assertTrue(activity["transmission_controls_pass"])

    def test_remeshing_preserves_attachment_identity_position_and_bonds(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            probe = ROOT/"tests/data/NeuralAttachmentProbe.java"
            compiled = subprocess.run([shutil.which("javac"),"--release","11","-cp",str(self.classes.classes),"-d",str(root),str(probe)],capture_output=True,text=True)
            self.assertEqual(compiled.returncode,0,compiled.stdout+compiled.stderr)
            for mode in ("split","merge","proximal-only","distal-only","tapered","rotated",
                         "passive-compress","passive-stretch","passive-branch","passive-soma",
                         "passive-invalid-length","passive-collapsed-length","passive-child-failure","active-growth"):
                with self.subTest(mode=mode):
                    result = subprocess.run([shutil.which("java"),"-ea","-Djava.awt.headless=true","-cp",
                        os.pathsep.join((str(root),str(self.classes.classes))),"NeuralAttachmentProbe",mode],capture_output=True,text=True,timeout=30)
                    self.assertEqual(result.returncode,0,result.stdout+result.stderr)


if __name__ == "__main__":
    unittest.main()

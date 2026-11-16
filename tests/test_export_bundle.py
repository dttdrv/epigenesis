from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
import zipfile

from brainc.development_bundle import compile_development
from brainc.export import ExportError, export_omex, export_rocrate
from brainc.source import FASTA_PROFILE, compile_source
from brainc.variant import compile_variants


ROOT = Path(__file__).resolve().parents[1]


class ExportBundleTests(unittest.TestCase):
    def _build_fixture(self, root: Path) -> tuple[Path, Path, Path]:
        fasta = ROOT / "tests" / "data" / "J02482.1.fasta"
        interp = ROOT / "examples" / "minimal-development" / "interpretation"
        src_dir = root / "source"
        dev_dir = root / "development"
        vc_path = root / "variants.json"

        compile_source(FASTA_PROFILE, {"sequence": fasta}, parameters={"wrapper": "identity"}).save(src_dir)
        # Reuse minimal-development chain with a matching request bound to src_dir?
        # Use minimal-development's own sequence.fasta so the committed chain binds cleanly.
        seq_fasta = ROOT / "examples" / "minimal-development" / "sequence.fasta"
        src_dir_min = root / "source_min"
        compile_source(FASTA_PROFILE, {"sequence": seq_fasta}, parameters={"wrapper": "identity"}).save(src_dir_min)
        req_obj = json.loads((interp / "request.json").read_text(encoding="utf-8"))
        src_desc = json.loads((src_dir_min / "source.json").read_text(encoding="utf-8"))
        req_obj["source"] = {
            "format": "brainc.source-descriptor",
            "version": 1,
            "artifact_sha256": src_desc["artifact_sha256"],
            "ir_sha256": src_desc["source_ir_sha256"],
        }
        from brainc._canonical import digest
        req_core = {k: v for k, v in req_obj.items() if k != "artifact_sha256"}
        req_obj["artifact_sha256"] = digest(req_core)
        req_path = root / "request.json"
        req_path.write_text(json.dumps(req_obj, indent=2), encoding="utf-8")

        resp_obj = json.loads((interp / "response.json").read_text(encoding="utf-8"))
        resp_obj["request_artifact_sha256"] = req_obj["artifact_sha256"]
        resp_core = {k: v for k, v in resp_obj.items() if k != "artifact_sha256"}
        resp_obj["artifact_sha256"] = digest(resp_core)
        resp_path = root / "response.json"
        resp_path.write_text(json.dumps(resp_obj, indent=2), encoding="utf-8")

        # Update manifest accepts to allow fasta-dna/v1
        man_obj = json.loads((interp / "manifest.json").read_text(encoding="utf-8"))
        man_obj["accepts"] = ["brainc.source-descriptor/v1;profile=fasta-dna/v1"]
        man_core = {k: v for k, v in man_obj.items() if k != "artifact_sha256"}
        man_obj["artifact_sha256"] = digest(man_core)
        man_path = root / "manifest.json"
        man_path.write_text(json.dumps(man_obj, indent=2), encoding="utf-8")

        # Re-bind request to new manifest
        req_obj["provider_manifest_sha256"] = man_obj["artifact_sha256"]
        req_core = {k: v for k, v in req_obj.items() if k != "artifact_sha256"}
        req_obj["artifact_sha256"] = digest(req_core)
        req_path.write_text(json.dumps(req_obj, indent=2), encoding="utf-8")

        resp_obj["request_artifact_sha256"] = req_obj["artifact_sha256"]
        resp_core = {k: v for k, v in resp_obj.items() if k != "artifact_sha256"}
        resp_obj["artifact_sha256"] = digest(resp_core)
        resp_path.write_text(json.dumps(resp_obj, indent=2), encoding="utf-8")

        compile_development(
            src_dir_min,
            man_path,
            req_path,
            resp_path,
            interp / "lowering-policy.json",
            interp / "target-contract.json",
        ).save(dev_dir)
        compile_variants(
            src_dir_min,
            [{"id": "ex-snv-0", "record_id": "example.1", "start": 0, "end": 1, "ref": "A", "alt": "T"}],
        ).save(vc_path)
        return src_dir_min, dev_dir, vc_path

    def test_rocrate_and_omex_are_deterministic_and_verified(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            src_dir, dev_dir, vc_path = self._build_fixture(root)

            crate_a = root / "crate_a"
            crate_b = root / "crate_b"
            sum_a = export_rocrate(dev_dir, crate_a, source_bundle=src_dir, variant_collection=vc_path)
            sum_b = export_rocrate(dev_dir, crate_b, source_bundle=src_dir, variant_collection=vc_path)
            self.assertEqual(sum_a["rocrate_metadata_sha256"], sum_b["rocrate_metadata_sha256"])
            self.assertEqual(sum_a["files"], sum_b["files"])

            meta = json.loads((crate_a / "ro-crate-metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(meta["@context"], "https://w3id.org/ro/crate/1.1/context")
            by_id = {node["@id"]: node for node in meta["@graph"]}
            self.assertIn("ro-crate-metadata.json", by_id)
            self.assertIn("./", by_id)
            self.assertIn("development/bundle.json", by_id)
            self.assertIn("source/source.json", by_id)
            self.assertIn("variants.json", by_id)
            self.assertEqual(
                by_id["variants.json"]["sha256"],
                hashlib.sha256((crate_a / "variants.json").read_bytes()).hexdigest(),
            )

            omex_a = root / "archive_a.omex"
            omex_b = root / "archive_b.omex"
            res_a = export_omex(dev_dir, omex_a, source_bundle=src_dir, variant_collection=vc_path)
            res_b = export_omex(dev_dir, omex_b, source_bundle=src_dir, variant_collection=vc_path)
            self.assertEqual(res_a["archive_sha256"], res_b["archive_sha256"])
            self.assertEqual(omex_a.read_bytes(), omex_b.read_bytes())

            with zipfile.ZipFile(omex_a, "r") as zf:
                names = set(zf.namelist())
                self.assertIn("manifest.xml", names)
                self.assertIn("ro-crate-metadata.json", names)
                self.assertIn("development/bundle.json", names)
                self.assertIn("source/source.json", names)
                self.assertIn("variants.json", names)
                manifest_xml = ET.fromstring(zf.read("manifest.xml"))
                locations = {el.attrib["location"] for el in manifest_xml}
                self.assertIn(".", locations)
                self.assertIn("./manifest.xml", locations)
                self.assertIn("./variants.json", locations)

            # CLI round-trip
            cli_omex = root / "cli.omex"
            cli_res = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "brainc",
                    "export-omex",
                    str(dev_dir),
                    "--source-bundle",
                    str(src_dir),
                    "--variants",
                    str(vc_path),
                    "-o",
                    str(cli_omex),
                ],
                check=False,
                capture_output=True,
                text=True,
                cwd=ROOT,
            )
            self.assertEqual(cli_res.returncode, 0, cli_res.stderr)
            self.assertEqual(cli_omex.read_bytes(), omex_a.read_bytes())

    def test_rejects_mismatched_variant_bundle_and_symlink_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            src_dir, dev_dir, vc_path = self._build_fixture(root)
            other_src = root / "other_src"
            other_vc = root / "other_vc.json"
            other_fasta = root / "other.fasta"
            other_fasta.write_text(">other\nACGTACGT\n", encoding="utf-8")
            compile_source(FASTA_PROFILE, {"sequence": other_fasta}, parameters={"wrapper": "identity"}).save(other_src)
            compile_variants(
                other_src,
                [{"id": "o-0", "record_id": "other", "start": 0, "end": 1, "ref": "A", "alt": "G"}],
            ).save(other_vc)
            with self.assertRaises(ExportError):
                export_omex(dev_dir, root / "bad.omex", source_bundle=src_dir, variant_collection=other_vc)
            # Even when source_bundle=None, variant collection must bind bundle's source
            with self.assertRaises(ExportError):
                export_omex(dev_dir, root / "bad_nosrc.omex", variant_collection=other_vc)

            # Tampered development_module must be rejected by export_rocrate / export_omex
            tampered_dev = root / "tampered_dev"
            import shutil
            shutil.copytree(dev_dir, tampered_dev)
            mod_path = tampered_dev / "development_module.json"
            mod_obj = json.loads(mod_path.read_text(encoding="utf-8"))
            mod_obj["module"]["tensors"][0]["storage"]["sha256"] = "0" * 64
            mod_path.write_text(json.dumps(mod_obj), encoding="utf-8")
            with self.assertRaises(ExportError):
                export_rocrate(tampered_dev, root / "bad_crate")

            if hasattr(os, "symlink"):
                sym = root / "sym.omex"
                os.symlink(root / "target.omex", sym)
                with self.assertRaises(ExportError):
                    export_omex(dev_dir, sym, source_bundle=src_dir, variant_collection=vc_path)


if __name__ == "__main__":
    unittest.main()

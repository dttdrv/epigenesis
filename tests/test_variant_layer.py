from __future__ import annotations

import copy
import gzip
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from brainc._canonical import digest
from brainc.source import (
    FASTA_PROFILE,
    GENBANK_PROFILE,
    GFF3_PROFILE,
    RAW_PROFILE,
    REFERENCE_FASTA_PROFILE,
    compile_source,
)
from brainc.variant import (
    VariantError,
    apply_variant_to_sequence,
    compile_variants,
    load_variant_collection,
    revert_variant,
    vrs_allele,
    vrs_allele_digest,
    vrs_location_digest,
)
from brainc.validator_variant import (
    validate_variant_paths,
)


SCHEMAS_DIR = ROOT / "brainc" / "schemas"
CONFORMANCE_DIR = ROOT / "tests" / "conformance"


def _assert_matches_schema(instance: object, schema_doc: dict, root_schema: dict | None = None, path: str = "$") -> None:
    """Minimal Draft 2020-12 structural validator using Python stdlib only."""
    if root_schema is None:
        root_schema = schema_doc
    if "$ref" in schema_doc:
        ref = schema_doc["$ref"]
        if not ref.startswith("#/$defs/"):
            raise AssertionError(f"Unsupported $ref {ref!r} at {path}")
        target = root_schema["$defs"][ref[len("#/$defs/") :]]
        _assert_matches_schema(instance, target, root_schema, path)
        return
    if "anyOf" in schema_doc:
        for branch in schema_doc["anyOf"]:
            try:
                _assert_matches_schema(instance, branch, root_schema, path)
                return
            except AssertionError:
                continue
        raise AssertionError(f"{path}: value {instance!r} matched no anyOf branch")
    if "const" in schema_doc:
        expected = schema_doc["const"]
        if type(instance) is not type(expected) or instance != expected:
            raise AssertionError(f"{path}: expected const {expected!r}, got {instance!r}")
    if "enum" in schema_doc:
        if instance not in schema_doc["enum"]:
            raise AssertionError(f"{path}: {instance!r} not in enum {schema_doc['enum']!r}")
    expected_type = schema_doc.get("type")
    if expected_type == "object":
        if type(instance) is not dict:
            raise AssertionError(f"{path}: expected object, got {type(instance).__name__}")
        req = schema_doc.get("required", [])
        for key in req:
            if key not in instance:
                raise AssertionError(f"{path}: missing required key {key!r}")
        if "minProperties" in schema_doc and len(instance) < schema_doc["minProperties"]:
            raise AssertionError(f"{path}: too few properties")
        if "maxProperties" in schema_doc and len(instance) > schema_doc["maxProperties"]:
            raise AssertionError(f"{path}: too many properties")
        props = schema_doc.get("properties", {})
        addl = schema_doc.get("additionalProperties", True)
        for key, val in instance.items():
            if key in props:
                _assert_matches_schema(val, props[key], root_schema, f"{path}.{key}")
            elif addl is False:
                raise AssertionError(f"{path}: unknown key {key!r}")
            elif isinstance(addl, dict):
                _assert_matches_schema(val, addl, root_schema, f"{path}.{key}")
    elif expected_type == "array":
        if type(instance) is not list:
            raise AssertionError(f"{path}: expected array, got {type(instance).__name__}")
        if "minItems" in schema_doc and len(instance) < schema_doc["minItems"]:
            raise AssertionError(f"{path}: fewer items than minItems")
        if "maxItems" in schema_doc and len(instance) > schema_doc["maxItems"]:
            raise AssertionError(f"{path}: more items than maxItems")
        if "items" in schema_doc:
            for idx, item in enumerate(instance):
                _assert_matches_schema(item, schema_doc["items"], root_schema, f"{path}[{idx}]")
    elif expected_type == "string":
        if type(instance) is not str:
            raise AssertionError(f"{path}: expected string, got {type(instance).__name__}")
        if "minLength" in schema_doc and len(instance) < schema_doc["minLength"]:
            raise AssertionError(f"{path}: shorter than minLength")
        if "maxLength" in schema_doc and len(instance) > schema_doc["maxLength"]:
            raise AssertionError(f"{path}: longer than maxLength")
        if "pattern" in schema_doc and re.fullmatch(schema_doc["pattern"], instance) is None:
            raise AssertionError(f"{path}: {instance!r} does not match pattern {schema_doc['pattern']!r}")
    elif expected_type == "integer":
        if type(instance) is not int:
            raise AssertionError(f"{path}: expected integer, got {type(instance).__name__}")
        if "minimum" in schema_doc and instance < schema_doc["minimum"]:
            raise AssertionError(f"{path}: {instance} < minimum {schema_doc['minimum']}")
        if "maximum" in schema_doc and instance > schema_doc["maximum"]:
            raise AssertionError(f"{path}: {instance} > maximum {schema_doc['maximum']}")
    elif expected_type == "boolean":
        if type(instance) is not bool:
            raise AssertionError(f"{path}: expected boolean, got {type(instance).__name__}")
    elif expected_type == "null":
        if instance is not None:
            raise AssertionError(f"{path}: expected null, got {instance!r}")


class VariantLayerTests(unittest.TestCase):
    def test_ga4gh_vrs_v2_conformance_and_pinned_j02482_vectors(self) -> None:
        vectors = json.loads((CONFORMANCE_DIR / "vrs-v2-vectors.json").read_text(encoding="utf-8"))
        apoe = vectors["ga4gh_apoe_vector"]
        loc_dig = vrs_location_digest(apoe["refgetAccession"], apoe["start"], apoe["end"])
        self.assertEqual(loc_dig, apoe["expected_location_digest"])
        allele_dig = vrs_allele_digest(loc_dig, apoe["state_sequence"])
        self.assertEqual(allele_dig, apoe["expected_allele_digest"])
        allele = vrs_allele(apoe["refgetAccession"], apoe["start"], apoe["end"], apoe["state_sequence"])
        self.assertEqual(allele["location"]["id"], apoe["expected_location_id"])
        self.assertEqual(allele["id"], apoe["expected_allele_id"])

        jvec = vectors["j02482_variant_collection"]
        bundle = compile_source(
            jvec["source_profile"],
            {"sequence": ROOT / jvec["source_path"]},
            parameters=jvec["parameters"],
        )
        vc = compile_variants(bundle, jvec["inputs"])
        payload = vc.to_dict()
        self.assertEqual(payload["variants_ir_sha256"], jvec["expected_variants_ir_sha256"])
        self.assertEqual(payload["artifact_sha256"], jvec["expected_artifact_sha256"])
        for expected_entry, actual_entry in zip(jvec["expected_entries"], payload["variants_ir"]["variants"]):
            self.assertEqual(actual_entry["id"], expected_entry["id"])
            self.assertEqual(actual_entry["edit_kind"], expected_entry["edit_kind"])
            self.assertEqual(actual_entry["vrs_location_id"], expected_entry["vrs_location_id"])
            self.assertEqual(actual_entry["vrs_allele_id"], expected_entry["vrs_allele_id"])

        vc_schema = json.loads((SCHEMAS_DIR / "variant-collection.v1.schema.json").read_text(encoding="utf-8"))
        sd_schema = json.loads((SCHEMAS_DIR / "source-descriptor.v1.schema.json").read_text(encoding="utf-8"))
        _assert_matches_schema(bundle.to_dict(), sd_schema)
        _assert_matches_schema(payload, vc_schema)

    def test_negative_conformance_vectors_rejected(self) -> None:
        neg = json.loads((CONFORMANCE_DIR / "variant-negative-vectors.json").read_text(encoding="utf-8"))
        bundle = compile_source(
            FASTA_PROFILE,
            {"sequence": ROOT / "tests" / "data" / "J02482.1.fasta"},
            parameters={"wrapper": "identity"},
        )
        for case in neg["invalid_compile_inputs"]:
            with self.subTest(case=case["name"]):
                with self.assertRaises(VariantError) as ctx:
                    compile_variants(bundle, [case["variant"]])
                self.assertIn(case["expected_code"], str(ctx.exception))

        for bad_edit in (
            "not-a-dict",
            {"start": 0, "end": 1, "ref": "A"},
            {"location": "bad", "ref": "A", "alt": "T"},
            {"start": 0, "end": 0, "ref": "A", "alt": "T"},
        ):
            with self.assertRaises(VariantError):
                apply_variant_to_sequence("ACGT", bad_edit)  # type: ignore[arg-type]
            with self.assertRaises(VariantError):
                revert_variant("ACGT", bad_edit)  # type: ignore[arg-type]

    def test_compile_and_validate_across_source_profiles_and_reversion(self) -> None:
        vc_schema = json.loads((SCHEMAS_DIR / "variant-collection.v1.schema.json").read_text(encoding="utf-8"))
        rep_schema = json.loads((SCHEMAS_DIR / "variant-validation-report.v1.schema.json").read_text(encoding="utf-8"))

        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)

            # 1. RAW IUPAC profile
            raw_path = tmp_dir / "raw.dna"
            raw_seq = "ACGTACGTACGT"
            raw_path.write_text(raw_seq, encoding="ascii")
            raw_bundle = compile_source(RAW_PROFILE, {"sequence": raw_path}, parameters={"record_id": "rec1"})
            raw_dir = tmp_dir / "src-raw"
            raw_bundle.save(raw_dir)
            edits = [
                {"id": "v-snv", "record_id": "rec1", "start": 1, "end": 2, "ref": "C", "alt": "T"},
                {"id": "v-mnv", "record_id": "rec1", "start": 4, "end": 6, "ref": "AC", "alt": "TG"},
                {"id": "v-ins", "record_id": "rec1", "start": 8, "end": 8, "ref": "", "alt": "GGG"},
                {"id": "v-del", "record_id": "rec1", "start": 8, "end": 11, "ref": "ACG", "alt": ""},
                {"id": "v-delins", "record_id": "rec1", "start": 0, "end": 3, "ref": "ACG", "alt": "T"},
            ]
            vc = compile_variants(raw_bundle, edits)
            vc_path = tmp_dir / "raw-variants.json"
            vc.save(vc_path)
            _assert_matches_schema(vc.to_dict(), vc_schema)
            loaded = load_variant_collection(vc_path)
            self.assertEqual(loaded.artifact_sha256, vc.artifact_sha256)

            for entry in vc.to_dict()["variants_ir"]["variants"]:
                mutated = apply_variant_to_sequence(raw_seq, entry)
                reverted = revert_variant(mutated, entry)
                self.assertEqual(reverted, raw_seq)

            report = validate_variant_paths(raw_dir, vc_path)
            self.assertTrue(report["valid"], report)
            _assert_matches_schema(report, rep_schema)

            # 2. fasta-reference-dna/v1 with gzip wrapper
            gz_fasta = tmp_dir / "ref.fa.gz"
            gz_fasta.write_bytes(gzip.compress(b">contigA\nGGGGAAAATTTTCCCC\n"))
            ref_bundle = compile_source(
                REFERENCE_FASTA_PROFILE,
                {"sequence": gz_fasta},
                parameters={"wrapper": "gzip"},
            )
            ref_dir = tmp_dir / "src-ref"
            ref_bundle.save(ref_dir)
            ref_vc = compile_variants(
                ref_dir,
                [{"id": "cA-1", "record_id": "contigA", "start": 4, "end": 8, "ref": "AAAA", "alt": "TTTT"}],
                source_inputs={"sequence": gz_fasta},
            )
            ref_vc_path = tmp_dir / "ref-variants.json"
            ref_vc.save(ref_vc_path)
            ref_rep_ok = validate_variant_paths(ref_dir, ref_vc_path, source_inputs={"sequence": gz_fasta})
            self.assertTrue(ref_rep_ok["valid"], ref_rep_ok)
            ref_rep_missing = validate_variant_paths(ref_dir, ref_vc_path)
            self.assertFalse(ref_rep_missing["valid"])
            self.assertEqual(ref_rep_missing["code"], "VARVAL001")

            # 3. GenBank profile (U49845.1)
            gb_bundle = compile_source(
                GENBANK_PROFILE,
                {"genbank": ROOT / "tests" / "data" / "U49845.1.gb"},
                parameters={},
            )
            gb_dir = tmp_dir / "src-gb"
            gb_bundle.save(gb_dir)
            gb_vc = compile_variants(
                gb_dir,
                [{"id": "u49845-0", "record_id": "U49845.1", "start": 0, "end": 1, "ref": "G", "alt": "A"}],
            )
            gb_vc_path = tmp_dir / "gb-variants.json"
            gb_vc.save(gb_vc_path)
            gb_rep = validate_variant_paths(gb_dir, gb_vc_path)
            self.assertTrue(gb_rep["valid"], gb_rep)

            # 4. GFF3 profile (J02482.1)
            gff_bundle = compile_source(
                GFF3_PROFILE,
                {
                    "sequence": ROOT / "tests" / "data" / "J02482.1.fasta",
                    "annotation": ROOT / "tests" / "data" / "J02482.1.gff3",
                },
                parameters={"sequence_profile": FASTA_PROFILE, "sequence_parameters": {"wrapper": "identity"}},
            )
            gff_dir = tmp_dir / "src-gff"
            gff_bundle.save(gff_dir)
            gff_vc = compile_variants(
                gff_dir,
                [{"id": "gff-0", "record_id": "J02482.1", "start": 0, "end": 1, "ref": "G", "alt": "C"}],
            )
            gff_vc_path = tmp_dir / "gff-variants.json"
            gff_vc.save(gff_vc_path)
            gff_rep = validate_variant_paths(gff_dir, gff_vc_path)
            self.assertTrue(gff_rep["valid"], gff_rep)

    def test_tamper_type_confusion_and_link_safety_attacks(self) -> None:
        rep_schema = json.loads((SCHEMAS_DIR / "variant-validation-report.v1.schema.json").read_text(encoding="utf-8"))
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            seq_path = tmp_dir / "seq.dna"
            seq_path.write_text("ACGTACGT", encoding="ascii")
            src_bundle = compile_source(RAW_PROFILE, {"sequence": seq_path}, parameters={"record_id": "r1"})
            src_dir = tmp_dir / "src"
            src_bundle.save(src_dir)
            vc = compile_variants(src_bundle, [{"id": "e1", "record_id": "r1", "start": 1, "end": 2, "ref": "C", "alt": "T"}])
            vc_path = tmp_dir / "vc.json"
            vc.save(vc_path)

            # Tamper vrs_allele_id and recompute all hashes -> independent validator must reject
            tampered_vrs = vc.to_dict()
            tampered_vrs["variants_ir"]["variants"][0]["vrs_allele_id"] = "ga4gh:VA." + ("A" * 32)
            tampered_vrs["variants_ir_sha256"] = digest(tampered_vrs["variants_ir"])
            tampered_vrs["artifact_sha256"] = digest({k: v for k, v in tampered_vrs.items() if k != "artifact_sha256"})
            tampered_path = tmp_dir / "tampered-vrs.json"
            tampered_path.write_text(json.dumps(tampered_vrs), encoding="utf-8")
            with self.assertRaises(VariantError):
                load_variant_collection(tampered_path)
            rep1 = validate_variant_paths(src_dir, tampered_path)
            self.assertFalse(rep1["valid"])
            self.assertEqual(rep1["code"], "VARVAL001")
            _assert_matches_schema(rep1, rep_schema)

            # JSON type confusion (True in place of start=1) even with matching hashes
            bool_tampered = copy.deepcopy(vc.to_dict())
            bool_tampered["variants_ir"]["variants"][0]["vrs_allele"]["location"]["start"] = True
            bool_tampered["variants_ir_sha256"] = digest(bool_tampered["variants_ir"])
            bool_tampered["artifact_sha256"] = digest({k: v for k, v in bool_tampered.items() if k != "artifact_sha256"})
            bool_path = tmp_dir / "bool-tampered.json"
            bool_path.write_text(json.dumps(bool_tampered), encoding="utf-8")
            with self.assertRaises(VariantError):
                load_variant_collection(bool_path)
            rep2 = validate_variant_paths(src_dir, bool_path)
            self.assertFalse(rep2["valid"])
            self.assertEqual(rep2["code"], "VARVAL001")

            # Symlink attack on variant collection file
            symlink_path = tmp_dir / "vc-symlink.json"
            symlink_path.symlink_to(vc_path)
            rep3 = validate_variant_paths(src_dir, symlink_path)
            self.assertFalse(rep3["valid"])
            self.assertEqual(rep3["code"], "VARVAL001")

            # Hardlink attack on posix
            if os.name == "posix":
                hardlink_path = tmp_dir / "vc-hardlink.json"
                os.link(vc_path, hardlink_path)
                try:
                    rep4 = validate_variant_paths(src_dir, hardlink_path)
                    self.assertFalse(rep4["valid"])
                    self.assertEqual(rep4["code"], "VARVAL001")
                finally:
                    hardlink_path.unlink(missing_ok=True)

    def test_cli_and_isolated_standalone_validator(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            src_dir = tmp_dir / "src"
            inputs_path = tmp_dir / "variants-input.json"
            vc_path = tmp_dir / "variants.json"
            report_path = tmp_dir / "report.json"

            subprocess.run(
                [
                    sys.executable, "-m", "brainc", "compile-source",
                    "--profile", "fasta-dna/v1",
                    "--sequence", str(ROOT / "tests" / "data" / "J02482.1.fasta"),
                    "--wrapper", "identity",
                    "--output", str(src_dir),
                ],
                cwd=ROOT, check=True, capture_output=True, text=True,
            )
            inputs_path.write_text(
                json.dumps([{"id": "v1", "record_id": "J02482.1", "start": 0, "end": 1, "ref": "G", "alt": "T"}]),
                encoding="utf-8",
            )
            subprocess.run(
                [
                    sys.executable, "-m", "brainc", "compile-variants",
                    str(src_dir), "--variants", str(inputs_path), "--output", str(vc_path),
                ],
                cwd=ROOT, check=True, capture_output=True, text=True,
            )
            val_proc = subprocess.run(
                [
                    sys.executable, "-m", "brainc", "validate-variants",
                    str(src_dir), str(vc_path), "--output", str(report_path),
                ],
                cwd=ROOT, check=True, capture_output=True, text=True,
            )
            self.assertEqual(val_proc.returncode, 0)
            standalone_proc = subprocess.run(
                [
                    sys.executable, "-I", str(ROOT / "brainc" / "validator_variant.py"),
                    str(src_dir), str(vc_path),
                ],
                cwd=tmp_dir, check=True, capture_output=True, text=True,
            )
            standalone_rep = json.loads(standalone_proc.stdout)
            saved_rep = json.loads(report_path.read_text(encoding="utf-8"))
            self.assertEqual(standalone_rep, saved_rep)
            self.assertTrue(standalone_rep["valid"])

    def test_lowercase_dna_and_surrogate_rejection(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            low_fa = tmp_dir / "lower.fasta"
            low_fa.write_text(">low1\nacgtACGTacgt\n", encoding="utf-8")
            src_dir = tmp_dir / "src_low"
            compile_source(FASTA_PROFILE, {"sequence": low_fa}, parameters={"wrapper": "identity"}).save(src_dir)
            vc = compile_variants(
                src_dir,
                [{"id": "low-snv", "record_id": "low1", "start": 0, "end": 1, "ref": "A", "alt": "T"}],
            )
            vc_path = tmp_dir / "low-vc.json"
            vc.save(vc_path)
            rep = validate_variant_paths(src_dir, vc_path)
            self.assertTrue(rep["valid"], rep)

            # Lone Unicode surrogate inside variant JSON fails cleanly with sealed VARVAL001
            bad_json = tmp_dir / "surrogate.json"
            bad_json.write_bytes(b'{"id":"\\ud800"}')
            bad_rep = validate_variant_paths(src_dir, bad_json)
            self.assertFalse(bad_rep["valid"])
            self.assertEqual(bad_rep["code"], "VARVAL001")


if __name__ == "__main__":
    unittest.main()

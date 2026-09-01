from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

from brainc._canonical import digest
from brainc.external_profile import (
    seal_profile_manifest,
    seal_source_descriptor,
    seal_validation_report,
)
from brainc.external_source import (
    ExternalSourceError,
    build_external_source_closure,
    build_external_source_closure_from_paths,
    input_references,
    profile_parameters,
    source_records,
    validate_external_source_closure,
)
from brainc.source import (
    EXTERNAL_PROFILE,
    SourceError,
    compile_external_source,
    compile_external_source_paths,
    load_source_bundle,
    validate_source_bundle,
)
from tests.test_external_profile import (
    FASTQ,
    NATIVE_INDEX,
    _closure,
    _profile_ir,
    _sha,
)


def _build_closure() -> tuple[dict, dict, dict, dict]:
    manifest, descriptor, report, _, _, _ = _closure()
    closure = build_external_source_closure(
        manifest,
        descriptor,
        report,
        original_payloads={"reads": FASTQ},
        native_artifact_payloads={"read-index": NATIVE_INDEX},
    )
    return manifest, descriptor, report, closure


def _reseal_source(value: dict) -> dict:
    value["source_ir_sha256"] = digest(value["source_ir"])
    value["artifact_sha256"] = digest(
        {key: member for key, member in value.items() if key != "artifact_sha256"}
    )
    return value


def _reseal_closure(value: dict) -> dict:
    value["closure_ir_sha256"] = digest(value["closure_ir"])
    value["artifact_sha256"] = digest(
        {key: member for key, member in value.items() if key != "artifact_sha256"}
    )
    return value


class ExternalSourceTests(unittest.TestCase):
    def test_path_admission_streams_originals_and_rejects_path_attacks(self) -> None:
        manifest, descriptor, report, closure = _build_closure()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "reads.fastq"
            native = root / "index.json"
            source.write_bytes(FASTQ)
            native.write_bytes(NATIVE_INDEX)
            observed = build_external_source_closure_from_paths(
                manifest,
                descriptor,
                report,
                original_paths={"reads": source},
                native_artifact_paths={"read-index": native},
            )
            self.assertEqual(observed, closure)
            bundle = compile_external_source_paths(
                manifest,
                descriptor,
                report,
                original_paths={"reads": source},
                native_artifact_paths={"read-index": native},
            )
            self.assertEqual(bundle.artifacts, {"external": closure})

            source.write_bytes(FASTQ + b"\n")
            with self.assertRaisesRegex(SourceError, "exact reference"):
                compile_external_source_paths(
                    manifest,
                    descriptor,
                    report,
                    original_paths={"reads": source},
                    native_artifact_paths={"read-index": native},
                )
            source.write_bytes(FASTQ)
            with self.assertRaises(SourceError):
                compile_external_source_paths(
                    manifest,
                    descriptor,
                    report,
                    original_paths={"wrong": source},
                    native_artifact_paths={"read-index": native},
                )
            if hasattr(Path, "symlink_to"):
                alias = root / "alias.fastq"
                alias.symlink_to(source)
                with self.assertRaises(SourceError):
                    compile_external_source_paths(
                        manifest,
                        descriptor,
                        report,
                        original_paths={"reads": alias},
                        native_artifact_paths={"read-index": native},
                    )

    def test_build_exposes_exact_profile_inputs_and_records(self) -> None:
        manifest, _, _, closure = _build_closure()
        self.assertEqual(validate_external_source_closure(closure), closure)
        self.assertEqual(
            profile_parameters(closure),
            {
                "profile_id": "org.example.fastq-sanger",
                "profile_version": 1,
                "profile_manifest_sha256": manifest["artifact_sha256"],
            },
        )
        self.assertEqual(
            input_references(closure),
            {
                "reads": {
                    "sha256": _sha(FASTQ),
                    "byte_length": len(FASTQ),
                }
            },
        )
        self.assertEqual(
            [(record["record_id"], record["bases"]) for record in source_records(closure)],
            [("read-1", 7), ("read-2", 5)],
        )

        bundle = compile_external_source(
            manifest,
            closure["closure_ir"]["source_descriptor"],
            closure["closure_ir"]["validation_report"],
            original_payloads={"reads": FASTQ},
            native_artifact_payloads={"read-index": NATIVE_INDEX},
        )
        source = bundle.to_dict()["source_ir"]
        self.assertEqual(source["profile"], EXTERNAL_PROFILE)
        self.assertEqual(source["parameters"], profile_parameters(closure))
        self.assertEqual(source["inputs"], input_references(closure))
        self.assertEqual(bundle.artifacts, {"external": closure})

    def test_bundle_save_load_round_trip_embeds_no_original_bytes(self) -> None:
        manifest, descriptor, report, closure = _build_closure()
        bundle = compile_external_source(
            manifest,
            descriptor,
            report,
            original_payloads={"reads": FASTQ},
            native_artifact_payloads={"read-index": NATIVE_INDEX},
        )
        wire = json.dumps(
            {"source": bundle.to_dict(), "artifacts": bundle.artifacts},
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=True,
        ).encode("utf-8")
        self.assertNotIn(FASTQ, wire)
        self.assertNotIn(NATIVE_INDEX, wire)
        self.assertNotIn(b"GATTACA", wire)
        self.assertNotIn(b"ACGTN", wire)
        self.assertEqual(
            closure["closure_ir"]["source_descriptor"]["frontend_ir"]["inputs"][0],
            {"role": "reads", "sha256": _sha(FASTQ), "byte_length": len(FASTQ)},
        )

        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / "external-source"
            paths = bundle.save(directory)
            self.assertEqual(set(paths), {"source.json", "external.json"})
            loaded = load_source_bundle(directory)
            validated = validate_source_bundle(loaded)
            self.assertEqual(validated.to_dict(), bundle.to_dict())
            self.assertEqual(validated.artifacts, bundle.artifacts)
            self.assertNotIn(FASTQ, (directory / "source.json").read_bytes())
            self.assertNotIn(FASTQ, (directory / "external.json").read_bytes())

    def test_coherent_closure_and_saved_profile_parameter_tampering_reject(self) -> None:
        manifest, descriptor, report, closure = _build_closure()
        attacked = deepcopy(closure)
        attacked["closure_ir"]["profile_manifest"]["profile_ir"]["profile"][
            "id"
        ] = "org.example.other-profile"
        with self.assertRaises(ExternalSourceError):
            validate_external_source_closure(_reseal_closure(attacked))

        bundle = compile_external_source(
            manifest,
            descriptor,
            report,
            original_payloads={"reads": FASTQ},
            native_artifact_payloads={"read-index": NATIVE_INDEX},
        )
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / "external-source"
            bundle.save(directory)
            source_path = directory / "source.json"
            source = json.loads(source_path.read_text(encoding="utf-8"))
            source["source_ir"]["parameters"]["profile_manifest_sha256"] = "0" * 64
            source_path.write_text(
                json.dumps(_reseal_source(source), ensure_ascii=False),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(SourceError, "profile does not match"):
                load_source_bundle(directory)

    def test_same_profile_name_and_version_cannot_substitute_manifest(self) -> None:
        manifest, descriptor, report, inputs, native, records = _closure()
        other_profile = _profile_ir()
        other_profile["grammar"]["authority"]["sha256"] = _sha(
            b"different pinned grammar authority"
        )
        other_manifest = seal_profile_manifest(other_profile)
        self.assertEqual(
            other_manifest["profile_ir"]["profile"],
            manifest["profile_ir"]["profile"],
        )
        self.assertNotEqual(
            other_manifest["artifact_sha256"], manifest["artifact_sha256"]
        )
        with self.assertRaises(ExternalSourceError):
            build_external_source_closure(
                other_manifest,
                descriptor,
                report,
                original_payloads={"reads": FASTQ},
                native_artifact_payloads={"read-index": NATIVE_INDEX},
            )

        other_descriptor = seal_source_descriptor(
            other_manifest,
            input_references=inputs,
            native_artifact_references=native,
            records=records,
        )
        other_report = seal_validation_report(
            other_manifest,
            other_descriptor,
            replayed_input_references=inputs,
            replayed_native_artifact_references=native,
            replayed_records=records,
        )
        with self.assertRaises(SourceError):
            compile_external_source(
                manifest,
                other_descriptor,
                other_report,
                original_payloads={"reads": FASTQ},
                native_artifact_payloads={"read-index": NATIVE_INDEX},
            )


if __name__ == "__main__":
    unittest.main()

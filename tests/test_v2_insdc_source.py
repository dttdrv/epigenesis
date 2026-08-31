from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

from brainc._canonical import digest
from brainc.insdc import FORMAT as GENBANK_FORMAT, VERSION as GENBANK_VERSION, GenBankCompiler
from brainc.v2 import (
    V2Error,
    compile_module,
    inline_storage,
    make_request,
    pack,
    policy_artifact,
    save,
    target_artifact,
)
from brainc.v2.limits import MAX_JSON_BYTES, MAX_SOURCE_JSON_BYTES
from brainc.v2.provider import SOURCE_FORMATS, load_manifest, load_source
from brainc.v2.target import DEV_DOMAIN, DEV_VERSION
from tests.test_v2 import _Chain


ROOT = Path(__file__).resolve().parents[1]
GENBANK = ROOT / "tests/data/U49845.1.gb"
RECORD_ID = "U49845.1"
RECORD_BASES = 5_028


def _seal(value: dict) -> dict:
    core = {key: item for key, item in value.items() if key != "artifact_sha256"}
    return {**core, "artifact_sha256": digest(core)}


def _write_padded_json(path: Path, value: dict, byte_length: int) -> None:
    raw = json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    if len(raw) > byte_length:
        raise AssertionError("test artifact exceeds requested boundary")
    padding = byte_length - len(raw)
    chunk = b" " * (1024 * 1024)
    with path.open("wb") as stream:
        while padding:
            written = min(padding, len(chunk))
            stream.write(chunk[:written])
            padding -= written
        stream.write(raw)


def _target() -> dict:
    return target_artifact(
        {
            "id": "org.example.genbank-source",
            "abi_major": 1,
            "opsets": [{"domain": DEV_DOMAIN, "version": DEV_VERSION}],
            "unit_schemas": [
                {
                    "id": "node",
                    "fields": [
                        {
                            "id": "value",
                            "type": {"dtype": "f64", "shape": []},
                            "unit": None,
                            "mutability": "state",
                            "numeric": {"kind": "float"},
                        }
                    ],
                }
            ],
            "edge_schemas": [],
            "ports": [],
            "rules": [],
        }
    )


class _GenBankChain:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.counter = 0
        self.source_path = root / "source.json"
        GenBankCompiler().compile_file(GENBANK).save(self.source_path)
        self.source, self.records = load_source(self.source_path)

        self.target = _target()
        self.target_path = self.write("target", self.target)
        axis = {
            "kind": "source-coordinate",
            "source_artifact_sha256": self.source["artifact_sha256"],
            "record_id": RECORD_ID,
            "start": RECORD_BASES - 1,
            "step": 1,
            "span": 1,
            "coordinate_system": "0-based-half-open",
        }
        self.output = {
            "id": "source.value",
            "type": {"dtype": "f64", "shape": [1]},
            "unit": None,
            "axes": [axis],
            "storage": inline_storage(pack("f64", [1.0])),
        }
        self.manifest = _seal(
            {
                "format": "brainc.provider-manifest",
                "version": 2,
                "provider": {"name": "org.example.provider", "version": "1"},
                "model_identity": {
                    "kind": "content-sha256",
                    "value": digest({"contract": "genbank-source-test", "version": 1}),
                },
                "accepts": [f"{GENBANK_FORMAT}/v{GENBANK_VERSION}"],
                "outputs": [
                    {
                        key: deepcopy(self.output[key])
                        for key in ("id", "type", "unit", "axes")
                    }
                ],
            }
        )
        self.manifest_path = self.write("manifest", self.manifest)
        self.request = make_request(
            self.source_path,
            self.manifest_path,
            [self.output["id"]],
        )
        self.request_path = self.write("request", self.request)
        self.response = _seal(
            {
                "format": "brainc.prediction-response",
                "version": 2,
                "request_artifact_sha256": self.request["artifact_sha256"],
                "provider": deepcopy(self.manifest["provider"]),
                "model_identity": deepcopy(self.manifest["model_identity"]),
                "outputs": [deepcopy(self.output)],
            }
        )
        self.response_path = self.write("response", self.response)
        target_binding = {
            "id": self.target["contract"]["id"],
            "abi_major": self.target["contract"]["abi_major"],
            "contract_sha256": self.target["contract_sha256"],
        }
        self.policy = policy_artifact(
            "org.example.genbank-source",
            target_binding,
            [{"id": self.output["id"], "from_output": self.output["id"]}],
            [],
        )
        self.policy_path = self.write("policy", self.policy)

    def write(self, label: str, value: dict) -> Path:
        self.counter += 1
        path = self.root / f"{self.counter:02d}-{label}.json"
        save(value, path)
        return path

    def compile(
        self,
        *,
        source_path: Path | None = None,
        request_path: Path | None = None,
    ) -> dict:
        return compile_module(
            source_path or self.source_path,
            self.manifest_path,
            request_path or self.request_path,
            self.response_path,
            self.policy_path,
            self.target_path,
        )


class GenBankSourceABITests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.chain = _GenBankChain(self.root)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_u49845_flows_through_request_and_module_with_exact_binding(self) -> None:
        self.assertEqual(
            SOURCE_FORMATS[(GENBANK_FORMAT, GENBANK_VERSION)],
            "bio_ir_sha256",
        )
        self.assertEqual(self.chain.records, {RECORD_ID: RECORD_BASES})
        expected_binding = {
            "format": GENBANK_FORMAT,
            "version": GENBANK_VERSION,
            "artifact_sha256": self.chain.source["artifact_sha256"],
            "ir_sha256": self.chain.source["bio_ir_sha256"],
        }
        self.assertEqual(self.chain.request["source"], expected_binding)
        first = self.chain.compile()
        second = self.chain.compile()
        self.assertEqual(first, second)
        self.assertEqual(first["sources"]["sequence"], expected_binding)

    def test_source_record_length_drives_coordinate_boundary(self) -> None:
        manifest = deepcopy(self.chain.manifest)
        manifest["outputs"][0]["axes"][0]["start"] = RECORD_BASES
        manifest_path = self.chain.write("past-record-end", _seal(manifest))
        with self.assertRaisesRegex(V2Error, "outside its sequence record"):
            make_request(
                self.chain.source_path,
                manifest_path,
                [self.chain.output["id"]],
            )

    def test_resealed_bio_ir_tamper_fails_genbank_replay(self) -> None:
        forged = deepcopy(self.chain.source)
        record = forged["bio_ir"]["records"][0]
        record["definition"] += " tampered"
        record["record_ir_sha256"] = digest(
            {key: value for key, value in record.items() if key != "record_ir_sha256"}
        )
        forged["bio_ir_sha256"] = digest(forged["bio_ir"])
        forged_path = self.chain.write("forged-source", _seal(forged))
        with self.assertRaisesRegex(V2Error, "does not match GenBank frontend replay"):
            load_source(forged_path)

    def test_resealed_request_cannot_substitute_the_bio_ir_binding(self) -> None:
        forged = deepcopy(self.chain.request)
        forged["source"]["ir_sha256"] = "0" * 64
        forged_path = self.chain.write("forged-request", _seal(forged))
        with self.assertRaisesRegex(V2Error, "not bound to the supplied sequence source"):
            self.chain.compile(request_path=forged_path)

    def test_source_accepts_exact_64_mib_and_rejects_plus_one(self) -> None:
        self.assertEqual(MAX_SOURCE_JSON_BYTES, 64 * 1024 * 1024)
        boundary = self.root / "source-boundary.json"
        _write_padded_json(boundary, self.chain.source, MAX_SOURCE_JSON_BYTES)
        source, records = load_source(boundary)
        self.assertEqual(source, self.chain.source)
        self.assertEqual(records, {RECORD_ID: RECORD_BASES})

        with boundary.open("ab") as stream:
            stream.write(b" ")
        with self.assertRaisesRegex(V2Error, f"exceeds byte limit {MAX_SOURCE_JSON_BYTES}"):
            load_source(boundary)

    def test_legacy_source_retains_exact_16_mib_boundary(self) -> None:
        legacy = _Chain(self.root / "legacy-boundary")
        boundary = self.root / "legacy-source-boundary.json"
        _write_padded_json(boundary, legacy.source, MAX_JSON_BYTES)
        source, records = load_source(boundary)
        self.assertEqual(source, legacy.source)
        self.assertEqual(records, legacy.records)

        with boundary.open("ab") as stream:
            stream.write(b" ")
        with self.assertRaisesRegex(V2Error, f"JSON byte limit {MAX_JSON_BYTES}"):
            load_source(boundary)

    def test_non_source_chain_artifacts_remain_at_16_mib(self) -> None:
        self.assertEqual(MAX_JSON_BYTES, 16 * 1024 * 1024)
        oversized = self.root / "manifest-oversized.json"
        _write_padded_json(oversized, self.chain.manifest, MAX_JSON_BYTES + 1)
        with self.assertRaisesRegex(V2Error, f"exceeds byte limit {MAX_JSON_BYTES}"):
            load_manifest(oversized)

    def test_legacy_source_request_and_module_hashes_are_unchanged(self) -> None:
        legacy = _Chain(self.root / "legacy")
        module = legacy.compile()
        self.assertEqual(
            legacy.source["artifact_sha256"],
            "117cb7f85fa9ab6fa27d0618d2666764015136394b179a524c4850ee47eacd91",
        )
        self.assertEqual(
            legacy.request["artifact_sha256"],
            "419cee703bebedc5dfad742fe6bf65733651ca63a2349fa21b6434e7fc4a1338",
        )
        self.assertEqual(
            module["artifact_sha256"],
            "4f786a7ab8b80f27c74dec36e635648952b390de6793b337daf8320842ec819c",
        )
        self.assertEqual(
            module["module_sha256"],
            "a71dc79b721db4eb8a89d74fbab1e072f89d4babff6bc05bd2fb8ccfb4ed0a8c",
        )


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import base64
from copy import deepcopy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import struct
import sys
import tempfile
import unittest
from unittest import mock
from typing import Any, Callable

from brainc._canonical import artifact_digest, digest
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
from brainc.v2._common import loads, pretty_bytes
from brainc.v2.limits import MAX_JSON_BYTES, MAX_STRING_BYTES
from brainc.v2.policy import (
    OP_EDGE_CREATE,
    OP_PORT_BIND,
    OP_RULE_ATTACH,
    OP_UNIT_CREATE,
    parse_policy,
)
from brainc.v2.provider import load_request, load_source
from brainc.v2.target import DEV_DOMAIN, parse_target_artifact
from brainc.v2.tensor import parse_storage, parse_type
from brainc.sequence_collection import SequenceCollectionCompiler


ROOT = Path(__file__).resolve().parents[1]
PHIX_FASTA = ROOT / "tests/data/J02482.1.fasta"
ZERO_SHA = "0" * 64
ONE_SHA = "1" * 64
UNIT = "org.example.signal"


def _seal(core: dict[str, Any]) -> dict[str, Any]:
    return {**core, "artifact_sha256": digest(core)}


def _reseal(artifact: dict[str, Any]) -> dict[str, Any]:
    return _seal({key: value for key, value in artifact.items() if key != "artifact_sha256"})


def _field(
    identifier: str,
    dtype: str,
    *,
    unit: str | None = None,
    mutability: str = "state",
) -> dict[str, Any]:
    kinds = {"bool": "boolean", "i64": "integer", "u64": "integer", "f64": "float"}
    return {
        "id": identifier,
        "type": {"dtype": dtype, "shape": []},
        "unit": unit,
        "mutability": mutability,
        "numeric": {"kind": kinds[dtype]},
    }


def _contract() -> dict[str, Any]:
    return {
        "id": "org.example.machine",
        "abi_major": 1,
        "opsets": [{"domain": DEV_DOMAIN, "version": 1}],
        "unit_schemas": [
            {
                "id": "node",
                "fields": [
                    _field("flag", "bool"),
                    _field("value", "f64", unit=UNIT),
                ],
            }
        ],
        "edge_schemas": [
            {"id": "link", "fields": [_field("weight", "i64")]},
        ],
        "ports": [
            {
                "id": "input.value",
                "direction": "input",
                "type": {"dtype": "f64", "shape": [2]},
                "unit": UNIT,
                "codec": {"id": "org.example.raw", "version": 1, "semantics_sha256": ONE_SHA},
            },
            {
                "id": "output.value",
                "direction": "output",
                "type": {"dtype": "f64", "shape": [2]},
                "unit": UNIT,
                "codec": {"id": "org.example.raw", "version": 1, "semantics_sha256": ONE_SHA},
            },
        ],
        "rules": [
            {
                "id": "org.example.adjust",
                "version": 1,
                "semantics_sha256": ZERO_SHA,
                "subject": "edge",
                "schema": "link",
                "phase": "UPDATE",
                "triggers": ["tick"],
                "parameters": [
                    {"id": "rate", "type": {"dtype": "f64", "shape": []}, "unit": None}
                ],
                "reads": [
                    {"scope": "source", "field": "value"},
                    {"scope": "subject", "field": "weight"},
                    {"scope": "target", "field": "value"},
                ],
                "writes": [{"scope": "subject", "field": "weight"}],
            }
        ],
    }


def _operations() -> list[dict[str, Any]]:
    return [
        {
            "id": "create.nodes",
            "op": OP_UNIT_CREATE,
            "version": 1,
            "schema": "node",
            "count": "t.count",
            "initializers": [
                {"field": "flag", "tensor": "t.flag"},
                {"field": "value", "tensor": "t.value"},
            ],
        },
        {
            "id": "create.links",
            "op": OP_EDGE_CREATE,
            "version": 1,
            "schema": "link",
            "sources": "create.nodes",
            "targets": "create.nodes",
            "pairs": "t.pairs",
            "initializers": [{"field": "weight", "tensor": "t.weight"}],
        },
        {
            "id": "attach.adjust",
            "op": OP_RULE_ATTACH,
            "version": 1,
            "rule": "org.example.adjust",
            "rule_version": 1,
            "subject": "create.links",
            "parameters": [{"id": "rate", "tensor": "t.rate"}],
        },
        {
            "id": "bind.input",
            "op": OP_PORT_BIND,
            "version": 1,
            "port": "input.value",
            "subject": "create.nodes",
            "field": "value",
            "indices": "t.indices.in",
        },
        {
            "id": "bind.output",
            "op": OP_PORT_BIND,
            "version": 1,
            "port": "output.value",
            "subject": "create.nodes",
            "field": "value",
            "indices": "t.indices.out",
        },
    ]


def _tensor(
    identifier: str,
    dtype: str,
    shape: list[int],
    values: list[Any],
    unit: str | None = None,
    axes: list[dict[str, Any] | None] | None = None,
) -> dict[str, Any]:
    return {
        "id": identifier,
        "type": {"dtype": dtype, "shape": shape},
        "unit": unit,
        "axes": [None] * len(shape) if axes is None else axes,
        "storage": inline_storage(pack(dtype, values)),
    }


class _Chain:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.source_path = root / "sequence.json"
        SequenceCollectionCompiler().compile_file(PHIX_FASTA).save(self.source_path)
        self.source, self.records = load_source(self.source_path)
        self.counter = 0

        self.target = target_artifact(_contract())
        self.target_path = self.write("target", self.target)
        identity = {
            "id": self.target["contract"]["id"],
            "abi_major": self.target["contract"]["abi_major"],
            "contract_sha256": self.target["contract_sha256"],
        }
        axis = {
            "kind": "source-coordinate",
            "source_artifact_sha256": self.source["artifact_sha256"],
            "record_id": "J02482.1",
            "start": 7,
            "step": 2,
            "span": 1,
            "coordinate_system": "0-based-half-open",
        }
        self.outputs = sorted(
            [
                _tensor("t.count", "u64", [], [2]),
                _tensor("t.flag", "bool", [], [False]),
                _tensor("t.indices.in", "u64", [2], [0, 1]),
                _tensor("t.indices.out", "u64", [2], [1, 0]),
                _tensor("t.pairs", "u64", [2, 2], [0, 1, 1, 0]),
                _tensor("t.rate", "f64", [], [0.5]),
                _tensor("t.value", "f64", [2], [1.25, 2.5], UNIT, [axis]),
                _tensor("t.weight", "i64", [2], [7, -3]),
            ],
            key=lambda item: item["id"],
        )
        self.output_ids = [item["id"] for item in self.outputs]
        self.manifest = _seal(
            {
                "format": "brainc.provider-manifest",
                "version": 2,
                "provider": {"name": "org.example.provider", "version": "1"},
                "model_identity": {
                    "kind": "content-sha256",
                    "value": digest({"algorithm": "fixture-independent-identity", "version": 1}),
                },
                "accepts": ["brain01.sequence-collection-ir/v1"],
                "outputs": [
                    {key: deepcopy(value[key]) for key in ("id", "type", "unit", "axes")}
                    for value in self.outputs
                ],
            }
        )
        self.manifest_path = self.write("manifest", self.manifest)
        self.request = make_request(self.source_path, self.manifest_path, self.output_ids)
        self.request_path = self.write("request", self.request)
        self.response = _seal(
            {
                "format": "brainc.prediction-response",
                "version": 2,
                "request_artifact_sha256": self.request["artifact_sha256"],
                "provider": deepcopy(self.manifest["provider"]),
                "model_identity": deepcopy(self.manifest["model_identity"]),
                "outputs": deepcopy(self.outputs),
            }
        )
        self.response_path = self.write("response", self.response)
        self.policy = policy_artifact(
            "org.example.lowering",
            identity,
            [{"id": identifier, "from_output": identifier} for identifier in self.output_ids],
            _operations(),
        )
        self.policy_path = self.write("policy", self.policy)

    def write(self, name: str, artifact: dict[str, Any]) -> Path:
        self.counter += 1
        path = self.root / f"{self.counter:02d}-{name}.json"
        save(artifact, path)
        return path

    def compile(self, **overrides: Path | None) -> dict[str, Any]:
        paths: dict[str, Path | None] = {
            "source_path": self.source_path,
            "manifest_path": self.manifest_path,
            "request_path": self.request_path,
            "response_path": self.response_path,
            "policy_path": self.policy_path,
            "target_path": self.target_path,
            "blob_root": None,
        }
        paths.update(overrides)
        return compile_module(**paths)  # type: ignore[arg-type]

    def changed_response(
        self, mutate: Callable[[dict[str, Any]], None]
    ) -> Path:
        response = deepcopy(self.response)
        mutate(response)
        return self.write("changed-response", _reseal(response))

    def changed_policy(self, mutate: Callable[[dict[str, Any]], None]) -> Path:
        policy = deepcopy(self.policy)
        mutate(policy)
        return self.write("changed-policy", _reseal(policy))


def _output(artifact: dict[str, Any], identifier: str) -> dict[str, Any]:
    return next(value for value in artifact["outputs"] if value["id"] == identifier)


class CompilerV2AdversarialTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.chain = _Chain(Path(self.temporary.name))

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def assertCompileRejects(self, **overrides: Path | None) -> None:
        with self.assertRaises(V2Error):
            self.chain.compile(**overrides)

    def test_deterministic_module_is_closed_exactly_bound_and_budgeted(self) -> None:
        first = self.chain.compile()
        second = self.chain.compile()
        self.assertEqual(first, second)
        self.assertEqual(artifact_digest(first), first["artifact_sha256"])
        self.assertEqual(digest(first["module"]), first["module_sha256"])
        self.assertEqual(first["format"], "brainc.development-module")
        self.assertEqual(first["version"], 1)
        self.assertEqual(
            first["module"]["budgets"],
            {"operations": 5, "tensor_bytes": 113, "units": 2, "edges": 2, "attachments": 1},
        )
        self.assertEqual(first["module"]["requirements"], [{"domain": DEV_DOMAIN, "version": 1}])
        expected_sources = {
            "sequence": self.chain.source["artifact_sha256"],
            "provider_manifest": self.chain.manifest["artifact_sha256"],
            "prediction_request": self.chain.request["artifact_sha256"],
            "prediction_response": self.chain.response["artifact_sha256"],
            "lowering_policy": self.chain.policy["artifact_sha256"],
            "target_contract": self.chain.target["artifact_sha256"],
        }
        for role, expected in expected_sources.items():
            self.assertEqual(first["sources"][role]["artifact_sha256"], expected)
        self.assertEqual(
            [tensor["id"] for tensor in first["module"]["tensors"]],
            self.chain.output_ids,
        )
        self.assertTrue(
            all(
                tensor["lineage"]
                == {"response_output": tensor["id"], "lowering_value": tensor["id"]}
                for tensor in first["module"]["tensors"]
            )
        )

    def test_every_artifact_layer_rejects_unknown_or_missing_keys(self) -> None:
        manifest = deepcopy(self.chain.manifest)
        manifest["surprise"] = None
        manifest_path = self.chain.write("bad-manifest", _reseal(manifest))

        request = deepcopy(self.chain.request)
        del request["requested_outputs"][0]["unit"]
        request_path = self.chain.write("bad-request", _reseal(request))

        response = deepcopy(self.chain.response)
        _output(response, "t.rate")["storage"]["surprise"] = 1
        response_path = self.chain.write("bad-response", _reseal(response))

        policy = deepcopy(self.chain.policy)
        policy["operations"][0]["surprise"] = 1
        policy_path = self.chain.write("bad-policy", _reseal(policy))

        target = deepcopy(self.chain.target)
        target["contract"]["unit_schemas"][0]["surprise"] = 1
        target["contract_sha256"] = digest(target["contract"])
        target_path = self.chain.write("bad-target", _reseal(target))

        for override in (
            {"manifest_path": manifest_path},
            {"request_path": request_path},
            {"response_path": response_path},
            {"policy_path": policy_path},
            {"target_path": target_path},
        ):
            with self.subTest(role=next(iter(override))):
                self.assertCompileRejects(**override)

    def test_exact_python_types_and_bool_substitution_are_rejected(self) -> None:
        contract = _contract()
        contract["abi_major"] = True
        with self.assertRaises(V2Error):
            target_artifact(contract)

        contract = _contract()
        contract["opsets"][0]["version"] = 1.0
        with self.assertRaises(V2Error):
            target_artifact(contract)

        policy = deepcopy(self.chain.policy)
        policy["operations"][0]["version"] = True
        with self.assertRaises(V2Error):
            parse_policy(_reseal(policy))

        response_path = self.chain.changed_response(
            lambda value: _output(value, "t.rate")["storage"].__setitem__("byte_length", True)
        )
        self.assertCompileRejects(response_path=response_path)

        for dtype, values in (("u64", [True]), ("bool", [1]), ("f64", [1])):
            with self.subTest(dtype=dtype):
                with self.assertRaises(V2Error):
                    pack(dtype, values)
        with self.assertRaises(V2Error):
            parse_type({"dtype": "u64", "shape": [True]}, "attack")
        with self.assertRaises(V2Error):
            make_request(self.chain.source_path, self.chain.manifest_path, tuple(self.chain.output_ids))  # type: ignore[arg-type]

        for field, value in (
            ("format", []),
            ("version", True),
            ("version", 1.0),
            ("version", []),
        ):
            with self.subTest(request_source_field=field, value=value):
                request = deepcopy(self.chain.request)
                request["source"][field] = value
                request_path = self.chain.write(
                    "bad-request-source-identity",
                    _reseal(request),
                )
                with self.assertRaisesRegex(V2Error, "format/version is unsupported"):
                    load_request(request_path)

    def test_json_boundary_rejects_duplicates_nonfinite_unsafe_and_oversize(self) -> None:
        for raw in (
            b'{"x":1,"x":2}',
            b'{"x":NaN}',
            b'{"x":Infinity}',
            b'{"x":9007199254740992}',
        ):
            with self.subTest(raw=raw):
                with self.assertRaises(V2Error):
                    loads(raw, "attack")
        with self.assertRaises(V2Error):
            loads(b" " * (MAX_JSON_BYTES + 1), "oversize")
        with self.assertRaises(V2Error):
            pretty_bytes({"x": "x" * (MAX_STRING_BYTES + 1)})

    def test_nonfinite_negative_zero_and_invalid_bool_tensor_bytes_are_rejected(self) -> None:
        attacks = (
            ("t.rate", struct.pack("<d", -0.0)),
            ("t.rate", struct.pack("<d", math.nan)),
            ("t.rate", struct.pack("<d", math.inf)),
            ("t.flag", b"\x02"),
        )
        for identifier, raw in attacks:
            def mutate(value: dict[str, Any], identifier: str = identifier, raw: bytes = raw) -> None:
                _output(value, identifier)["storage"] = inline_storage(raw)

            with self.subTest(identifier=identifier, raw=raw.hex()):
                self.assertCompileRejects(response_path=self.chain.changed_response(mutate))
        for value in (-0.0, math.nan, math.inf, -math.inf):
            with self.subTest(value=value):
                with self.assertRaises(V2Error):
                    pack("f64", [value])

    def test_inline_storage_is_length_digest_and_canonical_base64_bound(self) -> None:
        def wrong_digest(value: dict[str, Any]) -> None:
            _output(value, "t.rate")["storage"]["sha256"] = ZERO_SHA

        def wrong_length(value: dict[str, Any]) -> None:
            _output(value, "t.rate")["storage"]["byte_length"] = 7

        def whitespace(value: dict[str, Any]) -> None:
            storage = _output(value, "t.rate")["storage"]
            storage["data"] = storage["data"] + "\n"

        def noncanonical(value: dict[str, Any]) -> None:
            storage = _output(value, "t.rate")["storage"]
            storage["data"] = base64.b64encode(b"12345678").decode("ascii").rstrip("=")

        for name, attack in (
            ("digest", wrong_digest),
            ("length", wrong_length),
            ("whitespace", whitespace),
            ("base64", noncanonical),
        ):
            with self.subTest(name=name):
                self.assertCompileRejects(response_path=self.chain.changed_response(attack))

    def test_external_blob_requires_explicit_safe_content_addressed_storage(self) -> None:
        if os.name != "posix":
            self.skipTest("positive external-blob profile requires POSIX descriptors")
        raw = pack("i64", [7, -3])
        blob_sha = hashlib.sha256(raw).hexdigest()
        blob_root = Path(self.temporary.name) / "blobs"
        blob_root.mkdir()
        (blob_root / blob_sha).write_bytes(raw)

        def external(value: dict[str, Any]) -> None:
            _output(value, "t.weight")["storage"] = {
                "kind": "sha256-blob",
                "byte_length": len(raw),
                "sha256": blob_sha,
            }

        response_path = self.chain.changed_response(external)
        self.chain.compile(response_path=response_path, blob_root=blob_root)
        self.assertCompileRejects(response_path=response_path)

        alias = blob_root / "alias"
        alias.hardlink_to(blob_root / blob_sha)
        self.assertCompileRejects(response_path=response_path, blob_root=blob_root)

    def test_external_blob_fails_closed_without_posix_descriptor_guarantees(self) -> None:
        raw = pack("i64", [7, -3])
        storage = {
            "kind": "sha256-blob",
            "byte_length": len(raw),
            "sha256": hashlib.sha256(raw).hexdigest(),
        }
        with mock.patch("brainc.v2.tensor.os.name", "nt"):
            with self.assertRaisesRegex(V2Error, "requires a POSIX host"):
                parse_storage(storage, len(raw), "external", blob_root=self.temporary.name)

    def test_all_cross_artifact_bindings_reject_resealed_tampering(self) -> None:
        request_source = deepcopy(self.chain.request)
        request_source["source"]["artifact_sha256"] = ONE_SHA

        response_request = deepcopy(self.chain.response)
        response_request["request_artifact_sha256"] = ONE_SHA

        response_provider = deepcopy(self.chain.response)
        response_provider["provider"]["version"] = "2"

        policy_target = deepcopy(self.chain.policy)
        policy_target["target"]["contract_sha256"] = ONE_SHA

        policy_values = deepcopy(self.chain.policy)
        policy_values["values"].pop()

        attacks = (
            {"request_path": self.chain.write("request-source", _reseal(request_source))},
            {"response_path": self.chain.write("response-request", _reseal(response_request))},
            {"response_path": self.chain.write("response-provider", _reseal(response_provider))},
            {"policy_path": self.chain.write("policy-target", _reseal(policy_target))},
            {"policy_path": self.chain.write("policy-values", _reseal(policy_values))},
        )
        for attack in attacks:
            with self.subTest(role=next(iter(attack))):
                self.assertCompileRejects(**attack)

    def test_axes_are_checked_against_the_actual_source_after_contract_binding(self) -> None:
        base_axis = deepcopy(_output(self.chain.manifest, "t.value")["axes"][0])
        attacks = []
        wrong_source = deepcopy(base_axis)
        wrong_source["source_artifact_sha256"] = ONE_SHA
        attacks.append(wrong_source)
        unknown_record = deepcopy(base_axis)
        unknown_record["record_id"] = "not-a-record"
        attacks.append(unknown_record)
        overrun = deepcopy(base_axis)
        overrun.update({"start": 5385, "step": 2, "span": 1})
        attacks.append(overrun)

        for axis in attacks:
            manifest = deepcopy(self.chain.manifest)
            _output(manifest, "t.value")["axes"] = [axis]
            manifest = _reseal(manifest)
            manifest_path = self.chain.write("axis-manifest", manifest)
            with self.subTest(axis=axis):
                with self.assertRaises(V2Error):
                    make_request(
                        self.chain.source_path,
                        manifest_path,
                        self.chain.output_ids,
                    )

    def test_resealed_provider_values_still_face_operation_semantics(self) -> None:
        attacks = {
            "duplicate edge": ("t.pairs", "u64", [0, 1, 0, 1]),
            "edge endpoint": ("t.pairs", "u64", [0, 2, 1, 0]),
            "duplicate input index": ("t.indices.in", "u64", [0, 0]),
            "unsafe ABI integer": ("t.indices.out", "u64", [0, 2**53]),
        }
        for name, (identifier, dtype, values) in attacks.items():
            def mutate(
                value: dict[str, Any],
                identifier: str = identifier,
                dtype: str = dtype,
                values: list[int] = values,
            ) -> None:
                _output(value, identifier)["storage"] = inline_storage(pack(dtype, values))

            with self.subTest(name=name):
                self.assertCompileRejects(response_path=self.chain.changed_response(mutate))

    def test_invalid_operation_order_types_completeness_and_references_are_rejected(self) -> None:
        def missing_initializer(value: dict[str, Any]) -> None:
            value["operations"][0]["initializers"].pop()

        def edge_before_units(value: dict[str, Any]) -> None:
            value["operations"][0], value["operations"][1] = (
                value["operations"][1],
                value["operations"][0],
            )

        def wrong_rule_parameter(value: dict[str, Any]) -> None:
            value["operations"][2]["parameters"][0]["tensor"] = "t.count"

        def wrong_port_field(value: dict[str, Any]) -> None:
            value["operations"][4]["field"] = "flag"

        def unbound_port(value: dict[str, Any]) -> None:
            value["operations"].pop()

        def duplicate_attachment(value: dict[str, Any]) -> None:
            duplicate = deepcopy(value["operations"][2])
            duplicate["id"] = "attach.second"
            value["operations"].insert(3, duplicate)

        for name, attack in (
            ("initializer", missing_initializer),
            ("order", edge_before_units),
            ("parameter", wrong_rule_parameter),
            ("port", wrong_port_field),
            ("completeness", unbound_port),
            ("attachment", duplicate_attachment),
        ):
            with self.subTest(name=name):
                self.assertCompileRejects(policy_path=self.chain.changed_policy(attack))

    def test_target_semantics_reject_unsorted_ambiguous_or_unresolved_contracts(self) -> None:
        attacks: list[Callable[[dict[str, Any]], None]] = []

        def unsorted_fields(value: dict[str, Any]) -> None:
            value["unit_schemas"][0]["fields"].reverse()

        attacks.append(unsorted_fields)

        def constant_write(value: dict[str, Any]) -> None:
            value["edge_schemas"][0]["fields"][0]["mutability"] = "constant"

        attacks.append(constant_write)

        def unresolved_source(value: dict[str, Any]) -> None:
            value["rules"][0]["reads"][0]["field"] = "absent"

        attacks.append(unresolved_source)

        def duplicate_opset(value: dict[str, Any]) -> None:
            value["opsets"].append(deepcopy(value["opsets"][0]))

        attacks.append(duplicate_opset)

        for attack in attacks:
            contract = _contract()
            attack(contract)
            with self.subTest(attack=attack.__name__):
                with self.assertRaises(V2Error):
                    target_artifact(contract)

    def test_target_opset_and_policy_target_are_both_required(self) -> None:
        contract = _contract()
        contract["opsets"] = [{"domain": "org.example.other", "version": 1}]
        target = target_artifact(contract)
        target_path = self.chain.write("other-opset-target", target)
        policy = deepcopy(self.chain.policy)
        policy["target"] = {
            "id": target["contract"]["id"],
            "abi_major": target["contract"]["abi_major"],
            "contract_sha256": target["contract_sha256"],
        }
        policy_path = self.chain.write("other-opset-policy", _reseal(policy))
        self.assertCompileRejects(target_path=target_path, policy_path=policy_path)

    def test_output_helpers_reject_aliasing_and_noncanonical_values(self) -> None:
        contract = _contract()
        target = target_artifact(contract)
        contract["id"] = "org.example.mutated"
        self.assertEqual(target["contract"]["id"], "org.example.machine")

        values = [{"id": "x", "from_output": "x"}]
        operations: list[dict[str, Any]] = []
        policy = policy_artifact(
            "org.example.alias-test",
            {
                "id": target["contract"]["id"],
                "abi_major": 1,
                "contract_sha256": target["contract_sha256"],
            },
            values,
            operations,
        )
        values[0]["id"] = "changed"
        self.assertEqual(policy["values"][0]["id"], "x")
        with self.assertRaises(V2Error):
            inline_storage(bytearray(b"x"))  # type: ignore[arg-type]

    def test_v2_save_is_atomic_and_rejects_live_or_broken_symlinks(self) -> None:
        root = Path(self.temporary.name)
        victim = root / "victim.json"
        victim.write_text("preserve", encoding="utf-8")
        linked = root / "linked.json"
        linked.symlink_to(victim)
        with self.assertRaisesRegex(V2Error, "regular non-linked"):
            save(self.chain.target, linked)
        self.assertEqual(victim.read_text(encoding="utf-8"), "preserve")

        broken = root / "broken.json"
        broken.symlink_to(root / "missing.json")
        with self.assertRaisesRegex(V2Error, "regular non-linked"):
            save(self.chain.target, broken)


class RealPhiXRuntimeOracleTests(unittest.TestCase):
    def test_real_phix_compiles_and_develops_to_exact_body_and_state(self) -> None:
        runtime_root = ROOT.parent / "epigenesis-runtime"
        reference_path = runtime_root / "examples/reference_chain.py"
        if not reference_path.is_file():
            self.skipTest("sibling Epigenesis runtime source is not present")
        sys.path.insert(0, str(runtime_root))
        try:
            spec = importlib.util.spec_from_file_location("epigenesis_reference_chain", reference_path)
            if spec is None or spec.loader is None:
                self.fail("cannot load the reference-chain oracle")
            reference = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(reference)
            from epirun._canonical import load as runtime_load
            from epirun.contracts import load_development_module, load_target_contract
            from epirun.development import develop

            with tempfile.TemporaryDirectory() as temporary:
                chain = reference.build_chain(
                    Path(temporary) / "chain",
                    PHIX_FASTA,
                    observation_width=16,
                    action_count=8,
                )
                source, _ = runtime_load(chain["source"], "source oracle")
                target_raw, _ = runtime_load(chain["target"], "target oracle")
                module_raw, _ = runtime_load(chain["module"], "module oracle")
                target = load_target_contract(chain["target"])
                module = load_development_module(
                    chain["module"], target=target, blob_root=chain["blob_root"]
                )
                developed = develop(target, module)

            member = source["members"][0]
            self.assertEqual(member["record_id"], "J02482.1")
            self.assertEqual(
                len(member["sequence_artifact"]["sequence_ir"]["sequence"]), 5_386
            )
            self.assertEqual(
                target_raw["contract_sha256"],
                "03da2729deeeb4ab8457652d9be67e405880b17cde785111bb5522ff11d8a1ee",
            )
            self.assertEqual(
                module_raw["module_sha256"],
                "c80ed3dbca46519ea505575257af664efbd8ea502d67a1db2ca7633283b24873",
            )
            self.assertEqual(
                module_raw["module"]["budgets"],
                {"operations": 6, "tensor_bytes": 2330, "units": 24, "edges": 128, "attachments": 1},
            )
            self.assertEqual(developed.state.unit_count, 24)
            self.assertEqual(developed.state.edge_count, 128)
            self.assertEqual(len(developed.state.attachments), 1)
            self.assertEqual(len(developed.state.ports), 2)
            self.assertEqual(
                developed.final_state_sha256,
                "16a625b21f1ccc51cc3322ec203d99a0e35ada578f64356d14a7d6a908d7c8d4",
            )
        finally:
            try:
                sys.path.remove(str(runtime_root))
            except ValueError:
                pass


if __name__ == "__main__":
    unittest.main()

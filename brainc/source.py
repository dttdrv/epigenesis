"""DNA source admission through explicit profile contracts."""

from __future__ import annotations

from collections.abc import Mapping
import copy
from dataclasses import dataclass
import hashlib
import os
from pathlib import Path
import re
import stat
from typing import Any

from ._canonical import ContractError, SAFE_INTEGER, artifact_digest, digest, sha256
from ._io import (
    BoundedIOError,
    MAX_IDENTIFIER_BYTES,
    MAX_INPUT_BYTES,
    MAX_JSON_BYTES,
    MAX_JSON_DEPTH,
    MAX_JSON_MEMBERS,
    MAX_SOURCE_RECORDS,
    MAX_STRING_BYTES,
    load_json_object,
    pretty_json_bytes,
    read_regular_file,
)
from ._publish import PublicationError, publish_directory
from .bio import (
    FORMAT as GFF3_FORMAT,
    MAX_GFF3_BYTES,
    MAX_JSON_DEPTH as GFF3_JSON_DEPTH,
    MAX_JSON_MEMBERS as GFF3_JSON_MEMBERS,
    MAX_OUTPUT_BYTES as GFF3_OUTPUT_BYTES,
    VERSION as GFF3_VERSION,
    GFF3Compiler,
    GFF3Error,
    validate_gff3_artifact,
)
from .insdc import (
    FORMAT as GENBANK_FORMAT,
    MAX_GENBANK_BYTES,
    MAX_JSON_DEPTH as GENBANK_JSON_DEPTH,
    MAX_OUTPUT_BYTES as GENBANK_OUTPUT_BYTES,
    PROFILE as GENBANK_PROFILE,
    VERSION as GENBANK_VERSION,
    GenBankCompiler,
    GenBankError,
    validate_genbank_artifact,
)
from .external_source import (
    FORMAT as EXTERNAL_FORMAT,
    EXECUTABLE_VERSION as EXTERNAL_EXECUTABLE_VERSION,
    MAX_CLOSURE_BYTES as EXTERNAL_OUTPUT_BYTES,
    PROFILE as EXTERNAL_PROFILE,
    VERSION as EXTERNAL_VERSION,
    ExternalSourceError,
    build_external_source_closure,
    build_external_source_closure_from_paths,
    build_executable_external_source_closure,
    input_references as external_input_references,
    profile_parameters as external_profile_parameters,
    source_records as external_source_records,
    validate_external_source_closure,
)
from .sequence_collection import (
    FORMAT as COLLECTION_FORMAT,
    VERSION as COLLECTION_VERSION,
    SequenceCollectionArtifact,
    SequenceCollectionCompiler,
    SequenceCollectionError,
    _as_collection,
)
from .source_scale import (
    DEFAULT_LIMITS as DEFAULT_SCALE_LIMITS,
    FORMAT as REFERENCE_FORMAT,
    PROFILE as REFERENCE_FASTA_PROFILE,
    VERSION as REFERENCE_VERSION,
    ScaleLimits,
    SourceScaleError,
    compile_reference_fasta,
    validate_reference_catalog,
)


FORMAT = "brainc.source-descriptor"
VERSION = 1
RAW_PROFILE = "raw-iupac-dna/v1"
FASTA_PROFILE = "fasta-dna/v1"
GFF3_PROFILE = "gff3-external-sequence/v1"
PROFILES = frozenset(
    {
        RAW_PROFILE,
        FASTA_PROFILE,
        REFERENCE_FASTA_PROFILE,
        EXTERNAL_PROFILE,
        GENBANK_PROFILE,
        GFF3_PROFILE,
    }
)
PRODUCER = {
    "name": "brainc-source",
    "version": "1.0.0",
    "passes": [
        "select-source-profile",
        "compile-native-source",
        "bind-source-evidence",
        "emit-source-descriptor",
    ],
}

SOURCE_FILENAME = "source.json"
ROLE_FILENAMES = {
    "sequence": "sequence.json",
    "reference": "reference.json",
    "external": "external.json",
    "genbank": "genbank.json",
    "annotation": "annotation.json",
}
_REFGET_RE = re.compile(r"SQ\.[A-Za-z0-9_-]{32}\Z")
MAX_SOURCE_DESCRIPTOR_BYTES = GFF3_OUTPUT_BYTES


class SourceError(ContractError):
    """A source route, descriptor, or native source closure is invalid."""


@dataclass(frozen=True, slots=True)
class _Route:
    input_roles: tuple[str, ...]
    artifact_roles: tuple[str, ...]


def _fail(detail: str) -> SourceError:
    return SourceError(f"SOURCE001: {detail}")


def _mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise _fail(f"{label} must be a mapping")
    result = dict(value)
    if any(type(key) is not str for key in result):
        raise _fail(f"{label} keys must be strings")
    return result


def _keys(value: Any, expected: set[str], label: str) -> dict[str, Any]:
    if type(value) is not dict:
        raise _fail(f"{label} must be an object")
    missing = sorted(expected - value.keys())
    extra = sorted(value.keys() - expected)
    if missing or extra:
        raise _fail(
            f"{label} keys invalid; missing={missing or 'none'}, unknown={extra or 'none'}"
        )
    return value


def _integer(value: Any, label: str, *, minimum: int = 0) -> int:
    if type(value) is not int or not minimum <= value <= SAFE_INTEGER:
        raise _fail(f"{label} must be an integer in [{minimum}, {SAFE_INTEGER}]")
    return value


def _text(value: Any, label: str, *, identifier: bool = False) -> str:
    if type(value) is not str or not value or value != value.strip():
        raise _fail(f"{label} must be a nonempty trimmed string")
    if any(0xD800 <= ord(character) <= 0xDFFF for character in value):
        raise _fail(f"{label} contains a lone Unicode surrogate")
    ceiling = MAX_IDENTIFIER_BYTES if identifier else MAX_STRING_BYTES
    if len(value.encode("utf-8")) > ceiling:
        raise _fail(f"{label} exceeds {ceiling} UTF-8 bytes")
    if identifier and any(character.isspace() for character in value):
        raise _fail(f"{label} must be whitespace-free")
    return value


def _route(profile: Any, parameters_value: Any) -> tuple[_Route, dict[str, Any]]:
    if type(profile) is not str or profile not in PROFILES:
        raise _fail("source profile is unsupported")
    parameters = _mapping(parameters_value, "parameters")
    if profile == RAW_PROFILE:
        _keys(parameters, {"record_id"}, "parameters")
        _text(parameters["record_id"], "parameters.record_id", identifier=True)
        return _Route(("sequence",), ("sequence",)), copy.deepcopy(parameters)
    if profile in {FASTA_PROFILE, REFERENCE_FASTA_PROFILE}:
        _keys(parameters, {"wrapper"}, "parameters")
        if (
            type(parameters["wrapper"]) is not str
            or parameters["wrapper"] not in {"identity", "gzip"}
        ):
            raise _fail("parameters.wrapper must be identity or gzip")
        artifacts = (
            ("reference",)
            if profile == REFERENCE_FASTA_PROFILE
            else ("sequence",)
        )
        return _Route(("sequence",), artifacts), copy.deepcopy(parameters)
    if profile == GENBANK_PROFILE:
        _keys(parameters, set(), "parameters")
        return _Route(("genbank",), ("genbank",)), {}
    if profile == EXTERNAL_PROFILE:
        _keys(
            parameters,
            {"profile_id", "profile_version", "profile_manifest_sha256"},
            "parameters",
        )
        _text(parameters["profile_id"], "parameters.profile_id", identifier=True)
        _integer(parameters["profile_version"], "parameters.profile_version", minimum=1)
        sha256(
            parameters["profile_manifest_sha256"],
            "parameters.profile_manifest_sha256",
        )
        return _Route((), ("external",)), copy.deepcopy(parameters)

    _keys(parameters, {"sequence_profile", "sequence_parameters"}, "parameters")
    sequence_profile = parameters["sequence_profile"]
    if (
        type(sequence_profile) is not str
        or sequence_profile not in {RAW_PROFILE, FASTA_PROFILE}
    ):
        raise _fail("GFF3 sequence_profile must be raw IUPAC or FASTA")
    _, sequence_parameters = _route(sequence_profile, parameters["sequence_parameters"])
    normalized = {
        "sequence_profile": sequence_profile,
        "sequence_parameters": sequence_parameters,
    }
    return _Route(("sequence", "annotation"), ("sequence", "annotation")), normalized


def _input_paths(inputs_value: Any, expected: tuple[str, ...]) -> dict[str, Path]:
    inputs = _mapping(inputs_value, "inputs")
    _keys(inputs, set(expected), "inputs")
    result: dict[str, Path] = {}
    for role in expected:
        value = inputs[role]
        if not isinstance(value, (str, Path)) or not os.fspath(value):
            raise _fail(f"inputs.{role} must be a filesystem path")
        result[role] = Path(value)
    return result


def _read(path: Path, maximum: int, label: str) -> bytes:
    try:
        return read_regular_file(path, maximum_bytes=maximum, label=label)
    except BoundedIOError as failure:
        raise _fail(str(failure)) from failure


def _compile_sequence(
    profile: str,
    parameters: dict[str, Any],
    raw: bytes,
) -> SequenceCollectionArtifact:
    compiler = SequenceCollectionCompiler()
    try:
        if profile == RAW_PROFILE:
            artifact = compiler.compile_raw(raw, parameters["record_id"])
            expected_wrapper = None
            expected_kind = "raw-iupac"
        else:
            artifact = compiler.compile_fasta_bytes(raw)
            expected_wrapper = None if parameters["wrapper"] == "identity" else "gzip"
            expected_kind = "fasta"
    except SequenceCollectionError as failure:
        raise _fail(str(failure)) from failure
    payload = artifact.to_dict()
    if payload["inputs"]["kind"] != expected_kind or payload["inputs"]["wrapper"] != expected_wrapper:
        raise _fail("source bytes do not match the explicitly selected wrapper/profile")
    return artifact


def _input_reference(raw: bytes) -> dict[str, Any]:
    return {
        "sha256": hashlib.sha256(raw).hexdigest(),
        "byte_length": len(raw),
    }


def _artifact_reference(payload: dict[str, Any]) -> dict[str, Any]:
    artifact_format = payload.get("format")
    artifact_version = payload.get("version")
    if type(artifact_format) is not str or type(artifact_version) is not int:
        raise _fail("native source artifact format/version is unsupported")
    identity = (artifact_format, artifact_version)
    ir_fields = {
        (COLLECTION_FORMAT, COLLECTION_VERSION): "collection_ir_sha256",
        (REFERENCE_FORMAT, REFERENCE_VERSION): "catalog_sha256",
        (EXTERNAL_FORMAT, EXTERNAL_VERSION): "closure_ir_sha256",
        (EXTERNAL_FORMAT, EXTERNAL_EXECUTABLE_VERSION): "closure_ir_sha256",
        (GENBANK_FORMAT, GENBANK_VERSION): "bio_ir_sha256",
        (GFF3_FORMAT, GFF3_VERSION): "bio_ir_sha256",
    }
    field = ir_fields.get(identity)
    if field is None:
        raise _fail("native source artifact format/version is unsupported")
    try:
        artifact_sha = artifact_digest(payload)
    except (ContractError, KeyError) as failure:
        raise _fail(f"native source artifact digest is invalid: {failure}") from failure
    return {
        "format": identity[0],
        "version": identity[1],
        "artifact_sha256": artifact_sha,
        "ir_sha256": sha256(payload.get(field), f"native artifact {field}"),
    }


def _collection_records(payload: dict[str, Any]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for member in payload["members"]:
        sequence = member["sequence_artifact"]["sequence_ir"]
        records.append(
            {
                "record_id": sequence["record_id"],
                "bases": len(sequence["sequence"]),
                "sequence_sha256": sequence["canonical_sha256"],
                "refget_id": sequence["refget_id"],
            }
        )
    return records


def _genbank_records(payload: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "record_id": member["record_id"],
            "bases": member["sequence"]["bases"],
            "sequence_sha256": member["sequence"]["sha256"],
            "refget_id": member["sequence"]["refget_id"],
        }
        for member in payload["sequence_collection"]["members"]
    ]


def _descriptor_from_references(
    profile: str,
    parameters: dict[str, Any],
    input_references: dict[str, dict[str, Any]],
    artifacts: dict[str, dict[str, Any]],
    records: list[dict[str, Any]],
) -> dict[str, Any]:
    source_ir = {
        "profile": profile,
        "parameters": copy.deepcopy(parameters),
        "inputs": copy.deepcopy(input_references),
        "artifacts": {
            role: _artifact_reference(payload) for role, payload in artifacts.items()
        },
        "records": copy.deepcopy(records),
    }
    core = {
        "format": FORMAT,
        "version": VERSION,
        "producer": copy.deepcopy(PRODUCER),
        "source_ir": source_ir,
        "source_ir_sha256": digest(source_ir),
    }
    payload = {**core, "artifact_sha256": digest(core)}
    validate_source_descriptor(payload)
    _artifact_bytes(payload)
    return payload


def _descriptor(
    profile: str,
    parameters: dict[str, Any],
    inputs: dict[str, bytes],
    artifacts: dict[str, dict[str, Any]],
    records: list[dict[str, Any]],
) -> dict[str, Any]:
    return _descriptor_from_references(
        profile,
        parameters,
        {role: _input_reference(raw) for role, raw in inputs.items()},
        artifacts,
        records,
    )


def _artifact_bytes(payload: dict[str, Any]) -> bytes:
    identity = (payload.get("format"), payload.get("version"))
    if identity == (FORMAT, VERSION):
        maximum, depth, members = MAX_SOURCE_DESCRIPTOR_BYTES, MAX_JSON_DEPTH, MAX_JSON_MEMBERS
    elif identity == (GFF3_FORMAT, GFF3_VERSION):
        maximum, depth, members = GFF3_OUTPUT_BYTES, GFF3_JSON_DEPTH, GFF3_JSON_MEMBERS
    elif identity == (GENBANK_FORMAT, GENBANK_VERSION):
        maximum, depth, members = GENBANK_OUTPUT_BYTES, GENBANK_JSON_DEPTH, MAX_JSON_MEMBERS
    elif identity == (REFERENCE_FORMAT, REFERENCE_VERSION):
        maximum, depth, members = GFF3_OUTPUT_BYTES, MAX_JSON_DEPTH, MAX_JSON_MEMBERS
    elif identity in {
        (EXTERNAL_FORMAT, EXTERNAL_VERSION),
        (EXTERNAL_FORMAT, EXTERNAL_EXECUTABLE_VERSION),
    }:
        maximum, depth, members = EXTERNAL_OUTPUT_BYTES, MAX_JSON_DEPTH, MAX_JSON_MEMBERS
    else:
        maximum, depth, members = MAX_JSON_BYTES, MAX_JSON_DEPTH, MAX_JSON_MEMBERS
    try:
        return pretty_json_bytes(
            payload,
            ensure_ascii=False,
            maximum_bytes=maximum,
            maximum_depth=depth,
            maximum_members=members,
            maximum_string_bytes=MAX_STRING_BYTES,
        )
    except BoundedIOError as failure:
        raise _fail(f"artifact is outside output limits: {failure}") from failure


class SourceBundle:
    """Immutable-by-interface descriptor plus its exact native artifacts."""

    __slots__ = ("_descriptor", "_artifacts")

    def __init__(self, *_: Any, **__: Any) -> None:
        raise TypeError("SourceBundle instances are created by compile_source or load_source_bundle")

    @classmethod
    def _create(
        cls,
        descriptor: dict[str, Any],
        artifacts: dict[str, dict[str, Any]],
    ) -> SourceBundle:
        instance = object.__new__(cls)
        instance._descriptor = copy.deepcopy(descriptor)
        instance._artifacts = copy.deepcopy(artifacts)
        return instance

    @property
    def profile(self) -> str:
        return self._descriptor["source_ir"]["profile"]

    @property
    def artifacts(self) -> dict[str, dict[str, Any]]:
        return copy.deepcopy(self._artifacts)

    def to_dict(self) -> dict[str, Any]:
        return copy.deepcopy(self._descriptor)

    def save(self, directory: str | Path) -> dict[str, Path]:
        entries = {SOURCE_FILENAME: self._descriptor}
        entries.update(
            {
                ROLE_FILENAMES[role]: artifact
                for role, artifact in self._artifacts.items()
            }
        )
        try:
            return publish_directory(directory, entries, _artifact_bytes)
        except PublicationError as failure:
            state = " after commit" if failure.committed else ""
            raise _fail(
                f"cannot publish source bundle{state}: {failure.strerror}"
            ) from failure


def compile_source(
    profile: str,
    inputs: Mapping[str, str | Path],
    *,
    parameters: Mapping[str, Any],
    limits: ScaleLimits | None = None,
) -> SourceBundle:
    """Compile exactly one declared source profile; bytes never select a route."""

    route, normalized_parameters = _route(profile, parameters)
    if profile == EXTERNAL_PROFILE:
        raise _fail("external profiles must use compile_external_source")
    if limits is not None and type(limits) is not ScaleLimits:
        raise _fail("limits must be ScaleLimits")
    if limits is not None and profile != REFERENCE_FASTA_PROFILE:
        raise _fail("limits apply only to the reference FASTA profile")
    paths = _input_paths(inputs, route.input_roles)
    if profile == REFERENCE_FASTA_PROFILE:
        try:
            reference = compile_reference_fasta(
                paths["sequence"],
                wrapper=normalized_parameters["wrapper"],
                limits=DEFAULT_SCALE_LIMITS if limits is None else limits,
            )
        except SourceScaleError as failure:
            raise _fail(str(failure)) from failure
        artifacts = {"reference": reference}
        descriptor = _descriptor_from_references(
            profile,
            normalized_parameters,
            {"sequence": reference["inputs"]["source"]},
            artifacts,
            reference["sequence_catalog"]["records"],
        )
        return SourceBundle._create(descriptor, artifacts)

    snapshots = {
        role: _read(
            paths[role],
            MAX_GFF3_BYTES if role == "annotation" else (
                MAX_GENBANK_BYTES if role == "genbank" else MAX_INPUT_BYTES
            ),
            f"{role} source",
        )
        for role in route.input_roles
    }

    try:
        if profile in {RAW_PROFILE, FASTA_PROFILE}:
            sequence = _compile_sequence(profile, normalized_parameters, snapshots["sequence"])
            artifacts = {"sequence": sequence.to_dict()}
            records = _collection_records(artifacts["sequence"])
        elif profile == GENBANK_PROFILE:
            genbank = GenBankCompiler().compile_bytes(snapshots["genbank"]).to_dict()
            if genbank["profile"] != GENBANK_PROFILE:
                raise _fail("GenBank frontend emitted a different source profile")
            artifacts = {"genbank": genbank}
            records = _genbank_records(genbank)
        else:
            sequence_profile = normalized_parameters["sequence_profile"]
            sequence_parameters = normalized_parameters["sequence_parameters"]
            sequence = _compile_sequence(
                sequence_profile,
                sequence_parameters,
                snapshots["sequence"],
            )
            annotation = GFF3Compiler().compile_bytes(
                snapshots["annotation"], sequence
            ).to_dict()
            artifacts = {"sequence": sequence.to_dict(), "annotation": annotation}
            records = _collection_records(artifacts["sequence"])
    except (GenBankError, GFF3Error, SequenceCollectionError) as failure:
        raise _fail(str(failure)) from failure

    descriptor = _descriptor(
        profile,
        normalized_parameters,
        snapshots,
        artifacts,
        records,
    )
    return SourceBundle._create(descriptor, artifacts)


def compile_external_source(
    profile_manifest: dict[str, Any],
    source_descriptor: dict[str, Any],
    validation_report: dict[str, Any],
    *,
    original_payloads: Mapping[str, bytes],
    native_artifact_payloads: Mapping[str, bytes],
) -> SourceBundle:
    """Admit one explicitly selected, caller-attested external frontend."""

    try:
        external = build_external_source_closure(
            profile_manifest,
            source_descriptor,
            validation_report,
            original_payloads=original_payloads,
            native_artifact_payloads=native_artifact_payloads,
        )
    except ExternalSourceError as failure:
        raise _fail(str(failure)) from failure
    return _external_source_bundle(external)


def compile_external_source_paths(
    profile_manifest: dict[str, Any],
    source_descriptor: dict[str, Any],
    validation_report: dict[str, Any],
    *,
    original_paths: Mapping[str, str | Path],
    native_artifact_paths: Mapping[str, str | Path],
) -> SourceBundle:
    """Admit a caller-attested frontend while streaming original source paths."""

    try:
        external = build_external_source_closure_from_paths(
            profile_manifest,
            source_descriptor,
            validation_report,
            original_paths=original_paths,
            native_artifact_paths=native_artifact_paths,
        )
    except ExternalSourceError as failure:
        raise _fail(str(failure)) from failure
    return _external_source_bundle(external)


def compile_external_source_executable(
    profile_manifest: dict[str, Any],
    *,
    original_paths: Mapping[str, str | Path],
    frontend_executable: str | Path,
    validator_executable: str | Path,
) -> SourceBundle:
    """Compile any conforming DNA grammar through two pinned executables."""

    try:
        external = build_executable_external_source_closure(
            profile_manifest,
            original_paths=original_paths,
            frontend_executable=frontend_executable,
            validator_executable=validator_executable,
        )
    except ExternalSourceError as failure:
        raise _fail(str(failure)) from failure
    return _external_source_bundle(external)


def _external_source_bundle(external: dict[str, Any]) -> SourceBundle:
    try:
        parameters = external_profile_parameters(external)
        inputs = external_input_references(external)
        records = external_source_records(external)
    except ExternalSourceError as failure:
        raise _fail(str(failure)) from failure
    artifacts = {"external": external}
    descriptor = _descriptor_from_references(
        EXTERNAL_PROFILE,
        parameters,
        inputs,
        artifacts,
        records,
    )
    return SourceBundle._create(descriptor, artifacts)


def _validate_input_references(value: Any, roles: tuple[str, ...]) -> None:
    inputs = _keys(value, set(roles), "source_ir.inputs")
    for role in roles:
        item = _keys(
            inputs[role],
            {"sha256", "byte_length"},
            f"source_ir.inputs.{role}",
        )
        sha256(item["sha256"], f"source_ir.inputs.{role}.sha256")
        _integer(
            item["byte_length"],
            f"source_ir.inputs.{role}.byte_length",
        )


def _validate_external_input_references(value: Any) -> None:
    if type(value) is not dict or not value or len(value) > 64:
        raise _fail("source_ir.inputs must contain 1..64 external input roles")
    if any(type(role) is not str for role in value):
        raise _fail("source_ir.inputs external role names must be strings")
    for role, reference in value.items():
        _text(role, "source_ir.inputs external role", identifier=True)
        item = _keys(
            reference,
            {"sha256", "byte_length"},
            f"source_ir.inputs.{role}",
        )
        sha256(item["sha256"], f"source_ir.inputs.{role}.sha256")
        _integer(item["byte_length"], f"source_ir.inputs.{role}.byte_length")


def _validate_artifact_references(value: Any, roles: tuple[str, ...]) -> None:
    artifacts = _keys(value, set(roles), "source_ir.artifacts")
    for role in roles:
        item = _keys(
            artifacts[role],
            {"format", "version", "artifact_sha256", "ir_sha256"},
            f"source_ir.artifacts.{role}",
        )
        expected = {
            "sequence": (COLLECTION_FORMAT, COLLECTION_VERSION),
            "reference": (REFERENCE_FORMAT, REFERENCE_VERSION),
            "external": (EXTERNAL_FORMAT, EXTERNAL_VERSION),
            "genbank": (GENBANK_FORMAT, GENBANK_VERSION),
            "annotation": (GFF3_FORMAT, GFF3_VERSION),
        }[role]
        valid_identity = (
            item["format"] == expected[0]
            and type(item["version"]) is int
            and (
                item["version"] in {EXTERNAL_VERSION, EXTERNAL_EXECUTABLE_VERSION}
                if role == "external"
                else item["version"] == expected[1]
            )
        )
        if not valid_identity:
            raise _fail(f"source_ir.artifacts.{role} has the wrong format/version")
        sha256(item["artifact_sha256"], f"source_ir.artifacts.{role}.artifact_sha256")
        sha256(item["ir_sha256"], f"source_ir.artifacts.{role}.ir_sha256")


def validate_source_descriptor(payload: dict[str, Any]) -> dict[str, Any]:
    """Validate the closed descriptor without trusting native source content."""

    item = _keys(
        payload,
        {"format", "version", "producer", "source_ir", "source_ir_sha256", "artifact_sha256"},
        "source descriptor",
    )
    if item["format"] != FORMAT or type(item["version"]) is not int or item["version"] != VERSION:
        raise _fail("unsupported source descriptor format/version")
    if type(item["producer"]) is not dict or item["producer"] != PRODUCER:
        raise _fail("source descriptor producer is unsupported")
    source_ir = _keys(
        item["source_ir"],
        {"profile", "parameters", "inputs", "artifacts", "records"},
        "source_ir",
    )
    route, normalized_parameters = _route(source_ir["profile"], source_ir["parameters"])
    if source_ir["parameters"] != normalized_parameters:
        raise _fail("source_ir.parameters is not the exact normalized profile object")
    if source_ir["profile"] == EXTERNAL_PROFILE:
        _validate_external_input_references(source_ir["inputs"])
    else:
        _validate_input_references(source_ir["inputs"], route.input_roles)
    _validate_artifact_references(source_ir["artifacts"], route.artifact_roles)
    records = source_ir["records"]
    if type(records) is not list or not records or len(records) > MAX_SOURCE_RECORDS:
        raise _fail("source_ir.records must be a bounded nonempty array")
    seen: set[str] = set()
    for index, value in enumerate(records):
        record = _keys(
            value,
            {"record_id", "bases", "sequence_sha256", "refget_id"},
            f"source_ir.records[{index}]",
        )
        record_id = _text(record["record_id"], f"source_ir.records[{index}].record_id", identifier=True)
        if record_id in seen:
            raise _fail("source_ir.records contains duplicate record identifiers")
        seen.add(record_id)
        _integer(record["bases"], f"source_ir.records[{index}].bases", minimum=1)
        sha256(record["sequence_sha256"], f"source_ir.records[{index}].sequence_sha256")
        if type(record["refget_id"]) is not str or _REFGET_RE.fullmatch(record["refget_id"]) is None:
            raise _fail(f"source_ir.records[{index}].refget_id is invalid")
    if sha256(item["source_ir_sha256"], "source_ir_sha256") != digest(source_ir):
        raise _fail("source_ir_sha256 does not match source_ir")
    try:
        artifact_digest(item)
    except (ContractError, KeyError) as failure:
        raise _fail(f"source descriptor digest is invalid: {failure}") from failure
    _artifact_bytes(item)
    return copy.deepcopy(item)


def _read_descriptor_directory(directory: Path) -> dict[str, bytes]:
    limits = {
        SOURCE_FILENAME: MAX_SOURCE_DESCRIPTOR_BYTES,
        ROLE_FILENAMES["sequence"]: MAX_JSON_BYTES,
        ROLE_FILENAMES["reference"]: GFF3_OUTPUT_BYTES,
        ROLE_FILENAMES["external"]: EXTERNAL_OUTPUT_BYTES,
        ROLE_FILENAMES["genbank"]: GENBANK_OUTPUT_BYTES,
        ROLE_FILENAMES["annotation"]: GFF3_OUTPUT_BYTES,
    }
    inspected = directory.lstat()
    if stat.S_ISLNK(inspected.st_mode) or not stat.S_ISDIR(inspected.st_mode):
        raise _fail("source bundle must be a non-linked directory")
    if os.name != "posix":
        names = {entry.name for entry in directory.iterdir()}
        if not names <= limits.keys():
            raise _fail("source bundle contains an unknown filename")
        result: dict[str, bytes] = {}
        for name in names:
            child = directory / name
            before = child.lstat()
            if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode) or before.st_nlink != 1:
                raise _fail(f"source bundle child {name!r} is not one regular non-linked file")
            result[name] = _read(child, limits[name], f"source bundle child {name}")
            after = child.lstat()
            if (before.st_dev, before.st_ino, before.st_size) != (after.st_dev, after.st_ino, after.st_size):
                raise _fail(f"source bundle child {name!r} changed while reading")
        finished = directory.lstat()
        if (
            names != {entry.name for entry in directory.iterdir()}
            or (
                inspected.st_dev,
                inspected.st_ino,
                inspected.st_mtime_ns,
                inspected.st_ctime_ns,
            )
            != (
                finished.st_dev,
                finished.st_ino,
                finished.st_mtime_ns,
                finished.st_ctime_ns,
            )
        ):
            raise _fail("source bundle directory changed while reading")
        return result

    root_descriptor = -1
    child_descriptor = -1
    try:
        flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)
        root_descriptor = os.open(directory, flags)
        opened = os.fstat(root_descriptor)
        if not stat.S_ISDIR(opened.st_mode) or (opened.st_dev, opened.st_ino) != (inspected.st_dev, inspected.st_ino):
            raise _fail("source bundle directory changed while opening")
        names = set(os.listdir(root_descriptor))
        if not names <= limits.keys():
            raise _fail("source bundle contains an unknown filename")
        result = {}
        for name in names:
            before = os.stat(name, dir_fd=root_descriptor, follow_symlinks=False)
            if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode) or before.st_nlink != 1:
                raise _fail(f"source bundle child {name!r} is not one regular non-linked file")
            if before.st_size > limits[name]:
                raise _fail(f"source bundle child {name!r} exceeds its byte limit")
            child_descriptor = os.open(
                name,
                os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0),
                dir_fd=root_descriptor,
            )
            after_open = os.fstat(child_descriptor)
            identity = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns, before.st_ctime_ns)
            if (
                not stat.S_ISREG(after_open.st_mode)
                or after_open.st_nlink != 1
                or identity != (after_open.st_dev, after_open.st_ino, after_open.st_size, after_open.st_mtime_ns, after_open.st_ctime_ns)
            ):
                raise _fail(f"source bundle child {name!r} changed while opening")
            chunks: list[bytes] = []
            remaining = limits[name] + 1
            while remaining:
                chunk = os.read(child_descriptor, min(64 * 1024, remaining))
                if not chunk:
                    break
                chunks.append(chunk)
                remaining -= len(chunk)
            raw = b"".join(chunks)
            if len(raw) > limits[name] or len(raw) != before.st_size:
                raise _fail(f"source bundle child {name!r} changed or exceeds its byte limit")
            finished_child = os.fstat(child_descriptor)
            if identity != (finished_child.st_dev, finished_child.st_ino, finished_child.st_size, finished_child.st_mtime_ns, finished_child.st_ctime_ns):
                raise _fail(f"source bundle child {name!r} changed while reading")
            os.close(child_descriptor)
            child_descriptor = -1
            result[name] = raw
        finished_descriptor = os.fstat(root_descriptor)
        finished = directory.lstat()
        opened_identity = (
            opened.st_dev,
            opened.st_ino,
            opened.st_mtime_ns,
            opened.st_ctime_ns,
        )
        if (
            names != set(os.listdir(root_descriptor))
            or opened_identity
            != (
                finished_descriptor.st_dev,
                finished_descriptor.st_ino,
                finished_descriptor.st_mtime_ns,
                finished_descriptor.st_ctime_ns,
            )
            or (opened.st_dev, opened.st_ino) != (finished.st_dev, finished.st_ino)
        ):
            raise _fail("source bundle directory changed while reading")
        return result
    except SourceError:
        raise
    except OSError as failure:
        raise _fail(f"cannot read source bundle: {failure}") from failure
    finally:
        if child_descriptor >= 0:
            os.close(child_descriptor)
        if root_descriptor >= 0:
            os.close(root_descriptor)


def _loads(raw: bytes, label: str, *, maximum: int, depth: int, members: int) -> dict[str, Any]:
    try:
        return load_json_object(
            raw,
            label,
            maximum_bytes=maximum,
            maximum_depth=depth,
            maximum_members=members,
            maximum_string_bytes=MAX_STRING_BYTES,
        )
    except BoundedIOError as failure:
        raise _fail(f"invalid {label}: {failure}") from failure


def _validate_native(
    descriptor: dict[str, Any],
    raw_files: dict[str, bytes],
) -> dict[str, dict[str, Any]]:
    profile = descriptor["source_ir"]["profile"]
    route, _ = _route(profile, descriptor["source_ir"]["parameters"])
    expected_files = {SOURCE_FILENAME} | {ROLE_FILENAMES[role] for role in route.artifact_roles}
    if set(raw_files) != expected_files:
        raise _fail("source bundle filename closure does not match its profile")
    artifacts: dict[str, dict[str, Any]] = {}
    for role in route.artifact_roles:
        filename = ROLE_FILENAMES[role]
        if role == "annotation":
            payload = _loads(
                raw_files[filename], filename,
                maximum=GFF3_OUTPUT_BYTES, depth=GFF3_JSON_DEPTH, members=GFF3_JSON_MEMBERS,
            )
        elif role == "genbank":
            payload = _loads(
                raw_files[filename], filename,
                maximum=GENBANK_OUTPUT_BYTES, depth=GENBANK_JSON_DEPTH, members=MAX_JSON_MEMBERS,
            )
        elif role == "reference":
            payload = _loads(
                raw_files[filename], filename,
                maximum=GFF3_OUTPUT_BYTES, depth=MAX_JSON_DEPTH, members=MAX_JSON_MEMBERS,
            )
        elif role == "external":
            payload = _loads(
                raw_files[filename], filename,
                maximum=EXTERNAL_OUTPUT_BYTES, depth=MAX_JSON_DEPTH, members=MAX_JSON_MEMBERS,
            )
        else:
            payload = _loads(
                raw_files[filename], filename,
                maximum=MAX_JSON_BYTES, depth=MAX_JSON_DEPTH, members=MAX_JSON_MEMBERS,
            )
        artifacts[role] = payload

    return _validate_native_payloads(descriptor, artifacts)


def _validate_native_payloads(
    descriptor: dict[str, Any],
    artifacts: dict[str, dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    profile = descriptor["source_ir"]["profile"]
    route, _ = _route(profile, descriptor["source_ir"]["parameters"])
    if set(artifacts) != set(route.artifact_roles):
        raise _fail("native artifact role closure does not match its profile")
    collection: SequenceCollectionArtifact | None = None
    for role in route.artifact_roles:
        payload = artifacts[role]
        try:
            if role == "sequence":
                collection = _as_collection(payload)
            elif role == "reference":
                validate_reference_catalog(payload)
            elif role == "external":
                validate_external_source_closure(payload)
            elif role == "genbank":
                validate_genbank_artifact(payload)
            else:
                assert collection is not None
                validate_gff3_artifact(payload, collection)
        except (
            SequenceCollectionError,
            GenBankError,
            GFF3Error,
            ExternalSourceError,
        ) as failure:
            raise _fail(f"invalid {role} native artifact: {failure}") from failure

    expected_refs = {role: _artifact_reference(value) for role, value in artifacts.items()}
    if descriptor["source_ir"]["artifacts"] != expected_refs:
        raise _fail("source descriptor native artifact references do not match bundle children")
    if profile == GENBANK_PROFILE:
        records = _genbank_records(artifacts["genbank"])
    elif profile == REFERENCE_FASTA_PROFILE:
        records = artifacts["reference"]["sequence_catalog"]["records"]
    elif profile == EXTERNAL_PROFILE:
        records = external_source_records(artifacts["external"])
    else:
        records = _collection_records(artifacts["sequence"])
    if descriptor["source_ir"]["records"] != records:
        raise _fail("source descriptor record catalog does not match native source")

    descriptor_inputs = descriptor["source_ir"]["inputs"]
    if profile == GENBANK_PROFILE:
        native_source = artifacts["genbank"]["sequence_collection"]["inputs"]["source"]
        expected_input_digests = {"genbank": native_source["sha256"]}
        if descriptor_inputs["genbank"]["byte_length"] != native_source["byte_length"]:
            raise _fail(
                "source descriptor input byte length does not match native GenBank source"
            )
    elif profile == REFERENCE_FASTA_PROFILE:
        native_source = artifacts["reference"]["inputs"]["source"]
        if descriptor_inputs != {"sequence": native_source}:
            raise _fail(
                "source descriptor input identity does not match reference FASTA source"
            )
        expected_input_digests = {"sequence": native_source["sha256"]}
    elif profile == EXTERNAL_PROFILE:
        expected_inputs = external_input_references(artifacts["external"])
        if descriptor_inputs != expected_inputs:
            raise _fail(
                "source descriptor inputs do not match external frontend evidence"
            )
        expected_input_digests = {
            role: reference["sha256"]
            for role, reference in expected_inputs.items()
        }
    else:
        expected_input_digests = {
            "sequence": artifacts["sequence"]["inputs"]["raw_sha256"]
        }
        if profile == GFF3_PROFILE:
            expected_input_digests["annotation"] = artifacts["annotation"]["inputs"]["gff3_sha256"]
    if {
        role: value["sha256"] for role, value in descriptor_inputs.items()
    } != expected_input_digests:
        raise _fail("source descriptor input digests do not match native source inputs")

    if profile == RAW_PROFILE:
        sequence_inputs = artifacts["sequence"]["inputs"]
        if sequence_inputs["kind"] != "raw-iupac" or sequence_inputs["wrapper"] is not None:
            raise _fail("raw profile does not match its native collection")
    elif profile == FASTA_PROFILE:
        expected_wrapper = None if descriptor["source_ir"]["parameters"]["wrapper"] == "identity" else "gzip"
        sequence_inputs = artifacts["sequence"]["inputs"]
        if sequence_inputs["kind"] != "fasta" or sequence_inputs["wrapper"] != expected_wrapper:
            raise _fail("FASTA profile does not match its native collection")
    elif profile == REFERENCE_FASTA_PROFILE:
        if (
            artifacts["reference"]["inputs"]["profile"] != REFERENCE_FASTA_PROFILE
            or artifacts["reference"]["inputs"]["wrapper"]
            != descriptor["source_ir"]["parameters"]["wrapper"]
        ):
            raise _fail("reference FASTA profile does not match its native catalog")
    elif profile == EXTERNAL_PROFILE:
        if descriptor["source_ir"]["parameters"] != external_profile_parameters(
            artifacts["external"]
        ):
            raise _fail("external source profile does not match its evidence closure")
    elif profile == GENBANK_PROFILE and artifacts["genbank"]["profile"] != GENBANK_PROFILE:
        raise _fail("GenBank profile does not match its native artifact")
    elif profile == GFF3_PROFILE:
        nested = descriptor["source_ir"]["parameters"]
        sequence_inputs = artifacts["sequence"]["inputs"]
        expected_kind = "raw-iupac" if nested["sequence_profile"] == RAW_PROFILE else "fasta"
        expected_wrapper = (
            None
            if nested["sequence_profile"] == RAW_PROFILE or nested["sequence_parameters"]["wrapper"] == "identity"
            else "gzip"
        )
        if sequence_inputs["kind"] != expected_kind or sequence_inputs["wrapper"] != expected_wrapper:
            raise _fail("GFF3 sequence profile does not match its native collection")
    return copy.deepcopy(artifacts)


def validate_source_bundle(value: SourceBundle) -> SourceBundle:
    if not isinstance(value, SourceBundle):
        raise _fail("source bundle must be a SourceBundle")
    descriptor = validate_source_descriptor(value.to_dict())
    artifacts = _validate_native_payloads(descriptor, value.artifacts)
    return SourceBundle._create(descriptor, artifacts)


def load_source_bundle(directory: str | Path) -> SourceBundle:
    root = Path(directory)
    try:
        raw_files = _read_descriptor_directory(root)
    except FileNotFoundError as failure:
        raise _fail(f"source bundle does not exist: {failure}") from failure
    if SOURCE_FILENAME not in raw_files:
        raise _fail("source bundle is missing source.json")
    descriptor = validate_source_descriptor(
        _loads(
            raw_files[SOURCE_FILENAME], SOURCE_FILENAME,
            maximum=MAX_SOURCE_DESCRIPTOR_BYTES,
            depth=MAX_JSON_DEPTH,
            members=MAX_JSON_MEMBERS,
        )
    )
    artifacts = _validate_native(descriptor, raw_files)
    return SourceBundle._create(descriptor, artifacts)


__all__ = [
    "FASTA_PROFILE",
    "EXTERNAL_PROFILE",
    "FORMAT",
    "GFF3_PROFILE",
    "GENBANK_PROFILE",
    "PROFILES",
    "RAW_PROFILE",
    "REFERENCE_FASTA_PROFILE",
    "ROLE_FILENAMES",
    "SOURCE_FILENAME",
    "ScaleLimits",
    "SourceBundle",
    "SourceError",
    "VERSION",
    "compile_source",
    "compile_external_source",
    "compile_external_source_executable",
    "compile_external_source_paths",
    "load_source_bundle",
    "validate_source_bundle",
    "validate_source_descriptor",
]

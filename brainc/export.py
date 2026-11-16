"""Deterministic RO-Crate 1.1 and COMBINE/OMEX archive exporters for sealed Epigenesis bundles."""

from __future__ import annotations

import base64
import hashlib
import io
import os
from pathlib import Path
import secrets
import shutil
import stat
from typing import Any
from xml.sax.saxutils import escape as xml_escape
import zipfile

from ._canonical import ContractError, digest
from ._io import (
    BoundedIOError,
    MAX_JSON_BYTES,
    atomic_write_file,
    load_json_object,
    pretty_json_bytes,
    read_regular_file,
)
from ._publish import HAS_POSIX_DIRECTORY_PUBLICATION, rename_directory_noreplace
from .source import (
    EXTERNAL_PROFILE,
    FASTA_PROFILE,
    GENBANK_PROFILE,
    GFF3_PROFILE,
    RAW_PROFILE,
    REFERENCE_FASTA_PROFILE,
    validate_source_descriptor,
)
from .validator_development import (
    BUNDLE_FILENAME,
    CHILD_FILENAMES,
    MAX_GFF3_ARTIFACT_BYTES,
    MAX_SOURCE_DESCRIPTOR_BYTES,
    SOURCE_FILENAME,
    SOURCE_FILENAMES,
)
from .validator_variant import VariantValidationError, validate_variant_collection


RO_CRATE_CONTEXT = "https://w3id.org/ro/crate/1.1/context"
RO_CRATE_SPEC = "https://w3id.org/ro/crate/1.1"
OMEX_SPEC = "http://identifiers.org/combine.specifications/omex"
OMEX_MANIFEST_NS = "http://identifiers.org/combine.specifications/omex-manifest"
EPIGENESIS_SCHEMA_BY_FILE = {
    "development/bundle.json": "https://epigenesis.org/schemas/brainc.development-bundle.v1.schema.json",
    "source/source.json": "https://epigenesis.org/schemas/brainc.source-descriptor.v1.schema.json",
    "variants.json": "https://epigenesis.org/schemas/brainc.variant-collection.v1.schema.json",
}

_DEV_CHILD_META: dict[str, tuple[str, str, int]] = {
    "compilation_record": ("compilation_record.json", "brainc.development-compilation", 1),
    "development_module": ("development_module.json", "brainc.development-module", 1),
    "lowering_policy": ("lowering_policy.json", "brainc.lowering-policy", 2),
    "prediction_request": ("prediction_request.json", "brainc.prediction-request", 2),
    "prediction_response": ("prediction_response.json", "brainc.prediction-response", 2),
    "provider_manifest": ("provider_manifest.json", "brainc.provider-manifest", 2),
    "target_contract": ("target_contract.json", "brainc.target-contract", 1),
}

_SOURCE_EXPECTED_CLOSURE: dict[str, set[str]] = {
    RAW_PROFILE: {"source.json", "sequence.json"},
    FASTA_PROFILE: {"source.json", "sequence.json"},
    GFF3_PROFILE: {"source.json", "sequence.json", "annotation.json"},
    GENBANK_PROFILE: {"source.json", "genbank.json"},
    REFERENCE_FASTA_PROFILE: {"source.json", "reference.json"},
    EXTERNAL_PROFILE: {"source.json", "external.json"},
}

_HAS_DIR_FD = (
    os.name == "posix"
    and os.listdir in getattr(os, "supports_fd", set())
    and os.open in getattr(os, "supports_dir_fd", set())
    and getattr(os, "O_DIRECTORY", 0) != 0
    and getattr(os, "O_NOFOLLOW", 0) != 0
)


class ExportError(ContractError):
    """An RO-Crate or COMBINE/OMEX export failed validation or I/O safety checks."""


def _fail(detail: str) -> ExportError:
    return ExportError(f"EXPORT001: {detail}")


def _stat_identity(info: os.stat_result) -> tuple[Any, ...]:
    return (
        info.st_dev,
        info.st_ino,
        info.st_size,
        info.st_nlink,
        getattr(info, "st_mtime_ns", int(info.st_mtime * 1e9)),
        getattr(info, "st_ctime_ns", int(info.st_ctime * 1e9)),
    )


def _read_dir_snapshot(
    directory: Path,
    allowed_limits: dict[str, int],
    *,
    required: set[str],
    label: str,
) -> dict[str, bytes]:
    if not _HAS_DIR_FD:
        if directory.is_symlink() or not directory.is_dir():
            raise _fail(f"{label} must be a non-symlink directory")
        names = [entry.name for entry in directory.iterdir()]
        if not required <= set(names) or not set(names) <= set(allowed_limits):
            raise _fail(f"unexpected files in {label}: {sorted(names)}")
        return {
            name: read_regular_file(directory / name, maximum_bytes=allowed_limits[name], label=f"{label}/{name}")
            for name in sorted(names)
        }

    dir_fd = -1
    try:
        dir_fd = os.open(directory, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
        names = os.listdir(dir_fd)
        if not required <= set(names) or not set(names) <= set(allowed_limits):
            raise _fail(f"unexpected files in {label}: {sorted(names)}")
        flags = os.O_RDONLY | os.O_NOFOLLOW | getattr(os, "O_NONBLOCK", 0)
        out: dict[str, bytes] = {}
        for name in sorted(names):
            maximum = allowed_limits[name]
            cfd = os.open(name, flags, dir_fd=dir_fd)
            try:
                info = os.fstat(cfd)
                if not stat.S_ISREG(info.st_mode):
                    raise _fail(f"{label}/{name} is not a regular file")
                if info.st_nlink != 1:
                    raise _fail(f"{label}/{name} must have link count 1 (found {info.st_nlink})")
                if info.st_size > maximum:
                    raise _fail(f"{label}/{name} exceeds maximum byte size {maximum}")
                chunks: list[bytes] = []
                total = 0
                while True:
                    block = os.read(cfd, 65536)
                    if not block:
                        break
                    total += len(block)
                    if total > maximum:
                        raise _fail(f"{label}/{name} exceeded maximum byte size during read")
                    chunks.append(block)
                after = os.fstat(cfd)
                if _stat_identity(after) != _stat_identity(info) or total != info.st_size:
                    raise _fail(f"{label}/{name} changed during read")
                out[name] = b"".join(chunks)
            finally:
                os.close(cfd)
        return out
    except ExportError:
        raise
    except (OSError, BoundedIOError) as failure:
        raise _fail(f"cannot read {label}: {failure}") from failure
    finally:
        if dir_fd >= 0:
            os.close(dir_fd)


def _verify_seal(obj: Any, expected_format: str, expected_version: int, label: str) -> str:
    if type(obj) is not dict:
        raise _fail(f"{label} must be a JSON object")
    if obj.get("format") != expected_format or type(obj.get("version")) is not int or obj.get("version") != expected_version:
        raise _fail(f"{label} format or version mismatch")
    stored = obj.get("artifact_sha256")
    if type(stored) is not str:
        raise _fail(f"{label} missing artifact_sha256")
    core = {k: v for k, v in obj.items() if k != "artifact_sha256"}
    if digest(core) != stored:
        raise _fail(f"{label} artifact_sha256 seal mismatch")
    return stored


def _verify_development_snapshot(dev_bytes: dict[str, bytes]) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    bundle_obj = load_json_object(dev_bytes[BUNDLE_FILENAME], BUNDLE_FILENAME)
    _verify_seal(bundle_obj, "brainc.development-bundle", 1, BUNDLE_FILENAME)
    if bundle_obj.get("profile") != "development-module/v1":
        raise _fail("bundle.json profile must be 'development-module/v1'")
    if type(bundle_obj.get("source")) is not dict or type(bundle_obj["source"].get("descriptor")) is not dict:
        raise _fail("bundle.json missing source.descriptor")
    if type(bundle_obj.get("artifacts")) is not dict:
        raise _fail("bundle.json missing artifacts map")
    if set(bundle_obj["artifacts"]) != set(_DEV_CHILD_META):
        raise _fail("bundle.json artifacts roles mismatch")

    children: dict[str, dict[str, Any]] = {}
    child_shas: dict[str, str] = {}
    for role, (fname, exp_fmt, exp_ver) in _DEV_CHILD_META.items():
        obj = load_json_object(dev_bytes[fname], fname)
        sha = _verify_seal(obj, exp_fmt, exp_ver, fname)
        ref = bundle_obj["artifacts"].get(role)
        if (
            type(ref) is not dict
            or ref.get("format") != exp_fmt
            or ref.get("version") != exp_ver
            or ref.get("artifact_sha256") != sha
        ):
            raise _fail(f"bundle.json artifacts.{role} disagrees with {fname}")
        children[role] = obj
        child_shas[role] = sha

    # Verify development_module internal module_sha256 and inline-base64 digests
    dev_mod = children["development_module"]
    module_inner = dev_mod.get("module")
    if type(module_inner) is not dict or dev_mod.get("module_sha256") != digest(module_inner):
        raise _fail("development_module.module_sha256 mismatch")
    tensors = module_inner.get("tensors")
    if type(tensors) is not list:
        raise _fail("development_module.module.tensors must be an array")
    for idx, tensor in enumerate(tensors):
        if type(tensor) is not dict or type(tensor.get("storage")) is not dict:
            raise _fail(f"development_module tensor[{idx}] missing storage")
        storage = tensor["storage"]
        if storage.get("kind") == "inline-base64":
            data_b64 = storage.get("data")
            if type(data_b64) is not str:
                raise _fail(f"tensor[{idx}] inline-base64 missing data")
            try:
                raw_tensor = base64.b64decode(data_b64, validate=True)
            except ValueError as failure:
                raise _fail(f"tensor[{idx}] invalid base64: {failure}") from failure
            if (
                len(raw_tensor) != storage.get("byte_length")
                or hashlib.sha256(raw_tensor).hexdigest() != storage.get("sha256")
            ):
                raise _fail(f"tensor[{idx}] inline-base64 digest or byte_length mismatch")

    # Verify compilation_record and request/response/policy/target bindings
    comp = children["compilation_record"]
    src_desc_ref = bundle_obj["source"]["descriptor"]
    if comp.get("source") != src_desc_ref:
        raise _fail("compilation_record.source disagrees with bundle.json")
    comp_artifacts = comp.get("artifacts")
    if type(comp_artifacts) is not dict:
        raise _fail("compilation_record missing artifacts")
    for role in (
        "provider_manifest",
        "prediction_request",
        "prediction_response",
        "lowering_policy",
        "target_contract",
        "development_module",
    ):
        cref = comp_artifacts.get(role)
        if type(cref) is not dict or cref.get("artifact_sha256") != child_shas[role]:
            raise _fail(f"compilation_record.artifacts.{role} mismatch")
    if type(comp.get("result")) is not dict or comp["result"].get("module_sha256") != dev_mod["module_sha256"]:
        raise _fail("compilation_record.result.module_sha256 mismatch")

    req = children["prediction_request"]
    resp = children["prediction_response"]
    pol = children["lowering_policy"]
    tgt = children["target_contract"]
    if req.get("provider_manifest_sha256") != child_shas["provider_manifest"]:
        raise _fail("prediction_request.provider_manifest_sha256 mismatch")
    if req.get("source") != src_desc_ref:
        raise _fail("prediction_request.source disagrees with bundle.json")
    if resp.get("request_artifact_sha256") != child_shas["prediction_request"]:
        raise _fail("prediction_response.request_artifact_sha256 mismatch")
    if type(tgt.get("contract")) is not dict or type(pol.get("target")) is not dict:
        raise _fail("target_contract or lowering_policy missing contract/target")
    if pol["target"].get("contract_sha256") != digest(tgt["contract"]):
        raise _fail("lowering_policy.target.contract_sha256 disagrees with target_contract")

    return bundle_obj, children


def _verify_source_snapshot(src_bytes: dict[str, bytes]) -> tuple[dict[str, Any], dict[str, Any]]:
    src_desc = load_json_object(src_bytes[SOURCE_FILENAME], SOURCE_FILENAME)
    validate_source_descriptor(src_desc)
    source_ir = src_desc["source_ir"]
    profile = source_ir["profile"]
    if profile not in _SOURCE_EXPECTED_CLOSURE or set(src_bytes) != _SOURCE_EXPECTED_CLOSURE[profile]:
        raise _fail(f"source bundle files {sorted(src_bytes)} do not match profile {profile!r}")
    artifacts_meta = source_ir["artifacts"]
    expected_child_meta = {
        "sequence": ("sequence.json", "brain01.sequence-collection-ir", 1, "collection_ir_sha256"),
        "reference": ("reference.json", "brainc.reference-sequence-catalog", 1, "catalog_sha256"),
        "genbank": ("genbank.json", "brainc.bio.insdc-genbank-ir", 2, "bio_ir_sha256"),
        "annotation": ("annotation.json", "brainc.bio.gff3-ir", 1, "bio_ir_sha256"),
        "external": ("external.json", "brainc.external-source-closure", 1, "closure_ir_sha256"),
    }
    for role, meta in artifacts_meta.items():
        if role not in expected_child_meta:
            raise _fail(f"unexpected source artifact role: {role!r}")
        fname, exp_fmt, exp_ver, ir_key = expected_child_meta[role]
        obj = load_json_object(src_bytes[fname], fname)
        sha = _verify_seal(obj, exp_fmt, exp_ver, fname)
        if meta.get("artifact_sha256") != sha or obj.get(ir_key) != meta.get("ir_sha256"):
            raise _fail(f"source_ir.artifacts.{role} disagrees with {fname}")
    return src_desc, artifacts_meta


def _collect_export_entries(
    development_bundle: str | Path,
    *,
    source_bundle: str | Path | None = None,
    variant_collection: str | Path | None = None,
) -> tuple[dict[str, bytes], dict[str, Any]]:
    dev_dir = Path(development_bundle)
    dev_limits = {BUNDLE_FILENAME: MAX_JSON_BYTES, **{fname: MAX_JSON_BYTES for fname in CHILD_FILENAMES.values()}}
    try:
        dev_bytes = _read_dir_snapshot(
            dev_dir,
            dev_limits,
            required=set(dev_limits),
            label="development bundle",
        )
        bundle_obj, dev_children = _verify_development_snapshot(dev_bytes)
    except ExportError:
        raise
    except (ContractError, BoundedIOError, KeyError, TypeError, ValueError) as failure:
        raise _fail(f"invalid development bundle: {failure}") from failure

    entries: dict[str, bytes] = {f"development/{name}": raw for name, raw in dev_bytes.items()}
    bound_source_sha = bundle_obj["source"]["descriptor"]["artifact_sha256"]
    bound_source_ir_sha = bundle_obj["source"]["descriptor"].get("ir_sha256")

    source_artifact_sha: str | None = None
    src_bytes: dict[str, bytes] | None = None
    if source_bundle is not None:
        src_dir = Path(source_bundle)
        src_limits = {
            SOURCE_FILENAME: MAX_SOURCE_DESCRIPTOR_BYTES,
            **{
                fname: (MAX_GFF3_ARTIFACT_BYTES if fname == "annotation.json" else MAX_SOURCE_DESCRIPTOR_BYTES)
                for fname in SOURCE_FILENAMES.values()
            },
        }
        try:
            src_bytes = _read_dir_snapshot(
                src_dir,
                src_limits,
                required={SOURCE_FILENAME},
                label="source bundle",
            )
            src_desc, src_artifacts_meta = _verify_source_snapshot(src_bytes)
            source_artifact_sha = src_desc["artifact_sha256"]
            if source_artifact_sha != bound_source_sha or src_desc["source_ir_sha256"] != bound_source_ir_sha:
                raise _fail("development bundle does not bind the supplied source bundle")
            if bundle_obj["source"].get("native") != src_artifacts_meta:
                raise _fail("development bundle source.native disagrees with source bundle")
            for name, raw in src_bytes.items():
                entries[f"source/{name}"] = raw
        except ExportError:
            raise
        except (ContractError, BoundedIOError, KeyError, TypeError, ValueError) as failure:
            raise _fail(f"invalid source bundle: {failure}") from failure

    variant_artifact_sha: str | None = None
    if variant_collection is not None:
        var_path = Path(variant_collection)
        try:
            var_raw = read_regular_file(var_path, maximum_bytes=MAX_JSON_BYTES, label="variants.json")
            var_obj = load_json_object(var_raw, "variants.json")
            variant_artifact_sha = _verify_seal(var_obj, "brainc.variant-collection", 1, "variants.json")
            if type(var_obj.get("variants_ir")) is not dict or var_obj.get("variants_ir_sha256") != digest(var_obj["variants_ir"]):
                raise _fail("variants.json variants_ir_sha256 mismatch")
            var_src = var_obj.get("source")
            if type(var_src) is not dict or var_src.get("artifact_sha256") != bound_source_sha:
                raise _fail("variant collection source does not match development bundle source")
            if src_bytes is not None:
                profile = load_json_object(src_bytes[SOURCE_FILENAME], SOURCE_FILENAME)["source_ir"]["profile"]
                if profile in {RAW_PROFILE, FASTA_PROFILE, GFF3_PROFILE, GENBANK_PROFILE}:
                    rep = validate_variant_collection(src_bytes, var_obj)
                    if not rep["valid"]:
                        raise _fail(f"variant collection failed replay: {rep.get('detail')}")
            entries["variants.json"] = var_raw
        except ExportError:
            raise
        except (ContractError, BoundedIOError, VariantValidationError, KeyError, TypeError, ValueError) as failure:
            raise _fail(f"invalid variant collection: {failure}") from failure

    file_entities: list[dict[str, Any]] = []
    part_refs: list[dict[str, str]] = []
    for rel_path in sorted(entries):
        raw = entries[rel_path]
        part_refs.append({"@id": rel_path})
        entity: dict[str, Any] = {
            "@id": rel_path,
            "@type": "File",
            "name": rel_path,
            "encodingFormat": "application/json",
            "contentSize": len(raw),
            "sha256": hashlib.sha256(raw).hexdigest(),
        }
        if rel_path in EPIGENESIS_SCHEMA_BY_FILE:
            entity["conformsTo"] = {"@id": EPIGENESIS_SCHEMA_BY_FILE[rel_path]}
        file_entities.append(entity)

    root_dataset: dict[str, Any] = {
        "@id": "./",
        "@type": "Dataset",
        "name": "Epigenesis Sealed Provenance Bundle",
        "description": "Deterministic sealed provenance bundle for Development Module IR, source descriptor, and optional VRS 2.x variant collection.",
        "license": "https://spdx.org/licenses/Apache-2.0",
        "hasPart": part_refs,
        "epigenesisBundleSha256": bundle_obj["artifact_sha256"],
        "epigenesisModuleSha256": dev_children["development_module"]["artifact_sha256"],
    }
    if source_artifact_sha is not None:
        root_dataset["epigenesisSourceSha256"] = source_artifact_sha
    if variant_artifact_sha is not None:
        root_dataset["epigenesisVariantCollectionSha256"] = variant_artifact_sha

    crate_metadata = {
        "@context": RO_CRATE_CONTEXT,
        "@graph": [
            {
                "@id": "ro-crate-metadata.json",
                "@type": "CreativeWork",
                "conformsTo": {"@id": RO_CRATE_SPEC},
                "about": {"@id": "./"},
            },
            root_dataset,
            *file_entities,
        ],
    }
    meta_bytes = pretty_json_bytes(
        crate_metadata,
        ensure_ascii=False,
        maximum_bytes=MAX_JSON_BYTES,
        maximum_depth=64,
        maximum_members=1_000_000,
        maximum_string_bytes=1_048_576,
    )
    entries["ro-crate-metadata.json"] = meta_bytes
    summary = {
        "format": "epigenesis.export-summary",
        "version": 1,
        "rocrate_metadata_sha256": hashlib.sha256(meta_bytes).hexdigest(),
        "rocrate_canonical_sha256": digest(crate_metadata),
        "development_bundle_artifact_sha256": bundle_obj["artifact_sha256"],
        "source_artifact_sha256": source_artifact_sha,
        "variant_artifact_sha256": variant_artifact_sha,
        "files": {k: hashlib.sha256(v).hexdigest() for k, v in sorted(entries.items())},
    }
    return entries, summary


def export_rocrate(
    development_bundle: str | Path,
    output_dir: str | Path,
    *,
    source_bundle: str | Path | None = None,
    variant_collection: str | Path | None = None,
) -> dict[str, Any]:
    """Export a verified bundle closure as a deterministic RO-Crate 1.1 directory."""
    entries, summary = _collect_export_entries(
        development_bundle,
        source_bundle=source_bundle,
        variant_collection=variant_collection,
    )
    target = Path(output_dir)
    if not target.name or target.name in {".", ".."}:
        raise _fail(f"invalid RO-Crate output directory name: {target}")
    parent = target.parent
    if not parent.is_dir():
        raise _fail(f"RO-Crate output parent directory does not exist: {parent}")

    staging = parent / f".{target.name}.rocrate.{secrets.token_hex(8)}.tmp"
    pfd = -1
    committed = False
    try:
        if HAS_POSIX_DIRECTORY_PUBLICATION:
            pfd = os.open(parent, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
            os.mkdir(staging.name, 0o755, dir_fd=pfd)
        else:
            os.mkdir(staging, mode=0o755)
        for rel_path, raw in entries.items():
            dest = staging / rel_path
            dest.parent.mkdir(parents=True, exist_ok=True)
            atomic_write_file(dest, raw, maximum_bytes=MAX_SOURCE_DESCRIPTOR_BYTES, label=rel_path)
        if HAS_POSIX_DIRECTORY_PUBLICATION:
            rename_directory_noreplace(staging, target, parent_descriptor=pfd)
        else:
            if target.exists() or target.is_symlink():
                raise _fail(f"RO-Crate output directory already exists: {target}")
            staging.rename(target)
        committed = True
    except ExportError:
        raise
    except (OSError, BoundedIOError) as failure:
        raise _fail(f"cannot write RO-Crate directory: {failure}") from failure
    finally:
        if pfd >= 0:
            os.close(pfd)
        if not committed and staging.exists():
            shutil.rmtree(staging, ignore_errors=True)

    return {**summary, "mode": "ro-crate/1.1", "output": str(target)}


def _combine_manifest_xml(entries: dict[str, bytes]) -> bytes:
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<omexManifest xmlns="{OMEX_MANIFEST_NS}">',
        f'  <content location="." format="{OMEX_SPEC}"/>',
        f'  <content location="./manifest.xml" format="{OMEX_MANIFEST_NS}"/>',
    ]
    for rel_path in sorted(entries):
        escaped = xml_escape(f"./{rel_path}", {'"': "&quot;"})
        lines.append(f'  <content location="{escaped}" format="application/json"/>')
    lines.append("</omexManifest>\n")
    return "\n".join(lines).encode("utf-8")


def export_omex(
    development_bundle: str | Path,
    output_archive: str | Path,
    *,
    source_bundle: str | Path | None = None,
    variant_collection: str | Path | None = None,
) -> dict[str, Any]:
    """Export a verified bundle closure as a deterministic COMBINE/OMEX archive (.omex)."""
    entries, summary = _collect_export_entries(
        development_bundle,
        source_bundle=source_bundle,
        variant_collection=variant_collection,
    )
    manifest_bytes = _combine_manifest_xml(entries)
    archive_entries = {"manifest.xml": manifest_bytes, **entries}

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, mode="w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for rel_path in sorted(archive_entries):
            info = zipfile.ZipInfo(filename=rel_path, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            zf.writestr(info, archive_entries[rel_path])

    archive_bytes = buf.getvalue()
    try:
        atomic_write_file(
            output_archive,
            archive_bytes,
            maximum_bytes=MAX_SOURCE_DESCRIPTOR_BYTES,
            label="OMEX archive",
        )
    except BoundedIOError as failure:
        raise _fail(f"cannot write OMEX archive: {failure}") from failure

    return {
        **summary,
        "mode": "combine/omex",
        "output": str(output_archive),
        "archive_sha256": hashlib.sha256(archive_bytes).hexdigest(),
        "archive_bytes": len(archive_bytes),
    }


__all__ = [
    "ExportError",
    "export_omex",
    "export_rocrate",
]

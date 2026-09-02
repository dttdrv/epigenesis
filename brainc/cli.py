"""Command-line interface for the compiler-only distribution."""

from __future__ import annotations

import argparse
import json
import sys
from typing import Any

from ._canonical import ContractError, digest
from ._io import load_json_object, read_regular_file
from .bio import GFF3Compiler, load_gff3_artifact
from .bio_graph import compile_feature_graph
from .compiler import compile_program, load_policy, load_program, save as save_program
from .development_bundle import compile_development
from .external_profile import MAX_EVIDENCE_BYTES, MAX_PROFILE_BYTES
from .insdc import GenBankCompiler
from .insdc_graph import compile_insdc_graph
from .provider import (
    load_manifest,
    load_request,
    load_response,
    make_request,
    save as save_provider,
)
from .sequence import SequenceCompiler, SequenceCompilerError, load_sequence_artifact
from .sequence_collection import SequenceCollectionCompiler, SequenceCollectionError, load_sequence_collection
from .source import (
    EXTERNAL_PROFILE,
    FASTA_PROFILE,
    GFF3_PROFILE,
    GENBANK_PROFILE,
    PROFILES,
    RAW_PROFILE,
    REFERENCE_FASTA_PROFILE,
    ScaleLimits,
    compile_external_source_executable,
    compile_external_source_paths,
    compile_source,
)
from .validator import save_report, validate_chain
from .validator_bio import (
    MAX_BIO_ARTIFACT_BYTES,
    MAX_GFF3_BYTES as BIO_MAX_GFF3_BYTES,
    MAX_JSON_DEPTH as BIO_MAX_JSON_DEPTH,
    MAX_JSON_MEMBERS as BIO_MAX_JSON_MEMBERS,
    MAX_SOURCE_ARTIFACT_BYTES,
    MAX_SOURCE_INPUT_BYTES,
    MAX_TEXT_BYTES as BIO_MAX_TEXT_BYTES,
    validate_bio_chain,
)
from .validator_bio_graph import (
    validate_feature_graph_bundle as independently_validate_feature_graph_bundle,
)
from .validator_insdc import (
    INSDCValidationError,
    validate_genbank_paths as independently_validate_genbank_paths,
)
from .v2 import (
    compile_module,
    make_development_request,
    make_request as make_request_v2,
    save as save_v2,
)
from .validator_v2 import (
    load as load_v2,
    save_report as save_report_v2,
    source_record_lengths,
    validate_chain as validate_chain_v2,
    validate_manifest_artifact,
    validate_module_artifact,
    validate_policy_artifact,
    validate_request_artifact,
    validate_response_artifact,
    validate_source_artifact,
    validate_target_artifact,
)


def _role_path_argument(value: str) -> tuple[str, str]:
    role, separator, path = value.partition("=")
    if not separator or not role or not path:
        raise argparse.ArgumentTypeError("expected ROLE=PATH")
    return role, path


def _role_paths(values: list[tuple[str, str]], label: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for role, path in values:
        if role in result:
            raise ValueError(f"{label} role {role!r} was supplied more than once")
        result[role] = path
    return result


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="brainc",
        description="Deterministic compiler for content-bound biological sources",
    )
    parser.add_argument(
        "--diagnostics",
        choices=("text", "json"),
        default="text",
        help="render compilation failures as text or sealed JSON",
    )
    sub = parser.add_subparsers(dest="command", required=True)
    source = sub.add_parser(
        "compile-source",
        help="compile one explicitly selected DNA source profile to a source bundle",
    )
    source.add_argument(
        "--profile",
        choices=sorted(PROFILES - {EXTERNAL_PROFILE}),
        required=True,
    )
    source.add_argument("--sequence")
    source.add_argument("--genbank")
    source.add_argument("--annotation")
    source.add_argument("--record-id")
    source.add_argument("--wrapper", choices=("identity", "gzip"))
    for option in (
        "maximum-input-bytes",
        "maximum-logical-bytes",
        "maximum-records",
        "maximum-record-bases",
        "maximum-header-bytes",
        "maximum-gzip-members",
    ):
        source.add_argument(f"--{option}", type=int)
    source.add_argument(
        "--sequence-profile",
        choices=(RAW_PROFILE, FASTA_PROFILE),
        help="explicit external-sequence profile for GFF3",
    )
    source.add_argument("-o", "--output", required=True)
    external_source = sub.add_parser(
        "compile-external-source",
        help="compile a DNA grammar through pinned external tools",
    )
    external_source.add_argument("--profile-manifest", required=True)
    external_source.add_argument("--frontend-executable")
    external_source.add_argument("--validator-executable")
    external_source.add_argument("--source-descriptor")
    external_source.add_argument("--validation-report")
    external_source.add_argument(
        "--source-input",
        action="append",
        type=_role_path_argument,
        required=True,
        metavar="ROLE=PATH",
    )
    external_source.add_argument(
        "--native-artifact",
        action="append",
        type=_role_path_argument,
        default=[],
        metavar="ROLE=PATH",
    )
    external_source.add_argument("-o", "--output", required=True)
    sequence = sub.add_parser("compile-sequence", help="compile one FASTA record to Sequence IR")
    sequence.add_argument("input"); sequence.add_argument("--context"); sequence.add_argument("-o", "--output", required=True)
    collection = sub.add_parser("compile-collection", help="compile multi-FASTA or gzipped FASTA to collection IR")
    collection.add_argument("input"); collection.add_argument("-o", "--output", required=True)
    raw = sub.add_parser("compile-raw", help="compile one raw IUPAC DNA string file to collection IR")
    raw.add_argument("input"); raw.add_argument("--record-id", required=True); raw.add_argument("-o", "--output", required=True)
    genbank = sub.add_parser(
        "compile-genbank",
        help="compile GenBank physical-DNA records to content-bound structural IR",
    )
    genbank.add_argument("source")
    genbank.add_argument("-o", "--output", required=True)
    validate_genbank = sub.add_parser(
        "validate-genbank",
        help="independently replay GenBank bytes and validate a compiled artifact",
    )
    validate_genbank.add_argument("source")
    validate_genbank.add_argument("artifact")
    validate_genbank.add_argument("--report")
    insdc_graph = sub.add_parser(
        "compile-insdc-graph",
        help="lower compiled GenBank BioIR to deterministic feature state",
    )
    insdc_graph.add_argument("source_artifact")
    insdc_graph.add_argument("-o", "--output", required=True)
    gff3 = sub.add_parser(
        "compile-gff3",
        help="compile external-sequence GFF3 to a collection-bound BioIR artifact",
    )
    gff3.add_argument("gff3")
    gff3.add_argument("--sequence-collection", required=True)
    gff3.add_argument("-o", "--output", required=True)
    validate_gff3 = sub.add_parser(
        "validate-gff3",
        help="independently replay sequence and GFF3 sources into BioIR",
    )
    validate_gff3.add_argument("bio_ir")
    validate_gff3.add_argument("--sequence-collection", required=True)
    validate_gff3.add_argument("--sequence-input", required=True)
    validate_gff3.add_argument("--gff3-source", required=True)
    validate_gff3.add_argument("--source-metadata")
    validate_gff3.add_argument("--record-id")
    validate_gff3.add_argument("--report")
    feature_graph = sub.add_parser(
        "compile-feature-graph",
        help="mechanically lower collection-bound GFF3 BioIR to a structural graph bundle",
    )
    feature_graph.add_argument("bio_ir")
    feature_graph.add_argument("--sequence-collection", required=True)
    feature_graph.add_argument("--gff3-source")
    feature_graph.add_argument("-o", "--output", required=True)
    validate_feature_graph = sub.add_parser(
        "validate-feature-graph",
        help="independently reconstruct and validate a complete feature-graph bundle",
    )
    validate_feature_graph.add_argument("bundle")
    validate_feature_graph.add_argument("--sequence-collection", required=True)
    validate_feature_graph.add_argument("--bio-ir", required=True)
    validate_feature_graph.add_argument("--report")
    request = sub.add_parser("make-request", help="bind a sequence source to an external provider contract")
    request.add_argument("source"); request.add_argument("--manifest", required=True); request.add_argument("--output-id", action="append", required=True); request.add_argument("-o", "--output", required=True)
    request_v2 = sub.add_parser(
        "make-request-v2",
        help="bind a sequence source to a tensor-provider v2 contract",
    )
    request_v2.add_argument("source")
    request_v2.add_argument("--manifest", required=True)
    request_v2.add_argument("--output-id", action="append", required=True)
    request_v2.add_argument("-o", "--output", required=True)
    development_request = sub.add_parser(
        "make-development-request",
        help="bind a validated source bundle to a tensor-provider v2 contract",
    )
    development_request.add_argument("source_bundle")
    development_request.add_argument("--manifest", required=True)
    development_request.add_argument(
        "--output-id",
        action="append",
        required=True,
    )
    development_request.add_argument("-o", "--output", required=True)
    compile_cmd = sub.add_parser("compile", help="lower an externally predicted DNA source to state-program IR")
    compile_cmd.add_argument("source"); compile_cmd.add_argument("--manifest", required=True); compile_cmd.add_argument("--request", required=True)
    compile_cmd.add_argument("--response", required=True); compile_cmd.add_argument("--policy", required=True); compile_cmd.add_argument("-o", "--output", required=True)
    compile_v2 = sub.add_parser(
        "compile-v2",
        help="lower tensor-provider artifacts to a target-bound development module",
    )
    compile_v2.add_argument("source")
    compile_v2.add_argument("--manifest", required=True)
    compile_v2.add_argument("--request", required=True)
    compile_v2.add_argument("--response", required=True)
    compile_v2.add_argument("--policy", required=True)
    compile_v2.add_argument("--target", required=True)
    compile_v2.add_argument("--blob-root")
    compile_v2.add_argument("-o", "--output", required=True)
    development = sub.add_parser(
        "compile-development",
        help="compile a source bundle and caller-supplied interpretation to a development bundle",
    )
    development.add_argument("source_bundle")
    development.add_argument("--manifest", required=True)
    development.add_argument("--request", required=True)
    development.add_argument("--response", required=True)
    development.add_argument("--policy", required=True)
    development.add_argument("--target", required=True)
    development.add_argument("--blob-root")
    development.add_argument("-o", "--output", required=True)
    check = sub.add_parser("check", help="check one compiler artifact's closed schema and digest")
    check.add_argument("kind", choices=["sequence", "collection", "manifest", "request", "response", "policy", "program"]); check.add_argument("artifact")
    check_v2 = sub.add_parser(
        "check-v2",
        help="independently check one v2 artifact's closed contract and digest",
    )
    check_v2_kind = check_v2.add_subparsers(dest="kind", required=True)
    for kind in ("source", "manifest", "request", "policy", "target"):
        kind_parser = check_v2_kind.add_parser(kind)
        kind_parser.add_argument("artifact")
    response_check = check_v2_kind.add_parser("response")
    response_check.add_argument("artifact")
    response_check.add_argument(
        "--source",
        required=True,
        help="sequence source artifact used to verify tensor coordinate bindings",
    )
    response_check.add_argument("--blob-root")
    module_check = check_v2_kind.add_parser("module")
    module_check.add_argument("artifact")
    module_check.add_argument(
        "--target",
        required=True,
        help="target contract used to type-check module operations",
    )
    module_check.add_argument("--blob-root")
    validate = sub.add_parser("validate", help="independently replay a complete compiler chain")
    validate.add_argument("--source-input", "--fasta", dest="fasta", required=True)
    validate.add_argument("--context"); validate.add_argument("--record-id")
    validate.add_argument("--source-artifact", "--sequence", dest="sequence", required=True)
    validate.add_argument("--manifest", required=True); validate.add_argument("--request", required=True); validate.add_argument("--response", required=True)
    validate.add_argument("--policy", required=True); validate.add_argument("--program", required=True); validate.add_argument("--report")
    validate_v2 = sub.add_parser(
        "validate-v2",
        help="independently replay a complete tensor compiler chain",
    )
    validate_v2.add_argument(
        "--source-input", "--fasta", dest="fasta", required=True
    )
    validate_v2.add_argument("--context")
    validate_v2.add_argument("--record-id")
    validate_v2.add_argument(
        "--source-artifact", "--sequence", dest="sequence", required=True
    )
    validate_v2.add_argument("--manifest", required=True)
    validate_v2.add_argument("--request", required=True)
    validate_v2.add_argument("--response", required=True)
    validate_v2.add_argument("--policy", required=True)
    validate_v2.add_argument("--target", required=True)
    validate_v2.add_argument("--module", required=True)
    validate_v2.add_argument("--blob-root")
    validate_v2.add_argument("--report")
    return parser


def _check(kind: str, artifact: str) -> Any:
    return {
        "sequence": load_sequence_artifact,
        "collection": load_sequence_collection,
        "manifest": load_manifest,
        "request": load_request,
        "response": load_response,
        "policy": load_policy,
        "program": load_program,
    }[kind](artifact)


def _source_arguments(
    args: argparse.Namespace,
) -> tuple[dict[str, str], dict[str, Any]]:
    values = {
        name: getattr(args, name)
        for name in (
            "sequence",
            "genbank",
            "annotation",
            "record_id",
            "wrapper",
            "sequence_profile",
        )
    }

    def selected(*, required: set[str], allowed: set[str]) -> None:
        missing = sorted(name for name in required if values[name] is None)
        extra = sorted(
            name
            for name in values
            if name not in allowed and values[name] is not None
        )
        if missing or extra:
            raise ValueError(
                f"source profile flags invalid; missing={missing or 'none'}, "
                f"not-applicable={extra or 'none'}"
            )

    if args.profile == RAW_PROFILE:
        selected(required={"sequence", "record_id"}, allowed={"sequence", "record_id"})
        return {"sequence": args.sequence}, {"record_id": args.record_id}
    if args.profile in {FASTA_PROFILE, REFERENCE_FASTA_PROFILE}:
        selected(required={"sequence", "wrapper"}, allowed={"sequence", "wrapper"})
        return {"sequence": args.sequence}, {"wrapper": args.wrapper}
    if args.profile == GENBANK_PROFILE:
        selected(required={"genbank"}, allowed={"genbank"})
        return {"genbank": args.genbank}, {}

    selected(
        required={"sequence", "annotation", "sequence_profile"},
        allowed={"sequence", "annotation", "sequence_profile", "record_id", "wrapper"},
    )
    if args.sequence_profile == RAW_PROFILE:
        if args.record_id is None or args.wrapper is not None:
            raise ValueError(
                "raw GFF3 sequence profile requires --record-id and forbids --wrapper"
            )
        sequence_parameters = {"record_id": args.record_id}
    else:
        if args.wrapper is None or args.record_id is not None:
            raise ValueError(
                "FASTA GFF3 sequence profile requires --wrapper and forbids --record-id"
            )
        sequence_parameters = {"wrapper": args.wrapper}
    return (
        {"sequence": args.sequence, "annotation": args.annotation},
        {
            "sequence_profile": args.sequence_profile,
            "sequence_parameters": sequence_parameters,
        },
    )


def _source_limits(args: argparse.Namespace) -> ScaleLimits | None:
    names = (
        "maximum_input_bytes",
        "maximum_logical_bytes",
        "maximum_records",
        "maximum_record_bases",
        "maximum_header_bytes",
        "maximum_gzip_members",
    )
    overrides = {
        name: getattr(args, name)
        for name in names
        if getattr(args, name) is not None
    }
    if not overrides:
        return None
    if args.profile != REFERENCE_FASTA_PROFILE:
        raise ValueError(
            "source execution limits apply only to the reference FASTA profile"
        )
    defaults = ScaleLimits()
    values = {name: getattr(defaults, name) for name in names}
    values.update(overrides)
    return ScaleLimits(**values)


def _check_v2(args: argparse.Namespace) -> None:
    artifact, _ = load_v2(args.artifact, f"{args.kind} artifact")
    if args.kind == "source":
        validate_source_artifact(artifact)
    elif args.kind == "manifest":
        validate_manifest_artifact(artifact)
    elif args.kind == "request":
        validate_request_artifact(artifact)
    elif args.kind == "response":
        source, _ = load_v2(args.source, "sequence source artifact")
        records = source_record_lengths(source)
        validate_response_artifact(
            artifact,
            sequence_sha256=source["artifact_sha256"],
            records=records,
            blob_root=args.blob_root,
        )
    elif args.kind == "policy":
        validate_policy_artifact(artifact)
    elif args.kind == "target":
        validate_target_artifact(artifact)
    elif args.kind == "module":
        target, _ = load_v2(args.target, "target contract")
        validate_target_artifact(target)
        validate_module_artifact(
            artifact,
            target,
            blob_root=args.blob_root,
        )
    else:  # argparse keeps this unreachable; retain a closed internal boundary.
        raise ValueError(f"unsupported v2 artifact kind: {args.kind}")


def _bio_json(path: str, label: str, *, sequence_strings: bool) -> dict[str, Any]:
    maximum_bytes = (
        MAX_SOURCE_ARTIFACT_BYTES if sequence_strings else MAX_BIO_ARTIFACT_BYTES
    )
    raw = read_regular_file(
        path, maximum_bytes=maximum_bytes, label=label
    )
    return load_json_object(
        raw,
        label,
        maximum_bytes=maximum_bytes,
        maximum_depth=BIO_MAX_JSON_DEPTH,
        maximum_members=BIO_MAX_JSON_MEMBERS,
        maximum_string_bytes=BIO_MAX_TEXT_BYTES,
    )


def _external_json(path: str, label: str, maximum_bytes: int) -> dict[str, Any]:
    raw = read_regular_file(path, maximum_bytes=maximum_bytes, label=label)
    return load_json_object(raw, label, maximum_bytes=maximum_bytes)


def _emit_report(report: dict[str, Any], output: str | None) -> int:
    if output is not None:
        save_v2(report, output)
    print(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False))
    return 0


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "compile-source":
            inputs, parameters = _source_arguments(args)
            paths = compile_source(
                args.profile,
                inputs,
                parameters=parameters,
                limits=_source_limits(args),
            ).save(args.output)
            print(
                json.dumps(
                    {
                        "output": str(args.output),
                        "source": str(paths["source.json"]),
                    },
                    sort_keys=True,
                )
            )
        elif args.command == "compile-external-source":
            manifest = _external_json(
                args.profile_manifest,
                "external profile manifest",
                MAX_PROFILE_BYTES,
            )
            source_inputs = _role_paths(args.source_input, "source input")
            executable_mode = bool(
                args.frontend_executable or args.validator_executable
            )
            if executable_mode:
                if not args.frontend_executable or not args.validator_executable:
                    raise ValueError(
                        "executable mode requires both --frontend-executable and --validator-executable"
                    )
                if args.source_descriptor or args.validation_report or args.native_artifact:
                    raise ValueError(
                        "executable mode does not accept legacy descriptor, report, or native-artifact arguments"
                    )
                bundle = compile_external_source_executable(
                    manifest,
                    original_paths=source_inputs,
                    frontend_executable=args.frontend_executable,
                    validator_executable=args.validator_executable,
                )
            else:
                if not args.source_descriptor or not args.validation_report or not args.native_artifact:
                    raise ValueError(
                        "legacy mode requires --source-descriptor, --validation-report, and --native-artifact"
                    )
                bundle = compile_external_source_paths(
                    manifest,
                    _external_json(
                        args.source_descriptor,
                        "external source descriptor",
                        MAX_EVIDENCE_BYTES,
                    ),
                    _external_json(
                        args.validation_report,
                        "external validation report",
                        MAX_EVIDENCE_BYTES,
                    ),
                    original_paths=source_inputs,
                    native_artifact_paths=_role_paths(
                        args.native_artifact,
                        "native artifact",
                    ),
                )
            paths = bundle.save(args.output)
            print(
                json.dumps(
                    {
                        "output": str(args.output),
                        "source": str(paths["source.json"]),
                    },
                    sort_keys=True,
                )
            )
        elif args.command == "compile-sequence":
            SequenceCompiler().compile_file(args.input, args.context).save(args.output)
        elif args.command == "compile-collection":
            SequenceCollectionCompiler().compile_file(args.input).save(args.output)
        elif args.command == "compile-raw":
            SequenceCollectionCompiler().compile_raw_file(
                args.input, args.record_id
            ).save(args.output)
        elif args.command == "compile-genbank":
            GenBankCompiler().compile_file(args.source).save(args.output)
        elif args.command == "validate-genbank":
            try:
                report = independently_validate_genbank_paths(
                    args.artifact,
                    args.source,
                    report_path=args.report,
                )
            except INSDCValidationError as failure:
                print(f"brainc: validation failed: {failure}", file=sys.stderr)
                return 3
            print(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False))
            return 0
        elif args.command == "compile-insdc-graph":
            paths = compile_insdc_graph(args.source_artifact).save(args.output)
            print(
                json.dumps(
                    {"bundle": str(paths["bundle"]), "output": str(args.output)},
                    sort_keys=True,
                )
            )
        elif args.command == "compile-gff3":
            collection = load_sequence_collection(args.sequence_collection)
            GFF3Compiler().compile_file(args.gff3, collection).save(args.output)
        elif args.command == "validate-gff3":
            report = validate_bio_chain(
                _bio_json(args.bio_ir, "BioIR artifact", sequence_strings=False),
                _bio_json(
                    args.sequence_collection,
                    "sequence collection artifact",
                    sequence_strings=True,
                ),
                fasta_source=read_regular_file(
                    args.sequence_input,
                    maximum_bytes=MAX_SOURCE_INPUT_BYTES,
                    label="sequence source",
                ),
                gff3_source=read_regular_file(
                    args.gff3_source,
                    maximum_bytes=BIO_MAX_GFF3_BYTES,
                    label="GFF3 source",
                ),
                source_metadata=(
                    None
                    if args.source_metadata is None
                    else read_regular_file(
                        args.source_metadata,
                        maximum_bytes=BIO_MAX_TEXT_BYTES,
                        label="source metadata",
                    )
                ),
                record_id=args.record_id,
            )
            return _emit_report(report, args.report)
        elif args.command == "compile-feature-graph":
            collection = load_sequence_collection(args.sequence_collection)
            gff3_source = (
                None
                if args.gff3_source is None
                else read_regular_file(
                    args.gff3_source,
                    maximum_bytes=BIO_MAX_GFF3_BYTES,
                    label="GFF3 source",
                )
            )
            bio_ir = load_gff3_artifact(
                args.bio_ir, collection, gff3_source=gff3_source
            )
            bundle = compile_feature_graph(
                args.sequence_collection, bio_ir, gff3_source=gff3_source
            )
            paths = bundle.save(args.output)
            print(
                json.dumps(
                    {
                        "bundle": str(paths["bundle"]),
                        "output": str(args.output),
                    },
                    sort_keys=True,
                )
            )
        elif args.command == "validate-feature-graph":
            report = independently_validate_feature_graph_bundle(
                args.bundle, args.sequence_collection, args.bio_ir
            )
            return _emit_report(report, args.report)
        elif args.command == "make-request":
            save_provider(make_request(args.source, args.manifest, args.output_id), args.output)
        elif args.command == "make-request-v2":
            save_v2(
                make_request_v2(args.source, args.manifest, args.output_id),
                args.output,
            )
        elif args.command == "make-development-request":
            save_v2(
                make_development_request(
                    args.source_bundle,
                    args.manifest,
                    args.output_id,
                ),
                args.output,
            )
        elif args.command == "compile":
            save_program(compile_program(args.source, args.manifest, args.request, args.response, args.policy), args.output)
        elif args.command == "compile-v2":
            save_v2(
                compile_module(
                    args.source,
                    args.manifest,
                    args.request,
                    args.response,
                    args.policy,
                    args.target,
                    blob_root=args.blob_root,
                ),
                args.output,
            )
        elif args.command == "compile-development":
            paths = compile_development(
                args.source_bundle,
                args.manifest,
                args.request,
                args.response,
                args.policy,
                args.target,
                blob_root=args.blob_root,
            ).save(args.output)
            print(
                json.dumps(
                    {
                        "bundle": str(paths["bundle.json"]),
                        "output": str(args.output),
                    },
                    sort_keys=True,
                )
            )
        elif args.command == "check":
            _check(args.kind, args.artifact)
            print(json.dumps({"valid": True, "kind": args.kind, "artifact": args.artifact}, sort_keys=True))
        elif args.command == "check-v2":
            _check_v2(args)
            print(
                json.dumps(
                    {
                        "valid": True,
                        "version": 2,
                        "kind": args.kind,
                        "artifact": args.artifact,
                    },
                    sort_keys=True,
                )
            )
        elif args.command == "validate":
            report = validate_chain(fasta=args.fasta, context=args.context, sequence=args.sequence,
                                    manifest=args.manifest, request=args.request, response=args.response,
                                    policy=args.policy, program=args.program, record_id=args.record_id)
            if args.report: save_report(report, args.report)
            print(json.dumps(report, indent=2, sort_keys=True))
            return 0 if report["valid"] else 3
        elif args.command == "validate-v2":
            report = validate_chain_v2(
                fasta=args.fasta,
                context=args.context,
                sequence=args.sequence,
                manifest=args.manifest,
                request=args.request,
                response=args.response,
                policy=args.policy,
                target=args.target,
                module=args.module,
                blob_root=args.blob_root,
                record_id=args.record_id,
            )
            if args.report:
                save_report_v2(report, args.report)
            print(json.dumps(report, indent=2, sort_keys=True))
            return 0 if report["valid"] else 3
        return 0
    except (OSError, ContractError, SequenceCompilerError, SequenceCollectionError, ValueError) as failure:
        if args.diagnostics == "json":
            message = str(failure)
            candidate, separator, _ = message.partition(":")
            code = (
                candidate
                if (
                    separator
                    and candidate.isascii()
                    and candidate.isalnum()
                    and candidate.isupper()
                )
                else "CLI001"
            )
            core = {
                "format": "brainc.compiler-diagnostic",
                "version": 1,
                "valid": False,
                "command": args.command,
                "error": {"code": code, "message": message},
            }
            diagnostic = {**core, "diagnostic_sha256": digest(core)}
            print(
                json.dumps(
                    diagnostic,
                    sort_keys=True,
                    separators=(",", ":"),
                    ensure_ascii=False,
                ),
                file=sys.stderr,
            )
        else:
            print(f"brainc: {failure}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

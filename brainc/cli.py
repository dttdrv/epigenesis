"""Command-line interface for the compiler-only distribution."""

from __future__ import annotations

import argparse
import json
import sys
from typing import Any

from ._canonical import ContractError
from ._io import load_json_object, read_regular_file
from .bio import GFF3Compiler, load_gff3_artifact
from .bio_graph import compile_feature_graph
from .compiler import compile_program, load_policy, load_program, save as save_program
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
from .v2 import compile_module, make_request as make_request_v2, save as save_v2
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


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="brainc",
        description="Deterministic compiler for content-bound biological sources",
    )
    sub = parser.add_subparsers(dest="command", required=True)
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


def _emit_report(report: dict[str, Any], output: str | None) -> int:
    if output is not None:
        save_v2(report, output)
    print(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False))
    return 0


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "compile-sequence":
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
        print(f"brainc: {failure}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

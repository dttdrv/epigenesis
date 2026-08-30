"""Command-line interface for the compiler-only distribution."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

from ._canonical import ContractError
from .compiler import compile_program, load_policy, load_program, save as save_program
from .provider import load_manifest, load_request, load_response, make_request, save as save_provider
from .sequence import SequenceCompiler, SequenceCompilerError, load_sequence_artifact
from .sequence_collection import SequenceCollectionCompiler, SequenceCollectionError, load_sequence_collection
from .validator import save_report, validate_chain


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="brainc", description="Offline DNA-to-state-program compiler kernel")
    sub = parser.add_subparsers(dest="command", required=True)
    sequence = sub.add_parser("compile-sequence", help="compile one FASTA record to Sequence IR")
    sequence.add_argument("input"); sequence.add_argument("--context"); sequence.add_argument("-o", "--output", required=True)
    collection = sub.add_parser("compile-collection", help="compile multi-FASTA or gzipped FASTA to collection IR")
    collection.add_argument("input"); collection.add_argument("-o", "--output", required=True)
    raw = sub.add_parser("compile-raw", help="compile one raw IUPAC DNA string file to collection IR")
    raw.add_argument("input"); raw.add_argument("--record-id", required=True); raw.add_argument("-o", "--output", required=True)
    request = sub.add_parser("make-request", help="bind a sequence source to an external provider contract")
    request.add_argument("source"); request.add_argument("--manifest", required=True); request.add_argument("--output-id", action="append", required=True); request.add_argument("-o", "--output", required=True)
    compile_cmd = sub.add_parser("compile", help="lower an externally predicted DNA source to state-program IR")
    compile_cmd.add_argument("source"); compile_cmd.add_argument("--manifest", required=True); compile_cmd.add_argument("--request", required=True)
    compile_cmd.add_argument("--response", required=True); compile_cmd.add_argument("--policy", required=True); compile_cmd.add_argument("-o", "--output", required=True)
    check = sub.add_parser("check", help="check one compiler artifact's closed schema and digest")
    check.add_argument("kind", choices=["sequence", "collection", "manifest", "request", "response", "policy", "program"]); check.add_argument("artifact")
    validate = sub.add_parser("validate", help="independently replay a complete compiler chain")
    validate.add_argument("--source-input", "--fasta", dest="fasta", required=True)
    validate.add_argument("--context"); validate.add_argument("--record-id")
    validate.add_argument("--source-artifact", "--sequence", dest="sequence", required=True)
    validate.add_argument("--manifest", required=True); validate.add_argument("--request", required=True); validate.add_argument("--response", required=True)
    validate.add_argument("--policy", required=True); validate.add_argument("--program", required=True); validate.add_argument("--report")
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


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "compile-sequence":
            SequenceCompiler().compile_file(args.input, args.context).save(args.output)
        elif args.command == "compile-collection":
            SequenceCollectionCompiler().compile_file(args.input).save(args.output)
        elif args.command == "compile-raw":
            SequenceCollectionCompiler().compile_raw(Path(args.input).read_bytes(), args.record_id).save(args.output)
        elif args.command == "make-request":
            save_provider(make_request(args.source, args.manifest, args.output_id), args.output)
        elif args.command == "compile":
            save_program(compile_program(args.source, args.manifest, args.request, args.response, args.policy), args.output)
        elif args.command == "check":
            _check(args.kind, args.artifact)
            print(json.dumps({"valid": True, "kind": args.kind, "artifact": args.artifact}, sort_keys=True))
        elif args.command == "validate":
            report = validate_chain(fasta=args.fasta, context=args.context, sequence=args.sequence,
                                    manifest=args.manifest, request=args.request, response=args.response,
                                    policy=args.policy, program=args.program, record_id=args.record_id)
            if args.report: save_report(report, args.report)
            print(json.dumps(report, indent=2, sort_keys=True))
            return 0 if report["valid"] else 3
        return 0
    except (OSError, ContractError, SequenceCompilerError, SequenceCollectionError, ValueError) as failure:
        print(f"brainc: {failure}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

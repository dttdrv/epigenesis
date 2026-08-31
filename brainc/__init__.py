"""Epigenesis public API.

Public compiler symbols are loaded on first access. Keeping package import
side-effect free lets the independently implemented validators run without
loading compiler modules while preserving ``from brainc import ...``.
"""

from __future__ import annotations

from typing import Any


__version__ = "0.9.0"

__all__ = [
    "CompilerError",
    "GenBankArtifact",
    "GenBankCompiler",
    "GenBankError",
    "INSDCGraphBundle",
    "INSDCGraphError",
    "INSDCGraphValidationError",
    "INSDCValidationError",
    "ProviderError",
    "SequenceArtifact",
    "SequenceCollectionArtifact",
    "SequenceCollectionCompiler",
    "SequenceCollectionError",
    "SequenceCompiler",
    "SequenceCompilerError",
    "compile_program",
    "compile_insdc_graph",
    "load_manifest",
    "load_genbank_artifact",
    "load_insdc_graph_bundle_directory",
    "load_policy",
    "load_program",
    "load_request",
    "load_response",
    "load_sequence_artifact",
    "load_sequence_collection",
    "make_request",
    "validate_binding",
    "validate_chain",
    "validate_genbank",
    "validate_genbank_artifact",
    "validate_genbank_paths",
    "validate_genbank_report",
    "validate_insdc_graph_paths",
    "validate_insdc_graph_report",
]

def __getattr__(name: str) -> Any:
    if name in {"CompilerError", "compile_program", "load_policy", "load_program"}:
        from . import compiler as module
    elif name in {
        "GenBankArtifact",
        "GenBankCompiler",
        "GenBankError",
        "load_genbank_artifact",
        "validate_genbank_artifact",
    }:
        from . import insdc as module
    elif name in {
        "INSDCGraphBundle",
        "INSDCGraphError",
        "compile_insdc_graph",
    }:
        from . import insdc_graph as module
    elif name in {
        "INSDCGraphValidationError",
        "load_insdc_graph_bundle_directory",
        "validate_insdc_graph_paths",
        "validate_insdc_graph_report",
    }:
        from . import validator_insdc_graph as module
    elif name in {
        "INSDCValidationError",
        "validate_genbank",
        "validate_genbank_paths",
        "validate_genbank_report",
    }:
        from . import validator_insdc as module
    elif name in {
        "ProviderError",
        "load_manifest",
        "load_request",
        "load_response",
        "make_request",
        "validate_binding",
    }:
        from . import provider as module
    elif name in {
        "SequenceArtifact",
        "SequenceCompiler",
        "SequenceCompilerError",
        "load_sequence_artifact",
    }:
        from . import sequence as module
    elif name in {
        "SequenceCollectionArtifact",
        "SequenceCollectionCompiler",
        "SequenceCollectionError",
        "load_sequence_collection",
    }:
        from . import sequence_collection as module
    elif name == "validate_chain":
        from . import validator as module
    else:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(module, name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))

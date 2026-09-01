"""Epigenesis public API.

Public compiler symbols are loaded on first access. Keeping package import
side-effect free lets the independently implemented validators run without
loading compiler modules while preserving ``from brainc import ...``.
"""

from __future__ import annotations

from typing import Any


__version__ = "1.0.0"

__all__ = [
    "CompilerError",
    "DevelopmentBundle",
    "DevelopmentBundleError",
    "EXTERNAL_PROFILE",
    "FASTA_PROFILE",
    "GFF3_PROFILE",
    "GenBankArtifact",
    "GenBankCompiler",
    "GenBankError",
    "GENBANK_PROFILE",
    "INSDCGraphBundle",
    "INSDCGraphError",
    "INSDCGraphValidationError",
    "INSDCValidationError",
    "ProviderError",
    "PROFILES",
    "RAW_PROFILE",
    "REFERENCE_FASTA_PROFILE",
    "ScaleLimits",
    "SequenceArtifact",
    "SequenceCollectionArtifact",
    "SequenceCollectionCompiler",
    "SequenceCollectionError",
    "SequenceCompiler",
    "SequenceCompilerError",
    "SourceBundle",
    "SourceError",
    "compile_development",
    "compile_external_source",
    "compile_external_source_paths",
    "compile_program",
    "compile_insdc_graph",
    "compile_source",
    "load_manifest",
    "load_genbank_artifact",
    "load_insdc_graph_bundle_directory",
    "load_policy",
    "load_program",
    "load_request",
    "load_response",
    "load_sequence_artifact",
    "load_sequence_collection",
    "load_source_bundle",
    "make_development_request",
    "make_request",
    "validate_binding",
    "validate_chain",
    "validate_genbank",
    "validate_genbank_artifact",
    "validate_genbank_paths",
    "validate_genbank_report",
    "validate_insdc_graph_paths",
    "validate_insdc_graph_report",
    "validate_source_bundle",
    "validate_source_descriptor",
]

def __getattr__(name: str) -> Any:
    if name in {"CompilerError", "compile_program", "load_policy", "load_program"}:
        from . import compiler as module
    elif name in {
        "DevelopmentBundle",
        "DevelopmentBundleError",
        "compile_development",
    }:
        from . import development_bundle as module
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
    elif name in {
        "FASTA_PROFILE",
        "EXTERNAL_PROFILE",
        "GFF3_PROFILE",
        "GENBANK_PROFILE",
        "PROFILES",
        "RAW_PROFILE",
        "REFERENCE_FASTA_PROFILE",
        "ScaleLimits",
        "SourceBundle",
        "SourceError",
        "compile_source",
        "compile_external_source",
        "compile_external_source_paths",
        "load_source_bundle",
        "validate_source_bundle",
        "validate_source_descriptor",
    }:
        from . import source as module
    elif name == "validate_chain":
        from . import validator as module
    elif name == "make_development_request":
        from .v2 import provider as module
    else:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(module, name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))

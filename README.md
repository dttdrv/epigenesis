# Epigenesis

Epigenesis is a deterministic compiler kernel for DNA-bound state initialization.

It accepts raw IUPAC DNA, FASTA, multi-FASTA, and gzip-wrapped input; emits content-addressed Sequence IR; binds that exact source to an external interpretation provider; type-checks a caller-authored lowering policy; emits scalar state initializer IR; and independently replays every compiler-owned transformation and binding.

The compiler contains no training data, model weights, predictor, network client, simulator, world, runtime learner, checkpoint, or finalizer.

## Pipeline

    DNA source
      → Sequence or Collection IR
      → bound external provider request/response
      → typed lowering policy
      → scalar state initializer IR
      → independent replay validation

## Install

    python -m pip install .

## Use

    brainc compile-collection genome.fa -o source.json
    brainc make-request source.json --manifest provider.json \
      --output-id regulatory.score -o request.json
    brainc compile source.json --manifest provider.json \
      --request request.json --response response.json \
      --policy policy.json -o program.json
    brainc validate --source-input genome.fa --source-artifact source.json \
      --manifest provider.json --request request.json \
      --response response.json --policy policy.json \
      --program program.json

`--fasta` and `--sequence` remain aliases for the two validation source flags. Raw collection replay additionally requires `--record-id`.

## What “compiler” means here

DNA alone does not specify a computer program or a brain. The source package is the DNA, a versioned external interpretation, and an explicit lowering policy. Epigenesis checks and deterministically translates that package into typed IR. The current target is a validated scalar initializer, not executable machine code and not a biological simulation.

This is not the first system described as a DNA compiler. Attribute-grammar work compiled constrained synthetic DNA to dynamical models in 2009 ([Cai et al.](https://doi.org/10.1371/journal.pcbi.1000529)), while several established systems compile programs or circuits in the opposite direction, into DNA. A targeted primary-source review did not identify an exact precedent for this implementation’s combination of real-DNA ingress, an external interpretation boundary, typed lowering, content-addressed provenance, and a separate replay validator.

## Contracts and limits

Artifact hashes use [RFC 8785 JSON canonicalization](https://www.rfc-editor.org/rfc/rfc8785). Sequences and ordered collections expose [GA4GH refget](https://ga4gh.github.io/refget/sequences/) and [Sequence Collections](https://ga4gh.github.io/refget/seqcols/) identities. The provider ABI currently supports scalar number, integer, and boolean outputs.

FASTQ reads, VCF variants, annotations, tensors, target execution, biological correctness, and provider execution are outside v0.4. Inputs are materialized in memory, so format support is not yet a practical whole-assembly scaling claim. These boundaries are specified in [SPEC.md](SPEC.md).

## Evidence

The release gates compile the complete 5,386-base phiX174 genome ([NCBI `J02482.1`](https://www.ncbi.nlm.nih.gov/nuccore/J02482.1)) through collection IR, provider binding, lowering, and independent validation. The provider response in this gate is a synthetic contract fixture, so the test proves compiler integrity over real DNA and makes no biological interpretation claim. The gates also reject fully rehashed state-value and source-map forgeries, exercise RFC 8785 vectors and unsafe-number rejection, check the model-free boundary, and install/test the built wheel.

Run everything:

    python -m unittest discover -s tests -v
    python tests/verify.py real-dna
    python tests/verify.py tamper
    python tests/verify.py canonical
    python tests/verify.py contract-attacks
    python tests/verify.py boundary
    PIP_NO_INDEX=1 python tests/verify.py wheel

The NCBI fixture retains its accession, retrieval URL, and checksums in `tests/data/J02482.1.source.json`; it is public molecular data, not training data or project-authored source code.

Apache-2.0.

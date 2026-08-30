# Epigenesis

Epigenesis is a deterministic compiler kernel that lowers DNA and explicit biological interpretation into typed state-program IR.

It parses raw IUPAC DNA, FASTA, multi-FASTA, and gzipped FASTA; emits content-addressed sequence IR; exchanges bound files with an external prediction provider; type-checks and lowers a caller-authored policy; then independently replays the compilation.

The compiler contains no training data, model weights, predictor, network client, simulator, world, runtime learner, checkpoint, or finalizer.

## Pipeline

    DNA
      → Sequence IR
      → external provider request/response
      → typed lowering policy
      → state-program IR
      → independent validation

## Install

    python -m pip install .

## Use

    brainc compile-sequence input.fa -o sequence.json
    brainc make-request sequence.json --manifest provider.json \
      --output-id regulatory.score -o request.json
    brainc compile sequence.json --manifest provider.json \
      --request request.json --response response.json \
      --policy policy.json -o program.json
    brainc validate --fasta input.fa --sequence sequence.json \
      --manifest provider.json --request request.json \
      --response response.json --policy policy.json \
      --program program.json

## Boundary

Epigenesis validates provider identity, request binding, schemas, types, hashes, lowering arithmetic, and emitted target equality. It does not claim that an external provider's prediction is biologically correct.

The current provider ABI supports scalar number, integer, and boolean outputs. Track tensors, variant materialization, provider attestations, and execution backends remain separate extensions.

## Evidence

- Eight compiler-only regression tests
- Two unrelated external providers through one ABI
- Raw, single-record, multi-record, and gzipped DNA ingress
- Deterministic output with one-base sensitivity
- Rejection of fully rehashed target forgery
- Independent standard-library validator
- Reproducible, dependency-free wheel

No software license has been selected.

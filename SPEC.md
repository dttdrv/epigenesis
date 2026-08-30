# Epigenesis compiler contract 0.4

## Boundary

The compiler translates an explicit source package into replayable scalar state initializer IR:

1. one or more named IUPAC DNA strings;
2. a content-identified external provider contract and bound response;
3. a caller-authored typed lowering policy.

It proves syntactic validity, identity, binding, type correctness, deterministic arithmetic, and emitted-IR equality. It does not prove a provider’s biological claims.

## DNA source

Accepted logical inputs are one raw whitespace-free DNA string with an explicit record ID, or one or more FASTA records. Symbols are case-preserving `A C G T R Y S W K M B D H V N`; identity is computed over uppercase sequence. Gzip is a transport wrapper identified by magic bytes. Original and decompressed byte hashes are both retained.

FASTQ is a read-and-quality format, VCF is reference-relative variation, and GFF is annotation. They require additional typed semantics and are not accepted as aliases for DNA sequence.

Each member embeds local Sequence IR and a local source map. Collection members additionally carry source-global line maps. Duplicate record IDs are invalid. Sequence identity follows GA4GH refget Sequences v2; ordered collection identity follows refget Sequence Collections v1 with names and sequences as inherent attributes and lengths as non-inherent.

## Canonical values

Content digests are lowercase SHA-256 over RFC 8785 canonical JSON. Accepted values are the I-JSON subset represented by null, booleans, Unicode strings without lone surrogates, arrays, string-keyed objects, finite IEEE-754 binary64 numbers, and integers in `[-(2^53-1), 2^53-1]`. Duplicate JSON keys, non-finite numbers, and unsafe integers are invalid.

## External provider ABI

A manifest declares provider and model identity, accepted source tags, and ordered scalar outputs. A request binds the exact source and requested declared contracts. A response binds the exact request and must preserve provider/model identity and requested output order, ID, type, and unit. Output types are `number`, `integer`, and `boolean`.

Provider execution, model contents, training data, and biological truth remain outside the compiler and validator.

## Lowering

A policy declares target name/version, unique states, directed links, and ports. Numeric states apply finite `scale * input + offset` and an optional inclusive clamp. Integer states require nearest rounding and must resolve nonnegative; number states require no rounding. Boolean states pass through boolean outputs. Arithmetic overflow, unsafe integer results, duplicate graph elements, and unknown endpoints are invalid.

The emitted target has no runtime operational semantics in v0.4. It is a typed, content-addressed scalar initializer. A content-bound target contract and tensor ABI are later compiler extensions; tissue construction, learning, memory, simulation, and model finalization are downstream systems.

## Independent validation

The validator imports no compiler module. It reconstructs Sequence or Collection IR from the original bytes, rechecks every closed schema and digest, rebinds the provider chain, recomputes lowering arithmetic, and compares the complete expected program using canonical bytes. Rehashing modified semantics therefore does not make a forgery valid.

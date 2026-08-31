# Epigenesis compiler contract 0.9.0

## 1. Status

This document specifies package 0.9.0: the GenBank 0.8 frontend, Sequence
Collection v2, structural BioIR v2, independent validation report v2, and the
INSDC feature-state adapter v1. The GenBank compiler producer identity is 0.8.1;
the independent GenBank validator identity is 0.2.0.

The normative source profile is:

```text
genbank-273-traditional-dna-physical-structural/v2
```

It is pinned by authority-manifest SHA-256
`6dc3658a5f7d6a774ce9cd08ff6dc5327cbef8a9951c7679fd9a34650f912bf9`
to:

- NCBI GenBank Flat File Release 273.0, dated 15 August 2026; and
- INSDC Feature Table Definition 11.4, dated April 2026.

The words MUST, MUST NOT, REQUIRED, SHALL, SHALL NOT, SHOULD, SHOULD NOT,
RECOMMENDED, NOT RECOMMENDED, MAY, and OPTIONAL are interpreted as described by
[BCP 14](https://www.rfc-editor.org/info/bcp14) when capitalized.

## 2. Compilation model

The frontend accepts original GenBank bytes under the named profile and emits
one artifact containing:

1. an exact-source `brain01.sequence-collection-ir/v2`;
2. a collection-bound `brainc.bio.insdc-genbank-ir/v2` BioIR body;
3. a digest of that BioIR; and
4. a digest sealing the complete artifact.

Compilation consists of lexical analysis, fixed-field and feature-table
parsing, record binding, location parsing, coordinate and coverage checking,
normalization, deterministic emission, and sealing. The CLI and Python API are
invocation surfaces for this compiler contract.

The feature-state adapter accepts that closed GenBank artifact and emits a
reference-only index over nine sealed child artifacts. Its development module
contains one constant unit per ordered BioIR feature and eleven deterministic
structural columns. The bound GenBank BioIR remains the semantic authority.

Artifact behavior is selected only by explicit format, version, profile, and
source syntax. Molecular identifiers and content values are ordinary input
data.

## 3. Common JSON and digest rules

Every artifact is a closed JSON object. Missing and unknown members are errors
at every defined level. Decoders MUST reject:

- duplicate object keys;
- non-UTF-8 JSON;
- non-finite numbers;
- lone Unicode surrogates;
- integers outside `[-(2^53-1), 2^53-1]`; and
- booleans where an integer is required.

`artifact_sha256` is lowercase SHA-256 over the [RFC 8785 JSON Canonicalization
Scheme](https://www.rfc-editor.org/rfc/rfc8785) encoding of the complete object
with `artifact_sha256` omitted. `bio_ir_sha256`, `record_ir_sha256`, and
`report_sha256` apply the same rule to their specified bodies.

Array order is semantic unless a field is explicitly declared sorted. JSON
object member order and pretty-printing whitespace are not semantic.

## 4. GenBank source language

### 4.1 Source envelope

The source MUST be 1 to 16 MiB of ASCII bytes and contain at most 1,000,000
physical lines. Lines MAY use LF or CRLF and MUST NOT contain bare carriage
returns, tabs, control bytes, or non-ASCII bytes. One line is at most 1 MiB.

The source contains 1 to 40,000 ordered records. Each record begins with
`LOCUS` in its specified field and ends with exact `//`. At most one empty
physical separator line may follow a terminator. Leading source material,
nonempty separators, a missing terminator, and trailing material are invalid.

### 4.2 Physical DNA record

The LOCUS declaration MUST describe a positive `bp` length no greater than the
16 MiB aggregate-base ceiling. The admitted molecule forms are `DNA`, `ss-DNA`,
`ds-DNA`, and `ms-DNA`; topology is `linear` or `circular`; division is three
uppercase ASCII letters; and the date is a real Gregorian `DD-MON-YYYY` date.
This profile compiles physical ORIGIN records. CONTIG-only records are rejected.

The mandatory header order is:

1. `LOCUS`;
2. `DEFINITION`;
3. `ACCESSION`;
4. `VERSION`;
5. optional unique `DBLINK` and `PROJECT` blocks;
6. `KEYWORDS`;
7. optional singleton `SEGMENT` as one physical line containing `n of m`, where
   `1 <= n <= m` and `m >= 2`;
8. `SOURCE` and `ORGANISM`;
9. one or more ordered `REFERENCE` blocks, followed by optional unique
   `COMMENT` and `PRIMARY` blocks;
10. exact `FEATURES             Location/Qualifiers`; and
11. `ORIGIN`, sequence lines, and `//`.

Each REFERENCE is consecutively numbered from one and contains exactly one
JOURNAL subkeyword in its fixed field. Wrapped header values are joined in
source order and retain positional source maps.

ACCESSION values are unique uppercase accessions. VERSION contains the primary
accession plus a positive numeric version and MAY contain one legacy `GI:`
decimal. The primary VERSION identity is at most 64 ASCII bytes, equals the
first ACCESSION plus its version, and is unique in a multi-record input.

### 4.3 Features and qualifiers

Feature keys occupy the GenBank FEATURES key field. Locations and qualifiers
begin in their defined fixed columns. One feature row is at most 80 bytes; the
complete source remains subject to the general line ceiling. A qualifier is a
name with either no value or a bounded unquoted or quoted value. Doubled quotes
inside quoted values decode to one literal quote. Continuation text is joined
deterministically; translation continuations are joined without spaces.

The frontend retains feature order, qualifier order and duplicates, joined
logical location text, a location AST, normalized segments, unresolved remote
references, and positional source maps.

Every record contains at least one local `source` feature. The union of source
intervals MUST cover the complete record; overlaps are allowed and gaps are
not. Each source feature has exactly one nonempty `/organism`. `/mol_type` may
occur at most once on each source feature; when present in the record it is
uniformly present and equal across all source features.

### 4.4 Location language

The location parser admits these structural forms:

- one positive 1-based point;
- a closed span `a..b`, with `<` on the start and `>` on the end where used;
- a between-base boundary `a^b`;
- the archived uncertain-point form `a.b`;
- an exact-version remote prefix `ACCESSION.version:` on an atomic location;
- `complement(location)`;
- `join(location,location,...)`; and
- `order(location,location,...)`.

Compound locations contain at least two children. Nested compounds,
join/order mixing, double complement, leading-zero coordinates, invalid fuzz
placement, and coordinates above the I-JSON safe-integer ceiling are rejected.

Local coordinates are checked against the record length. A same-input remote
reference is checked against the referenced record length and topology. An
external exact-version remote reference is preserved with
`bounds_status: "unresolved"`; no accession lookup occurs during compilation.
Where bounds and topology are known, between-base coordinates are adjacent,
with `record_length^1` additionally valid for a circular record. An external
remote `a^1` form requires explicit reference length and circular topology;
without them it is rejected rather than guessed.

Each location emits its source AST plus ordered 0-based half-open segments.
Segments retain reference, orientation, interval kind, fuzz, and bounds status.

### 4.5 ORIGIN

ORIGIN sequence lines contain a right-aligned 1-based index followed by one to
six groups of lowercase IUPAC DNA. Groups contain ten bases except the final
group on a line; nonfinal sequence lines contain 60 bases. Indices are
continuous. The emitted sequence is uppercase and its length equals LOCUS.

The admitted alphabet is:

```text
ACGTRYSWKMBDHVN
```

## 5. Sequence Collection v2

The collection wire identity is:

```json
{"format":"brain01.sequence-collection-ir","version":2}
```

Its root members are exactly `format`, `version`, `compiler`, `inputs`,
`members`, `refget_seqcol`, and `artifact_sha256`.

`inputs` identifies the GenBank profile and stores the exact original bytes,
byte length, and SHA-256. Source and normalized sequence use the same storage
descriptor:

```json
{
  "kind": "chunked-inline-ascii",
  "version": 1,
  "chunk_bytes": 1048576,
  "chunks": ["..."]
}
```

Every nonfinal chunk is exactly 1 MiB and the final chunk is 1 byte to 1 MiB.
The source is stored once. Each normalized record sequence is stored once.

One ordered member corresponds to one GenBank record and has exactly:

- `record_id`: the record VERSION;
- `source_map`: the inclusive physical-line range of its ORIGIN sequence; and
- `sequence`: alphabet, uppercase normalization, base count, SHA-256, refget
  identifier, and chunk storage.

Members are unique and cover every LOCUS record exactly once in source order.
The source map MUST replay the stored sequence and bind the same VERSION and
LOCUS length.

Sequence identity follows [GA4GH refget Sequences
v2.0.0](https://ga4gh.github.io/refget/sequences/) using the `SQ.`
SHA-512/24URL identifier. Collection identity follows [GA4GH refget Sequence
Collections v1.0.0](https://ga4gh.github.io/refget/seqcols/). Names and sequence
identifiers are inherent; lengths are retained as a non-inherent collated
attribute.

## 6. Structural BioIR v2

The outer artifact wire identity is:

```json
{
  "format": "brainc.bio.insdc-genbank-ir",
  "version": 2,
  "profile": "genbank-273-traditional-dna-physical-structural/v2"
}
```

Its root contains exactly the identity above plus `authority`, `compiler`,
`sequence_collection`, `bio_ir`, `bio_ir_sha256`, and `artifact_sha256`.

The BioIR body contains the profile, authority, embedded Sequence Collection
artifact digest, coordinate-system declaration, and ordered records. Each
record contains exactly:

- record identity, accessions, version, and optional legacy GI;
- parsed LOCUS values;
- definition, keywords, source, organism name, and lineage;
- ordered normalized header blocks;
- ordered features;
- explicit relationships;
- record line range; and
- `record_ir_sha256`.

A feature contains `feature_id`, source ordinal, key, location, qualifiers, and
line range. `feature_id` is a compiler-private opaque content identity with the
prefix `io.github.dttdrv.epigenesis.insdc-feature.sha256.` followed by the
lowercase SHA-256 digest of record ID, ordinal, key, retained location text,
and ordered qualifier values. It is not an INSDC accession, URI/CURIE, or GA4GH
computed identifier.

Source maps in BioIR contain only physical line and column positions. Exact
source text is reconstructed from the single source store in Sequence
Collection v2. BioIR records do not duplicate raw source chunks or normalized
sequence strings.

## 7. Producer validation

`validate_genbank_artifact` validates the closed artifact, canonical seals,
embedded Sequence Collection, collection-to-BioIR binding, and complete replay
from the source bytes stored inside the artifact. When original bytes are
supplied, they MUST equal the embedded source exactly.

Producer validation is suitable for loading a compiler artifact within the
compiler package. Independent attestation uses section 8.

## 8. Independent validator

`brainc.validator_insdc` is a standard-library implementation that imports no
compiler module. Given the original GenBank bytes and compiled artifact it:

1. enforces independent byte, tree, and resource bounds;
2. parses the GenBank record grammar;
3. parses and normalizes INSDC locations;
4. reconstructs Sequence Collection v2 from exact ORIGIN lines;
5. recomputes SHA-256, refget Sequence, and refget Sequence Collection
   identities;
6. reconstructs position-only source maps, record IR, BioIR, and artifact;
7. verifies RFC 8785 seals and compares complete replay trees with exact JSON
   type equality; and
8. emits a sealed validation report only after exact equality succeeds.

The report wire identity is:

```json
{"format":"brainc.bio.insdc-genbank-validation-report","version":2}
```

It binds source, collection, BioIR, compiled artifact, replay summary, validator
identity, and report digest. A valid report always contains `valid: true`.
Invalid evidence raises an error and cannot be represented as a success report.

The integrated CLI returns 0 on success, 2 for compilation or argument errors,
and 3 for validation failure. The isolated validator returns 0 on success and 1
on validation failure.

## 9. INSDC feature-state adapter

### 9.1 Input and projection

The adapter accepts only a valid `brainc.bio.insdc-genbank-ir/v2` artifact with
profile `genbank-273-traditional-dna-physical-structural/v2` and authority
manifest SHA-256
`6dc3658a5f7d6a774ce9cd08ff6dc5327cbef8a9951c7679fd9a34650f912bf9`.
Its source binding contains the complete GenBank artifact SHA-256 and BioIR
SHA-256.

One feature-table entry becomes one `feature` unit. Unit order is BioIR record
array order followed by one-based feature ordinal, and unit identity is the
source `feature_id`. The target contract contains one `unit.create` operation,
constant scalar fields, and empty edge schemas, ports, and rules. The adapter
MUST emit no inferred relationships.

The development module is a structural projection. Exact location ASTs,
segment associations, qualifiers, sequence, and provenance are resolved from
the content-bound GenBank BioIR semantic authority.

### 9.2 Feature columns

Every column is a rank-one little-endian `u64` tensor with one value per unit.
The feature-count tensor is a unitless `u64` scalar.

| Column | Normative value |
|---|---|
| `fuzzy_boundary_count` | Count of non-null `start_fuzz` and `end_fuzz` values; archived within-position uncertainty is excluded |
| `key_code` | Zero-based code into distinct feature keys sorted by ascending UTF-8 bytes |
| `kind_mask` | Bitwise OR over segment kinds: interval `1`, between `2`, uncertain-point `4` |
| `ordinal` | One-based source feature ordinal |
| `orientation_mask` | Bitwise OR over segment orientations: forward `1`, reverse/complement `2` |
| `qualifier_count` | Length of the ordered qualifier array, including repeated qualifiers |
| `record_code` | Zero-based code into distinct record IDs sorted by ascending UTF-8 bytes |
| `remote_segment_count` | Count of segments with a non-null remote reference |
| `segment_count` | Length of the normalized segment array |
| `segment_extent_sum` | Sum of `end - start` across normalized 0-based half-open segments; between segments contribute zero, uncertain-point segments contribute their normalized envelope width, and overlapping extents count with multiplicity |
| `unresolved_segment_count` | Count of segments whose `bounds_status` is `unresolved` |

Dictionary codes are bundle-local and MAY renumber when dictionary membership
changes. Kind and orientation masks record presence; they do not retain
multiplicity, order, or segment association. `unresolved_segment_count` records
reference-resolution and bounds status rather than INSDC validity.

### 9.3 Artifact closure

The immutable backend is identified by:

```json
{
  "id": "io.github.dttdrv.epigenesis.bio.insdc-graph",
  "version": 1,
  "spec_sha256": "a8f2494dddd259f7f41604001f5dcae9e19b0e491eac7c13aec298f2f99a5da4"
}
```

The serialized `brainc.bio.insdc-graph-bundle/v1` object is a reference-only
index manifest. It binds the GenBank source dependency and exactly these nine
sealed children by `artifact_sha256`:

1. `backend_spec`;
2. `backend_semantics`;
3. `compilation_record`;
4. `development_module`;
5. `lowering_policy`;
6. `prediction_request`;
7. `prediction_response`;
8. `provider_manifest`; and
9. `target_contract`.

Directory publication writes `bundle.json` and one `<role>.json` file for every
child. The index establishes identity and integrity; validation receives the
bound GenBank artifact as the semantic source dependency.

### 9.4 Independent replay

A producer-independent validator MUST validate the compiled GenBank artifact,
reconstruct the dictionaries, unit order, all eleven columns, tensor bytes,
generic compiler chain, compilation record, child seals, and index. Exact replay
MUST reject direct or coherently resealed substitutions. A successful report is
`brainc.bio.insdc-graph-validation-report/v1` and binds the GenBank artifact,
BioIR, backend specification, bundle index, semantics, and development module.

The public replay interfaces are:

```python
validate_insdc_graph_bundle(
    bundle_index, artifacts_by_role, genbank_artifact, *, genbank_source=None
)
validate_insdc_graph_paths(
    genbank_path, genbank_artifact_path, bundle_directory, *, report_path=None
)
bundle_index, artifacts_by_role = load_insdc_graph_bundle_directory(path)
validate_insdc_graph_report(report)
```

The path loader accepts exactly `bundle.json` and the nine `<role>.json` child
files. The isolated command is:

```text
brainc-validate-insdc-graph ORIGINAL.gb GENBANK_ARTIFACT.json BUNDLE_DIRECTORY \
  [-o|--output|--report REPORT.json]
```

## 10. Resource profile

| Resource | Maximum |
|---|---:|
| Original GenBank source | 16 MiB |
| Emitted GenBank artifact | 64 MiB |
| Serialized validation report | 16 MiB |
| One feature-state child or index | 16 MiB |
| Total normalized bases | 16 MiB |
| Physical lines | 1,000,000 |
| One line / decoded string / location | 1 MiB |
| Records | 40,000 |
| Features across the source stream | 20,000 |
| Qualifiers across the source stream | 100,000 |
| Accessions per record | 10,000 |
| JSON nesting depth | 64 |
| JSON members | 1,000,000 |
| Location nesting depth | 32 |
| Location AST nodes | 300,005 |
| Compound location children | 100,000 |
| VERSION or remote identity | 64 ASCII bytes |
| Canonical storage chunk | 1 MiB |

These are independent admission ceilings and may bind in combination; meeting
one limit does not relax another. A host MAY apply tighter limits and MUST
reject before exceeding them.

## 11. Output publication

On descriptor-capable POSIX hosts, artifact and report publication uses bounded
same-parent exclusive staging, staged-file synchronization, descriptor-relative
final-component checks with `O_NOFOLLOW`, atomic replacement, published-file and
parent identity verification, and parent-directory synchronization.

The portable artifact fallback retains bounded same-parent staging,
staged-file synchronization, final-component link checks, and atomic
replacement. The portable validator fallback additionally rechecks the
published file and parent identities. Portable operation requires a
caller-controlled output parent that is not concurrently renamed or modified;
parent-directory synchronization is a POSIX-profile guarantee.

Feature-state publication writes `bundle.json` and the nine child files into a
same-parent staging directory, verifies their snapshots, synchronizes them, and
publishes the complete directory with one rename. The destination MUST NOT
already exist. On descriptor-capable POSIX hosts, the compiler pins the
non-linked parent and staging identities, creates regular files
descriptor-relatively with `O_EXCL` and `O_NOFOLLOW`, synchronizes both
directories, verifies the published identity and contents, and removes only
same-identity objects it created if publication fails.

## 12. Compatibility

Package 0.9.0 retains these independent versioned paths:

| Capability | Wire contracts | Commands |
|---|---|---|
| Sequence frontend | Sequence IR v2; Sequence Collection v1 | `compile-sequence`, `compile-collection`, `compile-raw` |
| Scalar lowering | 0.4 provider and state-program ABIs | `make-request`, `compile`, `check`, `validate` |
| Tensor/target lowering | 0.5 provider v2, target v1, lowering v2, development module v1 | `make-request-v2`, `compile-v2`, `check-v2`, `validate-v2` |
| External GFF3 | GFF3 BioIR v1 on Sequence Collection v1 | `compile-gff3`, `validate-gff3` |
| Feature graph | Feature-graph bundle v1 | `compile-feature-graph`, `validate-feature-graph` |
| INSDC feature state | Feature-state backend, semantics, compilation, bundle, and development-module v1 artifacts | `compile-insdc-graph`, `brainc-validate-insdc-graph` |

Sequence Collection v2 is a GenBank 0.8 contract. Existing tensor, GFF3, and
feature-graph paths continue to declare Sequence Collection v1 and MUST NOT
silently reinterpret v2.

## 13. Frontend extension contract

Additional biological languages enter as separately versioned frontends. A
publishable frontend defines:

1. a pinned primary authority;
2. a closed grammar and normalized IR;
3. exact source and cross-artifact identities;
4. explicit coordinate, reference, missing-value, and ordering semantics;
5. finite admission limits;
6. data-independent dispatch;
7. a declared target-lowering boundary; and
8. an independently implemented replay validator.

This contract permits heterogeneous DNA and biological model inputs while
retaining the semantics native to each source language.

Licensed under Apache-2.0.

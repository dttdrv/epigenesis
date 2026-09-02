# Epigenesis universal DNA translation compiler contract 1.0.0

## 1. Status and terminology

This document specifies Epigenesis 1.0, a universal DNA translation compiler:
profile-qualified biological source admission, executable external grammar
extension, exact source provenance, a common typed source boundary, Development
Module lowering, reference-only linking, and producer-isolated verification.

The interpretation manifest, request, and response are caller-supplied. This
contract does not specify biological prediction or derive an interpretation
from source DNA.

The words MUST, MUST NOT, REQUIRED, SHALL, SHALL NOT, SHOULD, SHOULD NOT,
RECOMMENDED, NOT RECOMMENDED, MAY, and OPTIONAL are interpreted as described by
[BCP 14](https://www.rfc-editor.org/info/bcp14) when capitalized.

An Epigenesis compilation has five inputs:

1. a source bundle produced under one exact biological source profile;
2. an interpretation manifest;
3. a source-bound interpretation request and response;
4. a lowering policy; and
5. a target contract.

Its primary output is `brainc.development-module/v1`. The linker publishes that
module with the complete interpretation chain and a compilation record under
`brainc.development-bundle/v1`. A producer-isolated implementation MUST reparse
built-in source profiles. For a version-2 executable external profile, it MUST
rerun the digest-pinned validator against fresh snapshots of the original
source and embedded profile-native IR.

## 2. Common wire rules

Every compiler contract artifact is a closed JSON object. Referenced source
bytes and external blobs are content-addressed byte resources. Missing or
unknown members are errors at every defined JSON level. A decoder MUST reject:

- duplicate object keys;
- non-UTF-8 JSON;
- non-finite numbers;
- lone Unicode surrogates;
- integers outside `[-(2^53-1), 2^53-1]`;
- booleans or floats where an integer is required; and
- values beyond the artifact's depth, member, string, or byte ceiling.

`artifact_sha256` is lowercase SHA-256 over the [RFC 8785 JSON Canonicalization
Scheme](https://www.rfc-editor.org/rfc/rfc8785) encoding of the complete object
with `artifact_sha256` omitted. Fields ending in `_ir_sha256`,
`catalog_sha256`, `module_sha256`, and `report_sha256` use the same canonical
encoding over their named body.

Array order is semantic unless a field explicitly declares a sorted set. JSON
object member order and presentation whitespace are not semantic. Validation
MUST compare exact JSON types; a matching canonical digest is not a type oracle.

## 3. Source compilation

### 3.1 Explicit dispatch

The source profile is a required caller input. Source content, path, extension,
record identifier, accession, organism, gene, digest, and record count MUST NOT
select or alter the frontend.

The built-in profiles are:

| Profile | Required input roles | Parameters | Native roles |
|---|---|---|---|
| `raw-iupac-dna/v1` | `sequence` | `record_id` | `sequence` |
| `fasta-dna/v1` | `sequence` | `wrapper` = `identity` or `gzip` | `sequence` |
| `fasta-reference-dna/v1` | `sequence` | `wrapper` = `identity` or `gzip` | `reference` |
| `genbank-273-traditional-dna-physical-structural/v2` | `genbank` | none | `genbank` |
| `gff3-external-sequence/v1` | `sequence`, `annotation` | nested raw or FASTA sequence profile | `sequence`, `annotation` |
| `external-dna-source/v1` | declared by manifest | exact external profile identity | `external` |

Unknown roles, parameters, wrappers, profiles, or native children MUST fail
before frontend use. The external profile is admitted through its dedicated
evidence API and command; it is not selected by sniffing a byte stream.

### 3.2 Source descriptor

The canonical typed source boundary is the composition of the common source
descriptor, its `source_ir.records` array, and the selected profile's
profile-native typed artifacts. The common array preserves ordered DNA record
identity. External source admission normalizes its
`brainc.sequence-record-catalog/v1` into this array. Native artifacts preserve
source-language semantics that cannot be represented without loss there.

Every accepted source emits `brainc.source-descriptor/v1`:

```text
format
version
producer
source_ir
  profile
  parameters
  inputs
  artifacts
  records
source_ir_sha256
artifact_sha256
```

Each `inputs` member contains exactly `sha256` and `byte_length`. Each
`artifacts` member contains exactly `format`, `version`, `artifact_sha256`, and
`ir_sha256`. Each ordered `records` member contains exactly `record_id`,
`bases`, `sequence_sha256`, and `refget_id`.

The descriptor MUST equal the native artifact closure:

- input roles and byte identities equal native provenance;
- native roles, formats, versions, artifact identities, and IR identities
  equal the bundle children;
- record order, identifiers, lengths, sequence digests, and refget identities
  equal the native record catalog; and
- normalized parameters equal the selected profile's exact parameter object.

A source bundle contains `source.json` plus only the filenames selected below:

| Native role | Filename |
|---|---|
| `sequence` | `sequence.json` |
| `reference` | `reference.json` |
| `external` | `external.json` |
| `genbank` | `genbank.json` |
| `annotation` | `annotation.json` |

The directory MUST contain no other file. Every child MUST be one stable,
regular, non-symbolic, single-link file and MUST remain unchanged for the
complete read transaction.

### 3.3 Raw and FASTA sequence profiles

Raw and FASTA sequence frontends accept the IUPAC DNA alphabet
`ACGTRYSWKMBDHVN` case-insensitively and normalize sequence letters to upper
case. FASTA identifiers are taken from the first nonempty defline token,
preserve record order, and MUST be unique. A selected wrapper MUST match the
physical source exactly.

The exact-source Sequence Collection IR binds the original bytes, normalized
sequence, record order, SHA-256, and GA4GH refget Sequence identities. The
legacy `fasta-dna/v1` route retains its established bounded in-memory contract.

### 3.4 Streaming reference FASTA

`fasta-reference-dna/v1` streams the physical source and emits
`brainc.reference-sequence-catalog/v1`. It MUST NOT place sequence strings in
the catalog or source descriptor.

The catalog binds:

- exact physical source SHA-256 and byte length;
- exact logical FASTA SHA-256 and byte length;
- explicit `identity` or `gzip` wrapper;
- ordered record identifiers, base counts, uppercase-sequence SHA-256, and
  GA4GH refget Sequence identifiers; and
- a GA4GH refget Sequence Collection identity over ordered names and sequences,
  with ordered lengths retained as a collated attribute.

Identity wrapping requires identical physical and logical references. Gzip
decoding accepts a bounded sequence of complete concatenated RFC 1952 members
and rejects leading, inter-member, and trailing bytes, padding, truncation,
corrupt checksums, and a member count beyond the execution budget.

LF, CRLF, bare CR, and mixed line endings have one declared grammar. A `>` byte
starts a record only at a logical line start; an interior `>` is invalid DNA.
Chunk boundaries MUST NOT affect the result.

Default execution budgets are 64 GiB physical and 64 GiB logical, 100,000
records, 2,147,483,647 bases per record, 1 MiB per defline, and 100,000 gzip
members. Callers MAY set tighter limits or raise physical/logical/member budgets
up to their published hard maxima. Execution limits do not alter accepted
artifact identity.

### 3.5 GenBank physical DNA and GFF3

The GenBank profile is pinned to NCBI GenBank Flat File Release 273.0 and INSDC
Feature Table Definition 11.4 by packaged authority manifest. It compiles
traditional physical-DNA records into exact Sequence Collection v2 and
`brainc.bio.insdc-genbank-ir/v2`, preserving ordered records, headers,
references, feature keys, qualifiers, location ASTs, normalized coordinates,
ORIGIN sequence, and source positions.

The GFF3 profile is a bounded structural profile derived from Sequence Ontology
GFF3 1.26 and compiled under an explicitly supplied raw or FASTA sequence
profile. It validates sequence-region bounds, feature coordinates, circular
origin crossing, phase, score, attributes, ordered relationships, reference
closure, and exact binding to the sequence collection. Feature types and
ontology claims remain opaque. Embedded FASTA is outside this profile.

## 4. Universal external frontend ABI

The external ABI admits any DNA-bearing source grammar through an explicitly
selected frontend and separate validator. Its current contracts are:

- `brainc.external-frontend-profile/v1` and `/v2`;
- `brainc.external-source-descriptor/v1`;
- `brainc.external-frontend-validation/v1`; and
- compiler closure `brainc.external-source-closure/v1` and `/v2`.

### 4.1 Profile manifest

Every manifest MUST declare:

- ABI `brainc.external-frontend-data/v1` for legacy data-only admission or
  `brainc.external-frontend-executable/v2` for executable admission;
- portable profile identifier and positive integer version;
- grammar identifier, version, authority URI, and authority SHA-256;
- one to 64 ordered input-role declarations with media type, wrapper, and
  maximum byte length;
- one to 64 ordered native JSON artifact declarations with media type, format,
  version, schema SHA-256, IR-digest field, and maximum byte length;
- record catalog schema `brainc.sequence-record-catalog/v1`;
- total input/native byte, record, base, and portable record-identifier limits,
  with record identifiers capped at 256 UTF-8 bytes; and
- validator protocol, distribution identity, distribution digest, executable
  identity, and executable digest; and
- for version 2, a frontend protocol and command plus `native` or `python`
  runtime declarations for both commands.

Version-2 frontend and validator executable digests MUST differ. Executable
paths are explicit caller arguments and MUST NOT be discovered from manifest
content. The implementation MUST verify a stable regular single-link program
snapshot against the declared SHA-256 and MUST invoke it without a shell.

The complete manifest digest is the external language identity. Every accepted
external source MUST provide a nonempty sequence-record catalog with positive
base counts, sequence SHA-256 values, and refget identities. Overlay formats
MUST enter through a composite profile that also binds their exact DNA-bearing
input. A provider accepting an external source MUST declare:

```text
brainc.source-descriptor/v1;profile=external-dna-source/v1;manifest=<sha256>
```

A generic external-profile acceptance tag is insufficient.

### 4.2 Executable frontend and replay

The frontend descriptor MUST bind the manifest, exact input roles and byte
identities, exact native roles/formats/versions/schemas/IR identities and byte
identities, and an ordered record catalog. Record identifiers MUST be unique;
ordinals MUST be contiguous from zero; every record MUST name a declared input
role and carry positive bases, sequence SHA-256, and refget identity.

For version 2, Epigenesis MUST snapshot every original source through a stable
file descriptor while enforcing role and cumulative byte limits. It passes
those private snapshots, exact references, and the sealed manifest to the
frontend over closed JSON on standard input. The frontend returns profile-native
JSON and ordered record entries. Epigenesis derives native byte references,
checks declared format/version and IR-digest fields, and seals the descriptor.

Epigenesis then invokes the separately digest-pinned validator with the same
snapshots, manifest, descriptor, and native JSON. The validator MUST reparse the
grammar, validate the native representation, and return exact input references,
native references, and records. Any difference fails admission. Execution MUST
have a finite wall timeout and bounded standard output and error. Nonzero exit,
standard-error output, malformed or duplicate-key JSON, unknown fields, or a
resource violation MUST fail before publication.

Version-2 closure artifacts embed the bounded profile-native JSON with the
manifest, descriptor, and replay report. Original source bytes remain external
and content-addressed. The final Development Module validator MUST repeat the
validator execution against fresh source snapshots and require the exact sealed
result.

Version-1 artifacts remain compatible. Their caller-attested data-only path
continues to verify exact original and separately supplied native bytes, but
does not establish grammar execution.

### 4.3 Executable protocol envelopes

Each command reads exactly one UTF-8, duplicate-free JSON object from standard
input and writes exactly one JSON object to standard output. It MUST
write nothing to standard error on success. Paths are ephemeral absolute paths
to compiler-owned source snapshots; a command MUST read them during that
invocation and MUST NOT persist them as source identities. The `inputs` array
is in manifest role order and has this exact member shape:

```text
{role, sha256, byte_length, path}
```

The frontend receives:

```text
format = brainc.external-frontend-execution-request
version = 1
profile_manifest = complete sealed v2 profile manifest
inputs = [{role, sha256, byte_length, path}, ...]
```

It returns:

```text
format = brainc.external-frontend-execution-result
version = 1
profile_manifest_sha256 = profile manifest artifact_sha256
native_artifacts = {declared-role: complete profile-native JSON, ...}
records = [{ordinal, input_role, record_id, bases,
            sequence_sha256, refget_id}, ...]
```

The compiler validates and seals that result, then gives the validator:

```text
format = brainc.external-validator-execution-request
version = 1
profile_manifest = complete sealed v2 profile manifest
source_descriptor = compiler-sealed external source descriptor
inputs = [{role, sha256, byte_length, path}, ...]
native_artifacts = {declared-role: complete profile-native JSON, ...}
```

The validator returns:

```text
format = brainc.external-validator-execution-result
version = 1
profile_manifest_sha256 = profile manifest artifact_sha256
source_descriptor_sha256 = source descriptor artifact_sha256
input_references = [{role, sha256, byte_length}, ...]
native_artifact_references =
  [{role, format, version, schema_sha256, sha256, byte_length, ir_sha256}, ...]
records = [{ordinal, input_role, record_id, bases,
            sequence_sha256, refget_id}, ...]
valid = true
```

Every protocol object shown above is closed. Profile-native object shape is
defined by its declared schema and enforced by the external validator. The
validator result MUST be exactly equal by JSON type and value to the references
and records derived and sealed by the compiler; `true` is not interchangeable
with `1`. A conforming frontend or validator MAY be implemented in any language
that can honor its manifest runtime and these process contracts.

## 5. Interpretation and target contracts

The interpretation chain uses these existing versioned artifacts:

| Role | Wire identity |
|---|---|
| Manifest | `brainc.provider-manifest/v2` |
| Request | `brainc.prediction-request/v2` |
| Response | `brainc.prediction-response/v2` |
| Lowering policy | `brainc.lowering-policy/v2` |
| Target contract | `brainc.target-contract/v1` |
| Result | `brainc.development-module/v1` |

The manifest MUST accept the exact source-language tag and declare ordered
output contracts. The request MUST bind the exact source and manifest. The
response MUST bind the request, provider identity, configuration identity, and
ordered output contracts.

Every tensor is checked for dtype, shape, axes, unit, storage kind, byte length,
digest, canonical inline encoding or explicit external blob reference, and
coordinate compatibility with the validated source records.

The target contract declares ABI major version, opsets, unit and edge schemas,
field types, numeric contracts, ports, and rules. The lowering policy MUST bind
that exact target, cover every required response output once, and emit only
operations authorized by the declared opsets and capabilities.

The compiler MUST validate operation order and references, schema/field
compatibility, tensor use, count sources, unit and edge endpoints, rule
subjects, phases, parameters, reads, writes, port closure, conflicts, and exact
operation/tensor/unit/edge/attachment budgets before emitting a module.

## 6. Development bundle and linker

`brainc.development-bundle/v1` is a flat reference-only index. Its directory
contains exactly:

```text
bundle.json
compilation_record.json
development_module.json
lowering_policy.json
prediction_request.json
prediction_response.json
provider_manifest.json
target_contract.json
```

The index binds the source descriptor, source-native artifact references,
seven child artifact references, and the sorted unique set of external tensor
blobs. It MUST NOT embed source DNA, source-native bodies, child bodies, or blob
bytes.

`brainc.development-compilation/v1` binds the exact source, all chain children,
the final `module_sha256`, and exact module budgets. The compilation is
deterministic: identical accepted inputs and execution-independent limits MUST
produce byte-identical files.

Publication is absent-to-complete and no-clobber. A destination that exists at
commit time MUST be preserved. A committed directory MUST never be removed due
to a later synchronization or verification error.

## 7. Producer-isolated validation

`brainc.development-validation-report/v1` is produced by an implementation
outside the producer module boundary. It MUST:

1. snapshot the exact source and development directory closures;
2. replay each built-in original source through a separately implemented
   frontend; for a version-2 external source, verify and execute the pinned
   validator over fresh source snapshots and embedded native JSON; for version
   1, verify the exact original/native paths and caller attestation;
3. validate descriptor/native/input/record equality;
4. replay manifest, request, response, tensors, target, policy, operations, and
   budgets;
5. reconstruct the exact Development Module and compilation record;
6. validate every unique external blob once and enforce cumulative budgets;
7. reconstruct the exact reference-only bundle; and
8. seal a deterministic success or failure report.

The validator MUST NOT import producer source, bundle, provider, compiler,
target, policy, tensor, canonicalization, or I/O modules. It MAY compose
separately implemented validators that share the same trust boundary.

Validation paths MUST reject symbolic links, multiply linked files, nonregular
files, unknown or missing children, duplicate JSON keys, path replacement,
directory replacement, oversized children, decompression overflow, allocation
failure, and coherently resealed semantic changes.

Success exits zero. After successful command-line parsing, semantic validation
failure exits one and emits a sealed report with `valid: false` and code
`DEVVAL001`. Command-line usage errors follow `argparse` and exit two.

The standalone command receives original roles as repeated
`--source-input ROLE=PATH` arguments. Version-1 external profiles additionally
receive every native role through `--native-input ROLE=PATH`; version-2
profiles receive the explicit validator through `--external-validator PATH`.

## 8. Diagnostics

Public compiler failures are normalized into stable domains:

| Code | Domain |
|---|---|
| `SOURCE001` | source dispatch, descriptor, closure, or publication |
| `SCALE001` | streaming reference FASTA admission |
| `EXT001` | external frontend manifest/descriptor/report |
| `EXTSRC001` | compiler-side external source closure |
| `DEVB001` | development compilation or bundle linking |
| `REFVAL001` | independent reference replay |
| `DEVVAL001` | independent complete-compilation replay |

Malformed user input MUST NOT escape as raw `TypeError`, `KeyError`,
`AssertionError`, `MemoryError`, or traceback. Internal programming errors are
not converted to successful diagnostics.

## 9. Resource and path invariants

All public input and output paths are bounded before untrusted allocation.
Accepted paths MUST resolve to one stable regular object of the required kind.
Secure directory publication and validation use descriptor-relative primitives
on Linux and macOS. Windows uses platform-specific identity checks for focused
file and publication operations, but complete secure bundle validation is not
part of the supported 1.0 host contract. A platform without the required
primitives MUST fail closed for that operation.

Every reader enforces per-role bytes, JSON depth, JSON member count, string
length, record count, and relevant decompressed/logical limits. Whole-chain
validation additionally enforces external blob count and cumulative bytes.
Allocation and recursion failures at public validation boundaries become stable
validation failures.

## 10. Compatibility

Package 1.0 retains these independently versioned paths:

| Capability | Commands |
|---|---|
| Sequence IR and Collection IR | `compile-sequence`, `compile-collection`, `compile-raw` |
| Scalar state-program lowering | `make-request`, `compile`, `check`, `validate` |
| Tensor/target lowering | `make-request-v2`, `compile-v2`, `check-v2`, `validate-v2` |
| External GFF3 BioIR | `compile-gff3`, `validate-gff3` |
| Structural feature graph | `compile-feature-graph`, `validate-feature-graph` |
| INSDC feature state | `compile-insdc-graph`, `brainc-validate-insdc-graph` |
| Final source/development compiler | `compile-source`, `compile-external-source`, `compile-development`, `brainc-validate-development` |

Package versions, producer versions, wire-format versions, frontend-profile
versions, and cited standards versions are separate namespaces. A wire-format
version changes only when that contract changes. A cited standards version is
not an Epigenesis release version.

## 11. Conformance

A conforming distribution MUST pass:

- complete unit and integration tests;
- universal translation convergence across every built-in profile and a real
  executable external FASTQ frontend/validator pair for equivalent DNA;
- rejection of malformed external grammar input, changed executable/source
  bytes, lying frontends, timeouts, output overflow, and malformed process JSON;
- generated unseen-source causality and no-content-dispatch checks;
- direct and coherently resealed cross-artifact attacks;
- path, resource, gzip, JSON, and allocation attacks;
- independent validator import-isolation checks;
- an installed-wheel real-DNA compilation and validation;
- an offline wheel smoke test; and
- an extracted-sdist self-test outside the checkout.

The Apache-2.0 license is the complete distribution license; no NOTICE file is
required by this project.

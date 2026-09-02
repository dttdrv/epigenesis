# Epigenesis

Epigenesis 1.0 is the universal DNA translation compiler: DNA-bearing source
formats enter one deterministic, verifiable pipeline and emerge as typed
Development Module IR.

```text
DNA source/profile
    -> common typed source boundary
    -> validated interpretation and target lowering
    -> Development Module IR
    -> independently verifiable compilation bundle
```

Built-in frontends cover raw IUPAC DNA, FASTA, streaming reference FASTA,
GenBank, and GFF3 with an exact sequence. The data-only external frontend ABI
extends the same pipeline to additional DNA grammars through explicit profile,
schema, provenance, and validation-evidence contracts.

Every accepted source keeps its exact byte identity, ordered DNA records,
profile-native structure, grammar authority, and content-derived sequence
identities. The common typed source boundary makes those formats available to
one target-neutral lowering system without erasing their native meaning.

The interpretation stage is an explicit compiler input. Epigenesis binds the
caller-supplied interpretation to the exact source, checks its tensors and
coordinates against a target ABI, emits only target-authorized operations, and
seals the complete result with RFC 8785 canonical JSON and SHA-256.

The compiler and its producer-isolated validators use only the Python 3.11+
standard library.

## Compile a DNA source

Install the package and compile an assembled genome without materializing its
sequence in the output:

```bash
python -m pip install .

brainc compile-source \
  --profile fasta-reference-dna/v1 \
  --sequence genome.fasta \
  --wrapper identity \
  --output genome.source
```

`genome.source` contains a sealed source descriptor and a compact reference
catalog. The catalog records exact physical and logical source identities,
per-record SHA-256 and GA4GH refget identifiers, ordered collection identity,
and base counts. Sequence strings stay in the original source file.

Compile the validated source and a caller-supplied typed interpretation into
the final target:

```bash
brainc compile-development genome.source \
  --manifest interpretation/manifest.json \
  --request interpretation/request.json \
  --response interpretation/response.json \
  --policy interpretation/lowering-policy.json \
  --target interpretation/target-contract.json \
  --output genome.development
```

Independently replay the complete compilation from the original DNA:

```bash
brainc-validate-development genome.source genome.development \
  --source-input sequence=genome.fasta \
  --output genome.validation.json
```

A successful report has `valid: true` and binds the original DNA, source
profile, native artifacts, interpretation chain, target contract, lowering,
Development Module, bundle closure, external blobs, and report itself.

## Source languages

Source selection is explicit. Bytes, filenames, accessions, organisms, genes,
and fixture identities never choose a frontend.

| Profile | Source program | Native result |
|---|---|---|
| `raw-iupac-dna/v1` | One raw IUPAC DNA sequence | Exact Sequence Collection IR |
| `fasta-dna/v1` | FASTA or concatenated gzip FASTA | Exact Sequence Collection IR |
| `fasta-reference-dna/v1` | Streaming chromosome/genome FASTA | Compact refget sequence catalog |
| `genbank-273-traditional-dna-physical-structural/v2` | Physical-DNA GenBank records | Sequence authority plus structural INSDC BioIR |
| `gff3-external-sequence/v1` | GFF3 plus an exact external sequence source | Collection-bound structural GFF3 BioIR |
| `external-dna-source/v1` | Manifest-declared DNA-bearing source plus validation evidence | Manifest-qualified external source closure |

The reference FASTA frontend streams identity or strict concatenated-gzip
input. Its default physical and logical budgets are 64 GiB, with explicit
caller overrides up to the I-JSON safe-integer ceiling:

```bash
brainc compile-source \
  --profile fasta-reference-dna/v1 \
  --sequence exceptional-assembly.fasta.gz \
  --wrapper gzip \
  --maximum-input-bytes 137438953472 \
  --maximum-logical-bytes 274877906944 \
  --output exceptional.source
```

## External frontend data ABI

Additional DNA-bearing grammars and composite DNA-plus-overlay profiles can
enter through a data-only frontend ABI. A frontend publishes three sealed
objects:

1. a profile manifest naming its grammar, authority binding, input roles,
   native artifact schema identities, validator identity, and finite limits;
2. a source descriptor containing exact input/native references and an ordered
   sequence-record catalog; and
3. a manifest-qualified validation report reproducing that descriptor.

The ABI separates frontend execution from compiler admission. Epigenesis
verifies the supplied report, streams and hashes the declared original inputs,
validates bounded native JSON identities, binds the exact profile-manifest
digest, and emits the standard source descriptor:

```bash
brainc compile-external-source \
  --profile-manifest fastq.profile.json \
  --source-descriptor reads.frontend.json \
  --validation-report reads.validation.json \
  --source-input reads=reads.fastq \
  --native-artifact read-index=reads.index.json \
  --output reads.source
```

After lowering that source into a development bundle, the final validator
rehashes the original/native files and verifies the same closure explicitly:

```bash
brainc-validate-development reads.source reads.development \
  --source-input reads=reads.fastq \
  --native-input read-index=reads.index.json \
  --output reads.validation.json
```

The ABI is manifest-extensible and sequence-catalog-bound. Every external source
provides at least one positive-length DNA record with SHA-256 and refget
identity; annotation, alignment, variation, and assembly-graph overlays require
an exact DNA-bearing input in the same composite profile. Role names, grammar
identity, wrappers, native IR formats, limits, and validator distribution
identities come from the manifest. Two otherwise matching profiles remain
distinct when their authority or schema changes because provider acceptance is
qualified by the complete manifest SHA-256.

## What the compiler checks

Epigenesis performs linked static validation across the whole compilation:

- built-in grammar semantics, roles, wrappers, record order, sequence alphabet,
  coordinates, references, and native IR closure;
- external manifest, validation-report, original-byte, native-JSON, and record-
  catalog closure through supplied evidence;
- source SHA-256, byte lengths, sequence digests, refget identities, canonical
  seals, source maps where emitted, and authority digests where declared;
- interpretation manifest, request, response, configuration identity, ordered
  outputs, tensor types, shapes, axes, storage, and exact source binding;
- target ABI, opsets, unit and edge schemas, fields, rules, ports,
  capabilities, operation order, references, and numeric contracts;
- lowering-policy coverage, tensor use, state operations, conflicts, exact
  unit/edge/attachment/tensor budgets, and deterministic module identity; and
- the flat bundle's exact child set, transitive references, unique external
  blobs, compilation record, and final validation report.

The result is `brainc.development-module/v1`: typed Development Module IR
containing unit, edge, rule-attachment, and port operations for a declared
target. The compiler's product is the validated IR and its complete
provenance/link bundle.

## Output contracts

`compile-development` publishes exactly eight files:

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

`bundle.json` is a reference-only digest index over source/native artifacts,
seven child contracts, and external tensor blobs. Publication is atomic,
no-clobber, and absent-to-complete on supported Linux and macOS hosts. Windows
has a focused publication primitive, but the complete secure bundle-validation
path is not advertised there. Unsupported operations fail closed.

| Contract | Identity |
|---|---|
| Source descriptor | `brainc.source-descriptor/v1` |
| Reference catalog | `brainc.reference-sequence-catalog/v1` |
| External frontend manifest | `brainc.external-frontend-profile/v1` |
| External source closure | `brainc.external-source-closure/v1` |
| Target contract | `brainc.target-contract/v1` |
| Development module | `brainc.development-module/v1` |
| Compilation bundle | `brainc.development-bundle/v1` |
| Independent result | `brainc.development-validation-report/v1` |

Existing Sequence IR, GenBank BioIR, GFF3 BioIR, feature-graph, and INSDC
feature-state commands remain versioned and compatible.

## Independent validation

The producer-isolated validator independently parses canonical JSON, snapshots
exact directories and paths, replays the built-in source frontends, verifies
external frontend closures against their original/native evidence, reconstructs
typed lowering, checks blob contents, and compares the complete result using
exact JSON types.

Dedicated validators are also available:

```bash
brainc-validate-genbank source.gb compiled.genbank.json -o report.json
brainc-validate-insdc-graph source.gb compiled.genbank.json graph.bundle -o report.json
brainc-validate-development source.bundle development.bundle \
  --source-input sequence=source.fasta -o report.json
```

After successful argument parsing, semantic validation failures return a sealed
failure report with a stable error code and exit status 1. Command-usage errors
use the standard `argparse` exit status 2.

## Reproducibility and hardening

```bash
python -m unittest discover -s tests -v
python tests/verify.py source-causality
python tests/verify.py universal-translation
python tests/verify.py boundary
python tests/verify.py contract-attacks
python tests/verify.py resource-path-safety
python tests/verify.py development-real-dna
PIP_NO_INDEX=1 python tests/verify.py wheel
PIP_NO_INDEX=1 python tests/verify.py sdist
```

The gates exercise generated unseen sources, actual NCBI DNA, coherent
resealing attacks, import isolation, duplicate JSON, JSON type confusion,
gzip member attacks, decompression budgets, symlinks, hardlinks, path races,
unknown bundle children, cumulative blob limits, wheel installation, and sdist
self-tests.

The pinned large-record gate uses NCBI `CP032762.1` (5,868,661 bp, 10,919
features). The installed final-compilation gate uses NCBI `J02482.1` (5,386 bp)
and requires byte-identical repeated compilation plus successful independent
replay.

## Scientific basis

The compiler maps biological source standards onto ordinary compiler concepts:
grammar profiles, typed IRs, symbol/reference binding, static semantics,
lowering, linking, reproducible builds, and independent validation.

- [INSDC Feature Table Definition](https://www.insdc.org/submitting-standards/feature-table/)
- [NCBI GenBank release format](https://www.ncbi.nlm.nih.gov/genbank/release/current)
- [Sequence Ontology GFF3 1.26 syntax basis](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
- [GA4GH refget Sequences](https://ga4gh.github.io/refget/sequences/)
- [GA4GH refget Sequence Collections](https://ga4gh.github.io/refget/seqcols/)
- [RFC 8785 JSON Canonicalization Scheme](https://www.rfc-editor.org/rfc/rfc8785)
- [Cai et al., 2009 DNA attribute-grammar compiler](https://doi.org/10.1371/journal.pcbi.1000529)

[SPEC.md](SPEC.md) is normative. [docs/STANDARDS.md](docs/STANDARDS.md) traces
the standards mapping, and [docs/PRIOR_ART.md](docs/PRIOR_ART.md) places the
compiler among biological compiler systems.

Apache-2.0.

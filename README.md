# Epigenesis

Epigenesis 0.9.0 is a deterministic compiler for biological source languages.
Its GenBank 0.8 frontend (producer 0.8.1) compiles traditional physical-DNA
flat files into two linked, content-addressed representations:

- `brain01.sequence-collection-ir/v2`, which binds the exact source bytes to
  normalized IUPAC DNA and GA4GH refget identities; and
- `brainc.bio.insdc-genbank-ir/v2`, a structural BioIR for records, headers,
  features, qualifiers, location syntax, normalized coordinates, and source
  positions.

Package 0.9.0 adds an independently replayable 11-column INSDC feature-state
structural projection: every ordered BioIR feature becomes one development unit
and the adapter emits no inferred relationships. The content-bound GenBank
BioIR remains the semantic authority. A reference-only index binds nine
independently sealed child artifacts, including the development module; the
compiled GenBank artifact remains an external, digest-bound source dependency.

The frontend is pinned to [NCBI GenBank Flat File Release
273.0](https://www.ncbi.nlm.nih.gov/genbank/release/current) and [INSDC Feature
Table Definition 11.4](https://www.insdc.org/submitting-standards/feature-table/).
Independent validators rebuild the GenBank and feature-state outputs from their
bound inputs.

## Quick start

Python 3.11 or newer is required. Runtime dependencies: Python standard library.

```bash
python -m pip install .

brainc compile-genbank CP032762.1.gb -o CP032762.1.epigenesis.json
brainc validate-genbank CP032762.1.gb CP032762.1.epigenesis.json \
  --report CP032762.1.validation.json
brainc compile-insdc-graph CP032762.1.epigenesis.json \
  -o CP032762.1.feature-state
brainc-validate-insdc-graph CP032762.1.gb \
  CP032762.1.epigenesis.json CP032762.1.feature-state \
  -o CP032762.1.feature-state.validation.json
```

GenBank validation prints the sealed report to standard output and optionally
writes the identical report. It also has an isolated entry point:

```bash
brainc-validate-genbank CP032762.1.gb CP032762.1.epigenesis.json \
  -o CP032762.1.validation.json
```

A successful report has `valid: true` and binds the original source, compiled
artifact, embedded sequence collection, replayed BioIR, and report itself by
digest.

`compile-insdc-graph` publishes `bundle.json` plus nine sealed JSON children in
a new output directory and prints the bundle path as JSON. The isolated graph
validator accepts the original GenBank bytes, compiled GenBank artifact, and
that directory; its successful report is also sealed.

Python API:

```python
from brainc.insdc import GenBankCompiler
from brainc.insdc_graph import compile_insdc_graph
from brainc.validator_insdc import validate_genbank
from brainc.validator_insdc_graph import validate_insdc_graph_bundle

source = open("CP032762.1.gb", "rb").read()
artifact = GenBankCompiler().compile_bytes(source)
report = validate_genbank(artifact.to_dict(), genbank_source=source)
assert report["valid"]
state = compile_insdc_graph(artifact)
graph_report = validate_insdc_graph_bundle(
    state.to_dict(), state.artifacts, artifact.to_dict(), genbank_source=source
)
assert graph_report["valid"]
state.save("CP032762.1.feature-state")
```

## Why this is a compiler

An API carries requests and responses. A format converter rewrites syntax. A
compiler accepts a defined source language, rejects programs outside that
language, checks static meaning, normalizes the accepted program, and emits a
defined target representation with deterministic identity.

The GenBank frontend performs those compiler operations: fixed-field parsing,
record and symbol binding, feature-location parsing, coordinate checking,
source-coverage checks, bounded normalization, content-addressed emission, and
semantic replay. The command line and Python calls are interfaces to the same
compiler.

The frontend artifact gives downstream targets exact sequence, structured
annotations, normalized intervals, and provenance without reparsing a flat file
or trusting an accession lookup.

The release contains compiler code, independent replay validators, pinned
standards metadata, conformance fixtures, and versioned downstream contracts.

## Compiler pipeline

### GenBank frontend 0.8

1. Read bounded ASCII GenBank bytes and verify the pinned authority manifest.
2. Parse physical DNA records, headers, references, FEATURES, qualifiers,
   locations, ORIGIN sequence, and exact record order.
3. Normalize uppercase IUPAC DNA and locations from 1-based closed source
   coordinates to 0-based half-open segments.
4. Emit Sequence Collection v2 with one exact source store, one normalized
   sequence store per record, SHA-256, refget Sequence v2 identifiers, and a
   refget Sequence Collections v1.0.0 identity.
5. Emit structural BioIR whose source maps contain positions rather than
   duplicated source text.
6. Seal the collection, each BioIR record, the BioIR body, and the outer
   artifact with RFC 8785 canonical JSON and SHA-256.
7. Independently replay the original bytes and seal the resulting validation
   report with the standard-library validator.

The profile accepts ordered multi-record files, LF and CRLF source lines,
traditional DNA LOCUS records, accession/version identity, legacy GI values,
wrapped fields, references, physical ORIGIN sequence, and the implemented INSDC
location grammar: points, spans, endpoint fuzz, between-base and archived
within positions, exact-version remote references, `complement`, `join`, and
`order`. References to another record in the same input are bounds-checked;
external references remain explicit and unresolved.

Every dispatch decision comes from the wire profile and source grammar. An
accession, organism, gene, sequence digest, or test-fixture identity never
selects compiler behavior.

### INSDC feature-state adapter

1. Validate the closed GenBank artifact, source profile, authority pin, and
   whole-artifact and BioIR digests.
2. Create one unit for every feature in record-array and feature-ordinal order.
3. Encode eleven constant structural columns as little-endian `u64` tensors.
4. Compile one `unit.create` operation against a closed target contract and
   emit no inferred relationships.
5. Publish a reference-only `bundle.json` index and nine sealed child artifacts.

| Column | Exact value |
|---|---|
| `fuzzy_boundary_count` | Count of non-null `<`/`>` endpoint markers; archived within-position uncertainty is excluded |
| `key_code` | Bundle-local code into the UTF-8-byte-sorted feature-key dictionary |
| `kind_mask` | Presence mask: interval `1`, between-base `2`, uncertain-point `4` |
| `ordinal` | One-based feature ordinal within its source record |
| `orientation_mask` | Presence mask: forward `1`, reverse/complement `2` |
| `qualifier_count` | Number of ordered qualifier entries, including repeats |
| `record_code` | Bundle-local code into the UTF-8-byte-sorted record dictionary |
| `remote_segment_count` | Segments carrying a remote record reference |
| `segment_count` | Number of normalized segments |
| `segment_extent_sum` | Sum of `end - start` for every normalized 0-based half-open segment; between sites contribute `0`, uncertain points contribute their envelope width, and overlaps count with multiplicity |
| `unresolved_segment_count` | Segments whose remote reference bounds remain unresolved |

Dictionary codes are bundle-local, and the masks record kind and orientation
presence; exact multiplicity, order, and segment association remain in the
bound BioIR.

## Real-genome evidence

The GenBank 0.8 frontend was exercised on the pinned NCBI GenBank record
[`CP032762.1`](https://www.ncbi.nlm.nih.gov/nuccore/CP032762.1), then replayed by
the independent validator.

| Measurement | Result |
|---|---:|
| GenBank source | 13,418,250 bytes |
| Sequence | 5,868,661 bp |
| Features / normalized segments | 10,919 / 10,919 |
| Compiled artifact | 37,192,363 bytes |
| Sequence chunks | 6 × at most 1 MiB |
| Unresolved references | 0 |

Content identities:

| Object | Identity |
|---|---|
| Original GenBank bytes | `c54efd1ee2a811c8efbe12e3f8ea9757013877ce7e5faf8c06d39ac326a91582` |
| Normalized sequence | `b2c317b26275d822aae0063819b0b4012f247404edf87982d1eed303416bde62` |
| refget Sequence v2 | `SQ.SFU4YXxY_mSKaqHE06_SmfWkLboBMB-N` |
| Sequence Collection artifact | `84f05893b51d644511327d97d9417c1cf7f13583f3076115ca73e7bc042a61f6` |
| refget Sequence Collection v1.0.0 | `gF82ixP2u8Jn7b814kR1kORCsMs3nNa7` |
| BioIR | `5367765fa8310e0491ec39e8e496750a18b2a51b5f52daaab5ad681b7f929035` |
| Compiled artifact | `97b13016b899da15a76b66c5d2a94bf839f225e42a3f92fdc20bb63e873768b1` |
| Validation report | `447f3c4e32714ae9efc637657ed777f80a8faaceec2302c38c7ba4b9ecee5c0c` |

The same artifact completed the feature-state lowering and independent graph
replay:

| Result | Value |
|---|---|
| Units / inferred relationships / operations | `10,919 / 0 / 1` |
| Logical tensor payload | `960,880` bytes |
| Bundle index | `b825447b352e8b0302b2649c42a9e6de59b7c2d5bf911a4f10ad55bcd700adf4` |
| Development module | `21642fdebbcc50a24d103e650fe5707313ab48463dd1ad5be28abed9a9db4a94` |
| Independent replay | `valid: true`; report `85bb27f7a6b5ba0d99155451a7085df4854ba9aeed5ddd53db01395021d3d93a` |

The repository stores the pinned provenance sidecar, not the 13.4 MB source.
Run the scale gate with an independently obtained matching record:

```bash
python tests/verify.py genbank-scale --source /absolute/path/CP032762.1.gb
```

The pinned sidecar records the source provenance, content hashes, artifact size,
and structural counts used by this gate.

## Contracts

| Layer | Wire identity | Primary role |
|---|---|---|
| Exact DNA collection | `brain01.sequence-collection-ir/v2` | Raw-byte authority, normalized sequence chunks, record order, source maps, SHA-256 and refget identities |
| Structural annotation | `brainc.bio.insdc-genbank-ir/v2` | GenBank record structure, feature-table AST, normalized segments, qualifiers and positional provenance |
| Independent result | `brainc.bio.insdc-genbank-validation-report/v2` | Original-source replay and digest closure |
| Feature-state backend | `brainc.bio.insdc-graph-backend-spec/v1` | Normative eleven-column projection, dictionaries, masks, target schema, and operation |
| Feature-state semantics | `brainc.bio.insdc-graph-semantics/v1` | Source binding, unit order, bundle-local dictionaries, and empty relationship order |
| Feature-state index | `brainc.bio.insdc-graph-bundle/v1` | References to the GenBank semantic authority and nine sealed children |
| Development target | `brainc.development-module/v1` | One constant feature unit per ordered source feature |
| Feature-state replay | `brainc.bio.insdc-graph-validation-report/v1` | Independent reconstruction of the bound source, index, and nine-child adapter closure |

All objects are closed JSON contracts. Unknown or missing members, duplicate
JSON keys, non-I-JSON values, malformed digests, invalid coordinates, and
coherently resealed semantic modifications fail validation.

## Existing frontends and lowerings

Package 0.9.0 retains the earlier versioned paths without reinterpreting them:

| Path | Commands | Contracts |
|---|---|---|
| Sequence | `compile-sequence`, `compile-collection`, `compile-raw` | Sequence IR v2 and Sequence Collection v1 |
| Scalar provider | `make-request`, `compile`, `check`, `validate` | 0.4 state-program chain |
| Tensor / target | `make-request-v2`, `compile-v2`, `check-v2`, `validate-v2` | 0.5 development-module chain |
| External GFF3 | `compile-gff3`, `validate-gff3` | GFF3 BioIR v1 bound to Sequence Collection v1 |
| Feature graph | `compile-feature-graph`, `validate-feature-graph` | 0.6 structural graph bundle |
| INSDC feature state | `compile-insdc-graph`, `brainc-validate-insdc-graph` | GenBank BioIR v2 to the eleven-column development bundle |

Sequence Collection v2 is the GenBank 0.8 sequence authority. The current
`compile-insdc-graph` command is its explicit feature-state lowering. The GFF3
and retained tensor paths continue to consume their declared Sequence
Collection v1 contracts.

## Verification

```bash
python -m unittest discover -s tests -v
python tests/verify.py real-dna
python tests/verify.py tamper
python tests/verify.py canonical
python tests/verify.py contract-attacks
python tests/verify.py boundary
python -m unittest tests.test_insdc_graph tests.test_cli_insdc_graph -v
python -m unittest tests.test_validator_insdc_graph \
  tests.test_validator_insdc_graph_paths \
  tests.test_validator_import_isolation -v
PIP_NO_INDEX=1 python tests/verify.py wheel
PIP_NO_INDEX=1 python tests/verify.py sdist
```

The release gates cover official GenBank source, multi-record compilation,
canonical chunks, exact source replay, rehashed tampering, resource ceilings,
fixture-independent production code, installed-wheel isolation, CLI behavior,
packaging of the pinned standards manifest, and exact feature-state lowering.

## Standards and extension path

The compiler grows through separately versioned source-language frontends, not
through accession-specific cases. Planned frontends include VCF 4.5 and GA4GH
VRS 2, SBOL 3.1, SBML Level 3, CellML 2.0, NeuroML 2, and SONATA. Each frontend
must bring a pinned authority, closed normalized IR, explicit links, finite
admission, target-bound lowering, and independent replay.

[SPEC.md](SPEC.md) defines the normative package 0.9.0 and GenBank frontend 0.8
contracts.
[docs/STANDARDS.md](docs/STANDARDS.md) traces the primary standards and prior
biological compilers.

Licensed under Apache-2.0.

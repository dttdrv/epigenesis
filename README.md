# Epigenesis

Epigenesis 1.5 is a universal DNA translation compiler targeting Development
Module IR, with experimental consumers for neural mechanisms, developmental
construction, coding consequences and DNA-binding preferences.

The [1.5 release notes](docs/RELEASE-1.5.md) describe the supported execution
paths and evidence. A validated natural-genome-to-brain model remains an open
research objective. The measured reporter and binding predictors failed their
scientific acceptance criteria; those results are preserved.

It accepts built-in and external DNA source languages, preserves exact source
identity and profile-native structure, normalizes sequence identity into a
common typed source boundary, validates a caller-supplied interpretation, and lowers it
against an explicit target contract.

```text
DNA source language
    -> grammar frontend
    -> typed source bundle + profile-native IR
    -> validated caller-supplied interpretation
    -> target-bound lowering
    -> Development Module IR
    -> sealed compilation bundle
```

Raw IUPAC, FASTA, streaming reference FASTA, GenBank, and GFF3 are built in.
Any other DNA-bearing grammar can enter through a version-2 executable external
frontend and a separate executable validator. Both programs are explicitly
selected and digest-pinned; the final Development Module validator reruns the
external validator against the original source.

The compiler and its producer-isolated validators use only Python 3.11+ and the
standard library.

## Complete source-to-module example

The repository includes every input for a reproducible compilation:

```bash
python3.11 -m pip install .

demo_dir=$(mktemp -d)
brainc compile-source \
  --profile fasta-reference-dna/v1 \
  --sequence examples/minimal-development/sequence.fasta \
  --wrapper identity \
  --output "$demo_dir/source"

brainc compile-development "$demo_dir/source" \
  --manifest examples/minimal-development/interpretation/manifest.json \
  --request examples/minimal-development/interpretation/request.json \
  --response examples/minimal-development/interpretation/response.json \
  --policy examples/minimal-development/interpretation/lowering-policy.json \
  --target examples/minimal-development/interpretation/target-contract.json \
  --output "$demo_dir/development"

brainc-validate-development "$demo_dir/source" "$demo_dir/development" \
  --source-input sequence=examples/minimal-development/sequence.fasta \
  --output "$demo_dir/validation.json"
```

The report has `valid: true`. Stable identities are recorded in
[`examples/minimal-development/expected.json`](examples/minimal-development/expected.json):

| Artifact | SHA-256 identity |
|---|---|
| Source descriptor | `7199f7b99358c3c2fa8b8df8f44a5431dac722aefc1d5ff18de9fc52d4ac4ced` |
| Development Module | `0d835308d23fe83a5f4cd0e28c26c678c510049a381752d8671c39e3805d2cea` |
| Development bundle | `9cc0b7c95027d777212beb75cdea04666c3333368ca5769c3203c6547d74b313` |
| Validation report | `3b8dafa503bcb376b50c798c0e4ceb2cf254505adeff83253d854fca61dc88d3` |

## Neural-fate mechanism experiment

The [neural-fate example](examples/neural-fate/README.md) compiles an artificial
512-base parameter recipe, validates its Development Module, and computes
Pax6/Olig2/Nkx2.2 dynamics from two published models:

```bash
demo_dir=$(mktemp -d)
python3.11 examples/neural-fate/experiment.py --output "$demo_dir/neural-fate"
```

It tests signal-history dependence, gene loss, and a steady-to-oscillating
parameter contrast. The output includes complete trajectories, sealed modules,
validation reports, and numerical refinement checks. Its caller-owned decoder
and target consumer are outside the compiler package. The recipe stores model
parameters; it does not infer a mechanism from natural DNA or construct a brain.
The [acceptance protocol](docs/NEURAL_MECHANISM_ACCEPTANCE.md) distinguishes
these executable results from biological validation.

## Measured neural regulatory sequences

The [neural-regulation example](examples/neural-regulation/README.md) fits a
sequence-only predictor to real enhancer perturbations measured during human
neural induction. It includes pinned source data, family-isolated evaluation,
competing baselines, and compilation of predicted reporter effects from two
exact DNA sequences. It is an external interpreter; the compiler retains its
caller-supplied interpretation boundary.

The first frozen evaluation reduced targeted-variant MSE by 20.97% versus
predicting no effect, but **failed its predictive acceptance criterion** because
superiority over a dinucleotide model was uncertain. It does not establish
endogenous expression, developmental dynamics, or brain construction.

## Local neural construction

The [construction example](examples/neural-construction/README.md) consumes a
validated artificial recipe, starts from one founder and grows cells, neurites
and physical contacts using the pinned CX3D source. An exact integrate-and-fire
assay tests transmission through the final contacts and its loss after lesions.
Recorded development and activity can be explored in an offline viewer.

```sh
python3.11 examples/neural-construction/construction.py --output /tmp/my-construction
```

It requires a JDK and `patch`. This demonstrates local construction and synthetic
transmission controls; natural DNA mapping, biological brain equivalence and
intelligence remain unvalidated. The compiler package still owns no simulator.

## Natural DNA coding consequences

The [coding example](examples/coding-consequences/README.md) reconstructs
published local Sox2 mutations on genuine genomic DNA, compiles the actual
coding nucleotides, validates their source binding, and translates the loaded
module with the NCBI standard genetic code. Tests distinguish real reference
haplotypes that encode identical normal proteins but different mutant peptides.
This predicts primary amino-acid sequence under an explicit coding-region
contract. Protein activity and its connection to neural development remain
unvalidated.

The [binding example](examples/neural-binding/README.md) assembles annotated
Nkx2-2 coding exons, checks its translated domain and computes relative DNA
preferences using the published FamilyCode single-residue model. Source and
arithmetic checks pass. Its frozen comparison uses all 92 published mutants
and 316 fitted matrices; the 66 clone-supported primary mutants fail all three
predictive criteria. A connection to neural-development parameters remains
unvalidated.

## Source languages

Source selection is explicit; Epigenesis never guesses a grammar from a
filename, accession, organism, or fixture.

| Profile | Input | Typed result |
|---|---|---|
| `raw-iupac-dna/v1` | One raw IUPAC DNA sequence | Exact Sequence Collection IR |
| `fasta-dna/v1` | FASTA or concatenated gzip FASTA | Exact Sequence Collection IR |
| `fasta-reference-dna/v1` | Streaming chromosome/genome FASTA | Compact refget sequence catalog |
| `genbank-273-traditional-dna-physical-structural/v2` | Physical-DNA GenBank records | Sequence authority plus structural INSDC BioIR |
| `gff3-external-sequence/v1` | GFF3 plus an exact sequence source | Collection-bound structural GFF3 BioIR |
| `external-dna-source/v1` | Any manifest-declared DNA grammar | Manifest-qualified native IR plus common sequence catalog |

Every route produces ordered records with record id, base count, SHA-256, and
GA4GH refget identity. Profile-specific information stays in the native IR
rather than being flattened away.

The reference FASTA frontend streams very large inputs without embedding
sequence strings:

```bash
brainc compile-source \
  --profile fasta-reference-dna/v1 \
  --sequence assembly.fasta.gz \
  --wrapper gzip \
  --maximum-input-bytes 137438953472 \
  --maximum-logical-bytes 274877906944 \
  --output assembly.source
```

## Add a DNA grammar

The executable external ABI makes the compiler format-extensible without a
core patch. A profile declares grammar authority, source roles, native JSON
schemas, record and byte budgets, and two command identities:

- the frontend parses stable snapshots and returns native typed JSON plus the
  common sequence-record catalog;
- the validator separately reparses those snapshots, validates native JSON,
  and reproduces the exact catalog and references.

The included FASTQ profile is a complete implementation:

```bash
python3.11 examples/external-fastq/profile.py \
  > /tmp/epigenesis-fastq-profile.json

brainc compile-external-source \
  --profile-manifest /tmp/epigenesis-fastq-profile.json \
  --frontend-executable examples/external-fastq/frontend.py \
  --validator-executable examples/external-fastq/validator.py \
  --source-input reads=examples/external-fastq/reads.fastq \
  --output /tmp/epigenesis-fastq-source
```

The compiler verifies both executable SHA-256 values, snapshots exact input
bytes, invokes no shell, bounds runtime and output, derives native references,
requires an exact validator replay, and atomically publishes the source bundle.
The bundle embeds its profile-native JSON and complete validation closure.

After compiling a Development Module from that source, validate it with the
same explicitly selected validator:

```bash
brainc-validate-development \
  /tmp/epigenesis-fastq-source \
  /tmp/epigenesis-fastq-development \
  --source-input reads=examples/external-fastq/reads.fastq \
  --external-validator examples/external-fastq/validator.py \
  --output /tmp/epigenesis-fastq-validation.json
```

See [`examples/external-fastq`](examples/external-fastq) for the protocol in a
small pair of independent programs.

Version-1 external manifests and closures remain readable for compatibility.
That legacy data-only path verifies caller-attested files and hashes. New
profiles should use the version-2 executable path, which performs the grammar
and semantic replay during compilation and final validation.

## Development Module lowering

The interpretation stage is explicit and content-bound. A provider response
supplies typed tensors for the exact source descriptor; a lowering policy maps
those tensors to operations authorized by the target ABI:

- `unit.create`
- `edge.create`
- `rule.attach`
- `port.bind`

Epigenesis validates tensor types, shapes, axes, storage, coordinates, source
binding, target schemas, operations, references, and budgets before emitting
`brainc.development-module/v1`.

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

Publication is atomic, no-clobber, and absent-to-complete on supported Linux
and macOS hosts.

## Contract identities

| Contract | Identity |
|---|---|
| Source descriptor | `brainc.source-descriptor/v1` |
| Reference catalog | `brainc.reference-sequence-catalog/v1` |
| External frontend manifest | `brainc.external-frontend-profile/v1` and `/v2` |
| External source closure | `brainc.external-source-closure/v1` and `/v2` |
| Target contract | `brainc.target-contract/v1` |
| Development Module | `brainc.development-module/v1` |
| Compilation bundle | `brainc.development-bundle/v1` |
| Validation result | `brainc.development-validation-report/v1` |

All JSON contracts are closed, duplicate-key rejecting I-JSON. Artifacts are
sealed with RFC 8785 canonicalization and SHA-256. Original bytes, native IR,
profile manifest, interpretation, lowering policy, target contract, module,
and validation result form one transitive content-addressed closure.

External commands are trusted code chosen by the caller. Digest pinning proves
which bytes ran; process limits are not a security sandbox.

## Validation and release gates

```bash
python3.11 -m unittest discover -s tests -v
python3.11 tests/verify.py source-causality
python3.11 tests/verify.py universal-translation
python3.11 tests/verify.py minimal-example
python3.11 tests/verify.py neural-mechanism-controls
python3.11 tests/verify.py neural-fate-experiment
python3.11 tests/verify.py neural-regulation
python3.11 tests/verify.py neural-construction
python3.11 tests/verify.py coding-consequences
python3.11 tests/verify.py neural-binding
python3.11 tests/verify.py boundary
python3.11 tests/verify.py contract-attacks
python3.11 tests/verify.py resource-path-safety
PIP_NO_INDEX=1 python3.11 tests/verify.py development-real-dna
PIP_NO_INDEX=1 python3.11 tests/verify.py wheel
PIP_NO_INDEX=1 python3.11 tests/verify.py sdist
```

The universal translation gate compiles all built-in profiles and the
executable FASTQ profile to the same typed record identity, lowers every route
to Development Module IR, and reruns the external validator. Security gates
cover coherent resealing, changed executables and sources, lying frontends,
malformed grammar input, duplicate JSON, type confusion, gzip bombs, resource
ceilings, symlinks, hardlinks, path races, unknown bundle children, and
cumulative blob limits.

The pinned large-record gate uses NCBI `CP032762.1`; the installed final-chain
gate uses NCBI `J02482.1` and requires byte-identical repeated compilation plus
producer-isolated replay.

## Standards

- [INSDC Feature Table Definition](https://www.insdc.org/submitting-standards/feature-table/)
- [NCBI GenBank release format](https://www.ncbi.nlm.nih.gov/genbank/release/current)
- [Sequence Ontology GFF3 1.26](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
- [GA4GH refget Sequences](https://ga4gh.github.io/refget/sequences/)
- [GA4GH refget Sequence Collections](https://ga4gh.github.io/refget/seqcols/)
- [RFC 8785 JSON Canonicalization Scheme](https://www.rfc-editor.org/rfc/rfc8785)

[`SPEC.md`](SPEC.md) is normative. [`docs/STANDARDS.md`](docs/STANDARDS.md)
traces standards mappings, and [`docs/PRIOR_ART.md`](docs/PRIOR_ART.md) places
the compiler among biological compiler systems.

Apache-2.0.

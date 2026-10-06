# Epigenesis

Epigenesis is a universal DNA translation compiler with reproducible
neural-mechanism experiments. Compare published measurements, explore regulatory
models, and trace each result back to its sequences, assumptions and source data.

Version 1.6 adds two offline neuroscience workflows:

- [Neural evidence](examples/neural-evidence/README.md): inspect a published
  CHE-1 depletion experiment, from individual measurements to intended DNA edits.
- [Neural patterning](examples/neural-patterning/README.md): explore how an
  enhancer substitution and signal history affect a four-gene developmental model.

[Release notes](docs/RELEASE-1.6.md) · [Specification](SPEC.md) ·
[Research roadmap](docs/NEURAL_PREDICTOR_ROADMAP.md)

## Start with a real experiment

From a repository checkout or extracted source archive, run with Python 3.11+.
Both demos use the standard library and bundled data; no installation, GPU or
network connection is needed. Choose an output directory that does not yet exist.

```sh
python3.11 examples/neural-evidence/evidence.py demo --output /tmp/neural-evidence
python3.11 examples/neural-evidence/evidence.py verify --output /tmp/neural-evidence
```

Open `/tmp/neural-evidence/index.html` in your browser.

The report reconstructs 109 supplied observations from
[Traets et al., Figure 6F](https://doi.org/10.7554/eLife.66955).
Select a measurement to inspect its original workbook cell, compare constructs,
inspect donor sequences, and export SVG, CSV or JSON. Print the report to PDF.

The practical question: does transferring a 12-base binding core preserve the
response obtained with a larger DNA fragment? The observed auxin/control mean
ratios are 8.70% for core transfer and 66.33% for fragment transfer. The viewer
shows how sequence extent, orientation and replacement geometry differ, helping
you design a comparison that separates these explanations. Sampling limits and
conflicting source labels remain attached to the observations.

## Explore a developmental mechanism

```sh
python3.11 examples/neural-patterning/patterning.py demo --output /tmp/neural-patterning
```

Open `/tmp/neural-patterning/index.html`. Move through developmental time,
compare an endogenous edit with a passive reporter, and examine signal history,
allele rescue, mediator rescue and phenocopy.

The experiment connects a real Nkx2.2 enhancer substitution, measured GLI binding
profiles and published Pax6/Olig2/Nkx2.2/Irx3 equations. The transfer from an
in-vitro binding-profile ratio to intracellular affinity is an explicit model
assumption. Simulated concentrations and measured embryo reporter counts are
presented separately. Numerical refinement and intervention checks accompany
every generated report.

See the [experiment documentation](examples/neural-patterning/README.md) for
source papers, units, assumptions and custom FASTA input.

## The compiler

Epigenesis preserves exact DNA source identity, validates a caller-supplied
interpretation, and lowers it against a target contract:

```text
DNA source → typed source bundle → explicit interpretation → Development Module
```

The interpretation supplies the biological or computational meaning. The
compiler checks types, shapes, coordinates, references and budgets, then seals
the compilation with content hashes. A separate validator replays the chain.

Raw IUPAC, FASTA, streaming reference FASTA, GenBank and GFF3 are built in.
Other grammars use a version-2 executable external frontend and a separate
executable validator. Both are digest-pinned and explicitly selected by the
caller. The [FASTQ example](examples/external-fastq/README.md) implements this
extension without changing the compiler.

## Complete source-to-module example

Install the compiler from the checkout into a virtual environment:

```sh
python3.11 -m venv .venv
. .venv/bin/activate
python -m pip install .

demo_dir=$(mktemp -d)
brainc compile-source \
  --profile fasta-reference-dna/v1 \
  --sequence examples/minimal-development/sequence.fasta \
  --wrapper identity --output "$demo_dir/source"

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

The validation report has `valid: true`. This small example compiles an explicit
base-count interpretation. The wheel installs the compiler and validators;
the checkout and source archive also contain the examples, data and tests.

<details>
<summary>Stable example identities</summary>

These identities are recorded in [expected.json](examples/minimal-development/expected.json)
and checked by the minimal-example release gate.

| Artifact | SHA-256 |
|---|---|
| Source descriptor | `7199f7b99358c3c2fa8b8df8f44a5431dac722aefc1d5ff18de9fc52d4ac4ced` |
| Development Module | `0d835308d23fe83a5f4cd0e28c26c678c510049a381752d8671c39e3805d2cea` |
| Development bundle | `9cc0b7c95027d777212beb75cdea04666c3333368ca5769c3203c6547d74b313` |
| Validation report | `3b8dafa503bcb376b50c798c0e4ceb2cf254505adeff83253d854fca61dc88d3` |

</details>

## Research objective

The long-term goal is to predict nervous-system development from a genome and
starting conditions, then connect the resulting neural structure to activity
and learning. The intended use is to compare genetic or environmental
interventions and identify experiments that distinguish competing mechanisms.

The current examples establish specific parts of that workflow. Useful
DNA-edit prediction is the next separate validation milestone: the candidate
must improve both effect prediction and experimental ranking against fixed
baselines on unseen sequence families. Current candidates have not passed both
requirements; reserved confirmation outcomes remain unopened. The
[roadmap](docs/NEURAL_PREDICTOR_ROADMAP.md) records the protocols, results and next
decisions, including negative results.

## More experiments

| Example | What you can examine |
|---|---|
| [Neural fate](examples/neural-fate/README.md) | An artificial DNA parameter recipe driving published regulatory dynamics |
| [Neural construction](examples/neural-construction/README.md) | A local development recipe producing cells, neurites, contacts and activity |
| [Coding consequences](examples/coding-consequences/README.md) | Natural Sox2 coding edits and independently checked translation consequences |
| [Neural regulation](examples/neural-regulation/README.md) | H1 neural-induction reporter data, evaluation and a failed predictive test |
| [Neural binding](examples/neural-binding/README.md) | Genomic Nkx2-2 sequence, a designed substitution and a failed binding benchmark |

Each example documents its inputs, runtime, scientific scope and reproduction
commands. Runtime requirements vary; neural construction also needs a JDK and `patch` (tested with Java 21).

## Verification

```sh
python3.11 -m unittest discover -s tests -v
python3.11 tests/verify.py minimal-example
python3.11 tests/verify.py universal-translation
python3.11 tests/verify.py boundary
PIP_NO_INDEX=1 python3.11 tests/verify.py wheel
PIP_NO_INDEX=1 python3.11 tests/verify.py sdist
```

The [CI workflow](.github/workflows/ci.yml) defines the complete release checks,
including independent source replay, contract attacks, resource and path safety,
real-DNA compilation, and both new demos. Package checks require the declared
build tooling; the offline commands use tools already installed on the host.

[SPEC.md](SPEC.md) defines the contracts. [Standards](docs/STANDARDS.md) maps
formats and identities to their specifications.
[Prior art](docs/PRIOR_ART.md) and the
[acceptance protocol](docs/NEURAL_MECHANISM_ACCEPTANCE.md) describe the research
context and criteria.

Compiler code: [Apache-2.0](LICENSE). Bundled scientific data and third-party
runtimes retain their documented licenses and attribution.

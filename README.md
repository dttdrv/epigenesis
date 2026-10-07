# Epigenesis

Epigenesis is a research platform for genome-driven neural development. Its goal
is to turn a genome and starting conditions into a system that develops neural
structure, produces activity, and learns.

```text
genome + starting conditions → development → neural structure → activity and learning
```

The central question is causal: how do genetic instructions and developmental
conditions produce a nervous system, and what changes when we intervene?
Epigenesis makes proposed mechanisms executable, connects them to biological
evidence, and preserves the source-to-result chain needed to examine them.

## Working components

Version 1.6 combines a universal DNA translation compiler with executable
experiments along this path:

| Component | What runs today |
|---|---|
| [Natural DNA to coding consequences](examples/coding-consequences/README.md) | Reconstruct Sox2 edits from genomic sequence, execute the declared coding annotation, and check the resulting peptides against independent reference sequences. |
| [Regulatory sequence to developmental dynamics](examples/neural-patterning/README.md) | Connect an Nkx2.2 enhancer substitution and measured GLI binding profiles to a four-gene model; test signal history, rescue and phenocopy under an explicit affinity-transfer assumption. |
| [Development to structure and activity](examples/neural-construction/README.md) | Grow cells, neurites and contacts from one founder and a shared artificial DNA recipe. Test signal transmission and contact lesions on the developed network. The final morphology and connections emerge during execution. |
| [Experimental evidence](examples/neural-evidence/README.md) | Reconstruct 109 published CHE-1 depletion observations, trace them to original workbook cells and intended DNA edits, and examine competing explanations for the measured response. |
| [Compilation and independent replay](SPEC.md) | Parse exact DNA sources, validate an explicit interpretation, compile a sealed Development Module, and replay its derivation with a separate validator. |

The research program is to connect these components through mechanisms that
survive independent biological tests. Each current experiment records its own
inputs, assumptions, interventions and verification results.

The two new 1.6 workflows, neural evidence and neural patterning, run offline
with Python's standard library. They are the quickest entry points below.

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

## Next proof: choose useful DNA edits

A researcher has candidate regulatory DNA edits and a limited experimental
budget. The next practical milestone is to predict each edit's direction and
size of effect, then select ten edits from distinct sequence families for
follow-up measurement. The current target is reporter activity in
WTC11-NGN2 induced neurons.

Success requires better effect prediction and better experimental ranking than
fixed sequence baselines on held-out families, with uncertainty reported for
both comparisons. This would establish a useful sequence-to-function link for
the larger development program and a concrete tool for planning experiments.

A compact [experimental predictor](examples/neural-predictor/README.md) now
ships with frozen weights and a CPU/Metal inference command. It combines a
local sequence baseline with two 11,368-parameter attention models using
[Interlace's rational gate](https://misul.org/interlace). You can score your own
270-base reference/edited pairs and inspect each component's contribution.

Current Salomon validation results use 7,390 edits in 3,029 families:

| Predictor or comparator | Effect-size MSE ↓ | Signed top-ten utility ↑ |
|---|---:|---:|
| Always predict zero | 0.02875171 | 0 (abstains) |
| Random selection, TRAIN-selected direction | n/a | 0.01674875 |
| k5 ridge, protocol's MSE-selected comparator | 0.02862511 | 0.03510890 |
| k3 ridge, highest utility among 16 baselines | 0.02869761 | 0.21515632 |
| Calibrated short model | 0.02863101 | 0.07609249 |
| **Bundled two-cycle mixture** | **0.02862935** | **0.07158510** |
| Subsequent precision-sampling trial | 0.02862093 | 0.00342057 |

The mixture reduces MSE by **0.426% relative to predicting zero**. Its MSE is
**0.01480% above k5**, calculated as `100 × (MSE_model / MSE_k5 − 1)`;
this percentage measures the gap to that comparator, not predictive accuracy.
MSE itself is in squared log2 reporter-effect units. Centered weighted R² is
−0.00555, so broad effect-size prediction remains weak.

Ranking utility is 2.04× the protocol's k5 comparator and 4.27× random selection,
while the k3 baseline ranks better. These are adaptive-validation point estimates
for ten selected families. Both original requirements remain fixed: MSE below
0.028625109823 and utility above 0.035108902127. Confirmation is unopened.

Independent inference reproduces **all 29,560 saved Metal values exactly**:
the baseline, both corrections and their mixture on every validation input.
This establishes numerical reproducibility; agreement with measured biology is
reported separately in the table. [Results and definitions](examples/neural-predictor/README.md#validation-results)
include the complete baseline panel and the underlying aggregate receipts.

The subsequent precision-sampling trial crosses the MSE limit but fails ranking.
Independent replay matches all 109,300 saved prediction values. It remains a
failed joint-acceptance result and is not the bundled predictor.

Interlace also improved distant-interaction and local-effect recovery in matched
synthetic tests, with 51.1% and 64.4% lower error than ordinary attention at equal
parameter count in one paired training seed. The biological experiment addresses the question of
whether that capacity transfers to measured neuronal reporter effects.

A completed audit recovers model-based measurement uncertainty for all 21,781
training edits. Uncalibrated neural corrections worsen error even in the most
precisely measured group. The completed sampling trial uses that uncertainty
without changing targets or acceptance limits. The research ledger records
287 biological fits, including failed experiments.
[Results and reproduction details](examples/neural-predictor/README.md)
include model identities, independent replay receipts and the uncertainty audit;
the [roadmap](docs/NEURAL_PREDICTOR_ROADMAP.md) preserves the full development history.

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

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

Today it helps a researcher answer a concrete question: which DNA intervention
would distinguish competing explanations of a neural regulatory response?
The working examples connect original measurements, reconstructed DNA edits,
binding models and developmental dynamics. A separate measurement pilot improves
omitted-repeat forecasts on four neuronal-assay components.

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
[Research roadmap](docs/NEURAL_PREDICTOR_ROADMAP.md) ·
[Neuroscience presentation](docs/NEUROSCIENCE_PRESENTATION.md)

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

## Measurement result

Before learning DNA effects, we need reliable measurements of those effects.
A [BIRD-based forecast](https://doi.org/10.1093/bioinformatics/btz545) predicts an
omitted RNA repeat from the other repeats in a fixed neuronal-assay pilot:

| Twelve omissions across four TRAIN components | Mean CRPS loss, lower is better |
|---|---:|
| Pooled Beta-binomial | 0.06219211 |
| BIRD repeat model | 0.05490724 |
| Relative loss reduction | 11.7% |

The gain's conditional 99% simulation interval is [0.00436489, 0.00975665].
All four component means improve; one individual omission favors pooling and
remains included. The interval describes numerical integration uncertainty for
this selected cohort. Broader biological generalization remains an open test.
BIRD and pooling also differ in priors and DNA use, so this comparison does not
isolate repeat variation as the cause of the gain.

[Replay the archived calculation](examples/neural-measurement/README.md) without
training, new random draws or access to private local paths. The package checks
saved forecast moments and score arithmetic; the original sampling run is
recorded separately. CRPS is a distributional forecast loss on RNA allele
fractions, not a percentage error in predicting a DNA edit.

Follow-up results remain visible in the
[research record](docs/NEURAL_PREDICTOR_ROADMAP.md#dna-normalized-effect-prediction).
A fresh DNA-normalized forecast reduces mean CRPS by 15.8% relative to pooling
and passes both fixed directional comparisons. A post hoc zero-effect forecast
scores better, keeping practical effect prediction unresolved. The separate
treatment-response comparison has a 14.22% point reduction but inconclusive
directional intervals. These are different targets and comparisons.

A new [temporal experiment](docs/NEURAL_PREDICTOR_ROADMAP.md#new-experiments-18-october-2026)
forecasts later enhancer responses from four early measurements. Across 351
excluded families, full-history ridge reduces error by 24% versus persistence
but matches a simpler averaging control. Timestral-5M and 209M also complete
the same test; neither improves on that control.

A subsequent [raw-count experiment](examples/neural-forecast/README.md) holds out
both sequence families and infection cultures. Its two-feature forecast has MSE
0.32132339 versus 0.32214740 for zero change, a 0.26% point improvement. All four
gain intervals include zero, and the early-mean control wins in one culture.
The completed test does not establish transferable forecasting superiority.
Its portable replay checks all 15 evaluation groups and preserves 164 forecasts
whose held-culture outcomes lack complete measurements. This is a different
observable and target from the preceding pooled-activity test.

## Experimental DNA-edit predictor

The next sequence-model milestone is to select ten useful regulatory edits from
distinct families for follow-up experiments. Success requires both better
prediction error and better selection utility than fixed sequence baselines.

The [bundled predictor](examples/neural-predictor/README.md) scores 270-base
reference/edited pairs in WTC11-NGN2 induced neurons. It combines a local sequence
baseline with two 11,368-parameter attention models using
[Interlace's rational gate](https://misul.org/interlace). Frozen weights and
CPU/Metal inference are included.

On 7,390 validation edits in 3,029 families, the mixture has MSE 0.02862935
versus the k5 baseline's 0.02862511. Its signed top-ten utility is 0.07158510
versus k5's 0.03510890; the k3 baseline ranks better at 0.21515632. It has not
passed the joint acceptance rule. MSE is in squared log2 reporter-effect units.

All 29,560 saved Metal prediction values reproduce exactly. This checks the
implementation; biological performance is measured by the comparisons above.
The [complete baseline panel and receipts](examples/neural-predictor/README.md#validation-results)
include the later failed trials and the corrected quadratic-readout tie.
The current research ledger records 445 fits, with confirmation unopened and the
bundled weights unchanged. The new temporal comparisons are complete and their
measured outcomes are retained. A useful general DNA-edit predictor remains an
open research problem.

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

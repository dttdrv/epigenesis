# Epigenesis

Epigenesis is a DNA compiler and simulation toolkit. It compiles a DNA source, read under an explicit interpretation, into a sealed Development Module that executable models consume, and a separate validator replays every compilation. Each simulated result can be traced to its source bytes and reproduced.

Experiments built on the compiler connect original measurements, reconstructed DNA edits, binding models and developmental dynamics. Every experiment records its inputs, assumptions, interventions and verification results, so a claim can be traced back to the source data. In the works reviewed in [Prior art](docs/reference/PRIOR_ART.md), no other system combines these boundaries; that is a scoped comparison of the cited corpus, not a benchmark against other simulators.

**Status.** Version 1.6 combines a universal DNA translation compiler with a set of executable experiments, all independently replayable. Its experimental edit predictor produces signed reporter-effect scores for the SK-N-SH assay (a preliminary 3.02% error reduction against a 5-mer baseline within that assay). Accuracy across studies is not established: thirteen frozen tests of sequence-based variant-effect prediction, run on 19 and 20 October 2026, found no reliable gain over simple baselines (mean correlation 0.049 on five external studies). [The technical paper](docs/paper/epigenesis-technical-paper.pdf) ([Bulgarian](docs/paper/epigenesis-technical-paper.bg.pdf)) states what is and is not established.

## Contents

| Path | What it holds |
|---|---|
| [`docs/paper/`](docs/paper/) | Technical paper as PDF (English, Bulgarian), with LaTeX sources and the build script |
| [`docs/`](docs/README.md) | Index of references, release contracts and research records |
| [`research/`](research/README.md) | Plans, hashes, scripts and results of the thirteen frozen tests, with a hash checker |
| [`examples/`](examples/) | One folder per experiment, each with its own README |
| [`studio/`](studio/README.md) | Epigenesis Studio: a desktop application (English and Bulgarian) to change a DNA site and simulate neural-development dynamics, scan and rank edits with the delivered predictor, and verify sealed results |
| `brainc/`, [`SPEC.md`](SPEC.md) | The DNA → Development Module compiler and its contract |
| `tests/` | 505 unit tests and the release gates in `tests/verify.py` |

## Quick start

Python 3.11 or later. The two offline demonstrations use only the standard library and bundled data; they need no installation, GPU or network. Choose an output directory that does not yet exist.

```sh
python3.11 examples/neural-evidence/evidence.py demo --output /tmp/neural-evidence
python3.11 examples/neural-evidence/evidence.py verify --output /tmp/neural-evidence
python3.11 examples/neural-patterning/patterning.py demo --output /tmp/neural-patterning
```

Open `index.html` in each output directory.

## The compiler

Epigenesis preserves exact DNA source identity, validates a caller-supplied interpretation and lowers it against a target contract:

```text
DNA source → typed source bundle → explicit interpretation → Development Module
```

The interpretation supplies the biological or computational meaning. The compiler checks types, shapes, coordinates, references and budgets, then seals the compilation with content hashes (SHA-256 over RFC 8785 canonical JSON). A separate validator replays the chain without importing the producer.

Raw IUPAC, FASTA, streaming reference FASTA, GenBank and GFF3 are built in. Other grammars use a version-2 executable external frontend and a separate executable validator. Both are digest-pinned and explicitly selected by the caller. The [FASTQ example](examples/external-fastq/README.md) implements this extension without changing the compiler. [SPEC.md](SPEC.md) defines the contracts; [Standards](docs/reference/STANDARDS.md) maps formats and identities to their specifications; [Prior art](docs/reference/PRIOR_ART.md) places the design among existing systems.

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

The validation report has `valid: true`. This small example compiles an explicit base-count interpretation. The wheel installs the compiler and validators; the checkout and source archive also contain the examples, data and tests.

<details>
<summary>Stable example identities</summary>

These identities are recorded in [expected.json](examples/minimal-development/expected.json) and checked by the minimal-example release gate.

| Artifact | SHA-256 |
|---|---|
| Source descriptor | `7199f7b99358c3c2fa8b8df8f44a5431dac722aefc1d5ff18de9fc52d4ac4ced` |
| Development Module | `0d835308d23fe83a5f4cd0e28c26c678c510049a381752d8671c39e3805d2cea` |
| Development bundle | `9cc0b7c95027d777212beb75cdea04666c3333368ca5769c3203c6547d74b313` |
| Validation report | `3b8dafa503bcb376b50c798c0e4ceb2cf254505adeff83253d854fca61dc88d3` |

</details>

## Executable experiments

| Experiment | What runs | Main source |
|---|---|---|
| [Coding consequences](examples/coding-consequences/README.md) | Reconstructs *Sox2* edits from genomic sequence, applies the declared coding annotation and checks the resulting peptides against independent reference sequences | Gong et al. (2020) |
| [Neural patterning](examples/neural-patterning/README.md) | Connects an *Nkx2.2* enhancer substitution and measured GLI binding profiles to a four-gene model; tests signal history, rescue and phenocopy under an explicit affinity-transfer assumption | Peterson et al. (2012); Cohen et al. (2014) |
| [Neural construction](examples/neural-construction/README.md) | Grows cells, neurites and contacts from one founder and a shared artificial DNA recipe; tests signal transmission and contact lesions on the developed network | Zubler et al. (2011, 2013) |
| [Neural evidence](examples/neural-evidence/README.md) | Reconstructs 109 published CHE-1 depletion observations, traces them to workbook cells and intended DNA edits, and compares competing explanations | Traets et al. (2021) |
| [Neural fate](examples/neural-fate/README.md) | An artificial DNA parameter recipe driving published regulatory dynamics | Balaskas et al. (2012) |
| [Neural measurement](examples/neural-measurement/README.md), [forecast](examples/neural-forecast/README.md) | Replay of repeat-forecast comparisons on neuronal reporter assays (measurement: BIRD model) | Majoros et al. (2020) |
| [Neural regulation](examples/neural-regulation/README.md) | H1 neural-induction reporter data and a failed predictive test | |
| [Neural binding](examples/neural-binding/README.md) | Genomic *Nkx2-2* sequence, a designed substitution and a failed binding benchmark | Liu et al. (2025) |
| [DNA-edit predictor](examples/neural-predictor/README.md) | Scores 270-base reference and edited pairs; see below | Siraj et al. (2026); Salomon et al. (2026) |

Each example documents its inputs, runtime, scientific scope and reproduction commands. Neural construction needs a JDK and `patch` (tested with Java 21).

## Case study: predicting variant effects

The question is whether sequence alone predicts what an unseen DNA alteration does in neural cells. We tested it with a frozen-protocol procedure: a plan with its decision gates is written and hashed before any outcome is read; the scorer is qualified on a synthetic positive and a synthetic null; each test is scored once; a label is earned only if every gate passes.

| Group | Tests | Result |
|---|---|---|
| Cross-study reporter effects (six neural studies) | 1–8 | Cross-study correlation of the Siraj-trained reference model about 0.047 (r² near 0.002). Models trained on other studies, a Sei ridge, a quadratic context head, an activity-transfer model, a convolutional network and a zero-shot binding model do not improve on it. The delivered Huber-linear model reaches 0.0490 against 0.0465, with a gain interval that includes zero. |
| Direction of eQTL effects (fetal and neural-derived cells) | 9–11 | Direction accuracy of the best model is 0.52–0.65 across six held-out sets, against 0.50–0.60 for the best baselines. The largest gain did not repeat in the next set. |
| Causal variant in a credible set | 12–13 | The nearest-TSS rule (0.460) beats every model. |

The records, with a script that verifies their hashes, are in [`research/`](research/README.md). Limits of the protocol, including that the hashes are local and that the tests were not corrected for multiplicity, are in Section 5 of the paper.

### Delivered predictor

The delivered model is the Huber-linear checkpoint. On 16,986 internal-tuning edits in 1,759 families, equal-family mean squared error falls from 0.29661 (k5) to 0.28767, a 3.02% reduction, and ten-edit signed utility rises from 0.857 to 1.483. Six of seven original requirements pass; the ten-edit gain lower bound against k5 is −0.150. These data were used throughout development. On five external neural reporter studies its correlation is 0.0490, and its mean squared error is 9.52% above predicting no effect. It is a research prototype for prioritizing edits in the Siraj assay, not a validated predictor. See [the predictor README](examples/neural-predictor/README.md). To use it interactively, install [Epigenesis Studio](studio/README.md) and run `epigenesis-studio`; Studio also runs the neural-patterning simulation on a chosen DNA site.

## Measurement result

Before learning DNA effects, the effects need reliable measurements. A [BIRD-based forecast](https://doi.org/10.1093/bioinformatics/btz545) predicts an omitted RNA repeat from the other repeats in a fixed neuronal-assay pilot:

| Twelve omissions across four TRAIN components | Mean CRPS loss, lower is better |
|---|---:|
| Pooled Beta-binomial | 0.06219211 |
| BIRD repeat model | 0.05490724 |
| Relative loss reduction | 11.7% |

The conditional 99% simulation interval for the gain is [0.00436489, 0.00975665]. All four component means improve; one omission favors pooling and remains included. The interval describes numerical integration uncertainty for this selected cohort, and BIRD and pooling also differ in priors and DNA use, so the comparison does not isolate repeat variation as the cause. [Replay the archived calculation](examples/neural-measurement/README.md) without training, new random draws or private paths.

Two follow-ups use different targets. A [temporal experiment](docs/research/NEURAL_PREDICTOR_ROADMAP.md#new-experiments-18-october-2026) forecasts later enhancer responses from four early measurements: across 351 excluded families, full-history ridge cuts error by 24% against persistence but matches a simpler averaging control. A [raw-count experiment](examples/neural-forecast/README.md) holds out sequence families and infection cultures and finds a 0.26% mean-squared-error improvement over zero change, with all four gain intervals including zero. Neither establishes transferable forecasting superiority.

## Verification

```sh
python3.11 -m unittest discover -s tests -v
python3.11 tests/verify.py minimal-example
python3.11 tests/verify.py universal-translation
python3.11 tests/verify.py boundary
PIP_NO_INDEX=1 python3.11 tests/verify.py wheel
PIP_NO_INDEX=1 python3.11 tests/verify.py sdist
python3 research/check_freeze.py
```

The [CI workflow](.github/workflows/ci.yml) defines the complete release checks, including independent source replay, contract attacks, resource and path safety, real-DNA compilation and both offline demonstrations. Package checks need the declared build tooling. The predictor and its feature extractors need NumPy and, for some, PyTorch; each example README lists its requirements.

## Foundations and sources

The platform builds on established models and methods. Its own contribution is the executable connection between exact source records, reconstructed experimental inputs, interventions and independent verification. Full citations for every work cited are in the [technical paper](docs/paper/epigenesis-technical-paper.pdf); the main ones follow.

**Compiling biological descriptions.**
[Cai et al. (2009)](https://doi.org/10.1371/journal.pcbi.1000529) compile synthetic DNA with an attribute grammar, the closest direct precedent.
[Bilitchenko et al. (2011)](https://doi.org/10.1371/journal.pone.0018882) (Eugene), [Beal et al. (2011)](https://doi.org/10.1371/journal.pone.0022490), [Nielsen et al. (2016)](https://doi.org/10.1126/science.aac7341) (Cello) and [Wilson et al. (2016)](https://doi.org/10.1021/acssynbio.5b00194) compile from design languages to DNA or networks.
[Bergmann et al. (2014)](https://doi.org/10.1186/s12859-014-0369-z) package models and data for reproduction (COMBINE).
Identities use [SHA-256](https://csrc.nist.gov/pubs/fips/180-4/upd1/final) over [RFC 8785](https://www.rfc-editor.org/rfc/rfc8785) canonical JSON.

**Simulated neural development.**
[Zubler and Douglas (2009)](https://doi.org/10.3389/neuro.10.025.2009), [Zubler et al. (2011)](https://doi.org/10.3389/fncom.2011.00057) and [Zubler et al. (2013)](https://doi.org/10.1371/journal.pcbi.1003173) build the CX3D substrate and the instruction language that neural construction uses.

**Neural-tube patterning and binding.**
[Balaskas et al. (2012)](https://doi.org/10.1016/j.cell.2011.10.047), [Cohen et al. (2014)](https://doi.org/10.1242/dev.112573) and [Herrera-Delgado et al. (2020)](https://doi.org/10.1103/PhysRevResearch.2.043069) give the Sonic Hedgehog readout models. [Peterson et al. (2012)](https://doi.org/10.1101/gad.207142.112) supply the *Nkx2.2* enhancer evidence and [Hallikas et al. (2006)](https://doi.org/10.1016/j.cell.2005.10.042) the binding affinities. [Traets et al. (2021)](https://doi.org/10.7554/eLife.66955) provide the CHE-1 observations, and [Gong et al. (2020)](https://doi.org/10.3389/fnmol.2020.00034) the *Sox2* mutant used for coding consequences.

**Reporter assays and their statistics.**
[Melnikov et al. (2012)](https://doi.org/10.1038/nbt.2137), [Tewhey et al. (2016)](https://doi.org/10.1016/j.cell.2016.04.027), [Kircher et al. (2019)](https://doi.org/10.1038/s41467-019-11526-w), [Kreimer et al. (2022)](https://doi.org/10.1038/s41467-022-28659-0), [Deng et al. (2024)](https://doi.org/10.1126/science.adh0559) and [Agarwal et al. (2025)](https://doi.org/10.1038/s41586-024-08430-9) established and extended massively parallel reporter assays.
[Ashuach et al. (2019)](https://doi.org/10.1186/s13059-019-1787-z), [Myint et al. (2019)](https://doi.org/10.1186/s12864-019-5556-x) and [Majoros et al. (2020)](https://doi.org/10.1093/bioinformatics/btz545) supply the statistical models.
[Gneiting and Raftery (2007)](https://doi.org/10.1198/016214506000001437) justify the proper scoring rules (CRPS).

**Data used in the predictive tests.**
Training: [Siraj et al. (2026)](https://doi.org/10.1038/s41586-026-10121-6) and [Salomon et al. (2026)](https://doi.org/10.64898/2026.07.16.738760).
Evaluation: [Lee et al. (2025)](https://doi.org/10.1016/j.cell.2024.12.022) (GSE244011), [McAfee et al. (2023)](https://doi.org/10.1016/j.xgen.2023.100404) (GSE211045), [Matoba et al. (2026)](https://doi.org/10.1038/s41593-026-02454-2), [Rummel et al. (2023)](https://doi.org/10.1016/j.cell.2023.09.015) (GSE225815), [Guo et al. (2023)](https://doi.org/10.1038/s41588-023-01533-5) (GSE221384) and [Weiss et al. (2021)](https://doi.org/10.7554/eLife.63713) (GSE152404).
Motifs: [JASPAR 2024](https://doi.org/10.1093/nar/gkad1059).

**Sequence-to-function models.**
[Zhou and Troyanskaya (2015)](https://doi.org/10.1038/nmeth.3547), [Kelley et al. (2016)](https://doi.org/10.1101/gr.200535.115), [Kelley et al. (2018)](https://doi.org/10.1101/gr.227819.117), [Avsec et al. (2021)](https://doi.org/10.1038/s41592-021-01252-x), [Vaishnav et al. (2022)](https://doi.org/10.1038/s41586-022-04506-6) and [Chen et al. (2022)](https://doi.org/10.1038/s41588-022-01102-2) (Sei) are background for the sequence-model tests; Sei supplied the chromatin-profile features used in tests 3 and 9–13. [Seitz et al. (2024)](https://doi.org/10.1038/s42256-024-00851-5) and [Tareen et al. (2022)](https://doi.org/10.1186/s13059-022-02661-7) motivated the latent-activity model.

**Natural variation in brain.**
[GTEx Consortium (2020)](https://doi.org/10.1126/science.aaz1776), [Kerimov et al. (2021)](https://doi.org/10.1038/s41588-021-00924-w) (eQTL Catalogue), [Wang et al. (2020)](https://doi.org/10.1111/rssb.12388) (SuSiE), [Walker et al. (2019)](https://doi.org/10.1016/j.cell.2019.09.021), [Aygün et al. (2021)](https://doi.org/10.1016/j.ajhg.2021.07.011), [Jerber et al. (2021)](https://doi.org/10.1038/s41588-021-00801-6), [Schwartzentruber et al. (2018)](https://doi.org/10.1038/s41588-017-0005-8) and [Young et al. (2021)](https://doi.org/10.1038/s41588-021-00875-2) supply fine-mapped eQTL sets for the direction and credible-set tests. Expression and chromatin context come from [Kang et al. (2011)](https://doi.org/10.1038/nature10523), [Miller et al. (2014)](https://doi.org/10.1038/nature13185), [Roadmap Epigenomics (2015)](https://doi.org/10.1038/nature14248) and [Pollard et al. (2010)](https://doi.org/10.1101/gr.097857.109) (phyloP).

## License and attribution

Compiler code is [Apache-2.0](LICENSE). Bundled scientific data and third-party runtimes keep their documented licenses and attributions (see each example's NOTICE). The author is responsible for the content; the [paper](docs/paper/epigenesis-technical-paper.pdf) lists code, data and licenses in its availability section.

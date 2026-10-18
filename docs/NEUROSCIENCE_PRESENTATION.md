# Epigenesis: from DNA edits to testable neural mechanisms

Presentation draft, 20 October 2026. Evidence snapshot: 18 October 2026.
The fixed measurement and temporal studies are complete. This document separates
their measured results from the remaining genome-driven development objective.

Epigenesis makes neural gene-regulation experiments executable and inspectable.
A researcher can trace a measurement to its source, reconstruct the intended DNA
edit, run a declared mechanism, and choose an intervention that distinguishes
competing explanations. The long-term goal is genome-driven neural development.

## The problem we can solve today

A DNA motif is only part of a regulatory mechanism. Experimental constructs can
also differ in orientation, surrounding sequence and spacing. Treating such
constructs as a clean test of one factor can lead to an ambiguous next experiment.
Epigenesis exposes these differences alongside the original measurements and
provides executable mechanisms for testing the consequences of assumptions.

The main demonstration reconstructs 109 published gcy-22 mRNA observations from
ASER neurons after CHE-1 depletion. The observed auxin/control ratios of group
means are 8.70% for transfer of the binding core and 66.33% for the larger
regulatory fragment. Traets et al. published this contrast. Epigenesis adds
inspectable source reconstruction, intended donor geometry and explicit
competing explanations. [Traets et al., 2021](https://doi.org/10.7554/eLife.66955)

## Demonstration, approximately eight minutes

Generate the reports before presenting. Both use Python 3.11+, the standard
library and bundled sources. Choose new output directories.

```sh
python3.11 examples/neural-evidence/evidence.py demo --output /tmp/epigenesis-evidence
python3.11 examples/neural-evidence/evidence.py verify --output /tmp/epigenesis-evidence
python3.11 examples/neural-patterning/patterning.py demo --output /tmp/epigenesis-patterning
```

1. Open `/tmp/epigenesis-evidence/index.html`. Show all three promoter groups.
   Select one dot and trace it to its original workbook cell. The endpoint is
   mRNA count; the relative display divides by the group's control mean.
2. Inspect the donor reconstruction. The short construct replaces 12 bases
   with 12, with its core opposite the target gene. The larger donor replaces
   160 bases with a 130-base fragment, with its core in the same direction.
   This reconstructs intended donor products; final clones were not sequenced
   by this project.
3. State the experimental decision: use matched constructs to vary core
   orientation, surrounding sequence and spacing separately. Keep the insertion
   site and the other factors fixed within each comparison. Measure independent
   biological preparations. This proposed follow-up has not been performed.
4. Open `/tmp/epigenesis-patterning/index.html`. Compare a passive reporter and
   an endogenous edit, then examine signal history, allele reversion, mediator
   rescue and phenocopy. The same binding change can propagate differently when
   the altered gene participates in feedback.
5. Show the quantitative measurement result below, then the development map.
   End with the specific experimental comparison a laboratory could falsify.

## What the mechanism demonstration computes

A specified Nkx2.2 enhancer substitution changes a measured GLI binding-profile
ratio. The model transfers that ratio to intracellular affinity, then integrates
published Pax6/Olig2/Nkx2.2/Irx3 regulatory dynamics. That transfer is an explicit
assumption. The model's relative concentrations and the paper's embryo reporter
counts are different observables.

Competitive Gli activation and repression explain why an affinity reduction can
raise expression in one signal context and lower it in another. This mechanism
was derived by Cohen et al.; Epigenesis reconstructs it and connects it to the
specified edit and intervention controls. [Cohen et al., 2014](https://doi.org/10.1242/dev.112573)
[Model inputs and full citations](../examples/neural-patterning/README.md)

## A completed measurement forecast

The BIRD-based pilot asks whether the other RNA repeats can forecast a repeat
that has been omitted from fitting. It uses four selected TRAIN genomic
components from Rummel's neuronal assay, with three omissions per component.

| Method | Mean continuous ranked probability score (CRPS) |
|---|---:|
| Pooled Beta-binomial | 0.06219211 |
| BIRD repeat model | 0.05490724 |
| Relative reduction in forecast loss | 11.7% |

The gain's conditional 99% simulation interval is [0.00436489, 0.00975665].
One omission favors pooling and remains included. The interval measures
integration uncertainty on this selected cohort. Admission used all repeats;
the result therefore needs a broader, independently admitted cohort before
generalization. The comparators differ in priors, DNA use and heterogeneity.

[Portable calculation replay](../examples/neural-measurement/README.md)
recomputes all twelve scores from saved moments without training or random draws.
This is a measurement forecast using assay counts. The sequence predictor's
separate acceptance requirements remain unmet. A fresh DNA-normalized forecast
passes both pooling comparisons with a 15.8% CRPS reduction, but a post hoc
zero-effect forecast scores better. The treatment-response comparison remains
inconclusive. Early-to-late neural reporter forecasting now runs across 351
excluded families: early measurements help, while a simple averaging control
matches full-history ridge and beats both Timestral sizes.

The subsequent raw-count test excludes both sequence families and infection
cultures. The candidate's MSE is 0.32132339 versus 0.32214740 for zero change,
a 0.26% point gain. All four gain intervals include zero, and a simpler model
wins in one culture. Independent replay confirms this failed acceptance result
across all 15 evaluation groups. It scores 11,696 complete curves and retains
164 eligible forecasts with incomplete outcomes. This is a different observable
and target from the pooled-data experiment, so the two losses cannot be compared
directly. [Replay the complete comparison](../examples/neural-forecast/README.md).
[Research record](NEURAL_PREDICTOR_ROADMAP.md#new-experiments-18-october-2026)

BIRD is the published statistical model. Our contribution here is its numerical
implementation and this specific evaluation.
[Majoros et al., 2020; online 2019](https://doi.org/10.1093/bioinformatics/btz545)
[Rummel et al., 2023](https://doi.org/10.1016/j.cell.2023.09.015)

## Development map

| Link in the long-term path | Current implementation | Next scientific test |
|---|---|---|
| Genome and starting conditions | Exact sequence records, coding annotations and explicit model inputs | Infer a useful regulatory mechanism from new sequences with calibrated uncertainty |
| Development | Four-gene patterning dynamics; separate artificial growth recipe | Predict measured time courses and interventions in independent preparations |
| Neural structure | Cells, neurites and contacts grow from one founder in CX3D | Validate geometry, cell types and connectivity against biological development |
| Neural activity | Signal transmission and contact-lesion tests on the developed graph | Establish biological electrophysiology and plasticity |
| Learning | Future research objective | Define and demonstrate a learning task in the developed system |

These are working experiments along the path. They are not yet one validated
natural-genome-to-brain model. The artificial construction recipe encodes fifteen
supplied parameters in 480 bases. Gene-like neural self-construction has direct
precedents, including Zubler et al.'s single-progenitor G-code framework.
[Zubler et al., 2011](https://doi.org/10.3389/fncom.2011.00057)

## Scientific contribution and next proof

The release's contribution is the inspectable connection between exact source
records, explicit mechanisms, interventions and evidence, together with the
reconstruction and numerical work needed to make those connections reliable.
The biological mechanisms, growth substrate and BIRD model are attributed to
their original authors. [Prior-art comparison](PRIOR_ART.md)

The next biological proof is a predeclared intervention whose measured outcome
distinguishes two competing mechanisms. The matched promoter comparisons above
are a concrete candidate. Before running them, a laboratory must specify exact
constructs, biological replication, endpoints and acceptance criteria. They can
produce useful mechanistic evidence without requiring a general DNA predictor
or a complete digital brain first.

The current sequence-model ledger records 445 fits. Frozen inference remains
available; its joint predictive acceptance rule is unmet and confirmation stays
unopened. The completed follow-up studies establish a bounded repeat-forecast
gain and identify where stronger controls remove apparent forecasting gains.
Reliable sequence-to-effect prediction and biological-development validation
remain open. The source-linked intervention workflows above are usable now.

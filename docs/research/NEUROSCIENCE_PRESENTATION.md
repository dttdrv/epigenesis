# Epigenesis: from DNA edits to testable neural mechanisms

Public demonstration guide, updated 20 October 2026. Use
[Epigenesis Studio](../../studio/README.md) for the live demonstration and the
technical paper in [English](../paper/epigenesis-technical-paper.pdf) or
[Bulgarian](../paper/epigenesis-technical-paper.bg.pdf) for the scientific account.
The order is **Simulate → Prediction → Receipts**.

The fixed measurement and temporal studies are complete. This document separates
their measured results from the remaining genome-driven development objective.

Epigenesis makes neural gene-regulation experiments executable and inspectable.
A researcher can trace a measurement to its source, reconstruct the intended DNA
edit, run a declared mechanism, and choose an intervention that distinguishes
competing explanations. The long-term goal is genome-driven neural development.

The decisive result is still ahead: predicting the downstream neural effect of an
unseen natural DNA alteration, better than simpler alternatives, on independently
collected data. No experiment here meets that bar. The platform's distinct
contribution is a compiler that combines several DNA source formats, explicit
interpretation, a typed Development Module and independent replay, with executable
experiments from measurement to activity built on it. In the prior work reviewed in
[the prior-art record](../reference/PRIOR_ART.md), no system combines those compiler boundaries.
That statement covers the cited papers, not all software.

## The problem we can solve today

A DNA motif is only part of a regulatory mechanism. Experimental constructs can
also differ in orientation, surrounding sequence and spacing. Treating such
constructs as a clean test of one factor can lead to an ambiguous next experiment.
Epigenesis exposes these differences alongside the original measurements and
provides executable mechanisms for testing the consequences of assumptions.

The primary demonstration changes a supported GLI binding site and computes
how a declared four-gene mechanism responds under different starting conditions
and signal histories. It includes the reference, variant, passive reporter and
reversion controls in the same run. This makes the effect of the assumed
binding change visible and reproducible.

A separate evidence demonstration reconstructs 109 published gcy-22 mRNA observations from
ASER neurons after CHE-1 depletion. The observed auxin/control ratios of group
means are 8.70% for transfer of the binding core and 66.33% for the larger
regulatory fragment. Traets et al. published this contrast. Epigenesis adds
inspectable source reconstruction, intended donor geometry and explicit
competing explanations. [Traets et al., 2021](https://doi.org/10.7554/eLife.66955)

## Demonstration, approximately eight minutes

Install Studio from the repository root before presenting (Python 3.11+ and
[uv](https://docs.astral.sh/uv/)):

```sh
uv tool install --editable studio --with-editable .
epigenesis-studio studio/epigenesis_studio/example-edits.json
```

1. Switch the interface to Bulgarian if desired and open **Simulate**. Compare
   reference `GACCACCCA` with variant `GAACACCCA`, using GLI2. Explain that the
   module stores the selected bases; the simulator reads them and derives the
   binding ratio from the measured table.
2. Run the signal protocol and inspect Pax6, Olig2, Nkx2.2 and Irx3 over time.
   Compare the network edit with the passive reporter. Change the signal
   history or starting state and rerun. Explain the affinity-transfer
   assumption and the meaning of model units.
3. Show the reversion control, download the simulation receipt and verify it
   in **Receipts**.
4. Open **Prediction → Edits**, use the supplied 64-edit table and set panel
   size to three. Inspect the signed scores and family-diverse shortlist,
   download that receipt and replay it. This is the separate Huber-linear
   reporter predictor; its scores do not drive the simulation.
5. Show the internal and external results below, then the development map.
   The next biological test measures interventions that distinguish mechanisms.

For the separate source-reconstruction demonstration, generate its offline
report with `python3.11 examples/neural-evidence/evidence.py demo --output /tmp/epigenesis-evidence`
(choose a new output directory). Follow the [evidence guide](../../examples/neural-evidence/README.md)
to trace individual observations to source cells and compare intended donor edits.
Matched constructs that vary orientation, surrounding sequence and spacing
separately are a proposed laboratory follow-up, not an experiment performed here.

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
[Model inputs and full citations](../../examples/neural-patterning/README.md)

## The current DNA-edit predictor

The delivered model is now **Huber-linear**: joint k5 and 532 compressed motif
potentials with a fixed robust loss. On the reused internal population of
16,986 edits and 1,759 families, equal-family MSE is 0.287667 versus k5's
0.296611, a 3.015% relative reduction. Utility at 50/200 selections is
0.986330/0.634112 versus 0.606767/0.508652. Six of seven original requirements
pass; the ten-edit lower-bound requirement remains open. External mean
correlation is 0.049 and MSE is 9.52% worse than predicting zero.

Use the [Studio workflow](../../studio/README.md#predictor-modes),
[standalone predictor](../../examples/neural-predictor/README.md#current-model-huber-linear)
and [frozen transfer result](../../research/frozen-tests/huber-transfer-20261020/RESULTS.md).
The operational example scores 64 edits from 11 families without
training; its three-edit panel illustrates the workflow, not new validation.
Siraj RESERVE remains sealed; the earlier Salomon confirmation is already spent.

### Earlier predictor experiments (historical sequence)

The lead-model descriptions below record the decisions at those earlier stages.
“Current” and “lead” in this historical section refer to those stages;
Huber-linear supersedes them in the delivered presentation. Local `build/`
paths name unpublished working records. Public results are indexed in
[`research/`](../../research/README.md).

The current research problem is to predict which DNA edits change a neural-cell
reporter's expression and prioritize a small experimental panel. The current
Siraj lead model adds motif context to a fixed k5 model and meets the predeclared
MSE-plus-50/200-selection development milestone. It is measured on the same
internal-TUNING population of 16,986 edits in 1,759 families used in earlier trials.

| Internal TUNING | k5 | Signed-log pair | Precision-weighted pair | Logarithmic potential | Motif context |
|---|---:|---:|---:|---:|---:|
| Equal-family MSE | 0.29661090 | 0.29144854 | 0.29257851 | 0.29178267 | **0.29130587** |
| Relative MSE reduction vs k5 | — | 1.740% | 1.359% | 1.628% | **1.789%** |
| Signed utility, 10 families | 0.85730 | 1.74277 | 1.74277 | 1.60123 | 1.28576 |
| Signed utility, 50 families | 0.60677 | 0.92245 | 0.92506 | 0.95974 | **0.99326** |
| Signed utility, 200 families | 0.50865 | 0.57264 | 0.63273 | 0.59432 | 0.58223 |

The context model's nominal gain lower bounds versus k5 are positive for MSE
(+0.000202), 50 selections (+0.153349) and 200 selections (+0.015840). Its
ten-selection lower bound remains negative (-0.581819), so the original joint
acceptance rule remains unmet. MSE improves against word and scrambled controls,
but the incremental gain over the earlier motif model or a duplicate-feature
control is not established. All eighteen fits and all 45 score comparisons
passed independent verification. This is a concrete internal development
milestone for prioritizing larger experimental panels.
Context protocol, result and verification (local record: `build/neural-predictor-reviewed-20261005/siraj-context-pwm-20261019/RESULTS.md`)

A fixed smooth-context follow-up reduces the extreme-feature error found in
FIT, but does not replace this result. Its TUNING MSE is 0.29204308 and its
10/50/200 selection utility is 1.601227/1.025679/0.586091. Both combined gates
fail. The twelve-fit comparison took 29.3 seconds on one CPU thread; all 51
score contrasts and saved predictions passed independent replay.
Follow-up and numerical checks (local record: `build/neural-predictor-reviewed-20261005/siraj-smooth-context-20261019/RESULTS.md`)

A fixed combination of the linear, context and smooth models has also been
tested and independently replayed. Its MSE is 0.29183895, behind context's
0.29130587; both combined gates fail. It used saved predictions and one CPU
thread, with no base retraining. Context remains the presentation's lead result.
Combination result (local record: `build/neural-predictor-reviewed-20261005/siraj-potential-stack-20261019/RESULTS.md`)

A fixed cross-motif feature test also leaves context as the lead. Its MSE is
0.29165827 versus context's 0.29130587; both combined gates fail, and added
predictive value over its redundant-feature control is unsupported. All twelve
saved-head predictions and 57 score contrasts passed independent replay. This records a tested
interaction hypothesis, rather than evidence that interactions are irrelevant.
Cross-profile result (local record: `build/neural-predictor-reviewed-20261005/siraj-crossprofile-pwm-20261019/RESULTS.md`)

Using measured activity to choose one motif interaction gives MSE 0.29129777,
almost unchanged from context's 0.29130587. It retains the additional 50/200-edit
result but misses the original ten-edit requirement; the incremental benefit
over context and the paired-label control is unsupported. All 24 fitting stages
and 33 score contrasts passed independent replay. Context remains the lead.
Activity-supervised interaction (local record: `build/neural-predictor-reviewed-20261005/siraj-supervised-crossprofile-20261019/RESULTS.md`)

A further fixed test adds motif-spacing information and a matched shuffled-position
control. Its MSE is 0.29269607, worse than the context lead; both acceptance
gates fail. All twelve fits and 57 score contrasts passed independent replay.
This particular spacing feature does not improve the completed result above.
Spacing test and controls (local record: `build/neural-predictor-reviewed-20261005/siraj-spatial-pwm-20261019/RESULTS.md`)

A same-cell activity-supervision follow-up uses 133,208 FIT-only endpoint
measurements to learn sequence filters. It reaches MSE 0.291247, 1.81% below
k5 on the reused TUNING population, but its incremental gain over context is
unsupported and shuffled activity labels give a slightly better point MSE.
The top-ten selections and utility are unchanged. The additional development
milestone passes; the original joint requirement remains unmet. Context stays
the lead. All 24 stages and 33 score contrasts pass independent verification.
Separately, the network reduces held-family activity-prediction MSE by 8.70%
against a fitting-mean baseline, versus 0.78% after label shuffling. It learns
activity signal, but this transfer has not established an improvement in
edit-effect prediction.
The 8.70% figure concerns activity levels, not our edit-effect acceptance target.
Measured allele-activity differences correlate 0.928 with the edit-effect
target; the network's excluded-family scalar differences correlate only 0.074.
Predicting activity levels and predicting the effect of a DNA edit therefore
need separate evaluation. The measured-difference comparison uses outcome
labels and is not a sequence-only prediction result.
Activity supervision and control (local record: `build/neural-predictor-reviewed-20261005/siraj-activity-filters-20261019/RESULTS.md`)

A fixed continuation then trains those filters directly on edit-effect errors
left by the context model. Real-label and shuffled-label MSE are essentially
tied: 0.2912454682 and 0.2912454885. The incremental gain over frozen filters
is unsupported, and the same top-ten edits remain selected. This controlled
test closes without replacing the context lead. All 24 stages complete in
9.23 minutes on one CPU thread; all 39 score contrasts pass independent replay.
Direct edit-effect training and verification (local record: `build/neural-predictor-reviewed-20261005/siraj-contrast-filters-20261019/RESULTS.md`)

A subsequent saved-checkpoint check finds that replacing the trained scalar
with ridge did not hide an aggregate FIT gain: scalar MSE is 0.34244662 versus
0.34243887 for ridge, with identical top-10/50/200 selections. Independent
replay confirms the comparison; no new fitting or TUNING access was needed.
Readout comparison (local record: `build/neural-predictor-reviewed-20261005/siraj-readout-diagnostic-20261019/RESULTS.md`)

The completed other-cell comparison adds 132,829 FIT-only HepG2/K562 effects
through separate training heads. Its MSE is 0.291321 versus AP-only training's
0.291245 and the shuffled control's 0.291257; the top ten remain unchanged.
It establishes no incremental benefit, so context remains the lead. The fixed
run took 7.80 minutes on one CPU thread, and independent replay verifies all
45 comparisons and 36 ridge equations. No GPU or new spending was used.
Auxiliary-cell comparison (local record: `build/neural-predictor-reviewed-20261005/siraj-multicell-effects-20261019/RESULTS.md`)

The three earlier motif models have a positive nominal family-bootstrap utility-gain
lower bound at 50 selections. Precision weighting and the logarithmic sequence
potential also have positive bounds at 200. None passes the required combined
MSE/top-ten rule. The signed-log pair model's MSE-gain lower
bound is -0.000131; its ten-selection utility-gain lower bound is -0.123633.
The population has supported adaptive development, so these intervals are not
fresh or selection-adjusted confirmation.

The useful technical finding is specific: rare, extremely large motif values
dominated errors in the raw representation; compressing them improved the
observed prediction and selection comparisons. Twelve reconstructed heads and
all 42 score comparisons passed independent checks. This verifies the reported
experiment and gives a concrete model-improvement lead. Independent families
must establish whether the gain generalizes.

The signed-log pair model predicts pairs of sequences and lacks a guaranteed common
activity potential: its predictions need not sum to zero around an edit cycle.
A separate bounded per-window model preserves that constraint but raises MSE
to 0.29820935, 0.539% worse than k5; its utility-gain bounds cross zero at all
three panel sizes. Independent replay confirms this failed comparison. The
completed logarithmic per-window model recovers the gain while preserving that
sequence potential. At 50 and 200 selections its utility-gain lower bounds are
positive against k5, word and scrambled-motif controls. Its MSE bound versus k5
is -0.000594 and its ten-selection bound is -0.147759. This is a useful lead for
prioritizing neural-cell reporter experiments, with independent confirmation
still required. Independent checks verified all twelve fit records, the two
final head equations and every score comparison.
Exact protocol, search history and replay (local record: `build/neural-predictor-reviewed-20261005/siraj-signedlog-replay-20261019/RESULTS.md`)
Bounded comparison (local record: `build/neural-predictor-reviewed-20261005/siraj-bounded-pwm-20261019/RESULTS.md`)
Logarithmic sequence-potential comparison (local record: `build/neural-predictor-reviewed-20261005/siraj-logpotential-pwm-20261019/RESULTS.md`)

A fixed uncertainty-weighted follow-up retains the positive 50/200-selection
bounds but raises MSE to 0.29302485. It does not replace the unweighted potential's
stronger error result. All twelve fits and independent replay are complete.
Follow-up comparison (local record: `build/neural-predictor-reviewed-20261005/siraj-logpotential-precision-20261019/RESULTS.md`)

An analysis restricted to the saved FIT predictions finds that all four
compressed recipes improve MSE in every fixed fold and every reported-SE
quartile. Their largest gains occur in the lowest-SE quartile. Compression also
reduces the largest adverse family's share of adverse error from 33.85% to
roughly 1.4–1.9%. Remaining losses are distributed across families, and the
precision model has little headroom from a global rescaling. These are adaptive
training diagnostics, not a measured noise ceiling or independent validation.
FIT-only diagnosis (local record: `build/neural-predictor-reviewed-20261005/siraj-fit-residual-diagnosis-20261019/REPORT.md`)

For the presentation: "We found and reproduced a concrete failure mode in
motif-based neural reporter prediction. A compact CPU model now meets our
declared internal MSE and 50/200-edit panel milestone. Ten-edit selection remains
unstable; the original acceptance rule and independent confirmation are still
ahead."

### Earlier independent confirmation, a different assay

The bundled activity-transfer model learns a short-word score from 30,800 TRAIN reporter
measurements, then transfers it to edit effects with one calibrated coefficient.
It passes both development requirements. The frozen model was then tested once
on 7,372 unseen edits across 2,938 sequence families.

| Confirmation | Family-mean MSE | Signed ten-family utility |
|---|---:|---:|
| Activity-transfer predictor | 0.02895676 | 0.16957315 |
| Fixed k5 baseline | 0.02899805 | 0.15278638 |

All seven point comparisons improve, but every paired family-bootstrap interval
crosses zero. The original acceptance rule remains unmet. This identifies a
useful direction for further predictor development without resolving statistical
superiority. The model and a ten-family experiment-selection command run locally;
the exported implementation reproduces all 14,762 saved development and
confirmation predictions exactly. [Run it and inspect the evidence](../../examples/neural-predictor/README.md#activity-transfer-predictor).

MSE uses squared log2 reporter-effect units. Selection utility is the average
observed effect multiplied by the model's predicted sign for its ten choices.
Neither number is a percentage error. The experiment measures reporter responses
in the specified WTC11-NGN2 context. [Salomon et al., 2026](https://doi.org/10.64898/2026.07.16.738760)

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

[Portable calculation replay](../../examples/neural-measurement/README.md)
recomputes all twelve scores from saved moments without training or random draws.
This is a measurement forecast using assay counts. The sequence predictor's
separate acceptance requirements remain unmet. A fresh DNA-normalized forecast
passes both pooling comparisons with a 15.8% CRPS reduction, but a post hoc
zero-effect forecast scores better. The treatment-response comparison remains
inconclusive.

A separate Kreimer H1 neural-induction experiment uses four early observations
at 0/3/6/12 hours to forecast 24/48/72 hours across 351 excluded sequence
families. Its equal-family MSE is 0.06015 for a two-feature ridge using the early
average and 12-hour value, 0.06034 for full-history ridge, 0.06534 for Timestral
5M, and 0.07504 for Timestral 209M. This is retrospective forecasting under
published all-time preprocessing. The simple control performs best in this
comparison; neither Timestral size passes the joint requirements.

The Timestral adapter preserves negative and zero measurements and masks the
missing 9-hour slot. Its normalization uses observed-context mean and sample
standard deviation, with a 1e-6 floor. Subsequent independent synthetic checks
verified zero and constant series, signed variation below/at/above the floor,
sample-variance calculations, masking and finite normalization roundtrips.
Those checks establish correct behavior on the tested boundaries; they do not
explain the larger model's forecasting error.
All 21 forecast quantiles were saved; the comparison above scores their median.
The completed temporal experiment has no DNA feature. Adding one requires a
separate assay-matched, nested experiment; Siraj SK-N-SH and Kreimer H1 have
different endpoints. Fixed temporal comparison (local record: `build/neural-predictor-reviewed-20261005/neural-temporal-forecast-20261018/RESULT.md`)
Normalization checks (local record: `build/neural-predictor-reviewed-20261005/neural-temporal-forecast-20261018/reviews/scaler-boundary-20261019/REPORT.md`)

The subsequent raw-count test excludes both sequence families and infection
cultures. The candidate's MSE is 0.32132339 versus 0.32214740 for zero change,
a 0.26% point gain. All four gain intervals include zero, and a simpler model
wins in one culture. Independent replay confirms this failed acceptance result
across all 15 evaluation groups. It scores 11,696 complete curves and retains
164 eligible forecasts with incomplete outcomes. This is a different observable
and target from the pooled-data experiment, so the two losses cannot be compared
directly. [Replay the complete comparison](../../examples/neural-forecast/README.md).
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
their original authors. [Prior-art comparison](../reference/PRIOR_ART.md)

The next biological proof is a predeclared intervention whose measured outcome
distinguishes two competing mechanisms. The matched promoter comparisons above
are a concrete candidate. Before running them, a laboratory must specify exact
constructs, biological replication, endpoints and acceptance criteria. They can
produce useful mechanistic evidence without requiring a general DNA predictor
or a complete digital brain first.

The latest motif comparisons are exploratory development on the Siraj assay.
The earlier activity-transfer confirmation remains unchanged and will not be
reused for tuning. The completed follow-up studies establish a bounded
repeat-forecast gain and identify where stronger controls remove apparent
forecasting gains. Fit accounting distinguishes the original model ledger,
the separately audited scratch search and numerical reproductions; a single
stale total must not stand in for that history.
Reliable sequence-to-effect prediction and biological-development validation
remain open. The source-linked intervention workflows above are usable now.

# Neural DNA-edit predictor roadmap

Updated 2026-10-18. The new fixed measurement and temporal comparisons are
complete. They resolve the earlier four-source pooling comparison and test
early-to-late forecasting, including Timestral-5M, Timestral-209M and raw-count
transfer across infection cultures. The latter test does not demonstrate
superiority over its controls. The bundled sequence predictor, its weights and
its unmet acceptance requirements remain unchanged.

## New experiments, 18 October 2026

A fresh DNA-normalized effect comparison completes all eight cross-barcode
transfers on the original four TRAIN components. BIRD CRPS is 0.51127630 versus
0.60712982 for pooling, a 15.79% reduction. Both predeclared directional gain
intervals are positive: [0.08694432, 0.11724122] and [0.05343141, 0.12160532].
These are conditional 99% Monte Carlo bounds, not biological confidence bounds.
Independent exact-moment replay reproduces every score. A post hoc zero-effect
forecast has lower CRPS, 0.42482818; BIRD central 80% empirical intervals cover
only two of eight observations. The fixed pooling comparison passes, while
practical effect-prediction superiority and calibration remain unestablished.
The original failed and unresolved attempts remain unchanged.

The temporal experiment uses measured effects at 0, 3, 6 and 12 hours to forecast
24, 48 and 72 hours in the Kreimer neural-induction assay. Five whole-family
outer folds cover 4,143 targeted perturbations in 351 existing TRAIN families.
Scales and tuning use only each training complement. Original validation and
test rows are not scored. The publication's pooled activity estimates used
all-timepoint processing, so this is a retrospective exploratory result.

| Forecast | Family-mean MSE | Within-family selection utility |
|---|---:|---:|
| Carry the 12-hour value forward | 0.07963470 | 0.19753080 |
| Early average and 12-hour value, ridge | 0.06015096 | 0.32444986 |
| All four early values, ridge | 0.06033572 | 0.32323640 |
| Timestral-5M median, zero-shot | 0.06533556 | 0.30991853 |
| Timestral-209M median, zero-shot | 0.07504111 | 0.30814769 |

Utility selects the predicted top 10% of changes within each family, rounded up,
and measures their observed RMS change from 12 hours. Random selection has
expected utility 0.19753080. Equal-score boundary ties receive fractional credit.
All seven error controls and four selection controls were fixed before the ridge
run. Neither full history nor either Timestral size passes all required controls.
The two-feature control matches the richer model. Early measurements have useful
predictive information here; greater model size has not added demonstrated value.
Both Timestral sizes use the same three-hour grid, missing 9-hour input, signed
effects and actual late-time offsets. Inference takes 0.67 seconds for 5M and
10.11 seconds for 209M on one CPU thread, with no GPU.

Independent saved-array scoring reproduces all three candidates and their
family bootstrap intervals. The temporal analysis adds 420 small nested ridge
fits and two zero-shot inference runs; the separate 445-fit sequence-model
ledger and bundled predictor weights remain unchanged. Full protocols, source
identities, controls and results are preserved locally under
`build/neural-predictor-reviewed-20261005/bird-effect-complement-followup` and
`build/neural-predictor-reviewed-20261005/neural-temporal-forecast-20261018`.

The subsequent raw-count experiment holds out both infection cultures and
sequence families. Early measurements from the other two cultures predict
changes from 12 hours to 24, 48 and 72 hours in the omitted culture. Its 12-hour
anchor is label-only. A new same-sample pooled RNA/DNA contrast avoids an
all-timepoint activity fit; upstream matrix selection remains unresolved.

| Raw-count change forecast | Family-balanced MSE |
|---|---:|
| Early mean and 12-hour value, ridge | 0.32132339 |
| Early-mean ridge | 0.32194217 |
| Latest-value ridge | 0.32211810 |
| Zero change | 0.32214740 |
| Training-population mean change | 0.32201750 |

All four simultaneous nominal 99% gain intervals contain zero; the early-mean
control also wins in culture 2. The fixed conditional acceptance rule is not met.
All 15 culture/fold groups are observed. Of 11,860 input-eligible forecasts,
11,696 have complete targets and 164 remain unscored. Loss averages variants
within family/culture cells, available cultures within families, then families.
The intervals resample whole families and are conditional on fitted models,
supported outcomes and these three cultures, not preparation-population bounds.

The [portable archive](../examples/neural-forecast/README.md) retains count
totals, support reasons, the frozen protocol, model coefficients, prediction
seals and all scores. Independent replay agrees within 1.78e-15. Before numerical
count conversion, 23 synthetic tests passed and 12 deliberate faults were
detected. The run adds 945 small ridge fits, separately from the 420 pooled-alpha
fits and 445 sequence-model fits. It completed in 8.99 seconds on one CPU thread.

This is one previously exposed H1 assay with three infected cultures and two
viral batches whose assignment is unavailable. It is not donor or viral-batch
validation. The pooled-alpha and raw-count experiments use different observables,
targets and exclusion schemes; their losses and relative gains are not a matched
before/after comparison. The new result does not establish why the earlier
pooled-data gain occurred.

## Earlier predictor and measurement experiments

The Matoba neural-progenitor experiment asks how a DNA edit changes the response
to Wnt stimulation. Its first bounded model comparison is complete. Rummel's
complementary neuron-stimulation experiment now has a fixed source partition and
TRAIN measurements. Its aggregate KCl interaction has weak RNA-view agreement
and closes at zero fits. A separate baseline-neuron effect screen passes both
fixed TRAIN comparisons: k5 slightly improves prediction error and improves
signed edit-selection utility. Independent replay confirms this exploratory
result. The subsequent independently tuned dinucleotide baseline beats k5 on
error; k5 retains its selection advantage and fails the combined progression
rule. The fixed-prediction error diagnosis is independently reproduced. The matched
anchor/residual comparison now has lower MSE but loses selection utility to
stronger direct k5. Ten further fits complete without a joint pass; independent
numerical replay and the root rerun reproduce the result. The subsequent fixed
full-insert motif comparison completes fifteen fits and fails both endpoints
against direct k5; independent saved-result replay and the root rerun agree.
A subsequent zero-fit audit finds that the strongest measured baseline edits
transfer positively between existing barcode/RNA views in both directions and
beat the fixed random-selection reference. Independent exact replay and the
root rerun reproduce both directions and all chromosome omissions.
A fixed crossover of the existing direct and anchored model selections now
locates the anchor's full-target utility loss entirely in changed selections;
both models retain identical signs on those picks. Direct k5's original picks
also exceed their model-sign random reference in both existing measurement views.
Independent exact replay and the root rerun agree. Three unavailable anchored
picks prevent a complete measurement-view comparison of the two models.
A five-fit joint k2/k5 comparison now improves error over direct k5 but reproduces
the anchored model's ten selections and signs. It fails the combined progression
rule; independent replay and the root rerun agree. This fixed recipe is closed.
A subsequent twenty-fit matched comparison adds 2,210 qualified same-study
training observations. Fixed k5 pooling worsens both error and edit-selection
utility against Dev-only training; all ten uncertainty bounds fail. Independent
saved-fit and score replay and both root reruns agree. The stronger Dev-only
model is retained. A subsequent source-defined comparison finds negative
cross-library agreement on 527 jointly eligible identical physical edits
(r = -0.1126), with both copy directions losing to zero. All 21 chromosome
omissions preserve these signs; an independent numerical replay agrees.
A fixed RNA/plasmid decomposition now reconstructs every saved label exactly.
Disagreement appears across the samples: eight of nine cross-library sample
correlations are negative, and the ninth is +0.0298. Independent high-precision
replay agrees. The author-method audit now confirms a concrete estimator
difference: the paper fits separate-library raw-count models after collapsing
barcodes, whereas this screen used fixed mean log ratios. Existing aggregate
counts are now exported and independently verified. Review of current MPRAnalyze
source rejects the proposed classic two-plasmid-observation setup before fitting.
A separate exact count-weighted contrast reduces mismatch by 6.7%, but agreement
remains negative and fails its fixed pooling-repair screen. Independent numerical
replay agrees. Original labels and the Dev-only predictor remain unchanged.
A fixed five-fold precision-weighting trial now raises signed selection utility
from 0.50658 to 0.84756 on the same 948 TRAIN development observations, a 67.31%
relative score increase. MSE is 0.01294% worse than Dev k5 and also loses to Dev
k2. Ten of twelve point gates and four of twelve paired interval gates pass;
the joint rule remains unmet. Independent numerical and score replay agree.
A fixed zero-fit audit now finds that this selection advantage also holds in
both existing measurement views, with utility gains of 0.43230 and 0.36292.
All ten selections have both views available. Two changed picks explain the
entire gain, with identical predicted directions on both lists. Independent
exact replay agrees. These related TRAIN measurements support the prioritization
result; the prior error loss and inconclusive utility interval remain open.
A subsequent nested combination of Dev k2 and precision k5 completes forty
inner sequence fits. Its MSE is 0.54542267 and utility 0.50184668, worsening both
against precision k5. Nine of fourteen point comparisons pass; none of the paired
interval lower bounds is positive. Independent saved-solution and score replay
and the root rerun agree. This combination closes without advancement.
A subsequent ten-fit quadratic-composition comparison improves both point
scores against all six prior controls: MSE 0.54323404 and utility 0.87841414.
Its constant-context control selects better (utility 1.18697653), with MSE
0.54331450. The candidate passes fifteen of sixteen point comparisons and four
paired interval gates, so the varying-context recipe closes without advancement.
Independent replay and the root rerun agree. A separate algebraic check shows
that the stronger selection control reduces to local k2/k5 edit deltas; its
result identifies useful readout structure without establishing a benefit from
varying sequence context. A fixed follow-up now attributes its entire selection
gain to one substituted edit. That edit lacks usable A/B views under the existing
QC, leaving incremental cross-view support unresolved. Independent exact replay
and the root rerun agree. Validation/test outcomes remain unopened.
A subsequent fixed count diagnostic rejects the common-RNA-fraction read model
for 438 of 948 eligible DevLib edits (46.20%), including the substituted pick.
Independent full replay agrees. Replicate variation must enter the measurement
uncertainty work before the pooled Swift sampler can serve as a default model.
The published BIRD replicate model now has a checked marginal density and an
independent-draw posterior sampler. Four fixed synthetic cases pass all 20
distribution comparisons against a separate generative reference. Subsequent
omitted-RNA forecasts improve CRPS by 11.7% in the fixed four-component pilot.
The separate DNA-effect comparison remains unresolved under its both-direction
rule. A full-count treatment-response extension now completes all 32 fixed
comparisons with 14.3% lower estimated mean CRPS than pooling. Its simulation
intervals remain too wide for the both-direction acceptance rule. Independent
replay verifies all 64 model results and the full aggregate. The original
CPU-limited attempt and the explicit completion amendment remain separate.
All 11,378 named Matoba design pairs have a
matching GRCh37 reference placement. Official component tables and assembly
alignments resolve 676 of the 680 alternate-reference cases to the same complete
primary interval and edit position. Four retain uncertainty and use conservative
genomic bounds for grouping. The one unnamed insert is explicitly excluded from
supervised use while remaining in the original design ledger.

The source-only partition contains 178 genomic families: 7,543 TRAIN pairs in
109 families, 1,645 validation pairs in 38 families, and 2,190 test pairs in 31
families. Independent replay reproduces every assignment, and the source
partition is frozen. Both approved count files are now fully imported and
verified. Fixed measurement QC retains 6,930 of the 7,543 TRAIN pairs across all
109 TRAIN families. The alternating-well correlations are 0.4971 for Vehicle
allele effects, 0.6315 for stimulated allele effects, and 0.1422 for their Wnt
interaction. These describe agreement conditional on shared barcode QC, not
independent donor replication or sequence-prediction accuracy. Validation and
test effects remain unexported. The first bounded TRAIN comparison is complete:
15 fit attempts across five excluded-family folds. Direct k5 ridge has MSE
0.01781801 and shared rank-one ridge 0.01776933, versus 0.01750798 for zero
effect. Sharing reduces error by 0.2732% relative to direct ridge but lowers
signed selection utility from 0.00071456 to -0.03394066. This sharing recipe
fails both progression requirements and is closed; independent result replay
passes and root reproduces it.

| Next step | Current state | Scientific purpose |
|---|---|---|
| Resolve designs and freeze genomic families/split | Complete: 11,378 source-eligible pairs in 178 families; one explicit source exclusion; independent replay matches every row | Keep related sequences and every observation of an edit together |
| Measure condition effects and their repeatability | Complete: 6,930 eligible TRAIN pairs; independent replay confirms all measurements and diagnostics; Wnt-interaction half agreement is weaker than either condition's allele agreement | Identify the measured signal and uncertainty before fitting |
| Train and compare a predictor | First bounded TRAIN comparison complete and independently replayed: shared ridge improves direct ridge's error slightly but worsens selection; both lose to zero on error | Beat simple baselines on both effect prediction and useful edit selection |
| Qualify complementary neural perturbation data | Rummel partition and TRAIN measurement complete: 948 eligible pairs in 800 families; the weak aggregate interaction route closes at zero fits | Test whether the assay supplies reproducible supervision before transferring or fitting it |
| Predict baseline-neuron regulatory edits | 130 Rummel predictor fits complete. Quadratic composition improves both point scores over six prior controls, but its constant-context control selects better; the full acceptance rule fails. Independent replay agrees | Improve prediction and useful selection together on the same complete population |
| Reuse full-insert motif features | Fixed fifteen-fit comparison independently replayed: full-insert motifs lose error to both matched controls and lose utility to scrambled motifs | Test sequence context that cancels out of raw k-mer deltas |
| Check whether strong measured edits repeat | Zero-fit audit independently replayed: positive transfer above random in both directions and all 22 chromosome omissions | Test measurement repeatability at the selection tail before another model recipe |
| Diagnose the model's selection loss | Zero-fit crossover independently replayed: anchor changes membership, with no sign changes on the selected lists; direct picks exceed model-sign random references in both available views | Identify the observed selection failure before choosing another learning objective |
| Qualify additional same-study training data | Measurement replay complete; 2,210 usable designs expand the training pool to 3,158 observations. Fixed pooling fails the matched predictor test | Source bindings and lossless count export pass. One count-weighting repair reduces mismatch but leaves negative agreement; keep the libraries separate |
| Confirm on untouched families | Pending a qualifying frozen candidate | Test generalization without another selection opportunity |
| Rank edits for follow-up experiments | Biological acceptance pending | Help allocate a limited experimental budget |

This addresses an early regulatory step in genome-driven development. Later
work must connect regulatory effects to developmental trajectories, neural
structure, activity and learning. The completed earlier predictor rounds below
remain part of the record, with their original outcomes and acceptance rules.

The earlier completed ledger is **300 biological fits**. Matoba adds **15 fit
attempts**, all completed: five complement means, five direct ridge and five
shared ridge. Rummel's baseline-neuron pilot adds **10 completed attempts**,
five k2 and five k5. Independent k2 tuning adds **15 completed attempts** and
reuses five prior fits. The matched direct/anchored comparison adds **10 completed
attempts** with five k2 reuses. The fixed full/overlap/permuted motif comparison
adds **15 completed attempts**. Joint k2/k5 adds **5 completed attempts**. The
matched added-data comparison adds **20 completed attempts**. Precision weighting
adds **5 completed predictor attempts**. Nested combination adds **40 completed
inner predictor attempts**. The varying/constant composition comparison adds
**10 completed attempts**, for **445 cumulatively**, plus 25 separately recorded
nuisance variance fits and five closed-form combination estimates. Zero prediction and independent
saved-result replay add no fits. The bundled
[two-cycle predictor](../examples/neural-predictor/README.md) is runnable on CPU
and Metal. Its validation utility is 0.07158510 and MSE is 0.02862935; utility
passes its comparator and MSE remains 0.01480% above its strict limit. The
TRAIN uncertainty audit and its seven-fit precision-sampling trial are complete.
The latter has MSE 0.02862093 and utility 0.00342057: it passes the error limit
but fails ranking, with no nomination. A subsequent six-fit blend of k3 and k5
has MSE 0.02863000 and utility −0.03630627, failing both gates. The subsequent
quadratic readout has MSE 0.02862839 and exactly the same signed selections as
k5, so it misses the MSE limit and ties ranking. Its raw floating-point utility
flag differs from the exact selection verdict; that discrepancy is retained.
Confirmation for the newer Salomon cohort remains unopened. The sections below retain earlier
rounds and their outcomes.

Read the [current baseline table](../examples/neural-predictor/README.md#validation-results)
with the latest result. The 0.01480% figure is a relative gap to k5 MSE, not a
prediction error rate. The bundled model improves MSE over zero effect by only
0.4256%; centered R² is −0.00555. Another baseline, k3 λ=0.01, has higher ranking
utility. Exact numerical replay and biological predictive accuracy are separate
results. The 196-fit ledger and MSE near 0.01554 below belong to the older Kosicki
rounds, whose measurements differ from the current Salomon cohort.

The original Kosicki reserved confirmation is complete and failed acceptance. The fixed
gapped-six predictor was evaluated on all 3,840 eligible edits in 50 families.
Its family-weighted MSE was 0.0169575704 versus 0.0169190102 for zero effect;
its ten-edit signed utility was -0.0379617141 versus +0.0225380178 for random
selection. All seven required improvement intervals crossed zero. No further
candidate selection or tuning can use this set as fresh confirmation.

Nine bounded development rounds preceded this test. The gapped-six model passed
both validation thresholds. The adopted, independently reviewed amendment first
retained eligible candidates passing both validation requirements, then selected
their minimum MSE. This post-validation nomination changed no weights, comparator
or confirmation threshold. The original minimum-MSE-first selection still fails
ranking and is preserved; the exact-word control remains ineligible. The latest
combination model reproduced gapped-six alone, adding no new predictor behavior.

A tenth development round is now complete: twelve compact motif models compared
local interactions and joint versus staged activity training across three fixed
seeds. None beats the required MSE baseline; no matched contrast improves both
endpoints in every seed. This brought the biological fitting ledger to 219. A
training-only diagnostic identifies strong task-gradient imbalance, and replay
identifies rounding-sensitive selection between mathematically tied local edits.
The full-length synthetic check is now complete: all nine runs pass every
recovery requirement. A separate inference repair resolves the recorded
rounding-sensitive structural ties. Neither adds a biological fit or changes
the failed predictive results.

An eleventh round compared static-balanced and paired-only objectives against
the existing joint checkpoints. All six new fits improve training MSE, while
all six worsen validation MSE and signed ranking. No candidate qualifies. The
biological fitting ledger is now 225. Independent replay verifies every saved
prediction and decision in this comparison.

A subsequent training-only diagnostic is complete. The exact-word and random
encoder controls have sharply different sequence-similarity structure despite
their similar historical validation errors. Their shared score does not justify
one common explanation or a new fit. The MPAC source audit recovered all ten
original design libraries and found exact matches to three validation-assigned
and five test-assigned families before measurement QC. Exact training membership
and upstream tuning exposure remain unresolved; MPAC is not admitted.
The outcome-free DNA-shape projection checks pass. Its frozen training-only
alignment diagnostic is now complete and independently reproduced: all three
progression screens failed. This pooled representation justified no new fit.

Independent replay verifies the completed rounds. Before the tenth round, the ledger contained 207 fits:
161 eligible candidate fits, 36 controls/baselines and ten auxiliary fold fits.
The ninth round's combination candidate duplicated the existing gapped model. The earlier
exact-word control reduces validation MSE by 0.467% against the required k-mer
baseline and has ranking utility 0.0251110 versus 0.0011730 for random selection.
These point estimates identify a follow-up lead. The earlier proposal to nominate
this control remains unadopted; the original eligibility rule is retained. Both required
sequence baselines now have independently checked training normal equations and
reproduced validation scores and rankings. All development rounds and their
accounting are preserved. Confirmation added no fits. The tenth round adds twelve
candidate fits, and the eleventh adds six: the ledger at that point is 179 candidate
fits, 36 controls/baselines and ten auxiliary fold fits, for 225 in total.

## Problem and deliverable

A researcher has candidate regulatory DNA edits and a limited experimental
budget. Predict each edit's effect on reporter activity in a specified human
neural preparation, then rank edits for follow-up measurements.

The deliverable accepts paired reference and edited sequences plus a supported
assay context. It returns predicted effects in declared units, a ranked table,
model and input identities, and measured validation performance. It must reject
unsupported sequence lengths or assay contexts. A confidence interval for average
benchmark performance is not a confidence interval for an individual prediction.

The existing target is H1 neural induction at 0, 3, 6, 12, 24, 48 and 72 hours.
The new experiment uses WTC11-NGN2 induced excitatory neurons, infected on day 7
and measured on day 14, from Kosicki et al. 2025. This establishes a separate,
explicitly named predictor if acceptance passes. It does not complete the H1
timecourse task.
Reporter activity, endogenous gene expression and neural cell identity are
separate endpoints.

## Starting evidence

The existing [neural-regulation experiment](../examples/neural-regulation/README.md)
has working data intake, sequence prediction, evaluation and compiler replay.
Its frozen model reduced mean squared error (MSE) by 20.97% against predicting
zero effect. It failed the required comparison against a dinucleotide model and
did not establish a benefit from time-specific predictions. The examined test
remains a historical result, not a fresh selection or confirmation set.

## Earlier Kosicki experiment roadmap

| Phase | Work and output | Required result before advancing | State |
|---|---|---|---|
| 1. Freeze the experiment | Select one accessible paired-edit dataset from the current shortlist. Write its assay contract, source hashes, sequence pairing, exclusions, split, endpoints and evaluation protocol. | Exact sequence and outcome formats are usable; held-out families can be isolated; evaluation outcomes have not informed model training or selection. | Complete; evidence below |
| 2. Build one successor | Use training and validation data to select a sequence model. Compare it with zero effect, the dinucleotide model and the existing k-mer approach under the same protocol. Save every attempted configuration. | Runnable predictor with fixed preprocessing, parameters and calibration; selection uses no held-out outcomes. | Nine rounds complete; existing gapped-six nominated through disclosed selection amendment |
| 3. Verify the implementation | Exercise data isolation, numerical calculations, input rejection, batch prediction and exact-source replay. Save the final candidate and baseline identities before scoring. | Independent reference calculations agree; deliberately broken controls fail; the independent review has no unresolved correctness or leakage finding. | Complete; statistics, inference, runner and final admission independently checked |
| 4. Evaluate once | Score the frozen candidate and baselines on the reserved data. Save all eligible predictions, exclusions, errors and uncertainty. | Both effect-prediction and ranking acceptance pass. A failure remains a failure; any subsequent model change requires a new confirmation set. | Complete; failed all seven comparison gates |
| 5. Deliver the experiment-selection tool | Expose batch prediction and ranking through a documented command, with a worked real-data example and reproducible report. | A clean run reproduces the accepted result and ranking; validation supports the stated assay and use. | Acceptance unmet; this candidate is not promoted |

The full-length positive control establishes trainability for its known function;
it does not settle the much larger biological task imbalance or representation
limits. Any further optimization comparison must test that specific biological
question prospectively. The repaired CPU reference inference path also needs
to be fixed before a new biological evaluation. Keep the desired single-model use:
predict effects and rank experiments from the same predictions. Development can
use exposed training/validation data; a new performance claim requires an unused
confirmation population and a frozen selection rule. Account for all exposure
to the completed test. The current evidence does not establish assay noise,
model capacity or gradient imbalance as the sole cause of predictive failure.

The experiment record below is chronological. Statements that outcomes were
unopened describe the corresponding earlier stage; the final confirmation
section records their subsequent evaluation.

The phase-1 literature search is closed. Kosicki was selected using experimental
design, access and split feasibility before downloading its measured outcomes.
A paper's reported absolute-activity correlation does
not establish accuracy for the difference between an edited and reference pair.

The frozen experiment and intake evidence are in
`build/neural-variants-20261005/`. The protocol was frozen before measurement
download. Exact sequences, strand-aware edits, connected families, previous H1
sequence/coordinate overlap and exposed published examples were checked. All
589 prior H1 intervals were mapped across genome assemblies before exclusion.

After measurement QC, training has 12,852 pairs in 157 families; validation has
2,635 in 41; the reserved test has 3,840 in 50. Only training and validation
effects were exported. Independent calculations verified every split and
eligibility decision and all 15,487 exported targets. Sixteen focused tests pass;
source replacement and resealed-receipt controls fail as intended.

The first model-selection round compared longer sequence words, a nonlinear
sequence-background model, and a combination of short and long words. Every fit
is recorded under `phase2-kmer-six/`, `phase2-context/` and `phase2-multiscale/`
in the experiment directory. The globally selected candidate has validation
MSE 0.0153321, compared with 0.0153782 for the best required sequence baseline.
Its signed ten-edit ranking utility is -0.0085593, below the random comparator's
0.0011730 and the required zero-utility comparator. It fails the development
requirements. These are validation results, not confirmation on reserved data.

The second round used the frozen protocol's permitted training-only
sequence-to-activity supervision. Of 33,679 design-admitted training reference
sequences, 29,981 in 13,555 families passed measurement quality checks. An
independent review verified every label and exclusion and reproduced the intake
byte for byte. A joint model fitted the paired effects and reference activity,
with an activity-only ablation and a same-supervision dinucleotide control.
The 80 declared fits are preserved in `phase2-auxiliary/`.

The best new candidate's validation MSE is 0.0157137 and ranking utility is
0.0701547. Better ranking does not offset its worse MSE. Global selection therefore
retains the earlier context model, which still fails ranking. The reserved test
remains unevaluated. Neither round will receive a further parameter-grid
extension. A subsequent model requires a specific new hypothesis and a bounded
plan; changing the required endpoints or choosing a runner-up by its ranking
score would not resolve this failure.

The third round added 532 admitted JASPAR2024 binding-profile scores to k4
features, using the same training-only reference-activity supervision. Its
12 declared candidate fits and one fixed, ineligible column-permutation
diagnostic completed. The best new candidate had validation MSE 0.0158121511
and signed ranking utility 0.0272682122. The earlier context model remains the
global MSE winner and still fails the ranking requirement. This brought the total
to 176 logged fits, including the diagnostic. Independent review reproduced all 13
prediction arrays, metrics and selected IDs; the reserved test stayed closed.
These results do not support a new predictive claim.

The completed diagnostic round found culture-pair edit-effect correlations of
0.484-0.581 on the exposed validation set. Larger measured effects are present:
143 of 150 rows with absolute mean effect at least 0.25 agree in direction across
all three cultures. The selected model's predictions have much less spread than
the observations, and its same ten chosen edits have negative signed utility in
one of the three cultures. These observations do not identify assay noise as the
limiting cause. Positive rescaling cannot repair the ranking because it leaves
both selection order and predicted signs unchanged.

Independent calculations reproduced the culture summaries, sequence diagnostics
and frozen selections. A prespecified training-only decomposition found that
62.6-65.6% of centered culture disagreement remains within parents when examining
parents with multiple edits. Removing parent means increases culture correlation
by only 0.0094-0.0108 in that population. This weakens the proposed explanation
that common parent offsets alone account for the disagreement; it does not
identify reference-measurement error.

The source audit verifies genomic-to-FASTA coordinates and all 24,942 labelled
single-base pairs. It does not establish the complete per-record link between
deposited sequences, synthesized oligo orientation and the reporter promoter.
The fourth experiment therefore used FASTA-relative strand and position features,
with that physical-orientation uncertainty retained. Two nested architectures
and four fixed penalties defined eight candidate fits, plus one ineligible
motif-structure control. Its pooled-feature ablation tested whether retaining
geometry helps. Eight implementation tests pass, and independent calculations
verify the scoring, channel scales, split isolation and candidate accounting.
The full label-free probe processed 15,896 unique sequences and 532 profiles
for both real and permuted motifs in 111 seconds, with 1.00 GiB measured peak
RSS on this Mac. Four synthetic solver calls measured engineering cost without
assay outcome exports. These checks passed independent review before fitting.
The original selection and both acceptance requirements remain unchanged; no
reserved outcomes were opened for these checks.

The fourth round completed once. Its best new candidate, the geometry model at
penalty 1, had validation MSE 0.0158400129 and signed-ranking utility 0.0287999876.
The earlier global context winner remains selected because its MSE is lower;
it still fails ranking. The conditional permuted-profile control had MSE
0.0156810798 and utility -0.0237892040 and remains ineligible. Independent replay
reconstructed all nine prediction vectors, weighted normal-equation residuals,
metrics, selected IDs and receipt identities. It reused the frozen low-level
motif scorer, which had separate independent likelihood and deliberate-defect
checks. No predictor qualifies from this round.

The history through the first four rounds has 185 logged fits: 151 eligible candidates and 34
baselines or controls. These four grids are closed. Further work is limited to
verifying distinct information/model hypotheses and their actual source,
context and access requirements before proposing another experiment. The
reserved test remains unread and the acceptance rules remain unchanged.

The next source candidate is [MPRALegNet](https://doi.org/10.1038/s41586-024-08430-9),
which adds external reporter measurements and learned nonlinear sequence
features. Its published mutation-effect examples are in HepG2 and K562; its
WTC11 model was trained in stem cells, not the selected day-14 neurons. The
official code and checkpoint listings are available, but the candidate is not
yet admitted to the benchmark. First establish a named checkpoint's supervised
training exposure, the construction of its 230 input bases and the normalization
of its training target. Full 270-base inference is computationally plausible;
its biological transfer remains an experiment. Do not crop the assayed DNA or
invent flanks to match the training length. This source review introduces no
new fit and does not change the evaluation population or acceptance rules.

Reviewed results, diagnostics, code and input identities are preserved locally
under `build/neural-predictor-reviewed-20261005/`. Its `snapshot.json` lists the
copied files and four omitted, regenerable numerical caches. This is an evidence
snapshot; the scripts retain their original research-workspace paths. The final
fourth-round report and all 62 manifest-pinned artifacts are in its
`geometry-final/` subdirectory; the earlier pre-fit snapshot is preserved
separately. The subsequent source audit, independent feasibility review and
bounded next-audit plan are preserved in `new-information/`; its manifest pins
47 files. Neither reviewed source has been admitted for another fit.

The subsequent bounded source projection is complete and preserved in
`input-contract/` with a 79-file manifest. It verifies 92,370 released WTC11
design rows, including 46,185 retained inputs assembled as 15 adapter bases,
200 core bases and 15 adapter bases. This resolves the released table's input
construction. Its relationship to the checkpoint's legacy training table,
target generation, complete supervised-family exposure and transfer to the
neuronal assay remain unresolved. No weights, predictions or new fits were
produced. The outcome-free Kosicki family inventory and independent graph
check are preserved in `exposure-inventory/`; they do not yet establish overlap
with the pretrained model's supervised data.

The subsequent exact-containment check in `exact-source-overlap/` found 807
Kosicki design records containing a full 200-base core from the released WTC11
table, including three families assigned to validation and three to the reserved
test before measurement quality filtering. Independent review verified all
positive matches and whole-family propagation. These are conservative exposure
flags, not proof of exact legacy-checkpoint membership. Negative containment
does not clear partial overlaps or other supervised sources. The existing
splits, evaluation population and 185-fit history remain unchanged.

A subsequent conservative coordinate screen used the public GRCh38 library
reference ENCFF508RRB. Independent replay verified 5,024 overlapping parents
and 3,268 flagged families, including seven validation and nine test-assigned
families before measurement filtering. All 807 exact sequence matches also map
to exact intervals in this reference. Missing element IDs still prevent a
complete checkpoint-membership map. This pretrained checkpoint is not admitted
to the current frozen benchmark; that source branch is closed for now.

The fifth round uses the original yeast-derived LegNet configuration at commit
`15b57b18855c449e5b8b467df6b84409ea7bda2c`, adapted to four DNA channels, full
270-base sequences and an unbounded scalar activity score. Its 2,124,913
parameters start from fresh weights. Two fixed candidates test paired-effect
supervision with and without the admitted training-only activity loss. Each
receives 1,000 Lion updates with an effective batch of 128 examples per task;
validation is scored once at the final update. The population, selection rule
and acceptance criteria remain unchanged.

The implementation passed 11 focused tests, four deliberate-defect controls,
independent review, CPU/Metal numerical comparison and exact synthetic GPU
checkpoint recovery. All 66 code, source, input and review artifacts are frozen
under `build/neural-predictor-reviewed-20261005/neural-network-round/freeze.json`
(SHA-256 `32e1c52769855683f84bb092a553ebd281c2d199231c1de2d9afca72f0c38830`).
Both biological fits completed in `neural-network-round/run-1/`. The joint model
has validation MSE 0.0158529157 and ranking utility 0.0095107575; the paired-only
model has MSE 0.0164350293 and utility 0.0035944047. The required k-mer baseline's
MSE is 0.0153782158. The earlier global winner remains selected. Independent
metric replay and bit-identical prediction replay from both saved models passed.
The ledger now contains 187 fits (153 eligible candidates, 34 controls/baselines).
The post hoc training-versus-validation diagnostic quantifies the generalization
gap without refitting or changing selection. The complete report is
`neural-network-round/RESULTS.md`; these two planned fits are closed.

The sixth round uses sequence-only hg38 pretraining from the pinned HyenaDNA
tiny checkpoint. Its frozen encoder supplies 128 strand-averaged feature
differences, added to the previous 462-feature base. Four regularized linear
predictors and one conditional random-encoder control were fitted. The best new
candidate's MSE is 0.0155413555 and utility is 0.0586510080. The random control's
MSE is 0.0153064118 and utility is -0.0163925187. The control is ineligible and
the historical global candidate remains selected. Eight focused tests, four
deliberate defects, independent metric verification and full encoder replay
passed. The ledger now contains 192 fits (157 eligible, 35 controls/baselines).
The report is `sequence-prior-round/RESULTS.md`; all five planned fits are closed.
The documented pretraining objective has no reporter outcomes, but genomic
sequence overlap is possible and disclosed. This is an assay-family holdout,
not a claim that pretraining never encountered the DNA.

The seventh round tests exact strand-averaged gapped counts, length 10 with six
retained positions, against an ineligible exact-ten-base-word control. Each
appends a block to the same base and uses training-only energy scaling and the
inherited penalty 0.1. The candidate's MSE is 0.0153543813 with utility 0.0455434550;
the control's MSE is 0.0153063780 with utility 0.0251110163. Both beat individual
required comparators. The control additionally beats the historical eligible
winner's MSE, but remains excluded under this round's prospective plan.

Eight focused tests, independent feature-kernel oracles, deliberate-defect tests
and independent saved-coefficient/metric/selection replay passed. A confirmed
snapshot-read integrity defect was fixed and tested before fitting. Both planned
fits are closed at 194 total (158 eligible, 36 controls/baselines). The full report
is `gapped-round/RESULTS.md`. `gapped-round/FOLLOWUP_DRAFT.md` specifies a proposed
validation-driven nomination of the exact existing control for confirmation.
The user subsequently preferred retaining the original eligibility, so that
draft was not adopted. The decision is recorded in `control-role-decision.json`
alongside the archived rounds. `baseline-readiness/review.json` verifies the
unchanged sequence baselines using training and validation only. No reserved
outcomes have been accessed.

The subsequent training-only diagnostic checks whether lower reference activity
is associated with smaller cross-culture edit-effect products. All 326 parents
join exactly to admitted auxiliary training measurements. Baseline activity from
one culture is correlated with the product of effects from the other two, using
the original equal-family pair weights and all three rotations. The correlations
are 0.042899, 0.014171 and 0.056212; their fixed mean is 0.037761. The mean's
2,000-draw family-bootstrap interval is [-0.016558, 0.096671], and every individual
interval also spans zero. All leave-one-family-out correlations remain positive.
The narrow directional hypothesis remains unresolved and does not justify a new
fitting round by itself. Five focused tests and three deliberate-fault checks
passed before the calculation. `context-response/RESULTS.md` records the result;
no new fits, validation calculations or reserved-outcome reads were added.

The eighth round compares two prospectively eligible shallow boosted-tree
predictors. Both use canonical k4 edit-count differences; one additionally uses
the symmetric midpoint of full-insert k3 counts. Fixed sign augmentation and
odd projection preserve zero, swap and reverse-complement behavior. Each model
received 100 iterations at depth 3, using training-only binning and mean-one
equal-family weights. No internal validation or post-result tuning was used.
The original exact-word control remains ineligible.

The delta-only model has validation MSE 0.0155772650 and signed ranking utility
-0.0468737604. The context model has MSE 0.0155400949 and utility 0.0451837081.
Both have 2,027 fitted scalar parameters under the prospective counting rule.
The context result passes the individual ranking threshold but fails the
required k-mer MSE comparison. Global minimum-MSE selection retains the historical
model. Both planned fits are closed at 196 total (160 eligible, 36 controls).
Nine focused tests, independent native-tree checks and deliberate defects passed
before fitting. Independent feature construction, scalar JSON inference and
native-model replay reproduce all 30,974 saved training/validation predictions
bit for bit; all new metrics, selected IDs and accounting agree. The report is
`paired-tree-round/RESULTS.md`. Reserved outcomes remain unopened.

The subsequent feature-block diagnostic reconstructs four saved joint ridge
models using training data only. The exact-word control places 99.69% of squared
kernel energy on the same example or within its family; the random Hyena control
places 98.13% between families. The gapped-six and pretrained Hyena blocks place
4.24% and 88.24% between families, respectively. These are matrix-entry energy
fractions, not accuracy fractions. Shared-column counts alone do not explain
this distinction. All four joint models improve training MSE through their
added blocks and interactions; changing their base coefficients alone increases
training error. These decompositions do not establish unseen-family performance
or justify promoting a control.

Six focused tests, deliberate-defect checks and independent row-level arithmetic
replay pass. Both Hyena feature/kernel decompositions were independently rebuilt
from saved caches; the result review did not independently rebuild the full
gapped kernels or base normal equations. The report and exact review limits are
in `block-diagnostic/RESULTS.md`. The ledger remains 196 fits, original eligibility
is unchanged, and no validation or reserved-outcome calculations were added.

The MPAC source audit is separate from those closed rounds. The author export
defines 200-base inserts with authentic reporter-vector flanks and a documented
sliding-window route for longer input. Applying it to all 270 assayed bases would
be a fragment-aggregate transfer experiment, with no invented genomic flanks.
Primary methods establish that the ten original design libraries cover the
Malinois training, validation and test pool. Their conservative union contains
929,712 source-library records, including unused oligos and repeated sequences.
All observed insert lengths and both orientations were compared with all 81,952
benchmark design sequences. The screen finds 1,017 matching records in 902
families: 13 training, three validation, five test and 881 outside the benchmark's
single-base set. These assignments precede measurement QC. Positives establish
possible source exposure, not exact retained checkpoint membership.

Independent full-string matching verifies every positive and whole-family
propagation. The complete report is `mpac-source/full-design-screen/RESULTS.md`.
Exact membership, original tuning exposure and partial genomic overlaps remain
unresolved; the checkpoint is not admitted. Author mapping code also distinguishes
native older-assembly variant IDs from hg38 analysis coordinates, so IDs must not
be silently treated as current insert intervals. No MPAC inference or new fit has
been performed. A subsequent bounded archive-prefix request retrieved opaque
compressed checkpoint bytes without interpreting or retaining tensor values.
Author source and input contract evidence are preserved in `mpac-source/`.

The source-only preprocessing notebook audit identifies why the published
performance table cannot establish original training retention. Training filters
DNA counts at 20 before aggregation; performance processing uses zero. Both final
schemas omit the count fields, and neither notebook exports a retained-row
manifest. A performance-table match therefore cannot prove passing the training
filter. The eight flagged families all match ordinary central variant designs
outside the original chromosomes 7/13 excluded from both fitting and tuning.
Original retained-row evidence is still required to resolve those exposure flags;
a larger performance-table download alone would not answer the question. The
archived processing contract and checks are in
`mpac-source/retention-contract/REPORT.md`; the full source audit is in
`mpac-source/admission-audit/FOLLOWUP.md`.

The original-release code confirms that training configuration is bundled with
model state in `artifacts/torch_checkpoint.pt`. Historical example commands use
different schemas and splits, so they cannot establish the final configuration.
A 65,536-byte archive-prefix probe stopped at an earlier checkpoint member and
recovered no configuration. The source trace and bounded result are preserved in
`mpac-source/artifact-trace/REPORT.md` and `mpac-source/prefix-probe/REPORT.md`.
This closes that bounded retrieval attempt; it does not resolve membership.

The next bounded representation is a four-coordinate DNA-shape projection from
the pinned DNAshapeR physical lookup. All 512 ordinary pentamer classes expand to
1,024 oriented words. Direct summation equals a fixed projection of pentamer
counts; this supplies coefficient constraints, not additional sequence
information. MGW and ProT use base-pair descriptors; Roll and HelT average the
two central steps of each pentamer. That pooled descriptor differs from a native
step-track average at the sequence boundaries. Every complete observed window
is used, including terminal edits, with no invented flanks.

The outcome-free check passes 2,430 synthetic 270-base substitutions, with exact
allele-swap antisymmetry and reverse-complement agreement within numerical
tolerance. Independent parsing reproduces all 6,144 field values exactly, and
independent Decimal calculations verify 2,478 substitutions including minimum
lengths. Wrong field indices, strand handling, omitted terminal windows and
source corruption are detected. This verifies the proposed representation;
the original R/C++ package was not executed. No assay data were read for these
checks. The frozen diagnostic compared this fixed projection with one tuple-permuted
control using 12,852 training rows in 157 families. Physical-minus-control
alignment differences were -0.00277396, -0.0000491640 and +0.00297658 across
the three culture pairs; their mean was +0.0000511509. The 10,000-draw
family-resampling descriptive interval was [-0.00680712, +0.00677215], and
37 of 157 family-deletion means were nonpositive. No resample or deletion was
undefined. All three progression screens failed. Independent reconstruction
reproduced the features, all resamples, deletions and decision. The pooled
hypothesis is closed without fitting or reseeding; this does not reject all
DNA-shape representations. The fitting ledger remains 196, and neither
validation nor reserved outcomes were read. See `dna-shape-source/CONTRACT.md`,
`dna-shape-source/RESULTS.md` and `dna-shape-alignment/RESULTS.md`.

The ninth round combined the existing context and gapped-six predictors using
five connected-family folds. Each fold rebuilt its basis and scales from the
remaining training families. Ten component fits supplied 25,704 excluded-family
predictions, from which one convex mixture weight was fitted. It chose the gapped
endpoint exactly. Independent replay reconstructed every fold prediction and
all three validation vectors, with matching metrics and selected IDs. This
round adds no new model behavior. Its result and checks are preserved in
`family-stack/RESULTS.md` and `family-stack/result-review/`.

The eligible gapped-six model has validation MSE 0.0153543813 and signed utility
0.0455434550. The respective thresholds are 0.0153782158 and 0.0011729805. The
historical context model has MSE 0.0153321180 and signed utility -0.0085593209.
These results explain the selection conflict: the original procedure selects
the lowest MSE before checking ranking. The reviewed proposal in
`joint-selection-proposal/AMENDMENT_DRAFT.md` instead requires both validation
gates before choosing minimum MSE. The user's contextual endorsement of the
single-model path adopted this exact proposal before confirmation. The approval
receipt preserves the actual message and interpretation. All 207 fits, the
original failed selection and excluded controls remain in the record.

Confirmation preparation has saved sequence-only predictions for all 3,840
reserved edits in 50 families, using the unchanged candidate and baseline
weights. Independent replay verifies every prediction and all 2,000 family
bootstrap samples. The new inference path also reproduces the 2,635 historical
validation predictions exactly. Statistical tests cover all seven required
comparisons, repeated family occurrences and strict acceptance. These steps add
no fits and used no reserved outcomes. Preparation and review evidence are in
`confirmation/`; their historical pending-state records remain unchanged.

The final confirmation used freeze
`2c308a08ab203f9814625eb9d6737d2bf15975e14e8e0ea4fe88c29a56edc141`.
Its one-use runner completed successfully, scored every eligible reserved pair
and preserved the scientific failure. The mean endpoint results are:

| Policy | Family-weighted MSE | Ten-edit signed utility |
|---|---:|---:|
| Gapped-six candidate | 0.0169575704 | -0.0379617141 |
| Dinucleotide baseline | 0.0170613548 | -0.0376025249 |
| k-mer baseline | 0.0172303711 | -0.0890435292 |
| Zero effect | 0.0169190102 | 0 |
| Random selection | — | +0.0225380178 |

All seven paired, whole-family 95% bootstrap intervals include zero, so no
required improvement comparison passes. The candidate's fixed ranking has
negative utility in each of the three cultures. Slightly better point MSE than
the sequence baselines does not resolve the failed zero comparison or ranking.
This result applies to the declared WTC11-NGN2 reporter assay and remains
conditional on the source's pooled preprocessing. It does not identify the
reason the validation result failed to generalize. See `confirmation/RESULTS.md`,
`confirmation/confirmation-1/result.json` and the frozen prediction/target files
for all intervals, selected edits and exact identities.

Independent replay expanded all 2,000 bootstrap draws and reproduced all seven
intervals, point metrics, culture results and selected IDs, with maximum numerical
disagreement 2.78e-17. It used the exported targets and frozen predictions without
reopening the raw measurement source or rerunning the confirmation. The review
is preserved in `confirmation/result-review/`.

The subsequent NT-v2 source audit is closed without checkpoint admission:
checkpoint-specific supervised selection exposure and overlapping-family
mapping remain unresolved. The culture-noise audit also adds no fit. Three
paired culture measurements alone do not identify reference-only measurement
error; abundance-aware alternatives require their underlying count data.
These audits are preserved in `nt-source-audit/` and
`measurement-noise-source/`.

The tenth round, `compact-transfer/`, uses a 4,993-parameter additive motif score
and a 26,505-parameter extension with local factorized interactions. Each score
uses the complete 270-base insert; effect is alternate minus reference. Joint
and staged training receive the same ordered 128,000 pair examples and 128,000
activity examples within each of three seeds. All twelve fixed final checkpoints
were scored once, with no grid extension. The original confirmation was not an
input. The design borrows from
[ExplaiNN](https://doi.org/10.1186/s13059-023-02985-y),
[BPNet](https://doi.org/10.1038/s41588-021-00782-6) and
[human MPRALegNet](https://doi.org/10.1038/s41586-024-08430-9), without importing
their pretrained weights or claiming architectural novelty.

Validation MSE ranges from 0.0158717786 to 0.0159574148, versus the required
0.0153782158. Training gains over zero are only 0.128–0.863%. No candidate meets
both gates, and none of the four matched contrasts improves both endpoints in
all three seeds. Full independent Metal replay exactly reproduces all 185,844
training and validation predictions and verifies 2,649 frozen identities.

CPU replay of the first three models exposes a separate ranking limitation:
nine tile selections change despite prediction differences below 1.94e-7.
The changed representatives have identical local model dependencies in both
orientations and are mathematically tied. Float32 subtraction of whole-sequence
scores breaks those ties differently by backend. Their different observed
effects change ranking utility. The original Metal metrics remain preserved;
portable ranking is unresolved, and no favorable CPU score replaces a result.

A fixed, training-only batch diagnostic finds activity gradients 2,813–3,230
times larger than paired-effect gradients at initialization, despite normalized
losses. Final joint-model ratios are 44–214 on that batch. Independent explicit
loss derivatives reproduce the gradient norms. This is a measured optimization
mismatch, not proof of a sole cause: AdamW rescales updates, staged training
removes the activity loss in its final phase, and trainable weights can compensate
for pooled local effects. The existing eight-base learning test does not verify
full-length sparse-edit trainability. Results, limitations, all twelve fits and
the next discriminating check are in `compact-transfer/RESULTS.md`.

The subsequent `trainability-control/` experiment uses an analytically
representable GC-based score on full 270-base sequences. Exact allele effects
include terminal attenuation and neutral substitutions. Joint, static-balanced
and pair-only training each ran with seeds 42–44 under the same 2,000-update
recipe. All nine runs pass all five recovery gates. Each predicts the direction
of all 2,160 nonzero held-out edits correctly, including all 16 terminal edits;
relative MSE ranges from 0.0000568 to 0.0002288 of the zero-effect MSE.
Independent Metal replay exactly reproduces all 174,960 effect/activity values.

Ordinary joint training also succeeds. This establishes basic full-length
trainability for this deliberately simple teacher, not the cause of the
biological failure. Its initial activity/pair gradient ratios are 4.59–72.56,
well below the inspected biological ratios. These nine synthetic optimizer runs
are recorded separately; the biological fitting ledger remains 219. Full
results and interpretation limits are in `trainability-control/RESULTS.md`.

The separate `stable-effects/` inference repair computes local allele differences
before pooling and shares one CPU float64 evaluation across identical contexts.
It changes no learned weights. All nine recorded unstable tie cases now receive
exactly equal effects within their request. On 7,911 fixed-checkpoint comparisons,
the maximum difference from full CPU float64 scoring is 2.25e-16. Independent
90-digit scalar-oracle tests detect fourteen deliberate faults. This repairs
the arithmetic and structural-tie handling; original experiment decisions are
preserved, and no biological accuracy gain or arbitrary cross-device ranking
guarantee is inferred. See `stable-effects/RESULTS.md`.

A primary-source audit identifies Rummel 2023 human induced-neuron MPRA as a
possible new confirmation source, with Lee 2025 and Matoba 2026 progenitor-cell
studies as alternatives. They have different insert lengths and assay contexts.
Exact design/overlap, normalization and population admission remain outstanding;
no per-variant outcome table was downloaded or evaluated. The source record is
`trainability-control/FRESH_CONFIRMATION_SOURCES.md`.

The eleventh round, `biological-optimization/`, fixes the repaired inference
path and compares three reused additive-joint checkpoints with six new fits:
static-balanced and paired-only objectives across seeds 42–44. Initialization,
the complete 192-sequence forward batch, ordered family-balanced examples,
paired coefficient and 2,000-update recipe remain matched. The balanced weights
come from initial training-only gradients and are fixed before fitting.

Every new fit lowers training MSE by 2.92–3.70% against its joint control, but
worsens both validation endpoints. None meets either development gate. Balanced
training slightly improves validation MSE against paired-only in every seed;
both methods select identical edit IDs and directions within each seed, so
their ranking utility is identical. This identifies an optimization effect
without a transferable prediction gain under the fixed protocol. It does not
establish a single cause for the remaining biological failure. The round adds
six biological fits, for 225 cumulatively; synthetic accounting remains nine.
Independent explicit full-sequence replay verifies 139,383 predictions within
1.50e-15 and agrees on all metrics, decisions and diagnostic-vector checks.
See `biological-optimization/RESULTS.md` for the table and independent review.

A subsequent training-only bound checks exact local-context collisions,
including allele-swap signs and strand symmetries. Its minimum error over free
context-group predictions is 0.00108891, or 6.408% of zero-prediction error.
The 1,082 nonsingleton groups contain 2,277 training rows. These exact collisions
alone leave substantial fitting headroom; the bound relaxes additional shared
motif constraints and does not identify an optimizer failure or a biological
noise ceiling. Independent alternate keys and weighted-variance calculations
reproduce the result and catch four faults. No fits or validation calculations
were added. See `alias-bound/RESULTS.md`.

The deeper source audit locates publisher descriptions of Rummel's separate
design table, but cannot verify its downloadable bytes/schema. Lee's public
code expects an external library-design file that was not found in its official
repositories or listed GEO files. Neither source is admitted. Exact sequences,
assay handling and whole-family overlap checks remain required; the records are
in `fresh-confirmation/`.

A newer neuronal source now has verified design and measurement intake. The July 2026 Salomon
preprint supplies a public author design and pair map for 46,374 primary
single-base edits across 18,572 reference constructs. Every pair passes complete
270-base sequence, orientation and coordinate-consistency checks. Independent
pairing and graph implementations reproduce the results. After prior-library
and known published-example exclusions, 45,041 design pairs in 17,546 connected
components remain before extended locus exclusions and measurement QC.

The prospective intake additionally excludes chromosome-17 components, covering
the discussed MAPT/KANSL1 region without guessing inversion boundaries. The
design split was frozen before the complete response download. The pinned author
coefficient table yields 21,781 training edits in 8,774 components, 7,390
validation edits in 3,029 components, and 7,372 reserved confirmation edits in
2,938 components. Independent parsing reproduces every export byte and verifies
all 301 frozen inputs. Confirmation values remain separate from development.
The target is the author's signed alternate-minus-reference log2 reporter
coefficient, with no significance or activity filter. Missing rows are recorded.
The context is `salomon2026-ngn2-author-a686f1d`; unresolved timing metadata is
preserved instead of assigning an unsupported collection day. The benchmark
conditions on the published processed coefficients.

This source can test a distinct explanation for the current transfer failure:
the breadth of directly measured edit responses. The proposed comparison holds
the learner and compute fixed while varying the breadth of direct-edit
supervision within the new study. The frozen protocol prescribes six fresh
LegNet-derived fits: 157 training components versus all eligible training
components, paired at seeds 42, 43 and 44. Sixteen sequence baselines are fitted
in the same target units. Both prediction error and ten-edit selection utility
must improve before any candidate can enter the single new confirmation.
The design and source contracts are in `salomon-source/`; the prospective
protocol, frozen split and measurement intake are in `salomon-intake/`.
The biological ledger contained 225 fits before this comparison. The completed
comparison used training freeze
`a356418d373435680f7134623a0ecbf54e5e03a77f0fa8fff25c0348a7c8b54d`.
All sixteen baselines are complete. Independent weighted-design SVD replay
reproduces their coefficients and predictions within 7.3e-16, with identical
metrics, rankings and baseline selections. The required validation MSE is
0.0286251098, and the highest comparator ranking utility is 0.0351089021.
The six neural fits completed the fixed 1,000-update budget and final-checkpoint
selection. Small-seed42 has validation MSE 0.0338784152
and ranking utility 0.1291521430; broad-seed42 has 0.0295765290 and
0.0468130343; small-seed43 has 0.0338055949 and 0.0497257430; broad-seed43 has
0.0293235670 and 0.0985900201; small-seed44 has 0.0316989500 and
0.0424524222; broad-seed44 has 0.0295586035 and 0.0689523746. All six exceed
the comparator ranking utility and miss the MSE requirement. The completed
runner nominates no candidate. Independent Metal replay exactly reproduces all
110,888 saved neural predictions, with matching metrics, rankings and nomination.
Broader
training improves MSE in all three seeds and utility in two, so the declared
requirement for improvement in both endpoints across all seeds fails. The
biological ledger now contains 247 fits. No Salomon confirmation outcome has
been opened by development.

The new scalar confirmation statistics and atomic one-use core pass independent
synthetic checks, including exact rational bootstrap oracles, each of the seven
strict gates, repeated-component ranking, concurrency and retained failures.
Prediction preparation and final command integration are implemented and pass
37 independent synthetic tests, with 35 deliberate faults detected across the
two reviews. The root reran both independent suites and all 13 local tests.
Review exposed and repaired extra manifest reads, redirected input paths, and
missing bindings between the consumed plan/results snapshots and their original
experiment identities. Preparation now admits an exact input catalogue and
requires the completed six-fit comparison plus matching independent replay.
The final command claims the fixed confirmation attempt before reading targets
and rechecks all inputs before publishing its result.
The reviewed implementation and evidence are frozen under
`667c1dd89e92aa90a0b5f1ca664389ae16d0240adf396d7e81e208779f926d4e`.
Actual prediction preparation and the one permitted confirmation remain
conditional on an eligible, independently replayed candidate. These checks are
in `salomon-confirmation/`; they add no biological fit or outcome evaluation.

A saved-prediction diagnostic is complete for the comparison.
It compares both models on their shared 157 training components, summarizes the
broad model's additional components separately, and decomposes validation MSE
into target magnitude, prediction magnitude and alignment. A fixed-list
cross-evaluation separates selection membership from predicted direction.
Six exact-arithmetic tests and seventeen synthetic integration tests pass;
thirteen deliberately broken arithmetic implementations are detected. The
reviewed diagnostic and evidence are frozen under
`4fafe39232114f4a56f55e03c0b5b5c0509b48a328ed2c44be6ed040056ab7af`.
Independent replay verifies all seventeen result cells and three fixed-list
tables. Small models remove 99.4–99.7% of the shared training population's
zero-prediction error; broad models remove 7.6–11.7%. At the same sample budget,
expected draws per component fall from 815.3 to 14.6. Lower prediction magnitude
accounts algebraically for 93.0–96.4% of the broad models' validation MSE
reductions. The paired ten-edit lists share no components. These observations
separate memorization, exposure and selection behavior without identifying one
cause of the biological error.

A subsequent exact-arithmetic bound tests whether positive global rescaling
alone could repair the validation error. For prediction p and target y, its
best possible MSE is E[y²] - max(E[py], 0)²/E[p²]. All six lower bounds remain
above the k5 comparator: 0.0286651285–0.0287465044 versus 0.0286251098.
Independent calculation from the saved arrays agrees. This rules out positive
rescaling alone; no calibration was fitted, and nomination and acceptance
requirements are unchanged.

An architecture review identifies a specific representation limit in the current
CNN. Its local receptive field spans 79 bases; global channel averages retain
some context but lose certain distant motif spacings. A synthetic pair with
motif separation of 100 versus 105 bases has identical neighborhood histograms
at all thirteen convolution depths, for both alleles and reverse complements.
Independent replay and integer neighborhood checks reproduce this witness;
boundary and nearby-interaction counterexamples break the equivalence as expected.
This establishes a constructed representation limit, not the cause of biological
error. A relative-position attention scorer is one candidate remedy; a larger
or dilated convolution can also retain the missing information.

The source review isolates Interlace's rational gate as a separate hypothesis.
It belongs inside a shared scalar sequence scorer: a proposed direct pair gate
fails exact three-sequence consistency despite preserving zero edits and swap
signs. The ordinary-gate comparator must retain the same attention, projections
and depth. Monodratic's default remote-block budget covers all eligible blocks
at 270 bases, so it supplies no sparse-routing comparison in this setting.
Evidence and the scalar-consistency counterexample are in
`build/neural-predictor-reviewed-20261005/architecture-diagnosis/`.
A four-fit synthetic experiment now tests a compact, 11,368-parameter
bidirectional attention scorer against the current CNN. Both learn the same
local control and a distance-dependent interaction teacher, using 26,208 edits
in disjoint layout groups. Independent rescanning verifies every label and
excludes shared reference or edited sequences across partitions, including
reverse complements. Exact neighborhood equivalence gives the CNN a validation
MSE lower bound of 10607/74420 on the interaction teacher; the new scorer must
beat that bound, remove at least 90% of zero-prediction error, and achieve at
least 95% sign accuracy on nonzero effects. Both models must pass the local
control. The four final checkpoints at 512 updates are fixed in advance.

The attention oracle, gradients and full 270-base input pass CPU checks and
CPU/Metal qualification. Independent runner review and fifteen deliberate
faults across the model, data and runner checks pass. Source, data and criteria
are frozen under
`db47d342af838ab6164fc38ab2814bce92a073033f621188f038e45b9dacd33d`.
All four fits completed once. Independent Metal replay exactly reproduces all
23,448 prediction and property values and verifies every teacher label, training
scale, checkpoint, metric and decision. Results are in `attention-capability/`:

| Model | Local validation MSE | Interaction validation MSE | Nonzero interaction sign accuracy |
|---|---:|---:|---:|
| Current CNN | 0.0000541633 | 0.1435870595 | 72.41% |
| Compact attention | 0.0004984695 | 0.1455605290 | 72.41% |
| Zero prediction | 0.5 | 0.1860655738 | 0% |

Both local controls pass. Attention misses all three interaction requirements;
the CNN respects its mathematical floor of 0.1425288901. This architecture and
training recipe therefore do not pass the capability test. These four
synthetic fits add no biological fits; that ledger remains 247. Model sizes
differ, so this experiment does not establish an efficiency comparison.

A training-only diagnostic was specified before the interaction outcomes were
opened. It compares fixed spacing and translation contrasts, isolates signed
and absolute position coefficients, checks their derivatives, and reconstructs
aligned per-token contributions. A source-only optimizer calculation also shows
that the fixed schedule limits each zero-initialized position coefficient to
about 0.12846 in real arithmetic. Dividing distance by 270 therefore limits the
direct absolute-distance logit change over 20 bases to about 0.00952. This is a
bound on one signal path, not on the full nonlinear model. The completed
diagnostic finds spacing contrasts below 0.000147 for a target of -2, with
nonzero derivatives and all four absolute-distance coefficients at 97.76–100%
of the schedule bound. An independent functional attention graph reproduces
all 16 fixed training contrasts, coefficient ablations, aligned token values
and eight derivative checks. These measurements motivate a conditioning test;
they do not identify a unique cause.

The frozen follow-up removes only the division by sequence length. It preserves
all initial tensors, the 512-update recipe, data and five acceptance checks,
and reuses the completed CNN controls. Both new checkpoints independently
reproduce all 11,724 saved prediction values exactly. Results in
`attention-position-units/` are:

| Attention distance units | Local MSE | Interaction MSE | Interaction sign accuracy |
|---|---:|---:|---:|
| Sequence-normalized, previous trial | 0.0004984695 | 0.1455605290 | 72.41% |
| Bases, follow-up | 0.0006014430 | 0.1449808052 | 72.41% |

The local control passes; all three interaction requirements still fail.
The 0.4% interaction-error reduction is insufficient. Changing these units
changes optimization conditioning, not the function family or parameter count,
and is not a new architecture. This second use of the synthetic development
split adds two synthetic fits and no biological fits. The same fixed training
diagnostic now finds maximum absolute interaction-spacing response 0.00488833
for a target of -2, and maximum translation response 0.0150248 for a target of
zero. Two of four spacing responses have the wrong direction. Larger position
sensitivity alone has not produced the required interaction response.

The diagnostic's original finite-difference step failed its existing float64
gradient tolerance after the unit change. That failed attempt is preserved.
A separately reviewed repeat derives the perturbation from the maximum
positional-feature magnitude, keeps derivatives in base-coefficient units and
retains both step sizes and tolerances. All eight gradient checks then pass.
Independent explicit-softmax replay agrees on all 16 records and four
ablations, with scalar discrepancies below 6e-15 and gradient discrepancies
below 1.5e-14. See `attention-position-units/mechanism-diagnostic-scaled-step/`.
This numerical repair changes verification resolution, not fitted predictions.
An unexecuted diagnostic combines the two admitted motif disruptions on the
same twelve training layouts. A double-mutant contrast separates partner
dependence from additive positional effects; a complementary contrast tests
whether a useful interaction is canceled by a partner-independent response.
This has not been executed and introduces no new acceptance rule. Its reviewed
definition is in `attention-position-units/NEXT_DIAGNOSTIC.md`.
The subsequent work prioritizes model changes and training directly. The exact
Interlace rational feed-forward gate has now been ported, preserving all initial
trainable tensors and the scalar edit interface. Its scaled implementation
passes independent value and derivative checks, including extreme inputs, and
CPU/Metal qualification. Two fixed 512-update fits completed, with all 11,724
prediction values independently reproduced. Local MSE is 0.0007130195;
interaction MSE is 0.1454870948 with 72.41% sign accuracy. This gate-only change
does not improve the interaction result. See `attention-interlace-gate/`.

A matched follow-up gives both ordinary and Interlace scorers 4,096 updates,
eight times the previous budget, with the corresponding longer OneCycle
schedule. Initial tensors, sampling, data, batches and thresholds remain fixed.
Both arms have completed and passed all three interaction requirements. Ordinary
attention reaches MSE 0.0011160242; Interlace reaches 0.0005454337, reducing error
from zero prediction by 99.40% and 99.71%, respectively. Every nonzero effect's
direction is correct in both models. Interlace has 51.13% lower MSE in this
matched seed42 comparison. Independent Metal replay reproduces all 11,724
prediction values exactly. Fresh local controls also pass: ordinary MSE
0.0000077815 and Interlace MSE 0.0000027699. Independent Metal replay reproduces
all 11,724 additional prediction values and confirms all five capability checks
for both architectures. Interlace lowers local-task MSE by 64.40% in this paired
comparison. See `attention-training-budget/` and
`attention-budget-local-confirmation/`.

A further candidate retains the Interlace gate and adds 64 coefficients that
carry signed and absolute base distances into attended values. This is a compact
adaptation of [Shaw et al. (2018), equation 3](https://aclanthology.org/N18-2074.pdf),
using a linear basis instead of the paper's clipped learned relative vectors.
The integrated full-length scorer passes independent CPU output and gradient
checks, including two deliberately broken implementations. Metal qualification
fails one gradient comparison at the unchanged float32 tolerance; no fit was
started. This optional candidate is deferred while the successful gate models
move to biological validation. Its failed check
and declared experiment remain in `attention-relative-values/`. These experiments
introduced no biological fits; their completed biological ledger was 247.

The qualified ordinary and Interlace scorers now enter a separately frozen
biological validation round in `attention-salomon-validation/`. Both start fresh
with seed42 and train on all 21,781 admitted edits in 8,774 components, using the
same 4,096-update recipe and the original common loss scale. They are evaluated
on 7,390 edits in 3,029 separate components. The existing k5 error and ranking
thresholds remain strict; both architectures are eligible, with minimum eligible
MSE and lexical architecture-name ties deciding nomination. Both fits completed
once, bringing the biological count to 249. Independent Metal replay reproduces
all 58,354 training, validation and pair-property values exactly:

| Scorer | Training MSE | Validation MSE | Signed top10 utility |
|---|---:|---:|---:|
| Ordinary attention | 0.0264657702 | 0.0291946709 | -0.0835244667 |
| Interlace gate | 0.0262335429 | 0.0289777495 | -0.0940942319 |
| Fixed k5 comparator | — | 0.0286251098 | 0.0351089021 |

Both candidates fail both biological requirements; no model is nominated and
confirmation remains sealed. A prospectively specified post-replay diagnostic
finds that even oracle positive rescaling cannot lower either validation MSE
below k5: the infima are 0.0287113703 and 0.0286757369. No calibrated candidate
was fitted or evaluated. The next source-backed step is separately admitted
element-activity supervision on training components, with its own output head
and physical sequence orientation preserved. It retains the signed paired
target and biological acceptance criteria. See the round's `RESULTS.md`,
`diagnostic/actual-1/RESULTS.json` and `PUBLISHED_METHODS.md`.

The activity-supervision trial is now complete in `attention-activity-training/`.
A whole-name intake admitted 30,800 physical activity measurements across the
existing 8,774 training families; 88 eligible unmeasured keys remain absent.
Independent source replay verified every admitted value, identity and physical
orientation. The first intake stopped before parsing coefficients because the
old source contract hashed the header without its trailing newline. A preserved
second revision corrected that header representation without changing the source
blob, parser, eligibility or missingness rules.

Each model adds a separate 33-parameter activity readout, totaling 11,401
parameters. Family-balanced activity targets use TRAIN-only centering/scaling,
physical orientation and a separate sampling stream. The paired target and its
loss denominator remain unchanged. Initial activity gradients were 81.34 and
108.59 times stronger than paired gradients in the ordinary and Interlace
trunks. A pooled gradient-energy rule fixed one common activity coefficient,
0.010327737397099176, before either fit. This is static initialization balance,
not dynamic GradNorm. Both fresh fits use the original seed42/4096 recipe;
activity supervision adds computation relative to historical paired-only runs.

| Joint scorer | Training edit MSE | Validation edit MSE | Signed top10 utility |
|---|---:|---:|---:|
| Ordinary attention | 0.0266293079 | 0.0290980677 | -0.0351043158 |
| Interlace gate | 0.0262956839 | 0.0287767719 | 0.0287056844 |

The two fits bring the biological ledger to 251. Interlace improves both
endpoints relative to its paired-only result, but its MSE remains 0.530% above
the strict threshold and utility remains 0.0064032177 below its threshold.
Neither model qualifies. Independent Metal replay exactly reproduces all
119,954 paired, activity and invariant-check values; confirmation remains sealed.
The activity training MSEs, 0.9665863352 and 0.9644827131 in standardized units,
describe training fit only. Source-based follow-up research isolates a second
shared attention/FFN pass with learned loop identifiers, matching the released
Interlace recurrence semantics. The ordinary control must receive the same
depth and identifiers. The current one-pass spacing result does not establish
an advantage from iterative composition. The alternative frozen-k5 residual
hypothesis remains separate; it retains 98.67% of TRAIN target squared magnitude
and has not been fitted. See `attention-activity-training/RESULTS.md` and
`followup-research/` for exact identities and scope.

The shared-depth capability study is complete in `attention-shared-depth/`.
Both models reuse their attention and feed-forward block for two passes, adding
a learned pass identifier before each pass and applying the final normalization
and readout once. The identifiers use a separate initialization stream, preserving
the original backbone and subsequent activity-head tensors. Each paired scorer
has 11,432 parameters. Four fresh final4,096 fits produce:

| Two-pass scorer | Local MSE | Interaction MSE | Nonzero interaction sign accuracy |
|---|---:|---:|---:|
| Ordinary attention | 0.0000037641 | 0.0014343516 | 99.8204% |
| Interlace gate | 0.0000042501 | 0.0006360701 | 100% |

Both pass all five predefined checks. Independent Metal replay reproduces all
23,448 final values exactly. Interlace interaction error is 55.65% lower than
the matched ordinary model. Historical one-pass interaction errors were lower
for both architectures, so this establishes retained capability rather than
a benefit from depth or iterative composition.

The biological follow-up in `attention-shared-depth-biology/` completed two
fresh final4,096 fits. Both assay heads use the same shared two-pass features,
totaling 11,465 parameters. The training populations, sampler streams, optimizer,
schedule, paired scale and historical activity coefficient remain fixed. No
gradient recalibration is performed. The original strict MSE and signed top10
utility thresholds still both apply.

| Two-pass joint scorer | Training edit MSE | Validation edit MSE | Signed top10 utility |
|---|---:|---:|---:|
| Ordinary attention | 0.0264659387 | 0.0292150736 | -0.0030449526 |
| Interlace gate | 0.0261816046 | 0.0291333610 | -0.0358172077 |

Neither model passes either requirement. Independent Metal replay reproduces
all 119,954 outputs exactly. The completed biological ledger is 253; confirmation
remains sealed. Interlace lowers TRAIN and validation error against its matched
ordinary control, but has worse selection utility. Both models have higher
validation MSE than their one-pass joint counterparts. This closes the tested
shared-depth recipe without a biological nomination.

The fixed-k5 residual trial is complete in `attention-fixed-prior/`. It uses the
qualified one-pass Interlace scorer, historical paired-only initialization and
recipe, original loss denominator, fixed baseline coefficients and both
acceptance requirements. Residuals are computed in float64 before the existing
training cast and included in checkpoint identity. Evaluation adds the baseline
and neural correction in float64 against the original targets. No depth
extension, activity objective or correction-weight search is included.

The one fresh final4,096 fit gives TRAIN MSE 0.0271146270, validation MSE
0.0288036124 and signed top10 utility 0.0394810097. Utility exceeds the fixed
comparator by 0.0043721075 (12.45%), while MSE remains 0.624% above its threshold.
The selection gate passes and the error gate fails; no candidate is nominated.
Independent Metal replay reproduces all 29,177 final values exactly, including
an independently reconstructed k5 contribution and residual checkpoint identity.
The completed biological ledger is 254; confirmation remains sealed. The
TRAIN-only component diagnostic favors increasing the correction rather than
shrinking it. Its best possible positive scalar gain lowers TRAIN MSE only to
0.0269622231; no calibrated predictor was produced.

A frozen-trunk diagnostic then computes the best linear readout on the existing
RC-averaged features, using float64 CPU inference and equal-family weighted SVD.
The numerical-rank TRAIN infimum is 0.0268690097 versus the current 0.0271146269
(0.906% lower). It retains 31 well-separated directions and discards the expected
LayerNorm null direction. The minimum-norm head has norm 34.62, beyond the
original schedule's fixed-trunk head bound of 6.47. This supports testing one
readout refit in `attention-linear-readout/`, preserving the existing trunk,
baseline and two acceptance criteria. The diagnostic optimizes the deployed
RC-averaged objective; earlier random-strand training also penalized disagreement
between orientations. It therefore does not isolate optimizer choice as the
cause of the gap. No validation predictions or fitted head were exported by
this diagnostic. Its inputs, tests and independent source review are retained
in `attention-fixed-prior/head-diagnosis/`.

The single readout refit in `attention-linear-readout/` is complete. It attains
TRAIN MSE 0.0268690098 in the original float32 Metal runtime, matching the
float64 optimum to 9.1e-11. Validation MSE rises to 0.0293084452, while signed
top10 utility rises to 0.0443042768 (26.19% above the comparator). Selection
passes and error fails. Independent feature reconstruction and least-squares
solution recover the deployed head bit-for-bit; all 12 non-readout tensors
remain identical and all 29,177 predictions replay exactly. This is one new
readout fit, bringing the completed ledger to 255, with no backbone or baseline
refits. No nomination is made. The improvement transfers to selection utility
but not to average effect-size error; a more powerful readout alone does not
close the biological acceptance test.

The TRAIN input-alias diagnostic finds no repeated physical input pairs. Reverse
complement folding creates 112 conflicting pairs of observations, covering
0.4262% of family-weighted mass and forcing at least 0.0000438587 MSE under that
invariance. This is 0.1554% of target energy; it does not explain most TRAIN error
or measure noise variance. Additional potential-cycle constraints are omitted.
The frozen calculation is in `architecture-diagnosis/training-aliases/completed/`.

Five-fold family calibration is complete in `attention-family-calibration/`.
Every fold refits its own k5 baseline, original-y loss scale and original
one-pass residual model without excluded-family labels. The excluded predictions
set one correction coefficient to 0.0334305183. Deployment reuses the original
full-TRAIN baseline and residual checkpoint; no full-data refit or validation
weight search occurs.

Validation MSE is 0.0286228342, just 0.00795% below the comparator, but signed
top10 utility falls to 0.0204151195, 41.85% below its requirement. Error passes
and selection fails; no nomination is made. Independent replay reconstructs
all five ridge solutions, fold checkpoints and excluded predictions, the
calibration arithmetic, and all final predictions. The eleven fits are complete
(five ridge, five neural, one scalar), bringing the ledger to 266. This closes
the declared calibration branch; changing bounds, weights or the head from
this validation result would be a new search.

The full-TRAIN gradient diagnostic in
`attention-fixed-prior/orientation-diagnosis/` is complete. At the original
residual checkpoint, the deployed-error and strand-disagreement gradients have
cosine -0.93442165. Their sum still gives local descent on deployed error. This
identifies opposing objectives at that checkpoint, not a training bug or a
generalization remedy. The original-target error agrees with the earlier
float64 readout diagnostic; no model weights change or new fit is counted.

The matched comparison is complete in `attention-strand-objective/`. Both
fresh Interlace models compute both orientations, with identical initialization,
family draws, updates and compute structure. `average` squares the error after
averaging the two edit scores; `expected` averages their two squared errors.
Both are eligible under the original strict gates, with minimum validation MSE
then lexical objective name determining nomination. Evidence for the alignment
hypothesis separately requires average to improve both metrics over expected.
`average` achieves TRAIN MSE 0.0231523775, validation MSE 0.0316932345 and utility
-0.0027484431. `expected` achieves TRAIN MSE 0.0226710484, validation MSE
0.0298840807 and utility -0.0435559398. Neither passes either gate. The first
improves utility relative to the second but worsens error, so the predeclared
alignment hypothesis is unsupported. Both fit TRAIN better than the original
random-strand residual model and generalize worse on validation. Independent
replay reproduces every saved prediction and verifies the completed ledger of
268. Close the objective comparison without a nomination or schedule extension.

The separate comparison of position-preserving real-motif versus
shuffled-motif frontend in `attention-position-motifs/`, using the original
random-strand objective, is complete. The prototype scans all 532 admitted JASPAR profiles
on both strands and projects their nonlinear tracks into the existing Interlace
block. Complete windows, source-only background moments, even-width strand
anchors and zero initialization pass independent CPU tests and deliberately
broken controls. Both arms have 45,416 trainable parameters. Native Metal qualification passes
all scores, paired losses and 14 parameter gradients at the real input size.
Under the 526-input freeze, real motifs reach TRAIN MSE 0.0220413128,
validation MSE 0.0303058035 and utility 0.0512856590. Utility is 46.1% above
the comparator, while error is 5.87% above its limit. The permuted control
reaches TRAIN MSE 0.0215161068, validation MSE 0.0300270864 and utility
-0.0241303565, failing both gates. Both arms were prospectively eligible;
neither is nominated. Motif-specific attribution required the real library
to improve both metrics over the control, so that hypothesis is unsupported.
Independent source/buffer reconstruction and native replay reproduce all
58,354 predictions exactly; the completed ledger is 270. Fixed motifs
plus attention have direct prior art in TIANA; the tested claim concerns this
paired neuronal-reporter endpoint, not invention of that general architecture.

The next bounded diagnostic in `attention-training-trajectory/` retains four
quarter-budget checkpoints from one exact repeat of the original calibration
fold-0 neural trajectory. It reuses the complement-trained k5 prior, preserves
the original 4096-update schedule and scores only TRAIN-complement and excluded
TRAIN families after fitting. Full final checkpoint state must match the
historical trajectory before interpretation. Component error and the identity
MSE(b+g)-MSE(b)=E[g²]-2E[g(y-b)] distinguish correction growth from transfer.
This diagnostic nominates no checkpoint and evaluates no outer validation or
confirmation. Native synthetic checks verify exact uninterrupted versus
four-segment training. Independent admission, label-isolation and failure tests
pass. The 612-input diagnostic completed its 4096 updates, but failed the
historical checkpoint equality check before inference. This attempt is counted
as fit 271 and has no prediction results or nominee. The first recorded training
objective mismatch is at update 576, before the first pause at 1024. Experiment
identity, scale, scheduler and all named random-generator states match; the
initiating cause is unresolved. Saved artifacts and an independent forensic
comparison are retained in `attention-training-trajectory/mismatch-review/`.
Full-shape synthetic gradient repetitions and two uninterrupted 640-update
synthetic fits match exactly. A CPU comparison also verifies identical packed
inputs and all 524,288 sampled IDs/strand choices, including augmented codes
and targets for the 73,728 pairs through the first divergent block. The
historical training discrepancy remains unresolved.

A separate secondary analysis in `attention-trajectory-secondary/` describes
the four fixed snapshots as their own recorded run. It preserves the failed
historical-recovery status and introduces no fit, calibration or nomination.
Both populations and all four checkpoints are fixed before scoring; independent
replay is required before interpretation. This development fold has previously
been used for calibration, so its learning curve cannot qualify a predictor or
establish an optimal stopping budget. Its purpose is to inform the next
prospectively specified experiment. The analysis and independent replay are now
complete: all 174,248 saved values agree exactly. Fitting MSE falls from
0.0274847033 at 1024 updates to 0.0256396979 at 4096 (6.7%); excluded-family
MSE rises from 0.0284318898 to 0.0290209407 (2.1%). No point beats that excluded
population's fixed-prior MSE of 0.0284220752. The correction's excluded-family
energy grows while its alignment with the target residual declines. This is
evidence of late deterioration in this one run and fold.

The sole candidate in `attention-prefix-budget/` is a fresh full-TRAIN
1024-update prefix of the original 4096-update training schedule. It uses the
original full-TRAIN baseline, residuals, loss scale, initialization and inference
path. The diagnostic snapshots remain ineligible. Fit 272 was counted before
training, and checkpoint progress and scheduler horizon were admitted before
validation conversion. It completed in 96.8 seconds including scoring: TRAIN MSE
0.0278152904, validation MSE 0.0286277662 and selection utility 0.0257376550.
MSE is 0.00928% above its limit and utility 26.69% below. Neither strict gate
passes; all 29,177 predictions replay exactly and no model is nominated.
Confirmation remains unopened. There is no continuation control, so the trial
tests candidate qualification rather than estimating a matched stopping-time
effect. The retained 4096-update OneCycle schedule is still rising at update
1024.

The subsequent `attention-short-cycle/` trial changes only the schedule horizon
to 1024 and keeps that same number of optimizer updates. It reuses the qualified
trial runtime and completes warmup and annealing before the fixed final endpoint.
The original OneCycle paper motivates a low-rate finishing phase; its reported
image-model results do not establish a gain for this Lion-trained sequence model.
The entire learning-rate path and cumulative exposure change, so this is a
schedule-allocation test, not an isolated cooldown comparison. See the trial's
`ASSESSMENT.md` and `SOURCES.json` for primary sources and exact rates.

Fit 273 completes in 101.4 seconds including scoring, with TRAIN MSE
0.0277553714, validation MSE 0.0286762398 and selection utility 0.0371846616.
Selection exceeds its requirement by 5.91%; MSE is still 0.17862% above its limit.
Independent replay reproduces all 29,177 predictions exactly. No candidate is
nominated and confirmation remains unopened. Compared with the stopped prefix,
utility improves and MSE worsens. Compared with the original 4096-update
correction model, MSE is lower and selection still passes, using one quarter of
the training updates. Both comparisons describe single runs. The bounded
short-budget trials are closed without a joint-gate success.

The completed trial, `attention-short-cycle-calibration/`, calibrates this
different correction using the existing five TRAIN-family folds. It reuses the
verified complement k5 coefficients and original-y scales; only five fresh
1024-update neural fits and one scalar fit are allowed. The coefficient uses
the existing clipped, pooled, equal-family least-squares rule and direct neural
corrections. An interior coefficient is applied to the unchanged full-TRAIN
checkpoint from fit 273 and tested once against both original gates. Exact
zero or one closes as a known failed endpoint without another validation score.
The earlier 4096-model calibration lost ranking utility. No nearby-coefficient
search or additional full-model fit is allowed.

The five fresh folds set alpha to 0.2372804265. Validation selection utility
reaches 0.0760924857, 116.73% above the comparator; MSE is 0.0286310126,
0.02062% above its strict limit. Independent replay reproduces all 109,300
prediction values exactly and verifies every fold, checkpoint, coefficient,
metric and accounting entry. Six fits bring the ledger to 279. The joint gate
remains unmet, no model is nominated, and confirmation remains unopened.
The utility gain is a point estimate after repeated validation-based development,
not a confidence bound. This calibration trial is closed with its original result.

The completed experiment, `attention-two-cycle-mixture/`, reuses the exact
excluded-family corrections from the 4096- and 1024-update models. One TRAIN-only
fit learns two nonnegative weights whose sum cannot exceed one; no neural or
ridge retraining is permitted. The same fixed coefficients combine their already
replayed full-TRAIN predictions for one validation evaluation. An optimum using
only one correction closes as an existing candidate before validation. This is
an ensemble construction, with both original gates unchanged and independent
coefficient, input-alignment and result checks required before advancement.

Fit 280 sets long and short weights to 0.0206557439 and 0.2223167744. MSE is
0.0286293465, 0.01480% above its limit; utility is 0.0715850998, 103.89% above
the comparator. The absolute MSE gap is 28.23% smaller than the single short
calibration's gap. Independent CPU replay verifies all 29,560 saved validation
values exactly and checks the constrained optimum with a separate profiled
solver and convex optimality conditions. No neural fitting or inference was
repeated. This fixed library closes without a nominee; confirmation remains
unopened. The subsequent TRAIN uncertainty diagnostic is completed below.

The separate `train-coefficient-uncertainty/` audit is complete. It pins limma
3.58.1 and filters original TRAIN identifiers before numeric conversion. All
21,781 coefficients in 8,774 families have recoverable moderated standard errors;
no values are unresolved. Equal-family weighting gives squared effect S =
0.028225847 and squared standard error N = 0.017436895, so N/S = 61.78%.
This is model-based coefficient uncertainty, with normalization and empirical-Bayes
moderation inherited from the published fit; it does not establish an empirical
noise floor. Baseline family-excluded MSE is 0.028175001; adding the long or short
unit correction increases it to 0.029025549 or 0.028212990. Both increase MSE in
every fixed uncertainty quartile. This finding concerns uncalibrated corrections,
not the fitted mixtures.

The independent checker reconstructs all 21,781 TRAIN rows exactly and agrees
on 53 summary values within 5.56e-17. No new fit, validation evaluation,
confirmation access, nomination or gate change occurs. The ledger remains 280.
This closes the diagnostic; any next predictor intervention requires its own
fixed protocol. See `architecture-diagnosis/TRAINING_NOISE_AUDIT.md` for the
initial source audit and `train-coefficient-uncertainty/RESULTS.md` for its result.

The repository now includes the [experimental inference example](../examples/neural-predictor/README.md),
with the exact attention source, tensor-only exports of both full-TRAIN models,
the k5 coefficients, the fixed mixture and aggregate evidence receipts. This
makes the fitted sequence-to-effect calculation available without the local
research directory. Training orchestration, raw assay tables and the reserved
confirmation are outside this inference package. No model is promoted by
packaging it.

The subsequent `attention-precision-sampling/` trial uses reported TRAIN
coefficient SEs in one prospectively fixed sampling rule. In each fitting
population, equal-family probabilities p become proportional to p/(SE²+vbar),
where vbar is that population's equal-family mean SE². Original targets,
complement priors, loss scales, initialization and the complete 1024-update
recipe stay unchanged. Six fresh neural fits and one equal-family scalar
calibration advance the ledger from 280 to 287; no ridge is refitted. Only the
final calibrated candidate receives a validation evaluation.

The fitted alpha is 0.4245090444. MSE is 0.0286209346, 0.01459% below the original
k5 limit, while utility is 0.0034205686, below random selection and 90.26% below
the required comparator. This specific precision rule fails joint acceptance.
All 109,300 saved prediction values replay exactly. Independent review verifies
the fit-local sampling, source joins, all checkpoints, folds, calibration,
metrics and seven-fit ledger. No confirmation is opened and no candidate is
nominated. The trial closes without a tempering grid, additional mixture or
metadata-causality claim. See `attention-precision-sampling/RESULTS.md`.

The subsequent TRAIN ranking diagnostic is complete, without additional fits.
It examines the unscaled long, short and precision-sampled corrections in the
five original excluded-family folds. Their pooled equal-family MSEs are
0.0290255495, 0.0282129898 and 0.0281834981, respectively, against the common
fold baseline's 0.0281750006. The long correction worsens MSE in all five folds;
short improves one fold and precision improves three. Correction energy and
residual alignment identify excess correction magnitude as part of the observed
error, without identifying a particular architectural cause.

Each correction improves signed top-ten utility in three folds and worsens it
in two. In nine of the fifteen comparisons, deleting one common family can
reverse the sign of the utility difference. All 295 family deletions include
the next eligible replacement. These are descriptive sensitivity checks on
fixed predictions, not confidence intervals or five independent confirmations.
Independent replay verifies all 20 cells, 15 comparisons and 295 deletions;
the maximum numerical difference is 5.56e-17. Fifteen synthetic tests and eight
deliberately broken controls qualify the diagnostic. It opens no validation or
confirmation outcomes. See `train-ranking-diagnostic/RESULTS.md`.

The next bounded trial, `kmer-resolution-mixture/`, tests complementary errors
between the fixed k3 λ=0.01 and k5 λ=0.1 predictors. Five new k3 fits use the
original family-fold complements and their own normalized family weights. A
single convex weight is then fitted to the assembled excluded-family predictions.
The k5 complement predictions and both full-TRAIN coefficient arrays are reused.
Selecting k3 after its earlier validation ranking makes this adaptive development;
the two original baseline roles and strict acceptance thresholds stay fixed.

The learned k3 weight is 0.3699603138. The one validation evaluation has MSE
0.0286300007 and signed top-ten utility −0.0363062730. MSE is 0.01709% above
k5, and utility is below the random comparator. Both gates fail. The six fits
advance the ledger from 287 to 293 without neural training. Both original
component prediction vectors reconstruct exactly. The trial closes without a
nominee, a weight sweep or a change to the bundled predictor.

Independent replay uses weighted-design SVD for the five new complement fits;
their largest coefficient difference is 5.80e-16. All 22,170 validation
prediction values and selected IDs reproduce exactly. Nine synthetic tests
and nine broken controls qualify the implementation. A root rerun of the same
checker agrees. These checks add no candidate or fitting choice. See
`kmer-resolution-mixture/RESULTS.md`.

A separate source audit identifies a pinned 83 MB barcode-count file with
three ordinal DNA/RNA assay pairs. A header-only probe confirms the columns;
the inspected workflow writes nullable integer counts, with blank cells
distinct from zero. A design-only check maps all 21,781 original TRAIN edits
across 8,774 families to 30,752 exact count names without split conflicts.

The TRAIN-only importer is implemented and frozen. Ten local tests, thirty
independent tests and nine deliberately broken controls exercise admission,
missing values, source identity, split isolation and atomic publication. The
review caught and repaired an integer-range check and a family-isolation check.
The importer retains original membership and skips non-TRAIN rows before count
conversion. The approved download and measurement run are complete: all 30,752
construct names and both alleles of all 21,781 TRAIN pairs are present. Exact-file
sample accessions, UMI merging and possible upstream pseudocounts remain unresolved.
This work adds no model fit or nomination. See `train-replicate-intake/PLAN.md`
and `train-replicate-intake/RESULTS.md` in the research work directory.

The corresponding TRAIN agreement analysis is also implemented. It compares
reference-to-alternate log count ratios across the three columns on common
barcode support, retains original family weights, and reports coverage,
centered correlation and direct-copy error against zero. Thirteen local tests,
ten independent test methods and nine deliberately broken controls qualify the
implementation. Separate synthetic cases show why correlated reporter activity,
constant effects and shared nonzero means cannot establish repeatable edit
effects. The fixed estimator and source limitations are in
`train-replicate-agreement/PLAN.md`.

The measured run admits 2,801,763 barcode rows and retains 1,910,615 (68.19%)
with all six DNA/RNA cells observed. All 21,781 edits in 8,774 families retain
positive aggregate counts and their original family weights. Allelic Pearson
correlations are 0.20311, 0.15860 and 0.18225 for column pairs 1/2, 1/3 and
2/3. Copying another column's raw effect increases MSE by 53.15%–72.48% versus
predicting zero. These are descriptive count-contrast results, not a sequence
prediction ceiling or a reconstruction of the published pooled BCalm targets.
Independent Decimal/Fraction replay checks all 65,343 effects, every integer
aggregate and coverage field, and all summary metrics; maximum effect and
summary differences are 2.22e-16 and 1.11e-16. See
`replicate-measurements-20261010/RESULTS.md`.

A subsequent fixed comparison uses the previously admitted family-excluded
k5 predictions and uncalibrated long, short and precision corrections. The
long and short corrections worsen pooled MSE against all three count columns;
precision improves column 2 and worsens columns 1 and 3. Short and precision
have small positive centered covariance with all three k5 residual columns,
but their correction energy exceeds twice the raw alignment in five of six
pooled comparisons. Long has negative residual covariance in column 1 and
positive covariance in columns 2 and 3. This distinguishes a magnitude problem
from agreement across columns without choosing a new scale or model.

The comparison retains all original pairs, weights and five excluded-family
folds. Eight tests and seven deliberate faults qualify it. The 90-file freeze
preserves the earlier input identities before joining predictions to effects.
Independent replay verifies all 72 model/column comparisons and 54 difference
covariances; all 1,494 numeric fields match within 2.22e-16.
The audit is adaptive descriptive work: original model training used pooled
targets derived from all three assay columns, so it is not an assay-held-out
prediction test. It adds no fit, new validation evaluation or nomination.
See `train-replicate-prediction-audit/RESULTS.md`.

The subsequent `attention-quadratic-readout/` trial adds one zero-initialized
coefficient per squared pooled feature to the Interlace scorer: 32 added
parameters, 11,400 total. Six models use the original complete 1024-update
recipe, equal-family sampling, targets and frozen k5 priors. The one
excluded-family TRAIN calibration chooses alpha 0.15097202949412603. Seven
begun and completed fits advance the ledger from 293 to 300. The complete
run takes 513.86 seconds on the local Metal path, including preparation and
scoring; this is one run duration, not a comparative performance benchmark.

Validation MSE is 0.028628390522333413, above 0.028625109823009157. The candidate
chooses the same ten ids and directions as k5, in a different order. Its
NumPy-mean utility is 0.03510890212666124, one step above the stored baseline
0.035108902126661234. Independent compensated sums and exact rational means
both give 0.035108902126661234 for both selections. Their exact utility gain
is zero. The original checker correctly stops on the differing gate flags;
retain that failure and the frozen runner record. Joint acceptance fails
under either arithmetic, with no nomination or confirmation. A separate
companion verifies all prediction arrays exactly and explicitly reports the
metric discrepancy; the original checker remains failed. The parent saved-array
reconstruction independently agrees. This trial closes; no nearby coefficient
or checkpoint search is admitted.

Qualification includes 18 independent CPU cases, nine deliberately broken
model/scoring variants, and repeated runner tests covering 12 orchestration
cases, 42 admission cases and nine faults. Two unnormalized CPU/Metal
gradient stress checks remain failed and documented. A separate predeclared
check of the actual normalized paired update, gradient accumulation, clipping
and inference agrees with float64 CPU under unchanged float32 tolerances.
Those checks establish the exercised arithmetic, not universal backend parity
or biological accuracy. See `attention-quadratic-readout/RESULTS.md`.

The saved-array utility correction is complete in `ranking-metric-repair/`.
It reuses the existing compensated-sum reference for both selected baselines
and six fixed candidates, while retaining the original MSE values and strict
thresholds. Only the quadratic utility flag changes; no joint acceptance
decision changes. The original runner fails the cancellation regression; the
reference passes five contract tests, including ranking permutations, signs,
family selection and strict improvement. Independent review reproduces all
eight results. No fit or inference runs during this correction.

A future runner must adopt only the utility correction and retain its admitted
input checks and original MSE calculation. Substituting the entire reference
metric also changes MSE reduction order and is outside this repair. Exact
rational means remain the independent oracle; compensated sums are stable
floating-point arithmetic, not a guarantee of correctly rounded rational
means for every input. The original baseline utility threshold is reproduced
exactly here. Frozen historical sources and records remain unchanged.

The assay source follow-on clarifies separate barcode-export and pooled-fit
processing: upstream counts may include configured offsets and filtering; the
pooled model performs its own normalization and weighting. Exact historical
settings and the author-column-to-culture mapping remain unresolved. Those
ordinal count ratios cannot be substituted for the original pooled target or
used as its noise ceiling. Both inspected public BCalm tags override the
inherited mpra response calculation and preserve missing barcode slots as NA;
this corrects an earlier default-accessor interpretation. Neither tag establishes
the historical `0.9.0` identity. The public history search and inspected publisher
supplement do not supply the missing execution record. See
`replicate-measurements-20261010/SOURCE_SEMANTICS.md`, `PROVENANCE_SEARCH.md`
and `SUPPLEMENT_READ.md` in that directory.

A renewed Matoba source audit locates and verifies the author repository's
150-base FASTA: 22,758 distinct sequences form 11,379 A1/A2 pairs, each differing
once at position 76. One pair has a missing SNP identifier and remains explicitly
unidentified. Independent parsing reproduces every pair. This is a candidate
design for testing DNA effects under vehicle versus Wnt stimulation in primary
neural progenitors. The source code uses the FASTA for motif analysis; its exact
connection to the full synthesized and barcode-mapped library remains unresolved.
The complete authorized barcode-map join now verifies 22,550 named FASTA
sequences by exact reverse complement, with no same-name sequence conflicts
or duplicate names. The map lacks 208 FASTA entries: 11,263 original pairs
have both alleles, 24 have only one and 92 have neither. Independent streaming
replay reproduces every projected row. These omissions remain explicit in
the original design population; biological orientation is still unresolved.
Official GEO metadata establishes **GSE289858** as the MPRA accession; the author
README's GSE289857 points to an unrelated leukemia study. Published allelic
results are activity-selected, and interaction results are further selected by
allelic significance, so those tables alone cannot support an all-variant
benchmark. No measurement file was opened. Source identity, sequence checks,
processing discrepancies and the remaining admission steps are recorded in
`matoba-source-admission/RESULTS.md`; this adds no fits or confirmation result.
Further source tracing establishes the intended condition-response endpoint as
`(A2−A1)Stim − (A2−A1)Vehicle`. Historical filtering remains unresolved: source
defaults use RNA count 1 while the paper specifies 5, and the retained-well
rules also differ. The prospective measurement contract is now fixed and
qualified on controlled counts. Coordinate metadata must additionally establish assembly,
strand and full insert overlap before the 205 preliminary lead-locus groups
can become independent evaluation families. These source checks precede any
new model or outcome-selected population.

The upstream coordinate-only projection is now complete. Its 9,331 annotation
rows each match one exact FASTA pair identifier, declared lead locus and both
allele bases; there are no duplicate or contradictory matches. **2,048 of the
11,379 design pairs have no annotation in this file.** All seven design fields
are retained separately from the thirteen association columns. The file supplies
no assembly or strand, so it cannot yet establish complete genomic families.
Seven whole locus groups and parts of 72 others lack annotations. The bounded
author-history review does not locate a generator or documented omission rule;
the missing pairs stay in the population accounting.
The isolated RDS reader is qualified against R-generated fixtures, with explicit
short-read and negative-length rejection. Source identity, every missing key,
checks and projection are retained in
`matoba-source-admission/coordinate-reader/`. The complete map/coordinate
intersection contains 9,236 pairs, with 2,143 original pairs lacking at least
one of those criteria. That intersection is not silently admitted as a smaller
population. The approved broader 11.96-MB source has now been downloaded and
inspected structurally. It contains 678,579 rows with SNP/allele/phenotype names
and GWAS statistics, but no chromosome, coordinate, assembly or strand column.
The statistics were not converted or evaluated. This file supplies no direct
coordinate table; the missing keys need a separate reference-variant route.

The existing family-graph helper can be reused unchanged after complete
placement admission. Four source-only controls pass, including 100 small
comparisons with an independent pairwise graph. Shared loci, chained and
one-base overlaps, chromosome boundaries and exact DNA/reverse complements
are covered. The caller must supply both allele records as a materialized
list and verify assembly, contigs and SNP consistency before grouping.
An official UCSC hg19 sequence query also resolves the first annotation-row
example: all 149 unchanged bases match under one-based, plus-strand placement,
and the full reference matches A2. Root independently reproduces that result;
the other origin/strand hypotheses differ at 105–121 flank positions.
This establishes a working reference-check route and a concrete A2-reference
example, not whole-library placement. No actual families or split have been
assigned. Evidence is in `matoba-source-admission/family-graph-check/`.

The reviewed batch reference checker has completed all 152 public GRCh37
batches. All 9,331 available-coordinate pairs have exactly one matching
origin/strand hypothesis, with all 149 unchanged bases matching and the reference
edit base equal to A1 or A2. The retained response bodies total 2,102,275 bytes.
Connection timeouts were resumed from verified complete prefixes; prior attempts
remain preserved. Independent replay from the original FASTA reproduces all
37,324 hypothesis records and the exact output bytes; root reran that replay.
All matching hypotheses use one-based coordinates and forward inserts. The
reference matches A1 for 4,465 pairs and A2 for 4,866. Evidence is in
`matoba-source-admission/reference-check/continued-4/` and
`complete-reference-review/`. This resolves placement
at supplied coordinates; global genomic uniqueness and the 2,048 missing
coordinates remain separate questions.

A separate official variant lookup recovers the first missing known-rs design,
rs12700001, with a complete A1 reference match. Root reproduces both strand
comparisons. The missing-ID projector and fixed-origin comparator now pass
seven independent controls and six root controls; eight deliberate faults
are detected. They preserve aliases, mapping multiplicity and unsupported
mappings. The complete lookup has now returned all 2,047 known missing IDs,
including 39 explicit aliases. All 3,836 returned mapping occurrences were
checked against public reference sequence in 44 batches. Resolving these
occurrences back to the original 2,048 design rows gives **1,367 unique,
680 ambiguous and one unknown-ID row**. Every named design has at least one
one-based placement with all 149 unchanged bases matching, the reference edit
base equal to A1 or A2, and both source alleles present in the database mapping.
The 680 ambiguous rows retain exact matches on chromosome 17 and alternate
reference representations; they are neither duplicated as examples nor assigned
an arbitrary winning location. These are placement results, not predictor scores.
The complete ledger is in `reference-check/variant-placement-resolution/`.
Independent input replay verifies all original owners, mapping occurrences,
aliases and request windows; root reproduced it. Original-row result review
also reconciles every row and all mapping checks. Ten independent controls pass
and five deliberate faults fail; root reran the controls. This result review
checks saved evidence and aggregation, not an independent replay of the actual
source DNA. A later coordinate-projection stage resolves alternate equivalence
as described below; it preserves this original strict placement ledger.
The one unknown-ID insert has no exact 149-flank match anywhere in its declared locus's bounded
119,923-base reference neighborhood, while all 98 known locus controls match.
Root reproduced that result from the original FASTA. The unknown pair stays
explicitly unresolved. The subsequent prospective source policy excludes it from
every supervised role before count access; its full identity and declared
lead-locus family remain recorded. A broader placement search can be separate
future work and does not create an independent singleton here.

The official Ensembl release-116 GRCh37 component tables bind `HSCHR17_1` and
`HG1146_PATCH` to accessioned scaffolds `GL000258.1` and `JH720447.1`. NCBI GFF3
alignments then project each complete insert and its edit position to primary
chromosome 17. Of the 680 originally ambiguous rows, 676 agree completely;
three have incomplete alignment and one has a 37-base discrepancy. Seven root
and five independent controls pass, including all strand combinations and
alignment gaps. Root reran the independent controls. An independent coordinate
replay reproduces all 805 alternate comparisons. It checks coordinate equivalence;
the earlier reference checks supply the separate DNA-match evidence.

`reference-check/ALTERNATE_PROJECTION_FREEZE.json` preserves this result.
The four unresolved rows retain their exact primary intervals and uncertainty.
Their official MAPT parent envelope is used only to join potentially related
families, reducing 179 primary-only components to 178. The existing admission
plan explicitly permits conservative bounds and a frozen source-only exclusion
for the unknown pair. `PARTITION_PLAN.md` adopts both options before assigning
the fixed hash-modulo-10 TRAIN/validation/test split. All 11,379 source rows
remain accounted for; 11,378 are eligible before measurement QC. Available
canonical IDs and aliases have no shared owners across this population.
Independent graph replay reproduces all assignments, hashes and partitions.
Root reran six independent code/join checks and the complete graph replay.
Eight deliberate code faults and four graph-replay faults fail their controls.
`PARTITION_FREEZE.json` records source admission before any count body was opened.
The fixed split has 7,543/1,645/2,190 TRAIN/validation/test source pairs in
109/38/31 families. The excluded unknown retains its source family but cannot
enter any supervised role. No family is reassigned after measurement QC.

Full count ingestion completed through the qualified adapter on one
low-priority CPU in 1,245.8 seconds, using 34.7 MB peak memory and no GPU.
It validated every row in both complete compressed streams and retained
1,221,929 Vehicle and 1,223,082 Stim rows with unique global barcode assignments.
Post-import reconciliation finds no mismatched barcode owner or undeclared
condition. The subsequent visible measurement work is TRAIN-only. A
selector review found and repaired a mutable-row isolation defect before use:
validated TRAIN execution identities are now captured before the first count
query. An independent regression fails on the prior code and passes on the
repair; root reran it and the three selector checks. The fixed diagnostic
compares alternating seven-well sets separately within each condition, retaining
the original full interaction target and giving eligible families equal weight.
Its agreement is conditional on shared QC; it is not independent replication
or a sequence-prediction noise ceiling. The atomic export retains every one of
the 7,543 TRAIN rows, including 613 explicitly rejected measurements. All 109
families retain eligible pairs. This measurement export preceded the bounded
model comparison reported below.

| TRAIN target | Half-well correlation | Zero-half MSE | Copy-other-half MSE |
|---|---:|---:|---:|
| Vehicle allele effect | 0.49713141 | 0.02964355 | 0.03014207 |
| Stim allele effect | 0.63152902 | 0.02703192 | 0.02015096 |
| Wnt interaction | 0.14221374 | 0.03069569 | 0.05274442 |

Every nonempty eligible family has equal weight. The full interaction's
zero-effect MSE is 0.01750798 in squared log2 reporter-effect units, not an
error percentage. Copying one interaction half to predict the other is worse
than zero; this does not determine a sequence-model noise ceiling or prove
that the interaction is unlearnable. The stronger within-condition allele
agreement motivates a controlled comparison of shared and direct interaction
learning, while retaining the Wnt interaction as the primary target. The actual
outputs and fixed method are in `matoba-measurement/train-measurements/` and
`TRAIN_DIAGNOSTIC_PLAN.md`. Independent replay reproduces every TRAIN row,
eligibility decision, reduced well/condition/full/half target and all family
contributions. Its 70,652 numerical comparisons pass; root reruns the checker
and reproduces the complete result. This checks arithmetic after barcode QC,
not an independent experiment.

The bounded sharing comparison is now implemented and has completed its 15
frozen attempts. It compares zero, a fold-complement mean, direct k5 ridge,
and rank-one joint prediction of the common allele effect and Wnt interaction.
The fixed k5 penalty is 0.1, with complement-only RMS scales and five complete
family folds. Independent linear condition heads would merely reproduce direct
interaction ridge; rank one introduces an actual sharing constraint.

Prospective review changed primary selection to two distinct-family edits
strictly within each excluded fold, ten choices across five batches. Pooled
out-of-fold top-ten selection would allow a family's labels to affect its own
selection through other folds' fitted models. An excluded-target perturbation
test now verifies unchanged coefficients, scales, predictions and selection
within that fold. Exact rational arithmetic on saved IEEE values determines
both progression endpoints. All 19 focused and independent synthetic checks
pass; deliberately broken solver and caller variants are rejected.

| TRAIN excluded-family model | Family-weighted Wnt MSE | Signed selection utility |
|---|---:|---:|
| Zero effect | 0.0175079789 | 0.0000000000 |
| Complement mean | 0.0175332147 | 0.0832594014 |
| Direct k5 ridge | 0.0178180113 | 0.0007145589 |
| Shared rank-one k5 ridge | 0.0177693295 | -0.0339406594 |

All 6,930 eligible pairs and 109 families remain in the comparison, with all
613 rejections retained. Sharing improves direct ridge's MSE by 0.2732%, but
its MSE is still 1.4928% above zero and selection worsens. Both progression
requirements fail. The mean control uses a fold-level offset and lexical tied
choices; its positive utility does not establish sequence discrimination.
This closes the fixed rank-one/k5/0.1 recipe without a parameter sweep or a
nomination. It does not identify measurement uncertainty as the sole cause.
Independent result replay verifies every source row, fold, coefficient,
prediction, exact metric, signed selection and ledger event. Root reproduces
all replay fields except runtime and peak memory. `matoba-model-pilot/PLAN.md`,
`FREEZE.json` and `experiment/complete/` retain the exact prospective contract,
coefficients, predictions and fit ledger. Validation and test labels stay closed.

A prospective count-to-response core is now implemented in
`matoba-measurement/`. It excludes wells 12/13 before QC, uses fixed paired
DNA/RNA thresholds and one-pass Tukey outliers, and preserves missing values.
Each condition requires at least ten shared physical wells for A1 and A2.
The target is the difference between condition means of paired-well A2−A1
log2 reporter ratios. It is an explicit new reanalysis, distinct from the
published mpralm normalization and weighted coefficients. Fourteen root and
nine independent controls pass; eight deliberate faults are detected. Review
repaired an equal-ratio rounding error that could change barcode eligibility.

The complete barcode index contains **2,880,805 unique assignments**. Its
global occurrence rule excludes 19,117 duplicated identifiers, including
same-owner repetitions. All 56 source variants with empty barcode fields
remain represented with zero assignments. Of the original 11,379 pairs,
11,234 have assignments for both alleles, 32 for only one, and 113 for neither;
every pair remains accounted for. The complete index took 23.92 seconds and
about 60 MB peak memory on one CPU. Seven root and seven final independent
controls pass, including empty-field preservation and failure cleanup.

All 64 sample-to-column identities are now verified directly from GEO
metadata. A bounded header-only probe confirms the two count matrices'
33-column schemas; their combined compressed size is 1,401,246,348 bytes.
The raw-count adapter connects synthetic map assignments and condition CSVs
to the qualified measurement core. Its transaction tests caught and repaired
partial imports under SQLite autocommit; all eight root controls now pass.
Independent adapter checks pass across all three transaction modes and detect
four deliberately wrong behaviors. The complete-file runner now passes five
root and eight independent controls; eight removed protections are detected.
It binds conditions to a separate source contract, verifies complete compressed
bytes and a private map snapshot, and publishes only a complete database.
Independent review reproduced and fixed a condition-swap defect in the original
caller; root replayed the final independent checks. Both approved full count
downloads are complete. Independent streamed hashes reproduce both downloaded
identities, exact byte counts and earlier source prefixes. `SOURCE_CONTRACT.json`
binds each condition to the verified public file, with source metadata and
approval provenance in `DOWNLOADS_VERIFIED.json`. The source-only partition,
complete count import, TRAIN measurement export and first bounded model
comparison are now completed as recorded above. Held-out evaluation remains ahead.
Evidence is in `matoba-measurement/barcode-index/` and `count-source/`.

The immediate scientific decision is whether the paired condition effects are
repeatable across measurements and predictable across independent sequence
families. More observations or a different architecture do not settle that
question. [Guo et al. 2023](https://www.nature.com/articles/s41588-023-01533-5)
is a complementary source candidate with measured allelic effects across neural
cell states. A source-only inventory identifies 23 deposited RNA files across
seven ES/neural/astrocyte contexts and a shared plasmid baseline. Exact insert
DNA, replicate-to-file correspondence and count processing still need resolution;
these files are not assumed to be 23 independent biological replicates. Its
145-base lentiviral assay needs a separate adapter from Matoba's 150-base,
paired-well Wnt assay. `guo-source-feasibility/` retains the primary-source review;
count bodies and outcome workbooks remain unopened.
Lee remains the previously recorded unresolved library-design route. These
datasets cannot be pooled as interchangeable Wnt-response labels.
[Nucleotide Transformer](https://pmc.ncbi.nlm.nih.gov/articles/PMC11810778/)
provides evidence for genome pretraining followed by supervised adaptation.
The [TACO preprint](https://arxiv.org/html/2503.07981v1) instead uses RL to design
sequences under learned fitness rewards and documents surrogate exploitation.
These are different tasks: pretraining is a testable representation choice;
RL would need an independently qualified reward and candidate evaluation.
This literature comparison adds no fits and changes no acceptance criteria.

The next representation investigation admits 347 representative composite
binding profiles from the pinned Xie2025 CAP-SELEX source. All differ exactly
from the prior 532-profile library, including reverse complements. A fixed
source-only alignment produces 312 independently composed single-factor
controls; 35 composites remain unavailable because a partner is missing or
too wide. Independent code reproduces every transformed matrix, selection,
placement and control exactly. This supplies a concrete comparison for joint
sequence preference beyond separately measured factor preferences. An exact
synthetic witness also distinguishes sequences with identical pooled k2–k5
counts, and independent high-precision scoring checks all 624 matrices. The
fixed no-fit diagnostic on original TRAIN fold 0 then conditions these features
against k5 and independent-factor features. It returns **STOP**: the three
cross-assay alignments are 0.000407700193, -0.000472964916 and -0.000905385925.
All 1,755 whole-family deletion means remain negative. Independent replay checks
every global and deletion term and reproduces STOP. The fixed representation
is closed without another fold or feature-selection retry. This was a
descriptive TRAIN screen, with no new predictor fit or confirmation access;
the fit ledger stays at 300. Source construction and the complete frozen result
are recorded in `composite-motif-source/RESULTS.md` and
`composite-train-diagnostic/RESULTS.md` in the research work directory.

Each successor uses the smallest model justified by training/validation evidence.
An available pretrained model is eligible only after its exact input, output,
training data, checkpoint, license and runtime are verified. Do not invent assay
flanks, convert scores into biochemical rates, or treat unspecified pretraining
overlap as independent validation. Exclude evaluation families whose outcomes
informed supervised pretraining, checkpoint selection, tuning or calibration.
If that exposure cannot be resolved, the model's result is exploratory only.
If the chosen dataset cannot support the declared experiment, return to phase 1
with the reason recorded.

## Rummel neuron-stimulation source qualification

The final publisher's Table S2 is now recovered and hashed. Its actual download
link and four sheet schemas resolve the earlier access limitation. Source-only
inspection finds 2,559 DevLib DNA pairs, each 160 bases long with one declared
single-base difference. Public GEO metadata establish three baseline and three
KCl RNA samples in N1 human neurons. The published treatment is five hours at
55 mM KCl on day 49. eQTLLib has three deposited human baseline samples but no
deposited human KCl GSM; it cannot supply a verified matched stimulation cohort.

The workbook also exposes specific input hazards. In eQTLLib, 2,243 sequence
pairs disagree with the order of their declared reference/alternate bases.
Across libraries, 1,542 identical DNA pairs have numeric coordinate labels offset
by two. After six fixed reference-genome witnesses, a complete DevLib check
matches all 2,559 original pairs locally to GRCh37. Every match is forward-strand,
with the full F/refSequence allele matching exactly and the edited base two
positions before the numeric name. This resolves the offset row by row; it
does not establish genome-wide unique mapping.

The source sequences retain restriction sites that the methods say were changed
before synthesis. DevLib has 152 pairs containing KpnI or XbaI sites; 15 overlap
the single edited base. The third listed enzyme, SfiI, adds 17 affected pairs,
making 169 in total. The remaining 2,390 pass this provisional sequence screen.
This is not a measured usable-pair count or an admitted evaluation population.
No substitution, allele swap or coordinate correction has been guessed.

Independent source parsing and root replay reproduce the sequence counts and
cross-library offsets. Independent SfiI scanning reproduces the three-enzyme
union, and literal negative controls reject incorrect projection, allele
position and motif boundaries. All 2,390 provisionally screened rows also pass
the complete reference check. The six initial witnesses replay locally.
Root also replays all 2,559 complete-library placements at their pinned origin.
The saved source, audits and current decision are in
`rummel-source-admission/ADMISSION.md` in the research work directory.

Source grouping now reproduces 2,059 strict overlap/exact-DNA components and
22 chromosome groups across all 2,559 pairs. Independent graph replay agrees on
every membership and hash. Whole-chromosome grouping is a conservative candidate
because the source supplies no complete lead-locus membership map. It also leaves
few independent evaluation units: chromosome 6 alone contains 947 of the 2,390
provisionally screened pairs. Nine exposed variant IDs now have verified GRCh37
placements, seven at exact DevLib edits. Their chromosome groups and unresolved
gene-only or prior-study exposure must be accounted for before a fresh split.
No evaluation roles have been assigned. Details and fault checks are in
`rummel-source-admission/family-feasibility/ASSESSMENT.md`.

The complete prior-source audit finds 280 distinct DevLib pairs overlapping
earlier libraries: 218 with Kosicki, 55 with Salomon and 10 with Matoba; these
counts overlap and must not be added. H1 has no detected interval or full-insert
containment match. Both alleles and all earlier physical source copies are
retained. Independent replay reproduces all 1,016 hit records. Cross-assembly
comparison uses every aligned chain block, verified against an independent
parser and 589 existing regions; one DevLib window retains a three-base mapping
gap. Source membership remains distinct from outcome exposure or fitted
membership. Under the conservative whole-source, whole-chromosome exclusion,
21 of 22 groups are touched, leaving five pairs on chromosome 13. This cannot
support fresh confirmation. The next design step must establish defensible
finer genomic families or retain this cohort for a separately declared
exploratory role. No acceptance requirement or evaluation assignment changed.
The evidence and limits are in
`rummel-source-admission/prior-overlap/ASSESSMENT.md`.

The next route is prospective, within-assay exploratory development. A contract
review confirms that whole-chromosome exclusion is an optional fallback, not a
general requirement of the sequence-family contract. Propagating prior-source
hits, known published-example matches and the mapping gap through local
sequence/interval components leaves 2,260 source pairs, including 2,105 passing
the provisional sequence screen, in 1,818 components. Independent replay
reproduces every count. These are candidate source counts, not measured
eligibility or independent biological replicates. Unresolved publication
ownership still prevents fresh-confirmation certification. Both numerical
acceptance requirements remain unchanged.

The approved DevLib download and design-only inspection are complete. The
39,678,063-byte archive passes the pinned gzip integrity check. All six
human-neuron RNA columns match their exact GEO owners. The reader explicitly
handles the deposited lowercase `eid` field; ambiguous aliases are rejected.
All 4,685 SNP IDs match original Ref/Alt ownership without allele, base or
chromosome conflicts. Of the original 2,559 pairs, 2,274 have both metadata IDs,
80 have Alt only, 57 have Ref only and 148 have neither. Every original pair
remains in the ledger; these counts do not establish measured eligibility.
The 222 control designs remain separate.

All raw SNP boundaries are one greater than the qualified source boundaries.
The observed offsets are retained without inferring a coordinate convention or
changing the source graph. The count file's `reference_name` supplies rsID
annotations, whereas exact element IDs join to the original coordinate names.
The count metadata does not resolve the 169 synthesis-sensitive sequence pairs.
Root and independent replay agree on every metadata record. Ten focused tests
pass, and three deliberately broken identity implementations are rejected.
The schema adaptation also passes independent fault checks.

The design-only phase established the inputs needed to freeze source admission,
the split, measurement and dependence sensitivity before TRAIN analysis. Its
evidence is preserved in
`rummel-source-admission/COUNT_DESIGN_ASSESSMENT.md` and
`rummel-source-admission/COUNT_DESIGN_CHECKPOINT.json`. Earlier preparation
records and their original file versions remain preserved.

This is a separate KCl-response target in mature neurons. It cannot be pooled
as an interchangeable Wnt-response label. The subsequent phase below fixes
source ownership, genomic families, exclusions and the split, then measures
shared-plasmid normalization and TRAIN repeatability. The Matoba paired-well
reducer is not applicable by analogy.
ValLib retains its result-selected ascertainment and is not fresh confirmation.
The original design-only phase is complete. The subsequent fixed TRAIN
measurement and diagnostic results are recorded below. Validation/test expression
values and the result workbook remain unopened; no new predictor fits were added.

### Rummel TRAIN measurement and diagnosis

The source-only partition retains all 2,559 original pairs: 1,209 TRAIN in 986
families, 416 validation in 342 families, 421 test in 352 families and 513 excluded
in 379 families. Exclusions propagate through the complete sequence/interval
graph, including unresolved synthesis changes and prior exposure. Availability
of count rows did not determine split membership. This is exploratory within-assay
development; unresolved broader locus and publication exposure prevent calling
it fresh external confirmation.

The prospective estimator sums each allele's barcode counts, then computes log2
RNA/DNA allele effects. The KCl-minus-baseline interaction cancels the shared
plasmid term algebraically. Its units and estimator differ from the author's
published MPRAnalyze coefficient. Every original TRAIN row is retained, including
261 unmeasurable pairs. The complete-case population contains 948 pairs in 800
families. An independent raw-count replay reproduces all 2,418 allele aggregates,
ownership and rejection records exactly. Validation/test numerical cells were
skipped before parsing.

The archive's `variantMapped` field contains literal Ref/Alt labels despite its
GEO description as Boolean. The first reader stopped before count parsing. The
corrected reader requires equality to both the pinned allele and element suffix;
all 429,217 retained TRAIN barcode rows satisfy this condition. Both attempts,
their protocols and checks are preserved. This correction changed no estimator,
source role or outcome-based eligibility rule.

| Fixed TRAIN diagnostic | Result |
|---|---:|
| Interaction correlation between one-sample and two-sample RNA views | -0.017139 |
| Interaction zero-effect loss on the full target | 0.359968 log2 squared |
| Mean zero-effect loss across the two diagnostic views | 0.883472 log2 squared |
| Copying one view to predict the other: loss | 1.796600 log2 squared |
| Baseline / KCl condition-effect view correlation | 0.497909 / 0.425463 |

These losses describe measurements, not sequence-model accuracy. Copying a view
loses to zero in both directions. Omitting each chromosome leaves interaction
correlation between -0.065879 and 0.017018. Condition-effect agreement includes
shared-plasmid error and is not independent biological replication.

A declared exploratory follow-up retains the full population and original family
weights. Interaction disagreement increases as RNA counts decrease, but its mean
is 8.28 times an independent-Poisson counting proxy. All four fixed count-scale
quarters have negative interaction correlation; the highest-coverage quarter is
-0.048499. All nine ways to divide the three samples per condition into one
versus two have correlations from -0.113492 to 0.082041. These checks support
measurement/data diagnosis before fitting this interaction. They identify
variation beyond simple counting but cannot separate its biological and technical
causes. No favorable count cutoff, sample split or replacement target was selected.

Evidence is in `rummel-measurement/PLAN.md`, `PRE_RUN.json`,
`train-measurements/`, `aggregate-replay/` and `COUNT_DIAGNOSIS.json` under the
reviewed build directory. The preceding roadmap is preserved as
`rummel-measurement/PRE_MEASUREMENT_ROADMAP.md`.

The subsequent matched barcode diagnostic is complete. SHA256 parity assigns
barcodes to two fixed halves; each view uses distinct RNA samples and barcode
counts. Its shared-versus-separate DNA comparison keeps the RNA numerators,
eligible pairs and original weights identical. All 1,209 original rows remain;
825 pairs in 704 families satisfy this diagnostic's support requirements,
retaining 86.8646% of the original eligible weight. This does not replace the
948-pair primary measurement population.

| Common allele-effect view correlation, same 825 pairs | Result |
|---|---:|
| Original full-barcode support | 0.671273 |
| Split barcode RNA with pooled DNA | 0.089945 |
| Same split barcode RNA with separate DNA | 0.226116 |

Separate DNA improves the matched comparison. These results therefore do not
attribute the stronger original full-support agreement to shared-DNA reuse
alone. The separate-DNA common correlation remains 0.192401–0.244275 across
chromosome omissions. This is modest conditional reporter agreement within the
same library/preparation, not independent biological replication. The matched
interaction remains exactly unchanged by DNA normalization (correlation
-0.082899), and all its chromosome-omission correlations are negative. The
fixed aggregate interaction route closes at zero new fits.

`BARCODE_DIAGNOSIS_PLAN.md`, `BARCODE_PRE_RUN.json`, `barcode-diagnosis/` and
`barcode-numeric-replay/` retain the follow-up. Independent arithmetic replay
matches every row, weight and chromosome omission within 3.11e-15. A separately
declared baseline-neuron allele-effect feasibility experiment can address the
roadmap's original edit-prioritization use case. It must retain all 948 original
eligible TRAIN pairs, keep the unsuccessful interaction route intact and pass
both effect-error and signed-selection comparisons before advancing.

### Rummel baseline-neuron sequence feasibility

The separate baseline-effect experiment retains the full 948 eligible pairs
in 800 families and all 261 original measurement rejections. Its five folds
derive from all 1,209 original TRAIN source pairs before measurement eligibility.
Raw canonical-RC delta features and ridge penalty 0.1 are fixed in advance.
Five k2 and five k5 fits use complement targets and complement family weights
only. Each excluded fold contributes two signed selections from distinct
families; no pooled ranking of OOF predictions is used.

| Fixed TRAIN arm | Family-weighted OOF MSE, log2 squared | Signed selection utility, log2 |
|---|---:|---:|
| Zero effect | 0.545893481102489 | 0 |
| k2 ridge | 0.5459176190298981 | -0.1324492889630082 |
| k5 ridge | 0.5456427270552119 | 0.4291458089483617 |

K5 strictly beats both comparators on both predefined point-estimate endpoints.
The MSE reduction against zero is only **0.0459346%**; this is relative loss
improvement, not prediction accuracy. Error improves over zero in two of five
folds and signed utility is positive in four. Independent replay reproduces all
predictions, exact scores and decisions without refitting; root reruns that
checker successfully. The run completes exactly ten starts on one CPU thread
with no GPU. Validation/test outcomes remain unopened.

This permitted the next declared TRAIN comparison. K5 stayed fixed while k2 was
independently tuned by minimum exact pooled MSE over the existing penalty grid.
All fifteen new starts completed; the five completed 0.1 fits were reused.
The selected k2 penalty is 1.0, with MSE **0.5441471964408893** and utility
**-0.17583489422728787**. Fixed k5 now loses on MSE by 0.274839% while retaining
its selection advantage, so the stronger joint progression rule fails. Independent
replay reproduces all fits, reuses, exact scores and decisions; no candidate is
nominated and no held-out outcomes are opened.

A declared no-fit diagnosis decomposes k5's excess MSE of 0.0014955306 into
extra prediction energy of 0.0120581392 minus additional target alignment of
0.0105626086. Chromosome 6 contributes +0.0051663086 to the total gap; only its
deletion reverses the saved-prediction comparison. Independent replay and the
root rerun reproduce all 170 exact diagnostic values. No chromosome is dropped or
reweighted. Independent source-only grouping finds 22 chromosome blocks, with
chromosome 6 holding 37–39% of each role's source families. Held-out chromosomes
also occur in TRAIN; subsequent uncertainty must respect this concentration and
retain the family-held-out scope.

This work assesses a possible useful sequence baseline for the assay. K5 remains
the prior k-mer comparator for later candidates; the result establishes no
improvement over that same mechanism. The uncertainty-aware acceptance rules
below remain in force. Evidence is in `rummel-baseline-pilot/RESULT.md`,
`FREEZE.json`, `experiment/complete/`, `result-review/` and `root-replay/` under
the reviewed build directory. `PRE_PILOT_ROADMAP.md` preserves the preceding
roadmap. The stronger comparison and its diagnosis are retained in
`rummel-baseline-tuning/RESULT.md`, `experiment/complete/`, `result-review/`,
`ERROR_DIAGNOSIS.json` and `chromosome-block-design/`.

### Fixed short-context anchor and matched control

The next fixed TRAIN comparison reuses each fold's exact tuned k2 lambda-1
coefficients. It fits direct k5 and a k5 correction to the complement-only k2
residual, both at lambda 1. All ten new starts complete; five k2 fits are reused.
The original 948 pairs, 800 families, 261 rejections and source folds remain.

| New model | Family-weighted excluded-fold MSE | Signed selection utility |
|---|---:|---:|
| Direct k5, lambda 1 | 0.5441923391392202 | 0.5065788215037964 |
| k2 anchor + k5 residual, lambda 1 | 0.5429088350809197 | 0.3585836973244997 |

The anchor improves MSE against all four required comparators: zero, tuned k2,
prior k5 and matched direct k5. It reduces MSE by 0.546745% relative to zero,
a modest relative loss reduction rather than an accuracy percentage. Selection
utility is below prior and direct k5, so the joint rule fails. The direct control
is not automatically promoted. This fixed recipe closes without nomination;
held-out outcomes remain unopened. Independent source qualification and root
synthetic tests pass. Independent numerical replay checks all reused/new normal
equations, predictions, scores and selected edits without refitting; root
reproduces every non-resource field. No existing comparator score drifts.

A source-only follow-up confirms a specific representational limit: unchanged
distant sequence cancels in raw k2/k5 deltas. The existing ROUND3 nonlinear
whole-sequence motif scorer distinguishes a synthetic pair of edits with the
same local context and different distant motifs. Full-sequence motif attention
also already exists. Reusing these on DevLib would test transfer to this assay;
it would not establish architectural novelty. No new biological fit follows
from that synthetic distinction. The current result and exact evidence are in
`rummel-anchored-pilot/RESULT.md`, `FREEZE.json`, `experiment/complete/`,
`source-review/`, `root-source-review/`, `result-review/`, `root-replay/` and
`CONTEXT_REPRESENTATION_CHECK.json`.

Source reuse is now checked separately in `rummel-motif-source/ADMISSION.json`
and `ADMISSION.md`: the exact 532-profile list and both matrix-source identities
match the existing admission. All profile widths are 5–19 bases and fit DevLib's
160-base inserts. Eligibility rests on declared binding-construction methods;
original pipelines and curation decisions are not exhaustively reconstructed.
Only fixed profiles/scoring helpers transfer, with no old reporter-trained
readout or 270-base outcome loader. This source qualification adds zero fits and
no downloads. The following comparison uses these admitted sources.

### Fixed full-insert, edit-overlap and permuted motifs

The frozen fifteen-fit comparison retains all 948 eligible pairs in 800 families,
all 261 rejections and the original source-family folds. Each arm fits ridge at
lambda 1 on raw k5 deltas plus 532 motif deltas. Each feature block's scale is
computed separately from its fitting complement's equal-family energy. The
same scale is applied to the excluded fold; targets and penalties stay fixed.
Real full-insert motifs are the sole candidate. The two controls use only windows
overlapping the edit, or full inserts with the inherited fixed column permutation.

| Model | Family-weighted excluded-fold MSE | Signed selection utility |
|---|---:|---:|
| Zero effect | 0.5458934811024890 | 0 |
| Tuned k2, lambda 1 | 0.5441471964408893 | -0.1758348942272879 |
| Prior k5, lambda 0.1 | 0.5456427270552119 | 0.4291458089483617 |
| Direct k5, lambda 1 | 0.5441923391392202 | 0.5065788215037964 |
| k2 anchor + k5 residual | 0.5429088350809197 | 0.3585836973244997 |
| Full-insert motifs + k5 | 0.5469254262339630 | 0.0841159425382546 |
| Edit-overlap motifs + k5 | 0.5445503752013403 | -0.0854566345787301 |
| Permuted motifs + k5 | 0.5426603027598709 | 0.3873616652724990 |

Full-insert motifs fail MSE against every comparator. Their utility exceeds zero,
k2 and the overlap control, but loses to the remaining four comparators. The
strict seven-comparator joint rule fails. Permuted motifs are an ineligible
control and also lose utility to direct k5. No control is promoted and this fixed
recipe closes without nomination, retries or a parameter sweep.

The result gives no predictive reason to prefer these real motifs in this assay.
It does not establish that distant context is useless: full-versus-overlap tests
this scorer's unchanged-window mass, while one permutation is a limited feature
control. Equal-energy scaling does not equate covariance or effective ridge
capacity. The repeated TRAIN-family screen supplies development evidence only.

Source qualification detected and repaired omitted-cache and omitted-source-pin
guards before biological fitting. Independent feature checks cover 255 rational
synthetic cases and 3,192 admitted-profile comparisons; four feature mutants fail.
Root tests cover a rational ridge solution, complement isolation, every strict
comparator endpoint, cache identity, failure accounting and output preservation;
four deliberate numerical/gate/fold faults fail. Fifteen starts and completions
are retained, adding no validation/test outcome access. Independent numerical
replay and the root rerun reproduce all eight exact score/selection reports and
all fifteen fitted solutions. Independent probability-product scoring reconstructs
all three feature arrays within the declared floating-point allowance; exact
motif identities, permutations and edit positions match. The checker has fourteen
passing synthetic tests and rejects five deliberate faults. Its initial setup
error looked for population files in the tuning cache; the corrected ownership
check is qualified and the failed attempt remains recorded. No biological fit
was repeated. The exact plan, source, inputs and outputs are under
`rummel-motif-pilot/`, including `FREEZE.json`, `experiment/complete/`,
`feature-review/`, `source-qualification-final/`, `root-source-review/`,
`result-review/` and `root-replay/`.

### Selection transfer between existing baseline measurement views

The next zero-fit diagnostic asks whether large measured effects repeat at the
ten-edit budget, even though the existing global separate-DNA correlation is
only 0.19091. It freezes both directions between the previously exported barcode/
RNA views. The 825-pair, 704-family availability mask is inherited, retaining
8339/9600 of the original pair-weight mass. All 1,209 source records, 948 primary
eligible pairs, 261 original exclusions and 123 additional unavailable records
stay in the ledger. The full-barcode model target and training population remain.

Within each original fold, select two distinct families by absolute discovering-
view effect, with lexical source-ID ties, then score the discovering sign times
the other view's effect on those exact edits. The random reference selects
available families uniformly and an available pair uniformly within each family,
using that pair's discovering sign. Exact expected utility and observed utility
are averaged equally across the five folds. This is a conditional selection
reference, separate from original family weights used for coverage.

| Direction | Transferred signed utility | Random expectation | Difference |
|---|---:|---:|---:|
| View A selects; view B measures | 0.4184015661 | 0.1257793756 | 0.2926221904 |
| View B selects; view A measures | 0.8046773388 | 0.1413481176 | 0.6633292212 |

Both directions exceed zero and their fixed random reference. All 22 declared
chromosome omissions retain both comparisons; no omission is selected. This is
observable repeatability of the measured selection tail on the existing matched
support. It supports investigating useful edit prioritization in this assay.
It does not show that DNA alone predicts those edits, estimate a biological noise
ceiling, supply confidence bounds or nominate a model. The views share library,
cultures and processing and have unequal one-versus-two RNA precision. Their
utilities are not directly comparable to the models' full-barcode utilities.

Four root tests and independent exact-arithmetic source qualification pass.
Four injected ranking/family/random-reference/gate faults are rejected after
adding a one-direction-only counterexample; the initially nondiscriminating
fault check is retained. The missing digest-owner pin was repaired before views
were inspected. The audit adds zero fits, downloads or held-out outcome access;
the cumulative fit ledger remains 365. Independent saved-result replay and the
root rerun agree on the complete source ledger, every exact value, all selected
edits and all 22 chromosome omissions. All 77 frozen source/evidence pins remain
unchanged. The next predictor decision can now target the repeatable selection
tail while retaining the original full-barcode error and utility requirements;
this diagnostic supplies no new model nomination or eligible subset.
Evidence is under `rummel-selection-repeatability/`: `PLAN.md`, `FREEZE.json`,
`audit/complete/`, `TEST.log`, `MUTATIONS.json`, `final-review/` and
`ROOT_REPLAY.json`.

### Fixed model selections and predicted directions

A zero-fit crossover keeps the direct-k5 and anchored models' original selections
over all 948 primary rows: two distinct families per original fold. It scores
each model's signs on both lists. The full-barcode target gives this exact table:

| Selected list | Direct model signs | Anchored model signs |
|---|---:|---:|
| Direct k5 | 0.5065788215037964 | 0.5065788215037964 |
| Anchored k2 + k5 residual | 0.3585836973244997 | 0.3585836973244997 |

Both fixed-sign membership contrasts are -0.14799512417929678; both fixed-list
sign contrasts are exactly zero. Thus changed selections explain the observed
utility loss. This is a mechanical decomposition of saved predictions, not an
identified biological cause or proof that a rank-only loss is better.

The direct model's original ten picks have both existing measurement views.
Their signed utility is 1.1961566121 versus a model-sign random expectation of
0.0541504951 in A, and 0.1524460024 versus 0.0128678675 in B. The reference samples
available families uniformly and a pair uniformly inside each family, using
the named DNA model's sign. It differs from the prior measured-view selector's
sign reference. These comparisons support continuing edit-prioritization research.
They are descriptive TRAIN results: the views share library, cultures and
processing, overlap the original training target, and have unequal precision.
Their utilities must not be compared directly to full-barcode model utilities.

Three anchored picks have unavailable views, so both complete anchored-view
endpoints and view crossover decompositions remain unavailable. No edit is
replaced or omitted from the selection budget. All 1,209 source records,
948 primary pairs, 800 families and original folds remain represented; the
matched-view support stays 825 pairs/704 families. No model is newly nominated.

Three root tests and independent exact hand checks pass; four deliberate faults
are caught. Independent replay and the root rerun reproduce every cell, random
reference, selected row, missing-value decision and contrast without fitting.
The 32 frozen pins remain unchanged. Evidence is in
`rummel-selection-diagnosis/{PLAN.md,FREEZE.json,RESULT.md,RESULTS.json,source-review/,ROOT_REPLAY.json}`.
The accompanying primary-source learning review retains supervised effect
regression: the ideal conditional mean optimizes both current endpoints.
Reliability weighting remains conditional on qualified assay-specific precision
estimates; guessed weights or a blind return to the failed prior precision trial
are not a next experiment. Cumulative biological fit starts remain 365.

### Joint short- and longer-word regression

One frozen candidate jointly fits the existing raw k2 and k5 feature blocks at
lambda 1, using the original targets, 948 eligible pairs, 800 families and five
source-family folds. Each complement supplies its own equal-family weights.
Eight saved comparator reports are reused without refitting. Joint fitting
changes coefficient allocation and regularization; it adds no sequence context.

| Model | Family-weighted excluded-fold MSE | Signed selection utility |
|---|---:|---:|
| Direct k5 | 0.5441923391392202 | 0.5065788215037964 |
| k2 anchor + k5 residual | 0.5429088350809197 | 0.3585836973244997 |
| Joint k2/k5 | 0.5429573411680731 | 0.3585836973244997 |

Joint fitting reduces MSE by 0.22694% relative to direct k5 but selects exactly
the same ten edits, with the same signs, as the anchored model. It loses MSE to
the anchor and permuted-motif control, and loses utility to direct k5, prior k5
and permuted motifs. Tying anchored utility also fails the strict rule. The
recipe closes without nomination, retries, a sweep or held-out outcome access.

All five starts complete. Four implementation tests and independent source
qualification pass; five deliberately broken implementations fail. The independent
checker has thirteen passing synthetic tests and detects three deliberate faults.
Independent replay and the root rerun verify all five saved joint solutions,
predictions, nine exact reports and sixteen comparisons without refitting.
Review caught inherited closed-comparator advancement fields in the runner and
checker before biological execution; regression checks preserve the repairs.
Evidence is under `rummel-joint-pilot/`, including `PLAN.md`, `FREEZE.json`,
`RESULT.md`, `experiment/complete/`, `source-review/`, `result-review/` and
`root-replay/`. The 124 pre-fit pins remain unchanged. Cumulative starts are 370.

The next direction is source-only qualification of the same study's eQTLLib
baseline measurements. The frozen design map and independent reconstruction agree
on all 4,325 eQTLLib plus 2,559 DevLib rows. Exactly 1,542 eQTLLib pairs overlap
DevLib by whole-allele DNA identity. Another 2,783 components have no DevLib owner
at this graph scope; genomic overlap and exposure may add dependencies. Original
DevLib roles remain fixed. Raw sequence order reverses declared allele order in
2,243 rows. The malformed source row 867 remains unresolved.

Independent literal reference matching and root replay now agree on every one of
the 4,325 source rows. All 4,324 valid pairs have a unique local GRCh37 match in
column F, forward strand, at the named coordinate. The existing cache supplied
313 rows; 70 bounded requests supplied the remaining envelopes in 986,915 response
bytes. A timeout after 66 batches and the explicit four-batch continuation remain
in the evidence. Row 867 has no guessed SNP envelope.

The augmented graph has 4,296 components over all 6,884 designs. Independent
reconstruction and root replay preserve all original owners and match every
output. Genomic overlap reveals a 13-row eQTLLib component linking a TRAIN family
to a test family: `dev:2254 -> eqtl:3828 -> dev:549`, with 66- and 56-base overlaps.
The original DevLib windows are separated by 38 bases. Keep their roles unchanged
and exclude new supervision from this component. Another eight-row eQTLLib group
touches a synthesis-excluded DevLib family and TRAIN.

The unchanged four-study prior audit covers 4,324 qualified pairs and identifies
392 with prior-source overlap. Independent replay agrees on all 1,445 physical
owner hits and every cross-assembly projection. The malformed isolated row 867
remains explicitly excluded, with incomplete physical prior screening.
Propagating available exposure, synthesis, mapping and owner exclusions retains
3,019 source candidates in 2,503 components: 2,279 rows in 1,929 components with
no DevLib owner, and 740 rows in 574 original-TRAIN components. These are not yet
admitted measurements. Additional publication leads, construct identity and
count eligibility remain open. The separate 28,368,772-byte eQTLLib archive was
downloaded after explicit approval under `rummel-eqtl-admission/COUNT_ACCESS_PLAN.md`.
Its design-only inspection is complete; numerical count analysis remains closed.
The publication follow-up maps the already recorded rs3764512 example to
GRCh37 chr18:52,495,324 and finds no overlap with current inserts. The two known
LRP1 identifiers were already included in the exclusion seeds. These checks leave
3,019 source candidates. A GATA1/OL43 benchmark context is identified, but the
original snippet provenance and exact construct membership remain unresolved.
An incidental EML1/OCR caption lacks an assembly in the retrieved excerpt and
remains an explicit exposure lead. `publication-followup/RESULT.md` records the
bounded checks; no gene name is used to invent a variant coordinate.
Independent literal reconstruction and the root rerun reproduce the point join,
all original owners and unchanged candidate counts. A positive-hit control catches
deliberately omitted seed propagation, which the real zero-hit case cannot test.
The subsequent source-only inner-fold audit finds one candidate component with
11 eQTLLib rows connecting original development folds 1 and 2. Preserving fold
membership requires whole-component fitting exclusions: fold 1 additionally
purges DevLib 2022; fold 2 additionally purges DevLib 2009, 2382 and 2383. Both
omit the 11 connected added rows. Future baseline arms need identical purges;
closed trial scores cannot stand in for those changed fits. The existing outer
TRAIN/test bridge also needs explicit treatment of its retained TRAIN owners
before a component-disjoint outer evaluation. `fold-admission/RESULT.md` and
`FOLD_GUARDS.json` retain the complete design-only constraints. No new fitting
protocol, role reassignment or measurement admission follows from this audit.
Independent reconstruction agrees with all 2,913 guard entries and the report;
the root rerun reproduces the independent audit. Both deliberately broken fold
guards fail the focused test.

The approved count archive passes its recorded length, SHA256 and complete gzip
integrity checks. The unchanged whitelist reader retains 8,133 distinct metadata
IDs and exposes no numerical counts or barcode support. The complete source join
preserves all 4,325 designs: 3,817 have both IDs, 121 only Ref, 153 only Alt and
234 neither. All 7,908 matched IDs resolve to their physical DNA columns:
3,938 Ref IDs to F and 3,970 Alt IDs to G, despite reversed workbook declarations.
There are no conflicting ID metadata or physical-base mismatches. All matched
intervals have a (+1,+1) offset from verified zero-based half-open placements;
the source convention is recorded as unresolved, without changing placements.
The 225 unmatched IDs remain explicit: 214 typed controls and eleven SNP-tagged
unknowns. Malformed source row 867 remains unresolved and has neither allele ID.

Among the unchanged 3,019 source candidates, 2,711 have both IDs, 81 only Ref,
88 only Alt and 139 neither. Identifier presence does not establish read-depth
eligibility. The three deposited human baseline columns map to GSM7056974–76;
the eleven mouse columns stay unselected. The `count` column is designated for
original-plasmid normalization, not an RNA replicate. Shared replicate names across
libraries do not establish paired preparations or independent donors.

Independent reconstruction and root replay reproduce the safe projection and
every source join record. Sample ownership is checked against retained GEO lines;
the physical-allele test passes and catches both deliberately broken mappings.
`COUNT_RESULT.md`, `COUNT_IDENTITY_JOIN.json`, `COUNT_DESIGN_SUMMARY.json` and the
two count-review directories retain the evidence. Full synthesized-construct
identity, remaining exposure decisions, count QC and the TRAIN-only measurement
contract remain open. The augmented inner-fold and outer-owner exclusions must
be fixed before a new evaluation. No numerical eQTLLib count analysis or new fits
were performed.
Existing TRAIN Ref/Alt aggregates could supply activity supervision, but they
reuse the same measurements, and earlier auxiliary-training trials failed the
combined gates. No new training recipe follows from that availability alone.
`rummel-eqtl-admission/RESEARCH.md` records this decision, alternative sources
and incidental published-example exposure leads. `RESULT.md`, `design-map/`,
`reference-complete/`, `genomic-map/`, `prior-adapter/`, `source-exclusions/`
and their review directories retain the completed source-map evidence.
This direction adds the one approved eQTLLib count archive and zero fits.

### Additional baseline-neuron measurements

The subsequent intended-insert admission and fixed TRAIN measurement are complete.
Screening the documented 190-base synthesis constructions identifies three new
KpnI adapter-junction risks, at eQ source rows 1075, 1809 and 2559. Their singleton
components have no DevLib owners and are excluded before count access. This is
a design-derived cloning risk, not observed truncation. The initial source ledger
is preserved; final admission retains 3,016 rows in 2,500 augmented components.
The intended 160-base core and physical F/G mapping define this exploratory
reporter target. Deposited orientation proportions and narrative exposure leads
remain unresolved, preventing an independent external-confirmation claim.

Fixed count QC yields 2,210 eligible designs in 1,892 components, with all 806
measurement rejections retained. Of these, 308 lack one or both element IDs and
498 fail positive RNA support. Equal-component first-RNA versus other-two-mean
correlation is 0.2397759825; directional copy gains are -0.1871447392 and
-0.4524115501. Full-label zero-effect MSE is 0.3520286241 squared log2 units.
All 22 chromosome omissions retain modest correlations, 0.23383–0.25011.
Shared DNA and unequal RNA precision limit the interpretation; these values are
not an independent noise ceiling. No outcome-dependent exclusions were added.

Independent replay and the root rerun verify all 6,032 allele aggregates, every
measurement and rejection, all component contributions and chromosome omissions.
Integer values match exactly; the largest floating difference is 8.99e-15 over
38,403 fields. Seven root tests and independent synthetic/mutation checks qualify
the access boundary, physical allele sign, plasmid term and component weights.
All 72 pre-run pins remain unchanged. Evidence is in `rummel-eqtl-measurement/`,
including `ADMISSION.md`, `PLAN.md`, `PRE_RUN.json`, `train-measurements/`,
`RESULT.md`, the retained reviews and `ROOT_NUMERIC_REPLAY.json`.

The reviewed next comparison preserves original folds and quarantines Dev TRAIN
684, 2254 and 2381 from fitting and scoring because their augmented components
reach excluded or held-out owners. Whole-component inner-fold purges apply to
every arm. The source-only union is 4,222 observations in 2,908 components before
measurement QC. Repeated cross-library pairs remain dependent observations.
The next fixed-recipe comparison will measure the effect of adding these data
against newly fitted, matched Dev-only and sequence baselines on the same Dev
outcomes, preserving both endpoint requirements. This phase adds zero fits;
cumulative starts remain 370 and validation/test outcomes remain unopened.

### Matched training-data augmentation

One fixed comparison fits raw canonical-RC k2 and k5 ridge at lambda 1 with
Dev-only and Dev-plus-eQ supervision. All 20 declared starts complete once in
4.19 seconds with one CPU/BLAS thread and no GPU. The sole candidate is augmented
k5; no model selection, penalty sweep or control promotion occurs.

| Model | Component-weighted Dev OOF MSE | Signed selection utility |
|---|---:|---:|
| Dev-only k5 | 0.5440146661 | 0.5065788215 |
| Augmented k5 | 0.5457496034 | 0.4012182070 |
| Dev-only k2 | 0.5439540581 | -0.1758348942 |
| Augmented k2 | 0.5461453315 | 0.1049741352 |
| Zero effect | 0.5457358786 | 0 |

Added supervision increases candidate MSE by 0.3189% relative to matched Dev-only
training and reduces utility by 0.1053606145 log2 units. It also loses MSE to zero.
Utility exceeds both k2 arms, zero and the two fixed random references, but the
combined criterion fails. All 2,000 paired chromosome draws are evaluable and all
ten improvement lower bounds are nonpositive. This closes the fixed pooling
recipe and retains the stronger Dev-only model.

All 948 original eligible Dev designs remain evaluation targets, grouped into
798 augmented components. The 2,210 added observations create a 3,158-record
training pool with 2,235 components. The three quarantined source TRAIN rows
already lacked eligible measurements. Original folds are preserved, component
purges use all source owners before QC, and a bridging component is counted once
in pooled loss. Selection remains two distinct components per original fold;
bootstrap copies are tagged resampling occurrences, not duplicate real edits.
No held-out outcome is opened.

Ten root tests and independent population, scoring and runner qualification pass.
Every source row, guard and fitting weight agrees with independent reconstruction.
Review repairs are supported by deliberate fault witnesses. `rummel-augmentation-pilot/`
retains the frozen 202-file input contract, all starts, coefficients, predictions,
draws and scores. Independent replay verifies all 20 normal equations, every
prediction and fit event, exact rational point reports, all 2,000 resamples and
22 chromosome omissions. Both root reruns agree. The largest scoring arithmetic
difference is 1.08e-15 across 8,297 floating fields; this measures replay agreement.
Cumulative biological fit starts are 390. The next source-defined diagnostic is agreement
between identical physical edits across libraries; pooling failure alone does
not identify measurement bias, distribution shift or a biological cause.


### Cross-library agreement at identical physical edits

A source-only literal ordered 160 nt Ref/Alt match identifies 695 physical units
in 574 augmented components, each with one source observation per library.
Swapped and reverse-complement alternatives match zero additional eQ rows. The
membership and analysis rule are fixed before the numeric join. Existing QC
leaves 527 complete pairs in 454 components; 37 are Dev-only, 27 eQ-only and 104
neither. Every source row and rejection remains in the ledger.

Component-weighted Dev/eQ correlation is **-0.11255317**. Copying measurements
between libraries gives MSE **0.79909687**, compared with zero-effect MSE
**0.43915958** on Dev and **0.28101260** on eQ. Both directional copy gains are
negative. All 21 source-chromosome omissions preserve negative correlation
(-0.1381 to -0.0633) and both negative gains. These are complete-case agreement
diagnostics and sensitivity checks, not independent biological confirmation.

The finding supports investigating library compatibility before treating the
added labels as interchangeable supervision. Separate barcode/plasmid pools,
unresolved culture/batch pairing and unknown orientation proportions remain
possible explanations. Copying noisy labels can fail even when estimating a
shared conditional mean is useful; this diagnostic establishes neither a causal
measurement defect nor impossibility of sequence prediction. It changes no
training eligibility, target or model-selection rule.

Three root tests, four independent hand tests and a required-input boundary
witness pass; seven deliberate statistical/ledger faults are detected across
the two test sets. Root repeats the independent tests and reconstructs all source
matches. Independent matrix arithmetic reconstructs the complete measurement
ledger and verifies 5,312 numeric fields, including component contributions and
all omissions, within 8.88e-16. `rummel-shared-pairs/` retains the fixed 44-file
contract, complete result and checks. This phase adds zero fits, downloads or
raw-count reads; cumulative biological fits remain **390**. The next narrow
question is whether RNA/plasmid terms account for the disagreement, under a new
fixed diagnostic on these same source-defined pairs.


### RNA/plasmid accounting on the shared edits

A fixed zero-fit diagnostic keeps the same 527 complete pairs in 454 components
and preserves all 695 source units. Exact sample/header bindings and 2,108 allele
IDs reassert; each normalized sample effect and original full label reconstructs
exactly from the sealed aggregates. No sign, sample or denominator is changed.

Weighted cross-library mismatch decomposes as **2.25333470** RNA-ratio mismatch
plus **1.61089629** plasmid-ratio mismatch minus **3.06513412** from their cross
term, giving the previous **0.79909687** full-effect mismatch. The four signed
covariance contributions are -0.07589098 (RNA/RNA), +0.04381258 and +0.01126957
(the two mixed terms), and -0.01867711 (plasmid/plasmid), totaling -0.03948593.
The existing normalization reduces mismatch relative to raw RNA ratios; correlated
terms do not identify causal fractions of noise or justify dropping normalization.

Eight of nine DNA-normalized cross-library sample correlations are negative
(-0.14298 to -0.00538); the ninth is +0.02979. Within-library correlations are
higher but share a plasmid denominator and do not establish independent
reliability. Every chromosome omission preserves negative raw-RNA and full-effect
covariance and lower mismatch after normalization. Individual cross-sample signs
vary, with five to nine negative correlations per omission. The result points
away from an obvious single-sample repair without excluding a library-wide
source-method problem.

Three root and five independent tests pass, as do six required-input boundary
cases and a source-binding review. Root repeats the independent checks. An
independent 60-digit Decimal/population-covariance replay verifies 7,514 fields
within 2.95e-15. `rummel-pair-decomposition/` retains the fixed 46-file contract,
complete ledger, six-sample matrices, all 21 omissions and witnesses. No new
fits, raw-archive reads, downloads or held-out outcomes; cumulative fits remain
**390**. The next bounded task compares the author's barcode QC, normalization
and estimation with this fixed aggregate-ratio target before a corrective model
or new pooling recipe is considered.

### Author count-model preparation

The source-method audit identifies a specific next comparison. Rummel's methods
collapse barcodes into allele-specific pEE counts, fit each library separately
with RNA replicate+status versus replicate, and retain DevLib measurements for
duplicate SNPs in a merged report. Our closed augmentation experiment tested a
different supervision hypothesis; the method difference is not a retroactive
implementation bug. The existing means of log ratios also differ from the
paper's count-likelihood estimator. [Cell methods](https://doi.org/10.1016/j.cell.2023.09.015).

The existing sealed allele aggregates now supply separate-library count matrices
without another assay download. Preparation retains the same eligible designs,
all rejected records and the component quarantine. Each pair has two original
plasmid allele observations and six baseline RNA observations, with explicit
allele/accession annotations. No shared DNA count becomes an independent DNA
replicate. Existing QC and labels remain intact. All 4,225 TRAIN records remain
in the ledger; 25,264 exported integer values for 948 Dev and 2,210 eQ observations
match exactly on independent readback. A Dev metadata-name mismatch stopped the
first export before publication; a repaired adapter and regression witness now
match the actual Dev/eQ schemas.

Exact historical reproduction remains unresolved: the author invocation, package
version, abundance-cutoff inequality and normalization universe were not recovered.
A prospective estimator must fix these choices explicitly and qualify its DNA/RNA
design, coefficient sign, units and convergence on synthetic cases using the actual
package. Current MPRAnalyze 1.30.0 source drops the allele DNA term from the proposed
two-observation `~status` design. Retaining both singleton DNA means instead creates
a saturated Gamma likelihood with no finite shape maximum. That proposed classic
setup closes before biological fitting. An isolated R runtime request is deferred;
the smaller diagnostic below uses the existing Python runtime. See the local
`build/neural-predictor-reviewed-20261005/rummel-count-model/` protocol and audits.

Method searches returned incidental published narrative. The named designs were
already present in the exclusion ledger; computed validation/test counts remain
closed. The retained source audit records that exposure explicitly.

### Fixed count-weighting diagnostic

One prospective contrast conditions on each RNA sample's allele total and treats
the original plasmid counts as fixed exposures. Its exact estimate is the log2
ratio of summed Alt/Ref RNA counts minus the plasmid log2 ratio. An independent
derivation and synthetic unequal-DNA, allele-swap and low-support sign controls
qualify this formula. It changes weighting, with no change to eligibility or
physical matching; it is not MPRAnalyze or an estimate of DNA uncertainty.

On the same 527 complete pairs in 454 components, correlation changes from
-0.11255317 to -0.08231955. Copy MSE falls from 0.79909687 to 0.74521677, a 6.7426%
reduction, but both copy gains remain negative (-0.31712427 toward Dev and
-0.48273905 toward eQ). All 21 chromosome omissions preserve negative correlation
and both negative gains. The fixed pooling-repair rule fails; no target or
predictor is replaced. Lower mismatch by itself would not pass, because the
alternative estimator can change scale.

The complete 695-unit ledger is preserved. Independent exact rational ratios,
60-digit Decimal logs and matrix statistics verify 7,548 numeric fields within
1.655e-14; the root reruns the independent synthetic driver. The phase computes
1,118 closed-form measurement estimates, with zero new predictor fits, raw-archive
reads, downloads, heldout outcomes or GPU work. Details and the fixed rule are in
`build/neural-predictor-reviewed-20261005/rummel-count-weighting/`.

### Fixed DevLib precision-weighting trial

One new candidate estimates an affine RNA variance trend inside each purged
training complement, then fits the existing k5 ridge model with normalized
inverse-variance weights. The original mean-log labels, 948 observations,
798 augmented genomic components, folds and evaluation weights remain fixed.
All five saved predictor controls and the same chromosome draws are reused.
The working variance trend is motivated by empirical MPRA precision weighting;
it does not estimate shared DNA error or reproduce mpralm's smoother.

Signed selection utility rises from 0.50657882 to 0.84756018, a 67.31% relative
score increase for the fixed ten-edit budget. MSE changes from 0.54401467 to
0.54408508, 0.01294% worse, and remains above the Dev k2 value 0.54395406.
Ten of twelve point gates pass. Four paired interval lower bounds are positive;
the utility improvement interval against Dev k5 spans [-0.25270804, 0.64449673].
All 2,000 chromosome draws are evaluable. The utility point gain survives 21 of
22 chromosome omissions, while its MSE gain survives only five. The fixed
recipe closes without advancement or further tuning.

Five sequence fits and five nuisance variance calibrations complete once in
4.49 seconds on one CPU thread. Independent Decimal variance calculations and
augmented least-squares replay reproduce weights and predictions, with prediction
differences below 1.354e-16. Independent saved-score replay and the root rerun
agree on all reports, gates, draws and omissions; 9,956 floating fields agree
within 1.087e-15. Controls remain unchanged. No GPU, raw archives, additional
downloads or validation/test outcomes are used.

The plan, result, fixed predictions and independent receipts are in
`build/neural-predictor-reviewed-20261005/rummel-dev-precision/`. The follow-up
below checks the new picks against existing separate-barcode/RNA views,
preserving all selections and missing observations.

### Precision-model selections across existing measurement views

The fixed zero-fit follow-up confirms a positive selection difference in both
existing views. View A utility rises from 1.19615661 to 1.62846156; view B rises
from 0.15244600 to 0.51536151. Both models exceed their respective conditional
model-sign random references in both views. Every selected edit has both views;
no missing edit is replaced or removed from the ten-edit budget.

Two substitutions account for the complete gain: `chr12_2314999` to
`chr6_28583394` in fold 2 and `chr9_7172500` to `chr5_109046212` in fold 4.
Each improves the signed effect on the full target and both views. Predicted
directions are identical on both fixed lists, so the exact crossover assigns
the gain entirely to selection membership.

The audit preserves the original 1,209-source ledger and its 261 exclusions.
It joins 825 available views to the current 948-observation, 798-component model
population, with 123 unavailable views retained explicitly. Available current
component weight is 8321/9576; the historical 800-family source weighting remains
separately reported. Independent replay and the root rerun reproduce every value
exactly, including pooled MSE, all twelve crossover cells and both paths.

These TRAIN views share library, cultures, processing and source measurements;
they support the observed prioritization gain without providing independent
biological replication. This audit changes neither the prior error result nor
its paired uncertainty gates. Precision k5 remains closed without advancement.
No new fits or heldout outcomes are used; the cumulative predictor count stays
395. The frozen plan, results and replay receipts are in
`build/neural-predictor-reviewed-20261005/rummel-precision-views/`.

### Nested combination of the error and selection models

One new candidate learns a convex combination of the matched Dev k2 and
precision-k5 predictions. Each outer fold's weight uses inner excluded-fold
predictions only, with complete pre-QC genomic ownership excluding both folds
from each learner and its RNA-variance calibration. Component-weighted MSE
determines one clipped analytic weight; no mixture grid or ranking adjustment
is introduced. Saved outer predictions and all six controls are reused.

The outer precision fractions are 1, 0.96609037, 0.98643489, 0 and 0.26208259.
Their resulting MSE is 0.54542267, versus 0.54395406 for Dev k2 and 0.54408508
for precision k5. Signed selection utility is 0.50184668, versus 0.50657882 for
Dev k5 and 0.84756018 for precision k5. Nine of fourteen point comparisons pass;
none of the fourteen paired interval lower bounds is positive. All 2,000 draws
are evaluable. All 22 chromosome omissions preserve the MSE loss to Dev k2 and
utility loss to precision k5. The fixed combination closes without advancement.

Forty inner sequence fits, twenty variance calibrations and five analytic
combination estimates finish once in 4.98 seconds on one CPU thread. Independent
replay and the root rerun verify all fitted equations, count-derived weights,
five combinations, source ownership, 130 fit events and complete scoring reports.
All exact point fields and decisions agree; 11,615 floating score fields agree
within 1.371e-15. This establishes execution consistency. The preceding precision
model's measured selection gain remains intact, while this attempted combination
fails to carry it into a better joint predictor.

No raw archives, additional downloads, GPU work or heldout outcomes are used.
Cumulative predictor starts are 435, with 25 variance calibrations and five
combination estimates separately recorded. Evidence is under
`build/neural-predictor-reviewed-20261005/rummel-nested-stack/`.

### Whole-insert composition and the constant-context control

The fixed ten-fit comparison adds 55 quadratic canonical-dinucleotide contrasts
to raw k5 edit deltas, with off-diagonal sqrt(2) scaling and normalization by the
159 source-defined dinucleotide windows. Its control replaces each row's count
midpoint with the fitting complement's precision-weighted mean. Both use the
same saved precision weights, lambda 1, 948 rows, 798 components and complete
pre-QC ownership purges. No nuisance model is refitted.

| Model | Component-weighted MSE | Signed selection utility |
|---|---:|---:|
| Varying composition | 0.5432340389 | 0.8784141358 |
| Constant composition, ineligible control | 0.5433145027 | 1.1869765263 |
| Previous precision k5 | 0.5440850847 | 0.8475601833 |
| Dev k5 | 0.5440146661 | 0.5065788215 |
| Dev k2 | 0.5439540581 | -0.1758348942 |
| Zero | 0.5457358786 | 0 |

The sole candidate improves both point metrics over all six prior controls.
Relative to precision k5, MSE falls 0.1564% and utility rises 3.6403%. Its MSE
improvement interval [-0.000029333, 0.001940829] and utility improvement interval
[-0.425327, 0.450478] include zero. Its one failed point comparison is utility
against the constant-context control: the loss is 0.30856239. All 22 chromosome
omissions retain that loss, alongside a small MSE advantage over the control.
Fifteen of sixteen point gates and four of sixteen interval gates pass; all
2,000 draws are evaluable. The complete rule fails and the recipe closes without
nomination or tuning. No control is promoted after observing its result.

The constant-context control has 40.0463% higher utility and 0.1416% lower MSE
than precision k5. A post-result, zero-fit algebraic witness reduces its readout
to fold-specific k2 coefficients beside k5: all 948 binary-coefficient rational
identities hold, rounded predictions differ by at most 1.388e-17, and every
selected ID and sign is retained. This verifies a local additive representation.
Added k2 directions, the induced ridge penalty and jointly refitted k5 weights
all contribute to the changed procedure; regularization alone is not isolated.
Whole-insert composition includes the edit neighborhood and discards position.
Neither the design nor this result identifies a causal biological interaction.

Ten sequence fits finish in 5.39 seconds on one CPU thread, raising the ledger
to 445 starts (130 Rummel), with no new variance calibrations. Independent source
and feature checks, saved-solution replay and the root rerun pass. All exact
point fields and decisions agree; the largest difference across 13,274 floating
score fields is 1.333e-15. No raw archive, GPU work or heldout outcome is used.
The fixed plan, full eight-model table, attempts, algebraic witness and independent
receipts are in
`build/neural-predictor-reviewed-20261005/rummel-context-interaction/`.

A parallel primary-source review identifies [BIRD](https://www.geneprediction.org/bird/)
as a separate future measurement option matching one DNA pair and three RNA
pairs. Its model propagates finite-read DNA uncertainty and RNA sample variation
with published priors. Its latent shrunken effect differs from the current mean
log-ratio label; runtime, posterior behavior and assay assumptions need separate
qualification. Full BIRD inference has not run, and no labels change. The cited
review is retained in that directory's `count-method-research/`.

### Count-uncertainty sampler qualification

The [historical Swift source](https://github.com/bmajoros/BIRD/blob/746ef6c49d3c915069beb82c2b392a757aadea83/swift.C)
provides a small two-stage alternative that runs with the installed NumPy.
Source history resolves a material discrepancy: a later change replaced the
documented constant concentration with an allele-asymmetric adaptive rule.
The new research kernel follows the older fixed-concentration distribution,
with explicit count orientation and the symmetric zero-count correction.
It returns stable signed log-effect draws and fails explicitly on invalid
inputs or numerical endpoints. It does not copy the author's defective median
routine or assume an independently seeded simulation reverses exactly.

Five root tests, two independent controlled checks and four independent distribution
and boundary tests pass. Root reruns reproduce both independent suites. The latter includes
joint DNA/RNA moments, so matching each marginal distribution alone is insufficient.
A concentration-input regression fails before repair and passes afterward;
deliberate formula faults are rejected. All fixtures are synthetic. This qualifies
the specified sampling kernel, not biological calibration or a new predictor.

Swift pools RNA replicates and uses a modular approximation in which RNA does
not update the DNA fraction. It therefore models finite-read uncertainty without
estimating RNA replicate heterogeneity or unobserved DNA preparation variation.
The existing labels and acceptance decisions remain fixed. A real measurement
experiment still requires a defined summary, calibration checks and its own
prospective comparison. No dependency installation, biological count read,
predictor fit or heldout outcome is added; the ledger remains 445 predictor starts.
Code, source history and qualification evidence are in
`build/neural-predictor-reviewed-20261005/swift-fixed-qualification/`.

### Measurement support for the constant-context selections

A fixed zero-fit comparison now traces the control's full-target gain over
precision k5 to one substituted edit: chr5_60696427 replaces chr11_130745047
in fold 1. The other nine picks remain unchanged, and both models predict the
same directions on all eleven edits in their union. Membership accounts for
the entire utility increase of 0.3394163430; direction changes contribute zero.

| Target | Precision k5 utility | Constant-context utility | Difference |
|---|---:|---:|---:|
| Full | 0.8475601833 | 1.1869765263 | +0.3394163430 |
| A | 1.6284615583 | unavailable | unavailable |
| B | 0.5153615077 | unavailable | unavailable |

The new edit has no usable A/B views under the existing joint measurement QC.
Its ledger records nonpositive or missing half-0 reference counts for DNA,
baseline RNA R3 and stimulated RNA R2. Both complete contrasts are therefore
unavailable. No missing edit is replaced or averaged away. These QC facts leave
support unresolved; they do not establish a cross-view loss or a false full effect.
The 40.05% full-target selection gain remains an observed control result.

Independent replay and the root rerun reproduce every value across 107 frozen
pins and three published outputs. Synthetic checks reject the old-model binding
and original-family remapping faults. No new fit, target, QC threshold or
acceptance opportunity is introduced. The control remains ineligible, the
composition recipe stays closed, and the ledger remains 445 predictor starts.
The fixed comparison, complete outputs and receipts are in
`build/neural-predictor-reviewed-20261005/rummel-constant-views/`.

### RNA replicate homogeneity

A fixed count diagnostic now tests whether the three baseline RNA samples share
one alternate fraction under independent binomial reads. Conditioning on their
totals and total alternate count removes the unknown fraction. Each of the same
948 eligible TRAIN edits receives 99,999 multivariate-hypergeometric null tables;
exact integer comparisons measure the predeclared alternate-fraction range.
The p-value is (b+1)/100000 with inclusive ties. Bonferroni correction uses all
1,209 original source pairs, assigning the 261 excluded pairs p=1 without reading
their counts. At alpha=0.05, rejection requires b<=3.

The working model rejects for 438 edits, 46.20% of the tested population.
The changed selection chr5_60696427 has RNA Ref/Alt counts 7/159, 27/90 and
24/97. Its alternate-fraction range is 407/2158, about 18.86 percentage points.
Two null draws meet that discrepancy, giving p=0.00003 and corrected p=0.03627.
This is a finite Monte Carlo test; no exactly enumerated tail is claimed and
no near-threshold row is rerun. The fraction of rejected edits is descriptive.

The result requires the next uncertainty model to account for replicate
variation. It cannot separate biological variation from technical effects,
barcode/PCR dependence or other violated read assumptions. It does not declare
the selected effect false or resolve missing A/B views. Original-universe
correction gives unconditional control under the specified count null, not a
biological discovery guarantee or a correction for earlier adaptive research.
Nonrejection does not establish calibrated pooled uncertainty.

Three root kernel tests, four independent kernel checks and two independent
driver tests pass, including root reruns and deliberate-fault witnesses. A
separately prepared checker, qualified by two tests and two rejected faults,
replays all 94,799,052 null tables once. Every scientific field in all three
outputs matches exactly, with all 64 frozen pins intact. Root verifies all
5,688 input count cells, fixed decisions and complete independent result equality.
No sequence model is fitted, relabeled or advanced; the ledger remains 445.
The protocol, complete results and primary-source references are in
`build/neural-predictor-reviewed-20261005/rummel-rna-homogeneity/`.

### Replicate-aware posterior inference

The BIRD author model now has a qualified marginal density and a bounded
rejection sampler for independent posterior draws. The density integrates each
RNA replicate's Beta fraction and retains the full shared-DNA posterior. An
independent rational/80-digit reference checks 430 comparisons; four root tests
and three further mathematical checks also pass. Eight deliberately introduced
mathematical faults fail. Review found and repaired loss of the RNA fraction at
extreme concentration values.

Four fixed synthetic cases compare 60,000 draws per method against a separately
implemented generative count-equality sampler. All 20 marginal CDF comparisons
pass: the largest distance is 0.00836667 against the fixed simultaneous limit
0.01939913. The total is 480,000 retained states. All eight sampler calls finish
within their proposal and CPU caps; combined sampling time is 1.341 seconds in
this run on one CPU thread. Independent deterministic replay verifies all 20
CDF distances and 120 quantile intervals from the eight saved arrays.

These results qualify inference on the fixed synthetic cases. The next
scientific test must compare prediction of an omitted RNA repeat under this
model and a pooled baseline, with a fixed TRAIN-only design. No target changes
or predictor advancement follow from the numerical pass. The ledger remains
445 predictor starts and heldout outcomes remain unopened. Code, protocols,
complete posterior draws and review evidence are in
`build/neural-predictor-reviewed-20261005/bird-marginal-qualification/` and
`build/neural-predictor-reviewed-20261005/bird-iid-qualification/`.

The subsequent real-count pilot selects four Dev TRAIN sources from four
distinct merged genomic components before reading their count values, with
three RNA omission roles each. It stops at the first source, 259, omission 0:
the prior rejection sampler retains only 610 of 20,000 required states within
2,000,000 proposals. No omission finishes and no aggregate score is produced.
The failed attempt is preserved without a retry or budget extension. This
identified assay-scale inference efficiency as the next implementation problem.
Four root tests, four independent driver tests and five rational score tests
pass, including omitted-count isolation and a repaired missing-input-pin guard.
The frozen plan, failure and evidence are in
`build/neural-predictor-reviewed-20261005/bird-predictive-pilot/`.

### Omitted RNA repeat prediction

A new weighted integration method preserves the BIRD author posterior while
retaining each integrated likelihood contribution. It uses a fast lifted-Gamma
ratio, a uniform weight bound from Beta log concavity, and fixed independent
prior particles. Exact rational and finite-product references qualify the
likelihood, bound, pooled Beta-binomial forecast and score calculations. Eleven
root tests and seventeen independent kernel/reference tests pass; eleven
deliberate numerical faults fail. Five further driver tests check original data
hashes, omitted-count isolation and failed-attempt retention. The qualification
is in `build/neural-predictor-reviewed-20261005/bird-importance-qualification/`.

The first weighted run completes the same four TRAIN sources and all twelve
omissions with 500,000 raw particles per job. Its gain estimate is 0.00730225,
but its conditional 99% simulation interval crosses zero. It stays closed as
unresolved. A separately frozen precision follow-up changes only the particle
count to 2,000,000 per job and uses fresh random streams. Its result alone is the
final precision inference; no particles or intervals are pooled between runs.

| Fixed four-component comparison | Mean CRPS loss |
| --- | ---: |
| BIRD repeat model | 0.05490723648280356 |
| Pooled Beta-binomial | 0.06219211366022381 |
| Pooled minus BIRD | 0.007284877177420252 |

The full conditional 99% Monte Carlo gain interval is
[0.004364886143923947, 0.009756654633382796]. The fixed aggregate comparison
passes, with 11.71% lower estimated loss. All four component mean gains are
positive. One omission, source 325 role 2, favors pooling and stays included:
gain -0.01653387, interval [-0.01970678, -0.01277616]. This is a comparison of
predicting RNA allele fractions from the other repeats, with CRPS on the
fraction scale. It does not measure sequence-predictor accuracy.

All twelve jobs finish within the 45 CPU-second cap, on one CPU thread with no
GPU. Their summed wall time is 24.581 seconds in this single execution, excluding
entry/exit checks. Independent replay verifies all 172 pins and twelve moment
archives, recovering the aggregate gain within 4.4e-16 and every score within
5.7e-15. A separate 720-state check on the actual training count pairs agrees
with the finite-product likelihood within 8.74e-11 log units.

The interval covers simulation error conditional on these four selected TRAIN
components; it is not a biological population confidence interval. Coverage of
the exact author posterior additionally assumes the qualified 1e-8 log-error
allowance holds over the prior. Admission already used the three RNA repeats.
The two methods differ in priors, DNA use and heterogeneity, so the comparison
does not isolate one mechanism. The next measurement test must address wider
coverage and DNA-normalized edit effects before changing predictor labels.
The predictor ledger stays at 445 and held-out outcomes remain unopened.
Complete plans, results and independent replay are in
`build/neural-predictor-reviewed-20261005/bird-weighted-predictive-pilot/` and
`build/neural-predictor-reviewed-20261005/bird-predictive-precision/`.

### DNA-normalized effect prediction

The next fixed comparison uses the same four TRAIN components and separate
barcode halves. A contains half0 DNA and baseline RNA R3; B contains half1 DNA
and baseline R4/R5. Each model forecasts the other view's DNA and RNA counts
jointly, conditional on its totals. A prospectively fixed symmetric half-count
defines a finite log2 RNA-odds minus DNA-odds effect. This is a new observable,
distinct from the original unsmoothed effect and latent BIRD theta. Inherited
eligibility still includes KCl; no unavailable row is replaced.

The first run stops after four of eight jobs because a rounded point CDF exceeds
one. A twenty-particle synthetic counterexample reproduces 1.0000000000000002;
the failed real histogram was not saved. The incomplete attempt stays sealed.
The repair caps the point CDF at its mathematical upper bound and saves raw
moments before scoring. Fourteen affected numerical tests and seven controlled
driver/repair tests pass. A separately frozen follow-up retains the model,
cohort, budget and acceptance rule, with fresh samples and no pooling of attempts.

All eight repaired comparisons complete, with 32 million proposals in sixteen
model-jobs. Mean CRPS is 0.5105277260763628 for BIRD and 0.6072364029592754 for
pooling, a 15.926% estimated reduction. The score is in log2-effect units; this
percentage is a relative loss reduction, not an error percentage.

| Gain: pooled minus BIRD | Estimate | Conditional 99% simulation interval |
| --- | ---: | --- |
| A to B | 0.10229415308172436 | [0.07301889919057131, 0.12283729312501623] |
| B to A | 0.09112320068410087 | [-0.016084221784680497, 0.13674844805547964] |
| Equal mean of eight comparisons | 0.09670867688291263 | [0.0284673387029454, 0.12979287059024794] |

All eight point gains are positive. The fixed primary rule requires positive
lower bounds in both directions, so its decision is **unresolved**. Source 259
B-to-A has the widest gain interval, [-0.3057030992967357, 0.11850868640566226],
and remains included at equal weight. The positive overall interval does not
override that directional requirement.

Independent replay verifies all 104 frozen inputs and all sixteen raw/scored
archive pairs. Exact rational point-CDF normalization reproduces the scored
results within 1.12e-15. The run uses one CPU thread, no GPU, and completes within
the 120 CPU-second cap; summed model-job wall time is 17.9154 seconds in this
single execution. No automatic precision extension follows this result.

This is technical prediction on four previously exposed TRAIN components, with
numerical assumptions and conditional simulation bounds. Shared cultures,
preparation and cross-half DNA exchangeability limit its biological scope.
CRPS improvement can reflect better uncertainty estimates without improving point
labels. Labels, predictor gates and the 445-fit ledger remain unchanged; validation
and confirmation stay closed. Artifacts and independent replay are in
`build/neural-predictor-reviewed-20261005/bird-effect-predictive-pilot/` and
`build/neural-predictor-reviewed-20261005/bird-effect-repair/`.

## Broader measurement coverage

A deterministic diagnosis attributes most of the earlier reverse-direction
simulation interval to source 259. A separately qualified CDF calculation now
uses both prefix and suffix importance moments, with the full simultaneous
probability budget. Independent exact-rational checks cover 202 interval corners;
fault tests reject the old normalization and unsafe squared-weight underflow.
The earlier four-source result remains closed and unresolved.

A new fixed sixteen-component comparison then stops before sampling: three
selected sources lack views under the inherited barcode rule. Sources 195 and
675 have supported integer-zero cells; source 850 also has an empty reference
half because its two reference barcodes share the same hash parity. Independent
closeout and root replay preserve all sixteen sources, the 32-direction ledger
and zero completed jobs. No unavailable source is replaced.

A separate experiment qualifies balanced barcode assignment: sort barcode
IDs by hash within each source and physical allele, then alternate halves.
Assignments are sealed before reading selected counts. The original full-count
baseline/KCl admission remains; new baseline views allow an observed zero allele
count at positive combined totals, while preserving missing values. The same
sixteen sources, finite half-count effect, models and both-direction CRPS rule
remain fixed. This changes the measurement views, not the sequence predictor.

The balanced construction recovers usable views for **15 of the 16** fixed
components. Source 850 now has one reference and one alternate barcode in each
half, but half1 baseline RNA R4 has zero combined reads. Its required view remains
unavailable under the declared positive-total rule. All sixteen components and
32 direction entries are retained; the comparison closes **incomplete with zero
simulation draws**. The source is neither replaced nor assigned an imputed
measurement. Independent closeout and root replay agree on all 257 frozen inputs,
assignments and ledger entries. This supplies a concrete support diagnosis, not
a score for either model.

Nine independent synthetic reader groups and ten deliberate faults verify the
assignment, zero/missing distinction, physical allele ownership and count-sum
reconciliation. Four preparation tests and four controlled driver tests cover
the frozen file identities, selected-only parsing, separate DNA denominators,
fixed cohort and unchanged directional gate. An initial preparation manifest
omitted the imported count parser; its source check stopped before archive
access. That failed attempt is preserved. The corrected manifest uses identical
code and the same unused random-stream recipe.

Evidence is retained under `bird-complement-qualification/`,
`bird-effect-coverage/`, `bird-coverage-admission-diagnosis/` and
`bird-balanced-effect/` and `bird-balanced-effect-run/` in the reviewed build
directory. Predictor fits remain 445, with validation and confirmation unopened.

## Full-count treatment-response forecasts

The measurement task now returns to the predictor's actual response target:
the KCl-minus-baseline difference in allele effects. The shared DNA log term
cancels exactly from that difference, as the original measurement contract
already specifies. Baseline DNA-effect forecasts are indirect evidence for this
target. The implemented two-condition extension shares one DNA posterior and
fits separate baseline/KCl repeat-variation distributions. It forecasts the
difference between mean half-count-smoothed RNA log odds, using the original
full counts. This finite observed target remains distinct from the original
unsmoothed predictor labels. Cultures in the two conditions are unpaired.

All sixteen fixed components have the required full-count views. The first
attempt saved one model's two million proposals, then stopped on squared-weight
underflow: some tiny positive weights had squares rounded to zero. A separately
qualified arithmetic repair retains that mass and accounts explicitly for the
rounding loss. The continuation reused the saved archive exactly and preserved
the models, cohort, random streams, confidence family and both-direction rule.

The continuation completed **31 of 32 paired comparisons** before the remaining
117 CPU seconds expired. All **62 saved model results** independently replay,
with a maximum numerical difference of 5.69e-13. The final source-direction has
an input record only; its interrupted proposals were not saved. The full
32-comparison ledger remains intact, and no overall average or directional
acceptance decision is reported for this incomplete run. The cap closes this
plan; the final comparison is not silently retried or dropped.

An explicit, independently reviewed operational amendment then permits one
15-CPU-second process to finish only the two missing original random streams.
It preserves all 62 completed model records and archives byte for byte, including
their scores. The final comparison completes in 3.92 wall seconds. All sixteen
components, 32 paired comparisons and 64 model jobs are now present. The earlier
incomplete runs remain immutable; the original CPU limit was not met.

The complete experiment has the following descriptive CRPS estimates; lower
loss is better. Gain is pooled loss minus BIRD loss. The intervals bound
simulation uncertainty conditional on the fixed data, model and declared
numerical assumptions.

| Training view | BIRD CRPS | Pooled CRPS | Gain interval |
|---|---:|---:|---:|
| One repeat per condition → other two | 0.74815 | 0.81868 | [-3.21985, 0.13245] |
| Two repeats per condition → other one | 0.64273 | 0.80490 | [-14.48194, 0.57873] |
| All 32 comparisons | 0.69544 | 0.81179 | [-8.85089, 0.35559] |

The estimated overall loss is **14.3% lower** than pooling, with positive point
gains in both directions. Both directional intervals cross zero, so the fixed
decision is **unresolved**. Precision is insufficient to establish the gain;
this result does not advance the predictor. The original simultaneous confidence
family applies to the complete experiment; it is not a new 99% guarantee
conditional on the 31 scores already seen when the completion was planned.

Independent replay verifies the 186 copied job files, the last input, the two
new model scores and all 32 aggregate contributions. Root replay agrees. The
comparison contains 128M unique scientific proposals; physical work also
includes an unknown duplicated prefix from the interrupted job. No proposal
count is inferred for those unsaved draws.

This test asks whether accounting for assay variation improves forecasts of
repeated treatment-response measurements. The next inference question is how
to obtain sufficiently precise estimates at low CPU cost without changing the
target or weakening the acceptance rule. Improved sampler efficiency must be
measured by valid interval width per CPU time, not effective-sample-size estimates
alone.

Predictor fits remain 445; labels and predictor acceptance requirements are
unchanged, and validation and confirmation remain unopened.

Evidence is retained under `bird-interaction-design/`, `bird-interaction-pilot/`,
`bird-moment-underflow/`, `bird-interaction-resume/` and
`bird-interaction-completion/` in the reviewed build directory. The original
numerical failure and the CPU-limited continuation are both preserved.

### Locating the inference bottleneck

A zero-draw diagnostic now replays the normalization bounds of all 64 models
exactly. Independent scalar arithmetic and the root replay agree within
1.39e-17. Five BIRD jobs have zero lower bounds: source 2227 in the first
direction, and sources 1528, 2227, 195 and 445 in the second. Their bounded-range
terms account for 87.6% to 99.0% of the normalization radius. Even setting the
empirical variance to zero leaves the diagnostic radius above each saved mean;
the worst ratio is 245.45. This diagnoses the current calculation at its saved
mean, rather than predicting how a new sample will behave.

Consequently, a proposal that merely lowers variance while retaining the same
loose weight bound is not yet a precision remedy. The simple half-prior defensive
mixture would halve the expected scaled denominator while keeping that bound,
doubling the relative range penalty. A useful proposal needs a demonstrably
better bound as well as useful samples. Changing the confidence allocation after
inspecting these draws would be a separate method analysis, not fresh
confirmation.

The diagnostic and independent replay are retained under
`bird-interaction-precision-replay/`. Its first attempt passed all numerical
comparisons but failed to serialize a NumPy Boolean; that attempt and partial
output remain under `bird-interaction-precision-diagnosis/`. The output repair
passes a regression that fails on the original code. No labels, predictor fits,
validation outcomes or acceptance decisions changed.

A deterministic follow-up intersects every saved CDF band with known full-count
support and monotonicity. This preserves coverage on every realization covered
by the original confidence event, with no new probability allocation. All 64
models replay independently and in the root check with exact agreement. Only
the outer endpoint bounds change; all model, directional and overall CRPS
intervals are unchanged. The result remains unresolved. This closes support
and monotonicity alone as a precision repair for these saved forecasts.
Four root checks and four deliberately broken variants qualify the change;
evidence is retained under `bird-cdf-refinement/`. No additional simulation was
needed.

A second follow-up tightens the global likelihood bound itself. For two
nonnegative RNA kernels under a common latent law, their integrated product is
at most one quarter of the squared maximum kernel sum. Exact rational cells
bound that sum over the full unit interval. This covers all latent Beta
parameters, including their tails, and leaves the saved weights unchanged.
The independent coverage review establishes a common minimal-range confidence
event while retaining the existing numerical assumptions and failure budget.

The bound decreases for 12 BIRD jobs and tightens 10 CRPS intervals. The
two-training-repeat gain interval narrows by 19.69%, from
`[-14.48194, 0.57873]` to `[-11.57646, 0.51875]`; the overall interval becomes
`[-7.39816, 0.32560]`. Point forecasts retain their 14.33% mean CRPS reduction
against pooling. Four normalization lower bounds remain zero and both
directional intervals cross zero, so the decision remains unresolved.

All 64 scores pass independent replay; the largest arithmetic difference is
`3.98e-13` in a squared-error interval. The retained checker reproduces the
independent reports byte for byte. Four tests detect three deliberate faults.
The scoring replay takes 9.68 seconds and uses no new draws or predictor fits.
Evidence is retained under `bird-amgm-envelope/`.

The next two source-only prototypes test stronger bounds. A compact
two-parameter Beta search retains the trivial bound one at its fixed 128-split
budget. A unimodal uniform-mixture construction with nine weights also fails
to improve the hardest example. Both negative results remain recorded under
`bird-compact-envelope/` and `bird-unimodal-envelope/`.

Refining that construction to 129 weights reduces the source-2227 two-repeat
joint bound by 27.57%, from `0.07004757` to `0.05073612`. Independent integer
arithmetic verifies all 4,096 kernel-cell bounds and both condition certificates.
The full 64-score replay changes only that job's normalization interval; its
lower endpoint remains zero. All point scores, CDF bands and CRPS intervals
remain exactly as in the accepted AM-GM result, and the decision remains
unresolved. This isolates the remaining calculation problem: a qualified
inference method or separately fixed sampling budget is needed to resolve the
weak normalization. Evidence is retained under `bird-unimodal-weights/` and
`bird-unimodal-replay/`. No new draws or predictor fits were used.

A finite prior-region proposal now addresses the range penalty itself. It
partitions the shared DNA state and each condition's scale, concentration and
RNA logit, retaining every tail. An exact 288-state oracle verifies the
importance-weight identity and detects lost shared-DNA and tail contributions.
A cheap likelihood cell bound passes 924 independent exact Beta comparisons.

On the already inspected source-2227 two-repeat training example, the fixed
table estimates a weight bound of `0.0010985474`, compared with `0.0507361212`,
a 46.18-fold reduction. The 60- and 80-digit tables agree below `4.83e-42`
relative error across 2,612 numeric values and finish in 0.70 seconds together.
These are proposal-design measurements. CDF arithmetic and truncated sampling
still require qualification before the estimate can become a certified weight
bound or improve forecast intervals. No new samples or fits were run, and the
acceptance result is unchanged. An independent 100-digit reconstruction checks
2,610 computational entries per table and all 1,018 file identities. The retained
root replay matches its scientific output exactly. Evidence is under
`bird-stratified-proposal/`.

The proposal's prior-probability block now passes on the fixed representable
grid. All 24 Beta/Gamma bin masses have positive directed-Decimal enclosures;
all tails and the original zero-target concentration mass are retained. The
fixed calculation takes 0.744 CPU seconds. Independent reconstruction verifies
every interval, its sums and widths, source ownership and all 1,092 file pins;
root replay matches. Six root test groups and 216 independent exact Beta cases
cover arithmetic, tail loss and source guards. Evidence is under
`bird-proposal-priors/`.

A separately checked sampler design uses explicit Normal density envelopes and
zero weights for missed cells, avoiding a conditional Normal CDF per particle.
Its 32-state exact oracle preserves the target and weight-range identity. The
new density-envelope table, truncated Beta/Gamma draw laws and full numerical
weight calculation remain to implement and qualify before a prospective run.
The earlier 46-fold table estimate does not certify this different proposal.
Forecast acceptance and predictor fit counts remain unchanged.

The explicit density-envelope table now has a certified exact-law range of
`0.0033101627`, 15.327 times smaller than the previous `0.0507361212` bound.
The fixed calculation takes 0.310 CPU seconds. Independent review checks all
2,048 proposal bounds, 512 likelihood bounds, 1,024 scheme decisions and 1,172
file identities. Exact rational composition verifies the range; an independent
directed-Decimal calculation reproduces every interval. Root replay matches.
Seven source test groups pass and five deliberate faults fail. Evidence is
under `bird-proposal-envelope/`.

The complete sampler now passes sixteen root test groups and independent
selection, prior-draw, mixture-weight and forecast-integration checks. The
130-digit full-weight oracle covers fifteen fixed cases across all eight DNA
cells and every proposal branch. Its maximum log-weight error is
`3.3998e-11`, within the unchanged `2.01e-8` allowance. Root replay reproduces
the independent numerical results. Evidence and regression failures are under
`bird-proposal-sampler/`. This qualification generated no scientific random
draws and fitted no predictor. A fixed two-million-attempt replacement for the
weak source-2227/direction-1 estimate is prepared under `bird-proposal-forecast/`;
it preserves the model, target and both-direction acceptance rule. Its confidence
family covers all 64 original streams and the one replacement.

That fixed run completed in 26.308 CPU seconds. At the same two-million-attempt
budget, effective particles rose from 22.960 to 912.779, and the normalization
lower bound rose from zero to `2.63989e-8`. The CRPS interval narrowed 44.633%.
These are improvements in sampling precision for the same model and target;
the replacement's point CRPS moved from 0.19603 to 0.22477 and remains in the
comparison. Independent and root replay reproduce all 65 scores, CDF hashes
and aggregates exactly. Directional gain intervals remain
`[-3.21990, 0.13247]` and `[-10.12452, 0.51882]`, so forecast acceptance is
unresolved. Three original normalization lower bounds remain zero
(2227/direction 0, 195/direction 1 and 445/direction 1); the new interval also
needs further precision. No predictor fits were run; validation and confirmation
remain closed.

## Evaluation contract

Freeze the following details in a versioned protocol before opening new held-out
outcomes. Record any already seen examples or published selected variant results.

1. Keep reference/edited pairs, overlapping loci and identical or
   reverse-complement sequences in the same connected family. Declare assembly,
   coordinates, orientation, assay flanks, replicates and exclusions. Fit target
   transformations and normalization on training data where the source permits;
   disclose preprocessing already pooled by the original study.
2. Define the effect and its units from the assay. Use family-weighted MSE for
   effect prediction and paired family-level uncertainty for baseline comparisons.
   Preserve zero-effect and independently tuned dinucleotide baselines; include
   the prior k-mer approach as an additional comparator. Fix candidate selection,
   bootstrap procedure and acceptance rule before final evaluation. Report
   preparation/batch uncertainty separately when independent replication permits.
3. Require a positive lower confidence bound for MSE improvement against every
   required baseline. For experiment selection, declare whether the objective is
   effect magnitude or a signed change, then freeze its utility metric and edit
   budget. Require a positive lower bound for ranking utility improvement against
   random selection and each eligible sequence baseline using paired,
   family-aware uncertainty. Better MSE alone does not establish useful ranking.
   Report prediction errors by edit type and sequence similarity; do not remove
   difficult rows using their measured effect.
4. Verify that identical input pairs predict zero effect and swapping the pair
   reverses the signed prediction when the declared effect definition requires
   it. Include shuffled-target and sequence-destroying controls with their exact
   interpretation. Do not assume random or scrambled DNA is biologically neutral.
5. Treat time dependence and transfer between cell types as separate hypotheses.
   A temporal claim must beat a time-averaged ablation. Held-out sequence families
   in one experiment establish that scope only; new donors, preparations or
   laboratories require corresponding independent evidence.

## Implementation boundaries and review

Use `examples/neural-regulation/regulation.py` and
`tests/test_neural_regulation.py` for the existing 171-base, seven-timepoint
contract. Preserve its negative result. If phase 1 selects a different assay,
give it a separate example and assay-specific tests; reuse only compatible
sequence, provenance and evaluation code. Keep experimental models outside the
compiler core. Save source receipts, protocols, fitted models and reports under
`build/` while developing; include the reproducible deliverable and source
attribution in the final example.

An independent reviewer checks each phase's evidence before advancement. The
review must identify concrete leakage, scope, baseline, numerical or claim
problems. Fix confirmed implementation defects and record scientific failures.
The reviewer does not fit the model or select it from final test performance.
Update this roadmap with each completed phase and its evidence path.

## Connection to the full goal

After a useful reporter predictor is validated, the next research milestone is
predicting an endogenous regulatory response to a DNA intervention in a declared
neural preparation. Subsequent milestones connect that response to development,
cell identity and structure, then activity and learning. Each connection needs
its own intervention evidence. The full
[genome-to-brain objective](../README.md#epigenesis) remains open.

## Sources under evaluation

- [Kreimer, Ashuach, Inoue et al., 2022](https://doi.org/10.1038/s41467-022-28659-0): the existing seven-timepoint H1 reporter experiment, whose test is already examined.
- [GSE152404](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152404): independent modern/archaic sequence-pair experiments; preparation and stage compatibility require an explicit contract.
- [Rummel et al., 2023](https://doi.org/10.1016/j.cell.2023.09.015), [GSE225815](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225815): DevLib baseline/KCl source partition and TRAIN count analysis complete; weak interaction repeatability directs further measurement diagnosis. This MPRA subseries sits within GSE225817.
- [Deng et al., 2024](https://doi.org/10.1126/science.adh0559): 270-base developing-cortex and organoid reporter experiments; usable pretrained checkpoints have not yet been verified.
- [Kosicki et al., 2025](https://doi.org/10.1038/s41467-025-60064-1): selected experiment; 270-base reporter sequences in WTC11-derived excitatory neurons, infected on day 7 and measured on day 14. Sequence and measurement intake are verified.
- [Ghandi et al., 2014](https://doi.org/10.1371/journal.pcbi.1003711): regulatory classification precedent for gapped sequence patterns; the seventh round adapts these features to paired neuronal reporter effects.
- [Gosai et al., 2024](https://doi.org/10.1038/s41586-024-08070-z): measured regulatory activity and original Malinois model/data used by the MPAC source under audit.
- [MPAC author export](https://github.com/Reilly-Lab-Yale/coda_mpac/tree/965fea279de8dace4580c81150b02136521fd890/hf_export): pinned model/input/strand contract; supervised-exposure admission remains under review.
- [Salomon et al., 2026 preprint](https://doi.org/10.64898/2026.07.16.738760): 270-base neuronal regulatory-variant experiment; the pinned primary-cohort design, overlap checks and processed-coefficient admission pass independent replay. [Pinned author code](https://github.com/kircherlab/80K-Analysis/tree/a686f1d552f97ebcc0224057b5a835002f7a3dc8).
- [Matoba et al., 2026](https://doi.org/10.1038/s41593-026-02454-2): 150-base cortical-variant reporter experiment in primary neural progenitors under vehicle and Wnt stimulation; source admission, count import, fixed TRAIN measurement export and first bounded model comparison complete; shared ridge fails the TRAIN progression screen. [GSE289858](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289858).

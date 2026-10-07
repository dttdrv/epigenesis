# Neural DNA-edit predictor roadmap

Updated 2026-10-07. Predictor development continues in parallel with the measured
neural-evidence workbench. A successful new predictor is not a prerequisite for
that practical deliverable. This document preserves the separate predictor
research plan and its independent acceptance requirements.

The latest completed ledger is **280 biological fits**. The bundled
[two-cycle predictor](../examples/neural-predictor/README.md) is runnable on CPU
and Metal. Its validation utility is 0.07158510 and MSE is 0.02862935; utility
passes its comparator and MSE remains 0.01480% above its strict limit. The
TRAIN uncertainty audit is complete. Confirmation for the newer Salomon cohort
remains unopened. The sections below retain earlier rounds and their outcomes.

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

## Roadmap

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

Each successor uses the smallest model justified by training/validation evidence.
An available pretrained model is eligible only after its exact input, output,
training data, checkpoint, license and runtime are verified. Do not invent assay
flanks, convert scores into biochemical rates, or treat unspecified pretraining
overlap as independent validation. Exclude evaluation families whose outcomes
informed supervised pretraining, checkpoint selection, tuning or calibration.
If that exposure cannot be resolved, the model's result is exploratory only.
If the chosen dataset cannot support the declared experiment, return to phase 1
with the reason recorded.

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
- [GSE225817](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225817): neural reporter variant experiments under metadata review.
- [Deng et al., 2024](https://doi.org/10.1126/science.adh0559): 270-base developing-cortex and organoid reporter experiments; usable pretrained checkpoints have not yet been verified.
- [Kosicki et al., 2025](https://doi.org/10.1038/s41467-025-60064-1): selected experiment; 270-base reporter sequences in WTC11-derived excitatory neurons, infected on day 7 and measured on day 14. Sequence and measurement intake are verified.
- [Ghandi et al., 2014](https://doi.org/10.1371/journal.pcbi.1003711): regulatory classification precedent for gapped sequence patterns; the seventh round adapts these features to paired neuronal reporter effects.
- [Gosai et al., 2024](https://doi.org/10.1038/s41586-024-08070-z): measured regulatory activity and original Malinois model/data used by the MPAC source under audit.
- [MPAC author export](https://github.com/Reilly-Lab-Yale/coda_mpac/tree/965fea279de8dace4580c81150b02136521fd890/hf_export): pinned model/input/strand contract; supervised-exposure admission remains under review.
- [Salomon et al., 2026 preprint](https://doi.org/10.64898/2026.07.16.738760): 270-base neuronal regulatory-variant experiment; the pinned primary-cohort design, overlap checks and processed-coefficient admission pass independent replay. [Pinned author code](https://github.com/kircherlab/80K-Analysis/tree/a686f1d552f97ebcc0224057b5a835002f7a3dc8).

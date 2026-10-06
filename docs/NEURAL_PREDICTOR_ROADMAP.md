# Neural DNA-edit predictor roadmap

Updated 2026-10-06. Predictor development continues in parallel with the measured
neural-evidence workbench. A successful new predictor is not a prerequisite for
that practical deliverable. This document preserves the separate predictor
research plan and its independent acceptance requirements.

Phase 1 is complete. Eight bounded model-selection rounds have completed. The
latest comparison fitted two shallow tree ensembles, with and without sequence
context. The context model improved this comparison's validation error and
signed ranking, but its MSE of 0.0155401 remains above the required k-mer
baseline's 0.0153782. The historical eligible winner therefore remains selected
and still fails ranking. Reserved confirmation is unevaluated.

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

Independent replay verifies both new models. The ledger contains 196 fits:
160 eligible candidates and 36 controls/baselines. The earlier
exact-word control reduces validation MSE by 0.467% against the required k-mer
baseline and has ranking utility 0.0251110 versus 0.0011730 for random selection.
These point estimates identify a follow-up lead. The user preferred keeping the
original eligibility rule after reviewing a proposed control nomination. That
proposal remains unadopted; the original control remains ineligible. Both required
sequence baselines now have independently checked training normal equations and
reproduced validation scores and rankings. The completed experiment and its
accounting are preserved. No reserved outcomes have been accessed.

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
| 2. Build one successor | Use training and validation data to select a sequence model. Compare it with zero effect, the dinucleotide model and the existing k-mer approach under the same protocol. Save every attempted configuration. | Runnable predictor with fixed preprocessing, parameters and calibration; selection uses no held-out outcomes. | Eight bounded rounds complete; original eligibility retained; paired trees do not qualify |
| 3. Verify the implementation | Exercise data isolation, numerical calculations, input rejection, batch prediction and exact-source replay. Save the final candidate and baseline identities before scoring. | Independent reference calculations agree; deliberately broken controls fail; the independent review has no unresolved correctness or leakage finding. | Per-experiment checks ongoing; final candidate review awaits qualification |
| 4. Evaluate once | Score the frozen candidate and baselines on the reserved data. Save all eligible predictions, exclusions, errors and uncertainty. | Both effect-prediction and ranking acceptance pass. A failure remains a failure; any subsequent model change requires a new confirmation set. | Not started |
| 5. Deliver the experiment-selection tool | Expose batch prediction and ranking through a documented command, with a worked real-data example and reproducible report. | A clean run reproduces the accepted result and ranking; validation supports the stated assay and use. | Not started |

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
[genome-to-brain objective](../README.md#research-objective) remains open.

## Sources under evaluation

- [Kreimer, Ashuach, Inoue et al., 2022](https://doi.org/10.1038/s41467-022-28659-0): the existing seven-timepoint H1 reporter experiment, whose test is already examined.
- [GSE152404](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152404): independent modern/archaic sequence-pair experiments; preparation and stage compatibility require an explicit contract.
- [GSE225817](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225817): neural reporter variant experiments under metadata review.
- [Deng et al., 2024](https://doi.org/10.1126/science.adh0559): 270-base developing-cortex and organoid reporter experiments; usable pretrained checkpoints have not yet been verified.
- [Kosicki et al., 2025](https://doi.org/10.1038/s41467-025-60064-1): selected experiment; 270-base reporter sequences in WTC11-derived excitatory neurons, infected on day 7 and measured on day 14. Sequence and measurement intake are verified.
- [Ghandi et al., 2014](https://doi.org/10.1371/journal.pcbi.1003711): regulatory classification precedent for gapped sequence patterns; the seventh round adapts these features to paired neuronal reporter effects.
- [Gosai et al., 2024](https://doi.org/10.1038/s41586-024-08070-z): measured regulatory activity and original Malinois model/data used by the MPAC source under audit.
- [MPAC author export](https://github.com/Reilly-Lab-Yale/coda_mpac/tree/965fea279de8dace4580c81150b02136521fd890/hf_export): pinned model/input/strand contract; supervised-exposure admission remains under review.

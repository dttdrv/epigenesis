# Cross-culture forecast of later enhancer-response changes

Preanalysis specification, 18 October 2026. Freeze this document, the executable, independent qualification, source identities and runtime before converting numerical count cells or fitting.

## Question and scope

Do early enhancer measurements from the other infected cultures predict later within-culture changes in an omitted culture, for excluded sequence families?

This is a retrospective assessment on the already exposed Kreimer et al. H1 neural-induction assay, DOI 10.1038/s41467-022-28659-0. It uses the deposited DNA/RNA count matrices instead of pooled MPRAnalyze alpha estimates. It does not alter the original DNA-only predictor, its acceptance rule, or its confirmation state. No claim of a new regulatory mechanism, independent donor validation, or population-level preparation generalization follows.

The three culture labels are rep1, rep2 and rep3. The paper describes three infections with two viral batches, but their batch assignment is absent from the retained annotations. Do not describe the three omissions as three independently held-out viral batches.

## Source admission before numerical conversion

Pin the downloaded DNA/RNA matrices, both annotation tables, original design FASTA/provenance, existing TRAIN identity/family list, recovered author reducer and all implementation/runtime dependencies. Retain the original 4,143 TRAIN targeted-row identities and existing source-qualified parent/locus/family assignments; raw measurement availability cannot redefine these families or introduce validation/test rows.

Each annotation table has 13,629 labels: seven hours [0,3,6,12,24,48,72], three culture labels and 649 slots. Its first field is the observation name; the remaining fields describe time, DNA/RNA type, replicate and slot. Join each count-matrix header explicitly to its own annotations. Align DNA and RNA by the complete canonical (culture,time,slot) key within the same construct row. Verify source row identities and dimensions. Do not rely on column order, or treat slot numbers as physical barcodes shared between different inserts.

Stream numerical conversion only for the fixed TRAIN targeted IDs and their required source-qualified WT IDs. Other source rows may be traversed to find identities, but their count fields are not converted or summarized. Require decimal integer TSV fields and reject negative, nonfinite, malformed or otherwise unsupported count representations. Do not silently coerce unknown sentinels, missing fields or parse failures into zero. A source/schema/integrity failure stops admission; a documented absent measurement remains explicitly missing.

The recovered author file `ratios_analysis2_DNA_RNA.m` has SHA-256 `ae5356a3fef76368b58cfd2c7604cc612466f837847260b38097c9cb8b4a6a66`. Lines 21-25 align slots locally; lines 39 and 45 select DNA>0 and RNA>0 within a single culture/time sample. The reducer does not require support across all times. Upstream matrix construction was not recovered, so earlier global selection or censoring remains unresolved. The assessment is conditional on those exported matrices and the sample-local positive-support rule.

## Fixed raw-count observable

For construct i, culture r and hour t, define S as its insert-local slots with both D>0 and R>0 in that same sample. Let D_i=sum(D over S) and R_i=sum(R over S), and retain |S|. If S is empty, this construct/sample is unmeasured. Do not use the author's NaN-to-zero replacement at line 53.

For a variant's parent locus, sum these construct-level R and D totals over all unambiguous, source-qualified WT IDs assigned to that locus. Count each matrix row once. Never duplicate a WT row because its FASTA record repeats. Preserve the source's oriented-parent identity checks; conflicting parent sequences are an admission failure.

If both variant and pooled WT have supported positive totals, define

```text
E(r,t) = [log(R_variant)-log(D_variant)-log(R_WT)+log(D_WT)] / log(2)
Y(r)   = [E(r,24)-E(r,12), E(r,48)-E(r,12), E(r,72)-E(r,12)]
```

Compute in log space; do not multiply counts into potentially overflowing products. There is no pseudocount, fitted asinh scale, library-size factor or all-time fit. Sample-depth factors cancel in the variant/WT contrast. Keep all positive low-count support; introduce no outcome-derived count cutoff or inverse-variance weighting.

This is a new pooled-count, DNA-normalized log2 contrast. It is not MPRAnalyze alpha, and it is not the author's mean of individual (1+RNA_CPM)/(1+DNA_CPM) ratios. The new ratio of summed counts weights slots differently. Its support rule conditions on positive matched RNA/DNA slots, so an all-zero-RNA or unsupported construct is missing, not established to have zero activity. Actual zero contrasts or changes calculated from supported measurements remain valid zero observations.

## Outer exclusion and prediction population

Use all three omitted-culture choices crossed with the existing five whole-family fold labels, SHA256('epigenesis-temporal-v1\n'+family) modulo five.

For omitted culture C and held family fold F, no C measurement enters fitting, tuning, feature preprocessing or training-label construction. No held-F family enters fitting or tuning in any culture. The other two cultures are A and B.

Prediction eligibility for a variant in F requires complete supported early E values at 0,3,6,12 hours in both A and B. It does not depend on any observation, support flag, count total or activity in C. Retain every fixed TRAIN row in the population ledger, including rows lacking one or both source prefixes and the exact missingness reason. Eligibility is specific to the culture omission and is identical for every method.

Train two ordered-direction observations for each eligible training variant when their required measurements exist: early A predicts changes Y(B), and early B predicts changes Y(A). Each observation requires its source prefix and its target's supported 12/24/48/72-hour measurements. Record missing training observations; all methods use the same training observations. Missingness in allowed training labels may determine whether a training observation is usable, but held-culture measurements cannot.

At evaluation, obtain one forecast from A's early context and one from B's, then average their three-dimensional predicted changes. The omitted culture's 12-hour measurement is never an input or anchor added to predictions; it is used only later to construct Y(C). This removes identical-measurement anchor noise between predictors and the change target.

## One candidate and four controls

The only candidate is ridge with features [mean(E0,E3,E6), E12], predicting the three changes directly, with an unpenalized intercept.

The four required controls are:

1. Ridge from mean(E0,E3,E6,E12) only.
2. Ridge from E12 only.
3. Zero change at every late horizon.
4. Training-population mean change at each horizon.

For learned models, use penalties [0.001,0.01,0.1,1,10]. Tune each model separately on the four remaining family-fold labels inside the outer training complement. Inner held families are excluded in both cultures and both directions. Center and standardize features from the current training complement only. Zero-variance coordinates contribute zero; do not estimate feature statistics from an evaluation family. Exact tuning ties choose the smallest penalty.

Give each training family equal total weight; within a family give each available target-culture direction equal weight, and within a direction give its usable variant rows equal weight. Use the same weights for centering, scaling, regression, population-change estimates and inner validation loss. Inner validation aggregates family means, rather than giving larger families or cultures with more supported rows extra weight.

All evaluated directions and models are fixed before results. Do not fit a full-data deployment model, add Timestral, search additional architectures, or change the observable after this result. There is no ranking endpoint or ranking gate in this first bounded raw-count assessment.

## Prediction seal and target availability

Use separate prediction and scoring stages for each outer (culture,fold) evaluation. Save and hash its complete prediction population, input-only eligibility/rejection ledger, fitted-model receipts and every method's predictions before its scorer computes held-target availability or scores. Once those predictions are sealed, that scorer may access C at 12/24/48/72 hours to determine target support and construct labels.

Cross-fitting means a culture's reduced counts can already exist in the cache because it is a training culture for another outer omission. This is not an untouched-data claim or a requirement never to read those bytes. Access is permitted only according to the predeclared ownership graph; the omitted culture's measurements must never enter the fit, eligibility or input preprocessing for the omission that holds it out. No scoring result or held-target support diagnostic from another omission can change eligibility, models or predictions.

For the primary score, a sealed, input-eligible curve is observed-complete only when its omitted-culture E12,E24,E48,E72 are all supported. Use this identical observed-complete set for all methods. Retain unsupported rows unscored in the ledger, with per-time, variant/WT and missing-row reasons. Do not impute zero, substitute a control value, use partial curves in the primary score, or quietly remove them from the stated prediction population.

## Loss, uncertainty and conditional gate

For each complete curve, average squared error over the three Y coordinates. Average over its observed-complete variants within each (family,omitted-culture) cell. For the overall family score, average its available omitted-culture cells equally, then average observed families equally. This explicitly defines a family-balanced risk conditional on supported complete measurements, not the risk of every originally eligible curve.

For each control, the paired gain is control loss minus candidate loss. Keep each family's culture cells together. Resample whole observed families 20,000 times with PCG64 seed 20261018. Report simultaneous nominal 99% two-sided percentile intervals for the four predeclared contrasts, tails 0.01/(2*4), using NumPy linear quantiles. Bootstrap intervals are conditional on the fitted cross-validation models, the observed-family set and these three cultures; they do not include fitting instability, missing-outcome uncertainty or population preparation variance.

Report each culture omission and outer family fold separately. The conditional forecast gate requires all four overall gain lower bounds strictly above zero and all four point gains strictly above zero within each of the three culture omissions. Each of the 15 (culture,outer-fold) strata must contain observable families; otherwise the planned assessment is incomplete. Nonfinite predictions, fit failure, missing required folds or integrity failures are incomplete/invalid results and cannot be hidden by dropping rows or methods.

There is no arbitrary coverage-percentage threshold. Missing endpoints do not alone invalidate an observed-case comparison. Report, for every culture/fold, the original fixed population, source-prefix eligibility, sealed predictions, complete targets, partial/unsupported targets, observed families and matched-slot/read support. The success label must be `conditional_on_supported_held_culture_measurements`. Never turn that conditional pass into a claim about unsupported rows or the general preparation population. If missingness prevents a stated stratum or contrast from being observed, report that limit rather than extrapolating its performance.

## Qualification and retained records

Before fitting, synthetic and source-identity checks must establish explicit header/annotation joins; no cross-insert slot mixing; exact WT summation without duplication; positive-support/missing-zero distinctions; log-ratio depth-factor cancellation; culture and family exclusions; inner-complement preprocessing and tuning; hierarchical weights; ordered-direction training; averaging two evaluation predictions; target-anchor isolation; and whole-family bootstrap ownership.

Mutation of any omitted-culture value or support flag must leave that outer fit, tuning, input preprocessing, prediction eligibility and prediction unchanged. Mutation of a held-family value after 12 hours in either source culture must likewise leave them unchanged. Source-culture held-family early values are legitimate prediction inputs and are allowed to change predictions. Synthetic missing targets must remain visible and must not be assigned zero loss or a zero biological response. Demonstrate relevant fault controls in isolated copies.

Retain exact file/runtime hashes, the matrix/source semantic audit, all count-reduction totals/support flags for admitted rows, source/family/culture/fold ownership, training-observation rejections, chosen penalties, predictions sealed before scoring, target availability, per-family/culture losses and paired contrasts. Run on one CPU/BLAS thread, without GPU, new data downloads, original validation/test count conversion, or changes to existing published predictor results.

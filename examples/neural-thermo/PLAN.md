# Thermodynamic transfer test — frozen plan

Written 2026-11-16, before any thermodynamic model was fitted and before any
test-family label was scored. This file is hashed at freeze time (see
`FREEZE.json`). Anything not stated here is exploratory.

## Question

Epigenesis' circuit scanner (`examples/neural-circuit`) rests on one explicit,
never-tested assumption: that a binding ratio computed from an *in-vitro* PWM
partition sum transfers to regulatory activity *in living cells*. The earlier
sequence-to-effect tests (research/) found that statistical learners do not
transfer. This test asks whether a physics-constrained enhancer model does.

> Does a thermodynamic occupancy model whose only sequence input is a PWM
> partition sum, with non-negative activator and repressor strengths per motif
> and timepoint, predict held-out-family motif-perturbation effects during
> human neural induction (GSE188264, 0–72 h) — and does that depend on the
> actual binding specificity of the PWMs, or only on having many features?

## Data and splits (unchanged)

`examples/neural-regulation/regulation.py::load_data`: 9,241 variants in 585
loci, 582 connected families (352 train / 112 validation / 118 test, assigned
by the frozen hash rule). Target: `regulation.target(row, scale)` (asinh/log2
units, variant minus parent mean, seven timepoints). Primary metric:
family-weighted MSE over `kind == "targeted"` rows (`regulation.family_losses`).
The 118 test families were examined by earlier experiments (README, Release
1.7), so this result is retrospective, not fresh confirmation; it is scored
once under the gates below.

## Model (frozen form)

For sequence `s`, motif `m` (532 bundled JASPAR 2024 log-odds PWMs):

```text
o_i      = (exp(f_i) + exp(r_i)) / 2          strand-averaged site odds
theta_i  = kappa o_i / (1 + kappa o_i)         saturating site occupancy
Theta_m  = sum_i theta_i                        partition input
F_t(s)   = log2(1 + sum_m A[m,t] Theta_m(s)) - log2(1 + sum_m R[m,t] Theta_m(s))
effect_t = F_t(variant) - F_t(parent)           A, R >= 0
```

`A` is activating and `R` quenching strength; both are free per motif and
timepoint (the time-0 basal term cancels in the difference). The loss is the
family-weighted mean squared error over the seven timepoints plus
`lambda (|A|^2 + |R|^2) + lambda rho * sum_t [(A_t - A_{t-1})^2 + (R_t - R_{t-1})^2]`.
Optimiser: `scipy.optimize.minimize` L-BFGS-B, bounds `[0, inf)`, start
`A = R = 0.01`, `maxiter = 500`, `ftol = 1e-10`; deterministic.

## Selection (train → validation only)

Grid: `kappa in {1e-4, 1e-3, 1e-2, 1e-1}`, `lambda in {1e-6, 1e-5, 1e-4, 1e-3,
1e-2}`, `rho in {0, 10, 100}` (60 fits). Fit on training rows (all kinds,
family-weighted, `scale = training_scale(train)`); choose the minimum
validation MSE on targeted validation rows (ties: larger `lambda`, then smaller
`kappa`). Development on train/validation (feature checks, grid scale) is
allowed and logged; the grid above is the search space and may not be widened
after the first test scoring.

## Final fit and the one test scoring

Refit the chosen configuration on train + validation (all kinds,
`scale = training_scale(fitting)`). Score the 118 test families once. Baselines
are refitted by the repository's own code and must reproduce its published
numbers (gate G0): zero effect (0.360870), dinucleotide (0.319358), `k5`
(0.285185), motif-temporal (0.252079).

## Controls

* **PWM-shuffle null** (decisive): 20 nulls, seeds 1–20. Every PWM's position
  columns are permuted (composition and information content preserved,
  specificity destroyed). Features are recomputed; the *selected* `kappa,
  lambda, rho` are reused without re-selection; same optimiser.
* **Linear same-feature control** (report-only): ridge on
  `Theta(variant) - Theta(parent)` at the selected `kappa`, penalty from
  `{0.01, 0.1, 1, 10, 100}` chosen on validation. Tests whether the saturating
  non-negative form adds anything over the same inputs used linearly.

## Statistics

Paired family bootstrap (resample the 118 families), 10,000 draws, seed
`20261116`, 95% percentile interval on the MSE reduction (baseline − model).

## Gates

| Gate | Condition |
|---|---|
| G0 integrity | Baselines reproduce the published values to `1e-5` |
| G1 | Gain vs zero effect: lower bound > 0 |
| G2 | Gain vs `k5`: lower bound > 0 |
| G3a | Real-PWM model beats the average-of-null per-family loss: lower bound > 0 |
| G3b | Real-PWM test MSE is lower than all 20 PWM-shuffle nulls |
| G4 | Non-inferior to motif-temporal: point gain ≥ −0.005 |

**Decision rule.** "Thermodynamic transfer supported on this dataset" iff
G0–G4 all pass. G3 is the claim that binding *specificity* carries the signal.
Passing G1–G2 without G3 means the physics form is a good regulariser but not
shown to be physics. Superiority to motif-temporal (lower bound > 0) is a
stretch result (G5), not required.

## Exploratory (not in the decision rule)

* **G6 dynamics sanity.** For the final model, the activating-strength centre
  of mass over time, `sum_t t A[m,t] / sum_t A[m,t]` (hours), for two fixed
  motif sets chosen before fitting from the bundled panel (it lacks core
  POU5F1/SOX2/NANOG/KLF4 motifs):
  early `{POU2F1::SOX2, KLF2, KLF17, TFAP2C, ESRRB, TBX3, FOXD3, SOX15}`,
  neural `{PAX6, ZIC1, ZIC4, SOX1, SOX21, SOX4, NEUROG1, NEUROG2, ASCL1,
  POU3F2, POU3F1, EMX2, FOXG1, LHX5, SIX3, OLIG2, HES5, RFX4, NR2F1, NFIA,
  NFIX}`. One-sided Mann–Whitney, early earlier than neural, alpha 0.05.
* Top motifs per timepoint, stratum breakdown (random, scrambled), and the
  `kappa/lambda/rho` surface.

## Conditional design stage

Only if G0–G4 pass: derive 171-bp enhancer designs with prescribed temporal
profiles from the model and score them with the independent frozen
motif-temporal and `k5` predictors. These are in-silico candidates for
experiments, not results.

## What this cannot show

One assay, one cell system, reporter constructs rather than endogenous loci,
published preprocessing, previously examined test families. A pass supports the
transfer assumption for GLI-independent TF binding in this system; it does not
validate the neural-tube circuit or any fate prediction.

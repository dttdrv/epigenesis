# Thermodynamic transfer test — results

Run 2026-11-16 under the frozen [plan](PLAN.md) ([hashes](FREEZE.json)). One
scoring of the 118 test families; every number below is in
[`result/report.json`](result/report.json). Runtime 436 s on four CPU cores.

## Outcome

**The pre-registered decision rule is not met: thermodynamic transfer is not
supported on this dataset.** Three of the six gates failed. The conditional
design stage was therefore not run.

| Gate | Result | Value |
|---|---|---|
| G0 baselines reproduce published numbers | pass | zero 0.360870, dinucleotide 0.319358, `k5` 0.285185, motif-temporal 0.252079 (exact to 6 digits) |
| G1 beats zero effect, lower bound > 0 | pass | gain 0.0923, 95% interval [0.0376, 0.1580] |
| G2 beats `k5`, lower bound > 0 | **fail** | gain 0.0166, [−0.0075, 0.0387] |
| G3a beats mean PWM-shuffle null, lower bound > 0 | **fail** | gain 0.0229, [−0.0079, 0.0527] |
| G3b beats all 20 PWM-shuffle nulls | pass | 20 of 20 |
| G4 non-inferior to motif-temporal (point gain ≥ −0.005) | **fail** | gain −0.0165, [−0.0411, +0.0053] |
| G5 (stretch) superior to motif-temporal | fail | same interval |

Test MSE on 1,415 targeted edits in 118 held-out families:

| Model | Family-weighted MSE |
|---|---:|
| Zero effect | 0.360870 |
| Dinucleotide | 0.319358 |
| `k5` ridge | 0.285185 |
| **Thermodynamic (this test)** | **0.268602** |
| Linear on the same site sums (control) | 0.261051 |
| Motif-temporal (frozen, 1.7.0) | 0.252079 |
| PWM-shuffle nulls (20) | mean 0.291487, range 0.273922–0.302086 |

Selected configuration (train → validation only): `kappa = 0.01`,
`lambda = 1e-3`, `rho = 10`; validation MSE 0.343482 against 0.340252 for
motif-temporal and 0.376377 for `k5` fitted the same way.

## What the result does and does not say

1. **The physics inputs carry signal; the physics form does not obviously add
   to it.** The real PWMs beat every one of 20 column-shuffled copies (the
   closest null, seed 4, is 0.0053 worse). Under exchangeability that is a
   one-sided rank p of 1/21 ≈ 0.048. But the pre-specified paired family test
   against the *average* null has an interval that includes zero, so G3a fails;
   the evidence that binding specificity matters is suggestive, not
   established under the plan. The same site sums used *linearly*
   (control, penalty 0.01, the edge of its grid) do better than the
   saturating non-negative model (0.2611 against 0.2686; gain interval
   [−0.033, +0.017]), so the thermodynamic nonlinearity and sign constraint
   are not shown to help.
2. **Not a good enough model to claim parity with the free-weight linear
   predictor.** It trails motif-temporal by 0.0165, a gap whose interval
   reaches +0.005 but is mostly on the wrong side of zero. Motif-temporal also
   uses `k5` words; this model does not.
3. **The recovered "TF dynamics" are essentially static.** The fitted
   activating strengths barely change across 0–72 h (for example
   `POU2F1::SOX2` 0.62 at 0 h, 0.69 at 48 h; the same four motifs lead at
   every timepoint, see `top_motifs_by_hour` in the report). The temporal smoothing and
   shrinkage selected on validation produced a nearly time-invariant program,
   so this fit does not recover regulator dynamics and the
   "inverse problem" ambition is not met.
4. **G6 could not be evaluated.** Only 2 of the 8 pre-specified early motifs
   kept non-zero activating strength (the rest were shrunk to exactly 0 by the
   non-negativity bound), below the 3-motif minimum that the code applies. That
   minimum is in `confirm.py`, not in `PLAN.md`. The medians in the report
   (early 15.3 h, neural 22.8 h) are descriptive only.
5. **Behaviour by stratum** (test, `report.json["categories"]`): on random
   mutation controls, where the true effect is near zero, the thermodynamic
   model is close to the zero-effect floor (0.1267 against 0.1223) and better
   than motif-temporal (0.1420) and `k5` (0.1494); on scrambled controls it is
   0.4722 against 0.4464 for motif-temporal.

## Caveats

* The test families were examined by earlier experiments, so this is
  retrospective; the gates were fixed before the model was fitted, and the test
  set was scored once.
* Selection used warm starts along the `lambda` path; a cold start of the chosen
  configuration reproduces its validation MSE (0.343471 against 0.343482).
  The final fit and nulls used cold starts, as the plan states.
* One assay, one cell system, reporter constructs, 118 families. The
  intervals are wide; a larger test could resolve G2 and G3a in either
  direction.
* The 532 bundled motifs lack core POU5F1, SOX2, NANOG and KLF4 profiles, which
  limits any dynamics check in this system.

## What this suggests next

A confirmatory test would need data not already used here, which this sandbox
cannot download. The most promising specification is the one the control points
to: PWM partition-sum features used linearly together with `k5` words, compared
with motif-temporal on a fresh held-out assay, with the same shuffled-PWM null.
It should be frozen and hashed before any new outcome is read.

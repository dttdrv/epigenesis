# Current delivery — Epigenesis 1.7.0

The interactive entry point is [Epigenesis Studio](../../studio/README.md) (`0.2.0`),
with English and Bulgarian interfaces. The compiler package is version **1.7.0**.

## DNA-site simulation and Enhancer Circuit

Simulate accepts a supported single substitution in the nine-base Nkx2.2 GLI
binding site, compiles its nucleotides into a validated Development Module,
and drives the Herrera-Delgado (2020) four-gene neural-tube ODE under a user
signal history. **Enhancer Circuit** (1.7) generalizes this to arbitrary
enhancers (`9` to `500` bp in Studio) via strand-symmetric Boltzmann GLI2/GLI3
partition sums, uncovering 7 cryptic cell-fate bifurcations (`|Δreporter| < 0.02`
with endogenous `Nkx2.2 -> Olig2` switch) across `chr15:90936956-90937127`.

These are conditional mechanistic experiments. The broader
genome → development → neural structure → activity → learning path remains
the research objective; its separate examples are not one validated biological model.

## Delivered predictors

Studio and `epigenesis.py predict` provide two NumPy CPU checkpoints (no
training occurs during use):

- **Siraj Huber-linear (default, `epigenesis-siraj-huber-linear-v1`)**: SK-N-SH
  24h reporter effect. On reused tuning data (`16,986` edits), equal-family MSE
  is `0.28766713` vs k5 `0.29661090` (3.02% reduction); on five external neural
  reporter studies, mean correlation is `0.0490` and MSE is 9.52% above zero.
- **Kreimer 0–72h (`--kreimer`, `epigenesis-kreimer-motif-temporal-v1`, 1.7)**:
  motif-temporal ridge model (5-mers + 532 JASPAR motifs) predicting H1 neural
  induction across `[0, 3, 6, 12, 24, 48, 72]` hours. Across 118 held-out
  GSE188264 families (`1,415` targeted edits), retrospective family-weighted MSE
  is `0.252079` vs k5 `0.285185` (95% bootstrap gain `[0.007810, 0.061363]`), and
  temporal-contrast MSE beats time-averaging (`0.031285` vs `0.032108`).
  Conditional on published MPRAnalyze preprocessing (`biological_development_accepted: false`).

The earlier activity-transfer and attention models remain available explicitly
for reproducing the separate Salomon assay results. Their confirmation is spent;
Siraj RESERVE remains sealed.

## Evidence and reproduction

- [1.7 release contract](../../docs/releases/RELEASE-1.7.md) and [Enhancer Circuit](../../examples/neural-circuit/README.md).
- [Predictor guide](../../examples/neural-predictor/README.md): Siraj inputs, archived models, and input contract.
- [Technical paper](../../docs/paper/epigenesis-technical-paper.pdf) ([Bulgarian](../../docs/paper/epigenesis-technical-paper.bg.pdf)): mechanisms, results and scientific qualifications.
- [Frozen research](../../research/README.md): thirteen cross-study and variant-prioritization tests, with plans and hash verification.
- [Roadmap](../../docs/research/NEURAL_PREDICTOR_ROADMAP.md): dated experiment history.

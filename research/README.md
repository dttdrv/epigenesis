# Research record

This folder holds the evidence behind the claims in the [technical paper](../docs/paper/epigenesis-technical-paper.pdf) ([Bulgarian](../docs/paper/epigenesis-technical-paper.bg.pdf)). It is a record, not a pipeline: the scripts reproduce what was scored, and the plans show what was fixed before outcomes were read.

## Frozen tests

Thirteen tests asked whether a sequence model predicts the effect of an unseen DNA alteration better than a simpler alternative. Each was planned, hashed and scored once. None earned its label.

| # | Test | Directory | Headline |
|---|---|---|---|
| 1 | Leave-one-study-out, three studies | `hnp-loso-20261019` | H2 r 0.015 vs S1 0.060 |
| 2 | Combined predictor, six studies | `hnp-loso6-20261019` | Cfit r 0.046 vs S2 0.050 |
| 3 | Sei ridge | `sei-features-20261019` | r −0.003; blend 0.041 vs S2 0.047 |
| 4 | Quadratic context head | `context-loso6-20261019` | S3 − S2 = −0.001 |
| 5 | Activity-transfer model | `activity-transfer-loso6-20261019` | r −0.010 vs S2 0.047 |
| 6 | Pooled convolutional network | `cnn-loso6-20261019` | r 0.014 vs S2 0.047 |
| 7 | Zero-shot expressed-TF binding | `mechanism-zeroshot-20261019` | r 0.0015 vs S2 0.0465 |
| 8 | Huber-linear transfer | `huber-transfer-20261020` | r 0.0490 vs S2 0.0465; MSE 9.5% above zero effect |
| 9 | eQTL direction, fetal neocortex | `eqtl-brain-20261019` | M5 0.566 vs B3 0.538 |
| 10 | eQTL direction, Jerber | `eqtl-brain-20261019` (`PLAN_CONFIRM.md`) | M5 0.603 vs B2 0.600 |
| 11 | eQTL direction, sensory neurons and microglia | `eqtl-brain-20261019` (`PLAN_REP.md`) | M5 0.518 vs B2 0.499; microglia 0.570 |
| 12 | Credible-set ranking | `eqtl-causal-20261019` (`PLAN.md`) | nearest TSS 0.460 beats all models |
| 13 | Credible-set ranking with accessibility | `eqtl-causal-20261019` (`PLAN2.md`) | 0.417 and 0.420, below nearest TSS |

`constraint-factor-20261019` is an exploratory, unfrozen analysis of evolutionary conservation (phyloP100way). It is kept for completeness and is not one of the thirteen.

Model labels (k5, S1, S2, Huber-linear, M4, M5 and the baselines) are defined in Section 5.1 of the paper.

## Contents of a test directory

`huber-transfer-20261020` (test 8) has a flatter layout (`run.py`, `metrics.py`, `prediction/`, `scoring/`) because it was exported separately; its freeze records pinned inputs rather than file digests, so the checker skips it.

| File | Role |
|---|---|
| `PLAN*.md` | The question, populations, models, baselines, metrics and decision gates, written before outcomes were read |
| `FREEZE*.json` | SHA-256 digests of the plan, scripts and inputs |
| `scripts/` | Feature, training and scoring code |
| `derived/` | Small result files: the scored result and the synthetic-positive and synthetic-null qualification runs. Large arrays stay in the untracked `build/` tree |
| `RESULTS.md` | The scored outcome and the label decision |

Sealed prediction files (`PREDICTIONS*.npz`), downloaded raw data and virtual environments are not in this repository. The Siraj RESERVE partition is still sealed.

## Checking the hashes

```bash
python3 research/check_freeze.py
```

The command recomputes each digest for the files held here and exits nonzero on a mismatch. Files that are not in this folder (for instance, large inputs) are reported as not present. To include them, point `--build` at the original working tree:

```bash
python3 research/check_freeze.py --build build/neural-predictor-reviewed-20261005
```

Some public copies remove editorial attribution or normalize file paths. [`PUBLICATION.json`](PUBLICATION.json) records both their original and published SHA-256 identities. The original freeze digests remain unchanged. The checker verifies these copies separately; it does not call an edited copy a byte-identical pre-evaluation artifact. Numerical results, model parameters and data are unchanged. The published context experiment directory is named `context-loso6-20261019`; external inputs must be supplied under the published directory names.

The hashes are local. A match to an original freeze records byte identity, not the time of freezing; an author who rewrites a plan and its hash together would not be caught. Publication-copy identity is a separate, later record. External pre-registration would close the timing gap. See Section 5 of the paper.

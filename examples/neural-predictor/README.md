# Experimental neuronal DNA-edit predictor

Score a single-base regulatory edit in a 270-base DNA sequence. The model
predicts the ALT minus REF log2 reporter-activity coefficient for the
WTC11-NGN2 neuronal assay studied by [Salomon et al.](https://doi.org/10.64898/2026.07.16.738760).
Its intended use is prioritizing edits for measurement and estimating their
direction and effect size.

The fitted predictor runs now. Its joint biological acceptance requirement
remains open: validation ranking utility improves, while prediction error is
still above the fixed limit. Final confirmation for this cohort is unopened.

## Run

Use Python 3.12 in a separate environment; the compiler itself needs neither
PyTorch nor NumPy.

```sh
python3.12 -m venv /tmp/epigenesis-predictor-env
. /tmp/epigenesis-predictor-env/bin/activate
python -m pip install -r examples/neural-predictor/requirements.txt
python examples/neural-predictor/predictor.py \
  --input examples/neural-predictor/edits.json \
  --output /tmp/edit-predictions.json
```

CPU is the default. On a Mac with supported Metal acceleration, add
`--device mps`. Choose an output path that does not already exist. Inference
uses bundled weights; it needs no network access after dependency installation.

`edits.json` contains **synthetic input sequences**, with an edit, an identical
pair, its swapped pair and its reverse complement. They demonstrate the input
contract and are not measured biological examples. Replace them with your own
JSON array; each record must contain exactly `id`, `reference` and `sequence`.
Identifiers must be unique. Both sequences must contain exactly 270 uppercase
ACGT bases and differ at no more than one position. Use the actual assay insert;
the command does not invent flanks or pad a shorter input.

The output gives the predicted `effect`, the local-word `baseline`, and the
unweighted `long_correction` and `short_correction`. A positive effect predicts
higher ALT reporter activity. Component values are statistical contributions,
not biochemical rates. The output also records the model-manifest hash, device
and dependency versions.

## Model

The predictor combines a reverse-complement canonical 5-mer ridge baseline
with two compact sequence scorers:

\[
\widehat{\Delta y}
= b + 0.020655743871331515\,g_{4096}
    + 0.2223167743737464\,g_{1024}.
\]

Each scorer has 11,368 parameters: a width-32 convolution, four-head attention
with signed and absolute distance biases, and a feed-forward block using
[Interlace's rational gate](https://misul.org/interlace),
\(a(1+b)/(1+a^2+b^2)\). The models differ in training schedule length: 4,096
and 1,024 updates. The baseline has 512 coefficients. Including the two mixture
weights, the combined predictor contains 23,250 fitted scalar parameters.
This is a two-model ensemble; it is not the full Interlace language architecture.

The neural alphabet is AGCT. Inference canonicalizes full sequences against
their reverse complements, sorts unique sequences, scores both strands in
batches of 64, and repeats the last sequence to fill a partial batch. Each
float32 strand score is converted to float64 before averaging. Subtracting the
reference score from the edited score gives the correction. Baseline counts,
matrix multiplication and the final mixture use float64. The training loss
scale is metadata only and is not multiplied into inference predictions.

The two full-TRAIN checkpoints and the baseline were exported without optimizer
state or cached measurements. Every tensor is unchanged. NPZ loading disables
pickles, checks file hashes, and validates names, shapes, dtypes and finiteness.
[The manifest](data/manifest.json) records source checkpoint identities, model
configuration, training recipes, original loss scales and the fitted mixture.

## Validation results

Training uses 21,781 edits in 8,774 sequence families. Validation uses 7,390 edits
in 3,029 families. Family assignment groups overlapping loci and identical or
reverse-complement sequences to reduce leakage. Each family has equal total
weight when calculating MSE. Mixture weights are learned from five-fold
excluded-family predictions within TRAIN, then applied unchanged to validation.

| Predictor | Family-weighted MSE ↓ | Signed top-ten utility ↑ |
|---|---:|---:|
| Always predict zero | 0.028751711602 | 0 (abstains) |
| Random selection, TRAIN-selected direction | n/a | 0.016748752672 |
| Dinucleotide ridge, k2 λ=1 | 0.028747144221 | 0.016577960175 |
| k5 ridge, λ=0.1; protocol comparator | 0.028625109823 | 0.035108902127 |
| k3 ridge, λ=0.01; highest baseline utility | 0.028697606252 | 0.215156319143 |
| Strict requirement | < 0.028625109823 | > 0.035108902127 |
| Long correction, calibrated | 0.028622834172 | 0.020415119456 |
| Short correction, calibrated | 0.028631012565 | 0.076092485687 |
| **Bundled two-cycle mixture** | **0.028629346486** | **0.071585099786** |
| Subsequent precision-sampling trial | 0.028620934594 | 0.003420568634 |

Ranking selects one edit per family, orders candidates by absolute predicted
effect, resolves ties by identifier, and takes ten families. Utility is the mean
of `sign(predicted effect) × measured effect` for those edits. It rewards correct
direction and measured magnitude; it is not a discovery rate or sign accuracy.
The protocol selects one dinucleotide model and one k-mer model by validation
MSE, then takes the strictest error and utility requirements among those selected
models, zero effect and random selection. Here k5 supplies both limits. It is
not the highest-utility model among the 16 baseline fits; k3 λ=0.01 ranks better
than the mixture. The [aggregate results](results/validation.json) preserve every
baseline, including this comparison. Random utility is the expected result of
uniform family selection, a uniform edit within each family, and the direction
chosen from the TRAIN mean. Predicting zero abstains and has zero signed utility.

The mixture improves utility by 103.89% over the protocol's k5 comparator and
327.41% over random selection. It misses the k5 error limit by 0.01480%, defined
as `100 × (MSE_model / MSE_k5 − 1)`. This is a relative baseline gap, not an error
percentage or closeness to perfect prediction. The short calibration improves
utility by 116.73% over k5 but has higher MSE. Both criteria must pass together.
These are point estimates from adaptive model development, not independent
confirmation or confidence bounds.
The research ledger contains 287 biological fits, including failed experiments;
the subsequent precision trial adds seven, while packaging adds none. The
original failed confirmation on an earlier cohort remains recorded in the
[roadmap](../../docs/NEURAL_PREDICTOR_ROADMAP.md).

MSE has units of squared log2 reporter effect. The mixture's RMSE is 0.16920
log2-effect units. Its reduction in MSE relative to always predicting zero,
`1 − MSE_model / MSE_zero`, is only 0.00425592, or **0.4256%**. This is zero-effect
skill, not centered R². Independently checked validation aggregates give a target
mean of 0.01674875 and centered variance of 0.02847119; conventional weighted
`R² = 1 − MSE_model / variance` is **−0.00555**. The validation-mean predictor in
that definition is a descriptive reference, not a deployable fitted baseline.
These summaries are derived from existing aggregate receipts without refitting.

No empirical between-culture noise ceiling has been established for this
cohort. The TRAIN standard-error audit below does not supply one. Culture
correlations and MSE values near 0.01554 in the historical roadmap concern the
earlier Kosicki dataset; they cannot be used to normalize these Salomon results.

## What the uncertainty audit found

The subsequent TRAIN-only audit recovered moderated coefficient standard errors
for all 21,781 training edits. Weighted squared standard error is 61.78% of the
observed effect second moment. This describes model-based reported uncertainty;
it does not measure an irreducible-noise floor. Source normalization and
empirical-Bayes moderation predate the benchmark's train/validation split.

| Reported uncertainty quartile | Baseline OOF MSE | Baseline + long correction | Baseline + short correction |
|---|---:|---:|---:|
| Lowest | 0.011644723 | 0.012121252 | 0.011664991 |
| Second | 0.018041001 | 0.018777220 | 0.018057714 |
| Third | 0.029375779 | 0.030212505 | 0.029466035 |
| Highest | 0.053640273 | 0.054993014 | 0.053665003 |

Both **uncalibrated unit corrections** worsen excluded-family MSE in every
quartile, including the most precisely measured group. The audit does not
evaluate the calibrated mixtures. It gives further model work a concrete
failure to explain and adds no fits or gate changes.

The subsequent fixed experiment tests precision-based sampling, with each
fitting population's probabilities proportional to
`family_probability / (SE² + E_family[SE²])`. Original targets, k5 priors,
architecture and 1,024-update recipe are unchanged. Six fresh neural fits
(full TRAIN and five family-fold complements) and one TRAIN-only scalar
calibration give alpha 0.42450904. Its MSE is 0.02862093, 0.01459% below k5,
but utility falls to 0.00342057, below the random comparator. Both original
criteria remain required, so this trial closes without nomination or a
replacement of the bundled weights. Independent replay matches all 109,300
saved values and verifies the seven-fit accounting. Confirmation stays unopened.

## Verification and scope

The portable inference path reproduces all 29,560 archived Metal values exactly:
baseline, both corrections and combined prediction for all 7,390 validation
inputs. These are matches to saved model outputs, not 29,560 exact matches to
biological measurements. The CPU path's largest combined-prediction difference
from those Metal values is 2.651e-8 on the checked Mac. This portability check reads input
sequences and saved predictions, without converting validation outcomes or
opening confirmation. Floating-point agreement on other hosts remains to be
measured. The archived biological scores above use the original Metal outputs.

Run the independent contract tests with the predictor dependencies installed:

```sh
python -m unittest discover -s examples/neural-predictor -p 'test_*.py' -v
```

The [results directory](results) includes aggregate validation values, independent
mixture and uncertainty review receipts, and the inference-export check.
The receipts identify the original local experiments; they are evidence records,
not a substitute for reproducing their training and assay intake. This package
supports inference reproduction and synthetic contract tests. The full local
training workflow, raw source tables, per-row measurement exports and reserved
confirmation are not included. No endogenous-development or cross-cell-type
prediction result follows from this reporter experiment.

## Sources

- [Salomon et al., 2026](https://doi.org/10.64898/2026.07.16.738760), neuronal regulatory-variant assay; [pinned analysis code](https://github.com/kircherlab/80K-Analysis/tree/a686f1d552f97ebcc0224057b5a835002f7a3dc8).
- [Interlace](https://github.com/MisulOrg/Interlace/tree/b7aeae7a3364be0f834ac2777acb887d7d2859cc), rational feed-forward gate adapted to the sequence scorer.
- [Breiman, 1996](https://statistics.berkeley.edu/sites/default/files/tech-reports/367.pdf), out-of-fold nonnegative regression for combining predictors. Our additional sum-at-most-one constraint is explicit regularization.
- [limma 3.58.1, empirical-Bayes implementation](https://github.com/bioconductor-source/limma/blob/27e81a0a1ac6a00c94ec680e1f37138e0445f41d/R/ebayes.R), the source of the moderated t-statistic relation used in the uncertainty audit.

Code and fitted weights are distributed under the repository's
[Apache-2.0 license](../../LICENSE). [Third-party attribution](NOTICE.md) retains
the upstream sources and their license notices.

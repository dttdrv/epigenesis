# Measured DNA effects during neural induction

This external interpreter learns how engineered changes to human enhancer
sequences affect a reporter during neural induction. It predicts seven values
from two actual DNA sequences, compiles them into a Development Module, and
checks the linked tensor against the saved model and exact source bytes.
Only the Python standard library and the existing compiler are required.

The first frozen evaluation **failed predictive acceptance**. It improved on
predicting no effect, but did not establish an advantage over a simpler sequence
model. This is a working experiment with a preserved negative result.

## Run from the repository or extracted source archive

```sh
python3.11 tests/verify.py neural-regulation
run_dir=$(mktemp -d)
python3.11 examples/neural-regulation/regulation.py benchmark --output "$run_dir/benchmark"
```

The software gate checks source admission, numerical fitting, family isolation,
and compiler integration. It does not accept the scientific model. The benchmark
writes `report.json`, both fitted models, `data-quality.json`, `split.csv`, and
every held-out prediction in `test-predictions.csv`. `model.json` also records
all validation choices, input provenance, and the fitting implementation hash.
Existing output directories are rejected. A completed benchmark exits zero even
if `report.json` has `predictive_acceptance: false`; invalid inputs or execution
failures exit nonzero. No download occurs during either command.

To compile a prediction, prepare `reference.dna` and `variant.dna`, each with
exactly 171 uppercase A/C/G/T bytes and no newline or assay flanks:

```sh
python3.11 examples/neural-regulation/regulation.py predict \
  --model "$run_dir/benchmark/model.json" \
  --reference reference.dna --variant variant.dna \
  --output "$run_dir/prediction"
```

The result contains the exact two-record FASTA, copied model, source closure,
caller interpretation, typed target, Development Module, and validation report.
Its constant `effect` tensor has seven `f64` values ordered by hours
`[0, 3, 6, 12, 24, 48, 72]`. The consumer verifies the compiler closure, model
identity, source layout, target, operation, count, and recomputed values.
Reading a saved prediction does not modify its directory. A valid compiler
bundle containing substituted effects is rejected by sequence-model replay.

## Source and exclusions

[Kreimer, Ashuach, Inoue et al. (2022)](https://doi.org/10.1038/s41467-022-28659-0),
*Massively parallel reporter perturbation assays uncover temporal regulatory
architecture during neural differentiation*, supplied the sequence design and
[GSE188264](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188264) activity
estimates. The study used a lentiviral reporter assay in H1 human embryonic stem
cells undergoing neural induction. These are engineered enhancer perturbations,
not measurements of endogenous target-gene expression.

`data/provenance.json` identifies the original URLs, lengths, and SHA-256 hashes.
The two bundled files retain their exact source bytes. Source material is
attributed to the authors under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/);
the analysis and software are separate. Each 201-base design has two verified
15-base assay flanks around its 171-base insert.

The intake records 10,041 physical design records, collapses eight exact
duplicate records, and excludes four IDs with conflicting sequences. It pairs
variants and controls to their exact parent locus and excludes 116 activity
rows lacking a measured parent. The resulting 9,241 variants cover 585 loci.
Overlapping loci and identical/reverse-complement insert sequences form 582
connected families: 352 training, 112 validation, and 118 test. Neither measured
significance nor motif activity labels select the rows. Zero activity estimates
are retained. Full exclusion IDs appear in `data-quality.json`.

## Analysis fixed before fitting

Each connected family's sorted locus list is hashed with the namespace
`epigenesis-neural-regulation-v1` followed by a newline. Digest modulo 10 assigns
0–5 to training, 6–7 to validation, and 8–9 to test. All variants, controls, and
timepoints from a family stay together. Sequence names and coordinates are used
for pairing and splitting only; they never enter the predictor.

For a training-defined positive scale `s`, transform each activity estimate `a`
as `g(a) = asinh(a/s)/ln(2)`. The scale is the median positive WT estimate across
training loci and times. The target is `g(variant) − mean(g(parent))` at each
timepoint. Duplicate WT labels are averaged only when their oriented sequences
and locus match. This transformed difference is defined at zero and is not
literally a log2 fold change. The saved model specifies its units through `s`.

Features are variant-minus-reference counts of reverse-complement-canonical
k-mers. A no-intercept ridge model minimizes mean squared error over families,
variants within each family, and seven timepoints, plus lambda times the
squared coefficient norm. Thus each family receives equal total fitting weight.
The seven output heads share sequence features but have separate coefficients;
this is a statistical sequence model, not a regulatory dynamics mechanism.

Validation on targeted variants chooses `k ∈ {3,4,5}` and
`lambda ∈ {0.01,0.1,1,10}`. Ties prefer smaller k and larger lambda. A k=2
baseline selects its own lambda using the same procedure. Both models are
refit on training plus validation, including a new scale from those families,
before the test predictions are evaluated once.

Primary acceptance requires the lower endpoint of a paired 95% family-bootstrap
MSE improvement interval to be positive against **both** zero effect and the
dinucleotide baseline on targeted perturbations. The bootstrap resamples whole
families 2,000 times with seed 20261004. RAND and SCRAM controls are scored
separately; neither is presumed biologically neutral. Time-centered errors and
a time-averaged prediction ablation are separate diagnostics. Shuffling whole
seven-value training target vectors with seeds 19, 73, and 211 provides three
diagnostic controls, not a calibrated significance test.

## First frozen result, 2026-10-04

The selected model used k=5 and lambda=0.01. Selection-training scale was
1.6251197955; refit scale was 1.61952825. On 1,415 targeted variants from 118
held-out families:

| Predictor | Family-weighted MSE | Improvement over this baseline, 95% interval |
|---|---:|---:|
| Sequence model | 0.285185 | — |
| Predict no effect | 0.360870 | [0.026133, 0.133503] |
| Dinucleotide model | 0.319358 | [−0.005632, 0.082658] |

The 20.97% MSE reduction versus zero effect does not satisfy the second
comparison. **Predictive acceptance is false.** Time-specific heads did not
establish a benefit over their time average: improvement interval
[−0.000736, 0.000796]. Their time-centered MSE was 0.032119, versus 0.032108 for
the averaged prediction. RAND point-estimate MSE was also worse than both
baselines: 0.149400 versus 0.122320 for zero effect and 0.129871 for dinucleotides.
Shuffled-target MSEs were 0.315996, 0.321234, and 0.325407; these do not establish
mechanistic validity. SCRAM MSE was 0.561121, versus 0.782811 and 0.635540.

The original fitted model SHA-256 is
`bd3151cd7d7d055a8e647dcaee0a88eb1bf663f0aa0a9e27b821bd245018586e`.
The pre-fit analysis document SHA-256 is
`76beb57293854c4aa702acff173159eddb1de2cc7bf1690a2526bcfc2cc3e5c6`.
Local evidence is in `build/neural-regulation-20261004/`, including the original
fitting source and protocol. Subsequent software fixes handle numeric overflow
and keep the loader read-only; they do not retune the model or change its score.
New runs record their own implementation and artifact identities.

## What this can establish

The published MPRAnalyze estimates were pooled and normalized before our split.
This is retrospective prediction for unseen enhancer families conditional on
published preprocessing. It is not a test in an untouched biological
preparation, donor, or batch; the intervals do not include those uncertainties.

Predictions apply to this reporter assay and engineered sequence distribution.
The coefficients do not identify a causal mediator. No independent rescue,
natural variant validation, calibrated biochemical rate, cell-fate prediction,
growth, or circuit construction is demonstrated. In particular, the assay's
MPRAnalyze `alpha` is not the `alpha` production parameter in the separate
[neural-fate ODE example](../neural-fate/README.md).

The current test set has now been examined. A successor model needs selection
without these outcomes and fresh evaluation data for a confirmatory claim.
The [neural-mechanism acceptance protocol](../../docs/NEURAL_MECHANISM_ACCEPTANCE.md)
lists the additional evidence required for a developmental mechanism.

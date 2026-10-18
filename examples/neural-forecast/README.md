# Early-to-late forecasts across neural infection cultures

This archive tests whether early reporter measurements in two infected H1
cultures forecast later changes in a third, for excluded sequence families.
It contains the completed comparison, including every failed control comparison
and unsupported prediction. The fixed acceptance rule was not met.

## Result

The candidate uses the mean response at 0, 3 and 6 hours and the response at
12 hours. It predicts changes from 12 hours to 24, 48 and 72 hours. Fitting
excludes the target culture and the evaluated sequence families. The target
culture's 12-hour measurement contributes only to the observed change.

| Method | Family-balanced mean squared error |
|---|---:|
| Two-feature ridge | 0.32132339 |
| Early-mean ridge | 0.32194217 |
| Latest-value ridge | 0.32211810 |
| Zero change | 0.32214740 |
| Training-population mean change | 0.32201750 |

MSE is in squared log2 reporter-contrast units. The candidate improves on zero
change by 0.26% in the point estimate. All four simultaneous nominal 99% gain
intervals include zero. The early-mean control also wins in culture 2.

All 351 families and all 15 culture/fold groups contribute to the comparison.
Of 11,860 input-eligible forecasts, 11,696 have complete held-culture targets;
164 remain explicitly unscored. Missing observations are never replaced by zero.

## Replay

Use Python 3.12+ and the pinned NumPy dependency:

```sh
python -m pip install -r examples/neural-forecast/requirements.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python examples/neural-forecast/replay.py
python -m unittest discover -s examples/neural-forecast -p 'test_*.py' -v
```

Replay verifies the archive hashes, reconstructs the observable from saved count
totals, checks culture/family ownership and applies the saved model coefficients.
It independently recomputes all losses and 20,000 whole-family bootstrap samples.
It uses no GPU, network, model fitting or new random model draws. Bootstrap
resampling uses the fixed seed. The original and replay runtimes are recorded.

Large JSON records are gzip-compressed without changing their decoded bytes.
[SUMMARY.json](SUMMARY.json) keeps the primary results readable on GitHub.
The archive can be copied outside the repository and replayed there. An optional
`--output /path/to/new-receipt.json` writes a receipt outside the archive. Run
without Python's `-O` option because replay uses checked assertions.

## Source and scope

The source is [Kreimer et al., 2022](https://doi.org/10.1038/s41467-022-28659-0),
with deposited matrices and annotations in
[GSE188264](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188264).
[sources.json](sources.json) records exact download identities. Raw matrices
and pretrained model weights are not needed for this calculation replay.

The observable is a new log2 contrast of pooled RNA/DNA counts against the
matched parent. Only same-sample barcode slots with positive DNA and RNA count
contribute. Sample-depth factors cancel. This differs from the authors' original
estimator and from the earlier pooled-alpha temporal experiment. Those MSEs and
relative improvements cannot be compared as a matched before/after test.

The three omissions hold out infection cultures, not donors. Two viral batches
were used; their assignment is unavailable. Original cohort selection used
published pooled activity availability, and upstream matrix selection remains
unresolved. This is a retrospective comparison conditional on supported
measurements, not prospective biological confirmation.

[PLAN.md](PLAN.md) preserves the specification frozen before numerical count
conversion. [original/](original/) contains the unchanged analysis code for
inspection; rerunning it needs the original source layout and pinned dependencies.
The supported portable command is the saved-calculation replay above. Per-row
support, individual parent totals, fitted coefficients, prediction seals and
losses are retained in [result/](result/).

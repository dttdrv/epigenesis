# Frozen Salomon-trained activity-transfer model on six neural reporter studies

Written 2026-10-19 before any activity-transfer prediction was compared with an outcome. The model is bundled `examples/neural-predictor` activity-transfer predictor, applied unchanged: a k5 baseline (512 coefficients) plus a linear correction from a 1,024-coefficient activity score learned on 30,800 Salomon WTC11-NGN2 reporter measurements. **Nothing is fitted on any of the six studies**; the test is whether a model trained on a different assay and cell state predicts allelic effects elsewhere.

Models: A4 = baseline + correction; Abase, Aact = its two parts; S2 = the existing Siraj-trained head; Ceq5 = z(A4)+z(S2), fixed equal weights (no tuning). Inputs are the same held-out-eligible variants and families as the six-study test; k-mer deltas are length-independent so the 150/200-base windows need no padding.

Scoring: identical to the six-study test (family-weighted r, 2,000-draw whole-family bootstrap, PCG64 20261019, 5-95% intervals, signed utility 50/200). Primary group is all six studies. Prespecified secondary group: the five studies without gse225815, because a count-depth check shows its allelic effect is mostly counting noise (83% of its variants have Poisson noise variance above 0.0625 log2^2; the other three count-based studies are below 0.01 at the median).

Checks (lower 5% bound > 0): G1 r(A4); G2 r(A4) - r(S2); G3 r(Ceq5) - r(S2). Label `activity_transfer_transfers_to_neural_reporter_assays` needs G1 and G2.

Limits: all six outcomes were seen in earlier tests (earlier three not blind). Scorer is qualified on a synthetic positive and null first; scored once.

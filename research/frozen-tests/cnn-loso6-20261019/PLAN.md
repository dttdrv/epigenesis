# Pooled-effect CNN, leave-one-study-out (frozen 2026-10-19, before any CNN prediction)

Question: does a nonlinear sequence model trained on all other effect measurements predict allelic effects in a held-out neural reporter study better than the Siraj-trained k-mer head S2?

**Model N1.** A scalar sequence scorer g (one-hot, 151-base window centred on the first reference/alternate difference, zero-padded; conv 4→64 k15, ReLU, conv 64→64 k9 dilation 2, ReLU, concatenated global max and mean pool, linear to 1; averaged over the sequence and its reverse complement). Predicted effect = g(alt) − g(ref). About 45k parameters.

**Training pool for held-out study H.** Siraj FIT edits (66,604; family weights; validation = Siraj TUNING, 16,986) plus every variant of the five other studies, with gse225815 variants kept only if their counting-noise variance V ≤ 0.0625 log2² (a count-depth rule, not an outcome rule). Outcomes are centred and scaled per dataset (family-weighted). Per-dataset total loss weight: Siraj 2, each other study 1. Loss: family-weighted squared error. AdamW lr 2e-3, weight decay 1e-2, batch 256, at most 10 epochs, two seeds (1, 2). The checkpoint with the lowest Siraj-TUNING family-weighted MSE is kept; the held-out study is never opened during training (asserted). N1 = mean of the two seeds' predictions.

**Comparators.** S2 (existing Siraj-trained k5 + motif head); Ceq6 = z(N1) + z(S2), fixed equal weights.

**Scoring.** Same held-out populations, families, family-weighted Pearson r, 2,000-draw whole-family bootstrap, PCG64 20261019, 5–95% intervals, signed utility at 50/200. **Primary group: the five studies without gse225815** (its effect is mostly counting noise); all six reported too. Checks (lower 5% bound > 0): G1 r(N1); G2 r(N1) − r(S2); G3 r(Ceq6) − r(S2). Label `pooled_cnn_transfers_to_neural_reporter_assays` needs G1 and G2.

**Limits.** All outcomes were seen in earlier tests (earlier three not blind; held-out outcomes are not used in training). Scorer qualified on a synthetic positive and null, then run once. Mac MPS/CPU only.

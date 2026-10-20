# Zero-shot expressed-TF binding model (scored once, 2026-10-19)

Plan `PLAN.md`, sealed in `FREEZE.json` before scoring. Nothing fitted: JASPAR 2024 motif occupancy change at the edited base, weighted by BrainSpan 8–13 pcw neocortex expression of the transcription factor, motif gain assumed to activate.

| Five studies without gse225815 | family-weighted r (5–95%) |
|---|---|
| Z1 expression-weighted binding | 0.0015 (lower bound below 0) |
| U1 same, all motifs equally weighted | −0.0245 |
| S2 Siraj-trained head | 0.0465 |
| Z1 with permuted weights (200 draws) | mean −0.0243, 95th percentile −0.0199 |

Gates: G1 (r > 0) **fails**; G2 (beat S2) fails (Z1 − S2 = −0.045); G3 (above permuted-weight null) passes only because the null is biased negative; G4 (beats unweighted) passes. Label `expression_weighted_binding_mechanism_predicts_reporter_effects` **not earned**. Reading: motif gain alone anti-correlates slightly with reporter effect; expression weighting removes that bias (+0.026) but leaves no predictive power. Per-study r: lee 0.008, mcafee −0.044, matoba 0.039, gse221384 −0.035, weiss 0.040, gse225815 0.001. Scorer qualified on a synthetic positive (label true) and null (false). Files: `PLAN.md`, `FREEZE.json`, `scripts/{features,score11}.py`, `derived/{PREDICTIONS.npz,RESULT.json,weights.json}`.

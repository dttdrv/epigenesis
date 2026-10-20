# Sei chromatin-profile features as a second model (scored once, 2026-10-19)

Question: does a model of a different kind (Sei sequence-to-chromatin predictions
on 4,096 bp of context, 40 sequence-class scores + 571 neural/brain profile
scores, ridge-trained on the other studies) fix what the k-mer/motif models miss?

Result: **no**. The Sei ridge (SeiR) has no signal on held-out variants
(r = -0.003, 5-95% interval -0.024 to 0.019) and is worse than the
Siraj-trained k-mer head S2 (r = 0.047; SeiR minus S2 = -0.050, interval
-0.078 to -0.021). The inner-chosen blend of Sei with S2 (Cfit2) gives
r = 0.041 and does not beat S2 alone (-0.006, interval -0.011 to -0.0005).
Label `sei_combination_transfers_conditional_on_neural_reporter_assays`
not earned (G2 false). Cfit2 beats SeiR (+0.044) and r > 0, but only because
SeiR is empty.

Per-study r (SeiR / S2): lee -0.011/0.019, mcafee -0.031/0.091, matoba -0.031/0.023,
gse225815 0.044/0.049, gse221384 -0.029/-0.014, weiss 0.041/0.115.
Noise ceilings are 0.60-0.98, so measurement noise is not the limit.

Scorer qualified first: synthetic positive passes all three comparisons, null
passes none. Limits: all six outcomes were used in earlier tests (earlier three
not blind); Sei forward strand only; context models not included.
Files: PLAN.md, FREEZE.json, derived/RESULT.json.

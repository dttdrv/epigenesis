# Salomon-trained activity-transfer model on six neural reporter studies (scored once, 2026-10-19)

Result: **no transfer**. The bundled model (k5 baseline + activity correction, trained only on Salomon WTC11-NGN2 data, applied unchanged) gives r = -0.010 (5-95% interval -0.029 to 0.011) against the Siraj-trained head's 0.047; A4 minus S2 = -0.057 (-0.084 to -0.029). Its activity part alone is -0.018. Averaging it with S2 (Ceq5, r = 0.024) is worse than S2 alone (-0.023, interval -0.036 to -0.008). Without gse225815: A4 -0.004, S2 0.047. Label not earned.

Utility at 50/200 edits: A4 -0.019/-0.018, S2 0.051/0.034.

## Correction to earlier reading of the six-study results
A count-depth check on raw counts found that gse225815's allelic effect is mostly counting noise: for 83% of its variants one allele has under ten RNA reads, so the log-ratio is set by the 0.5 pseudocount and DNA-count imbalance. That study's high replicate-half agreement (0.98) is the shared DNA count, not reproducible biology. On its 124 well-measured variants the Siraj-trained head scores r = -0.09; its apparent positive r comes from the noisy variants. The other three count-based studies have negligible counting noise (median Poisson variance under 0.002 log2^2). So "positive in 5 of 6 studies" overstated the evidence, and gse225815 should be dropped from future primary groups (secondary group `five_without_gse225815` here). The quoted noise ceilings for studies that share one DNA count across replicate halves are upper bounds, not proof that the models are the limit.

Where the models do carry signal: the Siraj-trained head's sign is right for 53-57% of the top-10% largest effects and 58-65% of the top 1% in lee, matoba and weiss (50% overall) — real but small, and exploratory (outcome-selected strata).
Files: PLAN.md, FREEZE.json, derived/RESULT.json; `precision-loso6-20261019/scripts/noise.py` computes the per-variant noise.

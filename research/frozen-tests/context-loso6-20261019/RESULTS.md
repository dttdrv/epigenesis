# motif-context model on six neural reporter studies (scored once, 2026-10-19)

Result: **no gain**. The Siraj-trained context head (S3) carries the same small signal as the existing head (S2): r = 0.046 vs 0.047 on held-out variants; S3 minus S2 = -0.001 (5-95% interval -0.005 to +0.002). Trained on the other six studies instead (H3), the same features give r = 0.020 (-0.005 to 0.043). The inner-chosen blend (Cfit3) is 0.044 and does not beat S3 (-0.002, interval -0.006 to +0.002). Neither label earned (G1 passes only because S3 matches S2; G2 and G3 fail).

Per-study r (S2 / S3 / H3): lee 0.019/0.013/0.023, mcafee 0.091/0.095/0.073, matoba 0.023/0.026/-0.012, gse225815 0.049/0.047/0.006, gse221384 -0.014/-0.020/0.017, weiss 0.115/0.112/0.011. Signed utility at 50/200 edits is also no better than S2 (S3 0.033/0.027 vs S2 0.051/0.034).

So 1.8-2.1% squared-error gain on Siraj edits (reproduced here at 2.1%) does not carry across to other neural reporter assays: the quadratic context term adds nothing out of assay. Noise ceilings are 0.60-0.98, so the models are the limit.

Scorer qualified first (synthetic positive passes all checks, null none). Features reproduce the reference cache to 1e-14. Limits: outcomes of all six studies were used in earlier tests (earlier three not blind); S3 reuses exploratory recipe; smooth-context model not included. Files: PLAN.md, FREEZE.json, derived/RESULT.json.

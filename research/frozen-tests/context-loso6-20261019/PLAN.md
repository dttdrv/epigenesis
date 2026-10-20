# motif-context model on the six neural reporter studies

Written 2026-10-19, before any context-model prediction exists. The context model is `siraj-context-pwm-20261019`: 532 log-potential motif deltas d = Q(alt) − Q(ref) plus 532 squared-context terms d·(Q(alt)+Q(ref)) (JASPAR 2024, scale .003). On Siraj's own held-out edits it beat the k-mer baseline by 1.8% in squared error. This tests whether it carries over to the six independent neural reporter studies (lee, mcafee, matoba, gse225815, gse221384, weiss), which every earlier head failed (best r ≈ 0.05).

## Models
- **S2** (existing): Siraj-trained k5 + signed-log motif head, predictions already sealed in `hnp-loso6-20261019`.
- **S3**: Siraj-trained, recipe: the same k5 head, plus a ridge on the k5 residual over the 1,064 context features (FIT+TUNING pool, equal family weights, penalty 3 as fixed in the original plan, no search). Features reproduce the reference cache to 1e-14. Nothing is tuned on any study.
- **H3**: ridge on [k5, context] trained on the other five studies (outcomes centred/scaled per study, penalty by inner leave-one-study-out over {0.1…1e5}).
- **Ceq3** = z(H3)+z(S3); **Cfit3** = w·z(H3)+(1−w)·z(S3), w from {0,.25,.5,.75,1} by inner leave-one-study-out.

## Test
Same held-out populations, families, family-weighted Pearson r, 2,000-draw whole-family bootstrap (PCG64 20261019), 5–95% intervals, signed utility at 50/200, and noise ceilings as the six-study test. Primary group is all six; study-equal mean.

Checks (lower 5% bound > 0): G1 r(S3); G2 r(S3) − r(S2) (does the context help beyond the existing head); G3 r(Cfit3) − r(S3).
Labels: `context_head_transfers_to_neural_reporter_assays` needs G1+G2; `context_combination_transfers_to_neural_reporter_assays` needs G1+G2+G3. Otherwise not earned.

## Limits
All six outcomes were used in earlier tests, so this is not a fresh blind test (earlier three are not blind; the first scored here for this model). Context features inherit adaptive history on Siraj (exploratory). Scorer is qualified on a synthetic positive and null before scoring; scored once; no reruns or tuning after seeing the result. smooth-context model is not included.

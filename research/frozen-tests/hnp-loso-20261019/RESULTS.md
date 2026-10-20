# Leave-one-study-out across human neural-progenitor reporter studies: result

19 October 2026. Frozen plan: `PLAN.md`; hashes in `FREEZE.json`. Predictions were sealed before any held-out outcome was
scored, and the scorer passed a synthetic positive and a synthetic null first. Scored once.

**Fail.** Models trained on two HNP variant studies do not predict the allelic effect of unseen variants in the third.
Not `cross_study_transfer_conditional_on_hnp_reporter_assays`.

Held-out variants (no variant of the other studies or Siraj within 2,000 bases): Lee 6,876 (2,588 families), McAfee 1,867
(736), Matoba TRAIN 4,055 (97 loci). Family-weighted Pearson r, McAfee and Matoba folds (the primary group):

| Model | r | 5-95% bootstrap |
|---|---:|---|
| H2 k5 + motif, trained on the two other studies | 0.015 | -0.024 to 0.054 |
| H1 k5, trained on the two other studies | -0.014 | -0.047 to 0.018 |
| S1 k5, trained on Siraj only | 0.060 | 0.025 to 0.095 |

G1 (H2 bound > 0) false. G2 (H2 minus S1 bound > 0) false; gain -0.045 (-0.098 to 0.005). G3 false. Per-study r for
H2: Lee -0.007, McAfee 0.045, Matoba -0.015. Signed utility at 50 selections: H2 -0.011, S1 +0.030.

Noise ceilings (Spearman-Brown): Lee 0.83, McAfee 0.60, Matoba 0.70. The outcomes are reproducible, so the models are the
limit. The inner check between training studies already showed r of about 0 to 0.05 for every penalty.

Observation, not a gate: the Siraj k5 model has a small positive r on McAfee and Matoba (0.060, bound above zero). It
explains well under 1% of the variance and was not a pre-specified claim.

Limits: the three studies share a cell type; endpoint is reporter activity; Lee outcomes were seen in the earlier
cross-lab test, so the Lee fold is not blind; the penalty grid was widened once, from 1e3 to 1e5, after training-study
results only (recorded in `FREEZE.json`); 30 McAfee and 29 Matoba rsIDs had no GRCh37 position and were dropped.

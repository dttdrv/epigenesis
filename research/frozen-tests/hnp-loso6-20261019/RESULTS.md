# Six-study leave-one-study-out test of a combined predictor: result

19 October 2026. Frozen plan: `PLAN.md`; hashes in `FREEZE.json`. Predictions were sealed before scoring; the scorer's
synthetic null failed every gate and its synthetic positive passed G1 and G3 (G2 was false there because the invented
truth was built from the combination, which already matched the Siraj-trained model). Scored once.

**Fail.** The combination does not beat the Siraj-trained model, so the label
`combination_transfers_conditional_on_neural_reporter_assays` is not earned.

Held-out variants: Lee 5,916, McAfee 1,682, Matoba 4,013, GSE225815 611, GSE221384 748, Weiss 9,247.
Family-weighted Pearson r, primary group (the three new studies, study-equal mean), 5-95% bootstrap:

| Model | r | interval |
|---|---:|---|
| Cfit, fitted blend of HNP-trained and Siraj-trained | 0.046 | 0.009 to 0.082 |
| S2, Siraj-trained k5 + motif | 0.050 | 0.015 to 0.085 |
| H2, trained on the other five neural studies | 0.005 | -0.035 to 0.043 |

G1 true. G2 false: Cfit minus S2 = -0.004 (-0.012 to 0.003). G3 true: Cfit minus H2 = 0.041 (0.006 to 0.077).
The inner leave-one-study-out choice gave the HNP-trained model weight 0 to 0.25 in every fold, so Cfit is close to S2.
Over all six studies: S1 (Siraj k5) r = 0.047 (0.026 to 0.069), H2 0.013, Cfit 0.045.

Per-study r for S2: Lee 0.019, McAfee 0.091, Matoba 0.023, GSE225815 0.049, GSE221384 -0.014, Weiss 0.115.
Noise ceilings (Spearman-Brown): Lee 0.83, McAfee 0.60, Matoba 0.71, GSE225815 0.98, GSE221384 0.82; Weiss has no replicates.

Reading: the Siraj-trained model carries a small, consistent cross-assay signal (positive in five of six studies,
r about 0.05, under 1% of variance); training on the other neural studies adds nothing. The measurements are reliable, so
the models are the limit.

Limits: reporter activity, not neural development; some studies share cell source or laboratory lineage; Weiss is
modern-versus-archaic alleles; the earlier-three folds are not blind; Salomon 2026 neuron data were excluded because the
variant design and tables are restricted on the IGVF portal; GSE225815 lost about 1,000 sites at the reference-base check.

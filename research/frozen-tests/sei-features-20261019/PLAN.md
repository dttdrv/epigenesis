# Sei chromatin-model features as a second model, six-study leave-one-study-out

Preanalysis plan, 19 October 2026. Freeze this file, the scripts and the sealed predictions before any outcome is scored.

## Question

The Siraj-trained k-mer/motif model (S2) carries a small cross-assay signal (r about 0.05) and the model trained on the
other neural reporter studies adds nothing. Does a model of a different kind, built on Sei (Chen et al. 2022, trained on
21,907 chromatin profiles including neural progenitors, neurons and fetal brain), predict allelic effects in the held-out
study, and does combining it with S2 beat S2 alone?

## Model

- Sei published weights (Zenodo 4906997), run unchanged on forward-strand 4,096-base hg19 contexts: the study's own allele
  sequence is written over the centre of the genomic context, variant at index 2048. No reverse-complement averaging.
- Features per variant (alt minus ref): the 40 histone-normalised sequence-class scores (sei-framework utils) and the
  probabilities of a neural/brain profile subset fixed by cell-type name (regex in `scripts/run_sei.py`, chosen from
  `target.names` only; Kidney and Bone_Marrow lines removed).
- SeiR: ridge on those features, trained on the other five studies (outcomes centred and scaled per study, equal study and
  family weights), penalty grid {0.1 ... 1e5} chosen by inner leave-one-study-out over the five training studies.
- S2: the Siraj-trained k5 + motif head, as in the six-study plan (no neural data).
- Cfit2 = w z(SeiR) + (1-w) z(S2), z = division by the SD over the held-out set; w in {0, .25, .5, .75, 1} chosen by inner
  leave-one-study-out r (SeiR refit on four studies per inner fold). Ceq2 = z(SeiR) + z(S2) as a fixed-weight reference.
- Held-out populations, families and training studies are exactly those of the six-study test (`hnp-loso6-20261019`).

## Scoring

Family-weighted Pearson r per fold; whole-family bootstrap, 2,000 draws, PCG64 seed 20261019, one-sided 95% bounds.
Primary group: all six studies, study-equal mean. Also reported: the three studies first scored in the six-study test, the
earlier three, per study, signed utility at 50 and 200.

- G1: lower bound of r(Cfit2) > 0.
- G2: lower bound of r(Cfit2) - r(S2) > 0 (Sei adds to the model it is combined with).
- G3: lower bound of r(Cfit2) - r(SeiR) > 0 (S2 adds to Sei).
Label `sei_combination_transfers_conditional_on_neural_reporter_assays` requires G1, G2, G3. Anything else is a failure.
The scorer is qualified on a synthetic positive and null first. Nothing is rerun.

## Limits

Reporter activity, not neural development. All six studies' outcomes were used as training targets or scored in earlier
tests, so no fold is blind to the data; the Sei features and the plan were fixed before any Sei prediction was compared with
an outcome. Sei is run on forward strand only, to fit the CPU budget. motif-context models are not included.

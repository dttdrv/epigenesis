# Predicting the direction of fine-mapped eQTLs in developing human brain (frozen 2026-10-19)

Question: from DNA sequence alone, can a model trained on **adult** brain eQTLs predict which way a natural variant moves gene expression in **fetal neocortex** (and in neural progenitors and neurons), better than simple alternatives? This is the natural-variant, independent-data form of the decisive goal; the target is gene-expression direction in developing brain, not neural structure.

## Data (eQTL Catalogue release 7 SuSiE credible sets, GRCh38; downloaded 2026-10-19)
- **Test:** Walker_2019 fetal neocortex (QTD000579) = primary (792 variant-gene pairs); Aygun_2021 neural progenitors (QTD000564) + neurons (QTD000569) = secondary (309 pairs).
- **Training:** adult brain gene-expression datasets: GTEx (13 brain regions), BrainSeq, CommonMind, ROSMAP, Braineac2 (14,625 rows).
- Variants: single-base, PIP ≥ 0.5, gene on the same chromosome with a GENCODE v39 TSS. Label = sign of beta (alt allele raises expression = 1). Nothing else of the test outcomes is used.
- **Split by chromosome:** a test pair on chromosome c is scored by a model trained only on training pairs from other chromosomes, so no locus is shared. Training rows are weighted 1/(number of rows of the same variant-gene pair).

## Models (score = predicted probability that the alt allele raises expression)
- **M4:** L2 logistic regression on Sei alt−ref scores (40 histone-normalised sequence classes + 571 neural/brain profiles; 4,096-base context, forward strand), each block repeated for three strand-aware TSS-distance bins (|d| ≤ 2 kb, 2–50 kb, > 50 kb). Features standardised on training rows; penalty from {0.01…10^4} by leave-chromosome-group-out CV (accuracy then log-loss tie-break) inside the training chromosomes.
- **M5:** M4 plus GC change, local k5 delta (150-base window) and the Siraj-trained reporter head S2 score.
- **Baselines:** B0 = chance (0.5); B1 = alt allele gains GC (G/C vs A/T) predicts up; B2 = L2 logistic on the 150-base k5 delta (same chromosome split); B3 = sign of the Siraj-trained reporter head S2 (a model trained on a different assay, no fitting here).

## Scoring (scored once)
Sign accuracy and AUROC on each test set; **5–95%** whole-gene bootstrap intervals (2,000 draws, PCG64 20261019; pairs sharing a gene are resampled together). Checks on Walker (lower 5% bound > 0): G1 accuracy(M4) − 0.5; G2 accuracy(M4) − accuracy(best of B1, B2, B3 by Walker point value — a deliberately strict comparator); G3 same for M5. Label `fetal_brain_eqtl_direction_predicted_from_sequence` needs G1 and G2 on Walker, and accuracy(M4) > 0.5 point value on both Aygun sets.

## Limits
Fine-mapped variants are mostly shared across brain tissues, which makes adult → fetal an easier transfer than a new cell type; shared causal biology is the point. SuSiE PIP ≥ 0.5 can still pick a tag of the causal variant. Direction of expression is a molecular readout, not neural structure. Adaptive: the design (features, gating) was chosen after seeing the table sizes only, not any accuracy. Scorer is qualified on a synthetic positive and null before scoring. No fitting, tuning or re-run after scoring.

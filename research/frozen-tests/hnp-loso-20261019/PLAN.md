# Leave-one-study-out prediction across human neural-progenitor reporter studies

Preanalysis plan, 19 October 2026. Freeze this file, the scripts and the sealed predictions before any held-out
outcome is opened.

## Question

Do models trained on the allelic effects of natural variants in two human neural-progenitor (HNP) reporter studies
predict the allelic effect of unseen natural variants in the third study, better than simpler alternatives, including
the Siraj-trained models that already failed on Lee?

Studies (all human neural progenitors, hg19, 150-base alleles):

- Lee 2025, GSE244011: AAV-MPRA, psychiatric-risk variants, 22 replicates.
- McAfee 2023, GSE211045: lentiMPRA, schizophrenia-associated variants, 10 replicates (rsID tagged; positions from Ensembl GRCh37).
- Matoba, TRAIN split only: Vehicle (unstimulated) allelic effect of brain-morphology GWAS variants, two half contrasts.
  Validation and test pairs stay closed.

This is a retrospective test on public data. The endpoint is reporter activity, not neural development. The three
studies share a cell type and, probably, laboratory lineage, so they are partially independent. Lee outcomes were seen
in the earlier cross-lab test; the Lee fold is therefore reported but is not part of the primary gate.

## Outcomes

Alt minus ref. Lee and McAfee: pooled-count log2((RNA+0.5)/(DNA+0.5)) contrast over all replicates. Matoba: mean of the two
Vehicle half contrasts (allele2 minus allele1, allele1 = reference). Odd/even replicate halves give a noise ceiling
(Matoba: the two halves). Units differ between studies, so every training outcome is centred and scaled to unit
family-weighted SD within its study before pooling, and the metric is scale-free.

## Population (inputs and DNA counts only)

Eligible: single-base REF/ALT windows at index 75; DNA total >= 1,000 per allele; Lee also drops the author exclude list;
variants whose rsID has no unique GRCh37 position are dropped (McAfee, Matoba). Families: 2,000-base connected components
(Lee, McAfee) or the Matoba genomic-locus hash.
Held-out set for a fold = eligible variants of that study with no variant of either other study or any Siraj design within
2,000 bases. The set is the same for every model. Training uses all eligible variants of the other two studies.

## Models (fixed before any held-out outcome)

- H0 zero effect.
- H1 k5: ridge on ALT-REF reverse-complement-collapsed 5-mer deltas, trained on the two other studies.
- H2 k5 + signed-log motif: ridge on [k5, slog(0.003 x JASPAR potential delta)] jointly, same training.
- S1, S2 the Siraj-trained k5 and k5+motif heads from the Lee transfer test (no HNP data).
- C exploratory blend: H2 and S2 each divided by their held-out SD, summed.
Ridge penalty for H1 and H2: grid {0.1, 1, 10, 100, 1e3, 1e4, 1e5} (the first run, grid to 1e3, selected the edge for most models; extended before any held-out outcome was opened), chosen by inner leave-one-study-out between the two training
studies (fit on one, family-weighted r on the other's held-out-eligible variants). Equal study weights, equal family
weights within a study. No held-out outcome tunes anything; no model is added after predictions are sealed.

## Scoring

Per fold: family-weighted Pearson r (families resampled whole; 2,000 draws; PCG64 seed 20261019; one-sided 95% bounds).
Primary group: McAfee and Matoba folds, study-equal mean r. Also reported: Lee fold alone, all three, per study, signed
utility (one edit per family chosen by |prediction|, top 50 and 200, mean of sign x effect) and noise ceilings.

- G1: lower bound of r(H2) > 0.
- G2: lower bound of r(H2) - r(S1) > 0 (beats the cross-assay model).
- G3: lower bound of r(H2) - r(H1) > 0 (motif features add to k5 under the same training).
Label `cross_study_transfer_conditional_on_hnp_reporter_assays` requires G1 and G2. If G3 also holds, the motif
features are credited; if not, k5 alone carries the result. Any other outcome is a failure and reported with equal
prominence. The scorer is qualified first on a synthetic positive and a synthetic null. Nothing is rerun after scoring.

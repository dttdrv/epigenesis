# Six-study leave-one-study-out test of a combined predictor

Preanalysis plan, 19 October 2026. Freeze this file, the scripts and the sealed predictions before any held-out
outcome of a new study is opened.

## Question

Does a combination of two models that fail in different places predict the allelic effect of unseen natural variants
in a neural-cell reporter study better than either model alone? The two models are an HNP-trained model (fitted on the
other neural studies) and the Siraj-trained model (fitted on SK-N-SH edits, no neural-progenitor data). The three-study
test showed neither transfers well alone (r about 0.0 to 0.06).

## Studies (hg19; neural progenitors unless stated)

Earlier three, outcomes already seen in the three-study test, so they are reported but not primary:
Lee 2025 (GSE244011), McAfee 2023 (GSE211045), Matoba TRAIN split (Vehicle condition).
New, primary (first scored here):

- GSE225815 (schizophrenia variants, iPSC-derived NPC columns of the barcode table; lentiviral/AAV libraries Dev, Val, eQTL pooled).
  SNPs with exact barcode match only; DNA >= 100 per allele; reference base must match GRCh37.
- GSE221384 (neuropsychiatric disease variants, anterior and posterior NPC columns, plasmid DNA >= 500 per allele;
  one pair per alternative allele).
- Weiss 2021, GSE152404 (modern against archaic alleles in NPCs, 200-base oligos, processed CPM from the authors:
  single-base differences, DNA CPM >= 30 per allele). These are fixed evolutionary differences, not disease variants.

Salomon 2026 (IGVF neuron lentiMPRA, 46,000 variants) was located and downloaded in part, but the variant design file and
variant-level tables are access-restricted or revoked on the portal, so it is excluded.

## Outcomes

Alt minus ref log2 activity (RNA over DNA, +0.5 pseudocount; Weiss CPM with +1). Lee, McAfee, Matoba as in the
three-study plan. Odd/even replicate columns give noise ceilings where replicates exist (not for Weiss). Training
outcomes are centred and scaled to unit family-weighted SD within their study before pooling.

## Population

Single-base variants with DNA support as above. Families: 2,000-base connected components, except Matoba (genomic
locus hash). A study's held-out set is its eligible variants with no variant of any other study or any Siraj design
within 2,000 bases; it is the same for every model. Training uses all eligible variants of the other five studies.

## Models (fixed before any new outcome is opened)

- H1 k5 ridge and H2 k5+motif ridge, trained on the other five studies, equal study weights and equal family weights.
  Penalty grid {0.1, 1, ..., 1e5} chosen by inner leave-one-study-out over the five training studies.
- S1, S2: the Siraj-trained k5 and k5+motif heads from the Lee transfer test.
- Ceq: z(H2) + z(S2), each divided by its SD over the held-out set.
- Cfit (the combination under test): w z(H2) + (1-w) z(S2), w in {0, .25, .5, .75, 1} chosen by inner
  leave-one-study-out r over the five training studies, with H2 refit on four studies per inner fold.
No new-study outcome tunes anything; no model is added after predictions are sealed.

## Scoring

Per fold family-weighted Pearson r; whole-family bootstrap, 2,000 draws, PCG64 seed 20261019, one-sided 95% bounds.
Primary group: the three new studies, study-equal mean. Also reported: the earlier three, all six, per study,
signed utility (top 50 and 200), noise ceilings.

- G1: lower bound of r(Cfit) > 0.
- G2: lower bound of r(Cfit) - r(S2) > 0 (the combination beats the Siraj-trained model).
- G3: lower bound of r(Cfit) - r(H2) > 0 (the combination beats the HNP-trained model).
Label `combination_transfers_conditional_on_neural_reporter_assays` requires G1, G2 and G3. Anything else is a failure,
reported with equal prominence. The scorer is first qualified on a synthetic positive and a synthetic null. Nothing is
rerun after scoring.

## Limits

Endpoint is reporter activity, not neural development. Several studies share cell source or laboratory lineage. The
Weiss set is evolutionary alleles. The earlier-three folds are not blind.

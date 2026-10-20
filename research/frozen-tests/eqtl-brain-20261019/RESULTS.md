# Direction of fine-mapped eQTLs in developing human brain (scored once, 2026-10-19)

Task: from DNA sequence, predict whether the alt allele raises or lowers gene expression, for fine-mapped natural variants (SuSiE PIP ≥ 0.5). Trained on adult brain (GTEx brain regions, BrainSeq, CommonMind, ROSMAP, Braineac2; 14,617 rows) with chromosome-held-out folds; tested on fetal neocortex (Walker 2019, 792 pairs) and neural progenitors / neurons (Aygun 2021, 229 / 80 pairs).

| Fetal neocortex (n = 792) | Sign accuracy (5–95%) | AUROC |
|---|---|---|
| M4: Sei chromatin scores × TSS-distance bin | 0.547 (0.517–0.575) | 0.582 |
| M5: M4 + GC change + local k5 + reporter head | 0.566 (0.537–0.594) | 0.606 |
| B3: Siraj reporter-trained head (no fitting) | 0.538 (0.508–0.567) | 0.580 |
| B1: alt gains GC | 0.525 (0.497–0.551) | 0.532 |
| B2: k5 logistic | 0.520 (0.492–0.550) | 0.519 |

M4 is above chance (lower bound 0.517) but not clearly above the best baseline: M4 − B3 = +0.009 (−0.032 to +0.050); M5 − B3 = +0.028 (−0.009 to +0.065). Label `fetal_brain_eqtl_direction_predicted_from_sequence` **not earned** (G2 and G3 fail).
Aygun progenitors: M4 0.528, M5 0.646 (M5 − best baseline B2 = +0.066, 0.013 to 0.117); Aygun neurons: M4 0.575, M5 0.600 (n = 80, intervals wide).

What is new: this is the first test where a model trained on other data beats chance on held-out natural variants in developing brain. The effect is small (55–57% direction; AUROC 0.58–0.61), and the reporter-assay head alone already reaches 54%, so part of the signal is shared regulatory sequence logic. Scorer qualified first (synthetic positive passes, nulls hover around 0.5; one null draw reached 0.547 by chance among nine seeds). Files: PLAN.md, FREEZE.json, derived/RESULT.json.

## Follow-up 1: Jerber differentiating neurons (scored once, frozen in `PLAN_CONFIRM.md`)

320 unique fine-mapped pairs (PIP ≥ 0.5), same models, same training. M5 60.3% (5–95%: 55.9–64.6, AUROC 0.62); M4 53.4% (48.8–58.1); B2 (k5 logistic) 60.0% (55.5–64.6); B1 (alt gains GC) 44.5%; B3 (reporter head) 49.1%. M5 is above chance but M5 − B2 = +0.3 points (−3.7 to +4.5). Label `eqtl_direction_confirmed_in_differentiating_neurons` **not earned** (C3 fails, C2 fails).

## Follow-up 2: Schwartzentruber sensory neurons and Young microglia (scored once, frozen in `PLAN_REP.md`)

Asked whether the k5 model B2, above chance in Aygun and Jerber, repeats in a fourth set. It does not: sensory neurons (PIP ≥ 0.1, n = 519) B2 49.9% (46.1–53.7), M5 51.8% (47.9–55.6), M4 53.8% (49.8–57.8); the PIP ≥ 0.5 subset has 45 pairs and cannot decide. Microglia (n = 272): B2 51.8%, M5 57.0% (52.0–61.9), reporter head 57.4%. Pooled descriptive result over the four neural-derived sets at PIP ≥ 0.5 (626 pairs): M5 59.9% (56.6–63.0), B2 56.7%, B3 50.3%, M4 52.9%, B1 45.4%. Label `kmer_eqtl_direction_replicates_in_neural_derived_cells` **not earned** (R1 fails).

Summary of the eQTL line: adult-brain-trained sequence models give 50–60% direction accuracy on neural-derived natural variants, varying by set (fetal neocortex 52–57%, Aygun progenitors 53–65%, Jerber 53–60%, sensory neurons 50–54%). No model is consistently above the simplest k-mer baseline, and the best set-level numbers did not repeat in the independent set. Not a predictor claim.

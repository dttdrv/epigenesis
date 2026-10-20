# Measured fetal-brain accessibility for causal-variant ranking (second look, frozen 2026-10-19, before any model output)

`PLAN.md` found no model above nearest-TSS in fetal neocortex and conservation at chance. Both were *predicted* or evolutionary features. This adds the **measured** regulatory landscape of the right tissue: Roadmap Epigenomics DNase-seq peaks of fetal brain (E081, E082) and of H1-derived neuronal progenitors (E007), lifted from hg19.

**Same population, labels, training sets, chromosome rule, conditional logit and scoring as `PLAN.md`.** The Walker and neural-derived test sets are therefore a second look with new features (adaptive; the first look is on record).

**New features per member.** For each track (E081, E082, E007): 1 if the variant lies in a peak (variant lifted from GRCh38 to GRCh37 with the UCSC chain), and log1p of the signalValue of the overlapping peak (0 if none). Plus log10(1 + distance to the nearest peak of any track), capped at 100 kb.

**Models.** B1 nearest the TSS, B2 distance bins, M4 distance + Sei (as before, recomputed in this run); **M6 distance + Sei + accessibility (primary)**; M7 distance + accessibility only. Conservation is dropped (no signal on training sets: pairwise AUC 0.495, a training-data diagnostic).

**Scoring (scored once).** Top-1 accuracy, pairwise AUC, 5–95% bootstrap over genes, PCG64 20261019. On Walker (lower 5% bound > 0): G1 top-1(M6) − top-1(B1); G2 top-1(M7) − top-1(B1); G3 top-1(M6) − top-1(M4). Label `fetal_accessibility_plus_sequence_model_prioritises_causal_variants` needs G1 and a positive point value of top-1(M6) − top-1(B1) on the neural-derived pool.

**Limits.** Roadmap fetal brain is bulk tissue; accessibility peaks are broad, so a modest gain is the realistic outcome. Lifting loses a few variants (they get no peak). Expression-level, not neural-structure, readout. Scorer qualified on a synthetic positive and null before scoring; no tuning or re-run after scoring.

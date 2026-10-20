# Which variant in a fine-mapped block changes expression in developing brain? (scored once, 2026-10-19)

Plan `PLAN.md`, sealed in `FREEZE.json` before scoring; scorer qualified on a synthetic positive (label true) and null (false). Conditional-logit ranking of credible-set members, trained on 5,836 adult-brain sets (27,875 members), leave-chromosome-out.

| Top-1 accuracy (5–95%) | Walker fetal neocortex (439 sets) | Neural-derived pool (420 sets) |
|---|---|---|
| random | 0.355 | 0.328 |
| B1 nearest the TSS | **0.460** (0.420–0.499) | 0.402 (0.362–0.442) |
| B2 distance bins only | 0.433 (0.403–0.463) | 0.381 (0.351–0.411) |
| B3 highest phyloP | 0.355 (0.320–0.392) | 0.333 (0.294–0.373) |
| M4 distance + Sei | 0.435 (0.396–0.474) | **0.426** (0.386–0.464) |
| M5 distance + Sei + phyloP | 0.419 (0.381–0.458) | 0.421 (0.381–0.461) |

Walker: no model beats the nearest-TSS rule; M5 − B2 = −0.014 (−0.051 to +0.023), M5 − B1 = −0.041 (−0.082 to +0.002). G1, G2, G3 fail. Neural-derived pool: M4 − B2 = +0.045 (+0.005 to +0.085), M5 − B2 = +0.040 (−0.001 to +0.081), but neither clearly beats B1 (M5 − B1 = +0.019, −0.029 to +0.067). Label `sequence_model_prioritises_causal_variants_beyond_distance` **not earned** (G1 fails on Walker).

Conservation carries no signal here: phyloP picks the top-PIP member no better than chance (pairwise AUC 0.495 on training sets, 0.503 on Walker), so the "matters from conservation" idea fails in natural-variant data as it did in reporter data. Sei features add a few points over distance in neural-derived cells but not in fetal neocortex, the same pattern as the direction tests. Training-fold CV top-1 (training data only): B2 0.381, M4 0.388, M5 0.392.

Files: `scripts/{build_cs,fetch_ctx_cs,run_sei_cs,phylop_cs,run_cs,score_cs}.py`, `derived/{cs_table.json,PREDICTIONS.npz,RESULT.json}`.

## Second look: measured fetal-brain accessibility (scored once, `PLAN2.md`, `FREEZE2.json`)

Added Roadmap DNase peaks of fetal brain (E081, E082) and H1-derived neuronal progenitors (E007), lifted to GRCh37, to the same ranking task. Walker (439 sets): nearest-TSS 46.0%, distance + Sei 43.5%, distance + Sei + accessibility (M6) 41.7%, distance + accessibility (M7) 42.0%; M6 − nearest-TSS = −4.3 points (−8.9 to 0.0). Neural-derived pool (420 sets): M7 43.1%, M4 42.6%, M6 40.2%, nearest-TSS 40.2%, distance bins 38.1%; M6 − nearest-TSS = 0.0 (−4.6 to +4.9). All four gates fail; label `fetal_accessibility_plus_sequence_model_prioritises_causal_variants` **not earned**. Training-fold CV (training data only) showed peaks adding about one point over distance (B2 0.381, M7 0.394), which did not carry to the test sets. Files: `scripts/{peaks_cs,run_cs2,score_cs2}.py`, `derived/{PREDICTIONS2.npz,RESULT2.json,peaks_cs.npz}`.

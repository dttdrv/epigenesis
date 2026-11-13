# Frozen Held-Out Test-Split Evaluation: V2F Tipping-Point Radar on GSE188264

Written BEFORE any unseen `test`-split locus or perturbation outcome (`r.alpha`) is read or scored.

## Purpose

Condition 5 of the PR #4 pre-merge audit requires writing the V2F Tipping-Point & Rescue Radar's predictions on held-out GSE188264 loci **before scoring them**, matching the hashed freeze protocol in `research/frozen-tests/`.

## Held-Out Cohort

- **Source dataset:** `examples/neural-regulation/data/design.fasta` and `alphas.csv.gz` (Kreimer et al. 2022 H1 neural-induction MPRA across `[0, 3, 6, 12, 24, 48, 72]` hours, loaded with `regulation.load_data()`).
- **Split discipline:** Use the exact deterministic family-grouped `train` / `validation` / `test` assignment from `regulation.family_groups`.
- **Unseen test cohort:** All `test`-split loci **excluding** the 5 curated benchmark loci in `circuit.CURATED_LOCI` (`chr1:175099416-175099587` falls in `test` and is excluded). This yields **118 unseen test loci** and **1,885 non-WT measured sequence perturbations**.

## Pre-Specified Model & Prediction Rules (No Outcome Access)

`predict.py` reads only `r.locus`, `r.split`, `r.id`, `r.kind`, `r.group`, `r.reference`, and `r.sequence` for the 118 unseen `test` loci (never touching `r.alpha` or `r.wt_alphas`):

1. **Circuit & Physics:** Default `neural-tube` circuit (`Pax6`, `Olig2`, `Nkx2.2`, `Irx3`), `GLI2` motif (`epsilon = 2e-5`), non-specific partition floor `Z_ns = 2e-4`, signaling regime `S ∈ [0.06, 0.85]` in `9` steps, passive reporter silence gate `max_S |Δreporter(S)| < 0.02`.
2. **Sequence Predictor Gate:** `kreimer_mpra_mean` (`epigenesis-kreimer-motif-temporal-v1`, trained only on `train+validation` with no `test` leakage) with silence threshold `|kreimer_mpra_mean| < 0.15`; `siraj_score` is logged separately as informational only and excluded from gating (`r = 0.049`).
3. **Per-Locus Scan (`118` loci × `3*L = 513` SNVs):** Predict WT critical signal `S_crit`, category counts (`cryptic_bifurcation`, `concordant_switch`, `canalized_buffered`, `neutral`), `doubly_silent_count`, and top cryptic edit per unseen locus.
4. **Per-Perturbation Variant Calls (`1,885` alleles):** Predict `affinity_ratio`, `max_abs_reporter_delta`, `regime_flip_fraction`, `max_abs_target_delta`, `delta_s_crit`, `signal_locked`, `category`, `kreimer_mpra_mean`, `siraj_score`, `sequence_silent`, and `doubly_silent`.

## Pre-Specified Checks & Decision Rules

`score.py` first verifies the SHA-256 digests in `FREEZE.json` (`PLAN.md`, `predict.py`, `score.py`, `predictions.json`), then scores `predictions.json` **once** against `r.alpha`:

- **Check H1 (Locus cryptic rate & null comparison on top unseen cryptic locus):** Compare the prevalence of `3*L` SNV cryptic switches across the 118 unseen test loci against curated loci, and run `n_shuffles = 100` dinucleotide shuffles + column-scrambled PWM controls on the top unseen cryptic locus (if any).
- **Check H2 (Passive-vs-measured reporter concordance across 1,885 unseen alleles):** Report Pearson and Spearman correlations of mechanistic `max_abs_reporter_delta` and `kreimer_mpra_mean` vs `|measured_alpha_mean|`, plus exact confusion counts across categories (`cryptic_bifurcation`, `concordant_switch`, `canalized_buffered`, `neutral`).
- **Check H3 (Cryptic enrichment vs concordant silence):** Test whether predicted `cryptic_bifurcation` variants on the 1,885 measured `test` alleles exhibit lower `|measured_alpha_mean|` (or higher fraction with `|measured_alpha_mean| < 0.15`) than predicted `concordant_switch` variants. Report `0` or any negative/null outcome verbatim.

All results carry `biological_acceptance = false`.

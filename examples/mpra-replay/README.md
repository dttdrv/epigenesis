# Replayable MPRA Variant-Effect Reanalysis (`GSE188264`)

This example addresses the MPRA reproducibility gap highlighted in Rosen et al. (2026) (*Genome Research* 36:1940–1951): pipelines standardize processing, while Epigenesis seals **exact source bytes, GA4GH VRS 2.x variant identities, declared calling parameters, and computed variant calls** into independently replayable bundles.

## Pipeline

1. **Source admission (`fasta-dna/v1`)**: Verifies the pinned SHA-256 digests of Kreimer et al. (2022) `GSE188264` `design.fasta` and `alphas.csv.gz` and compiles locus WT sequences into `source/` (`brainc.source-descriptor/v1`).
2. **Reference-bound variant layer (`brainc.variant-collection/v1`)**: Extracts minimal edit spans for all tested alleles and binds every allele to a `ga4gh:SL.<sha512t24u>` / `ga4gh:VA.<sha512t24u>` identifier in `variants.json`, verified via `brainc.validator_variant.validate_variant_paths`.
3. **Dual parameter-set compilation (`brainc.development-bundle/v1`)**:
   - `permissive_mean_fc`: All-timepoint (0–72h) mean $|\Delta\text{asinh\_log2}| \ge 0.50$, WT floor $0.0$.
   - `strict_late_qc`: Late-window (24–72h) mean $|\Delta\text{asinh\_log2}| \ge 1.00$, sign-consistent across 24h/48h/72h, WT floor $\ge 1.00$.
4. **Independent replay & export**: Both parameter sets compile into sealed `brainc.development-bundle/v1` modules verified by `validate_development_paths`, report the exact VRS 2.x variant calls that flip between parameterizations, and export deterministic RO-Crate 1.1 and COMBINE/OMEX archives.

## Usage

```bash
python3 examples/mpra-replay/mpra_replay.py --split test -o output/mpra-replay-test
```

# Epigenesis HelixLoop Benchmark & Audit Suite

Reproducible training, evaluation, zero-shot in vivo transfer, and pretraining-leakage auditing for **Epigenesis HelixLoop v4** (`studio/epigenesis_studio/helixloop.py` and `studio/epigenesis_studio/helixloop_engine.py`).

## Evaluation Tiers & Licenses

| Tier | Model | Pretrained Backbone | License | Shipped Checkpoint |
|---|---|---|---|---|
| **Tier 1 (From Scratch)** | `HelixLoop v4` (`3.46M` params) | None (`0` external weights) | **Apache-2.0** | `studio/epigenesis_studio/model/helixloop_v4_multicell_fp16.pt` (`6.9 MB`, SHA-256 verified) |
| **Tier 2 (Open Foundation-Hybrid)** | `Flashzoi + HelixLoop v4` (`32.1M` params/rep) | `johahi/flashzoi-replicate-{0..3}` | **MIT** (Flashzoi) / **Apache-2.0** (Borzoi & HelixLoop) | Reproducible via `benchmark_flashzoi.py` |
| **Published Closed/NC Reference** | `AlphaGenome-Encoder-FT` | AlphaGenome (Avsec et al., *Nature* 2026) | **CC-BY-NC 4.0** (non-commercial) | Published benchmark numbers |

## Key Scripts

- `benchmark_lentimpra.py` — Tier-1 single-cell ENCODE4 lentiMPRA (`HepG2`), 31-mer zero-overlap CAGI5 scrubbing, monotonic `T=1..8` loop scaling, and causal ablations (`helixloop_v4` vs `helixloop_v4_no_delta`, `helixloop_omega`, `legnet_official`, `legnet_large`, `dream_rnn`).
- `benchmark_multicell.py` — Tier-1 multi-lineage (`HepG2 + K562 + WTC11`, `281,170` scrubbed CREs, 0% overlap with Fold 1 or CAGI5) with Cell-State FiLM conditioning.
- `benchmark_flashzoi.py` — Tier-2 open MIT/Apache-2.0 `Flashzoi + HelixLoop v4` vs controlled `Flashzoi-Linear-FT` ablation across `johahi/flashzoi-replicate-{0..3}`.
- `eval_tier2_ensemble.py` — Exact 4-replicate ensemble evaluation on ENCODE4 Fold 1 and zero-shot CAGI5 (`39,170` SNVs / `29` elements).
- `evaluate_invivo.py` & `eval_invivo_4rep_and_folds.py` — Zero-shot in vivo human genetics benchmark (`ADASTRA` ChIP-seq ASB + `UDACHA` DNase/ATAC-seq, `1,974` verified GRCh38 SNPs across 6 panels) with exact Borzoi `fold3+fold4` unseen-locus stratification.
- `benchmark_gse188264.py` — Human 0–72 h neural-induction transfer (`GSE188264`, `1,895` held-out perturbations across `118` loci) + `51` strictly test-only double-mutant epistasis triplets.
- `audit_statistics_and_leakage.py` — Recomputes exact two-sided finite-sample Phipson–Smyth (2010) paired bootstrap $p$-values ($B = 2{,}000$) and Borzoi `fold3+fold4` stratification into `results/canonical_audit_summary.json`.

## Verified Artifacts (`results/`)

All raw benchmark logs and audited summaries are stored in `results/`:
- `benchmark_HepG2_omega_v1_s42.json`, `benchmark_HepG2_v4_hepg2_s42.json`, `benchmark_HepG2_v3_arch_std_s42.json`
- `benchmark_multicell_tier1_v4_mc_s42.json`
- `benchmark_tier2_flashzoi_rep{0,1,2,3}_s{42,43,44,45}.json`, `benchmark_tier2_4rep_ensemble.json`
- `benchmark_invivo_zeroshot.json`, `benchmark_invivo_4rep_audited.json`
- `benchmark_gse188264_v4.json`
- `canonical_audit_summary.json`

## Pretraining-Split Audit Disclosure

Per `calico/borzoi/README.md` and `borzoi_sequences_human.bed.gz` (`55,497` hg38 intervals), all four `johahi/flashzoi-replicate-{0..3}` checkpoints share the **identical** human pretraining split (`train = {fold0, fold1, fold2, fold5, fold6, fold7}`, `val = fold4`, `test = fold3`; zero MPRA, CAGI5, ADASTRA, UDACHA, or GSE188264 labels seen during pretraining). Ensembling `rep0..rep3` therefore performs pure replicate variance reduction with zero cross-fold leakage, leaving 25% of hg38 (`fold3 + fold4`) strictly unseen.

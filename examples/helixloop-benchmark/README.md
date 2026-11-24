# HelixLoop v4 Benchmark & Audit Suite

Reproducible training, evaluation, zero-shot *in vivo* transfer, epistasis testing, and pretraining-leakage auditing for **Epigenesis HelixLoop v4** (`studio/epigenesis_studio/helixloop.py`).

## Install Benchmark Dependencies

```sh
pip install -e . -e "studio[benchmarks]"
```

## Two-Tier Separation & Licenses

| Tier | Model | Pretrained Backbone | License | Shipped Artifact |
|---|---|---|---|---|
| **Tier 1 (From Scratch)** | `HelixLoop v4` (`3.46M`) | None (`0` external weights) | **Apache-2.0** | `studio/epigenesis_studio/model/helixloop_v4_multicell_fp16.pt` (`6.9 MB`) |
| **Tier 2 (Open Hybrid)** | `Flashzoi + HelixLoop v4` (`32.1M`/rep) | `johahi/flashzoi-replicate-{0..3}` | **MIT** + **Apache-2.0** | Reproducible via `benchmark_flashzoi.py` |
| **Published NC Reference** | `AlphaGenome-Encoder-FT` | AlphaGenome (*Nature* 2026) | *CC-BY-NC 4.0* | Published metrics only |

## Scripts

- `models.py` — Shared model registry (`HelixLoopV4Model`, `HelixLoopOmegaModel`, `OfficialMPRALegNet`, `LegNet-Large`, `OfficialDREAMRNN`, `MalinoisCNN`, `DeepSTARR`, JASPAR 2024 PWM bank).
- `benchmark_lentimpra.py` — Tier-1 single-cell ENCODE4 lentiMPRA (`HepG2`), 31-mer zero-overlap CAGI5 scrubbing, monotonic `T=1..8` loop scaling, and `BiGatedDeltaProduct1d` causal ablation.
- `benchmark_multicell.py` — Tier-1 multi-lineage (`HepG2 + K562 + WTC11`, `281,170` scrubbed CREs, 0% Fold-1 / CAGI5 overlap) with Cell-State FiLM.
- `benchmark_flashzoi.py` — Tier-2 `Flashzoi + HelixLoop v4` vs controlled `Flashzoi-Linear-FT` ablation across `johahi/flashzoi-replicate-{0..3}`.
- `eval_tier2_ensemble.py` — Exact 4-replicate ensemble evaluation on ENCODE4 Fold 1 and zero-shot CAGI5 (`39,170` SNVs).
- `evaluate_invivo.py` & `eval_invivo_4rep_and_folds.py` — Zero-shot *in vivo* `ADASTRA` (ChIP-seq ASB) + `UDACHA` (DNase/ATAC-seq) across 6 panels (`1,974` GRCh38 SNPs) with Borzoi `fold3+fold4` stratification.
- `benchmark_gse188264.py` — `GSE188264` 0–72 h neural induction (`1,895` held-out edits) + `51` strictly test-only double-mutant epistasis triplets.
- `audit_statistics_and_leakage.py` — Recomputes exact two-sided finite-sample Phipson–Smyth (2010) paired bootstrap $p$-values ($B=2{,}000$) and Borzoi `fold3+fold4` audits into `results/canonical_audit_summary.json`.

## Pretraining-Split Audit

Per `calico/borzoi/README.md` and `borzoi_sequences_human.bed.gz` (`55,497` hg38 intervals), all four `johahi/flashzoi-replicate-{0..3}` models share the **identical** split (`train = {fold0,1,2,5,6,7}`, `val = fold4`, `test = fold3`; zero MPRA, CAGI5, ADASTRA, UDACHA, or GSE188264 labels seen in pretraining). Ensembling `rep0..rep3` reduces replicate variance with zero cross-fold leakage, leaving 25% of hg38 (`fold3+fold4`) strictly unseen.

All verified JSON outputs live in [`results/`](results/).

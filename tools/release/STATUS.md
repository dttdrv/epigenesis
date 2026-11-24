# Current delivery — Epigenesis 1.7.0

The local workbench is [Epigenesis Studio](../../studio/README.md) (`0.4.0`, English and Bulgarian). The zero-dependency compiler package (`epigenesis-compiler`) is version **1.7.0**.

## Delivered predictors

Studio ships four SHA-256-sealed checkpoints inside `studio/epigenesis_studio/model/`:

- **HelixLoop v4 (`helixloop`, optional PyTorch, Apache-2.0):** 6.9 MB fp16 multi-lineage looped biophysical model (`HepG2`, `K562`, `WTC11`; zero external weights). Beats `OfficialMPRALegNet` from scratch (`0.7946` vs `0.7721` 3-cell ENCODE4 Fold-1 mean), generalizes zero-shot to *in vivo* `ADASTRA`/`UDACHA` SNPs (5/6 panels, pooled `+0.026`, $p=0.016$), and pairs with MIT Flashzoi in Tier 2 to clear non-commercial `AlphaGenome-Encoder-FT` (`0.8708` vs `0.8670`) and reach `0.7513` on non-additive `GSE188264` epistasis.
- **Biophysical Asymmetric 532-TF+ v2 (`biophys`, NumPy CPU):** Dual-scale Boltzmann thermodynamics + multi-scale $k$-mers predicting 0–72 h `GSE188264` neural-induction kinetics (`r = 0.5060`, targeted MSE `0.250815`, Top-100 repression `99.0%`).
- **Kreimer 0–72 h (`--kreimer`, NumPy CPU):** Symmetric 5-mer + 532-TF ridge (`GSE188264` targeted MSE `0.252079`).
- **Siraj Huber-linear (`siraj`, NumPy CPU):** SK-N-SH 24 h reporter model (internal MSE `-3.02%` vs 5-mer; external 5-study mean `r = 0.0490`).

## DNA-site simulation & Enhancer Circuit

Simulate and Enhancer Circuit compile regulatory edits into validated `brainc` Development Modules and drive the Herrera-Delgado neural-tube ODE (plus optional 3D CX3D morphogenesis) with exact allele reversion and closed-chromatin controls (`biological_acceptance = false`).

## Evidence and reproduction

- [HelixLoop v4 paper (.tex)](../../docs/paper/epigenesis-helixloop.tex) and [benchmark & audit suite](../../examples/helixloop-benchmark/README.md).
- [1.7 release contract](../../docs/releases/RELEASE-1.7.md) and [Platform technical paper](../../docs/paper/epigenesis-technical-paper.pdf) ([Bulgarian](../../docs/paper/epigenesis-technical-paper.bg.pdf)).
- [Frozen research ledger](../../research/README.md): 13 hash-verified cross-study checks (`research/check_freeze.py`).

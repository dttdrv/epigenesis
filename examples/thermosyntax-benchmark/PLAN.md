# Pre-Registered Benchmark Plan: Epigenesis ThermoSyntaxNet vs. Market SOTA Architectures

**Date:** 2026-11-23  
**Status:** Pre-registered before running held-out test evaluations.

## 1. Scientific Question

Does **Epigenesis ThermoSyntaxNet** — a hybrid differentiable biophysical-neural architecture combining:
1. a 1-bp-resolution **532-TF JASPAR 2024 Differentiable Thermodynamic Energy Field** (forward + reverse-complement Boltzmann occupancy $p_m^{\mathrm{sym}}(x), p_m^{\mathrm{ori}}(x)$ and global log-partition sufficient statistics $\ln Z_m, \max_x E_m(x)$),
2. a **Multi-Scale De Novo DNA Stem** ($k \in \{5, 11, 17\}$) fused with the biophysical field,
3. **Anti-Aliased Dual-Pooling (Avg + Max)** + **Helical-Phase Dilated Inverted-Residual SE Blocks**,
4. **B-DNA Helical-Phase (10.5 bp) + Log-Distance Relative-Position Multi-Head Self-Attention** over full cis-regulatory elements, and
5. **Multi-Query Regulatory Attention Pooling + Global Thermodynamic Partition Residual + Heteroscedastic Uncertainty**

strictly outperform published state-of-the-art short-sequence cis-regulatory architectures (**Official MPRALegNet** [Agarwal et al., *Nature* 2025], **Parameter-Matched MPRALegNet-Large**, **MalinoisCNN** [Gosai et al., *Nature* 2024], and **DeepSTARR** [de Almeida et al., *Nature Genetics* 2022]) under **identical training data, identical splits, identical data augmentation, and zero leakage**?

## 2. Benchmarks & Zero-Leakage Splits

### Benchmark A: ENCODE4 Human lentiMPRA (Agarwal et al., *Nature* 2025)
- **Cell lines:** `WTC11` (human iPSC, $n=92,370$), `HepG2` (hepatocyte, $n=245,852$), `K562` (erythroid, $n=393,328$).
- **CAGI5 Zero-Overlap Scrubbing:** Any sequence whose 200-bp insert shares any 31-mer (forward or reverse complement) with any of the 29 CAGI5/Kircher reference elements is **removed from training and validation** before training begins.
- **Official 10-Fold Split:**
  - **Train:** `folds 3..10` (with `rev == 0`, using online random RC + ±21 bp adapter-flanked shift augmentation, matching the official `human_legnet` protocol).
  - **Validation (model selection only):** `fold 2` (using `rev == 0`, averaging F/R insert pairs).
  - **Test (scored ONCE after checkpoint selection):** `fold 1` (strictly held-out genomic regions; evaluated on BOTH raw per-orientation sequences AND official CRE orientation-averaged targets `*_averaged.tsv`).

### Benchmark B: True Zero-Shot CAGI5 / Kircher 2019 Saturation Mutagenesis (39,170 GRCh38 SNVs)
- **Data:** All 39,170 GRCh38 SNVs across 29 human disease regulatory elements, 21 disjoint genomic loci, 10 cell types (`Kircher et al., Nature Communications 2019`).
- **Protocol:** **ZERO Kircher training labels.** Every SNV is embedded in its exact 200-bp GRCh38 genomic window centered on the variant and flanked by the 15-bp lentiMPRA adapters (`230 bp` total). Predicted effect is the exact orientation-averaged difference:
  $$\Delta \hat{y}(S_{\mathrm{ref}}, S_{\mathrm{alt}}) = \tfrac{1}{2}\bigl[f(S_{\mathrm{alt}}) - f(S_{\mathrm{ref}}) + f(\mathrm{RC}(S_{\mathrm{alt}})) - f(\mathrm{RC}(S_{\mathrm{ref}}))\bigr].$$
- **Evaluation:** Per-element Pearson $r$ and Spearman $\rho$ across:
  1. **Lineage-matched elements** (`HepG2` model on 7 liver elements: `F9, HNF4A, LDLR, LDLR.2, SORT1, SORT1.2, SORT1-flip`; `K562` model on 6 erythroid elements: `BCL11A, GP1BA, HBB, HBG1, PKLR-24h, PKLR-48h`; `WTC11` model on neural/pluripotent/unseen elements), and
  2. **All 29 elements (21 disjoint loci)** using the multi-lineage model trained strictly on scrubbed ENCODE4 lentiMPRA (`HepG2 + K562 + WTC11`).

### Benchmark C: Human Neural Induction 0–72h Perturbations & Test-Only Double-Mutant Epistasis (`GSE188264`)
- **Data:** Kreimer et al. (2022) `GSE188264`, 9,241 perturbations across 7 timepoints `[0, 3, 6, 12, 24, 48, 72]h`.
- **Splits:**
  1. Official family-grouped split (`train` / `validation` / `test` = 118 held-out loci, 1,895 edits), explicitly noting retrospective exposure of the 118 loci.
  2. **Strictly Test-Only Double-Mutant Triplets ($n = 51$ triplets in held-out test loci ONLY — 0% leakage from the 203 train/val triplets!)**: evaluating double-mutant prediction ($r, \rho$), rescue transition $\Delta(A \to A+B)$, AND non-additive epistasis $\epsilon = \Delta_{AB} - \Delta_A - \Delta_B$.

## 3. Decision Gates (Pre-Registered)

1. **G1 (Beat Market SOTA on Held-Out lentiMPRA):** Paired 1,000-draw bootstrap 95% CI lower bound of `Pearson r(ThermoSyntaxNet) - Pearson r(OfficialMPRALegNet)` $> 0$ AND `Pearson r(ThermoSyntaxNet) - Pearson r(MalinoisCNN)` $> 0$ on held-out `fold 1` (both raw and CRE-averaged).
2. **G2 (Beat Parameter-Matched SOTA):** `ThermoSyntaxNet` beats parameter-matched `OfficialMPRALegNet-Large` (proving gain is architectural, not parameter count).
3. **G3 (Causal Ablation Proof):** Full `ThermoSyntaxNet` beats both `ThermoSyntaxNet-NoBiophys` (removing the 532-TF JASPAR Thermodynamic Energy Field) and `ThermoSyntaxNet-NoAttn` (removing B-DNA Helical-Phase Relative Attention).
4. **G4 (Zero-Shot CAGI5 SNV Superiority):** On strictly zero-shot CAGI5 (39,170 SNVs, 0 Kircher labels, 0 overlapping 31-mers), `ThermoSyntaxNet` achieves higher locus-mean Pearson $r$ than `OfficialMPRALegNet`, `MalinoisCNN`, and `biophys v2` (`0.1273`), with paired 21-locus bootstrap 95% CI lower bound $> 0$.
5. **G5 (Calibration & Exact Symmetry Invariants):**
   - Exact reversion identity: $|\Delta\hat{y}(S, S)| = 0$.
   - Exact antisymmetry: $|\Delta\hat{y}(S_{\mathrm{alt}}, S_{\mathrm{ref}}) + \Delta\hat{y}(S_{\mathrm{ref}}, S_{\mathrm{alt}})| < 10^{-6}$.
   - Exact RC invariance: $|\Delta\hat{y}(\mathrm{RC}(S_{\mathrm{ref}}), \mathrm{RC}(S_{\mathrm{alt}})) - \Delta\hat{y}(S_{\mathrm{ref}}, S_{\mathrm{alt}})| < 10^{-6}$.
   - Uncertainty calibration: positive Spearman correlation between predicted $\hat{\sigma}$ and actual $|y - \hat{\mu}|$ on held-out test data.

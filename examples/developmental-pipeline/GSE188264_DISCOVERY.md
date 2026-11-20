# Out-of-Sample 532-TF Biophysical Discovery on Human Neural Induction (`GSE188264`)

Reproducible with:

```sh
python3 examples/developmental-pipeline/discover.py --output /tmp/epigenesis-gse188264-solve
# or via the unified CLI:
epigenesis discover --output /tmp/epigenesis-gse188264-solve
```

## Problem & Motivation

Earlier exploratory testing ([`examples/neural-v2f/heldout-freeze/RESULTS.md`](../neural-v2f/heldout-freeze/RESULTS.md)) showed that guessing a single transcription factor (`GLI2`) yields near-zero empirical concordance (`Pearson r = 0.0123`) across broad unseen human neural enhancer perturbations (`GSE188264`, Kreimer et al. 2022).

Using pure **strand-symmetric Boltzmann partition thermodynamics** $\Delta \ln(Z_m + Z_{\mathrm{ns}})$ ($Z_{\mathrm{ns}} = 2 \times 10^{-4}$, zero k-mers) across all **532 JASPAR 2024 human TFs**, `discover.py` fits exclusively on `train` (`5,559` edits), selects ridge regularization exclusively on `validation` (`1,787` edits; selected $\lambda = 10{,}000$, validation $r = 0.2686$), refits `train+validation`, and evaluates **once** on the strictly held-out `test` split (`1,895` unseen edits across `118` unseen human genomic loci).

## Verified Held-Out Results (`1,895` Unseen Human Neural Enhancer Edits)

| Metric (Unseen `test` Split, `118` Loci) | Pure 532-TF Boltzmann Model | Control / Prior Baseline |
|---|---:|---:|
| Mean-activity Pearson $r$ (all `1,895` test edits) | **`0.3439`** (`95%` family CI `[0.2016, 0.4949]`) | Label-permutation null: `0.0279` · Single-TF GLI2: `0.0123` |
| Targeted-edit Pearson $r$ (`1,415` edits) | **`0.3366`** | — |
| Top-10% strongest predicted edits direction accuracy | **`93.2%`** (`0.9316`) | Chance: `50.0%` |
| Family-weighted MSE reduction vs zero-effect | **`10.86%`** (`0.1937` vs `0.2173`) | Zero-effect: `0.0%` |
| Sparse top-25 TF biophysical model Pearson $r$ | **`0.2891`** | — |

## End-to-End Sealed Compilation of Top Corroborated Unseen Edit

- **Unseen test locus & perturbation:** `chr10:100206539-100206710` (`hit2diff_pert2_chr10:100206539-100206710_ELF1_known3`)
- **GA4GH VRS 2.x allele ID:** `ga4gh:VA.hgnZPbJRl0FdjaT06IpWZlz8apU8fsF9`
- **Predicted vs. measured mean $\Delta\operatorname{asinh}(\alpha)$:** `-0.6605` predicted vs. `-1.4182` measured (`0h -> 72h`: `-0.6244 -> -0.6968` predicted vs. `-1.3102 -> -1.1630` measured)
- **Independent verification:** Both Python (`brainc.validator_variant`, `brainc.validator_development`) and the compiled zero-dependency Rust validator (`epigenesis-validate-rs`) returned `variants_valid = true`, `development_valid = true`.
- **3D morphogenesis & spiking:** Compiled through `epigenesis develop` into CX3D + LIF (`cx3d_executed = true`, `cx3d_rescue_exact = true`) and exported as RO-Crate 1.1 and COMBINE/OMEX.

## Scientific Boundary

`biological_acceptance = false`. These results establish strictly out-of-sample biophysical reporter prediction (`r = 0.3439`, `93.2%` top-decile direction accuracy) and end-to-end cryptographically verified compilation within the `GSE188264` human neural-induction assay; cross-study transfer and *in vivo* embryonic phenotypes remain external validation problems.

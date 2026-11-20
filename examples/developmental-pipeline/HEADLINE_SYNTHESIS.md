# Inverse Developmental DNA Compiler (`epigenesis synthesize`) — Verified Headline Results

Reproducible in one command:

```sh
python3 epigenesis_cli.py synthesize benchmark --output /tmp/epigenesis-headline
# Or target any custom enhancer / circuit / breaking mutation:
python3 epigenesis_cli.py synthesize rescue --locus chr18:56475717-56475888 --breaking-edit 165:G>T
python3 epigenesis_cli.py synthesize reprogram --locus chr15:90936956-90937127 --target-fate Olig2
```

## 1. Proof-Carrying Epistatic Rescue of 3D Neural Morphogenesis (Primary Mutation Kept Intact)

Starting from real human neural enhancers in `GSE188264` where a single-base mutation looks nearly silent in a passive dish (`|Δreporter| ≈ 0.012`) yet collapses endogenous `Nkx2.2` (`-0.932`), flips cell identity (`Nkx2.2+ → Olig2+`), and perturbs 3D CX3D neuronal outgrowth (`+22.68` neurite-length drift), `epigenesis synthesize rescue` locks the broken base untouched and synthesizes a compensating 1-bp second-site DNA edit:

| Real Human Enhancer | Broken Edit (Intact) | Synthesized 1-bp Rescue | Double-Mutant $q$ | Max Gene $\Delta$ vs WT | 3D CX3D Neurite Recovery | Rust + VRS 2.x Replay |
|---|---|---|---|---:|---:|---|
| `chr18:56475717-56475888` | `165:G>T` (`ga4gh:VA.Tb8SIdM4...`) | **`163:T>C`** (`+2 bp`, `ga4gh:VA.DXH7NbxW...`) | **`1.000125`** | **`9.89 × 10⁻⁸`** | **`100.00%`** (`351.14` → `373.82` → **`351.14`**) | `variants_valid = true`, `development_valid = true` |
| `chr1:175099416-175099587` | `149:G>T` (`ga4gh:VA.Mjx89n4V...`) | **`147:A>T`** (`+2 bp`, `ga4gh:VA.9bUzqzHX...`) | **`0.412254`** | **`0.00111`** | **`91.13%`** (`351.14` → `373.83` → **`349.13`**) | `variants_valid = true`, `development_valid = true` |

Every synthesized design emits sealed GA4GH VRS 2.x identifiers, a `brainc.development-bundle/v1`, independent Python + zero-dependency Rust replay verification, LIF spiking transmission controls (`all_lif_transmission_controls_pass = true`), and portable RO-Crate 1.1 / COMBINE OMEX archives.

## 2. Strictly Out-of-Sample Wet-Lab Inverse Selection (`1,895` Unseen Human Neural Edits, `119` Unseen Loci)

Using pure 532-TF Boltzmann thermodynamics trained exclusively on `train+validation` (`GSE188264`, zero k-mers), selecting top-$K$ repressing vs. activating edits on `1,895` unseen human neural enhancer perturbations (`overall Pearson r = 0.3439`):

| Budget ($K$) | Top-$K$ Down: True Repression | Down Mean $\Delta\operatorname{asinh}$ | Top-$K$ Up: True Activation | Up Mean $\Delta\operatorname{asinh}$ | Up − Down Gap | Exact Permutation $p$ |
|---:|---:|---:|---:|---:|---:|---:|
| `10` | `80.0%` | `-0.498` | **`60.0%`** (vs `30.9%` chance) | `-0.261` | `+0.236` | `0.102` |
| `25` | **`88.0%`** | **`-0.548`** | `52.0%` | `-0.200` | **`+0.348`** | **`0.0040`** |
| `50` | **`90.0%`** | `-0.411` | `46.0%` | `-0.118` | **`+0.293`** | **`0.0020`** |
| `100` | **`95.0%`** | `-0.446` | `45.0%` | `-0.080` | **`+0.366`** | **`0.000999`** |

Random baseline across all `1,895` unseen edits: `69.1%` down / `30.9%` up / mean $\Delta = -0.222$.

## Honest Scientific Scope

`biological_acceptance = false`. The wet-lab Top-$K$ benchmark validates strictly out-of-sample human MPRA reporter reprogramming within `GSE188264`; the 3D CX3D morphogenesis and LIF spiking rescues validate mechanistic software compilation, parameter-uncertainty concordance, and cross-language cryptographic replay rather than *in vivo* embryonic efficacy.

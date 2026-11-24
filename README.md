<p align="center">
  <img src="desktop/macos/AppIcon.png" alt="Epigenesis App Icon" width="112" height="112" />
</p>

<h1 align="center">Epigenesis</h1>

<p align="center">
  <strong>Proof-Carrying DNA-to-Development Simulation, Inverse Regulatory Compiler &amp; Variant Screening Workbench</strong><br />
  Model how regulatory and coding DNA changes propagate through transcription-factor binding, chromatin state, gene-regulatory networks, and 3D neural morphogenesis — with deterministic cryptographic replay in Python and Rust.
</p>

<p align="center">
  <a href="https://github.com/dttdrv/epigenesis/releases"><strong>Download Desktop App</strong></a> &middot;
  <a href="docs/paper/epigenesis-technical-paper.pdf">Technical Paper (EN)</a> &middot;
  <a href="docs/paper/epigenesis-technical-paper.bg.pdf">Технически доклад (BG)</a> &middot;
  <a href="docs/paper/epigenesis-v2f-grand-challenge.pdf">V2F Grand Challenge Paper (EN)</a> &middot;
  <a href="docs/paper/epigenesis-v2f-grand-challenge.bg.pdf">V2F Доклад (BG)</a> &middot;
  <a href="docs/ARCHITECTURE.md">Architecture</a> &middot;
  <a href="SPEC.md">Compiler Spec</a>
</p>

---

## What Is Epigenesis?

Epigenesis is a **cross-platform scientific workbench and universal DNA translation compiler** for studying and designing genetic perturbations. Instead of burying biological assumptions inside one-off scripts, Epigenesis treats DNA sequences and GA4GH VRS 2.x genetic edits as typed, content-addressed inputs, compiles them under explicit caller-supplied interpretations into sealed Development Modules, and executes multi-scale biological simulations and sequence-to-function predictors whose results can be replayed independently in both Python and Rust.

### Core Capabilities

- **Multi-Locus Developmental Simulation (`Develop Pipeline`):** Couple cis-regulatory enhancers, chromatin accessibility gates, and coding-sequence translation into N-gene ordinary differential equation (ODE) networks with 90% parameter-uncertainty ensembles, exact reversion and closed-chromatin controls, and 3D CX3D + LIF neuronal morphogenesis.
- **Inverse Regulatory Compiler (`Inverse Compiler`):** Synthesize minimal 1-bp distal compensatory edits that rescue 3D morphogenesis while leaving an uneditable primary lesion intact, or compile single-base non-coding edits that switch cell fate without touching coding exons.
- **HelixLoop v4 Neural Regulatory Engine (`helixloop_engine` & `studio/epigenesis_studio/helixloop.py`):** 6.9 MB fp16 SHA-256-sealed Apache-2.0 looped biophysical DNA model (zero external weights) combining regularized RCPS, 6-channel DNA shape, 532 JASPAR 2024 TFs, and a shared-weight $T=6$ iterative core with exact parallel WY Householder state tracking (`BiGatedDeltaProduct1d`) and cell-state FiLM — guaranteeing bit-for-bit `0.0` reversion, antisymmetry, and reverse-complement invariants with calibrated uncertainty and per-TF attributions.
- **Asymmetric Biophysical 532-TF+ Variant Predictor (`Positional Scan` & `Batch Edits`):** Score saturation mutagenesis landscapes and arbitrary SNV/MNV/indel panels on CPU in milliseconds using dual-scale strand-symmetric Boltzmann thermodynamics across all 532 JASPAR 2024 human transcription factors, peak-site disruption, and reference-context gating — predicting 0–72 h human neural-induction kinetics with per-TF attributions.
- **Enhancer Circuit & Regulatory Site Analysis:** Scan arbitrary enhancers (`9–500` bp) for cryptic cell-fate bifurcations where passive reporters barely shift (`|Δ| < 0.02`) yet endogenous feedback flips lineage identity.
- **Cryptographic Provenance & Audit (`Verify Receipt`):** Every run exports an RFC 8785 SHA-256 sealed receipt (plus RO-Crate 1.1 and COMBINE OMEX archives) verified independently by producer-isolated Python and zero-dependency Rust validators (`epigenesis-validate-rs`).

---

## Desktop Application & Platform Availability

Epigenesis is designed for every platform. Today, the dedicated native desktop application is available for **macOS**, with native **Linux** and **Windows** applications coming soon™.

| Platform | Desktop Application | Status |
|---|---|---|
| **macOS 13+ (Apple Silicon & Intel)** | Native SwiftUI / AppKit application (**`Epigenesis-1.7.0-macOS.dmg`**) — direct local compute pipe (no browser wrapper, no localhost ports) | **Available now** |
| **Linux** | Native desktop installer | **Coming soon™** |
| **Windows** | Native desktop installer | **Coming soon™** |

### Installing the Desktop App (macOS)

Download **`Epigenesis-1.7.0-macOS.dmg`** from [GitHub Releases](https://github.com/dttdrv/epigenesis/releases), mount the disk image, and drag **`Epigenesis.app`** into **`Applications`**.

To build the `.dmg` locally from source (requires Python 3.11+, [uv](https://docs.astral.sh/uv/), and Xcode Command Line Tools):

```sh
uv tool install --editable studio --with-editable .
python3 desktop/build_desktop.py
open dist/Epigenesis-1.7.0-macOS.dmg
```

See [`desktop/README.md`](desktop/README.md) for workspace setup and shortcuts (`⌘O` Open, `⌘↩` Run, `⌘S` Save Receipt; full English & Bulgarian UI).

---

## Verified Scientific Benchmarks & Showcases

Every result below reproduces from bundled datasets (`ENCODE4` lentiMPRA, `CAGI5`, `ADASTRA`/`UDACHA`, `GSE188264`, Hallikas GLI profiles, Herrera-Delgado neural-tube ODE, CX3D 3D morphogenesis) with exact two-sided finite-sample Phipson–Smyth (2010) bootstrap tests ($B=2{,}000$) and Borzoi/Flashzoi `fold3+fold4` pretraining-leakage audits. Full details appear in [`docs/paper/epigenesis-helixloop.tex`](docs/paper/epigenesis-helixloop.tex) and [`examples/thermosyntax-benchmark/README.md`](examples/thermosyntax-benchmark/README.md).

### 1. Epigenesis HelixLoop v4 — Tier 1 (From Scratch, 6.9 MB fp16, Apache-2.0)

Trained strictly from scratch (`0` external weights) on `281,170` scrubbed ENCODE4 lentiMPRA CREs (`HepG2 + K562 + WTC11`, `0%` 31-mer overlap with Fold 1 or CAGI5). Bundled SHA-256-verified checkpoint: [`studio/epigenesis_studio/model/helixloop_v4_multicell_fp16.pt`](studio/epigenesis_studio/model/helixloop_v4_multicell_fp16.pt).

<p align="center">
  <img src="docs/figures/card1_tier1_from_scratch.png" alt="HelixLoop v4 Tier 1 From-Scratch Benchmark" width="960" />
</p>

| Model (From Scratch, Apache-2.0) | Params | HepG2 Fold-1 $r$ | K562 Fold-1 $r$ | WTC11 Fold-1 $r$ | 3-Lineage Mean | Zero-Shot CAGI5 (Match / Eryth) | Zero-Shot In Vivo Pooled (`n=1,974`) |
|---|---:|---:|---:|---:|---:|---:|---:|
| `OfficialMPRALegNet` (3-head) | `1.32M` | `0.7592` | `0.7620` | `0.7949` | `0.7721` | `0.2994` / `0.3907` | `r = 0.3700` / `AUC 0.7312` |
| **`HelixLoop v4` (FiLM, `6.9 MB` fp16)** | **`3.46M`** | **`0.7888`** (`p = 0.001`) | **`0.7900`** (`p = 0.001`) | **`0.8049`** (`p = 0.061` n.s.) | **`0.7946`** (`+0.0225`) | **`0.3079` / `0.4489`** (`+0.0582`) | **`r = 0.3792` / `AUC 0.7336`** (`+0.0260`, `p = 0.016`) |

- **Monotonic test-time loop scaling ($T=1 \to 6$):** HepG2 single-cell MSE drops `0.2808 → 0.2323 → 0.1795 → 0.1720` (`r: 0.8122 → 0.8191`, beating `OfficialDREAMRNN` `0.7361` by `+0.0830`, `p = 0.001`); multi-lineage mean improves monotonically `0.7859 → 0.7909 → 0.7938 → 0.7946`.
- **Causal ablation (`BiGatedDeltaProduct1d`):** Removing only the parallel WY Householder state tracker (`n_h = 2`) under identical seed/data drops zero-shot CAGI5 overall `-0.0175`, liver `-0.0487`, and neural `-0.0436`.
- **Exact invariants:** Bit-for-bit `0.0` reversion, antisymmetry, and reverse-complement invariance verified across CPU/GPU.

### 2. Tier 2 Open MIT/Apache-2.0 Foundation-Hybrid vs. Non-Commercial AlphaGenome

Coupling MIT-licensed `Flashzoi` (`johahi/flashzoi-replicate-0..3`; upstream Borzoi is Apache-2.0) with `HelixLoop v4` clears published non-commercial `AlphaGenome-Encoder-FT` (`CC-BY-NC 4.0`). **Zero cross-fold leakage audit:** all four Flashzoi replicates share the identical Borzoi hg38 pretraining split (`borzoi_sequences_human.bed.gz`: 75% train, 25% `fold3+fold4` strictly unseen; 0 MPRA/CAGI5/ADASTRA/UDACHA/GSE188264 labels).

<p align="center">
  <img src="docs/figures/card2_tier2_open_vs_alphagenome.png" alt="Tier 2 Open MIT/Apache-2.0 Flashzoi + HelixLoop v4 vs AlphaGenome" width="960" />
</p>

| Model | License | HepG2 $r$ | K562 $r$ | WTC11 $r$ | 3-Lineage Mean | CAGI5 Match / Liver / Eryth | In Vivo Unseen `fold3+4` (`n=453`) |
|---|---|---:|---:|---:|---:|---:|---:|
| `AlphaGenome-Encoder-FT` (published) | *CC-BY-NC 4.0* | **`0.888`** | `0.879` | `0.834` | `0.8670` | — | — |
| `Flashzoi-Linear-FT` (4-rep ablation) | **MIT** | `0.8789` | `0.8816` | **`0.8516`** | `0.8707` | **`0.3998`** / `0.5615` / **`0.5919`** | `r = 0.5521` / `AUC 0.8139` |
| **`Flashzoi + HelixLoop v4` (4-rep)** | **MIT + Apache-2.0** | `0.8796` | **`0.8822`** | `0.8505` | **`0.8708`** | `0.3966` / **`0.5624`** / `0.5877` | `r = 0.5466` / **`AUC 0.8250`** |

*Controlled ablation transparency:* On raw ENCODE4 Fold-1 CREs, 4-rep `Flashzoi + HelixLoop v4` (`0.8708`) and `Flashzoi-Linear-FT` (`0.8707`) are statistically tied (all per-lineage `p > 0.22`): MIT Flashzoi drives the primary `0.795 → 0.871` CRE jump, while `HelixLoop v4` preserves that frontier accuracy and adds calibrated uncertainty ($\sigma$), 532-TF Boltzmann attribution, `+0.0111` unseen-fold *in vivo* AUROC (`0.8250` vs `0.8139`), and `+0.1052` non-additive double-mutant epistasis (`0.7513` vs `0.6461`).

### 3. Zero-Shot In Vivo Genetics (`ADASTRA`/`UDACHA`), `GSE188264` & 51 Test-Only Epistasis Triplets

<p align="center">
  <img src="docs/figures/card3_zeroshot_invivo_epistasis.png" alt="Zero-Shot In Vivo Genetics, CAGI5 and Non-Additive Epistasis" width="960" />
</p>

Evaluated on `GSE188264` (`1,895` held-out 0–72 h human neural-induction perturbations across `118` loci; **note:** retrospective loci, 0 SNVs, multi-base motif ablations) plus **`51` strictly test-only double-mutant epistasis triplets** ($\varepsilon = \Delta_{AB} - \Delta_A - \Delta_B$, `0%` train/val locus overlap) and the `424` perturbations in strictly Borzoi-unseen `fold3+fold4`:

| Model | All `1,895` Edits $r$ | Targeted MSE (vs. 0) | Borzoi-Unseen `fold3+4` (`n=424`) | 51 Double-Mutants $r$ | **51 Non-Additive Epistasis $r$** |
|---|---:|---:|---:|---:|---:|
| Zero-effect / `k5` ridge | `0.0000` / — | `0.36087` / `0.28519` | — | — | `0.0000` |
| `Biophys Asymmetric 532-TF+ v2` (linear) | `0.5060` | `0.25082` (`-30.5%`) | `0.5276` | `0.7375` | `0.1509` |
| `OfficialMPRALegNet` (scratch / T1 xfer) | `0.3628` / `0.5536` | `0.30598` / `0.25422` | `0.3036` / `0.5196` | `0.2821` / `0.5162` | `0.1988` / `-0.2141` |
| **`HelixLoop v4` (scratch / T1 xfer)** | `0.3537` / **`0.6086`** | `0.29503` / **`0.24323`** (`-32.6%`) | `0.1858` / **`0.6646`** | `0.5625` / **`0.6619`** | `0.4283` / **`0.5645`** (`p = 0.017` vs LegNet) |
| `Flashzoi-Linear-FT` (T2 xfer) | **`0.7626`** | **`0.18013`** (`-50.1%`) | `0.7778` | **`0.8825`** | `0.6461` |
| **`Flashzoi + HelixLoop v4` (T2 xfer)** | `0.7553` | `0.18144` (`-49.7%`) | **`0.7985`** | `0.8194` | **`0.7513`** (`+0.1052`, `p = 0.053` vs Lin; `p = 0.001` vs v2) |

### 4. In-Silico Dynamical Rescue & Cross-Language Provenance Verification (`biological_acceptance = False`)
Inside the published Herrera-Delgado ventral neural-tube ODE coupled to 3D CX3D, single-base GLI2 motif perturbations can shift a passive promoter readout by $|\Delta| \approx 0.012$ while crossing an endogenous saddle-node bifurcation. Locking the primary substitution intact, Epigenesis compiles a 1-bp second-site substitution whose scalar GLI2 partition ratio $q \approx 1.0$ restores the deterministic ODE and CX3D trajectory and verifies the entire compilation byte-for-byte in Python and Rust (**note:** restoring scalar $q \to 1.0$ recovers deterministic ODE/CX3D output by construction, and neither human `GSE188264` testbed locus is a validated endogenous *Nkx2.2* enhancer *in vivo*):

| Human Enhancer (`GSE188264`) | Locked Primary Lesion | Synthesized 1-bp Rescue | Double-Mutant $q$ | Max Gene $\Delta$ vs WT | 3D CX3D Neurite Recovery | Python + Rust Replay |
|---|---|---|---:|---:|---:|---|
| `chr18:56475717-56475888` | `165:G>T` (`ga4gh:VA.Tb8SIdM4...`) | **`163:T>C`** (`+2 bp`) | **`1.000125`** | **`9.89 × 10⁻⁸`** | **`100.00%`** (`351.14 → 373.82 → 351.14`) | Verified |
| `chr1:175099416-175099587` | `149:G>T` (`ga4gh:VA.Mjx89n4V...`) | **`147:A>T`** (`+2 bp`) | **`0.412254`** | **`0.00111`** | **`91.13%`** (`351.14 → 373.83 → 349.13`) | Verified |

### 5. Single-Base Non-Coding Fate Reprogramming & Multi-Locus Chromatin Shielding
- **Non-coding fate switch (`chr15:90936956-90937127`):** Compiles **`61:C>G`** (`ga4gh:VA.RR2HJSiG...`, `q = 0.00030`), switching ventral `Nkx2.2 → Olig2` with **`100%` ensemble stability** (`n = 24`) and verified 3D CX3D arbor geometry (`351.14 → 377.41`, exact reversion & closed-chromatin shielding `Δ = 0.0`, OMEX archive `12,262` B).
- **Multi-locus chromatin gating (`tri-lineage-chromatin`):** Couples PAX6/BCL6 cis-enhancers through an explicit chromatin accessibility gate into an `Excitatory / Inhibitory / Progenitor` GRN (`Excitatory → Progenitor`, `66.7%` flip probability in open chromatin; `0.0` perturbation when chromatin closes), with mandatory `source_inputs` SHA-256 tamper rejection across Python and Rust.

*Validation scope:* Empirical prediction benchmarks evaluate human neural-induction reporter perturbations within `GSE188264`; mechanistic ODE and 3D CX3D morphogenesis runs evaluate published mathematical models under declared parameter assumptions. Cross-study transfer across independent reporter assays (thirteen frozen tests in [`research/`](research/README.md), mean external correlation `0.049`) and *in vivo* embryonic efficacy remain open experimental frontiers.

---

## Architecture

Epigenesis is organized into three strict layers ([`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md)):

1. **Provenance Compiler Core (`brainc/`, `SPEC.md`, `validators/rust/`):** Zero-dependency universal DNA translation compiler, GA4GH VRS 2.x reference-bound variant collections, RFC 8785 canonical seals, RO-Crate 1.1 / COMBINE OMEX export, and independent Python + Rust replay validators.
2. **Desktop Workbench & Unified Developmental Engine (`desktop/`, `studio/`, `examples/developmental-pipeline/`):** Native desktop interface driving multi-locus DNA + VRS → chromatin/coding → N-gene GRN uncertainty → inverse DNA compilation → biophysical 532-TF+ prediction → 3D CX3D + LIF morphogenesis.
3. **Frozen Research Record (`research/`, `examples/`):** Thirteen pre-registered, hash-verified cross-study predictive tests (`research/check_freeze.py`) alongside self-contained mechanistic building blocks (`neural-circuit`, `neural-patterning`, `coding-consequences`, `neural-construction`, `neural-evidence`, `neural-fate`).

In the works reviewed in [Prior art](docs/reference/PRIOR_ART.md), no other system combines these boundaries; that is a scoped comparison of the cited corpus, not a benchmark against other simulators.

## The Compiler Foundation

Epigenesis preserves exact DNA source identity, validates a caller-supplied interpretation and lowers it against a target contract:

```text
DNA source → typed source bundle → explicit interpretation → Development Module
```

The interpretation supplies the biological or computational meaning. The compiler checks types, shapes, coordinates, references and budgets, then seals the compilation with content hashes (SHA-256 over RFC 8785 canonical JSON). A separate validator (in Python and zero-dependency Rust) replays the chain without importing the producer.

Raw IUPAC, FASTA, streaming reference FASTA, GenBank and GFF3 are built in. Other grammars use a version-2 executable external frontend and a separate executable validator. Both are digest-pinned and explicitly selected by the caller. The [FASTQ example](examples/external-fastq/README.md) implements this extension without changing the compiler. [SPEC.md](SPEC.md) defines the contracts; [Standards](docs/reference/STANDARDS.md) maps formats and identities to their specifications.

<details>
<summary><strong>Scripting &amp; CI Command-Line Reference (Optional)</strong></summary>

```sh
# Reproduce showcases headlessly (Python 3.11+ stdlib; Java 21 for CX3D; NumPy for discover/predict)
python3.11 epigenesis_cli.py synthesize benchmark --output /tmp/epigenesis-headline
python3.11 epigenesis_cli.py synthesize reprogram --locus chr15:90936956-90937127 --target-fate Olig2
python3.11 epigenesis_cli.py discover --output /tmp/epigenesis-discover
python3.11 epigenesis_cli.py develop --preset tri-lineage-chromatin --skip-cx3d --output /tmp/epigenesis-tri

# Run full verification suite
python3.11 epigenesis_cli.py check --full
cargo test --offline --manifest-path validators/rust/Cargo.toml
python3.11 -m unittest discover -s tests -v
python3.11 research/check_freeze.py
```

</details>

## Complete source-to-module example

Install the compiler from the checkout into a virtual environment:

```sh
python3.11 -m venv .venv
. .venv/bin/activate
python -m pip install .

demo_dir=$(mktemp -d)
brainc compile-source \
  --profile fasta-reference-dna/v1 \
  --sequence examples/minimal-development/sequence.fasta \
  --wrapper identity --output "$demo_dir/source"

brainc compile-development "$demo_dir/source" \
  --manifest examples/minimal-development/interpretation/manifest.json \
  --request examples/minimal-development/interpretation/request.json \
  --response examples/minimal-development/interpretation/response.json \
  --policy examples/minimal-development/interpretation/lowering-policy.json \
  --target examples/minimal-development/interpretation/target-contract.json \
  --output "$demo_dir/development"

brainc-validate-development "$demo_dir/source" "$demo_dir/development" \
  --source-input sequence=examples/minimal-development/sequence.fasta \
  --output "$demo_dir/validation.json"
```

The validation report has `valid: true`.

<details>
<summary>Stable example identities</summary>

These identities are recorded in [expected.json](examples/minimal-development/expected.json) and checked by the minimal-example release gate.

| Artifact | SHA-256 |
|---|---|
| Source descriptor | `7199f7b99358c3c2fa8b8df8f44a5431dac722aefc1d5ff18de9fc52d4ac4ced` |
| Development Module | `0d835308d23fe83a5f4cd0e28c26c678c510049a381752d8671c39e3805d2cea` |
| Development bundle | `9cc0b7c95027d777212beb75cdea04666c3333368ca5769c3203c6547d74b313` |
| Validation report | `3b8dafa503bcb376b50c798c0e4ceb2cf254505adeff83253d854fca61dc88d3` |

</details>

## License and attribution

Code is [Apache-2.0](LICENSE). Bundled datasets and third-party components retain their documented licenses and attributions (JASPAR 2024 CC BY 4.0; CX3D 0.03; see each example's `NOTICE`). Full citations and availability statements appear in the [English](docs/paper/epigenesis-technical-paper.pdf) and [Bulgarian](docs/paper/epigenesis-technical-paper.bg.pdf) technical papers.

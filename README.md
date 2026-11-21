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

Every result below reproduces deterministically from bundled human datasets (`GSE188264`, Hallikas GLI profiles, Herrera-Delgado neural-tube ODE, CX3D 3D morphogenesis) and is documented in full in the [English](docs/paper/epigenesis-technical-paper.pdf) and [Bulgarian](docs/paper/epigenesis-technical-paper.bg.pdf) technical papers.

### 1. Out-of-Sample Human Neural-Induction Prediction (`GSE188264`, `1,895` Unseen Edits)
Trained strictly on `train + validation` families (hyperparameters selected exclusively on `validation`) and evaluated once on `118` strictly held-out human genomic loci (`1,895` total / `1,415` targeted edits across `[0, 3, 6, 12, 24, 48, 72]` h):

| Model | Held-Out Pearson $r$ (All `1,895`) | Targeted Pearson $r$ (`1,415`) | Targeted Family MSE (Lower Better) | Top-10% Direction Accuracy | Top-100 Wet-Lab Repression |
|---|---:|---:|---:|---:|---:|
| Zero-effect baseline | `0.0000` | `0.0000` | `0.360870` | `50.0%` | `69.1%` |
| Dinucleotide baseline | — | — | `0.319358` | — | — |
| Canonical 5-mer ridge (`k5`) | — | — | `0.285185` | — | — |
| Kreimer 0–72 h Motif-Temporal (v1) | — | — | `0.252079` | — | — |
| Pure 532-TF Boltzmann (zero $k$-mers) | `0.3439` | `0.3366` | — | `93.2%` | `95.0%` (`p = 0.000999`) |
| **Biophysical Asymmetric 532-TF+ (v2, default)** | **`0.5060`** | **`0.4987`** | **`0.250815`** | **`96.3%`** | **`99.0%`** |

### 2. Proof-Carrying Epistatic 3D Morphogenesis Rescue (Primary Mutation Kept Intact)
When a single-base enhancer mutation barely moves a passive reporter (`|Δreporter| ≈ 0.012`) yet collapses endogenous *Nkx2.2* (`-0.932`), flips ventral lineage identity (`Nkx2.2+ → Olig2+`), and distorts 3D CX3D neurite outgrowth (`+22.68` length drift), Epigenesis locks the primary lesion untouched and compiles a single 1-bp distal compensatory edit verified in both Python and Rust:

| Human Enhancer (`GSE188264`) | Locked Primary Lesion | Synthesized 1-bp Rescue | Double-Mutant $q$ | Max Gene $\Delta$ vs WT | 3D CX3D Neurite Recovery | Python + Rust Replay |
|---|---|---|---:|---:|---:|---|
| `chr18:56475717-56475888` | `165:G>T` (`ga4gh:VA.Tb8SIdM4...`) | **`163:T>C`** (`+2 bp`) | **`1.000125`** | **`9.89 × 10⁻⁸`** | **`100.00%`** (`351.14 → 373.82 → 351.14`) | Verified |
| `chr1:175099416-175099587` | `149:G>T` (`ga4gh:VA.Mjx89n4V...`) | **`147:A>T`** (`+2 bp`) | **`0.412254`** | **`0.00111`** | **`91.13%`** (`351.14 → 373.83 → 349.13`) | Verified |

### 3. Single-Base Non-Coding Fate Reprogramming & Multi-Locus Chromatin Shielding
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

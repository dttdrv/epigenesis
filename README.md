<p align="center">
  <img src="desktop/macos/AppIcon.png" alt="Epigenesis macOS App Icon" width="112" height="112" />
</p>

<h1 align="center">Epigenesis</h1>

<p align="center">
  <strong>Native macOS Desktop Workbench &amp; Proof-Carrying DNA-to-Development Compiler</strong><br />
  Simulate regulatory DNA, compile minimal epistatic rescues and cell-fate edits, and replay every result cryptographically in Python and Rust.
</p>

<p align="center">
  <a href="desktop/README.md"><strong>Download / Build macOS <code>.dmg</code></strong></a> &middot;
  <a href="docs/paper/epigenesis-technical-paper.pdf">Technical Paper (EN PDF)</a> &middot;
  <a href="docs/paper/epigenesis-technical-paper.bg.pdf">Технически доклад (BG PDF)</a> &middot;
  <a href="docs/ARCHITECTURE.md">Architecture</a> &middot;
  <a href="SPEC.md">Compiler Spec</a>
</p>

---

## Epigenesis for macOS (`.dmg`)

Epigenesis is built first and foremost as a **dedicated native macOS desktop program** (`Epigenesis.app`, distributed as **`Epigenesis-1.7.0-macOS.dmg`**) — a quiet, utilitarian laboratory workbench using native Swift, AppKit, and WKWebView with zero Electron overhead, automatic Light/Dark mode, and full English & Bulgarian UI.

| Platform | Status | How to run |
|---|---|---|
| **macOS 12+ (Apple Silicon & Intel)** | **Available now** | Native **`Epigenesis-1.7.0-macOS.dmg`** — mount & drag `Epigenesis.app` into `Applications` |
| **Linux** | **Coming soon™** | Native installer planned (local `epigenesis-studio` workbench available today) |
| **Windows** | **Coming soon™** | Native installer planned (local `epigenesis-studio` workbench available today) |

### Install & Open on macOS

Download `Epigenesis-1.7.0-macOS.dmg` from the [GitHub Releases](https://github.com/dttdrv/epigenesis/releases) page — or build it directly from this repository in seconds with Xcode Command Line Tools (`swiftc`, `iconutil`, `hdiutil`) and [uv](https://docs.astral.sh/uv/):

```sh
# 1. Install the local Studio runtime once (Python 3.11+)
uv tool install --editable studio --with-editable .

# 2. Build the native macOS .dmg installer (dist/Epigenesis-1.7.0-macOS.dmg)
python3 desktop/build_desktop.py

# 3. Open the disk image and drag Epigenesis.app into Applications
open dist/Epigenesis-1.7.0-macOS.dmg
```

Prefer installing straight to `~/Applications` without mounting the `.dmg`?

```sh
python3 desktop/build_desktop.py --install
open -a ~/Applications/Epigenesis.app
```

Everything runs **100% locally** on loopback (`127.0.0.1`), with standard macOS menus (`⌘O` Open, `⌘C`/`⌘V`, native `NSSavePanel` exports, Finder "Open With") and clean engine shutdown on `⌘Q`. See [`desktop/README.md`](desktop/README.md) and [`studio/README.md`](studio/README.md).

### Seven Integrated Laboratory Workbenches

| Tab | What you can do in the macOS app |
|---|---|
| **Develop Pipeline** | Compose multi-locus DNA + GA4GH VRS 2.x edits → chromatin accessibility & coding translation → N-gene GRN dynamics with 90% parameter-uncertainty intervals → optional 3D CX3D morphogenesis & LIF spiking, verified by Python + Rust. |
| **Inverse Compiler** | Synthesize minimal 1-bp epistatic rescues (keeping the primary lesion intact), compile single-base fate-reprogramming edits, or run out-of-sample 532-TF biophysical discovery on real human `GSE188264` enhancers. |
| **Simulate** | Perturb the 9-bp *Nkx2.2* GLI binding site with Hallikas GLI2/GLI3 profiles, inspect reporter-vs-endogenous divergence, and verify exact reversion controls. |
| **Enhancer Circuit** | Scan arbitrary enhancers (`9–500` bp) through Boltzmann partition physics across GLI2/GLI3 or 532 JASPAR 2024 motifs into `neural-tube`, `mutual-toggle`, or `tri-lineage-switch` ODEs to expose cryptic bifurcations. |
| **Scan & Edits** | Score all single-base substitutions or up to 5,000 arbitrary SNV/MNV/indel edit pairs with the Siraj Huber-linear or Kreimer 0–72 h temporal checkpoints and export family-diverse panels. |
| **Verify** | Drop any exported `.json` receipt to verify its RFC 8785 SHA-256 seal, model identity, and bit-exact deterministic replay. |

---

## Hard Problems Made Routine (Verified Showcases)

Traditional wet-lab screening and ad-hoc scripts struggle when regulatory edits interact nonlinearly across enhancers, chromatin state, gene-regulatory feedback, and 3D tissue growth — or when results cannot be independently reproduced. Every showcase below runs in seconds inside Epigenesis, emits GA4GH VRS 2.x + `brainc` canonical proofs, and is reproduced in full in the [English](docs/paper/epigenesis-technical-paper.pdf) and [Bulgarian](docs/paper/epigenesis-technical-paper.bg.pdf) technical papers.

### 1. Epistatic 3D Morphogenesis Rescue — With the Primary Mutation Left Intact
When a single-base enhancer lesion looks nearly silent in a passive reporter (`|Δreporter| ≈ 0.012`) yet collapses endogenous `Nkx2.2` (`-0.932`), flips ventral cell identity (`Nkx2.2+ → Olig2+`), and perturbs 3D neuronal outgrowth (`+22.68` neurite-length drift), reverting the original base is often blocked by PAM/structural constraints. Epigenesis locks the broken base untouched and synthesizes a single 1-bp distal compensatory edit that restores both gene-regulatory fate and 3D CX3D morphology:

| Human Enhancer (`GSE188264`) | Broken Edit (Kept Intact) | Synthesized 1-bp Rescue | Double-Mutant $q$ | Max Gene $\Delta$ vs WT | 3D CX3D Neurite Recovery | Python + Rust Replay |
|---|---|---|---:|---:|---:|---|
| `chr18:56475717-56475888` | `165:G>T` (`ga4gh:VA.Tb8SIdM4...`) | **`163:T>C`** (`+2 bp`, `ga4gh:VA.DXH7NbxW...`) | **`1.000125`** | **`9.89 × 10⁻⁸`** | **`100.00%`** (`351.14 → 373.82 → 351.14`) | `valid = true` |
| `chr1:175099416-175099587` | `149:G>T` (`ga4gh:VA.Mjx89n4V...`) | **`147:A>T`** (`+2 bp`, `ga4gh:VA.9bUzqzHX...`) | **`0.412254`** | **`0.00111`** | **`91.13%`** (`351.14 → 373.83 → 349.13`) | `valid = true` |

### 2. Zero-K-mer 532-TF Human Enhancer Discovery & Wet-Lab Inverse Design (`GSE188264`)
Local sequence regressions trained on hundreds of thousands of reporter measurements stall near `r = 0.049` across studies, and guessing a single transcription factor (`GLI2`) scores `r = 0.0123`. Using pure strand-symmetric **532-TF JASPAR 2024 Boltzmann thermodynamics** (`Δln(Z_m + Z_ns)`, zero k-mers) evaluated strictly out-of-sample on **`1,895` unseen human neural-induction edits** across `118–119` held-out loci:

- **Out-of-sample correlation:** **`r = 0.3439`** (`95%` family CI `[0.2016, 0.4949]`) vs label-permutation null (`0.0279`) and single-TF GLI2 (`0.0123`).
- **Top-decile direction accuracy:** **`93.2%`** (`0.9316`) vs `50.0%` chance.
- **Wet-lab inverse selection (Top-100 repressing designs):** **`95.0%` empirical wet-lab repression** (mean `Δasinh = -0.446`, Up − Down gap `+0.366`, exact permutation **`p = 0.000999`**; Top-25 achieves `88.0%`, `p = 0.0040`), with top corroborated unseen edit `chr10:100206539-100206710` (`ga4gh:VA.hgnZPbJR...`) compiled end-to-end through CX3D + LIF and verified in Python and Rust.

### 3. Non-Coding Single-Base Cell-Fate Reprogramming & Multi-Locus Chromatin Gating
- **Single-base fate switch without editing coding genes:** On human enhancer `chr15:90936956-90937127` (`513` SNVs scanned), Epigenesis compiles **`61:C>G`** (`ga4gh:VA.RR2HJSiG...`, `q = 0.00030`) switching ventral neural-tube `Nkx2.2` → `Olig2` with **`100%` target-fate ensemble concordance** (`24` parameter draws), shifting 3D CX3D arbor length (`351.14 → 377.41`), while exact reversion, closed-chromatin shielding, and passive-reporter isolation all hold to `0.0` error and pass Rust replay + RO-Crate 1.1 / COMBINE OMEX export (`12,262` B).
- **Multi-locus chromatin shielding (`tri-lineage-chromatin`):** Couples PAX6/BCL6 cis-enhancers through an explicit chromatin accessibility gate into a 3-lineage `Excitatory / Inhibitory / Progenitor` GRN (`Excitatory → Progenitor`, `66.7%` ensemble flip probability), where closing chromatin (`accessibility = 0.0`) completely silences the cis variant (`max_state_diff = 0.0`).

> **Honest Scientific Scope (`biological_acceptance = false`).** Software provenance, determinism, cross-language Python/Rust replay, and mechanistic causal controls pass; out-of-sample `GSE188264` wet-lab inverse selection validates strictly within the human neural-induction MPRA assay. *In vivo* embryonic phenotypes and cross-study sequence transfer remain unestablished (thirteen frozen cross-study tests in [`research/`](research/README.md) measured mean external correlation `0.049` and are preserved verbatim).

---

## Three-Layer Architecture

To keep biological ambition honest and core software rock-solid, Epigenesis separates three strict layers ([`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md)):

1. **Provenance Compiler Core (`brainc/`, `SPEC.md`, `validators/rust/`, `brainc/schemas/`):** Zero-dependency universal DNA translation compiler, GA4GH VRS 2.x reference-bound variant collections, RFC 8785 canonical seals, RO-Crate 1.1 / COMBINE OMEX export, and independent Python + Rust replay validators.
2. **Desktop Application & Unified Developmental Engine (`desktop/`, `studio/`, `examples/developmental-pipeline/`, `neural-circuit`, `neural-construction`):** Native macOS `.dmg` desktop app (`Studio 0.3.0`) driving multi-locus DNA + VRS → chromatin/coding → GRN uncertainty → inverse DNA compilation → 3D CX3D morphogenesis & LIF spiking.
3. **Exploratory Research & Honest Negative-Result Ledger (`research/`, `neural-v2f`, `mpra-replay`, `neural-predictor`, `neural-regulation`):** Clearly quarantined prototypes (`biological_acceptance = false`) and thirteen frozen cross-study predictive tests verified by `research/check_freeze.py`.

In the works reviewed in [Prior art](docs/reference/PRIOR_ART.md), no other system combines these boundaries; that is a scoped comparison of the cited corpus, not a benchmark against other simulators.

### Repository Map

| Path | Layer | Contents |
|---|---|---|
| [`desktop/`](desktop/README.md) | App | **Native macOS desktop app** (`EpigenesisApp.swift`, macOS HIG `AppIcon.png`, `build_desktop.py` → `Epigenesis-1.7.0-macOS.dmg`; Linux & Windows coming soon™) |
| [`studio/`](studio/README.md) | App | Epigenesis Studio `0.3.0` (English & Bulgarian) powering all seven desktop workbench modes and receipt replay |
| [`examples/developmental-pipeline/`](examples/developmental-pipeline/README.md) | Engine | Unified composable pipeline (`develop`), inverse DNA compiler (`synthesize`), and 532-TF biophysical discovery (`discover`) |
| `brainc/`, [`SPEC.md`](SPEC.md), `validators/rust/` | Core | Universal DNA compiler (`1.7.0`), GA4GH VRS 2.x variant layer, Draft 2020-12 schemas, RO-Crate/OMEX export, and Python + Rust replay validators |
| [`docs/`](docs/README.md), [`docs/paper/`](docs/paper/) | Docs | [Architecture](docs/ARCHITECTURE.md), standards trace, prior art, release contracts, and technical papers ([English](docs/paper/epigenesis-technical-paper.pdf), [Bulgarian](docs/paper/epigenesis-technical-paper.bg.pdf)) |
| [`research/`](research/README.md), [`examples/`](examples/) | Exploratory | Frozen negative-result ledger (13 tests + hash checker) and mechanistic/exploratory building blocks |
| [`epigenesis_cli.py`](epigenesis_cli.py), `tests/`, `studio/tests/` | CLI & QA | Optional standard-library scripting CLI, compiler/Rust/desktop/pipeline tests, Studio (`33`) tests, and `tests/verify.py` release gates |

---

## The Compiler Foundation

Epigenesis preserves exact DNA source identity, validates a caller-supplied interpretation and lowers it against a target contract:

```text
DNA source → typed source bundle → explicit interpretation → Development Module
```

The interpretation supplies the biological or computational meaning. The compiler checks types, shapes, coordinates, references and budgets, then seals the compilation with content hashes (SHA-256 over RFC 8785 canonical JSON). A separate validator (in Python and zero-dependency Rust) replays the chain without importing the producer.

Raw IUPAC, FASTA, streaming reference FASTA, GenBank and GFF3 are built in. Other grammars use a version-2 executable external frontend and a separate executable validator. Both are digest-pinned and explicitly selected by the caller. The [FASTQ example](examples/external-fastq/README.md) implements this extension without changing the compiler. [SPEC.md](SPEC.md) defines the contracts; [Standards](docs/reference/STANDARDS.md) maps formats and identities to their specifications; [Prior art](docs/reference/PRIOR_ART.md) places the design among existing systems.

---

<details>
<summary><strong>Advanced: Command-Line Workflows &amp; Automated Verification (Scripting / CI)</strong></summary>

For headless pipelines, batch sweeps, and CI automation (Python 3.11+ stdlib; add Java 21 for 3D CX3D morphogenesis and NumPy for `discover` / `predict`):

```sh
# Launch the native macOS desktop app from terminal
python3.11 epigenesis_cli.py app

# Reproduce the verified showcases from this README & the technical papers
python3.11 epigenesis_cli.py synthesize benchmark --output /tmp/epigenesis-headline
python3.11 epigenesis_cli.py synthesize rescue --locus chr18:56475717-56475888 --breaking-edit 165:G>T
python3.11 epigenesis_cli.py synthesize reprogram --locus chr15:90936956-90937127 --target-fate Olig2
python3.11 epigenesis_cli.py discover --output /tmp/epigenesis-discover
python3.11 epigenesis_cli.py develop --preset tri-lineage-chromatin --skip-cx3d --output /tmp/epigenesis-tri

# Full verification & release gates
python3.11 epigenesis_cli.py check --full
cargo test --offline --manifest-path validators/rust/Cargo.toml
python3.11 -m unittest discover -s tests -v
python3.11 tests/verify.py developmental-pipeline
PIP_NO_INDEX=1 python3.11 tests/verify.py wheel
PIP_NO_INDEX=1 python3.11 tests/verify.py sdist
python3.11 research/check_freeze.py
```

The [CI workflow](.github/workflows/ci.yml) runs all release checks, independent source replay, contract attacks, resource/path safety, real-DNA compilation, macOS `.app`/`.dmg` packaging tests, and the Studio suite (`33` tests in `studio/tests/`).

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

The validation report has `valid: true`. This small example compiles an explicit base-count interpretation. The wheel installs the compiler and validators; the checkout and source archive also contain the desktop app, Studio, examples, data, and tests.

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

## Executable workflows & research record

### Layer 2 — Developmental Engine & Mechanistic Building Blocks

| Module | What runs | Main source |
|---|---|---|
| [**Developmental pipeline & inverse compiler**](examples/developmental-pipeline/README.md) | Composes multi-locus DNA FASTA + GA4GH VRS 2.x edits → chromatin-gated TF binding & coding translation → `N`-gene GRN + uncertainty ensembles → inverse DNA synthesis (`rescue`/`reprogram`) & 532-TF discovery → 3D CX3D + LIF + RO-Crate/OMEX | `brainc`; Hallikas (2006); JASPAR 2024; Herrera-Delgado (2020); Zubler et al. (2011, 2013) |
| [Coding consequences](examples/coding-consequences/README.md) | Reconstructs *Sox2* edits from genomic sequence, applies the declared coding annotation and checks peptides against independent references | Gong et al. (2020) |
| [Neural patterning](examples/neural-patterning/README.md) | Connects an *Nkx2.2* enhancer substitution and measured GLI binding profiles to a four-gene model with signal history, rescue, and phenocopy controls | Peterson et al. (2012); Cohen et al. (2014) |
| [Neural circuit](examples/neural-circuit/README.md) | Scans arbitrary enhancers via Boltzmann TF partition physics (GLI2/GLI3, 532 JASPAR 2024 motifs, or custom PWMs) into $N$-gene ODEs (`neural-tube`, `mutual-toggle`, `tri-lineage-switch`), finds second-site rescues, and bridges to CX3D + LIF | Hallikas et al. (2006); Herrera-Delgado et al. (2020); JASPAR 2024 |
| [Neural construction](examples/neural-construction/README.md) | Grows 3D cells, neurites, and synaptic contacts from one founder and an artificial DNA recipe; tests LIF signal transmission and contact lesions | Zubler et al. (2011, 2013) |
| [Neural evidence](examples/neural-evidence/README.md) | Reconstructs 109 published CHE-1 depletion observations, traced to workbook cells and intended DNA edits | Traets et al. (2021) |
| [Neural fate](examples/neural-fate/README.md) | Artificial DNA parameter recipe driving published regulatory dynamics | Balaskas et al. (2012) |

### Layer 3 — Exploratory Research & Honest Negative-Result Ledger (`biological_acceptance = false`)

| Module | Scope & recorded result | Main source |
|---|---|---|
| [**Frozen transfer ledger (`research/`)**](research/README.md) | 13 pre-registered, hash-verified cross-study sequence tests (Siraj Huber-linear external mean `r = 0.049`, MSE `+9.52%` vs zero-effect; nearest-TSS `0.460` beats all models on credible sets). Verified via `research/check_freeze.py`. | Siraj et al. (2026); Lee (2025); McAfee (2023); Matoba (2026); Guo (2023); Weiss (2021) |
| [V2F radar](examples/neural-v2f/README.md) *(exploratory)* | Morphogen tipping-point sweep with 1,000-shuffle controls; frozen held-out test (`1,885` unseen `GSE188264` alleles) measured `r = 0.012`. | Altschul & Erickson (1985); Kreimer et al. (2022) |
| [MPRA replay](examples/mpra-replay/README.md) *(provenance demo)* | Seals bundled `GSE188264` `alphas.csv.gz` under two calling thresholds with VRS 2.x IDs and RO-Crate/OMEX export. | Kreimer et al. (2022); Rosen et al. (2026) |
| [Neural regulation](examples/neural-regulation/README.md), [predictor](examples/neural-predictor/README.md), [measurement](examples/neural-measurement/README.md), [forecast](examples/neural-forecast/README.md), [binding](examples/neural-binding/README.md) | Kreimer 0–72 h motif-temporal checkpoint (`test MSE 0.252079` vs `k5 0.285185` across 118 families), Siraj Huber-linear checkpoint, BIRD repeat forecast (`-11.7%` CRPS), and failed FamilyCode *Nkx2-2* binding check. | Kreimer et al. (2022); Siraj et al. (2026); Majoros et al. (2020); Liu et al. (2025) |

---

## License and attribution

Compiler, desktop app, and Studio code are [Apache-2.0](LICENSE). Bundled scientific data and third-party components keep their documented licenses and attributions (JASPAR 2024 CC BY 4.0; CX3D 0.03; see each example's `NOTICE`). Full citations and availability statements are in the [English](docs/paper/epigenesis-technical-paper.pdf) and [Bulgarian](docs/paper/epigenesis-technical-paper.bg.pdf) technical papers.

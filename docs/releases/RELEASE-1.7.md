# Epigenesis 1.7

Version **1.7.0** integrates two major experimental capabilities cleanly across
[`examples/`](../../examples/), [`studio/`](../../studio/README.md),
[`tools/release/epigenesis.py`](../../tools/release/epigenesis.py), and
[`tests/verify.py`](../../tests/verify.py):

1. **Kreimer H1 Neural-Induction 0–72h Motif-Temporal Predictor (`epigenesis-kreimer-motif-temporal-v1`)**
2. **Variant-to-Circuit Enhancer Scanner (`epigenesis-neural-circuit-v1`)**

Compiler version bumps cleanly to `1.7.0` while preserving all serialized
`brainc` v1/v2 contracts, producer identities (`1.0.0`, compiler `0.4.0` / `0.5.0`),
and zero runtime dependencies (`dependencies = []`) for the `epigenesis-compiler` wheel.

## What changed in 1.7.0

### 1. Kreimer 0–72h motif-temporal neural-induction predictor

- Built in [`examples/neural-regulation/motif_temporal.py`](../../examples/neural-regulation/motif_temporal.py)
  and packaged as a deterministic Studio checkpoint in
  [`studio/epigenesis_studio/kreimer.py`](../../studio/epigenesis_studio/kreimer.py)
  with [`kreimer_model.json`](../../studio/epigenesis_studio/model/kreimer_model.json)
  and [`kreimer_weights.npz`](../../studio/epigenesis_studio/model/kreimer_weights.npz)
  (`(1044, 7)` float64 weights combining 5-mer counts and 532 JASPAR 2024 motifs
  over `[0, 3, 6, 12, 24, 48, 72]` hour timepoints).
- **Retrospective evaluation on Kreimer et al. (2022) GSE188264** (`1,415` targeted variants across `118` held-out families; hyperparameters frozen on train->validation split):
  - Family-weighted targeted MSE: **`0.252079`** (`0.25207852`) vs `k5` (`0.285185`),
    dinucleotide (`0.319358`), and zero-effect (`0.360870`).
  - 95% family-bootstrap MSE improvement intervals: `[0.041258, 0.188584]` vs zero-effect,
    `[0.008174, 0.138345]` vs dinucleotide, and `[0.007810, 0.061363]` vs `k5`.
  - Temporal contrast MSE (beats time-averaged baseline): **`0.031285`** vs
    `0.032108` (`[0.000220, 0.001528]`).
  - **Scientific scope:** Retrospective comparison conditional on published MPRAnalyze
    preprocessing on the previously examined `118` test families; short-enhancer lentiviral
    reporter induction only (`biological_development_accepted: false`).
- **Studio & CLI integration:** Selectable via the model toggle (`Siraj` vs
  `Kreimer 0–72h`) in Studio **Scan** and **Edits** panels (with 7-timepoint trajectory strip,
  top contributing JASPAR 2024 motifs per edit via `explain_rows`, CSV columns, and full SHA-256 canonical receipt replay) and via
  `epigenesis predict --kreimer` or `python3 epigenesis.py predict --kreimer`.

### 2. Variant-to-Circuit Enhancer Scanner

- Built in [`examples/neural-circuit/README.md`](../../examples/neural-circuit/README.md)
  ([`circuit.py`](../../examples/neural-circuit/circuit.py)) and exposed live in
  Epigenesis Studio `0.2.0` under the new **Enhancer Circuit** tab
  ([`studio/epigenesis_studio/circuit.py`](../../studio/epigenesis_studio/circuit.py)).
- Connects arbitrary-length enhancers (`9` to `500` bp in Studio) via strand-symmetric
  Boltzmann partition sums (`q = Z_variant / Z_reference`) over Hallikas et al. (2006)
  `GLI2`/`GLI3` Table S1 profiles directly to the Herrera-Delgado et al. (2020)
  four-gene (`Pax6`, `Olig2`, `Nkx2.2`, `Irx3`) neural-tube regulatory ODE, compiling
  SNV, MNV, insertion, and deletion pairs plus exact reversions into sealed variable-length
  `brainc` v2 Development Modules (`--target Nkx2.2|Olig2|Both`).
- **chr15:90936956-90937127 scan (`171` bp, `513` SNVs at `S = 0.20`, `GLI2`):**
  - Exposes **`7` cryptic bifurcations** (`|Δreporter| < 0.02` yet coupled feedback flips fate),
    `20` concordant switches, `0` canalized buffered, and `486` neutral edits.
  - Featured edit **`60:C>A`**: `q = 0.042455`, passive reporter shift `−0.0124`,
    endogenous `Nkx2.2` change `−0.9318` (`Nkx2.2 -> Olig2`, `75.2×` amplification).
  - Spatial gradient P3 domain contracts (`max x = 0.25 -> 0.125`), and pulse-chase
    primed `[(5.0, 1.0), (15.0, 0.05)]` confirms competitive GliR derepression
    (`0.8985` vs `0.000136` activation-only ablation).
  - Exact 9-bp zero-epsilon parity with [`examples/neural-patterning`](../neural-patterning/README.md)
    (`0.04235549749669733`, `< 1e-15`), exact reversion (`max_state_diff = 0.0`), and RK4
    refinement error `9.10e-10 < 1e-4` all pass.
- **Multi-locus benchmark, second-site rescue & CX3D bridge:**
  - `multi_locus_benchmark` evaluates 5 curated GSE188264 human enhancers alongside measured
    variant alleles, identifying ODE-verified epistatic **second-site fate rescues** on
    suboptimal enhancers (`chr18:56475717-56475888` `165:G>T` rescued by `163:T>C`,
    `chr1:175099416-175099587` `149:G>T` rescued by `147:A>T`).
  - `bridge_to_construction` (`--with-construction` / `bridge`) deterministically lowers
    four-gene steady states into 480-bp [`examples/neural-construction`](../neural-construction/README.md)
    recipes, running CX3D morphogenesis and LIF transmission assays for reference, variant,
    and exact allele-rescue controls.
  - **Scientific scope:** Conditional mechanistic transfer from in vitro GLI ratios to ODE
    promoter affinity; uncalibrated model units; whole-organism/brain phenotypes remain unestablished.

### 3. Unified CLI & Release Ergonomics

- Added zero-dependency [`epigenesis_cli.py`](../../epigenesis_cli.py) (installed as
  `epigenesis = "epigenesis_cli:main"`) exposing `studio`, `circuit`, `scan`, `predict`,
  `evidence`, `patterning`, `construction`, and `check` under a single command with automatic
  repository discovery and collision-free `./epigenesis-out/` defaults.
- Expanded [`tools/release/epigenesis.py`](../../tools/release/epigenesis.py) with `scan`,
  `circuit` (`--locus`, `--target`, `--protein`, `--signal`, `--reference`/`--variant`,
  `--with-construction`, `--benchmark`), and `construction`.

## Compatibility

- `epigenesis-compiler` version is `1.7.0` (`dependencies = []`, Python `>=3.11`).
- `epigenesis-studio` bumps to `0.2.0` (`dependencies = ["numpy>=1.26", "epigenesis-compiler>=1.7"]`)
  with bilingual EN/BG UI coverage for **Enhancer Circuit** (locus presets, target selector,
  custom variant comparison, second-site rescue card) and **Kreimer 0–72h** (with top JASPAR motifs).
- Frozen research records ([`research/`](../../research/README.md)) and Siraj Huber-linear
  frozen weights remain untouched and verified via `research/check_freeze.py`.

## Verification gates

New and updated release checks in [`tests/verify.py`](../../tests/verify.py) and `tests/`:

- `python3.11 tests/verify.py neural-circuit`
- `neural-regulation-motif` (runs with NumPy, e.g. Studio Python environment)
- `tests/test_epigenesis_cli.py` (unified CLI & release launcher checks)
- `PIP_NO_INDEX=1 python3.11 tests/verify.py wheel` and `sdist` checking `1.7.0` metadata, `epigenesis` CLI script, and bundled files
- Studio integration suite (`"$(uv tool dir)/epigenesis-studio/bin/python" -m unittest discover -s studio/tests -v`)

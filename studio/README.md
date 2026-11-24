# Epigenesis Studio

Local workbench and predictor runtime (`0.4.0`) powering both the browser UI (`epigenesis-studio`) and the native macOS desktop app ([`desktop/README.md`](../desktop/README.md)). Supports English and Bulgarian.

## Install

Requires Python 3.11+ (with `pip` or [uv](https://docs.astral.sh/uv/)):

```sh
# Core Studio (NumPy CPU predictors + developmental simulator + receipts)
pip install -e . -e studio

# Optional: enable the bundled 6.9 MB HelixLoop v4 neural engine (PyTorch)
pip install -e . -e "studio[helixloop]"

# Or install globally as a tool via uv
uv tool install --editable studio --with-editable .
```

## Predictor Engines

Four SHA-256-sealed checkpoints ship inside `epigenesis_studio/model/`:

- **`helixloop` (`epigenesis-helixloop-v4-multicell`, optional PyTorch):** Bundled 6.9 MB fp16 multi-lineage looped biophysical model (`HepG2`, `K562`, `WTC11`; zero external weights; Apache-2.0). Guarantees exact `0.0` reversion, antisymmetry, and reverse-complement invariance, with calibrated uncertainty ($\sigma$), 532-TF JASPAR Boltzmann shifts, and non-additive epistasis decomposition (`helixloop_engine.decompose_epistasis`). Degrades cleanly when PyTorch is not installed.
- **`biophys` (`epigenesis-biophysical-temporal-v2`, NumPy CPU):** Asymmetric 532-TF Boltzmann + multi-scale $k$-mer ridge predicting 0–72 h `GSE188264` neural-induction kinetics (`r = 0.5060`, Top-100 repression `99.0%`).
- **`kreimer` (`epigenesis-kreimer-motif-temporal-v1`, NumPy CPU):** Symmetric 5-mer + 532-TF 0–72 h temporal model (`GSE188264` targeted MSE `0.252079`).
- **`siraj` (`epigenesis-siraj-huber-linear-v1`, NumPy CPU, default legacy endpoint):** SK-N-SH 24 h reporter effect model (internal MSE `-3.02%` vs 5-mer; external 5-study mean `r = 0.0490`).

## Usage

```sh
epigenesis-studio                  # launch local workbench (127.0.0.1:8765)
epigenesis-studio data/            # point at a specific folder
epigenesis-studio edits.csv        # open an edit table directly
epigenesis-studio verify receipt.json
```

| Mode | What it does |
|---|---|
| **Develop Pipeline** | Multi-locus DNA + VRS 2.x → chromatin/coding → N-gene GRN with 90% CIs → optional 3D CX3D + Python/Rust verification |
| **Inverse Compiler** | Synthesize 1-bp compensatory rescues or fate-reprogramming edits and run 532-TF discovery scans |
| **Simulate** | 9-bp *Nkx2.2* GLI site substitution in the Herrera-Delgado four-gene neural-tube ODE |
| **Enhancer Circuit** | Scan 9–500 bp enhancers for cryptic cell-fate bifurcations (`neural-tube`, `mutual-toggle`, `tri-lineage-switch`) |
| **Scan / Edits** | Positional saturation scan (up to 1,000 bp) or batch variant scoring (up to 5,000 edits) across `helixloop`, `biophys`, `kreimer`, or `siraj` |
| **Verify** | Offline RFC 8785 SHA-256 digest, model identity, and deterministic replay check (`1e-9` tolerance) |

## Tests

```sh
python3.11 -m unittest discover -s studio/tests -v
```

## License

Apache-2.0. Bundled JASPAR 2024 motif profiles are CC BY 4.0; manifests live in `epigenesis_studio/model/`.

# Epigenesis Desktop Application

Epigenesis is built for every platform. Today, the dedicated **100% native desktop application** is available for **macOS (`Epigenesis-1.7.0-macOS.dmg`)**, with native **Linux** and **Windows** applications coming soon™.

| Platform | Status | Distribution |
|---|---|---|
| **macOS 13+ (Apple Silicon & Intel)** | **Available now** | Native SwiftUI / AppKit **`Epigenesis-1.7.0-macOS.dmg`** (drag `Epigenesis.app` into `Applications`) |
| **Linux** | **Coming soon™** | Native desktop installer planned |
| **Windows** | **Coming soon™** | Native desktop installer planned |

## Architecture (Zero Electron, Zero Localhost)

Unlike web-wrapped desktop apps, Epigenesis for macOS is compiled directly to a **native SwiftUI + AppKit + CoreGraphics binary** with:
- **Zero `WKWebView`, zero browser wrapper, zero `localhost` HTTP server, and zero open TCP ports.**
- Direct communication with the local Python/Rust scientific engine strictly over private OS `stdin`/`stdout` pipes (`desktop/epigenesis_engine.py`).
- Standard macOS `NavigationSplitView`, unified `NSToolbar`, SF Symbols sidebar, native `NSOpenPanel` / `NSSavePanel`, full English & Bulgarian UI, automatic Light/Dark mode, and a macOS HIG-compliant minimalistic icon.

## Scientific Workbenches Included

1. **Develop Pipeline** — Multi-locus DNA + GA4GH VRS 2.x → chromatin accessibility & coding translation → N-gene GRN steady-state & 90% uncertainty plot → optional 3D CX3D + LIF assay with Python + Rust verification.
2. **Inverse Compiler** — Epistatic 3D morphogenesis rescues (keeping the primary lesion intact), non-coding single-base fate reprogramming, and out-of-sample 532-TF human enhancer discovery.
3. **Regulatory Site** — 9-bp *Nkx2.2* GLI site perturbation with SnapGene-style nucleotide track and interactive hover-crosshair 4-gene ODE trajectory plot.
4. **Enhancer Circuit** — Arbitrary enhancer Boltzmann partition scanning across `neural-tube`, `mutual-toggle`, and `tri-lineage-switch` for cryptic cell-fate bifurcations.
5. **Positional Scan** — Retina saturation mutagenesis heatmap with clickable base cells, 0–72 h neural-induction kinetics plot, and top JASPAR 2024 TF attributions — powered by **Biophysical Asymmetric 532-TF+ v2** (`r = 0.5060`, `99.0%` Top-100 repression), Kreimer 0–72 h, and Siraj.
6. **Batch Edits** — Multi-variant SNV/MNV/indel scoring and family-diverse panel selection.
7. **Verify Receipt** — Offline RFC 8785 SHA-256 digest and deterministic replay audit.

## Build & Install on macOS

```sh
# 1. Install the local Studio runtime once (Python 3.11+ and uv)
uv tool install --editable studio --with-editable .

# 2. Build Epigenesis.app and Epigenesis-1.7.0-macOS.dmg
python3 desktop/build_desktop.py
open dist/Epigenesis-1.7.0-macOS.dmg
```

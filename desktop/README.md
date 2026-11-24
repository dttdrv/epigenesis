# Epigenesis Desktop App

Native macOS desktop workbench (`Epigenesis-1.7.0-macOS.dmg`) built directly in SwiftUI, AppKit, and CoreGraphics — no Electron, no `WKWebView`, and no localhost HTTP ports. The UI talks to the local Python/Rust engine over private OS `stdin`/`stdout` pipes (`desktop/epigenesis_engine.py`).

Linux and Windows users can run the local Studio workbench (`epigenesis-studio`) and CLI (`epigenesis`) today; native desktop installers for Linux and Windows are planned.

## Included Tabs

1. **Develop Pipeline** — Multi-locus DNA + GA4GH VRS 2.x → chromatin/coding → N-gene GRN uncertainty plots → optional 3D CX3D + LIF assay with Python/Rust verification.
2. **Inverse Compiler** — 1-bp distal compensatory rescues, non-coding fate reprogramming, and 532-TF enhancer discovery.
3. **Regulatory Site** — 9-bp *Nkx2.2* GLI substitution with nucleotide track and interactive 4-gene ODE trajectories.
4. **Enhancer Circuit** — Boltzmann partition scanning for cryptic cell-fate bifurcations across `neural-tube`, `mutual-toggle`, and `tri-lineage-switch`.
5. **Positional Scan & Batch Edits** — Saturation mutagenesis heatmaps, 0–72 h kinetics, TF attributions, and family-diverse panel selection.
6. **Verify Receipt** — Offline SHA-256 and deterministic replay verification.

## Build on macOS

Requires macOS 13+, Python 3.11+, [uv](https://docs.astral.sh/uv/), and Xcode Command Line Tools:

```sh
uv tool install --editable studio --with-editable .
python3 desktop/build_desktop.py
open dist/Epigenesis-1.7.0-macOS.dmg
```

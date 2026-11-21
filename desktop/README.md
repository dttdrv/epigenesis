# Epigenesis for macOS (`.dmg`)

Epigenesis is designed to be used primarily as a **dedicated native macOS desktop program** (`Epigenesis.app`, packaged as **`Epigenesis-1.7.0-macOS.dmg`**) — built with native Swift, AppKit, and WKWebView with zero Electron bloat and a calm, utilitarian laboratory/office interface.

| Platform | Status | Distribution |
|---|---|---|
| **macOS 12+ (Apple Silicon & Intel)** | **Available now** | Native **`Epigenesis-1.7.0-macOS.dmg`** (drag `Epigenesis.app` into `Applications`) |
| **Linux** | **Coming soon™** | Native desktop installer planned (use local `epigenesis-studio` in the meantime) |
| **Windows** | **Coming soon™** | Native desktop installer planned (use local `epigenesis-studio` in the meantime) |

## Install & Open on macOS

### Option 1 — Build & Open the `.dmg` Installer (Recommended)

Requires macOS 12+ with Xcode Command Line Tools (`swiftc`, `iconutil`, `hdiutil`) and the local Studio Python environment (`uv tool install --editable studio --with-editable .`):

```sh
python3 desktop/build_desktop.py
open dist/Epigenesis-1.7.0-macOS.dmg
```

Drag **`Epigenesis.app`** onto **`Applications`** inside the mounted disk image, then open **Epigenesis** from Spotlight or Finder.

### Option 2 — Build & Install Directly to `~/Applications`

```sh
python3 desktop/build_desktop.py --install
open -a ~/Applications/Epigenesis.app
```

## What's Inside the Native macOS App

- **Native macOS chrome & shortcuts:** Standard macOS menu bar (`Epigenesis`, `File → Open ⌘O`, `Edit` Cut/Copy/Paste, `View` Zoom/Full Screen, `Window`), Retina `.icns` app icon, Finder "Open With" for `.fasta`/`.fa`/`.csv`/`.json` files, native `NSSavePanel` for saving CSVs and receipts, and automatic clean shutdown of the local `127.0.0.1` engine on Quit (`⌘Q`).
- **All Epigenesis capabilities in one window (English & Bulgarian, automatic Light/Dark mode):**
  1. **Develop Pipeline** — Multi-locus DNA + GA4GH VRS 2.x → chromatin/coding → GRN uncertainty ensembles → 3D CX3D + LIF with sealed Python + Rust replay verification.
  2. **Inverse Compiler** — Synthesize minimal epistatic rescue edits, fate-reprogramming edits, or run out-of-sample 532-TF biophysical discovery on `GSE188264`.
  3. **Simulate** — Interactive 9-bp GLI regulatory site patterning & reversion controls.
  4. **Enhancer Circuit** — Arbitrary enhancer Boltzmann TF scanning (GLI2/GLI3/532 JASPAR 2024 motifs), cryptic cell-fate bifurcations, and second-site rescues.
  5. **Scan & Edits** — Fast Siraj Huber-linear and Kreimer 0–72h temporal variant scoring.
  6. **Verify Receipts** — Offline cryptographic replay & tamper verification for saved `.json` receipts.

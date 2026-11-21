# Epigenesis Dedicated Desktop Application

Epigenesis ships a dedicated, utilitarian laboratory desktop program built **macOS-first** (native Swift/AppKit/WKWebView `.app` bundle) with clean **Windows** and **Linux** desktop launchers — zero Electron bloat and zero third-party GUI dependencies.

## Quick Launch (All Platforms)

From the repository checkout:

```sh
python3 epigenesis_cli.py app
# Or open a specific FASTA, edit table, or receipt directly:
python3 epigenesis_cli.py app studio/epigenesis_studio/example-edits.json
```

## Build & Install Dedicated Packages

```sh
# macOS (builds native dist/Epigenesis.app and installs into ~/Applications/Epigenesis.app):
python3 desktop/build_desktop.py --platform macos --install

# Build for all three platforms (macOS .app + Linux .desktop + Windows .bat/.vbs/.ico):
python3 desktop/build_desktop.py --platform all
```

### Platform Details

- **macOS (Primary — `dist/Epigenesis.app`):**
  - Native Swift + AppKit + WKWebView binary (`~150 KB`, arm64/Apple Silicon & Intel).
  - Standard macOS menu bar (`Epigenesis`, `File → Open ⌘O`, `Edit` Cut/Copy/Paste, `View` Zoom/Full Screen, `Window`), native macOS `NSSavePanel` for exporting CSVs and cryptographic receipts, and automatic clean shutdown of the local `127.0.0.1` engine on Quit (`⌘Q`).
  - Double-click `dist/Epigenesis.app` in Finder or install into `~/Applications` with `--install`.
- **Windows (`dist/windows/`):**
  - Double-click **`Epigenesis.vbs`** (silent windowless launcher) or **`Epigenesis.bat`**.
  - Starts the local loopback engine and opens a dedicated frameless desktop window via Edge/Chrome `--app=` mode (falling back cleanly to the default browser).
- **Linux (`dist/linux/`):**
  - Run **`epigenesis-desktop`** or copy **`epigenesis.desktop`** into `~/.local/share/applications/` for your desktop application menu.
  - Opens a dedicated Chromium/Chrome/Edge/Brave `--app=` desktop window and terminates the loopback engine on close.

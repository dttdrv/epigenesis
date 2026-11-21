#!/usr/bin/env python3
"""Build the native Epigenesis macOS desktop application (.app + .dmg installer).

Zero third-party build dependencies:
  - Compiles `desktop/macos/EpigenesisApp.swift` into a native AppKit/WKWebView
    `dist/Epigenesis.app` bundle with a generated Retina `.icns` icon and `runtime.json`.
  - Packages a standard drag-to-Applications disk image at
    `dist/Epigenesis-1.7.0-macOS.dmg` via Apple `hdiutil` and verifies it.

Native Linux and Windows desktop installers are coming soon (TM).
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import platform
import shutil
import struct
import subprocess
import sys
import tempfile
import zlib

REPO = Path(__file__).resolve().parents[1]
APP_VERSION = "1.7.0"
DMG_NAME = f"Epigenesis-{APP_VERSION}-macOS.dmg"

sys.path.insert(0, str(REPO / "desktop"))
import epigenesis_desktop as desk  # noqa: E402


# ---------------------------------------------------------------------------
# Pure-Python minimal laboratory icon generator (PNG + macOS ICNS)
# ---------------------------------------------------------------------------

def _png_chunk(tag: bytes, data: bytes) -> bytes:
    crc = zlib.crc32(tag + data) & 0xFFFFFFFF
    return struct.pack(">I", len(data)) + tag + data + struct.pack(">I", crc)


def render_icon_png(size: int) -> bytes:
    """Render a crisp, understated laboratory app icon (dark slate rounded square + minimalist helix)."""
    w = h = size
    rows = bytearray()
    radius = size * 0.22
    cx = (w - 1) / 2.0
    cy = (h - 1) / 2.0
    bg_r, bg_g, bg_b = 24, 27, 32
    fg_r, fg_g, fg_b = 240, 241, 243
    acc_r, acc_g, acc_b = 112, 165, 248

    stroke = max(1.6, size * 0.058)
    rung_h = max(1.2, size * 0.036)
    amp = size * 0.16
    top_y = size * 0.21
    bot_y = size * 0.79
    span_y = bot_y - top_y
    rungs = (0.18, 0.36, 0.64, 0.82)

    for y in range(h):
        rows.append(0)  # PNG filter type 0
        for x in range(w):
            dx = max(0.0, abs(x - cx) - (w / 2.0 - radius))
            dy = max(0.0, abs(y - cy) - (h / 2.0 - radius))
            dist_edge = math.hypot(dx, dy) - radius
            if dist_edge > 0.7:
                rows.extend((0, 0, 0, 0))
                continue
            alpha = 255 if dist_edge <= -0.7 else int(round(255 * (0.7 - dist_edge) / 1.4))

            g_shift = int(round((1.0 - y / max(1, h - 1)) * 12))
            r, g, b = bg_r + g_shift, bg_g + g_shift, bg_b + g_shift

            if top_y - stroke <= y <= bot_y + stroke:
                t = max(0.0, min(1.0, (y - top_y) / span_y))
                angle = t * math.pi * 2.0
                x1 = cx + amp * math.cos(angle)
                x2 = cx - amp * math.cos(angle)
                d1 = abs(x - x1)
                d2 = abs(x - x2)

                lo_x, hi_x = (x1, x2) if x1 < x2 else (x2, x1)
                for rt in rungs:
                    ry = top_y + rt * span_y
                    if abs(y - ry) <= rung_h * 0.55 and lo_x <= x <= hi_x:
                        r, g, b = acc_r, acc_g, acc_b
                        break

                d_min = min(d1, d2)
                if d_min <= stroke * 0.55:
                    r, g, b = fg_r, fg_g, fg_b
                elif d_min <= stroke * 0.85:
                    mix = (stroke * 0.85 - d_min) / (stroke * 0.30)
                    r = int(round(r * (1 - mix) + fg_r * mix))
                    g = int(round(g * (1 - mix) + fg_g * mix))
                    b = int(round(b * (1 - mix) + fg_b * mix))

            rows.extend((r, g, b, alpha))

    ihdr = struct.pack(">IIBBBBB", w, h, 8, 6, 0, 0, 0)
    idat = zlib.compress(bytes(rows), level=6)
    return b"\x89PNG\r\n\x1a\n" + _png_chunk(b"IHDR", ihdr) + _png_chunk(b"IDAT", idat) + _png_chunk(b"IEND", b"")


def build_macos_app(dist_dir: Path, *, install_dir: Path | None = None) -> Path:
    """Compile native Swift/AppKit/WKWebView Epigenesis.app bundle for macOS."""
    swiftc = shutil.which("swiftc")
    if swiftc is None:
        raise RuntimeError("macOS native build requires Xcode Command Line Tools (`swiftc`)")

    app_dir = dist_dir / "Epigenesis.app"
    contents = app_dir / "Contents"
    macos_dir = contents / "MacOS"
    res_dir = contents / "Resources"
    if app_dir.exists():
        shutil.rmtree(app_dir)
    macos_dir.mkdir(parents=True)
    res_dir.mkdir(parents=True)

    swift_src = REPO / "desktop" / "macos" / "EpigenesisApp.swift"
    bin_path = macos_dir / "Epigenesis"
    subprocess.run(
        [
            swiftc,
            "-O",
            str(swift_src),
            "-o",
            str(bin_path),
            "-framework",
            "AppKit",
            "-framework",
            "WebKit",
        ],
        check=True,
    )

    iconutil = shutil.which("iconutil")
    if iconutil is not None:
        with tempfile.TemporaryDirectory() as tmp:
            iconset = Path(tmp) / "AppIcon.iconset"
            iconset.mkdir()
            for sz, name in (
                (16, "icon_16x16.png"),
                (32, "icon_16x16@2x.png"),
                (32, "icon_32x32.png"),
                (64, "icon_32x32@2x.png"),
                (128, "icon_128x128.png"),
                (256, "icon_128x128@2x.png"),
                (256, "icon_256x256.png"),
                (512, "icon_256x256@2x.png"),
            ):
                (iconset / name).write_bytes(render_icon_png(sz))
            subprocess.run(
                [iconutil, "-c", "icns", str(iconset), "-o", str(res_dir / "AppIcon.icns")],
                check=True,
            )

    py_path = desk.find_studio_python(REPO)
    (res_dir / "runtime.json").write_text(
        json.dumps({"python": str(py_path), "repo": str(REPO)}, indent=2) + "\n",
        encoding="utf-8",
    )

    info_plist = f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>CFBundleName</key><string>Epigenesis</string>
  <key>CFBundleDisplayName</key><string>Epigenesis</string>
  <key>CFBundleIdentifier</key><string>org.epigenesis.studio</string>
  <key>CFBundleVersion</key><string>{APP_VERSION}</string>
  <key>CFBundleShortVersionString</key><string>{APP_VERSION}</string>
  <key>CFBundleExecutable</key><string>Epigenesis</string>
  <key>CFBundlePackageType</key><string>APPL</string>
  <key>CFBundleIconFile</key><string>AppIcon</string>
  <key>NSHighResolutionCapable</key><true/>
  <key>LSMinimumSystemVersion</key><string>12.0</string>
  <key>CFBundleDocumentTypes</key>
  <array>
    <dict>
      <key>CFBundleTypeName</key><string>Epigenesis DNA or Receipt File</string>
      <key>CFBundleTypeRole</key><string>Viewer</string>
      <key>CFBundleTypeExtensions</key>
      <array>
        <string>fasta</string><string>fa</string><string>fna</string>
        <string>json</string><string>csv</string><string>tsv</string>
      </array>
    </dict>
  </array>
</dict>
</plist>
"""
    (contents / "Info.plist").write_text(info_plist, encoding="utf-8")

    codesign = shutil.which("codesign")
    if codesign is not None:
        subprocess.run([codesign, "--force", "--deep", "--sign", "-", str(app_dir)], check=False, capture_output=True)

    if install_dir is not None:
        install_dir.mkdir(parents=True, exist_ok=True)
        dest = install_dir / "Epigenesis.app"
        if dest.exists():
            shutil.rmtree(dest)
        shutil.copytree(app_dir, dest)
    return app_dir


def build_macos_dmg(app_dir: Path, dist_dir: Path, *, version: str = APP_VERSION) -> Path:
    """Package `Epigenesis.app` into a compressed drag-to-Applications `.dmg` installer."""
    hdiutil = shutil.which("hdiutil")
    if hdiutil is None:
        raise RuntimeError("macOS `.dmg` creation requires Apple `hdiutil`")
    if not (app_dir / "Contents" / "MacOS" / "Epigenesis").is_file():
        raise FileNotFoundError(f"built app bundle not found at {app_dir}")

    dist_dir.mkdir(parents=True, exist_ok=True)
    dmg_path = dist_dir / f"Epigenesis-{version}-macOS.dmg"
    if dmg_path.exists():
        dmg_path.unlink()

    with tempfile.TemporaryDirectory(prefix="epigenesis-dmg-") as tmp:
        stage = Path(tmp) / "Epigenesis"
        stage.mkdir()
        shutil.copytree(app_dir, stage / "Epigenesis.app", symlinks=True)
        os.symlink("/Applications", stage / "Applications")
        subprocess.run(
            [
                hdiutil,
                "create",
                "-volname",
                "Epigenesis",
                "-srcfolder",
                str(stage),
                "-ov",
                "-format",
                "UDZO",
                str(dmg_path),
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
        )

    subprocess.run(
        [hdiutil, "verify", str(dmg_path)],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
    )
    return dmg_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--platform",
        choices=("macos", "linux", "windows"),
        default="macos",
        help="target desktop platform (default: macos; linux and windows coming soon)",
    )
    parser.add_argument("--dist", type=Path, default=REPO / "dist", help="output directory (default: dist/)")
    parser.add_argument("--no-dmg", action="store_true", help="build Epigenesis.app only and skip .dmg packaging")
    parser.add_argument("--install", action="store_true", help="also copy Epigenesis.app into ~/Applications")
    args = parser.parse_args(argv)

    if args.platform in ("linux", "windows"):
        print(
            f"Epigenesis native desktop installer for {args.platform.title()} is coming soon (TM).\n"
            "For now, macOS ships a native .dmg installer, and Linux/Windows users can run "
            "`epigenesis-studio` locally.",
            file=sys.stderr,
        )
        return 0

    if platform.system() != "Darwin":
        raise SystemExit("macOS .app / .dmg builds must be run on macOS (Darwin).")

    install_dir = Path.home() / "Applications" if args.install else None
    app_dir = build_macos_app(args.dist, install_dir=install_dir)
    built: dict[str, str] = {"macos_app": str(app_dir)}
    if install_dir is not None:
        built["installed_app"] = str(install_dir / "Epigenesis.app")
    if not args.no_dmg:
        built["macos_dmg"] = str(build_macos_dmg(app_dir, args.dist))

    print(json.dumps({"built": built}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

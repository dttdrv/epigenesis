#!/usr/bin/env python3
"""Build dedicated Epigenesis desktop application packages (macOS-first, plus Linux & Windows).

Zero third-party build dependencies:
  - macOS: Compiles `desktop/macos/EpigenesisApp.swift` into a native AppKit/WKWebView
    `dist/Epigenesis.app` bundle with a generated Retina `.icns` icon and `runtime.json`.
  - Linux: Generates `dist/linux/epigenesis-desktop`, `epigenesis.desktop`, and `epigenesis.png`.
  - Windows: Generates `dist/windows/Epigenesis.bat`, `Epigenesis.vbs` (silent launcher),
    `epigenesis_launcher.py`, and `epigenesis.ico`.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import platform
import shutil
import stat
import struct
import subprocess
import sys
import tempfile
import zlib

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "desktop"))
import epigenesis_desktop as desk  # noqa: E402


# ---------------------------------------------------------------------------
# Pure-Python minimal laboratory icon generator (PNG + ICO + macOS ICNS)
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
            # Rounded-rect signed distance
            dx = max(0.0, abs(x - cx) - (w / 2.0 - radius))
            dy = max(0.0, abs(y - cy) - (h / 2.0 - radius))
            dist_edge = math.hypot(dx, dy) - radius
            if dist_edge > 0.7:
                rows.extend((0, 0, 0, 0))
                continue
            alpha = 255 if dist_edge <= -0.7 else int(round(255 * (0.7 - dist_edge) / 1.4))

            # Subtle vertical lab-slate gradient inside the tile
            g_shift = int(round((1.0 - y / max(1, h - 1)) * 12))
            r, g, b = bg_r + g_shift, bg_g + g_shift, bg_b + g_shift

            if top_y - stroke <= y <= bot_y + stroke:
                t = max(0.0, min(1.0, (y - top_y) / span_y))
                angle = t * math.pi * 2.0
                x1 = cx + amp * math.cos(angle)
                x2 = cx - amp * math.cos(angle)
                d1 = abs(x - x1)
                d2 = abs(x - x2)

                # Base-pair rungs between strands
                lo_x, hi_x = (x1, x2) if x1 < x2 else (x2, x1)
                for rt in rungs:
                    ry = top_y + rt * span_y
                    if abs(y - ry) <= rung_h * 0.55 and lo_x <= x <= hi_x:
                        r, g, b = acc_r, acc_g, acc_b
                        break

                # Twin helix strands (antialiased)
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


def write_ico_from_png(png_bytes: bytes, dest: Path) -> None:
    """Wrap a 256x256 PNG inside a valid Windows .ico container."""
    header = struct.pack("<HHH", 0, 1, 1)
    entry = struct.pack("<BBBBHHII", 0, 0, 0, 0, 1, 32, len(png_bytes), 6 + 16)
    dest.write_bytes(header + entry + png_bytes)


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

    # Compile Swift source into native binary
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

    # Generate macOS .icns icon set
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

    # Record discovered Python + repo path so Finder/Spotlight launches work cleanly
    py_path = desk.find_studio_python(REPO)
    (res_dir / "runtime.json").write_text(
        json.dumps({"python": str(py_path), "repo": str(REPO)}, indent=2) + "\n",
        encoding="utf-8",
    )

    # Write macOS Info.plist
    info_plist = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>CFBundleName</key><string>Epigenesis</string>
  <key>CFBundleDisplayName</key><string>Epigenesis</string>
  <key>CFBundleIdentifier</key><string>org.epigenesis.studio</string>
  <key>CFBundleVersion</key><string>1.7.0</string>
  <key>CFBundleShortVersionString</key><string>1.7.0</string>
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

    # Ad-hoc sign locally on macOS so Gatekeeper/LaunchServices treat it cleanly
    codesign = shutil.which("codesign")
    if codesign is not None:
        subprocess.run([codesign, "--force", "--deep", "--sign", "-", str(app_dir)], check=False, capture_output=True)

    if install_dir is not None:
        install_dir.mkdir(parents=True, exist_ok=True)
        dest = install_dir / "Epigenesis.app"
        if dest.exists():
            shutil.rmtree(dest)
        shutil.copytree(app_dir, dest)
        return dest
    return app_dir


def build_linux_bundle(dist_dir: Path) -> Path:
    """Create Linux standalone launcher, XDG .desktop entry, and PNG icon."""
    out = dist_dir / "linux"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    icon_path = out / "epigenesis.png"
    icon_path.write_bytes(render_icon_png(256))

    launcher = out / "epigenesis-desktop"
    launcher.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail
REPO="${{EPIGENESIS_HOME:-{REPO}}}"
exec python3 "$REPO/desktop/epigenesis_desktop.py" --repo "$REPO" "$@"
""",
        encoding="utf-8",
    )
    launcher.chmod(launcher.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    desktop_entry = out / "epigenesis.desktop"
    desktop_entry.write_text(
        f"""[Desktop Entry]
Type=Application
Name=Epigenesis
GenericName=DNA Simulation & Analysis Workbench
Comment=DNA-to-Development Simulation, Inverse DNA Compiler & Cryptographic Replay
Exec={launcher} %f
Icon={icon_path}
Terminal=false
Categories=Science;Biology;Education;
MimeType=text/plain;application/json;text/csv;
""",
        encoding="utf-8",
    )
    return out


def build_windows_bundle(dist_dir: Path) -> Path:
    """Create Windows double-clickable .bat + silent .vbs launchers and .ico icon."""
    out = dist_dir / "windows"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    write_ico_from_png(render_icon_png(256), out / "epigenesis.ico")
    shutil.copyfile(REPO / "desktop" / "epigenesis_desktop.py", out / "epigenesis_launcher.py")

    (out / "Epigenesis.bat").write_text(
        "@echo off\r\n"
        "setlocal\r\n"
        "if \"%EPIGENESIS_HOME%\"==\"\" set \"EPIGENESIS_HOME=%~dp0..\\..\"\r\n"
        "py -3.11 \"%~dp0epigenesis_launcher.py\" --repo \"%EPIGENESIS_HOME%\" %*\r\n"
        "if errorlevel 1 python \"%~dp0epigenesis_launcher.py\" --repo \"%EPIGENESIS_HOME%\" %*\r\n",
        encoding="utf-8",
    )
    (out / "Epigenesis.vbs").write_text(
        'Set shell = CreateObject("WScript.Shell")\r\n'
        'dir = CreateObject("Scripting.FileSystemObject").GetParentFolderName(WScript.ScriptFullName)\r\n'
        'shell.Run """" & dir & "\\Epigenesis.bat""", 0, False\r\n',
        encoding="utf-8",
    )
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--platform", choices=("auto", "macos", "linux", "windows", "all"), default="auto")
    parser.add_argument("--dist", type=Path, default=REPO / "dist")
    parser.add_argument("--install", action="store_true", help="on macOS, also copy Epigenesis.app to ~/Applications")
    args = parser.parse_args(argv)

    built: dict[str, str] = {}
    target = args.platform
    if target == "auto":
        sys_name = platform.system()
        target = "macos" if sys_name == "Darwin" else ("windows" if sys_name == "Windows" else "linux")

    if target in ("macos", "all"):
        install_dir = Path.home() / "Applications" if args.install else None
        built["macos"] = str(build_macos_app(args.dist, install_dir=install_dir))
    if target in ("linux", "all"):
        built["linux"] = str(build_linux_bundle(args.dist))
    if target in ("windows", "all"):
        built["windows"] = str(build_windows_bundle(args.dist))

    print(json.dumps({"built": built}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

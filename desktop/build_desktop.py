#!/usr/bin/env python3
"""Build the native Epigenesis macOS desktop application (.app + .dmg installer).

Zero third-party build dependencies:
  - Renders a macOS Human Interface Guidelines (HIG) compliant minimalistic icon
    (1024x1024 canvas, 824x824 centered continuous superellipse squircle with 100px
    optical padding, subtle Dock drop shadow, specular rim, and full 10-size .icns).
  - Compiles `desktop/macos/EpigenesisApp.swift` into a native AppKit/WKWebView
    `dist/Epigenesis.app` bundle with `runtime.json` and ad-hoc code signature.
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
# macOS HIG-compliant minimalistic icon generator (1024x1024 squircle + .icns)
# ---------------------------------------------------------------------------

def _png_chunk(tag: bytes, data: bytes) -> bytes:
    crc = zlib.crc32(tag + data) & 0xFFFFFFFF
    return struct.pack(">I", len(data)) + tag + data + struct.pack(">I", crc)


def _encode_rgba_png(size: int, rgba_rows: bytes) -> bytes:
    ihdr = struct.pack(">IIBBBBB", size, size, 8, 6, 0, 0, 0)
    idat = zlib.compress(rgba_rows, level=6)
    return b"\x89PNG\r\n\x1a\n" + _png_chunk(b"IHDR", ihdr) + _png_chunk(b"IDAT", idat) + _png_chunk(b"IEND", b"")


def _render_icon_numpy(size: int) -> bytes:
    """Fast vectorized macOS HIG icon renderer (exact analytical antialiasing)."""
    import numpy as np

    coords = (np.arange(size, dtype=np.float64) + 0.5) / size  # 0..1
    xx, yy = np.meshgrid(coords, coords)
    # macOS HIG standard grid: 824x824 centered squircle inside 1024x1024 canvas
    half = 824.0 / 2048.0  # 0.40234375 (leaves exact 100px margin on 1024px canvas)
    p = 5.25  # Apple Big Sur+ macOS continuous-curvature squircle (flat sides, smooth corners)

    # Coordinates relative to canvas center
    ux = (xx - 0.5) / half
    uy = (yy - 0.5) / half
    # Superellipse implicit field F = (|ux|^p + |uy|^p)^(1/p); F=1 on boundary
    r_se = np.power(np.abs(ux) ** p + np.abs(uy) ** p, 1.0 / p)
    # Approximate Euclidean signed distance in pixels (negative inside, positive outside)
    grad_approx = np.maximum(1e-6, np.power(np.abs(ux) ** (2 * p - 2) + np.abs(uy) ** (2 * p - 2), 0.5) * (r_se ** (1.0 - p)))
    sdf_px = ((r_se - 1.0) / grad_approx) * (half * size)

    # Subtle macOS Dock drop shadow inside the 100px transparent padding
    uy_sh = (yy - (0.5 + 14.0 / 1024.0)) / half
    r_sh = np.power(np.abs(ux) ** p + np.abs(uy_sh) ** p, 1.0 / p)
    grad_sh = np.maximum(1e-6, np.power(np.abs(ux) ** (2 * p - 2) + np.abs(uy_sh) ** (2 * p - 2), 0.5) * (r_sh ** (1.0 - p)))
    sh_px = ((r_sh - 1.0) / grad_sh) * (half * size)
    sigma = max(1.0, size * (22.0 / 1024.0))
    sh_alpha = np.where(sh_px > 0.0, 0.28 * np.exp(-0.5 * (sh_px / sigma) ** 2), 0.28)
    # Fade shadow cleanly to 0 before canvas edge so there is never clipping
    edge_dist = np.minimum(np.minimum(xx, 1.0 - xx), np.minimum(yy, 1.0 - yy)) * size
    sh_alpha *= np.clip(edge_dist / max(1.0, size * 0.025), 0.0, 1.0)

    # Tile alpha (crisp 1.5px analytical antialiasing)
    tile_a = np.clip(0.5 - sdf_px / 1.4, 0.0, 1.0)

    # Understated macOS lab-slate vertical gradient (#2A3038 top -> #161A20 bottom)
    vy = np.clip((yy - (0.5 - half)) / (2.0 * half), 0.0, 1.0)
    bg_r = 42.0 * (1.0 - vy) + 21.0 * vy
    bg_g = 48.0 * (1.0 - vy) + 25.0 * vy
    bg_b = 56.0 * (1.0 - vy) + 32.0 * vy

    # Subtle radial top illumination (classic Apple macOS material polish)
    rad = np.hypot((xx - 0.5) * 1.15, (yy - 0.22) * 1.4)
    sheen = np.clip(1.0 - rad / 0.72, 0.0, 1.0) ** 2 * 10.0
    bg_r = np.minimum(255.0, bg_r + sheen)
    bg_g = np.minimum(255.0, bg_g + sheen)
    bg_b = np.minimum(255.0, bg_b + sheen * 1.1)

    # Crisp 1px top specular hairline + subtle bottom rim inside squircle
    rim_band = np.exp(-((sdf_px + max(1.0, size * 0.0028)) / max(0.7, size * 0.0022)) ** 2)
    top_weight = np.clip(1.0 - vy * 1.35, 0.0, 1.0)
    rim_boost = rim_band * (18.0 * top_weight + 4.0)
    bg_r = np.minimum(255.0, bg_r + rim_boost)
    bg_g = np.minimum(255.0, bg_g + rim_boost)
    bg_b = np.minimum(255.0, bg_b + rim_boost)

    # Minimalistic laboratory double-helix glyph centered in inner optical safe zone
    # Safe zone height: 0.25..0.75 (50% of canvas, ~62% of squircle)
    top_y, bot_y = 0.25, 0.75
    span_y = bot_y - top_y
    amp = 0.115
    stroke_px = max(1.5, size * 0.036)
    rung_half_px = max(1.0, size * 0.0145)

    # Sample helix polylines finely for exact Euclidean distance + true round endcaps
    ts = np.linspace(0.0, 1.0, 96)
    angles = (ts - 0.5) * (math.pi * 1.45)  # gracefully open at top & bottom, crisp central crossover
    ys_s = top_y + ts * span_y
    x1_s = 0.5 + amp * np.sin(angles)
    x2_s = 0.5 - amp * np.sin(angles)

    d1_min = np.full_like(xx, 1e9)
    d2_min = np.full_like(xx, 1e9)
    for k in range(len(ts) - 1):
        ax1, ay, bx1, by = x1_s[k], ys_s[k], x1_s[k + 1], ys_s[k + 1]
        vx1, vy1 = bx1 - ax1, by - ay
        l2_1 = vx1 * vx1 + vy1 * vy1
        u1 = np.clip(((xx - ax1) * vx1 + (yy - ay) * vy1) / l2_1, 0.0, 1.0)
        d1_min = np.minimum(d1_min, np.hypot(xx - (ax1 + u1 * vx1), yy - (ay + u1 * vy1)))

        ax2, bx2 = x2_s[k], x2_s[k + 1]
        vx2 = bx2 - ax2
        l2_2 = vx2 * vx2 + vy1 * vy1
        u2 = np.clip(((xx - ax2) * vx2 + (yy - ay) * vy1) / l2_2, 0.0, 1.0)
        d2_min = np.minimum(d2_min, np.hypot(xx - (ax2 + u2 * vx2), yy - (ay + u2 * vy1)))

    d1_px = d1_min * size
    d2_px = d2_min * size

    # Four balanced laboratory-cobalt base-pair rungs (two in upper loop, two in lower loop)
    rung_mask = np.zeros_like(xx)
    for rt in (0.12, 0.30, 0.70, 0.88):
        ry = top_y + rt * span_y
        ra = (rt - 0.5) * (math.pi * 1.45)
        rx1 = 0.5 + amp * math.sin(ra)
        rx2 = 0.5 - amp * math.sin(ra)
        lo_x, hi_x = min(rx1, rx2) + 0.010, max(rx1, rx2) - 0.010
        if hi_x > lo_x + 0.012:
            x_clamped = np.clip(xx, lo_x, hi_x)
            dr_px = np.hypot((xx - x_clamped) * size, (yy - ry) * size)
            rm = np.clip(0.5 - (dr_px - rung_half_px) / 1.3, 0.0, 1.0)
            rung_mask = np.maximum(rung_mask, rm)

    # Back strand (subtle depth) vs front strand (pure platinum)
    front_mask = np.clip(0.5 - (d1_px - stroke_px * 0.5) / 1.35, 0.0, 1.0)
    back_mask = np.clip(0.5 - (d2_px - stroke_px * 0.5) / 1.35, 0.0, 1.0)

    # Composite layers inside tile: bg -> rungs -> back strand -> front strand
    # Rungs: restrained laboratory cobalt (#5C9BFF)
    r = bg_r * (1.0 - rung_mask) + 92.0 * rung_mask
    g = bg_g * (1.0 - rung_mask) + 155.0 * rung_mask
    b = bg_b * (1.0 - rung_mask) + 255.0 * rung_mask

    # Back strand: soft cool silver (#B9C3D0)
    r = r * (1.0 - back_mask) + 185.0 * back_mask
    g = g * (1.0 - back_mask) + 195.0 * back_mask
    b = b * (1.0 - back_mask) + 208.0 * back_mask

    # Front strand: crisp off-white (#F3F5F8)
    r = r * (1.0 - front_mask) + 243.0 * front_mask
    g = g * (1.0 - front_mask) + 245.0 * front_mask
    b = b * (1.0 - front_mask) + 248.0 * front_mask

    # Alpha-composite tile over transparent drop shadow
    out_a = tile_a + sh_alpha * (1.0 - tile_a)
    safe_a = np.where(out_a > 1e-6, out_a, 1.0)
    out_r = np.where(out_a > 1e-6, (r * tile_a) / safe_a, 0.0)
    out_g = np.where(out_a > 1e-6, (g * tile_a) / safe_a, 0.0)
    out_b = np.where(out_a > 1e-6, (b * tile_a) / safe_a, 0.0)

    rgba = np.empty((size, size, 4), dtype=np.uint8)
    rgba[..., 0] = np.clip(np.rint(out_r), 0, 255).astype(np.uint8)
    rgba[..., 1] = np.clip(np.rint(out_g), 0, 255).astype(np.uint8)
    rgba[..., 2] = np.clip(np.rint(out_b), 0, 255).astype(np.uint8)
    rgba[..., 3] = np.clip(np.rint(out_a * 255.0), 0, 255).astype(np.uint8)

    # Prepend PNG filter byte 0 to each scanline
    buf = bytearray()
    for row in rgba:
        buf.append(0)
        buf.extend(row.tobytes())
    return _encode_rgba_png(size, bytes(buf))


def render_icon_png(size: int) -> bytes:
    """Render a macOS HIG-compliant minimalistic Epigenesis icon PNG of `size x size`."""
    try:
        import numpy  # noqa: F401
        return _render_icon_numpy(size)
    except ImportError:
        pass

    # Pure-Python fallback (identical macOS HIG 824/1024 squircle geometry)
    half = 824.0 / 2048.0
    p = 5.25
    top_y, bot_y = 0.25, 0.75
    span_y = bot_y - top_y
    amp = 0.115
    stroke_px = max(1.5, size * 0.036)
    rung_half_px = max(1.0, size * 0.0145)
    sigma = max(1.0, size * (22.0 / 1024.0))
    rows = bytearray()

    for yi in range(size):
        rows.append(0)
        y = (yi + 0.5) / size
        uy = (y - 0.5) / half
        uy_sh = (y - (0.5 + 14.0 / 1024.0)) / half
        vy = max(0.0, min(1.0, (y - (0.5 - half)) / (2.0 * half)))
        t = max(0.0, min(1.0, (y - top_y) / span_y))
        angle = (t - 0.5) * (math.pi * 1.45)
        x1 = 0.5 + amp * math.sin(angle)
        x2 = 0.5 - amp * math.sin(angle)
        slope_norm = math.hypot(1.0, (amp * (math.pi * 1.45) / span_y) * math.cos(angle))
        cap_pad = max(max(0.0, top_y - y), max(0.0, y - bot_y)) * size
        lo_x = min(x1, x2) + 0.010
        hi_x = max(x1, x2) - 0.010
        auy_p = abs(uy) ** p
        auy_sh_p = abs(uy_sh) ** p

        for xi in range(size):
            x = (xi + 0.5) / size
            ux = (x - 0.5) / half
            aux_p = abs(ux) ** p
            r_se = (aux_p + auy_p) ** (1.0 / p)
            grad = max(1e-6, ((abs(ux) ** (2 * p - 2) + abs(uy) ** (2 * p - 2)) ** 0.5) * (r_se ** (1.0 - p)))
            sdf_px = ((r_se - 1.0) / grad) * (half * size)
            tile_a = max(0.0, min(1.0, 0.5 - sdf_px / 1.4))

            r_sh = (aux_p + auy_sh_p) ** (1.0 / p)
            grad_sh = max(1e-6, ((abs(ux) ** (2 * p - 2) + abs(uy_sh) ** (2 * p - 2)) ** 0.5) * (r_sh ** (1.0 - p)))
            sh_px = ((r_sh - 1.0) / grad_sh) * (half * size)
            sh_a = 0.28 * math.exp(-0.5 * (sh_px / sigma) ** 2) if sh_px > 0.0 else 0.28
            edge_dist = min(x, 1.0 - x, y, 1.0 - y) * size
            sh_a *= max(0.0, min(1.0, edge_dist / max(1.0, size * 0.025)))

            out_a = tile_a + sh_a * (1.0 - tile_a)
            if out_a <= 1e-6:
                rows.extend((0, 0, 0, 0))
                continue

            r = 42.0 * (1.0 - vy) + 21.0 * vy
            g = 48.0 * (1.0 - vy) + 25.0 * vy
            b = 56.0 * (1.0 - vy) + 32.0 * vy

            if tile_a > 0.0:
                d1_px = math.hypot((abs(x - x1) / slope_norm) * size, cap_pad)
                d2_px = math.hypot((abs(x - x2) / slope_norm) * size, cap_pad)
                rm = 0.0
                if hi_x > lo_x + 0.015:
                    xc = max(lo_x, min(hi_x, x))
                    for rt in (0.19, 0.50, 0.81):
                        dr_px = math.hypot((x - xc) * size, (y - (top_y + rt * span_y)) * size)
                        rm = max(rm, max(0.0, min(1.0, 0.5 - (dr_px - rung_half_px) / 1.3)))
                r = r * (1.0 - rm) + 92.0 * rm
                g = g * (1.0 - rm) + 155.0 * rm
                b = b * (1.0 - rm) + 255.0 * rm
                bm = max(0.0, min(1.0, 0.5 - (d2_px - stroke_px * 0.5) / 1.35))
                r = r * (1.0 - bm) + 185.0 * bm
                g = g * (1.0 - bm) + 195.0 * bm
                b = b * (1.0 - bm) + 208.0 * bm
                fm = max(0.0, min(1.0, 0.5 - (d1_px - stroke_px * 0.5) / 1.35))
                r = r * (1.0 - fm) + 243.0 * fm
                g = g * (1.0 - fm) + 245.0 * fm
                b = b * (1.0 - fm) + 248.0 * fm

            rows.extend((
                int(round((r * tile_a) / out_a)),
                int(round((g * tile_a) / out_a)),
                int(round((b * tile_a) / out_a)),
                int(round(out_a * 255.0)),
            ))
    return _encode_rgba_png(size, bytes(rows))


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

    swift_sources = sorted((REPO / "desktop" / "macos").glob("*.swift"))
    bin_path = macos_dir / "Epigenesis"
    subprocess.run(
        [
            swiftc,
            "-O",
            *[str(src) for src in swift_sources],
            "-o",
            str(bin_path),
            "-framework",
            "SwiftUI",
            "-framework",
            "AppKit",
        ],
        check=True,
    )

    iconutil = shutil.which("iconutil")
    sips = shutil.which("sips")
    if iconutil is not None:
        with tempfile.TemporaryDirectory() as tmp:
            master_png = Path(tmp) / "master_1024.png"
            # Render full 1024x1024 macOS HIG master once (fast with NumPy; downsampled cleanly with sips)
            master_png.write_bytes(render_icon_png(1024))
            iconset = Path(tmp) / "AppIcon.iconset"
            iconset.mkdir()
            specs = (
                (16, "icon_16x16.png"),
                (32, "icon_16x16@2x.png"),
                (32, "icon_32x32.png"),
                (64, "icon_32x32@2x.png"),
                (128, "icon_128x128.png"),
                (256, "icon_128x128@2x.png"),
                (256, "icon_256x256.png"),
                (512, "icon_256x256@2x.png"),
                (512, "icon_512x512.png"),
                (1024, "icon_512x512@2x.png"),
            )
            for sz, name in specs:
                out_png = iconset / name
                if sz == 1024:
                    shutil.copyfile(master_png, out_png)
                elif sips is not None:
                    subprocess.run(
                        [sips, "-z", str(sz), str(sz), str(master_png), "--out", str(out_png)],
                        check=True,
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL,
                    )
                else:
                    out_png.write_bytes(render_icon_png(sz))
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
  <key>LSMinimumSystemVersion</key><string>13.0</string>
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
    parser.add_argument("--export-icon", type=Path, help="also write a 512x512 HIG PNG icon to this path")
    args = parser.parse_args(argv)

    if args.export_icon is not None:
        args.export_icon.parent.mkdir(parents=True, exist_ok=True)
        args.export_icon.write_bytes(render_icon_png(512))

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

from __future__ import annotations

import json
from pathlib import Path
import platform
import shutil
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "desktop"))

import build_desktop as bd  # noqa: E402


class DesktopBuildTests(unittest.TestCase):
    def test_icon_png_and_ico_generation(self) -> None:
        png = bd.render_icon_png(64)
        self.assertTrue(png.startswith(b"\x89PNG\r\n\x1a\n"))
        self.assertGreater(len(png), 200)
        with tempfile.TemporaryDirectory() as tmp:
            ico = Path(tmp) / "epigenesis.ico"
            bd.write_ico_from_png(png, ico)
            raw = ico.read_bytes()
            self.assertEqual(raw[:4], b"\x00\x00\x01\x00")
            self.assertIn(b"\x89PNG", raw)

    def test_linux_and_windows_bundles_build_cleanly(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            dist = Path(tmp) / "dist"
            linux_dir = bd.build_linux_bundle(dist)
            self.assertTrue((linux_dir / "epigenesis-desktop").is_file())
            self.assertTrue((linux_dir / "epigenesis.desktop").is_file())
            self.assertTrue((linux_dir / "epigenesis.png").is_file())
            self.assertIn("Type=Application", (linux_dir / "epigenesis.desktop").read_text(encoding="utf-8"))

            win_dir = bd.build_windows_bundle(dist)
            self.assertTrue((win_dir / "Epigenesis.bat").is_file())
            self.assertTrue((win_dir / "Epigenesis.vbs").is_file())
            self.assertTrue((win_dir / "epigenesis_launcher.py").is_file())
            self.assertTrue((win_dir / "epigenesis.ico").is_file())

    @unittest.skipUnless(platform.system() == "Darwin" and shutil.which("swiftc"), "macOS swiftc required")
    def test_macos_native_app_bundle_compiles_cleanly(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            dist = Path(tmp) / "dist"
            app = bd.build_macos_app(dist)
            self.assertTrue((app / "Contents" / "MacOS" / "Epigenesis").is_file())
            self.assertTrue((app / "Contents" / "Info.plist").is_file())
            self.assertTrue((app / "Contents" / "Resources" / "runtime.json").is_file())
            if shutil.which("iconutil"):
                self.assertTrue((app / "Contents" / "Resources" / "AppIcon.icns").is_file())
            cfg = json.loads((app / "Contents" / "Resources" / "runtime.json").read_text(encoding="utf-8"))
            self.assertEqual(Path(cfg["repo"]).resolve(), ROOT.resolve())


if __name__ == "__main__":
    unittest.main()

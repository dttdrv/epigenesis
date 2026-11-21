from __future__ import annotations

import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "desktop"))

import build_desktop as bd  # noqa: E402


class DesktopBuildTests(unittest.TestCase):
    def test_icon_png_generation(self) -> None:
        png = bd.render_icon_png(64)
        self.assertTrue(png.startswith(b"\x89PNG\r\n\x1a\n"))
        self.assertGreater(len(png), 200)

    def test_linux_and_windows_report_coming_soon_cleanly(self) -> None:
        for plat in ("linux", "windows"):
            self.assertEqual(bd.main(["--platform", plat]), 0)

    @unittest.skipUnless(platform.system() == "Darwin" and shutil.which("swiftc"), "macOS swiftc required")
    def test_macos_native_app_and_dmg_build_cleanly(self) -> None:
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

            if shutil.which("hdiutil"):
                dmg = bd.build_macos_dmg(app, dist)
                self.assertTrue(dmg.is_file())
                self.assertEqual(dmg.name, bd.DMG_NAME)
                self.assertGreater(dmg.stat().st_size, 10_000)
                ver = subprocess.run(["hdiutil", "verify", str(dmg)], capture_output=True, text=True)
                self.assertEqual(ver.returncode, 0, ver.stderr)


if __name__ == "__main__":
    unittest.main()

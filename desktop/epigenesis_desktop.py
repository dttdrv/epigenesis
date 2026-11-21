#!/usr/bin/env python3
"""Cross-platform dedicated desktop application launcher for Epigenesis Studio.

- macOS (primary): Launches the native Swift/AppKit/WKWebView `Epigenesis.app` bundle
  (building it automatically if `swiftc` is available and not yet built), with clean
  fallback to standalone app-mode window if needed.
- Windows & Linux: Starts the local loopback `epigenesis_studio` server on an ephemeral
  port and opens a dedicated standalone desktop window (Chrome/Edge/Chromium/Brave
  `--app=` mode with no browser tabs or toolbar, falling back to default browser) and
  terminates the local server cleanly when closed.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import webbrowser


REPO = Path(__file__).resolve().parents[1]


def find_studio_python(repo: Path) -> Path:
    """Locate a Python 3.11+ interpreter with NumPy installed for Epigenesis Studio."""
    candidates: list[Path] = []
    home = Path.home()
    if platform.system() == "Windows":
        candidates.extend([
            home / ".local/share/uv/tools/epigenesis-studio/Scripts/python.exe",
            repo / ".venv/Scripts/python.exe",
        ])
    else:
        candidates.extend([
            home / ".local/share/uv/tools/epigenesis-studio/bin/python",
            repo / ".venv/bin/python",
        ])
    candidates.append(Path(sys.executable))
    for cand in candidates:
        if not cand.is_file():
            continue
        probe = subprocess.run(
            [str(cand), "-c", "import sys, numpy; assert sys.version_info >= (3, 11)"],
            capture_output=True,
            check=False,
        )
        if probe.returncode == 0:
            return cand
    return Path(sys.executable)


def _find_chromium_app_browser() -> list[str] | None:
    system = platform.system()
    if system == "Windows":
        program_files = [
            Path(os.environ.get("PROGRAMFILES", r"C:\Program Files")),
            Path(os.environ.get("PROGRAMFILES(X86)", r"C:\Program Files (x86)")),
            Path(os.environ.get("LOCALAPPDATA", "")),
        ]
        rel_paths = [
            "Microsoft/Edge/Application/msedge.exe",
            "Google/Chrome/Application/chrome.exe",
            "BraveSoftware/Brave-Browser/Application/brave.exe",
        ]
        for base in program_files:
            for rel in rel_paths:
                cand = base / rel
                if cand.is_file():
                    return [str(cand)]
    elif system == "Linux":
        for name in ("chromium", "chromium-browser", "google-chrome", "microsoft-edge", "brave-browser"):
            found = shutil.which(name)
            if found:
                return [found]
    elif system == "Darwin":
        for app_bin in (
            "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
            "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge",
            "/Applications/Brave Browser.app/Contents/MacOS/Brave Browser",
            "/Applications/Chromium.app/Contents/MacOS/Chromium",
        ):
            if Path(app_bin).is_file():
                return [app_bin]
    return None


def launch_cross_platform_window(repo: Path, *, target_file: Path | None = None, port: int = 0) -> int:
    """Start loopback Epigenesis Studio and open a dedicated standalone window."""
    py = find_studio_python(repo)
    env = dict(os.environ)
    py_parts = [str(repo), str(repo / "studio")]
    if env.get("PYTHONPATH"):
        py_parts.append(env["PYTHONPATH"])
    env["PYTHONPATH"] = os.pathsep.join(py_parts)
    env["PYTHONUNBUFFERED"] = "1"

    cmd = [str(py), "-m", "epigenesis_studio.server", "--port", str(port), "--no-browser"]
    if target_file is not None:
        cmd.append(str(target_file))

    proc = subprocess.Popen(
        cmd,
        cwd=repo,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    url: str | None = None
    assert proc.stdout is not None
    try:
        while proc.poll() is None:
            line = proc.stdout.readline()
            if not line:
                break
            m = re.search(r"http://127\.0\.0\.1:\d+/", line)
            if m:
                url = m.group(0)
                break
        if url is None:
            err = proc.stderr.read() if proc.stderr else ""
            print(f"epigenesis-desktop: backend failed to start:\n{err}", file=sys.stderr)
            return 1

        browser_cmd = _find_chromium_app_browser()
        if browser_cmd is not None:
            win_proc = subprocess.Popen(
                [*browser_cmd, f"--app={url}", "--window-size=1280,860"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            win_proc.wait()
            return 0

        print(f"Epigenesis Studio running at {url} (Ctrl+C to stop)", flush=True)
        webbrowser.open(url)
        proc.wait()
        return 0
    except KeyboardInterrupt:
        return 0
    finally:
        if proc.poll() is None:
            proc.terminate()
            try:
                proc.wait(timeout=3)
            except subprocess.TimeoutExpired:
                proc.kill()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", nargs="?", type=Path, help="optional DNA FASTA, edit table, or receipt to open")
    parser.add_argument("--repo", type=Path, default=REPO, help="Epigenesis repository root")
    parser.add_argument("--port", type=int, default=0, help="loopback port (default: automatic free port)")
    parser.add_argument("--no-native", action="store_true", help="skip macOS native .app bundle and use standalone browser window")
    args = parser.parse_args(argv)
    repo = args.repo.resolve()

    if platform.system() == "Darwin" and not args.no_native:
        app_bundle = repo / "dist" / "Epigenesis.app"
        if not (app_bundle / "Contents" / "MacOS" / "Epigenesis").is_file() and shutil.which("swiftc"):
            build_script = repo / "desktop" / "build_desktop.py"
            if build_script.is_file():
                subprocess.run([sys.executable, str(build_script), "--platform", "macos"], cwd=repo, check=False)
        if (app_bundle / "Contents" / "MacOS" / "Epigenesis").is_file():
            open_cmd = ["open", "-a", str(app_bundle)]
            if args.file is not None:
                open_cmd.append(str(args.file.resolve()))
            return subprocess.call(open_cmd)

    return launch_cross_platform_window(repo, target_file=args.file, port=args.port)


if __name__ == "__main__":
    raise SystemExit(main())

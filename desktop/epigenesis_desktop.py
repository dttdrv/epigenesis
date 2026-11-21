#!/usr/bin/env python3
"""Launch the native Epigenesis desktop application (macOS-first; Linux & Windows coming soon).

- macOS (Primary): Opens the installed native Swift/AppKit/WKWebView `Epigenesis.app`
  (checking `/Applications`, `~/Applications`, and `dist/Epigenesis.app`, and building
  it automatically if Xcode Command Line Tools `swiftc` is present).
- Linux & Windows (Coming Soon™): Falls back cleanly to starting the local `127.0.0.1`
  Epigenesis Studio workbench in an app-mode window or default browser.
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


def _find_installed_macos_app(repo: Path) -> Path | None:
    for cand in (
        Path("/Applications/Epigenesis.app"),
        Path.home() / "Applications" / "Epigenesis.app",
        repo / "dist" / "Epigenesis.app",
    ):
        if (cand / "Contents" / "MacOS" / "Epigenesis").is_file():
            return cand
    return None


def launch_local_workbench(repo: Path, *, target_file: Path | None = None, port: int = 0) -> int:
    """Fallback local workbench launcher for non-macOS systems while native apps are coming soon."""
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

        print(
            f"Epigenesis Studio running locally at {url} "
            "(native Linux/Windows desktop installers coming soon; Ctrl+C to stop)",
            flush=True,
        )
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
    parser.add_argument("--port", type=int, default=0, help="loopback port for fallback mode")
    parser.add_argument("--no-native", action="store_true", help="skip macOS native .app and open in local browser")
    args = parser.parse_args(argv)
    repo = args.repo.resolve()

    if platform.system() == "Darwin" and not args.no_native:
        app_bundle = _find_installed_macos_app(repo)
        if app_bundle is None and shutil.which("swiftc"):
            build_script = repo / "desktop" / "build_desktop.py"
            if build_script.is_file():
                subprocess.run(
                    [sys.executable, str(build_script), "--platform", "macos", "--no-dmg"],
                    cwd=repo,
                    check=False,
                )
            app_bundle = _find_installed_macos_app(repo)
        if app_bundle is not None:
            open_cmd = ["open", "-a", str(app_bundle)]
            if args.file is not None:
                open_cmd.append(str(args.file.resolve()))
            return subprocess.call(open_cmd)

    return launch_local_workbench(repo, target_file=args.file, port=args.port)


if __name__ == "__main__":
    raise SystemExit(main())

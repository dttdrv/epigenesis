"""Build the English and Bulgarian paper PDFs from their LaTeX sources.

Uses Tectonic (a self-contained XeTeX engine; https://tectonic-typesetting.org). Set TECTONIC to its path,
or put it on PATH, or place it at build/tools/tectonic/tectonic. The sources use the macOS system fonts
Palatino and Menlo. Usage: python3 docs/paper/build.py
"""
import os, shutil, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
tectonic = os.environ.get("TECTONIC") or shutil.which("tectonic") or str(HERE.parents[1] / "build/tools/tectonic/tectonic")
if not Path(tectonic).is_file() and not shutil.which(tectonic):
    sys.exit("tectonic not found; set TECTONIC")
for name in ("epigenesis-technical-paper", "epigenesis-technical-paper.bg"):
    subprocess.run([tectonic, "-X", "compile", f"{name}.tex"], cwd=HERE, check=True)
    print(f"built {name}.pdf")

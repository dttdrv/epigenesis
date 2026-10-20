"""Build a clean, standalone release folder (default ../epigenesis-release) from the working tree.
It is the setuptools source archive (what MANIFEST.in publishes) plus tools/release/ (launcher, short README, status);
the full README follows the quick start in README.md. No build/, output/, caches or git data."""
import shutil, sys, tarfile, tempfile
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
dest = Path(sys.argv[1]) if sys.argv[1:] else ROOT.parent / 'epigenesis-release'
if dest.exists():
    sys.exit(f'{dest} exists; remove it first')
with tempfile.TemporaryDirectory() as tmp:
    subprocess.run([sys.executable, '-c', f'import setuptools.build_meta as b; b.build_sdist({tmp!r})'], check=True, cwd=ROOT, capture_output=True)
    (archive,) = Path(tmp).glob('*.tar.gz')
    with tarfile.open(archive) as t:
        t.extractall(tmp, filter='data')
    (src,) = [p for p in Path(tmp).iterdir() if p.is_dir()]
    for junk in src.glob('*.egg-info'):
        shutil.rmtree(junk)
    # quick start first; the full project README follows unchanged (a release test pins its contents)
    quick = (ROOT / 'tools/release/README.md').read_text().replace('](../../', '](').replace('](STATUS.md)', '](docs/STATUS.md)')
    full = (ROOT / 'README.md').read_text().replace('[CI workflow](.github/workflows/ci.yml)', 'CI workflow (in the source repository)')
    (src / 'README.md').write_text(quick.rstrip() + '\n\n---\n\n' + full)
    shutil.copy(ROOT / 'tools/release/epigenesis.py', src / 'epigenesis.py')
    status = (ROOT / 'tools/release/STATUS.md').read_text()
    (src / 'docs/STATUS.md').write_text(status.replace('](../../docs/', '](').replace('](../../', '](../'))
    shutil.copytree(src, dest)
print(dest)

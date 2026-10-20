"""Verify the SHA-256 digests in every FREEZE*.json under research/frozen-tests.

A freeze lists the plan, scripts and sealed inputs fixed before outcomes were scored. Files copied into this
folder are checked directly. Large sealed inputs (for example derived/PREDICTIONS.npz) stay in the untracked
build tree; pass --build to check those too. A mismatch exits nonzero; absent files are reported, not passed.
PUBLICATION.json records edited public copies separately from their original frozen identities.
"""
import argparse, hashlib, json, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent / "frozen-tests"
ap = argparse.ArgumentParser()
ap.add_argument("--build", type=Path, help="build/neural-predictor-reviewed-20261005, to resolve files not copied here")
args = ap.parse_args()

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

bad = absent = ok = copies = 0
publication = HERE.parent / "PUBLICATION.json"
exports = {}
if publication.is_file():
    record = json.loads(publication.read_text())
    if record.get("format") != "epigenesis-publication-copies-v1":
        sys.exit("unsupported publication format")
    exports = record["files"]
    for rel, entry in exports.items():
        path = HERE / rel
        if not path.is_file() or sha256(path) != entry["published_sha256"]:
            bad += 1
            print(f"MISMATCH publication copy {rel}")
for freeze in sorted(HERE.glob("*/FREEZE*.json")):
    test = freeze.parent.name
    files = json.loads(freeze.read_text()).get("files")
    if files is None:
        print(f"skipped {test}/{freeze.name}: no file digests in this freeze")
        continue
    for rel, want in files.items():
        for root in (freeze.parent, args.build / test if args.build else None):
            path = root / rel if root else None
            if path and path.is_file():
                break
        else:
            absent += 1
            continue
        actual = sha256(path)
        entry = exports.get(f"{test}/{rel}", {})
        if actual == want:
            ok += 1
        elif (path == freeze.parent / rel and entry.get("original_sha256") == want
              and entry.get("published_sha256") == actual):
            copies += 1
        else:
            bad += 1
            print(f"MISMATCH {test}/{rel}")
print(f"{ok} original digests verified, {copies} declared publication copies, "
      f"{bad} mismatched, {absent} not present locally")
sys.exit(1 if bad else 0)

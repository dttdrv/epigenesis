#!/usr/bin/env python3
"""Unified command-line interface for Epigenesis.

Single entry point across enhancer circuit scanning, variant-to-GRN-to-CX3D
bridging, Kreimer 0-72h / Siraj reporter prediction, Studio, and verification.
Uses Python standard library only.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys


def _find_repo(explicit: Path | None = None) -> Path:
    candidates: list[Path] = []
    if explicit is not None:
        candidates.append(explicit)
    env_home = os.environ.get("EPIGENESIS_HOME")
    if env_home:
        candidates.append(Path(env_home))
    candidates.extend([Path(__file__).resolve().parent, Path.cwd()])
    for root in candidates:
        resolved = root.resolve()
        if (resolved / "examples" / "neural-circuit" / "circuit.py").is_file():
            return resolved
    raise FileNotFoundError(
        "could not locate the Epigenesis repository checkout (set EPIGENESIS_HOME or pass --repo)"
    )


def _unique_out(base: Path) -> Path:
    if not base.exists():
        base.parent.mkdir(parents=True, exist_ok=True)
        return base
    for idx in range(1, 1000):
        cand = base.parent / f"{base.name}-{idx}"
        if not cand.exists():
            cand.parent.mkdir(parents=True, exist_ok=True)
            return cand
    raise FileExistsError(f"all default output directories under {base.parent} already exist")


def _resolve_out(explicit: Path | None, default_name: str) -> Path:
    if explicit is not None:
        return explicit
    return _unique_out(Path.cwd() / "epigenesis-out" / default_name)


def _run_py(cmd: list[str | Path], *, cwd: Path | None = None) -> int:
    return subprocess.call([sys.executable, *[str(part) for part in cmd]], cwd=cwd)


def _cmd_studio(args: argparse.Namespace) -> int:
    try:
        from epigenesis_studio import server
    except ImportError:
        repo = _find_repo(args.repo)
        studio_pkg = repo / "studio"
        if not (studio_pkg / "epigenesis_studio" / "server.py").is_file():
            print("epigenesis studio package is not installed; run `uv tool install -e studio`", file=sys.stderr)
            return 1
        env = dict(os.environ)
        py_parts = [str(repo), str(studio_pkg)]
        if env.get("PYTHONPATH"):
            py_parts.append(env["PYTHONPATH"])
        env["PYTHONPATH"] = os.pathsep.join(py_parts)
        cmd = [sys.executable, "-m", "epigenesis_studio.server"]
        if args.file is not None:
            cmd.append(str(args.file))
        if args.port is not None:
            cmd.extend(["--port", str(args.port)])
        if args.no_browser:
            cmd.append("--no-browser")
        return subprocess.call(cmd, env=env)
    argv = []
    if args.file is not None:
        argv.append(str(args.file))
    if args.port is not None:
        argv.extend(["--port", str(args.port)])
    if args.no_browser:
        argv.append("--no-browser")
    saved = sys.argv
    sys.argv = ["epigenesis-studio", *argv]
    try:
        rc = server.main()
        return 0 if rc is None else (int(rc) if isinstance(rc, int) else 1)
    finally:
        sys.argv = saved


def _cmd_circuit(args: argparse.Namespace) -> int:
    repo = _find_repo(args.repo)
    script = repo / "examples/neural-circuit/circuit.py"
    if args.benchmark:
        out = _resolve_out(args.output, "circuit-benchmark")
        cmd: list[str | Path] = [
            script, "benchmark", "--protein", args.protein, "--target", args.target,
            "--signal", str(args.signal), "--output", out,
        ]
        rc = _run_py(cmd)
        if rc == 0:
            print(f"\nBenchmark report: {out / 'benchmark.json'}")
        return rc
    if bool(args.reference) != bool(args.variant):
        print("circuit requires both --reference and --variant when comparing custom pairs", file=sys.stderr)
        return 2
    if args.reference and args.variant:
        if args.with_construction:
            out = _resolve_out(args.output, "circuit-bridge")
            cmd = [
                script, "bridge", "--reference", args.reference, "--variant", args.variant,
                "--protein", args.protein, "--target", args.target, "--signal", str(args.signal),
                "--duration", str(args.duration), "--output", out,
            ]
            rc = _run_py(cmd)
            if rc == 0:
                print(f"\nConstruction bridge report: {out / 'bridge.json'}")
            return rc
        out = _resolve_out(args.output, "circuit-compare")
        cmd = [
            script, "compare", "--reference", args.reference, "--variant", args.variant,
            "--protein", args.protein, "--target", args.target, "--signal", str(args.signal),
            "--output", out,
        ]
        rc = _run_py(cmd)
        if rc == 0:
            print(f"\nComparison report: {out / 'compare.json'}")
        return rc
    out = _resolve_out(args.output, "neural-circuit")
    cmd = [
        script, "demo", "--locus", args.locus, "--target", args.target,
        "--signal", str(args.signal), "--protein", args.protein, "--output", out,
    ]
    if args.with_construction:
        cmd.append("--with-construction")
    rc = _run_py(cmd)
    if rc == 0:
        print(f"\nOpen {out / 'index.html'} in your browser.")
    return rc


def _cmd_scan(args: argparse.Namespace) -> int:
    repo = _find_repo(args.repo)
    script = repo / "examples/neural-circuit/circuit.py"
    out = _resolve_out(args.output, "enhancer-scan")
    cmd: list[str | Path] = [
        script, "scan", "--fasta", args.fasta, "--protein", args.protein,
        "--target", args.target, "--signal", str(args.signal), "--output", out,
    ]
    rc = _run_py(cmd)
    if rc == 0:
        print(f"\nScan output: {out / 'scan.json'}")
    return rc


def _cmd_predict(args: argparse.Namespace) -> int:
    repo = _find_repo(args.repo)
    legacy = args.activity or args.attention
    if args.budget is not None and legacy:
        print("--budget belongs to the default Huber-linear and Kreimer models", file=sys.stderr)
        return 2
    if args.families and not legacy:
        print("--families requires --activity or --attention", file=sys.stderr)
        return 2
    if args.device is not None and not args.attention:
        print("--device requires --attention", file=sys.stderr)
        return 2
    if args.kreimer:
        script = repo / "studio/epigenesis_studio/kreimer.py"
    elif legacy:
        script = repo / "examples/neural-predictor" / ("predictor.py" if args.attention else "activity.py")
    else:
        script = repo / "studio/epigenesis_studio/model/predict.py"
    cmd: list[str | Path] = [script, "--input", args.input, "--output", args.output]
    if args.budget is not None:
        cmd.extend(["--budget", str(args.budget)])
    if args.families:
        cmd.extend(["--families", str(args.families)])
    if args.attention:
        cmd.extend(["--device", args.device or "cpu"])
    return _run_py(cmd)


def _cmd_evidence(args: argparse.Namespace) -> int:
    repo = _find_repo(args.repo)
    script = repo / "examples/neural-evidence/evidence.py"
    out = _resolve_out(args.output, "neural-evidence")
    rc = _run_py([script, "demo", "--output", out])
    if rc == 0:
        rc = _run_py([script, "verify", "--output", out])
    if rc == 0:
        print(f"\nOpen {out / 'index.html'} in your browser.")
    return rc


def _cmd_patterning(args: argparse.Namespace) -> int:
    repo = _find_repo(args.repo)
    script = repo / "examples/neural-patterning/patterning.py"
    out = _resolve_out(args.output, "neural-patterning")
    rc = _run_py([script, "demo", "--output", out])
    if rc == 0:
        print(f"\nOpen {out / 'index.html'} in your browser.")
    return rc


def _cmd_construction(args: argparse.Namespace) -> int:
    repo = _find_repo(args.repo)
    script = repo / "examples/neural-construction/construction.py"
    out = _resolve_out(args.output, "neural-construction")
    cmd: list[str | Path] = [script, "--output", out, "--duration", str(args.duration), "--seed", str(args.seed)]
    if args.recipe is not None:
        cmd.extend(["--recipe", args.recipe])
    rc = _run_py(cmd)
    if rc == 0:
        print(f"\nOpen {out / 'viewer.html'} in your browser.")
    return rc


def _cmd_check(args: argparse.Namespace) -> int:
    repo = _find_repo(args.repo)
    steps: list[list[str | Path]] = [
        ["-m", "unittest", "discover", "-s", repo / "tests", "-v"],
        [repo / "tests/verify.py", "boundary"],
        [repo / "tests/verify.py", "minimal-example"],
        [repo / "research/check_freeze.py"],
    ]
    if args.full:
        steps.extend([
            [repo / "tests/verify.py", "neural-circuit"],
            [repo / "tests/verify.py", "neural-construction"],
            [repo / "tests/verify.py", "contract-attacks"],
        ])
        try:
            import numpy  # noqa: F401
        except ImportError:
            pass
        else:
            steps.extend([
                [repo / "tests/verify.py", "neural-regulation-motif"],
                ["-m", "unittest", "discover", "-s", repo / "studio/tests", "-v"],
            ])
    for cmd in steps:
        rc = _run_py(cmd, cwd=repo)
        if rc != 0:
            return rc
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="epigenesis", description=__doc__)
    parser.add_argument("--repo", type=Path, help="optional path to Epigenesis repository root")
    sub = parser.add_subparsers(dest="command", required=True)

    st = sub.add_parser("studio", help="launch the local Epigenesis Studio web workbench")
    st.add_argument("file", nargs="?", type=Path, help="optional FASTA, edit table, receipt JSON, or folder")
    st.add_argument("--port", type=int, help="loopback port (default 8765)")
    st.add_argument("--no-browser", action="store_true", help="start server without opening a browser tab")

    cc = sub.add_parser("circuit", help="run Variant-to-Circuit analysis, comparison, benchmark, or CX3D bridge")
    cc.add_argument("--locus", default="chr15:90936956-90937127", help="GSE188264 reference enhancer locus for demo")
    cc.add_argument("--target", choices=("Nkx2.2", "Olig2", "Both"), default="Nkx2.2", help="downstream neural-tube GRN target")
    cc.add_argument("--protein", choices=("GLI2", "GLI3"), default="GLI2")
    cc.add_argument("--signal", type=float, default=0.20, help="Gli activator input S in [0, 1]")
    cc.add_argument("--reference", help="custom reference DNA sequence or FASTA path for pairwise comparison")
    cc.add_argument("--variant", help="custom variant DNA sequence or FASTA path (SNV, MNV, or indel)")
    cc.add_argument("--with-construction", action="store_true", help="chain GRN steady states into CX3D morphogenesis and LIF controls")
    cc.add_argument("--duration", type=float, default=3.0, help="CX3D simulation duration when --with-construction is set")
    cc.add_argument("--benchmark", action="store_true", help="run multi-locus GSE188264 cryptic bifurcation benchmark")
    cc.add_argument("--output", type=Path, help="output directory (auto-generated under ./epigenesis-out/ if omitted)")

    sc = sub.add_parser("scan", help="scan all 3*L single-base mutations of a FASTA enhancer through the circuit ODE")
    sc.add_argument("fasta", type=Path, help="FASTA file containing a single enhancer sequence")
    sc.add_argument("--protein", choices=("GLI2", "GLI3"), default="GLI2")
    sc.add_argument("--target", choices=("Nkx2.2", "Olig2", "Both"), default="Nkx2.2")
    sc.add_argument("--signal", type=float, default=0.20)
    sc.add_argument("--output", type=Path, help="output directory")

    pr = sub.add_parser("predict", help="score DNA edit pairs with Huber-linear (default) or Kreimer 0-72h")
    pr.add_argument("--input", type=Path, required=True, help="JSON list of {id, reference, sequence, optional family}")
    pr.add_argument("--output", type=Path, required=True, help="new JSON file to write")
    pr.add_argument("--budget", type=int, help="optional panel size for Huber-linear or Kreimer")
    group = pr.add_mutually_exclusive_group()
    group.add_argument("--kreimer", action="store_true", help="use the Kreimer H1 neural-induction 0-72h predictor")
    group.add_argument("--activity", action="store_true", help="use archived Salomon activity-transfer model")
    group.add_argument("--attention", action="store_true", help="use archived Salomon attention mixture")
    pr.add_argument("--families", help="legacy models only: JSON {edit id: family id}")
    pr.add_argument("--device", choices=("cpu", "mps"), help="attention mixture only")

    ev = sub.add_parser("evidence", help="replay the 109 published CHE-1 measurements traced to DNA edits")
    ev.add_argument("--output", type=Path, help="output directory")

    pt = sub.add_parser("patterning", help="run the 9-bp Nkx2.2 enhancer patterning control experiment")
    pt.add_argument("--output", type=Path, help="output directory")

    cx = sub.add_parser("construction", help="compile a 480-bp developmental recipe and run the CX3D + LIF assay")
    cx.add_argument("--recipe", type=Path, help="optional 480-bp recipe DNA file")
    cx.add_argument("--duration", type=float, default=3.0)
    cx.add_argument("--seed", type=int, default=7)
    cx.add_argument("--output", type=Path, help="output directory")

    ck = sub.add_parser("check", help="run self-checks (unit tests, compiler boundary, research freeze)")
    ck.add_argument("--full", action="store_true", help="also run neural-circuit, construction, attacks, and studio checks")

    return parser


def main(argv: list[str] | None = None) -> int:
    if sys.version_info < (3, 11):
        print("Python 3.11 or newer is required.", file=sys.stderr)
        return 1
    parser = build_parser()
    args = parser.parse_args(argv)
    dispatch = {
        "studio": _cmd_studio,
        "circuit": _cmd_circuit,
        "scan": _cmd_scan,
        "predict": _cmd_predict,
        "evidence": _cmd_evidence,
        "patterning": _cmd_patterning,
        "construction": _cmd_construction,
        "check": _cmd_check,
    }
    try:
        return dispatch[args.command](args)
    except (FileNotFoundError, FileExistsError, ValueError, OSError) as err:
        print(f"epigenesis: {err}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

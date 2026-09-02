"""Bounded execution support for explicitly selected external DNA tools."""

from __future__ import annotations

from collections.abc import Mapping
import hashlib
import os
from pathlib import Path
import stat
import subprocess
import sys
import tempfile
import threading
import time
from typing import Any

from ._canonical import ContractError, canonical_bytes, loads
from .external_profile import MAX_EVIDENCE_BYTES


MAX_STDERR_BYTES = 64 * 1024
MAX_EXECUTION_SECONDS = 60
MAX_EXECUTABLE_BYTES = 256 * 1024 * 1024
MAX_REQUEST_BYTES = 2 * 1024 * 1024 + 2 * MAX_EVIDENCE_BYTES
CHUNK_BYTES = 64 * 1024


class ExternalExecutionError(ContractError):
    """An external frontend or validator could not run within its contract."""


def _fail(detail: str) -> ExternalExecutionError:
    return ExternalExecutionError(f"EXTEXEC001: {detail}")


def _copy_regular(
    source: Path,
    target: Path,
    *,
    maximum_bytes: int,
    label: str,
    executable: bool = False,
) -> dict[str, Any]:
    source_fd = target_fd = -1
    try:
        before = source.lstat()
        if (
            stat.S_ISLNK(before.st_mode)
            or not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
        ):
            raise _fail(f"{label} must be one regular single-link file")
        if before.st_size > maximum_bytes:
            raise _fail(f"{label} exceeds byte ceiling {maximum_bytes}")
        flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
        flags |= getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
        source_fd = os.open(source, flags)
        opened = os.fstat(source_fd)
        stable = ("st_dev", "st_ino", "st_nlink", "st_size", "st_mtime_ns", "st_ctime_ns")
        if (
            not stat.S_ISREG(opened.st_mode)
            or opened.st_nlink != 1
            or any(getattr(before, field) != getattr(opened, field) for field in stable)
        ):
            raise _fail(f"{label} changed while opening")
        target_fd = os.open(target, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o700 if executable else 0o600)
        hasher = hashlib.sha256()
        length = 0
        while True:
            chunk = os.read(source_fd, min(CHUNK_BYTES, maximum_bytes + 1 - length))
            if not chunk:
                break
            length += len(chunk)
            if length > maximum_bytes:
                raise _fail(f"{label} exceeds byte ceiling {maximum_bytes}")
            hasher.update(chunk)
            view = memoryview(chunk)
            while view:
                written = os.write(target_fd, view)
                if written <= 0:
                    raise _fail(f"cannot snapshot {label}: short write")
                view = view[written:]
        os.fsync(target_fd)
        finished = os.fstat(source_fd)
        after = source.lstat()
        if (
            length != opened.st_size
            or stat.S_ISLNK(after.st_mode)
            or any(getattr(opened, field) != getattr(finished, field) for field in stable)
            or any(getattr(opened, field) != getattr(after, field) for field in stable)
        ):
            raise _fail(f"{label} changed while reading")
        return {"sha256": hasher.hexdigest(), "byte_length": length}
    except ExternalExecutionError:
        raise
    except (OSError, TypeError, ValueError) as failure:
        raise _fail(f"cannot snapshot {label}: {failure}") from failure
    finally:
        for descriptor in (target_fd, source_fd):
            if descriptor >= 0:
                try:
                    os.close(descriptor)
                except OSError:
                    pass


def snapshot_inputs(
    declarations: list[dict[str, Any]],
    paths: Mapping[str, str | os.PathLike[str]],
    root: str | os.PathLike[str],
    *,
    maximum_total_bytes: int,
) -> tuple[list[dict[str, Any]], dict[str, Path]]:
    """Copy exact inputs to stable private files and return ordered references."""

    expected = [member["role"] for member in declarations]
    if set(paths) != set(expected) or len(paths) != len(expected):
        raise _fail("original input role closure does not match the profile")
    directory = Path(root)
    references: list[dict[str, Any]] = []
    snapshots: dict[str, Path] = {}
    total = 0
    for index, declaration in enumerate(declarations):
        role = declaration["role"]
        snapshot = directory / f"input-{index}"
        observed = _copy_regular(
            Path(paths[role]),
            snapshot,
            maximum_bytes=declaration["maximum_byte_length"],
            label=f"original source {role}",
        )
        total += observed["byte_length"]
        if total > maximum_total_bytes:
            raise _fail("original inputs exceed their cumulative byte ceiling")
        references.append({"role": role, **observed})
        snapshots[role] = snapshot.resolve()
    return references, snapshots


def _drain(
    stream: Any,
    limit: int,
    chunks: list[bytes],
    overflow: threading.Event,
) -> None:
    total = 0
    while True:
        chunk = stream.read(CHUNK_BYTES)
        if not chunk:
            return
        total += len(chunk)
        if total > limit:
            overflow.set()
            return
        chunks.append(chunk)


def run_command(
    command: dict[str, Any],
    executable_path: str | os.PathLike[str],
    request: dict[str, Any],
    *,
    label: str,
    maximum_stdout_bytes: int = MAX_EVIDENCE_BYTES,
) -> dict[str, Any]:
    """Run one digest-pinned command with bounded JSON input and output."""

    try:
        request_bytes = canonical_bytes(request)
    except ContractError as failure:
        raise _fail(f"{label} request is not canonical I-JSON: {failure}") from failure
    if len(request_bytes) > MAX_REQUEST_BYTES:
        raise _fail(f"{label} request exceeds {MAX_REQUEST_BYTES} bytes")

    with tempfile.TemporaryDirectory(prefix="brainc-external-command-") as temporary:
        root = Path(temporary)
        suffix = Path(executable_path).suffix if command["runtime"] == "native" else ".py"
        executable = root / f"command{suffix}"
        observed = _copy_regular(
            Path(executable_path),
            executable,
            maximum_bytes=MAX_EXECUTABLE_BYTES,
            label=f"{label} executable",
            executable=True,
        )
        if observed["sha256"] != command["executable_sha256"]:
            raise _fail(f"{label} executable SHA-256 differs from the manifest")
        request_path = root / "request.json"
        request_path.write_bytes(request_bytes)
        argv = (
            [sys.executable, "-I", str(executable)]
            if command["runtime"] == "python"
            else [str(executable)]
        )
        environment = {
            "PATH": os.defpath,
            "PYTHONHASHSEED": "0",
            "TMPDIR": str(root),
        }
        if "SystemRoot" in os.environ:
            environment["SystemRoot"] = os.environ["SystemRoot"]
        stdout_chunks: list[bytes] = []
        stderr_chunks: list[bytes] = []
        overflow = threading.Event()
        try:
            with request_path.open("rb") as request_stream:
                process = subprocess.Popen(
                    argv,
                    stdin=request_stream,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    cwd=root,
                    env=environment,
                    shell=False,
                )
                assert process.stdout is not None and process.stderr is not None
                stdout_thread = threading.Thread(
                    target=_drain,
                    args=(process.stdout, maximum_stdout_bytes, stdout_chunks, overflow),
                    daemon=True,
                )
                stderr_thread = threading.Thread(
                    target=_drain,
                    args=(process.stderr, MAX_STDERR_BYTES, stderr_chunks, overflow),
                    daemon=True,
                )
                stdout_thread.start()
                stderr_thread.start()
                try:
                    deadline = time.monotonic() + MAX_EXECUTION_SECONDS
                    while process.poll() is None and not overflow.is_set():
                        if time.monotonic() >= deadline:
                            process.kill()
                            process.wait()
                            raise _fail(
                                f"{label} exceeded {MAX_EXECUTION_SECONDS} seconds"
                            )
                        time.sleep(0.01)
                    if overflow.is_set():
                        process.kill()
                        process.wait()
                        raise _fail(f"{label} exceeded its output byte ceiling")
                    process.wait()
                finally:
                    if process.poll() is None:
                        process.kill()
                        process.wait()
                    stdout_thread.join(timeout=1)
                    stderr_thread.join(timeout=1)
                    process.stdout.close()
                    process.stderr.close()
                if overflow.is_set():
                    raise _fail(f"{label} exceeded its output byte ceiling")
        except ExternalExecutionError:
            raise
        except (OSError, subprocess.SubprocessError) as failure:
            raise _fail(f"cannot execute {label}: {failure}") from failure
        stderr = b"".join(stderr_chunks).decode("utf-8", "replace").strip()
        if process.returncode != 0:
            detail = f": {stderr}" if stderr else ""
            raise _fail(f"{label} exited with status {process.returncode}{detail}")
        if stderr:
            raise _fail(f"{label} wrote to standard error: {stderr}")
        raw = b"".join(stdout_chunks)
        if not raw:
            raise _fail(f"{label} produced no JSON result")
        try:
            return loads(raw, f"{label} result")
        except ContractError as failure:
            raise _fail(f"{label} result is invalid: {failure}") from failure


__all__ = [
    "ExternalExecutionError",
    "MAX_EXECUTION_SECONDS",
    "MAX_EXECUTABLE_BYTES",
    "MAX_REQUEST_BYTES",
    "MAX_STDERR_BYTES",
    "run_command",
    "snapshot_inputs",
]

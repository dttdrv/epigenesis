"""Atomic, no-clobber publication of immutable artifact directories."""

from __future__ import annotations

from collections.abc import Callable, Mapping
import ctypes
import errno
import os
from pathlib import Path
import secrets
import stat
import sys
import tempfile
from typing import TypeVar


RENAME_NOREPLACE = 1
_T = TypeVar("_T")
_Snapshot = tuple[int, int, int, int, int]
_DIR_FD_OPERATIONS = (os.open, os.unlink, os.stat, os.mkdir, os.rmdir)
HAS_LINUX_DIRECTORY_PUBLICATION = (
    os.name == "posix"
    and sys.platform.startswith("linux")
    and all(
        operation in getattr(os, "supports_dir_fd", set())
        for operation in _DIR_FD_OPERATIONS
    )
)


class PublicationError(OSError):
    """A publication failure with an explicit native-commit state."""

    def __init__(self, detail: str, *, committed: bool, error: int = errno.EIO):
        super().__init__(error, detail)
        self.committed = committed


def _linux_renameat2() -> ctypes._CFuncPtr:
    try:
        renameat2 = ctypes.CDLL(None, use_errno=True).renameat2
    except (AttributeError, OSError) as failure:
        raise PublicationError(
            "renameat2 is unavailable for atomic directory publication",
            committed=False,
            error=errno.ENOTSUP,
        ) from failure
    renameat2.argtypes = (
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    )
    renameat2.restype = ctypes.c_int
    return renameat2


def _call_renameat2(
    renameat2: ctypes._CFuncPtr,
    parent_descriptor: int,
    source_name: str,
    destination_name: str,
) -> None:
    ctypes.set_errno(0)
    if renameat2(
        parent_descriptor,
        os.fsencode(source_name),
        parent_descriptor,
        os.fsencode(destination_name),
        RENAME_NOREPLACE,
    ) != 0:
        error = ctypes.get_errno() or errno.EIO
        raise OSError(error, os.strerror(error), destination_name)


def _directory_flags() -> int:
    flags = os.O_RDONLY
    flags |= getattr(os, "O_CLOEXEC", 0)
    flags |= getattr(os, "O_DIRECTORY", 0)
    flags |= getattr(os, "O_NOFOLLOW", 0)
    return flags


def _snapshot(metadata: os.stat_result) -> _Snapshot:
    return (
        metadata.st_dev,
        metadata.st_ino,
        metadata.st_size,
        metadata.st_mtime_ns,
        metadata.st_ctime_ns,
    )


def _same_directory(metadata: os.stat_result, identity: tuple[int, int]) -> bool:
    return (
        stat.S_ISDIR(metadata.st_mode)
        and not stat.S_ISLNK(metadata.st_mode)
        and (metadata.st_dev, metadata.st_ino) == identity
    )


def _validate_request(
    destination: str | Path,
    entries: Mapping[str, _T],
) -> tuple[Path, tuple[tuple[str, _T], ...]]:
    target = Path(destination)
    if target.name in {"", ".", ".."}:
        raise PublicationError(
            "output must name a new child directory",
            committed=False,
            error=errno.EINVAL,
        )
    if not isinstance(entries, Mapping) or not entries:
        raise PublicationError(
            "publication entries must be a nonempty mapping",
            committed=False,
            error=errno.EINVAL,
        )
    normalized: list[tuple[str, _T]] = []
    for filename, value in entries.items():
        if (
            type(filename) is not str
            or filename in {"", ".", ".."}
            or Path(filename).name != filename
            or "/" in filename
            or "\\" in filename
            or "\x00" in filename
        ):
            raise PublicationError(
                "publication entry names must be plain child filenames",
                committed=False,
                error=errno.EINVAL,
            )
        normalized.append((filename, value))
    normalized.sort(key=lambda item: item[0].encode("utf-8"))
    return target, tuple(normalized)


def _write_all(descriptor: int, raw: bytes) -> None:
    view = memoryview(raw)
    while view:
        written = os.write(descriptor, view)
        if written <= 0:
            raise OSError(errno.EIO, "artifact write made no progress")
        view = view[written:]


def _write_linux_child(
    directory_descriptor: int,
    filename: str,
    raw: bytes,
    owned_children: dict[str, tuple[int, int]],
) -> _Snapshot:
    descriptor = -1
    try:
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
        flags |= getattr(os, "O_CLOEXEC", 0)
        flags |= getattr(os, "O_NOFOLLOW", 0)
        descriptor = os.open(
            filename,
            flags,
            0o600,
            dir_fd=directory_descriptor,
        )
        opened = os.fstat(descriptor)
        if not stat.S_ISREG(opened.st_mode):
            raise OSError(errno.EINVAL, f"staged {filename} is not a regular file")
        identity = (opened.st_dev, opened.st_ino)
        owned_children[filename] = identity
        _write_all(descriptor, raw)
        os.fsync(descriptor)
        finished = os.fstat(descriptor)
        if (
            not stat.S_ISREG(finished.st_mode)
            or (finished.st_dev, finished.st_ino) != identity
            or finished.st_size != len(raw)
        ):
            raise OSError(errno.ESTALE, f"staged {filename} changed while writing")
        os.close(descriptor)
        descriptor = -1
        published = os.stat(
            filename,
            dir_fd=directory_descriptor,
            follow_symlinks=False,
        )
        if (
            not stat.S_ISREG(published.st_mode)
            or (published.st_dev, published.st_ino) != identity
            or published.st_size != len(raw)
        ):
            raise OSError(errno.ESTALE, f"staged {filename} changed after writing")
        return _snapshot(published)
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _require_linux_children(
    directory_descriptor: int,
    snapshots: Mapping[str, _Snapshot],
) -> None:
    if set(os.listdir(directory_descriptor)) != set(snapshots):
        raise OSError(errno.ESTALE, "staging directory contents changed")
    for filename, expected in snapshots.items():
        current = os.stat(
            filename,
            dir_fd=directory_descriptor,
            follow_symlinks=False,
        )
        if not stat.S_ISREG(current.st_mode) or _snapshot(current) != expected:
            raise OSError(errno.ESTALE, f"staged {filename} changed")


def _publish_linux(
    destination: Path,
    entries: tuple[tuple[str, _T], ...],
    encode: Callable[[_T], bytes],
    renameat2: ctypes._CFuncPtr,
) -> dict[str, Path]:
    parent_descriptor = -1
    staging_descriptor = -1
    staging_name: str | None = None
    staging_identity: tuple[int, int] | None = None
    parent_identity: tuple[int, int] | None = None
    owned_children: dict[str, tuple[int, int]] = {}
    snapshots: dict[str, _Snapshot] = {}
    committed = False
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        parent_inspected = destination.parent.lstat()
        if stat.S_ISLNK(parent_inspected.st_mode) or not stat.S_ISDIR(
            parent_inspected.st_mode
        ):
            raise OSError(errno.ENOTDIR, "output parent must be a non-linked directory")
        parent_descriptor = os.open(destination.parent, _directory_flags())
        parent_opened = os.fstat(parent_descriptor)
        parent_identity = (parent_opened.st_dev, parent_opened.st_ino)
        if not _same_directory(parent_opened, parent_identity) or parent_identity != (
            parent_inspected.st_dev,
            parent_inspected.st_ino,
        ):
            raise OSError(errno.ESTALE, "output parent changed while opening")
        try:
            os.stat(
                destination.name,
                dir_fd=parent_descriptor,
                follow_symlinks=False,
            )
        except FileNotFoundError:
            pass
        else:
            raise FileExistsError(errno.EEXIST, "output already exists", destination.name)

        for _ in range(128):
            candidate = f".{destination.name}.{secrets.token_hex(12)}"
            try:
                os.mkdir(candidate, 0o700, dir_fd=parent_descriptor)
            except FileExistsError:
                continue
            staging_name = candidate
            break
        if staging_name is None:
            raise OSError(errno.EEXIST, "cannot allocate a unique staging directory")
        staging_entry = os.stat(
            staging_name,
            dir_fd=parent_descriptor,
            follow_symlinks=False,
        )
        staging_identity = (staging_entry.st_dev, staging_entry.st_ino)
        if not _same_directory(staging_entry, staging_identity):
            raise OSError(errno.ESTALE, "staging entry is not a safe directory")
        staging_descriptor = os.open(
            staging_name,
            _directory_flags(),
            dir_fd=parent_descriptor,
        )
        staging_opened = os.fstat(staging_descriptor)
        if not _same_directory(staging_opened, staging_identity):
            raise OSError(errno.ESTALE, "staging directory changed while opening")

        for filename, value in entries:
            raw = encode(value)
            if type(raw) is not bytes:
                raise TypeError("publication encoder must return bytes")
            snapshots[filename] = _write_linux_child(
                staging_descriptor,
                filename,
                raw,
                owned_children,
            )
        _require_linux_children(staging_descriptor, snapshots)
        os.fsync(staging_descriptor)
        current_staging = os.stat(
            staging_name,
            dir_fd=parent_descriptor,
            follow_symlinks=False,
        )
        if not _same_directory(current_staging, staging_identity):
            raise OSError(errno.ESTALE, "staging directory changed before publication")
        _require_linux_children(staging_descriptor, snapshots)

        _call_renameat2(
            renameat2,
            parent_descriptor,
            staging_name,
            destination.name,
        )
        committed = True

        published = os.stat(
            destination.name,
            dir_fd=parent_descriptor,
            follow_symlinks=False,
        )
        if not _same_directory(published, staging_identity):
            raise OSError(errno.ESTALE, "published directory changed identity")
        _require_linux_children(staging_descriptor, snapshots)
        parent_finished = destination.parent.lstat()
        if not _same_directory(parent_finished, parent_identity):
            raise OSError(errno.ESTALE, "output parent changed during publication")
        os.fsync(parent_descriptor)
        final_entry = os.stat(
            destination.name,
            dir_fd=parent_descriptor,
            follow_symlinks=False,
        )
        if not _same_directory(final_entry, staging_identity):
            raise OSError(errno.ESTALE, "published directory changed after synchronization")
        _require_linux_children(staging_descriptor, snapshots)
        return {filename: destination / filename for filename, _ in entries}
    except PublicationError:
        raise
    except OSError as failure:
        detail = (
            "output path appeared during staging"
            if not committed and failure.errno in {errno.EEXIST, errno.ENOTEMPTY}
            else f"cannot publish directory: {failure}"
        )
        raise PublicationError(
            detail,
            committed=committed,
            error=failure.errno or errno.EIO,
        ) from failure
    finally:
        if not committed and staging_descriptor >= 0:
            for filename, identity in reversed(tuple(owned_children.items())):
                try:
                    current = os.stat(
                        filename,
                        dir_fd=staging_descriptor,
                        follow_symlinks=False,
                    )
                    if (current.st_dev, current.st_ino) == identity:
                        os.unlink(filename, dir_fd=staging_descriptor)
                except OSError:
                    pass
        if staging_descriptor >= 0:
            os.close(staging_descriptor)
        if (
            not committed
            and parent_descriptor >= 0
            and staging_name is not None
            and staging_identity is not None
        ):
            try:
                current = os.stat(
                    staging_name,
                    dir_fd=parent_descriptor,
                    follow_symlinks=False,
                )
                if _same_directory(current, staging_identity):
                    os.rmdir(staging_name, dir_fd=parent_descriptor)
            except OSError:
                pass
        if parent_descriptor >= 0:
            os.close(parent_descriptor)


def _write_windows_child(
    path: Path,
    raw: bytes,
    owned_children: dict[str, tuple[int, int]],
) -> _Snapshot:
    descriptor = -1
    try:
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
        flags |= getattr(os, "O_BINARY", 0)
        flags |= getattr(os, "O_NOINHERIT", 0)
        descriptor = os.open(path, flags, 0o600)
        opened = os.fstat(descriptor)
        if not stat.S_ISREG(opened.st_mode):
            raise OSError(errno.EINVAL, "staged output is not a regular file")
        identity = (opened.st_dev, opened.st_ino)
        owned_children[path.name] = identity
        _write_all(descriptor, raw)
        os.fsync(descriptor)
        finished = os.fstat(descriptor)
        if (
            not stat.S_ISREG(finished.st_mode)
            or (finished.st_dev, finished.st_ino) != identity
            or finished.st_size != len(raw)
        ):
            raise OSError(errno.ESTALE, "staged output changed while writing")
        os.close(descriptor)
        descriptor = -1
        visible = path.lstat()
        if (
            not stat.S_ISREG(visible.st_mode)
            or (visible.st_dev, visible.st_ino) != identity
            or visible.st_size != len(raw)
        ):
            raise OSError(errno.ESTALE, "staged output changed after writing")
        return _snapshot(visible)
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _require_windows_children(
    directory: Path,
    identity: tuple[int, int],
    snapshots: Mapping[str, _Snapshot],
) -> None:
    current_directory = directory.lstat()
    if not _same_directory(current_directory, identity):
        raise OSError(errno.ESTALE, "staging directory changed")
    if {entry.name for entry in directory.iterdir()} != set(snapshots):
        raise OSError(errno.ESTALE, "staging directory contents changed")
    for filename, expected in snapshots.items():
        current = (directory / filename).lstat()
        if not stat.S_ISREG(current.st_mode) or _snapshot(current) != expected:
            raise OSError(errno.ESTALE, f"staged {filename} changed")


def _publish_windows(
    destination: Path,
    entries: tuple[tuple[str, _T], ...],
    encode: Callable[[_T], bytes],
) -> dict[str, Path]:
    staging: Path | None = None
    staging_identity: tuple[int, int] | None = None
    parent_identity: tuple[int, int] | None = None
    owned_children: dict[str, tuple[int, int]] = {}
    snapshots: dict[str, _Snapshot] = {}
    committed = False
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        parent = destination.parent.lstat()
        parent_identity = (parent.st_dev, parent.st_ino)
        if not _same_directory(parent, parent_identity):
            raise OSError(errno.ENOTDIR, "output parent must be a non-linked directory")
        try:
            destination.lstat()
        except FileNotFoundError:
            pass
        else:
            raise FileExistsError(errno.EEXIST, "output already exists", str(destination))
        staging = Path(
            tempfile.mkdtemp(prefix=f".{destination.name}.", dir=destination.parent)
        )
        metadata = staging.lstat()
        staging_identity = (metadata.st_dev, metadata.st_ino)
        if not _same_directory(metadata, staging_identity):
            raise OSError(errno.ESTALE, "staging entry is not a safe directory")
        for filename, value in entries:
            raw = encode(value)
            if type(raw) is not bytes:
                raise TypeError("publication encoder must return bytes")
            snapshots[filename] = _write_windows_child(
                staging / filename,
                raw,
                owned_children,
            )
        _require_windows_children(staging, staging_identity, snapshots)
        os.rename(staging, destination)
        committed = True
        _require_windows_children(destination, staging_identity, snapshots)
        parent_finished = destination.parent.lstat()
        if not _same_directory(parent_finished, parent_identity):
            raise OSError(errno.ESTALE, "output parent changed during publication")
        return {filename: destination / filename for filename, _ in entries}
    except OSError as failure:
        detail = (
            "output path appeared during staging"
            if not committed and failure.errno in {errno.EEXIST, errno.ENOTEMPTY}
            else f"cannot publish directory: {failure}"
        )
        raise PublicationError(
            detail,
            committed=committed,
            error=failure.errno or errno.EIO,
        ) from failure
    finally:
        if not committed and staging is not None and staging_identity is not None:
            for filename, identity in reversed(tuple(owned_children.items())):
                try:
                    current = (staging / filename).lstat()
                    if (current.st_dev, current.st_ino) == identity:
                        (staging / filename).unlink()
                except OSError:
                    pass
            try:
                current = staging.lstat()
                if _same_directory(current, staging_identity):
                    staging.rmdir()
            except OSError:
                pass


def publish_directory(
    destination: str | Path,
    entries: Mapping[str, _T],
    encode: Callable[[_T], bytes],
) -> dict[str, Path]:
    """Publish exact encoded children through one absent-to-complete commit."""

    target, normalized = _validate_request(destination, entries)
    if not callable(encode):
        raise PublicationError(
            "publication encoder must be callable",
            committed=False,
            error=errno.EINVAL,
        )
    if os.name == "nt":
        return _publish_windows(target, normalized, encode)
    if not HAS_LINUX_DIRECTORY_PUBLICATION:
        raise PublicationError(
            "atomic directory publication is unsupported on this host",
            committed=False,
            error=errno.ENOTSUP,
        )
    renameat2 = _linux_renameat2()
    return _publish_linux(target, normalized, encode, renameat2)


def rename_directory_noreplace(
    source: Path,
    destination: Path,
    *,
    parent_descriptor: int = -1,
) -> None:
    """Low-level no-clobber rename retained for focused native probes."""

    source = Path(source)
    destination = Path(destination)
    if source.parent != destination.parent:
        raise OSError(errno.EXDEV, "publication paths must share one parent")
    if os.name == "nt":
        os.rename(source, destination)
        return
    if not HAS_LINUX_DIRECTORY_PUBLICATION or parent_descriptor < 0:
        raise OSError(errno.ENOTSUP, "descriptor-relative rename is unsupported")
    _call_renameat2(
        _linux_renameat2(),
        parent_descriptor,
        source.name,
        destination.name,
    )


__all__ = [
    "HAS_LINUX_DIRECTORY_PUBLICATION",
    "PublicationError",
    "publish_directory",
    "rename_directory_noreplace",
]

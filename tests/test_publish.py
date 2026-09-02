from __future__ import annotations

import errno
import os
from pathlib import Path
import stat
import tempfile
import unittest
from unittest import mock

from brainc import _publish as publication


class NativeDirectoryPublicationTests(unittest.TestCase):
    def test_native_commit_publishes_without_replacing_a_destination(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            destination = root / "published"
            publication.publish_directory(
                destination,
                {"artifact.json": b"{}\n"},
                lambda value: value,
            )
            self.assertEqual((destination / "artifact.json").read_bytes(), b"{}\n")

            with self.assertRaises(publication.PublicationError) as caught:
                publication.publish_directory(
                    destination,
                    {"artifact.json": b'{"changed":true}\n'},
                    lambda value: value,
                )
            self.assertFalse(caught.exception.committed)
            self.assertEqual((destination / "artifact.json").read_bytes(), b"{}\n")


@unittest.skipUnless(
    publication.HAS_POSIX_DIRECTORY_PUBLICATION,
    "requires POSIX descriptor-relative publication",
)
class DirectoryPublicationTests(unittest.TestCase):
    @staticmethod
    def _commit_hook():
        if publication.HAS_LINUX_DIRECTORY_PUBLICATION:
            return "_call_renameat2", publication._call_renameat2
        return "_call_renameatx_np", publication._call_renameatx_np

    def _entries(self) -> dict[str, bytes]:
        return {"a.json": b'{"a":1}\n', "b.json": b'{"b":2}\n'}

    @staticmethod
    def _encode(value: bytes) -> bytes:
        return value

    def test_commit_is_absent_to_complete(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            destination = root / "published"
            observed: list[tuple[bool, set[str]]] = []
            hook, native = self._commit_hook()

            def inspect_commit(renameat2, parent_fd, source_name, destination_name):
                observed.append((destination.exists(), set()))
                native(renameat2, parent_fd, source_name, destination_name)
                observed.append(
                    (destination.is_dir(), {entry.name for entry in destination.iterdir()})
                )

            with mock.patch.object(
                publication,
                hook,
                side_effect=inspect_commit,
            ):
                paths = publication.publish_directory(
                    destination, self._entries(), self._encode
                )
            self.assertEqual(
                observed,
                [(False, set()), (True, {"a.json", "b.json"})],
            )
            self.assertEqual(paths["a.json"].read_bytes(), b'{"a":1}\n')

    def test_existing_destinations_are_never_replaced(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            victim = root / "victim"
            victim.mkdir()
            cases = {
                "file": lambda path: path.write_text("preserve", encoding="utf-8"),
                "empty": lambda path: path.mkdir(),
                "nonempty": lambda path: (
                    path.mkdir(),
                    (path / "SENTINEL").write_text("preserve", encoding="utf-8"),
                ),
                "symlink": lambda path: path.symlink_to(victim, target_is_directory=True),
                "broken": lambda path: path.symlink_to(root / "missing"),
            }
            for name, create in cases.items():
                with self.subTest(name=name):
                    destination = root / name
                    create(destination)
                    before = destination.lstat()
                    with self.assertRaises(publication.PublicationError) as caught:
                        publication.publish_directory(
                            destination, self._entries(), self._encode
                        )
                    self.assertFalse(caught.exception.committed)
                    after = destination.lstat()
                    self.assertEqual(
                        (after.st_dev, after.st_ino),
                        (before.st_dev, before.st_ino),
                    )

    def test_destination_appearing_at_commit_is_preserved(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            destination = root / "appeared"
            hook, native = self._commit_hook()

            def race(renameat2, parent_fd, source_name, destination_name):
                destination.mkdir()
                (destination / "SENTINEL").write_text("preserve", encoding="utf-8")
                native(renameat2, parent_fd, source_name, destination_name)

            with mock.patch.object(
                publication,
                hook,
                side_effect=race,
            ), self.assertRaises(publication.PublicationError) as caught:
                publication.publish_directory(
                    destination, self._entries(), self._encode
                )
            self.assertFalse(caught.exception.committed)
            self.assertEqual(
                (destination / "SENTINEL").read_text(encoding="utf-8"),
                "preserve",
            )
            self.assertEqual(list(root.glob(".appeared.*")), [])

    def test_encode_failure_cleans_owned_staging_and_retry_succeeds(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            destination = root / "retry"
            calls = 0

            def fail_second(value: bytes) -> bytes:
                nonlocal calls
                calls += 1
                if calls == 2:
                    raise ValueError("injected encode failure")
                return value

            with self.assertRaisesRegex(ValueError, "injected encode failure"):
                publication.publish_directory(
                    destination, self._entries(), fail_second
                )
            self.assertFalse(destination.exists())
            self.assertEqual(list(root.glob(".retry.*")), [])
            publication.publish_directory(destination, self._entries(), self._encode)
            self.assertTrue((destination / "b.json").is_file())

    def test_file_and_staging_fsync_failures_are_precommit_and_retryable(self) -> None:
        for kind in ("file", "directory"):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                destination = root / "retry"
                native = os.fsync
                failed = False

                def inject(descriptor: int) -> None:
                    nonlocal failed
                    metadata = os.fstat(descriptor)
                    selected = (
                        stat.S_ISREG(metadata.st_mode)
                        if kind == "file"
                        else stat.S_ISDIR(metadata.st_mode)
                    )
                    if selected and not failed:
                        failed = True
                        raise OSError(errno.EINVAL, "injected fsync failure")
                    native(descriptor)

                with mock.patch.object(publication.os, "fsync", side_effect=inject):
                    with self.assertRaises(publication.PublicationError) as caught:
                        publication.publish_directory(
                            destination, self._entries(), self._encode
                        )
                self.assertTrue(failed)
                self.assertFalse(caught.exception.committed)
                self.assertFalse(destination.exists())
                self.assertEqual(list(root.glob(".retry.*")), [])
                publication.publish_directory(
                    destination, self._entries(), self._encode
                )

    def test_postcommit_child_swap_is_reported_without_cleanup(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            destination = root / "published"
            displaced = root / "displaced.json"
            hook, native = self._commit_hook()

            def swap_after_commit(renameat2, parent_fd, source_name, destination_name):
                native(renameat2, parent_fd, source_name, destination_name)
                (destination / "a.json").rename(displaced)
                (destination / "a.json").write_bytes(b'{"attacker":true}\n')

            with mock.patch.object(
                publication,
                hook,
                side_effect=swap_after_commit,
            ), self.assertRaises(publication.PublicationError) as caught:
                publication.publish_directory(
                    destination, self._entries(), self._encode
                )
            self.assertTrue(caught.exception.committed)
            self.assertEqual(
                (destination / "a.json").read_bytes(), b'{"attacker":true}\n'
            )
            self.assertEqual(displaced.read_bytes(), b'{"a":1}\n')
            self.assertTrue((destination / "b.json").is_file())

    def test_parent_fsync_failure_retains_complete_committed_tree(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            destination = root / "published"
            native = os.fsync
            directory_calls = 0

            def fail_parent_sync(descriptor: int) -> None:
                nonlocal directory_calls
                if stat.S_ISDIR(os.fstat(descriptor).st_mode):
                    directory_calls += 1
                    if directory_calls == 2:
                        raise OSError(errno.EIO, "injected parent fsync failure")
                native(descriptor)

            with mock.patch.object(
                publication.os,
                "fsync",
                side_effect=fail_parent_sync,
            ), self.assertRaises(publication.PublicationError) as caught:
                publication.publish_directory(
                    destination, self._entries(), self._encode
                )
            self.assertTrue(caught.exception.committed)
            self.assertEqual(
                {entry.name for entry in destination.iterdir()},
                {"a.json", "b.json"},
            )
            self.assertEqual((destination / "a.json").read_bytes(), b'{"a":1}\n')

    def test_missing_native_support_fails_before_filesystem_mutation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            parent = Path(temporary) / "not-created"
            destination = parent / "output"
            loader = (
                "_linux_renameat2"
                if publication.HAS_LINUX_DIRECTORY_PUBLICATION
                else "_darwin_renameatx_np"
            )
            unavailable = publication.PublicationError(
                "missing native rename",
                committed=False,
                error=errno.ENOTSUP,
            )
            with mock.patch.object(
                publication,
                loader,
                side_effect=unavailable,
            ), self.assertRaises(publication.PublicationError) as caught:
                publication.publish_directory(
                    destination, self._entries(), self._encode
                )
            self.assertFalse(caught.exception.committed)
            self.assertFalse(parent.exists())


if __name__ == "__main__":
    unittest.main()

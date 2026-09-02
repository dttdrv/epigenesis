from __future__ import annotations

from pathlib import Path
import tempfile
import unittest
from unittest import mock

from brainc import _io


class RegularFileTests(unittest.TestCase):
    def test_read_requests_binary_mode_when_the_platform_defines_it(self) -> None:
        fake_binary = 1 << 29
        observed: list[int] = []
        native_open = _io.os.open
        native_binary = getattr(_io.os, "O_BINARY", 0)

        def inspect_open(path, flags, *args, **kwargs):
            observed.append(flags)
            return native_open(path, flags & ~fake_binary, *args, **kwargs)

        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "input.txt"
            path.write_bytes(b"a\r\nb\r\n")
            with mock.patch.object(
                _io.os, "O_BINARY", native_binary | fake_binary, create=True
            ), mock.patch.object(_io.os, "open", side_effect=inspect_open):
                self.assertEqual(_io.read_regular_file(path), b"a\r\nb\r\n")

        self.assertTrue(observed)
        self.assertTrue(observed[0] & fake_binary)


if __name__ == "__main__":
    unittest.main()

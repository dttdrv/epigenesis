from __future__ import annotations

import hashlib
import math
import unittest
from unittest.mock import patch

from brainc._canonical import ContractError, artifact_digest, canonical_bytes, digest


class StreamingCanonicalDigestTests(unittest.TestCase):
    def test_streaming_digest_matches_canonical_bytes_without_calling_them(self) -> None:
        value = {
            "records": [
                {
                    "id": f"record-{index}",
                    "bases": index + 1,
                    "flags": [None, True, False],
                    "text": "ACGTRYSWKMBDHVN" * 8,
                }
                for index in range(10_000)
            ],
            "numbers": [0.0, -0.0, 1e-7, 1e-6, 1e20, 1e21],
        }
        expected = hashlib.sha256(canonical_bytes(value)).hexdigest()
        with patch("brainc._canonical.canonical_bytes", side_effect=AssertionError):
            self.assertEqual(digest(value), expected)

    def test_artifact_digest_uses_the_same_streaming_path(self) -> None:
        core = {"format": "example", "version": 1, "items": list(range(20_000))}
        sealed = {**core, "artifact_sha256": digest(core)}
        with patch("brainc._canonical.canonical_bytes", side_effect=AssertionError):
            self.assertEqual(artifact_digest(sealed), sealed["artifact_sha256"])

    def test_utf16_key_order_and_error_contract_are_unchanged(self) -> None:
        value = {"\uffff": 2, "\U0001f600": 1}
        self.assertEqual(canonical_bytes(value), '{"😀":1,"￿":2}'.encode())
        self.assertEqual(digest(value), hashlib.sha256(canonical_bytes(value)).hexdigest())

        invalid = [2**53, math.nan, "\ud800", (1, 2), {1: "non-string-key"}]
        for value in invalid:
            with self.subTest(value=repr(value)):
                with self.assertRaises(ContractError):
                    canonical_bytes(value)
                with self.assertRaises(ContractError):
                    digest(value)


if __name__ == "__main__":
    unittest.main()

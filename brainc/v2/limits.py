"""Finite resource ceilings for the tensor/development ABI.

These are compiler admission limits, not promises that every runtime can host
an artifact this large.  A runtime is expected to impose equal or tighter
limits before allocation.
"""

from __future__ import annotations


MAX_JSON_BYTES = 16 * 1024 * 1024
MAX_SOURCE_JSON_BYTES = 64 * 1024 * 1024
MAX_DECOMPRESSED_BYTES = 64 * 1024 * 1024
MAX_JSON_DEPTH = 64
MAX_JSON_MEMBERS = 1_000_000
MAX_STRING_BYTES = 1 * 1024 * 1024

MAX_TENSOR_RANK = 16
MAX_TENSOR_BYTES = 64 * 1024 * 1024
MAX_TOTAL_TENSOR_BYTES = 256 * 1024 * 1024
MAX_TENSORS = 100_000
MAX_OPERATIONS = 100_000
MAX_SOURCE_RECORDS = 100_000
MAX_UNITS = 10_000_000
MAX_EDGES = 10_000_000
MAX_ATTACHMENTS = 100_000
MAX_PORTS = 100_000
MAX_SCHEMAS = 100_000
MAX_RULES = 100_000
MAX_FIELDS_PER_SCHEMA = 100_000
MAX_IDENTIFIER_BYTES = 256

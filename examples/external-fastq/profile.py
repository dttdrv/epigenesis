#!/usr/bin/env python3
"""Emit the sealed profile for the adjacent FASTQ tools."""

import hashlib
import json
from pathlib import Path

from brainc.external_profile import (
    EXECUTABLE_ABI,
    FRONTEND_PROTOCOL,
    RECORD_CATALOG_SCHEMA,
    VALIDATOR_PROTOCOL,
    seal_profile_manifest,
)


ROOT = Path(__file__).resolve().parent


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def command(name, path):
    return {
        "id": f"org.epigenesis.fastq-{name}",
        "distribution": f"org.epigenesis.fastq-{name}-example",
        "version": "1.0.0",
        "runtime": "python",
        "distribution_sha256": sha(path.read_bytes()),
        "executable_sha256": sha(path.read_bytes()),
    }


profile_ir = {
    "abi": EXECUTABLE_ABI,
    "profile": {"id": "org.epigenesis.fastq-sanger", "version": 1},
    "grammar": {
        "id": "ncbi.sra.fastq-sanger",
        "version": "epigenesis-profile-1/phred33",
        "authority": {
            "uri": "https://github.com/dttdrv/epigenesis/blob/main/examples/external-fastq/GRAMMAR.md",
            "sha256": sha((ROOT / "GRAMMAR.md").read_bytes()),
        },
    },
    "inputs": [
        {
            "role": "reads",
            "media_type": "application/vnd.ncbi.fastq",
            "wrapper": "identity",
            "maximum_byte_length": 16 * 1024 * 1024,
        }
    ],
    "native_artifacts": [
        {
            "role": "fastq",
            "media_type": "application/json",
            "format": "org.epigenesis.fastq-ir",
            "version": 1,
            "schema_sha256": sha(b"org.epigenesis.fastq-ir.schema/v1"),
            "ir_digest_field": "fastq_ir_sha256",
            "maximum_byte_length": 16 * 1024 * 1024,
        }
    ],
    "catalog": {"schema": RECORD_CATALOG_SCHEMA},
    "limits": {
        "maximum_total_input_bytes": 16 * 1024 * 1024,
        "maximum_total_native_bytes": 16 * 1024 * 1024,
        "maximum_records": 100_000,
        "maximum_total_bases": 100_000_000,
        "maximum_record_id_bytes": 256,
    },
    "frontend": {
        "protocol": FRONTEND_PROTOCOL,
        "command": command("frontend", ROOT / "frontend.py"),
    },
    "validator": {
        "protocol": VALIDATOR_PROTOCOL,
        "command": command("validator", ROOT / "validator.py"),
    },
}

print(
    json.dumps(
        seal_profile_manifest(profile_ir, version=2),
        indent=2,
        sort_keys=True,
        ensure_ascii=False,
    )
)

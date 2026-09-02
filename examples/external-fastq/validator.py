#!/usr/bin/env python3
"""Separate FASTQ validator for the Epigenesis executable ABI."""

import base64
import hashlib
import json
import re
import sys


IUPAC = frozenset(b"ACGTRYSWKMBDHVN")
ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9._:/+@-]{0,255}\Z")


def canonical(value):
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"), sort_keys=True).encode()


def parse_fastq(raw):
    normalized = raw.replace(b"\r\n", b"\n")
    if b"\r" in normalized:
        raise ValueError("FASTQ contains a bare carriage return")
    lines = normalized.split(b"\n")
    if lines and lines[-1] == b"":
        lines.pop()
    if not lines or len(lines) % 4:
        raise ValueError("FASTQ must contain complete four-line records")
    records = []
    native_records = []
    ids = set()
    for ordinal in range(len(lines) // 4):
        header, sequence, separator, quality = lines[ordinal * 4 : ordinal * 4 + 4]
        if not header.startswith(b"@") or separator != b"+":
            raise ValueError("FASTQ header or separator is invalid")
        if any(byte < 32 or byte > 126 for byte in header[1:]):
            raise ValueError("FASTQ header must be printable ASCII")
        header_text = header[1:].decode("ascii")
        record_id = header_text.split(None, 1)[0]
        sequence = sequence.upper()
        if not record_id or ID.fullmatch(record_id) is None or record_id in ids:
            raise ValueError("FASTQ record id is invalid or duplicated")
        ids.add(record_id)
        if not sequence or any(base not in IUPAC for base in sequence):
            raise ValueError("FASTQ sequence contains non-IUPAC DNA")
        if len(sequence) != len(quality) or any(byte < 33 or byte > 126 for byte in quality):
            raise ValueError("FASTQ quality is invalid")
        records.append(
            {
                "ordinal": ordinal,
                "input_role": "reads",
                "record_id": record_id,
                "bases": len(sequence),
                "sequence_sha256": hashlib.sha256(sequence).hexdigest(),
                "refget_id": "SQ."
                + base64.urlsafe_b64encode(hashlib.sha512(sequence).digest()[:24]).decode(),
            }
        )
        native_records.append(
            {
                "ordinal": ordinal,
                "header": header_text,
                "record_id": record_id,
                "sequence": sequence.decode("ascii"),
                "quality": quality.decode("ascii"),
            }
        )
    return records, native_records


def main():
    request = json.load(sys.stdin)
    if set(request) != {
        "format",
        "version",
        "profile_manifest",
        "source_descriptor",
        "inputs",
        "native_artifacts",
    }:
        raise ValueError("request keys are invalid")
    if request["format"] != "brainc.external-validator-execution-request" or request["version"] != 1:
        raise ValueError("request identity is invalid")
    if len(request["inputs"]) != 1 or request["inputs"][0]["role"] != "reads":
        raise ValueError("FASTQ profile requires exactly the reads role")
    with open(request["inputs"][0]["path"], "rb") as stream:
        raw = stream.read()
    records, native_records = parse_fastq(raw)
    fastq_ir = {"quality_encoding": "phred33", "records": native_records}
    expected_native = {
        "format": "org.epigenesis.fastq-ir",
        "version": 1,
        "fastq_ir": fastq_ir,
        "fastq_ir_sha256": hashlib.sha256(canonical(fastq_ir)).hexdigest(),
    }
    if request["native_artifacts"] != {"fastq": expected_native}:
        raise ValueError("native FASTQ IR does not match independent replay")
    manifest = request["profile_manifest"]
    declaration = manifest["profile_ir"]["native_artifacts"][0]
    native_raw = canonical(expected_native)
    native_reference = {
        "role": "fastq",
        "format": declaration["format"],
        "version": declaration["version"],
        "schema_sha256": declaration["schema_sha256"],
        "sha256": hashlib.sha256(native_raw).hexdigest(),
        "byte_length": len(native_raw),
        "ir_sha256": expected_native[declaration["ir_digest_field"]],
    }
    input_reference = {
        "role": "reads",
        "sha256": hashlib.sha256(raw).hexdigest(),
        "byte_length": len(raw),
    }
    result = {
        "format": "brainc.external-validator-execution-result",
        "version": 1,
        "profile_manifest_sha256": manifest["artifact_sha256"],
        "source_descriptor_sha256": request["source_descriptor"]["artifact_sha256"],
        "input_references": [input_reference],
        "native_artifact_references": [native_reference],
        "records": records,
        "valid": True,
    }
    sys.stdout.buffer.write(canonical(result))


if __name__ == "__main__":
    try:
        main()
    except Exception as failure:
        print(f"FASTQ validator: {failure}", file=sys.stderr)
        raise SystemExit(2)

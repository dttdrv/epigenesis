#!/usr/bin/env python3
"""Conforming FASTQ frontend for the Epigenesis executable ABI."""

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
    record_ids = set()
    for ordinal in range(len(lines) // 4):
        header, sequence, separator, quality = lines[ordinal * 4 : ordinal * 4 + 4]
        if not header.startswith(b"@") or separator != b"+":
            raise ValueError("FASTQ header or separator is invalid")
        if any(byte < 32 or byte > 126 for byte in header[1:]):
            raise ValueError("FASTQ header must be printable ASCII")
        try:
            record_id = header[1:].split(None, 1)[0].decode("ascii")
            header_text = header[1:].decode("ascii")
        except UnicodeDecodeError as failure:
            raise ValueError("FASTQ record id must be ASCII") from failure
        sequence = sequence.upper()
        if not record_id or ID.fullmatch(record_id) is None:
            raise ValueError("FASTQ record id is not portable")
        if any(base not in IUPAC for base in sequence):
            raise ValueError("FASTQ sequence contains non-IUPAC DNA")
        if not sequence or len(sequence) != len(quality):
            raise ValueError("FASTQ sequence and quality lengths differ")
        if any(byte < 33 or byte > 126 for byte in quality):
            raise ValueError("FASTQ quality is outside Phred+33 bytes")
        if record_id in record_ids:
            raise ValueError("FASTQ record ids must be unique")
        record_ids.add(record_id)
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
    if set(request) != {"format", "version", "profile_manifest", "inputs"}:
        raise ValueError("request keys are invalid")
    if request["format"] != "brainc.external-frontend-execution-request" or request["version"] != 1:
        raise ValueError("request identity is invalid")
    if len(request["inputs"]) != 1 or request["inputs"][0]["role"] != "reads":
        raise ValueError("FASTQ profile requires exactly the reads role")
    with open(request["inputs"][0]["path"], "rb") as stream:
        records, native_records = parse_fastq(stream.read())
    fastq_ir = {"quality_encoding": "phred33", "records": native_records}
    native = {
        "format": "org.epigenesis.fastq-ir",
        "version": 1,
        "fastq_ir": fastq_ir,
        "fastq_ir_sha256": hashlib.sha256(canonical(fastq_ir)).hexdigest(),
    }
    result = {
        "format": "brainc.external-frontend-execution-result",
        "version": 1,
        "profile_manifest_sha256": request["profile_manifest"]["artifact_sha256"],
        "native_artifacts": {"fastq": native},
        "records": records,
    }
    sys.stdout.buffer.write(canonical(result))


if __name__ == "__main__":
    try:
        main()
    except Exception as failure:
        print(f"FASTQ frontend: {failure}", file=sys.stderr)
        raise SystemExit(2)

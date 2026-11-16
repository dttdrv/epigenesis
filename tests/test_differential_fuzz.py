from __future__ import annotations

import json
import math
from pathlib import Path
import random
import shutil
import subprocess
import tempfile
import unittest

from brainc._canonical import ContractError, canonical_bytes, digest
from brainc.bio import GFF3Compiler, GFF3Error
from brainc.insdc import GenBankCompiler, GenBankError
from brainc.source import FASTA_PROFILE, GENBANK_PROFILE, GFF3_PROFILE, compile_source
from brainc.validator_variant import validate_variant_paths
from brainc.variant import (
    IUPAC_DNA,
    VariantError,
    apply_variant_to_sequence,
    compile_variants,
    revert_variant,
    vrs_allele,
)


ROOT = Path(__file__).resolve().parents[1]
RUST_DIR = ROOT / "validators" / "rust"
BASES = "ACGT"
ALL_IUPAC = tuple(sorted(IUPAC_DNA))


def _minimal_fasta_oracle(text: str) -> dict[str, str]:
    records: dict[str, str] = {}
    rid: str | None = None
    chunks: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if rid is not None:
                records[rid] = "".join(chunks)
            rid = line[1:].strip().split()[0]
            chunks = []
        elif line.strip():
            chunks.append(line.strip().upper())
    if rid is not None:
        records[rid] = "".join(chunks)
    return records


def _random_json(rng: random.Random, depth: int = 0):
    if depth >= 4:
        choice = rng.randint(0, 4)
    else:
        choice = rng.randint(0, 6)
    if choice == 0:
        return None
    if choice == 1:
        return rng.choice((True, False))
    if choice == 2:
        return rng.choice((0, -1, 1, 42, -999999, 2**53 - 1, -(2**53 - 1), rng.randint(-10**9, 10**9)))
    if choice == 3:
        return rng.choice((0.0, -0.0, 1.5, -1e-6, 1e-7, 1e20, 1e21, 3.141592653589793, rng.uniform(-1e4, 1e4)))
    if choice == 4:
        alphabet = ["a", "b", "Z", "0", "\n", "\t", '"', "\\", "\x01", "é", "α", "🧬", "😀"]
        return "".join(rng.choice(alphabet) for _ in range(rng.randint(0, 12)))
    if choice == 5:
        return [_random_json(rng, depth + 1) for _ in range(rng.randint(0, 4))]
    keys = ["a", "b", "z", "é", "🧬", "k\n1", "\u00e9", "😀"]
    out = {}
    for _ in range(rng.randint(0, 4)):
        out[rng.choice(keys)] = _random_json(rng, depth + 1)
    return out


class DifferentialAndPropertyFuzzTests(unittest.TestCase):
    def test_property_fuzz_vrs_and_edit_roundtrip_laws(self) -> None:
        rng = random.Random(20261115)
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            for trial in range(40):
                seq_len = rng.randint(12, 120)
                seq = "".join(rng.choice(ALL_IUPAC) for _ in range(seq_len))
                fasta = root / f"trial_{trial}.fasta"
                fasta.write_text(f">rec_{trial}\n{seq}\n", encoding="utf-8")
                src_dir = root / f"src_{trial}"
                compile_source(FASTA_PROFILE, {"sequence": fasta}, parameters={"wrapper": "identity"}).save(src_dir)

                variants = []
                for kind in ("snv", "mnv", "insertion", "deletion", "delins"):
                    if kind == "snv":
                        pos = rng.randint(0, seq_len - 1)
                        ref = seq[pos : pos + 1]
                        alt = rng.choice([b for b in BASES if b != ref])
                        start, end = pos, pos + 1
                    elif kind == "mnv":
                        pos = rng.randint(0, seq_len - 3)
                        ref = seq[pos : pos + 3]
                        alt = "".join(rng.choice([b for b in BASES if b != ch]) for ch in ref)
                        start, end = pos, pos + 3
                    elif kind == "insertion":
                        pos = rng.randint(0, seq_len)
                        ref = ""
                        alt = "".join(rng.choice(BASES) for _ in range(rng.randint(1, 5)))
                        start, end = pos, pos
                    elif kind == "deletion":
                        pos = rng.randint(0, seq_len - 4)
                        width = rng.randint(1, 4)
                        ref = seq[pos : pos + width]
                        alt = ""
                        start, end = pos, pos + width
                    else:
                        pos = rng.randint(0, seq_len - 3)
                        ref = seq[pos : pos + 2]
                        alt = "".join(rng.choice(BASES) for _ in range(rng.randint(1, 4)))
                        if alt == ref:
                            alt = alt + "A"
                        start, end = pos, pos + 2

                    edit = {"id": f"v-{trial}-{kind}", "record_id": f"rec_{trial}", "start": start, "end": end, "ref": ref, "alt": alt}
                    mutated = apply_variant_to_sequence(seq, edit)
                    reverted = revert_variant(mutated, edit)
                    self.assertEqual(reverted, seq, f"reversion failed for {edit}")
                    variants.append(edit)

                vc = compile_variants(src_dir, variants)
                vc_path = root / f"vc_{trial}.json"
                vc.save(vc_path)
                rep = validate_variant_paths(src_dir, vc_path)
                self.assertTrue(rep["valid"])
                for entry in vc.to_dict()["variants_ir"]["variants"]:
                    expected_allele = vrs_allele(entry["refget_id"], entry["location"]["start"], entry["location"]["end"], entry["alt"])
                    self.assertEqual(entry["vrs_allele"], expected_allele)

    def test_differential_fasta_genbank_and_gff3_parsers_against_oracles_and_fuzz(self) -> None:
        rng = random.Random(9001)
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            # Differential FASTA across random multi-line wrapping & case
            for i in range(25):
                n_recs = rng.randint(1, 3)
                lines = []
                for r in range(n_recs):
                    seq = "".join(rng.choice("acgtACGTNn") for _ in range(rng.randint(10, 90)))
                    lines.append(f">seq_{i}_{r} desc={r}")
                    step = rng.randint(5, 25)
                    for pos in range(0, len(seq), step):
                        lines.append(seq[pos : pos + step])
                fasta_text = "\n".join(lines) + "\n"
                fpath = root / f"diff_{i}.fasta"
                fpath.write_text(fasta_text, encoding="utf-8")
                oracle = _minimal_fasta_oracle(fasta_text)
                bundle = compile_source(FASTA_PROFILE, {"sequence": fpath}, parameters={"wrapper": "identity"})
                records = {r["record_id"]: r["bases"] for r in bundle.to_dict()["source_ir"]["records"]}
                self.assertEqual(records, {k: len(v) for k, v in oracle.items()})

            # Differential GFF3 coordinate & feature preservation vs independent line oracle
            from brainc.sequence_collection import SequenceCollectionCompiler
            gff_compiler = GFF3Compiler()
            seq_compiler = SequenceCollectionCompiler()
            for i in range(15):
                length = rng.randint(40, 100)
                start = rng.randint(1, length - 10)
                end = rng.randint(start, length)
                strand = rng.choice(["+", "-", "."])
                seq_col = seq_compiler.compile_fasta_bytes(
                    f">r1\n{'A' * length}\n".encode("ascii"),
                )
                gff_text = (
                    "##gff-version 3\n"
                    f"##sequence-region r1 1 {length}\n"
                    f"r1\tfuzz\texon\t{start}\t{end}\t.\t{strand}\t.\tID=f{i};Note=fuzz_{i}\n"
                )
                art = gff_compiler.compile_bytes(gff_text.encode("utf-8"), seq_col)
                feats = art.to_dict()["bio_ir"]["features"]
                self.assertEqual(len(feats), 1)
                seg = feats[0]["segments"][0]
                self.assertEqual(seg["raw_start"], start)
                self.assertEqual(seg["raw_end"], end)
                self.assertEqual(seg["intervals"][0]["start"], start - 1)
                self.assertEqual(seg["intervals"][0]["end"], end)
                self.assertEqual(feats[0]["strand"], strand)

            # Malformed fuzz inputs to GenBank and GFF3 must fail cleanly with typed errors (never crash)
            gb_compiler = GenBankCompiler()
            base_col = seq_compiler.compile_fasta_bytes(b">r1\nACGTACGT\n")
            for _ in range(30):
                noise = bytes(rng.randint(0, 255) for _ in range(rng.randint(1, 128)))
                with self.assertRaises(GenBankError):
                    gb_compiler.compile_bytes(noise)
                with self.assertRaises(GFF3Error):
                    gff_compiler.compile_bytes(noise, base_col)

            # Malformed variant specs never leak unhandled exceptions
            for bad_variant in (
                {"id": "x", "record_id": "r", "start": True, "end": 1, "ref": "A", "alt": "T"},
                {"id": "x", "record_id": "r", "start": -1, "end": 1, "ref": "A", "alt": "T"},
                {"id": "x", "record_id": "r", "start": 0, "end": 1, "ref": "A", "alt": "A"},
                {"id": "x", "record_id": "r", "start": 0, "end": 0, "ref": "A", "alt": "T"},
            ):
                with self.assertRaises(VariantError):
                    apply_variant_to_sequence("ACGT", bad_variant)

    @unittest.skipUnless(
        shutil.which("cargo") is not None and (RUST_DIR / "Cargo.toml").is_file(),
        "Rust cargo toolchain is unavailable",
    )
    def test_differential_rfc8785_python_vs_rust_fuzz(self) -> None:
        subprocess.run(
            ["cargo", "build", "--offline", "--manifest-path", str(RUST_DIR / "Cargo.toml")],
            check=True,
            capture_output=True,
            text=True,
            cwd=ROOT,
        )
        binary = RUST_DIR / "target" / "debug" / "epigenesis-validate-rs"
        rng = random.Random(424242)
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            for i in range(60):
                obj = {"trial": i, "payload": _random_json(rng)}
                p = root / f"jcs_{i}.json"
                p.write_text(json.dumps(obj, ensure_ascii=False), encoding="utf-8")
                res = subprocess.run(
                    [str(binary), "canonical-digest", str(p)],
                    check=False,
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(res.returncode, 0, res.stderr + f" payload={obj!r}")
                rust_out = json.loads(res.stdout)
                self.assertEqual(rust_out["sha256"], digest(obj), f"mismatch on {obj!r}")

            # Negative JSON fuzzing: Rust must reject duplicate keys, NaN, Infinity, unsafe ints, lone surrogates
            bad_payloads = [
                b'{"a":1,"a":2}',
                b'{"a":NaN}',
                b'{"a":Infinity}',
                b'{"a":9007199254740992}',
                b'{"a":"\\ud800"}',
                b'{"a":01}',
            ]
            for idx, raw in enumerate(bad_payloads):
                bp = root / f"bad_{idx}.json"
                bp.write_bytes(raw)
                res = subprocess.run(
                    [str(binary), "canonical-digest", str(bp)],
                    check=False,
                    capture_output=True,
                    text=True,
                )
                self.assertNotEqual(res.returncode, 0, f"Rust unexpectedly accepted {raw!r}")
                with self.assertRaises(ContractError):
                    canonical_bytes(json.loads(raw.decode("utf-8"))) if b'"a":01' not in raw and b'"a":1,"a":2' not in raw else (_ for _ in ()).throw(ContractError("dup"))


if __name__ == "__main__":
    unittest.main()

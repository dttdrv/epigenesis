"""Asymmetric 532-TF Boltzmann + multi-scale sequence predictor (v2) for Epigenesis.

Combines canonical 5-mer and 3-mer sequence shifts, asymmetric motif-potential gain/loss,
dual-scale strand-symmetric Boltzmann log-partition shifts (Delta ln(Z + Z_ns)), peak-site
disruption/creation (Delta max score), and reference-context-gated motif disruption across
all 532 JASPAR 2024 vertebrate transcription factors, predicting 0-72h neural-induction
regulatory perturbation trajectories on CPU in milliseconds.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if (REPO / "brainc").is_dir() and str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from brainc import _canonical as canonical

MODEL_DIR = Path(__file__).resolve().parent / "model"
if str(MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(MODEL_DIR))

import features  # noqa: E402
import predict as siraj_predict  # noqa: E402

MANIFEST_PATH = MODEL_DIR / "biophys_model.json"
_CACHE = None


def available() -> bool:
    return MANIFEST_PATH.is_file() and (MODEL_DIR / "biophys_weights.npz").is_file()


def identity() -> dict | None:
    if not available():
        return None
    raw = MANIFEST_PATH.read_bytes()
    manifest = json.loads(raw)
    return dict(
        model=manifest["format"],
        manifest_sha256=canonical.bytes_digest(raw),
        endpoint=manifest["endpoint"],
    )


def load():
    global _CACHE
    if _CACHE is not None:
        return _CACHE
    if not available():
        raise ValueError("biophysical v2 checkpoint not installed")
    raw = MANIFEST_PATH.read_bytes()
    manifest = json.loads(raw)
    if manifest.get("format") != "epigenesis-biophysical-temporal-v2":
        raise ValueError("unsupported biophysical model format")
    for name, expected in manifest["sha256"].items():
        if Path(name).name != name or hashlib.sha256((MODEL_DIR / name).read_bytes()).hexdigest() != expected:
            raise ValueError("biophysical artifact identity mismatch: " + name)
    with np.load(MODEL_DIR / "biophys_weights.npz", allow_pickle=False) as ks:
        coefs = np.array(ks["coefficients"], dtype=np.float64)
    n_motifs = len(manifest["motif_ids"])
    with np.load(MODEL_DIR / "weights.npz", allow_pickle=False) as ms:
        matrices = [np.array(ms[f"motif_{i}"], dtype=np.float64) for i in range(n_motifs)]
    # Pre-group normalized matrices by width for fast vectorized einsum
    by_w: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    grouped: dict[int, list[tuple[int, np.ndarray]]] = {}
    for mi, m in enumerate(matrices):
        m_n = (m - m.max(axis=0, keepdims=True)).astype(np.float32)
        grouped.setdefault(m.shape[1], []).append((mi, m_n))
    for w, items in sorted(grouped.items()):
        by_w[w] = (
            np.array([x[0] for x in items], dtype=np.int64),
            np.stack([x[1] for x in items], axis=0),
        )
    _CACHE = (manifest, coefs, matrices, by_w)
    return _CACHE


def _boltzmann_stats(sequences: list[str], by_w: dict, n_motifs: int, z_ns1: float, z_ns2: float):
    lut = np.zeros(128, dtype=np.int8)
    lut[ord("A")] = 0
    lut[ord("C")] = 1
    lut[ord("G")] = 2
    lut[ord("T")] = 3
    n = len(sequences)
    logz1 = np.full((n, n_motifs), math.log(z_ns1), dtype=np.float32)
    logz2 = np.full((n, n_motifs), math.log(z_ns2), dtype=np.float32)
    maxs = np.full((n, n_motifs), -20.0, dtype=np.float32)
    by_len: dict[int, list[int]] = {}
    for i, s in enumerate(sequences):
        by_len.setdefault(len(s), []).append(i)
    for L, idxs in by_len.items():
        mn = len(idxs)
        arr = np.frombuffer("".join(sequences[i] for i in idxs).encode("ascii"), dtype=np.uint8).reshape(mn, L)
        si = lut[arr]
        rci = 3 - si[:, ::-1]
        ohf = np.zeros((mn, 4, L), dtype=np.float32)
        ohr = np.zeros((mn, 4, L), dtype=np.float32)
        ar = np.arange(mn)[:, None]
        ac = np.arange(L)[None, :]
        ohf[ar, si, ac] = 1.0
        ohr[ar, rci, ac] = 1.0
        dest = np.array(idxs, dtype=np.int64)
        for w, (mis, ws) in by_w.items():
            if w > L:
                continue
            wf = np.lib.stride_tricks.sliding_window_view(ohf, w, axis=2)
            wr = np.lib.stride_tricks.sliding_window_view(ohr, w, axis=2)
            for s in range(0, mn, 256):
                e = min(mn, s + 256)
                sf = np.einsum("nclw,mcw->nlm", wf[s:e], ws, optimize=True)
                sr = np.einsum("nclw,mcw->nlm", wr[s:e], ws, optimize=True)
                zf = np.exp(sf.astype(np.float64)).sum(axis=1)
                zr = np.exp(sr.astype(np.float64)).sum(axis=1)
                logz1[dest[s:e, None], mis[None, :]] = np.log(zf + zr + z_ns1)
                logz2[dest[s:e, None], mis[None, :]] = np.log(zf + zr + z_ns2)
                maxs[dest[s:e, None], mis[None, :]] = np.maximum(sf.max(axis=1), sr.max(axis=1))
    return logz1, logz2, maxs


def predict(rows: list[dict]) -> list[dict]:
    manifest, coefs, matrices, by_w = load()
    siraj_predict.validate(rows, manifest["k"])
    motif_ids = manifest["motif_ids"]
    nm = len(motif_ids)
    z_ns1 = float(manifest["z_ns1"])
    z_ns2 = float(manifest["z_ns2"])
    n_k5 = len(features.vocabulary(manifest["k"]))
    n_k3 = len(features.vocabulary(manifest["k_short"]))
    n_seq = n_k5 + n_k3

    # Deduplicate sequences for ultra-fast scan & batch scoring
    u_map: dict[str, int] = {}
    u_seqs: list[str] = []
    for r in rows:
        for s in (r["reference"], r["sequence"]):
            if s not in u_map:
                u_map[s] = len(u_seqs)
                u_seqs.append(s)

    logz1, logz2, maxs = _boltzmann_stats(u_seqs, by_w, nm, z_ns1, z_ns2)
    k5_mat = features.whole_counts(u_seqs, manifest["k"]).astype(np.float32)
    k3_mat = features.whole_counts(u_seqs, manifest["k_short"]).astype(np.float32)

    ri = np.array([u_map[r["reference"]] for r in rows], dtype=np.int64)
    ai = np.array([u_map[r["sequence"]] for r in rows], dtype=np.int64)
    dz1 = (logz1[ai] - logz1[ri]).astype(np.float64)
    dz2 = (logz2[ai] - logz2[ri]).astype(np.float64)
    dmax = (maxs[ai] - maxs[ri]).astype(np.float64)
    dk5 = (k5_mat[ai] - k5_mat[ri]).astype(np.float64)
    dk3 = (k3_mat[ai] - k3_mat[ri]).astype(np.float64)
    dz1_l = np.minimum(0.0, dz1)
    dz1_g = np.maximum(0.0, dz1)
    dmax_l = np.minimum(0.0, dmax)
    dmax_g = np.maximum(0.0, dmax)
    ref_z = np.clip(logz2[ri].astype(np.float64) - math.log(z_ns2), 0.0, 10.0)
    ctx_l = np.tanh(ref_z) * dz1_l
    agg = np.column_stack([
        dz1_l.sum(axis=1) / 100.0,
        dz1_g.sum(axis=1) / 100.0,
        dz1_l.min(axis=1),
        dz1_g.max(axis=1),
        dmax_l.min(axis=1),
        dmax_g.max(axis=1),
        np.abs(dz1).sum(axis=1) / 100.0,
        ref_z.mean(axis=1),
        np.ones(len(rows)),
    ])
    md = features.deltas(rows, matrices, float(manifest["scale"]))
    md_l = np.minimum(0.0, md)
    md_g = np.maximum(0.0, md)

    X_seq = np.hstack([dk5, dk3])
    X_tf = np.hstack([md_l, md_g, dz1_l, dz1_g, dz2, dmax_l, dmax_g, ctx_l, agg])
    word_mat = X_seq @ coefs[:n_seq, :]
    motif_mat = X_tf @ coefs[n_seq:, :]
    total_mat = word_mat + motif_mat
    if not np.isfinite(total_mat).all():
        raise ValueError("nonfinite biophysical predictions")

    # Per-motif attribution across all 8 motif channels (mean over 7 timepoints)
    tf_blocks = [md_l, md_g, dz1_l, dz1_g, dz2, dmax_l, dmax_g, ctx_l]
    per_motif_mean = np.zeros((len(rows), nm), dtype=np.float64)
    per_motif_h72 = np.zeros((len(rows), nm), dtype=np.float64)
    offset = n_seq
    for blk in tf_blocks:
        w_blk = coefs[offset : offset + nm, :]
        per_motif_mean += blk * w_blk.mean(axis=1)[None, :]
        per_motif_h72 += blk * w_blk[:, -1][None, :]
        offset += nm

    out = []
    for idx, (row, wm, mm, tm) in enumerate(zip(rows, word_mat, motif_mat, total_mat, strict=True)):
        contribs = per_motif_mean[idx]
        nonzero = [j for j in range(nm) if abs(contribs[j]) > 1e-6]
        nonzero.sort(key=lambda j: (-abs(float(contribs[j])), motif_ids[j]))
        top_motifs = [
            dict(
                id=motif_ids[j],
                delta=float(dz1[idx, j]),
                effect=float(contribs[j]),
                h72_effect=float(per_motif_h72[idx, j]),
            )
            for j in nonzero[:5]
        ]
        out.append(
            dict(
                id=row["id"],
                effect=float(tm.mean()),
                word_effect=float(wm.mean()),
                motif_effect=float(mm.mean()),
                temporal=[float(v) for v in tm],
                top_motifs=top_motifs,
            )
        )
    return out


def select(predictions, rows, budget):
    return siraj_predict.select(predictions, rows, budget)

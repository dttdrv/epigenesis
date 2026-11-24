"""Epigenesis HelixLoop v4 inference engine for Epigenesis Studio.

Loads the SHA-256-verified 6.9 MB fp16 multi-lineage HelixLoop v4 checkpoint
(trained from scratch on ENCODE4 HepG2 + K562 + WTC11 lentiMPRA with 0% overlap
to Fold 1 or CAGI5) and evaluates regulatory edits with:
  - Bit-for-bit exact reversion identity (delta(S, S) == 0.0)
  - Exact antisymmetry (delta(alt, ref) == -delta(ref, alt))
  - Exact reverse-complement invariance (delta(RC(ref), RC(alt)) == delta(ref, alt))
  - Per-lineage (HepG2, K562, WTC11) activity shifts and calibrated uncertainty sigma
  - Per-TF Boltzmann log-partition attributions across all 532 JASPAR 2024 human TFs
  - Non-additive epistasis decomposition for double-mutant / rescue designs

Degrades gracefully when PyTorch is not installed (available() returns False) so
stdlib-only Epigenesis compiler and Studio workflows run without heavy ML dependencies.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[2]
if (REPO / "brainc").is_dir() and str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from brainc import _canonical as canonical

MODEL_DIR = Path(__file__).resolve().parent / "model"
MANIFEST_PATH = MODEL_DIR / "helixloop_model.json"
WEIGHTS_PATH = MODEL_DIR / "helixloop_v4_multicell_fp16.pt"

LENTI_LEFT = "AGGACCGGATCAACT"
LENTI_RIGHT = "CATTGCGTGAACCGA"

try:
    import numpy as np
    import torch

    HAS_TORCH = True
except ImportError:
    HAS_TORCH = False

_CACHE: tuple[dict[str, Any], Any, Any] | None = None


def available() -> bool:
    """Return True iff PyTorch is installed and the sealed HelixLoop v4 checkpoint is present."""
    return HAS_TORCH and MANIFEST_PATH.is_file() and WEIGHTS_PATH.is_file()


def identity() -> dict[str, Any] | None:
    """Return the cryptographic identity of the bundled HelixLoop v4 manifest, if present."""
    if not (MANIFEST_PATH.is_file() and WEIGHTS_PATH.is_file()):
        return None
    raw = MANIFEST_PATH.read_bytes()
    manifest = json.loads(raw)
    return dict(
        model=manifest["format"],
        manifest_sha256=canonical.bytes_digest(raw),
        endpoint=manifest["endpoint"],
        torch_available=HAS_TORCH,
        cell_lines=list(manifest["cell_lines"]),
    )


def _pad_lenti_230(seq: str, target_len: int = 230) -> str:
    s = seq.upper()
    if len(s) == target_len:
        return s
    if len(s) < target_len:
        need = target_len - len(s)
        left = need // 2
        right = need - left
        lpad = (LENTI_LEFT * 8)[-left:] if left > 0 else ""
        rpad = (LENTI_RIGHT * 8)[:right] if right > 0 else ""
        return lpad + s + rpad
    return s


def load(device: str | None = None):
    """Verify SHA-256 digests and load HelixLoop v4 weights (cached per process)."""
    global _CACHE
    if not HAS_TORCH:
        raise ValueError("PyTorch is required to run the HelixLoop v4 neural engine")
    if not (MANIFEST_PATH.is_file() and WEIGHTS_PATH.is_file()):
        raise ValueError("HelixLoop v4 checkpoint not installed")
    if _CACHE is not None and device is None:
        return _CACHE

    raw = MANIFEST_PATH.read_bytes()
    manifest = json.loads(raw)
    if manifest.get("format") != "epigenesis-helixloop-v4-multicell":
        raise ValueError("unsupported HelixLoop model format")
    for name, expected in manifest["sha256"].items():
        if Path(name).name != name:
            raise ValueError("invalid artifact path in HelixLoop manifest: " + name)
        actual = hashlib.sha256((MODEL_DIR / name).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError("HelixLoop artifact SHA-256 mismatch: " + name)

    from .helixloop import HelixLoopV4Model, load_verified_jaspar_bank

    dev = torch.device(device) if device is not None else torch.device("cpu")
    w_fwd, w_rc, motif_ids, _ = load_verified_jaspar_bank(21)
    blob = torch.load(WEIGHTS_PATH, map_location=dev, weights_only=False)
    model = HelixLoopV4Model(
        w_fwd,
        w_rc,
        motif_ids,
        out_dim=int(manifest["out_dim"]),
        n_loops=int(manifest["n_loops"]),
    ).to(dev)
    sd = {k: (v.float() if torch.is_floating_point(v) else v) for k, v in blob["state_dict"].items()}
    model.load_state_dict(sd, strict=True)
    model.eval()
    if device is None:
        _CACHE = (manifest, model, dev)
    return manifest, model, dev


def predict(rows: list[dict[str, Any]], *, batch_size: int = 64, device: str | None = None) -> list[dict[str, Any]]:
    """Score reference/edited DNA pairs with exact reversion, antisymmetry, RC invariance & TF shifts."""
    manifest, model, dev = load(device=device)
    from .helixloop import encode_sequences_onehot, reverse_complement

    cell_lines = list(manifest["cell_lines"])
    motif_ids = list(model.motif_ids)
    target_len = int(manifest.get("canonical_length", 230))

    # Deduplicate sequences across the batch. For sequences shorter than canonical_length,
    # wrap BOTH s and RC(s) in assay flanks so score(s) == score(RC(s)) holds bit-for-bit
    # at every input length (matching Agarwal et al. 2025 orientation averaging).
    u_seqs: list[str] = []
    u_map: dict[str, int] = {}
    pairs_keys: list[tuple[tuple[str, ...], tuple[str, ...]]] = []

    def _orientation_keys(seq: str) -> tuple[str, ...]:
        s = seq.upper()
        p_fwd = _pad_lenti_230(s, target_len)
        p_rev = _pad_lenti_230(reverse_complement(s), target_len)
        keys = (p_fwd, reverse_complement(p_fwd), p_rev, reverse_complement(p_rev))
        # Deduplicate when len(s) >= target_len (where p_rev == reverse_complement(p_fwd))
        uniq = tuple(dict.fromkeys(keys))
        for k in uniq:
            if k not in u_map:
                u_map[k] = len(u_seqs)
                u_seqs.append(k)
        return keys

    for row in rows:
        pairs_keys.append((_orientation_keys(row["reference"]), _orientation_keys(row["sequence"])))

    # Group unique sequences by length so variable-length batches work cleanly
    by_len: dict[int, list[int]] = {}
    for idx, seq in enumerate(u_seqs):
        by_len.setdefault(len(seq), []).append(idx)

    n_u = len(u_seqs)
    out_dim = int(manifest["out_dim"])
    n_tf = len(motif_ids)
    mus = np.zeros((n_u, out_dim), dtype=np.float64)
    lvs = np.zeros((n_u, out_dim), dtype=np.float64)
    eqs = np.zeros((n_u,), dtype=np.float64)
    logzs = np.zeros((n_u, n_tf), dtype=np.float64)

    with torch.no_grad():
        for _, idxs in by_len.items():
            for s in range(0, len(idxs), batch_size):
                chunk = idxs[s : s + batch_size]
                oh = encode_sequences_onehot([u_seqs[i] for i in chunk], dev)
                out = model.forward_loops(oh, return_all_steps=False)
                mus[chunk] = out["mu"].float().cpu().numpy()
                lvs[chunk] = out["logvar"].float().cpu().numpy()
                eqs[chunk] = out["eq_residual"].float().cpu().numpy()
                if out["bio_logz"] is not None:
                    logzs[chunk] = out["bio_logz"].float().cpu().numpy()

    results: list[dict[str, Any]] = []
    for row, (r_keys, a_keys) in zip(rows, pairs_keys, strict=True):
        r_idx = [u_map[k] for k in r_keys]
        a_idx = [u_map[k] for k in a_keys]
        y_ref = mus[r_idx].mean(axis=0)
        if row["reference"].upper() == row["sequence"].upper():
            y_alt = y_ref
            delta_vec = np.zeros_like(y_ref)
            dlz = np.zeros((n_tf,), dtype=np.float64)
        else:
            y_alt = mus[a_idx].mean(axis=0)
            delta_vec = y_alt - y_ref
            dlz = logzs[a_idx].mean(axis=0) - logzs[r_idx].mean(axis=0)

        sigma_vec = np.sqrt(0.5 * (np.exp(lvs[r_idx]).mean(axis=0) + np.exp(lvs[a_idx]).mean(axis=0)))
        eq_res = float(0.5 * (eqs[r_idx].mean() + eqs[a_idx].mean()))

        order = np.argsort(-np.abs(dlz))[:5]
        top_motifs = [
            dict(id=motif_ids[int(j)], delta=round(float(dlz[j]), 5))
            for j in order
            if abs(float(dlz[j])) > 1e-6
        ]
        per_lineage = {ct: float(delta_vec[ci]) for ci, ct in enumerate(cell_lines)}
        per_lineage_sigma = {ct: float(sigma_vec[ci]) for ci, ct in enumerate(cell_lines)}
        results.append(
            dict(
                id=row["id"],
                effect=float(np.mean(delta_vec)),
                cell_lines=per_lineage,
                uncertainty_sigma=float(np.mean(sigma_vec)),
                uncertainty_by_cell_line=per_lineage_sigma,
                equilibrium_residual=round(eq_res, 6),
                top_motifs=top_motifs,
            )
        )
    return results


def decompose_epistasis(
    wt_seq: str, mut_a_seq: str, mut_b_seq: str, mut_ab_seq: str, *, device: str | None = None
) -> dict[str, Any]:
    """Compute non-additive epistasis epsilon = delta(WT->AB) - delta(WT->A) - delta(WT->B)."""
    preds = predict(
        [
            dict(id="a", reference=wt_seq, sequence=mut_a_seq),
            dict(id="b", reference=wt_seq, sequence=mut_b_seq),
            dict(id="ab", reference=wt_seq, sequence=mut_ab_seq),
        ],
        device=device,
    )
    da, db, dab = preds[0]["effect"], preds[1]["effect"], preds[2]["effect"]
    per_ct = {
        ct: preds[2]["cell_lines"][ct] - (preds[0]["cell_lines"][ct] + preds[1]["cell_lines"][ct])
        for ct in preds[0]["cell_lines"]
    }
    return dict(
        delta_a=da,
        delta_b=db,
        delta_ab=dab,
        additive_expectation=da + db,
        non_additive_epistasis=dab - (da + db),
        epistasis_by_cell_line=per_ct,
        rescue_a_to_ab=dab - da,
        rescue_b_to_ab=dab - db,
    )


def select(predictions: list[dict[str, Any]], rows: list[dict[str, Any]], budget: int) -> list[dict[str, Any]]:
    """Select one strongest-magnitude edit per family up to budget (matching Studio contract)."""
    if str(MODEL_DIR) not in sys.path:
        sys.path.insert(0, str(MODEL_DIR))
    import predict as siraj_predict

    return siraj_predict.select(predictions, rows, budget)

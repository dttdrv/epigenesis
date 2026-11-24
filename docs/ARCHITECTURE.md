# Epigenesis Architecture & Roadmap

Epigenesis is organized into three strict layers so exploratory ML experiments never compromise compiler provenance, and compliance plumbing never gets in the way of real biology.

## The Three Layers

```text
┌────────────────────────────────────────────────────────────────────────────┐
│ Layer 3: Benchmarks, Audits & Frozen Research Ledger                       │
│   examples/helixloop-benchmark/ (ENCODE4, CAGI5, ADASTRA/UDACHA, GSE188264 │
│   + Borzoi fold3+4 leakage audit), research/ (13 frozen predictor tests)   │
├────────────────────────────────────────────────────────────────────────────┤
│ Layer 2: HelixLoop v4 Engine, Studio (0.4.0) & Developmental Pipeline      │
│   studio/ (helixloop.py, helixloop_engine.py, 6.9 MB fp16 weights,         │
│   biophys/kreimer/siraj CPU models), desktop/ (native macOS .dmg),         │
│   examples/developmental-pipeline/ (ODE + 3D CX3D + inverse synthesis)     │
├────────────────────────────────────────────────────────────────────────────┤
│ Layer 1: Zero-Dependency Provenance Compiler Core (1.7.0)                  │
│   brainc/ (FASTA/GenBank/GFF3/VRS 2.x, RFC 8785 seals, RO-Crate/OMEX),     │
│   SPEC.md, validators/rust/ (independent zero-dependency Rust validator)   │
└────────────────────────────────────────────────────────────────────────────┘
```

### Boundary Rules
1. **Layer 1 (`brainc/`, `SPEC.md`, `validators/rust/`):** Zero third-party Python dependencies. Everything stays deterministic and independently replayable in Python and Rust.
2. **Layer 2 (`studio/`, `desktop/`, `examples/developmental-pipeline/`):** Composes biological predictors and mechanistic models with explicit controls, SHA-256 manifests, and graceful fallback when PyTorch is absent.
3. **Layer 3 (`examples/helixloop-benchmark/`, `research/`):** Pre-registered or audited benchmarks with exact two-sided finite-sample Phipson–Smyth $p$-values, explicit leakage audits, and negative results kept on record (`Omega v1` bypass diagnostic, frozen linear-model wall in `research/`).

## Current Status & Next Milestones

### Completed
- **Tier-1 HelixLoop v4 (Apache-2.0, `6.9 MB` fp16, zero external weights):** Beats `OfficialMPRALegNet` across multi-lineage ENCODE4 lentiMPRA (`0.7946` vs `0.7721`), scales monotonically with loop depth ($T=1\to 6$), transfers zero-shot to 1,974 *in vivo* `ADASTRA`/`UDACHA` SNPs (5/6 panels won, $p=0.016$), and improves `GSE188264` non-additive epistasis (`0.5645` vs `-0.2141`, $p=0.017$) with exact `0.0` symmetries.
- **Tier-2 Open Foundation-Hybrid (`MIT` Flashzoi + `Apache-2.0` HelixLoop v4):** 4-replicate ensemble clears non-commercial `AlphaGenome-Encoder-FT` (`0.8708` vs `0.8670`), holds `0.8250` AUROC on 453 strictly Borzoi-unseen `fold3+fold4` *in vivo* SNPs, and reaches `r = 0.7513` on strictly test-only non-additive double-mutant epistasis (`vs 0.6461` linear head).
- **Cryptographic provenance & replay:** Full Python + Rust validation across DNA sources, VRS 2.x edits, Development Modules, and Studio receipts.

### Still Open
1. **Long-range genomic context (524 kb–1 Mb):** HelixLoop v4 currently operates on 200–230 bp CRE windows (plus Flashzoi's convolutional stem). Extending the iterative state tracker across full-length Flashzoi transformer blocks for distal enhancer–gene linking is the next architectural step.
2. **Calibrating developmental ODEs against real perturbations:** Mechanistic GRN and CX3D parameters still use published literature units (`biological_acceptance = false`). Fitting and falsifying circuits directly against endogenous CRISPR/Perturb-seq trajectories remains open.
3. **Dynamic multi-cell morphogenesis:** Re-evaluating compiled cis-regulatory programs dynamically as cells divide and sense local morphogen gradients in 3D, rather than lowering a single steady-state snapshot.

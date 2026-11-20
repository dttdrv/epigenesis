# Epigenesis Architecture & Scientific Roadmap

Epigenesis is built around one clear objective: **general, honest, replayable DNA-to-development simulation software**—moving step by step from DNA sequence and genetic variation, through chromatin-gated gene regulation and multi-gene dynamics, to 3D neural morphogenesis and circuit activity.

To prevent exploratory experiments from compromising core software trust—and to prevent compliance plumbing from hijacking the biological goal—Epigenesis is organized into **three strict layers**.

## The Three Layers

```text
┌────────────────────────────────────────────────────────────────────────────┐
│ Layer 3: Exploratory Research & Negative-Result Ledger                     │
│   research/ (13 frozen predictor tests), examples/neural-v2f/ (V2F radar), │
│   examples/mpra-replay/ (GSE188264 provenance demo), neural-predictor/     │
│   Rule: biological_acceptance = false; negative results reported verbatim. │
├────────────────────────────────────────────────────────────────────────────┤
│ Layer 2: Unified Composable Developmental Engine                           │
│   examples/developmental-pipeline/ (`epigenesis develop`), neural-circuit, │
│   coding-consequences, neural-construction (CX3D + LIF), Studio            │
│   Rule: Composes DNA + VRS edits → chromatin/TF/coding → GRN uncertainty   │
│         → 3D neural morphogenesis & spiking through sealed brainc bundles. │
├────────────────────────────────────────────────────────────────────────────┤
│ Layer 1: Provenance Compiler Core (Zero-Dependency Foundation)             │
│   brainc/ (sources, VRS 2.x variants, Development Module IR, RO-Crate/     │
│   OMEX export), brainc/schemas/, SPEC.md, validators/rust/                 │
│   Rule: Neutral, deterministic, independently replayed in Python & Rust.   │
└────────────────────────────────────────────────────────────────────────────┘
```

### Boundary Rules
1. **If a change makes source/variant/module sealing or replay more exact**, it belongs in **Layer 1 (`brainc/`, `SPEC.md`, `validators/rust/`)**.
2. **If a change composes biological mechanisms into a general developmental workflow**, it belongs in **Layer 2 (`examples/developmental-pipeline/`, `studio/`, `epigenesis_cli.py`)**, with explicit assumptions, parameter-uncertainty ensembles, and built-in causal/negative controls.
3. **If a change tests an unproven scientific hypothesis or predictor prototype**, it belongs in **Layer 3 (`research/`, exploratory `examples/`)**, with `biological_acceptance = false` and pre-registered/frozen evaluation before any biological claim is made.

## Honest Scientific Status (What Works vs. What Is Still Open)

### Established Today (Software & Mechanistic Controls)
- **Cryptographic provenance & cross-language replay:** Sealed DNA sources, GA4GH VRS 2.x reference-bound edits, RFC 8785 Development Modules, Draft 2020-12 schemas, RO-Crate 1.1 / COMBINE OMEX export, and independent Python + zero-dependency Rust replay validators.
- **End-to-end composable developmental simulation (`epigenesis develop`):** Multi-locus cis-enhancer + coding ORF coupling through explicit chromatin accessibility/methylation gates into general `N`-gene GRNs (with deterministic parameter-uncertainty ensembles, data-to-mechanism falsification checks, exact allele reversion, and closed-chromatin null controls) and 3D CX3D morphogenesis + LIF spiking assays.

### Open Scientific Bottlenecks (Next Real Milestones)
1. **Breaking past the linear DNA-effect predictor wall:** Thirteen frozen tests in [`research/`](../research/README.md) proved that local ±18 bp linear 5-mer/PWM models hit a cross-study wall (`r ≈ 0.049`, dominated by GC composition). Genuine progress requires pretrained long-range sequence/chromatin backbones (e.g., Enformer/Borzoi/ChromBPNet embeddings) paired with cell-type chromatin context and count-uncertainty modeling (BIRD)—while keeping `Siraj RESERVE` sealed until external transfer beats baselines.
2. **Calibrating developmental mechanisms against endogenous perturbations:** Current GRN and CX3D parameters use uncalibrated literature/synthetic model units (`biological_acceptance = false`). The `evaluate_observations` falsification hook in `developmental-pipeline` is the entry point for fitting and falsifying regulatory networks directly against endogenous CRISPR/Perturb-seq or scMultiome differentiation trajectories.
3. **Scaling from single-cell GRNs to multi-cell tissue & cortical connectivity:** Extending the composed DNA → GRN → CX3D coupling so cells dynamically re-read local morphogen fields and re-evaluate their compiled cis-regulatory programs during 3D division and wiring, rather than lowering a single steady-state snapshot.

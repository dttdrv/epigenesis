# Variant-to-Circuit Enhancer Scanner

This experiment bridges arbitrary-length regulatory DNA enhancers directly to arbitrary $N$-gene regulatory circuits (default: `neural-tube` `Pax6`/`Olig2`/`Nkx2.2`/`Irx3` from Herrera-Delgado et al. 2020, plus `mutual-toggle`, tristable `tri-lineage-switch`, or custom 1–64 gene JSON `CircuitSpec`s) without restricting edits to a pre-annotated 9-base core.

Instead of stopping at passive promoter/MPRA log2 activity, it computes a strand-symmetric thermodynamic Boltzmann partition sum over all sliding windows using Hallikas et al. (2006) `GLI2`/`GLI3` profiles, any of the **532 bundled JASPAR 2024 motifs**, or **custom PWM JSON/dicts** (`q = Z_variant / Z_reference`) and integrates the coupled nonlinear ODE under endogenous feedback, spatial morphogen gradients, and Shh pulse-chase histories. Every featured variant pair compiles into a sealed `brainc` v2/v3 Development Module and verifies exact allele/mediator reversion.

## Key Result: Cryptic Bifurcations

On real human enhancer `chr15:90936956-90937127` (171 bp from GSE188264, scanning all `513` single-nucleotide substitutions at Shh signal `S = 0.20`):

- **7 cryptic bifurcations** (`|Δreporter| < 0.02`, yet endogenous `Nkx2.2 -> Olig2` cell-fate switch), alongside **20 concordant switches** and **486 neutral** variants.
- **Featured cryptic switch `60:C>A`:** `q = 0.042455`, passive reporter shift `−0.0124`, endogenous `Nkx2.2` change `−0.9318` (75.2× nonlinear amplification from the coupled Pax6/Olig2/Nkx2.2/Irx3 bistable switch).
- **Spatial domain contraction:** P3 (`Nkx2.2 > 0.5`) boundary shifts from `max x = 0.25` (reference) to `0.125` (variant) across 41 ventral-to-dorsal gradient positions.
- **Pulse-chase GliR derepression:** Primed `[(5.0, 1.0), (15.0, 0.05)]` ends at `Nkx2.2 = 0.8985` under competitive Gli binding versus `0.8638` reference, whereas `activation-only` ablation collapses to `0.000136`.
- **Multi-locus & second-site rescue:** Across 5 curated GSE188264 human enhancers (`benchmark`), suboptimal enhancers such as `chr18:56475717-56475888` (`165:G>T` rescued by `163:T>C`, `q = 1.0001`) and `chr1:175099416-175099587` (`149:G>T` rescued by `147:A>T`, `q = 0.4111`) support genuine second-site fate rescue back to `Nkx2.2`, whereas `chr15:90936956-90937127` relies on a single dominant core.
- **Generalized motifs, N-gene circuits, indels & CX3D bridge:** Supports wiring Hallikas `GLI2`/`GLI3`, any of the 532 JASPAR 2024 PWMs (e.g. `MA0069.1`, `MA0463.3`), or custom PWM JSONs into `neural-tube`, `mutual-toggle`, `tri-lineage-switch`, or custom 1–64 gene specs (`--circuit`, `--protein`, `--target`, `--mapping`), comparing arbitrary SNV/MNV/insertion/deletion pairs into sealed `brainc` v2/v3 Development Modules, and chaining any steady state into `neural-construction` CX3D + LIF controls (`bridge` / `--with-construction`).
- **Exact 9-bp parity & causal controls:** At `epsilon = 0.0`, 9-bp `GACCACCCA -> GAACACCCA` reproduces [`examples/neural-patterning`](../neural-patterning/README.md) `profile_ratio = 0.04235549749669733` (`< 1e-15`). Reversion (`max_state_diff = 0.0`), mediator rescue, phenocopy, passive-reporter isolation, and RK4 step-doubling refinement (`9.10e-10 < 1e-4`) all pass.

## Run

Python 3.11+ standard library only (no NumPy or network required for core physics and ODE; optional passive Kreimer/Siraj reporter scores attach automatically when NumPy is installed, and `bridge` / `--with-construction` uses Java 21 + `patch` via [`examples/neural-construction`](../neural-construction/README.md)).

```sh
python3.11 examples/neural-circuit/circuit.py list-circuits
python3.11 examples/neural-circuit/circuit.py list-motifs
python3.11 examples/neural-circuit/circuit.py demo --output /tmp/neural-circuit
python3.11 examples/neural-circuit/circuit.py scan --fasta my_enhancer.fasta --circuit tri-lineage-switch --target All --output /tmp/my-scan
python3.11 examples/neural-circuit/circuit.py compare --reference ref.fa --variant alt.fa --circuit mutual-toggle --protein MA0069.1 --output /tmp/my-compare
python3.11 examples/neural-circuit/circuit.py benchmark --output /tmp/circuit-benchmark
python3.11 examples/neural-circuit/circuit.py bridge --reference ref.fa --variant alt.fa --mapping map.json --output /tmp/circuit-bridge
python3.11 tests/verify.py neural-circuit
```

Also available via the unified `epigenesis` CLI (`epigenesis circuit`, `epigenesis scan`), live in [Epigenesis Studio](../../studio/README.md) under the **Enhancer Circuit** tab (`epigenesis-studio`), and via `python3 epigenesis.py circuit` in the release launcher.

## Scope and Limitations

- Conditional mechanistic model: transferring a measured in vitro PWM partition-sum ratio onto promoter affinity in ODE parameters is an explicit biological assumption.
- Concentrations, distances, and times are uncalibrated model units; nucleosome occupancy, co-factors, 3D chromatin looping, and embryo-level phenotype prediction are not established (`biological_acceptance = False`).

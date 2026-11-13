# Variant-to-Circuit Enhancer Scanner

This experiment bridges arbitrary-length regulatory DNA enhancers directly to the four-gene neural-tube regulatory circuit (`Pax6`, `Olig2`, `Nkx2.2`, `Irx3`; Herrera-Delgado et al. 2020) without restricting edits to a pre-annotated 9-base core.

Instead of stopping at passive promoter/MPRA log2 activity, it computes a strand-symmetric thermodynamic Boltzmann partition sum over all sliding 9-mer windows using Hallikas et al. (2006) `GLI2`/`GLI3` Table S1 binding ratios (`q = Z_variant / Z_reference`) and integrates the coupled nonlinear ODE under endogenous feedback, spatial morphogen gradients, and Shh pulse-chase histories. Every featured variant pair compiles into a sealed `brainc` Development Module and verifies exact allele/mediator reversion.

## Key Result: Cryptic Bifurcations

On real human enhancer `chr15:90936956-90937127` (171 bp from GSE188264, scanning all `513` single-nucleotide substitutions at Shh signal `S = 0.20`):

- **7 cryptic bifurcations** (`|Δreporter| < 0.02`, yet endogenous `Nkx2.2 -> Olig2` cell-fate switch), alongside **20 concordant switches** and **486 neutral** variants.
- **Featured cryptic switch `60:C>A`:** `q = 0.042455`, passive reporter shift `−0.0124`, endogenous `Nkx2.2` change `−0.9318` (75.2× nonlinear amplification from the coupled Pax6/Olig2/Nkx2.2/Irx3 bistable switch).
- **Spatial domain contraction:** P3 (`Nkx2.2 > 0.5`) boundary shifts from `max x = 0.25` (reference) to `0.125` (variant) across 41 ventral-to-dorsal gradient positions.
- **Pulse-chase GliR derepression:** Primed `[(5.0, 1.0), (15.0, 0.05)]` ends at `Nkx2.2 = 0.8985` under competitive Gli binding versus `0.8638` reference, whereas `activation-only` ablation collapses to `0.000136`.
- **Exact 9-bp parity & causal controls:** At `epsilon = 0.0`, 9-bp `GACCACCCA -> GAACACCCA` reproduces [`examples/neural-patterning`](../neural-patterning/README.md) `profile_ratio = 0.04235549749669733` (`< 1e-15`). Reversion (`max_state_diff = 0.0`), mediator rescue, phenocopy, passive-reporter isolation, and RK4 step-doubling refinement (`9.10e-10 < 1e-4`) all pass.

## Run

Python 3.11+ standard library only (no NumPy or network required for core physics and ODE; optional passive Kreimer/Siraj reporter scores attach automatically when NumPy is installed).

```sh
python3.11 examples/neural-circuit/circuit.py demo --output /tmp/neural-circuit
python3.11 examples/neural-circuit/circuit.py scan --fasta my_enhancer.fasta --protein GLI2 --signal 0.20 --output /tmp/my-scan
python3.11 tests/verify.py neural-circuit
```

Also integrated live into [Epigenesis Studio](../../studio/README.md) under the **Enhancer Circuit** tab (`epigenesis-studio`) and via `python3 epigenesis.py circuit` in the release launcher.

## Scope and Limitations

- Conditional mechanistic model: transferring a measured in vitro GLI2/GLI3 partition-sum ratio onto `Nkx2.2` promoter affinity in Herrera-Delgado et al. (2020) ODE parameters is an explicit biological assumption.
- Concentrations, distances, and times are uncalibrated model units; nucleosome occupancy, co-factors, 3D chromatin looping, and embryo-level phenotype prediction are not established (`biological_validation = False`).

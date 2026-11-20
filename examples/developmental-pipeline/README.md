# Unified Composable Developmental Pipeline

Connects Epigenesis's sealed core compiler (`brainc`) directly to its biological modeling layers in **one composable, end-to-end pipeline**:

```text
Multi-locus DNA FASTA
  → brainc.source-descriptor/v1
  → Reference-bound GA4GH VRS 2.x edits (brainc.variant-collection/v1)
  → Chromatin-gated multi-TF enhancer binding + NCBI coding ORF translation
  → Sealed brainc.development-bundle/v1 (+ optional Rust cross-check)
  → N-gene GRN dynamics + deterministic Parameter-Uncertainty Ensemble
  → Built-in causal controls (exact allele reversion, closed-chromatin cis null, passive reporter)
  → Optional empirical falsification against observed gene deltas
  → Optional 3D CX3D morphogenesis & LIF spiking transmission assay
  → Optional RO-Crate 1.1 & COMBINE/OMEX portable export
```

## Why This Exists

Earlier experiments in `examples/` lived in isolated silos (coding translation in `coding-consequences`, single-locus enhancers in `neural-circuit`, 3D morphogenesis in `neural-construction`). This pipeline composes them through a single declarative `epigenesis.developmental-program/v1` specification and seals every stage with `brainc`:

- **Multi-locus & multi-modal:** Wire multiple cis-regulatory enhancers (any Hallikas `GLI2`/`GLI3`, `532` JASPAR 2024 motifs, or custom PWMs) *and* coding ORFs into the same `N`-gene regulatory network.
- **Explicit chromatin context:** Every enhancer declares accessibility $a_\ell \in [0, 1]$ and methylation repression $m_\ell \in [0, 1]$ with log-linear gating $q_{\text{eff},\ell} = \exp(a_\ell (1 - m_\ell) \ln q_{\text{raw},\ell})$. Setting $a_\ell = 0$ acts as an exact closed-chromatin negative control ($q_{\text{eff}} = 1.0$).
- **Parameter-uncertainty ensembles:** Instead of brittle single-parameter claims, kinetic ODE parameters are perturbed deterministically across $K$ draws to report fate-transition probabilities and 90% delta intervals.
- **Data-to-mechanism falsification:** Supplying empirical observations (`--observations`) runs explicit residual & sign checks and marks `falsified: true` whenever the mechanism fails to match data.
- **3D neural morphogenesis & spiking:** Lowers reference, variant, and exact allele-reversion states into sealed 480-bp CX3D recipes and runs 3D cell growth + LIF transmission assays.

## Commands

```sh
# Fast pure-Python multi-locus demo (neural-tube enhancers + Pax6 coding ORF + uncertainty ensemble)
python3.11 epigenesis_cli.py develop --skip-cx3d --output /tmp/epigenesis-dev-fast

# Full end-to-end demo including 3D CX3D morphogenesis, LIF spiking, and RO-Crate/OMEX export (requires JDK 21)
python3.11 epigenesis_cli.py develop --export-bundles --output /tmp/epigenesis-dev-full

# Tri-lineage JASPAR chromatin demo
python3.11 epigenesis_cli.py develop --preset tri-lineage-chromatin --skip-cx3d --output /tmp/epigenesis-tri

# Custom multi-locus FASTA + program spec + variants + optional empirical observations
python3.11 epigenesis_cli.py develop \
  --fasta my_loci.fasta \
  --spec my_program.json \
  --variants my_edits.json \
  --observations my_obs.json \
  --output /tmp/epigenesis-custom

# Release verification gate
python3.11 tests/verify.py developmental-pipeline

# Out-of-sample 532-TF biophysical mechanism discovery on real human GSE188264 data (requires NumPy)
python3.11 epigenesis_cli.py discover --output /tmp/epigenesis-discover
```

## Scientific Scope

`biological_acceptance` is explicitly **`false`**. Passing software provenance, VRS 2.x verification, chromatin-gating controls, exact allele reversion, and deterministic CX3D+LIF replay establishes **mechanistic software integrity**, not *in vivo* embryonic or clinical phenotype validity.

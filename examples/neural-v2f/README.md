# Variant-to-Function (V2F) Tipping-Point & Rescue Radar

Detects **passive-reporter-silent, endogenous-circuit tipping** non-coding variants across full developmental signaling regimes—explaining why standard linear/MPRA constructs or single-dose assays can miss causal enhancers shifts, quantifying exact morphogen tipping-point shifts ($\Delta S_{\mathrm{crit}}$), proposing **dual cis/trans rescues**, and falsifying claims against **exact Altschul–Erickson dinucleotide shuffles** and **column-scrambled motif PWMs**.

## Plain-English Summary: What Problem Does This Solve?

In real genetics and cell biology, one of the hardest unsolved bottlenecks ("Variant-to-Function", V2F) is that a tiny DNA spelling change outside a gene can look completely silent in a fast dish-based **reporter assay** (or get scored near zero by linear sequence predictors) yet still push a developing cell over a **developmental tipping point** when feedback loops are active.

This radar automates four questions on any enhancer DNA sequence:
1. **Across the whole morphogen range (not just one dose), which mutations flip cell fate while keeping passive reporter readouts flat (`|Δreporter| < 0.02`)?**
   When optional NumPy sequence predictors (**Siraj Huber-linear** + **Kreimer 0–72h motif-temporal**) are available, it also flags **"doubly silent"** mutations where *both* the mechanistic reporter AND empirical MPRA sequence predictors stay near zero (`< 0.15`).
2. **How much does the mutation shift the developmental tipping point ($\Delta S_{\mathrm{crit}}$) or lock the fate out across $S \in [0, 1]$?**
3. **Can the defect be rescued?**
   - **Cis epistatic rescue:** Can a second DNA edit inside the same enhancer restore wild-type fate across the whole signaling regime (verified by re-sweeping the double mutant, not just at one point)?
   - **Trans signal rescue:** How much morphogen/pathway dose shift ($\Delta S_{\mathrm{rescue}}$) restores wild-type fate?
4. **Is this real motif-and-circuit structure or just fluff?**
   - **Matched dinucleotide shuffles:** Generates $N$ exact Altschul–Erickson Eulerian dinucleotide-preserving decoys (`Counter(dinucleotides)` verified identical) and compares cryptic bifurcation rates + empirical $p$-value.
   - **Scrambled PWM control:** Permutes motif matrix columns to test sequence-word specificity.
   - **Causal & cryptographic `brainc` verification:** Verifies exact allele reversion (`max_state_diff = 0.0`), passive-reporter isolation (`0.0`), and sealed SHA-256 `brainc` v2/v3 Development Modules for both the top candidate and cis-rescue pair.

## Commands

Core physics, regime scanning, dinucleotide shuffles, rescues, CSV/HTML reports, and `brainc` sealing require **Python 3.11+ standard library only**. Sequence predictors (`Siraj` + `Kreimer`) attach automatically when NumPy is installed.

```sh
# 1. Full radar on a curated human enhancer from GSE188264 (writes HTML, CSV, JSON, brainc bundles)
python3.11 epigenesis_cli.py radar --locus chr15:90936956-90937127 --output /tmp/v2f-demo

# 2. Scan any custom FASTA enhancer across any circuit & TF motif (Hallikas GLI2/GLI3, 532 JASPAR PWMs, or custom JSON)
python3.11 epigenesis_cli.py radar --fasta my_enhancer.fasta --circuit tri-lineage-switch --target All --output /tmp/v2f-fasta

# 3. Score a targeted variant table (JSON/CSV/TSV supporting SNVs, MNVs, insertions, deletions)
python3.11 epigenesis_cli.py radar --variants edits.csv --output /tmp/v2f-variants

# 4. Multi-locus benchmark across all 5 curated GSE188264 loci + real measured CRE perturbations + null controls
python3.11 epigenesis_cli.py radar --benchmark --output /tmp/v2f-bench

# 5. Direct script & release verification gate
python3.11 examples/neural-v2f/radar.py demo --output /tmp/v2f-direct
python3.11 tests/verify.py neural-v2f
```

## Outputs

- `radar_report.json` (`epigenesis.neural-v2f-radar/v1`): full regime sweep, tipping points ($\Delta S_{\mathrm{crit}}$), null-control statistics, causal checks, and `biological_acceptance: false`.
- `wetlab_candidates.csv` / `wetlab_candidates.json`: prioritized single-edit, cis-rescue double-edit, and trans-signal shift table.
- `index.html`: standalone report with plain-English summary, morphogen dose-response curves, passive-vs-endogenous scatter plot, candidate table, and null-control verification.
- `compiled/top_candidate` & `compiled/cis_rescue`: verified SHA-256 `brainc` v2/v3 Development Modules.

## Verified Multi-Locus Summary (`benchmark --shuffles 5`, NumPy attached)

- **5 curated GSE188264 loci** (`chr15`, `chr13`, `chr21`, `chr18`, `chr1`; `S ∈ [0.06, 0.85]`, 9 steps): **4 cryptic bifurcations** detected (`62:A>C`, `47:T>G`, `58:A>C`, `163:T>G`), and **all 4 are doubly silent** (`max |Δreporter| < 0.015`, `|Siraj| < 0.083`, `|Kreimer mean| < 0.053`, yet `max |ΔNkx2.2| ≈ 0.929` and $\Delta S_{\mathrm{crit}} \in [+0.0931, +0.0953]$).
- **Dual rescues:** Trans morphogen dose shifts (`ΔS ∈ [+0.002887, +0.005084]`) restore fate for **4 / 4**; full-regime second-site cis epistatic rescues succeed on **2 / 4** (`chr13` `118:T>A`, `chr18` `166:G>T`), while single-core loci (`chr15`, `chr21`) return `null`.
- **85 measured CRE perturbations in GSE188264:** `52` switch fate cleanly and `0` are cryptic (the designed MPRA CRE perturbations strongly disrupt whole cores rather than sitting on reporter-silent boundaries).
- **Null controls (`n = 5` dinucleotide shuffles per locus):** `null_dinucleotide_mean` `0.0–0.4` (`0.16` overall mean across 25 shuffled decoys vs `0.80` on real loci); column-scrambled PWM produces `0` cryptic hits on 4/5 loci and `2` ratio-only noise hits on `chr15` (reported verbatim alongside `30,295.8×` weaker baseline partition sum in `demo --shuffles 10`).

## Scope & Honest Limitations

- **Conditional mechanistic model (`biological_acceptance = False`):** Mapping PWM Boltzmann partition ratios onto promoter affinity inside developmental ODEs is a structural computational hypothesis; uncalibrated model concentrations/times, chromatin accessibility, 3D looping, and *in vivo* clinical/embryonic effects are **not** established.
- **Null-control candor:** When a dinucleotide-shuffled decoy or column-scrambled PWM happens to exhibit cryptic switches due to weak baseline $Z_{\mathrm{ref}}$ under self-normalized $q = Z_{\mathrm{mut}}/Z_{\mathrm{ref}}$, or when a single-core enhancer (`chr15:90936956-90937127`) admits no compensating second-site cis rescue, the radar reports those negative findings verbatim.

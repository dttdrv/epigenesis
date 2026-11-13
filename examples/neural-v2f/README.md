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

## Verified Multi-Locus Summary ($Z_{\mathrm{ns}} = 2 \times 10^{-4}$, NumPy attached)

- **Point-scan vs. regime-wide cryptic reconciliation (`7` / `33` / `4`):** A single-dose scan (`circuit.py` at $S = 0.20$) flags `7` cryptic SNVs on `chr15` and `33` total across the `5` curated loci (`chr15` 7, `chr13` 7, `chr21` 7, `chr18` 6, `chr1` 6), while `93` edits are reporter-silent at *at least one* step (`point_silent_switches`). Requiring passive reporter silence across the *entire* signaling grid ($S \in [0.06, 0.85]$, 9 steps; $\max_S |\Delta\mathrm{reporter}(S)| < 0.02$) filters out `29 / 33` edits whose reporter delta exceeds `0.02` at higher or lower morphogen doses, leaving **`4` strict regime-wide cryptic bifurcations**:

| Locus | Point Scan (`S = 0.20`) | Any-Step Point-Silent (`[0.06, 0.85]`) | Strict Regime-Wide Cryptic | Top Strict Edit |
|---|---:|---:|---:|---|
| `chr15:90936956-90937127` | 7 | 19 | 1 | `62:A>C` (`|Kreimer mean| = 0.000506`) |
| `chr13:26545315-26545486` | 7 | 19 | 1 | `47:T>G` (`|Kreimer mean| = 0.052383`) |
| `chr21:37477109-37477280` | 7 | 19 | 1 | `58:A>C` (`|Kreimer mean| = 0.018337`) |
| `chr18:56475717-56475888` | 6 | 18 | 1 | `163:T>G` (`|Kreimer mean| = 0.020890`) |
| `chr1:175099416-175099587` | 6 | 18 | 0 | `null` |
| **Total (`5` curated loci)** | **33** | **93** | **4** | **4 / 4 doubly silent** |

- **Kreimer-only `doubly_silent` gate & Siraj informational caveat:** `sequence_silent` / `doubly_silent` gate strictly on the Kreimer 0–72h motif-temporal predictor (`|kreimer_mpra_mean| < 0.15`); the bundled Siraj Huber-linear score is logged strictly as informational context because frozen cross-study validation measured mean correlation `r = 0.049`. All `4 / 4` curated regime-wide cryptic switches are doubly silent (`max |Δreporter| < 0.015`, `|Kreimer mean| < 0.053`, yet `max |ΔNkx2.2| ≈ 0.929` and $\Delta S_{\mathrm{crit}} \in [+0.0918, +0.0949]$).
- **Dual rescues:** Trans morphogen dose shifts (`ΔS ∈ [+0.001630, +0.004718]`) restore fate for **4 / 4**; full-regime second-site cis epistatic rescues succeed on **2 / 4** (`chr13` `118:T>A`, `chr18` `166:G>T`), while single-core loci (`chr15`, `chr21`) return `null`.
- **85 measured CRE perturbations in the 5 curated GSE188264 loci:** `52` switch fate cleanly and `0` are cryptic (the designed MPRA CRE perturbations strongly disrupt whole cores rather than sitting on reporter-silent boundaries).
- **Null controls (`n = 1,000` default dinucleotide shuffles; `benchmark --shuffles 5` across 5 loci):** Partition ratios use $q = (Z_{\mathrm{mut}} + Z_{\mathrm{ns}}) / (Z_{\mathrm{ref}} + Z_{\mathrm{ns}})$ with non-specific binding floor $Z_{\mathrm{ns}} = 2 \times 10^{-4}$ (auto-disabled when `--epsilon 0` for exact Hallikas parity). Under this default floor:
  - **1,000 dinucleotide shuffles on `chr15` (`DEFAULT_SHUFFLES = 1000`):** `null_cryptic_mean = 0.155` (`std = 1.3561`, `max = 19`, `empirical_p_ge = 0.04995`, **`6.45×`** enrichment over null, **`305.6×`** partition-sum separation), and column-scrambled PWM produces **`0` cryptic hits** (**`30,295.8×`** partition-sum separation).
  - **5-locus benchmark (`n = 5` shuffles per locus):** Column-scrambled PWMs produce `0` cryptic hits on **5 / 5** loci; per-locus `dinucleotide_null_mean` is `0.0` (`chr15`, `chr13`, `chr1`), `0.2` (`chr21`), and `0.4` (`chr18`) — averaging `0.12` vs `0.80` on real loci.
- **Frozen held-out evaluation on `118` unseen GSE188264 `test` loci (`1,885` unseen measured alleles; [`heldout-freeze/RESULTS.md`](heldout-freeze/RESULTS.md)):** Sealed SHA-256 predictions before single-shot scoring (`biological_acceptance = false`). `8 / 118` unseen loci harbor `24` strict SNV cryptic hits (`24 / 24` doubly silent; top unseen locus `chr3:20371482-20371653` has `14` cryptic SNVs vs `0 / 100` dinucleotide decoys, `empirical_p_ge = 0.009901`, `0` scrambled-PWM hits). Across `1,885` measured `test` alleles, GLI2 passive reporter delta does **not** predict empirical `|measured α|` (`Pearson r = 0.012258`, `Spearman ρ = 0.019435`), and although predicted cryptic CRE alleles (`7`) have lower mean `|measured α|` (`1.258827`) than concordant switches (`15`, `2.471281`), `0 / 7` are empirically silent (`< 0.15`).

## Scope & Honest Limitations

- **Conditional mechanistic model (`biological_acceptance = False`):** Mapping PWM Boltzmann partition ratios onto promoter affinity inside developmental ODEs is a structural computational hypothesis; uncalibrated model concentrations/times, chromatin accessibility, 3D looping, and *in vivo* clinical/embryonic effects are **not** established.
- **Null-control & held-out candor:** Even with the non-specific floor $Z_{\mathrm{ns}} = 2 \times 10^{-4}$, a dinucleotide-shuffled decoy can occasionally form a weak binding site (`chr15` 1,000-shuffle `max = 19`, `chr18` mean `0.4`, `chr21` mean `0.2`), single-core enhancers (`chr15:90936956-90937127`) admit no compensating second-site cis rescue, and single-TF passive reporter deltas do not predict multi-TF empirical MPRA outcomes across broad unseen GSE188264 CRE edits (`r = 0.012258`); the radar reports all such negative findings verbatim.

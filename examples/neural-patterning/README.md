# Sequence to neural cell state

This experiment connects a real regulatory DNA substitution to a published
neural-development model. It reads the selected nucleotides from a validated
Development Module, derives a measured GLI binding-profile ratio, and uses that
ratio in a declared thermodynamic model. It has an offline interactive viewer,
complete intervention trajectories, source receipts and falsification controls.

The use case is an enhancer-variant mechanism experiment: determine whether a
proposed binding change can alter neural regulatory dynamics, identify the
environments where its effect reverses, and choose interventions that distinguish
binding from competing explanations. The resulting predictions are conditional
on the assay-to-cell transfer. They are candidates for biological testing.

## Run

Python 3.11 or later is sufficient. No packages or network access are required.
Run from the repository root and select a new output directory:

```sh
python3.11 examples/neural-patterning/patterning.py demo --output /tmp/neural-patterning
```

Open `/tmp/neural-patterning/index.html` in a browser. Move the time slider or play
development. Switch between the hypothetical endogenous edit and a passive
reporter. Export the current spatial figure as SVG, save the page as PDF through
the print dialog, or use the JSON and CSV files directly.

The command exits nonzero if a mechanism check fails. Failed results remain in
the output directory. Existing directories are never replaced. Numerical
acceptance covers every spatial trajectory and every intervention, including
transients, at two output resolutions with adaptive integration.

## Scientific chain

1. Peterson et al. tested `TGGGTGGTC → TGGGTGTTC` in the Nkx2.2 −2 kb enhancer.
   The included FASTA files contain these nine-base cores, not complete enhancer
   constructs. In the Hallikas strand orientation they are
   `GACCACCCA → GAACACCCA`.
2. Hallikas et al. Table S1 gives measured single-substitution binding profiles.
   The GLI2 mutant/reference ratio is `0.04235549749669733`; the GLI3 ratio is
   `0.03529466123979041`. GLI2 is the primary choice and GLI3 is a sensitivity
   analysis. No parameter is fitted to the developmental result.
3. The model uses the four-factor equations and parameters in Herrera-Delgado
   et al. 2020, Appendix J. Pax6, Olig2, Nkx2.2 and Irx3 interact locally under
   `S = exp(-x / 0.15)`. Position, time and concentrations are model units.
4. A conditional binding intervention scales both activating and repressing
   Gli occupancy, following Cohen et al.'s competitive-site mechanism. It is
   applied to the Nkx2.2 regulatory context specified by the consumer.

For a target with published effective basal weight `w`, Gli coefficient `K`,
activation coefficient `c`, and cross-repression factor `D`, the probability of
the active promoter state is

```text
on  = w (1 + K) (1 + q c K S)
off = (1 + q K) D
occupancy = on / (on + off)
dX/dt = 2 occupancy - 2 X
```

`q` is the profile ratio transferred as an affinity multiplier. `GliA = S` and
`GliR = 1 − S`; the same `q` applies to both by assumption. For the reference
`q = 1`, the factor `1 + K` cancels and recovers the 2020 equations. Changing
only the activation term would omit the effect on competitive repression.

The resulting neutral signal is `S = 1/c = 0.1`. At a fixed endogenous state,
reducing affinity lowers occupancy above that signal and raises it below it.
An independent enumeration of promoter microstates tests the equation. The
input-only ablation must lose this sign reversal. That ablation is a numerical
control, not a biological GliR knockout.

## What the controls establish

Exact sequence reversion must restore the original compiled interpretation and
trajectory. A mediator rescue retains variant DNA but overrides its multiplier
with the reference value. A phenocopy retains reference DNA and imposes the
variant multiplier. The output records those overrides explicitly.

Reporter mode computes passive promoter occupancy on the unchanged reference
network. Endogenous mode changes the target equation and reruns all four genes.
Their separation prevents a reporter mutation from silently becoming an
endogenous genomic intervention. Signal-history controls show that identical
final environments can retain different regulatory states.

The spatial experiment starts from zero concentrations at each position. Each
position is an independent regulatory system; this is not a multicellular
mechanics or connectivity simulation. The display's expression-state threshold
of 0.5 is a visualization convention, not a measured cell-fate boundary.

RK4 step doubling controls local absolute error at `1e-10`. The output interval
is 0.1 model units, with smaller internal steps as needed. The complete
trajectories must agree within `1e-4` when the output interval is halved. Input
switches are integrated exactly; the occupancy reported at the switch timestamp
uses the preceding phase's signal (the left limit). Both endpoint state and
transient behavior are checked. Earlier fixed-step attempts failed transient
refinement and are not accepted as numerical references.

## Run another supported site

The consumer accepts the consensus or one measured, positive-weight
substitution from it. It rejects ambiguous bases, multiple substitutions and
zero-weight entries. A published zero does not establish zero physical affinity.
Both DNA strands and sites within multi-record FASTA files are supported.

Create a protocol file:

```json
{
  "phases": [[5, 1], [15, 0.05]],
  "initial_state": {"Pax6": 0.1, "Olig2": 0, "Nkx2.2": 0, "Irx3": 0.1}
}
```

Each pair is `[duration, Gli activator fraction]`. The fraction must be between
zero and one. The optional `initial_state` specifies all four genes, each between
zero and one in model-relative concentration units. These are supplied cellular
conditions, not raw fluorescence intensities or states inferred from the DNA.
Omitting the field starts every gene at zero. The report records the state used;
missing genes, extra genes, nonnumeric values and out-of-range values are rejected
before an output directory is created. Then run:

```sh
python3.11 examples/neural-patterning/patterning.py run \
  --fasta examples/neural-patterning/variant.fasta \
  --record nkx2_2_site --start 0 --strand -1 --protein GLI2 \
  --mode endogenous --protocol /tmp/protocol.json --output /tmp/my-site
```

`start` is a zero-based FASTA interval start. `strand` orients the selected
nine-base interval into the Hallikas consensus. The user supplies this site
annotation; the program does not infer enhancer identity from a nine-base motif.
This mode writes `trajectory.csv`, `report.json`, the compiled bundle and a file
hash manifest. Its source, selection, scientific data and interpreter identities
are checked before execution. A coherently resealed nucleotide substitution or
different interpreter is rejected against the external selection receipt and
the consumer's supported target.

## Biological comparison and remaining evidence

Peterson Figure 5D–E reports expression in 4/6 reference-reporter embryos and
0/13 variant-reporter embryos at E10.5. Those are expressing/transgenic counts,
not fluorescence intensities. The model has no measured reporter-detection
threshold, copy-number correction or embryo-level noise model, so those counts
are displayed as observations and never scored as reproduced predictions.
The four-gene model omits Foxa2 and cannot reproduce the paper's Foxa2 effector
experiment.

The binding-profile normalization remains a material uncertainty. The original
Cell methods equation is inconsistent at its consensus boundary. Later methods
describe an inverse-Kd normalization under which the within-column ratio cancels
the normalization, but the original processing pipeline was not fully recovered.
No absolute affinity, kinetic rate or expression fold change is inferred from
the table. Transferring the ratio to intracellular binding is explicit and
unvalidated. Complete enhancer sequences, chromatin, cofactors and measured
GliA/GliR concentrations are also missing.

The literal 2014 parameter caption and the 2018 supplement did not reproduce
the intended pattern in independent replays. This experiment uses the explicit
later 2020 parameters, not an inferred correction to those sources. Untuned
reference states from an independent transcription are included in `data.json`.
The working model establishes executable regulatory dynamics. It does not yet
establish natural-genome-to-brain development, neuronal connectivity or learning.

## Sources

- Hallikas O et al. (2006). Genome-wide prediction of mammalian enhancers based
  on analysis of transcription-factor binding affinity. Cell 124, 47–59.
  [DOI](https://doi.org/10.1016/j.cell.2005.10.042),
  [Table S1](https://ars.els-cdn.com/content/image/1-s2.0-S0092867405013954-mmc2.xls).
  The JSON retains all 72 GLI2/GLI3 core-profile values, original cell coordinates
  and the workbook SHA-256. No explicit dataset license was found; the source
  article is © 2006 Elsevier. The extracted numerical facts are attributed here.
- Peterson KA et al. (2012). Neural-specific Sox2 input and differential Gli-binding
  affinity provide context and positional information in Shh-directed neural
  patterning. Genes & Development 26, 2802–2816.
  [DOI](https://doi.org/10.1101/gad.207142.112).
- Cohen M et al. (2014). A theoretical framework for the regulation of Shh
  morphogen-controlled gene expression. Development 141, 3868–3878.
  [DOI](https://doi.org/10.1242/dev.112573). Competitive binding mechanism;
  the baseline parameters are taken from the later paper below.
- Herrera-Delgado E, Briscoe J, Sollich P (2020). Tractable nonlinear memory
  functions as a tool to capture and explain dynamical behaviors. Physical
  Review Research 2, 043069, Appendix J.
  [DOI](https://doi.org/10.1103/PhysRevResearch.2.043069). CC BY 4.0.
- Winklmayr M et al. (2010). Non-consensus GLI binding sites in Hedgehog target
  gene regulation. BMC Molecular Biology 11, 2.
  [DOI](https://doi.org/10.1186/1471-2199-11-2).
- Zhao L, Stoddard BL (2014). Rapid determination of homing endonuclease DNA
  binding specificity profile. Methods in Molecular Biology 1123, 127–134.
  [Primary methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC3969613/).
  Used only to inspect a later explicit affinity-normalization pipeline.

Run the focused behavioral tests with
`python3.11 -m unittest tests.test_neural_patterning -v`.

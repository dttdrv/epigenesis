# Epigenesis Studio

A local desktop application for working with DNA changes. Its workflows simulate a regulatory DNA binding site in a published four-gene neural-development model, scan an entire arbitrary enhancer sequence through Boltzmann GLI partition physics to expose cryptic cell-fate bifurcations, score edits with either the frozen Siraj Huber-linear checkpoint or the experimental Kreimer 0–72h neural-induction temporal predictor, and export receipts that anyone with the same installation can replay. The interface is in English and Bulgarian (the switch is in the sidebar).

## Simulate a DNA change

The Simulate mode runs the repository's own backend, [`examples/neural-patterning`](../examples/neural-patterning/README.md), live. A run does the following:

1. Takes a nine-base GLI binding site with one substitution from the Nkx2.2 reference. Supported substitutions have a positive measured binding ratio from Hallikas et al. (2006), GLI2 as the primary protein and GLI3 as a sensitivity check. The reference is included automatically in the comparison.
2. Compiles the site into a validated Development Module that stores the selected nucleotide bases. The simulator reads those bases back and derives the binding ratio from the measured table.
3. Applies the ratio, as an explicit and unvalidated assumption, to the Gli affinity of the Nkx2.2 promoter in the four-gene model of Herrera-Delgado et al. (2020): Pax6, Olig2, Nkx2.2 and Irx3.
4. Integrates the model under the signal history and starting state you set, for the reference, for the variant as an edit inside the network, and for the variant as a passive reporter. It also reverts the base and confirms that the reference trajectory returns exactly.

The page shows the conditions, the binding table, the trajectories (move the pointer to read all genes at one time), the reporter-versus-network comparison, the controls that passed and the six assumptions the result depends on. The simulation is conditional and mechanistic. It does not predict embryo outcomes, and only this site with one substitution is supported; a different genome is outside the supported path. Signal is the fraction of GLI in its activator form, and time and concentration are model units. A run is limited to 10 phases and 100 model units.

Simulate needs the repository checkout (the install below is editable, so it finds `examples/neural-patterning`). Without it the mode is disabled and the other modes still work.

## Enhancer Circuit mode (1.7)

The Enhancer Circuit tab connects arbitrary enhancers (9 to 500 bases) to the four-gene neural-tube ODE via [`examples/neural-circuit`](../examples/neural-circuit/README.md):

1. Computes a strand-symmetric Boltzmann partition ratio (`q = Z_variant / Z_reference`) over sliding 9-mer Hallikas GLI2/GLI3 binding profiles for every single-base substitution (with Laplace pseudocount `epsilon`, default `1e-4`), then classifies each edit as a cryptic bifurcation (`|Δreporter| < 0.02` with a fate flip), concordant switch, canalized buffered, or neutral.
2. Compiles the top cryptic switch (and exact sequence reversion) into `brainc` Development Modules, compares endogenous vs passive-reporter trajectories, checks the 21-point spatial Shh gradient P3 domain and pulse-chase competitive GliR derepression control, and displays Kreimer 0–72h and Siraj scores beside the circuit shift.

Requires the repository checkout; unavailable when installed standalone without `examples/`.

## Predictor modes

Two sequence-based predictor checkpoints are available in **Scan** and **Edits** via a model toggle:

- **Siraj Huber-linear (`epigenesis-siraj-huber-linear-v1`)**: scores SK-N-SH 24-hour log2-skew reporter effect (see [the predictor README](../examples/neural-predictor/README.md)). On 16,986 internal-tuning edits its family-weighted error is 3.02% below a 5-mer baseline, while on five external neural reporter studies correlation is 0.0490 and MSE is 9.52% above predicting no effect.
- **Kreimer 0–72h Neural Induction (`epigenesis-kreimer-motif-temporal-v1`, 1.7)**: outputs a 7-timepoint trajectory (`0, 3, 6, 12, 24, 48, 72` hours, plus its temporal mean) using 5-mer counts and 532 compressed JASPAR 2024 motifs fitted on Kreimer et al. (2022) GSE188264 train+validation families (`test MSE 0.252079` across 118 held-out families vs k5 `0.285185` and dinucleotide `0.319358`, see [`examples/neural-regulation`](../examples/neural-regulation/README.md)). Retrospective and conditional on MPRAnalyze preprocessing; not validated for endogenous loci or clinical use.

Studio uses the compiler's canonical JSON for receipts. Simulate and Enhancer Circuit compile their pairs into Development Modules; predictor Scan/Edits modes seal their outputs into receipts without passing through the Development Module pipeline.

## Install

Python 3.11 or later and [uv](https://docs.astral.sh/uv/). From the repository root:

```sh
uv tool install --editable studio --with-editable .
```

This puts `epigenesis-studio` on your `PATH`. The application depends on NumPy and on the compiler's canonical JSON module, which is why the compiler is installed alongside it.

## Use

```sh
epigenesis-studio                  # serves the app and opens the page; run it in any folder
epigenesis-studio data/            # looks in another folder
epigenesis-studio edits.csv        # opens one file straight away
epigenesis-studio locus.fa --port 9000 --no-browser
epigenesis-studio verify receipt.json
```

The command prints the address (`http://127.0.0.1:8765/`, or a free port when 8765 is taken, so several folders can run at once). The page lists the sequences (`.fa`, `.fasta`, `.fna`, `.fas`, `.seq`, `.txt`), edit tables (`.csv`, `.tsv`, `.json` with `reference` and `sequence` columns) and receipts found directly in that folder, recognised by content, and opens each in the right tab. **Open file** in the sidebar also imports a file from anywhere, and a file dropped anywhere on the window opens the same way. Only files in the list can be opened, so a request cannot reach outside the folder.

The server listens on 127.0.0.1 only and rejects requests whose Host header is not local. An explicit `--port` that is taken is an error; `--port 0` picks a free port.

| Mode | Input | Output |
|---|---|---|
| Simulate | A site substitution, a GLI protein, a signal protocol and an optional starting state | Trajectories for reference, variant and reversion, controls, assumptions, CSV, receipt |
| Enhancer Circuit | A FASTA or ACGT enhancer of 9 to 500 bases, GLI2/GLI3, pseudocount, signal phases | Cryptic bifurcation scatter, trajectory & spatial P3/GliR controls, brainc receipt |
| Scan | One FASTA record or plain sequence of 5 to 1,000 bases (ACGT only; several records are refused); Siraj or Kreimer 0–72h | Heatmap of all 3 changes at every position, temporal 0–72h profile (Kreimer), CSV, receipt |
| Edits | Up to 5,000 reference and edited pairs (CSV or JSON, columns `id`, `reference`, `sequence`, optional `family`; each sequence 5 to 2,000 bases, insertions and deletions allowed); Siraj or Kreimer 0–72h | Scores, 7-hour curves (Kreimer), panel of the requested size with one edit per family, CSV, receipt |
| Verify | A receipt file | Three checks: digest, identity of model and endpoint, replay of every stored result |

A single-substitution score in both Siraj and Kreimer modes depends only on the 18 bases either side of the change, so a scan scores each variant in a window of that size, with results identical to scoring whole sequences (tested to 1e-12). Edits mode scores the supplied pair, including insertions, deletions and multiple substitutions.

## Receipts

A receipt (`epigenesis-studio-receipt-v1`) holds the model name, its identity hash, the endpoint, the normalized request, the result and a SHA-256 over the RFC 8785 canonical JSON of all of it. Verification checks three things. The digest matches the content. The model name, identity hash and endpoint match this installation. Running the stored request again reproduces every field of the stored result, with numbers within 1e-9. A changed digest, a changed model or endpoint (even with the digest recomputed) and a changed result or input (even resealed) each fail the matching check. A receipt shows that a result came from this model and request. It says nothing about whether the prediction is right.

## Tests

```sh
$(uv tool dir)/epigenesis-studio/bin/python -m unittest discover -s studio/tests -v
```

The tests skip with a message when NumPy or the package is missing, and the Simulate / Enhancer Circuit tests skip without the repository checkout. They cover scan equals whole-sequence scoring (Siraj & Kreimer), the panel against the frozen model run and Kreimer checkpoint, cryptic bifurcations and brainc bundle verification in Enhancer Circuit, each receipt check on its own, input limits and refusals, the HTTP guards, the simulation against the backend's own command line and its published binding ratio, and that the English and Bulgarian interface texts have the same keys.

## Licenses

Apache-2.0, as the rest of the repository. The motif potentials derive from JASPAR 2024 (CC BY 4.0); the model files and their manifest are in `epigenesis_studio/model/`.

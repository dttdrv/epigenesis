# DNA-derived relative binding preferences

This example consumes a validated coding module, translates its nucleotide
tensor, checks an explicit protein-domain mapping, and predicts relative DNA
preferences with the published FamilyCode single-residue model. It uses the
original 414 homeodomain training motifs and the published BLOSUM62 table.
There is no additional runtime dependency.

The supplied natural example assembles the two annotated mouse **Nkx2-2**
coding exons from original NCBI genomic DNA. The reference translation matches
the independently retrieved 273-residue NP_035049.1 protein, including a codon
split across the exon junction. The P156L example is a **designed** single-base
substitution, not a measured Nkx2-2 mutant or validated neural intervention.

## Run

From the repository root, using Python 3.11 or newer and an absent output path:

```sh
python3.11 examples/coding-consequences/coding.py compile \
  --fasta examples/neural-binding/data/NC_000068.8_147025815_147028038.fasta \
  --contract examples/neural-binding/data/nkx2-2-p156l.json \
  --output /tmp/nkx2-2-p156l

python3.11 examples/neural-binding/binding.py /tmp/nkx2-2-p156l \
  --selection examples/neural-binding/data/nkx2-2-p156l-selection.json \
  --selection-sha256 7fec289fce11a847cdc6aee6c323dc929785355fc0d079830dfba3d69322267b
```

Use `nkx2-2-wt.json` and `nkx2-2-wt-selection.json` for the reference, with
selection receipt
`782bbc931d08141cec0a9ca9be2fa7fec0a037a0873a4edf7a277223ec76db25`.
Receipts are SHA-256 of the canonical JSON objects, not their formatted files.
The consumer revalidates the coding bundle; a stored `prediction.json` is never
an input. Reports identify the actual source, coding contract, model, domain
selection and implementation.

## Model and coordinates

The selection is a closed version-1 object binding the coding receipt, model
hash, reference alignment, peptide-position map and conditioning position.
`peptide_positions` contains zero-based translated-peptide indices, strictly
increasing, with `null` only for alignment gaps. Insertions absent from the
family alignment may be skipped explicitly. Every mapped residue except the
one being conditioned must match the selected reference alignment. The
conditioning position is **one-based within the alignment**; it must be mapped.
The example maps peptide positions 128–184 to the 57 alignment columns and
conditions column 28, corresponding to full-protein residue 156.

This reproduces the arithmetic of `makeSVDModel.singleMutation`, not the
generic iterative or two-position FamilyCode model. With one identical
residue feature in each of all three principal components, regressing a
component on its amino-acid group means gives that conditional mean. Reversing
the complete SVD therefore recovers the average of the training motifs after
each column is sum-normalized. The implementation computes this directly.
For an unseen amino acid it combines the group means using the published
`1 / (BLOSUM[a,a] - BLOSUM[a,observed])` weights. Groups are not weighted again
by their sample counts. The result is max-normalized and floored at 0.01.

`psam` has base rows A/C/G/T and the eight named `dna_positions`. Each column's
maximum is one; its sum need not be one. `centered_log_preference` is natural
log weight minus the column's mean log weight, with sign **−ΔΔG/RT**. The
report states whether interpolation was required and how many training rows
contained the conditioning residue. That count is not a confidence interval.

The single-residue approximation ignores protein-background epistasis: the
same conditioning residue gives the same prediction across backgrounds.
The model does not supply absolute affinity, protein abundance, transcription,
cell fate or neural connections. Separately normalized mutant and reference
motifs cannot determine an absolute change in Kd. The paper's tetrahedrally
regularized mutant statistic is different from a literal energy difference
and is not used here as a physical free-energy change.

## Evidence and limits

```sh
python3.11 tests/verify.py neural-binding
```

Tests cover independent numerical examples, interpolation, genuine joined
genomic translation, synonymous recoding, changed amino acids, domain maps,
source/model identities and externally retained receipts. A separate full
SVD/OLS calculation agrees across all 1,140 standard-residue/position queries
and 36,480 base preferences to within 1.9e-15. This numerical check uses an
independent NumPy calculation; original-R numerical parity is unverified.
The frozen experimental comparison below fails predictive acceptance.
Software checks do not establish biological validity; consumer reports keep
`biological_acceptance` false.

## Frozen mutant comparison

```sh
python3.11 examples/neural-binding/benchmark.py --output /tmp/familycode-benchmark.json
```

The output path must be absent. The evaluator verifies the frozen protocol,
mapping, experimental data, predictor and training-data identities before
running. It records every prediction, observation, residual, loss, baseline,
source correction and factor-level uncertainty calculation. No model is fit.

The [protocol](benchmark-data/protocol.json) was fixed before computing local
prediction errors; it was not a blinded preregistration. All 92 mutants and
316 released fitted matrices are used. Original clone peptide pairs support
66 primary mutants across 27 factors. The remaining 26 mutants have a mandatory
descriptive secondary analysis; their original clone sequences remain missing.
Fourteen primary mappings explicitly correct the original script's duplicate-
name row selection using original assayed clone peptides. This is not literal
original-script or original-R numerical parity.

Every available mutant and WT matrix contributes. The raw-preference controls
are the mean measured WT preference and the mean of all equally nearest
non-self training domains. The second endpoint compares the predicted mutation
effect with a zero-effect prediction in centered natural-log preference space.
Losses average over bases, mutant replicates, variants and finally factors;
each factor has equal weight. A paired factor bootstrap uses 20,000 resamples.

The primary result **fails all three required comparisons**. Positive values
below would favor the model; acceptance requires every lower interval bound
to exceed zero.

| Comparison | Baseline MSE minus model MSE | 95% factor-bootstrap interval |
|---|---:|---:|
| Raw preference versus measured WT | -0.05294 | [-0.07358, -0.03200] |
| Raw preference versus nearest non-self domain | -0.01226 | [-0.03615, 0.01160] |
| Mutation contrast versus zero effect | -0.11863 | [-0.18166, -0.06156] |

An independently written evaluator reproduces all observations, predictions,
residuals and aggregate statistics within 3.1e-16. Regression controls detect
discarded replicates, wrong log-averaging order, self-baseline leakage, lost
ties, variant pooling, relaxed acceptance and ignored input identities. These
checks validate this calculation, not the failed model's biological adequacy.

This comparison uses WT factors represented in training; it is not a new-factor
holdout. Released matrices summarize fitted PBM measurements, not absolute Kd
or independently established preparation counts. Assayed DNA haplotypes are
not supplied by the clone peptide table. The failure remains part of the
research record and does not change the frozen protocol.

The Nkx2-2 domain also matches human NKX2-2. Sequence identity alone does not
identify the species or construct of the original FamilyCode assay. The
original interval GenBank record is annotation evidence; its `REGION` header
is outside the compiler's supported GenBank profile. Original FASTA bytes are
admitted instead. The separately included NM record remains an mRNA oracle.

[Provenance](data/provenance.json) records source versions, hashes, conversion
and attribution. The original RDS inputs accompany their JSON export; two
independent readers agree exactly on the arrays. Duplicate gene names with
different sequences remain distinct rows. FamilyCode repository data retain
their MIT license; the FCpackage substitution table retains Artistic-2.0.
NCBI records retain their accessions and [data policy](https://www.ncbi.nlm.nih.gov/home/about/policies/).

Method: Liu et al. (2025), [doi:10.1093/nar/gkaf831](https://doi.org/10.1093/nar/gkaf831).

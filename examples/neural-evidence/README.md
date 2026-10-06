# Measured neural evidence

An offline workbench for a real experiment: compare gcy-22 expression in ASER
neurons after CHE-1 depletion, with a native promoter, a transferred CHE-1 binding
core, and a larger transferred CHE-1 regulatory fragment.

```sh
python3.11 examples/neural-evidence/evidence.py demo --output /tmp/neural-evidence
python3.11 examples/neural-evidence/evidence.py verify --output /tmp/neural-evidence
```

Open `/tmp/neural-evidence/index.html`. Python 3.11+ and its standard library are
sufficient. Generation and verification require no network. The destination must
be new; an existing report is never overwritten. Keep the output directory
intact to retain its source-download links. The HTML alone contains all plotted
measurements and can export CSV, JSON and SVG without a server. Print styles
support the browser's print/save-PDF dialog.

## The researcher’s task

Does moving the 12-base CHE-1 core reproduce the larger fragment's observed
expression retention after depletion? Inspect control and treated measurements,
select any dot to recover its exact workbook cell, compare donor geometry, and
check which sequence representations distinguish the constructs. Filter the
source table or export the observations and complete report for further analysis.

This provides a reproducible evidence comparison and helps specify the next
controlled experiment. It does not require a successful DNA-edit predictor.

## What the original data show

The source is [Traets et al. (2021), eLife 10:e66955, Figure 6F](https://doi.org/10.7554/eLife.66955).
The endpoint is gcy-22 mRNA counts measured by smFISH in ASER. The figure describes
early L3 animals after 24 hours of auxin; the analysis notebook specifies 1 mM.
All promoter constructs share the CHE-1::GFP::AID and TIR1 background. The native
promoter group is not unmodified N2.

| Promoter | Control: rows / sum / mean | Auxin: rows / sum / mean | Ratio of means |
|---|---|---|---|
| Core transfer, donor 3418 | 22 / 602 / 27.36 | 21 / 50 / 2.38 | 8.70% |
| Fragment transfer, donor 3415 | 12 / 316 / 26.33 | 30 / 524 / 17.47 | 66.33% |
| Native promoter | 14 / 358 / 25.57 | 10 / 30 / 3.00 | 11.73% |

All 109 supplied observations, including six zero counts and repeated labels,
are retained. Numeric blanks are omitted. Control and treated endpoints are
unpaired. Supplied row labels do not establish unique animals or independent
biological preparations, so the tool does not infer population confidence
intervals or hypothesis tests from them. Relative-scale dots divide each count
by that construct's control-group mean; they are not individual retention rates.

The donor/reference alignment reveals an important design distinction. The
short transfer replaces 12 bases with 12; its core runs opposite the target gene.
The larger donor replaces a 160-base interval with a 130-base CHE-1 fragment;
its core runs in the same direction as the target gene. Surrounding sequence,
orientation and spacing therefore vary together. These observations cannot
isolate a flank-only or cofactor mechanism. The reconstruction describes the
intended fully templated donor product, not verified final clone sequences.

The information comparison checks a finite-data constraint: if two feature
inputs are identical, any deterministic function of those inputs predicts a
single common value. Across the two observed auxin group means, with equal
weight, the smallest possible RMSE is half their absolute difference (7.54
molecules for the strand-normalized core). This is algebra, not a population
rejection or predictive-performance result. Distinct inputs remove that equality
constraint; they do not establish a working predictor.

## Sources and reconciliation

`case.json` declares the reviewed column selections, metadata, source hashes and
source URLs. Its identity is pinned in `evidence.py`. The script verifies source
bytes before parsing the original XLSX, DOCX and accession-version FASTA files.
It selects the same numeric columns as Figure 6 notebook cells 7–11, and uniquely
aligns the donor's two 35-base homology arms to both reference strands. Figure
6's original notebook and article XML are included for interpretation review.

The original paper has a caption conflict: Figure 6F prose names the reciprocal
che-1-promoter construct, whereas the main text, figure labels, notebook and
donor table support gcy-22 core-transfer allele gj2064. The workbook also has a
conflicting core label and a native-control L3/L4 label. Methods and Results
differ on ethanol concentration. These discrepancies remain visible in the
viewer and export. Our reconciliation is not an author-issued correction.

Source data and article: Traets et al., CC BY 4.0, as declared in the article.
NCBI reference accessions: NC_003283.11 and NC_003279.8; downloaded interval
coordinates and original retrieval URLs are recorded in `case.json`.
`data/figure6f.xlsx` and `data/figure6.ipynb` retain the original bytes of the
named Figure 6 source-archive members; other source files are likewise unedited.

`verify` regenerates the outputs and compares bytes using the installed source
and reviewed inputs. It is a reproducibility check, not an independent scientific
validator. Tests independently fix the original counts, sums, expected donor
coordinates and analytical fitting bound, and reject changed sources, invalid
counts, ambiguous alignments and formula cells.

```sh
python3.11 -m unittest discover -s tests -p test_neural_evidence.py -v
```

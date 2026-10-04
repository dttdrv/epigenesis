# Coding consequences from natural DNA

This consumer extracts explicitly selected coding segments from admitted
FASTA DNA, reconstructs a reference-bound local edit, compiles the resulting
nucleotides into a Development Module, validates it, and translates its loaded
tensor using the pinned NCBI standard genetic code. It predicts a peptide and
the first in-frame stop. The compiler itself still accepts a caller-supplied
interpretation.

From the repository root, with Python 3.11 or newer:

```sh
python3.11 examples/coding-consequences/coding.py compile \
  --fasta examples/coding-consequences/data/NC_007133.7_37347793_37350078.fasta \
  --contract examples/coding-consequences/data/sox2-minus4.json \
  --output /tmp/sox2-minus4
```

The output directory must not exist. The report includes `contract_sha256`.
Keep that receipt outside the candidate directory and supply it when replaying:

```sh
python3.11 examples/coding-consequences/coding.py consume /tmp/sox2-minus4 \
  --contract-sha256 RECEIPT_FROM_THE_COMPILE_REPORT
```

Use `sox2-wt.json` for the reference allele. The reference predicts 315 amino
acids; the displayed four-base deletion predicts 97, comprising 69 unchanged
residues and `NSGPSGNFCPRARSDHSSTKPNAFGLCT`, followed by a stop at codon 98.
These agree with the *predicted* mutant schematic in Gong et al. (2020),
[Supplementary Figure S1](data/sox2-pmc-image1.tiff). No mutant-protein abundance
measurement is asserted. Ordinary codon translation is established science;
this example makes its source, assumptions and execution auditable.

## Contract and execution

`contract.json` has a closed version-1 schema, exemplified by the supplied
contracts. `source_sha256` binds the exact original FASTA bytes. `record_id`
selects one record. `start` and `end` are zero-based, half-open coordinates in
that record's supplied strand; `strand` is `1` or `-1`. The consumer supports
only a complete reference CDS, unambiguous coding DNA, genetic code
`1` and canonical `ATG` initiation. The reference's first in-frame stop must be
its final codon. It does not infer gene locations, splice forms,
alternative initiation, RNA editing or translation efficiency.

Version 2 replaces `start` and `end` with `segments`, a nonempty list of
`{"start": integer, "end": integer}` objects in ascending genomic order.
Intervals are zero-based, half-open and nonoverlapping. Their bases are joined
before applying `strand` and translating, so a codon may span two exons.
This executes an explicit annotation; it does not predict splice sites.
Version 1 retains its contiguous-CDS schema.

An optional `edit` specifies a zero-based `start` in the **oriented reference
coding interval**, exact `deleted` bases and `inserted` bases. Reconstruction
checks the deleted sequence, applies the edit and preserves the selected
record's other bases. Reverse-strand edits are reverse-complemented back into
the genomic record. The derived FASTA contains only that selected record.
Unedited admission preserves the complete original FASTA bytes and records.
For version 2, edit coordinates refer to the oriented, joined CDS. A deleted
span must lie within one exon. Insertion exactly at an internal exon junction
is rejected, because its genomic location is ambiguous. Introns are preserved;
the consumer does not reinterpret a cross-junction transcript deletion as a
contiguous genomic deletion. Equivalent supported edit placements produce the same derived sequence but retain their
different provenance contracts.

The example accepts at most 1 MiB of source bytes and 10,000 template bases.
An edited template need not have a length divisible by three. Translation
stops at the first complete in-frame stop; `stop_offset` is zero-based in the
edited template. If no stop occurs, the status is `no-stop-in-template` and
`trailing_bases` reports any remaining one or two bases. This is a bounded
partial prediction, not a claim that the protein ends at the template boundary.
An edit that removes canonical initiation is rejected. No downstream genomic
flank is silently treated as transcribed sequence.

The Development Module has one unit, an ASCII-nucleotide `u64` tensor and a
translation rule pinned to this implementation, the code authority and the
exact caller contract. It contains no precomputed peptide or phenotype.
Consumption validates the source and module from the same in-memory snapshots,
reconstructs the complete template, checks every linked nucleotide, then
translates the loaded tensor. The external receipt detects replacement by an
internally consistent but different source or annotation contract. A receipt
does not establish that a caller's annotation is biologically correct.

## Source evidence and limits

The supplied genomic intervals are original NCBI FASTA responses:

| Reference | Assembly / strain | Local CDS, 1-based inclusive |
|---|---|---|
| NC_007133.7:37347793..37350078 | GRCz11 / Tuebingen | complement(1016..1963) |
| NC_141042.1:40458511..40460773 | GRCz13ab / M-AB | complement(993..1940) |

The FASTA sequence exactly matches each original GenBank ORIGIN. These GenBank
interval responses have `ACCESSION ... REGION:` headers outside the current
closed admission profile and fail with `INSDC026`; their original bytes are
retained as annotation evidence. The FASTA frontend imports no annotation.
`NM_213118.1.gb` is mRNA and supplies a separately identified protein oracle;
it is not admitted or relabeled as genomic DNA.

The four sequenced local alleles in Figure S1C are reconstructed on these
references. The original experimental clone's complete haplotype was not
sequenced here: unchanged flanks are assumed reference. The two references
encode the same normal protein, yet their identical 11-base deletions produce
241-residue mutant peptides differing at residues 116 and 230. This detects
shortcuts based on gene names, mutation sizes or protein lengths.

The in-frame three-base deletion and insertion are **not biological neutral
controls**: the paper describes similar phenotypes. Rescue is partial and
stage/context dependent. The model supplies no measured relationship from
peptide sequence to Sox2 activity, neural differentiation, neurite growth,
connections or behavior. `biological_acceptance` remains false. The broader
DNA-to-development objective is still open.

## Reproduction and attribution

```sh
python3.11 tests/verify.py coding-consequences
```

Tests cover all 64 codons, all three stops, partial templates, reverse strands,
flanks, alternate records, equivalent edits, same-size frameshifts, natural
reference differences, and coherent source/contract/tensor/rule substitutions.
The altered-tensor cases remain compiler-valid before consumer rejection,
including a change after the first stop that leaves the peptide unchanged.

[provenance.json](data/provenance.json) records exact sizes, hashes, URLs and
attribution. The unchanged Figure S1 TIFF was extracted from
`word/media/image1.tiff` of the original supplementary DOCX; that container's
hash and recovery route are recorded. The full article XML and extracted
supplement captions are included. These sources retain their original terms,
separately from the repository's Apache-2.0 code license.

Gong J et al. (2020), *The Requirement of Sox2 for the Spinal Cord Motor Neuron
Development of Zebrafish*, Front. Mol. Neurosci. 13:34,
[doi:10.3389/fnmol.2020.00034](https://doi.org/10.3389/fnmol.2020.00034).
Copyright © 2020 Gong, Hu, Huang, Hu, Wang, Zhao, Qian, Wang, Sheng, Lu, Wei and
Liu. The recovered article declares CC BY; no version is specified in its XML.
The NCBI sequence annotations and original
[genetic-code table](https://ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt) retain
their accession/version and original comments; see
[NCBI code documentation](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
and [NCBI data policy](https://www.ncbi.nlm.nih.gov/home/about/policies/).

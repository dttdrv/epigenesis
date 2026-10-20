# CRX binding-to-regulation bridge audit

Result, 2026-10-04: the recovered experiment does not identify the proposed
binding-to-regulation mechanism. No reporter model was fitted and no prediction
performance was calculated. This is a source and identifiability audit, not a
successful developmental model.

## Sources and genomic reconstruction

The original studies are Zheng et al., [eLife 87147](https://doi.org/10.7554/eLife.87147)
and [Genome Research, gr.279340.124](https://doi.org/10.1101/gr.279340.124).
Original reporter inputs are pinned to CRXHD_mpra commit
`cd0b681e6607454bcf611ac8c664647380342e16`; the mouse-model source is pinned to
`f68ab116e34ed29c4f9372fdf0a12f8e2b302f5c`. Eight original Spec-seq tables from
[GSE223658](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223658)
retain two lanes for WT, E80A, K88N and R90W.

The exact mouse genomic interval `NC_000073.6:15865847..15880055` contains three
reverse-strand coding exons. Their 900-base joined CDS produces the independently
retrieved 299-residue `NP_031796.1` sequence. Both original 190-base donors align
uniquely with one change. The E80A donor specifies `c.239A>C`, converting GAG to
GCG; the article's `c.239A>G` wording would instead produce E80G. K88N agrees
with `c.264G>T`, converting AAG to AAT. The discrepancy remains explicit.
Unchanged genomic flanks are assumed reference, not measured mutant haplotypes.

The existing coding consumer compiles and validates WT, E80A and K88N bundles.
Local exact sources, contracts, independent audit and replay are retained in
`build/neural-crx-genomic-20261004/`. This validates the coding reconstruction,
not protein activity, gene regulation or development.

## Coverage and independent units

Each binding table measures all 64 `TAANNN` and 64 `NNNTTA` rows: 127 distinct
literal words, with the palindromic `TAATTA` measurement duplicated after its
count correction. All 2,048 ratio/log calculations replay within 8.9e-16.
Within-lane reference ratios connect the two library gauges; they must not be
equated to one. Absolute affinity, active protein concentration and transfer
from fixed assay flanks to retinal CRE context remain unidentified.

The processed reporter matrix contains 4,368 constructs and 22 RNA preparation
columns across six genotypes. Each preparation pools three retinas. There are
233 entirely missing rows, including three scrambled controls; 415 observed
numeric zeros remain distinct from missing data. Shared preparation columns,
sequence backgrounds and normalization do not constitute independent animals.

An exhaustive sequence audit over 1,743 genomic backgrounds finds zero
nonidentical cis pairs with complete measured coverage of their changed binding
windows. A weaker condition, requiring supported intended monomer changes and
no measured-domain entry or exit, leaves 47 comparisons across 34 backgrounds;
41 comparisons across 32 backgrounds have reporter observations. Every such
pair still changes unmeasured words. Selection used sequence and availability,
not observed response magnitude.

## Why the proposed fit was stopped

An independent check of the actual CRE strings finds 37 single and nine double
palindromic-site losses, plus one loss with a reverse-complement orientation
flip. The measured transitions are `TAATTA` to `TAAGTA` or its reverse complement
`TACTTA`. For any strand-symmetric site-response function, the predicted change
is therefore the site-loss count times one fixed difference. A fitted gain per
genotype makes it exactly equivalent to a site-loss-count baseline, including
when the site response is saturating occupancy. Concentration and gain cannot
be separated with those contrasts.

Literal-orientation weighting adds the assay's fixed-flank asymmetry; it needs
an oriented-word-count baseline and cannot establish a binding mechanism by
beating an ordinary changed-base-count baseline. A constant unknown contribution
from unsupported windows cancels in the additive contrast, so it cannot be
calibrated either. Independent unknown effects per changed word could explain
arbitrary outcomes. No fit was run to disguise these limits.

Original sources retain separate terms: the 2023 article declares CC BY 4.0,
the 2025 article CC BY-NC 4.0, and the MPRA repository MIT. Local source archives
retain their notices. Source discrepancies in energy units, reporter ratio versus
difference terminology, missingness and published cohort filtering remain in
the preserved audit. The local evidence archive is
`build/neural-crx-audit-evidence-20261004/`; it is not a model-performance result.

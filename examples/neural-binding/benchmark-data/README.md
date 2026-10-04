# Frozen experimental inputs

`protocol.json` and `variant-mapping.json` were frozen and independently
reviewed on 2026-10-04 before local prediction/performance calculation.
Their source accounting and mapping checks contain no fitted prediction
results. `benchmark.py` pins all four JSON files by SHA-256.

`kock-panel.json` preserves the numeric contents and original SHA-256 identities
of all 316 released ProBound log-weight matrices: 222 mutant matrices and 94
WT matrices. These are fitted PBM summaries, not raw intensity measurements.
`input-manifest.json` lists 335 original or derived source inputs and their
conversion dependencies. The complete original freeze and verification files
are retained locally under `build/neural-binding-benchmark-freeze-20261004/`;
the source distribution contains the four JSON inputs needed to reproduce
the calculation, not every upstream file in that manifest.

Source: Liu et al., *Nucleic Acids Research* 53 (2025), gkaf831,
[doi:10.1093/nar/gkaf831](https://doi.org/10.1093/nar/gkaf831).
The original FamilyCode repository is pinned at
[`35726b009d7475b6783936530689a04c743a1900`](https://github.com/BussemakerLab/FamilyCode/tree/35726b009d7475b6783936530689a04c743a1900)
and retains its [MIT license](../data/FamilyCode-LICENSE).

Original experiment and clone supplement: Kock et al., *Nature Communications*
15 (2024), 3110, [doi:10.1038/s41467-024-47396-0](https://doi.org/10.1038/s41467-024-47396-0),
under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
These JSON files are local extractions and annotations: array conversion,
clone-to-domain mapping, provenance classifications and explicit source
discrepancies. They do not change upstream authorship or imply endorsement.

The published single-residue method and substitution table retain their
separate FCpackage provenance and Artistic-2.0 attribution in `../data/`.

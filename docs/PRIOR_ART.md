# Prior art and scientific position

Epigenesis treats biological records as source programs. An explicit source
profile defines their grammar and native semantics; a declared interpretation
provider supplies biological meaning; the compiler checks and lowers both into
a typed developmental IR; an exact linker closes the referenced artifacts; and
a separately implemented validator replays the compilation from the original
bytes.

The closest direct precedent is Cai et al. (2009). Their compiler lexed and
parsed synthetic DNA, evaluated an attribute grammar, emitted
molecular-interaction models, and simulated every construct in a small
genetic-parts design space. The paper explicitly compared DNA-to-model
translation with source-to-object-code compilation.
[Cai et al., 2009](https://doi.org/10.1371/journal.pcbi.1000529)

Epigenesis's contribution is a broader compiler architecture:

> Multiple explicit biological source languages and content-identified
> interpretation providers lower through one typed developmental-IR contract,
> producing exact reference-only bundles whose complete compilation can be
> replayed by an independent implementation.

A primary-source search through 2026-09-01 found no published compiler combining
that source/profile boundary, external interpretation-provider boundary, typed
developmental target, exact bundle/link closure, and independent compilation
replay in one system. The novelty claim is this combination, not the isolated
idea of compiling DNA or generating biological models.

## Closest systems

| Work | Compiled direction | Shared ground | Epigenesis contribution |
|---|---|---|---|
| [Cai attribute-grammar compiler](https://doi.org/10.1371/journal.pcbi.1000529) | Characterized synthetic DNA → molecular-interaction and SBML models | Direct DNA parsing, semantic actions, multipass compilation | Multiple real biological source profiles, provider-neutral interpretation, typed developmental IR, exact linking, and independent replay |
| [ModelSEED](https://doi.org/10.1038/nbt.1672), [CarveMe](https://doi.org/10.1093/nar/gky537), [gapseq](https://doi.org/10.1186/s13059-021-02295-1) | Microbial genome evidence → genome-scale metabolic model | Automated genome-conditioned model construction | The interpretation knowledge is an explicit, versioned provider input to a general compiler contract |
| [Karr whole-cell model](https://doi.org/10.1016/j.cell.2012.05.044) | Curated genotype and cellular evidence → whole-cell model | Genotype-linked executable model construction | Reusable source, provider, target, bundle, and replay interfaces rather than one organism-specific reconstruction |
| [BioCRNpyler](https://doi.org/10.1371/journal.pcbi.1009987) | Components, mechanisms, and context → chemical-reaction network and SBML | Biological model compilation with caller-supplied semantics | Original biological records, source maps, typed developmental lowering, and independently replayable linking |
| [Proto BioCompiler](https://doi.org/10.1371/journal.pone.0022490) | High-level spatial programs → genetic regulatory networks and models | Biological compiler pipeline and model targets | DNA/genome records are source languages rather than the generated implementation medium |
| [SBOL](https://sbolstandard.org/datamodel-specification/version-3.1.0/), [SBML](https://sbml.org/documents/specifications/level-3/version-2/), [CellML](https://www.cellml.org/cellml/2.0), [NeuroML](https://pmc.ncbi.nlm.nih.gov/articles/PMC11723582/), [SONATA](https://pmc.ncbi.nlm.nih.gov/articles/PMC7058350/) | Authored biological or neural descriptions → portable machine-readable models | Typed structures, dynamics, networks, and interoperable artifacts | Compilation from exact biological source and declared interpretation into one closed developmental target contract |
| [libRoadRunner](https://pmc.ncbi.nlm.nih.gov/articles/PMC4607739/), [NESTML](https://pmc.ncbi.nlm.nih.gov/articles/PMC12174165/), [RateML](https://pmc.ncbi.nlm.nih.gov/articles/PMC10013028/) | Model language → native or simulator-specific executable | Model compilation and code generation | The source-side biological language, interpretation, and artifact provenance remain bound through lowering |
| [GEC](https://pmc.ncbi.nlm.nih.gov/articles/PMC2843955/), [Eugene](https://doi.org/10.1371/journal.pone.0018882), [Cello](https://doi.org/10.1126/science.aac7341), [Genotype Specification Language](https://doi.org/10.1021/acssynbio.5b00194) | Logical program or design constraints → DNA construct | Compiler architecture, constraints, and characterized biological libraries | Opposite direction: Epigenesis starts from biological source and compiles toward a computer-native IR |

## Compiler boundary

Epigenesis separates five concerns that are often fused in genome-to-model
pipelines:

1. **Source profile.** Every accepted language names its grammar, version,
   authority, coordinate rules, reference semantics, and limits. Dispatch is
   explicit and independent of accessions, organisms, genes, fixture digests,
   and filenames.
2. **Source representation.** Parsing preserves exact source identity, source
   maps, normalized sequence and annotation structure, ordering, and unresolved
   references without silently inventing biological meaning.
3. **Interpretation provider.** A grammar, rule set, annotation database,
   parameter set, or analysis result is a caller-supplied artifact with a type,
   version, digest, configuration, and declared output.
4. **Typed lowering and linking.** Provider output lowers through one closed
   developmental IR. Units, edges, rules, and ports are type-checked; referenced
   artifacts are linked by digest; the final bundle proves exact transitive
   closure without duplicating large children.
5. **Independent replay.** A separate validator reparses the original bytes,
   rechecks the profile and provider bindings, recomputes lowering, verifies
   canonical identities and reference closure, and compares the complete
   compilation result without importing producer modules.

This division makes biological interpretation replaceable while keeping the
compiler deterministic. Different interpretation packages can be compared
against the same source and target ABI; a change in source, package,
configuration, or lowering produces a correspondingly different identity.

## Meaning of “universal”

Universal describes the architecture: every supported source language enters
through an explicit profile and every interpretation enters through the same
provider contract before lowering to the shared typed target. Adding a format or
scientific model extends a declared boundary instead of adding content-specific
dispatch or a parallel compiler.

The result is a reproducible translation layer between biological records and
computer-native developmental descriptions: one stable, inspectable artifact
with an exact source-to-target proof chain.

# Prior art and scientific position

Epigenesis treats biological records as source programs. Each built-in source
profile owns a finite grammar and native IR; the compiler then checks a declared
interpretation, target contract, and lowering policy before linking a typed
Development Module bundle. Producer-isolated validation reparses built-in
source bytes. External profiles instead supply a manifest-qualified caller
attestation, and Epigenesis verifies its exact manifest, original/native bytes,
sequence catalog, and digest closure.

The closest direct precedent is Cai et al. (2009). Their compiler lexed and
parsed synthetic DNA, evaluated an attribute grammar, emitted a molecular-
interaction representation, and compared the process with source-to-object-code
compilation. [Cai et al., 2009](https://doi.org/10.1371/journal.pcbi.1000529)

Epigenesis's contribution is the combination of:

> Multiple explicit DNA-bearing source profiles, content-identified
> interpretation inputs, one typed Development Module contract, exact
> reference-only linking, producer-isolated built-in reparse, and
> manifest-qualified external attestation verification.

In the primary works reviewed below, no system combines those boundaries in one
compiler. This is a scoped comparison of the cited corpus, not an exhaustive
claim over all published software.

## Closest compiler systems

| Work | Compiled direction | Relationship to Epigenesis |
|---|---|---|
| [Cai attribute-grammar compiler](https://doi.org/10.1371/journal.pcbi.1000529) | Characterized synthetic DNA to molecular-interaction representations | Direct DNA-as-source precedent using parsing, semantic actions, and multipass compilation |
| [Proto BioCompiler](https://doi.org/10.1371/journal.pone.0022490) | High-level spatial programs to genetic regulatory networks | Biological compiler pipeline, but its source is a programming language rather than exact biological records |
| [Eugene](https://doi.org/10.1371/journal.pone.0018882), [Cello](https://doi.org/10.1126/science.aac7341), and [Genotype Specification Language](https://doi.org/10.1021/acssynbio.5b00194) | Design languages and constraints to DNA constructs | Design-to-DNA compilation; Epigenesis is DNA-to-IR |
| [DNA Chisel](https://academic.oup.com/bioinformatics/article/36/16/4508/5869515) | DNA plus constraints and objectives to optimized DNA | Sequence optimization and constraint validation rather than source-to-target IR lowering |

## Compiler boundary

Epigenesis separates five compiler concerns:

1. **Source profile.** Dispatch is explicit and independent of accessions,
   organisms, genes, fixture digests, and filenames. Built-ins declare local
   grammar versions; an external manifest binds its grammar authority by URI and
   digest.
2. **Native representation.** Parsing binds exact source identity, ordering,
   normalized sequence or annotation structure, and source maps where the
   selected profile emits them.
3. **Executable external grammar boundary.** Epigenesis runs an explicitly
   selected, digest-pinned frontend and a separate validator over stable source
   snapshots, then preserves their native JSON and common sequence catalog.
4. **Typed lowering and linking.** Interpretation output lowers through one
   closed Development Module IR. Units, edges, rules, and ports are type-checked;
   the final bundle binds exact transitive reference closure without duplicating
   large children.
5. **Producer-isolated validation.** A separate implementation reparses
   built-ins, reruns version-2 external validators, recomputes lowering, checks
   canonical identities, and compares the complete compilation without
   importing producer modules.

This separation keeps source parsing, interpretation, target semantics, and
linking independently replaceable while preserving deterministic identities.

## Source-profile scope

The built-in profiles are enumerated in the compiler contract. Another
DNA-bearing grammar enters through a versioned profile when its frontend and
validator can produce and independently replay a closed native IR, a nonempty
sequence-record catalog, exact source and native identities, and finite limits.

The result is a reproducible compiler boundary between biological source
records and typed Development Module IR, published as a stable inspectable
bundle with a source-to-target evidence chain.

## Neural experiments and measurement inference

The experimental platform builds on established biological models and simulation
methods. Its contribution is the executable connection to exact source records,
the reconstruction of experimental inputs, and the accompanying intervention
and verification procedures.

| Work | Established result | Epigenesis implementation |
|---|---|---|
| [Zubler et al., 2011](https://doi.org/10.3389/fncom.2011.00057) | Gene-like instructions copied from a single progenitor guide replication, differentiation and neural self-construction in CX3D. | A separate artificial recipe with source provenance, geometry traces and transmission/lesion controls using the CX3D substrate. |
| [Cohen et al., 2014](https://doi.org/10.1242/dev.112573) | Competitive Gli activation and repression produce a neutral point: changing binding affinity has opposite effects on its two sides. | A specified Nkx2.2 substitution, measured binding profiles and published four-gene dynamics connected through an explicit affinity-transfer assumption. |
| [Traets et al., 2021](https://doi.org/10.7554/eLife.66955) | A larger CHE-1 promoter fragment transfers expression resilience that the binding core alone does not. | Reconstruction of 109 observations and intended donor products, with cell-level source links and explicit orientation, extent and spacing differences. |
| [Majoros et al., 2020; online 2019](https://doi.org/10.1093/bioinformatics/btz545) | BIRD models regulatory effects and variation between assay repeats. | A separately implemented and numerically checked forecast comparison on four fixed neuronal-assay TRAIN components. |
| [Bergmann et al., 2014](https://doi.org/10.1186/s12859-014-0369-z) | COMBINE archives package biological models, experiments, data and metadata for reproduction. | A narrower compiler contract that checks exact source identity, explicit interpretation, typed lowering and independent replay. |

Neural self-construction, the Gli affinity sign reversal and replicate-aware
assay inference therefore have direct precedents. The current release offers
inspectable reconstructions and a small new measurement evaluation. Establishing
a new biological mechanism requires a distinguishing intervention and independent
experimental measurements. The [presentation guide](../research/NEUROSCIENCE_PRESENTATION.md)
sets out a concrete experimental-design use case and the current evidence.

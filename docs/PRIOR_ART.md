# Prior art and scientific position

Epigenesis treats biological records as source programs. Each built-in source
profile owns a finite grammar and native IR; the compiler then checks a declared
interpretation, target contract, and lowering policy before linking a typed
Development Module bundle. Independent validation reparses built-in source
bytes. External profiles instead supply a manifest-qualified validation report,
and Epigenesis independently verifies its exact manifest, original/native bytes,
sequence catalog, and digest closure.

The closest direct precedent is Cai et al. (2009). Their compiler lexed and
parsed synthetic DNA, evaluated an attribute grammar, emitted a molecular-
interaction representation, and compared the process with source-to-object-code
compilation. [Cai et al., 2009](https://doi.org/10.1371/journal.pcbi.1000529)

Epigenesis's contribution is the combination of:

> Multiple explicit DNA-bearing source profiles, content-identified
> interpretation inputs, one typed Development Module contract, exact
> reference-only linking, independent built-in reparse, and manifest-qualified
> external evidence verification.

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
3. **External evidence boundary.** An external frontend runs outside Epigenesis
   and supplies a validation report. The core verifies the report and exact byte
   closure without loading frontend code or executing the declared grammar.
4. **Typed lowering and linking.** Interpretation output lowers through one
   closed Development Module IR. Units, edges, rules, and ports are type-checked;
   the final bundle binds exact transitive reference closure without duplicating
   large children.
5. **Independent validation.** A separate implementation reparses built-ins,
   verifies external closures, recomputes lowering, checks canonical identities,
   and compares the complete compilation without importing producer modules.

This separation keeps source parsing, interpretation, target semantics, and
linking independently replaceable while preserving deterministic identities.

## Meaning of “universal”

Universal describes the architecture, not a claim that every biological format
ships in 1.0. The built-in profiles are enumerated in the compiler contract.
Another DNA-bearing grammar enters through a versioned profile only when it can
produce a closed native IR, a nonempty sequence-record catalog, exact source and
native identities, finite limits, and manifest-qualified validation evidence.

The result is a reproducible compiler boundary between biological source
records and typed Development Module IR, published as a stable inspectable
bundle with a source-to-target evidence chain.

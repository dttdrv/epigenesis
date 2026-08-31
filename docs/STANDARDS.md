# Standards and prior-art trace

Epigenesis 0.9.0 treats each biological format as a source language with its own
authority, grammar, identity rules, coordinate system, and versioned IR. Its
GenBank 0.8 frontend accepts sequence and annotation together from a native
archival record; the package's feature-state adapter lowers that BioIR to
deterministic feature state.

## GenBank 0.8 authority chain

| Layer | Primary authority | Release used | Compiler treatment |
|---|---|---|---|
| GenBank flat file | [NCBI GenBank Release Notes](https://www.ncbi.nlm.nih.gov/genbank/release/current) | 273.0, 15 August 2026 | Traditional physical-DNA record envelope, fixed fields, record order, FEATURES, ORIGIN, and terminators |
| Feature table | [DDBJ/ENA/GenBank Feature Table Definition](https://www.insdc.org/submitting-standards/feature-table/) | 11.4, April 2026 | Feature keys and qualifiers as structural values; location AST and normalized segments; source-coverage semantics |
| Sequence identity | [GA4GH refget Sequences](https://ga4gh.github.io/refget/sequences/) | 2.0.0 | Uppercase IUPAC sequence identity using the `SQ.` SHA-512/24URL identifier |
| Collection identity | [GA4GH refget Sequence Collections](https://ga4gh.github.io/refget/seqcols/) | 1.0.0 | Ordered names and sequence IDs as inherent attributes; lengths retained as a collated non-inherent attribute |
| Canonical JSON | [RFC 8785](https://www.rfc-editor.org/rfc/rfc8785) | May 2020 | Canonical bytes for every artifact, BioIR, record, and report seal |
| Hashing | [FIPS 180-4](https://csrc.nist.gov/pubs/fips/180-4/upd1/final) | August 2015 | SHA-256 artifact and source identities; SHA-512 basis for refget identifiers |

The packaged authority manifest records the exact retrieved bytes and digests:

| Authority object | Byte length | SHA-256 |
|---|---:|---|
| NCBI `gbrel.txt` 273.0 | 1,235,842 | `9e26f1821a1846aba58d9b8bee167e89355c0711995cd2f21b716a6123e5bda0` |
| INSDC page API response | 206,717 | `724b3ca19099eb060830b05566f9c619dcd4d1da4abdeccdf24bae5a1b88a1a0` |
| Selected INSDC rendered definition | 196,690 | `40f7fbe1f778f5a2a45ede9567b5df97fd5eb5beb2cf2984c560527ae78940c3` |

The manifest itself is sealed as
`6dc3658a5f7d6a774ce9cd08ff6dc5327cbef8a9951c7679fd9a34650f912bf9`.
Compilation and independent validation both verify that pin before accepting the
profile.

## Implemented standard mapping

### GenBank record structure

The frontend follows the traditional GenBank entry organization described in
Release 273.0: keyword fields begin in the fixed header area, feature keys begin
in the FEATURES key field, sequence indices end in column 9, and `//` terminates
an entry. The compiler binds mandatory and optional blocks by grammar, including
wrapped ORGANISM text, optional singleton `SEGMENT n of m` between KEYWORDS and
SOURCE, consecutively numbered references, JOURNAL subkeywords, and optional
legacy VERSION GI values.

Physical DNA is authoritative in ORIGIN. LOCUS length, VERSION identity, record
order, source features, ORIGIN indices, IUPAC sequence, source maps, and emitted
collection members are checked as one linked compilation unit.

### INSDC feature locations

Feature Table 11.4 supplies the location model used by the compiler. The v0.8
profile represents points, closed spans, endpoint fuzz, between-base sites,
archived within positions, exact-version remote references, complements, joins,
and orders. It preserves the source expression as an AST and emits ordered
0-based half-open segments.

Same-input remote locations are resolved against explicit record lengths.
External remote locations remain named and carry unresolved bounds status. This
keeps the compiled output deterministic and makes later reference resolution an
explicit linking pass.

Feature keys and qualifier names and values are preserved structurally.

### Refget identity

Sequence Collection v2 retains three complementary identities:

1. SHA-256 of the exact original GenBank bytes;
2. SHA-256 and refget Sequence v2 identity of each normalized sequence; and
3. refget Sequence Collections v1.0.0 identity of the ordered record set.

Names and sequence identifiers determine the collection digest. Lengths are
stored and checked but do not alter that digest, matching the Sequence
Collections inherent/non-inherent attribute model.

### INSDC feature-state target

The adapter treats one ordered INSDC feature-table entry as one development
unit. Its eleven `u64` columns are a project-defined structural target ABI:
counts, bundle-local record and feature-key dictionary codes, kind and
orientation presence masks, and the exact normalized `segment_extent_sum`.
`segment_extent_sum` adds `end - start` for each 0-based half-open segment;
between sites contribute zero, uncertain-point locations contribute their
envelope width, and overlaps count with multiplicity.

The content-bound GenBank BioIR is the semantic authority for location ASTs,
segment associations, qualifiers, sequence, and provenance. The adapter emits
one ordered unit per feature and no inferred relationships. Its
reference-only index binds nine sealed children and the source artifact by
digest, allowing an independent implementation to replay the complete
projection.

Refget Sequences and Sequence Collections supply sequence and collection
identity algorithms in the GenBank frontend. Adapter artifact identities,
dictionaries, and masks remain in the Epigenesis ABI. Refget retrieval and
Sequence Collections service endpoints remain orthogonal transport services.

## Universal compiler architecture

“Universal” describes explicit profile dispatch: multiple registered source
languages compile into separately closed IRs that can be explicitly linked and
lowered. A source enters with a selected grammar, version, and authority, so its
native semantics remain mechanically checkable through compilation and replay.

A frontend qualifies for the architecture when it provides:

1. a named source language and pinned authority;
2. a finite grammar and closed normalized IR;
3. exact source identity and source maps;
4. explicit coordinate, reference, ordering, and missing-value semantics;
5. content-addressed cross-language links;
6. data-independent dispatch;
7. a target contract for every lowering; and
8. a separately implemented replay validator.

An API can transport or invoke these artifacts. The compiler is the language
implementation that performs acceptance, static checks, normalization,
lowering, emission, and rejection.

## Frontend map and roadmap

| Biological language | Primary specification | Compiler unit |
|---|---|---|
| Assembled sequence | [NCBI nucleotide FASTA](https://www.ncbi.nlm.nih.gov/genbank/fastaformat) | Existing FASTA, multi-FASTA, raw IUPAC, and gzip frontend on Sequence IR v2 / Sequence Collection v1; a future native v2 adapter can share chunked storage and refget collection identity |
| Reads with qualities | [Cock et al. FASTQ description](https://pmc.ncbi.nlm.nih.gov/articles/PMC2847217/) and [HTS format registry](https://samtools.github.io/hts-specs/) | A named FASTQ dialect must bind the selected quality-score convention, read names, pairing, sequence, and quality arrays instead of treating reads as assembled sequence |
| INSDC sibling records | [ENA EMBL flat file](https://ena-docs.readthedocs.io/en/latest/submit/fileprep/sequence-flatfile.html), [DDBJ flat file](https://www.ddbj.nig.ac.jp/ddbj/flat-file-e.html), and shared [INSDC Feature Table](https://www.insdc.org/submitting-standards/feature-table/) | EMBL and DDBJ outer-record frontends can lower their shared feature-table semantics into the same structural vocabulary while retaining each format's record grammar and source maps |
| Reads and alignments | [SAM/BAM 1.6 and CRAM 3.1](https://samtools.github.io/hts-specs/) | Read templates, flags, CIGAR, qualities, optional tags, reference bindings, coordinate conventions, compression dependencies, and indexes as a read/alignment IR |
| External feature annotation | [Sequence Ontology GFF3 1.26](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) and [NCBI GFF3](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/) | Existing bounded GFF3 BioIR v1 over Sequence Collection v1; future dialects add embedded FASTA, pragmas, and ontology-pinned semantics |
| Reference-relative variation | [VCF 4.5](https://samtools.github.io/hts-specs/VCFv4.5.pdf) | Header-defined fields, samples and genotypes, phasing, symbolic alleles, breakends, missingness, normalization state, and explicit reference binding |
| Computable variation | [GA4GH VRS 2.0.1](https://vrs.ga4gh.org/en/2.0/releases/2.0.html) | Release-pinned normalized VRS objects and identifiers linked to exact reference sequences |
| Synthetic design | [SBOL 3.1.0](https://sbolstandard.org/datamodel-specification/version-3.1.0/) | Components, sequences, features, constraints, interactions, models, and provenance as a design IR |
| Biochemical systems | [SBML Level 3 Version 2 Core, Release 2](https://sbml.org/specifications/sbml-level-3/version-2/core/release-2/sbml-level-3-version-2-release-2-core.pdf) | Species, compartments, reactions, units, rules, events, and package semantics |
| Mathematical cell models | [CellML 2.0.1](https://www.cellml.org/specifications/cellml_2.0/cellml_2_0_normative_specification.pdf) | Components, variables, units, equivalences, resets, MathML, and content-bound imports |
| Neural models | [NeuroML v2.3](https://docs.neuroml.org/Userdocs/NeuroMLv2.html) | Schema-bound neural structures and LEMS definitions as a model-language IR |
| Large neural networks | [SONATA developer guide](https://github.com/AllenInstitute/sonata/blob/master/docs/SONATA_DEVELOPER_GUIDE.md) | Node/edge tables, configuration, dynamics references, and index conventions as an external data bundle |

Identity and transport remain orthogonal to source-language parsing:

| Service standard | Version | Role in the architecture |
|---|---:|---|
| [GA4GH refget Sequences](https://ga4gh.github.io/refget/sequences/) | 2.0.0 | Content-derived identity and optional retrieval of one reference sequence |
| [GA4GH refget Sequence Collections](https://ga4gh.github.io/refget/seqcols/) | 1.0.0 | Content-derived identity and compatibility comparison for ordered reference collections |
| [GA4GH htsget](https://samtools.github.io/hts-specs/htsget.html) | 1.3.1 | Ticketed, range-selective transport of BAM/CRAM reads and VCF/BCF variants; retrieved bytes still enter their declared compiler frontend and evidence bundle |

Sequence Collection v2 currently belongs to the GenBank frontend. The retained
tensor, GFF3, and feature-graph paths declare Sequence Collection v1 and remain
separate until an explicit versioned adapter or native v2 lowering is specified.

## Prior biological compilers

Biological compilation has several established directions. The comparison below
locates Epigenesis by source language and target rather than by the shared word
“compiler.”

| Work | Compiled direction | Relationship |
|---|---|---|
| [Cai et al., 2009](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000529) | Constrained synthetic DNA parts to molecular-interaction and SBML models | Direct DNA-to-model precedent using an explicit attribute grammar and semantic actions |
| [Proto BioCompiler](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0022490) | High-level Proto programs to optimized genetic regulatory networks and simulation models | Biological programming language compiled toward networks and models |
| [Cello](https://pubmed.ncbi.nlm.nih.gov/27034378/) | Verilog logic to DNA using characterized genetic gates | Technology-mapped design-to-DNA compilation |
| [Eugene](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0018882) | Part/device DSL plus constraints to enumerated biological designs | Compiler-language and constraint system for synthetic design |
| [Genotype Specification Language](https://pubs.acs.org/doi/10.1021/acssynbio.5b00194) | High-level genotype specifications to large DNA constructs | Abstraction and automation for construct design |
| [BioCRNpyler](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009987) | Explicit components, mechanisms, and context to chemical reaction networks and SBML | Model compiler with biological semantics supplied in the source specification |
| [DNA Chisel](https://academic.oup.com/bioinformatics/article/36/16/4508/5869515) | DNA plus constraints and objectives to optimized DNA | Sequence optimization and constraint validation |

Epigenesis 0.9.0 takes a different implemented direction: real archival GenBank
bytes compile into exact sequence authority and structural BioIR, then lower to
an independently replayable feature-state development module.

## Release statement

Epigenesis 0.9.0 is a standards-pinned deterministic GenBank compiler and
feature-state lowering. It translates physical-DNA records into Sequence
Collection v2 and structural BioIR v2, then emits one constant development unit
per ordered feature through a sealed, independently replayable artifact closure.

Licensed under Apache-2.0.

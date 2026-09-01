# Standards and format trace

Epigenesis 1.0 treats each accepted biological format as an explicitly selected
source profile with defined grammar, identity rules, limits, and versioned IR.
Built-in raw DNA, FASTA, reference FASTA, GenBank, and bounded structural GFF3
frontends feed one source-descriptor contract. The GenBank profile carries a
packaged authority pin; raw/FASTA/GFF3 grammars are project-defined profiles.
The data-only external ABI binds caller-supplied grammar authorities and
manifest-qualified validation evidence for additional DNA-bearing profiles.

## GenBank v2 authority chain

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

Feature Table 11.4 supplies the location model used by the compiler. The
GenBank v2 profile, emitted with producer identity `brainc-insdc` 0.8.1,
represents points, closed spans, endpoint fuzz, between-base sites, archived
within positions, exact-version remote references, complements, joins, and
orders. It preserves the source expression as an AST and emits ordered 0-based
half-open segments.

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

“Universal” describes explicit profile dispatch: shipped built-ins and
manifest-qualified external DNA-bearing profiles compile into closed IRs that
can be linked and lowered through one target contract. Independent validation
reparses built-in sources. For an external profile, Epigenesis verifies the
manifest, supplied validation report, original/native bytes, sequence catalog,
and digest closure without executing the external grammar or validator.

A source frontend qualifies for this architecture when it provides:

1. a named DNA-bearing grammar and version, plus an authority URI and digest for
   an external profile;
2. a finite grammar and closed normalized IR;
3. exact source identity and source maps where the profile emits them;
4. explicit coordinate, reference, ordering, and missing-value semantics where
   applicable;
5. a nonempty ordered sequence catalog with lengths, SHA-256, and refget IDs;
6. data-independent dispatch and finite resource limits; and
7. deterministic native artifacts plus validation evidence independently
   checkable from exact bytes.

## Shipped and external-profile formats

| Source family | Primary authority or basis | Status in 1.0 | Compiler treatment |
|---|---|---|---|
| Raw IUPAC DNA | [INSDC nucleotide base codes](https://www.insdc.org/submitting-standards/feature-table/) | Built-in project profile | One named sequence lowered to exact Sequence Collection IR |
| FASTA | [NCBI nucleotide FASTA](https://www.ncbi.nlm.nih.gov/genbank/fastaformat) | Built-in Epigenesis grammar; NCBI is the format basis | Exact Sequence Collection route and streaming `fasta-reference-dna/v1` catalog route |
| Physical GenBank | [NCBI Release 273.0](https://www.ncbi.nlm.nih.gov/genbank/release/current/) and [INSDC 11.4](https://www.insdc.org/submitting-standards/feature-table/) | Built-in, packaged-authority-pinned | Traditional physical-DNA records to Sequence Collection v2 and structural INSDC BioIR |
| GFF3 plus sequence | [Sequence Ontology GFF3 1.26](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) | Built-in bounded structural profile | Syntax, coordinates, attributes, relationships, and exact sequence closure; ontology claims remain opaque |
| FASTQ reads | [Cock et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC2847217/) and [HTS registry](https://samtools.github.io/hts-specs/) | External profile required | Pin one FASTQ dialect and quality convention; the HTS registry notes that FASTQ has no formal definition and has incompatible variants |
| Physical EMBL/DDBJ records | [ENA EMBL flat file](https://ena-docs.readthedocs.io/en/latest/submit/fileprep/sequence-flatfile.html), [DDBJ flat file](https://www.ddbj.nig.ac.jp/ddbj/flat-file-e.html), and [INSDC 11.4](https://www.insdc.org/submitting-standards/feature-table/) | External profile required | Preserve each outer grammar while binding physical sequence and shared feature-table semantics |
| SAM/BAM 1.6 or CRAM 3.1 | [HTS specifications](https://samtools.github.io/hts-specs/) | External composite profile required | Bind read/alignment semantics and exact sequences; records with `SEQ=*` require a separately bound DNA source |
| VCF 4.5 or BCF 2.2 | [VCF 4.5/BCF 2.2](https://samtools.github.io/hts-specs/VCFv4.5.pdf) | External composite profile required | Bind the variant overlay to exact reference DNA and its declared normalization state |
| GFA 1.2 or GFA2 assembly graph | [GFA 1.2](https://gfa-spec.github.io/GFA-spec/GFA1.html) and [GFA2](https://gfa-spec.github.io/GFA-spec/GFA2.html) | External composite profile required | Bind graph topology and exact segment sequences; `*` sequence fields require a separately bound DNA source |

Identity and transport remain separate from source parsing:

| Service standard | Version | Role in the architecture |
|---|---:|---|
| [GA4GH refget Sequences](https://ga4gh.github.io/refget/sequences/) | 2.0.0 | Content-derived sequence identity and optional retrieval |
| [GA4GH refget Sequence Collections](https://ga4gh.github.io/refget/seqcols/) | 1.0.0 | Content-derived identity for ordered reference collections |
| [GA4GH htsget](https://samtools.github.io/hts-specs/htsget.html) | 1.3.1 | Transport for BAM/CRAM and VCF/BCF bytes that still enter a declared source profile |

Wire-version boundaries remain explicit: GenBank emits Sequence Collection v2;
GFF3 and retained feature-graph paths consume Sequence Collection v1.

## Release statement

Epigenesis 1.0 is a deterministic biological-source compiler. Its GenBank
profile is packaged-authority-pinned; other built-ins use project-defined
grammar versions; external profiles bind declared authorities and
manifest-qualified validation evidence. The independent validator reparses
built-in sources and verifies external source closures from exact evidence.

Licensed under Apache-2.0.

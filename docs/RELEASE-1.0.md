# Epigenesis 1.0 Release Contract

## Release boundary

Epigenesis 1.0 is the universal DNA translation compiler release. Its common
typed source boundary is the source descriptor, normalized `source_ir.records`
array, and profile-native typed artifacts. Its common translation target is
Development Module IR.

The release accepts biological interpretation as an explicit compiler input,
validates it against the exact DNA source, and lowers it to the selected target.
The external frontend ABI admits manifest-qualified data and evidence while
keeping frontend execution in the integrating toolchain.

## Version namespaces

The package and final source-to-development pipeline use release identity
`1.0.0`:

| Component | Producer identity |
|---|---|
| Source descriptor | `brainc-source/1.0.0` |
| Streaming reference FASTA | `brainc-reference-fasta/1.0.0` |
| External source closure | `brainc-external-source/1.0.0` |
| Development bundle | `brainc-development-bundle/1.0.0` |
| External source validator | `brainc-independent-external-source-validator/1.0.0` |

Compatibility commands retain their own implementation histories:

| Compatibility component | Producer identity |
|---|---|
| Sequence compiler | `brainc-dna/0.3.0` |
| Sequence collection compiler | `brainc-dna-collection/0.2.0` or `0.3.0`, selected by its wire path |
| GenBank frontend | `brainc-insdc-genbank/0.8.1` |
| GFF3 frontend | `brainc-bio-gff3/0.6.0` |
| GFF3 feature graph | `brainc-bio-feature-graph/0.6.0` |
| INSDC feature state | `brainc-insdc-feature-state/0.1.0` |
| Scalar state-program compiler | `brainc/0.4.0` |
| Typed tensor and target compiler | `brainc/0.5.0` |

These producer versions are not package versions and do not claim incomplete
1.0 source-to-development behavior. They identify the implementations that
remain available through compatibility commands.

Wire-format versions remain unchanged. A `v1` or `v2` suffix identifies a
contract schema and semantics, not the package release. Standards versions,
including GA4GH refget Sequence Collections 1.0.0, identify their cited
standards only.

## Compatibility and regeneration

This corrected 1.0 candidate retains the v1 wire contracts and final-pipeline
producer identities present at public `main` commit
`ac305318decc64a34afffee2bf42b048bc48ca81`. Existing artifacts from that commit
retain their identity when the same inputs and compiler contracts produce the
same payload.

The local 0.6 remediation was not committed, tagged, pushed, or published. Its
temporary `0.6.0` final-pipeline producer identities are not a supported wire
compatibility target. Regenerate artifacts created from that local draft with
the 1.0 candidate. The compiler does not silently reinterpret their producer
identity.

## Publication state

The repository worktree is a local candidate until all release gates pass and
an authorized commit, push, remote continuous integration run, exact remote
readback, tag, and package publication complete. Local test results do not
establish any of those states.

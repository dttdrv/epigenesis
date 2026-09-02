# External FASTQ frontend

This example adds FASTQ without adding a FASTQ parser to Epigenesis core. The
frontend and separately implemented validator communicate through the
version-2 executable external-profile ABI.

From the repository root:

```sh
python3.11 examples/external-fastq/profile.py > /tmp/epigenesis-fastq-profile.json
brainc compile-external-source \
  --profile-manifest /tmp/epigenesis-fastq-profile.json \
  --frontend-executable examples/external-fastq/frontend.py \
  --validator-executable examples/external-fastq/validator.py \
  --source-input reads=examples/external-fastq/reads.fastq \
  --output /tmp/epigenesis-fastq-source
```

The source bundle contains exact input identities, a common sequence-record
catalog, profile-native headers, sequences, qualities and record order, and the
validator replay closure.
External programs are explicitly selected and digest-checked; they are not
sandboxed. Use only programs you trust.

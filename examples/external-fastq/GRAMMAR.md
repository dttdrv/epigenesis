# Epigenesis FASTQ DNA profile 1

Input is one or more four-line FASTQ records. Line endings may be LF or CRLF;
bare carriage returns are invalid. The header starts with `@`, contains
printable ASCII, and its first token is the unique record id. The third line is
exactly `+`. Sequence is nonempty case-insensitive IUPAC DNA and is
canonicalized to uppercase. Quality has the same byte length as sequence and
every quality byte is in the printable Phred+33 range 33 through 126.

The profile-native FASTQ IR preserves the complete header, normalized sequence,
quality string, quality encoding, record id, and source order. The enclosing
source closure separately binds the exact physical bytes and line endings.

# Nothing fancy, just a simple store of useful(?) scripts

## Contents

### VCF

#### `subsample_vcf_evenly.sh`

Subsample an exact number of VCF records, distributing them as evenly as
possible across non-empty contigs and across each contig's record sequence.
The output retains the input header, uses natural contig order (for example,
`chr1`, `chr2`, ..., `chr10`, ..., `chrX`), is bgzip-compressed, and is supplied
with a tabix index.

```bash
VCF/subsample_vcf_evenly.sh INPUT.vcf.gz OUTPUT.vcf.gz RECORD_COUNT
```

For example:

```bash
VCF/subsample_vcf_evenly.sh \
  NA12877.vcf.gz \
  NA12877.subsample_100000.vcf.gz \
  100000
```

Use `-f` to replace an existing output and index. The input does not need to be
indexed. Requires `bcftools`, `awk`, `sort`, `cmp`, and `mktemp`. Compression
can use additional threads by setting `BCFTOOLS_THREADS`.

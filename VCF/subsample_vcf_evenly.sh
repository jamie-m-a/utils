#!/usr/bin/env bash

# Subsample an exact number of VCF records, distributing them as evenly as
# possible across non-empty contigs and across each contig's record sequence.
# The output is bgzip-compressed, naturally ordered by contig, and tabix-indexed.

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  subsample_vcf_evenly.sh [-f] INPUT.vcf[.gz] OUTPUT.vcf.gz RECORD_COUNT

Arguments:
  INPUT.vcf[.gz]  Source VCF; it does not need to be indexed.
  OUTPUT.vcf.gz   Destination bgzip-compressed VCF.
  RECORD_COUNT    Exact positive number of records to retain.

Options:
  -f              Replace OUTPUT.vcf.gz and its .tbi index if they exist.
  -h              Show this help text.

The script requires bcftools, awk, sort, cmp, and mktemp. Set BCFTOOLS_THREADS
to control compression threads (default: 1).
EOF
}

force=0
while getopts ':fh' option; do
  case "$option" in
    f) force=1 ;;
    h) usage; exit 0 ;;
    :) printf 'Error: -%s requires an argument\n' "$OPTARG" >&2; usage >&2; exit 2 ;;
    \?) printf 'Error: unknown option -%s\n' "$OPTARG" >&2; usage >&2; exit 2 ;;
  esac
done
shift $((OPTIND - 1))

if [[ $# -ne 3 ]]; then
  usage >&2
  exit 2
fi

input=$1
output=$2
target=$3
threads=${BCFTOOLS_THREADS:-1}

for command_name in bcftools awk sort cmp mktemp; do
  if ! command -v "$command_name" >/dev/null 2>&1; then
    printf 'Error: required command not found: %s\n' "$command_name" >&2
    exit 1
  fi
done

if [[ ! -r "$input" ]]; then
  printf 'Error: input is not readable: %s\n' "$input" >&2
  exit 1
fi
if [[ "$output" != *.vcf.gz ]]; then
  printf 'Error: output filename must end in .vcf.gz\n' >&2
  exit 2
fi
if [[ ! "$target" =~ ^[1-9][0-9]*$ ]]; then
  printf 'Error: RECORD_COUNT must be a positive integer\n' >&2
  exit 2
fi
if [[ ! "$threads" =~ ^[1-9][0-9]*$ ]]; then
  printf 'Error: BCFTOOLS_THREADS must be a positive integer\n' >&2
  exit 2
fi

output_dir=$(dirname "$output")
if [[ ! -d "$output_dir" || ! -w "$output_dir" ]]; then
  printf 'Error: output directory is not writable: %s\n' "$output_dir" >&2
  exit 1
fi

input_abs=$(cd "$(dirname "$input")" && pwd -P)/$(basename "$input")
output_abs=$(cd "$output_dir" && pwd -P)/$(basename "$output")
if [[ "$input_abs" == "$output_abs" ]]; then
  printf 'Error: input and output must be different files\n' >&2
  exit 2
fi

if (( force == 0 )) && [[ -e "$output" || -e "$output.tbi" ]]; then
  printf 'Error: output or index already exists (use -f to replace): %s\n' "$output" >&2
  exit 1
fi

tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/subsample_vcf_evenly.XXXXXX")
trap 'rm -rf -- "$tmpdir"' EXIT HUP INT TERM

counts_unsorted=$tmpdir/counts.unsorted.tsv
counts=$tmpdir/counts.tsv
quotas=$tmpdir/quotas.tsv
regions=$tmpdir/regions.tsv
selected=$tmpdir/selected.vcf.gz
final=$tmpdir/final.vcf.gz

printf 'Counting records by contig...\n' >&2
bcftools query -f '%CHROM\t%POS\n' "$input" \
  | awk 'BEGIN { OFS="\t" }
         { count[$1]++; if ($2 > max_pos[$1]) max_pos[$1]=$2 }
         END { for (chrom in count) print chrom, count[chrom], max_pos[chrom] }' \
  > "$counts_unsorted"

if [[ ! -s "$counts_unsorted" ]]; then
  printf 'Error: input contains no VCF records\n' >&2
  exit 1
fi

# Version sorting gives the expected human chromosome order: chr1, chr2, ...,
# chr10, ..., chrX. This order is also used for allocating remainder records.
LC_ALL=C sort -V -k1,1 "$counts_unsorted" > "$counts"

available=$(awk '{ total += $2 } END { print total + 0 }' "$counts")
if (( target > available )); then
  printf 'Error: requested %s records, but input contains only %s\n' \
    "$target" "$available" >&2
  exit 1
fi

# Capacity-aware water filling: quotas differ by at most one unless a small
# contig does not contain enough records, in which case its unused quota is
# redistributed among the remaining contigs.
awk -v target="$target" '
  BEGIN { OFS="\t" }
  {
    n++
    chrom[n]=$1
    capacity[n]=$2
    active[n]=1
  }
  END {
    remaining=target
    active_count=n
    while (active_count > 0) {
      share=int(remaining / active_count)
      removed=0
      for (i=1; i<=n; i++) {
        if (active[i] && capacity[i] <= share) {
          quota[i]=capacity[i]
          remaining-=quota[i]
          active[i]=0
          active_count--
          removed=1
        }
      }
      if (!removed) break
    }
    if (active_count > 0) {
      share=int(remaining / active_count)
      extra=remaining % active_count
      for (i=1; i<=n; i++) {
        if (active[i]) {
          quota[i]=share
          if (extra > 0 && capacity[i] > share) {
            quota[i]++
            extra--
          }
        }
      }
    }
    for (i=1; i<=n; i++) print chrom[i], capacity[i], quota[i]
  }
' "$counts" > "$quotas"

printf 'Selecting %s evenly spaced records...\n' "$target" >&2
bcftools view --no-version -Ov "$input" \
  | awk -v quota_file="$quotas" '
      BEGIN {
        while ((getline < quota_file) > 0) {
          available[$1]=$2
          quota[$1]=$3
        }
        close(quota_file)
      }
      /^#/ { print; next }
      {
        chrom=$1
        seen[chrom]++
        if (!(chrom in quota) || quota[chrom] == 0) next
        target_record=int(((chosen[chrom] + 0.5) * available[chrom]) / quota[chrom]) + 1
        if (seen[chrom] == target_record) {
          print
          chosen[chrom]++
        }
      }
    ' \
  | bcftools view --no-version --threads "$threads" -Oz -o "$selected"

bcftools index --threads "$threads" -f -t "$selected"

# Reading indexed contig regions in natural order fixes inputs whose records
# are grouped lexicographically (for example, chr1, chr10, ..., chr2).
awk 'BEGIN { OFS="\t" } { print $1, 1, $3 }' "$counts" > "$regions"
bcftools view --no-version --threads "$threads" -R "$regions" -Oz -o "$final" "$selected"
bcftools index --threads "$threads" -f -t "$final"

actual=$(bcftools index -n "$final")
if [[ "$actual" != "$target" ]]; then
  printf 'Error: validation failed: expected %s records, found %s\n' \
    "$target" "$actual" >&2
  exit 1
fi

bcftools view --no-version -h "$input" > "$tmpdir/input.header"
bcftools view --no-version -h "$final" > "$tmpdir/output.header"
if ! cmp -s "$tmpdir/input.header" "$tmpdir/output.header"; then
  printf 'Error: validation failed: output header differs from input header\n' >&2
  exit 1
fi

mv -f -- "$final" "$output"
mv -f -- "$final.tbi" "$output.tbi"

printf 'Created %s records in %s\n' "$actual" "$output" >&2
printf 'Created tabix index %s.tbi\n' "$output" >&2

#!/usr/bin/env bash
set -eu

separator='|'
rows=0
column=''

usage() {
    echo "Usage: $0 [-s separator] [-n rows] -c column file" >&2
    echo "Example: $0 -s '|' -n 10 -c SYMBOL variants.txt.gz" >&2
    exit 2
}

while getopts ':s:n:c:h' opt; do
    case "$opt" in
        s) separator=$OPTARG ;;
        n) rows=$OPTARG ;;
        c) column=$OPTARG ;;
        h) usage ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))

[[ -n "$column" && $# -eq 1 ]] || usage
[[ "$rows" =~ ^[0-9]+$ ]] || {
    echo "Row count must be a non-negative integer" >&2
    exit 2
}

file=$1
magic=$(LC_ALL=C od -An -tx1 -N2 "$file" | tr -d ' ')

if [[ "$magic" == "1f8b" ]]; then
    reader=(gzip -cd -- "$file")
else
    reader=(cat -- "$file")
fi

"${reader[@]}" |
awk -v FS="$separator" -v wanted="$column" -v limit="$rows" '
NR == 1 {
    sub(/\r$/, "", $NF)

    for (i = 1; i <= NF; i++) {
        if ($i == wanted) {
            column = i
            break
        }
    }

    if (!column) {
        printf "Column not found: %s\n", wanted > "/dev/stderr"
        exit 2
    }

    next
}

limit == 0 || emitted < limit {
    print $column
    emitted++
}
'

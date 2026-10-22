#!/bin/bash
# z3xz6_sweep -- probe: does the per-model awk extraction lose any label of a list?
# Written after a CRLF-terminated model list silently dropped its LAST label.
# probe: does the per-model awk extraction lose any of a given label list?
set -u
SRC=/mnt/f/physics/multi_geometry_campaign/data/literature/parr2020_full/extracted/Z3xZ6_1_1.txt
OUT=/tmp/extract_probe.txt
: > "$OUT"
for M in "$@"; do
  n=$(awk -v want="Label:$M" '
    {sub(/\r$/, "")}
    /^begin model/ {buf=$0"\n"; inb=1; keep=0; next}
    inb {buf=buf $0"\n"; if ($0 == want) keep=1}
    /^end model/ {if (inb && keep) {printf "%s", buf} inb=0; keep=0}
  ' "$SRC" | grep -c '^Label:')
  [ "$n" -eq 1 ] || echo "LOST or DUPLICATED: $M -> $n"
done
echo "probe done over $# labels"

#!/bin/bash
# z3xz6_sweep -- reproduce the coupling driver's own extraction loop verbatim and report
# any label that fails to land, so an extraction loss can never be silent.
# Reproduce the driver's own extraction loop verbatim and count what lands.
set -u
SRC=/mnt/f/physics/multi_geometry_campaign/data/literature/parr2020_full/extracted/Z3xZ6_1_1.txt
W=/tmp/exprobe; mkdir -p "$W"; : > "$W/models.txt"
for M in "$@"; do
  awk -v want="Label:$M" '
    {sub(/\r$/, "")}
    /^begin model/ {buf=$0"\n"; inb=1; keep=0; next}
    inb {buf=buf $0"\n"; if ($0 == want) keep=1}
    /^end model/ {if (inb && keep) {printf "%s", buf} inb=0; keep=0}
  ' "$SRC" >> "$W/models.txt"
done
echo "requested $# | landed $(grep -c '^Label:' "$W/models.txt")"
grep '^Label:' "$W/models.txt" | sed 's/Label://' | sort > "$W/got.txt"
printf '%s\n' "$@" | sort > "$W/want.txt"
echo "missing: $(comm -13 "$W/got.txt" "$W/want.txt" | tr '\n' ' ')"

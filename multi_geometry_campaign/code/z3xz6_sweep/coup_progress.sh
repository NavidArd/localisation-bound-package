#!/bin/bash
# z3xz6_sweep -- in-flight progress probe for the coupling segments.
for d in ~/z36sw/coup*/; do
  [ -f "$d/out.txt" ] || continue
  n=$(grep -c '^  W = ' "$d/out.txt" 2>/dev/null)
  m=$(grep -c '^Label:' "$d/models.txt" 2>/dev/null)
  age=$(( $(date +%s) - $(stat -c %Y "$d/cmd.txt") ))
  echo "$(basename $d) $n/$m  age=${age}s  bytes=$(stat -c %s "$d/out.txt")"
done

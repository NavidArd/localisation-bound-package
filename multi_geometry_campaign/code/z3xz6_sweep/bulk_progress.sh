#!/bin/bash
# z3xz6_sweep -- in-flight probe: completed / stopped model counts per bulk chunk.
# Reads orbifolder's own completion signal, not byte counts.
tot=0
for d in "$HOME"/z36sw/coupB*/; do
  [ -f "$d/out.txt" ] || continue
  n=$(grep -c ' done: "create coupling' "$d/out.txt" 2>/dev/null)
  s=$(grep -c ' stopped: "create coupling' "$d/out.txt" 2>/dev/null)
  tot=$((tot+n))
  echo "$(basename "$d") done=$n stopped=$s"
done
echo "TOTAL completed models in flight: $tot"

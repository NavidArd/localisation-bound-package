#!/bin/bash
# z3xz6_sweep -- in-flight probe of orbifolder's coupling COMPLETION SIGNAL.
#   `PID <n> done: "create coupling(F F F)" TIME: hh:mm:ss`     -> child finished, W valid
#   `PID <n> stopped: "create coupling(F F F)" TIME: hh:mm:ss`  -> child KILLED by wait(N);
#                                                                  a printed `W = 0` here
#                                                                  is a TRUNCATION
for d in ~/z36sw/coup*/; do
  [ -f "$d/out.txt" ] || continue
  nd=$(grep -c ' done: "create coupling' "$d/out.txt" 2>/dev/null)
  ns=$(grep -c ' stopped: "create coupling' "$d/out.txt" 2>/dev/null)
  echo "$(basename "$d") done=$nd stopped=$ns"
done
echo "--- creation TIME distribution over all live segments (done only):"
grep -h -o ' done: "create coupling(F F F)" TIME: [0-9:]*' ~/z36sw/coup*/out.txt 2>/dev/null \
  | awk -F'TIME: ' '{print $2}' | sort | uniq -c | sort -k2 | tail -20

#!/bin/bash
# z3xz6_sweep -- in-flight probe: orbifolder process states, per-chunk output size and
# age, and system load.  Used to distinguish a computing child (R, ~100% CPU) from a
# parent sleeping out its wait(N) (S, ~0%) from an unreaped zombie (Z).
echo "=== orbifolder processes (state, elapsed s, %cpu, rss MB) ==="
ps -eo stat,etimes,pcpu,rss,pid,ppid,comm | grep orbifolder | awk '{printf "%-5s %-7s %-6s %-8.0f %-7s %-7s\n",$1,$2,$3,$4/1024,$5,$6}' | head -20
echo "total: $(pgrep -c -x orbifolder)"
echo "=== per-chunk output ==="
for d in "$HOME"/z36sw/coupRP*/; do
  [ -f "$d/out.txt" ] || continue
  printf "%-14s bytes=%-9s age=%-6ss dots=%s\n" "$(basename "$d")" \
    "$(stat -c %s "$d/out.txt")" \
    "$(( $(date +%s) - $(stat -c %Y "$d/cmd.txt") ))" \
    "$(tail -c 200 "$d/out.txt" | grep -o '\.' | wc -l)"
done
echo "=== load ==="; uptime

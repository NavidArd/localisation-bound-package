#!/bin/bash
# z3xz6_sweep -- pause the shuffled bulk coupling run without touching the priority
# chunks, so a decisive batch is not starved by the background sweep.
# Pause the shuffled bulk coupling run without touching the priority (RP) chunks.
pkill -f rk_bulk_z36 2>/dev/null
for p in $(pgrep -f orbifolder 2>/dev/null); do
  d=$(tr '\0' ' ' < /proc/$p/cmdline 2>/dev/null)
  case "$d" in *coupB*) kill $p 2>/dev/null;; esac
done
sleep 2
echo "orbifolder processes left: $(pgrep -c -f orbifolder 2>/dev/null || echo 0)"

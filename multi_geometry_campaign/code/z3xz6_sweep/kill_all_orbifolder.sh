#!/bin/bash
# Hard-stop every orbifolder process.  Needed because the multi-model command stream lets
# `create coupling(F F F)` fork descendants that outlive the launcher shell: killing the
# launcher leaves ~40 sleeping PARENTS which keep spawning fresh children, so a single
# pkill pass never converges.  Loop until the count is zero.
for i in $(seq 1 12); do
  pkill -9 -x orbifolder 2>/dev/null
  pkill -9 -f 'src/orbifolder/orbifolder' 2>/dev/null
  sleep 1
  n=$(pgrep -c -x orbifolder 2>/dev/null || echo 0)
  echo "pass $i: $n left"
  [ "$n" -eq 0 ] && break
done

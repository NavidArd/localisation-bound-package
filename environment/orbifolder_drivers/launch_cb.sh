#!/bin/bash
# su5_construction -- parallel launcher for the CB-1 orbifolder census.
# Six concurrent segments (the PI's laptop has 6 cores).  Resumable: a segment whose
# .txt.gz already exists is skipped by the driver.
set -u
DRV=/mnt/f/physics/multi_geometry_campaign/code/su5_construction/run_orbifolder_cb.sh
N=${1:-2912}; CH=${2:-100}
S=0; I=0
while [ $S -lt $N ]; do
  E=$((S+CH-1)); [ $E -ge $N ] && E=$((N-1))
  bash "$DRV" $I $S $E > /tmp/cb_$I.log 2>&1 &
  I=$((I+1)); S=$((E+1))
  while [ "$(jobs -rp | wc -l)" -ge 6 ]; do sleep 3; done
done
wait
echo ALLDONE

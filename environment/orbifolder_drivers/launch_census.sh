#!/bin/bash
# z2xz4_sweep -- census LAUNCHER.  Stages the geometry ONCE (the shared-geometry race the
# Z3xZ6 pass diagnosed cost 16 models), then runs SEGMENTS with a bounded worker pool.
#
# Segments are 200 models wide, so no more than 200 model-units are ever at risk from a
# single lost worker, and a partial segment is re-run WHOLE (never repaired in place).
# A segment whose .gz already exists is skipped, so the launcher is RESUMABLE.
#
# Usage: wsl -d Ubuntu -- bash <this file> [NWORKERS] [SEGSIZE] [NMODELS]
set -u
BASE=/mnt/f/physics/multi_geometry_campaign
GEO=/mnt/f/physics/orbifolder-1.2.1/Geometry/Geometry_Z2xZ4_1_1.txt
RUN=/tmp/rc24.sh
N=${1:-6}; SZ=${2:-200}; TOT=${3:-11540}

mkdir -p ~/rcharge_build/Geometry $BASE/data/z2xz4_sweep/logs
cp -f "$GEO" ~/rcharge_build/Geometry/          # ONCE, before any worker starts
sed 's/\r$//' $BASE/code/z2xz4_sweep/run_census_z24.sh > $RUN

i=0
while [ $((i * SZ)) -lt $TOT ]; do
  A=$((i * SZ)); B=$((A + SZ - 1)); [ $B -ge $TOT ] && B=$((TOT - 1))
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done
  bash $RUN "$i" "$A" "$B" >> $BASE/data/z2xz4_sweep/logs/census_launch.log 2>&1 &
  i=$((i + 1))
done
wait
echo "ALL SEGMENTS DONE (n=$i)" | tee -a $BASE/data/z2xz4_sweep/logs/census_launch.log

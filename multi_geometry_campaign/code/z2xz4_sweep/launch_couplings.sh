#!/bin/bash
# z2xz4_sweep -- coupling LAUNCHER.  Stages the geometry ONCE, splits the model list into
# chunks, and runs them with a bounded worker pool.  Resumable: a chunk whose .gz exists is
# skipped and can be re-run WHOLE.
#
# Usage: wsl -d Ubuntu -- bash <this file> [NWORKERS] [CHUNKSIZE] [LISTFILE]
set -u
BASE=/mnt/f/physics/multi_geometry_campaign
GEO=/mnt/f/physics/orbifolder-1.2.1/Geometry/Geometry_Z2xZ4_1_1.txt
RUN=/tmp/rcc24.sh
N=${1:-6}; SZ=${2:-25}
LIST=${3:-$BASE/data/z2xz4_sweep/coupling_model_list.txt}

mkdir -p ~/rcharge_build/Geometry ~/z24sw/lists $BASE/data/z2xz4_sweep/logs
cp -f "$GEO" ~/rcharge_build/Geometry/
sed 's/\r$//' $BASE/code/z2xz4_sweep/run_couplings_z24.sh > $RUN
sed 's/\r$//' "$LIST" | grep -v '^$' > ~/z24sw/lists/all.txt
split -l $SZ -d -a 4 ~/z24sw/lists/all.txt ~/z24sw/lists/chunk_

for f in ~/z24sw/lists/chunk_*; do
  ID=${f##*chunk_}
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done
  bash $RUN "$ID" "$f" >> $BASE/data/z2xz4_sweep/logs/coup_launch.log 2>&1 &
done
wait
echo "ALL COUPLING CHUNKS DONE" | tee -a $BASE/data/z2xz4_sweep/logs/coup_launch.log

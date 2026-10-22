#!/bin/bash
# z3xz6_sweep -- launch the COUPLING runs.
#
# SIZING, from `pilot_timing_z36.sh` (recorded in the findings):
#   * `create coupling(F F F)` completes in ~12 s/model; the foundations pass's 601 s was
#     the wait(600) itself, not the computation -- at wait in {1,2,3,5,10} the printed
#     superpotential is BYTE-IDENTICAL to the wait(600) reference (908 terms, same md5).
#   * a 10-model spread segment took 123 s wall => 12.3 s/model.
#   * PROJECTION for the whole catalogue: 4,612 x 12.3 s = 15.8 CPU-hours; at 5 workers
#     ~3.2 h wall.  That is inside the 48 h ceiling, so the sweep runs the COUPLINGS ON
#     EVERY MODEL and the step-2 pre-filter is NOT load-bearing: it is stated, applied,
#     and then CHECKED against the full coupling data.  A filter that gated nothing cannot
#     have hidden a survivor.
#
# wait=30 is used in production -- a 2.5x margin over the measured completion, with the
# truncation guard (`pilot_waitguard_z36.sh`) re-running a spread at wait=120 and
# comparing the superpotential block byte for byte.
#
# Segment bookkeeping and resumability are as for the census: one gz per segment, an
# existing gz is SKIPPED, a partial segment is re-run whole and never repaired.
#
# Usage: wsl -d Ubuntu -- bash <this file> [SEGSZ] [WORKERS] [WAIT]
set -u
BASE=/mnt/f/physics/multi_geometry_campaign
SRC=$BASE/data/literature/parr2020_full/extracted/Z3xZ6_1_1.txt
SEGSZ=${1:-100}
WORKERS=${2:-5}
WAITN=${3:-30}
NMODELS=4612
RK=/tmp/rk_all_z36.sh
sed 's/\r$//' "$BASE/code/z3xz6_sweep/run_couplings_z36.sh" > $RK

# stage the geometry ONCE, before any worker starts (see run_*_z36.sh note)
GEO_STAGE=$BASE/data/z3xz6/patched_geometry/Geometry_Z3xZ6_1_1_Rsym_prime.txt
mkdir -p ~/rcharge_build/Geometry
cp -f "$GEO_STAGE" ~/rcharge_build/Geometry/

mkdir -p "$BASE/data/z3xz6_sweep/raw/coup" "$BASE/data/z3xz6_sweep/logs"
NSEG=$(( (NMODELS + SEGSZ - 1) / SEGSZ ))
echo "couplings: $NMODELS models, segsize $SEGSZ -> $NSEG segments, $WORKERS workers, wait=$WAITN"

seq 0 $((NSEG - 1)) | xargs -P "$WORKERS" -I{} bash -c \
  'S={}; SZ='"$SEGSZ"'; A=$((S*SZ)); B=$((S*SZ+SZ-1)); L=""
   for ((i=A;i<=B && i<'"$NMODELS"';i++)); do L="$L MSSM$i"; done
   bash '"$RK"' "$S" '"$WAITN"' $L >> '"$BASE"'/data/z3xz6_sweep/logs/coup_launch.log 2>&1'

echo "coupling launcher finished; segments present:"
ls -1 "$BASE/data/z3xz6_sweep/raw/coup" | grep -c '^coup_seg[0-9]'

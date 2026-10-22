#!/bin/bash
# z3xz6_sweep -- re-run the models that `coverage_z36.py` found missing from the
# published coupling segments, in small batches, and publish them as `coup_segGAPn`.
#
# Small batches (default 10) for two reasons: a batch that loses a label loses less, and
# the coverage gate is re-run afterwards so any residue is caught and re-queued.  A model
# still missing after two passes is reported, never silently dropped.
#
# Usage: wsl -d Ubuntu -- bash <this file> [BATCH] [WAIT] [WORKERS]
set -u
BASE=/mnt/f/physics/multi_geometry_campaign
LIST=$BASE/data/z3xz6_sweep/missing_models.txt
RK=/tmp/rk_gap_z36.sh
sed 's/\r$//' "$BASE/code/z3xz6_sweep/run_couplings_z36.sh" > $RK
GEO=$BASE/data/z3xz6/patched_geometry/Geometry_Z3xZ6_1_1_Rsym_prime.txt
mkdir -p ~/rcharge_build/Geometry; cp -f "$GEO" ~/rcharge_build/Geometry/

B=${1:-10}; WAITN=${2:-60}; P=${3:-6}
[ -s "$LIST" ] || { echo "no missing list at $LIST"; exit 0; }
tr ' ' '\n' < "$LIST" | sed '/^$/d' > /tmp/gap_all.txt
n=$(wc -l < /tmp/gap_all.txt)
echo "filling $n missing models, batch $B, wait $WAITN, workers $P"
split -l "$B" -d -a 4 /tmp/gap_all.txt /tmp/gapb_
ls /tmp/gapb_* | sed 's|.*/gapb_||' | xargs -P "$P" -I{} bash -c \
  'bash '"$RK"' "GAP{}" '"$WAITN"' $(cat /tmp/gapb_{}) >> '"$BASE"'/data/z3xz6_sweep/logs/coup_launch.log 2>&1'
rm -f /tmp/gapb_*
echo "gap fill finished"

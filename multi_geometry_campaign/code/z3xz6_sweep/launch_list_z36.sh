#!/bin/bash
# z3xz6_sweep -- run the coupling driver over an explicit MODEL LIST FILE, in parallel
# chunks.  Used for the decisive subsets (the Rule-P models) and for gap/straggler
# re-runs, where the model set is not a contiguous index range.
#
# Kept as a FILE rather than an inline shell string on purpose: every attempt to inline
# this through `wsl -- bash -c "..."` lost a level of quoting and silently produced EMPTY
# model lists (segments published with models=0).  A file has one level of quoting.
#
# Usage: wsl -d Ubuntu -- bash <this file> <LISTFILE> <TAG> <WAIT> <CHUNKS> <WORKERS>
set -u
BASE=/mnt/f/physics/multi_geometry_campaign
LIST=$1; TAG=$2; WAITN=${3:-200}; CH=${4:-6}; P=${5:-6}
RK=/tmp/rk_list_z36.sh
sed 's/\r$//' "$BASE/code/z3xz6_sweep/run_couplings_z36.sh" > $RK
GEO=$BASE/data/z3xz6/patched_geometry/Geometry_Z3xZ6_1_1_Rsym_prime.txt
mkdir -p ~/rcharge_build/Geometry; cp -f "$GEO" ~/rcharge_build/Geometry/

WD=/tmp/list_$TAG; rm -rf "$WD"; mkdir -p "$WD"
tr ' ' '\n' < "$LIST" | sed 's/\r$//' | sed '/^$/d' | sort -u > "$WD/all.txt"
N=$(wc -l < "$WD/all.txt")
[ "$N" -gt 0 ] || { echo "EMPTY LIST $LIST"; exit 1; }
echo "$TAG: $N models, $CH chunks, $P workers, wait=$WAITN"
split -n "l/$CH" -d -a 2 "$WD/all.txt" "$WD/b_"

for f in "$WD"/b_*; do echo "$f"; done | xargs -P "$P" -I{} bash -c \
  'f={}; s=$(basename "$f"); bash '"$RK"' "'"$TAG"'${s#b_}" '"$WAITN"' $(cat "$f") \
     >> '"$BASE"'/data/z3xz6_sweep/logs/coup_launch.log 2>&1'
echo "$TAG done"

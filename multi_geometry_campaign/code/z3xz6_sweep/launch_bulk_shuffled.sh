#!/bin/bash
# z3xz6_sweep -- the FULL-CATALOGUE coupling run, in RANDOMISED model order.
#
# WHY SHUFFLED.  Per-model coupling creation costs 20-70 s of real work, so the whole
# catalogue is a multi-hour job.  Running it in catalogue-index order would make any
# partial completion a PREFIX -- and catalogue index correlates with the PVW scan's own
# provenance order, so a prefix is a biased sample and supports no catalogue-level
# statement.  With a fixed-seed shuffle, whatever fraction completes is a UNIFORM RANDOM
# SAMPLE of the 4,612, and every partial statistic carries an honest binomial error bar
# instead of an unquantifiable selection effect.
#
# The seed is fixed and recorded so the ordering is reproducible.
#
# Usage: wsl -d Ubuntu -- bash <this file> [CHUNK] [WORKERS] [WAIT] [SEED]
set -u
BASE=/mnt/f/physics/multi_geometry_campaign
CH=${1:-25}; P=${2:-11}; WAITN=${3:-12}; SEED=${4:-20260818}
RK=/tmp/rk_bulk_z36.sh
sed 's/\r$//' "$BASE/code/z3xz6_sweep/run_couplings_z36.sh" > $RK
GEO=$BASE/data/z3xz6/patched_geometry/Geometry_Z3xZ6_1_1_Rsym_prime.txt
mkdir -p ~/rcharge_build/Geometry; cp -f "$GEO" ~/rcharge_build/Geometry/
mkdir -p "$BASE/data/z3xz6_sweep/logs"

WD=/tmp/bulk_z36; mkdir -p "$WD"
if [ ! -s "$WD/order.txt" ]; then
  python3 - "$SEED" > "$WD/order.txt" <<'PY'
import random, sys
r = random.Random(int(sys.argv[1]))
labs = ["MSSM%d" % i for i in range(4612)]
r.shuffle(labs)
print("\n".join(labs))
PY
fi
N=$(wc -l < "$WD/order.txt")
NCH=$(( (N + CH - 1) / CH ))
echo "bulk: $N models, chunk $CH -> $NCH chunks, $P workers, wait=$WAITN, seed=$SEED"
rm -f "$WD"/c_*
split -l "$CH" -d -a 4 "$WD/order.txt" "$WD/c_"

for f in "$WD"/c_*; do echo "$f"; done | xargs -P "$P" -I{} bash -c \
  'f={}; s=$(basename "$f"); bash '"$RK"' "B${s#c_}" '"$WAITN"' $(cat "$f") \
     >> '"$BASE"'/data/z3xz6_sweep/logs/coup_launch.log 2>&1'
echo "bulk done"

#!/bin/bash
# z3xz6_sweep -- launch the CENSUS over all 4,612 Z3xZ6_1_1 models.
#
# SEGMENT BOOKKEEPING (the campaign's environment rule).  The catalogue is cut into
# fixed index segments of SEGSZ models; each segment is one orbifolder invocation and
# one gzipped raw file `census_seg<N>.txt.gz`.  A segment whose gz already exists is
# SKIPPED by run_census_z36.sh, so the whole launcher is resumable: re-running it after
# an interruption completes only the missing segments.  Nothing is ever appended to an
# existing raw file, so a partial segment is re-run whole rather than repaired.
#
# Parallelism: WORKERS concurrent orbifolder processes (the machine has 6 physical
# cores; orbifolder is single-threaded and not BLAS-bound, so no thread pinning is
# needed here -- the CLAUDE.md BLAS rule engages only for the numpy stage, which this
# sweep does not use).
#
# Usage: wsl -d Ubuntu -- bash <this file> [SEGSZ] [WORKERS]
set -u
BASE=/mnt/f/physics/multi_geometry_campaign
SEGSZ=${1:-100}
WORKERS=${2:-6}
NMODELS=4612
RC=/tmp/rc_z36.sh
sed 's/\r$//' "$BASE/code/z3xz6_sweep/run_census_z36.sh" > $RC

# stage the geometry ONCE, before any worker starts (see run_*_z36.sh note)
GEO_STAGE=/mnt/f/physics/orbifolder-1.2.1/Geometry/Geometry_Z3xZ6_1_1.txt
mkdir -p ~/rcharge_build/Geometry
cp -f "$GEO_STAGE" ~/rcharge_build/Geometry/

mkdir -p "$BASE/data/z3xz6_sweep/raw/census" "$BASE/data/z3xz6_sweep/logs"
NSEG=$(( (NMODELS + SEGSZ - 1) / SEGSZ ))
echo "census: $NMODELS models, segsize $SEGSZ -> $NSEG segments, $WORKERS workers"

seq 0 $((NSEG - 1)) | xargs -P "$WORKERS" -I{} bash -c \
  'S={}; SZ='"$SEGSZ"'; bash '"$RC"' "$S" $((S*SZ)) $((S*SZ+SZ-1)) >> '"$BASE"'/data/z3xz6_sweep/logs/census_launch.log 2>&1'

echo "census launcher finished; segments present:"
ls -1 "$BASE/data/z3xz6_sweep/raw/census" | wc -l

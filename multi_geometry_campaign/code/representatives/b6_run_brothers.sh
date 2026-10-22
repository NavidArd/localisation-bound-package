#!/bin/bash
# Run orbifolder 1.2.1 on every segment decks_NNN.txt of the representative census
# (b5_build_brother_census.py + b5b_missing.py), stock geometry Z3xZ6_1_1 of this package,
# same command order as su5_construction/run_orbifolder_cb.sh (spectrum printed before any
# vev-configuration change).  Output gzipped per segment; segments present are skipped.
# Usage:  ORBIFOLDER_BIN=/path/to/orbifolder bash b6_run_brothers.sh     (6 parallel workers)
set -u
HERE=$(cd "$(dirname "$0")" && pwd); PKG=$(cd "$HERE/../../.." && pwd)
D=$PKG/multi_geometry_campaign/data/representatives/brothers
GEO=$PKG/orbifolder-1.2.1/Geometry/Geometry_Z3xZ6_1_1.txt
BIN=${ORBIFOLDER_BIN:-orbifolder}
TMPROOT=${TMPDIR:-/tmp}/rep_b6_$$
one() {
  S=$1; OUTG=$D/out_$S.txt.gz
  [ -s "$OUTG" ] && { echo "skip $S"; return; }
  W=$TMPROOT/$S; rm -rf "$W"; mkdir -p "$W/Geometry"; cp -f "$GEO" "$W/Geometry/"
  tr -d '\r' < "$D/decks_$S.txt" > "$W/models.txt"
  { for L in $(grep '^Label:' "$W/models.txt" | sed 's/Label://'); do
      echo "cd $L"; echo "cd spectrum"; echo "print(F) with internal information"; echo "cd .."
      echo "cd gauge group"; echo "print gauge group"; echo "cd .."; echo "cd .."
    done; echo exit; echo yes; } > "$W/cmd.txt"
  cd "$W"; T0=$(date +%s)
  timeout 7200 "$BIN" models.txt < cmd.txt > out.txt 2>&1; rc=$?
  gzip -c out.txt > "$OUTG.tmp" && mv "$OUTG.tmp" "$OUTG"
  echo "seg $S exit=$rc secs=$(( $(date +%s)-T0 ))" >> "$D/run.log"
  cd /; rm -rf "$W"
}
export -f one; export D GEO BIN TMPROOT
ls "$D" | sed -n 's/^decks_\([0-9]*\)\.txt$/\1/p' | xargs -P 6 -I{} bash -c 'one {}'
echo ALLDONE >> "$D/run.log"

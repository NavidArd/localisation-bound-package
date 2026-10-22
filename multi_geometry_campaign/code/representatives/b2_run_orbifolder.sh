#!/bin/bash
# Run orbifolder 1.2.1 on a deck file (the 19 representative decks b1_decks.txt of the model
# of section 7.1, or the nine V2 decks b11_decks.txt), with the STOCK geometry file
# Geometry_Z3xZ6_1_1.txt of this package and the command sequence of
# su5_construction/run_orbifolder_cb.sh (spectrum printed before any vev-configuration
# change).  Usage:  ORBIFOLDER_BIN=/path/to/orbifolder bash b2_run_orbifolder.sh <decks> <out>
set -u
SRC=$1; OUTF=$2
HERE=$(cd "$(dirname "$0")" && pwd); PKG=$(cd "$HERE/../../.." && pwd)
GEO=$PKG/orbifolder-1.2.1/Geometry/Geometry_Z3xZ6_1_1.txt
BIN=${ORBIFOLDER_BIN:-orbifolder}
W=${TMPDIR:-/tmp}/rep_b2_$$; rm -rf "$W"; mkdir -p "$W/Geometry"; cp -f "$GEO" "$W/Geometry/"
tr -d '\r' < "$SRC" > "$W/models.txt"
LABELS=$(grep '^Label:' "$W/models.txt" | sed 's/Label://')
{ for L in $LABELS; do
    echo "cd $L"; echo "cd spectrum"; echo "print(F) with internal information"; echo "cd .."
    echo "cd gauge group"; echo "print gauge group"; echo "print simple roots"; echo "cd .."
    echo "cd spectrum"; echo "print summary"; echo "cd .."; echo "cd .."
  done; echo exit; echo yes; } > "$W/cmd.txt"
cd "$W"
timeout 7200 "$BIN" models.txt < cmd.txt > out.txt 2>&1
echo "exit=$?"
cp -f out.txt "$OUTF"
cd /; rm -rf "$W"

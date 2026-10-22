#!/bin/bash
# gate_su5_construction -- fresh orbifolder run (V1 / V4).
# Same STOCK geometry and the same command order as the certified census driver, but the
# model deck is written by the GATE from the embedding data, not copied from the phase.
set -u
BASE=/mnt/f/physics/multi_geometry_campaign
SRC=$BASE/data/gate_su5_construction/${2:-gate_models.txt}
GEO=/mnt/f/physics/orbifolder-1.2.1/Geometry/Geometry_Z3xZ6_1_1.txt
DEST=$BASE/data/gate_su5_construction
BIN=~/rcharge_build/src/orbifolder/orbifolder
TAG=$1
W=~/gsc/$TAG
mkdir -p "$W" ~/rcharge_build/Geometry
[ -s ~/rcharge_build/Geometry/"$(basename "$GEO")" ] || cp -f "$GEO" ~/rcharge_build/Geometry/
LABELS=$(grep '^Label:' "$SRC" | sed 's/Label://')
{ for L in $LABELS; do
    echo "cd $L"; echo "cd spectrum"; echo "print(F) with internal information"; echo "cd .."
    echo "cd gauge group"; echo "print gauge group"; echo "print U1 generators"
    echo "print simple roots"; echo "print anomaly info"; echo "cd .."
    echo "cd vev-config"; echo "analyze config print SU(5) simple roots"
    echo "use config(SMConfig)"; echo "print gauge group"; echo "cd .."
    echo "cd spectrum"; echo "print summary"; echo "cd .."; echo "cd .."
  done; echo exit; echo yes; } > "$W/cmd.txt"
cd ~/rcharge_build
T0=$(date +%s)
timeout 20000 "$BIN" "$SRC" < "$W/cmd.txt" > "$W/out.txt" 2>&1
rc=$?
gzip -c "$W/out.txt" > "$DEST/gate_$TAG.txt.gz"
echo "gate_$TAG exit=$rc secs=$(( $(date +%s)-T0 )) bytes=$(wc -c < "$W/out.txt")"

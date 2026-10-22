#!/bin/bash
# z2xz4_sweep -- the print(bF) LEFT-CHIRAL-BASIS SPREAD PROBE.
#
# WHY THIS EXISTS.  `print(F)` CANNOT emit an `_r` state, so any "0 right-chiral states"
# check run on print(F) output is VACUOUS BY COMMAND (gate_z4xz4_foundations K8).  The
# census therefore cannot, by itself, establish which sectors host left-chiral fields --
# and this phase's CPT statement (which sectors are labelling images rather than physical
# absences) rests on exactly that.  `gate_z4xz4_sweep.md` attack item 5 charged the Z4xZ4
# phase with resting a 13-sector left-chiral basis on a FIVE-model print(bF) measurement.
# That charge is honoured up front here: this probe runs print(bF) on a SPREAD SAMPLE of
# 240 models (every 48th model of the catalogue), which is the evidence base the basis
# statement is then scoped to.
#
# Usage: wsl -d Ubuntu -- bash <this file> [STRIDE]
set -u
BASE=/mnt/f/physics/multi_geometry_campaign
SRC=$BASE/data/literature/parr2020_full/extracted/Z2xZ4_1_1.txt
GEO=/mnt/f/physics/orbifolder-1.2.1/Geometry/Geometry_Z2xZ4_1_1.txt
DEST=$BASE/data/z2xz4_sweep/raw
BIN=~/rcharge_build/src/orbifolder/orbifolder
STRIDE=${1:-48}

W=~/z24sw/bf
mkdir -p "$W" "$DEST" ~/rcharge_build/Geometry
cp -f "$GEO" ~/rcharge_build/Geometry/

awk -v s="$STRIDE" '
  {sub(/\r$/, "")}
  /^begin model/ {buf=$0"\n"; inb=1; keep=0; next}
  inb {buf=buf $0"\n"
       if ($0 ~ /^Label:MSSM/) {n=substr($0,11)+0; if (n % s == 0) keep=1}}
  /^end model/ {if (inb && keep) {printf "%s", buf} inb=0; keep=0}
' "$SRC" > "$W/models.txt"

LABELS=$(grep '^Label:' "$W/models.txt" | sed 's/Label://')
NM=$(echo $LABELS | wc -w)
echo "bF probe: $NM models (stride $STRIDE)"

{ for L in $LABELS; do
    echo "cd $L"; echo "cd spectrum"
    # `print summary` gives the FULL _l / _r listing (the gate's own k8 probe command);
    # `print(bF)` asks specifically for the right-chiral multiplets.  Both are issued so the
    # answer does not depend on one command's semantics.
    echo "print summary"
    echo "print(bF) with internal information"
    echo "cd .."; echo "cd .."
  done
  echo exit; echo yes; } > "$W/cmd.txt"

cd ~/rcharge_build
T0=$(date +%s)
timeout 7200 "$BIN" "$W/models.txt" < "$W/cmd.txt" > "$W/out.txt" 2>&1
T1=$(date +%s)
gzip -c "$W/out.txt" > "$DEST/bf_probe.txt.gz"
echo "bf probe exit=$? secs=$((T1-T0)) models=$NM bytes=$(wc -c < "$W/out.txt")"
rm -f "$W/out.txt"
echo DONE

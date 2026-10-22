#!/bin/bash
# su5_construction -- orbifolder driver for the CB-1 candidate embeddings.
#
# Same command sequence and the same STOCK geometry as the certified census driver
# code/z3xz6_sweep/run_census_z36.sh (spectrum prints BEFORE `analyze config`, because
# that command installs a new vev-configuration and would change the U(1) basis under
# which p_sh is reported).  One orbifolder invocation per segment.
#
# Usage: wsl -d Ubuntu -- bash <this file> <SEGID> <FIRST> <LAST>   (CB indices)
set -u
BASE=/mnt/f/physics/multi_geometry_campaign
SRC=$BASE/data/su5_construction/s5_models.txt
GEO=/mnt/f/physics/orbifolder-1.2.1/Geometry/Geometry_Z3xZ6_1_1.txt
DEST=$BASE/data/su5_construction/raw
LOGS=$BASE/data/su5_construction/logs
BIN=~/rcharge_build/src/orbifolder/orbifolder

SEG=$1; FIRST=$2; LAST=$3
W=~/cb1/seg$SEG
mkdir -p "$W" "$DEST" "$LOGS" ~/rcharge_build/Geometry
[ -s ~/rcharge_build/Geometry/"$(basename "$GEO")" ] || cp -f "$GEO" ~/rcharge_build/Geometry/

if [ -s "$DEST/cb_seg$SEG.txt.gz" ]; then
  echo "SKIP seg$SEG (already present)"; exit 0
fi

awk -v a="$FIRST" -v b="$LAST" '
  {sub(/\r$/, "")}
  /^begin model/ {buf=$0"\n"; inb=1; keep=0; next}
  inb {buf=buf $0"\n"
       if ($0 ~ /^Label:CB/) {n=substr($0,9)+0; if (n>=a && n<=b) keep=1}}
  /^end model/ {if (inb && keep) {printf "%s", buf} inb=0; keep=0}
' "$SRC" > "$W/models.txt"

LABELS=$(grep '^Label:' "$W/models.txt" | sed 's/Label://')
NM=$(echo $LABELS | wc -w)
echo "seg$SEG: models=$NM range=[$FIRST,$LAST]"
if [ "$NM" -eq 0 ]; then echo "EMPTY seg$SEG"; exit 1; fi

{ for L in $LABELS; do
    echo "cd $L"
    echo "cd spectrum"
    echo "print(F) with internal information"
    echo "cd .."
    echo "cd gauge group"
    echo "print gauge group"
    echo "print U1 generators"
    echo "print simple roots"
    echo "print anomaly info"
    echo "cd .."
    echo "cd vev-config"
    echo "analyze config print SU(5) simple roots"
    echo "use config(SMConfig)"
    echo "print gauge group"
    echo "cd .."
    echo "cd spectrum"
    echo "print summary"
    echo "cd .."
    echo "cd .."
  done
  echo exit; echo yes; } > "$W/cmd.txt"

cd ~/rcharge_build
T0=$(date +%s)
timeout 43200 "$BIN" "$W/models.txt" < "$W/cmd.txt" > "$W/out.txt" 2>&1
rc=$?
T1=$(date +%s)
gzip -c "$W/out.txt" > "$DEST/cb_seg$SEG.txt.gz"
echo "seg$SEG exit=$rc secs=$((T1-T0)) bytes=$(wc -c < "$W/out.txt") models=$NM" \
  | tee -a "$LOGS/cb_timing.log"
rm -f "$W/out.txt"
echo "DONE seg$SEG"

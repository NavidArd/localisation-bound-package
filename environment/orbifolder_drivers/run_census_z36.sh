#!/bin/bash
# z3xz6_sweep -- orbifolder CENSUS driver for the 4,612-model Z3xZ6_1_1 catalogue.
#
# ONE orbifolder invocation per SEGMENT of models (startup cost amortised), emitting
# per model, in this order and for the reason recorded in code/z6i/run_orbifolder_z6i.sh:
#
#   /spectrum>    print(F) with internal information     <- sector (k,l), fixed point,
#                                                           q_sh, p_sh, gamma, oscillators,
#                                                           representation, chirality
#   /gauge group> print gauge group / U1 generators / simple roots / anomaly info
#   /vev-config>  analyze config print SU(5) simple roots
#                 use config(SMConfig)   print gauge group
#   /spectrum>    print summary                          <- the SMConfig multiplicity table
#
# The spectrum prints come BEFORE `analyze config` because that command installs and
# selects a NEW vev-configuration and would change the U(1) basis under which p_sh is
# reported.  (Verbatim reason from the certified Z6-II / Z6-I drivers.)
#
# GEOMETRY: the STOCK Geometry_Z3xZ6_1_1.txt.  The massless spectrum and the gauge group
# are independent of declared discrete symmetries, so the census must NOT be run under a
# patched geometry -- the patch could otherwise shape its own evidence.  Coupling runs
# (run_couplings_z36.sh) use the _Rsym_prime patch; that split is the foundations pass's.
#
# Usage: wsl -d Ubuntu -- bash <this file> <SEGID> <FIRST> <LAST>
#        (FIRST/LAST are catalogue indices; labels are MSSM<idx>)
set -u
BASE=/mnt/f/physics/multi_geometry_campaign
SRC=$BASE/data/literature/parr2020_full/extracted/Z3xZ6_1_1.txt
GEO=/mnt/f/physics/orbifolder-1.2.1/Geometry/Geometry_Z3xZ6_1_1.txt
DEST=$BASE/data/z3xz6_sweep/raw/census
LOGS=$BASE/data/z3xz6_sweep/logs
BIN=~/rcharge_build/src/orbifolder/orbifolder

SEG=$1; FIRST=$2; LAST=$3
W=~/z36sw/seg$SEG
mkdir -p "$W" "$DEST" "$LOGS" ~/rcharge_build/Geometry
# NOTE: the geometry is staged ONCE by the launcher, never here.  Copying it
# per segment races sibling workers READING the same shared file: a worker can
# load a half-written geometry, print `"...txt" is corrupt` + `Orbifold "MSSMn"
# not known`, silently drop part of its model list, and still exit 0.  That cost
# 16 models of census segment 11 before it was caught by audit_raw_z36.sh.
[ -s ~/rcharge_build/Geometry/"$(basename "$GEO")" ] || cp -f "$GEO" ~/rcharge_build/Geometry/

if [ -s "$DEST/census_seg$SEG.txt.gz" ]; then
  echo "SKIP seg$SEG (already present)"; exit 0
fi

awk -v a="$FIRST" -v b="$LAST" '
  {sub(/\r$/, "")}
  /^begin model/ {buf=$0"\n"; inb=1; keep=0; next}
  inb {buf=buf $0"\n"
       if ($0 ~ /^Label:MSSM/) {n=substr($0,11)+0; if (n>=a && n<=b) keep=1}}
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
gzip -c "$W/out.txt" > "$DEST/census_seg$SEG.txt.gz"
echo "seg$SEG exit=$rc secs=$((T1-T0)) bytes=$(wc -c < "$W/out.txt") models=$NM" \
  | tee -a "$LOGS/census_timing.log"
rm -f "$W/out.txt"
echo "DONE seg$SEG"

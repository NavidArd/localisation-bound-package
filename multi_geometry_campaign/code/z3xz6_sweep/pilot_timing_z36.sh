#!/bin/bash
# z3xz6_sweep -- STEP-3 TIMING PILOT, and the wait()-sufficiency probe.
#
# `create coupling(F F F)` forks a child; `wait(N)` blocks the parent for up to N
# seconds.  The foundations pass used wait(600) and took 601 s wall -- i.e. the WAIT,
# not the computation, was the entire cost.  The pilot therefore has to establish two
# separate things, and it must never conflate them:
#
#   (a) the minimum wait at which the printed superpotential is BYTE-IDENTICAL to the
#       wait(600) reference (a sufficiency probe, not a speed measurement); and
#   (b) the per-model wall time of a segment at the chosen wait, over models drawn
#       from ACROSS the catalogue (not the first ten, which share provenance order).
#
# Silent truncation is the failure mode to fear: a too-short wait prints a PARTIAL W
# with no error.  So every production segment re-runs a checksum guard (see
# verify_coup_z36.py) and the chosen wait carries a stated margin over the measured
# completion time.
#
# Usage: wsl -d Ubuntu -- bash <this file>
set -u
BASE=/mnt/f/physics/multi_geometry_campaign
DEST=$BASE/data/z3xz6_sweep/raw/coup
RK=/tmp/rk_z36.sh
sed 's/\r$//' $BASE/code/z3xz6_sweep/run_couplings_z36.sh > $RK

echo "=== (a) wait-sufficiency probe on MSSM0 (reference: 908 terms at wait 600) ==="
for w in 1 2 3 5 10; do
  rm -f "$DEST/coup_segPW$w.txt.gz"
  T0=$(date +%s)
  bash $RK "PW$w" "$w" MSSM0 >/dev/null 2>&1
  T1=$(date +%s)
  n=$(zcat "$DEST/coup_segPW$w.txt.gz" | grep -A4000 '  W = ' | tr '+' '\n' | grep -c 'F_')
  h=$(zcat "$DEST/coup_segPW$w.txt.gz" | sed -n '/^  W = /,/^[[:space:]]*$/p' | md5sum | cut -c1-8)
  echo "  wait=$w  wall=$((T1-T0))s  terms=$n  md5=$h"
done

echo "=== (b) segment wall time, 10 models spread across the catalogue, wait=10 ==="
SPREAD="MSSM0 MSSM500 MSSM1000 MSSM1500 MSSM2000 MSSM2500 MSSM3000 MSSM3500 MSSM4000 MSSM4611"
rm -f "$DEST/coup_segPSPREAD.txt.gz"
T0=$(date +%s)
bash $RK PSPREAD 10 $SPREAD 2>&1 | tail -3
T1=$(date +%s)
echo "  spread segment wall=$((T1-T0))s for 10 models"
zcat "$DEST/coup_segPSPREAD.txt.gz" | grep -c '^  W = ' | sed 's/^/  W-blocks: /'
echo "DONE"

#!/bin/bash
# z3xz6_sweep -- the wait()-truncation guard.  Re-runs a 10-model spread segment at
# wait=120 and compares the superpotential block byte-for-byte against the wait=10 run.
# A too-short wait would print a PARTIAL W with no error; this is the only way to see it.
set -u
cd /mnt/f/physics/multi_geometry_campaign
D=data/z3xz6_sweep/raw/coup
SPREAD="MSSM0 MSSM500 MSSM1000 MSSM1500 MSSM2000 MSSM2500 MSSM3000 MSSM3500 MSSM4000 MSSM4611"
rm -f "$D/coup_segPSPREAD120.txt.gz"
T0=$(date +%s)
sed "s/\r$//" code/z3xz6_sweep/run_couplings_z36.sh > /tmp/rk_z36.sh; bash /tmp/rk_z36.sh PSPREAD120 120 $SPREAD 2>&1 | tail -2
T1=$(date +%s); echo "wall=$((T1-T0))s"
echo -n "wait10  W-md5: "; zcat "$D/coup_segPSPREAD.txt.gz"    | grep -A100000 "^  W = " | md5sum
echo -n "wait120 W-md5: "; zcat "$D/coup_segPSPREAD120.txt.gz" | grep -A100000 "^  W = " | md5sum

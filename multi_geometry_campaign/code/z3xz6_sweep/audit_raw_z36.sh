#!/bin/bash
# z3xz6_sweep -- integrity audit of the raw orbifolder segments.
#
# Catches the failure this sweep actually hit: every segment script copied the geometry
# file into the SHARED ~/rcharge_build/Geometry while sibling workers were READING it, so
# a worker could load a half-written file and reject part of its model list with
# `"...Geometry_Z3xZ6_1_1.txt" is corrupt` + `Orbifold "MSSMn" not known`.  The segment
# still exits 0 and still looks complete by byte count.  Only the loaded-model count and
# the corrupt/not-known strings reveal it.
set -u
D=${1:?dir}
PAT=${2:-census}
cd "$D"
for f in ${PAT}_seg*.txt.gz; do
  [ -e "$f" ] || continue
  t=$(zcat "$f" 2>/dev/null)
  n=$(printf '%s' "$t" | grep -m1 -oE '[0-9]+ orbifolds loaded' | grep -oE '^[0-9]+')
  c=$(printf '%s' "$t" | grep -c 'is corrupt')
  k=$(printf '%s' "$t" | grep -c 'not known')
  m=$(printf '%s' "$t" | grep -cE '^/> /MSSM[0-9]+>|/> /MSSM[0-9]+> ')
  echo "$f loaded=${n:-?} corrupt=$c notknown=$k blocks=$m"
done

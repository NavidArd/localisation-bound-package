#!/bin/bash
# z2xz4_sweep -- orbifolder COUPLING driver.
#
# ONE orbifolder invocation PER MODEL.  That is a correctness requirement, not a style
# choice: `create coupling` forks internally, and batching several models into one
# invocation makes the fork cascade emit MORE THAN ONE `W = ` print, one of which can be a
# spurious `W = 0`.  A model with != 1 `W = ` print is REJECTED by the parser rather than
# read, so the one-invocation-per-model discipline is what keeps that rejection rare.
#
# Commands per model:
#   /spectrum>            print(F) with internal information   <- so the ROLES can be
#                          re-derived inside the COUPLING run's own spectrum print and the
#                          agreement with the census becomes a CHECK, not an assumption
#   /couplings>           create coupling(F F F) allowed fields(<the model's Q, u^c and H_u
#                          role fields>) ; print superpotential
#
# THE `allowed fields(...)` RESTRICTION, and its scope, stated where it is used.
# Unrestricted `create coupling(F F F)` on these models tests C(n+2,3) ~ 2.2 million triples
# (n ~ 235 left-chiral fields) and costs ~540 s per model -- 236 core-hours over the 1,573
# models to be coupled, which the phase's scale discipline (guard 8) would have to stop and
# report.  orbifolder supports restricting the field pool: `allowed fields(...)` filters at
# the field-selection stage (`cyukawacouplings.cpp:779-830`) and every surviving triple is
# still subjected to the FULL gauge-invariance and space-group tests (`:1339-1421`).  The
# restriction used here is a ROLE restriction, NOT a block restriction: the pool is EVERY
# field the census gives the role Q, u^c or H_u in that model's physical orientation branch.
# Any up-type block on ANY channel, ANY leg assignment and ANY Higgs field therefore remains
# discoverable -- the completeness the sweep's guard 3 requires is preserved exactly.
# WHAT IS GIVEN UP, stated: the printed `W` is the superpotential RESTRICTED to that pool, so
# the term count is a restricted count and NO absence is claimed for couplings involving
# fields outside it (d^c, L, e^c, singlets, hidden matter).  Nothing downstream needs them.
#
# GEOMETRY: the STOCK Geometry_Z2xZ4_1_1.txt.
#   K6 SCOPE, stated where it is used: this geometry's plane orders are (4, 2, 4).  Plane 2
#   has order 2, which IS prime, so unlike Z4xZ4_1_1 a prime plane EXISTS here -- but plane 2
#   is never CAPABLE (the {3,4,6} letter-carrying theorem, selftest S2), so it carries no
#   classical exponent and the R-charge route it would license is not a route to a Yukawa
#   suppression statement.  Planes 1 and 3 are non-prime.  The stock geometry declares only
#   NON-R symmetries, so orbifolder's R-charge loop is INERT (CLAUDE.md 5.16 instance #2) and
#   its gamma rule is COMMENTED OUT upstream (instance #3); the coupling list orbifolder
#   returns is therefore a SUPERSET of the physical W, enforcing gauge invariance and the
#   SPACE-GROUP rule only.  Every absence statement derived downstream carries that scope.
#
# Usage: wsl -d Ubuntu -- bash <this file> <CHUNKID> <listfile: "LABEL<TAB>F_a F_b ...">
set -u
BASE=/mnt/f/physics/multi_geometry_campaign
SRC=$BASE/data/literature/parr2020_full/extracted/Z2xZ4_1_1.txt
GEO=/mnt/f/physics/orbifolder-1.2.1/Geometry/Geometry_Z2xZ4_1_1.txt
DEST=$BASE/data/z2xz4_sweep/raw/coup
LOGS=$BASE/data/z2xz4_sweep/logs
BIN=~/rcharge_build/src/orbifolder/orbifolder

CH=$1; LIST=$2
W=~/z24sw/coup$CH
mkdir -p "$W" "$DEST" "$LOGS" ~/rcharge_build/Geometry
[ -s ~/rcharge_build/Geometry/"$(basename "$GEO")" ] || cp -f "$GEO" ~/rcharge_build/Geometry/

if [ -s "$DEST/coup_seg$CH.txt.gz" ]; then
  echo "SKIP coup$CH (already present)"; exit 0
fi

: > "$W/all.txt"
T0=$(date +%s)
NM=0
while IFS=$'\t' read -r LAB FIELDS; do
  [ -z "$LAB" ] && continue
  awk -v want="Label:$LAB" '
    {sub(/\r$/, "")}
    /^begin model/ {buf=$0"\n"; inb=1; keep=0; next}
    inb {buf=buf $0"\n"; if ($0 == want) keep=1}
    /^end model/ {if (inb && keep) {printf "%s", buf} inb=0; keep=0}
  ' "$SRC" > "$W/m.txt"
  { echo "cd $LAB"
    echo "cd spectrum"
    echo "print(F) with internal information"
    echo "cd .."
    echo "cd couplings"
    if [ -n "${FIELDS:-}" ]; then
      echo "create coupling(F F F) allowed fields($FIELDS)"
    else
      echo "create coupling(F F F)"
    fi
    # `create coupling` FORKS (cprompt.cpp:5919-5975) and returns immediately.  Without an
    # explicit barrier the very next command runs in the PARENT while the child is still
    # computing, and `print superpotential` then prints `W = 0` for a model whose W is not
    # empty -- a silent false ABSENCE, the one error class this sweep cannot afford.  The
    # prompt's own `wait(n)` command (cprompt.cpp:616-650) blocks until every child PID
    # reports done; the parser additionally requires the `PID ... done:` marker and exactly
    # one `W = ` print per model, so a barrier failure is caught rather than absorbed.
    echo "wait(1)"
    echo "print superpotential"
    echo "cd .."
    echo "cd .."
    echo exit; echo yes; } > "$W/cmd.txt"
  cd ~/rcharge_build
  timeout 3600 "$BIN" "$W/m.txt" < "$W/cmd.txt" >> "$W/all.txt" 2>&1
  NM=$((NM + 1))
done < "$LIST"
T1=$(date +%s)
gzip -c "$W/all.txt" > "$DEST/coup_seg$CH.txt.gz"
echo "coup$CH models=$NM secs=$((T1-T0)) bytes=$(wc -c < "$W/all.txt")" \
  | tee -a "$LOGS/coup_timing.log"
rm -f "$W/all.txt"
echo "DONE coup$CH"

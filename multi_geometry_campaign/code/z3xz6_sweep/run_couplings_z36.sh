#!/bin/bash
# z3xz6_sweep -- orbifolder COUPLING driver for the `Z3xZ6_1_1` catalogue.
#
# GEOMETRY: data/z3xz6/patched_geometry/Geometry_Z3xZ6_1_1_Rsym_prime.txt, which declares
# ONLY the exact prime-plane law `sum (q_sh_2 + Osci_2) = -1 mod 3` (plane 2, the Z3
# plane -- the one plane of this geometry whose order is prime and whose R-charge law is
# therefore published, gamma-free and exact).  The two NON-prime planes 1 and 3 carry the
# unpublished two-generator gamma correction (RULE-GAP G-1) and are NOT declared here;
# they are applied OFFLINE with the three-branch gamma tag travelling.
#
# ONE ORBIFOLDER INVOCATION PER MODEL -- a correctness requirement, not a style choice.
# `create coupling(F F F)` FORKS a child that inherits the redirected stdin.  In a
# MULTI-model command stream that child goes on reading the SAME stream: it re-executes
# `wait`, re-prints the superpotential, and reaches the NEXT model's `create coupling`,
# forking again.  Three consequences, all measured on this catalogue:
#   (i)   every model printed its superpotential TWICE and the two prints were NOT always
#         equal -- 1 of 28 models in a wait=12 chunk printed `W = 0` first and the real
#         superpotential second, so a parser taking the FIRST print would have recorded a
#         SPURIOUS ABSENCE, the one error class this sweep cannot afford;
#   (ii)  the fork cascade reached ~7 live orbifolder processes per worker, parent and
#         descendants COMPETING for the same command stream;
#   (iii) throughput collapsed, because every surviving descendant also pays the `wait`.
# A single-model stream has no next `create coupling` for the child to reach, so the
# cascade cannot start and exactly one superpotential is printed.  (The foundations and
# gate specimens were single-model runs, which is why none of this was visible there, and
# why the multi-model shape inherited from the Z6-I driver looked validated.)
#
# The models of a chunk still run inside ONE segment script, appending to one raw file, so
# the output chunking is unchanged and one gz still covers ~50 models.  Only the
# orbifolder PROCESS boundary moved.
#
# SIZING FACT, measured rather than assumed: `wait(N)` costs N SECONDS ALWAYS -- it is a
# sleep, not a join.  What varies is only whether the forked child finished inside it,
# which orbifolder then reports as `done` (finished) or `stopped` (killed).  So the
# efficient protocol is a SHORT wait plus a re-run ladder over the `stopped` models, not a
# long wait for everyone: at wait=12 roughly half the catalogue reports `done` at ~50 s per
# model, and only the remainder needs a longer pass.
#
# Usage: wsl -d Ubuntu -- bash <this file> <SEGID> <WAIT> <LABEL...>
set -u
BASE=/mnt/f/physics/multi_geometry_campaign
SRC=$BASE/data/literature/parr2020_full/extracted/Z3xZ6_1_1.txt
GEO=$BASE/data/z3xz6/patched_geometry/Geometry_Z3xZ6_1_1_Rsym_prime.txt
DEST=$BASE/data/z3xz6_sweep/raw/coup
LOGS=$BASE/data/z3xz6_sweep/logs
BIN=~/rcharge_build/src/orbifolder/orbifolder

SEG=$1; WAITN=$2; shift 2
W=~/z36sw/coup$SEG
mkdir -p "$W" "$DEST" "$LOGS" ~/rcharge_build/Geometry
# NOTE: the geometry is staged ONCE by the launcher, never here.  Copying it per segment
# races sibling workers READING the same shared file: a worker can load a half-written
# geometry, print `"...txt" is corrupt` + `Orbifold "MSSMn" not known`, silently drop part
# of its model list, and still exit 0.  That cost 16 models of census segment 11 before an
# audit (audit_raw_z36.sh) caught it.
[ -s ~/rcharge_build/Geometry/"$(basename "$GEO")" ] || cp -f "$GEO" ~/rcharge_build/Geometry/

if [ -s "$DEST/coup_seg$SEG.txt.gz" ]; then
  echo "SKIP coup seg$SEG (already present)"; exit 0
fi

: > "$W/out.txt"
NM=0
cd ~/rcharge_build
T0=$(date +%s)
for M in "$@"; do
  NM=$((NM+1))
  awk -v want="Label:$M" '
    {sub(/\r$/, "")}
    /^begin model/ {buf=$0"\n"; inb=1; keep=0; next}
    inb {buf=buf $0"\n"; if ($0 == want) keep=1}
    /^end model/ {if (inb && keep) {printf "%s", buf} inb=0; keep=0}
  ' "$SRC" > "$W/model.txt"
  if [ ! -s "$W/model.txt" ]; then echo "EXTRACT-FAILED $M" >> "$W/out.txt"; continue; fi
  sed -i 's|Geometry_Z3xZ6_1_1\.txt|Geometry_Z3xZ6_1_1_Rsym_prime.txt|' "$W/model.txt"
  { echo "cd $M"; echo "cd spectrum"
    echo "print(F) with internal information"
    echo "cd .."; echo "cd couplings"
    echo "create coupling(F F F)"
    echo "wait($WAITN)"
    echo "print superpotential"
    echo "cd .."; echo "cd .."
    echo exit; echo yes; } > "$W/cmd.txt"
  timeout 1800 "$BIN" "$W/model.txt" < "$W/cmd.txt" >> "$W/out.txt" 2>&1
done
T1=$(date +%s)

NDONE=$(grep -o "/MSSM[0-9]*/couplings>" "$W/out.txt" 2>/dev/null | sort -u | wc -l)
NW=$(grep -c "^  W = " "$W/out.txt" 2>/dev/null || true)
# THE COMPLETION SIGNAL, and it supersedes any reasoning about `wait(N)`.
# orbifolder prints, when the forked coupling child finishes:
#     PID <n> done: "create coupling(F F F)" TIME: hh:mm:ss
#     <N> couplings created.
# A model WITHOUT that line whose `W = 0` is a TRUNCATION, not an absence -- measured:
# 4 of 10 probe models printed W = 0 at wait=12 with no `couplings created` line, while
# every model that DID carry the line had term count exactly equal to N.  A blind wait
# cannot distinguish the two; this line can, so the sweep gates on it rather than on the
# wait, and stragglers are re-run at a longer wait until they carry it.
NCREATED=$(grep -c "couplings created" "$W/out.txt" 2>/dev/null || true)
# COMPLETENESS GATE.  Published ONLY when every requested model reached `print
# superpotential` EXACTLY ONCE (NW == NDONE == NM).  NW != NM would mean the fork cascade
# is back.  An interrupted or partly-loaded run previously wrote a TRUNCATED gz that the
# resume logic then SKIPPED as if complete; a partial segment is now left unpublished and
# re-run whole, never repaired in place.
# Publication requires STRUCTURAL completeness (every requested model reached `print
# superpotential` exactly once).  It deliberately does NOT require every model to carry
# the completion signal: a 50-model chunk must not be discarded because one model was
# slow.  The per-model completion signal is enforced downstream by `coverage_z36.py`,
# which re-queues exactly the models that lack it -- so partial progress is kept and no
# unsignalled `W = 0` can reach a verdict.
# NM > 0 is part of the gate: an empty model list (a quoting accident in a caller)
# otherwise satisfies 0 == 0 == 0 and PUBLISHES an empty segment, which the resume
# logic then skips forever.  Observed once; guarded here.
if [ "$NM" -gt 0 ] && [ "$NDONE" -eq "$NM" ] && [ "$NW" -eq "$NM" ]; then
  gzip -c "$W/out.txt" > "$DEST/coup_seg$SEG.txt.gz.part"
  mv -f "$DEST/coup_seg$SEG.txt.gz.part" "$DEST/coup_seg$SEG.txt.gz"
  ST=PUBLISHED
else
  ST="INCOMPLETE-NOT-PUBLISHED"
fi
echo "coup seg$SEG exit=0 secs=$((T1-T0)) bytes=$(wc -c < "$W/out.txt") models=$NM modelsdone=$NDONE W-lines=$NW created=$NCREATED wait=$WAITN $ST" \
  | tee -a "$LOGS/coup_timing.log"
rm -f "$W/out.txt"
echo "DONE coup seg$SEG $ST"

#!/bin/bash
# z3xz6_sweep -- CALIBRATE the `wait(N)` needed for a SINGLE-model coupling run.
#
# With one model per orbifolder invocation there is no forked-child second print to act
# as a safety net, so `wait(N)` alone decides whether `print superpotential` sees a
# finished computation.  A too-short wait prints `W = 0` with no error -- a SPURIOUS
# ABSENCE.  MSSM255 is the known witness: it printed 0 at wait=12 and a real
# superpotential when given more time inside a multi-model stream.
#
# The calibration reports, per probe model, the term count as a function of the wait.  The
# production wait must sit safely above the point where every probe has STOPPED GROWING.
#
# Usage: wsl -d Ubuntu -- bash <this file> [WAITS...]
set -u
BASE=/mnt/f/physics/multi_geometry_campaign
RK=/tmp/rk_cal_z36.sh
sed 's/\r$//' "$BASE/code/z3xz6_sweep/run_couplings_z36.sh" > $RK
GEO=$BASE/data/z3xz6/patched_geometry/Geometry_Z3xZ6_1_1_Rsym_prime.txt
mkdir -p ~/rcharge_build/Geometry; cp -f "$GEO" ~/rcharge_build/Geometry/
D=$BASE/data/z3xz6_sweep/raw/coup

PROBES="MSSM255 MSSM1000 MSSM1500 MSSM2500 MSSM0 MSSM4569"
for WN in "${@:-15 30 60 120 240}"; do
  rm -f "$D/coup_segCAL$WN.txt.gz"
  T0=$(date +%s)
  bash $RK "CAL$WN" "$WN" $PROBES > /dev/null 2>&1
  T1=$(date +%s)
  echo -n "wait=$WN  wall=$((T1-T0))s  "
  python3 - "$D/coup_segCAL$WN.txt.gz" <<'PY'
import gzip, re, sys
t = gzip.open(sys.argv[1], "rt", errors="replace").read()
MARK = re.compile(r"/>\s*/(MSSM\d+)>")
m = [(x.start(), x.group(1)) for x in MARK.finditer(t)]
out = []
for i, (p, lab) in enumerate(m):
    b = t[p:(m[i+1][0] if i+1 < len(m) else len(t))]
    w = re.search(r"\n\s+W\s*=\s*(.*?)(?:\n\s*\n|\n/)", b, re.S)
    n = 0
    if w:
        n = sum(1 for x in " ".join(w.group(1).split()).split("+")
                if re.search(r"F_\d+", x))
    out.append("%s:%d" % (lab, n))
print(" ".join(out))
PY
done

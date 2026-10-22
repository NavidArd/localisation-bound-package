#!/bin/bash
# z3xz6_sweep -- PRODUCTION truncation guard for the coupling runs.
#
# `wait(N)` bounds how long the parent waits for the forked `create coupling` child; a
# too-short wait prints a PARTIAL superpotential with NO error, which is precisely the
# silent-loss shape this sweep must not tolerate.  The pilot established byte-identity at
# wait in {1,2,3,5,10} against the wait(600) reference on MSSM0 and at wait=10 vs wait=120
# on a 10-model spread.  Production ran at wait=12.
#
# This script re-runs a RANDOM SAMPLE of production models at wait=120 and compares their
# superpotential blocks byte for byte against what production published.  A single
# mismatch invalidates the wait and the affected segments must be re-run.
#
# Usage: wsl -d Ubuntu -- bash <this file> [NSAMPLE] [SEED]
set -u
BASE=/mnt/f/physics/multi_geometry_campaign
D=$BASE/data/z3xz6_sweep/raw/coup
N=${1:-24}
SEED=${2:-20260818}
RK=/tmp/rk_ver_z36.sh
sed 's/\r$//' "$BASE/code/z3xz6_sweep/run_couplings_z36.sh" > $RK
GEO=$BASE/data/z3xz6/patched_geometry/Geometry_Z3xZ6_1_1_Rsym_prime.txt
mkdir -p ~/rcharge_build/Geometry; cp -f "$GEO" ~/rcharge_build/Geometry/

MODELS=$(awk -v n="$N" -v s="$SEED" 'BEGIN{srand(s); for(i=0;i<n;i++) printf "MSSM%d\n", int(rand()*4612)}' | sort -u)
echo "verification sample: $(echo $MODELS | wc -w) models at wait=120"

rm -f "$D/coup_segVERIFY.txt.gz"
bash $RK VERIFY 120 $MODELS 2>&1 | tail -2

# production W blocks for the sampled models, extracted from the published segments
python3 - "$D" $MODELS <<'PY'
import gzip, os, re, sys, hashlib
D = sys.argv[1]; want = set(sys.argv[2:])
MARK = re.compile(r"/>\s*/(MSSM\d+)>")
def blocks(path):
    with gzip.open(path, "rt", errors="replace") as fh:
        t = fh.read()
    m = [(x.start(), x.group(1)) for x in MARK.finditer(t)]
    return {lab: t[p:(m[i+1][0] if i+1 < len(m) else len(t))]
            for i, (p, lab) in enumerate(m)}
def wpart(b):
    i = b.find("\n  W = ")
    return b[i:] if i >= 0 else ""
prod = {}
for fn in sorted(os.listdir(D)):
    if not (fn.startswith("coup_seg") and fn.endswith(".gz")) or "VERIFY" in fn:
        continue
    for lab, b in blocks(os.path.join(D, fn)).items():
        if lab in want:
            prod[lab] = wpart(b)
ver = {k: wpart(v) for k, v in blocks(os.path.join(D, "coup_segVERIFY.txt.gz")).items()
       if k in want}
bad = [k for k in sorted(ver) if prod.get(k) != ver[k]]
missing = sorted(want - set(prod))
print("compared %d models | identical %d | MISMATCH %d | not in production %d"
      % (len(ver), len(ver) - len(bad), len(bad), len(missing)))
if bad:
    print("MISMATCHING MODELS (the production wait is INSUFFICIENT):", bad)
if missing:
    print("models absent from the published segments:", missing)
PY

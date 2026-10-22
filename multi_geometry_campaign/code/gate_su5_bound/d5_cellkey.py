"""GATE gate_su5_bound D5 -- CELL-KEY COMPLETENESS (CLAUDE.md 5.22 first addendum (ii)).

The key used by verify.py and by the certified gate is `(sector (k,l), raw printed
`fixed point n_a`)` -- six lattice integers, verbatim.  Two failure modes:

  MERGE  two physically distinct fixed points share one key  -> FAKES multiplicity
         (a-fortiori SAFE for an absence claim: it can only create doublings, and none
          were seen)
  SPLIT  one physical fixed point is written under two keys   -> HIDES a doubling
         (a-fortiori UNSAFE: this is the direction that could kill the absence claim)

A SPLIT happens iff two DIFFERENT printed n_a give the same PHYSICAL localisation, i.e.
the same space-group conjugacy class.  The campaign's own certified reading of a printed
n_a is f_a = (1 - R_a)^{-1} n|_a mod 1 per plane, with cells identified up to the point
group ORBIT (CLAUDE.md 5.22: orbit form, NOT fixed-point equality).  This script:

  1. harvests EVERY distinct (sector, raw n_a) actually printed in the full raw census of
     both geometries, by its own regex over the gzipped text;
  2. maps each to the per-plane ORBIT tuple with the certified geometry machinery
     (z24sw_lib.leg_orbit_index / z36sw_lib.leg_orbit);
  3. tests INJECTIVITY of raw key -> orbit tuple, per sector.  A non-injective sector is a
     SPLIT and is listed in full;
  4. reports, for every sector, what is printed in the coordinates of an UNROTATED plane
     (the N=2 direction), which is where a spurious extra key would come from.

Outputs data/gate_su5_bound/d5_cellkey.json
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"

import glob                                                            # noqa: E402
import gzip                                                            # noqa: E402
import json                                                            # noqa: E402
import multiprocessing as mp                                           # noqa: E402
import re                                                              # noqa: E402
import sys                                                             # noqa: E402
from collections import Counter, defaultdict                           # noqa: E402
from fractions import Fraction as F                                    # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CAMP = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(CAMP, "data", "gate_su5_bound")

GEOM = {
    "Z2xZ4_1_1": dict(lib="z2xz4_sweep", mod="z24sw_lib", M=2, N=4,
                      v1=(F(0), F(1, 2), F(-1, 2)), v2=(F(1, 4), F(0), F(-1, 4)),
                      census=os.path.join(CAMP, "data", "z2xz4_sweep", "raw", "census")),
    "Z3xZ6_1_1": dict(lib="z3xz6_sweep", mod="z36sw_lib", M=3, N=6,
                      v1=(F(0), F(1, 3), F(-1, 3)), v2=(F(1, 6), F(0), F(-1, 6)),
                      census=os.path.join(CAMP, "data", "z3xz6_sweep", "raw", "census")),
}

_PAT = re.compile(r"sector \(k,l\)\s*:\s*\(([^)]*)\)\s*\n\s*"
                  r"fixed point n_a\s*:\s*\(([^)]*)\)")


def harvest(args):
    path, geom = args
    G = GEOM[geom]
    txt = gzip.open(path, "rt", errors="replace").read()
    c = Counter()
    for m in _PAT.finditer(txt):
        kl = tuple(int(x) for x in m.group(1).split(","))
        kl = (kl[0] % G["M"], kl[1] % G["N"])
        fp = tuple(int(x) for x in m.group(2).split(","))
        c[(kl, fp)] += 1
    return {"%s|%s" % (kl, fp): n for (kl, fp), n in c.items()}


def main():
    res = {}
    for geom in ("Z2xZ4_1_1", "Z3xZ6_1_1"):
        G = GEOM[geom]
        files = sorted(glob.glob(os.path.join(G["census"], "*.txt.gz")))
        with mp.Pool(6) as pool:
            parts = pool.map(harvest, [(f, geom) for f in files])
        tot = Counter()
        for p in parts:
            tot.update(p)
        sys.path.insert(0, os.path.join(CAMP, "code", G["lib"]))
        Z = __import__(G["mod"])

        bysec = defaultdict(dict)
        for k, n in tot.items():
            kl_s, fp_s = k.split("|")
            kl = eval(kl_s)                                            # noqa: S307
            fp = eval(fp_s)                                            # noqa: S307
            bysec[kl][fp] = n

        secs = {}
        splits_total = 0
        for kl in sorted(bysec):
            fps = bysec[kl]
            phi = [(kl[0] * G["v1"][i] + kl[1] * G["v2"][i]) % 1 for i in range(3)]
            unrot = [i + 1 for i in range(3) if phi[i] == 0]
            orb = {}
            err = None
            for fp in fps:
                try:
                    if geom == "Z2xZ4_1_1":
                        t = tuple(Z.leg_orbit_index(kl, fp, a) for a in (1, 2, 3))
                    else:
                        t = tuple((Z.leg_orbit(kl, fp, a) or (None,))[0]
                                  for a in (1, 2, 3))
                except Exception as e:                                 # noqa: BLE001
                    err = "%s: %s" % (type(e).__name__, e)
                    t = None
                orb[fp] = t
            inv = defaultdict(list)
            for fp, t in orb.items():
                inv[t].append(fp)
            collide = {str(t): [list(f) for f in v] for t, v in inv.items() if len(v) > 1}
            splits_total += sum(len(v) - 1 for v in inv.values() if len(v) > 1)
            # what is printed in the coordinates of an unrotated plane?
            unrot_vals = {}
            for a in unrot:
                unrot_vals["plane%d" % a] = sorted(
                    {(fp[2 * (a - 1)], fp[2 * (a - 1) + 1]) for fp in fps})
            secs[str(kl)] = dict(
                n_raw_keys=len(fps),
                n_orbit_tuples=len({t for t in orb.values()}),
                phi=[str(x) for x in phi],
                unrotated_planes=unrot,
                unrotated_plane_printed_values={k: [list(x) for x in v]
                                                for k, v in unrot_vals.items()},
                injective=(len(fps) == len({t for t in orb.values()})),
                colliding_orbit_tuples=collide,
                orbit_map={str(list(fp)): str(t) for fp, t in sorted(orb.items())},
                error=err,
                field_records=sum(fps.values()))
        res[geom] = dict(sectors=secs, total_raw_keys=sum(len(v) for v in bysec.values()),
                         split_excess_keys=splits_total,
                         any_split=splits_total > 0)
        print("=====", geom, "split_excess_keys=", splits_total)
        for s, d in secs.items():
            print("  ", s, "raw=%d orbits=%d inj=%s unrot=%s unrotvals=%s"
                  % (d["n_raw_keys"], d["n_orbit_tuples"], d["injective"],
                     d["unrotated_planes"], d["unrotated_plane_printed_values"]),
                  ("ERR " + d["error"]) if d["error"] else "")
            if d["colliding_orbit_tuples"]:
                print("      COLLISIONS:", json.dumps(d["colliding_orbit_tuples"])[:600])
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "d5_cellkey.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=str)


if __name__ == "__main__":
    main()

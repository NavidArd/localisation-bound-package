"""GATE su5_ten_bound -- D1/D2: the Delta = 0 hole.

The theorem's second line is "Delta = p - p'' is a NON-ZERO vector of Lambda".  The proof
justifies non-vanishing with: "The massless shell is partitioned among the printed 4D
irreps, so p'' != p."  That justification covers two states distinguished by their GAUGE
momentum.  It does NOT cover two states with the SAME p_sh distinguished by their
OSCILLATOR content -- and such states exist.

In a twisted sector with twist eta = (0, 1/6, 5/6) the available left-moving oscillators
are  alphabar^2 at mode 1/6, alpha^3 at mode 1/6, alpha^2 and alphabar^3 at mode 5/6.  The
two mode-1/6 oscillators contribute  (N_L - Nbar_L) = (0,-1,0)  and  (0,0,+1); dotted into
v = (0, 1/6, -1/6) BOTH give -1/6, and dotted into the other generator's twist
v = (0, 1/3, -1/3) BOTH give -1/3.  Their orbifold projection phases are therefore
IDENTICAL: if one survives, so does the other.  Two distinct massless left-chiral states,
same fixed point, same p_sh, same Standard-Model role, Delta = 0 -- and the theorem is
silent, because its inequality is vacuous when Delta vanishes.

The hole is not empty by budget arithmetic either.  A 10-role needs 2B >= |w|^2 = 6/5, and
on Z3xZ6_1_1 the top budget 31/36 minus the minimal oscillator mode 1/6 leaves B = 25/36,
i.e. R = 17/90 >= 0 -- such fields EXIST in the catalogue.  (On Z2xZ4_1_1 the minimal mode
is 1/4 and 13/16 - 1/4 = 9/16 < 3/5, so no oscillator 10-role can exist there at all: the
hole is empty on that geometry for arithmetic reasons, which is worth stating.)

This script decides the question empirically, on the full catalogues:

  T1  at every twisted cell, are there TWO OR MORE left-chiral records sharing one p_sh?
      If yes: their roles, oscillator strings, gamma phases and right-movers are printed.
  T2  restricted to records carrying a Standard-Model role, and separately to the 10-roles.
  T3  the arithmetic feasibility of an oscillator 10-role per geometry.

Output: data/gate_su5_bound/d2_pshdup_<geom>.json
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"

import glob                                                             # noqa: E402
import gzip                                                             # noqa: E402
import json                                                             # noqa: E402
import multiprocessing as mp                                            # noqa: E402
import sys                                                              # noqa: E402
from collections import Counter, defaultdict                            # noqa: E402
from fractions import Fraction as F                                     # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CAMP = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(CAMP, "code", "su5_ten_bound"))
OUT = os.path.join(CAMP, "data", "gate_su5_bound")

ROLES10 = {"Q", "uc", "ec", "Qbar", "ucbar", "ecbar"}
ROLES5 = {"dc", "L", "Hu", "dcbar"}

GEOM = {
    "Z2xZ4_1_1": dict(lib="z2xz4_sweep", mod="z24sw_lib", M=2, N=4,
                      census=os.path.join(CAMP, "data", "z2xz4_sweep", "raw", "census")),
    "Z3xZ6_1_1": dict(lib="z3xz6_sweep", mod="z36sw_lib", M=3, N=6,
                      census=os.path.join(CAMP, "data", "z3xz6_sweep", "raw", "census")),
}


def worker(args):
    path, geom = args
    G = GEOM[geom]
    sys.path.insert(0, os.path.join(CAMP, "code", G["lib"]))
    Z = __import__(G["mod"])
    acc = dict(models=0, cells=0, recs=0, skipped=Counter(),
               dup_groups=0, dup_groups_role=0, dup_groups_10=0,
               dup_size=Counter(), examples=[], examples10=[],
               osc10=Counter(), osc_by_fam=Counter())
    txt = gzip.open(path, "rt", errors="replace").read()
    for lab, blk in Z.split_models(txt):
        pre = blk.split("Analyzing")[0]
        try:
            if geom == "Z2xZ4_1_1":
                spec = Z.parse_spectrum(pre, M=G["M"], N=G["N"])
                factors = Z.parse_gauge_factors(blk)
                cfgs = Z.parse_smconfig(blk)
                if not cfgs or cfgs[0][0] is None or not spec:
                    acc["skipped"]["no_smconfig"] += 1
                    continue
                tY, tbl = cfgs[0]
                printed = defaultdict(int)
                for row in tbl:
                    printed[(row[0], row[1], row[2])] += row[3]
                pos = Z.observable_positions(spec, tY, factors, dict(printed))
                if pos is None:
                    acc["skipped"]["closure_failed"] += 1
                    continue
                roles, _ = Z.role_census(spec, tY, *pos)
            else:
                spec = Z.parse_spectrum(pre)
                gauge = Z.parse_gauge(blk)
                if gauge["tY"] is None or not spec:
                    acc["skipped"]["no_hypercharge"] += 1
                    continue
                pos = Z.observable_positions(spec, gauge["tY"], gauge)
                if pos is None:
                    acc["skipped"]["closure_failed"] += 1
                    continue
                a3, a2, sY0, tY, bi = pos
                roles, _ = Z.role_census(spec, tY, a3, a2, sY0)
        except Exception as e:                                          # noqa: BLE001
            acc["skipped"]["parse:%s" % type(e).__name__] += 1
            continue
        acc["models"] += 1
        cells = defaultdict(list)
        for l, r in spec.items():
            if r["p_sh"] is None or r["chir"] != "l":
                continue
            cells[(tuple(r["kl"]), tuple(r["fp"]) if r["fp"] else None)].append(r)
        for key, recs in cells.items():
            (kl, fp) = key
            if tuple(kl) == (0, 0):
                continue
            acc["cells"] += 1
            acc["recs"] += len(recs)
            byp = defaultdict(list)
            for r in recs:
                byp[tuple(r["p_sh"])].append(r)
                ro = roles.get(r["label"], {}).get("role")
                osc = (r.get("osc") or r.get("osc_raw") or "no").strip()
                if ro in ROLES10:
                    acc["osc10"][osc] += 1
                    acc["osc_by_fam"]["10|%s" % ("osc" if osc not in ("no", "") else "none")] += 1
                elif ro in ROLES5:
                    acc["osc_by_fam"]["5bar|%s" % ("osc" if osc not in ("no", "") else "none")] += 1
            for psh, rs in byp.items():
                if len(rs) < 2:
                    continue
                acc["dup_groups"] += 1
                acc["dup_size"][len(rs)] += 1
                ros = [roles.get(r["label"], {}).get("role") for r in rs]
                if any(x is not None for x in ros):
                    acc["dup_groups_role"] += 1
                if any(x in ROLES10 for x in ros):
                    acc["dup_groups_10"] += 1
                    if len(acc["examples10"]) < 20:
                        acc["examples10"].append(dict(
                            model=lab, kl=list(kl), fp=list(fp) if fp else None,
                            roles=ros, n=len(rs),
                            p_sh=[str(x) for x in psh],
                            osc=[(r.get("osc") or r.get("osc_raw")) for r in rs],
                            gamma=[[str(x) for x in (r.get("gamma") or ())] for r in rs],
                            qsh=[[str(x) for x in (r.get("qsh") or ())] for r in rs],
                            rep=[list(r["rep"]) for r in rs],
                            file=os.path.basename(path)))
                if len(acc["examples"]) < 20:
                    acc["examples"].append(dict(
                        model=lab, kl=list(kl), fp=list(fp) if fp else None,
                        roles=ros, n=len(rs),
                        osc=[(r.get("osc") or r.get("osc_raw")) for r in rs],
                        gamma=[[str(x) for x in (r.get("gamma") or ())] for r in rs],
                        rep=[list(r["rep"]) for r in rs],
                        file=os.path.basename(path)))
    for k in ("skipped", "dup_size", "osc10", "osc_by_fam"):
        acc[k] = dict(acc[k])
    return acc


def main():
    geom = sys.argv[1] if len(sys.argv) > 1 else "Z2xZ4_1_1"
    nseg = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    files = sorted(glob.glob(os.path.join(GEOM[geom]["census"], "*.txt.gz")))
    if nseg:
        files = files[:nseg]
    with mp.Pool(6) as pool:
        accs = pool.map(worker, [(f, geom) for f in files])
    out = dict(models=0, cells=0, recs=0, dup_groups=0, dup_groups_role=0,
               dup_groups_10=0)
    ctr = {k: Counter() for k in ("skipped", "dup_size", "osc10", "osc_by_fam")}
    ex, ex10 = [], []
    for a in accs:
        for k in out:
            out[k] += a[k]
        for k in ctr:
            ctr[k].update(a[k])
        ex += a["examples"][:5]
        ex10 += a["examples10"][:5]
    for k in ctr:
        out[k] = dict(ctr[k])
    out["examples"] = ex[:40]
    out["examples10"] = ex10[:40]
    out["geometry"] = geom
    with open(os.path.join(OUT, "d2_pshdup_%s.json" % geom), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print(json.dumps({k: v for k, v in out.items()
                      if k not in ("examples", "examples10")}, indent=1, default=str))
    print("\nDUP EXAMPLES (any):")
    for e in out["examples"][:5]:
        print("  ", e)
    print("DUP EXAMPLES (10-role):")
    for e in out["examples10"][:5]:
        print("  ", e)


if __name__ == "__main__":
    main()

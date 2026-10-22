"""root_lift QUESTION A -- does the centraliser of a constructing element see the whole
gauge-shift set?

For every twisted constructing element g of every campaign geometry, compute
Z_S(g) = {h : gh = hg}, the sublattice L_g subset Z^8 spanned by its (k,l,n) labels, and
decide whether Z^8 = L_g + ker_geo (convention C6 of rl_lib).  Index 1 means the
projection at g constrains Delta against EVERY gauge-shift generator, i.e.

    Delta a local root at g   =>   Delta a root of the unbroken 4d gauge group.

Index > 1 names exactly which shift/Wilson-line generators escape the centraliser, and the
lift is then CONDITIONAL on Delta pairing integrally with those.

Output: data/root_lift/questionA_centralizer.json
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CAMP = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)

import rl_lib as R                                                      # noqa: E402

OUT = os.path.join(CAMP, "data", "root_lift")
os.makedirs(OUT, exist_ok=True)


def main():
    guards = R.selftest(verbose=True)
    print()
    out = {"guards_1a_1b_1c": guards, "geometries": {}}
    for lab in R.CAMPAIGN_GEOMS:
        S = R.SpaceGroup(lab)
        g = S.G
        geo = {"N": S.N, "M": S.M, "wl_orders": list(g.wl_orders),
               "identical_wls": [list(x) for x in g.identical_wl_groups],
               "sectors": {}}
        for sec in g.sectors_ce:
            k, l = sec[0][0], sec[0][1]
            if (k, l) == (0, 0):
                continue
            key = "(%d, %d)" % (k, l)
            rows = []
            for ce in sec:
                v = S.lift_verdict(ce)
                rows.append(dict(ce=list(ce), index=v["index"],
                                 missing=v["missing_generators"],
                                 pg_image=["%d,%d" % t for t in v["centralizer_pointgroup"]],
                                 n_pure_transl_dim=v["n_pure_translations"]))
            idxs = sorted({r["index"] for r in rows})
            geo["sectors"][key] = dict(
                n_ce=len(sec),
                unrotated_planes=g.unrotated_planes(k, l),
                indices=idxs,
                lift_holds_all_ce=(idxs == [1]),
                lift_holds_n_ce=sum(1 for r in rows if r["index"] == 1),
                missing_union=sorted({m for r in rows for m in r["missing"]}),
                rows=rows)
        out["geometries"][lab] = geo
        nsec = len(geo["sectors"])
        nall = sum(1 for s in geo["sectors"].values() if s["lift_holds_all_ce"])
        nce = sum(s["n_ce"] for s in geo["sectors"].values())
        nce1 = sum(s["lift_holds_n_ce"] for s in geo["sectors"].values())
        print("%-12s  sectors %2d  lift-on-all-ce %2d   ce %4d  lift-ce %4d"
              % (lab, nsec, nall, nce, nce1))

    with open(os.path.join(OUT, "questionA_centralizer.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("\nwrote", os.path.join(OUT, "questionA_centralizer.json"))

    # ---- the paper's scope, printed explicitly
    print("\n=== PAPER A SCOPE: Z3xZ6_1_1 sectors T(0,1) and T(0,5) ===")
    S = R.SpaceGroup("Z3xZ6_1_1")
    for sec in S.G.sectors_ce:
        k, l = sec[0][0], sec[0][1]
        if (k, l) not in [(0, 1), (0, 5)]:
            continue
        for ce in sec:
            v = S.lift_verdict(ce)
            print("  ce=%s  index=%s  missing=%s  Z_S(g)->P image size=%d  pure-transl rank=%d"
                  % (list(ce), v["index"], v["missing_generators"],
                     len(v["centralizer_pointgroup"]), v["n_pure_translations"]))
            print("     kernel basis (pure translations in Z_S(g)):", v["kern_basis"])


if __name__ == "__main__":
    main()

"""V5/V8 -- gauge-group census over all 2,912 CB models.

The phase's funnel row reads "observable gauge group exactly SU(3) x SU(2) x U(1)_Y :
2,912 (by construction, A6)".  A6 kills only the twelve X/Y roots of the block-1 SU(5);
it says nothing about block-1 roots OUTSIDE that SU(5) (those involving coordinates
6,7,8 and the half-integer roots).  This script measures what orbifolder actually prints,
and checks the conditions that are actually load-bearing:

  (a) the first-E8 factor list as printed;
  (b) SU(3)_c = <e1-e2, e2-e3> and SU(2)_L = <e4-e5> are among the 4d simple roots;
  (c) no 4d root joins SU(3)_c or SU(2)_L to anything else (they are whole factors);
  (d) t_Y is orthogonal to EVERY printed simple root and lies in the span of the printed
      U(1) generators (so U(1)_Y really is a U(1) of the printed group);
  (e) no SU(5) among the 4d roots (the X/Y roots are dead).

Output: data/gate_su5_construction/g5b_gauge.json
"""
from __future__ import annotations

import glob
import gzip
import json
import os
import sys
from collections import Counter
from fractions import Fraction as F
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gs_lib as G
from g5_battery import solve_span

RAW = os.path.join(G.PHASE_DATA, "raw")
AC = [G.A3[0], G.A3[1], G.A2L]
XY = [r for r in G.e8xe8_roots()
      if all(r[k] == 0 for k in range(5, 16))
      and sum(1 for k in range(3) if r[k] != 0) == 1
      and sum(1 for k in (3, 4) if r[k] != 0) == 1]


def do(path):
    rows = []
    txt = gzip.open(path, "rt", errors="replace").read()
    for lab, blk in G.split_models(txt):
        gg = G.parse_gauge_group(blk)
        simple = G.parse_simple_roots(blk)
        u1 = G.parse_u1_generators(blk)
        if gg is None:
            continue
        first = gg.split(" and ")[0]
        rest = gg.split(" and ")[1:]
        sset = set(tuple(s) for s in simple)
        has_sm = all(tuple(a) in sset for a in AC)
        # (c) SU(3)_c / SU(2)_L are whole factors: no OTHER simple root is
        # non-orthogonal to them
        whole = True
        for a in AC:
            for s in simple:
                if tuple(s) in set(tuple(x) for x in AC):
                    continue
                if G.dot(a, s) != 0:
                    whole = False
        orth = all(G.dot(G.TY, s) == 0 for s in simple)
        span = solve_span(u1, G.TY) is not None if u1 else False
        xy_alive = sum(1 for r in XY if tuple(r) in sset)
        rows.append(dict(model=lab, gg=gg, first_e8=first,
                         n_U1=len(u1), n_simple=len(simple),
                         sm_roots_present=has_sm, sm_factors_whole=whole,
                         tY_orth=orth, tY_in_U1span=span,
                         xy_simple_alive=xy_alive))
    return rows


def main():
    paths = sorted(glob.glob(os.path.join(RAW, "*.txt.gz")))
    with Pool(6) as pool:
        rows = [r for c in pool.map(do, paths) for r in c]
    first = Counter(r["first_e8"] for r in rows)
    summary = dict(
        models=len(rows),
        first_E8_factor_histogram=dict(first),
        exactly_SU3xSU2_in_first_E8=sum(1 for r in rows
                                        if r["first_e8"] == "SU(3) x SU(2)"),
        more_than_SU3xSU2_in_first_E8=sum(1 for r in rows
                                          if r["first_e8"] != "SU(3) x SU(2)"),
        sm_roots_present=sum(1 for r in rows if r["sm_roots_present"]),
        sm_factors_whole=sum(1 for r in rows if r["sm_factors_whole"]),
        tY_orth_all_simple=sum(1 for r in rows if r["tY_orth"]),
        tY_in_printed_U1_span=sum(1 for r in rows if r["tY_in_U1span"]),
        any_XY_simple_root_alive=sum(1 for r in rows if r["xy_simple_alive"]),
        n_U1_hist=dict(Counter(r["n_U1"] for r in rows)),
    )
    print(json.dumps(summary, indent=1))
    json.dump(dict(summary=summary, rows=rows),
              open(os.path.join(G.OUT, "g5b_gauge.json"), "w"), indent=1)


if __name__ == "__main__":
    main()

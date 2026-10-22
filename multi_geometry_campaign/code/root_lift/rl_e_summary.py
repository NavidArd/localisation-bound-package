"""Totals over the per-sector lift-index table, so the headline is read and not transcribed.

`rl_a_centralizer.py` writes the lift index of every twisted constructing element of every
campaign geometry, sector by sector. The paper quotes two totals from it -- how many twisted
constructing elements there are, and at how many of them the local-to-four-dimensional root
lift is FORCED by the geometry -- and this script adds them up so that a claim row can point
at a stored number rather than at prose.

Input   ../../data/root_lift/questionA_centralizer.json
Output  ../../data/root_lift/questionA_summary.json

The per-geometry rows are kept, because the negative half of the statement is the important
half: on Z2xZ4_1_1, this paper's other catalogue, the lift is forced nowhere.
"""
from __future__ import annotations

import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "..", "data", "root_lift"))


def main() -> int:
    with open(os.path.join(DATA, "questionA_centralizer.json"), encoding="utf-8") as fh:
        doc = json.load(fh)

    per = {}
    n_ce = n_lift = 0
    for geom, g in doc["geometries"].items():
        ce = lift = 0
        sectors_all = []
        for sec, s in g["sectors"].items():
            ce += s["n_ce"]
            lift += s["lift_holds_n_ce"]
            if s.get("lift_holds_all_ce"):
                sectors_all.append(sec)
        per[geom] = {
            "twisted_constructing_elements": ce,
            "lift_forced_at": lift,
            "sectors_where_lift_holds_at_every_element": sectors_all,
        }
        n_ce += ce
        n_lift += lift

    out = {
        "per_geometry": per,
        "n_ce": n_ce,
        "n_lift": n_lift,
        "scope": ("The signature is claimed only where the lift is forced by the geometry "
                  "AND that element is the sector's only one: T(0,1) and T(0,5) of "
                  "Z3xZ6_1_1. On Z2xZ4_1_1 the lift is forced nowhere and no statement of "
                  "the signature's form may be made."),
    }

    path = os.path.join(DATA, "questionA_summary.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")

    for geom, row in sorted(per.items(), key=lambda kv: -kv[1]["lift_forced_at"]):
        print("%-12s  %3d elements, lift forced at %d  %s"
              % (geom, row["twisted_constructing_elements"], row["lift_forced_at"],
                 ",".join(row["sectors_where_lift_holds_at_every_element"])))
    print("TOTAL         %3d elements, lift forced at %d" % (n_ce, n_lift))
    print("wrote", path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

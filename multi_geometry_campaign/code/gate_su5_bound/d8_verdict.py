"""GATE d8d -- merge the gate's row-by-row verdicts on 7.1 into d8_classwide.json."""

from __future__ import annotations

# --- reproduction-package root shim -------------------------------------------------
# The line above originally read  F:\physics  -- the absolute path of the author's
# working tree.  This is the ONLY edit made to this file in the reproduction package;
# it points the same variable at the package root, which carries the same directory
# names.  Set PAPERA_PKG to relocate the package.  See
# environment/PATCHES_TO_PACKAGE_CODE.md.
# -------------------------------------------------------------------------------------
import os as _os
_PKG_ROOT = _os.environ.get('PAPERA_PKG') or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..'))

import json
import os
from fractions import Fraction as F

CAMP = _os.path.join(_PKG_ROOT, 'multi_geometry_campaign')
MINE = os.path.join(CAMP, "data", "gate_su5_bound", "d8_classwide.json")
PHASE = os.path.join(CAMP, "data", "su5_ten_bound", "classwide.json")

# the table AS PRINTED in results/su5_ten_bound_findings.md 7.1
PRINTED = {
    "Z2xZ4_1_1": (7, "13/16", "11/16", []),
    "Z4xZ4_1_1": (15, "13/16", "11/16", []),
    "Z3xZ3_1_1": (8, "7/9", "2/3", []),
    "Z6-I_1_1": (5, "3/4", "2/3", []),
    "Z6-II_1_1": (5, "7/9", "25/36", []),
    "Z4_1_1": (3, "3/4", "11/16", []),
    "Z3xZ6_1_1": (17, "31/36", "2/3", [[0, 1], [0, 5]]),
}


def main():
    mine = json.load(open(MINE, encoding="utf-8"))
    phase = json.load(open(PHASE, encoding="utf-8"))
    verd = {}
    for lab, (n, bmax, bmin, loop) in PRINTED.items():
        r = mine[lab]
        agree = (r["n_twisted"] == n and r["B_max"] == bmax and r["B_min"] == bmin
                 and sorted(r["ten_doubling_permitted"]) == sorted(loop))
        pj = phase[lab]
        agree_json = (pj["n_twisted"] == r["n_twisted"] and pj["B_max"] == r["B_max"]
                      and pj["B_min"] == r["B_min"]
                      and sorted(pj["ten_doubling_permitted"])
                      == sorted("(%d, %d)" % tuple(x) for x in r["ten_doubling_permitted"]))
        verd[lab] = dict(printed_md=dict(n_twisted=n, B_max=bmax, B_min=bmin,
                                         ten_doubling=loop),
                         gate=dict(n_twisted=r["n_twisted"], B_max=r["B_max"],
                                   B_min=r["B_min"],
                                   ten_doubling=r["ten_doubling_permitted"],
                                   B_multiset=r["B_multiset"]),
                         verdict_vs_markdown="AGREE" if agree else "DISAGREE",
                         verdict_vs_classwide_json="AGREE" if agree_json else "DISAGREE")
    mine["_gate_verdict_vs_phase_7_1"] = verd
    mine["_gate_notes"] = {
        "twist_provenance": "all seven read off PRIMARY orbifolder geometry files; "
                            "second route via campaign code/geomlib.py agrees 0 mismatches",
        "measured_cross_check": "per-sector MAX measured budget in "
                                "data/su5_ten_bound/verify_{Z2xZ4_1_1,Z3xZ6_1_1}.json equals "
                                "the recomputed N_osc=0 budget in every populated sector "
                                "(6/6 and 13/13)",
        "B_min_is_N_osc_0": "the table's B_min is min over SECTORS of 1-delta_c, NOT the "
                            "min over STATES of B; 4,682 role-carrying 5bar fields in the "
                            "phase's own data have B < 13/20 (Z2xZ4 4,038 at B in {9/16, "
                            "7/16}; Z3xZ6 644 at B in {4/9, 19/36, 7/12}) and are therefore "
                            "NOT permitted to double by the theorem",
        "unpopulated_sectors": {"Z2xZ4_1_1": [[1, 1]],
                                "Z3xZ6_1_1": [[1, 1], [1, 2], [1, 3], [2, 1]]},
        "susy_condition": "sum_a eta_a integral in 60/60 twisted sectors across all seven",
    }
    with open(MINE, "w", encoding="utf-8") as fh:
        json.dump(mine, fh, indent=1)
    for lab, v in verd.items():
        print("%-12s md:%-9s json:%-9s  gate=%s" %
              (lab, v["verdict_vs_markdown"], v["verdict_vs_classwide_json"],
               v["gate"]["B_multiset"]))


if __name__ == "__main__":
    main()

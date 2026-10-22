"""GATE d8c -- second, independent route to the same table.

Route A (d8_classwide.py): own regex parser of the primary geometry files.
Route B (here): the campaign's own primary-file parser, code/geomlib.py.

Any disagreement in the fractional twists or in delta_c is a parser bug on one side.
Also cross-checks the recomputed N_osc=0 budgets against the B values actually MEASURED
per sector in data/su5_ten_bound/verify_*.json (max over states = the N_osc=0 value),
and reports the sub-threshold budgets that oscillator-carrying states actually reach.
"""

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
import sys
from fractions import Fraction as F

CAMP = _os.path.join(_PKG_ROOT, 'multi_geometry_campaign')
sys.path.insert(0, os.path.join(CAMP, "code"))
import geomlib                                                        # noqa: E402

LABELS = ["Z2xZ4_1_1", "Z4xZ4_1_1", "Z3xZ3_1_1", "Z6-I_1_1",
          "Z6-II_1_1", "Z4_1_1", "Z3xZ6_1_1"]


def frac(x):
    return "%d/%d" % (x.numerator, x.denominator) if x.denominator != 1 else str(x.numerator)


def main():
    mine = json.load(open(os.path.join(CAMP, "data", "gate_su5_bound",
                                       "d8_classwide.json"), encoding="utf-8"))
    bad = 0
    for lab in LABELS:
        g = geomlib.Geometry(lab)
        tw_b = [[frac(x) for x in t] for t in g.twists]
        tw_a = mine[lab]["twist_vectors"]
        ok = tw_a == tw_b
        print("%-12s twists routeA==routeB: %s   %s" % (lab, ok, tw_b))
        if not ok:
            bad += 1
            print("   A:", tw_a)
        for s in mine[lab]["sectors"]:
            fr = g.frac_twist(s["k"], s["l"])
            dc = sum(e * (1 - e) for e in fr) / 2
            if [frac(e) for e in fr] != s["eta"] or frac(dc) != s["delta_c"]:
                bad += 1
                print("   MISMATCH (%d,%d): A eta=%s dc=%s | B eta=%s dc=%s"
                      % (s["k"], s["l"], s["eta"], s["delta_c"],
                         [frac(e) for e in fr], frac(dc)))
    print("\nroute-A/route-B mismatches:", bad)

    # ---- measured vs arithmetic, on the two geometries that carry data
    print("\nMEASURED per-sector budgets (verify_*.json) vs recomputed N_osc=0 budget:")
    for lab in ["Z2xZ4_1_1", "Z3xZ6_1_1"]:
        p = os.path.join(CAMP, "data", "su5_ten_bound", "verify_%s.json" % lab)
        d = json.load(open(p, encoding="utf-8"))
        byk = {tuple(int(x) for x in s.strip("()").split(",")): v
               for s, v in d["B_by_sector"].items()}
        rec = {(s["k"], s["l"]): s["B"] for s in mine[lab]["sectors"]}
        print(" ", lab)
        sub = 0
        for kl, cnt in sorted(byk.items()):
            Bs = sorted((F(x) for x in cnt), reverse=True)
            hit = "OK" if frac(Bs[0]) == rec[kl] else "!! MAX MISMATCH (rec %s)" % rec[kl]
            low = [frac(b) for b in Bs if b < F(13, 20)]
            sub += sum(cnt[frac(b)] for b in Bs if b < F(13, 20))
            print("    T%s measured B=%s  max=%s %s  below-13/20: %s"
                  % (kl, [frac(b) for b in Bs], frac(Bs[0]), hit, low or "-"))
        unpop = [kl for kl in rec if kl not in byk]
        print("    sectors with NO role-carrying field in the data:", sorted(unpop))
        print("    fields whose ACTUAL budget B < 13/20 (5bar-doubling NOT permitted "
              "for them):", sub)


if __name__ == "__main__":
    main()

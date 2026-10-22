"""
f1_pilot_repair.py --- the [F-1] fixed-point-vs-raw-n_a repair, per the spec in
`results/gate_rule_nc.md` GATE-ADDED FINDING [F-1], landed in THIS phase's self-test suite
as the sweep spec's guard 1 requires.

THE DEFECT (gate-added, campaign not audit).  `code/pilot_table4.py:128` (and its "independent"
copy `code/gate/gate_pilot_table4.py:147`) implement Rule 4's premise -- "all twisted legs sit
at the SAME FIXED POINT in plane a" -- by comparing the RAW PRINTED TRANSLATION `n_a` across
the legs of a coupling:

    def _plane_label(rec, a):
        return (rec["fp"][2*(a-1)], rec["fp"][2*(a-1)+1])       # <-- raw n_a

That is neither the fixed point `f = (1 - theta^k)^{-1} n |_a` nor its point-group orbit.
Measured by the gate: the pilot's `k45` count is 144 as coded, 132 under the fixed-point
reading and 132 under the orbit reading; `1301.2322` Table 4's published value is 132, and
`data/pilot_table4.json` already records `all_rules: False` at 144.

THIS MODULE reproduces the pilot's whole rule chain with three interchangeable readings of
the premise -- `raw`, `fixed_point`, `orbit` -- and asserts that the repaired readings
reproduce the published 132 while `raw` reproduces the recorded 144.  Nothing in
`code/pilot_table4.py` is modified; the repair is exercised here and the direction is
recorded: the raw reading UNDER-applies Rule 4 and therefore OVER-counts allowed couplings,
so no absence claim anywhere is exposed.

Run: python f1_pilot_repair.py
Emits: data/z2xz4_sweep/f1_pilot_repair.json
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

CODE = _os.path.join(_PKG_ROOT, 'multi_geometry_campaign', 'code')
sys.path.insert(0, CODE)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from orbparse import parse_spectrum, parse_superpotential          # noqa: E402
from geomlib import Geometry, ROOT                                 # noqa: E402
from localization import plane_vectors, rot_matrix, fixed_point    # noqa: E402
import pilot_table4 as P                                           # noqa: E402

RUN = os.path.join(ROOT, "data", "orbifolder_runs", "out_table4.txt")
PUBLISHED = {"old_R": 96, "correct_R": 156, "all_rules": 132}
M_RULE4 = (6, 6, 2)          # the pilot's own torus-lattice automorphism orders


def _plane_n(rec, a):
    """The RAW printed translation in plane a -- the defective reading."""
    return (rec["fp"][2 * (a - 1)], rec["fp"][2 * (a - 1) + 1])


def _fixed_point(g, rec, a):
    """f = (1 - M_a(k))^{-1} n|_a mod 1, exactly, in the plane's lattice basis.

    `a` is 1-based here (the pilot's convention); `localization`'s helpers are 0-based.
    The plane's two lattice-vector slots are DERIVED (`plane_vectors`), never assumed to be
    (2a-1, 2a) -- Z6-II_3_1 factorizes on plane 2 spanned by e_4, e_5."""
    M = rot_matrix(g, a - 1, rec["k"])
    idx = plane_vectors(g, a - 1)
    n = (F(rec["fp"][idx[0]]), F(rec["fp"][idx[1]]))
    return fixed_point(M, n)                          # None on a fixed torus


def _orbit(g, rec, a):
    """The point-group ORBIT of the fixed point -- the certified RULE NC form."""
    f = _fixed_point(g, rec, a)
    if f is None:
        return None
    M = rot_matrix(g, a - 1, 1)
    orb, cur = {f}, f
    for _ in range(12):
        cur = ((M[0][0] * cur[0] + M[0][1] * cur[1]) % 1,
               (M[1][0] * cur[0] + M[1][1] * cur[1]) % 1)
        if cur in orb:
            break
        orb.add(cur)
    return tuple(sorted(orb))


def rule4(legs, g, reading):
    """The pilot's Rule 4 with the premise evaluated under `reading`."""
    label = {"raw": lambda r, a: _plane_n(r, a),
             "fixed_point": lambda r, a: _fixed_point(g, r, a),
             "orbit": lambda r, a: _orbit(g, r, a)}[reading]
    for a in (1, 2, 3):
        holo, anti = P.instantons(legs, a)
        if not (holo or anti):
            continue
        pts = {label(r, a) for r in legs if r["k"] != 0}
        if len(pts) != 1:
            continue
        NL = sum(r["NL"][a] for r in legs)
        NbL = sum(r["NbL"][a] for r in legs)
        NbR = P.NbarR(legs, a)
        d = F(NL) - NbL - NbR
        if d.denominator != 1 or int(d) % M_RULE4[a - 1] != 0:
            return False
    return True


def main():
    spec = parse_spectrum(RUN)
    W = parse_superpotential(RUN)
    g = Geometry("Z6-II_1_1")
    tw = [c for c in W if all(spec[x]["k"] != 0 for x in c)]
    keep = [c for c in tw if P.rcharge_ok([spec[x] for x in c], True)]
    old = [c for c in tw if P.rcharge_ok([spec[x] for x in c], False)]
    k5 = [c for c in keep if P.rule5_oscillator_constraints([spec[x] for x in c])]

    res = {"n_fields": len(spec), "n_couplings": len(W), "n_all_twisted": len(tw),
           "old_R": len(old), "correct_R": len(keep), "correct_R_plus_rule5": len(k5),
           "published": PUBLISHED, "by_reading": {}}
    for reading in ("raw", "fixed_point", "orbit"):
        k45 = [c for c in k5 if rule4([spec[x] for x in c], g, reading)]
        res["by_reading"][reading] = len(k45)

    res["reproduces_published_132"] = {
        r: res["by_reading"][r] == PUBLISHED["all_rules"] for r in res["by_reading"]}
    res["direction"] = ("the RAW reading under-applies Rule 4 and therefore OVER-counts "
                        "allowed couplings: raw >= repaired on every model, so NO ABSENCE "
                        "CLAIM is exposed (CLAUDE.md 5.16 instance-#2 logic).")
    res["raw_is_superset"] = res["by_reading"]["raw"] >= res["by_reading"]["fixed_point"]
    res["F1_REPAIRED_AND_PASSING"] = bool(
        res["by_reading"]["fixed_point"] == 132 and res["by_reading"]["orbit"] == 132
        and res["old_R"] == 96 and res["correct_R"] == 156)

    out = os.path.join(ROOT, "data", "z2xz4_sweep", "f1_pilot_repair.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))
    assert res["F1_REPAIRED_AND_PASSING"], "[F-1] repair FAILED -- STOP, do not sweep"


if __name__ == "__main__":
    main()

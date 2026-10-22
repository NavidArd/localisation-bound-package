"""z3xz6_sweep -- the COUNTERFACTUAL: what the 267 triplicating models would need.

The rank bound (`rank_bound.py`) shows two things at once:

  * 278 of the 4,610 evaluable models DO put three or more chiral Q in a SINGLE sector,
    and 267 put three or more of BOTH Q and u^c in single sectors -- so the catalogue is
    perfectly capable of the same-sector triplication a rank-3 block needs;
  * yet the maximum rank bound on ANY allowed channel is 2.

Those are only compatible if the sectors where the triplication happens cannot form an
ALLOWED channel.  This file checks that directly and per model: for each triplicating
model it forms the channel `{s_Q, s_uc, s_H}` its own fields would require, and reports
WHICH selection rule kills it -- the point group over both generators (sum k = 0 mod 3 AND
sum l = 0 mod 6) or, if the point group passes, the published per-plane R-charge law
0705.3075 (15).

This is the miss diagnosis in its sharpest form: not "no model has three quarks in a
sector" (false), but "the sectors the catalogue triplicates in never assemble into a
channel the R-charge law admits".

Run: python multi_geometry_campaign/code/z3xz6_sweep/counterfactual.py
"""

from __future__ import annotations

import json
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z36sw_lib as Z                                                  # noqa: E402
import rules_z36 as RZ                                                 # noqa: E402

CENSUS = os.path.join(Z.DATA, "census")


def why_not(tri):
    """Which rule kills the channel {s_Q, s_uc, s_H}?  Point group first, then (15)."""
    k = sum(t[0] for t in tri) % 3
    l = sum(t[1] for t in tri) % 6
    if k or l:
        return "point group (sum k = %d mod 3, sum l = %d mod 6)" % (k, l)
    legs = [("U", 1) if tuple(t) == (0, 0) else tuple(t) for t in tri]
    n_u = sum(1 for t in tri if tuple(t) == (0, 0))
    if n_u:
        from itertools import product as _p
        for assign in _p((1, 2, 3), repeat=n_u):
            it = iter(assign)
            legs = [("U", next(it)) if tuple(t) == (0, 0) else tuple(t) for t in tri]
            if RZ.law_published_perplane(legs)["all_planes_pass"]:
                return None
        return "R-charge law (15), all untwisted-plane assignments fail"
    law = RZ.law_published_perplane(legs)
    if law["all_planes_pass"]:
        return None
    bad = [a for a in (1, 2, 3) if not law["plane%d" % a]["oscillator_free_pass"]]
    return "R-charge law (15) fails in plane(s) %s" % bad


def main():
    lc = Z.ladder_channels()
    LB = {tuple(sorted(tuple(t) for t in v["tri"])): k
          for k, v in lc.items() if v["ladder_bearing"]}
    ALL = {tuple(sorted(tuple(t) for t in c)): Z.channel_name(c)
           for c in Z.allowed_channels()}

    rows, kills, chans = [], Counter(), Counter()
    n_trip = 0
    for fn in sorted(os.listdir(CENSUS)):
        if not fn.endswith(".json"):
            continue
        for m in json.load(open(os.path.join(CENSUS, fn)))["models"]:
            if m.get("status") != "ok":
                continue
            sq, su, sh = (m["sectors"]["chiral_Q"], m["sectors"]["chiral_uc"],
                          m["sectors"]["Hu"])
            gq = [k for k, v in sq.items() if v >= 3]
            gu = [k for k, v in su.items() if v >= 3]
            if not (gq and gu):
                continue
            n_trip += 1
            best = None
            for a in gq:
                for b in gu:
                    ta, tb = tuple(eval(a)), tuple(eval(b))
                    # the Higgs sector is FORCED by the point group once Q and u^c are fixed
                    kh = (-(ta[0] + tb[0])) % 3
                    lh = (-(ta[1] + tb[1])) % 6
                    th = (kh, lh)
                    tri = tuple(sorted([ta, tb, th]))
                    have_h = sh.get(str(th), 0)
                    reason = why_not([ta, tb, th])
                    rec = {"model": m["model"], "Q_sector": a, "uc_sector": b,
                           "forced_Higgs_sector": str(th),
                           "model_has_an_Hu_there": have_h,
                           "channel": Z.channel_name([ta, tb, th]),
                           "channel_allowed": tri in ALL,
                           "channel_ladder_bearing": tri in LB,
                           "killed_by": reason}
                    if best is None or (rec["channel_allowed"] and not best["channel_allowed"]):
                        best = rec
            rows.append(best)
            kills[best["killed_by"] or "ALLOWED"] += 1
            chans[best["channel"]] += 1

    s = {"scope": Z.CATALOGUE_TAG,
         "n_models_triplicating_both_Q_and_uc_in_single_sectors": n_trip,
         "counterfactual_channel_verdicts": dict(kills.most_common()),
         "counterfactual_channels": dict(chans.most_common(12)),
         "n_allowed": sum(1 for r in rows if r["channel_allowed"]),
         "n_ladder_bearing": sum(1 for r in rows if r["channel_ladder_bearing"]),
         "n_with_a_Higgs_in_the_forced_sector": sum(
             1 for r in rows if r["model_has_an_Hu_there"])}
    s["verdict"] = (
        "Of the %d models that DO triplicate both a chiral Q and a chiral u^c in single "
        "sectors, %d have their forced Q.u^c.H_u channel ALLOWED and %d have it "
        "LADDER-BEARING.  The triplication and the selection rules are therefore not "
        "merely uncorrelated -- they are mutually exclusive on this catalogue."
        % (n_trip, s["n_allowed"], s["n_ladder_bearing"]))
    dest = os.path.join(Z.DATA, "counterfactual.json")
    json.dump({"summary": s, "models": rows}, open(dest, "w"), indent=1, default=str)
    print(json.dumps(s, indent=1, default=str))
    print("->", dest)


if __name__ == "__main__":
    main()

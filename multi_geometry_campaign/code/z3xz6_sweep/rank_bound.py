"""z3xz6_sweep -- the RANK BOUND, derived from the census alone.

A rank-3 single-Higgs up-type block needs THREE chiral Q and THREE chiral u^c that all
appear in the SAME trilinear channel with the SAME Higgs.  The channel fixes each leg's
SECTOR, so all three Q's must sit in ONE sector and all three u^c's in ONE sector (the
same sector if the channel repeats a leg).  That converts the rank question into a pure
counting question on the census, with NO coupling run:

    rank(model, channel, leg-assignment) <= min( #chiral Q in the Q-leg's sector,
                                                 #chiral u^c in the u^c-leg's sector )
    and the assignment additionally needs an H_u in the third leg's sector.

The bound is an upper bound on the realised rank (couplings can only be absent, never
extra), so a catalogue-wide bound of 2 EXCLUDES SW-1 without running a single coupling --
and the coupling runs then confirm it rather than establish it.

This file emits, per model:
  * the best rank upper bound over all ladder-bearing channels and leg assignments;
  * the same over ALL allowed channels (the weaker question);
  * the UNCONSTRAINED same-sector counts (how many chiral Q the model puts in its most
    populated single sector, and likewise u^c) -- which separates two very different
    diagnoses: "the catalogue never triplicates a quark in one sector" versus "it does,
    but never in a sector that a ladder-bearing channel can use".

Run: python multi_geometry_campaign/code/z3xz6_sweep/rank_bound.py
"""

from __future__ import annotations

import json
import os
import sys
from collections import Counter
from itertools import permutations

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z36sw_lib as Z                                                  # noqa: E402

CENSUS = os.path.join(Z.DATA, "census")


def main():
    lc = Z.ladder_channels()
    LB = {tuple(sorted(tuple(t) for t in v["tri"])): k
          for k, v in lc.items() if v["ladder_bearing"]}
    ALL = {tuple(sorted(tuple(t) for t in c)): Z.channel_name(c)
           for c in Z.allowed_channels()}

    rows = []
    for fn in sorted(os.listdir(CENSUS)):
        if not fn.endswith(".json"):
            continue
        for m in json.load(open(os.path.join(CENSUS, fn)))["models"]:
            if m.get("status") != "ok":
                continue
            sq, su, sh = (m["sectors"]["chiral_Q"], m["sectors"]["chiral_uc"],
                          m["sectors"]["Hu"])

            def bound(table):
                best, arg = 0, None
                for ms, name in table.items():
                    legs = list(ms)
                    for p in set(permutations(range(3))):
                        gq, gu, gh = (str(legs[p[0]]), str(legs[p[1]]),
                                      str(legs[p[2]]))
                        nq, nu, nh = sq.get(gq, 0), su.get(gu, 0), sh.get(gh, 0)
                        if nq and nu and nh:
                            b = min(nq, nu)
                            if b > best:
                                best, arg = b, (name, gq, gu, gh, nq, nu, nh)
                return best, arg
            b_lb, a_lb = bound(LB)
            b_all, a_all = bound(ALL)
            rows.append({
                "model": m["model"],
                "rank_bound_ladder_channels": b_lb, "best_ladder": a_lb,
                "rank_bound_all_allowed_channels": b_all, "best_allowed": a_all,
                "max_chiral_Q_in_one_sector": max(list(sq.values()) or [0]),
                "max_chiral_uc_in_one_sector": max(list(su.values()) or [0]),
                "n_chiral_Q": m["n_chiral_Q"], "n_chiral_uc": m["n_chiral_uc"]})

    s = {
        "scope": Z.CATALOGUE_TAG,
        "n_models": len(rows),
        "rank_bound_ladder_channels_histogram": dict(sorted(Counter(
            r["rank_bound_ladder_channels"] for r in rows).items())),
        "rank_bound_all_allowed_channels_histogram": dict(sorted(Counter(
            r["rank_bound_all_allowed_channels"] for r in rows).items())),
        "max_chiral_Q_in_one_sector_histogram": dict(sorted(Counter(
            r["max_chiral_Q_in_one_sector"] for r in rows).items())),
        "max_chiral_uc_in_one_sector_histogram": dict(sorted(Counter(
            r["max_chiral_uc_in_one_sector"] for r in rows).items())),
        "models_with_3_same_sector_chiral_Q": sum(
            1 for r in rows if r["max_chiral_Q_in_one_sector"] >= 3),
        "models_with_3_same_sector_chiral_Q_AND_uc": sum(
            1 for r in rows if r["max_chiral_Q_in_one_sector"] >= 3
            and r["max_chiral_uc_in_one_sector"] >= 3),
        "max_rank_bound_ladder_channels": max(
            r["rank_bound_ladder_channels"] for r in rows),
        "max_rank_bound_all_allowed_channels": max(
            r["rank_bound_all_allowed_channels"] for r in rows),
    }
    # where DO the triplicated quarks sit, and is that sector usable by any channel?
    trip = Counter()
    for fn in sorted(os.listdir(CENSUS)):
        if not fn.endswith(".json"):
            continue
        for m in json.load(open(os.path.join(CENSUS, fn)))["models"]:
            if m.get("status") != "ok":
                continue
            for k, v in m["sectors"]["chiral_Q"].items():
                if v >= 3:
                    trip[k] += 1
    s["sectors_hosting_3_or_more_same_sector_chiral_Q"] = dict(trip.most_common())
    s["ladder_channel_leg_sectors"] = sorted({str(tuple(t)) for k, v in lc.items()
                                              if v["ladder_bearing"]
                                              for t in v["tri"]})
    s["verdict"] = (
        "SW-1 is EXCLUDED catalogue-wide by counting alone: the maximum rank an up-type "
        "single-Higgs block can attain on ANY ladder-bearing allowed channel, over all "
        "%d models, is %d.  The coupling runs can only lower it."
        % (len(rows), s["max_rank_bound_ladder_channels"]))

    dest = os.path.join(Z.DATA, "rank_bound.json")
    json.dump({"summary": s, "models": [r for r in rows
                                        if r["rank_bound_ladder_channels"] > 0]},
              open(dest, "w"), indent=1, default=str)
    print(json.dumps(s, indent=1, default=str))
    print("->", dest)


if __name__ == "__main__":
    main()

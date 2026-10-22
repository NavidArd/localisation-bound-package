"""
alphabet_z36.py --- the suppression alphabet of `Z3xZ6_1_1`, recomputed from THIS phase's own
geometry engine and rule set, and compared entry by entry against the phase-2 survey payload.

Phase: multi_geometry_campaign / Z3xZ6 foundations (docs/z3xz6_foundations_task.md), item 5.

INDEPENDENCE, stated precisely (no overclaim).  This module shares NO code with
`alphabet_survey/code/phase2/{p2_geom,p2_norm}.py`.  It re-implements, from the source
equations, the three ingredients that produce a NORM:

  (a) Rule 5 capability, 1301.2322 (2.22)/(2.23) at L = 3   --- `rules_z36.capable`
  (b) the space-group rule for a trilinear,
          lambda_1 + R_1 lambda_2 + R_1 R_2 lambda_3 = 0 (mod Lambda),  lambda_a = (1 - R_a) P_a
      with leg 1 gauge-fixed to orbit representatives, leg 2 scanned over its orbit and a
      lattice box, and leg 3 SOLVED (never scanned) --- so no configuration is missed
  (c) the classical-action prefactor of 1301.2322 (2.30)/(2.31): with (n, m) the smallest pair
      such that n g_1 = m g_2 in Z_M x Z_N, and A = 1 - theta^{n g_1} in plane a,
          prefactor = rho(A u) / rho(u)      (the RATIO form --- the gate-certified fix; the
                                              value form is only valid on Grams where e_1, e_2,
                                              e_1+e_2 have equal norm)
      and  NORM = prefactor * rho_min,  required EQUAL across every admissible leg ordering.

Anti-claim: this reproduces the survey's GEOMETRY-level object.  It says nothing about models.

Outputs: data/z3xz6/alphabet_cross_check_z36.json
Run:     python multi_geometry_campaign/code/z3xz6/alphabet_z36.py
"""

from __future__ import annotations

import json
import os
import sys
from fractions import Fraction as F
from itertools import permutations, product
from typing import Dict, List, Optional, Tuple

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z36lib import Geom, mapply, mmul, msub_from_id, mdet, minv_rat, mapply_rat, sector_label, OUT_DIR, REPO  # noqa: E402
import rules_z36 as RZ  # noqa: E402

PHASE2 = os.path.join(REPO, "alphabet_survey", "data", "phase2", "Z3xZ6_1_1.json")
G = Geom()
ORDERS = list(permutations((0, 1, 2)))


def rho(a: int, u) -> F:
    g = G.gram[a - 1]
    return g[0][0] * u[0] * u[0] + 2 * g[0][1] * u[0] * u[1] + g[1][1] * u[1] * u[1]


def n_p(g1, g2) -> Optional[Tuple[int, int]]:
    """Smallest (n, m), n minimal then m, with n g_1 = m g_2 in Z_M x Z_N."""
    O = G.M * G.N
    for n in range(1, O + 1):
        for m in range(1, O + 1):
            if (n * g1[0] - m * g2[0]) % G.M == 0 and (n * g1[1] - m * g2[1]) % G.N == 0:
                return n, m
    return None


def prefactor(a: int, g1, g2):
    nm = n_p(g1, g2)
    if nm is None:
        return None, "no (n,m)"
    n = nm[0]
    ng = ((n * g1[0]) % G.M, (n * g1[1]) % G.N)
    R = G.rot(ng, a)
    if R is None:
        return None, "degenerate: 1 - theta^{n g} = 0"
    A = msub_from_id(R)
    if A == ((0, 0), (0, 0)):
        return None, "degenerate: 1 - theta^{n g} = 0"
    vals = set()
    for u in ((F(1), F(0)), (F(0), F(1)), (F(1), F(1))):
        vals.add(rho(a, mapply(A, u)) / rho(a, u))
    if len(vals) != 1:
        return None, "ANISOTROPIC prefactor ratio %s" % sorted(map(str, vals))
    return vals.pop(), "ok"


def orbits_of(g, a):
    obs = G.orbits_in_plane(g, a)
    return obs


def entry_norm(a: int, gs, orb_choice, box: int = 6):
    """NORM for one plane and one per-leg orbit-class assignment; None if not well defined."""
    per_ordering = {}
    for order in ORDERS:
        gg = [gs[i] for i in order]
        oo = [orb_choice[i] for i in order]
        pf, info = prefactor(a, gg[0], gg[1])
        if pf is None:
            per_ordering["".join(map(str, order))] = (None, info)
            continue
        Rs = [G.rot(x, a) for x in gg]
        if any(R is None for R in Rs):
            per_ordering["".join(map(str, order))] = (None, "fixed torus leg")
            continue
        A = [msub_from_id(R) for R in Rs]
        R1R2 = mmul(Rs[0], Rs[1])
        A3inv = minv_rat(A[2])
        Binv = minv_rat(mmul(Rs[0], Rs[1]))
        rhos = set()
        for P1 in oo[0]:
            for P2 in oo[1]:
                for n1 in range(-box, box + 1):
                    for n2 in range(-box, box + 1):
                        P2s = (P2[0] + n1, P2[1] + n2)
                        l1 = mapply(A[0], P1)
                        l2 = mapply(A[1], P2s)
                        rhs = (-(l1[0] + mapply(Rs[0], l2)[0]),
                               -(l1[1] + mapply(Rs[0], l2)[1]))
                        w = mapply_rat(Binv, rhs)          # = A3 P3
                        P3 = mapply_rat(A3inv, w)
                        P3m = (P3[0] % 1, P3[1] % 1)
                        if P3m not in oo[2]:
                            continue
                        rhos.add(rho(a, (P2s[0] - P1[0], P2s[1] - P1[1])))
        if not rhos:
            per_ordering["".join(map(str, order))] = (None, "space-group forbidden")
            continue
        per_ordering["".join(map(str, order))] = (pf * min(rhos), "ok")
    vals = {v for v, _ in per_ordering.values() if v is not None}
    if not vals:
        return None, "space-group forbidden / all orderings degenerate", per_ordering
    if len(vals) != 1:
        return None, "REJECTED non-invariant", per_ordering
    return vals.pop(), "ok", per_ordering


def channel_alphabet(gs) -> Dict[str, object]:
    out = {"gs": [list(g) for g in gs], "planes": {}}
    for a in (1, 2, 3):
        cap = RZ.capable(gs, a)
        if not cap["capable"] or not cap["all_legs_twisted_in_plane"]:
            out["planes"]["plane%d" % a] = {"status": "not capable" if not cap["capable"]
                                            else "capable but a leg is DELOCALIZED"}
            continue
        orbs = [orbits_of(g, a) for g in gs]
        entries, norms = [], set()
        for choice in product(*[range(len(o)) for o in orbs]):
            oc = [orbs[i][choice[i]] for i in range(3)]
            val, info, per = entry_norm(a, gs, oc)
            entries.append({"orbit_class_index": list(choice),
                            "orbit_sizes": [len(o) for o in oc],
                            "NORM": None if val is None else str(val), "outcome": info})
            if val is not None:
                norms.add(val)
        out["planes"]["plane%d" % a] = {"status": "computable", "entries": entries,
                                        "NORMs": sorted(str(x) for x in norms)}
    return out


def main() -> None:
    with open(PHASE2, "r", encoding="utf-8") as fh:
        p2 = json.load(fh)
    p2ch = {c["channel"]["multiset"]: c["channel"] for c in p2["channels"]}

    # the ladder-bearing maximal-block channel and its CPT partner, plus the two other
    # all-planes-capable channels found by the rule set --- the load-bearing subset
    targets = ["T(0,1)T(0,2)T(0,3)", "T(0,3)T(0,4)T(0,5)",
               "T(2,3)T(2,4)T(2,5)", "T(2,4)T(2,4)T(2,4)",
               "T(0,1)T(0,1)T(0,4)", "T(0,2)T(0,2)T(0,2)"]
    rep = {"geometry": "Z3xZ6_1_1", "scope": "geometry statement, not a model statement",
           "channels_compared": [], "summary": {}}
    n_plane_ok = n_plane = n_entry_ok = n_entry = 0
    for name in targets:
        if name not in p2ch:
            rep["channels_compared"].append({"channel": name, "status": "not in phase-2 payload"})
            continue
        gs = [tuple(x) for x in p2ch[name]["gs"]]
        mine = channel_alphabet(gs)
        theirs = p2ch[name]["planes"]
        row = {"channel": name, "planes": {}}
        for a in (1, 2, 3):
            key = "plane%d" % a
            m, t = mine["planes"][key], theirs[key]
            n_plane += 1
            if m["status"] != "computable" or t.get("status") != "computable":
                agree = (m["status"] == t.get("status"))
                row["planes"][key] = {"mine": m["status"], "theirs": t.get("status"),
                                      "agree": agree}
                n_plane_ok += 1 if agree else 0
                continue
            mine_by = {tuple(e["orbit_class_index"]): e["NORM"] for e in m["entries"]}
            th_by = {tuple(e["orbit_class_index"]): e["NORM"] for e in t["assignments"]}
            same_entries = mine_by == th_by
            same_alph = m["NORMs"] == sorted(t["NORMs"])
            n_entry += len(th_by)
            n_entry_ok += sum(1 for k in th_by if mine_by.get(k) == th_by[k])
            row["planes"][key] = {"mine_NORMs": m["NORMs"], "theirs_NORMs": sorted(t["NORMs"]),
                                  "entries_agree": same_entries, "alphabet_agrees": same_alph,
                                  "mine_entries": {str(list(k)): v for k, v in mine_by.items()},
                                  "theirs_entries": {str(list(k)): v for k, v in th_by.items()}}
            n_plane_ok += 1 if (same_entries and same_alph) else 0
        rep["channels_compared"].append(row)

    rep["summary"] = {"planes_compared": n_plane, "planes_agreeing": n_plane_ok,
                      "entries_compared": n_entry, "entries_agreeing": n_entry_ok,
                      "union_alphabet_phase2": p2["union_alphabet"],
                      "n_channels_mine": len(RZ.all_multisets()),
                      "n_channels_phase2": p2["n_channels"],
                      "n_capable_channels_mine": sum(
                          1 for tri in RZ.all_multisets()
                          if any(RZ.capable(tri, a)["capable"] for a in (1, 2, 3))),
                      "n_capable_channels_phase2": p2["n_channels_with_capacity"]}

    os.makedirs(OUT_DIR, exist_ok=True)
    with open(os.path.join(OUT_DIR, "alphabet_cross_check_z36.json"), "w", encoding="utf-8") as fh:
        json.dump(rep, fh, indent=1, default=str)

    s = rep["summary"]
    print("=" * 78)
    print("Z3xZ6_1_1 alphabet --- independent recomputation vs the phase-2 payload")
    print("=" * 78)
    print("channels   mine %s   phase-2 %s" % (s["n_channels_mine"], s["n_channels_phase2"]))
    print("capable    mine %s   phase-2 %s" % (s["n_capable_channels_mine"], s["n_capable_channels_phase2"]))
    print("planes agreeing  %d / %d" % (s["planes_agreeing"], s["planes_compared"]))
    print("entries agreeing %d / %d" % (s["entries_agreeing"], s["entries_compared"]))
    for row in rep["channels_compared"]:
        print("-", row["channel"])
        for k, v in row.get("planes", {}).items():
            if "mine_NORMs" in v:
                print("    %s mine %s theirs %s  entries_agree=%s"
                      % (k, v["mine_NORMs"], v["theirs_NORMs"], v["entries_agree"]))
            else:
                print("    %s %s / %s agree=%s" % (k, v["mine"], v["theirs"], v["agree"]))


if __name__ == "__main__":
    main()

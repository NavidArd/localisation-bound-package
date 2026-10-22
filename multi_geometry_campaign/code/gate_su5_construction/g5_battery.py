"""V5 + V1(E) -- the gate's own recount of the SM battery over all 2,912 CB models.

Three things are recomputed here by fresh code, and the phase's own numbers are
reproduced alongside so that any difference is localised:

  (1) Tr Y and Tr Y^3 over left-chiral STATES (weight = full printed dimension), and
      -- INDEPENDENTLY -- from orbifolder's OWN printed anomaly traces tr Q_i, tr Q_i^3,
      tr Q_i Q_j^2 after expressing t_Y in the printed (orthogonal) U(1) basis.
      This is the 5.21 "every defining condition" check done twice.
  (2) The generation count with TRIALITY-AWARE roles (gs_lib G5) beside the phase's
      triality-BLIND roles (su3_dim's sign is discarded in s8_battery.role_of).
      Convention-blindness is implemented as the ONLY legitimate freedom: an overall
      conjugation p -> -p, which flips BOTH the SU(3) triality AND Y.  Flipping Y alone
      (the phase's sY branch on a triality-blind table) is not a relabelling of the same
      spectrum.
  (3) The doubled-10-role test with triality-aware roles.

Also: t_Y . t_anom (Green-Schwarz orthogonality) per model from orbifolder's own U(1)
generators and traces -- the phase verified this on the flagship only (its scope limit 4).

Output: data/gate_su5_construction/g5_battery.json
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

RAW = os.path.join(G.PHASE_DATA, "raw")
SECT = {(0, 1), (0, 5)}
DELTA = [F(0)] * 8 + [F(1), F(1)] + [F(0)] * 6
TEN = ("Q", "uc", "ec")


def blind_role(t, d2, y):
    """the phase's table: SU(3) triality discarded (d3 = 3 for both 3 and 3bar)."""
    d3 = 1 if t == 0 else 3
    if d3 == 3 and d2 == 2 and y == F(1, 6):
        return "Q"
    if d3 == 3 and d2 == 2 and y == F(-1, 6):
        return "Qbar"
    if d3 == 3 and d2 == 1 and y == F(-2, 3):
        return "uc"
    if d3 == 3 and d2 == 1 and y == F(2, 3):
        return "ucbar"
    if d3 == 1 and d2 == 1 and y == F(1):
        return "ec"
    if d3 == 1 and d2 == 1 and y == F(-1):
        return "ecbar"
    if d3 == 3 and d2 == 1 and y == F(1, 3):
        return "dc"
    if d3 == 1 and d2 == 2 and y == F(-1, 2):
        return "L"
    if d3 == 1 and d2 == 2 and y == F(1, 2):
        return "Hu"
    return None


def solve_span(basis, target):
    n = len(basis)
    if n == 0:
        return None
    A = [[G.dot(basis[i], basis[j]) for j in range(n)] for i in range(n)]
    b = [G.dot(basis[i], target) for i in range(n)]
    M = [A[i][:] + [b[i]] for i in range(n)]
    r = 0
    for c in range(n):
        piv = next((i for i in range(r, n) if M[i][c] != 0), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        pv = M[r][c]
        M[r] = [x / pv for x in M[r]]
        for i in range(n):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        r += 1
    coef = [M[i][n] for i in range(n)]
    rec = [F(0)] * 16
    for ci, bi in zip(coef, basis):
        rec = G.add(rec, G.smul(ci, bi))
    return coef if rec == list(target) else None


def do_segment(path):
    rows = []
    txt = gzip.open(path, "rt", errors="replace").read()
    for lab, blk in G.split_models(txt):
        pre = blk.split("Analyzing")[0]
        fields = G.parse_fields(pre)
        left = [f for f in fields if f["chir"] == "l"]
        if not left:
            continue
        u1 = G.parse_u1_generators(blk)
        tr = G.parse_anomaly_traces(blk)
        trY = trY3 = F(0)
        awc, blc = Counter(), Counter()
        hits_aware, hits_blind = [], []
        lvl = 0
        for f in left:
            p = f["p_sh"]
            try:
                t = G.su3_triality(p)
                d2 = G.su2_dim(p)
            except ValueError:
                lvl += 1
                continue
            d3 = 1 if t == 0 else 3
            tot = 1
            for d in f["rep"]:
                tot *= abs(d)
            y = G.dot(G.TY, p)
            trY += tot * y
            trY3 += tot * y ** 3
            hid = int(F(tot, d3 * d2))
            ra = G.role_of(p, 1)
            if ra:
                awc[ra] += hid
            rb = blind_role(t, d2, y)
            if rb:
                blc[rb] += hid
            if f["kl"] in SECT and G.n2(p) == F(31, 18) \
                    and abs(G.dot(p, DELTA)) == 1:
                if ra in ("Q", "uc", "ec", "Qbar", "ucbar", "ecbar"):
                    hits_aware.append((f["label"], ra, tuple(f["kl"]), hid))
                rb2 = rb or blind_role(t, d2, -y)
                if rb2 in ("Q", "uc", "ec", "Qbar", "ucbar", "ecbar"):
                    hits_blind.append((f["label"], rb2, tuple(f["kl"]), hid))
        net_aware = {k: awc.get(k, 0) - awc.get(k + "bar", 0) for k in TEN}
        # the phase's convention-blind rule, reproduced exactly
        netQ_blind = {}
        for sY in (1, -1):
            c = Counter()
            for f in left:
                p = f["p_sh"]
                try:
                    t = G.su3_triality(p)
                    d2 = G.su2_dim(p)
                except ValueError:
                    continue
                d3 = 1 if t == 0 else 3
                tot = 1
                for d in f["rep"]:
                    tot *= abs(d)
                hid = int(F(tot, d3 * d2))
                r = blind_role(t, d2, sY * G.dot(G.TY, p))
                if r == "Q":
                    c["Q"] += hid
                elif r == "Qbar":
                    c["Q"] -= hid
            netQ_blind[sY] = c["Q"]
        # GS orthogonality from orbifolder's own data
        coef = solve_span(u1, G.TY) if u1 else None
        trY_orb = None
        trY3_orb = None
        tYdotanom = None
        if coef is not None and tr:
            nz = [i for i, c in enumerate(coef) if c != 0]
            s = F(0)
            okall = True
            for i in nz:
                v = tr.get("Q_%d" % (i + 1))
                if v is None:
                    okall = False
                    break
                s += coef[i] * v
            trY_orb = str(s) if okall else None
            if okall and len(nz) == 2:
                i, j = nz
                try:
                    t3 = (coef[i] ** 3 * tr["Q_%d^3" % (i + 1)]
                          + coef[j] ** 3 * tr["Q_%d^3" % (j + 1)]
                          + 3 * coef[i] * coef[j] ** 2
                          * tr["Q_%d Q_%d^2" % (i + 1, j + 1)]
                          + 3 * coef[j] * coef[i] ** 2
                          * tr["Q_%d Q_%d^2" % (j + 1, i + 1)])
                    trY3_orb = str(t3)
                except KeyError:
                    pass
            an = [F(0)] * 16
            for i, g in enumerate(u1):
                v = tr.get("Q_%d" % (i + 1))
                if v:
                    an = G.add(an, G.smul(v / G.n2(g), g))
            tYdotanom = str(G.dot(G.TY, an))
        rows.append(dict(model=lab, trY=str(trY), trY3=str(trY3),
                         trY_orb=trY_orb, trY3_orb=trY3_orb,
                         tY_dot_tanom=tYdotanom,
                         net_aware=net_aware, netQ_blind={str(k): v for k, v
                                                          in netQ_blind.items()},
                         hits_aware=len(hits_aware), hits_blind=len(hits_blind),
                         hit_roles_aware=[h[1] for h in hits_aware],
                         hit_roles_blind=[h[1] for h in hits_blind],
                         level1_skipped=lvl))
    return rows


def main():
    paths = sorted(glob.glob(os.path.join(RAW, "*.txt.gz")))
    with Pool(6) as pool:
        allrows = [r for chunk in pool.map(do_segment, paths) for r in chunk]
    print("models:", len(allrows))

    f = Counter()
    disagree_trY, disagree_gen = [], []
    aware_gen3, blind_gen3 = [], []
    aware_gen3_all, blind_gen3_all = [], []
    for r in allrows:
        trY0 = r["trY"] == "0"
        trY30 = r["trY3"] == "0"
        if r["trY_orb"] is not None and r["trY_orb"] != r["trY"]:
            disagree_trY.append(r["model"])
        # triality-aware, convention-blind = overall conjugation only
        na = r["net_aware"]
        awQ3 = abs(na["Q"]) == 3
        awall = any(all(na[k] == e * 3 for k in TEN) for e in (1, -1))
        # the phase's rule
        nb = r["netQ_blind"]
        blQ3 = (nb["1"] == 3 or nb["-1"] == 3)
        if awQ3 != blQ3:
            disagree_gen.append(r["model"])
        if awQ3:
            aware_gen3.append(r["model"])
        if blQ3:
            blind_gen3.append(r["model"])
        if awall:
            aware_gen3_all.append(r["model"])
        f[("trY0" if trY0 else "trYx",
           "trY30" if trY30 else "trY3x",
           "awQ3" if awQ3 else "-",
           "blQ3" if blQ3 else "-",
           "dblA" if r["hits_aware"] else "-",
           "dblB" if r["hits_blind"] else "-")] += 1

    summary = dict(
        models=len(allrows),
        trY0=sum(1 for r in allrows if r["trY"] == "0"),
        trY3_0=sum(1 for r in allrows if r["trY3"] == "0"),
        trY0_and_trY3_0=sum(1 for r in allrows
                            if r["trY"] == "0" and r["trY3"] == "0"),
        trY0_but_trY3_nonzero=sum(1 for r in allrows
                                  if r["trY"] == "0" and r["trY3"] != "0"),
        orbifolder_trY_route_available=sum(1 for r in allrows
                                           if r["trY_orb"] is not None),
        orbifolder_trY_disagreements=len(disagree_trY),
        orbifolder_trY3_route_available=sum(1 for r in allrows
                                            if r["trY3_orb"] is not None),
        orbifolder_trY3_disagreements=sum(
            1 for r in allrows
            if r["trY3_orb"] is not None and r["trY3_orb"] != r["trY3"]),
        gs_orth_checked=sum(1 for r in allrows if r["tY_dot_tanom"] is not None),
        gs_orth_nonzero=sum(1 for r in allrows
                            if r["tY_dot_tanom"] not in (None, "0")),
        netQ3_blind_PHASE=len(blind_gen3),
        netQ3_aware_GATE=len(aware_gen3),
        netQ3_disagreements=len(disagree_gen),
        gen3_all_three_aware=len(aware_gen3_all),
        doubled_blind_PHASE=sum(1 for r in allrows if r["hits_blind"]),
        doubled_aware_GATE=sum(1 for r in allrows if r["hits_aware"]),
        passing_phase_rule=sum(1 for r in allrows
                               if r["trY"] == "0"
                               and (r["netQ_blind"]["1"] == 3
                                    or r["netQ_blind"]["-1"] == 3)
                               and r["hits_blind"]),
        passing_gate_rule=sum(1 for r in allrows
                              if r["trY"] == "0" and abs(r["net_aware"]["Q"]) == 3
                              and r["hits_aware"]),
        level1_skipped_total=sum(r["level1_skipped"] for r in allrows),
        netQ_aware_hist=dict(Counter(r["net_aware"]["Q"] for r in allrows)),
        net_uc_aware_hist=dict(Counter(r["net_aware"]["uc"] for r in allrows)),
        net_ec_aware_hist=dict(Counter(r["net_aware"]["ec"] for r in allrows)),
        funnel={"|".join(k): v for k, v in sorted(f.items())},
    )
    print(json.dumps(summary, indent=1))
    json.dump(dict(summary=summary, rows=allrows,
                   aware_gen3=aware_gen3, blind_gen3=blind_gen3,
                   disagree_gen=disagree_gen, disagree_trY=disagree_trY),
              open(os.path.join(G.OUT, "g5_battery.json"), "w"), indent=1)
    print("wrote g5_battery.json")


if __name__ == "__main__":
    main()

"""V2 (final) -- the SM battery and the three-generation screen on the NEW V2 class.

Same battery as g5 (triality-aware roles, full-dimension trace weights, orbifolder's own
anomaly traces as the independent route), applied to the 1,432 case-B models.

THIS IS THE NUMBER THAT DECIDES WHETHER THE PHASE'S NEGATIVE IS CLASS-LOCAL.

Output: data/gate_su5_construction/g8_batteryB.json
"""
from __future__ import annotations

import gzip
import json
import os
import sys
from collections import Counter
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gs_lib as G
from g5_battery import solve_span

SECT = {(0, 1), (0, 5)}
TEN = ("Q", "uc", "ec")


def main():
    lift = json.load(open(os.path.join(G.OUT, "g4_lift_B.json")))
    p = [G.rat(t) for t in lift["V2"]]
    DELTA = [G.rat(t) for t in lift["Delta"]]
    print("case B  V2 =", ", ".join(str(t) for t in p))
    print("        Delta =", ", ".join(str(t) for t in DELTA), " |Delta|^2 =", G.n2(DELTA))
    assert G.n2(DELTA) == 2 and G.dot(p, DELTA) == 1
    assert all(G.dot(DELTA, c) == 0 for c in G.SM_CARTAN)

    txt = gzip.open(os.path.join(G.OUT, "gate_caseB.txt.gz"), "rt",
                    errors="replace").read()
    rows = []
    ggs = Counter()
    for lab, blk in G.split_models(txt, prefix="GB"):
        fields = G.parse_fields(blk.split("Analyzing")[0])
        left = [f for f in fields if f["chir"] == "l"]
        if not left:
            continue
        u1 = G.parse_u1_generators(blk)
        tr = G.parse_anomaly_traces(blk)
        simple = G.parse_simple_roots(blk)
        gg = G.parse_gauge_group(blk)
        ggs[gg.split(" and ")[0] if gg else "?"] += 1
        trY = trY3 = F(0)
        c = Counter()
        hits = []
        lvl = 0
        for f in left:
            q = f["p_sh"]
            try:
                t = G.su3_triality(q)
                d2 = G.su2_dim(q)
            except ValueError:
                lvl += 1
                continue
            d3 = 1 if t == 0 else 3
            tot = 1
            for d in f["rep"]:
                tot *= abs(d)
            y = G.dot(G.TY, q)
            trY += tot * y
            trY3 += tot * y ** 3
            hid = int(F(tot, d3 * d2))
            r = G.role_of(q, 1)
            if r:
                c[r] += hid
            if f["kl"] in SECT and G.n2(q) == F(31, 18) \
                    and r in G.TEN | G.TENBAR and abs(G.dot(q, DELTA)) == 1:
                pd = G.dot(q, DELTA)
                w = G.sub(q, G.smul(pd, DELTA))
                hits.append(dict(field=f["label"], kl=list(f["kl"]), fp=list(f["fp"]),
                                 role=r, rep=list(f["rep"]), osc=f["osc"], hid=hid,
                                 p2=str(G.n2(q)), w2=str(G.proj_wsm_norm2(q)),
                                 R=str(G.n2(q) - G.proj_wsm_norm2(q)), pdotD=str(pd),
                                 partner_norm2=str(G.n2(w)),
                                 partner_role=G.role_of(w, 1)))
        coef = solve_span(u1, G.TY) if u1 else None
        trY_orb = None
        tanom = None
        if coef is not None and tr:
            nz = [i for i, x in enumerate(coef) if x != 0]
            s = F(0)
            okall = True
            for i in nz:
                v = tr.get("Q_%d" % (i + 1))
                if v is None:
                    okall = False
                    break
                s += coef[i] * v
            trY_orb = str(s) if okall else None
            an = [F(0)] * 16
            for i, g in enumerate(u1):
                v = tr.get("Q_%d" % (i + 1))
                if v:
                    an = G.add(an, G.smul(v / G.n2(g), g))
            tanom = str(G.dot(G.TY, an))
        net = {k: c.get(k, 0) - c.get(k + "bar", 0) for k in TEN}
        rows.append(dict(model=lab, trY=str(trY), trY3=str(trY3), trY_orb=trY_orb,
                         tY_dot_tanom=tanom, net_aware=net, hits=len(hits),
                         hit_roles=[h["role"] for h in hits],
                         tY_orth=all(G.dot(G.TY, s) == 0 for s in simple),
                         level1_skipped=lvl,
                         dossier=hits[:2] if hits else []))
    print("models censused:", len(rows))
    trY0 = [r for r in rows if r["trY"] == "0"]
    dbl = [r for r in rows if r["hits"]]
    q3 = [r for r in rows if abs(r["net_aware"]["Q"]) == 3]
    allthree = [r for r in rows
                if any(all(r["net_aware"][k] == e * 3 for k in TEN) for e in (1, -1))]
    passing = [r for r in rows if r["trY"] == "0" and abs(r["net_aware"]["Q"]) == 3
               and r["hits"]]
    summary = dict(
        models=len(rows),
        first_E8_factor_histogram=dict(ggs),
        tY_orth_all=sum(1 for r in rows if r["tY_orth"]),
        level1_skipped_total=sum(r["level1_skipped"] for r in rows),
        trY0=len(trY0),
        trY0_and_trY3_0=sum(1 for r in rows if r["trY"] == "0" and r["trY3"] == "0"),
        orb_trY_disagreements=sum(1 for r in rows if r["trY_orb"] is not None
                                  and r["trY_orb"] != r["trY"]),
        gs_orth_nonzero=sum(1 for r in rows if r["tY_dot_tanom"] not in (None, "0")),
        doubled_10role=len(dbl),
        netQ3=len(q3),
        NET_3_IN_ALL_THREE=len(allthree),
        trY0_and_netQ3_and_doubled=len(passing),
        netQ_hist=dict(Counter(r["net_aware"]["Q"] for r in rows)),
        net_uc_hist=dict(Counter(r["net_aware"]["uc"] for r in rows)),
        net_ec_hist=dict(Counter(r["net_aware"]["ec"] for r in rows)),
        max_roles_at_pm3=dict(Counter(
            max(sum(1 for k in TEN if r["net_aware"][k] == e * 3) for e in (1, -1))
            for r in rows)),
        doubled_role_hist=dict(Counter(x for r in rows for x in r["hit_roles"])),
    )
    print(json.dumps(summary, indent=1))
    if dbl:
        print("\nexample case-B doubled-10 dossier:", json.dumps(dbl[0], indent=1)[:1400])
    json.dump(dict(summary=summary, rows=rows),
              open(os.path.join(G.OUT, "g8_batteryB.json"), "w"), indent=1)
    print("wrote g8_batteryB.json")


if __name__ == "__main__":
    main()

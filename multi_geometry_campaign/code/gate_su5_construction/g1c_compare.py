"""V1 -- compare the GATE's fresh orbifolder run against the phase's stored output.

Record-by-record on the printed spectrum (label, sector, fixed point, space-group charges,
representation, chirality, U(1) charges, p_sh, q_sh, oscillators) plus the printed gauge
group, U(1) generators, simple roots and anomaly traces.  Also dossiers the gate's
replacement flagship candidates with the same battery as g1.
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

RAW = os.path.join(G.PHASE_DATA, "raw")
DELTA = [F(0)] * 8 + [F(1), F(1)] + [F(0)] * 6


def stored_blocks(labels):
    want = set(labels)
    out = {}
    for fn in sorted(os.listdir(RAW)):
        t = gzip.open(os.path.join(RAW, fn), "rt", errors="replace").read()
        for lab, blk in G.split_models(t):
            if lab in want:
                out[lab] = blk
        if len(out) == len(want):
            break
    return out


def canon(fields):
    return [(f["label"], f["kl"], f["fp"], f["sgc"], f["rep"], f["chir"],
             tuple(str(x) for x in f["u1"]), tuple(str(x) for x in f["p_sh"]),
             tuple(str(x) for x in f["q_sh"]), f["osc"]) for f in fields]


def dossier(lab, blk):
    fields = G.parse_fields(blk.split("Analyzing")[0])
    left = [f for f in fields if f["chir"] == "l"]
    simple = G.parse_simple_roots(blk)
    res = dict(model=lab, gauge=G.parse_gauge_group(blk), n_records=len(fields))
    trY = trY3 = F(0)
    c = Counter()
    hits = []
    for f in left:
        p = f["p_sh"]
        t = G.su3_triality(p)
        d2 = G.su2_dim(p)
        d3 = 1 if t == 0 else 3
        tot = 1
        for d in f["rep"]:
            tot *= abs(d)
        y = G.dot(G.TY, p)
        trY += tot * y
        trY3 += tot * y ** 3
        hid = int(F(tot, d3 * d2))
        r = G.role_of(p, 1)
        if r:
            c[r] += hid
        if f["kl"] in ((0, 1), (0, 5)) and G.n2(p) == F(31, 18) \
                and r in G.TEN | G.TENBAR and abs(G.dot(p, DELTA)) == 1:
            pd = G.dot(p, DELTA)
            q = G.sub(p, G.smul(pd, DELTA))
            hits.append(dict(field=f["label"], kl=list(f["kl"]), fp=list(f["fp"]),
                             role=r, rep=list(f["rep"]), osc=f["osc"],
                             hid=hid, p2=str(G.n2(p)),
                             w2=str(G.proj_wsm_norm2(p)),
                             R=str(G.n2(p) - G.proj_wsm_norm2(p)),
                             pdotD=str(pd), q_norm2=str(G.n2(q)),
                             q_w2=str(G.proj_wsm_norm2(q)),
                             q_role=G.role_of(q, 1)))
    res["trY"] = str(trY)
    res["trY3"] = str(trY3)
    res["net_aware"] = {k: c.get(k, 0) - c.get(k + "bar", 0)
                        for k in ("Q", "uc", "ec", "dc", "L")}
    res["raw_counts"] = dict(c)
    res["hits"] = hits
    res["kY"] = str(2 * G.n2(G.TY))
    res["tY_orth_simple_roots"] = all(G.dot(G.TY, s) == 0 for s in simple)
    res["n_simple_roots"] = len(simple)
    return res


def main():
    fresh = gzip.open(os.path.join(G.OUT, "gate_v1.txt.gz"), "rt",
                      errors="replace").read()
    fb = dict(G.split_models(fresh))
    labels = list(fb)
    print("fresh run models:", labels)
    sb = stored_blocks(labels)
    out = dict(compare={}, dossiers={})
    for lab in labels:
        a = G.parse_fields(fb[lab].split("Analyzing")[0])
        cmp_ = dict(fresh_records=len(a))
        if lab in sb:
            b = G.parse_fields(sb[lab].split("Analyzing")[0])
            ca, cb = canon(a), canon(b)
            cmp_.update(stored_records=len(b), identical=(ca == cb),
                        first_diff=None)
            if ca != cb:
                for i, (x, y) in enumerate(zip(ca, cb)):
                    if x != y:
                        cmp_["first_diff"] = [i, str(x), str(y)]
                        break
            cmp_["gauge_fresh"] = G.parse_gauge_group(fb[lab])
            cmp_["gauge_stored"] = G.parse_gauge_group(sb[lab])
            cmp_["gauge_identical"] = cmp_["gauge_fresh"] == cmp_["gauge_stored"]
            cmp_["u1_identical"] = (G.parse_u1_generators(fb[lab])
                                    == G.parse_u1_generators(sb[lab]))
            cmp_["simple_identical"] = (G.parse_simple_roots(fb[lab])
                                        == G.parse_simple_roots(sb[lab]))
            cmp_["anomaly_identical"] = (G.parse_anomaly_traces(fb[lab])
                                         == G.parse_anomaly_traces(sb[lab]))
        else:
            cmp_["stored_records"] = None
        out["compare"][lab] = cmp_
        print(lab, cmp_)
        out["dossiers"][lab] = dossier(lab, fb[lab])
    for lab, d in out["dossiers"].items():
        print("\n--", lab, d["gauge"])
        print("   trY", d["trY"], "trY3", d["trY3"], "net", d["net_aware"])
        for h in d["hits"]:
            print("   hit", h)
    json.dump(out, open(os.path.join(G.OUT, "g1c_compare.json"), "w"), indent=1)
    print("\nwrote g1c_compare.json")


if __name__ == "__main__":
    main()

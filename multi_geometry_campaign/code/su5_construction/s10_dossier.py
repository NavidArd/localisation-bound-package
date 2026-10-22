"""S10 -- the two-10-role DOSSIER: the explicit constructed model, state by state.

The partner state is the WEYL REFLECTION  q = p - (p.Delta) Delta = sigma_Delta(p),
NOT p - Delta:  p.Delta = -1 occurs as often as +1 and the naive subtraction then leaves
the shell (|q|^2 = 103/18 instead of 31/18).  Caught by the |q|^2 invariant.
"""
import glob, gzip, json, os, sys
from collections import Counter
from fractions import Fraction as F
import sc_lib as S
sys.path.insert(0, os.path.join(S.CAMP, "code", "z3xz6_sweep"))
import z36sw_lib as SW
from s8_battery import P, DELTA, TY, A3, A2, split_models_cb, su3_dim, su2_dim, role_of, RAW, SECT
from s9_dossier import proj2


def main():
    build = {m["label"]: m for m in
             json.load(open(os.path.join(S.OUT, "s5_build.json")))["models"]}
    d8 = json.load(open(os.path.join(S.OUT, "s8_battery.json")))
    pm = {m["model"]: m for m in d8["per_model"]}
    want = set()
    for x in d8["doubled"]:
        v = x["netQ"].get("1", x["netQ"].get(1))
        if abs(int(v)) == 3:
            want.add(x["model"])
    print("candidates with |netQ| = 3 and a doubled 10-role:", len(want))

    best, alldoss = None, []
    for pth in sorted(glob.glob(os.path.join(RAW, "*.txt.gz"))):
        txt = gzip.open(pth, "rt", errors="replace").read()
        for lab, blk in split_models_cb(txt):
            if lab not in want:
                continue
            spec = SW.parse_spectrum(blk.split("Analyzing")[0])
            gauge = SW.parse_gauge(blk)
            left = [r for r in spec.values()
                    if r["chir"] == "l" and r["p_sh"] and r["rep"]]
            hits = []
            for r in left:
                if tuple(r["kl"]) not in SECT:
                    continue
                p = r["p_sh"]
                if S.dot(p, p) != F(31, 18):
                    continue
                d3, _ = su3_dim(p)
                d2 = su2_dim(p)
                y = S.dot(TY, p)
                ro = role_of(d3, d2, y, 1, 1) or role_of(d3, d2, y, -1, 1)
                if ro not in ("Q", "uc", "ec", "Qbar", "ucbar", "ecbar"):
                    continue
                pd = S.dot(p, DELTA)
                if abs(pd) != 1:
                    continue
                q = [a - pd * b for a, b in zip(p, DELTA)]          # sigma_Delta(p)
                assert S.dot(q, q) == S.dot(p, p)
                hits.append(dict(field=r["label"], kl=list(r["kl"]),
                                 fp=list(r["fp"]) if r["fp"] else None, role=ro,
                                 rep=list(r["rep"]), osc=str(r["osc_raw"]),
                                 qsh_right=[str(x) for x in r["qsh"]],
                                 p_sh=[str(x) for x in p],
                                 partner_p_sh=[str(x) for x in q],
                                 p2=str(S.dot(p, p)), q2=str(S.dot(q, q)),
                                 Delta_diff=[str(a - b) for a, b in zip(p, q)],
                                 p_dot_Delta=str(pd),
                                 w2_p=str(proj2(p)), w2_q=str(proj2(q)),
                                 R_p=str(S.dot(p, p) - proj2(p)),
                                 R_q=str(S.dot(q, q) - proj2(q)),
                                 Y_p=str(y), Y_q=str(S.dot(TY, q)),
                                 su3=(d3, su3_dim(q)[0]), su2=(d2, su2_dim(q))))
            if not hits:
                continue
            b = build[lab]
            V1 = [F(x) for x in b["V1"]]
            W3 = [F(x) for x in b["W3"]]
            Z = [F(0)] * 16
            ok, bad = S.mi_report(V1, P, [Z, Z, W3, W3, Z, Z])
            ent = dict(model=lab, V1=b["V1"], V2=[str(x) for x in P], W3=b["W3"],
                       ordV1=S.order_of(V1), ordV2=S.order_of(P), ordW3=S.order_of(W3),
                       MI_full_pass=ok, MI_failures=bad,
                       gauge_TestConfig=gauge.get("gauge_test"),
                       kY=str(2 * S.dot(TY, TY)), tY=[str(x) for x in TY],
                       Delta=[str(x) for x in DELTA],
                       Delta_norm2=str(S.dot(DELTA, DELTA)),
                       Delta_in_lattice=S.in_lattice(DELTA),
                       Delta_dot_V1=str(S.dot(DELTA, V1)),
                       Delta_dot_V2=str(S.dot(DELTA, P)),
                       Delta_dot_W3=str(S.dot(DELTA, W3)),
                       Delta_dot_SM=[str(S.dot(DELTA, v))
                                     for v in (A3[0], A3[1], A2[0], TY)],
                       netQ=pm[lab]["netQ"], n_hits=len(hits), hits=hits)
            alldoss.append(ent)
            if best is None or len(hits) > len(best["hits"]):
                best = ent
    json.dump(dict(flagship=best, n_candidates=len(alldoss),
                   all=alldoss[:50]),
              open(os.path.join(S.OUT, "s10_dossier.json"), "w"), indent=1)
    print("dossiers built:", len(alldoss))
    print(json.dumps({k: v for k, v in best.items() if k != "hits"}, indent=1))
    for h in best["hits"]:
        print("  ", h["field"], h["kl"], h["role"], "rep", h["rep"],
              "|p|^2", h["p2"], "|w|^2", h["w2_p"], "R", h["R_p"],
              "p.Delta", h["p_dot_Delta"], "osc", h["osc"])


if __name__ == "__main__":
    main()

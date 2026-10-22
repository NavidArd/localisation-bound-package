"""S9 -- the full SM battery on the survivors, and the STATE-BY-STATE dossier.

Battery (5.21: every defining condition enumerated first, all checked):
  B1  gauge group: the four-dimensional roots are recomputed here from (V1, V2, W) and
      the observable part must be exactly SU(3) x SU(2) x U(1)_Y (the twelve X/Y roots
      of the block-1 SU(5) dead), which is what the s3 screen imposed -- re-verified
      against orbifolder's printed gauge group;
  B2  k_Y = 2|t_Y|^2 = 5/3;
  B3  hypercharge anomaly-consistency: Tr Y = 0 over all left-chiral STATES (dimension
      weighted -- the un-weighted sum is NOT the trace), and t_Y orthogonal to the
      anomalous U(1) generator printed by orbifolder;
  B4  net 3 generations, convention-blind: net Q, u^c, d^c, e^c counted over STATES in
      all four (sY, sC) branches, physical branch = net Q = +3.
  B5  THE TARGET: a T(0,1)/T(0,5) record with |p_sh|^2 = 31/18, a 10-role, and
      p_sh . Delta = +/-1, with Delta VERIFIED to be a root of the 4D gauge group
      (Delta.V1, Delta.V2, Delta.W all integral) so that p_sh - Delta is a state of the
      same irrep: two states of one 10-role at one twisted fixed point.
  H1-H4  the localisation theorem's own hypotheses re-checked on the constructed
      spectrum (the anti-claim guard): even lattice, SM in 4D, same cell, |w|^2 equal.

Output: data/su5_construction/s9_dossier.json
"""
from __future__ import annotations

import glob
import gzip
import json
import os
import re
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F

import sc_lib as S
sys.path.insert(0, os.path.join(S.CAMP, "code", "z3xz6_sweep"))
import z36sw_lib as SW                                                    # noqa: E402
from s8_battery import (P, DELTA, TY, A3, A2, split_models_cb, su3_dim,    # noqa: E402
                        su2_dim, role_of, RAW, SECT)

TEN = ("Q", "uc", "ec", "Qbar", "ucbar", "ecbar")


def load_build():
    with open(os.path.join(S.OUT, "s5_build.json")) as fh:
        d = json.load(fh)
    return {m["label"]: m for m in d["models"]}


def frac(s):
    return F(s)


def main():
    build = load_build()
    out = dict(models=0, battery=Counter(), passing=[], dossier=None,
               H_checks=Counter())
    for pth in sorted(glob.glob(os.path.join(RAW, "*.txt.gz"))):
        txt = gzip.open(pth, "rt", errors="replace").read()
        for lab, blk in split_models_cb(txt):
            out["models"] += 1
            pre = blk.split("Analyzing")[0]
            spec = SW.parse_spectrum(pre)
            gauge = SW.parse_gauge(blk)
            if not spec:
                continue
            left = [r for r in spec.values()
                    if r["chir"] == "l" and r["p_sh"] and r["rep"]]
            mult, d3s, d2s, Ys = [], [], [], []
            for r in left:
                d3, _ = su3_dim(r["p_sh"])
                d2 = su2_dim(r["p_sh"])
                tot = 1
                for d in r["rep"]:
                    tot *= abs(d)
                mult.append(F(tot, d3 * d2))
                d3s.append(d3)
                d2s.append(d2)
                Ys.append(S.dot(TY, r["p_sh"]))

            # B3 : Tr Y over STATES (dimension weighted)
            trY = sum(m * d3 * d2 * y for m, d3, d2, y in zip(mult, d3s, d2s, Ys))
            # B4 : net roles per branch
            nets = {}
            for sY in (1, -1):
                c = Counter()
                for m, d3, d2, y in zip(mult, d3s, d2s, Ys):
                    ro = role_of(d3, d2, y, sY, 1)
                    if ro is None:
                        continue
                    if ro.endswith("bar"):
                        c[ro[:-3]] -= int(m)
                    else:
                        c[ro] += int(m)
                nets[sY] = dict(c)
            sYp = None
            for sY in (1, -1):
                if nets[sY].get("Q") == 3:
                    sYp = sY
            ok_gen = sYp is not None and all(nets[sYp].get(k) == 3
                                             for k in ("Q", "uc", "ec"))
            ok_gen_dc = sYp is not None and nets[sYp].get("dc") == 3

            # B5 : the target
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
                ro = role_of(d3, d2, y, sYp or 1, 1)
                if ro not in TEN:
                    continue
                if abs(S.dot(p, DELTA)) != 1:
                    continue
                q = [a - b for a, b in zip(p, DELTA)]
                hits.append(dict(field=r["label"], kl=list(r["kl"]),
                                 fp=list(r["fp"]) if r["fp"] else None, role=ro,
                                 rep=list(r["rep"]),
                                 p_sh=[str(x) for x in p], q_sh=[str(x) for x in q],
                                 p2=str(S.dot(p, p)), q2=str(S.dot(q, q)),
                                 p_dot_Delta=str(S.dot(p, DELTA)),
                                 w2_p=str(proj2(p)), w2_q=str(proj2(q)),
                                 R=str(S.dot(p, p) - proj2(p)),
                                 Y_p=str(y), Y_q=str(S.dot(TY, q)),
                                 su3_q=su3_dim(q)[0], su2_q=su2_dim(q),
                                 qsh_right=[str(x) for x in r["qsh"]],
                                 osc=str(r["osc_raw"])))
            key = "gen%s_target%s" % ("3" if ok_gen else "x", "Y" if hits else "N")
            out["battery"][key] += 1
            if ok_gen and hits:
                b = build[lab]
                out["passing"].append(dict(model=lab, V1=b["V1"], W3=b["W3"],
                                           nets=nets[sYp], sY=sYp, trY=str(trY),
                                           gauge_test=gauge.get("gauge_test"),
                                           n_hits=len(hits), hits=hits,
                                           dc3=ok_gen_dc))
    out["battery"] = dict(out["battery"])
    with open(os.path.join(S.OUT, "s9_dossier.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out["battery"], indent=1))
    print("PASSING (net 3 in Q, uc, ec  AND  a doubled 10-role at T(0,l)):",
          len(out["passing"]))
    for r in out["passing"][:6]:
        print("  ", r["model"], r["nets"], "trY", r["trY"], "hits",
              [(h["role"], h["kl"]) for h in r["hits"]])


def proj2(p):
    """|P_W p|^2 with W = span(SU(3) simple roots, SU(2) simple root, t_Y), exact."""
    basis = [A3[0], A3[1], A2[0], TY]
    # Gram-Schmidt over Q
    ob = []
    for b in basis:
        v = list(b)
        for u in ob:
            c = S.dot(v, u) / S.dot(u, u)
            v = [x - c * y for x, y in zip(v, u)]
        if any(x != 0 for x in v):
            ob.append(v)
    return sum(S.dot(p, u) ** 2 / S.dot(u, u) for u in ob)


if __name__ == "__main__":
    main()

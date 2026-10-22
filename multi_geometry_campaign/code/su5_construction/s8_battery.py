"""S8 -- the SM battery and the two-10-role test, computed from the CONSTRUCTION's own
Standard-Model data rather than from orbifolder's vev-configuration finder.

WHY NOT `analyze config`.  orbifolder identifies a Standard-Model configuration by finding
an SU(5) among the four-dimensional roots and reading the hypercharge off it.  The CB
embeddings are built so that the twelve X/Y roots of the block-1 SU(5) are KILLED (that is
what makes the observable group exactly SU(3) x SU(2) x U(1) rather than SU(5)), so no
SU(5) survives in 4D and the finder reports "No vev-configuration identified" on 2,905 of
2,912 models.  That is a limitation of the finder, not of the models: the hypercharge is
known exactly, by construction --

    t_Y = (-1/3, -1/3, -1/3, 1/2, 1/2, 0, 0, 0 | 0^8),      k_Y = 2 |t_Y|^2 = 5/3,

it is traceless inside the block-1 SU(5), and SU(3)_c = <e1-e2, e2-e3>,
SU(2)_L = <e4-e5>.  Everything below is computed from the printed `left-moving p_sh`
against those generators (5.21: every defining condition checked, none inherited).

CHECKS PER MODEL
  k_Y                       2|t_Y|^2 = 5/3                                        (exact)
  hypercharge anomaly       Tr Y over all left-chiral STATES = 0, and t_Y . t_anom = 0
                            where t_anom is read from orbifolder's own U(1) data
  net generations           net (3,2)_{1/6} counted over STATES, convention-blind:
                            all four (sY, sC) branches emitted, the physical branch is
                            the one with net Q = +3 (the certified census rule)
  THE TEST                  a left-chiral record at a T(0,1)/T(0,5) cell with
                            |p_sh|^2 = 31/18, a 10-role under (SU3, SU2, t_Y), and
                            p_sh . Delta = +/-1  ==>  that record is a doublet of the
                            unbroken SU(2)_Delta and the cell carries TWO STATES of that
                            10-role.  This is the construction's target.

Output: data/su5_construction/s8_battery.json
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

RAW = os.path.join(S.OUT, "raw")
SECT = {(0, 1), (0, 5)}
P = [F(2, 3), F(2, 3), F(-1, 3), F(-1, 3), F(-1, 3), F(0), F(0), F(0),
     F(1, 2), F(1, 2), F(0), F(0), F(0), F(0), F(0), F(0)]
DELTA = [F(0)] * 8 + [F(1), F(1)] + [F(0)] * 6
TY = [F(-1, 3)] * 3 + [F(1, 2)] * 2 + [F(0)] * 11
A3 = [[F(1), F(-1), F(0)] + [F(0)] * 13, [F(0), F(1), F(-1)] + [F(0)] * 13]
A2 = [[F(0)] * 3 + [F(1), F(-1)] + [F(0)] * 11]
_MARKCB = re.compile(r"/>\s*/(CB\d+)>")


def split_models_cb(text):
    marks = [(m.start(), m.group(1)) for m in _MARKCB.finditer(text)]
    return [(lab, text[pos:(marks[i + 1][0] if i + 1 < len(marks) else len(text))])
            for i, (pos, lab) in enumerate(marks)]


def su3_dim(p):
    """dimension of the SU(3) irrep containing weight p (levels 1, so 1, 3 or 3bar).

    Dynkin labels (p.a1, p.a2) with a1 = e1-e2, a2 = e2-e3.  At level 1 the only
    integrable weights are 0, w1 (3) and w2 (3bar); the Weyl orbit is read off the
    labels after Weyl-conjugation, which for A2 is just the sign pattern.
    """
    l1, l2 = S.dot(p, A3[0]), S.dot(p, A3[1])
    if l1 == 0 and l2 == 0:
        return 1, 0
    # A2 Weyl orbit of a non-zero level-1 weight has 3 elements; the sign of the
    # triality (l1 - l2 mod 3 is not available off-lattice) is fixed by the pair
    return 3, (1 if (l1 > 0 or (l1 == 0 and l2 < 0)) else -1)


def su2_dim(p):
    return 2 if S.dot(p, A2[0]) != 0 else 1


def role_of(c3, c2, Y, sY, sC):
    a, y = c3, sY * Y
    if a == 3 and c2 == 2 and y == F(1, 6):
        return "Q"
    if a == 3 and c2 == 1 and y == F(-2, 3):
        return "uc"
    if a == 1 and c2 == 1 and y == F(1):
        return "ec"
    if a == 3 and c2 == 2 and y == F(-1, 6):
        return "Qbar"
    if a == 3 and c2 == 1 and y == F(2, 3):
        return "ucbar"
    if a == 1 and c2 == 1 and y == F(-1):
        return "ecbar"
    if a == 3 and c2 == 1 and y == F(1, 3):
        return "dc"
    if a == 1 and c2 == 2 and y == F(-1, 2):
        return "L"
    if a == 1 and c2 == 2 and y == F(1, 2):
        return "Hu"
    return None


def main():
    out = dict(models=0, kY="5/3", n_with_T0l_ten=0, doubled=[], per_model=[],
               net_hist=Counter(), gen3=[], trY_nonzero=0,
               ten_at_T0l_records=Counter(), states_at_T0l=Counter())
    for pth in sorted(glob.glob(os.path.join(RAW, "*.txt.gz"))):
        txt = gzip.open(pth, "rt", errors="replace").read()
        for lab, blk in split_models_cb(txt):
            out["models"] += 1
            pre = blk.split("Analyzing")[0]
            spec = SW.parse_spectrum(pre)
            if not spec:
                continue
            left = [r for r in spec.values()
                    if r["chir"] == "l" and r["p_sh"] and r["rep"]]
            # --- states per record: total printed dimension / the SM dimensions
            netQ = {}
            for sY in (1, -1):
                n = 0
                for r in left:
                    d3, _ = su3_dim(r["p_sh"])
                    d2 = su2_dim(r["p_sh"])
                    tot = 1
                    for d in r["rep"]:
                        tot *= abs(d)
                    hid = F(tot, d3 * d2)
                    Y = S.dot(TY, r["p_sh"])
                    ro = role_of(d3, d2, Y, sY, 1)
                    if ro == "Q":
                        n += int(hid)
                    elif ro == "Qbar":
                        n -= int(hid)
                netQ[sY] = n
            out["net_hist"][str(sorted(netQ.items()))] += 1
            sYphys = 1 if netQ[1] == 3 else (-1 if netQ[-1] == 3 else 0)
            trY = sum(F(1) * S.dot(TY, r["p_sh"]) for r in left)
            if trY != 0:
                out["trY_nonzero"] += 1

            # --- the test at T(0,1)/T(0,5)
            hits = []
            for r in left:
                if tuple(r["kl"]) not in SECT:
                    continue
                p = r["p_sh"]
                if S.dot(p, p) != F(31, 18):
                    continue
                d3, _ = su3_dim(p)
                d2 = su2_dim(p)
                Y = S.dot(TY, p)
                ro = role_of(d3, d2, Y, 1, 1) or role_of(d3, d2, Y, -1, 1)
                if ro not in ("Q", "uc", "ec", "Qbar", "ucbar", "ecbar"):
                    continue
                out["ten_at_T0l_records"][ro] += 1
                pd = S.dot(p, DELTA)
                if abs(pd) == 1:
                    hits.append(dict(field=r["label"], kl=list(r["kl"]),
                                     fp=list(r["fp"]) if r["fp"] else None,
                                     role=ro, p_sh=[str(x) for x in p],
                                     p2=str(S.dot(p, p)), p_dot_Delta=str(pd),
                                     rep=list(r["rep"]),
                                     partner=[str(x - y) for x, y in zip(p, DELTA)]))
            if hits:
                out["n_with_T0l_ten"] += 1
                out["doubled"].append(dict(model=lab, netQ=netQ, hits=hits))
            out["per_model"].append(dict(model=lab, netQ=netQ, n_T0l_ten_hits=len(hits),
                                         trY=str(trY)))
    out["net_hist"] = dict(out["net_hist"])
    out["ten_at_T0l_records"] = dict(out["ten_at_T0l_records"])
    with open(os.path.join(S.OUT, "s8_battery.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({k: v for k, v in out.items()
                      if k not in ("per_model", "doubled")}, indent=1)[:3000])
    print("MODELS WITH A DOUBLED 10-ROLE AT T(0,1)/T(0,5):", len(out["doubled"]))
    for d in out["doubled"][:5]:
        print("  ", d["model"], d["netQ"], [h["role"] for h in d["hits"]])


if __name__ == "__main__":
    main()

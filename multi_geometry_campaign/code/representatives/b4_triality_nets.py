"""Referee #8, items 2.3 / 3: TRIALITY-AWARE, hidden-dimension-weighted net chiral indices
I_X = sum_M d_hidden(M) [1_{M~X} - 1_{M~Xbar}],  X in {Q, u^c, e^c, d^c, L},
computed from the printed left-moving p_sh of each left-chiral record, with the colour
representation read from the momentum (traceless part of coordinates 1-3 is a weight of
the 3 iff its entries are {2/3,-1/3,-1/3}, of the 3bar iff {-2/3,1/3,1/3}), the isospin
from p.(e4-e5) != 0, Y = t_Y.p with t_Y of eq. (7.1).  This is the rule the certified gate
(gate_su5_construction.md DEFECT F-1) adopted; s8_battery.role_of is triality-blind and is
NOT used for the role here.

Usable as a library: nets(spec_records) -> dict.  Run as a script on b2_out.txt it prints
the nets of every brother deck and writes b4_nets.json.
"""
import json
import os
_PKG = os.environ.get("PAPERA_PKG") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
import re
import sys
from fractions import Fraction as F

sys.path.insert(0, os.path.join(_PKG, "multi_geometry_campaign", "code", "su5_construction"))
sys.path.insert(0, os.path.join(_PKG, "multi_geometry_campaign", "code", "z3xz6_sweep"))
import sc_lib as S  # noqa: E402
import z36sw_lib as SW  # noqa: E402

D = os.path.join(_PKG, "multi_geometry_campaign", "data", "representatives")
TY = [F(-1, 3)] * 3 + [F(1, 2)] * 2 + [F(0)] * 11
ROLE = {(3, 2, F(1, 6)): "Q", (-3, 1, F(-2, 3)): "uc", (1, 1, F(1)): "ec",
        (-3, 1, F(1, 3)): "dc", (1, 2, F(-1, 2)): "L"}


def colour(p):
    m = (p[0] + p[1] + p[2]) / 3
    c = sorted([p[0] - m, p[1] - m, p[2] - m])
    if c == [0, 0, 0]:
        return 1
    if c == [F(-1, 3), F(-1, 3), F(2, 3)]:
        return 3
    if c == [F(-2, 3), F(1, 3), F(1, 3)]:
        return -3
    return None


def role(p):
    c3 = colour(p)
    c2 = 2 if (p[3] - p[4]) != 0 else 1
    Y = S.dot(TY, p)
    r = ROLE.get((c3, c2, Y))
    if r:
        return r, +1
    cc = -c3 if c3 in (3, -3) else c3
    r = ROLE.get((cc, c2, -Y))
    if r:
        return r, -1
    return None, 0


def nets(records):
    out = {k: 0 for k in ("Q", "uc", "ec", "dc", "L")}
    for r in records:
        if r["chir"] != "l" or not r["p_sh"] or not r["rep"]:
            continue
        p = r["p_sh"]
        tot = 1
        for d in r["rep"]:
            tot *= abs(d)
        c3 = colour(p)
        dsm = (3 if c3 in (3, -3) else 1) * (2 if (p[3] - p[4]) != 0 else 1)
        rl, sg = role(p)
        if rl:
            out[rl] += sg * (tot // dsm)
    return out


def split(text, pat=r"/> /(BR\d+)>"):
    marks = [(m.start(), m.group(1)) for m in re.finditer(pat, text)]
    return [(lab, text[pos:(marks[i + 1][0] if i + 1 < len(marks) else len(text))])
            for i, (pos, lab) in enumerate(marks)]


def main():
    text = open(os.path.join(D, "b2_out.txt"), errors="replace").read()
    res = {}
    for lab, blk in split(text):
        spec = SW.parse_spectrum(blk.split("print summary")[0])
        res[lab] = nets(spec.values())
        print(lab, res[lab])
    json.dump(res, open(os.path.join(D, "b4_nets.json"), "w"), indent=1)


if __name__ == "__main__":
    main()

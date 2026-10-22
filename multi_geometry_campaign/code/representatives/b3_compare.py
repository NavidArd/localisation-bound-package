"""Referee #8, item 2.3 (q.5): compare the orbifolder spectra of the brother decks
(b1_brothers.py, run by b2_run_orbifolder.sh) with the base model BR0 = CB1534.

Per model: number of left-chiral records; the multiset of (sector, fixed point, p_sh,
oscillators, right-mover, representation); the dimension-weighted net indices
I_X = sum_M d_hidden(M) [1_{M~X} - 1_{M~Xbar}] for X = Q, u^c, e^c, d^c, L (reviewer's
definition; the same rule as s8_battery.py of the certified census); the doubled T(0,1)
Q record (|p|^2 = 31/18, hidden dimension 2) of table 5 of paperA_v9d.

Output: review/data/paperA_referee8/b3_compare.json
"""
import json
import os
_PKG = os.environ.get("PAPERA_PKG") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
import re
import sys
from collections import Counter
from fractions import Fraction as F

sys.path.insert(0, os.path.join(_PKG, "multi_geometry_campaign", "code", "su5_construction"))
sys.path.insert(0, os.path.join(_PKG, "multi_geometry_campaign", "code", "z3xz6_sweep"))
import sc_lib as S  # noqa: E402
import z36sw_lib as SW  # noqa: E402
import s8_battery as B8  # noqa: E402

D = os.path.join(_PKG, "multi_geometry_campaign", "data", "representatives")
TY = B8.TY
ROLES = ("Q", "uc", "ec", "dc", "L")


def split(text):
    marks = [(m.start(), m.group(1)) for m in re.finditer(r"/> /(BR\d+)>", text)]
    return [(lab, text[pos:(marks[i + 1][0] if i + 1 < len(marks) else len(text))])
            for i, (pos, lab) in enumerate(marks)]


def analyse(blk):
    spec = SW.parse_spectrum(blk.split("print summary")[0])
    left = [r for r in spec.values() if r["chir"] == "l" and r["p_sh"] and r["rep"]]
    net = Counter()
    keyset = Counter()
    doubled = []
    for r in left:
        p = r["p_sh"]
        d3, _ = B8.su3_dim(p)
        d2 = B8.su2_dim(p)
        tot = 1
        for d in r["rep"]:
            tot *= abs(d)
        hid = F(tot, d3 * d2)
        Y = S.dot(TY, p)
        ro = B8.role_of(d3, d2, Y, 1, 1)
        if ro in ROLES:
            net[ro] += int(hid)
        elif ro and ro.endswith("bar") and ro[:-3] in ROLES:
            net[ro[:-3]] -= int(hid)
        keyset[(r["kl"], r["fp"], tuple(str(x) for x in p), r["osc_raw"],
                tuple(str(x) for x in r["qsh"]), r["rep"])] += 1
        if r["kl"] in ((0, 1), (0, 5)) and S.dot(p, p) == F(31, 18) and ro in (
                "Q", "uc", "ec", "Qbar", "ucbar", "ecbar") and hid == 2:
            doubled.append(dict(field=r["label"], kl=r["kl"], role=ro))
    gg = re.search(r'Gauge group in vev-configuration "TestConfig1": (.*)', blk)
    return dict(n_left=len(left), net={k: net[k] for k in ROLES}, keys=keyset,
                doubled=doubled, gauge=gg.group(1).strip() if gg else None)


def main():
    meta = json.load(open(os.path.join(D, "b1_decks.json")))
    text = open(os.path.join(D, "b2_out.txt"), errors="replace").read()
    res = {lab: analyse(b) for lab, b in split(text)}
    base = res["BR0"]["keys"]
    out = {}
    for lab, a in res.items():
        k = a["keys"]
        same = (k == base)
        common = sum((k & base).values())
        out[lab] = dict(meta=meta.get(lab), n_left=a["n_left"], gauge=a["gauge"],
                        net=a["net"], doubled_T0l=a["doubled"],
                        identical_to_BR0=same, records_shared_with_BR0=common,
                        records_only_here=sum((k - base).values()),
                        records_only_in_BR0=sum((base - k).values()))
    json.dump(out, open(os.path.join(D, "b3_compare.json"), "w"), indent=1, default=str)
    for lab, o in out.items():
        m = o["meta"] or {}
        print("%-5s %-4s lamV2=%-4s other=%-4s n_left=%3d same=%-5s only_here=%3d "
              "net=%s doubled=%s | %s" % (
                  lab, m.get("kind", ""), m.get("dot_V2_mod1", ""),
                  m.get("other_mod1", ""), o["n_left"], o["identical_to_BR0"],
                  o["records_only_here"], o["net"],
                  [(d["field"], d["role"]) for d in o["doubled_T0l"]], o["gauge"]))


if __name__ == "__main__":
    main()

"""gate_su5_bound D6 -- HAND verification of Z3xZ6_1_1 `010` / `001` cells.

Pulls 10 of the 90 cells out of the raw orbifolder census text, prints the actual printed
blocks (sector, fixed point n_a, left-moving p_sh, representation, U(1) charges) together
with the hypercharge generator and the selected (sY, sC) orientation branch, and
RE-DERIVES the Standard-Model role of every left-chiral field at the cell from
(colour dim, weak dim, Y = sY*sY0*tY.p_sh) without going through `role_census`'s output
dictionary.  The two readings are compared field by field.

Output: data/gate_su5_bound/d6_handcheck_Z3xZ6.json  (+ stdout transcript)
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"

import glob                                                             # noqa: E402
import gzip                                                             # noqa: E402
import json                                                             # noqa: E402
import re                                                               # noqa: E402
import sys                                                              # noqa: E402
from collections import defaultdict                                     # noqa: E402
from fractions import Fraction as F                                     # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CAMP = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(CAMP, "data", "gate_su5_bound")
sys.path.insert(0, os.path.join(CAMP, "code", "z3xz6_sweep"))
import z36sw_lib as Z                                                   # noqa: E402

CENSUS = os.path.join(CAMP, "data", "z3xz6_sweep", "raw", "census")
TEN = ("Q", "uc", "ec")

# independent role table, written out here rather than imported (5.21 (iii))
ROLE = {(3, 2, F(1, 6)): "Q", (-3, 1, F(-2, 3)): "uc", (1, 1, F(1)): "ec",
        (-3, 1, F(1, 3)): "dc", (1, 2, F(-1, 2)): "L", (1, 2, F(1, 2)): "Hu",
        (-3, 2, F(-1, 6)): "Qbar", (3, 1, F(2, 3)): "ucbar",
        (3, 1, F(-1, 3)): "dcbar", (1, 1, F(-1)): "ecbar"}


def raw_blocks(pre, kl, fp):
    """verbatim printed F_n blocks of the given (sector, fixed point)."""
    heads = [(m.start(), m.group(1))
             for m in re.finditer(r"(?m)^\s{2,8}(F_\d+)\s*$", pre)]
    out = []
    for i, (p, lab) in enumerate(heads):
        end = heads[i + 1][0] if i + 1 < len(heads) else min(len(pre), p + 3000)
        b = pre[p:end]
        ms = re.search(r"sector \(k,l\)\s*:\s*\(([^)]*)\)", b)
        mf = re.search(r"fixed point n_a\s*:\s*\(([^)]*)\)", b)
        if not ms or not mf:
            continue
        s = tuple(int(x) for x in ms.group(1).split(","))
        s = (s[0] % 3, s[1] % 6)
        f = tuple(int(x) for x in mf.group(1).split(","))
        if s == tuple(kl) and f == tuple(fp):
            out.append((lab, b.rstrip()))
    return out


def main(nwant=10, nseg=6):
    want = []
    report = dict(geometry="Z3xZ6_1_1", cells=[])
    for pth in sorted(glob.glob(os.path.join(CENSUS, "*.txt.gz")))[:nseg]:
        if len(want) >= nwant:
            break
        txt = gzip.open(pth, "rt", errors="replace").read()
        for lab, blk in Z.split_models(txt):
            if len(want) >= nwant:
                break
            pre = blk.split("Analyzing")[0]
            try:
                spec = Z.parse_spectrum(pre)
                gauge = Z.parse_gauge(blk)
                if gauge["tY"] is None or not spec:
                    continue
                pos = Z.observable_positions(spec, gauge["tY"], gauge)
                if pos is None:
                    continue
                a3, a2, sY0, tY, _bi = pos
                roles, diag = Z.role_census(spec, tY, a3, a2, sY0)
            except Exception:                                           # noqa: BLE001
                continue
            cells = defaultdict(list)
            for l, r in spec.items():
                if r["p_sh"] is None or r["chir"] != "l":
                    continue
                kl = tuple(r["kl"])
                if kl == (0, 0):
                    continue
                cells[(kl, tuple(r["fp"]) if r["fp"] else None)].append(r)
            for key, recs in cells.items():
                ros = set((roles.get(r["label"]) or {}).get("role") for r in recs)
                sig = "".join("1" if x in ros else "0" for x in TEN)
                if sig not in ("010", "001"):
                    continue
                sY, sC = diag["sY_selected"], diag["sC_selected"]
                fields = []
                agree = True
                for r in recs:
                    Y = sY * sY0 * sum(a * b for a, b in zip(tY, r["p_sh"]))
                    c3 = r["rep"][a3]
                    c3 = (sC * c3) if abs(c3) != 1 else 1
                    c2 = abs(r["rep"][a2])
                    mine = ROLE.get((c3, c2, Y))
                    theirs = (roles.get(r["label"]) or {}).get("role")
                    if mine != theirs:
                        agree = False
                    fields.append(dict(label=r["label"], rep=list(r["rep"]),
                                       c3=c3, c2=c2, Y=str(Y),
                                       p_sh=[str(x) for x in r["p_sh"]],
                                       qsh=[str(x) for x in (r["qsh"] or [])],
                                       osc=r.get("osc") or r.get("osc_raw"),
                                       role_independent=mine, role_code=theirs))
                report["cells"].append(dict(
                    model=lab, file=os.path.basename(pth), sector="T(%d,%d)" % key[0],
                    fixed_point_n_a=list(key[1]) if key[1] else None,
                    pattern=sig,
                    role_content=sorted(x for x in ros if x),
                    ten_roles_present=[x for x in TEN if x in ros],
                    Q_present=("Q" in ros), uc_present=("uc" in ros),
                    ec_present=("ec" in ros),
                    a3=a3, a2=a2, sY0=sY0, sY=sY, sC=sC,
                    tY=[str(x) for x in tY],
                    net_Q=diag["nets"]["sY=%d,sC=%d" % (sY, sC)],
                    independent_reading_agrees=agree,
                    fields=fields,
                    raw=[dict(label=l, block=b) for l, b in raw_blocks(pre, key[0], key[1])]))
                want.append(1)
                break                      # one cell per model, for spread
    with open(os.path.join(OUT, "d6_handcheck_Z3xZ6.json"), "w") as fh:
        json.dump(report, fh, indent=1, default=str)
    for c in report["cells"]:
        print("=" * 78)
        print("%s  %s  n_a=%s  pattern=%s  roles=%s  indep_agrees=%s"
              % (c["model"], c["sector"], c["fixed_point_n_a"], c["pattern"],
                 c["role_content"], c["independent_reading_agrees"]))
        for f in c["fields"]:
            print("   %-6s rep=%-22s (c3,c2,Y)=(%s,%s,%s)  role=%s / %s  osc=%s"
                  % (f["label"], f["rep"], f["c3"], f["c2"], f["Y"],
                     f["role_independent"], f["role_code"], f["osc"]))
    print("\ncells reported:", len(report["cells"]))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 10,
         int(sys.argv[2]) if len(sys.argv) > 2 else 6)

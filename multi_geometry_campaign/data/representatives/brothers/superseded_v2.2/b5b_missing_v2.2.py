"""Referee #8 item 2.3: the (model, a, b) brother slots that b5's small search (vectors
+-e_i +- e_j) left unbuilt.  Searches a larger lattice set (all E8 vectors with entries in
{-1,0,1} and up to four non-zero entries, plus the 128 spinor-type vectors (+-1/2)^8 of
even minus count, in either block) and, if still none passes the complete MI set, records
the obstruction.  New decks go to brothers/decks_9NN.txt for b6 to pick up.
Output: brothers/build_missing.json
"""
import json
import os
_PKG = os.environ.get("PAPERA_PKG") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
import sys
from itertools import combinations, product
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b5_build_brother_census as B5  # noqa: E402

S = B5.S
Z = B5.Z


def big_cands():
    out = []
    for blk in (0, 8):
        for k in (2, 4):
            for idx in combinations(range(8), k):
                for sg in product((1, -1), repeat=k):
                    v = [F(0)] * 16
                    for i, s in zip(idx, sg):
                        v[blk + i] = F(s)
                    out.append(v)
        for sg in product((1, -1), repeat=8):
            if sg.count(-1) % 2 == 0:
                v = [F(0)] * 16
                for i, s in enumerate(sg):
                    v[blk + i] = F(s, 2)
                out.append(v)
    return out


def main():
    bj = json.load(open(os.path.join(B5.OUT, "build.json")))
    miss = bj["missing"]
    A = dict(B5.read_decks(os.path.join(B5.CAMP, "su5_construction", "s5_models.txt")))
    Bd = dict(B5.read_decks(os.path.join(B5.CAMP, "gate_su5_construction",
                                         "caseB_models.txt")))
    C = big_cands()
    decks, meta, still = [], {}, []
    for lab, a, b in miss:
        rows = A.get(lab) or Bd.get(lab)
        V1, V2, W = rows[0], rows[1], rows[4]
        found = None
        for lam in C:
            t = B5.m1(S.dot(lam, V2)) * 3
            if t.denominator != 1 or int(t) != a:
                continue
            ok, bad = S.mi_report(B5.add(V1, lam), V2, [Z, Z, W, W, Z, Z])
            if ok:
                found = lam
                break
        if not found:
            # record why: which MI conditions fail for every lam with lam.V2 = a/3
            fails = set()
            for lam in C:
                t = B5.m1(S.dot(lam, V2)) * 3
                if t.denominator == 1 and int(t) == a:
                    fails |= set(S.mi_report(B5.add(V1, lam), V2, [Z, Z, W, W, Z, Z])[1])
            still.append(dict(model=lab, a=a, b=b, W_zero=all(x == 0 for x in W),
                              failing=sorted(fails)))
            continue
        nl = "X%sa%db%d" % (lab, a, b)
        decks.append("begin model\nLabel:%s\nSpaceGroup:Geometry/Geometry_Z3xZ6_1_1.txt"
                     "\nLattice:E8xE8\nShifts and Wilsonlines:\n%s\nend model\n"
                     % (nl, "\n".join(B5.fmt(r) for r in
                                      [B5.add(V1, found), V2, Z, Z, W, W, Z, Z])))
        meta[nl] = dict(parent=lab, cls="B", a=a, b=b, lam=[str(x) for x in found],
                        mu=[str(x) for x in Z])
    for k in range(0, len(decks), B5.SEG):
        with open(os.path.join(B5.OUT, "decks_9%02d.txt" % (k // B5.SEG)), "w",
                  newline="\n") as fh:
            fh.write("".join(decks[k:k + B5.SEG]))
    json.dump(dict(built=len(decks), meta=meta, still_missing=still),
              open(os.path.join(B5.OUT, "build_missing.json"), "w"), indent=1)
    print("missing", len(miss), "built", len(decks), "still", len(still), still[:3])


if __name__ == "__main__":
    main()

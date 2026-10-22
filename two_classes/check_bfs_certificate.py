"""Independent checker for paperA_two_classes_bfs_certificate.json.

Reads only the certificate (and nothing from paperA_two_classes_enum.py) and verifies, in exact
rational arithmetic, that:

  1. rows 0 and 1 are the two representatives printed in the paper,
         V2^A = (2/3,2/3,-1/3,-1/3,-1/3,0,0,0 | 1/2,1/2,0^6)
         V2^B = (2/3,2/3,-1/3,-1/3,-1/3,-1/2,-1/2,0 | 0^8);
  2. every row's momentum satisfies |p|^2 = 31/18, |P_W p|^2 = 6/5 (W = traceless part of
     coordinates 1..5 of the first E8), 6p in Lambda and 3p, 2p not in Lambda (order six);
  3. every generator is an automorphism of the declared kind: a root r of E8 (|r|^2 = 2, r in
     E8) acting on one block -- on the first block either orthogonal to W or one of the three
     Standard-Model simple roots e1-e2, e2-e3, e4-e5 -- or the global -1;
  4. for every non-seed row, g(parent momentum) - momentum lies in Lambda = E8 + E8, so the
     row's class is the image of its parent's class under g;
  5. following parents from any row ends at row 0 or row 1 (no cycles, no orphans);
  6. the momenta represent pairwise DISTINCT classes of Lambda (no two differ by a lattice
     vector), so the certificate covers exactly n_classes classes.

Usage:  python check_bfs_certificate.py      (about one minute)
"""
import json
import os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.path.join(HERE, "paperA_two_classes_bfs_certificate.json")


def in_e8(b):
    if not (all(x.denominator == 1 for x in b) or all(x.denominator == 2 for x in b)):
        return False
    return sum(b) % 2 == 0


def in_lat(v):
    return in_e8(v[:8]) and in_e8(v[8:])


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def proj_w2(p):
    f = p[:5]
    m = sum(f) / 5
    return sum((x - m) ** 2 for x in f)


def main():
    c = json.load(open(CERT))
    rows = c["rows"]
    P = [[F(x) for x in r["momentum"]] for r in rows]
    VA = [F(2, 3), F(2, 3), F(-1, 3), F(-1, 3), F(-1, 3), F(0), F(0), F(0),
          F(1, 2), F(1, 2)] + [F(0)] * 6
    VB = [F(2, 3), F(2, 3), F(-1, 3), F(-1, 3), F(-1, 3), F(-1, 2), F(-1, 2), F(0)] + [F(0)] * 8
    assert P[0] == VA and P[1] == VB, "rows 0/1 are not the printed representatives"
    assert rows[0]["parent"] is None and rows[1]["parent"] is None
    for p in P:
        assert dot(p, p) == F(31, 18) and proj_w2(p) == F(6, 5)
        assert in_lat([6 * x for x in p])
        assert not in_lat([3 * x for x in p]) and not in_lat([2 * x for x in p])
    n_check = 0
    for i, r in enumerate(rows):
        if r["parent"] is None:
            continue
        g = r["generator"]
        q = P[r["parent"]]
        if g["kind"] == "global -1":
            img = [-x for x in q]
        else:
            rt = [F(x) for x in g["root"]]
            assert len(rt) == 8 and dot(rt, rt) == 2 and in_e8(rt)
            if g["kind"] == "reflection, first E8":
                f = rt[:5]
                sm_simple = [[1, -1, 0, 0, 0, 0, 0, 0], [0, 1, -1, 0, 0, 0, 0, 0],
                             [0, 0, 0, 1, -1, 0, 0, 0]]
                is_sm = any(rt == [F(s * x) for x in a] for a in sm_simple for s in (1, -1))
                assert is_sm or all(x == f[0] for x in f), \
                    "first-block root neither perpendicular to W nor a Standard-Model simple root"
                cc = dot(q[:8], rt)
                img = [x - cc * y for x, y in zip(q[:8], rt)] + q[8:]
            elif g["kind"] == "reflection, second E8":
                cc = dot(q[8:], rt)
                img = q[:8] + [x - cc * y for x, y in zip(q[8:], rt)]
            else:
                raise AssertionError("unknown generator kind")
        assert in_lat([a - b for a, b in zip(img, P[i])]), "row %d: g(parent) not in class" % i
        n_check += 1
    for i in range(len(rows)):
        seen = set()
        j = i
        while rows[j]["parent"] is not None:
            assert j not in seen, "cycle at row %d" % i
            seen.add(j)
            j = rows[j]["parent"]
        assert j in (0, 1)
    # pairwise distinct classes: compare canonical residues mod Lambda via a key
    def key(v):
        out = []
        for blk in (v[:8], v[8:]):
            cands = []
            for t in (F(0), F(1, 2)):
                y = [x - t for x in blk]
                n = [x.numerator // x.denominator for x in y]
                cands.append((tuple(a - b for a, b in zip(y, n)), sum(n) % 2))
            out.append(min(cands))
        return tuple(out)
    keys = [key(p) for p in P]
    assert len(set(keys)) == len(P), "two rows represent the same class"
    # the key is exact: spot-confirm that equal keys <=> difference in Lambda on random pairs
    import random
    random.seed(1)
    for _ in range(2000):
        a, b = random.randrange(len(P)), random.randrange(len(P))
        assert (keys[a] == keys[b]) == in_lat([x - y for x, y in zip(P[a], P[b])])
    print("certificate verified: %d classes, %d parent links checked, all chains end at the "
          "two printed representatives, all classes distinct" % (len(P), n_check))


if __name__ == "__main__":
    main()

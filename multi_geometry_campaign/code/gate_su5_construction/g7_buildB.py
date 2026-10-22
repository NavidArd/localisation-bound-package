"""V2 (third half) -- build LITERAL MI-passing embeddings for the NEW V2 class (case B)
and emit an orbifolder deck, so the three-generation screen can be run at its true scope.

MI1 and MI6 are class-invariant; MI3, MI4 and MI5 are representative-dependent, so a
representative is searched for each class:  x' = x + 3*lambda, z' = z + 3*mu with
lambda, mu in Lambda.  Because y = 6*V2 is primitive in the unimodular Lambda there is
lambda0 with lambda0 . y = 1, and x . y = 1 mod 3 then guarantees a j with
(x + 3 j lambda0) . y = 1 mod 12 -- and x^2 mod 6 is unchanged.  MI4 is fixed last by a
search over a pool of lattice vectors orthogonal-enough to keep MI5.

Every emitted model is re-checked with the FULL MI battery before it is written.

Output: data/gate_su5_construction/g7_buildB.json, caseB_models.txt
"""
from __future__ import annotations

import json
import os
import random
import sys
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gs_lib as G
from g4_lift import E8B

Z16 = [F(0)] * 16


def lat_vec(coeffs, blockshift=0):
    v = [F(0)] * 16
    for k, c in enumerate(coeffs):
        for t in range(8):
            v[blockshift + t] += c * E8B[k][t]
    return v


def find_lambda0(y):
    """a lattice vector with lambda0 . y = 1 (exists: Lambda unimodular, y primitive)."""
    import itertools
    basis = [lat_vec([1 if k == i else 0 for k in range(8)], 0) for i in range(8)] + \
            [lat_vec([1 if k == i else 0 for k in range(8)], 8) for i in range(8)]
    vals = [G.dot(b, y) for b in basis]
    from math import gcd
    g = 0
    for v in vals:
        g = gcd(g, int(v))
    assert g == 1, ("y not primitive", vals)
    # extended gcd over the 16 basis dots
    cur, comb = 0, [0] * 16
    for i, v in enumerate(vals):
        v = int(v)
        if v == 0:
            continue
        if cur == 0:
            cur, comb = v, [1 if k == i else 0 for k in range(16)]
            continue
        a, b = cur, v
        x0, x1, y0, y1 = 1, 0, 0, 1
        while b:
            q = a // b
            a, b = b, a - q * b
            x0, x1 = x1, x0 - q * x1
            y0, y1 = y1, y0 - q * y1
        comb = [x0 * c for c in comb]
        comb[i] += y0
        cur = a
    assert cur == 1
    v = [F(0)] * 16
    for k, c in enumerate(comb):
        if c:
            v = G.add(v, G.smul(F(c), basis[k]))
    assert G.dot(v, y) == 1 and G.in_lattice(v)
    return v


def main():
    lift = json.load(open(os.path.join(G.OUT, "g4_lift_B.json")))
    p = [G.rat(t) for t in lift["V2"]]
    y = [6 * t for t in p]
    l0 = find_lambda0(y)
    rnd = random.Random(7)
    pool = []
    roots = G.e8xe8_roots()
    for _ in range(4000):
        v = [F(0)] * 16
        for _ in range(rnd.randrange(1, 4)):
            v = G.add(v, list(roots[rnd.randrange(len(roots))]))
        pool.append(v)

    ok, fail = [], []
    for m in lift["models"]:
        x = [G.rat(t) for t in m["x"]]
        z = [G.rat(t) for t in m["z"]]
        # MI3 : x . y = 1 mod 12
        for j in range(4):
            xx = G.add(x, G.smul(F(3 * j), l0))
            if int(G.dot(xx, y)) % 12 == 1:
                break
        else:
            fail.append((m["idx"], "MI3"))
            continue
        # MI5 : z . y = 0 mod 12
        for j in range(4):
            zz = G.add(z, G.smul(F(3 * j), l0))
            if int(G.dot(zz, y)) % 12 == 0:
                break
        else:
            fail.append((m["idx"], "MI5"))
            continue
        V1 = G.smul(F(1, 3), xx)
        W = G.smul(F(1, 3), zz)
        good, bad = G.mi_check(V1, p, [Z16, Z16, W, W, Z16, Z16])
        if not good:
            # fix MI4 (and anything else) by z' -> z' + 3*nu with nu . y = 0 mod 4
            done = False
            for nu in pool:
                if int(G.dot(nu, y)) % 4 != 0:
                    continue
                z2 = G.add(zz, G.smul(F(3), nu))
                W2 = G.smul(F(1, 3), z2)
                g2, b2 = G.mi_check(V1, p, [Z16, Z16, W2, W2, Z16, Z16])
                if g2:
                    W, good, bad, done = W2, True, {}, True
                    break
            if not done:
                fail.append((m["idx"], str(bad)))
                continue
        assert G.order_of(V1) in (1, 3) and G.order_of(W) in (1, 3)
        assert G.order_of(p) == 6
        ok.append(dict(label="GB%04d" % m["idx"], V1=[str(t) for t in V1],
                       V2=[str(t) for t in p], W3=[str(t) for t in W]))
    print("built:", len(ok), " failed:", len(fail), fail[:10])

    lines = []
    Zs = ", ".join(["0"] * 16)
    for mm in ok:
        rows = [", ".join(mm["V1"]), ", ".join(mm["V2"]), Zs, Zs,
                ", ".join(mm["W3"]), ", ".join(mm["W3"]), Zs, Zs]
        lines.append("begin model\nLabel:%s\nSpaceGroup:Geometry/Geometry_Z3xZ6_1_1.txt\n"
                     "Lattice:E8xE8\nShifts and Wilsonlines:\n%s\nend model\n"
                     % (mm["label"], "\n".join(rows)))
    open(os.path.join(G.OUT, "caseB_models.txt"), "w", newline="\n").write("".join(lines))
    json.dump(dict(n_built=len(ok), n_failed=len(fail), failures=fail, models=ok),
              open(os.path.join(G.OUT, "g7_buildB.json"), "w"), indent=1)
    print("wrote caseB_models.txt with", len(ok), "models")


if __name__ == "__main__":
    main()

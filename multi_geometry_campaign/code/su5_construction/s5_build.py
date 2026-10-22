"""S5 -- build literal, FULLY modular-invariant embeddings for every inequivalent
(V1, W) class pair over V2 = p, and emit the orbifolder model file.

The class-level screen (s3, s4) leaves MI3, MI4, MI5 to be fixed by choosing lattice
representatives; the representative freedom is x -> x + 3 L, z -> z + 3 M with
L, M in Lambda, which does not move any class, any gauge group, or any spectrum, but
does move the three representative-level modular-invariance quantities:

    MI3  x.y = 1 mod 12   (y = 6p)          -- moves by 3 (L.y)
    MI5  z.y = 0 mod 12                     -- moves by 3 (M.y)
    MI4  z.x = 0 mod 6                      -- moves by 3 (z.L + M.x + 3 M.L)

Since y is primitive in the unimodular Lambda, L.y runs over all of Z, so MI3 and MI5 are
always solvable given the class conditions x.y = 1 mod 3 and z.y = 0 mod 3; MI4 is then
fixed by a parity choice.  Rather than trust that argument, this script SEARCHES a pool of
lattice offsets and accepts a pair only when sc_lib.mi_report -- the complete orbifolder
condition set MI0..MI6 -- returns clean.  Failures are counted, not hidden.

Outputs: data/su5_construction/s5_models.txt   (orbifolder input, one block per model)
         data/su5_construction/s5_build.json   (index -> (classes, vectors), MI evidence)
"""
from __future__ import annotations

import json
import os
import random
import sys
from fractions import Fraction as F

import numpy as np

import sc_lib as S
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s3_screen import P, DELTA, TY, W_SM, Y6, e8_class_reps                # noqa: E402

random.seed(20260820)


def lattice_pool(n=400):
    """a pool of Lambda vectors: 0, the 480 roots, and random small combinations."""
    import s5lib as G5
    pool = [tuple([F(0)] * 16)]
    R = G5.e8xe8_roots()
    pool += [tuple(r) for r in R[:240]] + [tuple(r) for r in R[240:480]]
    for _ in range(n):
        v = [F(0)] * 16
        for _k in range(random.randint(1, 4)):
            r = R[random.randrange(len(R))]
            c = random.choice((-2, -1, 1, 2))
            v = [a + c * b for a, b in zip(v, r)]
        pool.append(tuple(v))
    return pool


def frac16(X2blk1, X2blk2):
    """class rep (2*x per block) -> V = x/3 as 16 Fractions."""
    return [F(int(c), 6) for c in X2blk1] + [F(int(c), 6) for c in X2blk2]


def fmt(v):
    return ", ".join(("%d" % x.numerator) if x.denominator == 1
                     else "%d/%d" % (x.numerator, x.denominator) for x in v)


def main():
    X2 = e8_class_reps()
    pairs = np.load(os.path.join(S.OUT, "s4_pairs.npy"))
    pool = lattice_pool()
    Z16 = [F(0)] * 16

    out, blocks = [], []
    nfail = 0
    for k, (a, b, i1, i2) in enumerate(pairs.tolist()):
        V1_0 = frac16(X2[a], X2[b])
        W_0 = frac16(X2[i1], X2[i2])
        got = None
        for L in pool:
            V1 = [u + v for u, v in zip(V1_0, L)]
            if S.order_of(V1) != 3 and S.order_of(V1) != 1:
                continue
            ok1, bad1 = S.mi_report(V1, P, [Z16, Z16, W_0, W_0, Z16, Z16])
            if "MI1" in bad1 or "MI3" in bad1 or "MI0_V1" in bad1:
                continue
            for M in pool:
                W = [u + v for u, v in zip(W_0, M)]
                if S.order_of(W) not in (1, 3):
                    continue
                ok, bad = S.mi_report(V1, P, [Z16, Z16, W, W, Z16, Z16])
                if ok:
                    got = (V1, W)
                    break
            if got:
                break
        if got is None:
            nfail += 1
            continue
        V1, W = got
        lab = "CB%04d" % len(blocks)
        blocks.append(dict(label=lab, classes=[int(a), int(b), int(i1), int(i2)],
                           V1=[str(x) for x in V1], W3=[str(x) for x in W],
                           ordV1=S.order_of(V1), ordW=S.order_of(W)))
        out.append("begin model\nLabel:%s\nSpaceGroup:Geometry/Geometry_Z3xZ6_1_1.txt\n"
                   "Lattice:E8xE8\nShifts and Wilsonlines:\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n"
                   "end model\n" % (lab, fmt(V1), fmt(P), fmt(Z16), fmt(Z16),
                                    fmt(W), fmt(W), fmt(Z16), fmt(Z16)))
        if len(blocks) % 200 == 0:
            print("  built %d / %d" % (len(blocks), len(pairs)), flush=True)

    with open(os.path.join(S.OUT, "s5_models.txt"), "w") as fh:
        fh.write("".join(out))
    with open(os.path.join(S.OUT, "s5_build.json"), "w") as fh:
        json.dump(dict(n_pairs=int(len(pairs)), n_built=len(blocks),
                       n_MI_representative_not_found=nfail,
                       V2_target=[str(x) for x in P], models=blocks), fh, indent=1)
    print("built %d of %d ; MI representative not found for %d"
          % (len(blocks), len(pairs), nfail))


if __name__ == "__main__":
    sys.path.insert(0, os.path.join(S.CAMP, "code", "su5_ten_bound"))
    main()

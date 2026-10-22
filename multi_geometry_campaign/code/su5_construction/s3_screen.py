"""S3 -- ROUTE B: the (V1, W) screen over the fixed target V2 = p.

STAGES (funnel counts reported at every stage)
---------------------------------------------
 0  class enumeration  : x = 3 V1 and z = 3 W in Lambda/3Lambda           3^16 each
 1  A1 SM-preserving   : the SM simple roots of the block-1 SU(5) integral on V1 and W
 2  A4 hidden-SU(2)    : Delta.x = Delta.z = 0 mod 3   (both doubled states survive)
 3  MI class conditions: x^2 = 2 mod 6 ; x.y = 1 mod 3 ; z^2 = 0 mod 6 ; z.y = 0 mod 3
                         (y = 6p; derived in the docstring below)
 4  V1 orbit reduction : quotient by the stabiliser  G0 = <sigma_alpha> for roots alpha
                         with alpha.p in Z, alpha _|_ W_SM, alpha _|_ Delta.
                         Block 1: the local SU(4) on coordinates 6,7,8 (12 roots, W(A3),
                         order 24).  Block 2: the local E7 (126 roots, W(E7)).
                         Both act blockwise, so orbits factorise.
 5  pair screen        : joint MI class condition z.x = 0 mod 3, and the 4D gauge group
                         H = { roots r : r.p, r.V1, r.W all in Z } analysed for
                         SU(3) x SU(2) x U(1)_Y x hidden with the SM the one in which p
                         is a 10-role.
 6  representative fix : explicit search for a literal (V1, W) satisfying the FULL
                         orbifolder condition set MI0..MI6 (sc_lib.mi_report).

WHY THE MI CLASS CONDITIONS ARE WHAT THEY ARE
  MI1  3(V1^2 - 2/9) in 2Z  <=>  x^2/3 - 2/3 in 2Z  <=>  x^2 = 2 mod 6.  Invariant under
       x -> x + 3L because 6 x.L in 6Z and 9 L^2 in 6Z (L^2 even).
  MI3  3(V1.V2 - 1/18) in 2Z  <=>  x.y = 1 mod 12  (y = 6p in Lambda, PRIMITIVE).
       Under x -> x + 3L the value moves by 3(L.y), and L.y ranges over all of Z because
       y is primitive in the self-dual Lambda; so a representative exists iff
       x.y = 1 mod 3, which IS class invariant.
  MI5  3(W.V2) in 2Z  <=>  z.y = 0 mod 12; representative exists iff z.y = 0 mod 3.
  MI4  3(W.V1) in 2Z  <=>  z.x = 0 mod 6; representative exists iff z.x = 0 mod 3.
  MI6  3 W^2 in 2Z    <=>  z^2 = 0 mod 6.  Invariant (same argument as MI1).

Outputs data/su5_construction/s3_screen.json (+ s3_candidates.json.gz)
"""
from __future__ import annotations

import gzip
import json
import os
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction as F

import numpy as np

import sc_lib as S

sys.path.insert(0, os.path.join(S.CAMP, "code", "su5_ten_bound"))
import s5lib as G5                                                       # noqa: E402

# ------------------------------------------------------------------ the target
P = [F(2, 3), F(2, 3), F(-1, 3), F(-1, 3), F(-1, 3), F(0), F(0), F(0),
     F(1, 2), F(1, 2), F(0), F(0), F(0), F(0), F(0), F(0)]
DELTA = [F(0)] * 8 + [F(1), F(1)] + [F(0)] * 6
TY = [F(-1, 3)] * 3 + [F(1, 2)] * 2 + [F(0)] * 11
SM_ROOTS = []
for (i, j) in ((0, 1), (1, 2), (3, 4)):
    r = [F(0)] * 16
    r[i], r[j] = F(1), F(-1)
    SM_ROOTS.append(r)
W_SM = SM_ROOTS + [TY]
Y6 = [6 * x for x in P]                       # y = 6p, in Lambda, primitive


def as_int(v, sc):
    return np.array([int(x * sc) for x in v], dtype=np.int64)


# ------------------------------------------------------------------ E8/3E8 classes
def e8_class_reps():
    """6561 representatives of E8/3E8, returned as 2*x (integer 8-vectors)."""
    A = np.array([
        [1, -1, -1, -1, -1, -1, -1, 1],
        [2, 2, 0, 0, 0, 0, 0, 0],
        [-2, 2, 0, 0, 0, 0, 0, 0],
        [0, -2, 2, 0, 0, 0, 0, 0],
        [0, 0, -2, 2, 0, 0, 0, 0],
        [0, 0, 0, -2, 2, 0, 0, 0],
        [0, 0, 0, 0, -2, 2, 0, 0],
        [0, 0, 0, 0, 0, -2, 2, 0],
    ], dtype=np.int64)
    coef = np.array(np.meshgrid(*([np.arange(3)] * 8), indexing="ij")).reshape(8, -1).T
    return coef @ A


def canon3(X2):
    """Canonical key of a class of E8/3E8, given as 2*x.

    2*(3 E8) = 6 D8  u  (6 D8 + (3,...,3)), so componentwise reduction mod 6 leaves
    exactly TWO candidate rows per class -- (2x mod 6) and ((2x + 3) mod 6) -- and the
    key is their lexicographic minimum.  (Reducing mod 6 alone is a FINER invariant than
    the class and silently splits classes; that bug was caught by a KeyError in the
    reflection permutation, i.e. by the group action failing to close.)
    """
    one = np.atleast_2d(X2)
    a = one % 6
    b = (one + 3) % 6
    take_b = np.zeros(len(a), bool)
    for k in range(a.shape[1]):
        undec = ~take_b & np.ones(len(a), bool)
        diff = a[:, k] != b[:, k]
        first = diff & ~np.any(a[:, :k] != b[:, :k], axis=1) if k else diff
        take_b |= first & (b[:, k] < a[:, k])
    out = np.where(take_b[:, None], b, a)
    return out[0] if np.ndim(X2) == 1 else out


def main():
    t0 = time.time()
    os.makedirs(S.OUT, exist_ok=True)
    rep = dict(stages={}, sizes={})

    X2 = e8_class_reps()                                     # 6561 x 8, = 2*x
    rep["sizes"]["classes_per_block"] = int(X2.shape[0])

    # ---------------- stage 1/2: per-block admissibility ------------------------
    # block 1 (SM-preserving): (x_i - x_j) = 0 mod 3 for (1,2),(2,3),(4,5)
    b1_ok = np.ones(len(X2), bool)
    for (i, j) in ((0, 1), (1, 2), (3, 4)):
        b1_ok &= (((X2[:, i] - X2[:, j]) // 2) % 3 == 0)
    # block 2 (Delta-integral): x_1 + x_2 = 0 mod 3
    b2_ok = (((X2[:, 0] + X2[:, 1]) // 2) % 3 == 0)
    rep["stages"]["block1_SM_preserving"] = int(b1_ok.sum())
    rep["stages"]["block2_Delta_integral"] = int(b2_ok.sum())

    # ---------------- stage 4: orbit reduction ---------------------------------
    roots8 = np.array([[int(2 * c) for c in r] for r in G5.e8_roots()], dtype=np.int64)

    def stab_roots_block1():
        out = []
        for r in roots8:                                    # r = 2*alpha
            al = r / 2.0
            # alpha . p in Z ; alpha _|_ W_SM
            if abs(sum(al[k] * float(P[k]) for k in range(8)) % 1.0) > 1e-9:
                continue
            good = True
            for w in W_SM:
                if abs(sum(al[k] * float(w[k]) for k in range(8))) > 1e-9:
                    good = False
                    break
            if good:
                out.append(r)
        return np.array(out, dtype=np.int64)

    def stab_roots_block2():
        out = []
        for r in roots8:
            al = r / 2.0
            if abs(sum(al[k] * float(P[8 + k]) for k in range(8)) % 1.0) > 1e-9:
                continue
            if abs(al[0] + al[1]) > 1e-9:                    # alpha _|_ Delta
                continue
            out.append(r)
        return np.array(out, dtype=np.int64)

    SR1, SR2 = stab_roots_block1(), stab_roots_block2()
    rep["sizes"]["stabiliser_roots_block1"] = int(len(SR1))
    rep["sizes"]["stabiliser_roots_block2"] = int(len(SR2))

    def orbits(mask, SR):
        """orbits of the reflection group <sigma_alpha : alpha in SR> on the classes
        selected by `mask`.  sigma_alpha(x) = x - (x.alpha) alpha ; in 2*x coordinates
        2x -> 2x - ((2x).(2alpha)/2) * (2 alpha) / 1 ... computed exactly below."""
        idx = np.nonzero(mask)[0]
        key = {tuple(canon3(X2[i])): i for i in idx}
        seen, orbs = set(), []
        for i in idx:
            k0 = tuple(canon3(X2[i]))
            if k0 in seen:
                continue
            stack, orb = [i], set()
            while stack:
                j = stack.pop()
                kj = tuple(canon3(X2[j]))
                if kj in orb:
                    continue
                orb.add(kj)
                xj = X2[j]
                for r in SR:
                    c = int(np.dot(xj, r))                   # = 2x . 2alpha = 4 x.alpha
                    assert c % 4 == 0 or c % 2 == 0
                    xa = c // 4 if c % 4 == 0 else None
                    if xa is None:                           # x.alpha half-integral?
                        xa = F(c, 4)
                        assert xa.denominator == 1
                        xa = int(xa)
                    y2 = xj - xa * r
                    ky = tuple(canon3(y2))
                    if ky not in orb and ky in key:
                        stack.append(key[ky])
                    elif ky not in orb:
                        # the reflected class must still be a class of E8/3E8 -- find it
                        kk = key.get(ky)
                        if kk is not None:
                            stack.append(kk)
            seen |= orb
            orbs.append(sorted(orb))
        return orbs

    o1 = orbits(b1_ok, SR1)
    o2 = orbits(b2_ok, SR2)
    rep["stages"]["block1_orbits"] = len(o1)
    rep["stages"]["block2_orbits"] = len(o2)
    rep["stages"]["V1_representatives"] = len(o1) * len(o2)
    print("block1 SM classes %d -> %d orbits ; block2 Delta classes %d -> %d orbits"
          % (b1_ok.sum(), len(o1), b2_ok.sum(), len(o2)), flush=True)


    # ================= stage 3/5: MI class filters and the pair screen ==========
    # exact integer arithmetic on X2 = 2*x
    def blk_dot(A2, B2):
        return (A2 @ B2) // 4 if False else None

    # per-block quantities, exact: x = X2/2
    # x.x  = (X2.X2)/4 ; x.v = (X2 . 2v)/4
    y1 = np.array([int(2 * c) for c in Y6[:8]], dtype=np.int64)   # 2*y block1
    y2 = np.array([int(2 * c) for c in Y6[8:]], dtype=np.int64)
    n1 = (X2 * X2).sum(1) // 4                     # x^2 per block-1 class
    d1 = (X2 @ y1) // 4                            # x.y  per block-1 class
    d2 = (X2 @ y2) // 4

    # index sets
    b1_idx = np.nonzero(b1_ok)[0]
    b2_idx = np.nonzero(b2_ok)[0]

    # V1 representatives = (orbit rep of block1) x (orbit rep of block2)
    key2i = {tuple(canon3(X2[i])): i for i in range(len(X2))}
    r1 = [key2i[o[0]] for o in o1]
    r2 = [key2i[o[0]] for o in o2]
    orb1_size = {key2i[o[0]]: len(o) for o in o1}
    orb2_size = {key2i[o[0]]: len(o) for o in o2}

    # --- V1 candidates: MI1 (x^2 = 2 mod 6) and MI3 (x.y = 1 mod 3)
    V1c = []
    for a in r1:
        for b in r2:
            if (n1[a] + n1[b]) % 6 != 2:
                continue
            if (d1[a] + d2[b]) % 3 != 1:
                continue
            V1c.append((a, b))
    rep["stages"]["V1_reps_after_MI1_MI3"] = len(V1c)

    # --- W candidates: MI6 (z^2 = 0 mod 6) and MI5 (z.y = 0 mod 3), full (no orbit cut)
    NB = len(b1_idx) * len(b2_idx)
    A = np.repeat(b1_idx, len(b2_idx))
    B = np.tile(b2_idx, len(b1_idx))
    okW = ((n1[A] + n1[B]) % 6 == 0) & ((d1[A] + d2[B]) % 3 == 0)
    A, B = A[okW], B[okW]
    rep["stages"]["W_space_A1_A4"] = int(NB)
    rep["stages"]["W_after_MI5_MI6"] = int(len(A))
    print("W candidates %d of %d" % (len(A), NB), flush=True)

    # --- bad roots: local roots at p with a Standard-Model component that are NOT
    #     one of the 8 SM roots.  Killing all of them <=> observable part is exactly
    #     SU(3) x SU(2) x U(1)_Y (and then k_Y = 5/3 and Y is SU(5)-embedded, hence
    #     automatically non-anomalous -- CLAUDE.md 5.21).
    allr = G5.e8xe8_roots()
    def frac_dot(u, v):
        return sum(a * b for a, b in zip(u, v))
    loc = [r for r in allr if frac_dot(r, P).denominator == 1]
    rep["sizes"]["local_roots_at_p"] = len(loc)
    smset = set()
    for (i, j) in ((0,1),(1,2),(3,4),(0,2)):
        for s in (1, -1):
            v = [F(0)]*16; v[i], v[j] = F(s), F(-s); smset.add(tuple(v))
    bad = []
    for r in loc:
        if any(frac_dot(r, w) != 0 for w in W_SM) and tuple(r) not in smset:
            bad.append(r)
    rep["sizes"]["SM_roots_kept"] = len(smset)
    rep["sizes"]["bad_XY_roots"] = len(bad)
    print("local roots %d ; SM roots %d ; bad X/Y roots %d" % (len(loc), len(smset), len(bad)), flush=True)

    # bad roots live in block 1 only (they are SU(5) roots); build 2*r block-1 arrays
    BR = np.array([[int(2 * c) for c in r[:8]] for r in bad], dtype=np.int64)
    # r . x = (2r . X2)/4
    RX = (X2 @ BR.T) // 4                          # 6561 x nbad, integer
    survives = (RX % 3 == 0)                       # bad root survives this block-1 class

    def mask_of(idx1):
        m = np.zeros(len(idx1), dtype=np.int64)
        s = survives[idx1]
        for k in range(len(bad)):
            m |= (s[:, k].astype(np.int64) << k)
        return m

    maskW = mask_of(A)
    maskV = mask_of(np.array([a for a, b in V1c], dtype=np.int64))

    # joint MI4 class condition:  z.x = 0 mod 3
    XdotX = None
    hits = []
    tot_pairs = 0
    for t, (a, b) in enumerate(V1c):
        mv = maskV[t]
        ok = (maskW & mv) == 0
        if not ok.any():
            continue
        # z.x mod 3 over the surviving W's
        zx = ((X2[A] @ X2[a]) // 4 + (X2[B] @ X2[b]) // 4) % 3
        ok &= (zx == 0)
        c = int(ok.sum())
        tot_pairs += c
        if c:
            hits.append(dict(v1_block1=int(a), v1_block2=int(b), n_W=c,
                             orb1=orb1_size[a], orb2=orb2_size[b]))
    rep["stages"]["pairs_SMexact_and_MI4"] = tot_pairs
    rep["stages"]["V1_reps_with_at_least_one_W"] = len(hits)
    print("PAIRS surviving the SM-exact + MI class screen: %d  (from %d V1 reps)"
          % (tot_pairs, len(hits)), flush=True)
    rep["hits_head"] = hits[:40]

    with open(os.path.join(S.OUT, "s3_screen.json"), "w") as fh:
        json.dump(rep, fh, indent=1)
    print(json.dumps(rep, indent=1))
    print("%.1fs" % (time.time() - t0))


if __name__ == "__main__":
    main()

"""V2 -- THE CLASS BOUNDARY.  Enumerate EVERY order-6 shift class V2 whose massless level
at B = 31/36 admits a doubled SU(5)-10 role -- not just the phase's target class.

=====================================================================================
DERIVATION (re-derived in this gate; nothing inherited)
=====================================================================================
T(0,l) has local twist omega = (1/6, 0, -1/6), so delta_c = (1/2)[eta(1-eta)] summed
= 5/36, the oscillator-free budget is B = 31/36 and |p_sh|^2 = 2B = 31/18.  The only
other level a 10-role could sit at is N = 1/6 (|p|^2 = 25/18); N = 1/3 gives
|p|^2 = 19/18 < |w|^2 = 6/5 and is empty.  The localisation theorem closes both mixed
channels, so a doubled 10 is oscillator-free on BOTH members (re-derived in g6).

L1 (verified directly against catalogue p_sh in g3) gives V_g(T(0,l)) = l*V2 exactly, so
the class of V2 IS the class of the massless weight p.  Hence

    V2 admits a doubled 10-role  <=>  V2 + Lambda contains p with |p|^2 = 31/18 whose
    SU(5) projection is a 10-weight, together with Delta in Lambda, |Delta|^2 = 2,
    Delta perp W_SM, p.Delta = +-1.

Decompose p = P_W p  +  (p.Delta/2) Delta  +  s  with s perp W_SM and perp Delta:
    |P_W p|^2 = 6/5 ,   R := 31/18 - 6/5 = 47/90 ,   |s|^2 = R - 1/2 = 1/45 .
Write p_{1..5} = m (1,1,1,1,1) + d, d the 10-weight deviation, and put Bv := 6 p_other,
so the five coordinates of 6p are Bv (three times) and Bv+6 (twice) and m = Bv/6 + 2/5.
The residual budget splits as

    R = 5 m^2  +  (p_6^2 + p_7^2 + p_8^2)  +  |p_B|^2  =  47/90.                (*)

TWO CLOSED-FORM LEMMAS THAT MAKE THE ENUMERATION FINITE AND TINY
----------------------------------------------------------------
LEMMA A (Delta in the SECOND E8).  Then p.Delta = p_B.Delta and Cauchy-Schwarz gives
1 = |p_B.Delta| <= |p_B| sqrt 2, so |p_B|^2 >= 1/2.  By (*) |p_B|^2 <= 47/90 < 1/2 + 1/36,
and 6 p_B lies in E8 so |6 p_B|^2 is an EVEN INTEGER in [18, 47/90*36 = 18.8] -- hence
|6 p_B|^2 = 18 exactly, Cauchy-Schwarz is saturated, and 6 p_B = +- 3 Delta.
Then |p_B|^2 = 1/2, so 5 m^2 + sum_{6,7,8} = 1/45, which forces Bv = -2 and p_6=p_7=p_8=0
(the only solution: 5m^2 = 1/45 already exhausts the budget).
==> the second-block case yields EXACTLY ONE class up to equivalence, and it is the
    phase's target p = (2/3,2/3,-1/3,-1/3,-1/3,0,0,0 | 1/2,1/2,0,0,0,0,0,0).

LEMMA B (Delta in the FIRST E8).  The 20 roots of Lambda_1 orthogonal to W_SM are the 12
of the form +-e_a +- e_b with a,b in {6,7,8} and the 8 of the form
(c,c,c,c,c, s_6/2, s_7/2, s_8/2) with c = +-1/2 (an A4' root system).  p.Delta involves
only the first block, so 1 = |p.Delta| <= sqrt(2) sqrt(5m^2 + sum_{6,7,8}), i.e.
5 m^2 + sum_{6,7,8} >= 1/2, and by (*) |p_B|^2 <= 47/90 - 1/2 = 1/45, so
|6 p_B|^2 <= 4/5 < 2 -- an even integer -- hence ZERO: p_B = 0 (mod nothing; exactly 0).
==> the first-block case has trivial second-block content, and 5m^2 + sum_{6,7,8} = 47/90
    exactly, which is an integer Diophantine problem in (Bv, c_6, c_7, c_8), enumerated
    exhaustively below.

The two cases are disjoint (|p_B|^2 = 1/2 vs 0).  The 10bar branch is the image of the 10
branch under the global -1 of W(E8) x W(E8), which fixes the SM data up to the (allowed)
conjugation of Y, and is therefore not enumerated separately.

EQUIVALENCE, declared once: G_1 = the reflections in the 8 SM roots and in the 20 roots
perpendicular to W_SM (this is exactly the subgroup of W(E8)_1 fixing the declared SM),
together with the global -1; G_2 = W(E8) on the second block.  The Delta-condition is
invariant under both.  Orbits factorise because the group acts blockwise.

Output: data/gate_su5_construction/g2_classes.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import product

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gs_lib as G

RTOT = F(47, 90)
NORM = F(31, 18)
W2 = F(6, 5)


def roots_perp_sm_block1():
    out = []
    for r in G.e8_roots():
        v = list(r) + [F(0)] * 8
        if all(G.dot(v, c) == 0 for c in G.SM_CARTAN):
            out.append(list(r))
    return out


PAIRS = [(i, j) for i in range(5) for j in range(i + 1, 5)]   # all ten 10-weights
# {i,j} in {0,1,2} -> u^c ;  one in {0,1,2} and one in {3,4} -> Q ;  {3,4} -> e^c.
# These are NOT related by the SM Weyl group S3 x S2, so all ten are enumerated.


def build_p(Bv, c, blk2, pair=(0, 1)):
    """p for the 10-weight with the given pair: 6p has Bv+6 at the pair, Bv elsewhere."""
    six = [(Bv + 6) if k in pair else Bv for k in range(5)] + [c[0], c[1], c[2]]
    return [x / 6 for x in six] + list(blk2)


def main():
    perp1 = roots_perp_sm_block1()
    assert len(perp1) == 20
    n_12 = sum(1 for r in perp1 if all(r[k] == 0 for k in range(5)))
    print("roots of Lambda_1 perpendicular to W_SM:", len(perp1),
          "( %d supported on coords 6,7,8 ; %d half-integer )" % (n_12, 20 - n_12))

    # ---------------------------------------------------------------- B values from (*)
    Bvals = []
    for num in range(-20, 1):
        Bv = F(num, 2)
        m = Bv / 6 + F(2, 5)
        rem = RTOT - 5 * m * m
        if rem >= 0:
            Bvals.append((Bv, m, rem))
    print("\nadmissible Bv = 6*p_other and residual budget rem = sum_{6,7,8} + |p_B|^2:")
    for Bv, m, rem in Bvals:
        print("   Bv = %-5s  m = %-7s  5m^2 = %-9s  rem = %s"
              % (Bv, m, 5 * m * m, rem))

    found = []

    # ---------------------------------------------------------------- LEMMA A
    print("\nLEMMA A  (Delta in the second E8):")
    print("   forced |6 p_B|^2 = 18 with 6 p_B = +-3 Delta, hence |p_B|^2 = 1/2,")
    print("   hence 5m^2 + sum_{6,7,8} = 47/90 - 1/2 = 1/45.")
    caseA = []
    for Bv, m, rem in Bvals:
        if rem != F(1, 2):
            continue                       # needs sum_{678} = 0 and 5m^2 = 1/45
        # sum_{6,7,8} must be exactly 0 -> c = 0 -> requires the integer branch
        if Bv.denominator != 1:
            continue
        c = (F(0), F(0), F(0))
        for pair in PAIRS:
            for r in G.e8_roots():
                blk2 = [x / 2 for x in r]   # p_B = +- Delta/2 ; sign absorbed by W(E8)
                p = build_p(Bv, c, blk2, pair)
                if G.n2(p) != NORM:
                    continue
                if not G.in_lattice([6 * x for x in p]):
                    continue
                caseA.append(p)
    print("   solutions (before class/orbit reduction):", len(caseA))
    for p in caseA:
        found.append(("A", p))

    # ---------------------------------------------------------------- LEMMA B
    print("\nLEMMA B  (Delta in the first E8):  p_B = 0 forced; solve")
    print("   5m^2 + (c6^2+c7^2+c8^2)/36 = 47/90 over c of the same integrality as Bv.")
    caseB = []
    for Bv, m, rem in Bvals:
        target = 36 * rem                    # = c6^2+c7^2+c8^2
        if target < 0:
            continue
        half = (Bv.denominator == 2)
        # 4*target must be a sum of three squares of integers of fixed parity
        t4 = 4 * target
        if t4.denominator != 1:
            print("   Bv = %-5s : needs (2c)^2-sum = %s  -> NOT an integer, no solution"
                  % (Bv, t4))
            continue
        t4 = int(t4)
        lim = int(t4 ** 0.5) + 1
        sols = []
        rng = range(-lim, lim + 1)
        for a in rng:
            for b in rng:
                for d in rng:
                    if a * a + b * b + d * d != t4:
                        continue
                    if half and not (a % 2 and b % 2 and d % 2):
                        continue
                    if (not half) and (a % 2 or b % 2 or d % 2):
                        continue
                    sols.append((F(a, 2), F(b, 2), F(d, 2)))
        nkept = 0
        for pair in PAIRS:
            for c in sols:
                p = build_p(Bv, c, [F(0)] * 8, pair)
                if G.n2(p) != NORM:
                    continue
                if not G.in_lattice([6 * x for x in p]):
                    continue
                hit = [r for r in perp1 if abs(G.dot(p[:8], r)) == 1]
                if not hit:
                    continue
                caseB.append((p, len(hit)))
                nkept += 1
        print("   Bv = %-5s : (2c)^2-sum = %-4d  %4d triples, %4d (pair,c) in Lambda "
              "with a valid Delta" % (Bv, t4, len(sols), nkept))
    for p, _ in caseB:
        found.append(("B", p))

    # ---------------------------------------------------------------- classes / orbits
    byclass = {}
    for tag, p in found:
        k = tuple(str(x) for x in G.lattice_class_key(p))
        byclass.setdefault(k, (tag, p))
    print("\ntotal solutions:", len(found), " distinct Lambda-classes:", len(byclass))

    # G_1 closure (block-1 reflections) and W(E8)_2 closure (block-2 reflections)
    gens1 = [list(r) for r in perp1]
    for a in (G.A3[0], G.A3[1], G.A2L):
        gens1.append(list(a[:8]))
    gens2 = [list(r) for r in G.e8_roots()]

    parent = {k: k for k in byclass}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def uni(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    keyset = set(byclass)
    for k, (tag, p) in list(byclass.items()):
        for r in gens1:
            c = G.dot(p[:8], r)
            if c == 0:
                continue
            q = [x - c * y for x, y in zip(p[:8], r)] + list(p[8:])
            kk = tuple(str(x) for x in G.lattice_class_key(q))
            if kk in keyset:
                uni(k, kk)
        for r in gens2:
            c = G.dot(p[8:], r)
            if c == 0:
                continue
            q = list(p[:8]) + [x - c * y for x, y in zip(p[8:], r)]
            kk = tuple(str(x) for x in G.lattice_class_key(q))
            if kk in keyset:
                uni(k, kk)
        # global -1
        q = [-x for x in p]
        kk = tuple(str(x) for x in G.lattice_class_key(q))
        if kk in keyset:
            uni(k, kk)
    orb = defaultdict(list)
    for k in byclass:
        orb[find(k)].append(k)
    print("inequivalent V2 classes (orbits under G_1 x W(E8)_2 x <-1>):", len(orb))

    phase_p = [F(2, 3), F(2, 3), F(-1, 3), F(-1, 3), F(-1, 3), F(0), F(0), F(0),
               F(1, 2), F(1, 2)] + [F(0)] * 6
    pk = tuple(str(x) for x in G.lattice_class_key(phase_p))
    where = find(pk) if pk in byclass else None
    print("phase target class present:", pk in byclass,
          "| its orbit size (classes):", len(orb[where]) if where else None)

    reps = []
    for root, members in orb.items():
        tag, p = byclass[root]
        # which Delta(s) work for this representative
        d1 = [r for r in perp1 if abs(G.dot(p[:8], r)) == 1]
        d2 = [r for r in G.e8_roots() if abs(G.dot(p[8:], r)) == 1]
        reps.append(dict(case=tag, n_classes_in_orbit=len(members),
                         p=[str(x) for x in p],
                         order=G.order_of(p), norm2=str(G.n2(p)),
                         w2=str(G.proj_wsm_norm2(p)),
                         role=G.role_of(p, 1),
                         MI2=str(G.N_ORD * (G.n2(p) - G.n2(G.V2_TWIST))),
                         n_Delta_block1=len(d1), n_Delta_block2=len(d2),
                         is_phase_target=(root == where),
                         second_block_zero=all(x == 0 for x in p[8:])))
    print("\n--- inequivalent V2 classes ---")
    for r in sorted(reps, key=lambda z: (z["case"], -z["n_classes_in_orbit"])):
        print("  case %s | %5d classes in orbit | role %-5s | ord %d | "
              "nDelta(1st,2nd)=(%d,%d) | phase=%s\n     p = %s"
              % (r["case"], r["n_classes_in_orbit"], r["role"], r["order"],
                 r["n_Delta_block1"], r["n_Delta_block2"], r["is_phase_target"],
                 ", ".join(r["p"])))

    out = dict(n_solutions=len(found), n_distinct_classes=len(byclass),
               n_inequivalent=len(orb), phase_target_present=pk in byclass,
               B_values=[[str(b), str(m), str(rr)] for b, m, rr in Bvals],
               n_caseA=len(caseA), n_caseB=len(caseB),
               representatives=reps,
               all_classes={"|".join(k): [byclass[k][0],
                                          [str(x) for x in byclass[k][1]]]
                            for k in byclass})
    json.dump(out, open(os.path.join(G.OUT, "g2_classes.json"), "w"), indent=1)
    print("\nwrote g2_classes.json")


if __name__ == "__main__":
    main()

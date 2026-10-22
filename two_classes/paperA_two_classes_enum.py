"""Fresh, independent enumeration of every order-six shift class V2 on Z3xZ6_1_1 whose
T(0,1) massless level admits a doubled SU(5)-10 role.

Written from the mathematics for paper15/review/paperA_two_classes.md.  Nothing is
imported from multi_geometry_campaign/code/{su5_construction,gate_su5_construction};
every object below (E8 membership, the 240 roots, the class key, the search box, the two
branches, the equivalence reduction) is built here from scratch, and the counts are
compared with the phase/gate numbers only at the very end.

CONVENTIONS, STATED ONCE (CLAUDE.md 5.19 rule (i))
--------------------------------------------------
C1  Lambda = E8 (+) E8 in the even coordinate system.  A block x in R^8 lies in E8 iff
    (all x_i in Z, or all x_i in Z + 1/2) AND sum(x) in 2Z.  Lambda is even; roots have
    norm^2 = 2.  16-vectors are written (block 1 | block 2) in the order orbifolder
    prints `left-moving p_sh`.
C2  Geometry Z3xZ6_1_1: twists v1 = (0, 1/3, -1/3) [theta, order 3] and
    v2 = (1/6, 0, -1/6) [omega, order 6], read off the primary geometry file at import.
    T(0,1) has delta_c = 5/36, oscillator-free budget B = 31/36, |p|^2 = 2B = 31/18.
C3  Standard Model, DECLARED (this is the scope hypothesis, not a derivation):
    SU(3)_c = <e1-e2, e2-e3>, SU(2)_L = <e4-e5>,
    t_Y = (-1/3,-1/3,-1/3, 1/2, 1/2, 0,0,0 | 0^8), k_Y = 2|t_Y|^2 = 5/3.
    W_SM = span{e1-e2, e2-e3, e4-e5, t_Y} is the traceless 4-space of coordinates 1..5,
    i.e. the Cartan subalgebra of the "block-1 SU(5)".  Its orthogonal complement inside
    coordinates 1..5 is the single direction (1,1,1,1,1).
C4  SIGN OF Delta.  For an ordered pair of states (p, p'') at one cell with the same
    W_SM projection we set Delta := p - p''.  Then |Delta|^2 = 2 and p.Delta = +1.  The
    unordered configuration is invariant under Delta -> -Delta, so the search condition
    is |p.Delta| = 1 and both signs are kept throughout.
C5  10-weights.  A 10-role's W_SM projection is w = (e_i + e_j) - (2/5)(1,1,1,1,1),
    i < j in {1..5}; |w|^2 = 6/5.  The ten pairs fall into three W_SM orbits
    ({i,j} in {1,2,3} -> u^c; one in {1,2,3} and one in {4,5} -> Q; {4,5} -> e^c) and all
    ten are enumerated, not one per orbit.  The 10bar branch is p -> -p and is handled by
    the global -1 in the equivalence group, not by a separate enumeration.

OUTPUT: paper15/review/code/paperA_two_classes_enum.json
"""
from __future__ import annotations

import json
import os
from collections import defaultdict
from fractions import Fraction as F
from itertools import product

HERE = os.path.dirname(os.path.abspath(__file__))
# --- reproduction-package root shim (see environment/PATCHES_TO_PACKAGE_CODE.md) ---
_PKG_ROOT = os.environ.get('PAPERA_PKG') or os.path.abspath(os.path.join(HERE, '..'))
GEOM = os.path.join(_PKG_ROOT, 'orbifolder-1.2.1', 'Geometry',
                    'Geometry_Z3xZ6_1_1.txt')

NORM = F(31, 18)          # |p|^2 at the oscillator-free level of T(0,1)
W2 = F(6, 5)              # |P_W p|^2 for every 10-role
RTOT = NORM - W2          # = 47/90, the residual budget
HALF = F(1, 2)


# ---------------------------------------------------------------- basic vector algebra
def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def n2(a):
    return dot(a, a)


def in_e8(blk):
    blk = list(blk)
    if all(x.denominator == 1 for x in blk):
        return sum(blk) % 2 == 0
    if all(x.denominator == 2 and x.numerator % 2 for x in blk):
        return sum(blk) % 2 == 0
    return False


def in_lattice(v):
    return in_e8(v[:8]) and in_e8(v[8:])


def order_of(v, cap=24):
    for n in range(1, cap + 1):
        if in_lattice([n * x for x in v]):
            return n
    return 0


# ---------------------------------------------------------------- E8 roots, built here
def e8_roots():
    R = []
    for i in range(8):
        for j in range(i + 1, 8):
            for si in (1, -1):
                for sj in (1, -1):
                    v = [F(0)] * 8
                    v[i], v[j] = F(si), F(sj)
                    R.append(tuple(v))
    for bits in range(256):
        s = [(-1) ** ((bits >> k) & 1) for k in range(8)]
        if s.count(-1) % 2 == 0:
            R.append(tuple(F(x, 2) for x in s))
    assert len(R) == 240
    for r in R:
        assert n2(r) == 2 and in_e8(r)
    return R


E8R = e8_roots()


def e8_shell(n2target):
    """Every x in E8 with |x|^2 = n2target (an even non-negative integer).  DFS over the
    two cosets Z^8 and (Z+1/2)^8, keeping only even coordinate sums."""
    assert n2target >= 0 and n2target % 2 == 0
    out = []
    for base in (F(0), F(1, 2)):
        cur = [F(0)] * 8
        top = int(n2target) + 1

        def rec(k, rem, twice_sum):
            if k == 8:
                if rem == 0 and twice_sum % 4 == 0:
                    out.append(tuple(cur))
                return
            for t in range(-2 * top, 2 * top + 1):
                val = base + t if base else F(t)
                if val * val > rem:
                    continue
                cur[k] = val
                rec(k + 1, rem - val * val, twice_sum + 2 * val)
            cur[k] = F(0)

        rec(0, F(n2target), F(0))
    return sorted(set(out))


# ---------------------------------------------------------------- Standard-Model data (C3)
A31 = [F(1), F(-1), F(0)] + [F(0)] * 13
A32 = [F(0), F(1), F(-1)] + [F(0)] * 13
A2L = [F(0)] * 3 + [F(1), F(-1)] + [F(0)] * 11
TY = [F(-1, 3)] * 3 + [F(1, 2)] * 2 + [F(0)] * 11
SM_CARTAN = [A31, A32, A2L, TY]
SM_ROOTS = [A31, A32, A2L]


def proj_w_norm2(p):
    """|P_{W_SM} p|^2 = |p_{1..5}|^2 - (sum p_{1..5})^2 / 5."""
    s = sum(p[:5])
    return sum(x * x for x in p[:5]) - s * s / 5


def roots_perp_sm_block1():
    out = []
    for r in E8R:
        v = list(r) + [F(0)] * 8
        if all(dot(v, c) == 0 for c in SM_CARTAN):
            out.append(list(r))
    return out


PERP1 = roots_perp_sm_block1()


# ---------------------------------------------------------------- exact class key mod Lambda
def _e8_class_key(blk):
    """Canonical label of blk + E8.  E8 = D8 u (D8 + s), s = (1/2)^8.

    Class mod D8 is (blk mod Z^8, parity of the integer part's coordinate sum); mod E8 we
    additionally identify that datum with the one obtained from blk - s.  This is EXACT:
    it does not quotient by e_i (which is in Z^8 but NOT in E8).
    """
    cands = []
    for t in (F(0), F(1, 2)):
        y = [x - t for x in blk]
        n = [F(v.numerator // v.denominator) for v in y]      # floor
        r = tuple(v - m for v, m in zip(y, n))
        cands.append((tuple(str(x) for x in r), int(sum(n)) % 2))
    return min(cands)


def class_key(v):
    return (_e8_class_key(v[:8]), _e8_class_key(v[8:]))


def coarse_key(v):
    """The gate's key (mod <Z^8, s> = D8*), reproduced for comparison only."""
    out = []
    for blk in (v[:8], v[8:]):
        r = tuple(x - F(x.numerator // x.denominator) for x in blk)
        h = tuple((y + F(1, 2)) - F((y + F(1, 2)).numerator // (y + F(1, 2)).denominator)
                  for y in r)
        out.append(tuple(str(z) for z in min(r, h)))
    return tuple(out)


# ---------------------------------------------------------------- self-tests (5.19 rule (i))
def selftests():
    rep = {}
    # (a) E8 theta series 1, 240, 2160, 6720, 17520, 30240, ... = 240 * sigma_3(n)
    theta = [1]
    for n in range(1, 4):
        sh = e8_shell(2 * n)
        theta.append(len(sh))
        for x in sh:
            assert in_e8(x) and n2(x) == 2 * n
    exp = [1, 240, 2160, 6720]
    assert theta == exp, (theta, exp)
    rep["e8_theta_0_to_10"] = theta
    # (b) the exact class key does NOT merge x and x + e1 (which the coarse key does)
    a = [F(0)] * 16
    b = [F(0)] * 16
    b[0] = F(1)
    assert class_key(a) != class_key(b)
    rep["exact_key_separates_e1"] = True
    rep["coarse_key_merges_e1"] = (coarse_key(a) == coarse_key(b))
    # (c) the exact class key DOES identify x and x + lambda for random lambda in Lambda
    import random
    random.seed(11)
    ok = 0
    for _ in range(200):
        lam = [F(0)] * 16
        for _ in range(3):
            r = random.choice(E8R)
            off = random.choice([0, 8])
            for k in range(8):
                lam[off + k] += r[k]
        x = [F(random.randint(-6, 6), 6) for _ in range(16)]
        if class_key(x) == class_key([u + w for u, w in zip(x, lam)]):
            ok += 1
    assert ok == 200
    rep["exact_key_lattice_invariance_200_of"] = ok
    # (d) W_SM is the traceless 4-space of coords 1..5
    assert len(PERP1) == 20
    rep["roots_perp_WSM_block1"] = len(PERP1)
    rep["  of which supported on coords 6,7,8"] = sum(
        1 for r in PERP1 if all(r[k] == 0 for k in range(5)))
    # (e) the geometry's twists, read off the primary file
    rows, on = [], False
    for line in open(GEOM, encoding="utf-8", errors="replace"):
        s = line.strip()
        if s == "begin twist":
            on = True
            continue
        if s == "end twist":
            break
        if on and s:
            rows.append([F(t) if "/" not in t else F(int(t.split("/")[0]),
                                                     int(t.split("/")[1]))
                         for t in s.split()])
    v1, v2 = tuple(rows[0][1:]), tuple(rows[1][1:])
    assert v1 == (F(0), F(1, 3), F(-1, 3)) and v2 == (F(1, 6), F(0), F(-1, 6)), (v1, v2)
    dc = sum(HALF * e * (1 - e) for e in [x % 1 for x in v2])
    assert dc == F(5, 36), dc
    assert 2 * (1 - dc) == NORM
    rep["v1"], rep["v2"], rep["delta_c(0,1)"] = str(v1), str(v2), str(dc)
    # (f) MI2 is automatic at |p|^2 = 31/18
    mi2 = 6 * (NORM - (F(1, 36) + F(1, 36)))
    assert mi2.denominator == 1 and mi2 % 2 == 0
    rep["MI2_value"] = str(mi2)
    return rep


# ---------------------------------------------------------------- the search box
def budget_box():
    """Every admissible Bv := 6 * p_other, with the box proved rather than assumed.

    p_{1..5} = m (1,1,1,1,1) + w  with w the 10-weight deviation, m = Bv/6 + 2/5.
    (*)  5 m^2 + (p6^2+p7^2+p8^2) + |p_B|^2 = R = 47/90, every term >= 0, so 5 m^2 <= R,
    i.e. |m| <= sqrt(47/450) = 0.32318..., i.e. Bv = 6(m - 2/5) in [-4.3391, -0.4609].
    6p in Lambda forces the eight block-1 coordinates of 6p to be all integers or all
    half-odd-integers, so Bv is in (1/2)Z.  The box therefore contains at most the
    half-integers of [-4.3391, -0.4609]; we scan the strictly larger range [-10, 0].
    """
    out = []
    for num in range(-20, 1):
        Bv = F(num, 2)
        m = Bv / 6 + F(2, 5)
        rem = RTOT - 5 * m * m
        if rem >= 0:
            out.append((Bv, m, rem))
    return out


def build_p(Bv, c, pB, pair):
    six = [(Bv + 6) if k in pair else Bv for k in range(5)] + list(c)
    return [x / 6 for x in six] + list(pB)


PAIRS = [(i, j) for i in range(5) for j in range(i + 1, 5)]
ROLE = {}
for (i, j) in PAIRS:
    ROLE[(i, j)] = ("uc" if j < 3 else ("ec" if i >= 3 else "Q"))


def role_of_weight(p):
    """Which 10-role the W_SM projection of p carries: the deviation of (p1..p5) from its
    own mean is a permutation of (3/5,3/5,-2/5,-2/5,-2/5) (a 10-weight, pair positions
    carrying 3/5) or its negative (a 10bar-weight)."""
    a = sum(p[:5]) / 5
    d = [x - a for x in p[:5]]
    for sgn, suf in ((1, ""), (-1, "bar")):
        hi = [k for k in range(5) if sgn * d[k] == F(3, 5)]
        lo = [k for k in range(5) if sgn * d[k] == F(-2, 5)]
        if len(hi) == 2 and len(lo) == 3:
            i, j = hi
            base = "uc" if j < 3 else ("ec" if i >= 3 else "Q")
            return base + suf
    return None


def c_triples(target, half):
    """All (c6,c7,c8) with c6^2+c7^2+c8^2 = target, of the parity fixed by Bv."""
    if target < 0:
        return []
    t4 = 4 * target
    if t4.denominator != 1:
        return []
    t4 = int(t4)
    lim = int(t4 ** 0.5) + 2
    out = []
    for a in range(-lim, lim + 1):
        for b in range(-lim, lim + 1):
            for d in range(-lim, lim + 1):
                if a * a + b * b + d * d != t4:
                    continue
                if half and not (a % 2 and b % 2 and d % 2):
                    continue
                if (not half) and (a % 2 or b % 2 or d % 2):
                    continue
                out.append((F(a, 2), F(b, 2), F(d, 2)))
    return out


def main():
    rep = {"selftests": selftests()}
    box = budget_box()
    rep["box"] = [[str(b), str(m), str(r)] for b, m, r in box]
    print("admissible Bv (%d):" % len(box))
    for Bv, m, r in box:
        print("   Bv=%-6s m=%-8s 5m^2=%-10s rem=%s" % (Bv, m, 5 * m * m, r))

    found = []          # (branch, p, Delta)

    # ---- BRANCH 2: Delta in the E8 factor NOT containing the Standard Model ----------
    # p.Delta = p_B . delta = +-1 with 6 p_B in E8.  Write x = 6 p_B, so x.delta = +-6.
    # Cauchy-Schwarz: 36 |p_B|^2 = |x|^2 >= (x.delta)^2 / 2 = 18, and (*) gives
    # |p_B|^2 <= 47/90, i.e. |x|^2 <= 18.8.  |x|^2 is an even integer, hence exactly 18,
    # C-S is SATURATED, and x = +-3 delta.  Then |p_B|^2 = 1/2 and rem must be >= 1/2.
    print("\nBRANCH 2 (Delta in the second E8):  |6 p_B|^2 = 18 forced, 6 p_B = +-3 delta")
    br2 = []
    for Bv, m, rem in box:
        if rem < HALF:
            print("   Bv=%-6s rem=%-8s < 1/2  -> no solution" % (Bv, rem))
            continue
        tgt = 36 * (rem - HALF)                      # = c6^2+c7^2+c8^2
        half = (Bv.denominator == 2)
        cs = c_triples(tgt, half)
        n_here = 0
        for pair in PAIRS:
            for c in cs:
                for delta in E8R:
                    pB = [z / 2 for z in delta]
                    p = build_p(Bv, c, pB, pair)
                    if n2(p) != NORM:
                        continue
                    if not in_lattice([6 * x for x in p]):
                        continue
                    D = [F(0)] * 8 + list(delta)
                    if abs(dot(p, D)) != 1:
                        continue
                    br2.append((p, D))
                    n_here += 1
        print("   Bv=%-6s rem=%-8s target c^2-sum=%-6s %2d triples -> %5d solutions"
              % (Bv, rem, tgt, len(cs), n_here))
    for p, D in br2:
        found.append(("2", p, D))

    # ---- BRANCH 1: Delta in the E8 factor CONTAINING the Standard Model ---------------
    # Delta perp W_SM, so p.Delta only sees the part of p_{1..8} orthogonal to W_SM,
    # of norm^2 5m^2 + sum_{6,7,8}.  C-S: 1 <= sqrt2 sqrt(5m^2+sum), so that is >= 1/2,
    # and (*) gives |p_B|^2 <= 47/90 - 1/2 = 1/45, i.e. |6 p_B|^2 <= 4/5 < 2: an even
    # non-negative integer below 2 is 0, so p_B = 0 exactly and sum_{6,7,8} = rem.
    print("\nBRANCH 1 (Delta in the first E8):  p_B = 0 forced")
    br1 = []
    for Bv, m, rem in box:
        tgt = 36 * rem
        half = (Bv.denominator == 2)
        cs = c_triples(tgt, half)
        n_here = 0
        for pair in PAIRS:
            for c in cs:
                p = build_p(Bv, c, [F(0)] * 8, pair)
                if n2(p) != NORM:
                    continue
                if not in_lattice([6 * x for x in p]):
                    continue
                for r in PERP1:
                    D = list(r) + [F(0)] * 8
                    if abs(dot(p, D)) != 1:
                        continue
                    br1.append((p, D))
                    n_here += 1
        print("   Bv=%-6s rem=%-8s target c^2-sum=%-8s %3d triples -> %4d (p,Delta)"
              % (Bv, rem, tgt, len(cs), n_here))
    for p, D in br1:
        found.append(("1", p, D))

    # The C5 conjugate branch: p is enumerated with a 10-weight W_SM projection; the
    # 10bar branch is exactly the image under the global -1 of Aut(Lambda).  Adding it
    # explicitly (rather than assuming it) makes the solution set closed under the
    # declared equivalence group, which is what licenses the union-find below.
    n_before = len(found)
    found += [(tag, [-x for x in p], [-y for y in D]) for tag, p, D in found]
    print("conjugate (10bar) branch added: %d -> %d" % (n_before, len(found)))
    rep["n_solutions_10_branch_only"] = n_before

    rep["n_branch2"], rep["n_branch1"] = len(br2), len(br1)
    rep["n_solutions"] = len(found)
    print("\nsolutions: branch 2 = %d, branch 1 = %d, total = %d"
          % (len(br2), len(br1), len(found)))

    # ---- sanity on every solution ----------------------------------------------------
    for tag, p, D in found:
        assert n2(p) == NORM
        assert proj_w_norm2(p) == W2
        assert in_lattice([6 * x for x in p])
        assert order_of(p) == 6
        assert n2(D) == 2 and in_lattice(D)
        assert all(dot(D, c) == 0 for c in SM_CARTAN)
        q = [a - b for a, b in zip(p, D)] if dot(p, D) == 1 else [a + b for a, b in
                                                                 zip(p, D)]
        assert n2(q) == NORM and proj_w_norm2(q) == W2
    print("all solutions pass: |p|^2 = 31/18, |P_W p|^2 = 6/5, 6p in Lambda, order 6,")
    print("Delta a root of Lambda perpendicular to W_SM, partner p -+ Delta on the shell")

    # ---- classes ---------------------------------------------------------------------
    byclass, bycoarse = {}, {}
    for tag, p, D in found:
        byclass.setdefault(class_key(p), (tag, p, D))
        bycoarse.setdefault(coarse_key(p), (tag, p))
    print("\ndistinct Lambda-classes (exact key): %d   (gate's coarser key: %d)"
          % (len(byclass), len(bycoarse)))
    rep["n_classes_exact"] = len(byclass)
    rep["n_classes_coarse"] = len(bycoarse)

    # ---- equivalence: G1 x W(E8)_2 x <-1> --------------------------------------------
    gens1 = [list(r) for r in PERP1] + [a[:8] for a in SM_ROOTS]
    gens2 = [list(r) for r in E8R]

    def orbits(keys, use_minus_one=True):
        parent = {k: k for k in keys}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def uni(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb

        ks = set(keys)
        for k in list(keys):
            p = byclass[k][1]
            for r in gens1:
                c = dot(p[:8], r)
                if c == 0:
                    continue
                q = [x - c * y for x, y in zip(p[:8], r)] + list(p[8:])
                kk = class_key(q)
                if kk in ks:
                    uni(k, kk)
            for r in gens2:
                c = dot(p[8:], r)
                if c == 0:
                    continue
                q = list(p[:8]) + [x - c * y for x, y in zip(p[8:], r)]
                kk = class_key(q)
                if kk in ks:
                    uni(k, kk)
            if use_minus_one:
                kk = class_key([-x for x in p])
                if kk in ks:
                    uni(k, kk)
        o = defaultdict(list)
        for k in keys:
            o[find(k)].append(k)
        return o

    # closure of the solution set under every generator (this is what makes the
    # union-find compute TRUE orbits, not merely connected components of a subgraph)
    ks = set(byclass)
    nclosed = 0
    for k in byclass:
        p = byclass[k][1]
        for r in [list(x) for x in PERP1] + [a[:8] for a in SM_ROOTS]:
            c = dot(p[:8], r)
            q = [x - c * y for x, y in zip(p[:8], r)] + list(p[8:])
            assert class_key(q) in ks, "solution set not closed under a block-1 generator"
            nclosed += 1
        for r in E8R:
            c = dot(p[8:], r)
            q = list(p[:8]) + [x - c * y for x, y in zip(p[8:], r)]
            assert class_key(q) in ks, "solution set not closed under a block-2 generator"
            nclosed += 1
        assert class_key([-x for x in p]) in ks
        nclosed += 1
    rep["closure_checks"] = nclosed
    print("solution set closed under all %d generator applications" % nclosed)

    orb = orbits(list(byclass), True)
    orb_no_m1 = orbits(list(byclass), False)
    print("inequivalent classes under G1 x W(E8)_2 x <-1> : %d" % len(orb))
    print("inequivalent classes under G1 x W(E8)_2        : %d" % len(orb_no_m1))
    rep["n_inequivalent"] = len(orb)
    rep["n_inequivalent_without_minus1"] = len(orb_no_m1)

    # ---- the distinguishing invariant ------------------------------------------------
    def min_norm_second_block(p):
        """min |y|^2 over y in p_B + E8, computed exactly (only 0 and 1/2 can occur here:
        6 p_B in E8 with |6 p_B|^2 in {0, 18}); computed by search over a small shell."""
        pB = list(p[8:])
        if all(x == 0 for x in pB):
            return F(0)
        best = n2(pB)
        for lam in [tuple([F(0)] * 8)] + [tuple(r) for r in E8R] + \
                   [tuple(-y for y in r) for r in E8R]:
            y = [a - b for a, b in zip(pB, lam)]
            if n2(y) < best:
                best = n2(y)
        return best

    reps = []
    for root, members in sorted(orb.items(), key=lambda kv: -len(kv[1])):
        tag, p, D = byclass[root]
        # prefer a member with the nicest printed form
        best = None
        for k in members:
            t, q, dd = byclass[k]
            score = (sum(abs(x.denominator) for x in q), [str(x) for x in q])
            if best is None or score < best[0]:
                best = (score, t, q, dd)
        _, tag, p, D = best
        d1 = [r for r in PERP1 if abs(dot(p[:8], r)) == 1]
        d2 = [r for r in E8R if abs(dot(p[8:], r)) == 1]
        reps.append(dict(
            branch=tag, n_classes_in_orbit=len(members),
            V2=[str(x) for x in p], Delta=[str(x) for x in D],
            order=order_of(p), norm2=str(n2(p)), w2=str(proj_w_norm2(p)),
            role=role_of_weight(p),
            MI2=str(6 * (n2(p) - F(1, 18))),
            n_Delta_in_SM_factor=len(d1), n_Delta_in_other_factor=len(d2),
            min_norm2_second_block=str(min_norm_second_block(p)),
            second_block_zero=all(x == 0 for x in p[8:])))
    rep["representatives"] = reps
    print("\n--- inequivalent classes ---")
    for r in reps:
        print("  branch %s | %d classes in orbit | order %d | MI2 %s | "
              "Delta in SM factor: %d, in other factor: %d | min|p_B|^2 = %s"
              % (r["branch"], r["n_classes_in_orbit"], r["order"], r["MI2"],
                 r["n_Delta_in_SM_factor"], r["n_Delta_in_other_factor"],
                 r["min_norm2_second_block"]))
        print("     V2    = (%s | %s)" % (", ".join(r["V2"][:8]), ", ".join(r["V2"][8:])))
        print("     Delta = (%s | %s)" % (", ".join(r["Delta"][:8]),
                                          ", ".join(r["Delta"][8:])))

    # ---- is the phase's target present, and in which orbit? --------------------------
    # ---- ORBIT CERTIFICATE: a breadth-first search from each representative, so that
    # every class carries an explicit word in the generators taking it to a
    # representative.  A referee can check any single row of this by hand.
    gen_list = ([("b1", i, list(r)) for i, r in enumerate(PERP1)]
                + [("b1", 100 + i, a[:8]) for i, a in enumerate(SM_ROOTS)]
                + [("b2", i, list(r)) for i, r in enumerate(E8R)]
                + [("m1", 0, None)])

    def apply_gen(p, g):
        kind, idx, r = g
        if kind == "b1":
            c = dot(p[:8], r)
            return [x - c * y for x, y in zip(p[:8], r)] + list(p[8:])
        if kind == "b2":
            c = dot(p[8:], r)
            return list(p[:8]) + [x - c * y for x, y in zip(p[8:], r)]
        return [-x for x in p]

    seeds = [sorted(orb.items(), key=lambda kv: -len(kv[1]))[i][0] for i in range(len(orb))]
    parent_word = {}
    frontier = list(seeds)
    for s in seeds:
        parent_word[s] = None
    while frontier:
        nxt = []
        for k in frontier:
            p = byclass[k][1]
            for g in gen_list:
                kk = class_key(apply_gen(p, g))
                if kk in byclass and kk not in parent_word:
                    parent_word[kk] = (k, g[0], g[1])
                    nxt.append(kk)
        frontier = nxt
    print("orbit certificate: %d of %d classes reached from the %d representatives"
          % (len(parent_word), len(byclass), len(seeds)))
    assert len(parent_word) == len(byclass)
    rep["orbit_certificate_classes_reached"] = len(parent_word)
    rep["orbit_certificate_generators"] = dict(
        block1_perp_roots=len(PERP1), block1_SM_simple_roots=len(SM_ROOTS),
        block2_roots=len(E8R), global_minus_one=1)

    # ---- named vectors: the paper-facing representatives and the phase/gate ones -----
    def which_orbit(v):
        k = class_key(v)
        if k not in byclass:
            return None
        for idx, (root, members) in enumerate(
                sorted(orb.items(), key=lambda kv: -len(kv[1]))):
            if k in members:
                return idx
        return None

    phase_p = [F(2, 3), F(2, 3), F(-1, 3), F(-1, 3), F(-1, 3), F(0), F(0), F(0),
               F(1, 2), F(1, 2)] + [F(0)] * 6
    gate_B = [F(7, 12), F(7, 12), F(-5, 12), F(-5, 12), F(-5, 12), F(-1, 12),
              F(-5, 12), F(7, 12)] + [F(0)] * 8
    nice_A = phase_p
    nice_B = [F(2, 3), F(2, 3), F(-1, 3), F(-1, 3), F(-1, 3), F(-1, 2), F(-1, 2),
              F(0)] + [F(0)] * 8
    named = {}
    for name, v in (("phase target (branch 2)", phase_p), ("gate case B", gate_B),
                    ("paper rep A", nice_A), ("paper rep B", nice_B)):
        w = which_orbit(v)
        named[name] = dict(present=class_key(v) in byclass, orbit=w,
                           norm2=str(n2(v)), w2=str(proj_w_norm2(v)),
                           order=order_of(v), role=role_of_weight(v),
                           MI2=str(6 * (n2(v) - F(1, 18))))
        print("\n%-24s |p|^2=%s  |P_W p|^2=%s  order=%d  role=%s  present=%s  orbit=%s"
              % (name, n2(v), proj_w_norm2(v), order_of(v), role_of_weight(v),
                 class_key(v) in byclass, w))
    rep["named_vectors"] = named

    # ---- BREADTH-FIRST CERTIFICATE, WRITTEN OUT (added 2026-10-22, package v2.1 amendment).
    # The search above asserts that every class is reached; this block stores the rows so
    # that the assertion can be checked without re-running anything.  It is seeded at the
    # two representatives the paper prints (eq. (reps) / table of representatives), and each
    # row records one class, its representative momentum, the class it was reached from
    # (its parent) and the generator g with  g(parent representative) in  class + Lambda.
    # Every generator is an involution, so g also carries the class back to its parent, and
    # following parents to a seed composes an automorphism carrying the class to a printed
    # representative.  Writes paperA_two_classes_bfs_certificate.json; the main JSON above
    # is not touched.  check_bfs_certificate.py verifies every row independently.
    cert_seeds = [class_key(nice_A), class_key(nice_B)]
    assert all(s in byclass for s in cert_seeds)
    cidx = {}
    rows = []
    frontier = []
    for si, (s, v) in enumerate(zip(cert_seeds, (nice_A, nice_B))):
        cidx[s] = len(rows)
        rows.append(dict(index=len(rows), momentum=[str(x) for x in v], parent=None,
                         generator=None, depth=0, seed=si))
        frontier.append((s, v))
    while frontier:
        nxt = []
        for k, pv in frontier:
            for g in gen_list:
                q = apply_gen(pv, g)
                kk = class_key(q)
                if kk in byclass and kk not in cidx:
                    cidx[kk] = len(rows)
                    gdesc = dict(kind={"b1": "reflection, first E8",
                                       "b2": "reflection, second E8",
                                       "m1": "global -1"}[g[0]],
                                 root=None if g[2] is None else [str(x) for x in g[2]])
                    rows.append(dict(index=len(rows), momentum=[str(x) for x in q],
                                     parent=cidx[k], generator=gdesc,
                                     depth=rows[cidx[k]]["depth"] + 1,
                                     seed=rows[cidx[k]]["seed"]))
                    nxt.append((kk, q))
        frontier = nxt
    assert len(rows) == len(byclass), (len(rows), len(byclass))
    cert = dict(
        note=("Breadth-first orbit certificate for the two-class theorem. Row i is one "
              "lattice class, represented by 'momentum'; 'parent' is the row it was reached "
              "from and 'generator' the involution g with g(momentum of parent) - momentum "
              "in Lambda = E8+E8. A 'reflection' with root r acts as x -> x - (x.r) r on the "
              "named E8 block; 'global -1' is x -> -x. Rows 0 and 1 are the two printed "
              "representatives V2^A and V2^B."),
        n_classes=len(rows),
        n_seeds=2,
        max_depth=max(r["depth"] for r in rows),
        rows_per_seed=[sum(1 for r in rows if r["seed"] == i) for i in (0, 1)],
        rows=rows)
    json.dump(cert, open(os.path.join(HERE, "paperA_two_classes_bfs_certificate.json"), "w"),
              indent=0)
    print("\nbreadth-first certificate written: %d classes, %d from V2^A and %d from V2^B, "
          "maximum depth %d" % (len(rows), cert["rows_per_seed"][0], cert["rows_per_seed"][1],
                                cert["max_depth"]))

    # ---- reconciliation with the gate's published 3,000 / 119 ------------------------
    # the gate's branch-A loop imposed rem == 1/2 (only Bv = -2) where the correct
    # condition is rem >= 1/2 (Bv in {-5/2, -2}); and it counted case B once per p
    # while counting case A once per (p, Delta).  Restricting to the gate's own subset
    # and applying the gate's own coarser class key must reproduce 3,000 and 119.
    subset = [p for tag, p, D in found
              if (tag == "2" and all(p[k] == F(0) for k in (5, 6, 7))) or tag == "1"]
    gate_coarse = {coarse_key(p) for p in subset}
    gate_exact = {class_key(p) for p in subset}
    n_gate_style = sum(1 for tag, p, D in found
                       if tag == "2" and all(p[k] == F(0) for k in (5, 6, 7))) \
        + len({tuple(str(x) for x in p) for tag, p, D in found if tag == "1"})
    print("\nreconciliation with the gate:")
    print("   gate-style solution count on the gate's own subset : %d  (gate: 3000)"
          % n_gate_style)
    print("   coarse-key classes on that subset                  : %d  (gate: 119)"
          % len(gate_coarse))
    print("   exact-key classes on that same subset              : %d" % len(gate_exact))
    rep["reconciliation"] = dict(gate_style_solutions=n_gate_style,
                                 coarse_classes_on_gate_subset=len(gate_coarse),
                                 exact_classes_on_gate_subset=len(gate_exact))

    # ---- counts in one convention ----------------------------------------------------
    distinct_p = {tuple(str(x) for x in p) for tag, p, D in found}
    n_pD = 0
    for q in distinct_p:
        pv = [F(x) for x in q]
        n_pD += sum(1 for r in PERP1 if abs(dot(pv[:8], r)) == 1) \
            + sum(1 for r in E8R if abs(dot(pv[8:], r)) == 1)
    print("\nin ONE convention: distinct p = %d ; (p, Delta) incidences = %d"
          % (len(distinct_p), n_pD))
    rep["n_distinct_p"] = len(distinct_p)
    rep["n_p_Delta_incidences"] = n_pD

    json.dump(rep, open(os.path.join(HERE, "paperA_two_classes_enum.json"), "w"),
              indent=1)
    print("\nwrote paperA_two_classes_enum.json")


if __name__ == "__main__":
    main()

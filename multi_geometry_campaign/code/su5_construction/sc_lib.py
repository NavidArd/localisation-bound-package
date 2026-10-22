"""su5_construction -- core library for the CB-1 global-embedding construction phase.

Spec: docs/su5_loophole_construction_task.md.  Geometry Z3xZ6_1_1 ONLY.

======================================================================================
CONVENTIONS, STATED ONCE (CLAUDE.md 5.19 (i))
======================================================================================
C1  Lattice.  Lambda = Lambda_E8 x Lambda_E8, even self-dual, EVEN COORDINATE SYSTEM:
    per block, x in Z^8 or x in (Z+1/2)^8, with sum(x) in 2Z.  Roots have norm^2 = 2.
    16-vectors are ordered (first E8 | second E8) exactly as orbifolder prints p_sh and
    exactly as the parr2020 catalogue rows are written.

C2  PLANES.  Read off the PRIMARY geometry file orbifolder-1.2.1/Geometry/
    Geometry_Z3xZ6_1_1.txt, `begin twist` (slot 0 is the 4d slot and is dropped):
        v1 = ( 0 , 1/3 , -1/3 )      order 3   (theta, the "M = 3" generator)
        v2 = ( 1/6,  0  , -1/6 )     order 6   (omega,  the "N = 6" generator)
    Hence: plane 1 is unrotated by THETA, plane 2 is unrotated by OMEGA, plane 3 is
    rotated by both.  Plane orders (6, 3, 6).
    THIS ORDERING IS LOAD-BEARING HERE and is NOT interchangeable with the plane-swapped
    version used by su5_ten_bound/measure.py: delta_c is permutation symmetric but the
    WILSON-LINE PLANE IS NOT.  The single independent Wilson line of this geometry is
    W_3 = W_4, order 3, in PLANE 2 (z3xz6_foundations_findings.md section 4, derived by
    Smith normal form; the geometry file's `allowed order of Wilson lines` block reads
    1 1 3 3 1 1 and `begin identical WLs` reads `3 4`).  Plane 2 is the plane OMEGA does
    not rotate --> the T(0,l) sectors are the N=2 sectors of plane 2.

C3  Catalogue rows.  parr2020 `Shifts and Wilsonlines:` block = 8 rows of 16 rationals:
    row 0 = V1 (order | 3), row 1 = V2 (order | 6), rows 2..7 = W_1..W_6.
    Verified by order measurement in s0_invariants.py, not assumed.

C4  Local shift.  V_g(k, l, n) = k V1 + l V2 + sum_a n_a W_a.

C5  MODULAR INVARIANCE -- the complete orbifolder-1.2.1 condition set, transcribed at
    SOURCE level in THIS phase (5.16 claim-in-source rule) from
      src/orbifolder/corbifoldgroup.cpp :: COrbifoldGroup::CheckModularInvariance, l.521-735
      src/orbifolder/cspacegroup.cpp l.1990-2016  (SP_factors: Orders = [M, N, w1..w6];
          SP[i][j] = Orders[i] if Orders[i]==Orders[j] else gcd(Orders[i],Orders[j]))
    For Z3xZ6_1_1  Orders = [3, 6, 1, 1, 3, 3, 1, 1]:
      MI0  ord(V1) divides 3 ; ord(V2) divides 6
      MI1  3 (V1.V1 - v1.v1)  == 0 mod 2       [prefactor M     = 3]
      MI2  6 (V2.V2 - v2.v2)  == 0 mod 2       [prefactor N     = 6]
      MI3  3 (V1.V2 - v1.v2)  == 0 mod 2       [SP[0][1] = gcd(3,6) = 3]
      MI4  3 (W_i.V1)         == 0 mod 2       [SP[i+2][0], i in {3,4} -> 3]
      MI5  3 (W_i.V2)         == 0 mod 2       [SP[i+2][1] = gcd(3,6) = 3]
      MI6  3 (W_i.W_j)        == 0 mod 2       [SP[i+2][j+2] = 3]
    (the source applies MI4/MI5/MI6 only to NON-ZERO Wilson lines, and MI6 for j >= i.)

C6  Records vs states (G-9).  A spectrum RECORD that is a non-singlet of a hidden factor
    is several massless STATES sharing the SM role and the cell.  Every count in this
    phase says which it is.  Z3xZ6's certified role census (z36sw_lib.role_census) sets
    role = None on hidden non-singlets; Z2xZ4's does not.  Both conventions are exercised
    where both geometries are used, and labelled.
"""

from __future__ import annotations

# --- reproduction-package root shim -------------------------------------------------
# The line above originally read  F:\physics  -- the absolute path of the author's
# working tree.  This is the ONLY edit made to this file in the reproduction package;
# it points the same variable at the package root, which carries the same directory
# names.  Set PAPERA_PKG to relocate the package.  See
# environment/PATCHES_TO_PACKAGE_CODE.md.
# -------------------------------------------------------------------------------------
import os as _os
_PKG_ROOT = _os.environ.get('PAPERA_PKG') or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..'))

import os
import re
from fractions import Fraction as F
from math import gcd

REPO = _PKG_ROOT
CAMP = os.path.join(REPO, "multi_geometry_campaign")
OUT = os.path.join(CAMP, "data", "su5_construction")
CAT_Z36 = os.path.join(CAMP, "data", "literature", "parr2020_full", "extracted",
                       "Z3xZ6_1_1.txt")
GEOM_FILE = os.path.join(REPO, "orbifolder-1.2.1", "Geometry", "Geometry_Z3xZ6_1_1.txt")

M_ORD, N_ORD = 3, 6
V1_TWIST = (F(0), F(1, 3), F(-1, 3))          # theta   (C2)
V2_TWIST = (F(1, 6), F(0), F(-1, 6))          # omega   (C2)
WL_ORDERS = (1, 1, 3, 3, 1, 1)                # geometry file, cross-checked in s0
SP = None                                     # filled by _make_sp()


def _make_sp():
    global SP
    orders = [M_ORD, N_ORD] + list(WL_ORDERS)
    SP = [[orders[i] if orders[i] == orders[j] else gcd(orders[i], orders[j])
           for j in range(8)] for i in range(8)]


_make_sp()


# ---------------------------------------------------------------- basic lattice
def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def in_e8(v):
    """even coordinate system membership test for one 8-block."""
    if all(x.denominator == 1 for x in v):
        return sum(v) % 2 == 0
    if all((x * 2).denominator == 1 and (x * 2).numerator % 2 == 1 for x in v):
        return sum(v) % 2 == 0
    return False


def in_lattice(v):
    return in_e8(v[:8]) and in_e8(v[8:])


def order_of(v, cap=64):
    for n in range(1, cap + 1):
        if in_lattice([n * x for x in v]):
            return n
    return 0


def reduce_mod_lattice(v):
    """canonical representative of v + Lambda with entries in [0,1) except that the
    half-integer class is respected; used only for hashing a shift CLASS."""
    out = []
    for blk in (v[:8], v[8:]):
        r = [x - (x.numerator // x.denominator) for x in blk]      # in [0,1)
        # the class is r + E8; E8 has index-1 in Z^8-union-(Z+1/2)^8 so we must also
        # allow subtracting the half-integer glue vector when that lands in E8
        out.append(r)
    a, b = out
    return tuple(a) + tuple(b)


def class_key(v):
    """A hashable key for v mod Lambda.  Two vectors have equal key iff they differ by a
    lattice vector.  Implemented by reducing coordinates mod 1 and then mod the E8 glue
    (which for E8 in the even coordinate system means: (r) and (r + (1/2)^8) are the same
    class iff (1/2)^8-shift lands back in E8 -- it does not, since (1/2)^8 has sum 4 which
    IS even, so (1/2)^8 IS in E8 and the two ARE the same class).  Hence per block we
    reduce mod 1 and then take the lexicographically smaller of r and r+(1/2)^8 mod 1."""
    key = []
    for blk in (v[:8], v[8:]):
        r = tuple(x - F(x.numerator // x.denominator) for x in blk)
        r2 = tuple((x + F(1, 2)) - F((x + F(1, 2)).numerator // (x + F(1, 2)).denominator)
                   for x in r)
        # r2 is in the same Lambda-class as r iff (1/2)^8 in E8 : sum = 4, even -> YES
        key.append(min(r, r2))
    return tuple(key[0]) + tuple(key[1])


_E8R = None
_E8xE8R = None


def e8_roots():
    global _E8R
    if _E8R is None:
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
        _E8R = R
    return _E8R


def e8xe8_roots():
    global _E8xE8R
    if _E8xE8R is None:
        z = (F(0),) * 8
        r8 = e8_roots()
        _E8xE8R = [tuple(r) + z for r in r8] + [z + tuple(r) for r in r8]
    return _E8xE8R


# ---------------------------------------------------------------- modular invariance
def mi_report(V1, V2, Ws):
    """The COMPLETE orbifolder condition set of C5.  Returns (ok, dict of failures)."""
    bad = {}
    o1, o2 = order_of(V1), order_of(V2)
    if o1 == 0 or M_ORD % o1 != 0:
        bad["MI0_V1"] = o1
    if o2 == 0 or N_ORD % o2 != 0:
        bad["MI0_V2"] = o2
    v11 = dot(V1_TWIST, V1_TWIST)
    v22 = dot(V2_TWIST, V2_TWIST)
    v12 = dot(V1_TWIST, V2_TWIST)
    s = M_ORD * (dot(V1, V1) - v11)
    if s.denominator != 1 or s % 2 != 0:
        bad["MI1"] = str(s)
    s = N_ORD * (dot(V2, V2) - v22)
    if s.denominator != 1 or s % 2 != 0:
        bad["MI2"] = str(s)
    s = SP[0][1] * (dot(V1, V2) - v12)
    if s.denominator != 1 or s % 2 != 0:
        bad["MI3"] = str(s)
    for i, W in enumerate(Ws):
        if all(x == 0 for x in W):
            continue
        if WL_ORDERS[i] and order_of(W) and WL_ORDERS[i] % order_of(W) != 0:
            bad["MI0_W%d" % (i + 1)] = order_of(W)
        s = SP[i + 2][0] * dot(W, V1)
        if s.denominator != 1 or s % 2 != 0:
            bad["MI4_W%d" % (i + 1)] = str(s)
        s = SP[i + 2][1] * dot(W, V2)
        if s.denominator != 1 or s % 2 != 0:
            bad["MI5_W%d" % (i + 1)] = str(s)
        for j in range(i, len(Ws)):
            Wj = Ws[j]
            if all(x == 0 for x in Wj):
                continue
            s = SP[i + 2][j + 2] * dot(W, Wj)
            if s.denominator != 1 or s % 2 != 0:
                bad["MI6_W%dW%d" % (i + 1, j + 1)] = str(s)
    return (not bad), bad


# ---------------------------------------------------------------- catalogue
def rat(t):
    t = t.strip()
    if "/" in t:
        a, b = t.split("/")
        return F(int(a), int(b))
    return F(int(t))


def read_catalogue(path=CAT_Z36):
    """yields (index, label, V1, V2, [W1..W6]) for every `begin model` block."""
    idx = -1
    rows, lab, reading, sg = None, None, False, None
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            s = raw.strip()
            if s == "begin model":
                idx += 1
                rows, lab, reading, sg = [], None, False, None
                continue
            if s == "end model":
                yield idx, lab, sg, rows
                continue
            if s.startswith("Label:"):
                lab = s[6:].strip()
            elif s.startswith("SpaceGroup:"):
                sg = s[11:].strip()
            elif s.startswith("Shifts and Wilsonlines"):
                reading = True
            elif reading and s and ":" not in s:
                f = [x.strip() for x in s.split(",")]
                if len(f) == 16:
                    rows.append([rat(x) for x in f])


# ---------------------------------------------------------------- shells
def shell(Vg, norm2, cap=None):
    """Complete list of p in Vg + Lambda with |p|^2 == norm2 (a Fraction).

    Enumerated block by block: for each split norm2 = n1 + n2 with n1 realised in
    block 1's coset and n2 in block 2's, take the product.  Per-block enumeration is a
    direct bounded search over the coset (Vg_blk + E8) using the even-coordinate
    description.  Self-tested in selftest() against the E8 theta series.
    """
    a = _block_shell(Vg[:8], norm2)
    b = _block_shell(Vg[8:], norm2)
    out = []
    for n1, L1 in a.items():
        n2 = norm2 - n1
        if n2 in b:
            for x in L1:
                for y in b[n2]:
                    out.append(tuple(x) + tuple(y))
                    if cap and len(out) > cap:
                        return out
    return out


_BLK_CACHE = {}


def _block_shell(V, maxnorm):
    """dict {norm2: [vectors]} for all v in V + E8 with |v|^2 <= maxnorm."""
    key = (tuple(V), maxnorm)
    if key in _BLK_CACHE:
        return _BLK_CACHE[key]
    res = {}
    # V + E8 = union of (V + Z^8 with sum in 2Z shifted) and (V + (Z+1/2)^8 ...)
    # enumerate over the two E8 cosets of D8: lambda in Z^8 with sum even, or
    # lambda in (Z+1/2)^8 with sum even.
    for half in (False, True):
        base = [x + (F(1, 2) if half else F(0)) for x in V]
        # each coordinate is base_i + integer; find integer ranges
        rng = []
        ok = True
        import math
        b = math.isqrt(int(maxnorm) + 2) + 3
        lim = F(maxnorm)
        for x in base:
            vals = [(x + n, n) for n in range(-b, b + 1) if (x + n) ** 2 <= lim]
            if not vals:
                ok = False
                break
            rng.append(vals)
        if not ok:
            continue
        # v = V (+ s if half) + n,  n in Z^8 with sum(n) even   [E8 = D8 u (D8 + s)]
        _rec(rng, 0, [], F(0), maxnorm, res, 0)
    _BLK_CACHE[key] = res
    return res


def _rec(rng, i, cur, acc, maxnorm, res, nsum):
    if acc > maxnorm:
        return
    if i == len(rng):
        if nsum % 2 == 0:
            res.setdefault(acc, []).append(tuple(cur))
        return
    for x, n in rng[i]:
        a2 = acc + x * x
        if a2 > maxnorm:
            continue
        cur.append(x)
        _rec(rng, i + 1, cur, a2, maxnorm, res, nsum + n)
        cur.pop()


# ---------------------------------------------------------------- FAST integer shell
def shell_int(Vg, num, den, scale=12):
    """Same as shell(), in pure integer arithmetic.

    Vg is a 16-vector of Fractions; the shell is {p in Vg + Lambda : |p|^2 = num/den}.
    Work with X = scale * p (integers).  Requires scale * Vg integral and
    scale^2 * num/den integral; both asserted.  Returns a list of 16-tuples of Fractions.
    """
    A = [x * scale for x in Vg]
    assert all(a.denominator == 1 for a in A), "scale too small for Vg"
    A = [int(a) for a in A]
    T = F(num, den) * scale * scale
    assert T.denominator == 1, "scale too small for the norm"
    T = int(T)
    b1 = _blk_int(A[:8], T, scale)
    b2 = _blk_int(A[8:], T, scale)
    out = []
    for n1, L1 in b1.items():
        n2 = T - n1
        L2 = b2.get(n2)
        if not L2:
            continue
        for x in L1:
            for y in L2:
                out.append(tuple(F(c, scale) for c in x) + tuple(F(c, scale) for c in y))
    return out


def _blk_int(A, T, scale):
    """{norm2*scale^2 : [integer 8-vectors]} for X in A + scale*E8 with |X|^2 <= T.

    scale*E8 = scale*D8  u  scale*(D8 + s), s = (1/2)^8, so the offsets are
    scale*n (sum n even) and scale*n + scale/2 * (1,..,1) (sum n even).  scale must be
    even for the second branch to be integral (12 is).
    """
    res = {}
    half = scale // 2
    for br in (0, 1):
        base = [a + (half if br else 0) for a in A]
        rng = []
        ok = True
        for x in base:
            vals = [(x + scale * n, n) for n in range(-4, 5)
                    if (x + scale * n) ** 2 <= T]
            if not vals:
                ok = False
                break
            rng.append(vals)
        if not ok:
            continue
        # remaining minimal squares for pruning
        mins = [0] * 9
        for i in range(7, -1, -1):
            mins[i] = mins[i + 1] + min(v * v for v, _ in rng[i])
        _rec_int(rng, mins, 0, [], 0, T, res, 0)
    return res


def _rec_int(rng, mins, i, cur, acc, T, res, nsum):
    if acc + mins[i] > T:
        return
    if i == 8:
        if nsum % 2 == 0:
            res.setdefault(acc, []).append(tuple(cur))
        return
    for x, n in rng[i]:
        a2 = acc + x * x
        if a2 + mins[i + 1] > T:
            continue
        cur.append(x)
        _rec_int(rng, mins, i + 1, cur, a2, T, res, nsum + n)
        cur.pop()


# ---------------------------------------------------------------- self-tests
def selftest(verbose=True):
    """5.19(i): one known invariant verified BEFORE any count is believed."""
    ok = True
    z = (F(0),) * 8
    # E8 theta coefficients 1, 240, 2160, 6720
    bs = _block_shell(z, F(6))
    got = [len(bs.get(F(n), [])) for n in (0, 2, 4, 6)]
    exp = [1, 240, 2160, 6720]
    ok &= got == exp
    if verbose:
        print("E8 theta 0,2,4,6 :", got, "expected", exp, "OK" if got == exp else "FAIL")
    # E8xE8 : 1, 480, 61920
    g2 = [len(shell((F(0),) * 16, F(n))) for n in (0, 2, 4)]
    e2 = [1, 480, 61920]
    ok &= g2 == e2
    if verbose:
        print("E8xE8 theta 0,2,4 :", g2, "expected", e2, "OK" if g2 == e2 else "FAIL")
    # coset representative must lie in its own shell
    p = [F(2, 3), F(2, 3), F(-1, 3), F(-1, 3), F(-1, 3), F(0), F(0), F(0),
         F(1, 2), F(1, 2), F(0), F(0), F(0), F(0), F(0), F(0)]
    sh = shell(p, dot(p, p))
    ok &= tuple(p) in set(sh)
    if verbose:
        print("coset rep in own shell:", tuple(p) in set(sh), " shell size", len(sh))
    # the FAST integer enumerator must agree with the exact-Fraction one, and reproduce
    # the E8xE8 theta coefficients too
    g3 = [len(shell_int((F(0),) * 16, n, 1)) for n in (0, 2, 4)]
    ok &= g3 == e2
    if verbose:
        print("E8xE8 theta via shell_int:", g3, "OK" if g3 == e2 else "FAIL")
    s2 = set(shell_int(p, 31, 18))
    ok &= s2 == set(sh)
    if verbose:
        print("shell_int == shell on the target coset:", s2 == set(sh), len(s2))
    # 480 roots, all norm 2, all in lattice
    R = e8xe8_roots()
    ok &= len(R) == 480 and all(dot(r, r) == 2 and in_lattice(r) for r in R)
    if verbose:
        print("roots:", len(R), "OK" if len(R) == 480 else "FAIL")
    return ok


if __name__ == "__main__":
    print("selftest:", "PASS" if selftest() else "FAIL")

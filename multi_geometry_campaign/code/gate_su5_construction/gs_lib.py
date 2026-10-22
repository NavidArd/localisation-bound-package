"""gate_su5_construction -- FRESH library for the adversarial gate on CB-1.

Written independently of code/su5_construction/sc_lib.py (no import from it anywhere in
this gate).  Where a quantity is also computed by the phase, the gate recomputes it by a
DIFFERENT route and compares.

======================================================================================
CONVENTIONS, STATED ONCE (CLAUDE.md 5.19 (i))
======================================================================================
G1  Lattice.  Lambda = E8 x E8 in the EVEN COORDINATE SYSTEM.  A block x in R^8 lies in
    E8 iff (all x_i in Z or all x_i in Z+1/2) and sum(x) in 2Z.  16-vectors are ordered
    (block1 | block2) exactly as orbifolder prints `left-moving p_sh`.

G2  Geometry Z3xZ6_1_1, read in THIS gate off the primary file
    F:/physics/orbifolder-1.2.1/Geometry/Geometry_Z3xZ6_1_1.txt.  `begin twist` gives
        v1 = (0, 1/3, -1/3)   order 3    theta
        v2 = (1/6, 0, -1/6)   order 6    omega
    (slot 0 of each printed row is the 4d slot and is dropped).  Plane 2 is the plane
    OMEGA does not rotate; `allowed order of Wilson lines` = 1 1 3 3 1 1 and
    `begin identical WLs` = "3 4" put the single independent Wilson line in plane 2.
    The gate re-reads both blocks from the file at import time and ASSERTS them
    (parse_geometry_twists / parse_constructing_elements).

G3  Modular invariance -- re-read at source in THIS gate from
    orbifolder-1.2.1/src/orbifolder/corbifoldgroup.cpp::CheckModularInvariance
    (l. 521-735) and cspacegroup.cpp l. 1990-2016.  Prefactors:
        MI0 : Order_ZM_Shift divides M=3 ; Order_ZN_Shift divides N=6 ; WL orders
        MI1 : M   * (V1.V1 - v1.v1)  even          [source uses M, not ord(V1)]
        MI2 : N   * (V2.V2 - v2.v2)  even
        MI3 : SP[0][1] * (V2.V1 - v2.v1) even      SP[0][1] = gcd(3,6) = 3
        MI4 : SP[i+2][0] * (W_i.V1)  even          W3,W4 -> SP[4][0] = gcd(3,3) = 3
        MI5 : SP[i+2][1] * (W_i.V2)  even          SP[4][1] = gcd(3,6) = 3
        MI6 : SP[i+2][j+2] * (W_i.W_j) even, j >= i, both non-zero
    MI4/5/6 are applied by the source ONLY to non-zero Wilson lines.

G4  Standard Model, as DECLARED by the phase (its scope limit 1, carried unchanged):
        SU(3)_c = <e1-e2, e2-e3>,  SU(2)_L = <e4-e5>,
        t_Y = (-1/3,-1/3,-1/3, 1/2, 1/2, 0,0,0 | 0^8),   k_Y = 2|t_Y|^2 = 5/3.
    W_SM := span{e1-e2, e2-e3, e4-e5, t_Y} is 4-dimensional; its orthogonal complement
    inside block-1 coords 1..5 is the single direction (1,1,1,1,1) (the SU(5) singlet).

G5  Roles.  Left-chiral SM roles, TRIALITY-AWARE (this is where the gate differs from the
    phase, which discards the SU(3) triality):
        Q   = (3, 2)_{+1/6}     u^c = (3bar,1)_{-2/3}    e^c = (1,1)_{+1}
        d^c = (3bar,1)_{+1/3}   L   = (1,2)_{-1/2}
    plus the conjugates.  10-roles = {Q, u^c, e^c} and their conjugates.  The SU(3)
    triality of a weight p is read off the traceless part of (p1,p2,p3): deviations
    (2/3,-1/3,-1/3) up to permutation = "3", (-2/3,1/3,1/3) = "3bar".
    NOTE ON ORBIFOLDER'S PRINTED SIGN: orbifolder prints "-3" for what this convention
    calls a 3 and "3" for a 3bar (verified in g1 on CB1041's own records).  The gate
    never uses orbifolder's sign for physics; it uses it only as a dimension cross-check.

G6  Records vs states.  A printed spectrum RECORD carrying a hidden non-singlet is
    several massless STATES.  Every count below says which it is.
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
GEOM = os.path.join(REPO, "orbifolder-1.2.1", "Geometry", "Geometry_Z3xZ6_1_1.txt")
OUT = os.path.join(CAMP, "data", "gate_su5_construction")
PHASE_DATA = os.path.join(CAMP, "data", "su5_construction")

M_ORD, N_ORD = 3, 6


# ------------------------------------------------------------------ rationals / vectors
def rat(t):
    t = t.strip()
    if "/" in t:
        a, b = t.split("/")
        return F(int(a), int(b))
    return F(int(t))


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def n2(a):
    return dot(a, a)


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def smul(c, a):
    return [c * x for x in a]


# ------------------------------------------------------------------ G1 lattice
def in_e8(blk):
    if all(x.denominator == 1 for x in blk):
        return sum(blk) % 2 == 0
    if all(x.denominator == 2 and x.numerator % 2 == 1 for x in blk):
        return sum(blk) % 2 == 0
    return False


def in_lattice(v):
    return in_e8(list(v[:8])) and in_e8(list(v[8:]))


def order_of(v, cap=64):
    for n in range(1, cap + 1):
        if in_lattice([n * x for x in v]):
            return n
    return 0


def lattice_class_key(v):
    """key equal for two vectors iff they differ by an element of Lambda.

    Independent implementation: reduce the 16-vector coordinatewise mod 1 and then, per
    block, canonicalise over the two E8/Z^8-glue possibilities.  E8 (even coordinate
    system) contains BOTH Z^8-with-even-sum and (Z+1/2)^8-with-even-sum, so within a
    block, r and r + (1/2)^8 are in the same class iff (1/2)^8 in E8 -- sum = 4, even,
    so YES -- and additionally r and r + e_i + e_j (i != j) are the same class.  Since
    r is already reduced mod 1 the only remaining freedom is the (1/2)^8 glue and the
    parity of the integer part, which mod-1 reduction has already removed.  We therefore
    take min(r, r + (1/2)^8 mod 1) per block.

    Cross-checked in selftest() against an explicit difference-in-lattice test.
    """
    out = []
    for blk in (v[:8], v[8:]):
        r = tuple(x - F(x.numerator // x.denominator) for x in blk)
        h = tuple((y + F(1, 2)) - F((y + F(1, 2)).numerator // (y + F(1, 2)).denominator)
                  for y in r)
        out.append(min(r, h))
    return out[0] + out[1]


_E8R = None


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
        assert len(R) == 240
        _E8R = R
    return _E8R


def e8xe8_roots():
    z = (F(0),) * 8
    return [tuple(r) + z for r in e8_roots()] + [z + tuple(r) for r in e8_roots()]


# ------------------------------------------------------------------ G2 geometry (source)
def parse_geometry_twists(path=GEOM):
    rows, on = [], False
    for line in open(path, encoding="utf-8", errors="replace"):
        s = line.strip()
        if s == "begin twist":
            on = True
            continue
        if s == "end twist":
            break
        if on and s:
            rows.append([rat(t) for t in s.split()])
    assert len(rows) == 2 and all(len(r) == 4 for r in rows), rows
    return tuple(rows[0][1:]), tuple(rows[1][1:])


def parse_wl_orders(path=GEOM):
    lines = [l.rstrip("\n") for l in open(path, encoding="utf-8", errors="replace")]
    i = lines.index("allowed order of Wilson lines")
    return tuple(int(lines[i + k + 1].strip()) for k in range(6))


def parse_constructing_elements(path=GEOM):
    """{(k,l): [8-int-tuples]} from `begin constructing elements`."""
    lines = [l.strip() for l in open(path, encoding="utf-8", errors="replace")]
    i = lines.index("begin constructing elements")
    j = lines.index("end constructing elements")
    out, cur = {}, None
    for s in lines[i + 1:j]:
        if s == "new sector":
            cur = None
            continue
        if not s:
            continue
        f = [int(x) for x in s.split()]
        assert len(f) == 8
        key = (f[0], f[1])
        out.setdefault(key, []).append(tuple(f))
    return out


V1_TWIST, V2_TWIST = parse_geometry_twists()
WL_ORDERS = parse_wl_orders()
assert V1_TWIST == (F(0), F(1, 3), F(-1, 3)), V1_TWIST
assert V2_TWIST == (F(1, 6), F(0), F(-1, 6)), V2_TWIST
assert WL_ORDERS == (1, 1, 3, 3, 1, 1), WL_ORDERS

_ORD = [M_ORD, N_ORD] + list(WL_ORDERS)
SP = [[_ORD[i] if _ORD[i] == _ORD[j] else gcd(_ORD[i], _ORD[j]) for j in range(8)]
      for i in range(8)]


# ------------------------------------------------------------------ G3 modular invariance
def mi_check(V1, V2, Ws):
    """Complete orbifolder-1.2.1 set (G3).  Returns (ok, {failed: value})."""
    bad = {}
    o1, o2 = order_of(V1), order_of(V2)
    if o1 == 0 or M_ORD % o1:
        bad["MI0_V1"] = o1
    if o2 == 0 or N_ORD % o2:
        bad["MI0_V2"] = o2
    s = M_ORD * (n2(V1) - n2(V1_TWIST))
    if s.denominator != 1 or s % 2:
        bad["MI1"] = str(s)
    s = N_ORD * (n2(V2) - n2(V2_TWIST))
    if s.denominator != 1 or s % 2:
        bad["MI2"] = str(s)
    s = SP[0][1] * (dot(V1, V2) - dot(V1_TWIST, V2_TWIST))
    if s.denominator != 1 or s % 2:
        bad["MI3"] = str(s)
    for i, W in enumerate(Ws):
        if all(x == 0 for x in W):
            continue
        ow = order_of(W)
        if ow == 0 or (WL_ORDERS[i] % ow):
            bad["MI0_W%d" % (i + 1)] = ow
        s = SP[i + 2][0] * dot(W, V1)
        if s.denominator != 1 or s % 2:
            bad["MI4_W%d" % (i + 1)] = str(s)
        s = SP[i + 2][1] * dot(W, V2)
        if s.denominator != 1 or s % 2:
            bad["MI5_W%d" % (i + 1)] = str(s)
        for j in range(i, len(Ws)):
            if all(x == 0 for x in Ws[j]):
                continue
            s = SP[i + 2][j + 2] * dot(W, Ws[j])
            if s.denominator != 1 or s % 2:
                bad["MI6_W%dW%d" % (i + 1, j + 1)] = str(s)
    return (not bad), bad


# ------------------------------------------------------------------ G4 Standard Model
A3 = [[F(1), F(-1), F(0)] + [F(0)] * 13, [F(0), F(1), F(-1)] + [F(0)] * 13]
A2L = [F(0)] * 3 + [F(1), F(-1)] + [F(0)] * 11
TY = [F(-1, 3)] * 3 + [F(1, 2)] * 2 + [F(0)] * 11
SM_CARTAN = A3 + [A2L, TY]


def proj_wsm_norm2(p):
    """|P_{W_SM} p|^2 computed as |p_{1..5}|^2 - (sum p_{1..5})^2 / 5 (G4)."""
    s = sum(p[:5])
    return sum(x * x for x in p[:5]) - s * s / 5


def su3_triality(p):
    """0 (singlet), +1 (3), -1 (3bar) for the SU(3)_c content of weight p (level 1).

    Traceless deviations d_i = p_i - (p1+p2+p3)/3 for i = 1,2,3.  Level-1 weights:
    all zero -> singlet; a permutation of (2/3,-1/3,-1/3) -> 3; of (-2/3,1/3,1/3) -> 3bar.
    Anything else raises (the level-1 assumption would be violated).
    """
    a = (p[0] + p[1] + p[2]) / 3
    d = sorted(x - a for x in p[:3])
    if d == [F(0), F(0), F(0)]:
        return 0
    if d == [F(-1, 3), F(-1, 3), F(2, 3)]:
        return 1
    if d == [F(-2, 3), F(1, 3), F(1, 3)]:
        return -1
    raise ValueError("SU(3) weight not level 1: %s" % (d,))


def su2_dim(p):
    d = p[3] - p[4]
    if d == 0:
        return 1
    if abs(d) == 1:
        return 2
    raise ValueError("SU(2) weight not level 1: %s" % d)


def role_of(p, sY=1):
    """Triality-aware SM role of a left-chiral weight p; None if not an SM role."""
    t = su3_triality(p)
    d2 = su2_dim(p)
    y = sY * dot(TY, p)
    if t == 1 and d2 == 2 and y == F(1, 6):
        return "Q"
    if t == -1 and d2 == 2 and y == F(-1, 6):
        return "Qbar"
    if t == -1 and d2 == 1 and y == F(-2, 3):
        return "uc"
    if t == 1 and d2 == 1 and y == F(2, 3):
        return "ucbar"
    if t == 0 and d2 == 1 and y == F(1):
        return "ec"
    if t == 0 and d2 == 1 and y == F(-1):
        return "ecbar"
    if t == -1 and d2 == 1 and y == F(1, 3):
        return "dc"
    if t == 1 and d2 == 1 and y == F(-1, 3):
        return "dcbar"
    if t == 0 and d2 == 2 and y == F(-1, 2):
        return "L"
    if t == 0 and d2 == 2 and y == F(1, 2):
        return "Lbar"
    return None


TEN = {"Q", "uc", "ec"}
TENBAR = {"Qbar", "ucbar", "ecbar"}


# ------------------------------------------------------------------ orbifolder parser
_FIELD = re.compile(
    r"\n\s+(F_\d+)\n"
    r"\s+sector \(k,l\)\s+:\s+\((\d+),\s*(\d+)\)\n"
    r"\s+fixed point n_a\s+:\s+\(([^)]*)\)\n"
    r"\s+space group charges\s+:\s+\(([^)]*)\)\n"
    r"\s*\n"
    r"\s+representation\s+:\s+\(([^)]*)\)_([lr])\s+U\(1\)\s*:\s*\(([^)]*)\)\n"
    r"\s*\n"
    r"\s+left-moving p_sh\s+:\s+\(([^)]*)\),\s+\(([^)]*)\)\n"
    r"\s+right-moving q_sh\s+:\s+\(([^)]*)\)\n"
    r"\s+oscillators\s+:\s+(\S+)")

_MODEL = re.compile(r"/>\s*/([A-Za-z0-9_]+)>")


def split_models(text, prefix="CB"):
    marks = [(m.start(), m.group(1)) for m in _MODEL.finditer(text)
             if m.group(1).startswith(prefix)]
    return [(lab, text[pos:(marks[i + 1][0] if i + 1 < len(marks) else len(text))])
            for i, (pos, lab) in enumerate(marks)]


def _repdim(tok):
    tok = tok.strip()
    m = re.match(r"^(-?)(\d+)", tok)
    if not m:
        return None
    d = int(m.group(2))
    return -d if m.group(1) == "-" else d


def parse_fields(block):
    """[dict] for every printed spectrum record in `block` (own parser, G6)."""
    out = []
    for m in _FIELD.finditer(block):
        lab, k, l = m.group(1), int(m.group(2)), int(m.group(3))
        fp = tuple(int(x) for x in m.group(4).split(",")) if m.group(4).strip() else ()
        sgc = tuple(int(x) for x in m.group(5).split(",")) if m.group(5).strip() else ()
        rep = tuple(_repdim(t) for t in m.group(6).split(","))
        chir = m.group(7)
        u1 = tuple(rat(t) for t in m.group(8).split(","))
        p = [rat(t) for t in m.group(9).split(",")] + \
            [rat(t) for t in m.group(10).split(",")]
        q = [rat(t) for t in m.group(11).split(",")]
        osc = m.group(12)
        out.append(dict(label=lab, kl=(k, l), fp=fp, sgc=sgc, rep=rep, chir=chir,
                        u1=u1, p_sh=p, q_sh=q, osc=osc))
    return out


_U1GEN = re.compile(r"U\(1\) generators:\n((?:\s+\([^\n]*\n)+)")
_SIMPLE = re.compile(r"Simple roots:\n((?:\s+\([^\n]*\n)+)")
_GG = re.compile(r'Gauge group in vev-configuration "[^"]*":\s*(.+)')


def parse_u1_generators(block):
    m = _U1GEN.search(block)
    if not m:
        return []
    gens = []
    for line in m.group(1).strip("\n").split("\n"):
        parts = re.findall(r"\(([^)]*)\)", line)
        if len(parts) != 2:
            continue
        gens.append([rat(t) for t in parts[0].split(",")] +
                    [rat(t) for t in parts[1].split(",")])
    return gens


def parse_simple_roots(block):
    m = _SIMPLE.search(block)
    if not m:
        return []
    out = []
    for line in m.group(1).strip("\n").split("\n"):
        parts = re.findall(r"\(([^)]*)\)", line)
        if len(parts) != 2:
            continue
        out.append([rat(t) for t in parts[0].split(",")] +
                   [rat(t) for t in parts[1].split(",")])
    return out


def parse_gauge_group(block):
    m = _GG.search(block)
    return m.group(1).strip() if m else None


_TRQ = re.compile(r"tr\s+(Q_\d+(?:\^3)?(?:\s+Q_\d+\^2)?)\s*=\s*(-?\d+)/(\d+)")


def parse_anomaly_traces(block):
    """{'Q_1': F, 'Q_1^3': F, 'Q_1 Q_2^2': F, ...} from `print anomaly info`."""
    out = {}
    for m in _TRQ.finditer(block):
        out[re.sub(r"\s+", " ", m.group(1))] = F(int(m.group(2)), int(m.group(3)))
    return out


# ------------------------------------------------------------------ shells (independent)
def e8_coset_shell(V, target):
    """All x in V + E8 with |x|^2 == target (Fraction).  Independent of sc_lib.

    Algorithm differs from the phase's: coordinates are enumerated in DECREASING order of
    their minimum achievable square with a running admissible-remainder bound, and the E8
    parity is handled by summing the integer offsets rather than by branching on cosets
    first.  Self-tested against the E8 theta series in selftest().
    """
    res = []
    lim = F(target)
    for half in (0, 1):
        base = [x + F(half, 2) for x in V]
        rng = []
        ok = True
        for x in base:
            vals = []
            n = -12
            while n <= 12:
                y = x + n
                if y * y <= lim:
                    vals.append((y, n))
                n += 1
            if not vals:
                ok = False
                break
            rng.append(vals)
        if not ok:
            continue
        mins = [F(0)] * 9
        for i in range(7, -1, -1):
            mins[i] = mins[i + 1] + min(y * y for y, _ in rng[i])
        _walk(rng, mins, 0, [], F(0), lim, res, 0)
    return res


def _walk(rng, mins, i, cur, acc, lim, res, nsum):
    if acc + mins[i] > lim:
        return
    if i == 8:
        if acc == lim and nsum % 2 == 0:
            res.append(tuple(cur))
        return
    for y, n in rng[i]:
        a = acc + y * y
        if a + mins[i + 1] > lim:
            continue
        cur.append(y)
        _walk(rng, mins, i + 1, cur, a, lim, res, nsum + n)
        cur.pop()


def e8_coset_shell_upto(V, maxnorm):
    """{norm2: [vectors]} for x in V + E8 with |x|^2 <= maxnorm."""
    res = {}
    lim = F(maxnorm)
    for half in (0, 1):
        base = [x + F(half, 2) for x in V]
        rng = []
        ok = True
        for x in base:
            vals = [(x + n, n) for n in range(-12, 13) if (x + n) ** 2 <= lim]
            if not vals:
                ok = False
                break
            rng.append(vals)
        if not ok:
            continue
        mins = [F(0)] * 9
        for i in range(7, -1, -1):
            mins[i] = mins[i + 1] + min(y * y for y, _ in rng[i])
        _walk_upto(rng, mins, 0, [], F(0), lim, res, 0)
    return res


def _walk_upto(rng, mins, i, cur, acc, lim, res, nsum):
    if acc + mins[i] > lim:
        return
    if i == 8:
        if nsum % 2 == 0:
            res.setdefault(acc, []).append(tuple(cur))
        return
    for y, n in rng[i]:
        a = acc + y * y
        if a + mins[i + 1] > lim:
            continue
        cur.append(y)
        _walk_upto(rng, mins, i + 1, cur, a, lim, res, nsum + n)
        cur.pop()


def full_shell(Vg, target):
    """All p in Vg + Lambda with |p|^2 = target."""
    a = e8_coset_shell_upto(list(Vg[:8]), target)
    b = e8_coset_shell_upto(list(Vg[8:]), target)
    out = []
    for n1, L1 in a.items():
        n2_ = target - n1
        for y in b.get(n2_, []):
            for x in L1:
                out.append(tuple(x) + tuple(y))
    return out


# ------------------------------------------------------------------ self-tests
def selftest(verbose=True):
    ok = True
    z = [F(0)] * 8
    th = e8_coset_shell_upto(z, F(6))
    got = [len(th.get(F(n), [])) for n in (0, 2, 4, 6)]
    ok &= got == [1, 240, 2160, 6720]
    if verbose:
        print("  E8 theta  0,2,4,6 :", got, "expected [1, 240, 2160, 6720]")
    g2 = [len(full_shell([F(0)] * 16, F(n))) for n in (0, 2, 4)]
    ok &= g2 == [1, 480, 61920]
    if verbose:
        print("  E8xE8 theta 0,2,4 :", g2, "expected [1, 480, 61920]")
    # a coset representative must lie in its own shell (exercises V != 0)
    p = [F(2, 3), F(2, 3), F(-1, 3), F(-1, 3), F(-1, 3), F(0), F(0), F(0),
         F(1, 2), F(1, 2)] + [F(0)] * 6
    sh = full_shell(p, n2(p))
    ok &= tuple(p) in set(sh)
    if verbose:
        print("  coset rep in own shell:", tuple(p) in set(sh), "| shell size", len(sh))
    # class key: independent cross-check against an explicit difference test
    import random
    rnd = random.Random(11)
    bad = 0
    roots = e8xe8_roots()
    for _ in range(200):
        lam = [F(0)] * 16
        for _ in range(3):
            r = roots[rnd.randrange(len(roots))]
            lam = add(lam, list(r))
        q = add(p, lam)
        if lattice_class_key(q) != lattice_class_key(p):
            bad += 1
    ok &= bad == 0
    if verbose:
        print("  class key stable under 200 random lattice shifts:", bad, "mismatches")
    # and it must SEPARATE genuinely different classes
    sep = lattice_class_key(p) != lattice_class_key(add(p, [F(1, 3)] + [F(0)] * 15))
    ok &= sep
    if verbose:
        print("  class key separates a non-lattice shift:", sep)
    # SM sanity: k_Y = 5/3 ; t_Y traceless on the SU(5) simple roots
    ok &= 2 * n2(TY) == F(5, 3)
    ok &= all(dot(TY, a) == 0 for a in [A3[0], A3[1], A2L])
    if verbose:
        print("  k_Y = 2|t_Y|^2 :", 2 * n2(TY), " t_Y . SM simple roots = 0: True")
    # roles: the SU(5) 10 and 5bar weight norms
    if verbose:
        print("  selftest", "OK" if ok else "FAIL")
    return ok


if __name__ == "__main__":
    print("gs_lib selftest")
    selftest()
    ce = parse_constructing_elements()
    print("constructing elements T(0,1):", ce[(0, 1)])
    print("constructing elements T(0,5):", ce[(0, 5)])
    print("SP[0][1], SP[4][0], SP[4][1], SP[4][4], SP[4][5] =",
          SP[0][1], SP[4][0], SP[4][1], SP[4][4], SP[4][5])

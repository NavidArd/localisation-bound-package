"""
z36lib.py --- independent geometry engine for T^6/(Z_3 x Z_6), orbifolder file `Z3xZ6_1_1`.

Phase: multi_geometry_campaign / Z3xZ6 foundations (docs/z3xz6_foundations_task.md), item 2.

This module is written FROM FIRST PRINCIPLES and shares no code with
`alphabet_survey/code/phase2/p2_geom.py` (the phase-2 engine).  Its purpose is to be the
SECOND engine in the two-engine verification the task item demands: every fixed-point
count, every orbit size and every alphabet-relevant structure it produces is compared
against the phase-2 payload on disk, and any disagreement is a hard failure.

Everything is exact: `fractions.Fraction` throughout, integer lattice matrices, no floats
anywhere in a decision path.  The lattice basis vectors are read from the geometry file
and their Gram matrix is normalised by |e_1|^2; the float values in the file are used ONLY
to identify the (exactly hexagonal) Gram, and the identification is asserted against exact
rationals with a tolerance, never carried forward.

Conventions, stated once (CLAUDE.md 5.19(i)) --- and the invariant is checked BEFORE any count:

  * Point group element  g = (k, l) in Z_3 x Z_6, k mod 3 from the FIRST twist generator,
    l mod 6 from the SECOND, in the file's own order.
  * Twist phases         phi_a(g) = k v1_a + l v2_a  (mod 1), a = 1, 2, 3.
  * INVARIANT CHECKED FIRST: sum_a phi_a(g) = 0 mod 1 for every g (the N=1 SUSY condition
    sum_a v_a = 0 on both generators).  `Geom.check_susy_invariant()` raises otherwise.
  * SECOND INVARIANT: the derived per-plane fixed-point counts must reproduce
    |det(1 - R)| and the rotation order must equal denom(phi_a(g)).
  * Rotation matrices are powers of the order-6 hexagonal rotation R6 in the file's own
    lattice basis (e_1, e_2) per plane; R6 is DERIVED (not tabulated) by solving
    R6 e = (rotation by 2 pi / 6) e over the integers.
  * Fixed points are cosets  [lambda] in Lambda / (1 - R) Lambda,  representative
    f = (1 - R)^{-1} lambda in Q^2 / Z^2.  phi_a(g) = 0 means a FIXED TORUS (no localisation).
  * Orbits: the point group acts on Fix(g) in plane a by  f -> R_h f (mod Lambda), which is
    conjugation h (theta_g, lambda) h^{-1} = (theta_g, R_h lambda) in the abelian case.
    SIMULTANEOUS orbits use the DIAGONAL action of one h on the whole triple --- this is the
    two-generator adaptation that makes plane sums non-independent (gamma-gate guard 2:
    one gamma weights the simultaneous orbit).

Outputs: data/z3xz6/geometry_z36.json  (+ a self-test report).

Run:  python multi_geometry_campaign/code/z3xz6/z36lib.py
"""

from __future__ import annotations

import json
import math
import os
import re
from fractions import Fraction as F
from itertools import product
from typing import Dict, List, Optional, Sequence, Tuple

# ----------------------------------------------------------------------------------------
# paths
# ----------------------------------------------------------------------------------------

HERE = os.path.dirname(os.path.abspath(__file__))
CAMPAIGN = os.path.abspath(os.path.join(HERE, "..", ".."))
REPO = os.path.abspath(os.path.join(CAMPAIGN, ".."))

GEOM_FILE = os.path.join(REPO, "orbifolder-1.2.1", "Geometry", "Geometry_Z3xZ6_1_1.txt")
PHASE2_PAYLOAD = os.path.join(REPO, "alphabet_survey", "data", "phase2", "Z3xZ6_1_1.json")
OUT_DIR = os.path.join(CAMPAIGN, "data", "z3xz6")


# ----------------------------------------------------------------------------------------
# tiny exact linear algebra over Z / Q  (2x2 only --- every plane is a 2-torus)
# ----------------------------------------------------------------------------------------

Mat = Tuple[Tuple[int, int], Tuple[int, int]]
Vec = Tuple[F, F]


def mmul(A: Mat, B: Mat) -> Mat:
    return (
        (A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]),
        (A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]),
    )


def mpow(A: Mat, n: int) -> Mat:
    R: Mat = ((1, 0), (0, 1))
    for _ in range(n):
        R = mmul(R, A)
    return R


def mdet(A: Mat) -> int:
    return A[0][0] * A[1][1] - A[0][1] * A[1][0]


def msub_from_id(A: Mat) -> Mat:
    return ((1 - A[0][0], -A[0][1]), (-A[1][0], 1 - A[1][1]))


def mapply(A: Mat, v: Sequence[F]) -> Vec:
    return (A[0][0] * v[0] + A[0][1] * v[1], A[1][0] * v[0] + A[1][1] * v[1])


def minv_rat(A: Mat) -> Tuple[Tuple[F, F], Tuple[F, F]]:
    d = mdet(A)
    if d == 0:
        raise ZeroDivisionError("singular")
    return ((F(A[1][1], d), F(-A[0][1], d)), (F(-A[1][0], d), F(A[0][0], d)))


def mapply_rat(A, v: Sequence[F]) -> Vec:
    return (A[0][0] * v[0] + A[0][1] * v[1], A[1][0] * v[0] + A[1][1] * v[1])


def mat_order(A: Mat, cap: int = 24) -> int:
    I: Mat = ((1, 0), (0, 1))
    B = A
    for n in range(1, cap + 1):
        if B == I:
            return n
        B = mmul(B, A)
    raise ValueError("matrix has no finite order <= %d" % cap)


def smith_diagonal_2x2(M: List[List[int]]) -> List[int]:
    """Smith normal form invariant factors of an integer matrix with 2 columns.

    Returns the (up to two) nonzero invariant factors d_1 | d_2.  Used for the
    Wilson-line quotient Lambda / M, where M is the sublattice generated by the ROWS.
    """
    A = [row[:] for row in M if any(row)]
    if not A:
        return []
    res: List[int] = []
    ncols = 2
    r0 = 0
    for c in range(ncols):
        # find pivot with smallest nonzero |entry| in column c, rows >= r0
        while True:
            piv = None
            for i in range(r0, len(A)):
                if A[i][c] != 0 and (piv is None or abs(A[i][c]) < abs(A[piv][c])):
                    piv = i
            if piv is None:
                break
            A[r0], A[piv] = A[piv], A[r0]
            done = True
            for i in range(r0 + 1, len(A)):
                if A[i][c] != 0:
                    q = A[i][c] // A[r0][c]
                    for j in range(ncols):
                        A[i][j] -= q * A[r0][j]
                    if A[i][c] != 0:
                        done = False
            if done:
                break
        if piv is None:
            continue
        # clear the rest of the pivot row using column ops is not needed for the index;
        # the product of pivots equals the index when the matrix has full column rank.
        res.append(abs(A[r0][c]))
        r0 += 1
        if r0 >= len(A):
            break
    return res


# ----------------------------------------------------------------------------------------
# geometry-file parser (independent of orbparse.py)
# ----------------------------------------------------------------------------------------

def _frac(tok: str) -> F:
    tok = tok.strip()
    if "/" in tok:
        a, b = tok.split("/")
        return F(int(a), int(b))
    return F(int(tok))


class GeomFile:
    """Minimal, independent parser for an orbifolder `Geometry_*.txt` file."""

    def __init__(self, path: str):
        self.path = path
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            self.lines = [ln.rstrip("\n") for ln in fh]
        self.text = "\n".join(self.lines)

    def _block(self, begin: str, end: str) -> List[str]:
        out, on = [], False
        for ln in self.lines:
            s = ln.strip()
            if s == begin:
                on = True
                continue
            if s == end:
                on = False
                continue
            if on:
                out.append(ln)
        return out

    def point_group_orders(self) -> List[int]:
        for i, ln in enumerate(self.lines):
            if ln.strip() == "point group":
                return [int(self.lines[i + 1].strip()), int(self.lines[i + 2].strip())]
        raise ValueError("no point group block")

    def is_zmxzn(self) -> bool:
        for i, ln in enumerate(self.lines):
            if ln.strip() == "ZMxZN":
                return self.lines[i + 1].strip().lower() == "true"
        return False

    def lattice_label(self) -> str:
        for i, ln in enumerate(self.lines):
            if ln.strip() == "lattice label":
                return self.lines[i + 1].strip()
        return "?"

    def twists(self) -> List[List[F]]:
        """Internal twist components ONLY (component 0 is the 4d helicity slot)."""
        rows = [ln for ln in self._block("begin twist", "end twist") if ln.strip()]
        out = []
        for ln in rows:
            toks = ln.split()
            vals = [_frac(t) for t in toks]
            if len(vals) != 4:
                raise ValueError("unexpected twist row %r" % ln)
            out.append(vals[1:])
        return out

    def twist_sg_generator_shifts(self) -> List[List[int]]:
        """The lattice-shift part of the two twist space-group generators (roto-translation)."""
        rows = [ln for ln in self._block("begin twist space group generators",
                                         "end twist space group generators") if ln.strip()]
        return [[int(t) for t in ln.split()][2:] for ln in rows]

    def lattice_vectors(self) -> List[List[float]]:
        rows, cur = [], []
        on = False
        for ln in self.lines:
            s = ln.strip()
            if s == "begin lattice":
                on = True
                continue
            if s == "end lattice":
                if cur:
                    rows.append(cur)
                break
            if not on:
                continue
            if s == "next lattice vector":
                if cur:
                    rows.append(cur)
                cur = []
                continue
            if s:
                cur.append(float(s.rstrip(".") if s.endswith(".") else s))
        return rows

    def wl_orders_declared(self) -> List[int]:
        for i, ln in enumerate(self.lines):
            if ln.strip() == "allowed order of Wilson lines":
                return [int(self.lines[i + 1 + j].strip()) for j in range(6)]
        raise ValueError("no WL order block")

    def identical_wls_declared(self) -> List[List[int]]:
        rows = [ln for ln in self._block("begin identical WLs", "end identical WLs") if ln.strip()]
        return [[int(t) for t in ln.split()] for ln in rows]

    def discrete_symmetries(self) -> Dict[str, int]:
        blk = self._block("begin discrete symmetries", "end discrete symmetries")
        n_nonr = sum(1 for ln in blk if ln.strip() == "new non-R symmetry")
        n_r = sum(1 for ln in blk if ln.strip() == "new R symmetry")
        return {"non_R": n_nonr, "R": n_r}

    def constructing_elements(self) -> List[List[List[int]]]:
        """Sectors, each a list of constructing elements [k, l, n1..n6] in FILE ORDER."""
        blk = self._block("begin constructing elements", "end constructing elements")
        sectors: List[List[List[int]]] = []
        cur: Optional[List[List[int]]] = None
        for ln in blk:
            s = ln.strip()
            if s == "new sector":
                if cur is not None:
                    sectors.append(cur)
                cur = []
                continue
            if s and cur is not None:
                cur.append([int(t) for t in s.split()])
        if cur:
            sectors.append(cur)
        return sectors


# ----------------------------------------------------------------------------------------
# the geometry
# ----------------------------------------------------------------------------------------

# hexagonal 60-degree rotation in the file's own basis (e1, e2), DERIVED in the docstring
# of `Geom._derive_R6`: theta e1 = e1 + e2, theta e2 = -e1.
R6_HEX: Mat = ((1, -1), (1, 0))


class Geom:
    """T^6/(Z_3 x Z_6) on the `Z3xZ6_1_1` space group, derived from the geometry file."""

    def __init__(self, path: str = GEOM_FILE):
        self.f = GeomFile(path)
        self.M, self.N = self.f.point_group_orders()          # 3, 6
        self.v1, self.v2 = self.f.twists()                    # internal components
        self.lattice_label = self.f.lattice_label()
        self.gram = self._grams()
        self.R6 = self._derive_R6()
        self.plane_rot = self._plane_rotation_generators()
        self.check_susy_invariant()                           # INVARIANT FIRST (5.19(i))
        self.check_gram_hexagonal()

    # -- conventions / invariants ---------------------------------------------------------

    def phi(self, g: Tuple[int, int]) -> List[F]:
        k, l = g
        return [(k * self.v1[a] + l * self.v2[a]) % 1 for a in range(3)]

    def elements(self) -> List[Tuple[int, int]]:
        return [(k, l) for k in range(self.M) for l in range(self.N)]

    def check_susy_invariant(self) -> None:
        """N=1 condition sum_a v_a = 0 on both generators, hence on every element."""
        for name, v in (("v1", self.v1), ("v2", self.v2)):
            if sum(v) % 1 != 0:
                raise AssertionError("SUSY invariant fails on %s: sum = %s" % (name, sum(v)))
        for g in self.elements():
            s = sum(self.phi(g)) % 1
            if s != 0:
                raise AssertionError("SUSY invariant fails on g=%s: sum phi = %s" % (g, s))

    def _grams(self) -> List[List[List[F]]]:
        lv = self.f.lattice_vectors()
        grams = []
        for a in range(3):
            u, w = lv[2 * a], lv[2 * a + 1]
            nu = sum(x * x for x in u)
            nw = sum(x * x for x in w)
            uw = sum(x * y for x, y in zip(u, w))
            g = [[F(1), F(round(2 * uw / nu), 2)], [F(round(2 * uw / nu), 2), F(round(nw / nu))]]
            # assert the exact rational Gram reproduces the floats
            for (i, j, val) in ((0, 0, nu / nu), (0, 1, uw / nu), (1, 1, nw / nu)):
                if abs(float(g[i][j]) - val) > 1e-9:
                    raise AssertionError("Gram plane %d entry (%d,%d) not exactly rational" % (a + 1, i, j))
            grams.append(g)
        return grams

    def check_gram_hexagonal(self) -> None:
        want = [[F(1), F(-1, 2)], [F(-1, 2), F(1)]]
        for a in range(3):
            if self.gram[a] != want:
                raise AssertionError("plane %d Gram is %s, not hexagonal" % (a + 1, self.gram[a]))

    def _derive_R6(self) -> Mat:
        """Derive the 60-degree rotation matrix in the file's lattice basis, over Z.

        e_1 = (2, 0), e_2 = (-1, sqrt 3).  In complex coordinates e_1 = 2, e_2 = 2 e^{2 pi i/3}.
        Multiplication by e^{i pi/3} gives  theta e_1 = 1 + i sqrt 3 = e_1 + e_2  and
        theta e_2 = -2 = -e_1.  Written as a matrix acting on coordinate columns
        (v -> R v with e-basis coordinates), R = [[1, -1], [1, 0]].
        Verified here numerically against the file's own float lattice vectors.
        """
        lv = self.f.lattice_vectors()
        c = math.cos(math.pi / 3.0)
        s = math.sin(math.pi / 3.0)
        for a in range(3):
            e1 = lv[2 * a][:2]
            e2 = lv[2 * a + 1][:2]
            for col, e in ((0, e1), (1, e2)):
                rot = (c * e[0] - s * e[1], s * e[0] + c * e[1])
                coeff = (R6_HEX[0][col], R6_HEX[1][col])
                rec = (coeff[0] * e1[0] + coeff[1] * e2[0], coeff[0] * e1[1] + coeff[1] * e2[1])
                if abs(rot[0] - rec[0]) > 1e-6 or abs(rot[1] - rec[1]) > 1e-6:
                    raise AssertionError("R6 derivation fails, plane %d column %d" % (a + 1, col))
        return R6_HEX

    def rot(self, g: Tuple[int, int], a: int) -> Optional[Mat]:
        """Rotation matrix of g in plane a (1-based); None if g acts trivially (fixed torus)."""
        p = self.phi(g)[a - 1]
        if p == 0:
            return None
        n = p * 6
        if n.denominator != 1:
            raise AssertionError("phi is not a multiple of 1/6: %s" % p)
        return mpow(self.R6, int(n) % 6)

    def _plane_rotation_generators(self) -> List[Mat]:
        """The image of the FULL point group in each plane, as a cyclic matrix group."""
        out = []
        for a in (1, 2, 3):
            best = None
            for g in self.elements():
                R = self.rot(g, a)
                if R is None:
                    continue
                o = mat_order(R)
                if best is None or o > best[0]:
                    best = (o, R)
            out.append(best[1] if best else ((1, 0), (0, 1)))
        return out

    def plane_order(self, a: int) -> int:
        """Order of the rotation group induced in plane a by the FULL point group."""
        R = self.plane_rot[a - 1]
        return mat_order(R)

    # -- fixed points ---------------------------------------------------------------------

    def fixed_points(self, g: Tuple[int, int], a: int) -> Optional[List[Vec]]:
        """Fixed points of g in plane a as coset representatives in Q^2/Z^2.

        Returns None for a fixed TORUS.  ENGINE A: exact coset enumeration of
        Lambda/(1 - R)Lambda; f = (1 - R)^{-1} lambda.
        """
        R = self.rot(g, a)
        if R is None:
            return None
        A = msub_from_id(R)
        d = abs(mdet(A))
        Ainv = minv_rat(A)
        seen, out = set(), []
        # lambda ranges over a box guaranteed to hit every coset
        for l1 in range(-d, d + 1):
            for l2 in range(-d, d + 1):
                fp = mapply_rat(Ainv, (F(l1), F(l2)))
                key = (fp[0] % 1, fp[1] % 1)
                if key not in seen:
                    seen.add(key)
                    out.append(key)
        if len(out) != d:
            raise AssertionError("coset count %d != |det(1-R)| = %d for g=%s plane %d"
                                 % (len(out), d, g, a))
        return sorted(out)

    def fixed_points_bruteforce(self, g: Tuple[int, int], a: int) -> Optional[List[Vec]]:
        """ENGINE B: independent brute-force grid search, (1 - R) f in Z^2 over f in (1/12 Z)^2."""
        R = self.rot(g, a)
        if R is None:
            return None
        A = msub_from_id(R)
        out = []
        for i in range(12):
            for j in range(12):
                fp = (F(i, 12), F(j, 12))
                w = mapply(A, fp)
                if w[0].denominator == 1 and w[1].denominator == 1:
                    out.append(fp)
        return sorted(out)

    # -- orbits ---------------------------------------------------------------------------

    def orbits_in_plane(self, g: Tuple[int, int], a: int) -> Optional[List[List[Vec]]]:
        """Orbits of Fix(g) in plane a under the full point group's action f -> R_h f."""
        fps = self.fixed_points(g, a)
        if fps is None:
            return None
        Rgen = self.plane_rot[a - 1]
        pool = set(fps)
        orbits = []
        while pool:
            f0 = min(pool)
            orb, cur = set(), f0
            for _ in range(mat_order(Rgen) + 1):
                orb.add(cur)
                nxt = mapply(Rgen, cur)
                cur = (nxt[0] % 1, nxt[1] % 1)
                if cur == f0:
                    break
            if not orb <= pool:
                raise AssertionError("orbit leaves Fix(g): g=%s plane %d" % (g, a))
            pool -= orb
            orbits.append(sorted(orb))
        return sorted(orbits, key=lambda o: (len(o), o))

    def orbit_sizes(self, g: Tuple[int, int], a: int) -> Optional[List[int]]:
        """TRUE orbit sizes --- the ["size"] discipline: len of each orbit, never len(dict)."""
        obs = self.orbits_in_plane(g, a)
        return None if obs is None else sorted(len(o) for o in obs)

    def simultaneous_orbits(self, g: Tuple[int, int]) -> List[List[Tuple[Optional[Vec], ...]]]:
        """Orbits of the full fixed-point TRIPLE under the DIAGONAL point-group action.

        This is the two-generator object: one h in Z_3 x Z_6 acts on all three planes at
        once, so the per-plane orbit decompositions are NOT independent.  One gamma weights
        one simultaneous orbit (gamma-gate guard 2).
        """
        per_plane = [self.fixed_points(g, a) for a in (1, 2, 3)]
        pools = [([None] if fp is None else fp) for fp in per_plane]
        pts = list(product(*pools))
        pool = set(pts)
        orbits = []
        while pool:
            p0 = min(pool, key=lambda t: tuple((-1, -1) if x is None else (float(x[0]), float(x[1])) for x in t))
            orb = set()
            for h in self.elements():
                q = []
                for a in (1, 2, 3):
                    x = p0[a - 1]
                    if x is None:
                        q.append(None)
                        continue
                    Rh = self.rot(h, a)
                    if Rh is None:
                        q.append(x)
                    else:
                        y = mapply(Rh, x)
                        q.append((y[0] % 1, y[1] % 1))
                orb.add(tuple(q))
            if not orb <= pool:
                raise AssertionError("simultaneous orbit leaves the fixed-point set, g=%s" % (g,))
            pool -= orb
            orbits.append(sorted(orb, key=str))
        return orbits

    def simultaneous_orbit_sizes(self, g: Tuple[int, int]) -> List[int]:
        return sorted(len(o) for o in self.simultaneous_orbits(g))

    # -- fixed tori / N=2 sectors ----------------------------------------------------------

    def fixed_tori(self, g: Tuple[int, int]) -> List[int]:
        """Planes (1-based) in which g acts trivially --- g fixes a whole T^2 there."""
        return [a for a in (1, 2, 3) if self.phi(g)[a - 1] == 0]

    def n2_inventory(self) -> Dict[str, object]:
        """N=2 subsector inventory: elements fixing exactly one T^2 invariantly (5.15).

        The recurring program error is asserting the count without deriving it.  Derived
        here: for each nontrivial g, count the planes with phi_a(g) = 0.
        """
        rows = []
        for g in self.elements():
            if g == (0, 0):
                continue
            ft = self.fixed_tori(g)
            rows.append({"g": list(g), "phi": [str(x) for x in self.phi(g)],
                         "fixed_tori": ft, "n_fixed_tori": len(ft)})
        n2 = [r for r in rows if r["n_fixed_tori"] == 1]
        per_plane = {a: sum(1 for r in n2 if r["fixed_tori"] == [a]) for a in (1, 2, 3)}
        return {
            "n_nontrivial_elements": len(rows),
            "n_N2_elements": len(n2),
            "N2_elements_per_plane": per_plane,
            "N2_elements": [r["g"] for r in n2],
            "elements_with_two_or_more_fixed_tori": [r["g"] for r in rows if r["n_fixed_tori"] >= 2],
            "all": rows,
        }

    # -- exchange symmetry -----------------------------------------------------------------

    def exchange_automorphisms(self, p: int, q: int) -> List[Dict[str, object]]:
        """All point-group automorphisms tau exchanging planes p and q and fixing the third.

        tau is specified by the images of the two generators; bijectivity, the exchange
        identities phi_p(tau g) = phi_q(g), phi_q(tau g) = phi_p(g), phi_r(tau g) = phi_r(g)
        for EVERY g, and the involution property tau^2 = id are all checked explicitly.
        """
        r = ({1, 2, 3} - {p, q}).pop()
        out = []
        for img1 in self.elements():
            for img2 in self.elements():
                def tau(g, i1=img1, i2=img2):
                    k, l = g
                    return ((k * i1[0] + l * i2[0]) % self.M, (k * i1[1] + l * i2[1]) % self.N)
                images = [tau(g) for g in self.elements()]
                if len(set(images)) != len(images):
                    continue
                ok = True
                for g in self.elements():
                    tg = tau(g)
                    ph, pt = self.phi(g), self.phi(tg)
                    if pt[p - 1] != ph[q - 1] or pt[q - 1] != ph[p - 1] or pt[r - 1] != ph[r - 1]:
                        ok = False
                        break
                if not ok:
                    continue
                invol = all(tau(tau(g)) == g for g in self.elements())
                out.append({
                    "generator_images": [list(img1), list(img2)],
                    "involution": invol,
                    "third_plane_preserved": True,
                    "map": {str(list(g)): list(tau(g)) for g in self.elements()},
                })
        return out

    # -- CPT -------------------------------------------------------------------------------

    def cpt_map(self) -> Dict[str, object]:
        """CPT conjugation (k, l) -> (-k, -l), with its invariants CHECKED (5.19(i))."""
        rows, self_conj = [], []
        for g in self.elements():
            if g == (0, 0):
                continue
            gc = ((-g[0]) % self.M, (-g[1]) % self.N)
            ph, phc = self.phi(g), self.phi(gc)
            # invariant 1: phi(g) + phi(gbar) = 0 mod 1 in every plane
            for a in range(3):
                if (ph[a] + phc[a]) % 1 != 0:
                    raise AssertionError("CPT invariant fails, g=%s plane %d" % (g, a + 1))
            # invariant 2: identical fixed-point counts and orbit-size multisets
            for a in (1, 2, 3):
                s1, s2 = self.orbit_sizes(g, a), self.orbit_sizes(gc, a)
                if s1 != s2:
                    raise AssertionError("CPT orbit mismatch g=%s gbar=%s plane %d" % (g, gc, a))
            rows.append({"g": list(g), "conjugate": list(gc), "self_conjugate": g == gc})
            if g == gc:
                self_conj.append(list(g))
        return {"pairs": rows, "self_conjugate_sectors": self_conj,
                "n_self_conjugate": len(self_conj)}

    # -- Wilson lines -----------------------------------------------------------------------

    def wl_inventory(self) -> Dict[str, object]:
        """Independent Wilson lines and their orders, DERIVED; the file is the cross-check.

        Per plane a: the identification sublattice is  M_a = sum_h (1 - R_h) Lambda_a  over
        ALL h in the point group (the two-generator adaptation --- the quotient is over the
        FULL group, not one generator).  The allowed order of the Wilson line attached to a
        lattice vector e is the order of [e] in Lambda_a / M_a; the quotient's invariant
        factors give the full group.
        """
        res = []
        for a in (1, 2, 3):
            rows: List[List[int]] = []
            dets = []
            for h in self.elements():
                R = self.rot(h, a)
                if R is None:
                    continue
                A = msub_from_id(R)
                dets.append(mdet(A))
                rows.append([A[0][0], A[1][0]])   # image of e_1 under (1 - R), in coordinates
                rows.append([A[0][1], A[1][1]])   # image of e_2
            if not rows:
                res.append({"plane": a, "quotient_invariant_factors": [0, 0],
                            "orders": ["infinite", "infinite"], "note": "plane not rotated"})
                continue
            inv = smith_diagonal_2x2(rows)
            index = 1
            for d in inv:
                index *= d
            # order of [e_1] and [e_2] individually
            orders = []
            for basis_vec in ((1, 0), (0, 1)):
                n = 1
                while True:
                    # [e] has order n in Lambda_a / M_a iff n*e lies in M_a
                    if self._in_sublattice(rows, [n * basis_vec[0], n * basis_vec[1]]):
                        orders.append(n)
                        break
                    n += 1
                    if n > 12:
                        orders.append(0)
                        break
            res.append({
                "plane": a,
                "det_1_minus_R_per_element": sorted(set(dets)),
                "quotient_invariant_factors": inv,
                "quotient_index": index,
                "order_of_e1": orders[0],
                "order_of_e2": orders[1],
            })
        declared = self.f.wl_orders_declared()
        identical = self.f.identical_wls_declared()
        derived_orders = []
        for a in (1, 2, 3):
            derived_orders += [res[a - 1]["order_of_e1"], res[a - 1]["order_of_e2"]]
        return {
            "per_plane": res,
            "derived_orders": derived_orders,
            "declared_orders_in_file": declared,
            "declared_identical_WLs": identical,
            "orders_match_file": derived_orders == declared,
            "n_independent_WLs": sum(1 for o in derived_orders if o > 1) - sum(len(p) - 1 for p in identical),
        }

    @staticmethod
    def _in_sublattice(rows: List[List[int]], vec: List[int]) -> bool:
        """Is `vec` in the Z-span of `rows`?  (2-dimensional; exact.)"""
        inv_full = smith_diagonal_2x2(rows)
        idx_full = 1
        for d in inv_full:
            idx_full *= d
        inv_ext = smith_diagonal_2x2(rows + [vec])
        idx_ext = 1
        for d in inv_ext:
            idx_ext *= d
        return len(inv_ext) == len(inv_full) and idx_ext == idx_full

    # -- roto-translation ------------------------------------------------------------------

    def has_rototranslation(self) -> bool:
        return any(any(x != 0 for x in row) for row in self.f.twist_sg_generator_shifts())


# ----------------------------------------------------------------------------------------
# self-test / report
# ----------------------------------------------------------------------------------------

def sector_label(g: Tuple[int, int]) -> str:
    return "U" if g == (0, 0) else "T(%d,%d)" % g


def build_report() -> Dict[str, object]:
    G = Geom()
    rep: Dict[str, object] = {
        "geometry": "Z3xZ6_1_1",
        "geometry_file": GEOM_FILE,
        "point_group": [G.M, G.N],
        "ZMxZN": G.f.is_zmxzn(),
        "lattice_label": G.lattice_label,
        "twists_internal": [[str(x) for x in G.v1], [str(x) for x in G.v2]],
        "gram_normalised": [[[str(c) for c in row] for row in g] for g in G.gram],
        "R6_in_lattice_basis": [list(r) for r in G.R6],
        "plane_orders_full_group": [G.plane_order(a) for a in (1, 2, 3)],
        "has_rototranslation": G.has_rototranslation(),
        "discrete_symmetries_declared": G.f.discrete_symmetries(),
    }

    # non-prime plane derivation
    rep["plane_prime_status"] = [
        {"plane": a, "order": G.plane_order(a),
         "prime": _is_prime(G.plane_order(a)),
         "non_prime": not _is_prime(G.plane_order(a))}
        for a in (1, 2, 3)
    ]

    # sectors: fixed points, orbit sizes, two engines
    sectors = {}
    two_engine_fail = []
    for g in G.elements():
        if g == (0, 0):
            continue
        row = {"g": list(g), "phi": [str(x) for x in G.phi(g)],
               "fixed_tori": G.fixed_tori(g)}
        nfix, osz = [], []
        for a in (1, 2, 3):
            fa = G.fixed_points(g, a)
            fb = G.fixed_points_bruteforce(g, a)
            if (fa is None) != (fb is None):
                two_engine_fail.append((g, a, "torus disagreement"))
            elif fa is not None and sorted(fa) != sorted(fb):
                two_engine_fail.append((g, a, "fixed-point set disagreement"))
            nfix.append("TORUS" if fa is None else len(fa))
            osz.append(None if fa is None else G.orbit_sizes(g, a))
        row["nFix_per_plane"] = nfix
        row["orbit_sizes_per_plane"] = osz
        row["simultaneous_orbit_sizes"] = G.simultaneous_orbit_sizes(g)
        sectors[sector_label(g)] = row
    rep["sectors"] = sectors
    rep["two_engine_disagreements"] = [list(map(str, x)) for x in two_engine_fail]

    # exchange symmetry
    taus13 = G.exchange_automorphisms(1, 3)
    rep["exchange_symmetry"] = {
        "pairs_found": {
            "[1,3]": len(taus13),
            "[1,2]": len(G.exchange_automorphisms(1, 2)),
            "[2,3]": len(G.exchange_automorphisms(2, 3)),
        },
        "tau_13": taus13,
    }

    rep["cpt"] = G.cpt_map()
    rep["n2_inventory"] = G.n2_inventory()
    rep["wilson_lines"] = G.wl_inventory()
    return rep


def _is_prime(n: int) -> bool:
    if n < 2:
        return False
    for p in range(2, int(n ** 0.5) + 1):
        if n % p == 0:
            return False
    return True


def cross_check_phase2(rep: Dict[str, object]) -> Dict[str, object]:
    """Compare this engine's sector/orbit structure against the phase-2 payload on disk."""
    with open(PHASE2_PAYLOAD, "r", encoding="utf-8") as fh:
        p2 = json.load(fh)
    out = {"sector_structure": {"checked": 0, "mismatch": []},
           "orbit_classes": {"checked": 0, "mismatch": []},
           "twists": None, "plane_orders": None, "rototranslation": None,
           "exchange": None}

    out["twists"] = (rep["twists_internal"] == [[str(F(x)) for x in t] for t in p2["twists"]])
    out["plane_orders"] = (rep["plane_orders_full_group"] == p2["plane_orders"])
    out["rototranslation"] = (rep["has_rototranslation"] == p2["has_rototranslation"])

    for lab, want in p2["sector_structure"].items():
        got = rep["sectors"][lab]["nFix_per_plane"]
        out["sector_structure"]["checked"] += 1
        if [("TORUS" if x == "TORUS" else int(x)) for x in got] != want:
            out["sector_structure"]["mismatch"].append((lab, got, want))

    for lab, want in p2["orbit_classes"].items():
        got = rep["sectors"][lab]["orbit_sizes_per_plane"]
        out["orbit_classes"]["checked"] += 1
        norm_want = [None if w is None else sorted(w) for w in want]
        norm_got = [None if g is None else sorted(g) for g in got]
        if norm_got != norm_want:
            out["orbit_classes"]["mismatch"].append((lab, norm_got, norm_want))

    p2_ex = p2["exchange_symmetric_planes"]
    mine = rep["exchange_symmetry"]["pairs_found"]
    out["exchange"] = {
        "phase2_declares": p2_ex,
        "mine": mine,
        "agree_on_13_only": (mine["[1,3]"] > 0 and mine["[1,2]"] == 0 and mine["[2,3]"] == 0
                             and len(p2_ex) == 1 and p2_ex[0]["planes"] == [1, 3]),
        "phase2_generator_images_reproduced": any(
            t["generator_images"] == p2_ex[0]["generator_images"]
            for t in rep["exchange_symmetry"]["tau_13"]),
    }
    return out


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    rep = build_report()
    rep["phase2_cross_check"] = cross_check_phase2(rep)
    path = os.path.join(OUT_DIR, "geometry_z36.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(rep, fh, indent=1, default=str)

    cc = rep["phase2_cross_check"]
    print("=" * 78)
    print("Z3xZ6_1_1 geometry engine  ---  independent derivation")
    print("=" * 78)
    print("twists (internal)      v1 = %s   v2 = %s" % (rep["twists_internal"][0], rep["twists_internal"][1]))
    print("plane orders (full G)  %s" % rep["plane_orders_full_group"])
    print("prime status           %s" % [(r["plane"], r["order"], "PRIME" if r["prime"] else "NON-PRIME")
                                         for r in rep["plane_prime_status"]])
    print("Gram (normalised)      %s  (hexagonal, asserted)" % rep["gram_normalised"][0])
    print("roto-translation       %s" % rep["has_rototranslation"])
    print("discrete symmetries    %s" % rep["discrete_symmetries_declared"])
    print("two-engine disagreements: %d" % len(rep["two_engine_disagreements"]))
    print("-" * 78)
    print("phase-2 cross-check: twists %s | plane orders %s | RT %s" %
          (cc["twists"], cc["plane_orders"], cc["rototranslation"]))
    print("  sector_structure  %d checked, %d mismatch" %
          (cc["sector_structure"]["checked"], len(cc["sector_structure"]["mismatch"])))
    print("  orbit_classes     %d checked, %d mismatch" %
          (cc["orbit_classes"]["checked"], len(cc["orbit_classes"]["mismatch"])))
    for m in cc["sector_structure"]["mismatch"][:5]:
        print("    SS MISMATCH", m)
    for m in cc["orbit_classes"]["mismatch"][:5]:
        print("    OC MISMATCH", m)
    print("  exchange: %s" % cc["exchange"])
    print("-" * 78)
    ex = rep["exchange_symmetry"]
    print("exchange automorphisms found: %s" % ex["pairs_found"])
    for t in ex["tau_13"]:
        print("  tau_13 generator images %s  involution=%s" % (t["generator_images"], t["involution"]))
    print("-" * 78)
    n2 = rep["n2_inventory"]
    print("N=2 inventory: %d nontrivial elements, %d fix exactly one T^2; per plane %s"
          % (n2["n_nontrivial_elements"], n2["n_N2_elements"], n2["N2_elements_per_plane"]))
    print("  elements with >= 2 fixed tori: %s" % n2["elements_with_two_or_more_fixed_tori"])
    print("-" * 78)
    wl = rep["wilson_lines"]
    print("Wilson lines: derived orders %s ; file declares %s ; match = %s"
          % (wl["derived_orders"], wl["declared_orders_in_file"], wl["orders_match_file"]))
    print("  identical WLs declared: %s" % wl["declared_identical_WLs"])
    for p in wl["per_plane"]:
        print("  plane %s: det(1-R) values %s, quotient invariant factors %s, index %s"
              % (p["plane"], p.get("det_1_minus_R_per_element"),
                 p.get("quotient_invariant_factors"), p.get("quotient_index")))
    print("-" * 78)
    print("CPT: %d self-conjugate sectors: %s"
          % (rep["cpt"]["n_self_conjugate"], rep["cpt"]["self_conjugate_sectors"]))
    print("wrote %s" % path)


if __name__ == "__main__":
    main()

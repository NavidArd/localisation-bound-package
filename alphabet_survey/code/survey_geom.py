"""ALPHABET SURVEY -- generalized orbifold geometry engine.

Direct generalization of the GATE-CERTIFIED engine
`hierarchy_analysis/code/gate/gate_geom.py` (see `results/gate_review.md`, VERDICT PASS)
from `Z6-II_1_1` to any of the 138 orbifolder geometry files.

Everything is PARSED or DERIVED, never assumed:
  * point group order(s), the ZMxZN flag and the twist vector(s) are parsed;
  * the six real-space lattice vectors are parsed;
  * factorizability is TESTED (each complex plane must be spanned by exactly two
    lattice vectors having no components outside that plane) -- a lattice failing the
    test is NOT-COMPUTABLE-YET, never forced;
  * each plane's Gram matrix is the TRUE per-plane Gram, exact and rational
    (gate defect G-5: plane 3 of Z6-II_1_1 is a SQUARE lattice, not hexagonal);
  * theta restricted to a plane is DERIVED by rotating the parsed basis vectors in
    real space and solving for integral coordinates -- it is never hard-coded;
  * the per-geometry INVARIANT CHECK is against orbifolder's own printed
    `begin constructing elements` block (the fixed-point data the file ships).

Scope boundary, pre-declared in `docs/criterion.md` sec.5 and enforced here:
this module refuses to build a `Geom` for a ZMxZN geometry, because the certified
source basis (1301.2322) states in its own footnote 1 that it does NOT treat
Z_N x Z_M.  `docs/task.md` scope rule (ii): flag rather than extend.
"""

import math
import os
from fractions import Fraction as F

I2 = ((F(1), F(0)), (F(0), F(1)))
Z2 = ((F(0), F(0)), (F(0), F(0)))


# ------------------------------------------------------------------ parsing
def parse_geometry(path):
    """Parse an orbifolder Geometry_*.txt into a raw dict. Nothing is interpreted."""
    lines = [l.rstrip("\n") for l in open(path)]
    st = [l.strip() for l in lines]
    out = {}

    i = st.index("point group")
    out["M"] = int(st[i + 1])
    out["N"] = int(st[i + 2])

    i = st.index("ZMxZN")
    out["ZMxZN"] = (st[i + 1].lower() == "true")

    i = st.index("lattice label")
    out["lattice_label"] = st[i + 1]

    # ---- twist vector(s): each line is "<time> v1 v2 v3" as rationals
    i = st.index("begin twist")
    j = st.index("end twist")
    twists = []
    for l in st[i + 1:j]:
        if not l:
            continue
        toks = l.split()
        v = []
        for tok in toks[1:]:                     # toks[0] is the 0/1 time component
            num, den = tok.split("/")
            v.append(F(int(num), int(den)))
        assert len(v) == 3, (path, l)
        twists.append(tuple(v))
    out["twists"] = twists

    # ---- lattice vectors
    i = st.index("begin lattice")
    j = st.index("end lattice")
    vecs, cur = [], None
    for l in st[i + 1:j]:
        if l == "next lattice vector":
            if cur is not None:
                vecs.append(cur)
            cur = []
        elif l:
            cur.append(float(l.rstrip(".")))
    if cur:
        vecs.append(cur)
    assert len(vecs) == 6 and all(len(x) == 6 for x in vecs), path
    out["vecs"] = vecs

    # ---- constructing elements, grouped by sector (orbifolder's printed fixed-point data)
    out["sectors"] = []
    if "begin constructing elements" in st:
        i = st.index("begin constructing elements")
        j = st.index("end constructing elements")
        cur = None
        for l in st[i + 1:j]:
            if l == "new sector":
                if cur is not None:
                    out["sectors"].append(cur)
                cur = []
            elif l:
                # entries may be integers or rationals (roto-translation lattices)
                toks = [F(x) for x in l.split()]
                assert len(toks) == 8, (path, l)
                cur.append((int(toks[0]), int(toks[1]), tuple(toks[2:])))
        if cur is not None:
            out["sectors"].append(cur)
    return out


def plane_pairs(vecs, diagnostic=False):
    """Which pair of lattice vectors spans each complex plane, or None if the
    lattice does not factorize over the three planes.

    With `diagnostic=True` returns instead the per-plane count of lattice vectors
    that live entirely inside that plane -- the measure of HOW FAR a lattice is from
    factorizing, which `docs/task.md` requires be reported for the
    NOT-COMPUTABLE-YET boundary ("listed with what closing them needs")."""
    counts, out, used = [], [], set()
    for a in range(3):
        c0, c1 = 2 * a, 2 * a + 1
        idx = [i for i, w in enumerate(vecs)
               if abs(w[c0]) + abs(w[c1]) > 1e-9
               and all(abs(w[c]) < 1e-9 for c in range(6) if c not in (c0, c1))]
        counts.append(len(idx))
        if len(idx) == 2 and not (used & set(idx)):
            used |= set(idx)
            out.append(tuple(idx))
        else:
            out.append(None)
    if diagnostic:
        return counts
    return None if any(p is None for p in out) else out


def gram(vecs, pair):
    """Exact rational Gram of the plane, normalised by |e_first|^2."""
    def dot(u, w):
        return sum(u[i] * w[i] for i in range(6))
    i, j = pair
    n = dot(vecs[i], vecs[i])
    raw = [[dot(vecs[i], vecs[i]) / n, dot(vecs[i], vecs[j]) / n],
           [dot(vecs[j], vecs[i]) / n, dot(vecs[j], vecs[j]) / n]]
    out = []
    for r in raw:
        rr = []
        for x in r:
            fr = F(x).limit_denominator(60)
            if abs(float(fr) - x) > 1e-7:
                return None                      # not a nice rational Gram
            rr.append(fr)
        out.append(tuple(rr))
    return (out[0], out[1])


def rot_in_plane(vecs, pair, plane_index, angle_turns):
    """theta restricted to a plane as an EXACT integer matrix on the (e_i, e_j) basis.

    DERIVED: rotate each basis vector by `angle_turns` full turns in real space,
    then solve for its integral coordinates.  Returns None if the image is not
    integral (which would mean the rotation is not a lattice automorphism).
    """
    th = 2 * math.pi * float(angle_turns)
    c, s = math.cos(th), math.sin(th)
    i, j = pair
    c0, c1 = 2 * plane_index, 2 * plane_index + 1
    ei = (vecs[i][c0], vecs[i][c1])
    ej = (vecs[j][c0], vecs[j][c1])
    det = ei[0] * ej[1] - ei[1] * ej[0]
    if abs(det) < 1e-9:
        return None
    cols = []
    for e in (ei, ej):
        r = (c * e[0] - s * e[1], s * e[0] + c * e[1])
        x = (r[0] * ej[1] - r[1] * ej[0]) / det
        y = (ei[0] * r[1] - ei[1] * r[0]) / det
        xi, yi = round(x), round(y)
        if abs(x - xi) > 1e-6 or abs(y - yi) > 1e-6:
            return None
        cols.append((F(xi), F(yi)))
    return ((cols[0][0], cols[1][0]), (cols[0][1], cols[1][1]))


# ------------------------------------------------------------------ 2x2 exact algebra
def mm(A, B):
    return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(2)) for j in range(2))
                 for i in range(2))


def mpow(A, n):
    R = I2
    for _ in range(n):
        R = mm(R, A)
    return R


def ap(A, v):
    return (A[0][0] * v[0] + A[0][1] * v[1], A[1][0] * v[0] + A[1][1] * v[1])


def minv(A):
    d = A[0][0] * A[1][1] - A[0][1] * A[1][0]
    return None if d == 0 else ((A[1][1] / d, -A[0][1] / d),
                                (-A[1][0] / d, A[0][0] / d))


def one_minus(A):
    return ((1 - A[0][0], -A[0][1]), (-A[1][0], 1 - A[1][1]))


def modl(f):
    return (f[0] - (f[0] // 1), f[1] - (f[1] // 1))


class NotComputable(Exception):
    """Raised when a geometry falls outside the pre-declared computable class."""
    def __init__(self, code, detail=""):
        Exception.__init__(self, "%s: %s" % (code, detail))
        self.code = code
        self.detail = detail


class Geom(object):
    """A factorizable Z_N orbifold geometry, with everything derived from the file."""

    def __init__(self, path, allow_ZMxZN=False):
        self.path = path
        self.name = os.path.basename(path)[len("Geometry_"):-len(".txt")]
        raw = parse_geometry(path)
        self.raw = raw
        self.M, self.N = raw["M"], raw["N"]
        self.ZMxZN = raw["ZMxZN"]

        if self.ZMxZN and not allow_ZMxZN:
            raise NotComputable(
                "NOT-COMPUTABLE-YET/source-basis",
                "Z_M x Z_N twist structure; 1301.2322 footnote 1 states it does NOT "
                "treat Z_N x Z_M ('we could make an analogous discussion'). "
                "criterion.md sec.5(2): flag, never extend by analogy.")

        self.vecs = raw["vecs"]
        self.pairs = plane_pairs(self.vecs)
        if self.pairs is None:
            cnt = plane_pairs(self.vecs, diagnostic=True)
            raise NotComputable(
                "NOT-COMPUTABLE-YET/non-factorizable",
                "lattice does not split into three orthogonal complex planes spanned "
                "by two lattice vectors each (per-plane in-plane lattice-vector counts "
                "%s, need [2,2,2]); the plane-wise machinery of the certified template "
                "does not apply (rules_z6ii.md sec.1.4 rule-gap class). Planes that DO "
                "factorize: %d of 3." % (cnt, sum(1 for c in cnt if c == 2)))

        self.gram = [gram(self.vecs, p) for p in self.pairs]
        if any(g is None for g in self.gram):
            raise NotComputable("NOT-COMPUTABLE-YET/irrational-gram",
                                "per-plane Gram is not rational at denominator <= 60")

        self.v = raw["twists"][0]
        self.theta = []
        for a in range(3):
            R = rot_in_plane(self.vecs, self.pairs[a], a, self.v[a])
            if R is None:
                raise NotComputable(
                    "NOT-COMPUTABLE-YET/theta-not-integral",
                    "theta does not act integrally on the parsed basis in plane %d" % (a + 1))
            self.theta.append(R)

        self._order = [self._ord(self.theta[a]) for a in range(3)]

        # theta must be an isometry of each plane's Gram -- CLAUDE.md 5.19(i)
        for a in range(3):
            for e in ((F(1), F(0)), (F(0), F(1)), (F(1), F(1))):
                if self.rho(a, ap(self.theta[a], e)) != self.rho(a, e):
                    raise NotComputable(
                        "NOT-COMPUTABLE-YET/theta-not-isometry",
                        "theta does not preserve the plane-%d Gram form" % (a + 1))

    # ---------------- basics
    @staticmethod
    def _ord(A):
        n, R = 1, A
        while R != I2:
            R = mm(R, A)
            n += 1
            if n > 24:
                raise NotComputable("NOT-COMPUTABLE-YET/theta-infinite-order", "")
        return n

    def order(self, a):
        return self._order[a]

    def rho(self, a, u):
        G = self.gram[a]
        return (G[0][0] * u[0] * u[0] + 2 * G[0][1] * u[0] * u[1]
                + G[1][1] * u[1] * u[1])

    def rot(self, a, k):
        return mpow(self.theta[a], k % self._order[a])

    def frac_twist(self, a, k):
        x = k * self.v[a]
        return x - (x // 1)

    def rigid_plane(self, a):
        """Is the plane's lattice SHAPE fixed by the point group?

        If theta restricted to the plane has order <= 2 the Gram is not determined by
        the point group -- a complex-structure modulus survives and the plane's NORM
        spectrum is not lattice-rigid.  Reported, never silently used.
        """
        return self._order[a] >= 3

    def isotropic_plane(self, a):
        """Is the plane's quadratic form proportional to the identity on the basis?"""
        G = self.gram[a]
        return G[0][0] == G[1][1]

    # ---------------- fixed points and orbits
    def fixed_points(self, a, k):
        """Fix(theta^k) in plane a, or None for a fixed torus / non-isolated case."""
        if self.frac_twist(a, k) == 0:
            return None                                   # fixed torus
        A = one_minus(self.rot(a, k))
        Ai = minv(A)
        if Ai is None:
            return None
        d = int(abs(A[0][0] * A[1][1] - A[0][1] * A[1][0]))
        out = set()
        for n1 in range(-d - 2, d + 3):
            for n2 in range(-d - 2, d + 3):
                out.add(modl(ap(Ai, (F(n1), F(n2)))))
        return sorted(out)

    def fp_of(self, a, k, n):
        A = one_minus(self.rot(a, k))
        Ai = minv(A)
        if Ai is None:
            return None
        return modl(ap(Ai, n))

    def orbits(self, a, k):
        """theta-orbits (under the FULL theta) of Fix(theta^k) in plane a."""
        fps = self.fixed_points(a, k)
        if fps is None:
            return None
        T, seen, out = self.theta[a], set(), []
        for f in fps:
            if f in seen:
                continue
            orb, cur = [], f
            for _ in range(self._order[a]):
                orb.append(cur)
                cur = modl(ap(T, cur))
                if cur == f:
                    break
            orb = sorted(set(orb))
            seen.update(orb)
            out.append(orb)
        return out

    def sim_orbits(self, k):
        """SIMULTANEOUS theta-orbits of the full 6d fixed-point set of sector k.

        A twisted state is not labelled by a fixed point per plane independently: it is
        labelled by the theta-orbit of its constructing element, and theta acts on all
        three planes AT ONCE (premise_note sec.3; gate_review sec.1.1).  The number of
        such orbits is therefore the number of physical states' geographic labels, and it
        is exactly what orbifolder prints as `constructing elements` for that sector --
        which makes `#sim_orbits == #constructing elements` a sharp per-geometry
        invariant.

        Returns a list of dicts: {members, class_tuple} where `class_tuple[a]` is the
        index of the plane-a theta-orbit the state sits in (None for a fixed torus).
        The class tuple is CONSTANT along a simultaneous orbit, because theta maps each
        plane's point inside its own plane-a orbit.
        """
        per_plane, orb_lists = [], []
        for a in range(3):
            fps = self.fixed_points(a, k)
            orbs = self.orbits(a, k)
            per_plane.append(fps)             # None for a fixed torus
            orb_lists.append(orbs)

        # enumerate the product of the isolated planes
        axes = [fps if fps is not None else [None] for fps in per_plane]
        pts = []
        for p0 in axes[0]:
            for p1 in axes[1]:
                for p2 in axes[2]:
                    pts.append((p0, p1, p2))

        def step(p):
            out = []
            for a in range(3):
                out.append(None if p[a] is None
                           else modl(ap(self.theta[a], p[a])))
            return tuple(out)

        seen, out = set(), []
        for p in pts:
            if p in seen:
                continue
            orb, cur = [], p
            for _ in range(max(self._order)):
                orb.append(cur)
                cur = step(cur)
                if cur == p:
                    break
            orb = sorted(set(orb), key=repr)
            seen.update(orb)
            ct = []
            for a in range(3):
                if orb[0][a] is None:
                    ct.append(None)
                else:
                    ct.append(next(i for i, o in enumerate(orb_lists[a])
                                   if orb[0][a] in o))
            # the class tuple must be constant along the orbit -- verified, not assumed
            for m in orb:
                for a in range(3):
                    if m[a] is None:
                        continue
                    assert m[a] in orb_lists[a][ct[a]], (k, a, m)
            out.append({"members": orb, "class_tuple": tuple(ct),
                        "size": len(orb)})
        return out

    # ---------------- the per-geometry invariant check
    def invariant_check(self):
        """Cross-check the derived fixed-point structure against orbifolder's OWN
        printed `constructing elements` block.  Returns (ok, report)."""
        rep = {"per_sector": [], "ok": True, "reasons": []}
        sectors = self.raw["sectors"]
        if not sectors:
            rep["ok"] = False
            rep["reasons"].append("no constructing-elements block in file")
            return False, rep

        for si, elems in enumerate(sectors):
            ks = set(e[0] for e in elems)
            if len(ks) != 1:
                rep["reasons"].append("sector %d mixes twist powers %s" % (si, sorted(ks)))
                rep["ok"] = False
                continue
            k = ks.pop()
            if any(e[1] != 0 for e in elems):
                rep["reasons"].append("sector %d has a nonzero second twist power" % si)
                rep["ok"] = False
                continue
            row = {"sector": "U" if k == 0 else "T%d" % k, "k": k,
                   "n_constructing_elements": len(elems)}
            if k == 0:
                rep["per_sector"].append(row)
                continue

            # TWO independent conditions, both required:
            #  (1) every listed element's fixed point lies in the DERIVED Fix(theta^k);
            #  (2) the listed elements touch exactly one representative of every derived
            #      theta-ORBIT.  orbifolder prints one constructing element per space-group
            #      CONJUGACY CLASS, and conjugation by theta moves a fixed point along its
            #      theta-orbit -- so (2) is a direct, per-geometry check of the certified
            #      theta-orbit convention (premise_note sec.3) against orbifolder's own output.
            covered = [set() for _ in range(3)]
            bad = 0
            for (_, _, n) in elems:
                for a in range(3):
                    if self.frac_twist(a, k) == 0:
                        continue                          # fixed torus: no constraint
                    i, j = self.pairs[a]
                    lam = (n[i], n[j])
                    f = self.fp_of(a, k, lam)
                    fps = self.fixed_points(a, k)
                    if f is None or f not in fps:
                        bad += 1
                    else:
                        covered[a].add(f)
            row["elements_with_a_fixed_point_off_the_derived_set"] = bad
            n_fix, n_orb, n_orb_touched, per_plane = [], [], [], []
            for a in range(3):
                fps = self.fixed_points(a, k)
                if fps is None:
                    per_plane.append("TORUS")
                    continue
                orbs = self.orbits(a, k)
                per_plane.append(len(fps))
                n_fix.append(len(fps))
                n_orb.append(len(orbs))
                touched = set()
                for p in covered[a]:
                    for oi, o in enumerate(orbs):
                        if p in o:
                            touched.add(oi)
                n_orb_touched.append(len(touched))
            row["derived_nFix_per_plane"] = per_plane
            row["derived_nOrbits_per_plane"] = n_orb
            row["orbits_touched_by_file_elements"] = n_orb_touched

            # THE sharp invariant: orbifolder prints exactly one constructing element
            # per space-group conjugacy class = per SIMULTANEOUS theta-orbit.
            nso = len(self.sim_orbits(k))
            row["derived_n_simultaneous_orbits"] = nso
            if nso != len(elems):
                rep["ok"] = False
                rep["reasons"].append(
                    "sector T%d: derived %d simultaneous theta-orbits but the file "
                    "prints %d constructing elements" % (k, nso, len(elems)))

            if bad:
                rep["ok"] = False
                rep["reasons"].append(
                    "sector T%d: %d listed elements have a fixed point outside the "
                    "derived Fix(theta^k)" % (k, bad))
            if n_orb != n_orb_touched:
                rep["ok"] = False
                rep["reasons"].append(
                    "sector T%d: derived theta-orbit counts %s not matched by the file's "
                    "own constructing elements %s" % (k, n_orb, n_orb_touched))
            rep["per_sector"].append(row)
        return rep["ok"], rep

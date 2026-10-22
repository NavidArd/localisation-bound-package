"""ALPHABET SURVEY PHASE 2 -- Z_M x Z_N orbifold geometry engine.

Extension of the phase-1 certified engine `alphabet_survey/code/survey_geom.py`
(gate: `alphabet_survey/results/gate_review.md`, VERDICT PASS -- A-1 CERTIFIED) from
factorizable Z_N geometries to the 24 factorizable Z_M x Z_N geometries that phase 1
recorded as blocked ONLY by the 1301.2322 footnote-1 scope boundary
(`results/findings.md` [D-3]: the 24-file list; the other 95 Z_M x Z_N files are ALSO
non-factorizable and stay NOT-COMPUTABLE-YET).

EVERY extension over the phase-1 engine is stated at equation level in
`results/phase2_findings.md` (E1-E8).  The ones this module carries:

  E1  sector label g = (k, l) in Z_M x Z_N; plane-a twist phase
          phi_a(g) = k v^1_a + l v^2_a  mod 1
      replacing the single-generator k v_a mod 1.  v^1, v^2 are the two twist rows
      the geometry file itself declares (never a paper heading -- the 1301.2322
      Z6-I lattice-heading trap, findings phase-1 Rule-6 refinement 3).
  E5  roto-translations: the file's `begin twist space group generators` block may
      attach FRACTIONAL lattice shifts s1, s2 to the two twist generators
      (16 of the 24 files do).  The translation coset of sector (k, l) is
          sigma(k,l) = sum_{i<k} theta^i s1  +  theta^k sum_{j<l} omega^j s2
      (space-group composition law (h,s)(g,t) = (hg, s + h t), applied per plane).
  E6  conjugacy classes: conjugating (g, lambda) by (h, s_h) gives
      (g, h lambda + (1-g) s_h), i.e. the fixed point maps AFFINELY
          f  |->  h f + s_h  (mod Lambda).
      Orbits are therefore taken under the affine action of the FULL two-generator
      point group, replacing phase 1's cyclic theta-orbits (to which this reduces
      exactly at s1 = s2 = 0, single generator).
  E7  fixed points with shifts:  F_a(g) = { (1 - R_a(g))^{-1} (n + sigma_a(g)) mod
      Lambda_a }.  A plane with phi_a(g) = 0 is a fixed torus; if additionally
      sigma_a(g) is fractional the sector acts FREELY in that plane (shifted torus)
      -- either way the leg is delocalized there and no coset vector exists, exactly
      the certified 'capable but a leg is DELOCALIZED' status.

CONVENTION DECLARED ONCE (CLAUDE.md 5.19(i)): the file's first twist row is the Z_M
generator theta (sector label k), the second the Z_N generator omega (label l); each
acts in plane a as the rotation by 2 pi phi_a.  The invariant that checks the
convention BEFORE any count: orbifolder's own printed `constructing elements` -- the
first two tokens of every element are (k, l) in exactly this order, and the derived
per-sector torus/fixed-point pattern and simultaneous-orbit count must reproduce the
file's own block sector by sector (`invariant_check`).  A geometry failing it is
FLAGGED and gets no alphabet (phase-1 scope rule iii).
"""

import os
import sys
from fractions import Fraction as F

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")))
from survey_geom import (I2, mm, mpow, ap, minv, one_minus, modl,          # noqa: E402
                         plane_pairs, gram, rot_in_plane, NotComputable)


def parse_geometry_mn(path):
    """Parse an orbifolder Geometry_*.txt including the Z_M x Z_N specifics.

    Superset of phase 1's `parse_geometry`: additionally reads the
    `begin twist space group generators` block (roto-translation shifts) and keeps
    BOTH constructing-element twist powers (k, l) instead of asserting l == 0.
    Nothing is interpreted here.
    """
    st = [l.strip() for l in open(path)]
    out = {}

    i = st.index("point group")
    out["M"] = int(st[i + 1])
    out["N"] = int(st[i + 2])

    i = st.index("ZMxZN")
    out["ZMxZN"] = (st[i + 1].lower() == "true")

    i = st.index("lattice label")
    out["lattice_label"] = st[i + 1]

    i = st.index("begin twist")
    j = st.index("end twist")
    twists = []
    for l in st[i + 1:j]:
        if not l:
            continue
        toks = l.split()
        v = []
        for tok in toks[1:]:
            num, den = tok.split("/")
            v.append(F(int(num), int(den)))
        assert len(v) == 3, (path, l)
        twists.append(tuple(v))
    out["twists"] = twists

    # ---- roto-translation shifts attached to the two twist generators (E5)
    i = st.index("begin twist space group generators")
    j = st.index("end twist space group generators")
    gens = []
    for l in st[i + 1:j]:
        if not l:
            continue
        toks = [F(x) for x in l.split()]
        assert len(toks) == 8, (path, l)
        gens.append({"kl": (int(toks[0]), int(toks[1])), "shift": tuple(toks[2:])})
    out["twist_generators"] = gens

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
                toks = [F(x) for x in l.split()]
                assert len(toks) == 8, (path, l)
                cur.append((int(toks[0]), int(toks[1]), tuple(toks[2:])))
        if cur is not None:
            out["sectors"].append(cur)
    return out


class GeomMN(object):
    """A factorizable Z_M x Z_N orbifold geometry, everything derived from the file."""

    def __init__(self, path):
        self.path = path
        self.name = os.path.basename(path)[len("Geometry_"):-len(".txt")]
        raw = parse_geometry_mn(path)
        self.raw = raw
        if not raw["ZMxZN"]:
            raise NotComputable("phase2-scope", "not a Z_M x Z_N geometry")
        self.M, self.N = raw["M"], raw["N"]

        self.vecs = raw["vecs"]
        self.pairs = plane_pairs(self.vecs)
        if self.pairs is None:
            cnt = plane_pairs(self.vecs, diagnostic=True)
            raise NotComputable(
                "NOT-COMPUTABLE-YET/non-factorizable",
                "in-plane lattice-vector counts %s, need [2,2,2]" % (cnt,))

        self.gram = [gram(self.vecs, p) for p in self.pairs]
        if any(g is None for g in self.gram):
            raise NotComputable("NOT-COMPUTABLE-YET/irrational-gram", "")

        # the two twist vectors, in the file's own order (convention block above)
        assert len(raw["twists"]) == 2, path
        self.v1, self.v2 = raw["twists"]

        # convention check: the file's twist-generator block must label the first
        # generator (1,0) and the second (0,1) -- parsed, not assumed
        gens = raw["twist_generators"]
        assert len(gens) == 2 and gens[0]["kl"] == (1, 0) and gens[1]["kl"] == (0, 1), \
            (path, [g["kl"] for g in gens])
        # per-plane components of the two generator shifts, in the plane's own basis
        self.s1 = [self._plane_comp(gens[0]["shift"], a) for a in range(3)]
        self.s2 = [self._plane_comp(gens[1]["shift"], a) for a in range(3)]
        self.has_rototranslation = any(x != (F(0), F(0)) for x in self.s1 + self.s2)

        # per-plane rotations of the two generators, derived from real space (E1)
        self.theta, self.omega = [], []
        for a in range(3):
            Ra = rot_in_plane(self.vecs, self.pairs[a], a, self.v1[a])
            Rb = rot_in_plane(self.vecs, self.pairs[a], a, self.v2[a])
            if Ra is None or Rb is None:
                raise NotComputable("NOT-COMPUTABLE-YET/theta-not-integral",
                                    "plane %d" % (a + 1))
            self.theta.append(Ra)
            self.omega.append(Rb)

        # both generators must be exact isometries of every plane form (5.19(i))
        for a in range(3):
            for Rg in (self.theta[a], self.omega[a]):
                for e in ((F(1), F(0)), (F(0), F(1)), (F(1), F(1))):
                    if self.rho(a, ap(Rg, e)) != self.rho(a, e):
                        raise NotComputable("NOT-COMPUTABLE-YET/theta-not-isometry",
                                            "plane %d" % (a + 1))

        # commutativity in every plane (abelian point group), verified not assumed
        for a in range(3):
            assert mm(self.theta[a], self.omega[a]) == mm(self.omega[a], self.theta[a]), path

        self.sectors_all = [(k, l) for k in range(self.M) for l in range(self.N)]
        self._plane_order = [self._ord2(a) for a in range(3)]
        self._sim_cache = {}

    # ---------------- basics
    def _plane_comp(self, shift6, a):
        i, j = self.pairs[a]
        return (shift6[i], shift6[j])

    def _ord2(self, a):
        """Order of the point-group image in plane a (cyclic, generated by both)."""
        seen, frontier = {I2}, [I2]
        while frontier:
            cur = frontier.pop()
            for Rg in (self.theta[a], self.omega[a]):
                nx = mm(cur, Rg)
                if nx not in seen:
                    seen.add(nx)
                    frontier.append(nx)
        return len(seen)

    def plane_order(self, a):
        return self._plane_order[a]

    def rho(self, a, u):
        G = self.gram[a]
        return (G[0][0] * u[0] * u[0] + 2 * G[0][1] * u[0] * u[1]
                + G[1][1] * u[1] * u[1])

    def phi(self, a, g):
        """Plane-a twist phase of sector g = (k, l):  k v1_a + l v2_a mod 1  (E1)."""
        k, l = g
        x = k * self.v1[a] + l * self.v2[a]
        return x - (x // 1)

    def rot(self, a, g):
        """Plane-a rotation of sector g:  theta_a^k omega_a^l  (k < M, l < N: direct)."""
        k, l = g
        return mm(mpow(self.theta[a], k), mpow(self.omega[a], l))

    def sigma(self, a, g):
        """Plane-a translation coset offset of sector g (E5):
        sigma(k,l) = sum_{i<k} theta^i s1 + theta^k sum_{j<l} omega^j s2."""
        k, l = g
        out = (F(0), F(0))
        for i in range(k):
            t = ap(mpow(self.theta[a], i), self.s1[a])
            out = (out[0] + t[0], out[1] + t[1])
        Tk = mpow(self.theta[a], k)
        for j in range(l):
            t = ap(Tk, ap(mpow(self.omega[a], j), self.s2[a]))
            out = (out[0] + t[0], out[1] + t[1])
        return out

    def gen_affine(self, a):
        """The two generators' affine actions in plane a: (R, s) pairs (E6)."""
        return ((self.theta[a], self.s1[a]), (self.omega[a], self.s2[a]))

    def rigid_plane(self, a):
        """Plane shape fixed by the point group iff its image has order >= 3."""
        return self._plane_order[a] >= 3

    # ---------------- fixed points and orbits
    def freely_acting(self, a, g):
        """phi_a(g) = 0 but sigma_a(g) fractional: pure shifted translation (E7)."""
        if self.phi(a, g) != 0:
            return False
        s = self.sigma(a, g)
        return not (s[0].denominator == 1 and s[1].denominator == 1)

    def fixed_points(self, a, g):
        """F_a(g), or None for a fixed torus / freely-acting plane (E7)."""
        if self.phi(a, g) == 0:
            return None
        A = one_minus(self.rot(a, g))
        Ai = minv(A)
        if Ai is None:
            return None
        d = int(abs(A[0][0] * A[1][1] - A[0][1] * A[1][0]))
        s = self.sigma(a, g)
        out = set()
        for n1 in range(-d - 2, d + 3):
            for n2 in range(-d - 2, d + 3):
                out.add(modl(ap(Ai, (F(n1) + s[0], F(n2) + s[1]))))
        out = sorted(out)
        assert len(out) == d, (self.name, a, g, len(out), d)
        return out

    def fp_of(self, a, g, lam):
        A = one_minus(self.rot(a, g))
        Ai = minv(A)
        if Ai is None:
            return None
        return modl(ap(Ai, lam))

    def orbits(self, a, g):
        """Orbits of F_a(g) under the affine action of the FULL point group (E6)."""
        fps = self.fixed_points(a, g)
        if fps is None:
            return None
        gens = self.gen_affine(a)
        seen, out = set(), []
        for f in fps:
            if f in seen:
                continue
            orb, frontier = {f}, [f]
            while frontier:
                cur = frontier.pop()
                for (Rg, s) in gens:
                    t = ap(Rg, cur)
                    nx = modl((t[0] + s[0], t[1] + s[1]))
                    if nx not in orb:
                        orb.add(nx)
                        frontier.append(nx)
            orb = sorted(orb)
            assert all(x in fps for x in orb), (self.name, a, g)
            seen.update(orb)
            out.append(orb)
        return out

    def sim_orbits(self, g):
        """Simultaneous affine point-group orbits of the 6d fixed-point set of
        sector g -- the physical geographic labels; #sim_orbits must equal the
        number of constructing elements the file prints for that sector."""
        if g in self._sim_cache:
            return self._sim_cache[g]
        per_plane, orb_lists = [], []
        for a in range(3):
            per_plane.append(self.fixed_points(a, g))
            orb_lists.append(self.orbits(a, g))

        axes = [fps if fps is not None else [None] for fps in per_plane]
        pts = [(p0, p1, p2) for p0 in axes[0] for p1 in axes[1] for p2 in axes[2]]

        def images(p):
            out = []
            for gi in range(2):
                q = []
                for a in range(3):
                    if p[a] is None:
                        q.append(None)
                    else:
                        (Rg, s) = self.gen_affine(a)[gi]
                        t = ap(Rg, p[a])
                        q.append(modl((t[0] + s[0], t[1] + s[1])))
                out.append(tuple(q))
            return out

        seen, out = set(), []
        for p in pts:
            if p in seen:
                continue
            orb, frontier = {p}, [p]
            while frontier:
                cur = frontier.pop()
                for nx in images(cur):
                    if nx not in orb:
                        orb.add(nx)
                        frontier.append(nx)
            orb = sorted(orb, key=repr)
            seen.update(orb)
            ct = []
            for a in range(3):
                if orb[0][a] is None:
                    ct.append(None)
                else:
                    ct.append(next(i for i, o in enumerate(orb_lists[a])
                                   if orb[0][a] in o))
            # class tuple constant along the orbit -- verified, not assumed
            for m in orb:
                for a in range(3):
                    if m[a] is None:
                        continue
                    assert m[a] in orb_lists[a][ct[a]], (g, a, m)
            out.append({"members": orb, "class_tuple": tuple(ct), "size": len(orb)})
        self._sim_cache[g] = out
        return out

    # ---------------- the per-geometry invariant check
    def invariant_check(self):
        """Cross-check the derived structure against orbifolder's OWN printed
        `constructing elements`, sector by sector.  Same two conditions as the
        certified phase-1 check: (1) every listed element's fixed point lies in the
        derived F_a(g); (2) #derived simultaneous orbits == #printed elements, and
        the printed elements touch every derived per-plane orbit."""
        rep = {"per_sector": [], "ok": True, "reasons": []}
        sectors = self.raw["sectors"]
        if not sectors:
            rep["ok"] = False
            rep["reasons"].append("no constructing-elements block in file")
            return False, rep

        for si, elems in enumerate(sectors):
            kls = set((e[0], e[1]) for e in elems)
            if len(kls) != 1:
                rep["reasons"].append("sector %d mixes twist labels %s" % (si, sorted(kls)))
                rep["ok"] = False
                continue
            g = kls.pop()
            row = {"sector": "U" if g == (0, 0) else "T(%d,%d)" % g,
                   "g": list(g), "n_constructing_elements": len(elems)}
            if g == (0, 0):
                rep["per_sector"].append(row)
                continue

            covered = [set() for _ in range(3)]
            bad = 0
            for (_, _, n) in elems:
                for a in range(3):
                    if self.phi(a, g) == 0:
                        continue
                    i, j = self.pairs[a]
                    lam = (n[i], n[j])
                    f = self.fp_of(a, g, lam)
                    fps = self.fixed_points(a, g)
                    if f is None or f not in fps:
                        bad += 1
                    else:
                        covered[a].add(f)
            row["elements_with_a_fixed_point_off_the_derived_set"] = bad
            n_orb, n_orb_touched, per_plane = [], [], []
            for a in range(3):
                fps = self.fixed_points(a, g)
                if fps is None:
                    per_plane.append("TORUS(free)" if self.freely_acting(a, g) else "TORUS")
                    continue
                orbs = self.orbits(a, g)
                per_plane.append(len(fps))
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

            nso = len(self.sim_orbits(g))
            row["derived_n_simultaneous_orbits"] = nso
            if nso != len(elems):
                rep["ok"] = False
                rep["reasons"].append(
                    "sector T(%d,%d): derived %d simultaneous orbits but the file "
                    "prints %d constructing elements" % (g[0], g[1], nso, len(elems)))
            if bad:
                rep["ok"] = False
                rep["reasons"].append(
                    "sector T(%d,%d): %d listed elements have a fixed point outside "
                    "the derived set" % (g[0], g[1], bad))
            if n_orb != n_orb_touched:
                rep["ok"] = False
                rep["reasons"].append(
                    "sector T(%d,%d): derived per-plane orbit counts %s not matched "
                    "by the file's own elements %s" % (g[0], g[1], n_orb, n_orb_touched))
            rep["per_sector"].append(row)
        return rep["ok"], rep

    # ---------------- exchange symmetry (quality axis, reported only)
    def automorphisms(self):
        """All automorphisms of Z_M x Z_N as generator-image pairs (a, b)."""
        G = self.sectors_all
        out = []
        for ga in G:
            for gb in G:
                imgs = set()
                for (k, l) in G:
                    imgs.add(((k * ga[0] + l * gb[0]) % self.M,
                              (k * ga[1] + l * gb[1]) % self.N))
                if len(imgs) == len(G):
                    out.append((ga, gb))
        return out

    def exchange_symmetric_planes(self):
        """Unordered plane pairs (p, q) whose exchange extends to a point-group
        automorphism preserving the twist assignment (twist+lattice level; a
        symmetry-point equal-eps locus exists on such a pair).  Roto-translation
        shifts are NOT tested here -- recorded as a caveat where they exist."""
        auts = self.automorphisms()
        out = []
        for p in range(3):
            for q in range(p + 1, 3):
                if self.gram[p] != self.gram[q]:
                    continue
                r = 3 - p - q
                found = None
                for (ga, gb) in auts:
                    def tau(g, ga=ga, gb=gb):
                        return ((g[0] * ga[0] + g[1] * gb[0]) % self.M,
                                (g[0] * ga[1] + g[1] * gb[1]) % self.N)
                    if all(self.phi(p, tau(g)) == self.phi(q, g)
                           and self.phi(q, tau(g)) == self.phi(p, g)
                           and self.phi(r, tau(g)) == self.phi(r, g)
                           for g in self.sectors_all):
                        found = (ga, gb)
                        break
                if found:
                    out.append({"planes": [p + 1, q + 1],
                                "generator_images": [list(found[0]), list(found[1])]})
        return out

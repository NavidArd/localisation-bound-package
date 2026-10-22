"""multi_geometry_campaign -- exact-rational geometry library for the classified
heterotic orbifolds shipped with orbifolder 1.2.1, specialised to Z6-II.

Everything here is DERIVED from the primary data on disk:
    orbifolder-1.2.1/Geometry/Geometry_<label>.txt      (the space-group definition)
    orbifolder-1.2.1/All138SpaceGroups.txt              (the trivial-model templates)
No number in this module is quoted from the literature.

CONVENTIONS, STATED ONCE (CLAUDE.md 5.19(i)) -- and checked against invariants in
`selftest()` BEFORE any count is reported:

  * The geometry file's `begin twist` block carries FOUR components; component 0 is
    the 4d helicity slot and is always 0/1.  The internal twist is components 1..3:
        v = (v_1, v_2, v_3).
    For every Z6-II lattice orbifolder ships v = (1/6, 1/3, -1/2).
  * Local twist of the k-th sector: v_k = k*v, and its plane-a "fractional twist" is
        f_a(k) = frac(k*v_a)  in  [0,1).
    Plane a is UNROTATED in sector k iff f_a(k) == 0; the sector then has a FIXED
    TORUS in plane a (an N=2 subsector if exactly one plane is unrotated, N=4 if all
    three are, i.e. the untwisted sector).
  * The right-moving shifted momentum of a TWISTED left-chiral multiplet, in
    orbifolder's storage convention, is
        q_sh_a = f_a(k) - 1        when f_a(k) != 0
    and, in an unrotated plane, q_sh_a in {0, -1} -- the two values being the two
    N=2 "chiralities".  This is the Z6-II generalisation of the rcharge_census gate
    G4 identity, and it is re-verified here against orbifolder's own printed
    q_sh for every field of every model actually run (see code/verify_qsh.py).
  * The literature H-momentum is r_a = -q_sh_a.  Orbifolder's discrete R-charge is
        R^orb_a = q_sh_a + Osci_a ,   Osci_a = Nbar_a - N_a ,
    and R^orb = -R^lit exactly (rcharge_census/docs/rule_derivation.md sec.0).
  * A space-group element is stored by orbifolder as 8 integers
        (k, l, n_1, ..., n_6)
    meaning g = (theta^k omega^l, sum_i n_i e_i).  For Z_N (ZMxZN false) l is always 0.

  * PLANE ROTATION ORDER  N_a := order of the point-group rotation acting on plane a
        N_a = denominator(v_a)      (v_a in lowest terms)
    For Z6-II: N = (6, 3, 2).  Plane a is called NON-PRIME iff N_a is not prime.
    Z6-II therefore has EXACTLY ONE non-prime plane (plane 1, order 6).

Nothing in this module is Z6-II-specific by construction; the Z6-II specialisations
are asserted, not assumed, in selftest().
"""

import os
import re
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
PROG = os.path.abspath(os.path.join(ROOT, ".."))
GEOMDIR = os.path.join(PROG, "orbifolder-1.2.1", "Geometry")
ALL138 = os.path.join(PROG, "orbifolder-1.2.1", "All138SpaceGroups.txt")


# --------------------------------------------------------------------------- parsing

def _rat(tok):
    tok = tok.strip()
    if "/" in tok:
        n, d = tok.split("/")
        return F(int(n), int(d))
    return F(tok)


def _block(text, name):
    m = re.search(r"^begin %s\s*$(.*?)^end %s\s*$" % (re.escape(name), re.escape(name)),
                  text, re.S | re.M)
    return m.group(1) if m else None


class Geometry(object):
    """One orbifolder geometry file, parsed exactly."""

    def __init__(self, label):
        self.label = label
        self.path = os.path.join(GEOMDIR, "Geometry_%s.txt" % label)
        with open(self.path) as fh:
            self.text = fh.read()
        t = self.text

        m = re.search(r"^point group\s*$\s*(\d+)\s*$\s*(\d+)\s*$", t, re.M)
        self.N, self.M = int(m.group(1)), int(m.group(2))
        self.is_ZMxZN = re.search(r"^ZMxZN\s*$\s*(\w+)", t, re.M).group(1) == "true"
        self.lattice_label = re.search(r"^lattice label\s*$\s*(\S+)", t, re.M).group(1)

        tw = [_rat(x) for x in _block(t, "twist").split()]
        # 4 components for Z_N, 8 for Z_NxZ_M
        if len(tw) == 4:
            self.twists = [tuple(tw[1:4])]
        elif len(tw) == 8:
            self.twists = [tuple(tw[1:4]), tuple(tw[5:8])]
        else:
            raise ValueError("%s: unexpected twist block length %d" % (label, len(tw)))

        self.sg_twist_gens = [tuple(int(x) for x in ln.split())
                              for ln in _block(t, "twist space group generators").strip().splitlines()]
        self.sg_shift_gens = [tuple(int(x) for x in ln.split())
                              for ln in _block(t, "shift space group generators").strip().splitlines()]

        # allowed Wilson-line orders (6 numbers, one per lattice vector)
        m = re.search(r"^allowed order of Wilson lines\s*$((?:\s*\d+\s*$){6})", t, re.M)
        self.wl_orders = [int(x) for x in m.group(1).split()]
        idw = _block(t, "identical WLs")
        # NOTE: the block may carry SEVERAL lines, each one an independent set of
        # identified Wilson lines (e.g. Z6-II_3_1 has "1 2 3" and "4 5").  Flattening
        # them into one list is WRONG and was caught by gate S3 of parse_parr_z6ii.py.
        self.identical_wl_groups = ([[int(x) for x in ln.split()]
                                     for ln in idw.strip().splitlines() if ln.strip()]
                                    if idw else [])
        self.identical_wls = ([g[0] for g in self.identical_wl_groups]
                              if False else sorted({i for gr in self.identical_wl_groups
                                                    for i in gr}))

        # lattice vectors (floats -- used only for reporting the plane type)
        lat = _block(t, "lattice")
        vecs, cur = [], []
        for ln in lat.strip().splitlines():
            ln = ln.strip()
            if ln.startswith("next lattice vector"):
                if cur:
                    vecs.append([float(x) for x in cur])
                cur = []
            elif ln:
                cur.append(ln)
        if cur:
            vecs.append([float(x) for x in cur])
        self.lattice_vectors = vecs

        # constructing elements, grouped by sector, in file order
        ce = _block(t, "constructing elements")
        self.sectors_ce = []
        for chunk in ce.split("new sector"):
            chunk = chunk.strip()
            if not chunk:
                continue
            self.sectors_ce.append([tuple(int(x) for x in ln.split())
                                    for ln in chunk.splitlines() if ln.strip()])

        # declared discrete symmetries
        self.nonR, self.Rsym = _parse_symmetries(_block(t, "discrete symmetries"))

    # ---------------------------------------------------------------- derived

    @property
    def v(self):
        """The single internal twist vector (Z_N geometries only)."""
        if len(self.twists) != 1:
            raise ValueError("%s is Z_MxZ_N; use .twists" % self.label)
        return self.twists[0]

    def local_twist(self, k, l=0):
        if len(self.twists) == 1:
            return tuple(k * a for a in self.twists[0])
        return tuple(k * a + l * b for a, b in zip(*self.twists))

    def frac_twist(self, k, l=0):
        return tuple(a - (a // 1) for a in self.local_twist(k, l))   # frac in [0,1)

    def unrotated_planes(self, k, l=0):
        return [a for a, f in enumerate(self.frac_twist(k, l)) if f == 0]

    def susy_of_sector(self, k, l=0):
        """N=1 (no unrotated plane), N=2 (one), N=4 (three, i.e. untwisted)."""
        n = len(self.unrotated_planes(k, l))
        return {0: 1, 1: 2, 3: 4}[n]

    @property
    def plane_orders(self):
        """N_a = order of the point-group rotation acting on plane a."""
        if len(self.twists) == 1:
            return tuple(a.denominator for a in self.twists[0])
        # Z_MxZ_N: order of the cyclic group generated by all (k v_theta + l v_omega)_a
        out = []
        for a in range(3):
            ords = set()
            for k in range(self.N):
                for l in range(self.M):
                    x = k * self.twists[0][a] + l * self.twists[1][a]
                    ords.add((x - (x // 1)).denominator)
            out.append(max(ords))
        return tuple(out)

    @property
    def plane_factorization(self):
        """Which complex planes the LATTICE factorizes on.

        Returns (factorized_planes, mixed_blocks).  A plane a is `factorized` iff
        every lattice vector with support in plane a has support in NO other plane,
        and plane a is spanned by exactly two of them.  Planes that are not
        factorized are grouped into connected blocks.  This is the criterion that
        decides which row of arXiv:1301.2322 Table 3 applies to the geometry
        (factorizable / partially factorizable / non-factorizable).
        """
        sup = []
        for v in self.lattice_vectors:
            sup.append(frozenset(j // 2 for j, x in enumerate(v) if abs(x) > 1e-9))
        # union-find over planes connected by a shared lattice vector
        parent = {a: a for a in range(3)}

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        for s_ in sup:
            ss = sorted(s_)
            for b in ss[1:]:
                ra, rb = find(ss[0]), find(b)
                if ra != rb:
                    parent[rb] = ra
        blocks = {}
        for a in range(3):
            blocks.setdefault(find(a), []).append(a)
        blocks = [sorted(v) for v in blocks.values()]
        fac = [b[0] for b in blocks if len(b) == 1]
        mixed = [b for b in blocks if len(b) > 1]
        return sorted(fac), sorted(mixed)

    @property
    def factorizability(self):
        f, m = self.plane_factorization
        if not m:
            return "factorizable"
        if not f:
            return "non-factorizable"
        return "partially factorizable"

    @property
    def nonprime_planes(self):
        return [a for a, n in enumerate(self.plane_orders) if not _is_prime(n)]

    def sector_summary(self):
        rows = []
        for k, ce in enumerate(self.sectors_ce):
            rows.append(dict(k=k,
                             n_constructing_elements=len(ce),
                             local_twist=[str(x) for x in self.local_twist(k)],
                             frac_twist=[str(x) for x in self.frac_twist(k)],
                             unrotated_planes=[a + 1 for a in self.unrotated_planes(k)],
                             susy=self.susy_of_sector(k),
                             fixed_point_indices_per_plane=_fp_spread(ce)))
        return rows


def _is_prime(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def _fp_spread(ce):
    """How many distinct values each of the six lattice-translation slots takes."""
    return [len({e[2 + i] for e in ce}) for i in range(6)]


def _parse_symmetries(block):
    """Return (list of non-R, list of R) symmetry records."""
    if block is None:
        return [], []
    nonR, Rs = [], []
    lines = [ln.strip() for ln in block.strip().splitlines() if ln.strip()]
    i = 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith("new non-R symmetry") or ln.startswith("new R symmetry"):
            isR = ln.startswith("new R symmetry")
            order = int(lines[i + 1])
            wcharge = _rat(lines[i + 2])
            ops, j = [], i + 3
            while j < len(lines) and not lines[j].startswith("new "):
                ops.append([_rat(x) for x in lines[j].split()])
                j += 1
            rec = dict(order=order, superpotential_charge=str(wcharge),
                       operators=[[str(x) for x in o] for o in ops])
            (Rs if isR else nonR).append(rec)
            i = j
        else:
            i += 1
    return nonR, Rs


# --------------------------------------------------------------------------- selftest

def selftest(verbose=True):
    """Invariant checks that must pass before any Z6-II count is believed."""
    out = []

    def chk(name, cond, detail=""):
        out.append((name, bool(cond), detail))
        if verbose:
            print("%-46s %s  %s" % (name, "PASS" if cond else "FAIL", detail))
        return cond

    # I1 -- the Z3xZ3 geometry the programme already trusts reproduces its known facts
    g3 = Geometry("Z3xZ3_1_1")
    chk("I1a Z3xZ3 twists",
        g3.twists == [(F(0), F(1, 3), F(-1, 3)), (F(1, 3), F(0), F(-1, 3))],
        str(g3.twists))
    chk("I1b Z3xZ3 plane orders = (3,3,3)", g3.plane_orders == (3, 3, 3), str(g3.plane_orders))
    chk("I1c Z3xZ3 has NO non-prime plane", g3.nonprime_planes == [], str(g3.nonprime_planes))
    chk("I1d Z3xZ3 has 9 sectors", len(g3.sectors_ce) == 9, str(len(g3.sectors_ce)))
    # I1e -- CLAUDE.md 5.15, the programme's own verified ground truth:
    # "6 of the 9 Gamma-elements fix exactly one T^2_a invariantly (two per factor)".
    # This is the invariant that licenses reading fixed-torus structure off this module.
    n2 = [(k, l) for k in range(3) for l in range(3) if g3.susy_of_sector(k, l) == 2]
    per_plane = {}
    for k, l in n2:
        per_plane.setdefault(g3.unrotated_planes(k, l)[0], []).append((k, l))
    chk("I1e Z3xZ3: 6 N=2 sectors, two per plane (CLAUDE 5.15)",
        len(n2) == 6 and sorted(len(v) for v in per_plane.values()) == [2, 2, 2],
        "%s" % ({a + 1: v for a, v in sorted(per_plane.items())},))
    n1 = [(k, l) for k in range(3) for l in range(3) if g3.susy_of_sector(k, l) == 1]
    chk("I1f Z3xZ3: only T(1,1),T(2,2) are N=1", sorted(n1) == [(1, 1), (2, 2)], str(n1))

    # I2 -- Z6-II structure
    for lab in ["Z6-II_1_1", "Z6-II_2_1", "Z6-II_3_1", "Z6-II_4_1"]:
        g = Geometry(lab)
        chk("I2 %s twist = (1/6,1/3,-1/2)" % lab,
            g.v == (F(1, 6), F(1, 3), F(-1, 2)), str(g.v))
        chk("I2 %s plane orders = (6,3,2)" % lab, g.plane_orders == (6, 3, 2),
            str(g.plane_orders))
        chk("I2 %s exactly one non-prime plane (plane 1)" % lab,
            g.nonprime_planes == [0], str(g.nonprime_planes))
        susy = [g.susy_of_sector(k) for k in range(6)]
        chk("I2 %s SUSY per sector = [4,1,2,2,2,1]" % lab,
            susy == [4, 1, 2, 2, 2, 1], str(susy))
        chk("I2 %s zero R symmetries declared" % lab, g.Rsym == [],
            "%d non-R, %d R" % (len(g.nonR), len(g.Rsym)))

    return out


if __name__ == "__main__":
    selftest()

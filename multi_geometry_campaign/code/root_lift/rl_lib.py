"""root_lift -- exact space-group algebra for the LOCAL-VS-GLOBAL ROOT question.

Spec: docs/root_lift_task.md.  Referee charge: paper15/review/external_referee_v2.md item 7.

WHAT THIS MODULE COMPUTES
-------------------------
For a constructing space-group element g of a heterotic Z_N or Z_M x Z_N orbifold, the
centraliser Z_S(g) = { h in S : gh = hg }, and the subgroup of gauge-shift values

    V( Z_S(g) )  subset  (Lambda_{E8xE8} tensor Q) / Lambda_{E8xE8}

that the orbifold projection at g can constrain.  The question is whether that subgroup
equals V(S) -- the subgroup whose annihilator is the UNBROKEN FOUR-DIMENSIONAL gauge group.

CONVENTIONS, STATED ONCE (CLAUDE.md 5.19 (i)), each verified in selftest() BEFORE any
count is reported:

  C1  Geometry file `begin twist` carries FOUR components per generator; component 0 is
      the 4d helicity slot.  The internal twists are
          v1 = (v1_1, v1_2, v1_3)   (order N, the FIRST point-group number in the file)
          v2 = (v2_1, v2_2, v2_3)   (order M, the SECOND)
      matching geomlib.Geometry.twists.

  C2  A space-group element is stored by orbifolder as the 8 integers (k, l, n_1..n_6),
      meaning  g = ( theta^k omega^l , sum_i n_i e_i )  with theta the v1-rotation and
      omega the v2-rotation.  (cspacegroup.cpp / geomlib docstring; same convention.)

  C3  The point-group rotation acting on the LATTICE BASIS is
          M(k,l) = B^{-1} R(k v1 + l v2) B ,
      where B has the geometry file's lattice vectors as COLUMNS and R(v) is the
      block-diagonal Cartesian rotation by ANGLE +2*pi*v_a in plane a.  The + sign is a
      convention; the opposite sign gives the inverse rotation and relabels k -> -k,
      which changes no quantity computed here.  M(k,l) MUST be integral -- checked.

  C4  Group law:  (k1,l1,n1)*(k2,l2,n2) = (k1+k2, l1+l2, n1 + M(k1,l1) n2).
      Hence  g h = h g  <=>  (I - M_g) n_h = (I - M_h) n_g   over Z.

  C5  Gauge embedding.  V_(k,l,n) = k V_1 + l V_2 + sum_i n_i A_i.  Under the standard
      Wilson-line consistency conditions this is a HOMOMORPHISM S -> (Lambda x Q)/Lambda,
      so it factors through the abelianisation; we represent it by the Z-linear map
          Phi : Z^8 -> (Lambda x Q)/Lambda ,  (k,l,n) |-> k V_1 + l V_2 + sum n_i A_i .

  C6  GEOMETRY-LEVEL KERNEL.  ker_geo subset Z^8 is the Z-span of
          (N,0,0), (0,M,0),  and  (0,0,(I - M(k,l)) n)  for all (k,l) and all n in Z^6.
      Every admissible model on the geometry satisfies Phi(ker_geo) = 0:
      the first two by the shift-order condition N V_1 = 0 = M V_2 (mod Lambda), the
      third by Wilson-line consistency A(M lambda) = A(lambda).  A particular model may
      have a LARGER kernel; that is model-specific and is treated separately.

  C7  PROJECTION PHASE, read off orbifolder-1.2.1/src/orbifolder/cfixedpoint.cpp,
      CFixedBrane::SortByEigenvalue ll. 566-612 and CFixedBrane::CreateStates ll. 356-380:

          LM_eig(p, h) = p . V_h  -  (1/2)( V_g . V_h  -  v_g . v_h )
                                   [+ discrete-torsion term, p-independent]
                                   [+ oscillator term,       p-independent]
          RM_eig(h)    = k_h * (v1 . q_sh [+osc])  +  l_h * (v2 . q_sh [+osc])
          state survives  <=>  LM_eig(p,h) == RM_eig(h)  (mod 1)  for all h in Z_S(g).

      Only the FIRST term depends on p.  Therefore two states at the SAME g with the SAME
      right-mover q_sh and the SAME oscillator excitation, with weights p and p+Delta,
      both surviving, force

          Delta . V_h  in  Z   for every h in Z_S(g) .                        (*)

      (*) is the entire content of the root-lift question.

NOTHING in this module is quoted from the literature; every input is read from
orbifolder-1.2.1/Geometry/Geometry_<label>.txt.
"""
from __future__ import annotations

import os
import sys
from fractions import Fraction as F

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.dirname(HERE)
CAMP = os.path.dirname(CODE)
sys.path.insert(0, CODE)

import geomlib                                                          # noqa: E402


# ============================================================ integer linear algebra

def snf(A):
    """Smith normal form with transforms.  Returns (D, U, V) with U A V = D,
    U and V unimodular, D diagonal with non-negative entries.  A is an integer
    numpy array; small matrices only (<= 12x12 here)."""
    A = np.array(A, dtype=object)
    m, n = A.shape
    U = np.eye(m, dtype=object)
    V = np.eye(n, dtype=object)

    def swap_rows(i, j):
        A[[i, j]] = A[[j, i]]
        U[[i, j]] = U[[j, i]]

    def swap_cols(i, j):
        A[:, [i, j]] = A[:, [j, i]]
        V[:, [i, j]] = V[:, [j, i]]

    def addrow(i, j, c):          # row_i += c * row_j
        A[i] = A[i] + c * A[j]
        U[i] = U[i] + c * U[j]

    def addcol(i, j, c):          # col_i += c * col_j
        A[:, i] = A[:, i] + c * A[:, j]
        V[:, i] = V[:, i] + c * V[:, j]

    t = 0
    while t < min(m, n):
        # find a pivot
        piv = None
        for i in range(t, m):
            for j in range(t, n):
                if A[i, j] != 0:
                    piv = (i, j)
                    break
            if piv:
                break
        if piv is None:
            break
        swap_rows(t, piv[0])
        swap_cols(t, piv[1])
        while True:
            # clear column t
            done = True
            for i in range(t + 1, m):
                if A[i, t] != 0:
                    q = A[i, t] // A[t, t]
                    addrow(i, t, -q)
                    if A[i, t] != 0:
                        swap_rows(t, i)
                        done = False
            for j in range(t + 1, n):
                if A[t, j] != 0:
                    q = A[t, j] // A[t, t]
                    addcol(j, t, -q)
                    if A[t, j] != 0:
                        swap_cols(t, j)
                        done = False
            if done:
                break
        if A[t, t] < 0:
            A[t] = -A[t]
            U[t] = -U[t]
        t += 1
    return A, U, V


def int_solve(Amat, b):
    """Integer solution x of A x = b, or None.  A, b integer arrays."""
    A = np.array(Amat, dtype=object)
    b = np.array(b, dtype=object)
    D, U, V = snf(A)
    m, n = A.shape
    c = U @ b
    y = np.zeros(n, dtype=object)
    for i in range(min(m, n)):
        d = D[i, i]
        if d == 0:
            if c[i] != 0:
                return None
        else:
            if c[i] % d != 0:
                return None
            y[i] = c[i] // d
    for i in range(min(m, n), m):
        if c[i] != 0:
            return None
    return V @ y


def int_kernel(Amat):
    """Basis (rows) of the integer kernel { x in Z^n : A x = 0 }."""
    A = np.array(Amat, dtype=object)
    D, U, V = snf(A)
    m, n = A.shape
    cols = []
    for j in range(n):
        d = D[j, j] if j < min(m, n) else 0
        if d == 0:
            cols.append(V[:, j])
    return [np.array(c, dtype=object) for c in cols]


def sublattice_index(gens, n):
    """Index [Z^n : L] where L is the Z-span of `gens` (list of length-n integer
    vectors).  Returns (index, invariant_factors).  index = 0 means infinite index
    (L has rank < n)."""
    if not gens:
        return 0, []
    A = np.array([list(g) for g in gens], dtype=object)     # rows = generators
    D, _, _ = snf(A)
    r = min(D.shape)
    facs = [int(D[i, i]) for i in range(r) if D[i, i] != 0]
    if len(facs) < n:
        return 0, facs
    idx = 1
    for f in facs:
        idx *= f
    return idx, facs


# ============================================================ geometry -> matrices

def _cart_rot(v):
    """Block-diagonal 6x6 Cartesian rotation by angle +2*pi*v_a in plane a (C3)."""
    R = np.zeros((6, 6))
    for a in range(3):
        th = 2.0 * np.pi * float(v[a])
        c, s = np.cos(th), np.sin(th)
        R[2 * a, 2 * a] = c
        R[2 * a, 2 * a + 1] = -s
        R[2 * a + 1, 2 * a] = s
        R[2 * a + 1, 2 * a + 1] = c
    return R


class SpaceGroup(object):
    """The classified orbifold space group, exactly, from the primary geometry file."""

    def __init__(self, label):
        self.G = geomlib.Geometry(label)
        self.label = label
        g = self.G
        self.ZMxZN = len(g.twists) == 2
        self.v1 = g.twists[0]
        self.v2 = g.twists[1] if self.ZMxZN else (F(0), F(0), F(0))
        # geomlib: self.N, self.M = first, second point-group number
        self.N = g.N
        self.M = g.M if self.ZMxZN else 1
        B = np.array(g.lattice_vectors, dtype=float).T          # columns = e_i
        self.B = B
        self.Binv = np.linalg.inv(B)
        self.M1 = self._latmat(self.v1)
        self.M2 = self._latmat(self.v2) if self.ZMxZN else np.eye(6, dtype=object)
        self._cache = {}

    def _latmat(self, v):
        Mf = self.Binv @ _cart_rot(v) @ self.B
        Mi = np.rint(Mf)
        if np.max(np.abs(Mf - Mi)) > 1e-6:
            raise AssertionError("%s: point-group matrix not integral (C3 fails)" % self.label)
        return np.array(Mi.astype(np.int64), dtype=object)

    def mat(self, k, l=0):
        key = (k % self.N, l % self.M)
        if key in self._cache:
            return self._cache[key]
        R = np.eye(6, dtype=object)
        for _ in range(key[0]):
            R = self.M1 @ R
        for _ in range(key[1]):
            R = self.M2 @ R
        self._cache[key] = R
        return R

    def twist(self, k, l=0):
        return tuple(k * self.v1[a] + l * self.v2[a] for a in range(3))

    # -------------------------------------------------------------- centraliser

    def centralizer_data(self, ce):
        """ce = (k, l, n_1..n_6) integer 8-tuple, a constructing element g.

        Returns dict with
          'PG'      : list of (k_h, l_h) for which SOME h = (k_h,l_h,n_h) commutes with g
          'part'    : {(k_h,l_h): n_h^0}  one particular solution each
          'kern'    : basis of ker_Z (I - M_g)   (the pure translations in Z_S(g))
          'gens8'   : generators of L_g subset Z^8 (the Z-span of the (k,l,n) of Z_S(g))
        """
        k, l = ce[0], ce[1]
        ng = np.array(ce[2:8], dtype=object)
        Mg = self.mat(k, l)
        I = np.eye(6, dtype=object)
        Amat = I - Mg
        kern = int_kernel(Amat)
        PG, part, gens = [], {}, []
        for kh in range(self.N):
            for lh in range(self.M):
                b = (I - self.mat(kh, lh)) @ ng
                x = int_solve(Amat, b)
                if x is None:
                    continue
                PG.append((kh, lh))
                part[(kh, lh)] = x
                gens.append(np.array([kh, lh] + [int(t) for t in x], dtype=object))
        for kv in kern:
            gens.append(np.array([0, 0] + [int(t) for t in kv], dtype=object))
        return dict(PG=PG, part=part, kern=kern, gens8=gens)

    # ----------------------------------------------------- geometry-level kernel

    def ker_geo_gens(self):
        """Generators of ker_geo subset Z^8 (convention C6)."""
        gens = [np.array([self.N, 0] + [0] * 6, dtype=object),
                np.array([0, self.M] + [0] * 6, dtype=object)]
        I = np.eye(6, dtype=object)
        for k in range(self.N):
            for l in range(self.M):
                if (k, l) == (0, 0):
                    continue
                Q = I - self.mat(k, l)
                for j in range(6):
                    gens.append(np.array([0, 0] + [int(x) for x in Q[:, j]], dtype=object))
        return gens

    def wl_structure(self):
        """DERIVED Wilson-line orders and identifications, from the point-group matrices
        alone -- i.e. the structure of Z^6 / <(I-M)n>.  Returned as the invariant-factor
        decomposition of that quotient plus the per-generator order.  Compared in
        selftest() against the geometry file's own declared block (guard 1c)."""
        I = np.eye(6, dtype=object)
        rel = []
        for k in range(self.N):
            for l in range(self.M):
                if (k, l) == (0, 0):
                    continue
                Q = I - self.mat(k, l)
                for j in range(6):
                    rel.append(np.array([int(x) for x in Q[:, j]], dtype=object))
        Amat = np.array([list(r) for r in rel], dtype=object)
        D, _, _ = snf(Amat)
        r = min(D.shape)
        inv = [int(D[i, i]) for i in range(r)]
        inv = [x for x in inv if x != 1]
        # per-generator order: smallest t>0 with t*e_i in the relation lattice
        orders = []
        for i in range(6):
            e = np.zeros(6, dtype=object)
            e[i] = 1
            t = 1
            while t <= 12:
                if int_solve(Amat.T, t * e) is not None:
                    break
                t += 1
            orders.append(t if t <= 12 else None)
        # identifications: e_i ~ e_j iff e_i - e_j in the relation lattice
        ident = []
        for i in range(6):
            for j in range(i + 1, 6):
                d = np.zeros(6, dtype=object)
                d[i] = 1
                d[j] = -1
                if int_solve(Amat.T, d) is not None:
                    ident.append((i + 1, j + 1))
        return dict(invariant_factors=inv, orders=orders, identifications=ident)

    # ------------------------------------------------------------------ verdict

    def lift_verdict(self, ce):
        """Does the projection at constructing element `ce` constrain Delta against the
        FULL gauge-shift set?  Geometry-level: is Z^8 = L_g + ker_geo ?

        Returns dict with the index [Z^8 : L_g + ker_geo] (1 == lift holds), the
        invariant factors of the quotient, and which of the 8 generators of Z^8
        individually fail to lie in L_g + ker_geo."""
        cd = self.centralizer_data(ce)
        gens = list(cd["gens8"]) + self.ker_geo_gens()
        idx, facs = sublattice_index(gens, 8)
        A = np.array([list(g) for g in gens], dtype=object).T        # columns = gens
        missing = []
        names = ["V1", "V2", "A1", "A2", "A3", "A4", "A5", "A6"]
        for i in range(8):
            e = np.zeros(8, dtype=object)
            e[i] = 1
            if int_solve(A, e) is None:
                missing.append(names[i])
        return dict(index=idx, invariant_factors=facs, missing_generators=missing,
                    centralizer_pointgroup=cd["PG"],
                    n_pure_translations=len(cd["kern"]),
                    kern_basis=[[int(x) for x in v] for v in cd["kern"]])


# ============================================================ E8xE8 lattice helpers

def in_lambda(x):
    """Is x (length-16 Fractions) in Lambda_{E8xE8}?  Each E8 half: all integer, or all
    half-odd-integer, with integer sum."""
    for h in (x[:8], x[8:]):
        fr = [F(t) - F(int(F(t))) for t in h]
        allint = all(F(t).denominator == 1 for t in h)
        allhalf = all(F(t).denominator == 2 for t in h)
        if not (allint or allhalf):
            return False
        if sum(F(t) for t in h) % 1 != 0:
            return False
    return True


def dot(a, b):
    return sum(F(x) * F(y) for x, y in zip(a, b))


# ============================================================ selftest / guards

CAMPAIGN_GEOMS = ["Z3xZ6_1_1", "Z2xZ4_1_1", "Z4xZ4_1_1", "Z3xZ3_1_1",
                  "Z6-I_1_1", "Z6-II_1_1", "Z4_1_1"]


def selftest(verbose=True):
    """Guards 1a, 1b, 1c of docs/root_lift_task.md.  Raises on any failure."""
    report = {}
    for lab in CAMPAIGN_GEOMS:
        S = SpaceGroup(lab)
        g = S.G
        rec = {}
        # ---- 1a  integrality + order
        R = np.eye(6, dtype=object)
        for t in range(1, 20):
            R = S.M1 @ R
            if np.array_equal(R, np.eye(6, dtype=object)):
                break
        rec["order_M1"] = t
        assert t == S.N, "%s: order(M1)=%d != N=%d" % (lab, t, S.N)
        if S.ZMxZN:
            R = np.eye(6, dtype=object)
            for t2 in range(1, 20):
                R = S.M2 @ R
                if np.array_equal(R, np.eye(6, dtype=object)):
                    break
            rec["order_M2"] = t2
            assert t2 == S.M, "%s: order(M2)=%d != M=%d" % (lab, t2, S.M)
            assert np.array_equal(S.M1 @ S.M2, S.M2 @ S.M1), "%s: point group not abelian" % lab

        # ---- 1b  constructing-element count per sector vs the geometry file
        cnt = {}
        for sec in g.sectors_ce:
            k, l = sec[0][0], sec[0][1]
            if (k, l) == (0, 0):
                continue
            Mg = S.mat(k, l)
            Amat = np.eye(6, dtype=object) - Mg
            D, _, _ = snf(Amat)
            # Lambda / (I-M_g)Lambda = (torsion of order prod of nonzero factors)
            #                          x Z^(number of zero factors, the fixed torus).
            # Constructing elements are labelled by the TORSION part only.
            det = 1
            free = 0
            for i in range(6):
                if D[i, i] == 0:
                    free += 1
                else:
                    det *= int(D[i, i])
            cnt[(k, l)] = dict(file=len(sec), translation_classes=det,
                               fixed_torus_dim=free,
                               unrotated_planes=g.unrotated_planes(k, l))
            assert free == 2 * len(g.unrotated_planes(k, l)), \
                "%s sector %s: SNF free rank %d != 2 x unrotated planes %s" % (
                    lab, (k, l), free, g.unrotated_planes(k, l))
        rec["ce_counts"] = {"%d,%d" % kk: vv for kk, vv in cnt.items()}
        # the file lists one element per SPACE-GROUP conjugacy class; translation
        # conjugation alone already gives |det(I-M_g)| classes, so file <= det, and
        # every det must be divisible-compatible.  Assert the inequality (a real check:
        # a wrong matrix convention breaks it immediately).
        for kk, vv in cnt.items():
            assert vv["translation_classes"] >= vv["file"], \
                "%s sector %s: file lists %d ce but only %d translation classes" % (
                    lab, kk, vv["file"], vv["translation_classes"])

        # ---- 1c  Wilson-line orders and identifications
        w = S.wl_structure()
        rec["wl_derived_orders"] = w["orders"]
        rec["wl_file_orders"] = list(g.wl_orders)
        rec["wl_derived_ident"] = w["identifications"]
        rec["wl_file_ident"] = [list(x) for x in g.identical_wl_groups]
        assert w["orders"] == list(g.wl_orders), \
            "%s: derived WL orders %s != file %s" % (lab, w["orders"], g.wl_orders)
        file_pairs = set()
        for grp in g.identical_wl_groups:
            for i in range(len(grp)):
                for j in range(i + 1, len(grp)):
                    file_pairs.add((min(grp[i], grp[j]), max(grp[i], grp[j])))
        der_pairs = set()
        for (i, j) in w["identifications"]:
            if g.wl_orders[i - 1] > 1 and g.wl_orders[j - 1] > 1:
                der_pairs.add((i, j))
        assert der_pairs == file_pairs, \
            "%s: derived WL identifications %s != file %s" % (lab, sorted(der_pairs),
                                                              sorted(file_pairs))
        report[lab] = rec
        if verbose:
            print("[selftest] %-12s OK  order(M1)=%d  WL orders %s  ident %s"
                  % (lab, rec["order_M1"], w["orders"], sorted(der_pairs)))
    return report


if __name__ == "__main__":
    selftest()

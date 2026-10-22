"""ALPHABET SURVEY PHASE 2 -- classical instanton NORM and channel alphabets on
Z_M x Z_N geometries.

Generalization of the phase-1 `alphabet_survey/code/survey_norm.py` (itself the
gate-certified generalization of `hierarchy_analysis/code/gate/gate_norm.py`).
Equation-level extensions (stated in full in `results/phase2_findings.md` E1-E8):

  E2  channel (point-group) rule: a trilinear sector multiset {g1, g2, g3},
      g_alpha = (k_alpha, l_alpha), is a channel iff
          sum k_alpha = 0 mod M   AND   sum l_alpha = 0 mod N
      (replaces sum k = 0 mod N; this is the rotation part of the space-group
      closure g1 g2 g3 = (1, 0), from the group law and nothing else).
  E3  capability (1301.2322 (2.22)/(2.23), L = 3): the coefficient counting is
      LOCAL -- it depends on the participating twist fields only through their
      plane-i monodromy exponents.  For a Z_M x Z_N sector these are
          k^i_alpha = phi_i(g_alpha) = (k_alpha v^1_i + l_alpha v^2_i) mod 1,
      substituted into the SAME formulas on the SAME domains
      (holomorphic 0 < k <= 1, antiholomorphic 0 <= k < 1):
          M_i^holo = sum_alpha [k^i_alpha or 1],   holo capable iff 2 - M^holo > 0
          M_i^anti = sum_alpha  k^i_alpha,         anti capable iff M^anti - 1 > 0.
      This substitution is exactly the "analogous discussion" 1301.2322 footnote 1
      declines to make; it is the PROGRAM'S OWN extension, owned as such.
  E4  loop integers (the sentence after (2.29)): smallest positive (n_p, m_p) with
          n_p g_p = m_p g_{p+1}   in Z_M x Z_N
      (componentwise mod (M, N); replaces n_p k_p = m_p k_{p+1} mod N; same
      lexicographic smallest-n-then-smallest-m convention as the certified code).
      Prefactor |1 - R_a(n_p g_p)|^2 with R_a the plane rotation of E1.
  E8  prefactor ISOTROPY: the FIXED ratio form (gate-certified,
      `multi_geometry_campaign/results/gate_z4_foundations.md` Z1):
          pf_a = rho(A u) / rho(u),  A = 1 - R_a(n_1 g_1),
      which equals |1 - e^{i phi}|^2 = 2 - 2 cos(2 pi n_1 phi_a(g_1)) exactly and is
      u-independent for any rotation.  Evaluated on u in {e1, e2, e1+e2} and required
      EXACTLY equal (rational arithmetic).  The phase-1 equality-of-values check is
      valid only when the test vectors have equal Gram norm (hexagonal/square with
      equal basis norms); the ratio form is valid on every Gram.  This is the fix
      the Z4-foundations gate certified with zero contamination.

Space-group rule, lambda part (identical structure to the certified engine, with
sector rotations R_a(g) and roto-translation cosets carried by the fixed points
themselves, p2_geom E5/E7):  lambda_1 + R_1 lambda_2 + R_1 R_2 lambda_3 = 0 with
lambda_a = (1 - R_a) P_a; leg 1 gauge-fixed to orbit representatives, leg 2 scanned
over orbit x lattice box, leg 3 SOLVED and accepted iff its fixed point lands in
leg 3's orbit.  Box adaptivity as phase 1: every minimum recomputed at box+2 and
rejected unless stable.
"""

import itertools
import os
import sys
from fractions import Fraction as F

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")))
from survey_geom import ap, minv, mm, one_minus, modl               # noqa: E402

ORDERS = [(0, 1, 2), (1, 2, 0), (2, 0, 1), (1, 0, 2), (0, 2, 1), (2, 1, 0)]


def sector_name(g):
    return "U" if g == (0, 0) else "T(%d,%d)" % g


# ---------------------------------------------------------------- Rule 5 (capability, E3)
def capable(g, gs, a):
    """(holomorphic?, antiholomorphic?) instanton coefficients in plane a, L = 3."""
    kk = [g.phi(a, s) for s in gs]
    M_holo = sum(x if x != 0 else F(1) for x in kk)
    M_anti = sum(kk)
    return (2 - M_holo > 0), (M_anti - 1 > 0)


def all_legs_localized(g, gs, a):
    """Every leg needs an isolated fixed point in plane a (untwisted-in-plane legs,
    fixed tori and freely-acting shifted tori are all delocalized -- E7)."""
    return all(s != (0, 0) and g.phi(a, s) != 0 for s in gs)


# ---------------------------------------------------------------- the NORM
def n_p(g1, g2, M, N):
    """Smallest (n, m) with n g1 = m g2 in Z_M x Z_N (E4)."""
    O = 1
    for x in (M, N):
        O = O * x // __import__("math").gcd(O, x)
    for n in range(1, O + 1):
        for m in range(1, O + 1):
            if ((n * g1[0] - m * g2[0]) % M == 0 and
                    (n * g1[1] - m * g2[1]) % N == 0):
                return n, m
    raise ValueError((g1, g2, M, N))


def prefactor(g, a, g1, g2):
    """|1 - R_a(n_1 g_1)|^2, FIXED ratio form (E8).

    Returns (value, 'ok'), (0, 'degenerate'), or (None, 'ANISOTROPIC') if the three
    exact ratios disagree (impossible for a true rotation; kept as a hard check)."""
    n, _ = n_p(g1, g2, g.M, g.N)
    ng = ((n * g1[0]) % g.M, (n * g1[1]) % g.N)
    A = one_minus(g.rot(a, ng))
    if A == ((F(0), F(0)), (F(0), F(0))):
        return F(0), "degenerate"
    ratios = set()
    for u in ((F(1), F(0)), (F(0), F(1)), (F(1), F(1))):
        ratios.add(g.rho(a, ap(A, u)) / g.rho(a, u))
    if len(ratios) != 1:
        return None, "ANISOTROPIC"
    return ratios.pop(), "ok"


_seps_cache = {}


def allowed_seps(g, a, gs, orbs, order, box):
    """Multiset of rho(P_2 - P_1) over space-group-allowed configurations.

    Identical to the certified engine with R_a(g) in place of theta^k; cached on the
    ordered (sector, orbit) configuration so repeated-leg orderings are not
    recomputed."""
    g1, g2, g3 = [gs[i] for i in order]
    o1, o2, o3 = [orbs[i] for i in order]
    key = (id(g), a, g1, g2, g3, tuple(o1), tuple(o2), tuple(o3), box)
    if key in _seps_cache:
        return _seps_cache[key]
    R1, R2, R3 = [g.rot(a, s) for s in (g1, g2, g3)]
    A1, A2, A3 = [one_minus(R) for R in (R1, R2, R3)]
    A3i = minv(A3)
    if A3i is None:
        _seps_cache[key] = None
        return None
    R12i = minv(mm(R1, R2))
    if R12i is None:
        _seps_cache[key] = None
        return None
    B = mm(A3i, R12i)                       # P3 = -B (A1 P1 + R1 A2 P2)
    R1A2 = mm(R1, A2)
    o3m = {modl(f) for f in o3}
    G = g.gram[a]
    du_i = ap(R1A2, (F(1), F(0)))
    du_j = ap(R1A2, (F(0), F(1)))
    seps = {}
    for f1 in o1:
        t = ap(A1, f1)
        for f2 in o2:
            u0 = ap(R1A2, f2)
            base = (t[0] + u0[0], t[1] + u0[1])
            v0 = (f2[0] - f1[0], f2[1] - f1[1])
            # rho(v0 + (i,j)) expanded once
            c0 = g.rho(a, v0)
            cx = 2 * (G[0][0] * v0[0] + G[0][1] * v0[1])
            cy = 2 * (G[0][1] * v0[0] + G[1][1] * v0[1])
            for i in range(-box, box + 1):
                bi = (base[0] + i * du_i[0], base[1] + i * du_i[1])
                for j in range(-box, box + 1):
                    need = (bi[0] + j * du_j[0], bi[1] + j * du_j[1])
                    P3 = ap(B, (-need[0], -need[1]))
                    if modl(P3) not in o3m:
                        continue
                    r = (c0 + i * cx + j * cy
                         + G[0][0] * i * i + 2 * G[0][1] * i * j + G[1][1] * j * j)
                    seps[r] = seps.get(r, 0) + 1
    _seps_cache[key] = seps
    return seps


def entry_norm(g, a, gs, orbs, box=6):
    """Ordering-invariant NORM in plane a (structure identical to the certified
    engine: 6 leg orderings, exact invariance demanded, box stability demanded;
    'space-group forbidden' and 'rejected' are SEPARATE outcomes -- the phase-1
    D-6(a) None-conflation is not reproduced here)."""
    info = {"degenerate_orderings": [], "per_ordering": {}, "sep_spectrum": {},
            "box_stable": True, "anisotropic_orderings": []}
    vals = {}
    for od in ORDERS:
        tag = "".join(map(str, od))
        pf, st = prefactor(g, a, gs[od[0]], gs[od[1]])
        if st == "degenerate":
            info["degenerate_orderings"].append(tag)
            continue
        if st == "ANISOTROPIC":
            info["anisotropic_orderings"].append(tag)
            continue
        seps = allowed_seps(g, a, gs, orbs, od, box)
        if not seps:
            continue
        rmin = min(seps)
        seps2 = allowed_seps(g, a, gs, orbs, od, box + 2)
        if not seps2 or min(seps2) != rmin:
            info["box_stable"] = False
        info["per_ordering"][tag] = (str(pf), str(rmin))
        info["sep_spectrum"][tag] = sorted((str(x) for x in seps), key=lambda s: F(s))
        vals[tag] = pf * rmin
    info["n_admissible_orderings"] = len(vals)
    if not vals:
        info["outcome"] = "space-group forbidden"
        return None, info
    uniq = set(vals.values())
    if len(uniq) != 1:
        info["NON_INVARIANT"] = {k: str(v) for k, v in vals.items()}
        info["outcome"] = "REJECTED non-invariant"
        return None, info
    if not info["box_stable"]:
        info["REJECTED"] = "box-unstable minimum"
        info["outcome"] = "REJECTED box-unstable"
        return None, info
    info["outcome"] = "ok"
    return uniq.pop(), info


# ---------------------------------------------------------------- channel alphabet
def leg_class_tuples(g, s, planes):
    """{class tuple restricted to `planes`: #simultaneous orbits sharing it} --
    built from the affine simultaneous orbits (p2_geom E6), never a plane product."""
    out = {}
    for so in g.sim_orbits(s):
        key = tuple(so["class_tuple"][a - 1] for a in planes)
        out[key] = out.get(key, 0) + 1
    return out


def channel_alphabet(g, gs, box=6):
    """Full alphabet of a trilinear channel on a Z_M x Z_N geometry."""
    res = {"multiset": "".join(sector_name(s) for s in gs),
           "gs": [list(s) for s in gs], "planes": {}, "capable_planes": [],
           "computable_planes": [], "notes": []}

    for a in range(3):
        h, an = capable(g, list(gs), a)
        pl = {"rule5_holomorphic": h, "rule5_antiholomorphic": an,
              "capable": bool(h or an),
              "all_legs_localized": all_legs_localized(g, list(gs), a),
              "plane_shape_rigid": g.rigid_plane(a),
              "gram": [[str(y) for y in r] for r in g.gram[a]]}
        if pl["capable"]:
            res["capable_planes"].append(a + 1)
        if not pl["capable"] or not pl["all_legs_localized"]:
            pl["status"] = ("capable but a leg is DELOCALIZED" if pl["capable"]
                            else "not capable")
            res["planes"]["plane%d" % (a + 1)] = pl
            continue

        orbs = [g.orbits(a, s) for s in gs]
        pl["orbit_classes_per_leg"] = [len(o) for o in orbs]
        pl["orbit_sizes_per_leg"] = [[len(x) for x in o] for o in orbs]
        rows, norms = [], set()
        for combo in itertools.product(*[range(len(o)) for o in orbs]):
            choice = [orbs[i][combo[i]] for i in range(3)]
            Nv, info = entry_norm(g, a, list(gs), choice, box=box)
            row = {"orbit_class_index": list(combo),
                   "orbit_sizes": [len(c) for c in choice],
                   "NORM": (None if Nv is None else str(Nv)),
                   "outcome": info["outcome"],
                   "degenerate_orderings": info["degenerate_orderings"],
                   "n_admissible_orderings": info["n_admissible_orderings"],
                   "per_ordering_pf_rhomin": info["per_ordering"]}
            for key in ("NON_INVARIANT", "REJECTED", "anisotropic_orderings"):
                if info.get(key):
                    row[key] = info[key]
            if Nv is not None:
                norms.add(Nv)
            rows.append(row)
        pl["assignments"] = rows
        pl["NORMs"] = sorted((str(x) for x in norms), key=lambda s: F(s))
        pl["status"] = "computable"
        res["computable_planes"].append(a + 1)
        res["planes"]["plane%d" % (a + 1)] = pl

    planes = res["computable_planes"]
    if not planes:
        res["alphabet_tuples"] = [["0", "0", "0"]]
        res["alphabet_total"] = ["0"]
        res["leg_class_tuples"] = None
        res["entries"] = None
        res["note_all_unsuppressed"] = (
            "no computable plane: every classical factor is exactly 1")
        return res

    norm_tab = {}
    for a in planes:
        pl = res["planes"]["plane%d" % a]
        norm_tab[a] = {tuple(row["orbit_class_index"]): row["NORM"]
                       for row in pl["assignments"]}

    legs = [leg_class_tuples(g, gs[i], planes) for i in range(3)]
    res["leg_class_tuples"] = [
        {"sector": sector_name(gs[i]),
         "n_distinct_class_tuples": len(legs[i]),
         "n_simultaneous_orbits": sum(legs[i].values()),
         "class_tuples": {repr(k_): v for k_, v in sorted(legs[i].items(), key=repr)}}
        for i in range(3)]

    entries, tuples = [], set()
    for c1 in legs[0]:
        for c2 in legs[1]:
            for c3 in legs[2]:
                t = [F(0), F(0), F(0)]
                ok = True
                for pi, a in enumerate(planes):
                    key = (c1[pi], c2[pi], c3[pi])
                    v = norm_tab[a].get(key)
                    if v is None:
                        ok = False
                        break
                    t[a - 1] = F(v)
                if not ok:
                    continue
                entries.append({"legs": [list(c1), list(c2), list(c3)],
                                "exponent_tuple": [str(x) for x in t],
                                "total": str(sum(t))})
                tuples.add(tuple(str(x) for x in t))
    res["entries"] = entries
    res["alphabet_tuples"] = sorted((list(t) for t in tuples), key=repr)
    res["alphabet_total"] = sorted({e["total"] for e in entries}, key=F)
    return res


def all_multisets_mn(M, N):
    """Every trilinear sector multiset {g1, g2, g3} with sum = (0, 0) mod (M, N) (E2)."""
    G = [(k, l) for k in range(M) for l in range(N)]
    seen, out = set(), []
    for g1 in G:
        for g2 in G:
            g3 = ((-g1[0] - g2[0]) % M, (-g1[1] - g2[1]) % N)
            key = tuple(sorted((g1, g2, g3)))
            if key in seen:
                continue
            seen.add(key)
            out.append(key)
    return out

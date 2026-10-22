"""
z24sw_lib.py --- the Z2xZ4_1_1 sweep library (the square arc's last catalogue).

CONVENTIONS AND PROVENANCE, STATED ONCE, AT THE TOP (CLAUDE.md 5.19(i)).

  * GEOMETRY.  Everything is derived from the STOCK orbifolder geometry file
    `orbifolder-1.2.1/Geometry/Geometry_Z2xZ4_1_1.txt` through the CERTIFIED
    alphabet-survey phase-2 engine `alphabet_survey/code/phase2/p2_geom.py::GeomMN`
    (gate: `alphabet_survey/results/gate_phase_two.md`; its `invariant_check` reproduces
    the file's own printed `constructing elements` sector by sector).  The file's FIRST
    twist row is the Z_M generator theta (sector label k), the second the Z_N generator
    omega (label l) -- the engine asserts this against the file's own
    `begin twist space group generators` block.
  * NORMS / EXPONENTS.  The classical suppression exponents come from the CERTIFIED
    engine `alphabet_survey/code/phase2/p2_norm.py` (the same code that produced the
    certified payload `alphabet_survey/data/phase2/Z2xZ4_1_1.json`), never a
    reimplementation.  This module's own table is cross-checked cell-for-cell against
    that payload in the STEP-0 self-test.
  * THE CAP IS QUOTED, NEVER RE-TESTED (guard 2).  `gate_z4xz4_foundations.md` K2 certified
    the single-plane cap of 2 as SQUARE-GENERIC *by computation on this very geometry*:
    Z2xZ4_1_1, 12 capable channel-planes, [Lambda:L] = 2 on 12/12, classes-per-leg
    [2,2,3] on 12/12, allowed-classes-per-bucket histogram {1: 168, 2: 24}, CAP 2.
    `assert_cap_quoted_not_retested` RAISES on any observed excess: a single-plane rank-3
    column is a parser-bug signal first and an escalation second.
  * RULE NC IS IMPOSED FROM THE START (guard 1; CLAUDE.md 5.22; certified statement in
    `results/gate_rule_nc.md` N10).  On any plane p on which the lattice factorizes, if a
    coupling satisfies the space-group rule and at least one leg is DELOCALISED in p, then
    all legs LOCALISABLE in p carry the SAME point-group ORBIT of plane-p fixed points.
    Necessary, not sufficient; vacuous with fewer than two localisable legs.  The ORBIT
    form is the certified one -- strict fixed-point equality is REFUTED campaign-wide
    (860,992 failures).  Every index-level table this module builds carries NC.
  * THE {3,4,6} LETTER-CARRYING THEOREM, applied to the order-2 generator's plane.
    On Z2xZ4_1_1, phi_2(k,l) = k/2, so a leg is localised in plane 2 iff k is odd; the
    point-group rule forces sum k = 0 mod 2, so plane 2 carries 0 or exactly 2 localisable
    legs and is NEVER capable.  Plane 2 therefore carries NO LETTER -- and is exactly where
    RULE NC bites.  Derived in `plane2_theorem()`, asserted in the self-test.

BRANCH TAGS that travel with every statement this module produces:

  [D-3] degenerate-prefactor : a cell may be space-group ALLOWED yet have a degenerate
                               classical prefactor, so (2.30)/(2.31) gives it no norm.
                               Branch 'suppressed' reads such a cell as absent; branch
                               'live' reads it as PRESENT WITH UNDETERMINED EXPONENT.
                               NEVER as an absence.  This geometry's own degenerate-cell
                               census is DERIVED (`degenerate_census()`), not inherited.
  [D-5] K reading            : Rule 6's K_sector vs K_orbit lift, carried per channel.
  gamma branch               : gamma is a PAIR (gamma^theta, gamma^omega) on a Z_M x Z_N
                               geometry; the character value on h = (k,l) is
                               k*gamma^theta + l*gamma^omega.  Branches 'correlator'
                               (sgn +1), 'LEEFT' (sgn -1) and the 'blind' control.
  K6 R-charge scope          : plane orders here are (4, 2, 4).  Plane 2 has order 2, which
                               IS prime -- so unlike Z4xZ4_1_1 this geometry HAS a prime
                               plane and the R-charge route is available there in principle;
                               planes 1 and 3 are non-prime.  Stated at every point of use.
  sY / sC orientation        : hypercharge sign and colour orientation are INDEPENDENT
                               per-model conventions; all four branches computed (5.21).
  protection                 : PROTECTED / IMPOSED / UNDETERMINED via the CORRECTED
                               full-coset machinery (gate_z4xz4_sweep S5); UNDETERMINED is
                               a reported category, never silently IMPOSED.

SCOPE.  Every model-level number produced here is a CATALOGUE statement about the 11,540
PVW `Z2xZ4_1_1` models.  Geometry statements are marked as such and never conflated.
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

import gzip
import json
import os
import re
import sys
from collections import defaultdict, Counter
from fractions import Fraction as F
from itertools import product

REPO = _PKG_ROOT
CAMP = os.path.join(REPO, "multi_geometry_campaign")
DSW = os.path.join(CAMP, "data", "z2xz4_sweep")
GEOM_FILE = os.path.join(REPO, "orbifolder-1.2.1", "Geometry", "Geometry_Z2xZ4_1_1.txt")
CATALOGUE = os.path.join(CAMP, "data", "literature", "parr2020_full", "extracted",
                         "Z2xZ4_1_1.txt")
CATALOGUE_SHA256 = "6c73de64aeddb5688a50789d8afcf58ff4cb1abe598ad12d49a7c4a39a6d2b7d"
CATALOGUE_BYTES = 7560329
N_MODELS = 11540
CERTIFIED_ALPH = os.path.join(REPO, "alphabet_survey", "data", "phase2", "Z2xZ4_1_1.json")

# the certified engines (imported, never reimplemented)
sys.path.insert(0, os.path.join(REPO, "alphabet_survey", "code"))
sys.path.insert(0, os.path.join(REPO, "alphabet_survey", "code", "phase2"))
from p2_geom import GeomMN                                              # noqa: E402
import p2_norm as PN                                                    # noqa: E402
from survey_geom import ap, mm, minv, one_minus, modl, I2               # noqa: E402

CERTIFIED_CAP = 2            # gate_z4xz4_foundations.md K2, measured ON Z2xZ4_1_1
PLANE_ORDERS = (4, 2, 4)     # asserted against the geometry in the self-test
GAMMA_BRANCHES = ("blind", "correlator", "LEEFT")

_G = None


def geom():
    global _G
    if _G is None:
        _G = GeomMN(GEOM_FILE)
    return _G


def certified_alphabet():
    with open(CERTIFIED_ALPH, "r", encoding="utf-8") as fh:
        return json.load(fh)


def sector_name(g):
    return "U" if tuple(g) == (0, 0) else "T(%d,%d)" % tuple(g)


def channel_name(gs):
    return "".join(sector_name(g) for g in gs)


def assert_cap_quoted_not_retested(observed_max, where):
    """Guard 2.  The cap is CERTIFIED (K2, on this geometry) and QUOTED.  A parser that
    appears to exceed it is a parser bug until proven otherwise -- raise loudly."""
    if observed_max > CERTIFIED_CAP:
        raise AssertionError(
            "CAP VIOLATION at %s: %d single-plane classes allowed in one bucket; the "
            "certified square-generic cap is %d (gate_z4xz4_foundations.md K2, measured on "
            "Z2xZ4_1_1 itself: 12/12 capable channel-planes, [Lambda:L]=2, bucket histogram "
            "{1:168, 2:24}).  PARSER BUG UNTIL PROVEN OTHERWISE -- ESCALATE, do not publish."
            % (where, observed_max, CERTIFIED_CAP))


# =======================================================================================
#  channels, capability, computability
# =======================================================================================

def all_channels():
    """Every trilinear sector multiset with sum k = 0 mod M and sum l = 0 mod N (E2).

    INCLUDES multisets with untwisted legs: a U leg is delocalised in every plane, so such
    a channel carries no letter, but it is a legitimate superpotential channel and the
    catalogue-wide block scan must not scope it out (guard 3, the S1-F completeness lesson).
    """
    g = geom()
    return [tuple(t) for t in PN.all_multisets_mn(g.M, g.N)]


def capable_planes(gs):
    """Planes on which a classical exponent is defined: Rule 5 coefficient condition AND
    every leg localised (the certified `p2_norm.channel_alphabet` criterion)."""
    g = geom()
    out = []
    for a in range(3):
        h, an = PN.capable(g, list(gs), a)
        if (h or an) and PN.all_legs_localized(g, list(gs), a):
            out.append(a + 1)
    return out


def localisable_legs(gs, a):
    """Indices of legs with an isolated fixed point in plane a (1-based plane)."""
    g = geom()
    return [i for i, s in enumerate(gs)
            if tuple(s) != (0, 0) and g.phi(a - 1, tuple(s)) != 0]


def nc_planes(gs):
    """Planes on which RULE NC has a premise for this channel: at least one leg
    DELOCALISED and at least two legs LOCALISABLE.  Returns {plane: [leg indices]}."""
    out = {}
    for a in (1, 2, 3):
        loc = localisable_legs(gs, a)
        if len(loc) >= 2 and len(loc) < 3:
            out[a] = loc
    return out


def plane2_theorem():
    """THE {3,4,6} LETTER-CARRYING THEOREM, applied here to the ORDER-2 generator's plane.

    Derived, not asserted:
      phi_a(k,l) = k v1_a + l v2_a.  On Z2xZ4_1_1 the file's twists are
      v1 = (0, 1/2, -1/2) and v2 = (1/4, 0, -1/4), so phi_2(k,l) = k/2 exactly.
      A leg is localised in plane 2 iff k is odd.  The point-group (channel) rule forces
      sum_alpha k_alpha = 0 mod 2, so the number of odd-k legs in any channel is EVEN:
      0 or 2.  A plane is capable only when ALL THREE legs are localised, so plane 2 is
      NEVER capable -- it carries NO LETTER -- and whenever it carries any localised leg
      it carries EXACTLY TWO, which is precisely RULE NC's premise.

    Returns the measured evidence for the self-test."""
    g = geom()
    rows = []
    for gs in all_channels():
        n_odd = sum(1 for s in gs if s[0] % 2 == 1)
        loc = localisable_legs(gs, 2)
        rows.append({"channel": channel_name(gs), "n_odd_k": n_odd,
                     "n_localisable_plane2": len(loc),
                     "plane2_capable": 2 in capable_planes(gs)})
    return {"phi2_is_k_over_2": all(g.phi(1, (k, l)) == F(k % 2, 2)
                                    for k in range(g.M) for l in range(g.N)),
            "n_localisable_plane2_histogram":
                dict(Counter(r["n_localisable_plane2"] for r in rows)),
            "plane2_capable_anywhere": any(r["plane2_capable"] for r in rows),
            "plane_orders": [g.plane_order(a) for a in range(3)],
            "rows": rows}


# =======================================================================================
#  RULE NC  (guard 1 -- imposed in EVERY table this module builds)
# =======================================================================================

def nc_orbit_compatibility(gs, a):
    """The orbit-vs-equality check RULE NC's index comparison requires (5.22 / gate N1(b)).

    NC says the localisable legs share a point-group ORBIT.  Comparing orbit INDICES across
    two different sectors is only meaningful if the two sectors' plane-a fixed-point SETS
    and ORBIT PARTITIONS coincide.  They do, because the two legs' plane-a rotations are
    mutually inverse (SO(2) is abelian and the third leg is trivial there) and
    fix(g^-1) = fix(g).  VERIFIED here, per channel per plane, never assumed."""
    g = geom()
    loc = localisable_legs(gs, a)
    # NC's premise: at least one leg DELOCALISED and at least two LOCALISABLE.  With three
    # legs that is exactly |loc| == 2.  |loc| == 3 is a capable plane (no NC premise);
    # |loc| <= 1 is vacuous.
    if len(loc) != 2:
        return {"vacuous": True, "n_localisable": len(loc)}
    fps = [g.fixed_points(a - 1, tuple(gs[i])) for i in loc]
    orbs = [g.orbits(a - 1, tuple(gs[i])) for i in loc]
    same_fps = all(set(map(tuple, fps[0])) == set(map(tuple, f)) for f in fps[1:])
    canon = [sorted(tuple(sorted(map(tuple, o))) for o in ob) for ob in orbs]
    same_part = all(canon[0] == c for c in canon[1:])
    # the derivation's own premise: the localisable legs' rotations multiply to 1
    prod = I2
    for i in loc:
        prod = mm(prod, g.rot(a - 1, tuple(gs[i])))
    return {"vacuous": False, "legs": loc, "same_fixed_point_set": same_fps,
            "same_orbit_partition": same_part, "rotations_multiply_to_identity":
                prod == I2, "n_classes": len(orbs[0])}


def nc_ok(gs, cells):
    """RULE NC on a full cell.  `cells[i]` is leg i's class tuple over the THREE planes
    (None where the leg is delocalised).  Returns True iff, in every plane with a
    delocalised leg and >= 2 localisable legs, all localisable legs share the class index."""
    for a in (1, 2, 3):
        loc = localisable_legs(gs, a)
        if len(loc) < 2 or len(loc) == 3:
            continue
        vals = {cells[i][a - 1] for i in loc}
        if len(vals) != 1:
            return False
    return True


# =======================================================================================
#  the index-level exponent table  (NC-corrected, [D-3]-tagged)
# =======================================================================================

_TABLE = None


def _norm_table(gs, a, box=6):
    """{(c0,c1,c2) in plane a: (NORM or None, outcome)} using the CERTIFIED entry_norm."""
    g = geom()
    orbs = [g.orbits(a - 1, tuple(s)) for s in gs]
    out = {}
    for combo in product(*[range(len(o)) for o in orbs]):
        choice = [orbs[i][combo[i]] for i in range(3)]
        Nv, info = PN.entry_norm(g, a - 1, list(gs), choice, box=box)
        out[combo] = (None if Nv is None else F(Nv), info["outcome"])
    return out


def allowed_exact(gs, a, classes):
    """EXACT space-group allowedness of one plane-a class triple; no box, no gauge fix.

    L := A1.Lambda + R1 A2.Lambda + R1 R2 A3.Lambda   (A_i = 1 - R_i)
    allowed  <=>  A1 P1 + R1 A2 P2 + R1 R2 A3 P3  in  L  for SOME points in the orbits.

    This is the same criterion the foundations gate certified (k1_cap_independent /
    k2_genericity), transcribed here for a general Z_M x Z_N plane.  It is INDEPENDENT of
    the prefactor, so the difference between it and NORM-presence is exactly the [D-3]
    degenerate class."""
    g = geom()
    Rs = [g.rot(a - 1, tuple(s)) for s in gs]
    if any(g.phi(a - 1, tuple(s)) == 0 for s in gs):
        return None
    A = [one_minus(R) for R in Rs]
    M1 = A[0]
    M2 = mm(Rs[0], A[1])
    M3 = mm(mm(Rs[0], Rs[1]), A[2])
    L = _lattice_hnf([M1, M2, M3])
    orbs = [g.orbits(a - 1, tuple(gs[i]))[classes[i]] for i in range(3)]
    for P1 in orbs[0]:
        for P2 in orbs[1]:
            for P3 in orbs[2]:
                v = (ap(M1, P1)[0] + ap(M2, P2)[0] + ap(M3, P3)[0],
                     ap(M1, P1)[1] + ap(M2, P2)[1] + ap(M3, P3)[1])
                if _in_lat(v, L):
                    return True
    return False


def _lattice_hnf(mats):
    from math import gcd
    cols = []
    for M in mats:
        cols.append([int(M[0][0]), int(M[1][0])])
        cols.append([int(M[0][1]), int(M[1][1])])
    cols = [c for c in cols if c != [0, 0]]
    if not cols:
        return ((0, 0), (0, 0))
    while True:
        nz = [c for c in cols if c[0] != 0]
        if len(nz) <= 1:
            break
        nz.sort(key=lambda c: abs(c[0]))
        p = nz[0]
        moved = False
        for c in cols:
            if c is p or c[0] == 0:
                continue
            q = c[0] // p[0]
            if q:
                c[0] -= q * p[0]
                c[1] -= q * p[1]
                moved = True
            elif abs(c[0]) < abs(p[0]):
                moved = True
        if not moved:
            break
    piv = [c for c in cols if c[0] != 0]
    rest = [c for c in cols if c[0] == 0]
    if not piv:
        a11, a21 = 0, 0
    else:
        piv.sort(key=lambda c: abs(c[0]))
        a11, a21 = piv[0]
        for c in piv[1:]:
            q = c[0] // a11
            rest.append([0, c[1] - q * a21])
    g2 = 0
    for c in rest:
        g2 = gcd(g2, abs(c[1]))
    return ((a11, 0), (a21, g2))


def _in_lat(v, L):
    a, b, gg = L[0][0], L[1][0], L[1][1]
    if a == 0:
        return v[0] == 0 and (v[1] == 0 if gg == 0 else (F(v[1]) / gg).denominator == 1)
    x = F(v[0]) / a
    if x.denominator != 1:
        return False
    r = v[1] - x * b
    if gg == 0:
        return r == 0
    return (F(r) / gg).denominator == 1


def leg_class_tuples(s):
    """{full 3-plane class tuple: n_simultaneous_orbits} for a sector, from the certified
    affine simultaneous orbits (p2_geom E6).  None on a delocalised plane."""
    g = geom()
    if tuple(s) == (0, 0):
        return {(None, None, None): 1}
    out = {}
    for so in g.sim_orbits(tuple(s)):
        out[tuple(so["class_tuple"])] = out.get(tuple(so["class_tuple"]), 0) + 1
    return out


def build_table(box=6):
    """The phase's index-level exponent table for EVERY channel, NC IMPOSED FROM THE START.

    A CELL is a triple of FULL 3-plane class tuples, one per leg.  For each cell:
      * RULE NC is checked in every plane with a delocalised leg and >= 2 localisable legs;
        a cell violating it is `nc_forbidden` and is NEVER counted as present.
      * in every COMPUTABLE plane the certified NORM must exist; if the plane is space-group
        allowed but prefactor-degenerate the cell is tagged [D-3] and carried in the `live`
        branch with an UNDETERMINED exponent (never read as an absence).
      * `total` = sum of per-plane norms over the computable planes.

    Returns {channel: {...}} and is cached on disk."""
    global _TABLE
    if _TABLE is not None:
        return _TABLE
    g = geom()
    out = {}
    for gs in all_channels():
        nm = channel_name(gs)
        comp = capable_planes(gs)
        ncp = nc_planes(gs)
        legs = [leg_class_tuples(s) for s in gs]
        ntab = {a: _norm_table(gs, a, box=box) for a in comp}
        atab = {a: {k: allowed_exact(gs, a, list(k)) for k in ntab[a]} for a in comp}
        entries, d3_cells, nc_killed = [], [], 0
        for c0 in legs[0]:
            for c1 in legs[1]:
                for c2 in legs[2]:
                    cells = (c0, c1, c2)
                    if not nc_ok(gs, cells):
                        nc_killed += 1
                        continue
                    t = [F(0), F(0), F(0)]
                    ok, degen = True, []
                    for a in comp:
                        key = (c0[a - 1], c1[a - 1], c2[a - 1])
                        nv, outc = ntab[a][key]
                        if nv is None:
                            if atab[a][key]:
                                degen.append(a)      # [D-3]: allowed, no norm
                            else:
                                ok = False
                                break
                        else:
                            t[a - 1] = nv
                    if not ok:
                        continue
                    e = {"legs": [list(c0), list(c1), list(c2)],
                         "exponent_tuple": [str(x) for x in t],
                         "total": str(sum(t)),
                         "D3_degenerate_planes": degen}
                    if degen:
                        d3_cells.append(e)
                    else:
                        entries.append(e)
        out[nm] = {"channel": nm, "gs": [list(x) for x in gs],
                   "capable_planes": comp, "nc_planes": {str(k): v for k, v in ncp.items()},
                   "n_cells_before_NC": int(_prod(len(x) for x in legs)),
                   "n_cells_killed_by_NC": nc_killed,
                   "entries": entries,
                   "D3_entries": d3_cells,
                   "alphabet_total": sorted({e["total"] for e in entries}, key=F),
                   "alphabet_tuples": sorted({tuple(e["exponent_tuple"]) for e in entries})}
    _TABLE = out
    return out


def _prod(it):
    r = 1
    for x in it:
        r *= x
    return r


def degenerate_census():
    """This geometry's OWN [D-3] census (guard 5), derived not inherited: per capable
    channel-plane, cells space-group ALLOWED by the exact criterion vs cells given a NORM."""
    rows = []
    for gs in all_channels():
        nm = channel_name(gs)
        for a in capable_planes(gs):
            nt = _norm_table(gs, a)
            n_allowed = sum(1 for k in nt if allowed_exact(gs, a, list(k)))
            n_norm = sum(1 for k, v in nt.items() if v[0] is not None)
            rows.append({"channel": nm, "plane": a, "n_cells": len(nt),
                         "n_space_group_allowed_exact": n_allowed,
                         "n_with_norm": n_norm, "n_degenerate": n_allowed - n_norm})
    return {"rows": rows,
            "n_capable_channel_planes": len(rows),
            "n_cells": sum(r["n_cells"] for r in rows),
            "n_allowed": sum(r["n_space_group_allowed_exact"] for r in rows),
            "n_with_norm": sum(r["n_with_norm"] for r in rows),
            "n_degenerate": sum(r["n_degenerate"] for r in rows),
            "degenerate_channel_planes": [(r["channel"], r["plane"]) for r in rows
                                          if r["n_degenerate"]]}


def bucket_cap_scan():
    """Guard 2 incidental check: allowed classes per (channel, plane, leg, other-two)
    bucket, by the exact SPACE-GROUP criterion alone.  Quoted against CERTIFIED_CAP; an
    excess RAISES."""
    g = geom()
    hist, mx, worst = Counter(), 0, []
    for gs in all_channels():
        for a in capable_planes(gs):
            orbs = [g.orbits(a - 1, tuple(s)) for s in gs]
            ncl = [len(o) for o in orbs]
            for leg in range(3):
                oth = [i for i in range(3) if i != leg]
                for o1 in range(ncl[oth[0]]):
                    for o2 in range(ncl[oth[1]]):
                        hits = []
                        for v in range(ncl[leg]):
                            cl = [0, 0, 0]
                            cl[leg] = v
                            cl[oth[0]] = o1
                            cl[oth[1]] = o2
                            if allowed_exact(gs, a, cl):
                                hits.append(v)
                        if not hits:
                            continue
                        hist[len(hits)] += 1
                        if len(hits) > mx:
                            mx = len(hits)
                        if len(hits) >= 3:
                            worst.append({"channel": channel_name(gs), "plane": a,
                                          "leg": leg, "classes": hits})
    assert_cap_quoted_not_retested(mx, "bucket_cap_scan")
    return {"histogram": dict(sorted(hist.items())), "max": mx,
            "buckets_with_3_or_more": worst,
            "CERTIFIED_CAP": CERTIFIED_CAP,
            "note": "QUOTED from gate_z4xz4_foundations.md K2 (measured on Z2xZ4_1_1); "
                    "this scan is an incidental consistency check, not a re-derivation."}


# =======================================================================================
#  the exchange automorphism sigma, DERIVED from the geometry  (guard 3 / protection)
# =======================================================================================

def _lattice_autos():
    """All signed-permutation-of-plane-basis isometries of the 6d lattice that map planes
    onto planes.  On this geometry every Gram is the identity and the three planes are
    spanned by (e1,e2), (e3,e4), (e5,e6), so the full isometry group of the (orthonormal)
    lattice is the signed permutation group of the 6 coordinates -- of which we need those
    permuting the plane BLOCKS.  Enumerated exactly; the count is checked against the
    certified tau-premise record (n_aut = 46080 for this lattice class)."""
    from itertools import permutations as perms
    g = geom()
    ident = ((F(1), F(0)), (F(0), F(1)))
    assert all(tuple(tuple(r) for r in g.gram[a]) == ident for a in range(3)), g.gram
    idx = [g.pairs[a] for a in range(3)]
    out = []
    for pperm in perms(range(3)):
        for inner in product((0, 1), repeat=3):
            for signs in product((1, -1), repeat=6):
                M = [[0] * 6 for _ in range(6)]
                ok = True
                for a in range(3):
                    src = idx[a]
                    dst = idx[pperm[a]]
                    d = (dst[1], dst[0]) if inner[a] else dst
                    M[src[0]][d[0]] = signs[2 * a]
                    M[src[1]][d[1]] = signs[2 * a + 1]
                if ok:
                    out.append((tuple(pperm), tuple(tuple(r) for r in M)))
    return out


def exchange_sigmas(pair=(1, 3)):
    """Every lattice isometry sigma that TRANSPOSES `pair`, fixes the third plane setwise,
    and NORMALISES the point group -- together with the induced point-group automorphism
    tau and the induced Wilson-line permutation-with-signs.

    Conditions imposed (the tau-premise memo's link L1-L3, re-derived here):
      * sigma is a lattice isometry (integer, orthogonal for the identity Gram);
      * sigma maps plane a onto plane sigma(a) with {a,b} transposed and c fixed;
      * sigma normalises the point group: sigma R_g sigma^{-1} = R_{tau(g)} in every plane
        simultaneously, for a single tau: Z_M x Z_N -> Z_M x Z_N;
      * tau is an automorphism of the point group and (checked, not assumed) an involution.

    Returns a list of distinct (tau, WL condition) pairs after exact deduplication.  The
    per-model lift condition is then EXISTENTIAL over this list (any admissible sigma that
    lifts protects the model) -- stated wherever protection is reported."""
    g = geom()
    a, b = pair
    c = ({1, 2, 3} - {a, b}).pop()
    want = {a - 1: b - 1, b - 1: a - 1, c - 1: c - 1}
    R6 = {}
    for kl in [(k, l) for k in range(g.M) for l in range(g.N)]:
        M = [[0] * 6 for _ in range(6)]
        for p in range(3):
            R = g.rot(p, kl)
            i, j = g.pairs[p]
            M[i][i], M[i][j] = int(R[0][0]), int(R[1][0])
            M[j][i], M[j][j] = int(R[0][1]), int(R[1][1])
        R6[kl] = tuple(tuple(r) for r in M)
    sectors = sorted(R6)

    def m6(A, B):
        return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(6)) for j in range(6))
                     for i in range(6))

    def inv6(A):                      # orthogonal integer matrix: inverse = transpose
        return tuple(tuple(A[j][i] for j in range(6)) for i in range(6))

    found = {}
    for pperm, S in _lattice_autos():
        if {p: pperm[p] for p in range(3)} != want:
            continue
        Si = inv6(S)
        tau = {}
        good = True
        for kl in sectors:
            conj = m6(m6(Si, R6[kl]), S)
            hit = [x for x in sectors if R6[x] == conj]
            if len(hit) != 1:
                good = False
                break
            tau[kl] = hit[0]
        if not good:
            continue
        # tau must be a group automorphism and an involution (checked, not assumed)
        if any(tau[tau[kl]] != kl for kl in sectors):
            continue
        if tau[(0, 0)] != (0, 0):
            continue
        addok = True
        for x in sectors:
            for y in sectors:
                z = ((x[0] + y[0]) % g.M, (x[1] + y[1]) % g.N)
                tz = ((tau[x][0] + tau[y][0]) % g.M, (tau[x][1] + tau[y][1]) % g.N)
                if tau[z] != tz:
                    addok = False
                    break
            if not addok:
                break
        if not addok:
            continue
        # induced Wilson-line map: sigma(e_i) = sum_j S[i][j] e_j  =>  W_i -> sum_j S[i][j] W_j.
        # EVERY Wilson line here has order dividing 2 (the file declares 2,2,2,2,2,2 and the
        # catalogue confirms it on 11,540/11,540), so -W = W mod Lambda and the condition
        # depends only on S mod 2.  Deduplicate on (tau, S mod 2) -- exactly the Z4xZ4 phase's
        # discipline, transcribed.
        wl = tuple(tuple(abs(S[i][j]) % 2 for j in range(6)) for i in range(6))
        key = (tuple(sorted(tau.items())), wl)
        if key not in found:
            found[key] = {"tau": {str(k): list(v) for k, v in sorted(tau.items())},
                          "generator_images": [list(tau[(1, 0)]), list(tau[(0, 1)])],
                          "wl_matrix": [list(r) for r in wl],
                          "sigma_matrix": [list(r) for r in S]}
    return list(found.values())


# =======================================================================================
#  catalogue parsing -- row identities DERIVED, never assumed
# =======================================================================================

_ROWPAT = re.compile(r"^\s*-?\d")


def parse_catalogue(path=CATALOGUE, limit=None):
    """{label: {'rows': [8 x 16 Fractions]}}.

    ROW IDENTITY IS DERIVED, NOT ASSUMED (the S2-inversion lesson).  `derive_row_identity`
    below decides, from the geometry file's own declarations, which rows are the two shift
    vectors and which are the six Wilson lines, by checking the ORDER conditions the file
    declares (`allowed order of Wilson lines`: 2,2,2,2,2,2) and the identical-WL block
    (`1 2`, `5 6`) against the catalogue itself on all 11,540 models."""
    out = {}
    lab, rows, inb = None, [], False
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith("begin model"):
                lab, rows, inb = None, [], True
                continue
            if line.startswith("end model"):
                if lab and len(rows) == 8:
                    out[lab] = {"rows": rows}
                inb = False
                if limit and len(out) >= limit:
                    break
                continue
            if not inb:
                continue
            if line.startswith("Label:"):
                lab = line.split(":", 1)[1].strip()
            elif "," in line and _ROWPAT.match(line):
                v = [F(x.strip()) for x in line.split(",")]
                if len(v) == 16:
                    rows.append(v)
    return out


def geometry_declarations():
    """The file's own `allowed order of Wilson lines` and `identical WLs` blocks."""
    st = [l.strip() for l in open(GEOM_FILE, encoding="utf-8")]
    i = st.index("allowed order of Wilson lines")
    orders = [int(st[i + 1 + k]) for k in range(6)]
    ident = []
    if "begin identical WLs" in st:
        i = st.index("begin identical WLs")
        j = st.index("end identical WLs")
        for l in st[i + 1:j]:
            if l:
                ident.append([int(x) for x in l.split()])
    return {"wl_orders": orders, "identical_WLs": ident}


def derive_row_identity(cat):
    """DERIVE which catalogue rows are (V1, V2) and which are (W1..W6).

    Evidence used, all from the geometry file itself and the catalogue, never from a
    convention:
      (a) the file declares `allowed order of Wilson lines` = 2 on all six, so every Wilson
          -line row must have order dividing 2 in (Lambda x Q)/Lambda;
      (b) the file declares `identical WLs` {1,2} and {5,6}, so the corresponding rows must
          be EQUAL on every model;
      (c) the twist orders are M = 2 and N = 4, so V1 must have order dividing 2 and V2
          order dividing 4 -- and V2 must genuinely reach order 4 somewhere in the
          catalogue, otherwise the assignment is degenerate.
    Returns the measured evidence plus the decided assignment."""
    ev = {"n_models": len(cat)}
    ordhist = {r: Counter() for r in range(8)}
    eq = Counter()
    for m in cat.values():
        rows = m["rows"]
        for r in range(8):
            ordhist[r][order_mod_lattice(rows[r])] += 1
        for (i, j) in ((2, 3), (4, 5), (6, 7), (2, 4), (4, 6)):
            if rows[i] == rows[j]:
                eq["%d==%d" % (i, j)] += 1
    ev["order_histogram_per_row"] = {str(r): dict(ordhist[r]) for r in range(8)}
    ev["row_equality_counts"] = dict(eq)
    ev["geometry_declarations"] = geometry_declarations()
    ok_wl = all(max(ordhist[r]) <= 2 for r in range(2, 8))
    ev["rows_2_to_7_all_order_le_2"] = ok_wl
    ev["rows_2_3_equal_on_all"] = eq.get("2==3", 0) == len(cat)
    ev["rows_6_7_equal_on_all"] = eq.get("6==7", 0) == len(cat)
    ev["rows_4_5_equal_on_all"] = eq.get("4==5", 0) == len(cat)
    ev["row0_order_le_2"] = max(ordhist[0]) <= 2
    ev["row1_order_divides_4"] = all(o in (1, 2, 4) for o in ordhist[1])
    ev["row1_reaches_order_4"] = ordhist[1].get(4, 0) > 0
    ev["assignment"] = {"V1": 0, "V2": 1, "W1": 2, "W2": 3, "W3": 4, "W4": 5,
                        "W5": 6, "W6": 7}
    ev["DERIVED"] = bool(ok_wl and ev["rows_2_3_equal_on_all"]
                         and ev["rows_6_7_equal_on_all"] and ev["row0_order_le_2"]
                         and ev["row1_order_divides_4"] and ev["row1_reaches_order_4"])
    return ev


# ---- E8xE8 lattice arithmetic (x32 integer units; from the CORRECTED gate machinery) ---

def to32(v):
    out = []
    for x in v:
        n = x * 32
        assert n.denominator == 1, "component %s has denominator not dividing 32" % x
        out.append(int(n))
    return out


def in_E8_32(v):
    r = v[0] % 32
    if r not in (0, 16):
        return False
    for x in v:
        if x % 32 != r:
            return False
    return sum(v) % 64 == 0


def in_E8xE8_32(v):
    return in_E8_32(v[:8]) and in_E8_32(v[8:])


def order32(v):
    for n in range(1, 33):
        if in_E8xE8_32([n * x for x in v]):
            return n
    return None


def order_mod_lattice(v):
    return order32(to32(v))


# =======================================================================================
#  orbifolder census parsing   (ported verbatim from the certified z44sw_lib / z36sw_lib,
#  with the geometry-dependent sector modulus taken from the geometry, not hard-coded)
# =======================================================================================

SM_ROLES = {(3, 2, F(1, 6)): "Q", (-3, 2, F(-1, 6)): "Qbar",
            (-3, 1, F(-2, 3)): "uc", (3, 1, F(2, 3)): "ucbar",
            (-3, 1, F(1, 3)): "dc", (3, 1, F(-1, 3)): "dcbar",
            (1, 2, F(-1, 2)): "L", (1, 2, F(1, 2)): "Hu",
            (1, 1, F(1)): "ec", (1, 1, F(-1)): "ecbar"}


def _rat(t):
    t = t.strip()
    return F(t) if t else F(0)


def _vec_pairs(s):
    out = []
    for grp in re.findall(r"\(([^)]*)\)", s):
        out += [_rat(x) for x in grp.split(",")]
    return out


def split_models(text):
    """[(label, block)] from a census segment.  orbifolder emits prompts INLINE, so the
    boundary regex must not be anchored at start-of-line (the Z4xZ4 lesson)."""
    idx = [(m.start(), m.group(1)) for m in re.finditer(r"/> /(MSSM\d+)> /", text)]
    out = []
    for i, (p, lab) in enumerate(idx):
        end = idx[i + 1][0] if i + 1 < len(idx) else len(text)
        out.append((lab, text[p:end]))
    return out


_OSC_TOK = re.compile(r"(-?)a_\(\s*(-?\d+(?:/\d+)?)\s*\)\^(?:\((-\d)\)|(\d))")


def parse_oscillators(s):
    """`no` or `a_(mode)^(+/-n)` tokens -> (N_L, Nbar_L) per plane.

    WHICH tower is the unbarred N_L is a LABELLING CHOICE (the 5.19 conjugation family) and
    is carried as a tagged branch, never assumed.  NOTE (gate_z4xz4_sweep S7.2): the MODE
    NUMBER is captured and retained here as well as the occurrence count, because counting
    occurrences where the law wants the mode-number sum is the 5.21(iv) 'right name,
    neighbouring quantity' failure that produced the Z4xZ4 R-charge anomalies."""
    NL, NbL = {1: 0, 2: 0, 3: 0}, {1: 0, 2: 0, 3: 0}
    ML, MbL = {1: F(0), 2: F(0), 3: F(0)}, {1: F(0), 2: F(0), 3: F(0)}
    s = (s or "no").strip()
    if s in ("no", ""):
        return NL, NbL, ML, MbL, s
    hits = _OSC_TOK.findall(s)
    if not hits:
        raise ValueError("unparsable oscillator string %r" % s)
    for _sgn, mode, neg, pos in hits:
        a = abs(int(neg if neg else pos))
        if a == 0 or a > 3:
            raise ValueError("oscillator plane index %d in %r" % (a, s))
        (NbL if neg else NL)[a] += 1
        (MbL if neg else ML)[a] += abs(F(mode))
    return NL, NbL, ML, MbL, s


def parse_spectrum(block, M=2, N=4):
    """{F_n: record} from `print(F) with internal information`."""
    heads = [(m.start(), m.group(1))
             for m in re.finditer(r"(?m)^\s{2,8}(F_\d+)\s*$", block)]
    out = {}
    for i, (pos, lab) in enumerate(heads):
        end = heads[i + 1][0] if i + 1 < len(heads) else min(len(block), pos + 3000)
        b = block[pos:end]

        def grab(pat):
            m = re.search(pat, b)
            return m.group(1) if m else None

        sec = grab(r"sector \(k,l\)\s*:\s*\(([^)]*)\)")
        if sec is None:
            continue
        kl = [int(x) for x in sec.split(",")]
        kl = (kl[0] % M, kl[1] % N)
        fps = grab(r"fixed point n_a\s*:\s*\(([^)]*)\)")
        fp = tuple(int(x) for x in fps.split(",")) if fps else None
        psh_s = grab(r"left-moving p_sh\s*:\s*(.*)")
        psh = tuple(_vec_pairs(psh_s)) if psh_s else None
        repm = re.search(
            r"representation\s*:\s*\(([^)]*)\)_([lr])\s+U\(1\):\s*\(([^)]*)\)", b)
        rep = chir = None
        if repm:
            rep = tuple(int(re.match(r"^(-?\d+)", x.strip()).group(1))
                        for x in repm.group(1).split(","))
            chir = repm.group(2)
        osc = grab(r"oscillators\s*:\s*(.*)")
        NL, NbL, ML, MbL, osc_raw = parse_oscillators(osc)
        qs = grab(r"right-moving q_sh\s*:\s*\(([^)]*)\)")
        gam = grab(r"gamma phases\s*:\s*\(([^)]*)\)")
        out[lab] = dict(label=lab, kl=kl, fp=fp, p_sh=psh, rep=rep, chir=chir,
                        osc=osc_raw,
                        qsh=tuple(_rat(x) for x in qs.split(",")) if qs else None,
                        gamma=tuple(_rat(x) for x in gam.split(",")) if gam else None,
                        NL=(NL[1], NL[2], NL[3]), NbL=(NbL[1], NbL[2], NbL[3]),
                        ML=(str(ML[1]), str(ML[2]), str(ML[3])),
                        MbL=(str(MbL[1]), str(MbL[2]), str(MbL[3])))
    return out


_ROW = re.compile(
    r"(?m)^\s*(\d+)\s+\(\s*(-?\d+),\s*(-?\d+)\)_([lr])\s+U\(1\)\s*:\s*\(\s*([^)]*)\)")


def parse_smconfig(block):
    """[(tY, sm_table)] -- ONE entry per printed SM vev-configuration; the caller SELECTS
    BY THE CLOSURE (never the first generator, never a merged table)."""
    cuts = [m.start() for m in re.finditer(r"Hypercharge generator:", block)]
    out = []
    for i, p in enumerate(cuts):
        end = cuts[i + 1] if i + 1 < len(cuts) else len(block)
        seg = block[p:end]
        m = re.search(r"Hypercharge generator:\s*\n(.*)", seg)
        tY = _vec_pairs(m.group(1)) if m else None
        tbl = [(int(mm_.group(2)), int(mm_.group(3)), _rat(mm_.group(5)),
                int(mm_.group(1)), mm_.group(4)) for mm_ in _ROW.finditer(seg)]
        out.append((tY, tbl))
    if not out:
        out = [(None, [(int(mm_.group(2)), int(mm_.group(3)), _rat(mm_.group(5)),
                        int(mm_.group(1)), mm_.group(4)) for mm_ in _ROW.finditer(block)])]
    return out


def parse_gauge_factors(block):
    m = re.search(r'Gauge group in vev-configuration "TestConfig1":\s*(.*)', block)
    if not m:
        return []
    s = m.group(1).split("and U(1)")[0]
    return [(t, int(n)) for t, n in re.findall(r"(SU|SO|E)\((\d+)\)", s)]


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def hidden_singlet(rep, a3, a2):
    return all(abs(d) == 1 for k, d in enumerate(rep) if k not in (a3, a2))


def observable_positions(spec, tY, factors, printed):
    """(a3, a2, sY) identified by REBUILDING orbifolder's own printed SMConfig table from
    the full left-chiral spectrum -- never by row position (the Z6-I closure)."""
    if not printed or tY is None:
        return None
    left = [r for r in spec.values() if r["rep"] and r["chir"] == "l" and r["p_sh"]]
    for a3 in range(len(factors)):
        if factors[a3] != ("SU", 3):
            continue
        for a2 in range(len(factors)):
            if a2 == a3 or factors[a2] != ("SU", 2):
                continue
            for sY in (1, -1):
                tbl = defaultdict(int)
                bad = False
                for r in left:
                    if a3 >= len(r["rep"]) or a2 >= len(r["rep"]):
                        bad = True
                        break
                    hid = 1
                    for k, dd in enumerate(r["rep"]):
                        if k not in (a3, a2):
                            hid *= abs(dd)
                    c3 = r["rep"][a3] if abs(r["rep"][a3]) != 1 else 1
                    tbl[(c3, abs(r["rep"][a2]), sY * dot(tY, r["p_sh"]))] += hid
                if not bad and dict(tbl) == dict(printed):
                    return (a3, a2, sY)
    return None


def role_census(spec, tY, a3, a2, sY0):
    """{label: {role, Y, kl, fp, mult}} + orientation diagnostics.  sY and sC are
    INDEPENDENT per-model conventions (5.21); ALL FOUR branches are computed and the
    physical one selected by net Q = +3.  Hidden multiplicity is carried, exactly as
    orbifolder's own SMConfig table carries it."""
    y0cache = {lab: (dot(tY, r["p_sh"]) if (r["rep"] and r["chir"] == "l" and r["p_sh"])
                     else None) for lab, r in spec.items()}
    hidcache = {}
    branches, nets = {}, {}
    for sY in (1, -1):
        for sC in (1, -1):
            out = {}
            for lab, r in spec.items():
                if r["rep"] is None or r["chir"] != "l" or r["p_sh"] is None:
                    continue
                if a3 >= len(r["rep"]) or a2 >= len(r["rep"]):
                    continue
                Y = sY * sY0 * y0cache[lab]
                hid = hidcache.get(lab)
                if hid is None:
                    hid = 1
                    for k, dd in enumerate(r["rep"]):
                        if k not in (a3, a2):
                            hid *= abs(dd)
                    hidcache[lab] = hid
                c3 = r["rep"][a3]
                c3 = (sC * c3) if abs(c3) != 1 else 1
                c2 = abs(r["rep"][a2])
                out[lab] = dict(role=SM_ROLES.get((c3, c2, Y)), Y=Y, kl=r["kl"],
                                fp=r["fp"], mult=hid,
                                hidden_singlet=hidden_singlet(r["rep"], a3, a2))
            branches[(sY, sC)] = out

    def net(o, pos, neg):
        return (sum(v["mult"] for v in o.values() if v["role"] == pos)
                - sum(v["mult"] for v in o.values() if v["role"] == neg))
    for key, o in branches.items():
        nets[key] = {"Q": net(o, "Q", "Qbar"), "uc": net(o, "uc", "ucbar"),
                     "dc": net(o, "dc", "dcbar"), "ec": net(o, "ec", "ecbar")}
    best = max(nets, key=lambda k: (nets[k]["Q"], nets[k]["uc"]))
    ties = [k for k in nets if (nets[k]["Q"], nets[k]["uc"])
            == (nets[best]["Q"], nets[best]["uc"])]
    diag = dict(sY_selected=best[0], sC_selected=best[1], sY_from_closure=sY0,
                nets={"sY=%d,sC=%d" % k: v for k, v in nets.items()},
                n_argmax_ties=len(ties),
                net_Q_is_plus_three=(nets[best]["Q"] == 3))
    return branches[best], diag


def leg_orbit_index(kl, fp, a):
    """Orbit-class index of a field's plane-a fixed point; None on a delocalised plane.

    orbifolder prints the space-group element's lattice part n_a in the SIX lattice-vector
    coordinates; the localisation in plane a is f = (1 - R_a)^{-1}(n|_a + sigma_a) mod 1 --
    the same construction the certified alphabet uses to label its orbit classes.  This
    geometry has NO roto-translation (has_rototranslation = False, certified), so
    sigma_a = 0; the term is nevertheless carried so the code is correct if it is ever
    pointed at a roto-translation geometry."""
    g = geom()
    kl = tuple(kl)
    if kl == (0, 0) or g.phi(a - 1, kl) == 0 or fp is None:
        return None
    A = one_minus(g.rot(a - 1, kl))
    Ai = minv(A)
    if Ai is None:
        return None
    i1, i2 = g.pairs[a - 1]
    s = g.sigma(a - 1, kl)
    f = modl(ap(Ai, (F(fp[i1]) + s[0], F(fp[i2]) + s[1])))
    for i, o in enumerate(g.orbits(a - 1, kl)):
        if f in o:
            return i
    return "OFF_LATTICE"


def read_gz(path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        return fh.read()


def write_json(name, obj):
    os.makedirs(DSW, exist_ok=True)
    with open(os.path.join(DSW, name), "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=1, default=str)

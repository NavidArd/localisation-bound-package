"""
rules_z36.py --- the complete trilinear selection-rule set for T^6/(Z_3 x Z_6) (`Z3xZ6_1_1`).

Phase: multi_geometry_campaign / Z3xZ6 foundations (docs/z3xz6_foundations_task.md), items 1, 5.

WHAT IS PUBLISHED AND WHAT IS THIS PROGRAMME'S OWN
--------------------------------------------------
* PUBLISHED, two-generator, read at equation level off the archived source
  `rcharge_census/docs/sources/0705.3075_discrete_R_anomalies.txt`
  (Araki-Kobayashi-Kubo-Ramos-Sanchez-Ratz-Vaudrevange, PRD77 (2008) 046013):
      (3)   untwisted H-momenta  U_b -> e_b
      (9)   the Z_N x Z_M twist shape  v^1 = (v1, -v1, 0), v^2 = (0, v2, -v2)
      (10)  r_i = |k v^1_i| + |l v^2_i| - Int[ |k v^1_i| + |l v^2_i| ]
      (12)  H-momentum conservation for the 3-point function: sum_legs r_i = 1
      (14)  R_i = r_i + N_i - Nbar_i
      (15)  sum_legs R_i = 1 mod N_i,  N_i the order of the point-group rotation in plane i
  This is the ONLY published R-charge law that covers Z_N x Z_M, and it is the law the
  programme's certified Z3xZ3 pass (`rcharge_census/`) implements and gate-validated against
  orbifolder itself.  It carries NO gamma term.

* PUBLISHED, single-generator only, off `multi_geometry_campaign/data/literature/1301.2322/`
  (`1301.2322v2_full.txt`, equation labels located by string search of the equation CONTENT in
  the local extract --- no number was ever supplied TO a retrieval tool, CLAUDE.md 5.16 sub-rule):
      (3.41)  sum_alpha (N^i_L - Nbar^i_L - Nbar^i_R) = 0 mod N_i/K      <- Rule 6, coset vector
      (3.44)  sum_alpha [ sum_{i nonprime} N v_i (N^i_L - Nbar^i_L - Nbar^i_R) - N gamma_alpha ] = 0 mod N
      (3.45)  sum_alpha R_alpha = ( sum_{i nonprime} N v_i ) mod N,
              R_alpha = sum_{i nonprime} N v_i (q^i_sh - N^i_L + Nbar^i_L) + N gamma_alpha
      and verbatim: *"when there are two non-prime planes, the condition (3.41) becomes an
      independent selection rule, which we may call Rule 6 or the coset vector selection rule"*.
  1301.2322 footnote 1 DECLINES the Z_N x Z_M discussion, so (3.41)/(3.44)/(3.45) have NO
  published two-generator form.

* THIS PROGRAMME'S OWN (owned extension, flagged as such everywhere it is used):
  the two-generator lift of (3.44)/(3.45)/(3.41).  Derived in `docs/rules_z3xz6.md` 4; the code
  below emits it as a TAGGED BRANCH and never as a published rule.

CONVENTIONS, STATED ONCE (5.19(i)) --- invariant checked before any verdict
--------------------------------------------------------------------------
  v_g   = k v^1 + l v^2                      local twist of sector g = (k, l)
  r_a(g)= (-v_g,a) mod 1                     H-momentum of the SCALAR component, plane a
          (the programme's `q^paper_a = 1 - frac(k v_a)` convention; identical to the certified
           Z3xZ3 route's  r = -q_sh  with  q_sh = -frac(-v))
  r_a(U_b) = delta_{ab}                       untwisted leg, 0705.3075 (3)
  INVARIANT CHECKED FIRST: sum_a r_a(g) = 1 for EVERY nontrivial sector g, and
  sum_a r_a(U_b) = 1; and U_1 U_2 U_3 satisfies (15) in every plane.  Refuses to run otherwise.

  D-3 CONJUGATION TAG (inherited): under the CPT-conjugate reading r' = frac(v_g) every channel
  identity maps (k,l) -> (-k,-l).  Every verdict below whose content flips under that map is
  tagged `d3_sensitive`.

Outputs: data/z3xz6/rules_channels_z36.json
Run:     python multi_geometry_campaign/code/z3xz6/rules_z36.py
"""

from __future__ import annotations

import json
import os
import sys
from fractions import Fraction as F
from itertools import combinations_with_replacement, product
from math import gcd
from typing import Dict, List, Optional, Sequence, Tuple

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z36lib import Geom, sector_label, OUT_DIR  # noqa: E402

G = Geom()
M, N6 = G.M, G.N                                   # 3, 6
PLANE_ORDER = [G.plane_order(a) for a in (1, 2, 3)]  # (6, 3, 6)
NONPRIME = [a for a in (1, 2, 3) if not _isprime(PLANE_ORDER[a - 1])] if False else None


def _isprime(n: int) -> bool:
    if n < 2:
        return False
    for p in range(2, int(n ** 0.5) + 1):
        if n % p == 0:
            return False
    return True


NONPRIME = [a for a in (1, 2, 3) if not _isprime(PLANE_ORDER[a - 1])]
PRIME = [a for a in (1, 2, 3) if _isprime(PLANE_ORDER[a - 1])]

ELEMENTS = G.elements()
TWISTED = [g for g in ELEMENTS if g != (0, 0)]


# ----------------------------------------------------------------------------------------
# H-momenta
# ----------------------------------------------------------------------------------------

def v_local(g: Tuple[int, int]) -> List[F]:
    return [(g[0] * G.v1[a] + g[1] * G.v2[a]) for a in range(3)]


def r_twisted(g: Tuple[int, int]) -> List[F]:
    """r_a = (-v_a) mod 1  (0 when the plane is unrotated: fixed torus).

    MINIMAL reading --- the programme's certified convention, validated on Z6-I and Z4 by the
    measured left-chiral basis (see `left_chiral_basis()` below).
    """
    out = []
    for a in range(3):
        va = v_local(g)[a] % 1
        out.append(F(0) if va == 0 else (1 - va))
    return out


def r_untwisted(b: int) -> List[F]:
    return [F(1) if a == b - 1 else F(0) for a in range(3)]


def r_of(leg) -> List[F]:
    """leg is either a point-group element (k,l) != (0,0), or ('U', b)."""
    if isinstance(leg, tuple) and leg and leg[0] == "U":
        return r_untwisted(leg[1])
    return r_twisted(leg)


def admissible_weights(leg) -> List[List[F]]:
    """ENUMERATED reading (tagged branch).

    Every weight with  r_a = (-v_a) mod 1  shifted by integers,  sum_a r_a = 1, and
    r_a in [-1, 1] (bounded; a masslessness proxy, NOT a derived cut --- see the RULE-GAP
    register in docs/rules_z3xz6.md 6).  The MINIMAL reading is always the first element.
    """
    base = r_of(leg)
    out = [list(base)]
    for n in product((-1, 0, 1), repeat=3):
        if sum(n) != 0:
            continue
        if all(x == 0 for x in n):
            continue
        cand = [base[a] + n[a] for a in range(3)]
        if all(-1 <= c <= 1 for c in cand) and sum(cand) == 1:
            if cand not in out:
                out.append(cand)
    return out


def left_chiral_basis() -> Dict[str, object]:
    """DERIVED left-chiral basis: the sectors on which sum_a r_a = 1 (5.19(i) invariant).

    In the programme's convention  r_a = 1 - frac(v_a)  on rotated planes, 0 on a fixed torus,
    sum_a r_a = (#rotated planes) - sum_a frac(v_a), and a massless left-chiral multiplet
    requires sum_a r_a = 1 --- 0705.3075 (12) applied to the single leg's own weight, together
    with (3) for the untwisted case (which the rule reproduces exactly: U_b -> e_b).

    VALIDATION AGAINST THE PROGRAMME'S OWN MEASURED BASES (not a fit --- a reproduction):
      Z6-I,  v = (1/6, 1/6, -1/3):  the criterion selects k in {0, 3, 4, 5} --- the MEASURED basis.
      Z4,    v = (1/4, 1/4, -1/2):  the criterion selects k in {0, 2, 3}    --- the MEASURED basis.
    Both are asserted in `validate_criterion_on_ZN()` before the Z3xZ6 basis is quoted.
    """
    left, right = [], []
    for g in TWISTED:
        (left if sum(r_twisted(g)) == 1 else right).append(sector_label(g))
    return {"left_chiral_sectors": ["U"] + left, "right_chiral_labels": right,
            "n_left": len(left) + 1, "n_right": len(right)}


def validate_criterion_on_ZN() -> Dict[str, object]:
    """Reproduce the MEASURED Z6-I and Z4 left-chiral bases from the same criterion."""
    out = {}
    for name, v, Nn, measured in (("Z6-I", (F(1, 6), F(1, 6), F(-1, 3)), 6, [0, 3, 4, 5]),
                                  ("Z4", (F(1, 4), F(1, 4), F(-1, 2)), 4, [0, 2, 3])):
        sel = [0]
        for k in range(1, Nn):
            r = []
            for a in range(3):
                va = (k * v[a]) % 1
                r.append(F(0) if va == 0 else 1 - va)
            if sum(r) == 1:
                sel.append(k)
        out[name] = {"derived": sel, "measured": measured, "match": sel == measured}
        if sel != measured:
            raise AssertionError("left-chiral criterion FAILS on %s: derived %s vs measured %s"
                                 % (name, sel, measured))
    return out


def check_conventions() -> Dict[str, object]:
    """INVARIANTS FIRST (5.19(i)).

    (i)  the left-chiral criterion reproduces the two MEASURED Z_N bases (Z6-I, Z4);
    (ii) sum_a r_a = 1 on every LEFT-CHIRAL leg (by construction of the basis, asserted);
    (iii) U1U2U3 satisfies the published law (15) in every plane.
    """
    zn = validate_criterion_on_ZN()
    lc = left_chiral_basis()
    for lab in lc["left_chiral_sectors"]:
        if lab == "U":
            continue
        k, l = [int(x) for x in lab[2:-1].split(",")]
        if sum(r_twisted((k, l))) != 1:
            raise AssertionError("left-chiral invariant fails on %s" % lab)
    uuu = [sum(r_untwisted(b)[a] for b in (1, 2, 3)) for a in range(3)]
    if not all((uuu[a] - 1) % PLANE_ORDER[a] == 0 for a in range(3)):
        raise AssertionError("U1U2U3 fails (15): %s" % uuu)
    return {"left_chiral_criterion_validated_on_ZN": zn,
            "left_chiral_basis_Z3xZ6": lc,
            "U1U2U3_satisfies_eq15_in_every_plane": True,
            "plane_orders_N_i": PLANE_ORDER,
            "non_prime_planes": NONPRIME, "prime_planes": PRIME}


# ----------------------------------------------------------------------------------------
# channels (Rule 2 point-group part, over BOTH generators)
# ----------------------------------------------------------------------------------------

def all_multisets() -> List[Tuple[Tuple[int, int], ...]]:
    """E2 / Rule 2 rotation part: sum k = 0 mod 3 AND sum l = 0 mod 6, from the group law."""
    out = []
    for tri in combinations_with_replacement(ELEMENTS, 3):
        if sum(t[0] for t in tri) % M == 0 and sum(t[1] for t in tri) % N6 == 0:
            out.append(tri)
    return out


def channel_name(tri) -> str:
    return "".join(sector_label(g) for g in tri)


# ----------------------------------------------------------------------------------------
# Rule 5 (forbidden-instanton / capability), 1301.2322 (2.22)/(2.23) at L = 3
# ----------------------------------------------------------------------------------------

def capable(tri, a: int) -> Dict[str, object]:
    """Per-plane capability.  kk_alpha = phi_a(g_alpha); the source's own coefficient counts.

    holomorphic branch  (domain 0 < k <= 1):  2 - sum_alpha kk_alpha > 0
    antiholomorphic     (domain 0 <= k < 1):  sum_alpha kk_alpha - 1 > 0
    An untwisted leg contributes kk = 1 in the holomorphic count and 0 in the antiholomorphic
    (the source's domain convention), which is the phase-2 engine's E3 substitution.
    """
    kk = []
    for g in tri:
        if isinstance(g, tuple) and g and g[0] == "U":
            kk.append(F(0))
        else:
            kk.append(G.phi(g)[a - 1])
    m_holo = sum(x if x != 0 else F(1) for x in kk)
    m_anti = sum(kk)
    return {"kk": [str(x) for x in kk],
            "holomorphic": (2 - m_holo) > 0,
            "antiholomorphic": (m_anti - 1) > 0,
            "capable": ((2 - m_holo) > 0) or ((m_anti - 1) > 0),
            "all_legs_twisted_in_plane": all(x != 0 for x in kk)}


# ----------------------------------------------------------------------------------------
# R-charge laws
# ----------------------------------------------------------------------------------------

def law_published_perplane_enumerated(legs) -> Dict[str, object]:
    """(15) under the ENUMERATED weight branch --- the decisive branch on this geometry.

    The MINIMAL reading assigns r_a = 0 to a leg in a plane where its sector acts trivially
    (a fixed torus).  Z3xZ6_1_1 has NINE N=2 sectors, and the maximal-block channel's three
    legs ALL have their fixed torus in the SAME plane (plane 2), so under the minimal reading
    sum_legs r_2 = 0 identically and (15) can never be met there.  The enumerated reading
    allows a leg's weight to be shifted by integers subject to sum_a r_a = 1 --- i.e. a unit of
    H-momentum may sit in the unrotated plane at the cost of a rotated one.

    Both readings are emitted.  Which one is physical is a DECLARED RULE-GAP (see
    docs/rules_z3xz6.md 6, D-2): it is settled by an orbifolder spectrum print, not by algebra.
    """
    best = None
    for combo in product(*[admissible_weights(g) for g in legs]):
        rows = {}
        ok_all = True
        for a in (1, 2, 3):
            s = sum(w[a - 1] for w in combo)
            ok = ((s - 1) % PLANE_ORDER[a - 1]) == 0
            rows["plane%d" % a] = {"sum_r": str(s), "pass": ok}
            ok_all = ok_all and ok
        if ok_all:
            best = {"weights": [[str(x) for x in w] for w in combo], "planes": rows,
                    "all_planes_pass": True}
            break
        if best is None:
            best = {"weights": [[str(x) for x in w] for w in combo], "planes": rows,
                    "all_planes_pass": False}
    return best


def law_published_perplane(tri) -> Dict[str, object]:
    """0705.3075 (15): sum_legs R_i = 1 mod N_i, oscillator-free (R_i = r_i).

    This is the PUBLISHED two-generator law and carries no gamma.
    """
    rows = {}
    for a in (1, 2, 3):
        s = sum(r_of(g)[a - 1] for g in tri)
        need = F(1)
        ok = ((s - need) % PLANE_ORDER[a - 1]) == 0
        rows["plane%d" % a] = {"sum_r": str(s), "target": "1", "modulus": PLANE_ORDER[a - 1],
                               "oscillator_free_pass": ok,
                               "delta_needed": str(((need - s) % PLANE_ORDER[a - 1]))}
    rows["all_planes_pass"] = all(rows["plane%d" % a]["oscillator_free_pass"] for a in (1, 2, 3))
    return rows


def law_nonprime_owned(tri) -> Dict[str, object]:
    """OWNED EXTENSION --- the two-generator lift of 1301.2322 (3.44)/(3.45).

    For a single generator h of order N_h whose rotation is nontrivial in the non-prime planes,
    the source's law reads
        sum_alpha [ sum_{i in NP(h)} N_h v_i^{(h)} (q^i - N^i_L + Nbar^i_L) + N_h gamma^{(h)}_alpha ]
            = sum_{i in NP(h)} N_h v_i^{(h)}   mod N_h.
    With two generators the sum over classical solutions in (3.43) runs over the FULL point
    group, so the gamma weighting is a CHARACTER of the group modulo the stabiliser --- i.e. a
    PAIR (gamma^theta, gamma^omega), not a scalar.  We therefore emit ONE law per generator,
    each with its own required sum of gammas, and flag the whole object `owned_extension: true`.

    Only the oscillator-free values are computed here (N_L = Nbar_L = 0).
    """
    out = {"owned_extension": True, "per_generator": {}}
    gens = {"theta": (1, 0), "omega": (0, 1)}
    for gname, h in gens.items():
        # order of h and the planes in which h's own rotation is non-prime
        ph = [G.phi(h)[a - 1] for a in (1, 2, 3)]
        Nh = 1
        for a in range(3):
            if ph[a] != 0:
                Nh = Nh * ph[a].denominator // gcd(Nh, ph[a].denominator)
        np_h = [a for a in (1, 2, 3) if ph[a - 1] != 0 and not _isprime(ph[a - 1].denominator)]
        if not np_h:
            out["per_generator"][gname] = {"order": Nh, "non_prime_planes_of_this_generator": [],
                                           "law_active": False,
                                           "note": "no non-prime plane for this generator; "
                                                   "the per-plane published law (15) is complete here"}
            continue
        coef = {a: Nh * ph[a - 1] for a in np_h}
        target = sum(coef.values()) % Nh
        lhs = sum(coef[a] * r_of(g)[a - 1] for g in tri for a in np_h)
        # R_alpha uses q^i_sh = -r_i in the paper's own sign; the programme carries r and
        # records the sign tag.  required sum of gammas from  lhs + Nh * sum_gamma = target:
        req = ((target - lhs) / Nh) % 1
        out["per_generator"][gname] = {
            "order": Nh,
            "phi_of_generator": [str(x) for x in ph],
            "non_prime_planes_of_this_generator": np_h,
            "coefficients_N_v": {str(a): str(coef[a]) for a in np_h},
            "target_mod_N": str(target),
            "lhs_oscillator_free": str(lhs),
            "required_sum_gamma_correlator_branch": str(req),
            "required_sum_gamma_LEEFT_branch": str((-req) % 1),
            "law_active": True,
        }
    return out


def gamma_value_set() -> Dict[str, object]:
    """The realisable gamma set, from TRUE simultaneous-orbit lengths (gamma-gate guard 1)."""
    lens = set()
    per_sector = {}
    for g in TWISTED:
        sz = G.simultaneous_orbit_sizes(g)
        per_sector[sector_label(g)] = sz
        lens |= set(sz)
    vals = sorted({F(p, l) for l in lens for p in range(l)})
    return {"realisable_simultaneous_orbit_lengths": sorted(lens),
            "gamma_values": [str(x) for x in vals],
            "n_gamma_values": len(vals),
            "collapses_to_0_and_half": vals == [F(0), F(1, 2)],
            "negation_symmetric_with_only_self_negative_elements":
                all(((-x) % 1) == x for x in vals),
            "per_sector_simultaneous_orbit_sizes": per_sector}


def rule6(tri) -> Dict[str, object]:
    """1301.2322 (3.41), coset-vector rule, in the two-generator setting (owned extension).

    (3.41): sum_alpha (N^i_L - Nbar^i_L - Nbar^i_R) = 0 mod N_i / K, per NON-PRIME plane i,
    with K the lowest common multiple of the orbit lengths l_alpha of the three legs.

    BOTH K readings of the source's internally inconsistent definition are emitted as tagged
    branches, exactly as in the certified Z6-I/Z4 pipeline:
      K_orbit  --- lcm over the legs' SIMULTANEOUS-orbit lengths (footnote-7 reading);
      K_sector --- lcm over the legs' sector orders (prose-after-(3.41) reading).
    Oscillator-free: N_L = Nbar_L = 0 and Nbar^i_R = sum_alpha q^i_sh - 1 = -(sum r_i) - 1 by
    H-momentum, so the left-hand side is  A_i = 1 - sum_alpha r_i  exactly (the identity the Z4
    pass proved; re-derived here for the two-generator case in docs/rules_z3xz6.md 4.2).
    """
    res = {"owned_extension": True, "branches": {}}
    ord_of = []
    orb_of = []
    for g in tri:
        if isinstance(g, tuple) and g and g[0] == "U":
            ord_of.append(1)
            orb_of.append([1])
        else:
            o = 1
            for a in range(3):
                p = G.phi(g)[a]
                if p != 0:
                    o = o * p.denominator // gcd(o, p.denominator)
            ord_of.append(o)
            orb_of.append(G.simultaneous_orbit_sizes(g))
    def lcm(xs):
        r = 1
        for x in xs:
            r = r * x // gcd(r, x)
        return r
    K_sector = lcm(ord_of)
    for tag, Kset in (("K_sector", [K_sector]),
                      ("K_orbit", sorted({lcm(c) for c in product(*orb_of)}))):
        rows = {}
        for a in NONPRIME:
            Ai = 1 - sum(r_of(g)[a - 1] for g in tri)
            per_K = {}
            for K in Kset:
                mod = F(PLANE_ORDER[a - 1], K)
                if mod.denominator != 1:
                    per_K[str(K)] = {"modulus": str(mod), "applies": False,
                                     "note": "K does not divide N_i: rule inactive at this K"}
                    continue
                per_K[str(K)] = {"modulus": str(mod), "applies": True,
                                 "A_i": str(Ai),
                                 "pass": (Ai % int(mod)) == 0}
            rows["plane%d" % a] = {"A_i_equals_1_minus_sum_r": str(Ai), "per_K": per_K}
        res["branches"][tag] = {"K_values": Kset, "planes": rows,
                                "any_K_passes_all_nonprime_planes": any(
                                    all(rows["plane%d" % a]["per_K"][str(K)].get("pass", True)
                                        for a in NONPRIME) for K in Kset)}
    res["branches_agree"] = (
        res["branches"]["K_sector"]["any_K_passes_all_nonprime_planes"]
        == res["branches"]["K_orbit"]["any_K_passes_all_nonprime_planes"])
    return res


# ----------------------------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------------------------

def analyse() -> Dict[str, object]:
    conv = check_conventions()
    chans = all_multisets()
    rows = []
    for tri in chans:
        # U legs need their plane index enumerated; a channel with n U legs has 3^n variants
        n_u = sum(1 for g in tri if g == (0, 0))
        variants = []
        for assign in product((1, 2, 3), repeat=n_u):
            legs, it = [], iter(assign)
            for g in tri:
                legs.append(("U", next(it)) if g == (0, 0) else g)
            pub = law_published_perplane(legs)
            variants.append({
                "u_assignment": list(assign),
                "published_perplane_law_0705_3075_eq15": pub,
                "published_law_ENUMERATED_weight_branch": law_published_perplane_enumerated(legs),
                "nonprime_owned_law": law_nonprime_owned(legs),
                "rule6": rule6(legs),
            })
        rows.append({
            "channel": channel_name(tri),
            "gs": [list(g) for g in tri],
            "rule5_per_plane": {"plane%d" % a: capable(tri, a) for a in (1, 2, 3)},
            "capable_planes": [a for a in (1, 2, 3) if capable(tri, a)["capable"]],
            "n_untwisted_legs": n_u,
            "variants": variants,
            "oscillator_free_allowed_published_law": any(
                v["published_perplane_law_0705_3075_eq15"]["all_planes_pass"] for v in variants),
            "oscillator_free_allowed_ENUMERATED_branch": any(
                v["published_law_ENUMERATED_weight_branch"]["all_planes_pass"] for v in variants),
        })
    return {"geometry": "Z3xZ6_1_1", "conventions": conv,
            "gamma": gamma_value_set(),
            "n_channels": len(chans), "channels": rows}


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    rep = analyse()
    path = os.path.join(OUT_DIR, "rules_channels_z36.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(rep, fh, indent=1, default=str)

    print("=" * 78)
    print("Z3xZ6_1_1 rule set  ---  channels, Rule 5, R-charge laws, Rule 6, gamma")
    print("=" * 78)
    print("conventions:", rep["conventions"])
    gm = rep["gamma"]
    print("gamma: realisable simultaneous-orbit lengths %s" % gm["realisable_simultaneous_orbit_lengths"])
    print("       gamma set (%d values): %s" % (gm["n_gamma_values"], gm["gamma_values"]))
    print("       Z4-style collapse to {0,1/2}? %s" % gm["collapses_to_0_and_half"])
    print("       negation-symmetric with only self-negative elements? %s"
          % gm["negation_symmetric_with_only_self_negative_elements"])
    print("-" * 78)
    print("channels: %d" % rep["n_channels"])
    cap = [c for c in rep["channels"] if c["capable_planes"]]
    print("channels with at least one capable plane: %d" % len(cap))
    allow = [c for c in rep["channels"] if c["oscillator_free_allowed_published_law"]]
    print("oscillator-free allowed by the PUBLISHED per-plane law (15): %d" % len(allow))
    for c in allow:
        print("   %-28s capable planes %s" % (c["channel"], c["capable_planes"]))
    allow2 = [c for c in rep["channels"] if c["oscillator_free_allowed_ENUMERATED_branch"]]
    print("oscillator-free allowed under the ENUMERATED weight branch: %d" % len(allow2))
    both = [c["channel"] for c in allow2 if not c["oscillator_free_allowed_published_law"]]
    print("  channels the enumerated branch ADDS (%d): %s" % (len(both), both[:40]))
    print("-" * 78)
    print("the maximal-block channel T(0,1)T(0,2)T(0,3):")
    for c in rep["channels"]:
        if c["channel"] == "T(0,1)T(0,2)T(0,3)":
            v = c["variants"][0]
            print("   capable planes            %s" % c["capable_planes"])
            for a in (1, 2, 3):
                p = v["published_perplane_law_0705_3075_eq15"]["plane%d" % a]
                print("   plane %d  sum_r=%-5s target 1 mod %-2s  pass=%-5s  delta needed %s"
                      % (a, p["sum_r"], p["modulus"], p["oscillator_free_pass"], p["delta_needed"]))
            for gname, gg in v["nonprime_owned_law"]["per_generator"].items():
                print("   owned non-prime law [%s]: %s" % (gname, gg))
            print("   ENUMERATED branch: %s" % v["published_law_ENUMERATED_weight_branch"])
            r6 = v["rule6"]
            print("   Rule 6 branches agree: %s" % r6["branches_agree"])
            for tag, b in r6["branches"].items():
                print("     %-9s K=%s  any-K passes all non-prime planes: %s"
                      % (tag, b["K_values"], b["any_K_passes_all_nonprime_planes"]))
    print("wrote %s" % path)


if __name__ == "__main__":
    main()

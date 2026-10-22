"""z3xz6_sweep -- STEP 1: the OWED hypercharge dual-route verification.

Gate requirement 4 of `results/gate_z3xz6_foundations.md`: *"Hypercharge dual-route on
>= 5 models, with the anomaly conditions, before any model-level flavour verdict -- still
OWED"*.  The foundations pass could not run it (no orbifolder run); this file discharges
it, and if it fails on ANY model the sweep STOPS, because the catalogue's whole
defining-condition battery would then be compromised.

CLAUDE.md 5.21 is the binding lesson, verbatim in effect: *"closing a DATA-GAP requires
verifying EVERY defining condition of the object, enumerated first, not the convenient
subset"*.  The Z3xZ3 failure it records was a fitted hypercharge reproducing exact SM
values on every gauge structure -- check PASS -- while being GS-ANOMALOUS
(`t_Y . t_anom = -3`, `Tr Y = -36`).  The conditions are therefore enumerated FIRST:

  C1  t_Y lies in the rational span of the model's own printed SU(5) origin (LITERAL
      membership, exact Gaussian elimination -- not the statement that the C1+C2 solution
      space is one-dimensional, which is kept as a SEPARATE named check; that conflation
      was gate ZF-1 defect D-7 on Z6-I);
  C2  t_Y is orthogonal to every SU(3)_C / SU(2)_L simple root;
  C3  the GUT normalisation k_Y = 2 t_Y . t_Y = 5/3;
  C4  t_Y . t_anom = 0 -- GS-orthogonality, the condition the Z3xZ3 fit silently broke;
  C5  Tr Y = 0 over the full SMConfig spectrum;
  C6  Tr Y^3 = 0 (the cubic hypercharge anomaly);
  C7  every SMConfig representation carries a Y in the SM set;
  C8  net generations: the chiral excess of the five SM multiplets is exactly three --
      an INDEPENDENT defining condition, being a property of the spectrum rather than of
      the generator.

ROUTE 1 (tool).     orbifolder's own printed `Hypercharge generator`.
ROUTE 2 (offline).  t_Y RE-DERIVED over Q from the printed SU(5)-origin roots by
                    C1 + C2 + C3 alone -- exact linear algebra, no fit, no least squares,
                    no tolerance.  5.21(ii): a fitted representation of an object with a
                    known defining constraint must impose or verify that constraint IN
                    THE FIT.  There is no fit here; the solution space is one-dimensional
                    and the constraint set fixes t_Y UP TO SIGN.
The routes are compared VECTOR BY VECTOR, exactly.

SIGN (gate ZF-1 defect D-6, inherited).  C1+C2+C3 fix t_Y only up to an overall sign, so
Route 2 cannot produce a signed vector from its own constraints.  Its sign is fixed by an
INTERNAL rule (first non-zero component positive) that NEVER looks at Route 1, and two
numbers are reported: `routes_agree_up_to_sign` (the physics claim -- Y -> -Y is
particle <-> antiparticle) and `routes_agree_including_sign` (a strictly weaker
bookkeeping coincidence, NOT expected to be universal).  No claim here rests on the second.

C4's PROVENANCE -- an improvement on the Z6-I pass, which RECORDED this as unrepaired
(defect D-8).  There, `t_anom` was TAKEN to be the first printed U(1).  Here it is
IDENTIFIED from orbifolder's own printed `tr Q_i` block: the anomalous direction is the
U(1) whose gravitational trace is non-zero, and the check is reported as
non-discriminating whenever `t_Y` and `t_anom` have disjoint E8 support (in which case
C4 is automatic and proves nothing).

Run: python multi_geometry_campaign/code/z3xz6_sweep/hypercharge_z36.py [N_MODELS]
"""

from __future__ import annotations

import json
import math
import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z36sw_lib as Z                                                  # noqa: E402

SM_SET = {F(1, 6), F(-1, 6), F(2, 3), F(-2, 3), F(1, 3), F(-1, 3),
          F(1, 2), F(-1, 2), F(1), F(-1), F(0)}
CORE_KEYS = ("C1_pass", "C2_orthogonal_to_SM_roots", "C3_kY_is_5_3", "C4_pass",
             "C5_pass", "C6_pass", "C7_values_in_SM_set", "C8_pass")


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def nullspace_in_span(basis, constraints):
    """x = sum c_a basis[a] with x.r = 0 for all r in constraints; exact over Q."""
    n = len(basis)
    rows = [[dot(b, r) for b in basis] for r in constraints]
    piv, where = 0, [-1] * n
    for col in range(n):
        sel = None
        for i in range(piv, len(rows)):
            if rows[i][col] != 0:
                sel = i
                break
        if sel is None:
            continue
        rows[piv], rows[sel] = rows[sel], rows[piv]
        f = rows[piv][col]
        rows[piv] = [x / f for x in rows[piv]]
        for i in range(len(rows)):
            if i != piv and rows[i][col] != 0:
                g = rows[i][col]
                rows[i] = [a - g * b for a, b in zip(rows[i], rows[piv])]
        where[col] = piv
        piv += 1
    free = [c for c in range(n) if where[c] == -1]
    if len(free) != 1:
        return None, len(free)
    fc = free[0]
    coef = [F(0)] * n
    coef[fc] = F(1)
    for c in range(n):
        if where[c] != -1:
            coef[c] = -rows[where[c]][fc]
    x = [F(0)] * 16
    for a, c in enumerate(coef):
        if c:
            x = [xi + c * bi for xi, bi in zip(x, basis[a])]
    return x, 1


def in_span(x, basis):
    """C1 LITERALLY: is x in the rational span of `basis`?  Exact elimination."""
    if not basis:
        return None
    n, dim = len(basis), len(x)
    rows = [[basis[a][i] for a in range(n)] + [x[i]] for i in range(dim)]
    piv = 0
    for col in range(n):
        sel = None
        for i in range(piv, dim):
            if rows[i][col] != 0:
                sel = i
                break
        if sel is None:
            continue
        rows[piv], rows[sel] = rows[sel], rows[piv]
        f = rows[piv][col]
        rows[piv] = [v / f for v in rows[piv]]
        for i in range(dim):
            if i != piv and rows[i][col] != 0:
                g = rows[i][col]
                rows[i] = [a - g * b for a, b in zip(rows[i], rows[piv])]
        piv += 1
    for i in range(dim):
        if all(v == 0 for v in rows[i][:n]) and rows[i][n] != 0:
            return False
    return True


def fix_sign_internal(v):
    """Route 2's OWN sign rule; depends on nothing outside `v` (defect D-6)."""
    for x in v:
        if x != 0:
            return [-y for y in v] if x < 0 else list(v)
    return list(v)


def net_generations(sm):
    """C8.  Chiral excess of the five SM multiplets, computed on CONJUGATE PAIRS so the
    model's overall Y-sign convention cannot change the answer.  The pair itself is the
    physical one (writing d^c as (3,1,+1/3) names a DIFFERENT pair and returns 0 on every
    model -- a 5.19 conjugation slip caught by the census reading 0 where the others
    read 3)."""
    tot = {}
    for r in sm:
        key = (r["d3"], r["d2"], r["Y"])
        tot[key] = tot.get(key, 0) + r["mult"]
    pairnet = {}
    for (d3, d2, Y), n in tot.items():
        conj = (-d3 if abs(d3) != 1 else 1, d2, -Y)
        if conj == (d3, d2, Y):
            continue
        key = tuple(sorted([(d3, d2, Y), conj], key=repr))
        pairnet[key] = abs(n - tot.get(conj, 0))
    want = {"q": (3, 2, F(1, 6)), "u^c": (-3, 1, F(-2, 3)), "d^c": (-3, 1, F(1, 3)),
            "l": (1, 2, F(-1, 2)), "e^c": (1, 1, F(1))}

    def census(s):
        o = {}
        for name, (a3, a2, aY) in want.items():
            aY = s * aY
            k = tuple(sorted([(a3, a2, aY), (-a3 if abs(a3) != 1 else 1, a2, -aY)],
                             key=repr))
            o[name] = pairnet.get(k, 0)
        return o
    cp, cm = census(1), census(-1)
    assert cp["l"] == cm["l"] and cp["e^c"] == cm["e^c"], (cp, cm)
    s = 1 if cp["q"] >= cm["q"] else -1
    out = dict(cp if s == 1 else cm)
    out["_Y_sign_convention"] = s
    return out


def identify_anomalous(gauge):
    """C4's t_anom, IDENTIFIED from the printed `tr Q_i` block rather than assumed to be
    the first U(1) (the Z6-I pass RECORDED this gap as defect D-8; it is closed here).
    Returns (index_1based, t_anom, provenance)."""
    tr = gauge.get("trQ") or {}
    nz = [i for i, v in tr.items() if v != 0]
    if len(nz) == 1 and nz[0] - 1 < len(gauge["u1"]):
        return nz[0], gauge["u1"][nz[0] - 1], "identified from tr Q_i (unique non-zero)"
    if not gauge["u1"]:
        return None, None, "no U(1) generators printed"
    if not nz:
        return None, None, "anomaly-free: every printed tr Q_i vanishes"
    return (nz[0], gauge["u1"][nz[0] - 1],
            "several non-zero tr Q_i %s; first taken -- FLAGGED" % nz)


def analyse_model(label, block):
    g = Z.parse_gauge(block)
    spec = Z.parse_spectrum(block)
    r = {"model": label, "gauge_TestConfig": g["gauge_test"],
         "gauge_SMConfig": g["gauge_sm"], "n_u1": len(g["u1"]),
         "n_simple_roots": len(g["roots"]), "n_su5_roots": len(g["su5"]),
         "n_hypercharge_blocks": g["n_hyper_blocks"], "n_fields": len(spec),
         "n_SMConfig_rows": len(g["sm"])}
    tY1, su5 = g["tY"], g["su5"]
    smroots = [a for a in g["roots"] if a in su5]
    r["n_sm_roots_identified"] = len(smroots)
    tY2, nfree = (nullspace_in_span(su5, smroots) if (su5 and smroots) else (None, None))
    if tY2 is not None:
        n2 = dot(tY2, tY2)
        if n2 != 0:
            s2 = F(5, 6) / n2
            rn, rd = math.isqrt(s2.numerator), math.isqrt(s2.denominator)
            if rn * rn == s2.numerator and rd * rd == s2.denominator:
                tY2 = fix_sign_internal([F(rn, rd) * x for x in tY2])
            else:
                r["C3_scale_irrational"] = str(s2)
    r["route1_tY"] = [str(x) for x in tY1] if tY1 else None
    r["route2_tY"] = [str(x) for x in tY2] if tY2 else None
    r["route2_sign_convention"] = "internal: first non-zero component positive"
    upto = bool(tY1 and tY2 and (tY1 == tY2 or tY1 == [-x for x in tY2]))
    r["routes_agree_up_to_sign"] = upto
    r["routes_agree_including_sign"] = bool(tY1 and tY2 and tY1 == tY2)

    c = {}
    if tY1:
        c["C1_pass"] = in_span(tY1, su5) if su5 else None
        c["C1_solution_space_dim_1"] = (nfree == 1) if (su5 and smroots) else None
        c["C2_orthogonal_to_SM_roots"] = all(dot(tY1, a) == 0 for a in smroots)
        c["C3_kY"] = str(2 * dot(tY1, tY1))
        c["C3_kY_is_5_3"] = (2 * dot(tY1, tY1) == F(5, 3))
        ai, tan, prov = identify_anomalous(g)
        c["C4_tanom_index"] = ai
        c["C4_tanom_provenance"] = prov
        c["C4_tY_dot_tanom"] = str(dot(tY1, tan)) if tan else None
        c["C4_pass"] = (dot(tY1, tan) == 0) if tan else None
        if tan:
            sy = {0 if i < 8 else 1 for i, x in enumerate(tY1) if x != 0}
            sa = {0 if i < 8 else 1 for i, x in enumerate(tan) if x != 0}
            c["C4_disjoint_E8_support"] = not (sy & sa)
    if g["sm"]:
        trY = sum(F(x["mult"]) * abs(x["d3"]) * x["d2"] * x["Y"] for x in g["sm"])
        trY3 = sum(F(x["mult"]) * abs(x["d3"]) * x["d2"] * x["Y"] ** 3 for x in g["sm"])
        c["C5_TrY"] = str(trY)
        c["C5_pass"] = trY == 0
        c["C6_TrY3"] = str(trY3)
        c["C6_pass"] = trY3 == 0
        c["C7_values_in_SM_set"] = all(x["Y"] in SM_SET for x in g["sm"])
        gens = net_generations(g["sm"])
        c["C8_net_generations"] = {k: v for k, v in gens.items() if k[0] != "_"}
        c["C8_Y_sign_convention"] = gens["_Y_sign_convention"]
        c["C8_pass"] = all(v == 3 for k, v in gens.items() if k[0] != "_")
    r["checks"] = c
    assert len(CORE_KEYS) == 8, "the 'all eight' aggregate must contain eight"
    r["conditions_aggregated"] = list(CORE_KEYS)
    r["all_conditions_pass"] = all(c.get(k) is True for k in CORE_KEYS)
    return r


def pick_models(n):
    """A SPREAD across the catalogue plus the three specimens the foundations and the gate
    used -- never the first n, whose provenance order is shared."""
    base = ["MSSM0", "MSSM24", "MSSM26"]
    step = max(1, 4612 // max(1, n - len(base)))
    out = list(base)
    for i in range(0, 4612, step):
        lab = "MSSM%d" % i
        if lab not in out:
            out.append(lab)
        if len(out) >= n:
            break
    return out


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 25
    want = set(pick_models(n))
    rep = {"phase": "z3xz6_sweep STEP 1 (gate requirement 4)",
           "scope": Z.CATALOGUE_TAG,
           "conditions": ["C1 in SU(5) Cartan (literal membership)",
                          "C2 orthogonal to SM simple roots",
                          "C3 k_Y = 2 t_Y.t_Y = 5/3", "C4 t_Y.t_anom = 0 (GS)",
                          "C5 Tr Y = 0", "C6 Tr Y^3 = 0", "C7 SM values",
                          "C8 net generations = 3"],
           "models": []}
    seen = set()
    for fn in sorted(os.listdir(Z.RAW_CENSUS)):
        if not fn.endswith(".gz") or "PILOT" in fn:
            continue
        txt = Z.read_raw(os.path.join(Z.RAW_CENSUS, fn))
        for lab, block in Z.split_models(txt):
            if lab in want and lab not in seen:
                seen.add(lab)
                rep["models"].append(analyse_model(lab, block))
        if seen >= want:
            break
    rep["models"].sort(key=lambda r: int(r["model"][4:]))
    m = rep["models"]
    rep["summary"] = {
        "n_models": len(m),
        "routes_agree_up_to_sign": sum(1 for r in m if r["routes_agree_up_to_sign"]),
        "routes_agree_including_sign":
            sum(1 for r in m if r["routes_agree_including_sign"]),
        "all_eight_conditions_pass": sum(1 for r in m if r["all_conditions_pass"]),
        "three_generations": sum(1 for r in m if r["checks"].get("C8_pass")),
        "C4_identified_from_trQ":
            sum(1 for r in m if "identified" in
                (r["checks"].get("C4_tanom_provenance") or "")),
        "C4_non_discriminating_disjoint_support":
            sum(1 for r in m if r["checks"].get("C4_disjoint_E8_support")),
        "route2_sign_convention":
            "internal: first non-zero component positive (Route 2 never reads Route 1)"}
    dest = os.path.join(Z.DATA, "hypercharge_dual_route_z36.json")
    json.dump(rep, open(dest, "w"), indent=1, default=str)

    print("Z3xZ6_1_1 hypercharge -- dual route + full defining-condition battery (5.21)")
    print("%-10s %-7s %-6s %-9s %-11s %-6s %-12s %s"
          % ("model", "routes", "k_Y", "tY.tanom", "TrY/TrY^3", "SMval", "generations",
             "ALL"))
    print("-" * 92)
    for r in m:
        c = r["checks"]
        g = c.get("C8_net_generations", {})
        print("%-10s %-7s %-6s %-9s %-11s %-6s %-12s %s"
              % (r["model"], "AGREE" if r["routes_agree_up_to_sign"] else "DIFFER",
                 c.get("C3_kY"), c.get("C4_tY_dot_tanom"),
                 "%s/%s" % (c.get("C5_TrY"), c.get("C6_TrY3")),
                 "ok" if c.get("C7_values_in_SM_set") else "NO",
                 ",".join(str(g.get(k, "?"))
                          for k in ("q", "u^c", "d^c", "l", "e^c")),
                 "PASS" if r["all_conditions_pass"] else "FAIL"))
    s = rep["summary"]
    print("-" * 92)
    print("models %d | routes agree up to sign %d | sign also matches %d | all eight "
          "C1..C8 %d | 3 generations %d" %
          (s["n_models"], s["routes_agree_up_to_sign"],
           s["routes_agree_including_sign"], s["all_eight_conditions_pass"],
           s["three_generations"]))
    print("C4: t_anom identified from tr Q_i on %d/%d models; non-discriminating "
          "(disjoint E8 support) on %d" %
          (s["C4_identified_from_trQ"], s["n_models"],
           s["C4_non_discriminating_disjoint_support"]))
    print("->", dest)
    assert s["n_models"] >= 5, "gate requirement 4 needs >= 5 models"
    if s["all_eight_conditions_pass"] != s["n_models"]:
        bad = [r["model"] for r in m if not r["all_conditions_pass"]]
        raise SystemExit("STOP: the defining-condition battery FAILS on %s -- the "
                         "catalogue's condition battery is compromised and no verdict "
                         "row may be produced (gate requirement 4)." % bad)
    print("\nSTEP 1 PASS: every defining condition holds on every model tested.")


if __name__ == "__main__":
    main()

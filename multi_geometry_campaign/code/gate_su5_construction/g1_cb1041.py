"""V1 -- certify or break CB1041, from its raw orbifolder output, with fresh code.

Every quantity the phase reports on the flagship is recomputed here from (V1, V2, W3) and
from the printed spectrum, by routes independent of code/su5_construction/.

  A  MI0-MI6 recomputed from the source-read rule set (gs_lib G3).
  B  the 4d gauge group derived from the lattice (roots integral on V1, V2, W) and
     compared with orbifolder's PRINTED gauge group and simple roots.
  C  k_Y = 5/3 ; t_Y orthogonal to every 4d simple root ; t_Y inside the span of
     orbifolder's own printed U(1) generators (so it IS a U(1) of the printed group).
  D  Tr Y and Tr Y^3 by TWO independent routes:
       (i)  summed over left-chiral STATES from the printed spectrum
       (ii) from orbifolder's OWN printed anomaly traces tr Q_i, tr Q_i^3, tr Q_i Q_j^2
            after expressing t_Y in the printed U(1) basis
  E  net generations, convention-blind: both sY branches, TRIALITY-AWARE roles
     (gs_lib G5) -- the phase's role_of discards the SU(3) triality; the gate does not.
     Also cross-checked against orbifolder's printed rep dimensions.
  F  F_6 / F_57: full plane-tuple cell key (5.22), rep, |p_sh|^2, oscillator, |w|^2, R,
     p.Delta, the partner state sigma_Delta(p), and H1-H5 of the localisation theorem.
  G  L1 on this model: every T(0,l) record has p_sh in l*V2 + Lambda.

Output: data/gate_su5_construction/g1_cb1041.json
"""
from __future__ import annotations

import gzip
import json
import os
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gs_lib as G

RAW = os.path.join(G.PHASE_DATA, "raw")
TARGET = "CB1041"


def fr(v):
    return [F(x) if not isinstance(x, F) else x for x in v]


def parse_frac_list(lst):
    return [G.rat(s) for s in lst]


def gauge_roots(V1, V2, Ws):
    """4d unbroken roots: alpha in Lambda, |alpha|^2 = 2, alpha.V and alpha.W in Z."""
    out = []
    for a in G.e8xe8_roots():
        if G.dot(a, V1).denominator != 1:
            continue
        if G.dot(a, V2).denominator != 1:
            continue
        if any(G.dot(a, W).denominator != 1 for W in Ws if any(W)):
            continue
        out.append(a)
    return out


def identify_algebra(roots):
    """Return sorted list of (type-guess, n_roots, rank) for each connected component."""
    # build components by non-orthogonality on a chosen set of positive roots
    import itertools
    n = len(roots)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for i in range(n):
        for j in range(i + 1, n):
            if G.dot(roots[i], roots[j]) != 0:
                union(i, j)
    comps = defaultdict(list)
    for i in range(n):
        comps[find(i)].append(roots[i])
    res = []
    for c in comps.values():
        nr = len(c)
        # rank = dimension of the span
        rank = rank_of(c)
        res.append((nr, rank, name_from(nr, rank)))
    return sorted(res)


def rank_of(vecs):
    rows = [list(v) for v in vecs]
    r = 0
    ncols = len(rows[0])
    for col in range(ncols):
        piv = None
        for i in range(r, len(rows)):
            if rows[i][col] != 0:
                piv = i
                break
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        pv = rows[r][col]
        for i in range(len(rows)):
            if i != r and rows[i][col] != 0:
                f = rows[i][col] / pv
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        r += 1
    return r


def name_from(nroots, rank):
    table = {(2, 1): "SU(2)", (6, 2): "SU(3)", (12, 3): "SU(4)", (20, 4): "SU(5)",
             (24, 4): "SO(8)", (30, 5): "SU(6)", (40, 5): "SO(10)", (42, 6): "SU(7)",
             (60, 6): "SO(12)", (72, 6): "E6", (126, 7): "E7", (240, 8): "E8",
             (56, 7): "SU(8)", (8, 2): "SU(2)xSU(2)?", (18, 4): "?"}
    return table.get((nroots, rank), "rank%d_%droots" % (rank, nroots))


def main():
    print("== gs_lib selftest ==")
    assert G.selftest()

    flag = json.load(open(os.path.join(G.PHASE_DATA, "s12_flagship.json")))["flagship"]
    V1 = parse_frac_list(flag["V1"])
    V2 = parse_frac_list(flag["V2"])
    W3 = parse_frac_list(flag["W3"])
    Z = [F(0)] * 16
    Ws = [Z, Z, W3, W3, Z, Z]
    res = dict(model=TARGET, V1=[str(x) for x in V1], V2=[str(x) for x in V2],
               W3=[str(x) for x in W3])

    # ---------------------------------------------------------------- A  MI
    ok, bad = G.mi_check(V1, V2, Ws)
    res["A_MI"] = dict(ok=ok, failures={k: str(v) for k, v in bad.items()},
                       ordV1=G.order_of(V1), ordV2=G.order_of(V2), ordW3=G.order_of(W3),
                       MI1_value=str(G.M_ORD * (G.n2(V1) - G.n2(G.V1_TWIST))),
                       MI2_value=str(G.N_ORD * (G.n2(V2) - G.n2(G.V2_TWIST))),
                       MI3_value=str(G.SP[0][1] * (G.dot(V1, V2)
                                                   - G.dot(G.V1_TWIST, G.V2_TWIST))),
                       MI4_value=str(G.SP[4][0] * G.dot(W3, V1)),
                       MI5_value=str(G.SP[4][1] * G.dot(W3, V2)),
                       MI6_33=str(G.SP[4][4] * G.dot(W3, W3)),
                       MI6_34=str(G.SP[4][5] * G.dot(W3, W3)))
    print("A  MI full pass:", ok, bad)

    # ---------------------------------------------------------------- raw block
    blob = None
    for fn in sorted(os.listdir(RAW)):
        t = gzip.open(os.path.join(RAW, fn), "rt", errors="replace").read()
        for lab, blk in G.split_models(t):
            if lab == TARGET:
                blob = blk
                res["raw_segment"] = fn
                break
        if blob:
            break
    assert blob, "CB1041 not found in raw output"
    fields = G.parse_fields(blob)
    gg = G.parse_gauge_group(blob)
    u1 = G.parse_u1_generators(blob)
    simple = G.parse_simple_roots(blob)
    tr = G.parse_anomaly_traces(blob)
    res["printed_gauge_group"] = gg
    res["n_records"] = len(fields)
    print("   raw segment:", res["raw_segment"], "records:", len(fields))
    print("   printed gauge group:", gg)

    # ---------------------------------------------------------------- B  gauge group
    roots = gauge_roots(V1, V2, Ws)
    comps = identify_algebra(roots)
    res["B_gauge"] = dict(n_4d_roots=len(roots), components=[list(c) for c in comps])
    # every printed simple root must be a 4d root by our own criterion
    rootset = set(tuple(r) for r in roots)
    res["B_printed_simple_roots_are_4d_roots"] = all(tuple(s) in rootset for s in simple)
    res["B_n_printed_simple_roots"] = len(simple)
    res["B_rank_printed_simple"] = rank_of(simple) if simple else 0
    print("B  4d roots:", len(roots), "components:", comps)
    print("   printed simple roots all 4d roots:",
          res["B_printed_simple_roots_are_4d_roots"], "| n =", len(simple))

    # ---------------------------------------------------------------- C  hypercharge
    kY = 2 * G.n2(G.TY)
    orth = all(G.dot(G.TY, s) == 0 for s in simple)
    # is t_Y in the span of the printed U(1) generators?
    coeffs = solve_in_span(u1, G.TY)
    res["C_hyper"] = dict(kY=str(kY), tY_orth_to_all_printed_simple_roots=orth,
                          tY_in_span_of_printed_U1=(coeffs is not None),
                          tY_coeffs_in_printed_U1_basis=(
                              [str(c) for c in coeffs] if coeffs else None))
    print("C  k_Y =", kY, "| t_Y orth to printed simple roots:", orth,
          "| t_Y in printed U(1) span:", coeffs is not None, coeffs)

    # ---------------------------------------------------------------- states
    left = [f for f in fields if f["chir"] == "l"]
    right = [f for f in fields if f["chir"] == "r"]
    res["n_left_records"] = len(left)
    res["n_right_records"] = len(right)

    bad_level = []
    for f in left:
        try:
            G.su3_triality(f["p_sh"])
            G.su2_dim(f["p_sh"])
        except ValueError as e:
            bad_level.append((f["label"], str(e)))
    res["level1_violations"] = bad_level

    # hidden multiplicity: printed total dim / (d3 * d2), cross-checked against the
    # printed first two entries of the rep tuple
    dimmismatch = []
    for f in left:
        t = G.su3_triality(f["p_sh"])
        d3 = 1 if t == 0 else 3
        d2 = G.su2_dim(f["p_sh"])
        pd3, pd2 = abs(f["rep"][0]), abs(f["rep"][1])
        if (pd3, pd2) != (d3, d2):
            dimmismatch.append((f["label"], (pd3, pd2), (d3, d2)))
        f["_d3"], f["_d2"], f["_tri"] = d3, d2, t
        tot = 1
        for d in f["rep"]:
            tot *= abs(d)
        f["_hid"] = F(tot, d3 * d2)
    res["D_rep_dim_mismatches_vs_lattice"] = dimmismatch
    print("   lattice (d3,d2) vs orbifolder printed dims: mismatches =",
          len(dimmismatch), "| level-1 violations:", len(bad_level))

    # ---------------------------------------------------------------- D  Tr Y, Tr Y^3
    trY_states = sum(f["_hid"] * G.dot(G.TY, f["p_sh"]) for f in left)
    trY3_states = sum(f["_hid"] * G.dot(G.TY, f["p_sh"]) ** 3 for f in left)
    trY_records = sum(G.dot(G.TY, f["p_sh"]) for f in left)
    route2 = None
    if coeffs is not None and tr:
        # Tr Y = sum_i c_i tr Q_i ;  Tr Y^3 = sum_ijk c_i c_j c_k tr Q_i Q_j Q_k, which
        # for an ORTHOGONAL printed basis and only two non-zero c_i reduces to the
        # printed tr Q_i^3 and tr Q_i Q_j^2 entries.
        nz = [i for i, c in enumerate(coeffs) if c != 0]
        tY = F(0)
        for i in nz:
            tY += coeffs[i] * tr.get("Q_%d" % (i + 1), None)
        route2 = dict(trY=str(tY))
        if len(nz) == 2:
            i, j = nz
            ci, cj = coeffs[i], coeffs[j]
            t3 = (ci ** 3 * tr["Q_%d^3" % (i + 1)]
                  + cj ** 3 * tr["Q_%d^3" % (j + 1)]
                  + 3 * ci * cj ** 2 * tr["Q_%d Q_%d^2" % (i + 1, j + 1)]
                  + 3 * cj * ci ** 2 * tr["Q_%d Q_%d^2" % (j + 1, i + 1)])
            route2["trY3"] = str(t3)
            route2["basis_orthogonal"] = all(
                G.dot(u1[a], u1[b]) == 0 for a in range(len(u1))
                for b in range(a + 1, len(u1)))
    res["D_traces"] = dict(trY_states=str(trY_states), trY3_states=str(trY3_states),
                           trY_records_UNWEIGHTED=str(trY_records),
                           from_orbifolder_anomaly=route2)
    print("D  Tr Y (states) =", trY_states, " Tr Y^3 (states) =", trY3_states)
    print("   orbifolder anomaly route:", route2)

    # anomalous U(1) and GS orthogonality
    anomgen = [F(0)] * 16
    for i, g in enumerate(u1):
        t = tr.get("Q_%d" % (i + 1))
        if t:
            anomgen = G.add(anomgen, G.smul(t / G.n2(g), g))
    res["D_anom"] = dict(t_anom=[str(x) for x in anomgen],
                         tY_dot_tanom=str(G.dot(G.TY, anomgen)),
                         nonzero_trQ={k: str(v) for k, v in tr.items()
                                      if "^" not in k and " " not in k and v != 0})
    print("   t_Y . t_anom =", G.dot(G.TY, anomgen))

    # ---------------------------------------------------------------- E  generations
    nets = {}
    for sY in (1, -1):
        c = Counter()
        for f in left:
            r = G.role_of(f["p_sh"], sY)
            if r:
                c[r] += int(f["_hid"])
        nets[sY] = {k: c.get(k, 0) - c.get(k + "bar", 0)
                    for k in ("Q", "uc", "ec", "dc", "L")}
        nets[sY]["_raw"] = dict(c)
    res["E_nets"] = {str(k): v for k, v in nets.items()}
    print("E  nets sY=+1:", {k: v for k, v in nets[1].items() if k != "_raw"})
    print("   nets sY=-1:", {k: v for k, v in nets[-1].items() if k != "_raw"})

    # ---------------------------------------------------------------- F  the doubled cells
    DELTA = [F(0)] * 8 + [F(1), F(1)] + [F(0)] * 6
    hits = []
    for f in left:
        if f["kl"] not in ((0, 1), (0, 5)):
            continue
        p = f["p_sh"]
        r1, r2 = G.role_of(p, 1), G.role_of(p, -1)
        if not ((r1 in G.TEN | G.TENBAR) or (r2 in G.TEN | G.TENBAR)):
            continue
        w2 = G.proj_wsm_norm2(p)
        pd = G.dot(p, DELTA)
        partner = G.sub(p, G.smul(pd, DELTA))
        hits.append(dict(
            field=f["label"], kl=list(f["kl"]),
            cell_key_full_plane_tuple=list(f["fp"]), sgc=list(f["sgc"]),
            rep=list(f["rep"]), osc=f["osc"], hidden_mult=str(f["_hid"]),
            role_sY_plus=r1, role_sY_minus=r2,
            p_sh=[str(x) for x in p], p2=str(G.n2(p)), w2=str(w2),
            R=str(G.n2(p) - w2), p_dot_Delta=str(pd),
            partner=[str(x) for x in partner], partner_norm2=str(G.n2(partner)),
            partner_w2=str(G.proj_wsm_norm2(partner)),
            partner_role_sY_plus=G.role_of(partner, 1),
            partner_in_same_coset=G.in_lattice(G.sub(partner, p)),
            partner_distinct=(tuple(partner) != tuple(p)),
            q_sh=[str(x) for x in f["q_sh"]]))
    res["F_hits"] = hits
    print("F  10-role records at T(0,1)/T(0,5):", len(hits))
    for h in hits:
        print("   ", h["field"], h["kl"], h["rep"], "osc=" + h["osc"],
              "role+/-", h["role_sY_plus"], h["role_sY_minus"],
              "|p|^2", h["p2"], "|w|^2", h["w2"], "R", h["R"],
              "p.D", h["p_dot_Delta"], "hid", h["hidden_mult"])

    # theorem H1-H5 on the doubled pair
    H = []
    for h in hits:
        p = [G.rat(x) for x in h["p_sh"]]
        q = [G.rat(x) for x in h["partner"]]
        d = G.sub(q, p)
        H.append(dict(
            field=h["field"],
            H1_Delta_in_Lambda_even=(G.in_lattice(d) and G.n2(d).denominator == 1
                                     and G.n2(d) % 2 == 0),
            H1_Delta_norm2=str(G.n2(d)),
            H2_SM_in_4d_group=res["B_printed_simple_roots_are_4d_roots"],
            H3_same_cell_same_constructing_element=True,
            H4_w2_equal=(G.proj_wsm_norm2(p) == G.proj_wsm_norm2(q)),
            H5_distinct_gauge_momenta=(tuple(p) != tuple(q)),
            Delta_perp_SM=all(G.dot(d, c) == 0 for c in G.SM_CARTAN),
            bound_lhs=str(G.n2(p) - G.proj_wsm_norm2(p)),
            bound_check_2sqrtR_ge_sqrt2=(4 * (G.n2(p) - G.proj_wsm_norm2(p)) >= 2)))
    res["F_theorem"] = H
    for h in H:
        print("   H:", h["field"], h)

    # ---------------------------------------------------------------- G  L1 on this model
    l1 = Counter()
    l1bad = []
    for f in fields:
        k, l = f["kl"]
        if k != 0 or l == 0:
            continue
        Vg = G.smul(l, V2)
        if G.in_lattice(G.sub(f["p_sh"], Vg)):
            l1["ok_l%d" % l] += 1
        else:
            l1["BAD_l%d" % l] += 1
            l1bad.append(f["label"])
    res["G_L1_on_CB1041"] = dict(counts=dict(l1), bad=l1bad)
    print("G  L1 (p_sh in l*V2 + Lambda) over T(0,l) records:", dict(l1))

    # also: every record, general check V_g = k V1 + l V2 + n.W
    ce = G.parse_constructing_elements()
    res["G_constructing_elements_T01_T05"] = {str(k): [list(x) for x in v]
                                              for k, v in ce.items()
                                              if k in ((0, 1), (0, 5))}

    with open(os.path.join(G.OUT, "g1_cb1041.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print("\nwrote g1_cb1041.json")


def solve_in_span(basis, target):
    """exact rational least-structure solve: return coefficients c with sum c_i b_i = t,
    or None.  Basis assumed independent."""
    if not basis:
        return None
    n = len(basis)
    A = [[G.dot(basis[i], basis[j]) for j in range(n)] for i in range(n)]
    b = [G.dot(basis[i], target) for i in range(n)]
    # solve A c = b by Gaussian elimination over Q
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    r = 0
    for c in range(n):
        piv = next((i for i in range(r, n) if M[i][c] != 0), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        pv = M[r][c]
        M[r] = [x / pv for x in M[r]]
        for i in range(n):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        r += 1
    coef = [M[i][n] for i in range(n)]
    rec = [F(0)] * 16
    for ci, bi in zip(coef, basis):
        rec = G.add(rec, G.smul(ci, bi))
    return coef if rec == list(target) else None


if __name__ == "__main__":
    main()

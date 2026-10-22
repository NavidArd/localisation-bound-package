"""GATE su5_ten_bound -- D7/D9: is the loophole STRUCTURALLY closed, or genuinely open?

The phase reports the loophole "open in principle and unrealised in the catalogue" on the
strength of an exact shell count at 40 cells.  That is an absence of counterexample, not a
proof of possibility.  The letter adjudication (DL-2 vs DL-1-with-caveat) turns on which:

  * if the three loophole conditions are unsatisfiable in principle, the honest letter is
    DL-1 with a non-universality caveat;
  * if a lattice-consistent configuration EXISTS, DL-2 stands and CB-1 gets a real target.

This script settles it by CONSTRUCTION, from scratch, with no model data.

Construction.  Put the Standard Model inside an SU(5) whose roots are e_i - e_j, i<j<=5,
inside the first E8 (these are D8 roots, hence E8 roots).  Then
    W = SM Cartan span = {x in span(e_1..e_5) : sum x_i = 0},   dim 4,
    the 10-weights are  w = e_i + e_j - (2/5)u,  u = e_1+...+e_5,  |w|^2 = 6/5,
    the 5bar-weights are -f_i = -(e_i - u/5),                     |w|^2 = 4/5,
so the hypercharge is automatically GUT-normalised (k_Y = 5/3) -- the same normalisation
the catalogue measures.

Required, for a SECOND state of the same 10-role at one T(0,1)/T(0,5) fixed point of
Z3xZ6_1_1:
    |p|^2 = 2B = 31/18                       (budget, no oscillator)
    P_W p = w, a 10-weight                   (the role)
    a root Delta with Delta . W = 0 and p . Delta = 1
                                             (then p and p - Delta are two massless states
                                              of the SAME role at the SAME cell)
    6 p in Lambda                            (the constructing element of T(0,1) has order
                                              6, so 6 V_g in Lambda and V_g == p mod Lambda)

The script solves this exactly, verifies every condition independently, enumerates the
complete massless shell of the resulting local shift, and checks the leading modular-
invariance condition N(V^2 - v^2) = 0 mod 2 for the Z6 generator.

Output: data/gate_su5_bound/d7_construction.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
CAMP = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(CAMP, "code", "su5_ten_bound"))
sys.path.insert(0, HERE)
import s5lib as S                                                       # noqa: E402
from d7_shell import e8_coset_shell, shell                              # noqa: E402

OUT = os.path.join(CAMP, "data", "gate_su5_bound")


def vec(*xs):
    return tuple(F(x) if not isinstance(x, F) else x for x in xs)


def main():
    rep = {}
    # ---------------------------------------------------------------- the SM inside SU(5)
    e = [tuple(F(1) if k == i else F(0) for k in range(16)) for i in range(16)]
    u = tuple(sum(e[i][k] for i in range(5)) for k in range(16))          # e1+..+e5
    su5_roots = [tuple(a - b for a, b in zip(e[i], e[j]))
                 for i in range(5) for j in range(5) if i != j]
    assert all(S.dot(r, r) == 2 for r in su5_roots) and len(su5_roots) == 20
    ROOTS = set(S.e8xe8_roots())
    assert all(r in ROOTS for r in su5_roots), "SU(5) roots must be E8 roots"
    rep["su5_roots_are_E8_roots"] = True

    # colour = indices 0,1,2 ; weak = 3,4
    a3 = [tuple(a - b for a, b in zip(e[0], e[1])),
          tuple(a - b for a, b in zip(e[1], e[2]))]
    a2 = [tuple(a - b for a, b in zip(e[3], e[4]))]
    # t_Y: the SU(5) Cartan direction orthogonal to SU(3) and SU(2) roots
    tY = tuple(F(-1, 3) if k < 3 else (F(1, 2) if k < 5 else F(0)) for k in range(16))
    assert all(S.dot(tY, r) == 0 for r in a3 + a2)
    rep["k_Y"] = str(2 * S.dot(tY, tY))                                  # must be 5/3
    smbasis = a3 + a2 + [tY]

    # ---------------------------------------------------------------- the 10-weight
    w = tuple(a + b - F(2, 5) * c for a, b, c in zip(e[0], e[1], u))     # e1+e2-(2/5)u
    rep["w2_of_10weight"] = str(S.dot(w, w))                             # must be 6/5
    f0 = tuple(a - F(1, 5) * c for a, c in zip(e[0], u))
    rep["w2_of_5weight"] = str(S.dot(f0, f0))                            # must be 4/5

    # ---------------------------------------------------------------- solve for p
    # p = (e1 + e2 + tau*u) (+) (1/2, 1/2, 0^6),  Delta = e_9 + e_10 in the SECOND E8.
    # |p|^2 = 31/18 fixes tau; 6p in Lambda fixes tau to a rational with denominator 3.
    sols = []
    for m in range(-40, 41):                       # 6*tau = m, integrality of 6p
        tau = F(m, 6)
        z = tuple(a + b + tau * c for a, b, c in zip(e[0], e[1], u))[:8]
        y = (F(1, 2), F(1, 2)) + (F(0),) * 6
        p = z + y
        if S.dot(p, p) != F(31, 18):
            continue
        sols.append((m, p))
    rep["n_solutions_scanned"] = len(sols)
    assert sols, "no p with |p|^2 = 31/18 in this family"
    m, p = sols[0]
    rep["tau"] = str(F(m, 6))
    rep["p"] = [str(x) for x in p]
    rep["p2"] = str(S.dot(p, p))

    Delta = tuple(F(1) if k in (8, 9) else F(0) for k in range(16))
    rep["Delta"] = [str(x) for x in Delta]
    rep["Delta_is_root"] = Delta in ROOTS
    rep["Delta_perp_W"] = all(S.dot(Delta, b) == 0 for b in smbasis)
    rep["p_dot_Delta"] = str(S.dot(p, Delta))
    q = tuple(a - b for a, b in zip(p, Delta))
    rep["q2"] = str(S.dot(q, q))

    # ---- both states carry the SAME SM role: identical SU(3)/SU(2) labels and Y
    def smq(x):
        return ([str(S.dot(x, r)) for r in a3], [str(S.dot(x, r)) for r in a2],
                str(S.dot(tY, x)), str(S.proj_norm2(x, smbasis)))
    rep["sm_p"] = smq(p)
    rep["sm_q"] = smq(q)
    rep["same_role"] = smq(p) == smq(q)
    rep["PW_p_equals_10weight"] = str(S.proj_norm2(p, smbasis)) == "6/5"

    # ---- 6p in Lambda (order-6 local shift), and the order is exactly 6
    rep["6p_in_Lambda"] = S.in_lattice(tuple(6 * x for x in p))
    rep["order_of_Vg"] = next(n for n in range(1, 13)
                              if S.in_lattice(tuple(n * x for x in p)))

    # ---- the SM group is UNBROKEN at this local shift (H2 holds locally)
    rep["SM_roots_integral_on_Vg"] = all(S.dot(r, p).denominator == 1
                                         for r in su5_roots)
    g = S.analyse_group(list(p))
    rep["local_group"] = g["name"]
    rep["local_n_roots"] = g["n_roots"]

    # ---- leading modular-invariance condition for the Z6 generator of Z3xZ6:
    #      v = (0, 1/6, -1/6);  N (V^2 - v^2) = 0 mod 2 with N = 6.
    v2 = F(0) ** 2 + F(1, 6) ** 2 + F(1, 6) ** 2
    mi = 6 * (S.dot(p, p) - v2)
    rep["v2"] = str(v2)
    rep["N(V2-v2)"] = str(mi)
    rep["modular_leading_ok"] = (mi.denominator == 1 and int(mi) % 2 == 0)

    # ---- delta_c and the budget, recomputed
    eta = [F(0), F(1, 6), F(5, 6)]
    dc = sum(x * (1 - x) for x in eta) / 2
    rep["delta_c"] = str(dc)
    rep["B"] = str(1 - dc)
    rep["2B"] = str(2 * (1 - dc))
    rep["budget_matches_p2"] = (2 * (1 - dc) == S.dot(p, p))
    rep["R"] = str(2 * (1 - dc) - F(6, 5))
    rep["s2_required"] = str(2 * (1 - dc) - F(6, 5) - F(1, 2))

    # ---- COMPLETE massless shell of this local shift, and the role multiplicities
    sts = shell(list(p), F(31, 18))
    rep["shell_size"] = len(sts)
    from collections import Counter
    cnt = Counter()
    for x in sts:
        l3 = tuple(S.dot(x, r) for r in a3)
        l2 = tuple(S.dot(x, r) for r in a2)
        if any(t.denominator != 1 for t in l3 + l2):
            cnt["NONINT"] += 1
            continue
        key = (tuple(int(t) for t in l3), tuple(int(t) for t in l2),
               str(S.dot(tY, x)), str(S.proj_norm2(x, smbasis)))
        cnt[str(key)] += 1
    rep["shell_printed_in_shell"] = (tuple(p) in set(sts) and tuple(q) in set(sts))

    def shellkey(x):
        return str((tuple(int(S.dot(x, r)) for r in a3),
                    tuple(int(S.dot(x, r)) for r in a2),
                    str(S.dot(tY, x)), str(S.proj_norm2(x, smbasis))))
    key_p = shellkey(p)
    rep["key_p"] = key_p
    rep["key_q"] = shellkey(q)
    rep["shell_multiplicity_of_the_10_role"] = cnt.get(key_p, 0)
    rep["shell_keys_with_multiplicity_ge2"] = {k: c for k, c in cnt.items() if c >= 2}
    rep["shell_max_multiplicity"] = max(cnt.values())

    ok = (rep["Delta_is_root"] and rep["Delta_perp_W"] and rep["same_role"]
          and rep["6p_in_Lambda"] and rep["p2"] == "31/18" and rep["q2"] == "31/18"
          and rep["k_Y"] == "5/3" and rep["w2_of_10weight"] == "6/5"
          and rep["SM_roots_integral_on_Vg"]
          and rep["shell_multiplicity_of_the_10_role"] >= 2)
    rep["VERDICT_loophole_lattice_realisable"] = bool(ok)

    with open(os.path.join(OUT, "d7_construction.json"), "w") as fh:
        json.dump(rep, fh, indent=1, default=str)
    print(json.dumps(rep, indent=1, default=str))


if __name__ == "__main__":
    main()

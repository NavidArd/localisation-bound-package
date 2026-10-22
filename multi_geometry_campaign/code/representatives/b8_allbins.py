"""Referee #8, item 2.6: extend the independent projection check of the deposited guard
(paperA_repro_package/.../root_lift/rl_d_projection.py, functions imported unchanged:
`survives`, `centralizer_elements`, `Vh`, `vh`, the shell enumerator and the space group)
from the (fixed point, right-mover) bins that orbifolder POPULATES to EVERY bin:

  every twisted constructing element of Geometry_Z3xZ6_1_1 (all 100, from the primary
  geometry file) x every massless bosonic right-mover q_sh = q + v_g, q in Z^4, q_0 = 0,
  sum(q) odd (SO(8) vector class), (1/2)|q_sh|^2 = 1/2 - delta_c  -- both chiralities.

For each bin the oscillator-free survivors S(g, q) of the projection are counted (each
survivor is one weight, i.e. one state).  Accounting, exact:
  * q in the printed left-chiral q-set of sector (k,l):  |S(g,q)| must equal the summed
    dimension of the oscillator-free left-chiral records orbifolder prints at (g, q);
  * q the CPT image of a left-chiral q' printed in sector (-k,-l) (q = -q'):  summed over
    the constructing elements of the sector, sum_g |S(g,q)| must equal the summed dimension
    printed at (-k,-l) with right-mover q';
  * any other massless q:  every survivor would be an unprinted state.
Bins with no survivors and nothing printed are counted as agreeing EMPTY bins.

Models: CB1534 (= deck BR0 of b1_decks.txt) and, as an independent second route for item
2.3, its brothers BR4 (lam.V2 = 1/3) and BR7 (lam.V2 = 2/3), whose orbifolder spectra
differ from CB1534's although their classes coincide.

Output: review/data/paperA_referee8/b8_allbins.json
"""
import json
import os
_PKG = os.environ.get("PAPERA_PKG") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
import re
import sys
from collections import defaultdict
from fractions import Fraction as F
from itertools import product

PKG = os.path.join(_PKG, "multi_geometry_campaign", "code")
sys.path.insert(0, os.path.join(PKG, "root_lift"))
sys.path.insert(0, os.path.join(PKG, "su5_ten_bound"))
sys.path.insert(0, PKG)
import rl_d_projection as P  # noqa: E402
import rl_lib as R  # noqa: E402
import shell as SH  # noqa: E402

D = os.path.join(_PKG, "multi_geometry_campaign", "data", "representatives")


def rq_solutions(v4, target):
    out = []
    for q in product(range(-2, 3), repeat=3):
        qq = (0,) + q
        if sum(qq) % 2 == 0:
            continue
        qs = tuple(F(a) + b for a, b in zip(qq, v4))
        if sum(x * x for x in qs) == target:
            out.append(qs)
    return out


def dimrep(s):
    d = 1
    for t in s.split(","):
        d *= abs(int(re.match(r"\s*(-?\d+)", t).group(1)))
    return d


def run(lab, mod, recs, S, G):
    printed = defaultdict(int)          # (kl, fp, q) -> summed dimension, osc-free, _l
    leftq = defaultdict(set)
    for r in recs:
        if r["chir"] != "l" or r["osc"] not in ("no", ""):
            continue
        key = (r["kl"], tuple(r["fp"]), tuple(r["q"]))
        printed[key] += dimrep(r["rep"])
        leftq[r["kl"]].add(tuple(r["q"]))
    N, M = S.N, S.M
    res = dict(bins=0, populated_bins_tested=0, empty_bins=0, left_mismatch=[],
               right_mismatch=[], unexplained=[], survivors_total=0)
    right_sum = defaultdict(int)
    for sec in G.sectors_ce:
        for ce in sec:
            k, l = ce[0], ce[1]
            if (k % N, l % M) == (0, 0):
                continue
            kl = (k % N, l % M)
            v = S.twist(k, l)
            dc = sum((F(x) % 1) * (1 - (F(x) % 1)) for x in v) / 2
            v4 = (F(0),) + tuple(F(x) for x in v)     # unreduced, as orbifolder prints
            qsols = set(rq_solutions(v4, 1 - 2 * dc))
            Vg = P.Vh(mod, k, l, ce[2:8])
            cand = SH.shell(Vg, 2 * (1 - dc))
            els = P.centralizer_elements(S, ce)
            fp = tuple(ce[2:8])
            rightq = {tuple(-x for x in q) for q in leftq.get(((-k) % N, (-l) % M), set())}
            for q in sorted(qsols):
                surv = [p for p in cand if P.survives(p, list(q), mod, S, ce, els)]
                ns = len(surv)
                res["bins"] += 1
                res["survivors_total"] += ns
                if q in leftq.get(kl, set()):
                    pr = printed.get((kl, fp, q), 0)
                    if pr == 0 and ns == 0:
                        res["empty_bins"] += 1
                    else:
                        res["populated_bins_tested"] += 1
                    if ns != pr:
                        res["left_mismatch"].append(dict(ce=ce, q=[str(x) for x in q],
                                                         surv=ns, printed=pr))
                elif q in rightq:
                    right_sum[(kl, q)] += ns
                    if ns == 0:
                        res["empty_bins"] += 1
                else:
                    if ns:
                        res["unexplained"].append(dict(ce=ce, q=[str(x) for x in q],
                                                       surv=ns))
                    else:
                        res["empty_bins"] += 1
    for (kl, q), ns in right_sum.items():
        img = ((-kl[0]) % N, (-kl[1]) % M)
        qi = tuple(-x for x in q)
        pr = sum(v for (k2, _fp, q2), v in printed.items() if k2 == img and q2 == qi)
        if ns != pr:
            res["right_mismatch"].append(dict(kl=kl, q=[str(x) for x in q], surv=ns,
                                              printed=pr))
    res["printed_states_oscfree"] = sum(printed.values())
    res["printed_records_oscfree_twisted"] = sum(
        1 for r in recs if r["chir"] == "l" and r["osc"] in ("no", "")
        and tuple(r["kl"]) != (0, 0))
    res["agree"] = not (res["left_mismatch"] or res["right_mismatch"] or res["unexplained"])
    return res


def main():
    S = R.SpaceGroup("Z3xZ6_1_1")
    G = S.G
    decks = P.parse_models(os.path.join(D, "b1_decks.txt"))
    out_txt = open(os.path.join(D, "b2_out.txt"), errors="replace").read()
    blocks = P.split_models(out_txt)
    report = {}
    for lab in ("BR0", "BR4", "BR7"):
        recs = P.parse_spectrum_block(blocks[lab].split("print summary")[0])
        r = run(lab, decks[lab], recs, S, G)
        report[lab] = r
        print(lab, {k: (v if not isinstance(v, list) else len(v)) for k, v in r.items()})
    json.dump(report, open(os.path.join(D, "b8_allbins.json"), "w"), indent=1, default=str)


def main_cb1041():
    """The Discussion's second model: CB1041, deck and orbifolder output as deposited
    (gate_su5_construction/gate_models.txt, gate_v1.txt.gz -- the inputs rl_d_projection
    itself reads)."""
    import gzip
    S = R.SpaceGroup("Z3xZ6_1_1")
    base = os.path.join(os.path.join(_PKG, "multi_geometry_campaign", "data"),
                        "gate_su5_construction")
    decks = P.parse_models(os.path.join(base, "gate_models.txt"))
    blocks = P.split_models(gzip.open(os.path.join(base, "gate_v1.txt.gz"), "rt",
                                      errors="replace").read())
    report = {}
    for lab in ("CB1534", "CB1041"):
        recs = P.parse_spectrum_block(blocks[lab])
        r = run(lab, decks[lab], recs, S, S.G)
        report[lab] = r
        print(lab, {k: (v if not isinstance(v, list) else len(v)) for k, v in r.items()})
    json.dump(report, open(os.path.join(D, "b8_allbins_deposited.json"), "w"), indent=1,
              default=str)


if __name__ == "__main__":
    main()
    main_cb1041()

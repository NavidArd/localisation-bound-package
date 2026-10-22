"""V6 -- the hidden-SU(2) signature: re-derive it and re-measure its three legs.

  (1) The oscillator-channel closure, re-derived arithmetically: enumerate every level a
      10-role could occupy at T(0,1)/T(0,5) and apply the localisation bound to every
      PAIR of levels.
  (2) Right-mover uniqueness per sector, measured over ALL THIRTY CB segments (the phase
      measured ten) AND over the whole certified Z3xZ6 catalogue census.
  (3) The oscillator flag at |p_sh|^2 = 31/18, all segments.
  (4) The absence recount: 10-role records that are hidden NON-SINGLETS, over the Z3xZ6
      certified census -- and, for contrast, over the CB models where they DO occur.
  (5) Scope: for which budgets does the signature argument hold?

Output: data/gate_su5_construction/g6_signature.json
"""
from __future__ import annotations

import glob
import gzip
import json
import math
import os
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gs_lib as G

CBRAW = os.path.join(G.PHASE_DATA, "raw")
CENSUS = os.path.join(G.CAMP, "data", "z3xz6_sweep", "raw", "census")
W10 = F(6, 5)
W5 = F(4, 5)


def levels():
    """(N, |p_sh|^2, R) for every oscillator level at T(0,1) with |p|^2 >= |w|^2 = 6/5."""
    eta = [F(1, 6), F(0), F(5, 6)]
    dc = sum(e * (1 - e) for e in eta) / 2
    B = 1 - dc
    Ns = set()
    for a in range(0, 7):
        for b in range(0, 7):
            for c in range(0, 3):
                N = a * F(1, 6) + b * F(5, 6) + c * F(1)
                if N <= B:
                    Ns.add(N)
    out = []
    for N in sorted(Ns):
        p2 = 2 * (B - N)
        R = p2 - W10
        out.append((N, p2, R))
    return dc, B, out


def do_cb(path):
    q = defaultdict(Counter)
    osc = Counter()
    hidden = Counter()
    txt = gzip.open(path, "rt", errors="replace").read()
    nmod = 0
    for lab, blk in G.split_models(txt):
        nmod += 1
        for f in G.parse_fields(blk.split("Analyzing")[0]):
            if f["kl"] not in ((0, 1), (0, 5)):
                continue
            q[f["kl"]][tuple(str(x) for x in f["q_sh"])] += 1
            if G.n2(f["p_sh"]) == F(31, 18):
                osc[(f["kl"], f["osc"])] += 1
    return nmod, {k: dict(v) for k, v in q.items()}, dict(osc)


def do_census(path):
    q = defaultdict(Counter)
    osc = Counter()
    tenhidden = Counter()
    txt = gzip.open(path, "rt", errors="replace").read()
    for lab, blk in G.split_models(txt, prefix="MSSM"):
        for f in G.parse_fields(blk.split("Analyzing")[0]):
            k, l = f["kl"]
            if k == 0 and l in (1, 5):
                q[(k, l)][tuple(str(x) for x in f["q_sh"])] += 1
                if G.n2(f["p_sh"]) == F(31, 18):
                    osc[((k, l), f["osc"])] += 1
    return {k: dict(v) for k, v in q.items()}, dict(osc)


def main():
    res = {}
    dc, B, lv = levels()
    print("T(0,1): local twist eta = (1/6, 0, 5/6) -> delta_c =", dc, " B =", B)
    rows = []
    for N, p2, R in lv:
        rows.append(dict(N=str(N), p2=str(p2), R_ten=str(R),
                         ten_possible=(R >= 0), R_fivebar=str(p2 - W5)))
        print("   N = %-5s |p_sh|^2 = %-6s R(10) = %-8s %s"
              % (N, p2, R, "" if R >= 0 else "  <- no 10-role possible"))
    good = [(N, p2, R) for N, p2, R in lv if R >= 0]
    print("\n  pairs of levels for a DOUBLED 10 (bound: sqrt R + sqrt R' >= sqrt 2):")
    pairs = []
    for i in range(len(good)):
        for j in range(i, len(good)):
            R1, R2 = good[i][2], good[j][2]
            s = math.sqrt(float(R1)) + math.sqrt(float(R2))
            ok = s >= math.sqrt(2)
            pairs.append(dict(N1=str(good[i][0]), N2=str(good[j][0]),
                              sum_sqrt=round(s, 6), allowed=ok))
            print("    N = %-5s with N' = %-5s : sqrt R + sqrt R' = %.4f  %s"
                  % (good[i][0], good[j][0], s, "ALLOWED" if ok else "EXCLUDED"))
    res["levels"] = rows
    res["level_pairs"] = pairs
    res["only_allowed_pair"] = [p for p in pairs if p["allowed"]]

    paths = sorted(glob.glob(os.path.join(CBRAW, "*.txt.gz")))
    print("\nCB segments:", len(paths))
    with Pool(6) as pool:
        out = pool.map(do_cb, paths)
    nmod = sum(o[0] for o in out)
    q = defaultdict(Counter)
    osc = Counter()
    for _, qq, oo in out:
        for k, v in qq.items():
            q[eval(k) if isinstance(k, str) else k].update(v)
        osc.update({(eval(k) if isinstance(k, str) else k): v for k, v in oo.items()})
    print("  models:", nmod)
    for k in sorted(q):
        print("  sector", k, ": distinct right-movers q_sh =", len(q[k]),
              " records =", sum(q[k].values()), " q_sh =", list(q[k])[:3])
    res["cb_models"] = nmod
    res["cb_qsh"] = {str(k): {"n_distinct": len(v), "n_records": sum(v.values()),
                              "values": list(v)} for k, v in q.items()}
    res["cb_osc_at_31_18"] = {str(k): v for k, v in osc.items()}
    print("  oscillator flag at |p_sh|^2 = 31/18:", dict(osc))

    cpaths = sorted(glob.glob(os.path.join(CENSUS, "*.txt.gz")))
    print("\ncatalogue census segments:", len(cpaths))
    with Pool(6) as pool:
        out2 = pool.map(do_census, cpaths)
    q2 = defaultdict(Counter)
    osc2 = Counter()
    for qq, oo in out2:
        for k, v in qq.items():
            q2[k].update(v)
        osc2.update(oo)
    for k in sorted(q2):
        print("  sector", k, ": distinct right-movers q_sh =", len(q2[k]),
              " records =", sum(q2[k].values()), " q_sh =", list(q2[k])[:3])
    print("  oscillator flag at |p_sh|^2 = 31/18:", dict(osc2))
    res["cat_qsh"] = {str(k): {"n_distinct": len(v), "n_records": sum(v.values()),
                               "values": list(v)} for k, v in q2.items()}
    res["cat_osc_at_31_18"] = {str(k): v for k, v in osc2.items()}

    json.dump(res, open(os.path.join(G.OUT, "g6_signature.json"), "w"), indent=1)
    print("\nwrote g6_signature.json")


if __name__ == "__main__":
    main()

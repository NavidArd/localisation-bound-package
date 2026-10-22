"""S6 -- THE DOUBLING MECHANISM, MEASURED (spec step S2).

The certified per-sector multiplicity table (data/su5_ten_bound/verify_Z3xZ6_1_1.json,
`mult_by_fam_sector`) shows that at the very sectors that host the loophole --
T(0,1) and T(0,5), the two single-fixed-point N=2 cells -- the 5bar roles DO reach
record multiplicity 2 (5 cells in T(0,1), 4 cells in T(0,5)), while the 10 roles never
exceed 1.  So record-level doubling of one role at one of THESE cells is a realised
phenomenon and NOT forbidden by the orbifold projection.

This script extracts those cells and measures, for each same-role pair (p, q):

    Delta = p - sigma(q)   for sigma in W_SM aligning the SM weights,
    |Delta|^2,  p.Delta,  Delta.V1,  Delta.V2,  Delta.W3          (V1, V2, W3 from the
                                                                   model's own catalogue
                                                                   row -- s0 verified)
    and the two records' right-moving q_sh.

WHY THIS MATTERS FOR THE CONSTRUCTION.  If Delta.V1 and Delta.W3 are integral in every
realised pair, the two states are a doublet of an unbroken SU(2) and the doubling could
only ever appear as ONE record with hidden multiplicity 2; the construction must then
impose that.  If they are NOT, the projection admits two independent records, and the
mechanism that lets it (different right-movers in the N=2 sector) must be identified.
Either way the answer is measured, not assumed.

Output: data/su5_construction/s6_mechanism.json
"""
from __future__ import annotations

import glob
import gzip
import json
import os
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F

import sc_lib as S

sys.path.insert(0, os.path.join(S.CAMP, "code", "z3xz6_sweep"))
import z36sw_lib as SW                                                    # noqa: E402

CENSUS = os.path.join(S.CAMP, "data", "z3xz6_sweep", "raw", "census")
SECT = {(0, 1), (0, 5)}


def catalogue_index():
    out = {}
    for idx, lab, sg, rows in S.read_catalogue():
        out[lab] = rows
    return out


def main():
    cat = catalogue_index()
    res = dict(cells=[], n_cells=0, roles=Counter(), verdict={})
    nseg = 0
    for pth in sorted(glob.glob(os.path.join(CENSUS, "*.txt.gz"))):
        nseg += 1
        txt = gzip.open(pth, "rt", errors="replace").read()
        for lab, blk in SW.split_models(txt):
            pre = blk.split("Analyzing")[0]
            spec = SW.parse_spectrum(pre)
            gauge = SW.parse_gauge(blk)
            if not spec or gauge.get("tY") is None:
                continue
            pos = SW.observable_positions(spec, gauge["tY"], gauge)
            if pos is None:
                continue
            a3, a2, sY, tv, bi = pos
            roles, diag = SW.role_census(spec, tv, a3, a2, sY)
            cells = defaultdict(list)
            for l, r in spec.items():
                if r["p_sh"] is None or r["chir"] != "l":
                    continue
                if tuple(r["kl"]) not in SECT:
                    continue
                ro = roles.get(l, {}).get("role")
                if ro is None:
                    continue
                cells[(tuple(r["kl"]), tuple(r["fp"]) if r["fp"] else None, ro)].append(r)
            for key, recs in cells.items():
                if len(recs) < 2:
                    continue
                (kl, fp, ro) = key
                V1, V2, Ws = cat[lab][0], cat[lab][1], cat[lab][2:]
                ent = dict(model=lab, kl=list(kl), fp=list(fp) if fp else None, role=ro,
                           n_records=len(recs), pairs=[])
                for i in range(len(recs)):
                    for j in range(i + 1, len(recs)):
                        p, q = recs[i]["p_sh"], recs[j]["p_sh"]
                        d = [a - b for a, b in zip(p, q)]
                        ent["pairs"].append(dict(
                            p2=str(S.dot(p, p)), q2=str(S.dot(q, q)),
                            d2=str(S.dot(d, d)), pdotd=str(S.dot(p, d)),
                            d_in_lattice=S.in_lattice(d),
                            dV1=str(S.dot(d, V1)), dV2=str(S.dot(d, V2)),
                            dW3=str(S.dot(d, Ws[2])),
                            dV1_int=S.dot(d, V1).denominator == 1,
                            dW3_int=S.dot(d, Ws[2]).denominator == 1,
                            qsh_i=[str(x) for x in recs[i]["qsh"]],
                            qsh_j=[str(x) for x in recs[j]["qsh"]],
                            qsh_equal=recs[i]["qsh"] == recs[j]["qsh"],
                            osc_i=str(recs[i].get("osc_raw")), osc_j=str(recs[j].get("osc_raw")),
                            NL_i=str(recs[i].get("NL")), NL_j=str(recs[j].get("NL")),
                            rep_i=str(recs[i]["rep"]), rep_j=str(recs[j]["rep"])))
                res["cells"].append(ent)
                res["roles"][ro] += 1
        print("  segment %d done, cells so far %d" % (nseg, len(res["cells"])), flush=True)

    res["n_cells"] = len(res["cells"])
    res["roles"] = dict(res["roles"])
    # verdict
    allp = [pr for c in res["cells"] for pr in c["pairs"]]
    res["verdict"] = dict(
        n_pairs=len(allp),
        d2_histogram=dict(Counter(pr["d2"] for pr in allp)),
        pdotd_histogram=dict(Counter(pr["pdotd"] for pr in allp)),
        dV1_integral=sum(1 for pr in allp if pr["dV1_int"]),
        dW3_integral=sum(1 for pr in allp if pr["dW3_int"]),
        both_integral=sum(1 for pr in allp if pr["dV1_int"] and pr["dW3_int"]),
        qsh_equal=sum(1 for pr in allp if pr["qsh_equal"]),
    )
    with open(os.path.join(S.OUT, "s6_mechanism.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res["verdict"], indent=1))


if __name__ == "__main__":
    main()

"""S0 -- conventions, parser and MODULAR-INVARIANCE rule-set validation.

The known-model invariant of this phase (CLAUDE.md 5.19(i)): the complete orbifolder
modular-invariance condition set of sc_lib.C5 must pass on ALL 4,612 catalogue models of
Z3xZ6_1_1.  If it does not, the rule set is wrong -- not the catalogue.  Only after it
passes may the same code be used to REJECT a candidate embedding in S4.

Also re-derives (does not assume) the row->object assignment C3 by measuring the order of
every row, and reproduces the seven defining conditions C1..C7 of the foundations phase.

Output: data/su5_construction/s0_invariants.json
"""
from __future__ import annotations

import json
import os
from collections import Counter
from fractions import Fraction as F

import sc_lib as S


def main():
    os.makedirs(S.OUT, exist_ok=True)
    n = 0
    ordhist = Counter()
    mi_fail = Counter()
    mi_fail_ex = []
    c_fail = Counter()
    wl_nontrivial = 0
    wl_vectors = set()
    sg_ok = 0
    for idx, lab, sg, rows in S.read_catalogue():
        n += 1
        if len(rows) != 8 or any(len(r) != 16 for r in rows):
            c_fail["C1_block"] += 1
            continue
        if sg and sg.endswith("Geometry_Z3xZ6_1_1.txt"):
            sg_ok += 1
        else:
            c_fail["C2_spacegroup"] += 1
        V1, V2, Ws = rows[0], rows[1], rows[2:]
        o1, o2 = S.order_of(V1), S.order_of(V2)
        ordhist[(o1, o2)] += 1
        if not (o1 and 3 % o1 == 0):
            c_fail["C3_ordV1"] += 1
        if not (o2 and 6 % o2 == 0):
            c_fail["C4_ordV2"] += 1
        for i, W in enumerate(Ws):
            if S.WL_ORDERS[i] == 1 and any(x != 0 for x in W):
                c_fail["C7_zeroWL_slot%d" % (i + 1)] += 1
            if any(x != 0 for x in W):
                ow = S.order_of(W)
                if not (ow and S.WL_ORDERS[i] % ow == 0):
                    c_fail["C5_ordW%d" % (i + 1)] += 1
        if Ws[2] != Ws[3]:
            c_fail["C6_W3neqW4"] += 1
        if any(x != 0 for x in Ws[2]):
            wl_nontrivial += 1
            wl_vectors.add(tuple(Ws[2]))
        ok, bad = S.mi_report(V1, V2, Ws)
        if not ok:
            for k in bad:
                mi_fail[k] += 1
            if len(mi_fail_ex) < 5:
                mi_fail_ex.append(dict(idx=idx, label=lab, bad={k: str(v)
                                                               for k, v in bad.items()}))

    res = dict(
        models=n,
        spacegroup_ok=sg_ok,
        defining_condition_failures=dict(c_fail),
        order_histogram={"%d,%d" % k: v for k, v in sorted(ordhist.items())},
        MI_failures=dict(mi_fail),
        MI_failure_examples=mi_fail_ex,
        MI_pass=n - sum(1 for _ in ()) if not mi_fail else None,
        wl_nontrivial=wl_nontrivial,
        wl_distinct_vectors=len(wl_vectors),
        rule_set="orbifolder-1.2.1 CheckModularInvariance, SP_factors from cspacegroup.cpp",
        SP_matrix=S.SP,
    )
    with open(os.path.join(S.OUT, "s0_invariants.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "SP_matrix"}, indent=1))


if __name__ == "__main__":
    print("selftest:", "PASS" if S.selftest(verbose=False) else "FAIL")
    main()

"""S1 -- L1 (the local shift at T(0,l)) and ROUTE A (the complete catalogue screen).

L1, DERIVED AT SOURCE, not assumed
----------------------------------
orbifolder-1.2.1/Geometry/Geometry_Z3xZ6_1_1.txt, `begin constructing elements`:
sector (0,1) lists exactly ONE constructing element, `0 1 0 0 0 0 0 0`, and sector (0,5)
exactly one, `0 5 0 0 0 0 0 0` -- both with ZERO translation.  (Geometric reason: omega
leaves plane 2 unrotated -- twist row 2 = (1/6, 0, -1/6) -- so an element (omega^l, n)
with n having a plane-2 component has NO fixed point; and this geometry's only Wilson
line, W_3 = W_4 of order 3, lives in plane 2.)  Therefore

    V_g( T(0,1) ) = V2      and     V_g( T(0,5) ) = 5 V2 = -V2  (mod Lambda)

EXACTLY, for every model, with no Wilson-line contribution and no fixed-point dependence.
This script re-parses that block mechanically and asserts it.

ROUTE A
-------
Because of L1 the loophole condition depends on V2 ALONE.  Per the theorem at B = 31/36,
4R = 94/45 < 4 forces |Delta|^2 = 2 and then p.Delta = 1; so a doubled 10-role at T(0,l)
requires the massless shell {p : |p|^2 = 31/18, p in V_g + Lambda} to contain a PAIR
p, q with p - q a ROOT.  That is a necessary condition INDEPENDENT of the model's
Standard-Model identification (it is implied by, and weaker than, "same role"), so a model
failing it cannot host the loophole.  This screens all 4,612 catalogue models exactly.

Outputs: data/su5_construction/s1_liftspace.json
"""
from __future__ import annotations

import json
import os
import re
import sys
import time
from collections import Counter
from fractions import Fraction as F

import numpy as np

import sc_lib as S

TARGET_NORM = F(31, 18)


def parse_constructing_elements(path=S.GEOM_FILE):
    """{(k,l): [ (n1..n6), ... ]} straight off the primary geometry file."""
    out, cur = {}, None
    on = False
    for raw in open(path, "r", encoding="utf-8", errors="replace"):
        s = raw.strip()
        if s == "begin constructing elements":
            on = True
            continue
        if s == "end constructing elements":
            break
        if not on or s == "new sector":
            continue
        f = [int(x) for x in s.split()]
        if len(f) != 8:
            continue
        out.setdefault((f[0], f[1]), []).append(tuple(f[2:]))
    return out


def root_pairs(sh):
    """Indices (i<j) of shell members with |p_i - p_j|^2 = 2.

    On this shell every member has |p|^2 = 31/18, so
        |p - q|^2 = 2*(31/18) - 2 p.q = 2   <=>   p.q = 13/18,
    and then p.(p-q) = 31/18 - 13/18 = 1 automatically.  In the integer scaling X = 12p
    this is X.Y = 144 * 13/18 = 104.  Computed as a Gram matrix (numpy), which also
    removes the 480-root loop entirely.
    """
    if not sh:
        return []
    X = np.array([[int(c * 12) for c in p] for p in sh], dtype=np.int64)
    G = X @ X.T
    iu = np.triu_indices(len(sh), 1)
    m = G[iu] == 104
    return list(zip(iu[0][m].tolist(), iu[1][m].tolist()))


def main():
    os.makedirs(S.OUT, exist_ok=True)
    ce = parse_constructing_elements()
    l1 = {}
    for l in (1, 5):
        els = ce.get((0, l), [])
        l1["T(0,%d)" % l] = dict(n_constructing_elements=len(els),
                                 translations=[list(e) for e in els],
                                 all_zero=all(all(x == 0 for x in e) for e in els))
    print("L1:", json.dumps(l1))
    assert all(v["all_zero"] and v["n_constructing_elements"] == 1 for v in l1.values()), \
        "L1 FAILED -- T(0,l) constructing elements are not unique / not translation-free"

    t0 = time.time()
    n = 0
    hist = Counter()
    hits = []
    shell_sizes = Counter()
    for idx, lab, sg, rows in S.read_catalogue():
        n += 1
        V1, V2, Ws = rows[0], rows[1], rows[2:]
        for l in (1, 5):
            Vg = [l * x for x in V2]
            sh = S.shell_int(Vg, 31, 18)
            shell_sizes[(l, len(sh))] += 1
            if not sh:
                hist[(l, "empty_shell")] += 1
                continue
            pr = root_pairs(sh)
            hist[(l, "pairs=%d" % len(pr))] += 1
            if pr:
                hits.append(dict(idx=idx, label=lab, l=l, shell=len(sh), n_pairs=len(pr)))
        if n % 500 == 0:
            print("  %d models  %.1fs" % (n, time.time() - t0), flush=True)

    res = dict(
        L1=l1,
        models=n,
        target_norm2="31/18",
        histogram={"T(0,%d) %s" % k: v for k, v in sorted(hist.items(), key=str)},
        shell_size_histogram={"T(0,%d) size %d" % k: v
                              for k, v in sorted(shell_sizes.items())},
        n_models_with_root_pair=len(hits),
        hits=hits[:200],
        seconds=round(time.time() - t0, 1),
    )
    with open(os.path.join(S.OUT, "s1_liftspace.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "hits"}, indent=1)[:4000])
    print("models with a root-pair in the T(0,l) shell:", len(hits))


if __name__ == "__main__":
    assert S.selftest(verbose=False)
    main()

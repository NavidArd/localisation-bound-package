"""V4 -- does the declared equivalence really close?

Two tests, both settled by orbifolder rather than by the gate's own algebra:

  MERGE   take 20 of the phase's models, apply a generator of the declared equivalence
          group (a reflection sigma_alpha with alpha perp W_SM, alpha perp Delta,
          alpha . V2 = 0) simultaneously to V1, V2 and W3, and check that orbifolder
          returns the SAME spectrum.  The fingerprint is basis-independent: the multiset
          of (sector, fixed point, |p_sh|^2, sorted |rep dims|, oscillator flag) over all
          printed records, plus the gauge-group string.
  SEPARATE take 10 pairs of models the reduction declares INEQUIVALENT and check that
          their fingerprints differ -- otherwise the reduction under-merges and 2,912
          over-counts.

Output: data/gate_su5_construction/g9_equiv.json  (+ one orbifolder run)
"""
from __future__ import annotations

import gzip
import json
import os
import random
import sys
from collections import Counter
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gs_lib as G

Z = [F(0)] * 16


def fingerprint(blk):
    fp = Counter()
    for f in G.parse_fields(blk.split("Analyzing")[0]):
        fp[(f["kl"], f["fp"], str(G.n2(f["p_sh"])),
            tuple(sorted(abs(d) for d in f["rep"])), f["chir"], f["osc"])] += 1
    return fp


def block(label, V1, V2, W3):
    Zs = ", ".join(["0"] * 16)
    rows = [", ".join(str(x) for x in V1), ", ".join(str(x) for x in V2), Zs, Zs,
            ", ".join(str(x) for x in W3), ", ".join(str(x) for x in W3), Zs, Zs]
    return ("begin model\nLabel:%s\nSpaceGroup:Geometry/Geometry_Z3xZ6_1_1.txt\n"
            "Lattice:E8xE8\nShifts and Wilsonlines:\n%s\nend model\n"
            % (label, "\n".join(rows)))


def main():
    build = json.load(open(os.path.join(G.PHASE_DATA, "s5_build.json")))
    V2 = [G.rat(t) for t in build["V2_target"]]
    DELTA = [F(0)] * 8 + [F(1), F(1)] + [F(0)] * 6
    gens = []
    for r in G.e8xe8_roots():
        v = list(r)
        if all(G.dot(v, c) == 0 for c in G.SM_CARTAN) and G.dot(v, DELTA) == 0 \
                and G.dot(v, V2) == 0:
            gens.append(v)
    print("equivalence generators (alpha perp SM, perp Delta, perp V2):", len(gens),
          "(block1 %d, block2 %d)" % (sum(1 for g in gens if any(g[:8])),
                                      sum(1 for g in gens if any(g[8:]))))

    rnd = random.Random(3)
    models = [m for m in build["models"]]
    pick = rnd.sample(models, 20)
    deck, pairs = [], []
    for n, m in enumerate(pick):
        V1 = [G.rat(t) for t in m["V1"]]
        W3 = [G.rat(t) for t in m["W3"]]
        al = gens[rnd.randrange(len(gens))]
        def refl(v):
            c = G.dot(v, al)
            return [x - c * y for x, y in zip(v, al)]
        V1b, V2b, W3b = refl(V1), refl(V2), refl(W3)
        assert G.in_lattice(G.sub(V2b, V2)), "sigma must fix the V2 class"
        ok, bad = G.mi_check(V1b, V2b, [Z, Z, W3b, W3b, Z, Z])
        assert ok, (m["label"], bad)
        deck.append(block("EQA%02d" % n, V1, V2, W3))
        deck.append(block("EQB%02d" % n, V1b, V2b, W3b))
        pairs.append(("EQA%02d" % n, "EQB%02d" % n, m["label"], "merge"))

    sep = rnd.sample(models, 20)
    for n in range(10):
        a, b = sep[2 * n], sep[2 * n + 1]
        deck.append(block("SPA%02d" % n, [G.rat(t) for t in a["V1"]], V2,
                          [G.rat(t) for t in a["W3"]]))
        deck.append(block("SPB%02d" % n, [G.rat(t) for t in b["V1"]], V2,
                          [G.rat(t) for t in b["W3"]]))
        pairs.append(("SPA%02d" % n, "SPB%02d" % n, a["label"] + "/" + b["label"],
                      "separate"))
    p = os.path.join(G.OUT, "equiv_models.txt")
    open(p, "w", newline="\n").write("".join(deck))
    print("wrote", p, "with", len(deck), "models")
    json.dump(dict(pairs=pairs), open(os.path.join(G.OUT, "g9_pairs.json"), "w"), indent=1)


def analyse():
    pairs = json.load(open(os.path.join(G.OUT, "g9_pairs.json")))["pairs"]
    txt = gzip.open(os.path.join(G.OUT, "gate_equiv.txt.gz"), "rt",
                    errors="replace").read()
    blocks = {}
    for lab, blk in G.split_models(txt, prefix="EQ"):
        blocks[lab] = blk
    for lab, blk in G.split_models(txt, prefix="SP"):
        blocks[lab] = blk
    res = []
    nm = ns = 0
    for a, b, src, kind in pairs:
        if a not in blocks or b not in blocks:
            res.append(dict(a=a, b=b, kind=kind, status="MISSING"))
            continue
        fa, fb = fingerprint(blocks[a]), fingerprint(blocks[b])
        ga, gb = G.parse_gauge_group(blocks[a]), G.parse_gauge_group(blocks[b])
        same = (fa == fb) and (ga == gb)
        res.append(dict(a=a, b=b, src=src, kind=kind, identical=same,
                        n_records=(sum(fa.values()), sum(fb.values())),
                        gauge=(ga, gb)))
        if kind == "merge" and same:
            nm += 1
        if kind == "separate" and not same:
            ns += 1
    print("MERGE pairs identical   :", nm, "/", sum(1 for r in res
                                                    if r.get("kind") == "merge"))
    print("SEPARATE pairs distinct :", ns, "/", sum(1 for r in res
                                                    if r.get("kind") == "separate"))
    for r in res:
        if r.get("kind") == "merge" and not r.get("identical"):
            print("  MERGE FAILURE:", r)
        if r.get("kind") == "separate" and r.get("identical"):
            print("  SEPARATION FAILURE (reduction under-merges):", r)
    json.dump(dict(results=res, merge_identical=nm, separate_distinct=ns),
              open(os.path.join(G.OUT, "g9_equiv.json"), "w"), indent=1)
    print("wrote g9_equiv.json")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "analyse":
        analyse()
    else:
        main()

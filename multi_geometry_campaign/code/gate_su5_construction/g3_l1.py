"""V3 -- L1's pinning, attacked directly against catalogue p_sh data.

The phase ASSERTED L1 ("V_g(T(0,l)) = l*V2 exactly, no Wilson-line contribution") from
the geometry file, and its own attack list #1 says the direct p_sh comparison was not run.
The gate runs it.

  (a) the primary geometry file: T(0,1) and T(0,5) each have exactly ONE constructing
      element, both translation-free (re-read here);
  (b) plane 2 is the omega-unrotated plane and is where the single Wilson line lives
      (allowed orders 1 1 3 3 1 1, identical WLs "3 4"), so no (omega^l, n) with a
      plane-2 translation can have a fixed point -- checked mechanically by asking the
      geometry file whether ANY listed constructing element of a T(0,l) sector carries a
      non-zero translation;
  (c) THE DIRECT TEST: over the certified 4,612-model Z3xZ6 catalogue and its own census
      output, every printed T(0,l) record must satisfy  p_sh - l*V2 in Lambda, INCLUDING
      the models whose W3 is non-zero (which is where a Wilson-line contribution would
      show up);
  (d) the counterfactual: for those same models, p_sh - (l*V2 + n*W3) is NOT in Lambda for
      n = 1, 2 -- so the test discriminates.

Output: data/gate_su5_construction/g3_l1.json
"""
from __future__ import annotations

import glob
import gzip
import json
import os
import sys
from collections import Counter
from fractions import Fraction as F
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gs_lib as G

CAT = os.path.join(G.CAMP, "data", "literature", "parr2020_full", "extracted",
                   "Z3xZ6_1_1.txt")
CENSUS = os.path.join(G.CAMP, "data", "z3xz6_sweep", "raw", "census")


def read_catalogue(path=CAT):
    out = {}
    lab, rows, reading = None, None, False
    for raw in open(path, encoding="utf-8", errors="replace"):
        s = raw.strip()
        if s == "begin model":
            lab, rows, reading = None, [], False
        elif s == "end model":
            if lab:
                out[lab] = rows
        elif s.startswith("Label:"):
            lab = s[6:].strip()
        elif s.startswith("Shifts and Wilsonlines"):
            reading = True
        elif reading and s and ":" not in s:
            f = [x.strip() for x in s.split(",")]
            if len(f) == 16:
                rows.append([G.rat(x) for x in f])
    return out


CATALOGUE = None


def _init():
    global CATALOGUE
    CATALOGUE = read_catalogue()


def do(path):
    cnt = Counter()
    bad = []
    cf = Counter()
    txt = gzip.open(path, "rt", errors="replace").read()
    for lab, blk in G.split_models(txt, prefix="MSSM"):
        rows = CATALOGUE.get(lab)
        if not rows or len(rows) < 8:
            cnt["no_catalogue_row"] += 1
            continue
        V1, V2 = rows[0], rows[1]
        Ws = rows[2:8]
        W3 = Ws[2]
        nzw = any(x != 0 for x in W3)
        for f in G.parse_fields(blk.split("Analyzing")[0]):
            k, l = f["kl"]
            if k != 0 or l == 0:
                continue
            Vg = G.smul(l, V2)
            ok = G.in_lattice(G.sub(f["p_sh"], Vg))
            tag = "W3nz" if nzw else "W3z"
            cnt["%s_l%d_%s" % (tag, l, "ok" if ok else "BAD")] += 1
            if not ok:
                bad.append((lab, f["label"], l))
            if nzw and l in (1, 5):
                for n in (1, 2):
                    alt = G.add(Vg, G.smul(n, W3))
                    cf["l%d_plus%dW3_%s" % (l, n,
                                            "inLambda" if G.in_lattice(
                                                G.sub(f["p_sh"], alt)) else "not")] += 1
    return cnt, bad, cf


def main():
    global CATALOGUE
    CATALOGUE = read_catalogue()
    print("catalogue models:", len(CATALOGUE))

    ce = G.parse_constructing_elements()
    tr = {}
    for kl, els in sorted(ce.items()):
        if kl[0] != 0:
            continue
        tr[str(kl)] = dict(n_elements=len(els),
                           translation_free=[e for e in els if all(x == 0 for x in e[2:])],
                           with_translation=[e for e in els if any(x != 0 for x in e[2:])])
    print("(a) T(0,l) constructing elements:")
    for k, v in tr.items():
        print("   ", k, "n =", v["n_elements"],
              "translation-free:", len(v["translation_free"]),
              "with translation:", len(v["with_translation"]))

    # (b) which T(0,l) sectors carry a plane-2 translation at all?
    plane2 = {k: [e for e in v["with_translation"] if e[4] or e[5]] for k, v in tr.items()}
    print("(b) T(0,l) constructing elements with a PLANE-2 translation (slots 3,4 of n):",
          {k: len(v) for k, v in plane2.items()})

    paths = sorted(glob.glob(os.path.join(CENSUS, "*.txt.gz")))
    print("census segments:", len(paths))
    with Pool(6, initializer=_init) as pool:
        res = pool.map(do, paths)
    cnt = Counter()
    bad = []
    cf = Counter()
    for c, b, f in res:
        cnt.update(c)
        bad.extend(b)
        cf.update(f)
    print("(c) direct p_sh test over the census:")
    for k in sorted(cnt):
        print("   ", k, cnt[k])
    print("    BAD total:", len(bad), bad[:10])
    print("(d) counterfactual  p_sh - (l V2 + n W3):")
    for k in sorted(cf):
        print("   ", k, cf[k])

    out = dict(constructing_elements=tr,
               plane2_translations={k: len(v) for k, v in plane2.items()},
               census_counts=dict(cnt), bad=bad[:200], n_bad=len(bad),
               counterfactual=dict(cf))
    json.dump(out, open(os.path.join(G.OUT, "g3_l1.json"), "w"), indent=1)
    print("wrote g3_l1.json")


if __name__ == "__main__":
    main()

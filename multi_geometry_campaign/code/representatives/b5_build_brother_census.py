"""Referee #8, item 2.3 (q.5): build the BROTHER census of the 4,032 Standard-Model-intact
S1-S4 embeddings of paperA_v9d section 7.2.

For every SM-intact model (2,912 class-A decks su5_construction/s5_models.txt, all SM-intact;
1,120 class-B decks gate_su5_construction/caseB_models.txt listed in g8b_smvalid.json) and
every (a, b) in {0,1,2}^2 we emit the deck with
    V1 -> V1 + lam,   W3 = W4 -> W3 + mu,     lam, mu in Lambda,
lam.V2 = a/3 (mod 1), mu.V2 = b/3 (mod 1), chosen from the vectors +-e_i +- e_j of either
E8 block such that the COMPLETE orbifolder modular-invariance set (sc_lib.mi_report,
transcribed at source level) passes for the combined data.  (a, b) = (0, 0) is the census
representative itself (re-run as a pipeline control).  The classes V1+Lambda, V2+Lambda,
W3+Lambda are unchanged by construction, so every class-level condition of section 7.2
(SM roots integral, the twelve X/Y roots non-integral, Delta.V1, Delta.W3 integral) and the
four-dimensional gauge group are unchanged.

Output: review/data/paperA_referee8/brothers/decks_XXX.txt (200 models per segment) and
build.json (label -> parent, (a,b), lam, mu; and any (model, a, b) with no MI-valid choice).
"""
import json
import os
_PKG = os.environ.get("PAPERA_PKG") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
import sys
from itertools import combinations
from fractions import Fraction as F

sys.path.insert(0, os.path.join(_PKG, "multi_geometry_campaign", "code", "su5_construction"))
import sc_lib as S  # noqa: E402

CAMP = os.path.join(_PKG, "multi_geometry_campaign", "data")
OUT = os.path.join(_PKG, "multi_geometry_campaign", "data", "representatives", "brothers")
SEG = 200
Z = [F(0)] * 16


def read_decks(path):
    res, cur, lab = [], None, None
    for raw in open(path, encoding="utf-8", errors="replace"):
        s = raw.strip()
        if s == "begin model":
            cur, lab = [], None
        elif s.startswith("Label:"):
            lab = s[6:]
        elif s == "end model":
            res.append((lab, cur))
            cur = None
        elif cur is not None and "," in s:
            cur.append([S.rat(t) for t in s.split(",")])
    return res


def cands():
    out = []
    for blk in (0, 8):
        for i, j in combinations(range(8), 2):
            for si in (1, -1):
                for sj in (1, -1):
                    v = [F(0)] * 16
                    v[blk + i], v[blk + j] = F(si), F(sj)
                    out.append(v)
    return out


def m1(x):
    return x - (x.numerator // x.denominator)


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def fmt(v):
    return ", ".join(str(x) for x in v)


def main():
    os.makedirs(OUT, exist_ok=True)
    A = read_decks(os.path.join(CAMP, "su5_construction", "s5_models.txt"))
    B = read_decks(os.path.join(CAMP, "gate_su5_construction", "caseB_models.txt"))
    validB = set(json.load(open(os.path.join(CAMP, "gate_su5_construction",
                                             "g8b_smvalid.json")))["valid"])
    models = [("A", l, r) for l, r in A] + [("B", l, r) for l, r in B if l in validB]
    assert len(models) == 4032, len(models)
    C = cands()
    meta, missing, decks = {}, [], []
    for cls, lab, rows in models:
        V1, V2, Ws = rows[0], rows[1], rows[2:8]
        assert Ws[2] == Ws[3] and all(w == Z for k, w in enumerate(Ws) if k not in (2, 3))
        W = Ws[2]
        byV = {}
        for lam in C:
            byV.setdefault(int(m1(S.dot(lam, V2)) * 3) if (m1(S.dot(lam, V2)) * 3).denominator == 1 else None, []).append(lam)
        for a in range(3):
            for b in range(3):
                found = None
                lams = [Z] if a == 0 else byV.get(a, [])
                mus = [Z] if b == 0 else byV.get(b, [])
                for lam in lams:
                    nv1 = add(V1, lam)
                    for mu in mus:
                        nw = add(W, mu)
                        ok, _ = S.mi_report(nv1, V2, [Z, Z, nw, nw, Z, Z])
                        if ok:
                            found = (lam, mu, nv1, nw)
                            break
                    if found:
                        break
                if not found:
                    missing.append([lab, a, b])
                    continue
                lam, mu, nv1, nw = found
                nl = "X%s%s" % (lab, "a%db%d" % (a, b))
                rowsout = [nv1, V2, Z, Z, nw, nw, Z, Z]
                decks.append("begin model\nLabel:%s\nSpaceGroup:Geometry/Geometry_Z3xZ6_1_1.txt"
                             "\nLattice:E8xE8\nShifts and Wilsonlines:\n%s\nend model\n"
                             % (nl, "\n".join(fmt(r) for r in rowsout)))
                meta[nl] = dict(parent=lab, cls=cls, a=a, b=b,
                                lam=[str(x) for x in lam], mu=[str(x) for x in mu])
    for k in range(0, len(decks), SEG):
        with open(os.path.join(OUT, "decks_%03d.txt" % (k // SEG)), "w", newline="\n") as fh:
            fh.write("".join(decks[k:k + SEG]))
    json.dump(dict(n_models=len(models), n_decks=len(decks), missing=missing, meta=meta),
              open(os.path.join(OUT, "build.json"), "w"))
    print("models", len(models), "decks", len(decks), "segments",
          (len(decks) + SEG - 1) // SEG, "missing", len(missing), missing[:10])


if __name__ == "__main__":
    main()

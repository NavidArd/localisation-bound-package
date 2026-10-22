"""Referee #8, item 2.3 (q.5): are the spectra of the S1-S4 embeddings functions of the
class pair (V1+Lambda, W3+Lambda), or do different modular-invariant representatives of
the SAME classes (zero discrete torsion) give different spectra ("brother models",
Ploger--Ramos-Sanchez--Ratz--Vaudrevange, mirage torsion)?

Builds, for the paper's model CB1534 (eq. (7.1) of paperA_v9d; deck
paperA_repro_package/construction/CB1534/phase_deck_s5_models_block.txt), representatives
V1 + lam and W3 = W4 -> W3 + mu with lam, mu in Lambda (E8xE8, even coordinates), keeps
only those that pass the COMPLETE orbifolder modular-invariance set (sc_lib.mi_report,
transcribed at source level from corbifoldgroup.cpp), and classifies each by the
projection-phase shift it induces, lam.V2 mod 1 and mu.V2 mod 1.

Output: review/data/paperA_referee8/b1_decks.txt (orbifolder input, stock geometry
Z3xZ6_1_1) and b1_decks.json (which lam/mu each label carries).
"""
import json
import os
_PKG = os.environ.get("PAPERA_PKG") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
import sys
from fractions import Fraction as F
from itertools import combinations

sys.path.insert(0, os.path.join(_PKG, "multi_geometry_campaign", "code", "su5_construction"))
import sc_lib as S  # noqa: E402

OUT = os.path.join(_PKG, "multi_geometry_campaign", "data", "representatives")
q = F
V1 = [q(1, 2), q(3, 2), q(1, 2), q(1, 6), q(7, 6), q(1, 2), q(1, 2), q(1, 2)] + [q(0)] * 8
V2 = [q(2, 3), q(2, 3), q(-1, 3), q(-1, 3), q(-1, 3), q(0), q(0), q(0),
      q(1, 2), q(1, 2)] + [q(0)] * 6
W3 = [q(2, 3), q(-4, 3), q(2, 3)] + [q(0)] * 5 + [q(1, 3), q(-1, 3), q(2, 3)] + [q(0)] * 5
Z = [q(0)] * 16


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def frac1(x):
    return x - (x.numerator // x.denominator)


def small_lattice_vectors():
    """integer vectors +-e_i +- e_j (roots of D8 inside each E8) in either block."""
    out = []
    for blk in (0, 8):
        for i, j in combinations(range(8), 2):
            for si in (1, -1):
                for sj in (1, -1):
                    v = [q(0)] * 16
                    v[blk + i], v[blk + j] = q(si), q(sj)
                    out.append(v)
    return out


def fmt(v):
    return ", ".join(str(x) for x in v)


def deck(label, v1, v2, w):
    rows = [v1, v2, Z, Z, w, w, Z, Z]
    s = ["begin model", "Label:%s" % label,
         "SpaceGroup:Geometry/Geometry_Z3xZ6_1_1.txt", "Lattice:E8xE8",
         "Shifts and Wilsonlines:"]
    s += [fmt(r) for r in rows]
    s.append("end model")
    return "\n".join(s) + "\n"


def main():
    ok, bad = S.mi_report(V1, V2, [Z, Z, W3, W3, Z, Z])
    assert ok, bad
    picks = {}
    for lam in small_lattice_vectors():
        nv1 = add(V1, lam)
        ok, _ = S.mi_report(nv1, V2, [Z, Z, W3, W3, Z, Z])
        if not ok:
            continue
        key = ("V1", str(frac1(S.dot(lam, V2))), str(frac1(S.dot(lam, W3))))
        picks.setdefault(key, lam)
    for mu in small_lattice_vectors():
        nw = add(W3, mu)
        ok, _ = S.mi_report(V1, V2, [Z, Z, nw, nw, Z, Z])
        if not ok:
            continue
        key = ("W", str(frac1(S.dot(mu, V2))), str(frac1(S.dot(mu, V1))))
        picks.setdefault(key, mu)
    decks, meta = [deck("BR0", V1, V2, W3)], {"BR0": dict(kind="base")}
    n = 0
    for key, v in sorted(picks.items()):
        n += 1
        lab = "BR%d" % n
        if key[0] == "V1":
            decks.append(deck(lab, add(V1, v), V2, W3))
        else:
            decks.append(deck(lab, V1, V2, add(W3, v)))
        meta[lab] = dict(kind=key[0], shift=[str(x) for x in v],
                         dot_V2_mod1=key[1], other_mod1=key[2])
    with open(os.path.join(OUT, "b1_decks.txt"), "w", newline="\n") as fh:
        fh.write("".join(decks))
    with open(os.path.join(OUT, "b1_decks.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    print(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()

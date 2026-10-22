"""GATE d8 -- independent recomputation of the class-wide twisted-sector budget table.

Adversarial re-derivation of `results/su5_ten_bound_findings.md` 3.4 / 7.1.

Twist vectors are read OFF the PRIMARY orbifolder geometry files
(F:/physics/orbifolder-1.2.1/Geometry/Geometry_<label>.txt, `begin twist` block) --
NOT taken from the phase's hard-coded GEOM dicts.  The parser is deliberately dumb:
it reads the 4-component `a/b` rows verbatim and drops the leading (0th) entry.

For each twisted sector (k,l) != (0,0):
    v      = k*v1 + l*v2                   (component-wise, exact Fractions)
    eta_a  = v_a mod 1
    delta_c= 1/2 sum_a eta_a (1 - eta_a)
    B      = 1 - delta_c                   (N_osc = 0)
SUSY check: sum_a eta_a must be an integer (N=1 SUSY, v in the SU(3) lattice).

Also re-derives the N=2 threshold claim of 7.1 and searches exhaustively for
three-plane-rotated twists with delta_c <= 3/20.

Writes data/gate_su5_bound/d8_classwide.json.
"""

from __future__ import annotations

# --- reproduction-package root shim -------------------------------------------------
# The line above originally read  F:\physics  -- the absolute path of the author's
# working tree.  This is the ONLY edit made to this file in the reproduction package;
# it points the same variable at the package root, which carries the same directory
# names.  Set PAPERA_PKG to relocate the package.  See
# environment/PATCHES_TO_PACKAGE_CODE.md.
# -------------------------------------------------------------------------------------
import os as _os
_PKG_ROOT = _os.environ.get('PAPERA_PKG') or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..'))

import json
import os
import re
from collections import Counter
from fractions import Fraction as F
from math import gcd, isqrt

GEODIR = _os.path.join(_PKG_ROOT, 'orbifolder-1.2.1', 'Geometry')
OUT = _os.path.join(_PKG_ROOT, 'multi_geometry_campaign', 'data', 'gate_su5_bound', 'd8_classwide.json')

LABELS = ["Z2xZ4_1_1", "Z4xZ4_1_1", "Z3xZ3_1_1", "Z6-I_1_1",
          "Z6-II_1_1", "Z4_1_1", "Z3xZ6_1_1"]

THR10 = F(17, 20)
THR5 = F(13, 20)


def read_geometry(label):
    """Return (M, N, [twist vectors as 3-tuples of Fraction], provenance lines)."""
    path = os.path.join(GEODIR, "Geometry_%s.txt" % label)
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    # point group orders: line index 0 == 'point group', then two ints
    assert lines[0].strip().lower() == "point group", lines[0]
    M = int(lines[1].strip())
    N = int(lines[2].strip())
    i0 = next(i for i, s in enumerate(lines) if s.strip() == "begin twist")
    i1 = next(i for i, s in enumerate(lines) if s.strip() == "end twist")
    tw, prov = [], []
    for i in range(i0 + 1, i1):
        raw = lines[i].strip()
        if not raw:
            continue
        parts = re.findall(r"-?\d+/\d+|-?\d+", raw)
        assert len(parts) == 4, (label, raw, parts)
        vec = [F(p) for p in parts]
        assert vec[0] == 0, (label, raw)          # 0th (time-like) entry must vanish
        tw.append(tuple(vec[1:]))
        prov.append(dict(file=path, line=i + 1, text=raw))
    return M, N, tw, prov


def delta_c(v):
    return sum((x % 1) * (1 - (x % 1)) for x in v) / 2


def susy_sum(v):
    return sum(x % 1 for x in v)


def sectors(M, N, tw):
    """Yield (k, l, v) for all (k,l) != (0,0), k mod M, l mod N."""
    if len(tw) == 1:
        v1, v2 = tw[0], (F(0), F(0), F(0))
    else:
        v1, v2 = tw[0], tw[1]
    for k in range(M):
        for l in range(N):
            if (k, l) == (0, 0):
                continue
            v = tuple(k * v1[a] + l * v2[a] for a in range(3))
            yield k, l, v


def frac(x):
    return "%d/%d" % (x.numerator, x.denominator) if x.denominator != 1 else str(x.numerator)


def main():
    res = {"_meta": {"source": "primary orbifolder geometry files",
                     "geometry_dir": GEODIR,
                     "formula": "delta_c = 1/2 sum_a eta_a (1-eta_a); B = 1 - delta_c (N_osc=0)",
                     "thresholds": {"ten": "17/20", "fivebar": "13/20"}}}

    for lab in LABELS:
        M, N, tw, prov = read_geometry(lab)
        secs = []
        susy_bad = []
        bcount = Counter()
        for k, l, v in sectors(M, N, tw):
            dc = delta_c(v)
            B = 1 - dc
            s = susy_sum(v)
            if s.denominator != 1:
                susy_bad.append(dict(k=k, l=l, v=[frac(x % 1) for x in v], sum=frac(s)))
            secs.append(dict(k=k, l=l,
                             eta=[frac(x % 1) for x in v],
                             delta_c=frac(dc), B=frac(B),
                             susy_sum=frac(s),
                             ten_doubling_permitted=bool(B >= THR10),
                             fivebar_doubling_permitted=bool(B >= THR5)))
            bcount[B] += 1
        Bs = [F(s["B"]) for s in secs]
        res[lab] = dict(
            point_group=[M, N],
            twist_vectors=[[frac(x) for x in t] for t in tw],
            provenance=prov,
            n_twisted=len(secs),
            n_distinct_B=len(set(Bs)),
            B_multiset={frac(b): c for b, c in sorted(bcount.items(), reverse=True)},
            B_max=frac(max(Bs)), B_min=frac(min(Bs)),
            ten_doubling_permitted=[[s["k"], s["l"]] for s in secs if s["ten_doubling_permitted"]],
            fivebar_doubling_forbidden=[[s["k"], s["l"]] for s in secs
                                        if not s["fivebar_doubling_permitted"]],
            susy_violations=susy_bad,
            sectors=secs,
        )

    # ---------------- N=2 threshold re-derivation ----------------
    # (x, 0, -x): delta_c = 1/2 [x(1-x) + 0 + (1-x)x] = x(1-x)
    # x(1-x) <= 3/20  <=>  x <= (1 - sqrt(2/5))/2   (lower root branch)
    import decimal
    decimal.getcontext().prec = 40
    D = decimal.Decimal
    xstar = (1 - (D(2) / D(5)).sqrt()) / 2
    n2 = {"delta_c_formula": "x(1-x)",
          "x_star_exact": "(1 - sqrt(2/5))/2",
          "x_star_decimal": str(+xstar)[:20],
          "check_1_over_6": {"x": "1/6", "delta_c": frac(F(1, 6) * F(5, 6)),
                             "delta_c_le_3_20": F(1, 6) * F(5, 6) <= F(3, 20)},
          "check_1_over_5": {"x": "1/5", "delta_c": frac(F(1, 5) * F(4, 5)),
                             "delta_c_le_3_20": F(1, 5) * F(4, 5) <= F(3, 20)},
          "check_1_over_4": {"x": "1/4", "delta_c": frac(F(1, 4) * F(3, 4)),
                             "delta_c_le_3_20": F(1, 4) * F(3, 4) <= F(3, 20)},
          }
    # all k/N with N<=30 satisfying x(1-x) <= 3/20, x in (0,1), gcd(k,N)=1
    ok = []
    for Nn in range(2, 31):
        for k in range(1, Nn):
            if gcd(k, Nn) != 1:
                continue
            x = F(k, Nn)
            if x * (1 - x) <= F(3, 20):
                ok.append("%d/%d" % (k, Nn))
    n2["all_primitive_k_over_N_passing_upto_N30"] = ok
    res["_N2_threshold"] = n2

    # ---------------- three-plane-rotated exhaustive search ----------------
    # all eta in (0,1) with common denominator N <= 60, sum eta_a integer,
    # minimal delta_c with ALL THREE eta_a nonzero.
    best = None
    hits = []
    for Nn in range(2, 61):
        for a in range(1, Nn):
            for b in range(1, Nn):
                c = (-a - b) % Nn
                if c == 0:
                    continue
                eta = (F(a, Nn), F(b, Nn), F(c, Nn))
                dc = sum(e * (1 - e) for e in eta) / 2
                if best is None or dc < best[0]:
                    best = (dc, Nn, (a, b, c))
                if dc <= F(3, 20):
                    hits.append(dict(N=Nn, eta=["%d/%d" % (a, Nn), "%d/%d" % (b, Nn),
                                                "%d/%d" % (c, Nn)], delta_c=frac(dc)))
    res["_three_plane_rotated"] = {
        "search": "all eta_a = a/N in (0,1), a=1..N-1, N<=60, sum eta_a in Z, all three nonzero",
        "min_delta_c": frac(best[0]), "at_N": best[1], "at": list(best[2]),
        "n_hits_delta_c_le_3_20": len(hits),
        "first_hits": hits[:12],
        "smallest_N_with_hit": min([h["N"] for h in hits]) if hits else None,
    }

    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1)

    # ---------------- console report ----------------
    print("%-12s %4s %8s %8s  %-28s %s" % ("geometry", "#tw", "B_max", "B_min",
                                           "B multiset (value x mult)", "10-doubling"))
    for lab in LABELS:
        r = res[lab]
        ms = ", ".join("%sx%d" % (k, v) for k, v in r["B_multiset"].items())
        td = r["ten_doubling_permitted"] or "NONE"
        print("%-12s %4d %8s %8s  %-28s %s" % (lab, r["n_twisted"], r["B_max"],
                                               r["B_min"], ms, td))
        if r["susy_violations"]:
            print("   !! SUSY VIOLATIONS:", r["susy_violations"])
        if r["fivebar_doubling_forbidden"]:
            print("   !! 5bar forbidden in:", r["fivebar_doubling_forbidden"])
    print()
    print("twists read:", {l: res[l]["twist_vectors"] for l in LABELS})
    print()
    print("N=2 threshold:", json.dumps(n2, indent=1))
    print()
    print("three-plane:", json.dumps(res["_three_plane_rotated"], indent=1)[:1500])


if __name__ == "__main__":
    main()

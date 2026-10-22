"""su5_ten_bound STEP 4 -- FULL-DATA VERIFICATION of the localisation bound.

THE THEOREM UNDER TEST (docs/su5_ten_bound_task.md; derivation in
results/su5_ten_bound_findings.md).  At one twisted fixed point, let two DISTINCT massless
left-chiral states p, p' carry the SAME Standard-Model role.  Then

    Delta = p - p' is a NON-ZERO vector of the EVEN lattice Lambda_{E8xE8}   =>  |Delta|^2 >= 2
    Delta is ORTHOGONAL to the SM Cartan span W  (same SM weight)
    |p|^2 = |p'|^2 = 2B                                (same masslessness budget)

so, writing R = |P_{W^perp} p|^2 = 2B - |w|^2,

    |Delta|^2 = 2B + 2B' - 2|w|^2 - 2 (p_perp . p'_perp) <= (sqrt(R) + sqrt(R'))^2
    =>  sqrt(R) + sqrt(R') >= sqrt(2),     R = 2B - |w|^2, R' = 2B' - |w|^2.

When both states are oscillator-free (B = B', R = R') this reduces to R >= 1/2.  The
GENERAL form is the one tested here: testing min(R) >= 1/2 is the equal-R special case and
is wrong for a pair straddling two oscillator levels.

This script measures R for EVERY SM-role field at EVERY twisted cell of the on-disk
catalogues and tests the prediction

    per-cell multiplicity >= 2   =>   some pair has sqrt(R_i) + sqrt(R_j) >= sqrt(2)

with |w|^2 obtained by EXPLICIT projection onto the SM Cartan span (never assumed).
It also records every cell where a 10-role has R >= 1/2 -- those are exactly the cells the
theorem does NOT protect, i.e. the DL-2 loophole targets.

Output: data/su5_ten_bound/verify_<geom>.json
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"

import glob                                                             # noqa: E402
import gzip                                                             # noqa: E402
import json                                                             # noqa: E402
import multiprocessing as mp                                            # noqa: E402
import sys                                                              # noqa: E402
from collections import Counter, defaultdict                            # noqa: E402
from fractions import Fraction as F                                     # noqa: E402

import re                                                             # noqa: E402

import numpy as np                                                      # noqa: E402

# oscillator token: `a_(mode)^(+/-n)`.  The MODE NUMBER is the quantity the mass formula
# wants, not the occurrence count (CLAUDE.md 5.21 (iv)).  The Z2xZ4 parser already returns
# it as ML/MbL; the Z3xZ6 parser returns only occurrence counts, so it is re-read here from
# the raw string with the SAME regex, and the two routes are cross-checked where both exist.
_OSC = re.compile(r"(-?)a_\(\s*(-?\d+(?:/\d+)?)\s*\)\^(?:\((-\d)\)|(\d))")


def mode_sum(rec):
    if "ML" in rec:
        return sum(F(x) for x in rec["ML"]) + sum(F(x) for x in rec["MbL"])
    s = (rec.get("osc") or rec.get("osc_raw") or "no").strip()
    if s in ("no", ""):
        return F(0)
    return sum(abs(F(m)) for _sg, m, _n, _p in _OSC.findall(s))

HERE = os.path.dirname(os.path.abspath(__file__))
CAMP = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
OUT = os.path.join(CAMP, "data", "su5_ten_bound")

ROLES10 = {"Q", "uc", "ec", "Qbar", "ucbar", "ecbar"}
ROLES5 = {"dc", "L", "Hu", "dcbar"}

GEOM = {
    "Z2xZ4_1_1": dict(lib="z2xz4_sweep", mod="z24sw_lib", M=2, N=4,
                      census=os.path.join(CAMP, "data", "z2xz4_sweep", "raw", "census"),
                      v1=(F(0), F(1, 2), F(-1, 2)), v2=(F(1, 4), F(0), F(-1, 4))),
    "Z3xZ6_1_1": dict(lib="z3xz6_sweep", mod="z36sw_lib", M=3, N=6,
                      census=os.path.join(CAMP, "data", "z3xz6_sweep", "raw", "census"),
                      v1=(F(0), F(1, 3), F(-1, 3)), v2=(F(1, 6), F(0), F(-1, 6))),
}

_S = None
_R2 = None


def _init():
    global _S, _R2
    import s5lib
    _S = s5lib
    _R2 = np.array([[int(2 * x) for x in a] for a in _S.e8xe8_roots()], dtype=np.int64)


def mask_of(V):
    D = 1
    for x in V:
        D = D * x.denominator // np.gcd(D, x.denominator)
    Vi = np.array([int(x * D) for x in V], dtype=np.int64)
    return ((_R2 @ Vi) % (2 * D)) == 0


def pair_ok(a, b):
    """(sqrt(a) + sqrt(b))^2 >= 2, exactly over Q."""
    t = 2 - a - b
    if t <= 0:
        return True
    return 4 * a * b >= t * t


def worker(args):
    path, geom = args
    if _S is None:
        _init()
    G = GEOM[geom]
    sys.path.insert(0, os.path.join(CAMP, "code", G["lib"]))
    Z = __import__(G["mod"])
    S = _S
    allr = S.e8xe8_roots()

    acc = dict(models=0, skipped=Counter(), kY=Counter(),
               cells=0, fields=0,
               R_by_fam=defaultdict(Counter),
               w2_by_fam=defaultdict(Counter),
               mult_by_fam=defaultdict(Counter),
               mult_by_fam_sector=defaultdict(Counter),
               B_by_sector=defaultdict(Counter),
               violations=[], viol_n=0,
               loophole_cells=[], loophole_n=0,
               unprotected_5=0,
               I1_bad=0, I1_checked=0,
               smbasis_fail=0)
    txt = gzip.open(path, "rt", errors="replace").read()
    for lab, blk in Z.split_models(txt):
        pre = blk.split("Analyzing")[0]
        try:
            if geom == "Z2xZ4_1_1":
                spec = Z.parse_spectrum(pre, M=G["M"], N=G["N"])
                factors = Z.parse_gauge_factors(blk)
                cfgs = Z.parse_smconfig(blk)
                if not cfgs or cfgs[0][0] is None or not spec:
                    acc["skipped"]["no_smconfig"] += 1
                    continue
                tY, tbl = cfgs[0]
                printed = defaultdict(int)
                for row in tbl:
                    printed[(row[0], row[1], row[2])] += row[3]
                pos = Z.observable_positions(spec, tY, factors, dict(printed))
                if pos is None:
                    acc["skipped"]["closure_failed"] += 1
                    continue
                roles, diag = Z.role_census(spec, tY, *pos)
            else:
                spec = Z.parse_spectrum(pre)
                gauge = Z.parse_gauge(blk)
                if gauge["tY"] is None or not spec:
                    acc["skipped"]["no_hypercharge"] += 1
                    continue
                pos = Z.observable_positions(spec, gauge["tY"], gauge)
                if pos is None:
                    acc["skipped"]["closure_failed"] += 1
                    continue
                a3, a2, sY0, tY, bi = pos
                roles, diag = Z.role_census(spec, tY, a3, a2, sY0)
        except Exception as e:                                        # noqa: BLE001
            acc["skipped"]["parse_error:%s" % type(e).__name__] += 1
            continue
        acc["models"] += 1
        acc["kY"][str(2 * S.dot(tY, tY))] += 1

        cells = defaultdict(list)
        for l, r in spec.items():
            if r["p_sh"] is None or r["chir"] != "l":
                continue
            cells[(tuple(r["kl"]), tuple(r["fp"]) if r["fp"] else None)].append(r)

        # ---- 4D gauge group -> SM Cartan span
        inter = None
        for key, recs in cells.items():
            m = mask_of(recs[0]["p_sh"])
            inter = m if inter is None else (inter & m)
        R4 = [allr[i] for i in np.nonzero(inter)[0]] if inter is not None else []
        Sp4 = S.simple_roots(R4) if R4 else []
        smbasis = None
        if Sp4:
            comps4, adj4 = S._components(Sp4)
            qrec = [r for r in spec.values()
                    if roles.get(r["label"], {}).get("role") in ("Q", "Qbar")]
            drec = [r for r in spec.values()
                    if roles.get(r["label"], {}).get("role") in ("dc", "dcbar")]
            lrec = [r for r in spec.values()
                    if roles.get(r["label"], {}).get("role") in ("L", "Hu")]
            f3, f2 = None, None
            for c in comps4:
                kind, rk = S._classify(c, adj4, Sp4)
                so = [Sp4[i] for i in c]
                if kind == "A" and rk == 2 and (qrec or drec):
                    if all(any(S.dot(x["p_sh"], a) != 0 for a in so) for x in qrec + drec):
                        f3 = so if f3 is None else "AMBIG"
                if kind == "A" and rk == 1 and (qrec or lrec):
                    if all(any(S.dot(x["p_sh"], a) != 0 for a in so) for x in qrec + lrec):
                        f2 = so if f2 is None else "AMBIG"
            if isinstance(f3, list) and isinstance(f2, list):
                smbasis = list(f3) + list(f2) + [tuple(tY)]
        if smbasis is None:
            acc["smbasis_fail"] += 1
            continue

        for key, recs in cells.items():
            (kl, fp) = key
            if tuple(kl) == (0, 0):
                continue
            v = [(kl[0] * G["v1"][i] + kl[1] * G["v2"][i]) for i in range(3)]
            dc = sum((x % 1) * (1 - (x % 1)) for x in v) / 2
            rmult = Counter()
            info = {}
            for r in recs:
                ro = roles.get(r["label"], {}).get("role")
                if ro is None:
                    continue
                rmult[ro] += 1
                B = S.dot(r["p_sh"], r["p_sh"]) / 2
                Mosc = mode_sum(r)
                acc["I1_checked"] += 1
                if B != 1 - dc - Mosc:
                    acc["I1_bad"] += 1
                w2 = S.proj_norm2(r["p_sh"], smbasis)
                Rv = 2 * B - w2
                fam = "10" if ro in ROLES10 else ("5bar" if ro in ROLES5 else "?")
                acc["fields"] += 1
                acc["R_by_fam"][fam][str(Rv)] += 1
                acc["w2_by_fam"][fam][str(w2)] += 1
                acc["B_by_sector"][str(tuple(kl))][str(B)] += 1
                info.setdefault(ro, []).append((Rv, B, w2, r["label"],
                                               r.get("osc") or r.get("osc_raw")))
                if fam == "10" and pair_ok(Rv, Rv) and acc["loophole_n"] < 400:
                    acc["loophole_n"] += 1
                    acc["loophole_cells"].append(
                        dict(model=lab, kl=list(kl), fp=list(fp) if fp else None,
                             role=ro, R=str(Rv), B=str(B), w2=str(w2), file=os.path.basename(path)))
                elif fam == "10" and pair_ok(Rv, Rv):
                    acc["loophole_n"] += 1
            acc["cells"] += 1
            for ro, mm in rmult.items():
                fam = "10" if ro in ROLES10 else ("5bar" if ro in ROLES5 else "?")
                acc["mult_by_fam"][fam][mm] += 1
                acc["mult_by_fam_sector"]["%s|%s" % (tuple(kl), fam)][mm] += 1
                if mm >= 2:
                    # THE THEOREM, general form: two states with the same role at one
                    # cell require  sqrt(R_i) + sqrt(R_j) >= sqrt(2)  for SOME pair
                    # (equivalently R >= 1/2 when both are oscillator-free, R_i = R_j).
                    # Testing min(R) >= 1/2 is the EQUAL-R special case and is WRONG when
                    # the two fields sit at different oscillator levels.
                    Rs = sorted((x[0] for x in info[ro]), reverse=True)
                    ok = any(pair_ok(Rs[i], Rs[j])
                             for i in range(len(Rs)) for j in range(i + 1, len(Rs)))
                    Rv = Rs[0]
                    if not ok:
                        acc["viol_n"] += 1
                        if len(acc["violations"]) < 60:
                            acc["violations"].append(
                                dict(model=lab, kl=list(kl), fp=list(fp) if fp else None,
                                     role=ro, mult=mm, R=str(Rv),
                                     fields=[(str(a), str(b), str(c), d, e)
                                             for a, b, c, d, e in info[ro]],
                                     file=os.path.basename(path)))
    for k in ("R_by_fam", "w2_by_fam", "mult_by_fam", "mult_by_fam_sector",
              "B_by_sector"):
        acc[k] = {a: dict(b) for a, b in acc[k].items()}
    acc["skipped"] = dict(acc["skipped"])
    acc["kY"] = dict(acc["kY"])
    return acc


def merge(accs):
    out = dict(models=0, skipped=Counter(), kY=Counter(), cells=0, fields=0,
               R_by_fam=defaultdict(Counter), w2_by_fam=defaultdict(Counter),
               mult_by_fam=defaultdict(Counter),
               mult_by_fam_sector=defaultdict(Counter),
               B_by_sector=defaultdict(Counter),
               violations=[], viol_n=0, loophole_cells=[], loophole_n=0,
               I1_bad=0, I1_checked=0, smbasis_fail=0)
    for a in accs:
        for k in ("models", "cells", "fields", "viol_n", "loophole_n", "I1_bad",
                  "I1_checked", "smbasis_fail"):
            out[k] += a[k]
        for k in ("skipped", "kY"):
            out[k].update(a[k])
        for k in ("R_by_fam", "w2_by_fam", "mult_by_fam", "mult_by_fam_sector",
                  "B_by_sector"):
            for x, d in a[k].items():
                out[k][x].update({str(kk): vv for kk, vv in d.items()})
        out["violations"] += a["violations"][:20]
        out["loophole_cells"] += a["loophole_cells"][:40]
    for k in ("skipped", "kY"):
        out[k] = dict(out[k])
    for k in ("R_by_fam", "w2_by_fam", "mult_by_fam", "mult_by_fam_sector",
              "B_by_sector"):
        out[k] = {a: dict(b) for a, b in out[k].items()}
    return out


def main():
    geom = sys.argv[1] if len(sys.argv) > 1 else "Z2xZ4_1_1"
    nseg = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    files = sorted(glob.glob(os.path.join(GEOM[geom]["census"], "*.txt.gz")))
    if nseg:
        files = files[:nseg]
    tasks = [(f, geom) for f in files]
    with mp.Pool(6, initializer=_init) as pool:
        accs = pool.map(worker, tasks)
    out = merge(accs)
    out["geometry"] = geom
    out["n_segments"] = len(files)
    with open(os.path.join(OUT, "verify_%s.json" % geom), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print(json.dumps({k: v for k, v in out.items()
                      if k not in ("violations", "loophole_cells")},
                     indent=1, default=str)[:8000])
    print("VIOLATIONS:", out["viol_n"], " LOOPHOLE-CELLS(10-role,R>=1/2):", out["loophole_n"])
    for v in out["violations"][:6]:
        print("   ", v)
    for v in out["loophole_cells"][:6]:
        print("  LH", v)


if __name__ == "__main__":
    main()

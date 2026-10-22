"""root_lift QUESTION B -- the data test.

On both on-disk catalogues (Z3xZ6_1_1, Z2xZ4_1_1): for every twisted cell keyed by the
FULL (k, l, n_1..n_6) tuple, for every pair of distinct printed records sharing a
Standard-Model role, compute Delta = p - p' and test whether Delta is a root of the
UNBROKEN FOUR-DIMENSIONAL gauge group -- i.e. Delta . V_1, Delta . V_2 and Delta . A_i in Z
for ALL i = 1..6 (CLAUDE.md 5.21 (i): every defining condition, not a convenient subset).

TWO INDEPENDENT ROUTES to the 4d condition, computed and compared per model:
  route G ("generators"): recover V_1, V_2, A_1..A_6 modulo Lambda from the printed p_sh
      of the zero-translation and unit-translation cells (V_g = p_sh mod Lambda, premise P1
      of su5_ten_bound), then evaluate all eight pairings.  This is the complete condition.
  route I ("intersection"): the su5_ten_bound route -- alpha is a 4d root iff
      alpha . V_g in Z at EVERY populated cell of the model.  This is a NECESSARY condition
      that becomes sufficient when the populated cells generate V(S).
Route G is the load-bearing one; route I is the cross-check and covers models where a
generator cell is unpopulated.

Each pair is tagged with the Question A verdict for its constructing element (does the
centraliser see the whole gauge-shift set?), with whether the two records share the
right-mover q_sh, and with whether both are oscillator-free.

Output: data/root_lift/questionB_<geom>.json
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

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.dirname(HERE)
CAMP = os.path.dirname(CODE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(CODE, "su5_ten_bound"))

OUT = os.path.join(CAMP, "data", "root_lift")
os.makedirs(OUT, exist_ok=True)

ROLES10 = {"Q", "uc", "ec", "Qbar", "ucbar", "ecbar"}
ROLES5 = {"dc", "L", "Hu", "dcbar"}

GEOM = {
    "Z2xZ4_1_1": dict(lib="z2xz4_sweep", mod="z24sw_lib", M=2, N=4,
                      census=os.path.join(CAMP, "data", "z2xz4_sweep", "raw", "census")),
    "Z3xZ6_1_1": dict(lib="z3xz6_sweep", mod="z36sw_lib", M=3, N=6,
                      census=os.path.join(CAMP, "data", "z3xz6_sweep", "raw", "census")),
}

# filled per process
_S5 = None
_LIFT = None            # {(k,l,n...) : bool}
_GENCELL = None         # {name: (kl, fp)}  cells whose V_g gives a generator
_WLORD = None
_WLIDENT = None
_R2 = None
_RL = None


def _mask_of(V):
    """Boolean mask over e8xe8_roots(): alpha . V in Z.  Integer arithmetic (verify.py)."""
    import numpy as np
    D = 1
    for x in V:
        D = D * F(x).denominator // np.gcd(D, F(x).denominator)
    Vi = np.array([int(F(x) * D) for x in V], dtype=np.int64)
    return ((_R2 @ Vi) % (2 * D)) == 0


def _init(geom):
    global _S5, _LIFT, _GENCELL, _WLORD, _R2, _RL
    import numpy as np
    import s5lib
    import rl_lib as R
    _S5 = s5lib
    _RL = R
    _R2 = np.array([[int(2 * x) for x in a] for a in s5lib.e8xe8_roots()], dtype=np.int64)
    S = R.SpaceGroup(geom)
    _LIFT = {}
    for sec in S.G.sectors_ce:
        if (sec[0][0], sec[0][1]) == (0, 0):
            continue
        for ce in sec:
            _LIFT[tuple(ce)] = (S.lift_verdict(ce)["index"] == 1)
    # cells from which the eight gauge-shift generators can be recovered:
    #   V_1  <- ce (1,0,0...0)          V_2 <- ce (0,1,0...0)
    #   A_i  <- ce (k,l,n) with n = e_i, minus the (k,l,0) cell
    gen = {"V1": None, "V2": None}
    for i in range(6):
        gen["A%d" % (i + 1)] = None
    for sec in S.G.sectors_ce:
        k, l = sec[0][0], sec[0][1]
        if (k, l) == (0, 0):
            continue
        zero = tuple([k, l] + [0] * 6)
        if zero not in [tuple(c) for c in sec]:
            continue
        for ce in sec:
            n = list(ce[2:8])
            nz = [j for j in range(6) if n[j] != 0]
            if len(nz) == 1 and n[nz[0]] == 1 and gen["A%d" % (nz[0] + 1)] is None:
                gen["A%d" % (nz[0] + 1)] = (("base", (k, l), tuple([0] * 6)),
                                            ("plus", (k, l), tuple(n)))
    # V1, V2 from the pure-generator sectors
    for sec in S.G.sectors_ce:
        k, l = sec[0][0], sec[0][1]
        zero = tuple([k, l] + [0] * 6)
        if zero not in [tuple(c) for c in sec]:
            continue
        if (k, l) == (1, 0) and gen["V1"] is None:
            gen["V1"] = (("one", (k, l), tuple([0] * 6)),)
        if (k, l) == (0, 1) and gen["V2"] is None:
            gen["V2"] = (("one", (k, l), tuple([0] * 6)),)
    _GENCELL = gen
    _WLORD = list(S.G.wl_orders)
    global _WLIDENT
    _WLIDENT = [list(x) for x in S.G.identical_wl_groups]


def _reduce(v):
    """Representative of V_g.  DELIBERATELY the printed p_sh itself, unreduced.

    DEFECT CAUGHT IN THIS PHASE (2026-09-04, by the route-G / route-I cross-check):
    an earlier version reduced each component into [0,1).  That shifts the vector by an
    element of Z^16, which is NOT contained in Lambda_{E8xE8} (the integral part of the E8
    lattice is the EVEN-sum sublattice D8).  Against a SPINORIAL root -- half-integer
    components -- an odd-sum integer shift moves the pairing by 1/2, and 26 of 287 pairs
    in one census segment flipped verdict.  p_sh is already an exact representative of
    V_g + Lambda, so no reduction is needed or legitimate."""
    return [F(x) for x in v]


def worker(args):
    path, geom = args
    if _S5 is None:
        _init(geom)
    G = GEOM[geom]
    sys.path.insert(0, os.path.join(CODE, G["lib"]))
    Z = __import__(G["mod"])
    S5 = _S5
    allr = S5.e8xe8_roots()

    acc = dict(models=0, skipped=Counter(), pairs=0,
               by=defaultdict(Counter), examples=[],
               gen_coverage=Counter(), route_disagree=0, route_disagree_ex=[],
               d2=Counter(), notlambda=0, notperp=0,
               liftsector_distinct_record_sameq=0,
               liftsector_distinct_record_sameq_ex=[],
               mult_at_lift_cells=Counter())
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
                _roles, _diag = Z.role_census(spec, tY, *pos)
                a3, a2, sY0 = pos[0], pos[1], pos[2]
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
                a3, a2, sY0, tY, _bi = pos
                _roles, _diag = Z.role_census(spec, tY, a3, a2, sY0)
        except Exception as e:                                          # noqa: BLE001
            acc["skipped"]["parse_error:%s" % type(e).__name__] += 1
            continue
        acc["models"] += 1

        # ---- OWN role assignment, in the branch the library selected, carrying the
        # HIDDEN MULTIPLICITY.  z36sw_lib.role_census drops non-hidden-singlets (role
        # None); that is exactly the configuration under review here, so it cannot be
        # used.  z24sw_lib carries `mult`; this reproduces both uniformly.
        sY, sC = _diag["sY_selected"], _diag["sC_selected"]
        roles = {}
        for lab_f, r in spec.items():
            if r["rep"] is None or r["chir"] != "l" or r["p_sh"] is None:
                continue
            if a3 >= len(r["rep"]) or a2 >= len(r["rep"]):
                continue
            Y = sY * sY0 * S5.dot(tY, r["p_sh"])
            c3 = r["rep"][a3]
            c3 = (sC * c3) if abs(c3) != 1 else 1
            c2 = abs(r["rep"][a2])
            hid = 1
            for kk, dd in enumerate(r["rep"]):
                if kk not in (a3, a2):
                    hid *= abs(dd)
            roles[lab_f] = dict(role=Z.SM_ROLES.get((c3, c2, Y)), hid=hid)

        cells = defaultdict(list)
        for _l, r in spec.items():
            if r["p_sh"] is None or r["chir"] != "l":
                continue
            cells[(tuple(r["kl"]), tuple(r["fp"]) if r["fp"] else tuple([0] * 6))].append(r)

        # ---- route G: recover the eight gauge-shift generators mod Lambda
        Vg_of = {k: _reduce(v[0]["p_sh"]) for k, v in cells.items()}
        gens = {}
        for name, spec_g in _GENCELL.items():
            if spec_g is None:
                continue
            if len(spec_g) == 1:
                _t, kl, fp = spec_g[0]
                if (kl, fp) in Vg_of:
                    gens[name] = Vg_of[(kl, fp)]
            else:
                (_t0, kl0, fp0), (_t1, kl1, fp1) = spec_g
                if (kl0, fp0) in Vg_of and (kl1, fp1) in Vg_of:
                    gens[name] = [a - b for a, b in zip(Vg_of[(kl1, fp1)],
                                                        Vg_of[(kl0, fp0)])]
        # geometry-trivial Wilson lines are exactly zero
        R = _RL
        for i, o in enumerate(_WLORD):
            if o == 1:
                gens["A%d" % (i + 1)] = [F(0)] * 16
        # geometry-identified Wilson lines are equal (declared in the geometry file AND
        # re-derived from the point-group matrices in rl_lib.selftest, guard 1c)
        for grp in _WLIDENT:
            known = [j for j in grp if ("A%d" % j) in gens]
            if known:
                v = gens["A%d" % known[0]]
                for j in grp:
                    gens["A%d" % j] = v
        have_all = all(("V1" in gens, "V2" in gens) +
                       tuple("A%d" % (i + 1) in gens for i in range(6)))
        acc["gen_coverage"]["complete" if have_all else "partial"] += 1

        # ---- route I: intersection over populated cells
        import numpy as np
        inter = None
        for _k, v in Vg_of.items():
            m = _mask_of(v)
            inter = m if inter is None else (inter & m)
        R4_I = {tuple(allr[i]) for i in np.nonzero(inter)[0]} if inter is not None else set()

        # ---- SM Cartan span, for the orthogonality test
        smb = None
        try:
            smb = [tuple(x) for x in (list(a3) + list(a2) + [tuple(tY)])]
        except Exception:                                               # noqa: BLE001
            smb = None

        for key, recs in cells.items():
            (kl, fp) = key
            if tuple(kl) == (0, 0):
                continue
            ce = tuple(list(kl) + list(fp))
            lifts = _LIFT.get(ce)
            byrole = defaultdict(list)
            for r in recs:
                ro = roles.get(r["label"], {}).get("role")
                if ro in ROLES10 or ro in ROLES5:
                    byrole[ro].append(r)
            if lifts:
                for ro, rs in byrole.items():
                    fam0 = "10" if ro in ROLES10 else "5bar"
                    acc["mult_at_lift_cells"]["%s:records=%d" % (fam0, len(rs))] += 1
                    for r in rs:
                        acc["mult_at_lift_cells"][
                            "%s:hidmult=%d" % (fam0, roles[r["label"]]["hid"])] += 1
            for ro, rs in byrole.items():
                fam0 = "10" if ro in ROLES10 else "5bar"
                for r in rs:
                    acc["by"]["hidmult_all"][
                        "%s|lift=%s|hid=%d" % (fam0, lifts,
                                               roles[r["label"]]["hid"])] += 1
            for ro, rs in byrole.items():
                if len(rs) < 2:
                    continue
                fam = "10" if ro in ROLES10 else "5bar"
                for i in range(len(rs)):
                    for j in range(i + 1, len(rs)):
                        a, b = rs[i], rs[j]
                        D = [F(x) - F(y) for x, y in zip(a["p_sh"], b["p_sh"])]
                        acc["pairs"] += 1
                        if not R.in_lambda(D):
                            acc["notlambda"] += 1
                            continue
                        d2 = S5.dot(D, D)
                        acc["d2"][str(d2)] += 1
                        perp = True
                        if smb:
                            perp = all(S5.dot(D, x) == 0 for x in smb)
                        if not perp:
                            acc["notperp"] += 1
                        sameq = (tuple(a.get("qsh") or ()) == tuple(b.get("qsh") or ()))

                        def _of(rr):
                            t = (rr.get("osc") or rr.get("osc_raw") or "no")
                            return str(t).strip() in ("no", "", "None")
                        oscfree = _of(a) and _of(b)

                        def _os(rr):
                            return str(rr.get("osc") or rr.get("osc_raw") or "no").strip()
                        sameosc = _os(a) == _os(b)
                        okI = tuple(D) in R4_I if d2 == 2 else None
                        okG = None
                        if have_all:
                            okG = all(S5.dot(D, gens[n]) % 1 == 0
                                      for n in ["V1", "V2"] + ["A%d" % (t + 1)
                                                               for t in range(6)])
                        if okG is not None and okI is not None and (okG != okI):
                            acc["route_disagree"] += 1
                            if len(acc["route_disagree_ex"]) < 10:
                                acc["route_disagree_ex"].append(
                                    dict(model=lab, kl=list(kl), fp=list(fp), role=ro,
                                         okG=okG, okI=okI))
                        tag = "%s|lift=%s|sameq=%s|sameosc=%s|4droot=%s" % (
                            fam, lifts, sameq, sameosc,
                            okG if okG is not None else ("I:%s" % okI))
                        acc["by"]["pairs"][tag] += 1
                        acc["by"]["sector"]["%s|%s|lift=%s" % (fam, str(tuple(kl)),
                                                              lifts)] += 1
                        if lifts and sameq and sameosc:
                            acc["liftsector_distinct_record_sameq"] += 1
                            if len(acc["liftsector_distinct_record_sameq_ex"]) < 20:
                                acc["liftsector_distinct_record_sameq_ex"].append(
                                    dict(model=lab, kl=list(kl), fp=list(fp), role=ro,
                                         d2=str(d2), okG=okG, okI=okI,
                                         oscfree=oscfree,
                                         file=os.path.basename(path)))
                        if okG is False and len(acc["examples"]) < 40:
                            acc["examples"].append(
                                dict(model=lab, kl=list(kl), fp=list(fp), role=ro,
                                     d2=str(d2), sameq=sameq, oscfree=oscfree,
                                     sameosc=sameosc,
                                     lifts=lifts, file=os.path.basename(path),
                                     dots={n: str(S5.dot(D, gens[n]))
                                           for n in ["V1", "V2"] + ["A%d" % (t + 1)
                                                                    for t in range(6)]}))
    # make picklable
    acc["by"] = {k: dict(v) for k, v in acc["by"].items()}
    acc["skipped"] = dict(acc["skipped"])
    acc["gen_coverage"] = dict(acc["gen_coverage"])
    acc["d2"] = dict(acc["d2"])
    acc["mult_at_lift_cells"] = dict(acc["mult_at_lift_cells"])
    return acc


def merge(accs):
    out = dict(models=0, skipped=Counter(), pairs=0, by=defaultdict(Counter),
               examples=[], gen_coverage=Counter(), route_disagree=0,
               route_disagree_ex=[], d2=Counter(), notlambda=0, notperp=0,
               liftsector_distinct_record_sameq=0,
               liftsector_distinct_record_sameq_ex=[],
               mult_at_lift_cells=Counter())
    for a in accs:
        for k in ("models", "pairs", "route_disagree", "notlambda", "notperp",
                  "liftsector_distinct_record_sameq"):
            out[k] += a[k]
        for k in ("skipped", "gen_coverage", "d2", "mult_at_lift_cells"):
            out[k].update(a[k])
        for k, v in a["by"].items():
            out["by"][k].update(v)
        out["examples"].extend(a["examples"][:40])
        out["route_disagree_ex"].extend(a["route_disagree_ex"][:10])
        out["liftsector_distinct_record_sameq_ex"].extend(
            a["liftsector_distinct_record_sameq_ex"][:20])
    out["by"] = {k: dict(v) for k, v in out["by"].items()}
    for k in ("skipped", "gen_coverage", "d2", "mult_at_lift_cells"):
        out[k] = dict(out[k])
    out["examples"] = out["examples"][:60]
    return out


def main():
    for geom, G in GEOM.items():
        files = sorted(glob.glob(os.path.join(G["census"], "*.txt.gz")))
        print("=== %s : %d census files" % (geom, len(files)))
        with mp.Pool(6, initializer=_init, initargs=(geom,)) as pool:
            accs = pool.map(worker, [(f, geom) for f in files])
        res = merge(accs)
        with open(os.path.join(OUT, "questionB_%s.json" % geom), "w") as fh:
            json.dump(res, fh, indent=1, default=str)
        print("  models=%d  same-role pairs=%d  Delta not in Lambda=%d  not perp W=%d"
              % (res["models"], res["pairs"], res["notlambda"], res["notperp"]))
        print("  |Delta|^2 histogram:", res["d2"])
        print("  generator coverage:", res["gen_coverage"])
        print("  route G vs route I disagreements:", res["route_disagree"])
        print("  distinct-record same-role pairs at LIFT cells with same q_sh, "
              "both oscillator-free:", res["liftsector_distinct_record_sameq"])
        print("  pair tags:")
        for k, v in sorted(res["by"].get("pairs", {}).items()):
            print("     %-52s %d" % (k, v))
        print("  multiplicity at lift-sector cells:", res["mult_at_lift_cells"])
        print()


if __name__ == "__main__":
    main()

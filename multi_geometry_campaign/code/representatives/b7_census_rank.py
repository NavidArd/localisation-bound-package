"""Referee #8, items 2.4 and 2.5, on the complete catalogue census (both geometries).

Re-runs the deposited UNIFORM census (paperA_repro_package/.../uniform_census.py; its
parsers, Standard-Model identification, role assignment, hidden dimensions and W-basis are
imported unchanged) with `census_model` replaced by an instrument that records, for every
(populated twisted fixed point, Standard-Model representation) group of the `uniform`
convention (the one tables 1-4 of paperA_v9d print):

  * the state list (R_i, d_i, N_osc,i) -- d = hidden dimension;
  * item 2.4: the reviewer's mean-budget bound  n_st <= 1/(1 - Rbar),  Rbar the state-
    weighted mean of R (a theorem under H1-H5: sum_{i<j}|s_i - s_j|^2 >= n(n-1));
    any violation is printed (there must be none); and the below-threshold states
    (R < 1/2) that DO have a same-representation partner at their fixed point;
  * item 2.5: two rank ceilings for the A_{n-1} of proposition 4.1:
      r_lat = largest n with A_{n-1} inside the roots of Lambda orthogonal to W, taken per
              E8 factor (an A_{n-1} is connected, so it lies in one factor) -- valid for
              every group with all R_i < 1 (parts (a)-(c));
      r_loc = the same inside the roots of the LOCAL group H(V_g) orthogonal to W -- valid
              when all oscillator numbers in the group differ by integers (part (d)).
    Largest n with A_{n-1} in a simple component: A_k -> k+1, D_k -> k, E6 -> 6, E7 -> 8,
    E8 -> 9.  H(V_g) is read off any record p of the cell (alpha.V_g = alpha.p mod 1).

Output: review/data/paperA_referee8/b7_<geom>.json  (per realised (family, R_max) row:
largest multiplicity in multiplets and states, floor(1/(1-R)), max over cells of the two
rank ceilings and of min(norm, rank), mean-bound violations, mixed-group compositions).

E_6 banner parse (v2.2 repair, 2026-10-22).  The paper's tables are produced by
`uniform_census_e6.py`, which rebinds the two banner parsers so that the four models whose
E_6 factor is printed before the Standard-Model factors (MSSM967, MSSM1409 on Z2xZ4_1_1;
MSSM825, MSSM2669 on Z3xZ6_1_1) are parsed.  As first deposited, this script imported
`uniform_census` WITHOUT that patch and so covered 16,134 of the 16,138 models.  Every
worker now calls `uniform_census_e6._patch()` (imported unchanged, idempotent) before
`uniform_census.worker`, so the instrument covers the same 11,526 + 4,612 models as the
paper's tables.
"""
import os
_PKG = os.environ.get("PAPERA_PKG") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"
import glob  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402
from collections import Counter, defaultdict  # noqa: E402
from fractions import Fraction as F  # noqa: E402

import numpy as np  # noqa: E402

PKG = os.path.join(_PKG, "multi_geometry_campaign", "code", "uniform_census")
sys.path.insert(0, PKG)
import uniform_census as UC  # noqa: E402

OUT = os.path.join(_PKG, "multi_geometry_campaign", "data", "representatives")
DATA = os.path.join(_PKG, "multi_geometry_campaign", "data")
# The census is read from the package's own raw census (uniform_census.GEOM).

_TAG = [0]
_orig_newacc = UC.newacc


def tagged_newacc():
    a = _orig_newacc()
    a["_tag"] = _TAG[0] % 3          # worker builds (curr, asrun, unif) in this order
    _TAG[0] += 1
    a["_mine"] = []
    return a


UC.newacc = tagged_newacc
_CACHE = {}


def maxA(roots):
    """largest n with A_{n-1} in the root system `roots` (list of 16-tuples of F)."""
    if not roots:
        return 1
    S5 = UC._S5
    sp = S5.simple_roots(roots)
    comps, adj = S5._components(sp)
    best = 1
    for c in comps:
        kind, rk = S5._classify(c, adj, sp)
        n = {"A": rk + 1, "D": rk}.get(kind)
        if kind == "E":
            n = {6: 6, 7: 8, 8: 9}[rk]
        best = max(best, n)
    return best


def perp_mask(smb):
    D = 1
    for v in smb:
        for x in v:
            D = D * F(x).denominator // np.gcd(D, F(x).denominator)
    Bm = np.array([[int(F(x) * D) for x in v] for v in smb], dtype=np.int64).T
    return np.all((UC._R2 @ Bm) == 0, axis=1)          # _R2 = 2*roots, exact


def rank_bound(mask):
    key = mask.tobytes()
    if key not in _CACHE:
        allr = UC._S5.e8xe8_roots()
        idx = np.nonzero(mask)[0]
        r1 = [allr[i] for i in idx if i < 240]
        r2 = [allr[i] for i in idx if i >= 240]
        _CACHE[key] = (maxA(r1), maxA(r2))
    return _CACHE[key]


def my_census_model(geom, lab, spec, roles, hid, smbasis, acc, memo):
    if acc.get("_tag") != 2:
        return
    G = UC.GEOM[geom]
    S5 = UC._S5
    assert UC._R2.shape[0] == 480 and all(
        (i < 240) == (all(x == 0 for x in a[8:])) for i, a in enumerate(S5.e8xe8_roots()))
    pm = perp_mask(smbasis)
    lat = rank_bound(pm)
    cells = defaultdict(list)
    for r in spec.values():
        if r["chir"] != "l" or r["p_sh"] is None or tuple(r["kl"]) == (0, 0):
            continue
        cells[(tuple(r["kl"]), tuple(r["fp"]) if r["fp"] else None)].append(r)
    for key, recs in cells.items():
        kl = key[0]
        byrole = defaultdict(list)
        for r in recs:
            ro = roles.get(r["label"])
            if ro in UC.ROLES10 or ro in UC.ROLES5:
                byrole[ro].append(r)
        if not byrole:
            continue
        loc = rank_bound(pm & UC.mask_of(recs[0]["p_sh"]))
        for ro, rs in byrole.items():
            fam = "10" if ro in UC.ROLES10 else "5bar"
            st = []
            for r in rs:
                B = S5.dot(r["p_sh"], r["p_sh"]) / 2
                w2 = S5.proj_norm2(r["p_sh"], smbasis)
                st.append((2 * B - w2, hid.get(r["label"], 1), UC.mode_sum(r)))
            acc["_mine"].append(dict(geom=geom, model=lab, kl=list(kl), role=ro, fam=fam,
                                     st=[(str(a), d, str(o)) for a, d, o in st],
                                     lat=list(lat), loc=list(loc)))


UC.census_model = my_census_model


def job(args):
    import uniform_census_e6 as UE      # v2.2 repair: the E_6-aware banner parsers
    UE._patch()
    _TAG[0] = 0
    c, a, u = UC.worker(args)
    return u["_mine"]


def summarise(groups):
    rows = defaultdict(lambda: dict(n_groups=0, max_rec=0, max_st=0, norm=None,
                                    lat_max=0, loc_max=0, min_norm_lat=0, min_norm_loc=0,
                                    lat_tighter_cells=0, loc_tighter_cells=0,
                                    mixed=Counter()))
    viol, below_partnered, below_states = [], Counter(), Counter()
    for g in groups:
        Rs = [F(a) for a, _d, _o in g["st"]]
        Ds = [d for _a, d, _o in g["st"]]
        Os = [F(o) for _a, _d, o in g["st"]]
        n_rec, n_st = len(Rs), sum(Ds)
        Rmax = max(Rs)
        row = rows["%s|%s|%s" % (g["geom"], g["fam"], Rmax)]
        row["n_groups"] += 1
        row["max_rec"] = max(row["max_rec"], n_rec)
        row["max_st"] = max(row["max_st"], n_st)
        norm = int(1 / (1 - Rmax)) if Rmax < 1 else None
        row["norm"] = norm
        lat = max(g["lat"])
        osc_int = all((o - Os[0]).denominator == 1 for o in Os)
        loc = max(g["loc"]) if osc_int else lat
        row["lat_max"] = max(row["lat_max"], lat)
        row["loc_max"] = max(row["loc_max"], loc)
        if norm is not None:
            row["min_norm_lat"] = max(row["min_norm_lat"], min(norm, lat))
            row["min_norm_loc"] = max(row["min_norm_loc"], min(norm, loc))
            row["lat_tighter_cells"] += int(lat < norm)
            row["loc_tighter_cells"] += int(loc < norm)
        if len(set(Rs)) > 1:
            row["mixed"][str(sorted(((str(r), d) for r, d in zip(Rs, Ds)),
                                    reverse=True))] += 1
        Rbar = sum(r * d for r, d in zip(Rs, Ds)) / n_st
        if Rbar < 1 and n_st * (1 - Rbar) > 1:
            viol.append(g)
        for r, d in zip(Rs, Ds):
            if r < F(1, 2):
                below_states[g["geom"]] += d
                if n_st >= 2:
                    below_partnered["%s|%s|%s" % (g["geom"], r, Rmax)] += d
        if n_st > (lat if Rmax < 1 else 10 ** 9):
            viol.append(dict(kind="rank_lat", **g))
        if osc_int and n_st > loc and Rmax < 1:
            viol.append(dict(kind="rank_loc", **g))
    for r in rows.values():
        r["mixed"] = dict(r["mixed"])
    return dict(rows=dict(sorted(rows.items())), violations=viol[:20],
                n_violations=len(viol), below_threshold_states=dict(below_states),
                below_threshold_states_with_partner=dict(below_partnered),
                n_groups=len(groups))


def main():
    geom = sys.argv[1]
    files = sorted(glob.glob(os.path.join(UC.GEOM[geom]["census"], "*.txt.gz")))
    # duplicate-label map exactly as the deposited census (its pass 1)
    with mp.Pool(4) as pool:
        labmap = dict(pool.map(UC.labels_of, [(f, geom) for f in files]))
    seen, skipmap = set(), {}
    for f in files:
        sk = []
        for lab in labmap[f]:
            if lab in seen:
                sk.append(lab)
            else:
                seen.add(lab)
        skipmap[f] = sk
    with mp.Pool(4) as pool:
        res = pool.map(job, [(f, geom, skipmap[f]) for f in files])
    groups = [g for chunk in res for g in chunk]
    models = len({g["model"] for g in groups})
    out = summarise(groups)
    out["models_with_groups"] = models
    json.dump(out, open(os.path.join(OUT, "b7_%s.json" % geom), "w"), indent=1,
              default=str)
    print(geom, "groups", len(groups), "models", models, "violations", out["n_violations"])
    for k, v in out["rows"].items():
        print(k, {x: y for x, y in v.items() if x != "mixed"}, v["mixed"])
    print("below-threshold states with a partner:", out["below_threshold_states_with_partner"])


if __name__ == "__main__":
    main()

"""su5_ten_bound STEP 2 -- THE LOCAL DATA.

Per twisted fixed-point cell, on a sample spanning Wilson-line classes of BOTH campaign
geometries that carry the measurement (Z2xZ4_1_1, Z3xZ6_1_1):

  * the local shift V_g (= p_sh mod Lambda, premise P1),
  * the local unbroken gauge group H(V_g) and its simple factors,
  * the local irrep of every SM-role field under every factor,
  * the conformal-weight decomposition  h_tot = |p_sh|^2/2 = h_SM + h_hid,
    with h_SM obtained by EXPLICIT PROJECTION onto the SM Cartan span
    (SU(3)_c simple roots, SU(2)_L simple root, t_Y) -- never assumed from k_Y,
  * the per-cell per-role multiplicity, on the certified census cell key.

Pre-registered invariants (docs/su5_ten_bound_task.md 3.2), checked here:
  I1  |p_sh|^2/2 == 1 - delta_c(k,l) - N_osc                (budget framework)
  I2  every SM-role state is a LEVEL-1 weight of every local factor
  P1  p_sh mod Lambda constant per cell
  P2  intersection of local groups over cells == printed 4D gauge group

Output: data/su5_ten_bound/local_cells_<geom>.json(.gz)
"""
from __future__ import annotations

import glob
import gzip
import json
import os
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CAMP = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
import s5lib as S                                                       # noqa: E402

OUT = os.path.join(CAMP, "data", "su5_ten_bound")
os.makedirs(OUT, exist_ok=True)

ROLES10 = ("Q", "uc", "ec")
ROLES10B = ("Qbar", "ucbar", "ecbar")
ROLES5 = ("dc", "L", "Hu")
ROLES5B = ("dcbar",)

GEOMS = {
    "Z2xZ4_1_1": dict(lib="z2xz4_sweep", mod="z24sw_lib", M=2, N=4,
                      census=os.path.join(CAMP, "data", "z2xz4_sweep", "raw", "census"),
                      v1=(F(0), F(1, 2), F(-1, 2)), v2=(F(1, 4), F(0), F(-1, 4))),
    "Z3xZ6_1_1": dict(lib="z3xz6_sweep", mod="z36sw_lib", M=3, N=6,
                      census=os.path.join(CAMP, "data", "z3xz6_sweep", "raw", "census"),
                      v1=(F(1, 3), F(0), F(-1, 3)), v2=(F(0), F(1, 6), F(-1, 6))),
}

# ---------------------------------------------------------------- fast root machinery
_R2 = None          # 480 x 16 integer matrix, roots x 2


def roots2():
    global _R2
    if _R2 is None:
        _R2 = np.array([[int(2 * x) for x in a] for a in S.e8xe8_roots()], dtype=np.int64)
    return _R2


def local_root_mask(V):
    """Boolean mask over the 480 roots: alpha . V in Z."""
    D = 1
    for x in V:
        D = D * x.denominator // np.gcd(D, x.denominator)
    Vi = np.array([int(x * D) for x in V], dtype=np.int64)
    d = roots2() @ Vi                      # = 2D * (alpha.V)
    return (d % (2 * D)) == 0


_GCACHE = {}


def local_group(V):
    """Cached local-group analysis keyed on V mod Lambda."""
    key = S.reduce_mod_lattice(V)
    g = _GCACHE.get(key)
    if g is None:
        mask = local_root_mask(V)
        allr = S.e8xe8_roots()
        LR = [allr[i] for i in np.nonzero(mask)[0]]
        if not LR:
            g = dict(n_roots=0, rank_ss=0, factors=[], simple=[], name="U(1)^16")
        else:
            Sp = S.simple_roots(LR)
            comps, adj = S._components(Sp)
            facs = []
            for c in comps:
                kind, r = S._classify(c, adj, Sp)
                _, _, order = S._order_component(c, adj, Sp)
                facs.append(dict(name=S.group_name(kind, r), kind=kind, rank=r,
                                 idx=c, order=order))
            nm = " x ".join(sorted(f["name"] for f in facs))
            if len(Sp) < 16:
                nm += " x U(1)^%d" % (16 - len(Sp))
            g = dict(n_roots=len(LR), rank_ss=len(Sp), factors=facs, simple=Sp, name=nm)
        _GCACHE[key] = g
    return g


# ---------------------------------------------------------------- SM Cartan projection
def proj_norm2(p, basis):
    """|P_W p|^2 where W = span(basis), exact over Q."""
    n = len(basis)
    if n == 0:
        return F(0)
    G = [[S.dot(basis[i], basis[j]) for j in range(n)] for i in range(n)]
    b = [S.dot(basis[i], p) for i in range(n)]
    # solve G x = b by Gaussian elimination over Q
    A = [row[:] + [b[i]] for i, row in enumerate(G)]
    piv = []
    r = 0
    for c in range(n):
        pr = next((i for i in range(r, n) if A[i][c] != 0), None)
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(n):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c)
        r += 1
    x = [F(0)] * n
    for i, c in enumerate(piv):
        x[c] = A[i][n]
    return sum(x[i] * b[i] for i in range(n))


def delta_c(v):
    """1/2 sum_a v_a (1 - v_a), v reduced to [0,1)."""
    return sum((x % 1) * (1 - (x % 1)) for x in v) / 2


# ---------------------------------------------------------------- main
def run(geom, nmodels, tag):
    G = GEOMS[geom]
    sys.path.insert(0, os.path.join(CAMP, "code", G["lib"]))
    SW = __import__(G["mod"])
    M, N = G["M"], G["N"]

    rep = dict(geometry=geom, n_models=0, I1_checked=0, I1_bad=0, I1_bad_examples=[],
               I2_checked=0, I2_bad=0, P1_cells=0, P1_bad=0,
               P2_match=0, P2_total=0, P2_mismatch_examples=[],
               kY=Counter(), local_groups=Counter(),
               sm_factor_shape=Counter(),
               role_localrep=defaultdict(Counter),
               hSM=defaultdict(Counter), hhid=defaultdict(Counter),
               mult_vs_hhid=defaultdict(Counter),
               cells=[])
    done = 0
    for pth in sorted(glob.glob(os.path.join(G["census"], "*.txt.gz"))):
        if done >= nmodels:
            break
        txt = gzip.open(pth, "rt", errors="replace").read()
        for lab, blk in SW.split_models(txt):
            if done >= nmodels:
                break
            pre = blk.split("Analyzing")[0]
            spec = SW.parse_spectrum(pre, M=M, N=N)
            factors = SW.parse_gauge_factors(blk)
            cfgs = SW.parse_smconfig(blk)
            if not cfgs or cfgs[0][0] is None or not spec:
                continue
            tY, tbl = cfgs[0]
            printed = defaultdict(int)
            for row in tbl:
                printed[(row[0], row[1], row[2])] += row[3]
            pos = SW.observable_positions(spec, tY, factors, dict(printed))
            if pos is None:
                continue
            roles, diag = SW.role_census(spec, tY, *pos)
            done += 1
            rep["n_models"] += 1
            rep["kY"][str(2 * S.dot(tY, tY))] += 1

            # ---- cells
            cells = defaultdict(list)
            for l, r in spec.items():
                if r["p_sh"] is None or r["chir"] != "l":
                    continue
                cells[(tuple(r["kl"]), tuple(r["fp"]) if r["fp"] else None)].append(r)

            # ---- P1 + P2
            inter = None
            for key, recs in cells.items():
                rep["P1_cells"] += 1
                if len({S.reduce_mod_lattice(r["p_sh"]) for r in recs}) != 1:
                    rep["P1_bad"] += 1
                m = local_root_mask(recs[0]["p_sh"])
                inter = m if inter is None else (inter & m)
            rep["P2_total"] += 1
            allr = S.e8xe8_roots()
            R4 = [allr[i] for i in np.nonzero(inter)[0]] if inter is not None else []
            Sp4 = S.simple_roots(R4) if R4 else []
            comps4, adj4 = S._components(Sp4) if Sp4 else ([], {})
            f4 = []
            for c in comps4:
                kind, r_ = S._classify(c, adj4, Sp4)
                _, _, order = S._order_component(c, adj4, Sp4)
                f4.append(dict(name=S.group_name(kind, r_), kind=kind, rank=r_,
                               order=order))
            got = Counter(f["name"] for f in f4)
            want = Counter("%s(%d)" % (t, n) for t, n in factors)
            if got == want:
                rep["P2_match"] += 1
            elif len(rep["P2_mismatch_examples"]) < 5:
                rep["P2_mismatch_examples"].append(
                    dict(model=lab, computed=dict(got), printed=dict(want)))

            # ---- SM Cartan basis: the 4D SU(3) under which Q is a triplet, the 4D SU(2)
            #      under which Q is a doublet, and t_Y.  Identified by the ROLE data, not
            #      by orbifolder's factor ordering (5.21).
            qrec = [r for r in spec.values()
                    if roles.get(r["label"], {}).get("role") in ("Q", "Qbar")]
            drec = [r for r in spec.values()
                    if roles.get(r["label"], {}).get("role") in ("dc", "dcbar")]
            lrec = [r for r in spec.values()
                    if roles.get(r["label"], {}).get("role") in ("L", "Hu")]
            f3 = f2 = None
            for f in f4:
                if f["name"] == "SU(3)" and qrec and drec:
                    so = [Sp4[i] for i in f["order"]]
                    if all(any(S.dot(x["p_sh"], a) != 0 for a in so) for x in qrec + drec):
                        f3 = so
                if f["name"] == "SU(2)" and qrec and lrec:
                    so = [Sp4[i] for i in f["order"]]
                    if all(any(S.dot(x["p_sh"], a) != 0 for a in so)
                           for x in qrec + lrec):
                        if f2 is None:
                            f2 = so
                        else:
                            f2 = "AMBIG"
            smbasis = None
            if f3 and isinstance(f2, list):
                smbasis = list(f3) + list(f2) + [tuple(tY)]
                rep["sm_factor_shape"]["SU(3)+SU(2)+tY"] += 1
            else:
                rep["sm_factor_shape"]["UNRESOLVED(f3=%s,f2=%s)" %
                                       (bool(f3), "AMBIG" if f2 == "AMBIG" else bool(f2))] += 1

            # ---- per-cell
            for key, recs in sorted(cells.items(), key=lambda kv: str(kv[0])):
                (kl, fp) = key
                if tuple(kl) == (0, 0):
                    continue
                lroles = {r["label"]: roles.get(r["label"], {}).get("role") for r in recs}
                if not any(lroles.values()):
                    continue
                V = recs[0]["p_sh"]
                g = local_group(V)
                rep["local_groups"][g["name"]] += 1
                v = [k_ * a + l_ * b for k_, l_, a, b in
                     [(kl[0], kl[1], G["v1"][i], G["v2"][i]) for i in range(3)]]
                dc = delta_c(v)
                ent = dict(model=lab, kl=list(kl), fp=list(fp) if fp else None,
                           group=g["name"], delta_c=str(dc), fields=[])
                rmult = Counter()
                for r in recs:
                    ro = lroles.get(r["label"])
                    if ro is None:
                        continue
                    rmult[ro] += 1
                    h_tot = S.dot(r["p_sh"], r["p_sh"]) / 2
                    Nosc = sum(r["NL"]) + sum(r["NbL"])
                    # I1
                    rep["I1_checked"] += 1
                    # N_osc contributes its MODE numbers, not its occurrence count (5.21 iv)
                    Mosc = sum(F(x) for x in r["ML"]) + sum(F(x) for x in r["MbL"])
                    if h_tot != 1 - dc - Mosc:
                        rep["I1_bad"] += 1
                        if len(rep["I1_bad_examples"]) < 8:
                            rep["I1_bad_examples"].append(
                                dict(model=lab, kl=list(kl), role=ro, h=str(h_tot),
                                     dc=str(dc), Mosc=str(Mosc), osc=r["osc"]))
                    tags = []
                    for f in g["factors"]:
                        rr = S.rep_of(r["p_sh"], f, g["simple"])
                        rep["I2_checked"] += 1
                        if not rr["ok"]:
                            tags.append((f["name"], "NONINT"))
                            continue
                        if not rr["level1_ok"] and rr["level"] != 2:
                            rep["I2_bad"] += 1
                        if any(rr["hw"]):
                            tags.append((f["name"], rr.get("name") or str(rr["hw"])))
                    hSM = proj_norm2(r["p_sh"], smbasis) / 2 if smbasis else None
                    hhid = (h_tot - hSM) if hSM is not None else None
                    fam = ("10" if ro in ROLES10 + ROLES10B else
                           "5bar" if ro in ROLES5 + ROLES5B else "?")
                    ent["fields"].append(dict(role=ro, fam=fam, h=str(h_tot),
                                              hSM=str(hSM), hhid=str(hhid),
                                              Mosc=str(Mosc), osc=r["osc"],
                                              localrep=tags))
                    rep["role_localrep"][ro][
                        " + ".join("%s:%s" % t for t in tags) or "all-singlet"] += 1
                    if hSM is not None:
                        rep["hSM"][fam][str(hSM)] += 1
                        rep["hhid"][fam][str(hhid)] += 1
                for ro, mm in rmult.items():
                    fam = ("10" if ro in ROLES10 + ROLES10B else
                           "5bar" if ro in ROLES5 + ROLES5B else "?")
                    rep["mult_vs_hhid"][fam][mm] += 1
                ent["mult"] = dict(rmult)
                rep["cells"].append(ent)

    for k in ("kY", "local_groups", "sm_factor_shape"):
        rep[k] = dict(rep[k].most_common(40))
    for k in ("role_localrep", "hSM", "hhid", "mult_vs_hhid"):
        rep[k] = {a: dict(b.most_common(30)) for a, b in rep[k].items()}
    with gzip.open(os.path.join(OUT, "local_cells_%s.json.gz" % tag), "wt") as fh:
        json.dump(rep, fh, indent=1, default=str)
    summ = {k: v for k, v in rep.items() if k != "cells"}
    summ["n_cells"] = len(rep["cells"])
    print(json.dumps(summ, indent=1, default=str)[:9000])
    return rep


if __name__ == "__main__":
    geom = sys.argv[1] if len(sys.argv) > 1 else "Z2xZ4_1_1"
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    run(geom, n, sys.argv[3] if len(sys.argv) > 3 else geom)

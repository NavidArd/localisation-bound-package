"""gate_su5_bound D6 -- DIRECT test of the Section 5 weight relation on real data.

Section 5 asserts, for the three SU(5)-10 roles at one cell,

    phi(u^c, e^c) = 2 phi(Q, e^c)  (mod 1),

which is the statement that

    mu := p(u^c) + p(e^c) - 2 p(Q)   lies in the ROOT LATTICE of SU(3)_c x SU(2)_L,

since only then does mu.V_h vanish mod 1 for every projecting shift V_h (the roots of the
4d-unbroken group pair integrally with every shift -- hypothesis H2).  That in turn is
TRUE ONLY IF the three roles are weights of ONE local SU(5) 10 at the cell; if they come
from two different local multiplets, mu acquires a component ORTHOGONAL to the
Standard-Model Cartan span and the relation fails.  The document does not test this.

This script tests it three ways on the printed representatives, all exactly over Q:

  T1  |p(Q) - p(e^c)|^2 == 2   and  |p(u^c) - p(e^c)|^2 == 4
      -- representative-INDEPENDENT signatures of a single SU(5) 10 (the weights of the
         10 are e_i + e_j with e_i.e_j = delta_ij - 1/5, so those two norms are fixed
         whichever of the 6 Q-weights and 3 u^c-weights orbifolder chose to print).
  T2  mu.t_Y == 0, mu in Lambda_{E8xE8}, |mu|^2 an even integer  (necessary for T3).
  T3  P_{W-perp} mu == 0, i.e. |mu|^2 == |P_W mu|^2, with W the SM Cartan 4-plane built
      exactly as verify.py builds it (4d root system -> the SU(3) and SU(2) acting on the
      roles, plus t_Y).  This is the load-bearing condition itself.

Output: data/gate_su5_bound/d6_weightrel_<geom>.json
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

import numpy as np                                                      # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CAMP = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(CAMP, "data", "gate_su5_bound")
sys.path.insert(0, os.path.join(CAMP, "code", "su5_ten_bound"))

GEOM = {
    "Z2xZ4_1_1": dict(lib="z2xz4_sweep", mod="z24sw_lib", M=2, N=4,
                      census=os.path.join(CAMP, "data", "z2xz4_sweep", "raw", "census")),
    "Z3xZ6_1_1": dict(lib="z3xz6_sweep", mod="z36sw_lib", M=3, N=6,
                      census=os.path.join(CAMP, "data", "z3xz6_sweep", "raw", "census")),
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


def worker(args):
    path, geom, do_T3 = args
    if _S is None:
        _init()
    S = _S
    allr = S.e8xe8_roots()
    G = GEOM[geom]
    sys.path.insert(0, os.path.join(CAMP, "code", G["lib"]))
    Z = __import__(G["mod"])
    acc = dict(models=0, cells_Qe=0, cells_ue=0, cells_111=0,
               d_Qe=Counter(), d_ue=Counter(), d_Qu=Counter(),
               mu_norm2=Counter(), mu_tY_nonzero=0, mu_not_in_lattice=0,
               mu_norm2_odd=0,
               T3_checked=0, T3_pass=0, T3_fail=0, T3_nobasis=0,
               fail_examples=[])
    txt = gzip.open(path, "rt", errors="replace").read()
    for lab, blk in Z.split_models(txt):
        pre = blk.split("Analyzing")[0]
        try:
            if geom == "Z2xZ4_1_1":
                spec = Z.parse_spectrum(pre, M=G["M"], N=G["N"])
                factors = Z.parse_gauge_factors(blk)
                cfgs = Z.parse_smconfig(blk)
                if not cfgs or cfgs[0][0] is None or not spec:
                    continue
                tY, tbl = cfgs[0]
                printed = defaultdict(int)
                for row in tbl:
                    printed[(row[0], row[1], row[2])] += row[3]
                pos = Z.observable_positions(spec, tY, factors, dict(printed))
                if pos is None:
                    continue
                roles, _d = Z.role_census(spec, tY, *pos)
            else:
                spec = Z.parse_spectrum(pre)
                gauge = Z.parse_gauge(blk)
                if gauge["tY"] is None or not spec:
                    continue
                pos = Z.observable_positions(spec, gauge["tY"], gauge)
                if pos is None:
                    continue
                a3, a2, sY0, tY, _bi = pos
                roles, _d = Z.role_census(spec, tY, a3, a2, sY0)
        except Exception:                                               # noqa: BLE001
            continue
        acc["models"] += 1

        cells = defaultdict(list)
        for l, r in spec.items():
            if r["p_sh"] is None or r["chir"] != "l":
                continue
            kl = tuple(r["kl"])
            if kl == (0, 0):
                continue
            cells[(kl, tuple(r["fp"]) if r["fp"] else None)].append(r)

        smbasis = None
        if do_T3:
            inter = None
            for _k, recs in cells.items():
                m = mask_of(recs[0]["p_sh"])
                inter = m if inter is None else (inter & m)
            R4 = [allr[i] for i in np.nonzero(inter)[0]] if inter is not None else []
            Sp4 = S.simple_roots(R4) if R4 else []
            if Sp4:
                comps4, adj4 = S._components(Sp4)
                qrec = [r for r in spec.values()
                        if (roles.get(r["label"]) or {}).get("role") in ("Q", "Qbar")]
                drec = [r for r in spec.values()
                        if (roles.get(r["label"]) or {}).get("role") in ("dc", "dcbar")]
                lrec = [r for r in spec.values()
                        if (roles.get(r["label"]) or {}).get("role") in ("L", "Hu")]
                f3, f2 = None, None
                for c in comps4:
                    kind, rk = S._classify(c, adj4, Sp4)
                    so = [Sp4[i] for i in c]
                    if kind == "A" and rk == 2 and (qrec or drec):
                        if all(any(S.dot(x["p_sh"], a) != 0 for a in so)
                               for x in qrec + drec):
                            f3 = so if f3 is None else "AMBIG"
                    if kind == "A" and rk == 1 and (qrec or lrec):
                        if all(any(S.dot(x["p_sh"], a) != 0 for a in so)
                               for x in qrec + lrec):
                            f2 = so if f2 is None else "AMBIG"
                if isinstance(f3, list) and isinstance(f2, list):
                    smbasis = list(f3) + list(f2) + [tuple(tY)]

        for key, recs in cells.items():
            byrole = {}
            for r in recs:
                ro = (roles.get(r["label"]) or {}).get("role")
                if ro in ("Q", "uc", "ec"):
                    byrole.setdefault(ro, r)
            if len(byrole) < 2:
                continue
            pQ = byrole.get("Q", {}).get("p_sh")
            pu = byrole.get("uc", {}).get("p_sh")
            pe = byrole.get("ec", {}).get("p_sh")

            def n2(a, b):
                d = tuple(x - y for x, y in zip(a, b))
                return S.dot(d, d), d
            if pQ and pe:
                v, _ = n2(pQ, pe)
                acc["d_Qe"][str(v)] += 1
                acc["cells_Qe"] += 1
            if pu and pe:
                v, _ = n2(pu, pe)
                acc["d_ue"][str(v)] += 1
                acc["cells_ue"] += 1
            if pQ and pu:
                v, _ = n2(pQ, pu)
                acc["d_Qu"][str(v)] += 1
            if pQ and pu and pe:
                acc["cells_111"] += 1
                mu = tuple(a + b - 2 * c for a, b, c in zip(pu, pe, pQ))
                m2 = S.dot(mu, mu)
                acc["mu_norm2"][str(m2)] += 1
                if S.dot(mu, tuple(tY)) != 0:
                    acc["mu_tY_nonzero"] += 1
                if not S.in_lattice(mu):
                    acc["mu_not_in_lattice"] += 1
                if m2.denominator != 1 or int(m2) % 2 != 0:
                    acc["mu_norm2_odd"] += 1
                if do_T3:
                    if smbasis is None:
                        acc["T3_nobasis"] += 1
                    else:
                        acc["T3_checked"] += 1
                        if S.proj_norm2(mu, smbasis) == m2:
                            acc["T3_pass"] += 1
                        else:
                            acc["T3_fail"] += 1
                            if len(acc["fail_examples"]) < 10:
                                acc["fail_examples"].append(
                                    dict(model=lab, key=str(key), mu2=str(m2),
                                         pW2=str(S.proj_norm2(mu, smbasis))))
    for k in ("d_Qe", "d_ue", "d_Qu", "mu_norm2"):
        acc[k] = dict(acc[k])
    return acc


def main():
    geom = sys.argv[1] if len(sys.argv) > 1 else "Z2xZ4_1_1"
    nT3 = int(sys.argv[2]) if len(sys.argv) > 2 else 8      # segments to run T3 on
    files = sorted(glob.glob(os.path.join(GEOM[geom]["census"], "*.txt.gz")))
    args = [(f, geom, i < nT3) for i, f in enumerate(files)]
    with mp.Pool(6) as pool:
        accs = pool.map(worker, args)
    out = dict(geometry=geom, n_files=len(files), models=0, cells_Qe=0, cells_ue=0,
               cells_111=0, mu_tY_nonzero=0, mu_not_in_lattice=0, mu_norm2_odd=0,
               T3_checked=0, T3_pass=0, T3_fail=0, T3_nobasis=0,
               d_Qe=Counter(), d_ue=Counter(), d_Qu=Counter(), mu_norm2=Counter(),
               fail_examples=[])
    for a in accs:
        for k in ("models", "cells_Qe", "cells_ue", "cells_111", "mu_tY_nonzero",
                  "mu_not_in_lattice", "mu_norm2_odd", "T3_checked", "T3_pass",
                  "T3_fail", "T3_nobasis"):
            out[k] += a[k]
        for k in ("d_Qe", "d_ue", "d_Qu", "mu_norm2"):
            out[k].update(a[k])
        out["fail_examples"] += a["fail_examples"][:3]
    for k in ("d_Qe", "d_ue", "d_Qu", "mu_norm2"):
        out[k] = dict(out[k])
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "d6_weightrel_%s.json" % geom), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print(json.dumps(out, indent=1, default=str)[:4000])


if __name__ == "__main__":
    main()

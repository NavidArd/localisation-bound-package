"""su5_ten_bound STEP 2b -- LOCAL-GUT STRUCTURE and the u^c/e^c CO-LOCATION.

For a sample of models of BOTH geometries this records, per twisted fixed-point cell:

  * the local unbroken group H(V_g) and its simple factors,
  * the local irrep of every SM-role field (level-1, Weyl-conjugated highest weight),
    reported BOTH by name and by CONVENTION-FREE CLASS.

CONJUGATION GUARD (5.19).  The two spinor nodes of D_n are exchanged by the diagram's
outer automorphism, and this code's Bourbaki node ordering picks one of them ARBITRARILY.
"16" versus "16bar" is therefore NOT a stable label here.  What IS stable is the
CONGRUENCE CLASS -- vectorial / spinorial / adjoint -- which is what the derivation uses,
and it is what this script reports as `cls`.  The named label is carried alongside, tagged,
and no conclusion rests on it.

Also measures the (c) reproduction target directly: the per-model cell SETS of Q, u^c, e^c.

Output: data/su5_ten_bound/localgut_<geom>.json
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"

import glob                                                             # noqa: E402
import gzip                                                             # noqa: E402
import json                                                             # noqa: E402
import sys                                                              # noqa: E402
from collections import Counter, defaultdict                            # noqa: E402
from fractions import Fraction as F                                     # noqa: E402

import numpy as np                                                      # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CAMP = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
import s5lib as S                                                       # noqa: E402
from verify import GEOM, ROLES10, ROLES5                                # noqa: E402

OUT = os.path.join(CAMP, "data", "su5_ten_bound")
_R2 = np.array([[int(2 * x) for x in a] for a in S.e8xe8_roots()], dtype=np.int64)
_GC = {}


def mask_of(V):
    D = 1
    for x in V:
        D = D * x.denominator // np.gcd(D, x.denominator)
    Vi = np.array([int(x * D) for x in V], dtype=np.int64)
    return ((_R2 @ Vi) % (2 * D)) == 0


def lgroup(V):
    key = S.reduce_mod_lattice(V)
    g = _GC.get(key)
    if g is None:
        allr = S.e8xe8_roots()
        LR = [allr[i] for i in np.nonzero(mask_of(V))[0]]
        if not LR:
            g = dict(factors=[], simple=[], name="U(1)^16", rank_ss=0)
        else:
            Sp = S.simple_roots(LR)
            comps, adj = S._components(Sp)
            facs = []
            for c in comps:
                kind, r = S._classify(c, adj, Sp)
                _, _, order = S._order_component(c, adj, Sp)
                facs.append(dict(name=S.group_name(kind, r), kind=kind, rank=r,
                                 order=order))
            nm = " x ".join(sorted(f["name"] for f in facs))
            if len(Sp) < 16:
                nm += " x U(1)^%d" % (16 - len(Sp))
            g = dict(factors=facs, simple=Sp, name=nm, rank_ss=len(Sp))
        _GC[key] = g
    return g


def congruence_class(kind, rank, hw):
    """Convention-free congruence class of a level-1 highest weight.

    A_n : (sum_i i*lambda_i) mod (n+1)  -- but reported as 'k-fold antisymmetric' since
          the direction (k vs n+1-k) is the conjugation ambiguity; we report min(k, n+1-k)
          and tag whether it is self-conjugate.
    D_n : 'v' (vector), 's' (spinorial, EITHER chirality -- the outer automorphism makes
          the two indistinguishable without an extra convention), 'adj', '1'.
    E6  : '27-class' (either 27 or 27bar), '1', 'adj'.  E7: '56', '1', 'adj'.
    """
    if not any(hw):
        return "1"
    if kind == "A":
        n1 = rank + 1
        k = sum((i + 1) * hw[i] for i in range(rank)) % n1
        if sum(hw) > 1:
            return "adj" if sum(hw) == 2 and hw[0] == 1 and hw[-1] == 1 else "high"
        return "L%d" % min(k, n1 - k)
    if kind == "D":
        if sum(hw) > 1:
            return "adj/high"
        if hw[0] == 1:
            return "v"
        if hw[rank - 1] == 1 or hw[rank - 2] == 1:
            return "s"
        return "other"
    if kind == "E":
        if sum(hw) > 1:
            return "adj/high"
        if rank == 6:
            return "27class"
        if rank == 7:
            return "56"
        return "other"
    return "?"


def run(geom, nmodels, nseg=6):
    G = GEOM[geom]
    sys.path.insert(0, os.path.join(CAMP, "code", G["lib"]))
    Z = __import__(G["mod"])
    allr = S.e8xe8_roots()
    rep = dict(geometry=geom, models=0, cells=0,
               local_groups=Counter(),
               role_cls=defaultdict(Counter), role_named=defaultdict(Counter),
               ten_same_irrep=Counter(),
               coloc=Counter(), coloc_models=0,
               sm_factor_for_10=Counter(),
               cells_sample=[])
    done = 0
    for pth in sorted(glob.glob(os.path.join(G["census"], "*.txt.gz")))[:nseg]:
        if done >= nmodels:
            break
        txt = gzip.open(pth, "rt", errors="replace").read()
        for lab, blk in Z.split_models(txt):
            if done >= nmodels:
                break
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
                    roles, _ = Z.role_census(spec, tY, *pos)
                else:
                    spec = Z.parse_spectrum(pre)
                    gauge = Z.parse_gauge(blk)
                    if gauge["tY"] is None or not spec:
                        continue
                    pos = Z.observable_positions(spec, gauge["tY"], gauge)
                    if pos is None:
                        continue
                    a3, a2, sY0, tY, bi = pos
                    roles, _ = Z.role_census(spec, tY, a3, a2, sY0)
            except Exception:                                          # noqa: BLE001
                continue
            done += 1
            rep["models"] += 1

            cells = defaultdict(list)
            for l, r in spec.items():
                if r["p_sh"] is None or r["chir"] != "l":
                    continue
                cells[(tuple(r["kl"]), tuple(r["fp"]) if r["fp"] else None)].append(r)

            # ---- (c) co-location of the three 10-roles, twisted cells only
            sets = {ro: set() for ro in ("Q", "uc", "ec")}
            for key, recs in cells.items():
                if tuple(key[0]) == (0, 0):
                    continue
                for r in recs:
                    ro = roles.get(r["label"], {}).get("role")
                    if ro in sets:
                        sets[ro].add(key)
            rep["coloc_models"] += 1
            rep["coloc"]["uc_cells"] += len(sets["uc"])
            rep["coloc"]["ec_cells"] += len(sets["ec"])
            rep["coloc"]["Q_cells"] += len(sets["Q"])
            rep["coloc"]["uc_eq_ec"] += int(sets["uc"] == sets["ec"])
            rep["coloc"]["uc_ec_shared"] += len(sets["uc"] & sets["ec"])
            rep["coloc"]["uc_only_vs_ec"] += len(sets["uc"] - sets["ec"])
            rep["coloc"]["ec_only_vs_uc"] += len(sets["ec"] - sets["uc"])
            rep["coloc"]["Q_eq_uc"] += int(sets["Q"] == sets["uc"])
            rep["coloc"]["Q_subset_uc"] += int(sets["Q"] <= sets["uc"])
            rep["coloc"]["uc_minus_Q"] += len(sets["uc"] - sets["Q"])
            rep["coloc"]["Q_minus_uc"] += len(sets["Q"] - sets["uc"])
            rep["coloc"]["all_three_equal"] += int(
                sets["Q"] == sets["uc"] == sets["ec"])

            # ---- local groups / reps
            for key, recs in cells.items():
                (kl, fp) = key
                if tuple(kl) == (0, 0):
                    continue
                have = [(r, roles.get(r["label"], {}).get("role")) for r in recs]
                have = [(r, ro) for r, ro in have if ro]
                if not have:
                    continue
                g = lgroup(recs[0]["p_sh"])
                rep["local_groups"][g["name"]] += 1
                rep["cells"] += 1
                ent = dict(model=lab, kl=list(kl), fp=list(fp) if fp else None,
                           group=g["name"], fields=[])
                ten_sig = []
                for r, ro in have:
                    fam = "10" if ro in ROLES10 else ("5bar" if ro in ROLES5 else "?")
                    tags, cls = [], []
                    for f in g["factors"]:
                        rr = S.rep_of(r["p_sh"], f, g["simple"])
                        if not rr["ok"]:
                            tags.append((f["name"], "NONINT"))
                            continue
                        c = congruence_class(f["kind"], f["rank"], rr["hw"])
                        if c != "1":
                            tags.append((f["name"], rr.get("name") or str(rr["hw"])))
                            cls.append((f["name"], c))
                    sig = " + ".join("%s:%s" % t for t in cls) or "all-singlet"
                    rep["role_cls"][ro][sig] += 1
                    rep["role_named"][ro][
                        " + ".join("%s:%s" % t for t in tags) or "all-singlet"] += 1
                    ent["fields"].append(dict(role=ro, fam=fam, cls=sig,
                                              named=tags,
                                              osc=r.get("osc") or r.get("osc_raw")))
                    if fam == "10":
                        ten_sig.append(sig)
                        for nm, c in cls:
                            rep["sm_factor_for_10"]["%s:%s" % (nm, c)] += 1
                if ten_sig:
                    rep["ten_same_irrep"][
                        "n=%d all_same=%s" % (len(ten_sig),
                                              len(set(ten_sig)) == 1)] += 1
                    if len(rep["cells_sample"]) < 40:
                        rep["cells_sample"].append(ent)
    rep["local_groups"] = dict(rep["local_groups"].most_common(30))
    rep["coloc"] = dict(rep["coloc"])
    rep["ten_same_irrep"] = dict(rep["ten_same_irrep"])
    rep["sm_factor_for_10"] = dict(rep["sm_factor_for_10"].most_common(20))
    for k in ("role_cls", "role_named"):
        rep[k] = {a: dict(b.most_common(12)) for a, b in rep[k].items()}
    with open(os.path.join(OUT, "localgut_%s.json" % geom), "w") as fh:
        json.dump(rep, fh, indent=1, default=str)
    print(json.dumps({k: v for k, v in rep.items() if k != "cells_sample"},
                     indent=1, default=str)[:9000])


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else "Z2xZ4_1_1",
        int(sys.argv[2]) if len(sys.argv) > 2 else 60,
        int(sys.argv[3]) if len(sys.argv) > 3 else 6)

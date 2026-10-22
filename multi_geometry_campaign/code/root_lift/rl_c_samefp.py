"""root_lift QUESTION C -- what "the same fixed point" must mean in hypothesis H3.

Three checks, all from the primary geometry file plus the on-disk catalogues:

  C-i   V_g is a CLASS FUNCTION.  For h = (k_h, l_h, n_h),
            h g h^{-1} = ( theta_g , (I - M_g) n_h + M_h n_g ) ,
        so V_{hgh^{-1}} - V_g = A((I - M_g) n_h) + A((M_h - I) n_g), and both terms vanish
        modulo Lambda by Wilson-line consistency.  Verified numerically on every listed
        constructing element of every campaign geometry, against the DERIVED Wilson-line
        relation lattice (rl_lib.wl_structure), not against an assumed one.

  C-ii  The constructing elements orbifolder lists per sector are pairwise NON-conjugate,
        so "same printed (k, l, n_1..n_6)" and "same conjugacy class" coincide on the
        catalogues.  Checked by explicit integer solvability of the conjugation equation.

  C-iii The converse fails: DISTINCT constructing elements of one sector can carry the
        SAME local shift V_g modulo Lambda.  Measured on the catalogues.  This is why H3
        must be phrased as "the same constructing element", never "the same local shift".

Output: data/root_lift/questionC_samefp.json
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"

import glob                                                             # noqa: E402
import gzip                                                             # noqa: E402
import itertools                                                        # noqa: E402
import json                                                             # noqa: E402
import multiprocessing as mp                                            # noqa: E402
import sys                                                              # noqa: E402
from collections import Counter                                         # noqa: E402
from fractions import Fraction as F                                     # noqa: E402

import numpy as np                                                      # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.dirname(HERE)
CAMP = os.path.dirname(CODE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(CODE, "su5_ten_bound"))

import rl_lib as R                                                      # noqa: E402

OUT = os.path.join(CAMP, "data", "root_lift")
os.makedirs(OUT, exist_ok=True)

GEOM = {
    "Z2xZ4_1_1": dict(lib="z2xz4_sweep", mod="z24sw_lib", M=2, N=4,
                      census=os.path.join(CAMP, "data", "z2xz4_sweep", "raw", "census")),
    "Z3xZ6_1_1": dict(lib="z3xz6_sweep", mod="z36sw_lib", M=3, N=6,
                      census=os.path.join(CAMP, "data", "z3xz6_sweep", "raw", "census")),
}


def wl_relation_matrix(S):
    """Rows spanning the Wilson-line relation lattice { n : A(n) == 0 mod Lambda for every
    admissible model } = Z-span of the columns of (I - M(k,l))."""
    I = np.eye(6, dtype=object)
    rel = []
    for k in range(S.N):
        for l in range(S.M):
            if (k, l) == (0, 0):
                continue
            Q = I - S.mat(k, l)
            for j in range(6):
                rel.append([int(x) for x in Q[:, j]])
    return np.array(rel, dtype=object)


def geometry_checks():
    out = {}
    for lab in R.CAMPAIGN_GEOMS:
        S = R.SpaceGroup(lab)
        R_rel = wl_relation_matrix(S).T          # columns = relations
        I = np.eye(6, dtype=object)
        n_cls, n_bad = 0, 0
        conj_pairs = []
        for sec in S.G.sectors_ce:
            k, l = sec[0][0], sec[0][1]
            if (k, l) == (0, 0):
                continue
            Mg = S.mat(k, l)
            # ---- C-i : conjugation moves V_g by a relation vector
            for ce in sec:
                ng = np.array(ce[2:8], dtype=object)
                for kh in range(S.N):
                    for lh in range(S.M):
                        Mh = S.mat(kh, lh)
                        for j in range(6):
                            nh = np.zeros(6, dtype=object)
                            nh[j] = 1
                            d = (I - Mg) @ nh + (Mh - I) @ ng
                            n_cls += 1
                            if R.int_solve(R_rel, d) is None:
                                n_bad += 1
            # ---- C-ii : listed constructing elements pairwise non-conjugate
            for a, b in itertools.combinations(sec, 2):
                na = np.array(a[2:8], dtype=object)
                nb = np.array(b[2:8], dtype=object)
                for kh in range(S.N):
                    for lh in range(S.M):
                        rhs = nb - S.mat(kh, lh) @ na
                        if R.int_solve(I - Mg, rhs) is not None:
                            conj_pairs.append([list(a), list(b), kh, lh])
                            break
                    else:
                        continue
                    break
        out[lab] = dict(Ci_checks=n_cls, Ci_violations=n_bad,
                        Cii_conjugate_pairs_found=len(conj_pairs),
                        Cii_examples=conj_pairs[:5])
        print("%-12s  C-i checks %6d violations %d   |   C-ii conjugate listed pairs %d"
              % (lab, n_cls, n_bad, len(conj_pairs)))
    return out


# ------------------------------------------------------------------ C-iii on the data

def worker(args):
    path, geom = args
    G = GEOM[geom]
    sys.path.insert(0, os.path.join(CODE, G["lib"]))
    Z = __import__(G["mod"])
    acc = dict(models=0, cells=0, sector_groups=0, shift_collisions=0,
               collision_ex=[], sectors_with_collision=Counter(),
               distinct_cells_per_sector=Counter())
    txt = gzip.open(path, "rt", errors="replace").read()
    for lab, blk in Z.split_models(txt):
        pre = blk.split("Analyzing")[0]
        try:
            spec = (Z.parse_spectrum(pre, M=G["M"], N=G["N"])
                    if geom == "Z2xZ4_1_1" else Z.parse_spectrum(pre))
        except Exception:                                               # noqa: BLE001
            continue
        if not spec:
            continue
        acc["models"] += 1
        cells = {}
        for _l, r in spec.items():
            if r["p_sh"] is None or r["chir"] != "l":
                continue
            key = (tuple(r["kl"]), tuple(r["fp"]) if r["fp"] else tuple([0] * 6))
            cells.setdefault(key, list(r["p_sh"]))
        acc["cells"] += len(cells)
        bysec = {}
        for (kl, fp), p in cells.items():
            if tuple(kl) == (0, 0):
                continue
            bysec.setdefault(tuple(kl), []).append((fp, p))
        for kl, lst in bysec.items():
            acc["sector_groups"] += 1
            acc["distinct_cells_per_sector"]["%s:%d" % (str(kl), len(lst))] += 1
            for (fa, pa), (fb, pb) in itertools.combinations(lst, 2):
                d = [F(x) - F(y) for x, y in zip(pa, pb)]
                if R.in_lambda(d):
                    acc["shift_collisions"] += 1
                    acc["sectors_with_collision"][str(kl)] += 1
                    if len(acc["collision_ex"]) < 15:
                        acc["collision_ex"].append(
                            dict(model=lab, kl=list(kl), fpA=list(fa), fpB=list(fb),
                                 file=os.path.basename(path)))
    acc["sectors_with_collision"] = dict(acc["sectors_with_collision"])
    acc["distinct_cells_per_sector"] = dict(acc["distinct_cells_per_sector"])
    return acc


def main():
    R.selftest(verbose=False)
    print("[guards 1a-1c] PASS\n")
    rep = dict(geometry=geometry_checks(), data={})
    print()
    for geom, G in GEOM.items():
        files = sorted(glob.glob(os.path.join(G["census"], "*.txt.gz")))
        with mp.Pool(6) as pool:
            accs = pool.map(worker, [(f, geom) for f in files])
        tot = dict(models=0, cells=0, sector_groups=0, shift_collisions=0,
                   collision_ex=[], sectors_with_collision=Counter())
        for a in accs:
            for k in ("models", "cells", "sector_groups", "shift_collisions"):
                tot[k] += a[k]
            tot["sectors_with_collision"].update(a["sectors_with_collision"])
            tot["collision_ex"].extend(a["collision_ex"][:3])
        tot["sectors_with_collision"] = dict(tot["sectors_with_collision"])
        tot["collision_ex"] = tot["collision_ex"][:20]
        rep["data"][geom] = tot
        print("%-12s models %5d  populated cells %7d  (sector,model) groups %6d"
              % (geom, tot["models"], tot["cells"], tot["sector_groups"]))
        print("             distinct-fixed-point pairs in one sector sharing V_g mod "
              "Lambda: %d" % tot["shift_collisions"])
        print("             by sector:", tot["sectors_with_collision"])
    with open(os.path.join(OUT, "questionC_samefp.json"), "w") as fh:
        json.dump(rep, fh, indent=1, default=str)
    print("\nwrote", os.path.join(OUT, "questionC_samefp.json"))


if __name__ == "__main__":
    main()

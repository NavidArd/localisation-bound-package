"""gate_su5_bound D6 -- DIRECT MEASUREMENT of the Section 5 projection phase phi.

orbifolder's own source fixes what the printed `gamma phases` are.  `corbifold.cpp:238-271`
builds `Gamma_Centralizer` = the space-group TWIST generators (theta, omega) followed by the
six SHIFT generators (lattice translations), i.e. exactly the 8 printed slots; and
`cstate.cpp:100-112` computes, for each of them,

    gamma_h(x) = q_sh(x) . v_h  -  p_sh(x) . V_h   (mod 1).

Therefore, for two states x, y at the SAME cell with the SAME right-mover q_sh (measured
identical at 22,375/22,375 and 12,285/12,285 multi-10-role cells),

    phi_h(x, y) := (p_x - p_y) . V_h  =  gamma_h(y) - gamma_h(x)   (mod 1),

so Section 5's phi is DIRECTLY READABLE off the printed spectrum for these eight elements.
This script

  (1) tests the Section 5 identity   phi(u^c,e^c) = 2 phi(Q,e^c) (mod 1)   slot by slot on
      every cell hosting all three 10-roles -- the identity itself, on real data;
  (2) histograms phi(Q,e^c) per geometry, which is the DECLARED ATTACK #5 question: is
      phi in {0, 1/2} on Z2xZ4_1_1?
  (3) histograms every printed gamma component, to show the denominators actually realised.

Note the semantics: for the Gamma_Centralizer elements the phase is ABSORBED into the
gamma-twisted linear combination over the fixed-point orbit -- it labels the state, it does
not remove it.  A non-zero phi_h here therefore does NOT mean one of the pair is projected
out; the projection proper is by the centralizer elements that FIX the point.  The identity
(1) is nevertheless the same linear-algebra identity Section 5 asserts, so testing it here
tests the algebra exactly.

Output: data/gate_su5_bound/d6_phi_<geom>.json
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
CAMP = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(CAMP, "data", "gate_su5_bound")

GEOM = {
    "Z2xZ4_1_1": dict(lib="z2xz4_sweep", mod="z24sw_lib", M=2, N=4,
                      census=os.path.join(CAMP, "data", "z2xz4_sweep", "raw", "census")),
    "Z3xZ6_1_1": dict(lib="z3xz6_sweep", mod="z36sw_lib", M=3, N=6,
                      census=os.path.join(CAMP, "data", "z3xz6_sweep", "raw", "census")),
}
TEN = ("Q", "uc", "ec")


def worker(args):
    path, geom = args
    G = GEOM[geom]
    sys.path.insert(0, os.path.join(CAMP, "code", G["lib"]))
    Z = __import__(G["mod"])
    acc = dict(models=0, cells_all3=0, cells_Qe=0, cells_ue=0,
               identity_checked=0, identity_pass=0, identity_fail=0,
               qsh_mismatch_cells=0,
               phi_Qe=Counter(), phi_ue=Counter(), phi_Qu=Counter(),
               gamma_components=Counter(), fail_examples=[])
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
        cells = defaultdict(dict)
        for l, r in spec.items():
            if r["p_sh"] is None or r["chir"] != "l" or tuple(r["kl"]) == (0, 0):
                continue
            ro = (roles.get(l) or {}).get("role")
            if ro in TEN and r.get("gamma"):
                cells[(tuple(r["kl"]), tuple(r["fp"]) if r["fp"] else None)]\
                    .setdefault(ro, r)
        for key, byrole in cells.items():
            for r in byrole.values():
                for x in r["gamma"]:
                    acc["gamma_components"][str(F(x) % 1)] += 1
            if len(byrole) < 2:
                continue
            if len(set(str(r.get("qsh")) for r in byrole.values())) != 1:
                acc["qsh_mismatch_cells"] += 1
                continue          # phi is not gamma-difference then; excluded, counted

            def phi(x, y):        # (p_x - p_y).V_h = gamma_h(y) - gamma_h(x)
                return tuple((F(b) - F(a)) % 1
                             for a, b in zip(byrole[x]["gamma"], byrole[y]["gamma"]))

            def key8(t):
                return ",".join(str(x) for x in t)
            if "Q" in byrole and "ec" in byrole:
                acc["cells_Qe"] += 1
                acc["phi_Qe"][key8(phi("Q", "ec"))] += 1
            if "uc" in byrole and "ec" in byrole:
                acc["cells_ue"] += 1
                acc["phi_ue"][key8(phi("uc", "ec"))] += 1
            if "Q" in byrole and "uc" in byrole:
                acc["phi_Qu"][key8(phi("Q", "uc"))] += 1
            if len(byrole) == 3:
                acc["cells_all3"] += 1
                a = phi("Q", "ec")
                b = phi("uc", "ec")
                acc["identity_checked"] += 1
                if all((2 * x - y) % 1 == 0 for x, y in zip(a, b)):
                    acc["identity_pass"] += 1
                else:
                    acc["identity_fail"] += 1
                    if len(acc["fail_examples"]) < 10:
                        acc["fail_examples"].append(dict(
                            model=lab, key=str(key),
                            phi_Qe=[str(x) for x in a], phi_ue=[str(x) for x in b]))
    for k in ("phi_Qe", "phi_ue", "phi_Qu", "gamma_components"):
        acc[k] = dict(acc[k])
    return acc


def main():
    geom = sys.argv[1] if len(sys.argv) > 1 else "Z2xZ4_1_1"
    files = sorted(glob.glob(os.path.join(GEOM[geom]["census"], "*.txt.gz")))
    with mp.Pool(6) as pool:
        accs = pool.map(worker, [(f, geom) for f in files])
    out = dict(geometry=geom, models=0, cells_all3=0, cells_Qe=0, cells_ue=0,
               identity_checked=0, identity_pass=0, identity_fail=0,
               qsh_mismatch_cells=0,
               phi_Qe=Counter(), phi_ue=Counter(), phi_Qu=Counter(),
               gamma_components=Counter(), fail_examples=[])
    for a in accs:
        for k in ("models", "cells_all3", "cells_Qe", "cells_ue", "identity_checked",
                  "identity_pass", "identity_fail", "qsh_mismatch_cells"):
            out[k] += a[k]
        for k in ("phi_Qe", "phi_ue", "phi_Qu", "gamma_components"):
            out[k].update(a[k])
        out["fail_examples"] += a["fail_examples"][:3]
    # denominators actually realised in phi(Q,e^c)
    den, denu = Counter(), Counter()
    for src, dst in ((out["phi_Qe"], den), (out["phi_ue"], denu)):
        for s, n in src.items():
            for tok in s.split(","):
                tok = tok.strip()
                if tok:
                    dst[F(tok).denominator] += n
    out["phi_Qe_denominators"] = dict(den)
    out["phi_ue_denominators"] = dict(denu)
    for k in ("phi_Qe", "phi_ue", "phi_Qu", "gamma_components"):
        out[k] = dict(Counter(out[k]).most_common(20))
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "d6_phi_%s.json" % geom), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print(json.dumps(out, indent=1, default=str)[:5000])


if __name__ == "__main__":
    main()

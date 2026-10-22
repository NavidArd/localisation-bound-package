"""gate_su5_bound D6 -- is the per-state `gamma phases` slot CONSTANT among the 10-roles
of one cell?

The Section 5 argument needs the whole non-`p.V_h` part of the orbifold projection phase to
CANCEL between two 10-roles at the same cell.  The vacuum phase depends only on (g, h) and
cancels unconditionally; the right-mover term cancels iff `R = q_sh + N - Nbar` agrees; and
orbifolder additionally prints a PER-STATE `gamma phases` vector (the stabiliser character
the state carries).  If two 10-roles at a cell carry different gamma, the difference of
their projection phases is NOT `Delta.V_h` and the relation is not available at that cell.

Light-weight: a fixed number of leading census segments per geometry, single process.

Output: data/gate_su5_bound/d6_gamma_homog.json
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

HERE = os.path.dirname(os.path.abspath(__file__))
CAMP = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(CAMP, "data", "gate_su5_bound")
TEN = ("Q", "uc", "ec")

GEOM = {
    "Z2xZ4_1_1": dict(lib="z2xz4_sweep", mod="z24sw_lib", M=2, N=4,
                      census=os.path.join(CAMP, "data", "z2xz4_sweep", "raw", "census")),
    "Z3xZ6_1_1": dict(lib="z3xz6_sweep", mod="z36sw_lib", M=3, N=6,
                      census=os.path.join(CAMP, "data", "z3xz6_sweep", "raw", "census")),
}


def run(geom, nseg=6):
    G = GEOM[geom]
    sys.path.insert(0, os.path.join(CAMP, "code", G["lib"]))
    Z = __import__(G["mod"])
    rep = dict(geometry=geom, nseg=nseg, models=0, cells_multi10=0,
               gamma_same=0, gamma_diff=0, qsh_same=0, qsh_diff=0,
               gamma_diff_examples=[], gamma_values=Counter())
    for pth in sorted(glob.glob(os.path.join(G["census"], "*.txt.gz")))[:nseg]:
        txt = gzip.open(pth, "rt", errors="replace").read()
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
                    roles, _ = Z.role_census(spec, tY, *pos)
                else:
                    spec = Z.parse_spectrum(pre)
                    gauge = Z.parse_gauge(blk)
                    if gauge["tY"] is None or not spec:
                        continue
                    pos = Z.observable_positions(spec, gauge["tY"], gauge)
                    if pos is None:
                        continue
                    a3, a2, sY0, tY, _bi = pos
                    roles, _ = Z.role_census(spec, tY, a3, a2, sY0)
            except Exception:                                           # noqa: BLE001
                continue
            rep["models"] += 1
            cells = defaultdict(list)
            for l, r in spec.items():
                if r["p_sh"] is None or r["chir"] != "l":
                    continue
                if tuple(r["kl"]) == (0, 0):
                    continue
                ro = (roles.get(l) or {}).get("role")
                if ro in TEN:
                    cells[(tuple(r["kl"]), tuple(r["fp"]) if r["fp"] else None)].append(
                        (ro, l, r))
            for key, recs in cells.items():
                for _ro, _l, r in recs:
                    rep["gamma_values"][str(r.get("gamma"))] += 1
                if len(recs) < 2:
                    continue
                rep["cells_multi10"] += 1
                gs = set(str(r.get("gamma")) for _a, _b, r in recs)
                qs = set(str(r.get("qsh")) for _a, _b, r in recs)
                if len(gs) == 1:
                    rep["gamma_same"] += 1
                else:
                    rep["gamma_diff"] += 1
                    if len(rep["gamma_diff_examples"]) < 10:
                        rep["gamma_diff_examples"].append(dict(
                            model=lab, key=str(key),
                            g=[(a, str(r.get("gamma"))) for a, _b, r in recs]))
                rep["qsh_same" if len(qs) == 1 else "qsh_diff"] += 1
    rep["gamma_values"] = dict(rep["gamma_values"].most_common(8))
    return rep


if __name__ == "__main__":
    out = {}
    for g in ("Z2xZ4_1_1", "Z3xZ6_1_1"):
        out[g] = run(g, int(sys.argv[1]) if len(sys.argv) > 1 else 6)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "d6_gamma_homog.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print(json.dumps(out, indent=1, default=str)[:4000])

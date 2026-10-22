"""su5_ten_bound -- probe of the (c) reproduction target's RESIDUAL.

`localgut.py` reproduces the gate's exact u^c / e^c co-location (0 u^c-only, 0 e^c-only)
but shows Q is NOT co-located with them at a minority of cells.  This script asks WHAT IS
AT those cells: is the 10-role's local irrep still a single one, is the missing Q present
as a Qbar, or is it absent from the printed spectrum altogether?

The distinction matters for the honesty of the (c) claim: the derivation says the massless
SHELL at a cell is a union of complete local-group irreps; the PRINTED spectrum is that
shell after the orbifold projection.  If the missing Q is a projection effect the
derivation's scope is exactly stated; if it is something else the claim must be weakened.

Output: data/su5_ten_bound/coloc_probe.json
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

import numpy as np                                                      # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CAMP = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
import s5lib as S                                                       # noqa: E402
from localgut import lgroup, congruence_class                           # noqa: E402
from verify import GEOM                                                 # noqa: E402

OUT = os.path.join(CAMP, "data", "su5_ten_bound")
TEN = ("Q", "uc", "ec")
TENB = ("Qbar", "ucbar", "ecbar")


def run(geom="Z2xZ4_1_1", nmodels=400, nseg=14):
    G = GEOM[geom]
    sys.path.insert(0, os.path.join(CAMP, "code", G["lib"]))
    Z = __import__(G["mod"])
    rep = dict(geometry=geom, models=0,
               pattern_at_10cells=Counter(),
               uc_no_Q_group=Counter(), uc_no_Q_content=Counter(),
               Q_no_uc_group=Counter(), Q_no_uc_content=Counter(),
               examples=[])
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
            for key, recs in cells.items():
                if tuple(key[0]) == (0, 0):
                    continue
                pres = {}
                for r in recs:
                    ro = roles.get(r["label"], {}).get("role")
                    if ro:
                        pres.setdefault(ro, []).append(r)
                if not any(ro in pres for ro in TEN + TENB):
                    continue
                sig = "".join("1" if ro in pres else "0" for ro in TEN)
                sigb = "".join("1" if ro in pres else "0" for ro in TENB)
                rep["pattern_at_10cells"]["Q,uc,ec=%s | bars=%s" % (sig, sigb)] += 1
                if sig in ("011", "111", "100"):
                    pass
                if sig == "011" or sig == "100":
                    g = lgroup(recs[0]["p_sh"])
                    content = sorted(set(
                        roles.get(r["label"], {}).get("role") for r in recs
                        if roles.get(r["label"], {}).get("role")))
                    tgt = "uc_no_Q" if sig == "011" else "Q_no_uc"
                    rep[tgt + "_group"][g["name"]] += 1
                    rep[tgt + "_content"][",".join(content)] += 1
                    if len(rep["examples"]) < 25:
                        det = []
                        for r in recs:
                            ro = roles.get(r["label"], {}).get("role")
                            if not ro:
                                continue
                            tags = []
                            for f in g["factors"]:
                                rr = S.rep_of(r["p_sh"], f, g["simple"])
                                if rr["ok"] and any(rr["hw"]):
                                    tags.append((f["name"],
                                                 congruence_class(f["kind"], f["rank"],
                                                                  rr["hw"])))
                            det.append(dict(role=ro, osc=r["osc"],
                                            h=str(S.dot(r["p_sh"], r["p_sh"]) / 2),
                                            cls=tags))
                        rep["examples"].append(
                            dict(kind=tgt, model=lab, kl=list(key[0]),
                                 fp=list(key[1]) if key[1] else None,
                                 group=g["name"], fields=det))
    for k in ("pattern_at_10cells", "uc_no_Q_group", "uc_no_Q_content",
              "Q_no_uc_group", "Q_no_uc_content"):
        rep[k] = dict(rep[k].most_common(25))
    with open(os.path.join(OUT, "coloc_probe_%s.json" % geom), "w") as fh:
        json.dump(rep, fh, indent=1, default=str)
    print(json.dumps({k: v for k, v in rep.items() if k != "examples"},
                     indent=1, default=str)[:5000])
    print("--- examples ---")
    print(json.dumps(rep["examples"][:4], indent=1, default=str)[:3500])


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else "Z2xZ4_1_1",
        int(sys.argv[2]) if len(sys.argv) > 2 else 400,
        int(sys.argv[3]) if len(sys.argv) > 3 else 14)

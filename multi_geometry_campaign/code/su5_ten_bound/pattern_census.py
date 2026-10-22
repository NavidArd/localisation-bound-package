"""su5_ten_bound -- FULL-CATALOGUE census of the (Q, u^c, e^c) occupancy pattern per cell.

The derivation predicts two of the eight patterns to be IMPOSSIBLE, from a relation that
uses nothing but the unbrokenness of SU(3)_c x SU(2)_L in four dimensions:

  the orbifold projection phase separating u^c from e^c is exactly TWICE the phase
  separating Q from e^c (mod 1), because Q -> u^c and Q -> e^c differ by broken SU(5)
  roots that are equal modulo UNBROKEN SU(3)_c / SU(2)_L roots.

  =>  pattern 110 (Q, u^c, no e^c)  is impossible
  =>  pattern 101 (Q, e^c, no u^c)  is impossible
  (patterns 010 and 001 -- u^c and e^c separating from EACH OTHER -- are NOT excluded by
   this relation; their absence is measured, not derived.)

This script counts all eight patterns over every twisted cell of both catalogues.
Output: data/su5_ten_bound/patterns_<geom>.json
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

HERE = os.path.dirname(os.path.abspath(__file__))
CAMP = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
from verify import GEOM                                                 # noqa: E402

OUT = os.path.join(CAMP, "data", "su5_ten_bound")
TEN = ("Q", "uc", "ec")
TENB = ("Qbar", "ucbar", "ecbar")


def worker(args):
    path, geom = args
    G = GEOM[geom]
    sys.path.insert(0, os.path.join(CAMP, "code", G["lib"]))
    Z = __import__(G["mod"])
    acc = dict(models=0, cells10=0, pat=Counter(), patbar=Counter(),
               pat_by_sector=Counter(), examples=[])
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
        except Exception:                                              # noqa: BLE001
            continue
        acc["models"] += 1
        cells = defaultdict(set)
        for l, r in spec.items():
            if r["p_sh"] is None or r["chir"] != "l":
                continue
            ro = roles.get(l, {}).get("role")
            if ro:
                cells[(tuple(r["kl"]), tuple(r["fp"]) if r["fp"] else None)].add(ro)
        for key, ros in cells.items():
            if tuple(key[0]) == (0, 0):
                continue
            sig = "".join("1" if x in ros else "0" for x in TEN)
            sigb = "".join("1" if x in ros else "0" for x in TENB)
            if sig != "000":
                acc["cells10"] += 1
                acc["pat"][sig] += 1
                acc["pat_by_sector"]["%s|%s" % (tuple(key[0]), sig)] += 1
                if sig in ("110", "101", "010", "001") and len(acc["examples"]) < 30:
                    acc["examples"].append(dict(model=lab, kl=list(key[0]),
                                                fp=list(key[1]) if key[1] else None,
                                                pat=sig, roles=sorted(ros),
                                                file=os.path.basename(path)))
            if sigb != "000":
                acc["patbar"][sigb] += 1
    for k in ("pat", "patbar", "pat_by_sector"):
        acc[k] = dict(acc[k])
    return acc


def main():
    geom = sys.argv[1] if len(sys.argv) > 1 else "Z2xZ4_1_1"
    files = sorted(glob.glob(os.path.join(GEOM[geom]["census"], "*.txt.gz")))
    with mp.Pool(6) as pool:
        accs = pool.map(worker, [(f, geom) for f in files])
    out = dict(geometry=geom, models=0, cells10=0, pat=Counter(), patbar=Counter(),
               pat_by_sector=Counter(), examples=[])
    for a in accs:
        out["models"] += a["models"]
        out["cells10"] += a["cells10"]
        for k in ("pat", "patbar", "pat_by_sector"):
            out[k].update(a[k])
        out["examples"] += a["examples"][:10]
    for k in ("pat", "patbar", "pat_by_sector"):
        out[k] = dict(out[k])
    out["DERIVED_IMPOSSIBLE"] = {p: out["pat"].get(p, 0) for p in ("110", "101")}
    out["NOT_DERIVED_measured_zero"] = {p: out["pat"].get(p, 0) for p in ("010", "001")}
    with open(os.path.join(OUT, "patterns_%s.json" % geom), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print(json.dumps({k: v for k, v in out.items() if k != "pat_by_sector"},
                     indent=1, default=str)[:5000])


if __name__ == "__main__":
    main()

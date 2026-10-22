"""GATE gate_su5_bound D5 -- is the (sector, raw n_a) key ever a SPLIT of one physical cell?

d5_cellkey.py found raw-key pairs that land in the SAME PER-PLANE orbit tuple.  A per-plane
reading is COARSER than the space-group conjugacy class (the conjugating element acts on all
three planes at once), so a per-plane collision is not by itself a split.  This script
decides the question with the phase's own premise P1: the local shift

        V_g  =  p_sh  mod  Lambda_{E8xE8}

is a property of the CELL.  If two raw keys carry DIFFERENT V_g classes they are different
physical cells and the key does not split anything.  If they carry the same V_g class in
every model, the keys are indistinguishable by local data and the split cannot be excluded.

We also record, for each key, the multiset of |p_sh|^2 and of representation tuples, so the
comparison does not rest on V_g alone.

Outputs data/gate_su5_bound/d5_keysplit.json
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"

import glob                                                            # noqa: E402
import gzip                                                            # noqa: E402
import json                                                            # noqa: E402
import sys                                                             # noqa: E402
from collections import Counter, defaultdict                           # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CAMP = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(CAMP, "data", "gate_su5_bound")
sys.path.insert(0, os.path.join(CAMP, "code", "su5_ten_bound"))
import s5lib as S                                                      # noqa: E402

PAIRS = {
    "Z2xZ4_1_1": [((0, 2), (1, 0, 0, 0, 1, 0), (0, 1, 0, 0, 1, 0))],
    "Z3xZ6_1_1": [((0, 2), (0, 1, 0, 0, 0, 1), (0, 2, 0, 0, 0, 1)),
                  ((0, 4), (0, 1, 0, 0, 0, 1), (0, 2, 0, 0, 0, 1)),
                  ((2, 4), (0, 1, 0, 0, 0, 1), (0, 2, 0, 0, 0, 1)),
                  ((2, 4), (0, 1, 0, 1, 0, 1), (0, 2, 0, 1, 0, 1)),
                  ((2, 4), (0, 1, 0, 2, 0, 1), (0, 2, 0, 2, 0, 1))],
}
GEOM = {"Z2xZ4_1_1": dict(lib="z2xz4_sweep", mod="z24sw_lib", M=2, N=4,
                          census=os.path.join(CAMP, "data", "z2xz4_sweep",
                                              "raw", "census")),
        "Z3xZ6_1_1": dict(lib="z3xz6_sweep", mod="z36sw_lib", M=3, N=6,
                          census=os.path.join(CAMP, "data", "z3xz6_sweep",
                                              "raw", "census"))}


def vg(p):
    return tuple(str(x) for x in S.reduce_mod_lattice(list(p)))


def main():
    res = {}
    for geom, pairs in PAIRS.items():
        G = GEOM[geom]
        sys.path.insert(0, os.path.join(CAMP, "code", G["lib"]))
        Z = __import__(G["mod"])
        out = []
        for (kl, fpA, fpB) in pairs:
            rec = dict(sector=list(kl), keyA=list(fpA), keyB=list(fpB),
                       models=0, models_both=0,
                       Vg_A=Counter(), Vg_B=Counter(),
                       same_Vg_models=0, diff_Vg_models=0,
                       same_p2_multiset=0, diff_p2_multiset=0,
                       nfields_A=0, nfields_B=0, example=None)
            out.append(rec)
        files = sorted(glob.glob(os.path.join(G["census"], "*.txt.gz")))[:3]
        for path in files:
            txt = gzip.open(path, "rt", errors="replace").read()
            for lab, blk in Z.split_models(txt):
                pre = blk.split("Analyzing")[0]
                try:
                    spec = Z.parse_spectrum(pre, M=2, N=4) if geom == "Z2xZ4_1_1" \
                        else Z.parse_spectrum(pre)
                except Exception:                                      # noqa: BLE001
                    continue
                buckets = defaultdict(list)
                for l, r in spec.items():
                    if r["p_sh"] is None or r["fp"] is None:
                        continue
                    buckets[(tuple(r["kl"]), tuple(r["fp"]))].append(r)
                for rec, (kl, fpA, fpB) in zip(out, pairs):
                    A = buckets.get((kl, fpA), [])
                    B = buckets.get((kl, fpB), [])
                    rec["models"] += 1
                    if not A or not B:
                        continue
                    rec["models_both"] += 1
                    rec["nfields_A"] += len(A)
                    rec["nfields_B"] += len(B)
                    va = {vg(r["p_sh"]) for r in A}
                    vb = {vg(r["p_sh"]) for r in B}
                    rec["Vg_A"][str(sorted(va))] += 1
                    rec["Vg_B"][str(sorted(vb))] += 1
                    if va == vb:
                        rec["same_Vg_models"] += 1
                    else:
                        rec["diff_Vg_models"] += 1
                        if rec["example"] is None:
                            rec["example"] = dict(model=lab,
                                                  file=os.path.basename(path),
                                                  Vg_A=sorted(va)[:2],
                                                  Vg_B=sorted(vb)[:2])
                    pa = Counter(str(sum(x * x for x in r["p_sh"])) for r in A)
                    pb = Counter(str(sum(x * x for x in r["p_sh"])) for r in B)
                    if pa == pb:
                        rec["same_p2_multiset"] += 1
                    else:
                        rec["diff_p2_multiset"] += 1
        for rec in out:
            rec["Vg_A_distinct"] = len(rec["Vg_A"])
            rec["Vg_B_distinct"] = len(rec["Vg_B"])
            rec["Vg_A"] = dict(list(rec["Vg_A"].items())[:3])
            rec["Vg_B"] = dict(list(rec["Vg_B"].items())[:3])
        res[geom] = out
        print("=====", geom)
        for rec in out:
            print(json.dumps({k: v for k, v in rec.items()
                              if k not in ("Vg_A", "Vg_B")}, default=str))
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "d5_keysplit.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=str)


if __name__ == "__main__":
    main()

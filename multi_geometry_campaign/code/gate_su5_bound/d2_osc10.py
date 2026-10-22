"""GATE su5_ten_bound -- D1/D2 follow-up: the OSCILLATOR 10-roles on Z3xZ6_1_1.

d2_pshdup.py found that two distinct printed left-chiral fields CAN share one p_sh at one
cell (601 groups / 600 models), differing only in oscillator content -- so the proof's
justification for `p'' != p` ("the massless shell is partitioned among the printed 4D
irreps") is false as written.  In the sample none of those groups carried a Standard-Model
role.  This script closes the question on the FULL catalogue for the dangerous case:

  * every left-chiral 10-role field carrying an oscillator, with its cell, p_sh, mode sum,
    budget B, R, oscillator string and gamma phases;
  * whether any two of them sit at the SAME cell in the SAME role (which would be a
    doubling the theorem forbids at R = 17/90 -- a live counterexample if Delta != 0, and a
    hole in the theorem if Delta == 0);
  * the same for 5bar-roles, as a control.

Output: data/gate_su5_bound/d2_osc10_<geom>.json
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"

import glob                                                             # noqa: E402
import gzip                                                             # noqa: E402
import json                                                             # noqa: E402
import multiprocessing as mp                                            # noqa: E402
import re                                                               # noqa: E402
import sys                                                              # noqa: E402
from collections import Counter, defaultdict                            # noqa: E402
from fractions import Fraction as F                                     # noqa: E402

import numpy as np                                                      # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CAMP = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(CAMP, "code", "su5_ten_bound"))
OUT = os.path.join(CAMP, "data", "gate_su5_bound")

ROLES10 = {"Q", "uc", "ec", "Qbar", "ucbar", "ecbar"}
ROLES5 = {"dc", "L", "Hu", "dcbar"}
CENSUS = os.path.join(CAMP, "data", "z3xz6_sweep", "raw", "census")
V1 = (F(0), F(1, 3), F(-1, 3))
V2 = (F(1, 6), F(0), F(-1, 6))
_OSC = re.compile(r"(-?)a_\(\s*(-?\d+(?:/\d+)?)\s*\)\^(?:\((-\d)\)|(\d))")

_S = None
_R2 = None


def _init():
    global _S, _R2
    import s5lib
    _S = s5lib
    _R2 = np.array([[int(2 * x) for x in a] for a in _S.e8xe8_roots()], dtype=np.int64)


def worker(path):
    if _S is None:
        _init()
    sys.path.insert(0, os.path.join(CAMP, "code", "z3xz6_sweep"))
    Z = __import__("z36sw_lib")
    S = _S
    allr = S.e8xe8_roots()
    acc = dict(models=0, osc10=[], n_osc10=0, n_osc5=0, same_cell_role_osc10=[],
               skipped=Counter())
    txt = gzip.open(path, "rt", errors="replace").read()
    for lab, blk in Z.split_models(txt):
        pre = blk.split("Analyzing")[0]
        try:
            spec = Z.parse_spectrum(pre)
            gauge = Z.parse_gauge(blk)
            if gauge["tY"] is None or not spec:
                acc["skipped"]["no_hypercharge"] += 1
                continue
            pos = Z.observable_positions(spec, gauge["tY"], gauge)
            if pos is None:
                acc["skipped"]["closure_failed"] += 1
                continue
            a3i, a2i, sY0, tY, bi = pos
            roles, _d = Z.role_census(spec, tY, a3i, a2i, sY0)
        except Exception as e:                                          # noqa: BLE001
            acc["skipped"]["parse:%s" % type(e).__name__] += 1
            continue
        acc["models"] += 1
        cells = defaultdict(list)
        for l, r in spec.items():
            if r["p_sh"] is None or r["chir"] != "l":
                continue
            cells[(tuple(r["kl"]), tuple(r["fp"]) if r["fp"] else None)].append(r)
        # SM Cartan span
        inter = None
        for key, recs in cells.items():
            V = recs[0]["p_sh"]
            D = 1
            for x in V:
                D = D * x.denominator // np.gcd(D, x.denominator)
            m = ((_R2 @ np.array([int(x * D) for x in V], dtype=np.int64))
                 % (2 * D)) == 0
            inter = m if inter is None else (inter & m)
        R4 = [allr[i] for i in np.nonzero(inter)[0]] if inter is not None else []
        if not R4:
            continue
        Sp4 = S.simple_roots(R4)
        comps4, adj4 = S._components(Sp4)
        qrec = [r for r in spec.values()
                if roles.get(r["label"], {}).get("role") in ("Q", "Qbar")]
        drec = [r for r in spec.values()
                if roles.get(r["label"], {}).get("role") in ("dc", "dcbar")]
        lrec = [r for r in spec.values()
                if roles.get(r["label"], {}).get("role") in ("L", "Hu")]
        f3 = f2 = None
        for c_ in comps4:
            kind, rk = S._classify(c_, adj4, Sp4)
            so = [Sp4[i] for i in c_]
            if kind == "A" and rk == 2 and (qrec or drec):
                if all(any(S.dot(x["p_sh"], a) != 0 for a in so) for x in qrec + drec):
                    f3 = so if f3 is None else "AMBIG"
            if kind == "A" and rk == 1 and (qrec or lrec):
                if all(any(S.dot(x["p_sh"], a) != 0 for a in so) for x in qrec + lrec):
                    f2 = so if f2 is None else "AMBIG"
        if not (isinstance(f3, list) and isinstance(f2, list)):
            continue
        smbasis = list(f3) + list(f2) + [tuple(tY)]

        for key, recs in cells.items():
            (kl, fp) = key
            if tuple(kl) == (0, 0):
                continue
            v = [(kl[0] * V1[i] + kl[1] * V2[i]) for i in range(3)]
            dc = sum((x % 1) * (1 - (x % 1)) for x in v) / 2
            byrole = defaultdict(list)
            for r in recs:
                ro = roles.get(r["label"], {}).get("role")
                if ro is None:
                    continue
                osc = (r.get("osc") or r.get("osc_raw") or "no").strip()
                has = osc not in ("no", "")
                if ro in ROLES5 and has:
                    acc["n_osc5"] += 1
                if ro not in ROLES10 or not has:
                    if ro in ROLES10:
                        byrole[ro].append(r)
                    continue
                acc["n_osc10"] += 1
                byrole[ro].append(r)
                ms = sum(abs(F(m)) for _s, m, _n, _p in _OSC.findall(osc))
                B = S.dot(r["p_sh"], r["p_sh"]) / 2
                w2 = S.proj_norm2(r["p_sh"], smbasis)
                if len(acc["osc10"]) < 200:
                    acc["osc10"].append(dict(
                        model=lab, kl=list(kl), fp=list(fp) if fp else None, role=ro,
                        osc=osc, mode_sum=str(ms), B=str(B), dc=str(dc),
                        I1_ok=(B == 1 - dc - ms), w2=str(w2), R=str(2 * B - w2),
                        p_sh=[str(x) for x in r["p_sh"]],
                        gamma=[str(x) for x in (r.get("gamma") or ())],
                        n_at_cell_same_role=None, file=os.path.basename(path)))
            # any (cell, role) with >= 2 records where at least one has an oscillator
            for ro, rs in byrole.items():
                if len(rs) < 2:
                    continue
                oscs = [(r.get("osc") or r.get("osc_raw") or "no").strip() for r in rs]
                acc["same_cell_role_osc10"].append(dict(
                    model=lab, kl=list(kl), fp=list(fp) if fp else None, role=ro,
                    n=len(rs), osc=oscs,
                    p_sh=[[str(x) for x in r["p_sh"]] for r in rs],
                    same_psh=(len({tuple(r["p_sh"]) for r in rs}) < len(rs)),
                    file=os.path.basename(path)))
    acc["skipped"] = dict(acc["skipped"])
    return acc


def main():
    files = sorted(glob.glob(os.path.join(CENSUS, "*.txt.gz")))
    if len(sys.argv) > 1:
        files = files[:int(sys.argv[1])]
    with mp.Pool(6, initializer=_init) as pool:
        accs = pool.map(worker, files)
    out = dict(models=0, n_osc10=0, n_osc5=0, osc10=[], same_cell_role_osc10=[],
               skipped=Counter())
    for a in accs:
        for k in ("models", "n_osc10", "n_osc5"):
            out[k] += a[k]
        out["osc10"] += a["osc10"][:40]
        out["same_cell_role_osc10"] += a["same_cell_role_osc10"][:40]
        out["skipped"].update(a["skipped"])
    out["skipped"] = dict(out["skipped"])
    out["n_10role_doubled_cells_found"] = len(out["same_cell_role_osc10"])
    out["osc10_summary"] = dict(
        by_osc=dict(Counter(e["osc"] for e in out["osc10"])),
        by_sector=dict(Counter(str(tuple(e["kl"])) for e in out["osc10"])),
        by_role=dict(Counter(e["role"] for e in out["osc10"])),
        by_R=dict(Counter(e["R"] for e in out["osc10"])),
        by_B=dict(Counter(e["B"] for e in out["osc10"])),
        I1_ok=dict(Counter(e["I1_ok"] for e in out["osc10"])),
        distinct_cells=len({(e["model"], tuple(e["kl"]), tuple(e["fp"] or ()))
                            for e in out["osc10"]}),
    )
    with open(os.path.join(OUT, "d2_osc10_Z3xZ6_1_1.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print(json.dumps({k: v for k, v in out.items()
                      if k not in ("osc10", "same_cell_role_osc10")},
                     indent=1, default=str))
    print("\n10-ROLE DOUBLED CELLS:", len(out["same_cell_role_osc10"]))
    for e in out["same_cell_role_osc10"][:5]:
        print("  ", {k: e[k] for k in ("model", "kl", "fp", "role", "n", "osc", "same_psh")})
    print("\nOSC 10-ROLE FIELDS (first 6):")
    for e in out["osc10"][:6]:
        print("  ", {k: e[k] for k in ("model", "kl", "fp", "role", "osc", "mode_sum",
                                       "B", "R", "I1_ok")})


if __name__ == "__main__":
    main()

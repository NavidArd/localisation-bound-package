"""GATE su5_ten_bound -- D7 completeness: EVERY unprotected cell, not a capped sample.

The phase's loophole list is truncated twice: verify.py records at most 400 loophole
entries per worker AND merge() keeps only the first 40 per worker accumulator; shell.py
then runs with maxcells=40.  The findings report 65 unprotected FIELDS but shell-enumerate
40 cells.  Completeness by search, not by scoping, is the campaign's recurring lesson
(CLAUDE.md 5.22), so this script enumerates ALL of them with no cap, and runs an
INDEPENDENT shell enumerator (d7_shell.py, self-tested against the E8 theta series).

For every twisted cell of Z3xZ6_1_1 carrying a 10-role field with R = 2B - |w|^2 >= 1/2
(i.e. a field the theorem does NOT protect) it records the cell, the local group, whether
the local group contains an SU(2) in the Standard-Model commutant, and the EXACT
multiplicity of that field's SM key in the complete pre-projection massless shell.

It also verifies the phase's design decision to skip oscillator states: with the minimal
oscillator mode 1/6 the budget drops to 25/36, giving R = 17/90 < 1/2, so an oscillator
cell can never be unprotected.  That is checked, not assumed.

Output: data/gate_su5_bound/d7_loophole_full.json
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

import numpy as np                                                      # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CAMP = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(CAMP, "code", "su5_ten_bound"))
sys.path.insert(0, HERE)
OUT = os.path.join(CAMP, "data", "gate_su5_bound")

ROLES10 = {"Q", "uc", "ec", "Qbar", "ucbar", "ecbar"}
CENSUS = os.path.join(CAMP, "data", "z3xz6_sweep", "raw", "census")
V1 = (F(0), F(1, 3), F(-1, 3))
V2 = (F(1, 6), F(0), F(-1, 6))

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
    import d7_shell as SH
    sys.path.insert(0, os.path.join(CAMP, "code", "z3xz6_sweep"))
    Z = __import__("z36sw_lib")
    S = _S
    allr = S.e8xe8_roots()
    acc = dict(models=0, cells=0, ten_fields=0, unprotected=[], skipped=Counter(),
               R10=Counter(), osc_unprotected=0)
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
            acc["skipped"]["no_4d_roots"] += 1
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
            acc["skipped"]["smbasis_fail"] += 1
            continue
        smbasis = list(f3) + list(f2) + [tuple(tY)]

        for key, recs in cells.items():
            (kl, fp) = key
            if tuple(kl) == (0, 0):
                continue
            acc["cells"] += 1
            v = [(kl[0] * V1[i] + kl[1] * V2[i]) for i in range(3)]
            dc = sum((x % 1) * (1 - (x % 1)) for x in v) / 2
            for r in recs:
                ro = roles.get(r["label"], {}).get("role")
                if ro not in ROLES10:
                    continue
                acc["ten_fields"] += 1
                B = S.dot(r["p_sh"], r["p_sh"]) / 2
                w2 = S.proj_norm2(r["p_sh"], smbasis)
                Rv = 2 * B - w2
                acc["R10"][str(Rv)] += 1
                if Rv < F(1, 2):
                    continue
                osc = (r.get("osc") or r.get("osc_raw") or "no").strip()
                if osc not in ("no", ""):
                    acc["osc_unprotected"] += 1
                # ---------- the loophole test, exactly and independently
                V = list(recs[0]["p_sh"])
                c = S.dot(r["p_sh"], r["p_sh"])
                sts = SH.shell(V, c)
                stset = set(sts)
                inshell = tuple(r["p_sh"]) in stset
                cnt = Counter()
                cnt10 = Counter()
                for x in sts:
                    l3 = tuple(S.dot(x, a) for a in f3)
                    l2 = tuple(S.dot(x, a) for a in f2)
                    if any(t.denominator != 1 for t in l3 + l2):
                        continue
                    kk = str((tuple(int(t) for t in l3), tuple(int(t) for t in l2),
                              str(S.dot(tY, x))))
                    cnt[kk] += 1
                    # restrict to keys that ARE 10-roles: |P_W x|^2 = 6/5
                    if S.proj_norm2(x, smbasis) == F(6, 5):
                        cnt10[kk] += 1
                mykey = str((tuple(int(S.dot(r["p_sh"], a)) for a in f3),
                             tuple(int(S.dot(r["p_sh"], a)) for a in f2),
                             str(S.dot(tY, r["p_sh"]))))
                g = S.analyse_group(V)
                # a commutant SU(2): an SU(2) factor whose root is orthogonal to the SM span
                comsu2 = False
                for fdesc in g["factors"]:
                    if fdesc["name"] != "SU(2)":
                        continue
                    a = g["simple"][fdesc["idx"][0]]
                    if all(S.dot(a, b) == 0 for b in smbasis):
                        comsu2 = True
                # any root at all orthogonal to W and integral on V_g (the real condition 2)
                anyperp = any(all(S.dot(a, b) == 0 for b in smbasis)
                              for a in S.local_roots(V, allr))
                acc["unprotected"].append(dict(
                    model=lab, kl=list(kl), fp=list(fp) if fp else None, role=ro,
                    R=str(Rv), B=str(B), osc=osc, dc=str(dc),
                    group=g["name"], shell_size=len(sts),
                    printed_in_own_shell=bool(inshell),
                    shell_mult=cnt.get(mykey, 0),
                    shell_max_mult=(max(cnt.values()) if cnt else 0),
                    shell_max_mult_10role=(max(cnt10.values()) if cnt10 else 0),
                    n_10role_states_in_shell=sum(cnt10.values()),
                    commutant_SU2=bool(comsu2),
                    any_root_perp_W_local=bool(anyperp),
                    file=os.path.basename(path)))
    acc["skipped"] = dict(acc["skipped"])
    acc["R10"] = dict(acc["R10"])
    return acc


def main():
    files = sorted(glob.glob(os.path.join(CENSUS, "*.txt.gz")))
    if len(sys.argv) > 1:
        files = files[:int(sys.argv[1])]
    with mp.Pool(6, initializer=_init) as pool:
        accs = pool.map(worker, files)
    out = dict(models=0, cells=0, ten_fields=0, unprotected=[], skipped=Counter(),
               R10=Counter(), osc_unprotected=0)
    for a in accs:
        for k in ("models", "cells", "ten_fields", "osc_unprotected"):
            out[k] += a[k]
        out["unprotected"] += a["unprotected"]
        out["skipped"].update(a["skipped"])
        out["R10"].update(a["R10"])
    out["skipped"] = dict(out["skipped"])
    out["R10"] = dict(out["R10"])
    up = out["unprotected"]
    out["n_unprotected_fields"] = len(up)
    out["n_unprotected_cells"] = len({(u["model"], tuple(u["kl"]),
                                       tuple(u["fp"] or ())) for u in up})
    out["summary"] = dict(
        by_sector=dict(Counter(str(tuple(u["kl"])) for u in up)),
        by_role=dict(Counter(u["role"] for u in up)),
        shell_mult=dict(Counter(u["shell_mult"] for u in up)),
        shell_max_mult=dict(Counter(u["shell_max_mult"] for u in up)),
        shell_max_mult_10role=dict(Counter(u["shell_max_mult_10role"] for u in up)),
        n_10role_states_in_shell=dict(Counter(u["n_10role_states_in_shell"] for u in up)),
        commutant_SU2=dict(Counter(u["commutant_SU2"] for u in up)),
        any_root_perp_W=dict(Counter(u["any_root_perp_W_local"] for u in up)),
        printed_in_own_shell=dict(Counter(u["printed_in_own_shell"] for u in up)),
        groups=dict(Counter(u["group"] for u in up)),
        R=dict(Counter(u["R"] for u in up)),
        osc=dict(Counter(u["osc"] for u in up)),
    )
    with open(os.path.join(OUT, "d7_loophole_full.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print(json.dumps({k: v for k, v in out.items() if k != "unprotected"},
                     indent=1, default=str)[:6000])


if __name__ == "__main__":
    main()

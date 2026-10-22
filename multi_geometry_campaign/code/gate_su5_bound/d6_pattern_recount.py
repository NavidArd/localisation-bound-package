"""gate_su5_bound D6 -- INDEPENDENT recount of the (Q, u^c, e^c) occupancy pattern census.

Re-counts, from the raw orbifolder census text, the 7-pattern table printed in
`results/su5_ten_bound_findings.md` Section 6.3.  Written from scratch against the
certified parsers (`z24sw_lib`, `z36sw_lib`); it does NOT import `pattern_census.py`.

Three independent cell KEYINGS are counted simultaneously, because CLAUDE.md 5.22 first
addendum (ii) makes the key the load-bearing choice:

  key F  (kl, full printed `fixed point n_a` 6-tuple)      -- what the phase used
  key O  (kl, per-plane ORBIT-CLASS index)                 -- the physical fixed point:
         two different printed n_a in the SAME space-group orbit are the SAME cell.
         If the phase's key OVER-SPLITS, key O merges cells and can MANUFACTURE the
         two patterns the derivation forbids (110 = 100 u 010, 101 = 100 u 001).
  key S  (kl only)                                         -- the deliberately INCOMPLETE
         key, counted only to size the merging hazard.  Never a physics statement.

It additionally tests, per cell, the load-bearing assumption of the Section 5 phase
argument that the document never states: that the 10-role fields at a cell share ONE
right-mover `q_sh` and ONE oscillator content, so that the vacuum/gamma/right-mover part
of the orbifold projection phase CANCELS in the difference between two of them.

Usage:  python d6_pattern_recount.py <geom> [n_files]
Output: data/gate_su5_bound/d6_patterns_<geom>.json
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"

import glob                                                             # noqa: E402
import gzip                                                             # noqa: E402
import json                                                             # noqa: E402
import multiprocessing as mp                                            # noqa: E402
import random                                                           # noqa: E402
import sys                                                              # noqa: E402
from collections import Counter, defaultdict                            # noqa: E402

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
TENB = ("Qbar", "ucbar", "ecbar")
PATS = ("111", "110", "101", "011", "100", "010", "001")


def _orbit_key(Z, geom, kl, fp, cache):
    """Per-plane orbit-class index of the printed fixed point; the PHYSICAL cell label."""
    k = (tuple(kl), tuple(fp) if fp else None)
    if k in cache:
        return cache[k]
    try:
        if geom == "Z2xZ4_1_1":
            o = tuple(Z.leg_orbit_index(kl, fp, a) for a in (1, 2, 3))
        else:
            r = [Z.leg_orbit(kl, fp, a) for a in (1, 2, 3)]
            o = tuple(None if x is None else x[0] for x in r)
    except Exception as e:                                              # noqa: BLE001
        o = ("ERR", str(e)[:40])
    o = (tuple(kl), o)
    cache[k] = o
    return o


def worker(args):
    path, geom, want_orbit = args
    G = GEOM[geom]
    sys.path.insert(0, os.path.join(CAMP, "code", G["lib"]))
    Z = __import__(G["mod"])
    acc = dict(models=0, models_skipped=0,
               patF=Counter(), patO=Counter(), patS=Counter(), patbarF=Counter(),
               cellsF=0, cellsO=0, cellsS=0,
               by_sector=Counter(), orbit_err=0,
               n_fp_per_kl_collapsed=0,           # cells that key O merged
               qsh_same=0, qsh_diff=0, osc_same=0, osc_diff=0,
               qsh_diff_examples=[], hits=[])
    ocache = {}
    txt = gzip.open(path, "rt", errors="replace").read()
    for lab, blk in Z.split_models(txt):
        pre = blk.split("Analyzing")[0]
        try:
            if geom == "Z2xZ4_1_1":
                spec = Z.parse_spectrum(pre, M=G["M"], N=G["N"])
                factors = Z.parse_gauge_factors(blk)
                cfgs = Z.parse_smconfig(blk)
                if not cfgs or cfgs[0][0] is None or not spec:
                    acc["models_skipped"] += 1
                    continue
                tY, tbl = cfgs[0]
                printed = defaultdict(int)
                for row in tbl:
                    printed[(row[0], row[1], row[2])] += row[3]
                pos = Z.observable_positions(spec, tY, factors, dict(printed))
                if pos is None:
                    acc["models_skipped"] += 1
                    continue
                roles, _d = Z.role_census(spec, tY, *pos)
            else:
                spec = Z.parse_spectrum(pre)
                gauge = Z.parse_gauge(blk)
                if gauge["tY"] is None or not spec:
                    acc["models_skipped"] += 1
                    continue
                pos = Z.observable_positions(spec, gauge["tY"], gauge)
                if pos is None:
                    acc["models_skipped"] += 1
                    continue
                a3, a2, sY0, tY, _bi = pos
                roles, _d = Z.role_census(spec, tY, a3, a2, sY0)
        except Exception:                                               # noqa: BLE001
            acc["models_skipped"] += 1
            continue
        acc["models"] += 1

        # ---- gather role-carrying LEFT-chiral fields, per cell, under three keyings
        cF = defaultdict(list)
        cS = defaultdict(list)
        cO = defaultdict(list)
        for l, r in spec.items():
            if r["p_sh"] is None or r["chir"] != "l":
                continue
            ro = (roles.get(l) or {}).get("role")
            if not ro:
                continue
            kl = tuple(r["kl"])
            if kl == (0, 0):
                continue
            fp = tuple(r["fp"]) if r["fp"] else None
            rec = (ro, l, fp, r.get("qsh"), r.get("osc") or r.get("osc_raw"))
            cF[(kl, fp)].append(rec)
            cS[kl].append(rec)
            if want_orbit:
                cO[_orbit_key(Z, geom, kl, fp, ocache)].append(rec)

        def tally(d, patc, cellc_key, bar=None, sector=False, collect=False):
            n = 0
            for key, recs in d.items():
                ros = set(x[0] for x in recs)
                sig = "".join("1" if x in ros else "0" for x in TEN)
                if sig != "000":
                    n += 1
                    patc[sig] += 1
                    if sector:
                        acc["by_sector"]["%s|%s" % (key[0], sig)] += 1
                    if collect and sig in ("110", "101", "010", "001"):
                        tens = [x for x in recs if x[0] in TEN]
                        acc["hits"].append(dict(
                            model=lab, key=str(key), pat=sig,
                            roles=sorted(ros),
                            ten_fields=[dict(role=x[0], label=x[1], fp=list(x[2]) if x[2] else None,
                                             qsh=[str(y) for y in x[3]] if x[3] else None,
                                             osc=x[4]) for x in tens],
                            file=os.path.basename(path)))
                if bar is not None:
                    sb = "".join("1" if x in ros else "0" for x in TENB)
                    if sb != "000":
                        bar[sb] += 1
                # right-mover / oscillator homogeneity among the 10-roles at this cell
                if collect:
                    tens = [x for x in recs if x[0] in TEN]
                    if len(tens) >= 2:
                        if len(set(str(x[3]) for x in tens)) == 1:
                            acc["qsh_same"] += 1
                        else:
                            acc["qsh_diff"] += 1
                            if len(acc["qsh_diff_examples"]) < 8:
                                acc["qsh_diff_examples"].append(dict(
                                    model=lab, key=str(key), pat=sig,
                                    q=[(x[0], str(x[3]), x[4]) for x in tens]))
                        if len(set((x[4] or "no") for x in tens)) == 1:
                            acc["osc_same"] += 1
                        else:
                            acc["osc_diff"] += 1
            return n

        acc["cellsF"] += tally(cF, acc["patF"], None, bar=acc["patbarF"],
                               sector=True, collect=True)
        acc["cellsS"] += tally(cS, acc["patS"], None)
        if want_orbit:
            acc["cellsO"] += tally(cO, acc["patO"], None)
            acc["n_fp_per_kl_collapsed"] += max(0, len(cF) - len(cO))
    for k in ("patF", "patO", "patS", "patbarF", "by_sector"):
        acc[k] = dict(acc[k])
    return acc


def main():
    geom = sys.argv[1] if len(sys.argv) > 1 else "Z2xZ4_1_1"
    nfiles = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    want_orbit = os.environ.get("D6_ORBIT", "1") == "1"
    files = sorted(glob.glob(os.path.join(GEOM[geom]["census"], "*.txt.gz")))
    if nfiles and nfiles < len(files):
        random.Random(20260820).shuffle(files)
        files = sorted(files[:nfiles])
    with mp.Pool(6) as pool:
        accs = pool.map(worker, [(f, geom, want_orbit) for f in files])
    out = dict(geometry=geom, n_files=len(files), files=[os.path.basename(f) for f in files],
               models=0, models_skipped=0,
               cellsF=0, cellsO=0, cellsS=0, qsh_same=0, qsh_diff=0,
               osc_same=0, osc_diff=0, n_fp_per_kl_collapsed=0,
               patF=Counter(), patO=Counter(), patS=Counter(), patbarF=Counter(),
               by_sector=Counter(), hits=[], qsh_diff_examples=[])
    for a in accs:
        for k in ("models", "models_skipped", "cellsF", "cellsO", "cellsS",
                  "qsh_same", "qsh_diff", "osc_same", "osc_diff",
                  "n_fp_per_kl_collapsed"):
            out[k] += a[k]
        for k in ("patF", "patO", "patS", "patbarF", "by_sector"):
            out[k].update(a[k])
        out["hits"] += a["hits"]
        out["qsh_diff_examples"] += a["qsh_diff_examples"][:4]
    for k in ("patF", "patO", "patS", "patbarF", "by_sector"):
        out[k] = dict(out[k])
    out["TABLE"] = {p: out["patF"].get(p, 0) for p in PATS}
    out["TABLE_orbitkey"] = {p: out["patO"].get(p, 0) for p in PATS}
    out["TABLE_sectoronlykey"] = {p: out["patS"].get(p, 0) for p in PATS}
    out["TABLE_conj"] = {p: out["patbarF"].get(p, 0) for p in PATS}
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "d6_patterns_%s.json" % geom), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print(json.dumps({k: v for k, v in out.items()
                      if k not in ("by_sector", "hits", "files", "patF", "patO",
                                   "patS", "patbarF", "qsh_diff_examples")},
                     indent=1, default=str))
    print("n_hits(010/001/110/101) =", len(out["hits"]))


if __name__ == "__main__":
    main()

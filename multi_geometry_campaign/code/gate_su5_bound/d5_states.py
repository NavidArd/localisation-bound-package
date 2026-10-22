"""GATE gate_su5_bound D5 -- (a) RECORDS vs STATES, (b) the 599 recoverable excluded models.

TWO ATTACKS.

(a) THE COUNTER COUNTS RECORDS, NOT STATES.  verify.py does `rmult[ro] += 1` per spectrum
    RECORD; the certified gate (g1d_crosstab.py) does `nf[(kl, fp, role)] += 1` likewise.
    A record whose representation is a NON-SINGLET of a hidden factor is several distinct
    massless left-chiral states -- orbifolder's own SMConfig table counts them with that
    multiplicity (`printed[...] += row[3]`).  Those states share the SM role and the cell
    and differ by a hidden root, so they are exactly the configuration the theorem forbids:
    Delta = p - p' in Lambda, Delta orthogonal to W, |Delta|^2 >= 2  =>  R >= 1/2, while the
    measured R for 10-roles on Z2xZ4_1_1 is at most 17/40 < 1/2.  If ANY twisted 10-role
    record carries hidden multiplicity >= 2, the theorem is FALSIFIED and the census cannot
    see it.  This script measures the hidden-multiplicity distribution directly.

(b) THE EXCLUDED MODELS.  d5_excluded.py showed 599 of the 615 Z2xZ4 "closure_failed"
    models have a LATER printed vev-configuration that DOES close; the census takes cfgs[0]
    only.  Those models are therefore censusable and were never censused.  This script
    censuses them with the first closing configuration and repeats the per-cell 10-role
    multiplicity measurement on them.

Outputs data/gate_su5_bound/d5_states.json
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"

import glob                                                            # noqa: E402
import gzip                                                            # noqa: E402
import json                                                            # noqa: E402
import multiprocessing as mp                                           # noqa: E402
import sys                                                             # noqa: E402
from collections import Counter, defaultdict                           # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CAMP = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(CAMP, "data", "gate_su5_bound")

ROLES10 = {"Q", "uc", "ec", "Qbar", "ucbar", "ecbar"}
ROLES5 = {"dc", "L", "Hu", "dcbar"}

GEOM = {
    "Z2xZ4_1_1": dict(lib="z2xz4_sweep", mod="z24sw_lib", M=2, N=4,
                      census=os.path.join(CAMP, "data", "z2xz4_sweep", "raw", "census")),
    "Z3xZ6_1_1": dict(lib="z3xz6_sweep", mod="z36sw_lib", M=3, N=6,
                      census=os.path.join(CAMP, "data", "z3xz6_sweep", "raw", "census")),
}


def _roles_for(Z, geom, blk, pre, cfg_index):
    """(roles, spec) using the cfg_index-th printed vev-configuration, or None."""
    if geom == "Z2xZ4_1_1":
        spec = Z.parse_spectrum(pre, M=2, N=4)
        factors = Z.parse_gauge_factors(blk)
        cfgs = Z.parse_smconfig(blk)
        if not cfgs or not spec or cfg_index >= len(cfgs):
            return None
        tY, tbl = cfgs[cfg_index]
        if tY is None:
            return None
        printed = defaultdict(int)
        for row in tbl:
            printed[(row[0], row[1], row[2])] += row[3]
        pos = Z.observable_positions(spec, tY, factors, dict(printed))
        if pos is None:
            return None
        roles, _ = Z.role_census(spec, tY, *pos)
        return roles, spec
    spec = Z.parse_spectrum(pre)
    gauge = Z.parse_gauge(blk)
    if gauge["tY"] is None or not spec or cfg_index != 0:
        return None
    pos = Z.observable_positions(spec, gauge["tY"], gauge)
    if pos is None:
        return None
    a3, a2, sY0, tY, bi = pos
    roles, _ = Z.role_census(spec, tY, a3, a2, sY0)
    return roles, spec


def _tally(acc, tag, roles, spec):
    cells_rec = defaultdict(Counter)
    cells_state = defaultdict(Counter)
    for lab, r in spec.items():
        if r["p_sh"] is None or r["chir"] != "l":
            continue
        rr = roles.get(lab)
        if not rr or rr.get("role") is None:
            continue
        kl = tuple(r["kl"])
        if kl == (0, 0):
            continue
        key = (kl, tuple(r["fp"]) if r["fp"] else None)
        ro = rr["role"]
        fam = "10" if ro in ROLES10 else ("5bar" if ro in ROLES5 else "?")
        mult = rr.get("mult", 1)
        cells_rec[key][ro] += 1
        cells_state[key][ro] += mult
        acc[tag + "hidmult"][fam][str(mult)] += 1
        if fam == "10" and mult >= 2 and len(acc[tag + "hidmult_ex"]) < 40:
            acc[tag + "hidmult_ex"].append(
                dict(role=ro, kl=list(kl), fp=list(r["fp"]) if r["fp"] else None,
                     mult=mult, rep=list(r["rep"]), label=lab))
    for key, c in cells_rec.items():
        for ro, m in c.items():
            fam = "10" if ro in ROLES10 else ("5bar" if ro in ROLES5 else "?")
            acc[tag + "mult_rec"][fam][str(m)] += 1
    for key, c in cells_state.items():
        for ro, m in c.items():
            fam = "10" if ro in ROLES10 else ("5bar" if ro in ROLES5 else "?")
            acc[tag + "mult_state"][fam][str(m)] += 1
    # distinct cells hosting at least one member of each family
    for key, c in cells_rec.items():
        fams = {"10" if ro in ROLES10 else ("5bar" if ro in ROLES5 else "?") for ro in c}
        for f in fams:
            acc[tag + "distinct_cells"][f] += 1


def _blank():
    return dict(models_main=0, models_recovered=0, recover_cfg=Counter(),
                main_mult_rec=defaultdict(Counter), main_mult_state=defaultdict(Counter),
                main_hidmult=defaultdict(Counter), main_hidmult_ex=[],
                main_distinct_cells=Counter(),
                rec_mult_rec=defaultdict(Counter), rec_mult_state=defaultdict(Counter),
                rec_hidmult=defaultdict(Counter), rec_hidmult_ex=[],
                rec_distinct_cells=Counter(),
                unrecoverable=0)


def worker(args):
    path, geom = args
    G = GEOM[geom]
    sys.path.insert(0, os.path.join(CAMP, "code", G["lib"]))
    Z = __import__(G["mod"])
    acc = _blank()
    txt = gzip.open(path, "rt", errors="replace").read()
    for lab, blk in Z.split_models(txt):
        pre = blk.split("Analyzing")[0]
        try:
            got = _roles_for(Z, geom, blk, pre, 0)
        except Exception:                                              # noqa: BLE001
            got = None
        if got is not None:
            acc["models_main"] += 1
            _tally(acc, "main_", *got)
            continue
        # ---- excluded by the census: try later configurations
        rec = None
        if geom == "Z2xZ4_1_1":
            try:
                n = len(Z.parse_smconfig(blk))
            except Exception:                                          # noqa: BLE001
                n = 0
            for i in range(1, n):
                try:
                    rec = _roles_for(Z, geom, blk, pre, i)
                except Exception:                                      # noqa: BLE001
                    rec = None
                if rec is not None:
                    acc["recover_cfg"][str(i)] += 1
                    break
        if rec is None:
            acc["unrecoverable"] += 1
        else:
            acc["models_recovered"] += 1
            _tally(acc, "rec_", *rec)
    for k in ("main_mult_rec", "main_mult_state", "main_hidmult",
              "rec_mult_rec", "rec_mult_state", "rec_hidmult"):
        acc[k] = {a: dict(b) for a, b in acc[k].items()}
    for k in ("recover_cfg", "main_distinct_cells", "rec_distinct_cells"):
        acc[k] = dict(acc[k])
    return acc


def main():
    res = {}
    for geom in ("Z2xZ4_1_1", "Z3xZ6_1_1"):
        files = sorted(glob.glob(os.path.join(GEOM[geom]["census"], "*.txt.gz")))
        with mp.Pool(6) as pool:
            accs = pool.map(worker, [(f, geom) for f in files])
        out = _blank()
        for a in accs:
            out["models_main"] += a["models_main"]
            out["models_recovered"] += a["models_recovered"]
            out["unrecoverable"] += a["unrecoverable"]
            out["recover_cfg"].update(a["recover_cfg"])
            for k in ("main_mult_rec", "main_mult_state", "main_hidmult",
                      "rec_mult_rec", "rec_mult_state", "rec_hidmult"):
                for f, d in a[k].items():
                    out[k][f].update(d)
            for k in ("main_distinct_cells", "rec_distinct_cells"):
                out[k].update(a[k])
            out["main_hidmult_ex"] += a["main_hidmult_ex"][:5]
            out["rec_hidmult_ex"] += a["rec_hidmult_ex"][:5]
        for k in ("main_mult_rec", "main_mult_state", "main_hidmult",
                  "rec_mult_rec", "rec_mult_state", "rec_hidmult"):
            out[k] = {a_: dict(b) for a_, b in out[k].items()}
        for k in ("recover_cfg", "main_distinct_cells", "rec_distinct_cells"):
            out[k] = dict(out[k])
        out["main_hidmult_ex"] = out["main_hidmult_ex"][:25]
        out["rec_hidmult_ex"] = out["rec_hidmult_ex"][:25]
        res[geom] = out
        print("=====", geom)
        print(json.dumps(out, indent=1, default=str)[:4000])
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "d5_states.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=str)


if __name__ == "__main__":
    main()

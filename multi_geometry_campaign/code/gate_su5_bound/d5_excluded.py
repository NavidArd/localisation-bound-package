"""GATE gate_su5_bound D5 -- the EXCLUDED MODELS (attack #3) and the model accounting.

verify.py reports `closure_failed: 615` on Z2xZ4_1_1 and `2` on Z3xZ6_1_1 and censuses the
rest.  Two questions this script answers OFF THE RAW TEXT, with its own parsing:

  (1) COMPLETENESS.  Is (censused + skipped) equal to the number of model blocks actually
      present in the raw census, and to orbifolder's own "N orbifolds loaded" declaration?
      Anything missing is a SILENT drop.
  (2) SAFETY.  Can the excluded models be given a hypercharge normalisation k_Y = 2|t_Y|^2
      at all?  If any excluded model carries k_Y > 5/3 the e^c threshold moves toward the
      permissive side and the exclusion is not a-fortiori safe.  We read EVERY printed
      `Hypercharge generator:` of every excluded model and compute k_Y exactly over Q.

  (3) WHY they fail.  verify.py (and the certified gate g1d_crosstab.py before it) takes
      cfgs[0] -- the FIRST printed vev-configuration -- and declares closure failure if
      that one does not rebuild.  We count how many configurations each excluded model
      prints, so the reader can see whether "closure failed" means "no SM config" or
      "the SM config was not the first one printed".

Outputs data/gate_su5_bound/d5_excluded.json
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"

import glob                                                            # noqa: E402
import gzip                                                            # noqa: E402
import json                                                            # noqa: E402
import multiprocessing as mp                                           # noqa: E402
import re                                                              # noqa: E402
import sys                                                             # noqa: E402
from collections import Counter                                        # noqa: E402
from fractions import Fraction as F                                    # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CAMP = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(CAMP, "data", "gate_su5_bound")

GEOM = {
    "Z2xZ4_1_1": dict(lib="z2xz4_sweep", mod="z24sw_lib", M=2, N=4,
                      census=os.path.join(CAMP, "data", "z2xz4_sweep", "raw", "census")),
    "Z3xZ6_1_1": dict(lib="z3xz6_sweep", mod="z36sw_lib", M=3, N=6,
                      census=os.path.join(CAMP, "data", "z3xz6_sweep", "raw", "census")),
}

# ---- own parsing (deliberately not imported from the census libs) ----------------
_BLK = re.compile(r"/> /(MSSM\d+)> /")
_LOADED = re.compile(r"(\d+) orbifolds loaded")
_HYP = re.compile(r"Hypercharge generator:\s*\n(.*)")
_ROW = re.compile(
    r"(?m)^\s*(\d+)\s+\(\s*(-?\d+),\s*(-?\d+)\)_([lr])\s+U\(1\)\s*:\s*\(\s*([^)]*)\)")


def own_split(text):
    idx = [(m.start(), m.group(1)) for m in _BLK.finditer(text)]
    return [(lab, text[p:(idx[i + 1][0] if i + 1 < len(idx) else len(text))])
            for i, (p, lab) in enumerate(idx)]


def own_vec(s):
    out = []
    for grp in re.findall(r"\(([^)]*)\)", s):
        out += [F(x.strip()) if x.strip() else F(0) for x in grp.split(",")]
    return out


def own_generators(blk):
    """every printed hypercharge generator, in print order, with its k_Y = 2|t_Y|^2."""
    cuts = [m.start() for m in re.finditer(r"Hypercharge generator:", blk)]
    res = []
    for i, p in enumerate(cuts):
        seg = blk[p:(cuts[i + 1] if i + 1 < len(cuts) else len(blk))]
        m = _HYP.search(seg)
        if not m:
            res.append((None, None, 0))
            continue
        tY = own_vec(m.group(1))
        kY = 2 * sum(x * x for x in tY)
        res.append((tY, kY, len(_ROW.findall(seg))))
    return res


def worker(args):
    path, geom = args
    G = GEOM[geom]
    sys.path.insert(0, os.path.join(CAMP, "code", G["lib"]))
    Z = __import__(G["mod"])
    txt = gzip.open(path, "rt", errors="replace").read()
    acc = dict(loaded=sum(int(x) for x in _LOADED.findall(txt)),
               own_blocks=0, lib_blocks=0,
               ok=0, skip=Counter(), labels=Counter(),
               excl_kY=Counter(), excl_ncfg=Counter(),
               excl_kY_by_ncfg=Counter(),
               excl_examples=[], ok_kY=Counter(),
               excl_any_cfg_closes=0, excl_cfg0_no_generator=0)
    own = own_split(txt)
    acc["own_blocks"] = len(own)
    acc["lib_blocks"] = len(Z.split_models(txt))
    for lab, blk in own:
        acc["labels"][lab] += 1
        pre = blk.split("Analyzing")[0]
        gens = own_generators(blk)
        try:
            if geom == "Z2xZ4_1_1":
                spec = Z.parse_spectrum(pre, M=G["M"], N=G["N"])
                factors = Z.parse_gauge_factors(blk)
                cfgs = Z.parse_smconfig(blk)
                if not cfgs or cfgs[0][0] is None or not spec:
                    acc["skip"]["no_smconfig"] += 1
                    continue
                from collections import defaultdict as _dd
                closes = []
                for (tY, tbl) in cfgs:
                    if tY is None:
                        closes.append(False)
                        continue
                    printed = _dd(int)
                    for row in tbl:
                        printed[(row[0], row[1], row[2])] += row[3]
                    closes.append(
                        Z.observable_positions(spec, tY, factors, dict(printed)) is not None)
                good = closes[0]
            else:
                spec = Z.parse_spectrum(pre)
                gauge = Z.parse_gauge(blk)
                if gauge["tY"] is None or not spec:
                    acc["skip"]["no_hypercharge"] += 1
                    continue
                closes = [Z.observable_positions(spec, gauge["tY"], gauge) is not None]
                good = closes[0]
        except Exception as e:                                          # noqa: BLE001
            acc["skip"]["parse_error:%s" % type(e).__name__] += 1
            continue
        kYs = sorted({str(k) for _t, k, _n in gens if k is not None})
        if good:
            acc["ok"] += 1
            for k in kYs:
                acc["ok_kY"][k] += 1
        else:
            acc["skip"]["closure_failed"] += 1
            acc["excl_ncfg"][len(gens)] += 1
            for k in kYs:
                acc["excl_kY"][k] += 1
            acc["excl_kY_by_ncfg"]["ncfg=%d|%s" % (len(gens), ",".join(kYs))] += 1
            if any(closes[1:]):
                acc["excl_any_cfg_closes"] += 1
            if gens and gens[0][0] is None:
                acc["excl_cfg0_no_generator"] += 1
            if len(acc["excl_examples"]) < 8:
                acc["excl_examples"].append(
                    dict(model=lab, file=os.path.basename(path), n_cfg=len(gens),
                         kY=[str(k) for _t, k, _n in gens],
                         rows=[n for _t, _k, n in gens],
                         later_cfg_closes=bool(any(closes[1:]))))
    for k in ("skip", "labels", "excl_kY", "excl_ncfg", "excl_kY_by_ncfg", "ok_kY"):
        acc[k] = dict(acc[k])
    return acc


def main():
    res = {}
    for geom in ("Z2xZ4_1_1", "Z3xZ6_1_1"):
        files = sorted(glob.glob(os.path.join(GEOM[geom]["census"], "*.txt.gz")))
        with mp.Pool(6) as pool:
            accs = pool.map(worker, [(f, geom) for f in files])
        tot = dict(loaded=0, own_blocks=0, lib_blocks=0, ok=0,
                   excl_any_cfg_closes=0, excl_cfg0_no_generator=0)
        cs = {k: Counter() for k in ("skip", "labels", "excl_kY", "excl_ncfg",
                                     "excl_kY_by_ncfg", "ok_kY")}
        ex = []
        for a in accs:
            for k in tot:
                tot[k] += a[k]
            for k in cs:
                cs[k].update(a[k])
            ex += a["excl_examples"][:2]
        dup = {k: v for k, v in cs["labels"].items() if v > 1}
        res[geom] = dict(tot, skip={k: v for k, v in cs["skip"].items()},
                         excl_kY=dict(cs["excl_kY"]), excl_ncfg=dict(cs["excl_ncfg"]),
                         excl_kY_by_ncfg=dict(cs["excl_kY_by_ncfg"]),
                         ok_kY=dict(cs["ok_kY"]),
                         n_distinct_labels=len(cs["labels"]),
                         duplicate_labels=len(dup),
                         duplicate_label_examples=sorted(dup)[:12],
                         excl_examples=ex[:20])
        res[geom]["accounting_ok"] = (
            tot["own_blocks"] == tot["lib_blocks"] == tot["loaded"]
            == tot["ok"] + sum(cs["skip"].values()))
        print(geom, json.dumps({k: v for k, v in res[geom].items()
                                if k != "excl_examples"}, indent=1)[:3000])
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "d5_excluded.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=str)


if __name__ == "__main__":
    main()

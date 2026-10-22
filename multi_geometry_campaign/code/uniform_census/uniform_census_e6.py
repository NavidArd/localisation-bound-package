r"""Paper A v6 --- the uniform-convention, state-level census WITH the `E_6` parse repaired.

Why this script exists
----------------------
Two independent corrections to Paper A's Table 1 are outstanding and they overlap:

* `paper15/review/paperA_uniform_census.md` (2026-09-05) replaces the table's two role
  conventions and its as-run model count by one uniform convention (ungated roles, first
  CLOSING embedding, each catalogue model counted once).  Its engine is
  `paper15/review/code_referee_round3/uniform_census.py`.
* `paper15/review/exposure_gamma_e6.md` (2026-09-07) repairs the `E_6` hole in the two
  certified banner parsers, which silently drops four models whose `E_6` factor is printed
  BEFORE the Standard-Model factors (`MSSM967`, `MSSM1409` on `Z2xZ4_1_1`; `MSSM825`,
  `MSSM2669` on `Z3xZ6_1_1`).

The uniform census predates the `E_6` finding and therefore does NOT contain those four
models (it reports 11,524 of 11,540 and 4,610 of 4,612).  Neither memo gives the combined
figures, so they are computed here rather than inferred: `uniform_census.worker` is run
unmodified, with the two defective parsers rebound AT RUNTIME ONLY.  No file under
`multi_geometry_campaign/code/` is edited (CLAUDE.md 8.6); the certified libraries are
imported verbatim and patched in the child process' memory.

The two rebound functions
-------------------------
* `z2xz4_sweep/z24sw_lib.py:889 parse_gauge_factors` --- `r"(SU|SO|E)\((\d+)\)"` never
  matches `E_6`, so every factor after it is shifted by one position.
* `z3xz6_sweep/z36sw_lib.py:392 parse_banner` --- `r"(SU|SO|SP|E)[_]?\((\d+)\)"` makes the
  underscore optional but the parentheses mandatory, so `E_6` is unmatched there too.

Both replacements are the ones `paper15/review/code_exposure/ex2_recompute.py:48` and
`ex4_z3xz6_recompute.py:41` verified: `r"(SU|SO|SP|E)(?:\((\d+)\)|_(\d+))"`, order
preserved.  No model on either catalogue prints an `SP` factor, so on every model without
an exceptional factor the corrected parse is bit-identical to the library one.

Usage
-----
    python uniform_census_e6.py Z2xZ4_1_1
    python uniform_census_e6.py Z3xZ6_1_1

Output: `paper15/review/data_paperA_v6/uniform_e6_<geom>.json`, the same accumulator
schema as `uniform_census.py` writes, plus a `patched_models` record naming the models
whose parse actually changed.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"                                                # noqa: E402

import glob                                                             # noqa: E402
import json                                                             # noqa: E402
import multiprocessing as mp                                            # noqa: E402
import re                                                               # noqa: E402
import sys                                                              # noqa: E402
import time                                                             # noqa: E402
from collections import Counter                                         # noqa: E402

# Package copy: paths derived from this file's own location.
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "..", ".."))
R3 = os.path.dirname(os.path.abspath(__file__))
CAMP = os.path.join(REPO, "multi_geometry_campaign")
OUT = os.path.join(REPO, "results_deposited", "uniform_census")

sys.path.insert(0, R3)

FIXRE = re.compile(r"(SU|SO|SP|E)(?:\((\d+)\)|_(\d+))")

_PATCHED = False


def _patch():
    """Rebind the two defective banner parsers in the child process, once."""
    global _PATCHED
    if _PATCHED:
        return
    for lib in ("z2xz4_sweep", "z3xz6_sweep"):
        p = os.path.join(CAMP, "code", lib)
        if p not in sys.path:
            sys.path.insert(0, p)

    import z24sw_lib                                                    # noqa: E402
    import z36sw_lib                                                    # noqa: E402

    def corrected_pgf(block):
        """`E_n`-aware, position-preserving `z24sw_lib.parse_gauge_factors`."""
        m = re.search(r'Gauge group in vev-configuration "TestConfig1":\s*(.*)', block)
        if not m:
            return []
        s = m.group(1).split("and U(1)")[0]
        return [(t, int(a or b)) for t, a, b in FIXRE.findall(s)]

    def corrected_banner(b):
        """`E_n`-aware, position-preserving `z36sw_lib.parse_banner`."""
        out = []
        if not b:
            return out
        head = b.split(" and U(1)")[0]
        for t, a, c in FIXRE.findall(head):
            n = a or c
            if t == "SU" and n == "1":
                continue
            out.append((t, int(n)))
        return out

    z24sw_lib.parse_gauge_factors = corrected_pgf
    z36sw_lib.parse_banner = corrected_banner
    _PATCHED = True


def worker_e6(args):
    _patch()
    import uniform_census as uc
    return uc.worker(args)


def scan_e6(args):
    """Which blocks of this segment carry an exceptional factor, and where."""
    import gzip
    path, geom = args
    _patch()
    import uniform_census as uc
    sys.path.insert(0, os.path.join(CAMP, "code", uc.GEOM[geom]["lib"]))
    Z = __import__(uc.GEOM[geom]["mod"])
    txt = gzip.open(path, "rt", errors="replace").read()
    hits = []
    for lab, blk in Z.split_models(txt):
        if geom == "Z2xZ4_1_1":
            m = re.search(r'Gauge group in vev-configuration "TestConfig1":\s*(.*)', blk)
            head = m.group(1).split("and U(1)")[0] if m else ""
        else:
            g = Z.parse_gauge(blk)
            head = (g.get("gauge_test") or "").split(" and U(1)")[0]
        fixed = [(t, int(a or c)) for t, a, c in FIXRE.findall(head)]
        if geom == "Z3xZ6_1_1":
            fixed = [x for x in fixed if not (x[0] == "SU" and x[1] == 1)]
        exc = [i for i, (t, _n) in enumerate(fixed) if t in ("E",)]
        if exc:
            hits.append(dict(model=lab, position=exc[0], factors=len(fixed)))
    return hits


def main():
    geom = sys.argv[1]
    import uniform_census as uc
    files = sorted(glob.glob(os.path.join(uc.GEOM[geom]["census"], "*.txt.gz")))
    if len(sys.argv) > 2:
        files = files[:int(sys.argv[2])]
    t0 = time.time()

    with mp.Pool(6) as pool:
        labmap = dict(pool.map(uc.labels_of, [(f, geom) for f in files]))
    seen, skipmap, dups = set(), {}, Counter()
    for f in files:
        sk = []
        for lab in labmap[f]:
            if lab in seen:
                sk.append(lab)
                dups[lab] += 1
            else:
                seen.add(lab)
        skipmap[f] = sk
    n_blocks = sum(len(v) for v in labmap.values())
    print("pass1: %d blocks, %d distinct labels, %d duplicate occurrences in %d labels"
          % (n_blocks, len(seen), sum(dups.values()), len(dups)), flush=True)

    with mp.Pool(6) as pool:
        exc = pool.map(scan_e6, [(f, geom) for f in files])
    excmodels = [h for hs in exc for h in hs]
    print("exceptional-factor models: %d (position 0: %d)"
          % (len(excmodels), sum(1 for h in excmodels if h["position"] == 0)), flush=True)

    with mp.Pool(6) as pool:
        res = pool.map(worker_e6, [(f, geom, skipmap[f]) for f in files])
    C, A, U = uc.newacc(), uc.newacc(), uc.newacc()
    for c, a, u in res:
        uc.merge_into(C, c)
        uc.merge_into(A, a)
        uc.merge_into(U, u)
    if geom == "Z2xZ4_1_1":
        C = A
    out = dict(geometry=geom, n_segments=len(files),
               seconds=round(time.time() - t0, 1),
               parse="E6-corrected (runtime rebind; no library file edited)",
               blocks=n_blocks, distinct_labels=len(seen),
               duplicate_occurrences=sum(dups.values()),
               duplicate_labels=sorted(dups)[:40],
               exceptional_models=excmodels,
               current=uc.jsonable(C), asrun=uc.jsonable(A), uniform=uc.jsonable(U))
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "uniform_e6_%s.json" % geom), "w") as fh:
        json.dump(out, fh, indent=1)
    for nm, acc in (("CURRENT", C), ("ASRUN", A), ("UNIFORM", U)):
        print("--- %s %s: models %d (recovered %d) fields10 rec %d st %d | "
              "pairs10 rec %d | fields5 rec %d st %d | pairs5 rec %d | viol %d/%d"
              % (nm, geom, acc["models"], acc["models_recovered"],
                 acc["f_rec_by_fam"]["10"], acc["f_st_by_fam"]["10"],
                 acc["pairs_rec_by_fam"]["10"],
                 acc["f_rec_by_fam"]["5bar"], acc["f_st_by_fam"]["5bar"],
                 acc["pairs_rec_by_fam"]["5bar"], acc["viol_rec"], acc["viol_st"]),
              flush=True)
        print("    same-role record pairs total %d ; w2 %s ; fp_len %s ; I1_bad %d ; "
              "notperp %d ; notlambda %d ; d2 %s ; br_disagree %d"
              % (acc["pairs_rec_by_fam"]["10"] + acc["pairs_rec_by_fam"]["5bar"],
                 {k: dict(v) for k, v in acc["w2_by_fam"].items()},
                 dict(acc["fp_len"]), acc["I1_bad"], acc["d_notperp"],
                 acc["d_notlambda"], dict(acc["d_d2"]), acc["br_disagree"]),
              flush=True)


if __name__ == "__main__":
    main()

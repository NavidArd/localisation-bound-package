r"""Paper A v6 --- the four-dimensional-root test, re-run on the UNIFORM pair population.

Why
---
Section 8's negative half prints the number of same-role pairs of distinct records at which
the difference fails to be a root of the unbroken four-dimensional gauge group.  That test
lives in `multi_geometry_campaign/code/root_lift/rl_b_data.py` and was run under the
convention `paperA_v5.tex` prints: first PRINTED embedding on `Z2xZ4_1_1`, duplicate
`Z3xZ6_1_1` labels counted twice, defective `E_6` parse.  Its population is 24,406 + 7,888
= 32,294 pairs.

Table 1 of v6 counts the same-role pairs under the uniform convention (first CLOSING
embedding, each model once, `E_6` repaired).  That population is larger, so leaving the
Section 8 sentence at its old population would either print two pair counts under two
conventions --- which is exactly what the referee charged --- or assert the test over pairs
it never saw.  This script re-runs the test over the uniform population instead.

How
---
`rl_b_data.worker` is imported and run unmodified.  Three things are rebound at runtime, in
the child process only; no file under `multi_geometry_campaign/code/` is edited:

1. the two defective banner regexes (`z24sw_lib.parse_gauge_factors`,
   `z36sw_lib.parse_banner`), as in `uniform_census_e6.py`;
2. `z24sw_lib.parse_smconfig`, rotated so that the FIRST CLOSING configuration is
   `cfgs[0]` --- `rl_b_data` takes `cfgs[0]` outright, so this is how the first-closing
   rule reaches it;
3. `split_models`, filtered against the duplicate-label map measured in pass 1, so that
   the ten `Z3xZ6_1_1` models the census output carries twice are counted once.

What is NOT claimed here
------------------------
`rl_b_data`'s own `notperp` counter is a dead branch (`paperA_uniform_census.md` F-1): it
builds its Standard-Model basis from integer factor positions inside a bare `try/except`,
so the test never runs and the counter can never increment.  Nothing in this script repairs
it and no orthogonality number is taken from it.  The orthogonality of every uniform pair
is measured independently, against a real basis, by `uniform_census_e6.py`.

Usage:  python rl_b_uniform.py Z2xZ4_1_1
        python rl_b_uniform.py Z3xZ6_1_1
Output: `paper15/review/data_paperA_v6/rootlift_uniform_<geom>.json`
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"                                                # noqa: E402

import glob                                                             # noqa: E402
import gzip                                                             # noqa: E402
import json                                                             # noqa: E402
import multiprocessing as mp                                            # noqa: E402
import re                                                               # noqa: E402
import sys                                                              # noqa: E402
import time                                                             # noqa: E402
from collections import Counter, defaultdict                            # noqa: E402

# Package copy: paths derived from this file's own location.
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "..", ".."))
CAMP = os.path.join(REPO, "multi_geometry_campaign")
RL = os.path.join(CAMP, "code", "root_lift")
OUT = os.path.join(REPO, "results_deposited", "uniform_census")

sys.path.insert(0, RL)

FIXRE = re.compile(r"(SU|SO|SP|E)(?:\((\d+)\)|_(\d+))")

_SKIP: set[str] = set()
_DONE = False


def _corrected_pgf(block):
    m = re.search(r'Gauge group in vev-configuration "TestConfig1":\s*(.*)', block)
    if not m:
        return []
    s = m.group(1).split("and U(1)")[0]
    return [(t, int(a or b)) for t, a, b in FIXRE.findall(s)]


def _corrected_banner(b):
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


def _patch(geom, skip):
    """Rebind the three things, in this process, once."""
    global _DONE, _SKIP
    _SKIP = set(skip)
    if _DONE:
        return
    for lib in ("z2xz4_sweep", "z3xz6_sweep"):
        p = os.path.join(CAMP, "code", lib)
        if p not in sys.path:
            sys.path.insert(0, p)
    import z24sw_lib                                                    # noqa: E402
    import z36sw_lib                                                    # noqa: E402

    z24sw_lib.parse_gauge_factors = _corrected_pgf
    z36sw_lib.parse_banner = _corrected_banner

    orig_smcfg = z24sw_lib.parse_smconfig

    def first_closing_smconfig(block):
        """Rotate the printed configurations so the first CLOSING one is `cfgs[0]`."""
        cfgs = orig_smcfg(block)
        if not cfgs:
            return cfgs
        pre = block.split("Analyzing")[0]
        try:
            spec = z24sw_lib.parse_spectrum(pre, M=2, N=4)
            factors = _corrected_pgf(block)
        except Exception:                                               # noqa: BLE001
            return cfgs
        if not spec:
            return cfgs
        for i, (tY, tbl) in enumerate(cfgs):
            if tY is None:
                continue
            printed = defaultdict(int)
            for row in tbl:
                printed[(row[0], row[1], row[2])] += row[3]
            try:
                if z24sw_lib.observable_positions(spec, tY, factors,
                                                  dict(printed)) is not None:
                    return cfgs[i:] + cfgs[:i]
            except Exception:                                           # noqa: BLE001
                continue
        return cfgs

    z24sw_lib.parse_smconfig = first_closing_smconfig

    for mod in (z24sw_lib, z36sw_lib):
        orig_split = mod.split_models

        def filtered(txt, _o=orig_split):
            return [(lab, blk) for lab, blk in _o(txt) if lab not in _SKIP]

        mod.split_models = filtered

    _DONE = True


def labels_of(args):
    path, geom = args
    import rl_b_data as RB
    sys.path.insert(0, os.path.join(CAMP, "code", RB.GEOM[geom]["lib"]))
    Z = __import__(RB.GEOM[geom]["mod"])
    txt = gzip.open(path, "rt", errors="replace").read()
    return path, [lab for lab, _b in Z.split_models(txt)]


def worker_u(args):
    path, geom, skip = args
    _patch(geom, skip)
    import rl_b_data as RB
    return RB.worker((path, geom))


def main():
    geom = sys.argv[1]
    import rl_b_data as RB
    files = sorted(glob.glob(os.path.join(RB.GEOM[geom]["census"], "*.txt.gz")))
    if len(sys.argv) > 2:
        files = files[:int(sys.argv[2])]
    t0 = time.time()

    with mp.Pool(6) as pool:
        labmap = dict(pool.map(labels_of, [(f, geom) for f in files]))
    seen, skipmap = set(), {}
    for f in files:
        sk = []
        for lab in labmap[f]:
            if lab in seen:
                sk.append(lab)
            else:
                seen.add(lab)
        skipmap[f] = sk
    print("pass1: %d blocks, %d distinct labels" %
          (sum(len(v) for v in labmap.values()), len(seen)), flush=True)

    with mp.Pool(6) as pool:
        accs = pool.map(worker_u, [(f, geom, skipmap[f]) for f in files])
    res = RB.merge(accs)

    tags = res["by"].get("pairs", {})
    root_true = sum(v for k, v in tags.items() if k.endswith("4droot=True")
                    or k.endswith("4droot=I:True"))
    root_false = sum(v for k, v in tags.items() if k.endswith("4droot=False")
                     or k.endswith("4droot=I:False"))
    res["uniform_convention"] = ("first CLOSING embedding; duplicate labels dropped; "
                                 "E_6 parse repaired; runtime rebind only")
    res["pairs_4droot_true"] = root_true
    res["pairs_4droot_false"] = root_false
    res["seconds"] = round(time.time() - t0, 1)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "rootlift_uniform_%s.json" % geom), "w") as fh:
        json.dump(res, fh, indent=1, default=str)
    print("models=%d  same-role pairs=%d  Delta not in Lambda=%d  d2=%s"
          % (res["models"], res["pairs"], res["notlambda"], res["d2"]), flush=True)
    print("4d-root TRUE=%d  FALSE=%d  (sum %d)  route disagreements=%d"
          % (root_true, root_false, root_true + root_false,
             res["route_disagree"]), flush=True)
    print("NOTE: res['notperp']=%s is rl_b_data's dead counter and is not a measurement."
          % res["notperp"], flush=True)


if __name__ == "__main__":
    main()

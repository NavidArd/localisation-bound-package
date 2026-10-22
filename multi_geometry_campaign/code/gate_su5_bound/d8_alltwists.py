"""GATE d8b -- extend the budget scan to EVERY point group shipped with orbifolder 1.2.1.

Tests the findings' 7.1 prediction ("Z6 is the first order at which the bound can fail,
and it fails only for the smallest rotation at that order") against the complete list of
6d abelian toroidal orbifold point groups, read off the primary geometry files.

For each distinct point group (M,N,twists) we enumerate the twisted sectors and report
every sector with B = 1 - delta_c >= 17/20 (i.e. delta_c <= 3/20), classified as
N=2-type (one unrotated plane) or fully rotated.

Writes data/gate_su5_bound/d8_alltwists.json.
"""

from __future__ import annotations

# --- reproduction-package root shim -------------------------------------------------
# The line above originally read  F:\physics  -- the absolute path of the author's
# working tree.  This is the ONLY edit made to this file in the reproduction package;
# it points the same variable at the package root, which carries the same directory
# names.  Set PAPERA_PKG to relocate the package.  See
# environment/PATCHES_TO_PACKAGE_CODE.md.
# -------------------------------------------------------------------------------------
import os as _os
_PKG_ROOT = _os.environ.get('PAPERA_PKG') or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..'))

import glob
import json
import os
import re
from fractions import Fraction as F

GEODIR = _os.path.join(_PKG_ROOT, 'orbifolder-1.2.1', 'Geometry')
OUT = _os.path.join(_PKG_ROOT, 'multi_geometry_campaign', 'data', 'gate_su5_bound', 'd8_alltwists.json')


def read_geometry(path):
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    if lines[0].strip().lower() != "point group":
        return None
    M, N = int(lines[1].strip()), int(lines[2].strip())
    i0 = lines.index("begin twist")
    i1 = lines.index("end twist")
    tw = []
    for i in range(i0 + 1, i1):
        raw = lines[i].strip()
        if not raw:
            continue
        parts = re.findall(r"-?\d+/\d+|-?\d+", raw)
        vec = [F(p) for p in parts]
        tw.append(tuple(vec[1:]))
    return M, N, tuple(tw)


def frac(x):
    return "%d/%d" % (x.numerator, x.denominator) if x.denominator != 1 else str(x.numerator)


def main():
    seen = {}
    for path in sorted(glob.glob(os.path.join(GEODIR, "Geometry_*.txt"))):
        g = read_geometry(path)
        if g is None:
            continue
        M, N, tw = g
        key = (M, N, tw)
        seen.setdefault(key, []).append(os.path.basename(path))

    out = {"_meta": {"geometry_dir": GEODIR, "n_files": sum(len(v) for v in seen.values()),
                     "n_distinct_point_groups": len(seen)}}
    rows = []
    for (M, N, tw), files in sorted(seen.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        v1 = tw[0]
        v2 = tw[1] if len(tw) > 1 else (F(0), F(0), F(0))
        fails = []
        for k in range(M):
            for l in range(N):
                if (k, l) == (0, 0):
                    continue
                eta = tuple((k * v1[a] + l * v2[a]) % 1 for a in range(3))
                dc = sum(e * (1 - e) for e in eta) / 2
                if dc <= F(3, 20):
                    fails.append(dict(k=k, l=l, eta=[frac(e) for e in eta],
                                      delta_c=frac(dc), B=frac(1 - dc),
                                      n_unrotated_planes=sum(1 for e in eta if e == 0)))
        rows.append(dict(point_group=[M, N],
                         twists=[[frac(x) for x in t] for t in tw],
                         example_file=files[0], n_files=len(files),
                         n_twisted=M * N - 1,
                         n_sectors_permitting_10_doubling=len(fails),
                         failing_sectors=fails))
    out["point_groups"] = rows

    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)

    print("distinct point groups:", len(seen), " files:", out["_meta"]["n_files"])
    print()
    print("%-10s %-40s %5s %s" % ("PG", "twists", "#tw", "sectors with B>=17/20"))
    for r in rows:
        tag = "" if r["n_sectors_permitting_10_doubling"] == 0 else "  <== LOOPHOLE"
        print("%-10s %-40s %5d %d%s" % (
            "Z%dxZ%d" % tuple(r["point_group"]),
            ";".join(",".join(t) for t in r["twists"]),
            r["n_twisted"], r["n_sectors_permitting_10_doubling"], tag))
        for f in r["failing_sectors"]:
            print("        (%d,%d) eta=%s  dc=%s B=%s  unrot=%d"
                  % (f["k"], f["l"], f["eta"], f["delta_c"], f["B"], f["n_unrotated_planes"]))


if __name__ == "__main__":
    main()

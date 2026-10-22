"""GATE gate_su5_bound D5 -- does the raw `fixed point n_a` key SPLIT a physical cell?

The decisive test.  The physical cells of a twisted sector are the SIMULTANEOUS affine
point-group orbits of the sector's 6d fixed-point set -- the "constructing elements" the
geometry file itself prints.  The certified engine
`alphabet_survey/code/phase2/p2_geom.py::GeomMN.sim_orbits` computes them.

If  #distinct printed n_a  ==  #sim_orbits  for every sector, the printed key is in
bijection with the physical cells and CANNOT split one cell into two (the direction that
would hide a same-cell doubling).  If  #printed > #sim_orbits, it splits.
If  #printed < #sim_orbits, it merges (the a-fortiori-safe direction).

d5_cellkey.py's PER-PLANE orbit test reported apparent collisions; a per-plane reading is
strictly coarser than the simultaneous one (the conjugating element acts on all three
planes at once), so that test cannot decide the question and this one replaces it.

Outputs data/gate_su5_bound/d5_keycount.json
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

import os
import sys

REPO = _PKG_ROOT
sys.path.insert(0, os.path.join(REPO, "alphabet_survey", "code"))
sys.path.insert(0, os.path.join(REPO, "alphabet_survey", "code", "phase2"))

import glob                                                            # noqa: E402
import gzip                                                            # noqa: E402
import json                                                            # noqa: E402
import re                                                              # noqa: E402
from collections import Counter, defaultdict                           # noqa: E402

from p2_geom import GeomMN                                             # noqa: E402

CAMP = os.path.join(REPO, "multi_geometry_campaign")
OUT = os.path.join(CAMP, "data", "gate_su5_bound")

GEOM = {
    "Z2xZ4_1_1": dict(file=os.path.join(REPO, "orbifolder-1.2.1", "Geometry",
                                        "Geometry_Z2xZ4_1_1.txt"), M=2, N=4,
                      census=os.path.join(CAMP, "data", "z2xz4_sweep", "raw", "census")),
    "Z3xZ6_1_1": dict(file=os.path.join(REPO, "orbifolder-1.2.1", "Geometry",
                                        "Geometry_Z3xZ6_1_1.txt"), M=3, N=6,
                      census=os.path.join(CAMP, "data", "z3xz6_sweep", "raw", "census")),
}

_PAT = re.compile(r"sector \(k,l\)\s*:\s*\(([^)]*)\)\s*\n\s*"
                  r"fixed point n_a\s*:\s*\(([^)]*)\)")


def main():
    res = {}
    for geom, G in GEOM.items():
        g = GeomMN(G["file"])
        printed = defaultdict(Counter)
        for path in sorted(glob.glob(os.path.join(G["census"], "*.txt.gz"))):
            txt = gzip.open(path, "rt", errors="replace").read()
            for m in _PAT.finditer(txt):
                kl = tuple(int(x) for x in m.group(1).split(","))
                kl = (kl[0] % G["M"], kl[1] % G["N"])
                printed[kl][tuple(int(x) for x in m.group(2).split(","))] += 1
        rows = {}
        bad = []
        for kl in sorted(set(list(printed) + list(g.sectors_all))):
            try:
                so = g.sim_orbits(kl)
                nso = len(so) if so is not None else None
            except Exception as e:                                     # noqa: BLE001
                nso = "ERR:%s" % type(e).__name__
            npr = len(printed.get(kl, {}))
            verdict = ("no sector data" if npr == 0 else
                       "BIJECTIVE" if npr == nso else
                       "SPLIT (printed > physical)" if isinstance(nso, int) and npr > nso
                       else "MERGE (printed < physical)")
            rows[str(kl)] = dict(printed_keys=npr, sim_orbits=nso, verdict=verdict,
                                 records=sum(printed.get(kl, {}).values()))
            if npr and isinstance(nso, int) and npr != nso:
                bad.append((str(kl), npr, nso))
        res[geom] = dict(sectors=rows, mismatches=bad,
                         any_split=any(p > s for _k, p, s in bad),
                         any_merge=any(p < s for _k, p, s in bad))
        print("=====", geom)
        for k, v in rows.items():
            print("  %-8s printed=%-4s sim_orbits=%-6s %s  (records %d)"
                  % (k, v["printed_keys"], v["sim_orbits"], v["verdict"], v["records"]))
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "d5_keycount.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=str)


if __name__ == "__main__":
    main()

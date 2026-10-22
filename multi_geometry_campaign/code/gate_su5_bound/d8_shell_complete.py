"""GATE d8e -- complete the loophole shell enumeration the phase truncated.

`code/su5_ten_bound/shell.py` is invoked with `run(only=only, maxcells=40)`.  There are
**65** unprotected (loophole) cells in `verify_Z3xZ6_1_1.json`, so the enumeration stopped
at the cap after 40 and 25 cells -- a census-segment-ordered TAIL, not a random sample --
were never enumerated.  results/su5_ten_bound_findings.md 7.3 nonetheless states
"at every populated instance the pre-projection shell already contains exactly one state
of that role" and 0 states the multiplicity "is established by an exact shell count".

This script re-runs the SAME enumerator with the cap lifted and writes its own output
into data/gate_su5_bound/ (the phase's file is not touched).
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

import json
import os
import sys
from collections import Counter

CAMP = _os.path.join(_PKG_ROOT, 'multi_geometry_campaign')
GOUT = os.path.join(CAMP, "data", "gate_su5_bound")
sys.path.insert(0, os.path.join(CAMP, "code", "su5_ten_bound"))

import shell                                                            # noqa: E402

# redirect the enumerator's output into the gate's own data directory
shell.OUT = GOUT


def main():
    d = json.load(open(os.path.join(CAMP, "data", "su5_ten_bound",
                                    "verify_Z3xZ6_1_1.json"), encoding="utf-8"))
    only = {c["model"] for c in d["loophole_cells"]}
    only_files = {c["file"] for c in d["loophole_cells"]}
    shell.only_files = only_files
    print("targeting %d loophole models in %d files (phase ran maxcells=40)"
          % (len(only), len(only_files)))
    shell.run(only=only, maxcells=10 ** 6)

    r = json.load(open(os.path.join(GOUT, "shell_Z3xZ6_1_1.json"), encoding="utf-8"))
    cells = r["cells"]
    print("\n=== GATE RESULT ===")
    print("cells enumerated:", len(cells))
    print("sector split:", Counter(tuple(c["kl"]) for c in cells))
    print("roles:", Counter(k for c in cells for k in c["printed_10"]))
    print("printed state found in own shell:",
          Counter(c["printed_states_found_in_shell"] for c in cells))
    print("shell multiplicity of printed 10-role:",
          Counter(v for c in cells for v in c["shell_multiplicity_of_printed_10"].values()))
    print("commutant SU(2) present:", Counter(c["has_commutant_SU2"] for c in cells))
    print("shell sizes:", min(c["n_shell"] for c in cells), "-",
          max(c["n_shell"] for c in cells))
    print("local groups:", Counter(c["group"] for c in cells))
    os.replace(os.path.join(GOUT, "shell_Z3xZ6_1_1.json"),
               os.path.join(GOUT, "d8_shell_complete_Z3xZ6_1_1.json"))


if __name__ == "__main__":
    main()

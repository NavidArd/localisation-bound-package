"""z3xz6_sweep -- CROSS-PHASE CHECK (CLAUDE.md 8.6, raised the turn it was found).

The Z3xZ6 sweep discovered that `create coupling(F F F)` forks a child inheriting the
redirected stdin, so in a MULTI-model orbifolder stream every model prints its
superpotential TWICE and the two prints are not always equal -- a model whose creation had
not finished when the first `print superpotential` executed prints `W = 0` first and the
real superpotential second.

`code/z6i_sweep/zsw_lib.py::parse_W` uses `_WPAT.search`, i.e. it takes the **FIRST**
print, and `code/z6i_sweep/run_couplings_z6i.sh` runs ALL models in ONE stream.  The Z6-I
sweep's headline was **SW-3 EMPTY** -- an ABSENCE claim -- which is exactly what a
truncated first print manufactures.

This script measures whether the exposure is REALISED on the Z6-I raw logs on disk: for
every model it extracts all `W = ` prints and reports any whose prints disagree, and
whether the FIRST print is the shorter one.

Run: python multi_geometry_campaign/code/z3xz6_sweep/check_z6i_exposure.py
"""

from __future__ import annotations

import glob
import gzip
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MGC = os.path.abspath(os.path.join(HERE, "..", ".."))
RAW = os.path.join(MGC, "data", "z6i_sweep", "raw")
MARK = re.compile(r"/>\s*/(MSSM\d+)>")
WPAT = re.compile(r"\n\s+W\s*=\s*(.*?)(?:\n\s*\n|\n/)", re.S)


def terms(txt):
    return tuple(sorted(tuple(re.findall(r"b?F_\d+", t))
                        for t in " ".join(txt.split()).split("+")
                        if re.findall(r"F_\d+", t)))


def main():
    rep = {"files": [], "models_with_multiple_prints": 0,
           "models_with_disagreeing_prints": [], "models_where_first_is_shorter": []}
    seen = {}
    for f in sorted(glob.glob(os.path.join(RAW, "*.gz"))):
        with gzip.open(f, "rt", errors="replace") as fh:
            t = fh.read()
        marks = [(m.start(), m.group(1)) for m in MARK.finditer(t)]
        rows = []
        for i, (p, lab) in enumerate(marks):
            b = t[p:(marks[i + 1][0] if i + 1 < len(marks) else len(t))]
            prints = [terms(m.group(1)) for m in WPAT.finditer(b)]
            if not prints:
                continue
            rows.append((lab, [len(x) for x in prints]))
            if len(prints) > 1:
                rep["models_with_multiple_prints"] += 1
                if len({x for x in prints}) > 1:
                    rep["models_with_disagreeing_prints"].append(
                        {"file": os.path.basename(f), "model": lab,
                         "print_sizes": [len(x) for x in prints]})
                    if len(prints[0]) < max(len(x) for x in prints):
                        rep["models_where_first_is_shorter"].append(lab)
            seen.setdefault(lab, []).append((os.path.basename(f),
                                             [len(x) for x in prints]))
        rep["files"].append({"file": os.path.basename(f), "models": len(rows),
                             "print_sizes": dict(rows)})
    rep["n_models_seen"] = len(seen)
    rep["verdict"] = (
        "EXPOSURE REALISED -- the Z6-I sweep's first-print parser read a TRUNCATED "
        "superpotential on at least one model, and its SW-3 EMPTY headline rests on "
        "absence claims that must be re-derived from the longest print."
        if rep["models_where_first_is_shorter"] else
        "EXPOSURE NOT REALISED on the logs on disk -- every model's prints agree, so the "
        "first-print parser happened to read the complete superpotential everywhere.  The "
        "CODE is still wrong in principle (it would silently truncate on a slower model) "
        "and should be changed to take the longest print, but no Z6-I number moves.")
    dest = os.path.join(MGC, "data", "z3xz6_sweep", "z6i_first_print_exposure.json")
    json.dump(rep, open(dest, "w"), indent=1)
    print("Z6-I first-print exposure check")
    print("  models seen                    :", rep["n_models_seen"])
    print("  models printing W more than once:", rep["models_with_multiple_prints"])
    print("  models with DISAGREEING prints  :", len(rep["models_with_disagreeing_prints"]))
    print("  models where FIRST is shorter   :", rep["models_where_first_is_shorter"])
    print("\n", rep["verdict"])
    print("->", dest)


if __name__ == "__main__":
    main()

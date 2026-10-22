"""z3xz6_sweep -- COVERAGE GATE over the published coupling segments.

Two independent extraction hazards were observed while running this sweep, both silent
and both of the "still exits 0, still looks complete" kind:

  (a) a shared-geometry race -- several segment scripts copying the geometry file into the
      SHARED orbifolder Geometry directory while sibling workers were reading it, so a
      worker loaded a half-written file and dropped part of its model list
      (`Orbifold "MSSMn" not known`);  cost: 16 census models, caught by audit_raw_z36.sh;
  (b) a per-batch label loss in the awk extraction (one label of 44 absent from the
      generated model file, a different one on two different runs).

Neither can be prevented by care alone, so the sweep does not rely on care: this gate
enumerates the models actually PRESENT in the published segments, compares against the
4,612 catalogue labels, and emits the missing list for `fill_gaps_z36.sh` to re-run.  It
is the gate that makes the individual extraction hazards non-fatal.

Run: python coverage_z36.py [--emit-missing]
"""

from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z36sw_lib as Z                                                  # noqa: E402

N = 4612


def present_models():
    have = {}
    for fn in sorted(os.listdir(Z.RAW_COUP)):
        if not (fn.startswith("coup_seg") and fn.endswith(".gz")):
            continue
        if fn.startswith("coup_segP"):                 # pilot / probe segments
            continue
        txt = Z.read_raw(os.path.join(Z.RAW_COUP, fn))
        for lab, block in Z.split_models(txt):
            # a model counts as covered only if it actually reached `print superpotential`
            if "\n  W = " in block:
                have.setdefault(lab, []).append(fn)
    return have


def unsignalled():
    """Models PRESENT in a segment but WITHOUT orbifolder's completion signal.

    These are the models whose forked coupling child was KILLED by `wait(N)`
    (`PID <n> stopped: ...`).  Their printed `W = 0` is a TRUNCATION, not an absence, and
    they must be re-run at a larger wait before any statement about them is used.
    """
    bad = []
    for fn_ in sorted(os.listdir(Z.RAW_COUP)):
        if not (fn_.startswith("coup_seg") and fn_.endswith(".gz")):
            continue
        if fn_.startswith("coup_segP"):
            continue
        for lab, block in Z.split_models(Z.read_raw(os.path.join(Z.RAW_COUP, fn_))):
            if ' done: "create coupling' not in block:
                bad.append(lab)
    return sorted(set(bad), key=lambda x: int(x[4:]))


def main():
    have = present_models()
    want = ["MSSM%d" % i for i in range(N)]
    missing = [m for m in want if m not in have]
    unsig = unsignalled()
    dup = {k: v for k, v in have.items() if len(v) > 1}
    rep = {"n_catalogue": N, "n_covered": len(have), "n_missing": len(missing),
           "missing": missing,
           "n_present_but_unsignalled": len(unsig),
           "present_but_unsignalled": unsig[:200],
           "n_models_in_more_than_one_segment": len(dup),
           "duplicated_examples": {k: dup[k] for k in list(dup)[:5]},
           "segments": sorted(set(f for v in have.values() for f in v))}
    dest = os.path.join(Z.DATA, "coverage.json")
    json.dump(rep, open(dest, "w"), indent=1)
    print("coupling coverage: %d / %d models   missing %d   in >1 segment %d"
          % (len(have), N, len(missing), len(dup)))
    if missing:
        print("missing (first 30):", missing[:30])
    print("->", dest)
    print("present but WITHOUT the completion signal (must be re-run): %d" % len(unsig))
    if "--emit-missing" in sys.argv:
        p = os.path.join(Z.DATA, "missing_models.txt")
        todo = missing + unsig            # both classes need a (longer-wait) re-run
        with open(p, "w", newline=chr(10)) as fh:
            fh.write(" ".join(sorted(set(todo), key=lambda x: int(x[4:]))) + chr(10))
        print("re-run list (missing + unsignalled) -> %s  (%d models)"
              % (p, len(set(todo))))


if __name__ == "__main__":
    main()

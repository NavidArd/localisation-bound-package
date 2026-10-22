"""z3xz6_sweep -- run the whole downstream pipeline and print the findings numbers.

Order matters and is enforced here:
  1. coverage_z36  -- which models actually have a COMPLETE coupling run (the completion
                      signal, not merely a printed `W`);
  2. sweep         -- per-model ladder analysis over the published segments;
  3. aggregate     -- catalogue roll-up, gates X1/X2/X3, C(S_min) branch statistics.

Census-derived results (rank bound, counterfactual, three-layer diagnosis) do NOT depend
on any of this and are produced separately; they are re-printed here so the findings can
be assembled from one place.

Run: python multi_geometry_campaign/code/z3xz6_sweep/finalize.py [--workers N]
"""

from __future__ import annotations

import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z36sw_lib as Z                                                  # noqa: E402

PY = sys.executable


def run(script, *args):
    print("\n" + "=" * 88)
    print("### %s %s" % (script, " ".join(args)))
    print("=" * 88)
    r = subprocess.run([PY, os.path.join(HERE, script), *args],
                       capture_output=True, text=True)
    sys.stdout.write(r.stdout[-8000:])
    if r.returncode:
        sys.stdout.write(r.stderr[-3000:])
    return r.returncode


def main():
    workers = "6"
    if "--workers" in sys.argv:
        workers = sys.argv[sys.argv.index("--workers") + 1]
    run("coverage_z36.py", "--emit-missing")
    run("sweep.py", "--workers", workers)
    run("aggregate.py")
    print("\n" + "=" * 88)
    print("### census-derived results (coverage-independent)")
    print("=" * 88)
    for f in ("rank_bound.json", "counterfactual.json", "gamma_on_legs.json"):
        p = os.path.join(Z.DATA, f)
        if os.path.exists(p):
            d = json.load(open(p))
            print("--- %s" % f)
            print(json.dumps(d.get("summary", d), indent=1, default=str)[:2500])


if __name__ == "__main__":
    main()

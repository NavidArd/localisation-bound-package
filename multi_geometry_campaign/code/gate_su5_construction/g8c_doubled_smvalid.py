"""How many of the 4,032 Standard-Model-valid loophole models carry a doubled norm-6/5 role.

Why this script exists
----------------------
The construction's headline negative is stated over the **Standard-Model-valid** models of
the complete two-class loophole space -- 2,912 of class A plus 1,120 of class B, 4,032 in
all -- and the count of those that carry a doubled norm-6/5 role is quoted alongside it.
That count was carried in prose without a producing script, and the figure in circulation
(1,482) does not follow from the deposited tables. This script derives it from them.

  class A   every one of the 2,912 built embeddings is Standard-Model-valid, so the
            doubled count is `s8_battery.json::n_with_T0l_ten` unrestricted.

  class B   only 1,120 of the 1,432 built lifts are Standard-Model-valid
            (`g8b_smvalid.json::valid`), so the doubled count must be taken over that
            subset, not over all 1,432. `g8_batteryB.json::summary.doubled_10role` = 400
            is the unrestricted figure and is NOT the one the headline needs.

The restriction matters: 400 unrestricted against 362 on the valid subset.

Inputs  (all deposited, all read-only)
    ../../data/su5_construction/s8_battery.json
    ../../data/gate_su5_construction/g8_batteryB.json
    ../../data/gate_su5_construction/g8b_smvalid.json
Output
    ../../data/gate_su5_construction/g8c_doubled_smvalid.json

Note on scope: this script re-uses `g8b_smvalid.json`, which is itself the reproduction gap
R-1 recorded in `construction/README.md` section 3 -- no script on disk produces it. That gap
is unchanged by this script; what changes is that the doubled count is no longer prose.
"""
from __future__ import annotations

import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "..", "data"))
GSC = os.path.join(DATA, "gate_su5_construction")
SC = os.path.join(DATA, "su5_construction")


def load(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def main() -> int:
    a = load(os.path.join(SC, "s8_battery.json"))
    b = load(os.path.join(GSC, "g8_batteryB.json"))
    smv = load(os.path.join(GSC, "g8b_smvalid.json"))

    valid = set(smv["valid"])
    rows = b["rows"]

    a_doubled = a["n_with_T0l_ten"]
    a_models = a["n_models"] if "n_models" in a else 2912

    b_built = len(rows)
    b_doubled_all = sum(1 for r in rows if r["hits"] > 0)
    b_doubled_valid = sum(1 for r in rows if r["hits"] > 0 and r["model"] in valid)

    out = {
        "classA": {
            "models": a_models,
            "all_sm_valid": True,
            "doubled": a_doubled,
            "source": "su5_construction/s8_battery.json::n_with_T0l_ten",
        },
        "classB": {
            "built": b_built,
            "sm_valid": len(valid),
            "doubled_over_all_built": b_doubled_all,
            "doubled_over_sm_valid": b_doubled_valid,
            "source": "gate_su5_construction/g8_batteryB.json::rows[hits>0]"
                      " intersected with g8b_smvalid.json::valid",
        },
        "sm_valid_models_total": a_models + len(valid),
        "doubled_over_sm_valid_total": a_doubled + b_doubled_valid,
        "doubled_over_all_built_total": a_doubled + b_doubled_all,
        "consistency": {
            "classB_summary_doubled_10role": b["summary"]["doubled_10role"],
            "agrees_with_recount_over_all_built":
                b["summary"]["doubled_10role"] == b_doubled_all,
        },
    }

    path = os.path.join(GSC, "g8c_doubled_smvalid.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")

    print("class A  : %d models, all SM-valid, %d doubled" % (a_models, a_doubled))
    print("class B  : %d built, %d SM-valid; doubled %d over all built, %d over SM-valid"
          % (b_built, len(valid), b_doubled_all, b_doubled_valid))
    print("TOTAL    : %d SM-valid models, %d carry a doubled norm-6/5 role"
          % (out["sm_valid_models_total"], out["doubled_over_sm_valid_total"]))
    print("(over all %d BUILT models, unrestricted, the figure is %d)"
          % (a_models + b_built, out["doubled_over_all_built_total"]))
    print("wrote", path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

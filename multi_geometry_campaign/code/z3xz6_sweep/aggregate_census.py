"""z3xz6_sweep -- catalogue-level roll-up of the STEP-2 role census.

Emits `data/z3xz6_sweep/census_summary.json`: the role statistics over all 4,612 models,
the Rule-P classification, the carrier-channel occupancy distribution, the orientation
(sY, sC) branch statistics, and the D-2 uncovered-element tagging.

Every number here is a CATALOGUE statement (`z36sw_lib.CATALOGUE_TAG`), never a geometry
statement: channel ALLOWEDNESS is a property of the space group; channel OCCUPANCY is a
property of the model.
"""

from __future__ import annotations

import json
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z36sw_lib as Z                                                  # noqa: E402

CENSUS = os.path.join(Z.DATA, "census")


def main():
    rows = []
    for fn in sorted(os.listdir(CENSUS)):
        if fn.endswith(".json"):
            rows.extend(json.load(open(os.path.join(CENSUS, fn)))["models"])
    rows.sort(key=lambda r: int(r["model"][4:]))
    n = len(rows)
    ok = [r for r in rows if r.get("status") == "ok"]

    s = {"scope": Z.CATALOGUE_TAG,
         "n_models": n, "n_status_ok": len(ok),
         "status_histogram": dict(Counter(r.get("status") for r in rows))}

    s["net_Q_is_plus_three"] = sum(
        1 for r in ok if r["orientation"]["net_Q_is_plus_three"])
    s["net_all_four_species_three"] = sum(
        1 for r in ok if r["orientation"]["net_all_three"])
    s["orientation_branch_selected"] = dict(Counter(
        "sY=%d,sC=%d" % (r["orientation"]["sY_selected"],
                         r["orientation"]["sC_selected"]) for r in ok))
    s["sC_flipped_relative_to_naive"] = sum(
        1 for r in ok if r["orientation"]["sC_selected"] == -1)

    s["chiral_Q_distribution"] = dict(sorted(Counter(r["n_chiral_Q"] for r in ok).items()))
    s["chiral_uc_distribution"] = dict(sorted(Counter(r["n_chiral_uc"] for r in ok).items()))
    s["Hu_distribution"] = dict(sorted(Counter(r["n_Hu"] for r in ok).items()))
    s["models_with_3plus_chiral_Q_and_uc_and_a_Hu"] = sum(
        1 for r in ok if r["n_chiral_Q"] >= 3 and r["n_chiral_uc"] >= 3
        and r["n_Hu"] >= 1)

    sec_q, sec_u, sec_h = Counter(), Counter(), Counter()
    for r in ok:
        for k, v in r["sectors"]["chiral_Q"].items():
            sec_q[k] += v
        for k, v in r["sectors"]["chiral_uc"].items():
            sec_u[k] += v
        for k, v in r["sectors"]["Hu"].items():
            sec_h[k] += v
    s["sector_occupancy_chiral_Q"] = dict(sec_q.most_common())
    s["sector_occupancy_chiral_uc"] = dict(sec_u.most_common())
    s["sector_occupancy_Hu"] = dict(sec_h.most_common())
    s["models_with_a_chiral_Q_in_a_right_chiral_label"] = sum(
        1 for r in ok for k in r["sectors"]["chiral_Q"]
        if tuple(eval(k)) in Z.RIGHT_CHIRAL_LABELS)

    s["ruleP_pass"] = sum(1 for r in ok if r.get("ruleP_pass"))
    s["ruleP_fail"] = len(ok) - s["ruleP_pass"]
    chp = Counter()
    for r in ok:
        for k, v in (r.get("ruleP_ladder_channel_triples") or {}).items():
            chp[k] += 1
    s["ruleP_models_per_ladder_channel"] = dict(chp.most_common())
    s["ruleP_carrier_models"] = sum(
        1 for r in ok if r.get("ruleP_carrier_triples"))
    ach = Counter()
    for r in ok:
        for k in (r.get("allowed_channel_triples") or {}):
            ach[k] += 1
    s["models_with_a_QucHu_triple_on_ANY_allowed_channel"] = dict(ach.most_common())
    s["n_models_with_any_allowed_QucHu_triple"] = sum(
        1 for r in ok if r.get("allowed_channel_triples"))

    s["D2_uncovered_tag_models"] = sum(1 for r in ok if r.get("D2_uncovered_tag"))
    s["gamma_denominators_union"] = sorted({d for r in ok
                                            for d in r.get("gamma_denominators", [])})
    s["oscillator_fields_distribution"] = dict(sorted(Counter(
        r["n_fields_with_oscillators"] for r in ok).items()))
    s["models_with_a_gamma_pair_witness"] = sum(
        1 for r in ok if r.get("n_fields_both_twist_gammas_nonzero", 0) > 0)

    s["ladder_bearing_channels"] = {
        k: {"total_exponents": v["total_exponents"],
            "computable_planes": v["computable_planes"],
            "tau_invariant": v.get("tau_invariant")}
        for k, v in Z.ladder_channels().items() if v["ladder_bearing"]}

    # Also emit the Rule-P model list as a plain, LF-terminated file for the shell
    # drivers.  A Windows text-mode write turns the newline into CR+LF, and the
    # shell's $(cat file) then hands awk a label with a trailing CR -- which
    # silently loses exactly the LAST label in the list.  That cost MSSM4569 from
    # the priority batch before an extraction probe caught it.  The explicit
    # newline argument below is load-bearing, not cosmetic.
    with open(os.path.join(Z.DATA, "rulep_models.txt"), "w", newline=chr(10)) as fh:
        fh.write(" ".join(r["model"] for r in ok if r.get("ruleP_pass")) + chr(10))

    dest = os.path.join(Z.DATA, "census_summary.json")
    json.dump({"summary": s,
               "ruleP_pass_models": [r["model"] for r in ok if r.get("ruleP_pass")],
               "models_with_any_allowed_triple":
                   [r["model"] for r in ok if r.get("allowed_channel_triples")]},
              open(dest, "w"), indent=1, default=str)
    print(json.dumps(s, indent=1, default=str))
    print("->", dest)


if __name__ == "__main__":
    main()

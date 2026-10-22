"""z3xz6_sweep -- STEP 5: the structural miss diagnosis.

Answers, from the CENSUS ALONE and at the level of a proof rather than a tally, WHY a
given ladder-bearing channel is or is not reachable by an up-type Yukawa across the
catalogue.  The argument has the same three-layer shape the Z6-I sweep's diagnosis used:

  LAYER 1 (geometry, model-independent).  Which channels are ALLOWED by the published
          per-plane R-charge law, which are LADDER-BEARING on the TOTAL exponent, and
          which are EMPTY BY GEOMETRY (the all-T(0,l) family).  No model is consulted.
  LAYER 2 (catalogue occupancy, per leg).  For each ladder-bearing channel and each leg,
          how many models place a CHIRAL Q / a CHIRAL u^c / an H_u on that leg's sector.
          A leg that NO model can supply with a chiral Q kills the Q assignment at that
          leg for the whole catalogue.
  LAYER 3 (joint occupancy).  Whether the required assignment occurs in a SINGLE model:
          the per-leg counts can all be non-zero while their intersection is empty, and
          that is a strictly stronger statement than any per-leg tally.

The output separates GEOMETRY statements from CATALOGUE statements at every line, and
never promotes a catalogue emptiness into a geometry claim.

Run: python multi_geometry_campaign/code/z3xz6_sweep/diagnosis.py
"""

from __future__ import annotations

import json
import os
import sys
from collections import Counter, defaultdict
from itertools import permutations

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z36sw_lib as Z                                                  # noqa: E402

CENSUS = os.path.join(Z.DATA, "census")


def load_rows():
    rows = []
    for fn in sorted(os.listdir(CENSUS)):
        if fn.endswith(".json"):
            rows.extend(json.load(open(os.path.join(CENSUS, fn)))["models"])
    return [r for r in rows if r.get("status") == "ok"]


def main():
    rows = load_rows()
    lc = Z.ladder_channels()

    # per-sector catalogue occupancy (models, not fields)
    occ = {"chiral_Q": Counter(), "chiral_uc": Counter(), "Hu": Counter()}
    have = {"chiral_Q": defaultdict(set), "chiral_uc": defaultdict(set),
            "Hu": defaultdict(set)}
    for r in rows:
        for role in occ:
            for k in r["sectors"][role]:
                occ[role][k] += 1
                have[role][k].add(r["model"])

    rep = {"scope_geometry": Z.GEOMETRY_TAG, "scope_catalogue": Z.CATALOGUE_TAG,
           "n_models_ok": len(rows),
           "layer1_geometry": {
               "n_channels": 66,
               "n_allowed_published_law": len(Z.allowed_channels()),
               "n_allowed_all_legs_twisted": len(Z.allowed_channels_twisted()),
               "n_ladder_bearing": sum(1 for v in lc.values() if v["ladder_bearing"]),
               "empty_by_geometry_channels":
                   sorted({Z.channel_name(c) for c in Z.EMPTY_BY_GEOMETRY}),
               "empty_by_geometry_reason": Z.EMPTY_BY_GEOMETRY_REASON,
               "nine_factor_channels": sorted(
                   k for k, v in lc.items() if v["n_total"] == 9),
               "tau_invariant_ladder_channels": sorted(
                   k for k, v in lc.items()
                   if v["ladder_bearing"] and v.get("tau_invariant"))},
           "layer2_per_leg_occupancy": {},
           "layer3_joint": {}}

    for name, v in sorted(lc.items()):
        if not v["ladder_bearing"]:
            continue
        legs = [tuple(t) for t in v["tri"]]
        rep["layer2_per_leg_occupancy"][name] = {
            str(g): {"models_with_a_chiral_Q_here": occ["chiral_Q"][str(g)],
                     "models_with_a_chiral_uc_here": occ["chiral_uc"][str(g)],
                     "models_with_an_Hu_here": occ["Hu"][str(g)]}
            for g in legs}

        # LAYER 3: which leg->role assignments are even conceivable, and which occur
        assigns = []
        for p in set(permutations(range(3))):
            gq, gu, gh = legs[p[0]], legs[p[1]], legs[p[2]]
            nq, nu, nh = (occ["chiral_Q"][str(gq)], occ["chiral_uc"][str(gu)],
                          occ["Hu"][str(gh)])
            joint = (have["chiral_Q"][str(gq)] & have["chiral_uc"][str(gu)]
                     & have["Hu"][str(gh)])
            assigns.append({"Q_leg": str(gq), "uc_leg": str(gu), "Hu_leg": str(gh),
                            "models_per_leg": [nq, nu, nh],
                            "models_with_ALL_THREE": len(joint),
                            "example_models": sorted(joint)[:5]})
        best = max(a["models_with_ALL_THREE"] for a in assigns)
        killer = [a for a in assigns if 0 in a["models_per_leg"]]
        rep["layer3_joint"][name] = {
            "assignments": assigns,
            "max_models_with_a_complete_assignment": best,
            "n_assignments_killed_at_LAYER_2": len(killer),
            "n_assignments_total": len(assigns),
            "verdict": ("reachable in %d model(s)" % best if best else
                        ("UNREACHABLE catalogue-wide -- every leg->role assignment is "
                         "killed at LAYER 2 (some leg's sector hosts a chiral Q / u^c / "
                         "H_u in NO model)" if len(killer) == len(assigns) else
                         "UNREACHABLE catalogue-wide -- every per-leg count is non-zero "
                         "but NO SINGLE MODEL supplies all three simultaneously "
                         "(LAYER-3 emptiness, strictly stronger than any per-leg tally)"))}

    # sectors that host NO chiral Q anywhere in the catalogue
    allsec = {str((k, l)) for k in range(3) for l in range(6)}
    rep["sectors_hosting_no_chiral_Q_in_any_model"] = sorted(
        s for s in allsec if occ["chiral_Q"][s] == 0)
    rep["sectors_hosting_no_chiral_uc_in_any_model"] = sorted(
        s for s in allsec if occ["chiral_uc"][s] == 0)
    rep["sectors_hosting_no_Hu_in_any_model"] = sorted(
        s for s in allsec if occ["Hu"][s] == 0)
    rep["note_right_chiral_labels"] = (
        "T(1,1), T(1,2), T(1,3), T(2,1) print only as CPT conjugates on this geometry "
        "(measured left-chiral basis, foundations addendum A4).  A zero count on those "
        "four labels is LABELLING, NOT ABSENCE, and no absence claim rests on them.")
    rep["right_chiral_labels"] = [str(g) for g in Z.RIGHT_CHIRAL_LABELS]

    dest = os.path.join(Z.DATA, "diagnosis.json")
    json.dump(rep, open(dest, "w"), indent=1, default=str)

    print("=" * 92)
    print("Z3xZ6_1_1 sweep -- structural miss diagnosis (three layers)")
    print("=" * 92)
    L1 = rep["layer1_geometry"]
    print("LAYER 1 (GEOMETRY): 66 channels | %d allowed | %d allowed & all legs twisted "
          "| %d ladder-bearing" % (L1["n_allowed_published_law"],
                                   L1["n_allowed_all_legs_twisted"],
                                   L1["n_ladder_bearing"]))
    print("  nine-factor channels: %s" % L1["nine_factor_channels"])
    print("  tau-invariant ladder channels: %s" % L1["tau_invariant_ladder_channels"])
    print("  empty by geometry (PRE-DECLARED, parser check only): %s"
          % L1["empty_by_geometry_channels"])
    print("\nLAYER 2/3 (CATALOGUE, %d models):" % rep["n_models_ok"])
    print("%-22s %-38s %s" % ("channel", "per-leg (chiralQ / chiral uc / Hu)", "verdict"))
    for name in sorted(rep["layer3_joint"]):
        occs = rep["layer2_per_leg_occupancy"][name]
        s = " ".join("%s:%d/%d/%d" % (g, o["models_with_a_chiral_Q_here"],
                                      o["models_with_a_chiral_uc_here"],
                                      o["models_with_an_Hu_here"])
                     for g, o in occs.items())
        print("%-22s %-38s %s" % (name, s[:38], rep["layer3_joint"][name]["verdict"]))
    print("\nsectors hosting NO chiral Q in ANY of the %d models: %s"
          % (rep["n_models_ok"], rep["sectors_hosting_no_chiral_Q_in_any_model"]))
    print("sectors hosting NO chiral u^c in ANY model:            %s"
          % rep["sectors_hosting_no_chiral_uc_in_any_model"])
    print("(right-chiral LABELS, absence is labelling not physics: %s)"
          % rep["right_chiral_labels"])
    print("->", dest)


if __name__ == "__main__":
    main()

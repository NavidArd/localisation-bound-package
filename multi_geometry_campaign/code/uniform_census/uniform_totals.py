r"""Reduce the two uniform-convention accumulators to the both-catalogue totals the paper prints.

Reads   results_deposited/uniform_census/uniform_e6_<geom>.json  (written by
        uniform_census_e6.py)
Writes  results_deposited/uniform_census/uniform_totals.json

Every entry is a sum over the two geometries of one accumulator key, or a ratio of two
such sums; the key each one came from is recorded beside it in the `from` block, so no
number here is typed by hand. `claims.json` resolves the Table 1 rows against this file.

Usage:  python uniform_totals.py
"""
from __future__ import annotations

import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
D = os.path.join(REPO, "results_deposited", "uniform_census")
GEOMS = ("Z2xZ4_1_1", "Z3xZ6_1_1")

# name -> path inside the `uniform` accumulator
KEYS = {
    "models": ("models",),
    "twisted_model_cells": ("mcells_tw",),
    "fields_10_rec": ("f_rec_by_fam", "10"),
    "fields_5bar_rec": ("f_rec_by_fam", "5bar"),
    "fields_10_st": ("f_st_by_fam", "10"),
    "fields_5bar_st": ("f_st_by_fam", "5bar"),
    "slots_10": ("slots_by_fam", "10"),
    "slots_5bar": ("slots_by_fam", "5bar"),
    "below_slots_5bar": ("below_slots", "5bar"),
    "cells_with_a_10_role": ("occ", "ten_any"),
    "cells_with_a_5bar_role": ("occ", "five_any"),
    "violations_record": ("viol_rec",),
    "violations_state": ("viol_st",),
    "unprotected_10_fields": ("loophole_f_rec",),
    "unprotected_10_fields_state": ("loophole_f_st",),
    "pairs_rec_10": ("pairs_rec_by_fam", "10"),
    "pairs_rec_5bar": ("pairs_rec_by_fam", "5bar"),
    "pairs_st_10": ("pairs_st_by_fam", "10"),
    "pairs_st_5bar": ("pairs_st_by_fam", "5bar"),
    "pairs_st_cross_5bar": ("pairs_st_cross_by_fam", "5bar"),
    "delta_not_in_lambda": ("d_notlambda",),
    "delta_not_perp_W": ("d_notperp",),
    "delta_pairs_tested": ("d_pairs",),
    "budget_identity_checked": ("I1_checked",),
    "budget_identity_bad": ("I1_bad",),
    "smbasis_fail": ("smbasis_fail",),
    "dup_groups": ("dup_groups",),
    "dup_cells": ("dup_cells",),
    "dup_groups_carrying_a_role": ("dup_groups_role",),
    "ceiling_overruns": ("ceil_over",),
    "kY_5_over_3": ("kY", "5/3"),
    "models_recovered_at_a_later_embedding": ("models_recovered",),
    "constructing_elements_twisted": ("geo_cells_tw",),
}


def dig(acc, path):
    v = acc
    for k in path:
        v = v[k]
    return int(v)


def main():
    acc = {g: json.load(open(os.path.join(D, "uniform_e6_%s.json" % g)))["uniform"]
           for g in GEOMS}
    out = {"note": ("Both-catalogue totals of the uniform-convention census with the E_6 "
                    "parse repaired. Each value is the sum over the two geometries of the "
                    "accumulator key named in `from`; the two ratios are quotients of two "
                    "such sums."),
           "geometries": list(GEOMS), "from": {k: list(v) for k, v in KEYS.items()},
           "per_geometry": {}, "totals": {}}
    for name, path in KEYS.items():
        per = {g: dig(acc[g], path) for g in GEOMS}
        out["per_geometry"][name] = per
        out["totals"][name] = sum(per.values())

    t = out["totals"]
    t["role_carrying_fields_rec"] = t["fields_10_rec"] + t["fields_5bar_rec"]
    t["role_carrying_fields_st"] = t["fields_10_st"] + t["fields_5bar_st"]
    t["pairs_rec"] = t["pairs_rec_10"] + t["pairs_rec_5bar"]
    t["pairs_st"] = t["pairs_st_10"] + t["pairs_st_5bar"]
    t["pairs_st_cross"] = t["pairs_st_cross_5bar"] + 0
    t["ratio_10_pairs_per_cell"] = round(t["slots_10"] / t["cells_with_a_10_role"], 3)
    t["ratio_5bar_pairs_per_cell"] = round(t["slots_5bar"] / t["cells_with_a_5bar_role"], 3)

    # mixed-budget (cell, role) groups -- the numeral of the ceiling annexe's caption.
    # Keyed "family|R_max" for the groups and "family|R_max|R" for the fields inside them.
    out["mixed_budget"] = {g: {"groups": acc[g]["mixed_groups"],
                               "fields": acc[g]["mixed_fields"]} for g in GEOMS}

    # the two multiplicity histograms, per geometry, as deposited
    out["histograms"] = {g: {"mult_records": acc[g]["mult_rec_by_fam"],
                             "mult_states": acc[g]["mult_st_by_fam"],
                             "hidden_dimension": acc[g]["hid_by_fam"],
                             "R_by_family": acc[g]["R_rec_by_fam"]} for g in GEOMS}

    with open(os.path.join(D, "uniform_totals.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for k in sorted(t):
        print("%-40s %s" % (k, t[k]))


if __name__ == "__main__":
    main()

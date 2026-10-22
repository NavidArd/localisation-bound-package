"""
finalize.py --- assembles the phase's funnel and the ladder verdict from the on-disk payloads.

Reads only what other modules wrote; computes no physics of its own.  Its purpose is that the
findings' funnel and letter are produced from the JSON on disk rather than transcribed from a
console log (the campaign's orchestrator discipline: verify deliverables on disk).

Emits: data/z2xz4_sweep/funnel.json
"""

from __future__ import annotations

import gzip
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z24sw_lib as L                                                   # noqa: E402


def j(name):
    p = os.path.join(L.DSW, name)
    return json.load(open(p, encoding="utf-8")) if os.path.exists(p) else {}


def main():
    st = j("selftest.json")
    cen = j("census_summary.json")
    lev = j("levels.json")
    blk = j("blocks_summary.json")
    swp = j("sweep_summary.json")
    prot = j("protection.json")
    diag = j("diagnosis.json")
    mh = j("multihiggs.json")
    bf = j("bf_basis.json")
    pil = j("pilot_timing.json")

    funnel = {
        "catalogue (SHA-256 %s, %d bytes)" % (L.CATALOGUE_SHA256[:8], L.CATALOGUE_BYTES):
            L.N_MODELS,
        "orbifolder census, 58 segments, corrupt-geometry lines": 0,
        "models seen by the census": cen.get("n_models_seen"),
        "evaluable (SMConfig closure reproduced)":
            cen.get("status_histogram", {}).get("OK"),
        "  excluded: no SM vev-configuration at all":
            cen.get("status_histogram", {}).get("NO_SMCONFIG", 0),
        "  excluded: SM_ID DATA-GAP":
            cen.get("status_histogram", {}).get("SM_ID_DATA_GAP", 0),
        "net Q = +3 on the selected orientation branch": cen.get("n_net_Q_plus_three"),
        "argmax ties on the orientation branch": cen.get("n_argmax_ties_total"),
        "off-lattice fixed points": cen.get("n_off_lattice_fixed_points"),
        "INDEX LEVEL (all 15 channels x all 6 leg assignments x every H_u field)": "",
        "  models assembling any block": blk.get("n_models_with_any_block"),
        "  index-level blocks": blk.get("n_blocks"),
        "  models reaching generic rank 3": blk.get("n_models_index_rank_3"),
        "  models carrying a 0:n:2n ladder": blk.get("n_models_index_ladder_0_n_2n"),
        "  models carrying the pre-declared 0:2:4": blk.get("n_models_index_ladder_0_2_4"),
        "COUPLING LEVEL": "",
        "  models coupled": swp.get("n_models_analysed"),
        "  status histogram": swp.get("status_histogram"),
        "  label mismatches vs the census": swp.get("label_mismatches"),
        "  terms != couplings created": swp.get("n_terms_vs_couplings_created_mismatches"),
        "  realised blocks": swp.get("n_blocks"),
        "  models with a realised 0:n:2n ladder": swp.get("n_models_any_ladder"),
        "  models with the pre-declared 0:2:4": swp.get("n_models_SW1_0_2_4"),
        "SW-1": swp.get("n_models_SW1_0_2_4"),
    }

    letter = "SW-3 EMPTY"
    if swp.get("n_models_SW1_0_2_4"):
        letter = "SW-1 REALIZED"
    elif swp.get("n_models_any_ladder"):
        letter = "SW-1 REALIZED (n != 2 ladder; report the n)"
    elif (blk.get("n_models_with_any_block") or 0) > 0:
        letter = "SW-2 PARTIAL or SW-3 EMPTY -- decided by whether any up-structure assembles"

    rep = {"geometry": "Z2xZ4_1_1",
           "scope": "CATALOGUE statement about the 11,540 PVW Z2xZ4_1_1 models.  The "
                    "GEOMETRY-level admissibility is levels.json and is never conflated "
                    "with this.",
           "funnel": funnel,
           "ladder_letter_computed_from_disk": letter,
           "geometry_admits": {
               "channels_admitting_0_2_4":
                   lev.get("L_B_admissibility", {}).get("channels_admitting_0_2_4"),
               "n_blocks_admitting_0_2_4":
                   lev.get("L_B_admissibility", {}).get("n_blocks_admitting_0_2_4"),
               "NC_corrected_three_rung_buckets":
                   lev.get("L_A_buckets_NC", {}).get("n_ladder_buckets"),
               "buckets_without_NC":
                   lev.get("L_A_buckets_no_NC", {}).get("n_ladder_buckets")},
           "protection": {
               "exchange_pair": prot.get("exchange_pair"),
               "n_sigma_conditions": prot.get("n_sigma_conditions"),
               "PROTECTED": prot.get("n_protected_any_sigma"),
               "IMPOSED": prot.get("n_imposed_any_sigma"),
               "UNDETERMINED": prot.get("n_undetermined")},
           "csmin": {"annihilations_per_branch": swp.get("csmin_annihilations_per_branch"),
                     "capacity": swp.get("csmin_capacity")},
           "multi_higgs": {"groups": mh.get("n_multi_higgs_groups"),
                           "uniform_vev_ladders": mh.get("n_uniform_vev_ladders"),
                           "purchasable_ladders": mh.get("n_purchasable_ladders"),
                           "purchase_cost": mh.get("purchase_cost_histogram")},
           "left_chiral_basis": {"probed": bf.get("n_models_probed"),
                                 "right_chiral_states": bf.get(
                                     "total_right_chiral_states_printed"),
                                 "dead_sectors": bf.get(
                                     "sectors_with_no_left_chiral_field_in_any_of_the_11524_models")},
           "cpt_effective_denominator":
               diag.get("cpt_audit", {}).get("honest_effective_denominator"),
           "selftest": st.get("VERDICT"),
           "pilot": pil}
    L.write_json("funnel.json", rep)
    print(json.dumps(rep, indent=1, default=str)[:6000])


if __name__ == "__main__":
    main()

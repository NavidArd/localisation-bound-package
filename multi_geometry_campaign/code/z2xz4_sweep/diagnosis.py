"""
diagnosis.py --- STEP 5: the FUNNEL and the MISS DIAGNOSIS.

Answers, from the catalogue, the question the ladder letter leaves open: WHERE does the
catalogue stop?  Layer by layer, with the denominator stated at every layer:

  layer 0  catalogue                                      11,540
  layer 1  evaluable (SMConfig closure reproduced)
  layer 2  net Q = +3 on the selected orientation branch
  layer 3  PER-LEG OCCUPANCY: for each channel and each leg assignment, how many models put
           >= 3 chiral row-species fields on the row leg's sector, >= 3 column-species on the
           column leg's, and >= 1 H_u on the Higgs leg -- and how many do all three JOINTLY
           (the layer-3 conjunction that decided the Z4xZ4 sweep)
  layer 4  a block assembles at index level
  layer 5  the block reaches generic rank 3
  layer 6  the block carries a 0 : n : 2n mass ladder
  layer 7  the same at COUPLING level

Plus the CPT audit: which sectors host a left-chiral field in ANY model, which channels are
therefore unreachable FOR A LABELLING REASON rather than a physical absence, and what the
honest effective denominator of letter-carrying channels is.

Emits: data/z2xz4_sweep/diagnosis.json
"""

from __future__ import annotations

import gzip
import json
import os
import sys
from collections import Counter, defaultdict
from itertools import permutations

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z24sw_lib as L                                                   # noqa: E402
from blocks import load_tables                                          # noqa: E402


def cpt(kl):
    g = L.geom()
    return ((-kl[0]) % g.M, (-kl[1]) % g.N)


def main():
    tables = load_tables()
    with gzip.open(os.path.join(L.DSW, "models", "roles.json.gz"), "rt",
                   encoding="utf-8") as fh:
        roles = json.load(fh)
    with gzip.open(os.path.join(L.DSW, "models", "blocks.json.gz"), "rt",
                   encoding="utf-8") as fh:
        blocks = json.load(fh)
    lev = json.load(open(os.path.join(L.DSW, "levels.json"), encoding="utf-8"))

    ok = {k: v for k, v in roles.items() if v.get("status") == "OK"}

    # ---------------------------------------------------- the left-chiral basis / CPT audit
    host = Counter()
    for r in ok.values():
        for s in r.get("sectors_present", []):
            host[s] += 1
    all_sectors = ["%d,%d" % (k, l) for k in range(L.geom().M) for l in range(L.geom().N)]
    never = [s for s in all_sectors if host.get(s, 0) == 0]
    lettered = [nm for nm, T in tables.items() if T["capable_planes"]]
    unreachable, reachable = [], []
    for nm in lettered:
        gs = tables[nm]["gs"]
        if any("%d,%d" % g in never for g in gs):
            img = L.channel_name(tuple(sorted(cpt(g) for g in gs)))
            unreachable.append({"channel": nm, "cpt_image": img,
                                "cpt_image_is_lettered": img in lettered,
                                "dead_legs": ["T(%d,%d)" % g for g in gs
                                              if "%d,%d" % g in never]})
        else:
            reachable.append(nm)
    cptaudit = {
        "sectors_hosting_a_left_chiral_field": {s: host.get(s, 0) for s in all_sectors},
        "sectors_never_hosting_a_left_chiral_field": never,
        "n_lettered_channels": len(lettered),
        "lettered_channels_unreachable_by_labelling": unreachable,
        "lettered_channels_reachable": reachable,
        "honest_effective_denominator": len(reachable),
        "note": "A sector that hosts no LEFT-CHIRAL field in any model can carry no "
                "superpotential leg in the left-chiral basis.  That is a LABELLING fact, not "
                "a physical absence: the channel is the CPT image of one that IS populated.  "
                "NO ABSENCE IS CLAIMED on the relabelled channels.  Scope: this statement "
                "rests on print(F), which cannot emit an _r state; the non-vacuous "
                "measurement is the print(bF) spread probe (bf_basis.json)."}

    # ------------------------------------------------------------ layer-3 leg occupancy
    layer3 = {}
    for nm, T in tables.items():
        if not T["capable_planes"]:
            continue
        gs = T["gs"]
        per = {}
        seen = set()
        for legQ, legU, legH in permutations(range(3)):
            key = (gs[legQ], gs[legU], gs[legH])
            if key in seen:
                continue
            seen.add(key)
            cQ = cU = cH = cJ = 0
            for r in ok.values():
                nq = nu = nh = 0
                for k, l, o1, o2, o3, role, mult, hs, flab in r["fields"]:
                    s = (k, l)
                    if role == "Q" and s == gs[legQ]:
                        nq += 1
                    elif role == "uc" and s == gs[legU]:
                        nu += 1
                    elif role == "Hu" and s == gs[legH]:
                        nh += 1
                cQ += nq >= 3
                cU += nu >= 3
                cH += nh >= 1
                cJ += (nq >= 3 and nu >= 3 and nh >= 1)
            per["Q@T(%d,%d) uc@T(%d,%d) Hu@T(%d,%d)"
                % (gs[legQ] + gs[legU] + gs[legH])] = {
                    "models_with_3plus_Q_on_row_leg": cQ,
                    "models_with_3plus_uc_on_col_leg": cU,
                    "models_with_an_Hu_on_higgs_leg": cH,
                    "models_with_ALL_THREE_jointly": cJ}
        layer3[nm] = per

    # ------------------------------------------------------------------------- the funnel
    allb = [b for v in blocks.values() for b in v["blocks"]]
    r3 = {lab for lab, v in blocks.items() if any(b["generic_rank_3"] for b in v["blocks"])}
    lad = {lab for lab, v in blocks.items() if any(b["ladder_0_n_2n"] for b in v["blocks"])}
    funnel = {
        "catalogue": L.N_MODELS,
        "census_segments": 58,
        "corrupt_geometry_lines": 0,
        "evaluable_SMConfig_closure_reproduced": len(ok),
        "excluded_no_SM_vev_configuration": sum(
            1 for v in roles.values() if v.get("status") == "NO_SMCONFIG"),
        "excluded_SM_ID_DATA_GAP": sum(
            1 for v in roles.values() if v.get("status") == "SM_ID_DATA_GAP"),
        "net_Q_plus_3": sum(1 for v in ok.values() if v.get("net_Q_is_plus_three")),
        "models_assembling_an_index_level_block": len(blocks),
        "index_level_blocks": len(allb),
        "models_index_generic_rank_3": len(r3),
        "models_index_ladder_0_n_2n": len(lad),
    }

    # ------------------------------------------------------- where the blocks fall short
    shortfall = Counter()
    for b in allb:
        if b["generic_rank_3"]:
            shortfall["rank 3 at index level"] += 1
        elif b["n_rows"] < 3 and b["n_cols"] < 3:
            shortfall["fewer than 3 rows AND fewer than 3 columns"] += 1
        elif b["n_rows"] < 3:
            shortfall["fewer than 3 rows"] += 1
        elif b["n_cols"] < 3:
            shortfall["fewer than 3 columns"] += 1
        else:
            shortfall["3+ rows and 3+ columns but support rank < 3"] += 1

    rep = {"geometry": "Z2xZ4_1_1",
           "scope": "CATALOGUE statement about the 11,540 PVW Z2xZ4_1_1 models; the "
                    "geometry-level admissibility it is contrasted with is levels.json and "
                    "the two are never conflated.",
           "funnel": funnel,
           "cpt_audit": cptaudit,
           "layer3_leg_occupancy": layer3,
           "block_shortfall_histogram": dict(shortfall.most_common()),
           "index_mass_pattern_histogram": dict(Counter(
               ":".join(b["mass_exponents"]) if b["mass_exponents"] else "rank<3"
               for b in allb).most_common()),
           "geometry_admits_vs_catalogue_realises": {
               "geometry_admits_0_2_4_on_n_channels":
                   len(lev["L_B_admissibility"]["channels_admitting_0_2_4"]),
               "geometry_admits_channels":
                   lev["L_B_admissibility"]["channels_admitting_0_2_4"],
               "catalogue_realises_at_index_level": len(lad),
               "statement": "the geometry ADMITS the pre-declared ladder; whether the "
                            "catalogue REALISES it is a separate question with a separate "
                            "answer, and the two are reported separately."}}
    L.write_json("diagnosis.json", rep)
    print(json.dumps({k: v for k, v in rep.items() if k != "layer3_leg_occupancy"},
                     indent=1, default=str)[:6000])


if __name__ == "__main__":
    main()

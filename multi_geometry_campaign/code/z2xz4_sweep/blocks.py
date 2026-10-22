"""
blocks.py --- STEP 2, catalogue half: the CATALOGUE-WIDE INDEX-LEVEL BLOCK SCAN.

COMPLETENESS BY SEARCH, NOT BY SCOPING (guard 3 -- the S1-F lesson of gate_z4xz4_sweep,
where a funnel that searched 12 of 51 channels and 1 of 6 leg assignments missed 24 blocks
carrying the pre-declared pattern).  This scan runs, for EVERY one of the 11,540 models:

    ALL 15 certified channels  x  ALL SIX assignments of the three legs to (Q, u^c, H_u)
    x  EVERY H_u-role field the model puts on the Higgs leg.

Nothing is pre-filtered.  The `interesting` flag computed at the end may ORDER the coupling
runs; it never scopes a discovery claim, and the coupling stage additionally couples a
fixed-seed RANDOM CONTROL sample drawn from the whole catalogue.

WHAT A BLOCK IS.  Rows = the model's own ROW-species fields in the row leg's sector, columns
= its COLUMN-species fields in the column leg's sector, at ONE fixed Higgs FIELD.  Two fields
at the SAME generation index are two independent rows: they carry the same classical exponent
and independent O(1) coefficients.  The exponent of a cell is read off the phase's
NC-CORRECTED index table (guard 1); a cell absent from that table is a structural zero, not a
small coefficient.

The verdict per block is the BASIS-INDEPENDENT mass ladder d_k (minimal exponent sum over
k x k minors with a generic non-vanishing determinant), so no statement depends on a
row/column labelling.

Emits: data/z2xz4_sweep/models/blocks.json.gz, data/z2xz4_sweep/blocks_summary.json,
       data/z2xz4_sweep/coupling_model_list.txt
"""

from __future__ import annotations

import gzip
import json
import os
import random
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import permutations

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z24sw_lib as L                                                   # noqa: E402
from levels import minors33, term_rank                                  # noqa: E402

ROW_SPECIES = "Q"
COL_SPECIES = "uc"
HIGGS = "Hu"
CONTROL_N = 600
CONTROL_SEED = 20260819


def load_tables():
    """{channel: {(idxQ, idxU, idxH) keyed by leg position: total exponent}} from the phase's
    NC-corrected index table."""
    with gzip.open(os.path.join(L.DSW, "table_z24.json.gz"), "rt", encoding="utf-8") as fh:
        tab = json.load(fh)
    out = {}
    for nm, v in tab.items():
        idx = {}
        for e in v["entries"]:
            key = tuple(tuple(None if x is None else int(x) for x in leg)
                        for leg in e["legs"])
            idx[key] = F(e["total"])
        out[nm] = {"gs": [tuple(x) for x in v["gs"]], "idx": idx,
                   "capable_planes": v["capable_planes"],
                   "nc_planes": v["nc_planes"]}
    return out


def minors_of(cells, rows, cols):
    """(d1, d2, d3) for an arbitrary rows x cols support, by min-cost matching over the
    ACTUAL fields (no occupancy freedom -- this is the model's own block)."""
    from levels import minors_by_flow
    return minors_by_flow(cells, rows, cols, 1, 1, 3)


def analyse_model(rec, tables):
    """Every index-level block this model carries, over all channels and assignments."""
    by_role = defaultdict(list)
    for k, l, o1, o2, o3, role, mult, hs, flab in rec["fields"]:
        by_role[role].append({"kl": (k, l), "idx": (o1, o2, o3), "field": flab,
                              "mult": mult})
    if not (by_role[ROW_SPECIES] and by_role[COL_SPECIES] and by_role[HIGGS]):
        return []
    bysec = {r: defaultdict(list) for r in (ROW_SPECIES, COL_SPECIES, HIGGS)}
    for r in (ROW_SPECIES, COL_SPECIES, HIGGS):
        for f in by_role[r]:
            bysec[r][f["kl"]].append(f)

    out = []
    for nm, T in tables.items():
        gs = T["gs"]
        if not T["capable_planes"]:
            continue                        # no letter: cannot carry a graded block
        seen_assign = set()
        for legQ, legU, legH in permutations(range(3)):
            key = (gs[legQ], gs[legU], gs[legH])
            if key in seen_assign:
                continue
            seen_assign.add(key)
            sQ, sU, sH = gs[legQ], gs[legU], gs[legH]
            if not (bysec[ROW_SPECIES].get(sQ) and bysec[COL_SPECIES].get(sU)
                    and bysec[HIGGS].get(sH)):
                continue
            for hf in bysec[HIGGS][sH]:
                hidx = hf["idx"]
                cells, rowf, colf = {}, {}, {}
                for qf in bysec[ROW_SPECIES][sQ]:
                    for uf in bysec[COL_SPECIES][sU]:
                        k = [None, None, None]
                        k[legQ] = qf["idx"]
                        k[legU] = uf["idx"]
                        k[legH] = hidx
                        val = T["idx"].get(tuple(k))
                        if val is None:
                            continue
                        cells[(qf["field"], uf["field"])] = val
                        rowf[qf["field"]] = qf
                        colf[uf["field"]] = uf
                if len(rowf) < 2 or len(colf) < 2:
                    continue
                rws = sorted(rowf)
                cls = sorted(colf)
                sup = {k: True for k in cells}
                d = minors_of(cells, rws, cls)
                steps = None
                if d[2] is not None:
                    steps = [d[0], d[1] - d[0], d[2] - d[1]]
                colmap, rowmap = defaultdict(set), defaultdict(set)
                for (r, c), v in cells.items():
                    colmap[c].add(v)
                    rowmap[r].add(v)
                col_det = all(len(s) == 1 for s in colmap.values())
                row_det = all(len(s) == 1 for s in rowmap.values())
                charges = sorted({next(iter(s)) for s in colmap.values()}) if col_det else []
                rcharges = sorted({next(iter(s)) for s in rowmap.values()}) if row_det else []
                out.append({
                    "channel": nm, "legQ": legQ, "legU": legU, "legH": legH,
                    "sector_Q": list(sQ), "sector_U": list(sU), "sector_H": list(sH),
                    "higgs_field": hf["field"], "higgs_index": list(hidx),
                    "row_fields": {rws[i]: list(rowf[rws[i]]["idx"]) for i in range(len(rws))},
                    "col_fields": {cls[i]: list(colf[cls[i]]["idx"]) for i in range(len(cls))},
                    "n_rows": len(rws), "n_cols": len(cls), "n_cells": len(cells),
                    "n_cells_full": len(rws) * len(cls),
                    "support_rank": term_rank(rws, cls, sup),
                    "column_determined": col_det,
                    "column_charges": [str(x) for x in charges],
                    "row_determined": row_det,
                    "row_charges": [str(x) for x in rcharges],
                    "minors_d": [str(x) if x is not None else "inf" for x in d],
                    "mass_exponents": [str(x) for x in steps] if steps else None,
                    "generic_rank_3": steps is not None,
                    "ladder_0_n_2n": bool(steps and steps[0] == 0 and steps[1] > 0
                                          and steps[2] == 2 * steps[1]),
                    "ladder_0_2_4": bool(steps and [int(x) for x in steps] == [0, 2, 4]),
                    "cells": {"%s|%s" % k: str(v) for k, v in sorted(cells.items())},
                })
    return out


def main():
    tables = load_tables()
    with gzip.open(os.path.join(L.DSW, "models", "roles.json.gz"), "rt",
                   encoding="utf-8") as fh:
        roles = json.load(fh)
    print("models with a census record:", len(roles), flush=True)

    allb, per_model = {}, {}
    for lab, rec in roles.items():
        if rec.get("status") != "OK":
            continue
        bl = analyse_model(rec, tables)
        if bl:
            allb[lab] = {"label": lab, "blocks": bl}
        per_model[lab] = len(bl)

    blocks = [b for v in allb.values() for b in v["blocks"]]
    ok = [r for r in roles.values() if r.get("status") == "OK"]
    lad = sorted({lab for lab, v in allb.items()
                  if any(b["ladder_0_n_2n"] for b in v["blocks"])})
    r3 = sorted({lab for lab, v in allb.items()
                 if any(b["generic_rank_3"] for b in v["blocks"])})
    coldet3 = sorted({lab for lab, v in allb.items()
                      if any(b["column_determined"] and len(b["column_charges"]) >= 3
                             for b in v["blocks"])})
    rowdet3 = sorted({lab for lab, v in allb.items()
                      if any(b["row_determined"] and len(b["row_charges"]) >= 3
                             for b in v["blocks"])})

    rep = {"geometry": "Z2xZ4_1_1",
           "scope": "CATALOGUE statement about the 11,540 PVW Z2xZ4_1_1 models",
           "search_space": "ALL 15 channels x ALL 6 leg assignments x EVERY H_u field "
                           "(completeness by search, guard 3 -- no pre-filter scopes any "
                           "discovery claim here)",
           "n_models_evaluable": len(ok),
           "n_models_with_any_block": len(allb),
           "n_blocks": len(blocks),
           "n_models_index_rank_3": len(r3),
           "n_models_index_ladder_0_n_2n": len(lad),
           "n_models_index_ladder_0_2_4": len(sorted({lab for lab, v in allb.items()
                                                      if any(b["ladder_0_2_4"]
                                                             for b in v["blocks"])})),
           "n_models_column_determined_3_charges": len(coldet3),
           "n_models_row_determined_3_charges": len(rowdet3),
           "mass_pattern_histogram": dict(Counter(
               ":".join(b["mass_exponents"]) if b["mass_exponents"] else "rank<3"
               for b in blocks).most_common()),
           "support_rank_histogram": dict(Counter(b["support_rank"] for b in blocks)),
           "channel_histogram": dict(Counter(b["channel"] for b in blocks).most_common()),
           "n_rows_histogram": dict(Counter(b["n_rows"] for b in blocks)),
           "n_cols_histogram": dict(Counter(b["n_cols"] for b in blocks)),
           "ladder_models": lad,
           "index_rank3_models": r3[:400],
           "column_determined_3_charge_models": coldet3[:400],
           "row_determined_3_charge_models": rowdet3[:400]}
    L.write_json("blocks_summary.json", rep)
    with gzip.open(os.path.join(L.DSW, "models", "blocks.json.gz"), "wt",
                   encoding="utf-8") as fh:
        json.dump(allb, fh, default=str)

    # ---- the coupling model list.
    # EVERY model that assembles ANY index-level block is coupled, not merely those whose
    # index-level block already looks like a ladder.  The reason is directional and was
    # checked, not assumed: orbifolder's realised support is a SUBSET of the index-level
    # support, and DELETING cells can RAISE a minor -- so a coupling-level ladder can appear
    # from an index-level non-ladder (e.g. 0:0:4 -> 0:2:4 if the second unsuppressed cell is
    # not generated).  Restricting the coupling runs to index-level ladder candidates would
    # therefore scope the discovery claim, which guard 3 forbids.  A fixed-seed RANDOM
    # CONTROL from the rest of the catalogue is coupled alongside.
    interesting = sorted(set(allb) | set(lad) | set(r3) | set(coldet3) | set(rowdet3))
    rnd = random.Random(CONTROL_SEED)
    pool = sorted(lab for lab, r in roles.items() if r.get("status") == "OK"
                  and lab not in interesting)
    control = rnd.sample(pool, min(CONTROL_N, len(pool)))
    lst = sorted(set(interesting) | set(control), key=lambda s: int(s[4:]))
    with open(os.path.join(L.DSW, "coupling_model_list.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lst) + "\n")
    rep2 = {"n_interesting": len(interesting), "n_control": len(control),
            "n_to_couple": len(lst), "control_seed": CONTROL_SEED}
    L.write_json("coupling_list_summary.json", rep2)
    print(json.dumps({k: v for k, v in rep.items()
                      if k not in ("ladder_models", "index_rank3_models",
                                   "column_determined_3_charge_models",
                                   "row_determined_3_charge_models")},
                     indent=1, default=str))
    print(json.dumps(rep2, indent=1))


if __name__ == "__main__":
    main()

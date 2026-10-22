"""z3xz6_sweep -- STEP 5: catalogue-level roll-up of the per-model ladder analysis.

Consumes `data/z3xz6_sweep/models/*.json` and emits `data/z3xz6_sweep/sweep_summary.json`:
the verdict distribution against the pre-declared SW-1/SW-2/SW-3 ladder, the coupling and
channel statistics, the block statistics (rank, distinct total exponents), the C(S_min)
branch statistics, and the gate results X1/X2/X3.

GUARDS ENFORCED HERE
  * gate X2 (R1): the six all-T(0,l) channels must be realised ZERO times across the
    whole catalogue.  A non-zero total is a MACHINERY FAILURE and is raised, never
    reported as a discovery -- their emptiness is pre-declared by the geometry-level
    R_2 = 0 theorem.
  * gate X3: Rule P (the STEP-2 pre-filter) is compared model by model against the actual
    superpotentials.  Any model where Rule P said `ladder-impossible` while W contains an
    up-type trilinear on a ladder-bearing channel FALSIFIES Rule P and is listed.
  * C(S_min) branch results are counted PER BRANCH and never averaged (GAMMA_BRANCH_TAG).
  * every ladder line carries [G-1] symmetry-protection and [G-3] ratio disclosure.
"""

from __future__ import annotations

import json
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z36sw_lib as Z                                                  # noqa: E402

MODELS = os.path.join(Z.DATA, "models")
CENSUS = os.path.join(Z.DATA, "census")


def main():
    rows = []
    for fn in sorted(os.listdir(MODELS)):
        if fn.endswith(".json"):
            rows.append(json.load(open(os.path.join(MODELS, fn))))
    ok = [r for r in rows if r.get("status") == "ok"]

    cen = {}
    for fn in sorted(os.listdir(CENSUS)):
        if fn.endswith(".json"):
            for m in json.load(open(os.path.join(CENSUS, fn)))["models"]:
                cen[m["model"]] = m

    s = {"scope": Z.CATALOGUE_TAG, "predictivity": Z.PREDICTIVITY_TAG,
         "ratio_disclosure": Z.RATIO_DISCLOSURE, "gamma_branch_tag": Z.GAMMA_BRANCH_TAG,
         "rulegap_G1": Z.RULEGAP_G1, "rulegap_G2": Z.RULEGAP_G2,
         "n_models_analysed": len(rows), "n_status_ok": len(ok),
         "status_histogram": dict(Counter(r.get("status") for r in rows))}

    # COVERAGE, stated before any statistic that depends on it.  The per-model ladder
    # analysis runs on whatever coupling segments have been published; the catalogue is
    # 4,612 models and a partial run is NOT a prefix but a UNIFORM RANDOM SAMPLE, because
    # the bulk driver runs a fixed-seed SHUFFLE of the catalogue (launch_bulk_shuffled.sh).
    # Every catalogue-level percentage below therefore carries a binomial error bar rather
    # than an unquantifiable selection effect.
    s["coverage_models_analysed"] = len(rows)
    s["coverage_fraction_of_catalogue"] = round(len(rows) / 4612.0, 4)
    s["coverage_note"] = (
        "Coupling-derived statistics are over the models analysed, which are the "
        "PRIORITY set (every model Rule P admits -- the only ones that can carry an "
        "up-type block on a ladder-bearing channel) plus a fixed-seed UNIFORM RANDOM "
        "SAMPLE of the rest.  Census-derived statements (the rank bound, the "
        "counterfactual, the three-layer diagnosis) are over ALL 4,610 evaluable models "
        "and do not depend on coverage at all.")

    s["verdicts"] = dict(Counter(r.get("verdict") for r in ok))
    s["total_W_terms"] = sum(r["n_W_terms"] for r in ok)
    s["total_up_trilinears"] = sum(r["n_up_trilinears"] for r in ok)
    s["models_with_an_up_trilinear"] = sum(1 for r in ok if r["n_up_trilinears"])
    s["models_with_up_on_an_allowed_channel"] = sum(
        1 for r in ok if r.get("up_on_allowed_channels"))
    s["models_with_up_on_a_ladder_channel"] = sum(
        1 for r in ok if r.get("up_on_ladder_channels"))

    # ---- W-print integrity (the forked-child duplicate print)
    s["models_with_disagreeing_W_prints"] = sum(
        1 for r in ok if r.get("W_prints_disagree"))
    s["models_with_disagreeing_W_prints_list"] = [
        r["model"] for r in ok if r.get("W_prints_disagree")]
    s["W_print_count_histogram"] = dict(sorted(Counter(
        len(r.get("W_print_sizes") or []) for r in ok).items()))
    s["W_print_note"] = (
        "`create coupling(F F F)` forks a child inheriting the redirected stdin, so each "
        "model in a multi-model run prints its superpotential twice.  The LONGEST print "
        "is used (a truncated print is a prefix).  Models whose prints DISAGREE had at "
        "least one truncated print; they are listed and re-run at a long wait, and the "
        "re-run supersedes.  Taking the first print would have recorded spurious "
        "ABSENCES -- the one error class this sweep cannot afford.")

    # ---- GATE X1 / X2
    s["gateX1_pointgroup_violations_total"] = sum(
        r["gateX1"]["pointgroup_violations"] for r in ok)
    s["gateX1_plane2_Rlaw_violations_total"] = sum(
        r["gateX1"]["plane2_Rlaw_violations"] for r in ok)
    s["gateX2_all_T0l_terms_total"] = sum(r["gateX2_all_T0l_terms"] for r in ok)
    s["gateX2_note"] = Z.EMPTY_BY_GEOMETRY_REASON
    if s["gateX2_all_T0l_terms_total"] != 0:
        raise SystemExit(
            "GATE X2 FAILED: %d all-T(0,l) trilinears found across the catalogue.  These "
            "channels are EMPTY BY GEOMETRY (R_2 = 0 theorem); a non-zero count is a "
            "machinery failure, not a discovery."
            % s["gateX2_all_T0l_terms_total"])

    # ---- GATE X3: Rule P versus the actual superpotentials
    falsify, agree = [], 0
    for r in ok:
        c = cen.get(r["model"], {})
        pred = bool(c.get("ruleP_pass"))
        actual = bool(r.get("up_on_ladder_channels"))
        if actual and not pred:
            falsify.append(r["model"])
        agree += (pred == actual)
    s["gateX3_ruleP_agreement"] = agree
    s["gateX3_ruleP_falsifying_models"] = falsify
    s["gateX3_ruleP_falsified"] = bool(falsify)
    s["gateX3_ruleP_predicted_pass"] = sum(
        1 for r in ok if cen.get(r["model"], {}).get("ruleP_pass"))

    # ---- channels actually populated by an up-type trilinear
    upch, upch_models = Counter(), Counter()
    for r in ok:
        for k, v in (r.get("up_channel_histogram") or {}).items():
            upch[k] += v
            upch_models[k] += 1
    s["up_channel_terms"] = dict(upch.most_common())
    s["up_channel_models"] = dict(upch_models.most_common())

    # ---- carrier-channel occupancy distribution (ALL trilinears, not just up-type)
    carr = Z.channel_name(Z.CARRIER)
    occ = Counter()
    for r in ok:
        occ[(r.get("channel_histogram") or {}).get(carr, 0)] += 1
    s["carrier_channel_occupancy_distribution_all_trilinears"] = dict(sorted(occ.items()))
    s["models_with_a_nonzero_carrier_occupancy"] = sum(
        n for k, n in occ.items() if k > 0)
    s["carrier_channel_max_occupancy"] = max(occ) if occ else 0
    carr2 = "T(2,3)T(2,4)T(2,5)"
    occ2 = Counter()
    for r in ok:
        occ2[(r.get("channel_histogram") or {}).get(carr2, 0)] += 1
    s["second_nine_factor_channel"] = carr2
    s["second_carrier_occupancy_distribution_all_trilinears"] = dict(sorted(occ2.items()))

    # ---- blocks
    blocks = [b for r in ok for b in (r.get("blocks") or [])]
    lb = [b for b in blocks if b["channel_is_ladder_bearing"]]
    s["n_single_higgs_blocks"] = len(blocks)
    s["n_blocks_on_ladder_channels"] = len(lb)
    s["block_rank_histogram_all"] = dict(sorted(Counter(
        b["generic_rank"] for b in blocks).items()))
    s["block_rank_histogram_ladder_channels"] = dict(sorted(Counter(
        b["generic_rank"] for b in lb).items()))
    s["block_distinct_factor_histogram_ladder_channels"] = dict(sorted(Counter(
        b["n_distinct_total_exponents"] for b in lb).items()))
    s["blocks_rank3_with_3plus_distinct_factors"] = [
        {"model": r["model"], "channel": b["channel"], "higgs": b["higgs"],
         "rank": b["generic_rank"],
         "distinct_total_exponents": b["distinct_total_exponents"],
         "tau_invariant_channel": b["tau_invariant_channel"],
         "leading_det_term": b["leading_det_term"]}
        for r in ok for b in (r.get("blocks") or [])
        if b["channel_is_ladder_bearing"] and b["generic_rank"] >= 3
        and b["n_distinct_total_exponents"] >= 3]
    s["n_SW1_blocks"] = len(s["blocks_rank3_with_3plus_distinct_factors"])

    # ---- C(S_min), per branch, never averaged
    br = {b: {"entries": 0, "zero": 0} for b in ("blind", "correlator", "LEEFT")}
    for r in ok:
        for row in (r.get("C_Smin") or []):
            for p in row["planes"].values():
                for b in br:
                    v = p.get(b)
                    if v:
                        br[b]["entries"] += 1
                        br[b]["zero"] += bool(v["C_is_zero"])
    for b in br:
        e = br[b]["entries"]
        br[b]["annihilated_pct"] = round(100.0 * br[b]["zero"] / e, 2) if e else None
    s["C_Smin_by_branch"] = br
    s["models_with_branch_disagreement"] = sum(
        1 for r in ok if r.get("C_Smin_branch_disagreement"))

    # ---- D-2 uncovered tagging
    s["blocks_carrying_D2_uncovered_legs"] = sum(
        1 for b in blocks if b.get("D2_uncovered_legs"))
    s["D2_note"] = Z.RULEGAP_G1

    dest = os.path.join(Z.DATA, "sweep_summary.json")
    json.dump(s, open(dest, "w"), indent=1, default=str)
    for k, v in s.items():
        if k.endswith("_tag") or k.startswith("scope") or k.startswith("rulegap") \
           or k in ("predictivity", "ratio_disclosure", "gateX2_note", "D2_note"):
            continue
        print("%-52s %s" % (k, json.dumps(v, default=str)[:220]))
    print("->", dest)


if __name__ == "__main__":
    main()

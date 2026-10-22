"""z3xz6_sweep -- guard R3: the gamma PAIRS actually carried by candidate ladder legs.

`C(S_min)` is weighted by `exp(-2 pi i * sgn * sum_alpha gamma_alpha(h_alpha))` with
`gamma_alpha(h) = k gamma^theta_alpha + l gamma^omega_alpha` the CHARACTER value of the
leg's stabiliser -- read off orbifolder's printed `gamma phases` slots 0 and 1, never
modelled.  The three-branch tag (blind / correlator `-gamma` / LEEFT `+gamma`) travels
with every `C(S_min)` statement.

But the branch question is only LIVE where the legs actually carry non-zero gamma: on
`MSSM0`'s 60 realisations of the carrier channel every leg has `(gamma^theta, gamma^omega)
= (0,0)`, so all three branches coincide THERE BY THE DATA -- a physical fact about that
model, not evidence that the machinery is branch-insensitive (selftest S9 settles the
latter separately, by injecting a character and watching `|C|` go 6.0 -> 0.0).

This file measures, over the models Rule P admits, the gamma pairs on their chiral
`Q` / `u^c` / `H_u` candidates, so the findings can say how load-bearing the branch
question is on THIS catalogue rather than assuming either way.

Run: python multi_geometry_campaign/code/z3xz6_sweep/gamma_on_legs.py
"""

from __future__ import annotations

import json
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z36sw_lib as Z                                                  # noqa: E402

CENSUS = os.path.join(Z.DATA, "census")


def main():
    want = set(json.load(open(os.path.join(Z.DATA, "census_summary.json")))
               ["ruleP_pass_models"])
    lc = Z.ladder_channels()
    LB_SECTORS = {tuple(t) for k, v in lc.items() if v["ladder_bearing"]
                  for t in v["tri"]}

    pairs = Counter()
    per_model = {}
    nonzero_models = 0
    for fn in sorted(os.listdir(Z.RAW_CENSUS)):
        if not fn.endswith(".gz") or "PILOT" in fn:
            continue
        txt = Z.read_raw(os.path.join(Z.RAW_CENSUS, fn))
        for lab, block in Z.split_models(txt):
            if lab not in want:
                continue
            spec = Z.parse_spectrum(block)
            gauge = Z.parse_gauge(block)
            pos = Z.observable_positions(spec, gauge["tY"], gauge) if gauge["tY"] else None
            if pos is None:
                continue
            a3, a2, sY0, tY, _bi = pos
            cen, _d = Z.role_census(spec, tY, a3, a2, sY0)
            loc = Counter()
            for l, v in cen.items():
                if v["role"] in ("Q", "uc", "Hu") and tuple(v["kl"]) in LB_SECTORS:
                    gp = Z.gamma_pair(spec[l])
                    key = "(%s,%s)" % (gp[0], gp[1])
                    pairs[key] += 1
                    loc[key] += 1
            per_model[lab] = dict(loc)
            if any(k != "(0,0)" for k in loc):
                nonzero_models += 1

    rep = {"scope": Z.CATALOGUE_TAG, "branch_tag": Z.GAMMA_BRANCH_TAG,
           "n_models": len(per_model),
           "gamma_pairs_on_ladder_sector_Q_uc_Hu": dict(pairs.most_common()),
           "models_with_a_nonzero_gamma_on_such_a_leg": nonzero_models,
           "per_model": per_model}
    rep["verdict"] = (
        "The three-branch gamma question is LIVE on %d of the %d Rule-P models: their "
        "candidate ladder legs carry a non-zero character pair, so the blind / correlator "
        "/ LEEFT branches can differ there and must be reported separately, never "
        "averaged." % (nonzero_models, len(per_model))
        if nonzero_models else
        "The three-branch gamma question is INERT on every Rule-P model: every candidate "
        "ladder leg carries (gamma^theta, gamma^omega) = (0,0), so the three branches "
        "coincide BY THE DATA.  That is a measured property of this catalogue's candidate "
        "legs, NOT a statement that the machinery is branch-insensitive (selftest S9 "
        "shows it is not) and NOT a resolution of RULE-GAP G-1.")
    dest = os.path.join(Z.DATA, "gamma_on_legs.json")
    json.dump(rep, open(dest, "w"), indent=1, default=str)
    print("gamma pairs on candidate ladder legs (Q / u^c / H_u), %d Rule-P models:"
          % len(per_model))
    for k, v in pairs.most_common():
        print("   %-12s %d" % (k, v))
    print("\nmodels with a NON-ZERO gamma on such a leg: %d / %d"
          % (nonzero_models, len(per_model)))
    print("\n" + rep["verdict"])
    print("->", dest)


if __name__ == "__main__":
    main()

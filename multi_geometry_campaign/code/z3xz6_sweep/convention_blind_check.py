"""z3xz6_sweep -- CONVENTION-BLIND check of the layer-2 absence (self-attack #3).

The headline rests on: *the sectors `T(0,3)` and `T(2,4)` host a chiral `Q` or `u^c` in NO
model of the catalogue*.  That statement is produced by the role census, which selects one
of four `(sY, sC)` sign branches by net `Q = +3` -- and 49% of the catalogue selects
`sC = -1`, so a branch slip would move exactly this kind of number.

This file re-derives the same absence WITHOUT any sign convention, using the
`gate_z6i_sweep/gza_census.py` pattern: `|c3|`, `|c2|` and `|Y|` are invariant under every
`(sY, sC)` choice, so

    a quark doublet in sector s  requires a field in s with  |c3| = 3 and |c2| = 2
    an up-type singlet in s      requires a field in s with  |c3| = 3, |c2| = 1, |Y| = 2/3

and no missed branch can hide one.  The check is deliberately GENEROUS: it ignores
chirality, ignores whether the field is a hidden-sector singlet, and ignores the
hypercharge for the doublet case.  An absence that survives a generous convention-blind
scan is far stronger than one selected by a branch.

Run: python multi_geometry_campaign/code/z3xz6_sweep/convention_blind_check.py
"""

from __future__ import annotations

import json
import os
import sys
from collections import Counter
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z36sw_lib as Z                                                  # noqa: E402


def main():
    tri_doublet = Counter()     # sector -> models with ANY |c3|=3, |c2|=2 field
    tri_singlet = Counter()     # sector -> models with ANY |c3|=3, |c2|=1, |Y|=2/3 field
    n_models = 0
    for fn in sorted(os.listdir(Z.RAW_CENSUS)):
        if not fn.endswith(".gz") or "PILOT" in fn:
            continue
        txt = Z.read_raw(os.path.join(Z.RAW_CENSUS, fn))
        for lab, block in Z.split_models(txt):
            spec = Z.parse_spectrum(block)
            gauge = Z.parse_gauge(block)
            if not gauge["tY"]:
                continue
            fac = Z.parse_banner(gauge["gauge_test"])
            i3 = [i for i, f in enumerate(fac) if f == ("SU", 3)]
            i2 = [i for i, f in enumerate(fac) if f == ("SU", 2)]
            if not i3 or not i2:
                continue
            n_models += 1
            d_here, s_here = set(), set()
            for r in spec.values():
                if not r["rep"] or not r["p_sh"]:
                    continue
                for a3 in i3:
                    if a3 >= len(r["rep"]) or abs(r["rep"][a3]) != 3:
                        continue
                    absY = abs(Z.dot(gauge["tY"], r["p_sh"]))
                    for a2 in i2:
                        if a2 >= len(r["rep"]):
                            continue
                        if abs(r["rep"][a2]) == 2:
                            d_here.add(str(r["kl"]))
                    if all(abs(r["rep"][a2]) == 1 for a2 in i2 if a2 < len(r["rep"])) \
                       and absY == F(2, 3):
                        s_here.add(str(r["kl"]))
            for s in d_here:
                tri_doublet[s] += 1
            for s in s_here:
                tri_singlet[s] += 1

    allsec = [str((k, l)) for k in range(3) for l in range(6)]
    rep = {"scope": Z.CATALOGUE_TAG, "n_models_scanned": n_models,
           "method": ("convention-blind: |c3|, |c2| and |Y| are invariant under both the "
                      "sY (overall Y sign) and sC (3 vs 3bar) conventions, so no missed "
                      "branch can hide a candidate.  Deliberately generous: chirality, "
                      "hidden-singlet-ness and (for doublets) hypercharge are IGNORED."),
           "models_with_a_colour_triplet_SU2_doublet_per_sector":
               {s: tri_doublet.get(s, 0) for s in allsec},
           "models_with_a_colour_triplet_SU2_singlet_absY_2_3_per_sector":
               {s: tri_singlet.get(s, 0) for s in allsec},
           "sectors_with_NO_triplet_doublet_anywhere":
               [s for s in allsec if tri_doublet.get(s, 0) == 0],
           "sectors_with_NO_triplet_singlet_absY_2_3_anywhere":
               [s for s in allsec if tri_singlet.get(s, 0) == 0],
           "right_chiral_labels": [str(g) for g in Z.RIGHT_CHIRAL_LABELS]}
    key = ["(0, 3)", "(2, 4)"]
    rep["carrier_leg_sectors_checked"] = {
        s: {"models_with_a_triplet_doublet": tri_doublet.get(s, 0),
            "models_with_a_triplet_singlet_absY_2_3": tri_singlet.get(s, 0)} for s in key}
    rep["verdict"] = (
        "CONFIRMED convention-blind: %s host NO colour-triplet SU(2)-doublet and NO "
        "colour-triplet SU(2)-singlet with |Y| = 2/3 in ANY of the %d models scanned, so "
        "the layer-2 absence that kills the carrier channel does NOT depend on the "
        "(sY, sC) branch selection."
        % (" and ".join(key), n_models)
        if all(tri_doublet.get(s, 0) == 0 and tri_singlet.get(s, 0) == 0 for s in key)
        else "REFUTED -- a candidate exists convention-blind; the layer-2 claim must be "
             "re-derived.")
    dest = os.path.join(Z.DATA, "convention_blind_check.json")
    json.dump(rep, open(dest, "w"), indent=1, default=str)
    print("convention-blind scan over %d models" % n_models)
    print("  models with a colour-triplet x SU(2)-doublet, per sector:")
    for s in allsec:
        print("     %-8s %d" % (s, tri_doublet.get(s, 0)))
    print("  sectors with NO triplet-doublet anywhere:",
          rep["sectors_with_NO_triplet_doublet_anywhere"])
    print("  sectors with NO triplet-singlet |Y|=2/3 anywhere:",
          rep["sectors_with_NO_triplet_singlet_absY_2_3_anywhere"])
    print("\n" + rep["verdict"])
    print("->", dest)


if __name__ == "__main__":
    main()

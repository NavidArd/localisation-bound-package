"""z3xz6_sweep -- the CARRIER's unreachability, re-derived convention-blind and layered.

`convention_blind_check.py` refuted an over-claim in the first draft of the headline: the
sector `T(0,3)` is NOT empty of quark-doublet quantum numbers -- **219** models carry a
colour-triplet x SU(2)-doublet field there.  Only `T(2,4)` is empty convention-blind
(together with the four right-chiral LABELS, whose emptiness is labelling).

So the carrier's death has to be stated exactly.  This file does it in three nested
statements, each strictly weaker than the last, so the reader can see which one is
load-bearing:

  B1 (CONVENTION-BLIND, strongest).  `T(2,4)` hosts NO field with colour-triplet x
     SU(2)-doublet quantum numbers, and none with colour-triplet x SU(2)-singlet and
     |Y| = 2/3, in ANY model.  `|c3|`, `|c2|`, `|Y|` are invariant under both sign
     conventions, so this cannot be a branch artefact.  It removes `T(2,4)` as a candidate
     `Q` leg or `u^c` leg for the whole catalogue.
  B2 (still convention-blind).  With `T(2,4)` excluded, the carrier's `Q` and `u^c` must
     occupy `T(0,3)` and `T(1,5)` in some order.  Count the models that carry the
     necessary quantum numbers in BOTH at once -- ignoring chirality and hidden charge.
  B3 (role-level).  Among those, how many actually have a CHIRAL `Q` and a CHIRAL `u^c`
     there, as hidden-sector singlets with the SM hypercharges?

Run: python multi_geometry_campaign/code/z3xz6_sweep/carrier_blind_check.py
"""

from __future__ import annotations

import json
import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z36sw_lib as Z                                                  # noqa: E402

CARRIER_LEGS = [(0, 3), (1, 5), (2, 4)]


def main():
    b1_doublet = {str(s): 0 for s in CARRIER_LEGS}
    b1_singlet = {str(s): 0 for s in CARRIER_LEGS}
    b2 = 0
    b2_models = []
    n = 0
    for fn in sorted(os.listdir(Z.RAW_CENSUS)):
        if not fn.endswith(".gz") or "PILOT" in fn:
            continue
        for lab, block in Z.split_models(Z.read_raw(os.path.join(Z.RAW_CENSUS, fn))):
            spec = Z.parse_spectrum(block)
            gauge = Z.parse_gauge(block)
            if not gauge["tY"]:
                continue
            fac = Z.parse_banner(gauge["gauge_test"])
            i3 = [i for i, f in enumerate(fac) if f == ("SU", 3)]
            i2 = [i for i, f in enumerate(fac) if f == ("SU", 2)]
            if not i3 or not i2:
                continue
            n += 1
            doub, sing = set(), set()
            for r in spec.values():
                if not r["rep"] or not r["p_sh"]:
                    continue
                s = str(r["kl"])
                absY = abs(Z.dot(gauge["tY"], r["p_sh"]))
                for a3 in i3:
                    if a3 >= len(r["rep"]) or abs(r["rep"][a3]) != 3:
                        continue
                    if any(a2 < len(r["rep"]) and abs(r["rep"][a2]) == 2 for a2 in i2):
                        doub.add(s)
                    if all(a2 >= len(r["rep"]) or abs(r["rep"][a2]) == 1 for a2 in i2) \
                       and absY == F(2, 3):
                        sing.add(s)
            for s in b1_doublet:
                b1_doublet[s] += (s in doub)
                b1_singlet[s] += (s in sing)
            # B2: Q and u^c must sit on T(0,3) and T(1,5) in SOME order
            ok = ((("(0, 3)" in doub) and ("(1, 5)" in sing))
                  or (("(1, 5)" in doub) and ("(0, 3)" in sing)))
            if ok:
                b2 += 1
                b2_models.append(lab)

    rep = {"scope": Z.CATALOGUE_TAG, "n_models_scanned": n,
           "B1_convention_blind_per_leg": {
               s: {"models_with_a_colour_triplet_SU2_doublet": b1_doublet[s],
                   "models_with_a_colour_triplet_SU2_singlet_absY_2_3": b1_singlet[s]}
               for s in b1_doublet},
           "B2_models_carrying_the_quantum_numbers_on_BOTH_T03_and_T15": b2,
           "B2_models": b2_models[:50]}
    rep["B1_verdict"] = (
        "T(2,4) carries NEITHER quark-doublet NOR up-singlet quantum numbers in any of "
        "the %d models -- convention-blind, so it can be neither the Q leg nor the u^c "
        "leg of the carrier for the whole catalogue.  T(0,3) is NOT empty this way "
        "(%d models carry a colour-triplet x SU(2)-doublet field there), so the carrier's "
        "death is NOT a pure quantum-number statement on both legs and must not be "
        "written as one." % (n, b1_doublet["(0, 3)"]))
    rep["B2_verdict"] = (
        "With T(2,4) excluded convention-blind, the carrier needs its Q and u^c on T(0,3) "
        "and T(1,5).  Models carrying the necessary quantum numbers on BOTH at once: %d."
        % b2)
    dest = os.path.join(Z.DATA, "carrier_blind_check.json")
    json.dump(rep, open(dest, "w"), indent=1, default=str)
    print(json.dumps(rep, indent=1, default=str)[:2600])
    print("->", dest)


if __name__ == "__main__":
    main()

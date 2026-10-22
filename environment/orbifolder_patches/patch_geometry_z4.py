"""R-symmetry patches for the Z4_1_1 geometry file (task item 3).

`Geometry_Z4_1_1.txt` declares FIVE non-R symmetries and ZERO R symmetries (verified
by z6ilib.parse_geometry), so orbifolder's R-charge loop in `CreatePreCouplings` is
inert on every stock Z4_1_1 run -- CLAUDE.md sec.5.16 instance #2 ("right file, real
code, wrong input data"), distribution-wide per the R-charge census generalisation.

WHAT IS EXPRESSIBLE ON Z4 -- MORE THAN ON EITHER Z6 GEOMETRY.

  The declaration mechanism is `GetDiscreteCharge = -R^paper` (orbifolder
  `cfield.cpp:118-128` with `OsciContribution = +N_L - Nbar_L`,
  `cfixedpoint.cpp:267-284`, and `q_sh = -q^paper`; both signs flip together, so the
  mechanism is exact INCLUDING oscillator states -- the source reading verified in
  the Z6-I foundations, `docs/rules_z6i.md` sec.6.1).  A geometry-file R symmetry can
  express any GAMMA-FREE law of the form sum_a (charge-op . R_a) = target mod order.

  Z4's printed rule set (1301.2322 (3.49)-(3.52)) contains THREE gamma-free laws:

    (3.50)  plane 3:  sum R^3 = 1 mod 2     [the prime-plane law]
    (3.51)  plane 1:  sum R^1 = 1 mod 2     [the (theta_1)^2 extra law]
    (3.52)  plane 2:  sum R^2 = 1 mod 2     [the (theta_2)^2 extra law]

  and ONE gamma-carrying law, (3.49): sum R = 2 mod 4 with R_a carrying + 4 gamma_a.
  The three mod-2 laws are EXACT as declarations (no gamma anywhere in them);
  the mod-4 law's gamma term is NOT expressible (GetDiscreteCharge never touches
  gamma_phases), so its gamma-free truncation is KNOWN-TOO-STRONG: with
  sum(gamma) = 1/2 the truncation demands sum R = 2 mod 4 where the true law demands
  sum R = 0 mod 4, and vice versa -- it would kill every sum-gamma = 1/2 coupling.

  Contrast: Z6-II could express two of three laws exactly; Z6-I one of two.  Z4_1_1
  expresses THREE of its FOUR printed laws exactly -- the largest expressible
  fraction of any geometry in the campaign so far.  Only the mod-4 refinement (the
  law that distinguishes sum-gamma = 0 from 1/2, i.e. the one that carries the
  orbit-class information) needs the offline route (`code/z4/rules_z4.py`).

Emitted patches (each re-parsed after writing; engine re-loaded; theta, orbit counts
and orbit-size histograms asserted unchanged):

  Geometry_Z4_1_1_Rsym_exact.txt -- (3.50) + (3.51) + (3.52), all exact.  This is
      the patch a Z4 run should use, with (3.49) applied offline.
  Geometry_Z4_1_1_Rsym_full4.txt -- additionally the gamma-free truncation of
      (3.49) (order 4, two-slot operator (0,1,1,0), target 2).  KNOWN-TOO-STRONG;
      emitted only to quantify what the gamma term buys on this geometry.  It may
      not carry a physics verdict.

Reproduce:  python multi_geometry_campaign/code/z4/patch_geometry_z4.py
Emits:      data/z4/patched_geometry/*.txt
"""

import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "z6i"))
from z6ilib import Geom, parse_geometry  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(HERE))
PHYS = os.path.dirname(ROOT)
GEOMDIR = os.path.join(PHYS, "orbifolder-1.2.1", "Geometry")
OUT = os.path.join(ROOT, "data", "z4", "patched_geometry")


def rsym_block(order, slots, target):
    """One `new R symmetry` block.  slots: 1-based plane indices with weight 1 in the
    charge operator (slot 0 is the 4d helicity slot).  target: the residue the
    PAPER's law requires; the declaration carries -target because
    GetDiscreteCharge = -R^paper."""
    op = ["0", "0", "0", "0"]
    for p in slots:
        op[p] = "1"
    return "new R symmetry\n%d\n%d/1\n%s\n" % (order, -target, " ".join(op))


def main():
    os.makedirs(OUT, exist_ok=True)
    src_path = os.path.join(GEOMDIR, "Geometry_Z4_1_1.txt")
    g = Geom(src_path)
    d = g.d
    nR = [s for s in d["discrete_symmetries"] if s["kind"] == "R"]
    print("Z4_1_1: plane orders %s, declared R symmetries %d, non-R %d"
          % (g.plane_order, len(nR), len(d["discrete_symmetries"]) - len(nR)))
    assert not nR, "already declares R symmetries"

    # the three exact mod-2 laws: (3.51) plane 1, (3.52) plane 2, (3.50) plane 3
    blocks_exact = (rsym_block(2, [1], 1) + rsym_block(2, [2], 1)
                    + rsym_block(2, [3], 1))
    # the gamma-free truncation of (3.49): order 4, planes 1+2 combined, target 2
    tgt = sum(g.N * g.v[a] for a in g.nonprime_planes)
    assert tgt.denominator == 1 and int(tgt) == 2, tgt
    blocks_full4 = blocks_exact + rsym_block(g.N, [1, 2], int(tgt))

    src = open(src_path).read()
    for tag, blocks, n_expect in (("exact", blocks_exact, 3),
                                  ("full4", blocks_full4, 4)):
        new = src.replace("end discrete symmetries",
                          blocks + "end discrete symmetries")
        assert new != src
        path = os.path.join(OUT, "Geometry_Z4_1_1_Rsym_%s.txt" % tag)
        open(path, "w").write(new)
        chk = parse_geometry(path)
        got = [s for s in chk["discrete_symmetries"] if s["kind"] == "R"]
        assert len(got) == n_expect, (tag, got)
        # patched file must still load; theta and orbits unchanged
        tmp = os.path.join(GEOMDIR, "Geometry__z4_tmp.txt")
        shutil.copy(path, tmp)
        try:
            gg = Geom(tmp)
            assert gg.plane_order == g.plane_order
            assert all(gg.orbit_size_histogram(k) == g.orbit_size_histogram(k)
                       for k in range(1, g.N))
        finally:
            os.remove(tmp)
        print("   %-5s -> %s  (%d R symmetries: %s)"
              % (tag, os.path.basename(path), len(got),
                 [(s["order"], s["superpotential_charge"], s["rows"][0])
                  for s in got]))
    print()
    print("SCOPE: *_Rsym_exact is exact ((3.50)+(3.51)+(3.52), all gamma-free).")
    print("*_Rsym_full4 adds the gamma-free truncation of (3.49): KNOWN-TOO-STRONG,")
    print("it may not carry a physics verdict.  Rules 4/5/6 are ABSENT from")
    print("orbifolder entirely and are applied offline by code/z4/rules_z4.py.")


if __name__ == "__main__":
    main()

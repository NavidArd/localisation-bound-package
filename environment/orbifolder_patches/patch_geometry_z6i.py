"""R-symmetry patches for the Z6-I geometry files.

`Geometry_Z6-I_1_1.txt` declares TWO non-R symmetries and ZERO R symmetries (verified by
`z6ilib.parse_geometry`; gamma-gate guard 4 flagged this in advance and it is confirmed
here).  So orbifolder's R-charge loop in `CreatePreCouplings` is inert on every stock
Z6-I run -- CLAUDE.md sec.5.16 instance #2, "right file, real code, wrong input data".

WHAT IS EXPRESSIBLE HERE, AND WHY THE PRIME/NON-PRIME SPLIT RUNS THE OTHER WAY THAN Z6-II.

  Z6-II has ONE non-prime plane and TWO prime ones: two of its three laws are exactly
  expressible and only the third needs the offline route.  Z6-I is the mirror image:
  TWO non-prime planes served by ONE combined law, and only ONE prime plane.

  * plane 3 (prime, N^3 = 3), 1301.2322 (3.53):
        sum_a R^3_a = 1 mod 3 ,  R^3_a = q^3_sh a - N^3_L a + Nbar^3_L a
    EXPRESSIBLE, exactly, as an order-3 R symmetry with charge operator e_3.

  * planes 1 and 2 (both non-prime, N^i = 6), 1301.2322 (3.54):
        sum_a R_a = 2 mod 6 ,
        R_a = sum_{i=1,2} (q^i_sh a - N^i_L a + Nbar^i_L a) + 6 gamma_a
    The gamma-free part IS expressible -- and, unlike on Z6-II, it needs a charge
    operator with TWO non-zero slots, (0,1,1,0), because the law is a SINGLE condition
    summed over the two non-prime planes and NOT two per-plane conditions (the one-theta,
    one-gamma-slot discipline: gamma-gate guard 2).  The `+ 6 gamma_a` term is NOT
    expressible: `CField::GetDiscreteCharge` (cfield.cpp:118-128) builds the charge as
    `sum_k (q_sh[k] + OsciContribution[k]) * ChargeOperator[k]` and never touches
    `gamma_phases`.  Declaring it therefore applies the OLD, TOO-STRONG law.

SIGN, verified against the source rather than assumed (this corrects a prose slip in
`rules_z6ii.md` sec.3.1, flagged per CLAUDE.md sec.8.6 -- the conclusion there is right,
the stated reason is sign-inverted):
  `OsciContribution` (cfixedpoint.cpp:267-284) is `+1` per COMPLEX (unbarred) and `-1`
  per non-complex (barred) oscillator, i.e. `+N_L - Nbar_L`, NOT `-N_L + Nbar_L`.
  Together with orbifolder's `q_sh = -q^paper` this gives
        GetDiscreteCharge = -(q^paper - N_L + Nbar_L) = -R^paper ,
  so the declaration must carry superpotential charge  -(target).  Both signs flip
  together, which is why the mechanism is exact INCLUDING oscillator states -- a point
  that matters much more on Z6-I than on Z6-II, because three of Z6-I's six capacity
  channels (T1T1T4, T2T2T2, T2T5T5) admit no oscillator-free coupling at all.

Emitted patches (both re-parsed and checked after writing):
  *_Rsym_prime.txt  -- plane 3 only.  Everything it enforces is CORRECT; this is the
                       patch a Z6-I sweep should use, with (3.54) applied offline.
  *_Rsym_old.txt    -- plane 3 plus the gamma-free (3.54).  KNOWN-INCOMPLETE / TOO
                       STRONG; emitted only to quantify what the gamma term buys on this
                       geometry (the Z6-I analogue of Table 4's 96 -> 156).
"""

import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z6ilib import Geom, parse_geometry  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(HERE))))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
GEOMDIR = os.path.join(ROOT, "orbifolder-1.2.1", "Geometry")
OUT = os.path.join(ROOT, "multi_geometry_campaign", "data", "z6i", "patched_geometry")

LABELS = ["Z6-I_1_1", "Z6-I_2_1"]


def rsym_block(order, slots, target):
    """slots: list of 1-based plane indices carrying weight 1 in the charge operator.
    target: the residue the PAPER's law requires; the declaration carries -target."""
    op = ["0", "0", "0", "0"]
    for p in slots:
        op[p] = "1"
    return "new R symmetry\n%d\n%d/1\n%s\n" % (order, -target, " ".join(op))


def main():
    os.makedirs(OUT, exist_ok=True)
    for lab in LABELS:
        src_path = os.path.join(GEOMDIR, "Geometry_%s.txt" % lab)
        if not os.path.exists(src_path):
            print("skip (not shipped):", lab)
            continue
        try:
            g = Geom(src_path)
        except Exception as e:
            # Z6-I_2_1 is NOT factorizable (its lattice does not split into the three
            # complex planes).  1301.2322's Z6-I rule set (3.53)-(3.58) is derived for
            # the FACTORIZABLE lattice and does not apply there; the applicable law is
            # (3.65)/(3.66).  Recorded as a RULE-GAP, not patched.
            print("%s: NOT FACTORIZABLE (%s) -- RULE-GAP, no patch emitted" % (lab, e))
            continue
        d = g.d
        nR = [s for s in d["discrete_symmetries"] if s["kind"] == "R"]
        print("%s: plane orders %s, non-prime %s, prime %s, declared R symmetries %d, "
              "non-R %d" % (lab, g.plane_order, [a + 1 for a in g.nonprime_planes],
                            [a + 1 for a in g.prime_planes], len(nR),
                            len(d["discrete_symmetries"]) - len(nR)))
        if nR:
            print("   !! already declares R symmetries -- not patching")
            continue
        prime = [a + 1 for a in g.prime_planes]
        nonprime = [a + 1 for a in g.nonprime_planes]
        # targets from 1301.2322: prime plane -> 1 mod N^i;  non-prime combined ->
        # (sum_{i nonprime} N v^i) mod N
        blocks_prime = "".join(rsym_block(g.plane_order[p - 1], [p], 1) for p in prime)
        tgt = sum(g.N * g.v[a] for a in g.nonprime_planes)
        assert tgt.denominator == 1, tgt
        blocks_old = blocks_prime + rsym_block(g.N, nonprime, int(tgt))

        src = open(src_path).read()
        for tag, blocks in (("prime", blocks_prime), ("old", blocks_old)):
            new = src.replace("end discrete symmetries",
                              blocks + "end discrete symmetries")
            assert new != src
            path = os.path.join(OUT, "Geometry_%s_Rsym_%s.txt" % (lab, tag))
            open(path, "w").write(new)
            chk = parse_geometry(path)
            got = [s for s in chk["discrete_symmetries"] if s["kind"] == "R"]
            n_expect = len(prime) + (1 if tag == "old" else 0)
            assert len(got) == n_expect, (tag, got)
            # the patched file must still load in the engine (theta, orbits unchanged)
            tmp = os.path.join(GEOMDIR, "Geometry__z6i_tmp.txt")
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
    print("SCOPE: *_Rsym_prime is exact.  *_Rsym_old is TOO STRONG (it imposes (3.54)")
    print("without the +6 gamma term) and may not carry a physics verdict.  Neither")
    print("expresses Rule 4, Rule 5 or Rule 6 -- those are ABSENT from orbifolder")
    print("entirely and are applied offline by code/z6i/rules_z6i.py.")


if __name__ == "__main__":
    main()

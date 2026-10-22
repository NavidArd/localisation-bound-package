"""multi_geometry_campaign -- R-symmetry patches for the four Z6-II geometry files.

CLAUDE.md §5.16 instance #2 generalises: ALL 138 geometry files shipped with orbifolder
1.2.1 declare ZERO R symmetries, so `CreatePreCouplings`' R-charge loop is inert on every
stock run.  `rcharge_census` fixed this for Z3xZ3 by patching the geometry file to declare
the three order-3 R symmetries.  This script does the Z6-II analogue -- BUT the Z6-II case
is NOT a transcription of the Z3xZ3 one, and the difference is the whole point:

  * Planes 2 and 3 have PRIME rotation order (3 and 2).  arXiv:1301.2322 eqs. (3.60),
    (3.61) give, for these, exactly the traditional law with NO gamma term:
        sum_a (q^i_sh,a - N^i_L,a + Nbar^i_L,a) = 1 mod N^i ,  N^2 = 3, N^3 = 2.
    In orbifolder's sign convention that is  sum_a (q_sh + Osci) = -1 mod N^i, which IS
    expressible as an `SDiscreteSymmetry` (order N^i, superpotential charge -1, charge
    operator e_i).  These two patches are EXACT.

  * Plane 1 has NON-PRIME rotation order 6.  Its law, arXiv:1301.2322 eq. (3.59), is
        sum_a (q^1_sh,a - N^1_L,a + Nbar^1_L,a + 6 gamma_a) = 1 mod 6 ,
    and the `+6 gamma_a` term CANNOT be expressed as an orbifolder discrete symmetry:
    `CField::GetDiscreteCharge` (cfield.cpp:118-128) builds the R-charge as
    `sum_k (q_sh[k] + OsciContribution[k]) * ChargeOperator[k]` and has no access to
    `gamma_phases` at all.  So a plane-1 R-symmetry declaration would make orbifolder
    apply the OLD, TOO-STRONG rule.  On the published Z6-II benchmark of 1301.2322
    Table 4 that rule keeps 96 trilinears where the correct rule keeps 156 -- it kills
    60 couplings that genuinely exist.

Therefore TWO patches are emitted per geometry, and they are not interchangeable:

  *_Rsym_prime.txt   declares ONLY the two prime-plane R symmetries (orders 3 and 2).
                     Everything it enforces is CORRECT.  This is the patch a Z6-II
                     coupling sweep should use; the plane-1 law is then applied OFFLINE,
                     gamma-aware, by code/zrules.py.
  *_Rsym_old.txt     declares all three (orders 6, 3, 2), i.e. the traditional
                     pre-1301.2322 rule set.  It is KNOWN-INCOMPLETE and is emitted only
                     to reproduce the "old R-charges" column of Table 4 inside orbifolder.
                     DO NOT use it for a physics verdict.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from geomlib import Geometry, GEOMDIR, ROOT                          # noqa: E402

OUT = os.path.join(ROOT, "data", "patched_geometry")
LABELS = ["Z6-II_1_1", "Z6-II_2_1", "Z6-II_3_1", "Z6-II_4_1"]


def rsym_block(order, plane):
    op = ["0", "0", "0", "0"]
    op[plane] = "1"
    return "new R symmetry\n%d\n-1/1\n%s\n" % (order, " ".join(op))


def main():
    os.makedirs(OUT, exist_ok=True)
    for lab in LABELS:
        g = Geometry(lab)
        N = g.plane_orders
        assert N == (6, 3, 2), "%s: unexpected plane orders %s" % (lab, N)
        assert g.Rsym == [], "%s already declares R symmetries" % lab
        src = open(os.path.join(GEOMDIR, "Geometry_%s.txt" % lab)).read()
        for tag, planes in (("prime", [2, 3]), ("old", [1, 2, 3])):
            blocks = "".join(rsym_block(N[p - 1], p) for p in planes)
            new = src.replace("end discrete symmetries", blocks + "end discrete symmetries")
            assert new != src
            path = os.path.join(OUT, "Geometry_%s_Rsym_%s.txt" % (lab, tag))
            with open(path, "w") as fh:
                fh.write(new)
            # re-parse the patch and check it declares what we intended
            import shutil
            tmp = os.path.join(GEOMDIR, "Geometry__tmp_check.txt")
            shutil.copy(path, tmp)
            try:
                gg = Geometry("_tmp_check")
                got = [(r["order"], r["superpotential_charge"], r["operators"])
                       for r in gg.Rsym]
                assert len(got) == len(planes), got
                for (o, w, ops), p in zip(got, planes):
                    assert o == N[p - 1] and w == "-1" and ops[0][p] == "1", (o, w, ops, p)
            finally:
                os.remove(tmp)
            print("%-14s %-6s -> %s   (orders %s)" %
                  (lab, tag, os.path.basename(path), [N[p - 1] for p in planes]))


if __name__ == "__main__":
    main()

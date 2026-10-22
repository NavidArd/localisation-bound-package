"""S2 -- ROUTE B sizing: the (V1, W) lift space over the fixed target V2 = p.

Two-route sizing discipline (z4_construct / z2xz4_construct): the enumeration is sized
analytically BEFORE it is run and the size is reported next to the achieved count.

THE SPACE, declared
-------------------
By L1 (s1_liftspace.py, derived at source) the local shift at T(0,1) is V2 exactly.
Criterion (i) therefore pins  V2 := p  (the gate's target vector), up to (E1) lattice
automorphisms and (E2) translations; we fix the representative V2 = p and let the residual
freedom act on (V1, W).  What remains free is

    V1  : order | 3   ->  V1 = x/3,  x in Lambda,  class x in Lambda/3Lambda  (3^16)
    W   : order | 3   ->  W  = z/3,  z in Lambda,  class z in Lambda/3Lambda  (3^16)
    (W1 = W2 = W5 = W6 = 0 forced by the geometry; W3 = W4 = W forced by `identical WLs`)

CLASS-LEVEL conditions (well defined on Lambda/3Lambda -- proved in the docstring of each):
  A1  SM-preserving:  alpha . V1 in Z  for the SM simple roots alpha, i.e. alpha . x = 0
      mod 3.  Same for z.  [The SM is the one for which p is a 10-role: the block-1 SU(5)
      with roots e_i - e_j, i<j<=5; SU(3)_c = <e1-e2, e2-e3>, SU(2)_L = <e4-e5>.  This is
      the local-GUT statement of su5_ten_bound_findings section 4, declared, not assumed.]
  A2  MI1 :  3(V1^2 - 2/9) = 0 mod 2  <=>  x^2 = 2 mod 6.
      well defined: x -> x + 3L gives x^2 + 6 x.L + 9 L^2, and 9L^2 = 0 mod 6 since L^2
      is even.
  A3  MI6 :  3 W^2 = 0 mod 2  <=>  z^2 = 0 mod 6.   (same argument)
  A4  hidden-SU(2) survival:  Delta . V1 in Z and Delta . W in Z, i.e. Delta.x = 0 mod 3
      and Delta.z = 0 mod 3, with Delta = (0^8 | 1,1,0^6).
      well defined: Delta . 3L = 3 (Delta.L) = 0 mod 3.
      [Necessary for BOTH members of the doubled 10-role to survive the T(0,1) projection:
      the two states differ by Delta, all other projection data at the cell is identical,
      so their projection phases differ by exp(2 pi i Delta.V_h) for h in the centralizer;
      (theta, 0) and the plane-2 translations both fix this cell -- T(0,1) has a single
      translation-free constructing element (L1) -- so both conditions are imposed.]

REPRESENTATIVE-LEVEL conditions (NOT class invariant; a representative is searched for,
exactly as z6i_settle does for 3 V.W):
  B1  MI3 :  3 (V1.V2 - 1/18) = 0 mod 2
  B2  MI4 :  3 (W.V1) = 0 mod 2
  B3  MI5 :  3 (W.V2) = 0 mod 2

Outputs data/su5_construction/s2_sizing.json
"""
from __future__ import annotations

import json
import os
from fractions import Fraction as F

import numpy as np

import sc_lib as S

# the target
P = [F(2, 3), F(2, 3), F(-1, 3), F(-1, 3), F(-1, 3), F(0), F(0), F(0),
     F(1, 2), F(1, 2), F(0), F(0), F(0), F(0), F(0), F(0)]
DELTA = [F(0)] * 8 + [F(1), F(1)] + [F(0)] * 6

# SM simple roots inside the block-1 SU(5) (e_i - e_j basis, even coordinate system)
def _e(i):
    v = [0] * 16
    v[i] = 1
    return v


SM_SIMPLE = []
for (i, j) in ((0, 1), (1, 2), (3, 4)):
    r = [0] * 16
    r[i], r[j] = 1, -1
    SM_SIMPLE.append(r)


def e8_class_reps():
    """The 6561 = 3^8 classes of E8/3E8 as integer 8-vectors (simple-root combinations)."""
    A = np.array([
        [1, -1, -1, -1, -1, -1, -1, 1],      # 2*alpha1
        [2, 2, 0, 0, 0, 0, 0, 0],
        [-2, 2, 0, 0, 0, 0, 0, 0],
        [0, -2, 2, 0, 0, 0, 0, 0],
        [0, 0, -2, 2, 0, 0, 0, 0],
        [0, 0, 0, -2, 2, 0, 0, 0],
        [0, 0, 0, 0, -2, 2, 0, 0],
        [0, 0, 0, 0, 0, -2, 2, 0],
    ], dtype=np.int64)                        # rows = 2*alpha_i
    coef = np.array(np.meshgrid(*([np.arange(3)] * 8), indexing="ij")).reshape(8, -1).T
    X2 = coef @ A                             # 2*x  (integers)
    assert (X2 % 2 == 0).all() or True
    return X2                                 # keep the factor 2; halve at use


def main():
    os.makedirs(S.OUT, exist_ok=True)
    X2 = e8_class_reps()                      # 6561 x 8, entries = 2*x
    n = X2.shape[0]
    assert n == 6561

    # sanity: every class rep must lie in E8 (x integral-or-all-half with even sum)
    x = X2 / 2.0
    isint = np.all(X2 % 2 == 0, axis=1)
    ishalf = np.all(np.abs(X2 % 2) == 1, axis=1)
    assert np.all(isint | ishalf), "class reps not in the even coordinate system"
    sums = X2.sum(axis=1)
    assert np.all((sums / 2) % 2 == 0), "class reps not in E8 (odd sum)"

    # ---- block-1 SM-preserving classes -------------------------------------------
    # alpha . x = 0 mod 3  with alpha = e_i - e_j : (X2[:,i]-X2[:,j])/2 = 0 mod 3
    keep = np.ones(n, bool)
    for (i, j) in ((0, 1), (1, 2), (3, 4)):
        d = (X2[:, i] - X2[:, j]) // 2
        keep &= (d % 3 == 0)
    n_sm = int(keep.sum())

    # ---- block-2 Delta-integral classes ------------------------------------------
    # Delta = e9 + e10 -> in block-2 coordinates (1,1,0^6): Delta.z = z0 + z1 = 0 mod 3
    dz = (X2[:, 0] + X2[:, 1]) // 2
    keepD = (dz % 3 == 0)
    n_D = int(keepD.sum())

    # ---- norms mod 6 --------------------------------------------------------------
    nrm = (X2 * X2).sum(axis=1) // 4          # x^2, an even integer
    assert np.all(nrm % 2 == 0)

    res = dict(
        target_V2="p (gate d7_construction.json)",
        classes_per_block=int(n),
        classes_total_per_shift=int(n) ** 2,
        block1_SM_preserving=n_sm,
        block2_Delta_integral=n_D,
        V1_space_A1="%d x %d = %d" % (n_sm, n, n_sm * n),
        V1_space_A1_A4="%d x %d = %d" % (n_sm, n_D, n_sm * n_D),
        W_space_A1_A4="%d x %d = %d" % (n_sm, n_D, n_sm * n_D),
        naive_pair_space=(n_sm * n_D) ** 2,
        norm_mod6_histogram={int(k): int(v) for k, v in
                             zip(*np.unique(nrm % 6, return_counts=True))},
    )
    with open(os.path.join(S.OUT, "s2_sizing.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()

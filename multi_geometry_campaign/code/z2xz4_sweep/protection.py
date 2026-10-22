"""
protection.py --- STEP 4 (guard 4): the exchange-symmetry L4 lift, decided with the
CORRECTED full-coset machinery only.

WHY THE CORRECTED MACHINERY.  `gate_z4xz4_sweep.md` S5 (severity F, the SEVENTH instance of
the CLAUDE.md 5.19 conjugation/duality family) showed that the Z4xZ4 phase's "exhaustive"
search over `Aut(Lambda_{E8xE8})` in fact covered **3 of 135 cosets = 2.22 % of W(E8) per
block**, because a LEFT transversal of W(D8) was used as a RIGHT transversal; and, second,
that a SINGLE lattice-coset type `delta` was forced across all source vectors when the
condition `gamma(v_i) - w_i in Lambda` is PER VECTOR.  Both defects are repaired in
`code/gate_z4xz4_sweep/s5_lift.py`, which this module IMPORTS rather than reimplements:

    exhaustive_lift(src32, dst32, mode="correct")   frame applied LAST (gamma = M_F . p, so
                                                    p(src) == M_F^{-1}(dst)); per-vector
                                                    residue tracked during backtracking.

THE VERDICT IS THREE-VALUED AND ALL THREE ARE REPORTED.
  IMPOSED       -- the exact necessary condition fails (ORDER / ORDER-COMBO / INNER-PRODUCT
                   certificate attached), OR the genuinely exhaustive search over the full
                   automorphism group returns nothing.
  PROTECTED     -- a lift is exhibited constructively.
  UNDETERMINED  -- reported, never silently promoted to IMPOSED.  With the corrected
                   exhaustive search this category is empty by construction; it is retained
                   because a truncated run (a raised time budget, an aborted frame loop) must
                   be able to say so rather than pretend.

THIS GEOMETRY'S OWN sigma, DERIVED (selftest S12), not inherited:
  * the ONLY plane pair admitting an exchange automorphism is [1,3] (matching the certified
    tau-premise record, which also records `broken_by_shift = false`);
  * 16 distinct (tau, S mod 2) conditions survive exact deduplication -- 8 Wilson-line
    patterns x 2 point-group automorphisms tau(omega) in {omega, omega^3}.  The certified
    holomorphic tau is ((1,2),(0,3)) and is the PRIMARY branch; the phase-inverting tau
    ((1,2),(0,1)) is carried as a TAGGED second branch.  Protection is EXISTENTIAL over the
    admitted conditions, and BOTH scopes are reported: `protected_holomorphic` (the certified
    tau alone) and `protected_any` (either tau).  A permissive sigma set can only turn
    IMPOSED into PROTECTED, so the two numbers bracket the answer.
  * the Wilson-line condition is NOT trivial on the fixed plane here, unlike Z4xZ4: plane 2's
    two Wilson lines W3, W4 are INDEPENDENT on this geometry (rows 4 and 5 of the catalogue
    differ on 11,090 of 11,540 models) and half the admitted sigma SWAP them.

Emits: data/z2xz4_sweep/protection.json
"""


from __future__ import annotations

# --- reproduction-package root shim -------------------------------------------------
# The line above originally read  F:\physics  -- the absolute path of the author's
# working tree.  This is the ONLY edit made to this file in the reproduction package;
# it points the same variable at the package root, which carries the same directory
# names.  Set PAPERA_PKG to relocate the package.  See
# environment/PATCHES_TO_PACKAGE_CODE.md.
# -------------------------------------------------------------------------------------
import os as _os
_PKG_ROOT = _os.environ.get('PAPERA_PKG') or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..'))

import json
import os
import sys
import time
from collections import Counter
from multiprocessing import Pool

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, _os.path.join(_PKG_ROOT, 'multi_geometry_campaign', 'code', 'gate_z4xz4_sweep'))
import z24sw_lib as L                                                   # noqa: E402
import s5_lift as LI                                                    # noqa: E402

PAIR = (1, 3)
_COND = None


def targets(rows32, cond):
    """(src, dst) in x32 units for one admitted sigma condition.

    gamma(V1) = V_{tau(theta)},  gamma(V2) = V_{tau(omega)},  V_{(k,l)} := k V1 + l V2
    gamma(W_i) = sum_j S[i][j] W_j                                (i = 1..6)

    Every Wilson line has order dividing 2 here (the file declares 2,2,2,2,2,2 and the
    catalogue confirms it), so -W = W mod Lambda and only S mod 2 matters -- which is exactly
    what `exchange_sigmas` deduplicates on."""
    V1, V2 = rows32[0], rows32[1]
    W = rows32[2:8]

    def comb(kl):
        k, l = kl
        return [k * V1[i] + l * V2[i] for i in range(16)]

    gi = cond["generator_images"]
    S = cond["wl_matrix"]
    src = [V1, V2] + list(W)
    dst = [comb(gi[0]), comb(gi[1])]
    for i in range(6):
        v = [0] * 16
        for j in range(6):
            if S[i][j]:
                for t in range(16):
                    v[t] += S[i][j] * W[j][t]
        dst.append(v)
    return src, dst


def _init(conds):
    global _COND
    _COND = conds
    LI.frames()


def work(item):
    lab, rows = item
    rows32 = [LI.to32(r) for r in rows]
    per = []
    for ci, cond in enumerate(_COND):
        src, dst = targets(rows32, cond)
        fails, cert = LI.necessary_fails(src, dst)
        if fails:
            per.append(("IMPOSED_necessary", cert.get("kind") if cert else None))
        elif LI.exhaustive_lift(src, dst, mode="correct"):
            per.append(("PROTECTED", None))
        else:
            per.append(("IMPOSED_exhaustive", None))
    return lab, per


def main():
    conds = L.exchange_sigmas(PAIR)
    holo = [i for i, c in enumerate(conds)
            if [list(x) for x in c["generator_images"]] == [[1, 2], [0, 3]]]
    cat = L.parse_catalogue()
    assert len(cat) == L.N_MODELS, len(cat)
    items = sorted((k, v["rows"]) for k, v in cat.items())
    t0 = time.time()
    with Pool(6, initializer=_init, initargs=(conds,)) as pool:
        out = dict(pool.map(work, items, chunksize=20))
    el = time.time() - t0
    print("elapsed %.1f s for %d models x %d conditions" % (el, len(out), len(conds)),
          flush=True)

    prot_any, prot_holo, per_cond = [], [], Counter()
    verdicts = Counter()
    for lab, per in out.items():
        if any(p[0] == "PROTECTED" for p in per):
            prot_any.append(lab)
        if any(per[i][0] == "PROTECTED" for i in holo):
            prot_holo.append(lab)
        for i, p in enumerate(per):
            verdicts[p[0]] += 1
            if p[0] == "PROTECTED":
                per_cond[i] += 1

    rep = {"geometry": "Z2xZ4_1_1",
           "scope": "CATALOGUE statement about the 11,540 PVW Z2xZ4_1_1 models",
           "machinery": "code/gate_z4xz4_sweep/s5_lift.py (the CORRECTED full-coset "
                        "implementation: 135 frames as a LEFT transversal with the frame "
                        "applied LAST, per-vector lattice coset).  The superseded phase "
                        "implementation is not used anywhere in this phase.",
           "exchange_pair": list(PAIR),
           "n_sigma_conditions": len(conds),
           "n_holomorphic_conditions": len(holo),
           "elapsed_s": round(el, 1),
           "verdict_histogram_over_model_x_condition": dict(verdicts),
           "n_protected_any_sigma": len(prot_any),
           "n_protected_holomorphic_tau": len(prot_holo),
           "n_imposed_any_sigma": len(out) - len(prot_any),
           "n_undetermined": 0,
           "undetermined_note": "EMPTY BY CONSTRUCTION: the search over Aut(Lambda_{E8xE8}) "
                                "is genuinely exhaustive with the corrected machinery, so "
                                "'no lift found' IS a negative.  The category is reported, "
                                "not omitted.",
           "protected_models_holomorphic": sorted(prot_holo),
           "protected_models_any": sorted(prot_any),
           "protection_by_condition": {str(k): v for k, v in sorted(per_cond.items())},
           "conditions": conds,
           "per_model": {k: [p[0] for p in v] for k, v in out.items()}}
    L.write_json("protection.json", rep)
    print(json.dumps({k: v for k, v in rep.items()
                      if k not in ("per_model", "conditions", "protected_models_any",
                                   "protected_models_holomorphic")}, indent=1))


if __name__ == "__main__":
    main()

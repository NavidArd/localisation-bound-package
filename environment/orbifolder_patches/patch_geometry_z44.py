"""
patch_geometry_z44.py --- R-symmetry patches for `Geometry_Z4xZ4_1_1.txt`.

Phase: multi_geometry_campaign / Z4xZ4 foundations (docs/z4xz4_foundations_task.md), item 3.

WHY A PATCH IS NEEDED (the three-surface check, surface C).
Orbifolder's `cyukawacouplings.cpp::CreatePreCouplings` DOES iterate the geometry file's
`DiscreteRSymmetries` --- surfaces A (code path exists) and B (reachable) pass --- but ALL 138
shipped geometry files declare ZERO R symmetries, and `Geometry_Z4xZ4_1_1.txt` is one of them
(5 non-R, 0 R; re-verified fresh on THIS file by `z44lib.GeomFile.discrete_symmetries`).
Surface C therefore FAILS on the stock file: the rule has no data to act on (CLAUDE.md 5.16
instance #2 + #3).  The patch supplies what can be supplied.

WHAT IS EXPRESSIBLE --- AND ON THIS GEOMETRY THE ANSWER IS "NOTHING, EXACTLY".
`CField::GetDiscreteCharge` computes  sum_k (q_sh[k] + Osci[k]) * ChargeOperator[k]  and has NO
access to the gamma phases, so only GAMMA-FREE laws are declarable.  The published two-generator
law 0705.3075 (15) is per-plane and gamma-free AS PRINTED, but it is complete only where the
plane is PRIME; in a non-prime plane the published Z_N treatment (1301.2322 (3.44)/(3.45))
replaces the naive per-plane law by a COMBINED law carrying `+ N gamma`.

  Z4xZ4_1_1 plane orders are (4, 4, 4).  ALL THREE PLANES ARE NON-PRIME.  There is therefore
  NO plane on this geometry whose R-charge law is both exactly expressible and known-correct.
  The defensible patch set is EMPTY --- reported as EMPTY, never as a vacuous pass.

Compare, so the scoping is visible: Z4_1_1 had 3 of 4 printed laws exactly expressible (the
largest fraction in the campaign); Z3xZ6_1_1 had its one prime plane expressible; Z4xZ4_1_1 has
none.  This is the sharpest form the campaign's central RULE-GAP has taken.

Two artefacts are emitted anyway, both TAGGED:
  Geometry_Z4xZ4_1_1_Rsym_all3.txt   --- the naive per-plane law in all three planes.
                                         KNOWN-POTENTIALLY-TOO-STRONG in all three.
                                         BRANCH MEASUREMENT ONLY --- no physics verdict may
                                         rest on a run using it.
  Geometry_Z4xZ4_1_1_Rsym_none.txt   --- the stock file copied verbatim, so that a run pair
                                         (stock vs all3) measures exactly the effect of the
                                         over-strong branch and nothing else.

Block format (read off the certified two-generator precedent
`rcharge_census/data/patched_geometry/Geometry_Z3xZ3_1_1_Rsym.txt`, whose R blocks orbifolder
consumed successfully in that phase's gate-validated run):

    new R symmetry
    <order>
    <superpotential charge>
    <4 operator entries: slot 0 = 4d helicity, slots 1..3 = the three complex planes>

The superpotential charge is written `-1/1` because orbifolder's own charge is
`GetDiscreteCharge = -R^paper` (the sign convention recorded once in `rules_z4xz4.md` 0).

Verification (run every time, temp file removed in `finally`): the patched file is re-parsed by
this phase's own independent engine and must reproduce, unchanged, the twists, the plane orders,
the per-sector fixed-point counts and the simultaneous-orbit size multisets.

Outputs: data/z4xz4/patched_geometry/*.txt  (+ data/z4xz4/patches_z44.json)
Run:     python multi_geometry_campaign/code/z4xz4/patch_geometry_z44.py
"""

from __future__ import annotations

import json
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z44lib import Geom, GEOM_FILE, OUT_DIR, _is_prime  # noqa: E402

PATCH_DIR = os.path.join(OUT_DIR, "patched_geometry")
GEOMDIR = os.path.dirname(GEOM_FILE)


def rsym_block(order: int, plane: int, w_charge: str = "-1/1") -> str:
    op = ["0", "0", "0", "0"]
    op[plane] = "1"                       # slot 0 is the 4d helicity slot; planes are 1..3
    return "new R symmetry\n%d\n%s\n%s\n" % (order, w_charge, " ".join(op))


def build(src: str, planes) -> str:
    blocks = "".join(rsym_block(o, p) for p, o in planes)
    if "end discrete symmetries" not in src:
        raise AssertionError("geometry file has no discrete-symmetry block to extend")
    return src.replace("end discrete symmetries", blocks + "end discrete symmetries")


def verify(path: str, ref: Geom) -> dict:
    """Re-parse the patched file with this phase's own engine; nothing geometric may move."""
    tmp = os.path.join(GEOMDIR, "Geometry__z44_tmp_check.txt")
    try:
        shutil.copyfile(path, tmp)
        g = Geom(tmp)
        out = {
            "twists_unchanged": ([str(x) for x in g.v1] == [str(x) for x in ref.v1]
                                 and [str(x) for x in g.v2] == [str(x) for x in ref.v2]),
            "plane_orders_unchanged": ([g.plane_order(a) for a in (1, 2, 3)]
                                       == [ref.plane_order(a) for a in (1, 2, 3)]),
            "gram_unchanged": g.gram == ref.gram,
            "R_symmetries_declared": g.f.discrete_symmetries()["R"],
            "non_R_symmetries_declared": g.f.discrete_symmetries()["non_R"],
            "WL_block_unchanged": (g.f.wl_orders_declared() == ref.f.wl_orders_declared()
                                   and g.f.identical_wls_declared()
                                   == ref.f.identical_wls_declared()),
        }
        fp_ok, orb_ok = True, True
        for el in g.elements():
            if el == (0, 0):
                continue
            for a in (1, 2, 3):
                fa, fb = g.fixed_points(el, a), ref.fixed_points(el, a)
                if (fa is None) != (fb is None) or (fa is not None and len(fa) != len(fb)):
                    fp_ok = False
            if g.simultaneous_orbit_sizes(el) != ref.simultaneous_orbit_sizes(el):
                orb_ok = False
        out["fixed_point_counts_unchanged"] = fp_ok
        out["simultaneous_orbit_sizes_unchanged"] = orb_ok
        out["all_pass"] = all([out["twists_unchanged"], out["plane_orders_unchanged"],
                               out["gram_unchanged"], out["WL_block_unchanged"], fp_ok, orb_ok])
        return out
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def main() -> None:
    ref = Geom()
    os.makedirs(PATCH_DIR, exist_ok=True)
    with open(GEOM_FILE, "r", encoding="utf-8") as fh:
        src = fh.read()

    orders = {a: ref.plane_order(a) for a in (1, 2, 3)}
    prime = [a for a in (1, 2, 3) if _is_prime(orders[a])]
    nonprime = [a for a in (1, 2, 3) if a not in prime]

    specs = [
        ("Geometry_Z4xZ4_1_1_Rsym_none.txt", [],
         "the stock file, copied verbatim.  The DEFENSIBLE patch set on this geometry is "
         "EMPTY-BY-DERIVATION: there is no prime plane, so no printed R-charge law is both "
         "gamma-free-as-printed and known-complete.  Reported as EMPTY, not as a pass."),
        ("Geometry_Z4xZ4_1_1_Rsym_all3.txt", [(a, orders[a]) for a in (1, 2, 3)],
         "KNOWN-POTENTIALLY-TOO-STRONG in ALL THREE planes %s: every plane is non-prime, so in "
         "every plane the published Z_N treatment replaces the naive per-plane law by a "
         "combined law carrying +N*gamma, which GetDiscreteCharge cannot express.  BRANCH "
         "MEASUREMENT ONLY --- no physics verdict may rest on a run using this file."
         % nonprime),
    ]

    rep = {"geometry": "Z4xZ4_1_1", "plane_orders": orders,
           "prime_planes": prime, "non_prime_planes": nonprime,
           "defensible_patch_set_is_empty": len(prime) == 0,
           "stock_file_declares": ref.f.discrete_symmetries(),
           "comparison_across_the_campaign": {
               "Z4_1_1": "2 non-prime + 1 prime; 3 of 4 printed laws exactly expressible",
               "Z3xZ6_1_1": "2 non-prime + 1 prime; the prime-plane law exactly expressible",
               "Z4xZ4_1_1": "3 non-prime + 0 prime; NOTHING exactly expressible"},
           "three_surface_check": {
               "A_code_path_exists": "CreatePreCouplings iterates DiscreteRSymmetries "
                                     "(source level, verified by the certified pipeline)",
               "B_reachable": "yes --- the loop runs on every model",
               "C_input_data_populated": "NO on the stock file: 0 R symmetries declared "
                                         "(5 non-R, 0 R, re-verified fresh on this file)"},
           "patches": []}

    for name, planes, note in specs:
        path = os.path.join(PATCH_DIR, name)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(build(src, planes) if planes else src)
        v = verify(path, ref)
        rep["patches"].append({"file": path, "planes_declared": [p for p, _ in planes],
                               "orders": [o for _, o in planes], "note": note,
                               "verification": v})
        print("%-40s planes %s orders %s  verification all_pass=%s (R declared: %s)"
              % (name, [p for p, _ in planes], [o for _, o in planes],
                 v["all_pass"], v["R_symmetries_declared"]))
        print("    %s" % note)

    with open(os.path.join(OUT_DIR, "patches_z44.json"), "w", encoding="utf-8") as fh:
        json.dump(rep, fh, indent=1, default=str)
    print("prime planes: %s   non-prime: %s   DEFENSIBLE PATCH SET EMPTY: %s"
          % (prime, nonprime, rep["defensible_patch_set_is_empty"]))
    print("wrote", os.path.join(OUT_DIR, "patches_z44.json"))


if __name__ == "__main__":
    main()

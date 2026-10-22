"""
patch_geometry_z36.py --- R-symmetry patches for `Geometry_Z3xZ6_1_1.txt`.

Phase: multi_geometry_campaign / Z3xZ6 foundations (docs/z3xz6_foundations_task.md), item 3.

WHY A PATCH IS NEEDED (the three-surface check, surface C).
Orbifolder's `cyukawacouplings.cpp::CreatePreCouplings` DOES iterate the geometry file's
`DiscreteRSymmetries` --- surfaces A (code path exists) and B (reachable) pass --- but ALL 138
shipped geometry files declare ZERO R symmetries, and `Geometry_Z3xZ6_1_1.txt` is one of them
(3 non-R, 0 R; re-verified fresh on THIS file by `z36lib.GeomFile.discrete_symmetries`).
Surface C therefore FAILS on the stock file: the rule has no data to act on (CLAUDE.md 5.16
instance #2 + #3).  The patch supplies it.

WHAT IS EXPRESSIBLE.
`CField::GetDiscreteCharge` computes  sum_k (q_sh[k] + Osci[k]) * ChargeOperator[k]  and has NO
access to the gamma phases.  Only GAMMA-FREE laws are declarable.  On Z3xZ6_1_1:

  * plane 2 is PRIME (order 3): its per-plane law is gamma-free and EXACTLY expressible.
  * planes 1 and 3 are NON-PRIME (order 6): the published Z_N treatment (1301.2322 (3.44)/(3.45))
    replaces the naive per-plane law there by a COMBINED law carrying `+ N gamma`.  Declaring the
    naive per-plane law in a non-prime plane is therefore KNOWN-POTENTIALLY-TOO-STRONG --- the
    exact shape of the Z4 pass's `_Rsym_full4` caveat, inherited here and NOT resolved.

Two patches are emitted, and the difference between them IS the phase's central rule question:
  Geometry_Z3xZ6_1_1_Rsym_prime.txt  --- plane 2 only.  DEFENSIBLE; use this one for any run
                                         whose output carries a physics verdict.
  Geometry_Z3xZ6_1_1_Rsym_all3.txt   --- all three planes.  TAGGED KNOWN-POTENTIALLY-TOO-STRONG
                                         in planes 1 and 3; for branch measurement only.

Block format (read off the certified two-generator precedent
`rcharge_census/data/patched_geometry/Geometry_Z3xZ3_1_1_Rsym.txt`, whose R blocks orbifolder
consumed successfully in that phase's gate-validated run):

    new R symmetry
    <order>
    <superpotential charge>
    <4 operator entries: slot 0 = 4d helicity, slots 1..3 = the three complex planes>

The superpotential charge is written `-1/1` because orbifolder's own charge is
`GetDiscreteCharge = -R^paper` (the sign convention recorded once in `rules_z3xz6.md` 0).

Verification (run every time, temp file removed in `finally`): the patched file is re-parsed by
this phase's own independent engine and must reproduce, unchanged, the twists, the plane orders,
the per-sector fixed-point counts and the simultaneous-orbit size multisets.

Outputs: data/z3xz6/patched_geometry/*.txt  (+ data/z3xz6/patches_z36.json)
Run:     python multi_geometry_campaign/code/z3xz6/patch_geometry_z36.py
"""

from __future__ import annotations

import json
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from z36lib import Geom, GEOM_FILE, OUT_DIR, sector_label  # noqa: E402

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
    tmp = os.path.join(GEOMDIR, "Geometry__z36_tmp_check.txt")
    try:
        shutil.copyfile(path, tmp)
        g = Geom(tmp)
        out = {
            "twists_unchanged": ([str(x) for x in g.v1] == [str(x) for x in ref.v1]
                                 and [str(x) for x in g.v2] == [str(x) for x in ref.v2]),
            "plane_orders_unchanged": ([g.plane_order(a) for a in (1, 2, 3)]
                                       == [ref.plane_order(a) for a in (1, 2, 3)]),
            "R_symmetries_declared": g.f.discrete_symmetries()["R"],
            "non_R_symmetries_declared": g.f.discrete_symmetries()["non_R"],
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
                               fp_ok, orb_ok])
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
    prime = [a for a in (1, 2, 3) if orders[a] in (2, 3, 5, 7, 11)]
    nonprime = [a for a in (1, 2, 3) if a not in prime]

    specs = [
        ("Geometry_Z3xZ6_1_1_Rsym_prime.txt", [(a, orders[a]) for a in prime],
         "DEFENSIBLE: prime plane(s) only; gamma-free and exactly expressible."),
        ("Geometry_Z3xZ6_1_1_Rsym_all3.txt", [(a, orders[a]) for a in (1, 2, 3)],
         "KNOWN-POTENTIALLY-TOO-STRONG in the non-prime planes %s: the published Z_N treatment "
         "replaces the naive per-plane law there by a combined law carrying +N*gamma, which "
         "GetDiscreteCharge cannot express.  Branch measurement only." % nonprime),
    ]

    rep = {"geometry": "Z3xZ6_1_1", "plane_orders": orders,
           "prime_planes": prime, "non_prime_planes": nonprime,
           "stock_file_declares": ref.f.discrete_symmetries(),
           "three_surface_check": {
               "A_code_path_exists": "CreatePreCouplings iterates DiscreteRSymmetries (source level, "
                                     "verified by the certified pipeline)",
               "B_reachable": "yes --- the loop runs on every model",
               "C_input_data_populated": "NO on the stock file: 0 R symmetries declared "
                                         "(3 non-R, 0 R, re-verified fresh on this geometry file)"},
           "patches": []}

    for name, planes, note in specs:
        path = os.path.join(PATCH_DIR, name)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(build(src, planes))
        v = verify(path, ref)
        rep["patches"].append({"file": path, "planes_declared": [p for p, _ in planes],
                              "orders": [o for _, o in planes], "note": note, "verification": v})
        print("%-42s planes %s orders %s   verification all_pass=%s (R declared: %s)"
              % (name, [p for p, _ in planes], [o for _, o in planes],
                 v["all_pass"], v["R_symmetries_declared"]))
        print("    %s" % note)

    with open(os.path.join(OUT_DIR, "patches_z36.json"), "w", encoding="utf-8") as fh:
        json.dump(rep, fh, indent=1, default=str)
    print("wrote", os.path.join(OUT_DIR, "patches_z36.json"))


if __name__ == "__main__":
    main()

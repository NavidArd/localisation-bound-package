"""GATE gate_su5_bound D5 -- STATE-LEVEL theorem test + the Z3xZ6 hidden-non-singlet DROP.

TWO DEFECTS FOUND UPSTREAM, TESTED HERE.

(A) RECORDS vs STATES.  A spectrum record that is a NON-SINGLET of a hidden factor is
    several distinct massless left-chiral states sharing the SM role, the cell and |p|^2,
    differing by a hidden root (Delta in Lambda, Delta orthogonal to W, |Delta|^2 >= 2).
    The theorem therefore demands R >= 1/2 for such a record.  verify.py counts records,
    so it cannot see this.  Here every record is expanded to its states and the general
    pair test  sqrt(R_i) + sqrt(R_j) >= sqrt(2)  is applied to the expanded list.

(B) THE TWO GEOMETRIES DO NOT USE THE SAME DEFINITION.
      z24sw_lib.role_census  assigns a role to hidden NON-singlets (carrying `mult`);
      z36sw_lib.role_census  sets role = None for every hidden non-singlet.
    So Z3xZ6_1_1's 37,334 / 125,146 silently EXCLUDE a population that Z2xZ4_1_1's
    67,555 / 221,005 include.  Excluded 10-role fields are the dangerous direction: one of
    them at a cell that already carries a 10-role is a HIDDEN DOUBLING.  This script
    re-runs the Z3xZ6 role assignment with the gate removed and counts what comes back.

R is computed here as R = |p_sh|^2 - |w|^2 straight off the PRINTED p_sh, with
|w|^2 = 6/5 (10-roles) and 4/5 (5bar-roles) -- the GUT-normalised values forced by
k_Y = 5/3, which the phase verified by explicit projection on all 481,728 fields.  No twist
vector and no delta-c enters, so this test is independent of the hard-coded GEOM dict.

Outputs data/gate_su5_bound/d5_hidden.json
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"

import glob                                                            # noqa: E402
import gzip                                                            # noqa: E402
import json                                                            # noqa: E402
import multiprocessing as mp                                           # noqa: E402
import sys                                                             # noqa: E402
from collections import Counter, defaultdict                           # noqa: E402
from fractions import Fraction as F                                    # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CAMP = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(CAMP, "data", "gate_su5_bound")

ROLES10 = {"Q", "uc", "ec", "Qbar", "ucbar", "ecbar"}
ROLES5 = {"dc", "L", "Hu", "dcbar"}
W2 = {"10": F(6, 5), "5bar": F(4, 5)}

SM_ROLES = {(3, 2, F(1, 6)): "Q", (-3, 2, F(-1, 6)): "Qbar",
            (-3, 1, F(-2, 3)): "uc", (3, 1, F(2, 3)): "ucbar",
            (-3, 1, F(1, 3)): "dc", (3, 1, F(-1, 3)): "dcbar",
            (1, 2, F(-1, 2)): "L", (1, 2, F(1, 2)): "Hu",
            (1, 1, F(1)): "ec", (1, 1, F(-1)): "ecbar"}

GEOM = {
    "Z2xZ4_1_1": dict(lib="z2xz4_sweep", mod="z24sw_lib", M=2, N=4,
                      census=os.path.join(CAMP, "data", "z2xz4_sweep", "raw", "census")),
    "Z3xZ6_1_1": dict(lib="z3xz6_sweep", mod="z36sw_lib", M=3, N=6,
                      census=os.path.join(CAMP, "data", "z3xz6_sweep", "raw", "census")),
}


def pair_ok(a, b):
    t = 2 - a - b
    if t <= 0:
        return True
    return 4 * a * b >= t * t


def fam_of(ro):
    return "10" if ro in ROLES10 else ("5bar" if ro in ROLES5 else None)


def worker(args):
    path, geom = args
    G = GEOM[geom]
    sys.path.insert(0, os.path.join(CAMP, "code", G["lib"]))
    Z = __import__(G["mod"])
    acc = dict(models=0,
               hid_x_R=defaultdict(Counter),        # fam -> "hid|R" -> n records
               state_viol=0, state_viol_ex=[],
               rec_viol=0,
               mult_state=defaultdict(Counter),
               ungated_extra=defaultdict(Counter),  # fam -> "hid" -> n records recovered
               ungated_extra_R=defaultdict(Counter),
               ungated_new_doubling=0, ungated_new_doubling_ex=[],
               gated_roles=0, ungated_roles=0,
               ungated_branch_same=0, ungated_branch_diff=0)
    txt = gzip.open(path, "rt", errors="replace").read()
    for lab, blk in Z.split_models(txt):
        pre = blk.split("Analyzing")[0]
        try:
            if geom == "Z2xZ4_1_1":
                spec = Z.parse_spectrum(pre, M=2, N=4)
                factors = Z.parse_gauge_factors(blk)
                cfgs = Z.parse_smconfig(blk)
                if not cfgs or cfgs[0][0] is None or not spec:
                    continue
                tY, tbl = cfgs[0]
                printed = defaultdict(int)
                for row in tbl:
                    printed[(row[0], row[1], row[2])] += row[3]
                pos = Z.observable_positions(spec, tY, factors, dict(printed))
                if pos is None:
                    continue
                a3, a2, sY0 = pos
                roles, diag = Z.role_census(spec, tY, a3, a2, sY0)
                sC = diag["sC_selected"]
                sY = diag["sY_selected"]
            else:
                spec = Z.parse_spectrum(pre)
                gauge = Z.parse_gauge(blk)
                if gauge["tY"] is None or not spec:
                    continue
                pos = Z.observable_positions(spec, gauge["tY"], gauge)
                if pos is None:
                    continue
                a3, a2, sY0, tY, bi = pos
                roles, diag = Z.role_census(spec, tY, a3, a2, sY0)
                sC = diag["sC_selected"]
                sY = diag["sY_selected"]
        except Exception:                                              # noqa: BLE001
            continue
        acc["models"] += 1

        # ---- ungated re-assignment (same slots, same branch, no hidden_singlet gate)
        ung = {}
        for l, r in spec.items():
            if r["rep"] is None or r["chir"] != "l" or r["p_sh"] is None:
                continue
            if a3 >= len(r["rep"]) or a2 >= len(r["rep"]):
                continue
            Y = sY * sY0 * Z.dot(tY, r["p_sh"])
            c3 = r["rep"][a3]
            c3 = (sC * c3) if abs(c3) != 1 else 1
            c2 = abs(r["rep"][a2])
            hid = 1
            for k, dd in enumerate(r["rep"]):
                if k not in (a3, a2):
                    hid *= abs(dd)
            ung[l] = (SM_ROLES.get((c3, c2, Y)), hid)

        cells = defaultdict(lambda: defaultdict(list))
        cells_u = defaultdict(lambda: defaultdict(list))
        for l, r in spec.items():
            if r["p_sh"] is None or r["chir"] != "l":
                continue
            kl = tuple(r["kl"])
            if kl == (0, 0):
                continue
            key = (kl, tuple(r["fp"]) if r["fp"] else None)
            p2 = sum(x * x for x in r["p_sh"])
            gro = (roles.get(l) or {}).get("role")
            uro, uhid = ung.get(l, (None, 1))
            if gro is not None:
                fam = fam_of(gro)
                if fam:
                    acc["gated_roles"] += 1
                    Rv = p2 - W2[fam]
                    hid = (roles.get(l) or {}).get("mult", uhid)
                    acc["hid_x_R"][fam]["%d|%s" % (hid, Rv)] += 1
                    cells[key][gro] += [Rv] * hid
            if uro is not None:
                fam = fam_of(uro)
                if fam:
                    acc["ungated_roles"] += 1
                    Rv = p2 - W2[fam]
                    cells_u[key][uro] += [Rv] * uhid
                    if gro is None:
                        acc["ungated_extra"][fam][str(uhid)] += 1
                        acc["ungated_extra_R"][fam][str(Rv)] += 1
        # ---- state-level theorem test on the GATED assignment
        for key, d in cells.items():
            for ro, Rs in d.items():
                fam = fam_of(ro)
                acc["mult_state"][fam][str(len(Rs))] += 1
                if len(Rs) < 2:
                    continue
                Rs = sorted(Rs, reverse=True)
                ok = any(pair_ok(Rs[i], Rs[j])
                         for i in range(len(Rs)) for j in range(i + 1, len(Rs)))
                if not ok:
                    acc["state_viol"] += 1
                    if len(acc["state_viol_ex"]) < 20:
                        acc["state_viol_ex"].append(
                            dict(model=lab, geom=geom, kl=list(key[0]),
                                 fp=list(key[1]) if key[1] else None, role=ro,
                                 R=[str(x) for x in Rs],
                                 file=os.path.basename(path)))
        # ---- does removing the gate create a 10-role doubling the census cannot see?
        for key, d in cells_u.items():
            for ro, Rs in d.items():
                if fam_of(ro) != "10":
                    continue
                nrec_gated = len(cells[key].get(ro, []))
                if len(Rs) >= 2 and nrec_gated < len(Rs):
                    acc["ungated_new_doubling"] += 1
                    if len(acc["ungated_new_doubling_ex"]) < 20:
                        acc["ungated_new_doubling_ex"].append(
                            dict(model=lab, geom=geom, kl=list(key[0]),
                                 fp=list(key[1]) if key[1] else None, role=ro,
                                 n_states_ungated=len(Rs), n_states_gated=nrec_gated,
                                 R=[str(x) for x in sorted(Rs, reverse=True)],
                                 file=os.path.basename(path)))
    for k in ("hid_x_R", "mult_state", "ungated_extra", "ungated_extra_R"):
        acc[k] = {a: dict(b) for a, b in acc[k].items()}
    return acc


def main():
    res = {}
    for geom in ("Z2xZ4_1_1", "Z3xZ6_1_1"):
        files = sorted(glob.glob(os.path.join(GEOM[geom]["census"], "*.txt.gz")))
        with mp.Pool(6) as pool:
            accs = pool.map(worker, [(f, geom) for f in files])
        out = dict(models=0, state_viol=0, rec_viol=0, gated_roles=0, ungated_roles=0,
                   ungated_new_doubling=0, state_viol_ex=[], ungated_new_doubling_ex=[])
        cs = {k: defaultdict(Counter) for k in
              ("hid_x_R", "mult_state", "ungated_extra", "ungated_extra_R")}
        for a in accs:
            for k in ("models", "state_viol", "rec_viol", "gated_roles", "ungated_roles",
                      "ungated_new_doubling"):
                out[k] += a[k]
            for k in cs:
                for f, d in a[k].items():
                    cs[k][f].update(d)
            out["state_viol_ex"] += a["state_viol_ex"][:3]
            out["ungated_new_doubling_ex"] += a["ungated_new_doubling_ex"][:3]
        for k in cs:
            out[k] = {a_: dict(b) for a_, b in cs[k].items()}
        out["state_viol_ex"] = out["state_viol_ex"][:25]
        out["ungated_new_doubling_ex"] = out["ungated_new_doubling_ex"][:25]
        # hidden multiplicity >= 2 with R < 1/2 -> theorem-fatal records
        fatal = {}
        for fam, d in out["hid_x_R"].items():
            for k, n in d.items():
                h, R = k.split("|")
                if int(h) >= 2 and F(R) < F(1, 2):
                    fatal["%s|%s" % (fam, k)] = n
        out["records_hidden_nonsinglet_with_R_below_half"] = fatal
        res[geom] = out
        print("=====", geom)
        print(json.dumps({k: v for k, v in out.items()
                          if k not in ("state_viol_ex", "ungated_new_doubling_ex")},
                         indent=1, default=str)[:5000])
        for e in out["state_viol_ex"][:5]:
            print("  STATE-VIOL", e)
        for e in out["ungated_new_doubling_ex"][:5]:
            print("  UNGATED-DOUBLING", e)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "d5_hidden.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=str)


if __name__ == "__main__":
    main()

"""
selftest.py --- STEP 0: THE FOUNDATIONS-LITE BATTERY (sweep guard 9).

`Z2xZ4_1_1` inherits certified square-plane machinery (gate_z4xz4_foundations.md K2 was
measured ON THIS GEOMETRY) but has never had a foundations phase of its own.  This battery
is therefore the gate before the sweep, with HARD ASSERTS throughout.  ANY FAILURE = STOP,
report a Z24F-gap, DO NOT SWEEP.

  S1  geometry: twists, planes, Gram, plane orders, roto-translation status -- all DERIVED
      from the stock geometry file, none assumed.
  S2  the {3,4,6} LETTER-CARRYING THEOREM applied to the ORDER-2 generator's plane:
      plane 2 is never capable and carries no letter; it carries 0 or exactly 2 localisable
      legs, which is exactly RULE NC's premise.
  S3  Wilson-line inventory: SNF of (1 - R_a) per plane -> the WL order group; the
      identical-WL identifications DERIVED from the point-group action on the lattice basis
      and checked against the file's own `identical WLs` block.
  S4  the ALPHABET reproduced against the CERTIFIED phase-2 payload
      `alphabet_survey/data/phase2/Z2xZ4_1_1.json`, channel by channel, cell by cell.
  S5  SECTOR INVENTORY, TWO ENGINES: the certified engine's simultaneous-orbit count vs the
      geometry file's own printed `constructing elements`, plus an independent
      determinant/brute-force fixed-point count per sector-plane.
  S6  PARSER ROW IDENTITIES DERIVED from the catalogue itself plus the file's declarations
      (the S2-inversion lesson: measure, do not assume).
  S7  the [F-1] fixed-point-vs-raw-n_a pilot check REPAIRED and PASSING (132 = the published
      1301.2322 Table 4 value), per `results/gate_rule_nc.md`.
  S8  RULE NC: the orbit-vs-equality premise verified on every NC channel-plane; NC's bite
      and NON-VACUITY measured against a null model.
  S9  the CAP QUOTED, not re-tested: the incidental space-group-only bucket scan must
      reproduce K2's published Z2xZ4_1_1 histogram {1: 168, 2: 24} with max 2.
  S10 this geometry's OWN [D-3] degenerate-prefactor census, derived.
  S11 catalogue integrity: SHA-256, model count, 8 rows per model.
  S12 the exchange automorphism sigma DERIVED; pairs, tau, and the induced Wilson-line
      condition, checked against the certified tau-premise and phase-2 records.

Emits: data/z2xz4_sweep/selftest.json
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from collections import Counter
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z24sw_lib as L                                                   # noqa: E402
from survey_geom import mm, one_minus, I2, ap                           # noqa: E402

FAIL = []


def chk(name, cond, detail=""):
    ok = bool(cond)
    print("  [%s] %-46s %s" % ("PASS" if ok else "FAIL", name, detail))
    if not ok:
        FAIL.append(name)
    return ok


def snf2(M):
    """Smith normal form invariant factors (d1 | d2) of an integer 2x2 matrix."""
    from math import gcd
    a, b, c, d = int(M[0][0]), int(M[0][1]), int(M[1][0]), int(M[1][1])
    g1 = gcd(gcd(abs(a), abs(b)), gcd(abs(c), abs(d)))
    det = abs(a * d - b * c)
    if g1 == 0:
        return (0, 0)
    if det == 0:
        return (g1, 0)
    return (g1, det // g1)


def main():
    rep = {"phase": "z2xz4_sweep", "step": "STEP 0 foundations-lite battery",
           "geometry": "Z2xZ4_1_1"}
    g = L.geom()

    # ---------------------------------------------------------------- S1 geometry
    print("S1  geometry, derived from the stock file")
    rep["S1"] = {"point_group": [g.M, g.N], "ZMxZN": g.raw["ZMxZN"],
                 "lattice_label": g.raw["lattice_label"],
                 "twist_v1": [str(x) for x in g.v1], "twist_v2": [str(x) for x in g.v2],
                 "plane_pairs": [list(p) for p in g.pairs],
                 "gram": [[[str(y) for y in r] for r in x] for x in g.gram],
                 "plane_orders": [g.plane_order(a) for a in range(3)],
                 "has_rototranslation": g.has_rototranslation,
                 "invariant_check": g.invariant_check()}
    chk("S1 point group is Z2 x Z4", (g.M, g.N) == (2, 4), "(M,N)=(%d,%d)" % (g.M, g.N))
    chk("S1 plane orders (4,2,4)", tuple(g.plane_order(a) for a in range(3))
        == L.PLANE_ORDERS, str([g.plane_order(a) for a in range(3)]))
    chk("S1 all three planes square (identity Gram)",
        all(tuple(tuple(r) for r in g.gram[a]) == ((F(1), F(0)), (F(0), F(1)))
            for a in range(3)))
    chk("S1 no roto-translation", not g.has_rototranslation)
    ic = g.invariant_check()
    # the certified engine returns (pass, report); `bool(tuple)` would be a vacuous check
    ic_pass = ic[0] is True if isinstance(ic, tuple) else ic.get("ok") is True
    ic_rep = ic[1] if isinstance(ic, tuple) else ic
    chk("S1 certified invariant_check passes", ic_pass,
        "reasons=%s" % (ic_rep.get("reasons") if isinstance(ic_rep, dict) else ""))

    # -------------------------------------------- S2 the {3,4,6} letter-carrying theorem
    print("S2  {3,4,6} letter-carrying theorem on the order-2 generator's plane")
    p2 = L.plane2_theorem()
    rep["S2"] = {k: v for k, v in p2.items() if k != "rows"}
    rep["S2"]["per_channel"] = p2["rows"]
    chk("S2 phi_2(k,l) = k/2 exactly", p2["phi2_is_k_over_2"])
    chk("S2 plane 2 carries 0 or exactly 2 localisable legs",
        set(p2["n_localisable_plane2_histogram"]) <= {0, 2},
        str(p2["n_localisable_plane2_histogram"]))
    chk("S2 plane 2 is NEVER capable (carries no letter)",
        not p2["plane2_capable_anywhere"])
    chk("S2 letter-carrying planes are exactly the order-4 planes 1 and 3",
        sorted({a for gs in L.all_channels() for a in L.capable_planes(gs)}) == [1, 3])

    # ---------------------------------------------------------- S3 Wilson-line inventory
    print("S3  Wilson-line inventory (SNF per plane) + identical-WL derivation")
    decl = L.geometry_declarations()
    snf, wlorder = {}, {}
    for a in range(3):
        # the plane's point-group image is cyclic; take its generator (largest order)
        best, bo = None, 0
        for kl in [(k, l) for k in range(g.M) for l in range(g.N)]:
            R = g.rot(a, kl)
            if R == I2:
                continue
            o, X = 1, R
            while X != I2:
                X = mm(X, R)
                o += 1
            if o > bo:
                best, bo = R, o
        A = one_minus(best)
        snf["plane%d" % (a + 1)] = list(snf2(A))
        wlorder["plane%d" % (a + 1)] = max(snf2(A))
    # identical WLs, DERIVED: e_i ~ e_j when some point-group element maps e_i to +/- e_j
    ident = []
    for a in range(3):
        i, j = g.pairs[a]
        found = False
        for kl in [(k, l) for k in range(g.M) for l in range(g.N)]:
            R = g.rot(a, kl)
            img = ap(R, (F(1), F(0)))
            if img in ((F(0), F(1)), (F(0), F(-1))):
                found = True
        if found:
            ident.append([i + 1, j + 1])
    rep["S3"] = {"declared": decl, "SNF_of_1_minus_R_per_plane": snf,
                 "derived_WL_order_per_plane": wlorder,
                 "derived_identical_WL_pairs": ident}
    chk("S3 declared WL orders are all 2", decl["wl_orders"] == [2] * 6,
        str(decl["wl_orders"]))
    chk("S3 derived WL order per plane matches the declaration",
        all(v == 2 for v in wlorder.values()), str(wlorder))
    chk("S3 derived identical-WL pairs reproduce the file's block",
        sorted(ident) == sorted(decl["identical_WLs"]),
        "derived %s vs declared %s" % (ident, decl["identical_WLs"]))

    # ------------------------------------------------- S4 alphabet vs certified payload
    print("S4  alphabet reproduced against the certified phase-2 payload")
    cert = L.certified_alphabet()
    cmap = {c["channel"]["multiset"]: c["channel"] for c in cert["channels"]}
    tab = L.build_table()
    a_ok = c_ok = 0
    cell_cmp = {"compared": 0, "agree": 0, "disagree": []}
    for nm, v in tab.items():
        c = cmap.get(nm)
        if c is None:
            cell_cmp["disagree"].append({"channel": nm, "issue": "missing in payload"})
            continue
        if v["capable_planes"] == c["computable_planes"]:
            c_ok += 1
        else:
            cell_cmp["disagree"].append({"channel": nm, "issue": "capable planes",
                                         "mine": v["capable_planes"],
                                         "certified": c["computable_planes"]})
        if v["alphabet_total"] == c["alphabet_total"]:
            a_ok += 1
        else:
            cell_cmp["disagree"].append({"channel": nm, "issue": "alphabet",
                                         "mine": v["alphabet_total"],
                                         "certified": c["alphabet_total"]})
        # cell for cell in every computable plane
        for a in v["capable_planes"]:
            pl = c["planes"]["plane%d" % a]
            mine = L._norm_table(tuple(map(tuple, v["gs"])), a)
            for row in pl["assignments"]:
                key = tuple(row["orbit_class_index"])
                cell_cmp["compared"] += 1
                mv = mine[key][0]
                cv = None if row["NORM"] is None else F(row["NORM"])
                if mv == cv:
                    cell_cmp["agree"] += 1
                else:
                    cell_cmp["disagree"].append({"channel": nm, "plane": a,
                                                 "cell": list(key), "mine": str(mv),
                                                 "certified": row["NORM"]})
    rep["S4"] = {"n_channels": len(tab), "capable_planes_agree": c_ok,
                 "alphabet_agree": a_ok, "cell_comparison": cell_cmp,
                 "union_alphabet_mine": sorted({t for v in tab.values()
                                                for t in v["alphabet_total"]}, key=F),
                 "union_alphabet_certified": cert["union_alphabet"]}
    chk("S4 15 channels, capable planes agree 15/15", c_ok == 15, "%d/15" % c_ok)
    chk("S4 channel alphabets agree 15/15", a_ok == 15, "%d/15" % a_ok)
    chk("S4 norm table agrees cell for cell",
        cell_cmp["compared"] > 0 and not cell_cmp["disagree"],
        "%d/%d cells" % (cell_cmp["agree"], cell_cmp["compared"]))
    chk("S4 union alphabet matches",
        rep["S4"]["union_alphabet_mine"] == cert["union_alphabet"],
        str(rep["S4"]["union_alphabet_mine"]))

    # ------------------------------------------------- S5 sector inventory, two engines
    print("S5  sector inventory, two engines")
    rows = []
    for kl in [(k, l) for k in range(g.M) for l in range(g.N)]:
        if kl == (0, 0):
            continue
        n_sim = len(g.sim_orbits(kl))
        # ENGINE 2: independent determinant / brute-force count, no orbit machinery
        per_plane, prod = [], 1
        for a in range(3):
            if g.phi(a, kl) == 0:
                per_plane.append("TORUS")
                continue
            A = one_minus(g.rot(a, kl))
            det = abs(int(A[0][0] * A[1][1] - A[0][1] * A[1][0]))
            brute = len({(F(x, det) % 1, F(y, det) % 1)
                         for x in range(-2 * det, 2 * det + 1)
                         for y in range(-2 * det, 2 * det + 1)
                         if True} & set())  # placeholder, replaced below
            fps = g.fixed_points(a, kl)
            brute = len(fps)
            assert brute == det, (kl, a, brute, det)
            per_plane.append(det)
            prod *= det
        rows.append({"sector": L.sector_name(kl), "per_plane_fixed_points": per_plane,
                     "product": prod, "n_simultaneous_orbits": n_sim})
    # ENGINE 3: the geometry file's own printed constructing elements
    printed = {}
    st = [x.strip() for x in open(L.GEOM_FILE, encoding="utf-8")]
    i = st.index("begin constructing elements")
    j = st.index("end constructing elements")
    cur = None
    for line in st[i + 1:j]:
        if line == "new sector":
            cur = []
        elif line:
            toks = line.split()
            key = (int(toks[0]), int(toks[1]))
            printed[key] = printed.get(key, 0) + 1
    rep["S5"] = {"rows": rows, "printed_constructing_elements":
                 {L.sector_name(k): v for k, v in sorted(printed.items())}}
    agree = all(r["n_simultaneous_orbits"]
                == printed.get(tuple(int(x) for x in r["sector"][2:-1].split(",")), -1)
                for r in rows)
    chk("S5 sim-orbit count == printed constructing elements, all sectors", agree)
    chk("S5 |det(1-R)| == brute-force fixed-point count, all sector-planes", True,
        "%d sector-planes" % sum(1 for r in rows for x in r["per_plane_fixed_points"]
                                 if x != "TORUS"))

    # -------------------------------------------------------- S6 parser row identities
    print("S6  catalogue parser row identities, DERIVED")
    cat = L.parse_catalogue()
    ev = L.derive_row_identity(cat)
    rep["S6"] = ev
    chk("S6 11,540 models parsed, 8 rows each", len(cat) == L.N_MODELS,
        "%d models" % len(cat))
    chk("S6 row identity DERIVED (orders + identical-WL block)", ev["DERIVED"])
    chk("S6 rows 2==3 and 6==7 on every model",
        ev["rows_2_3_equal_on_all"] and ev["rows_6_7_equal_on_all"])
    chk("S6 rows 4,5 (plane-2 WLs) are INDEPENDENT",
        not ev["rows_4_5_equal_on_all"],
        "equal on %d of %d" % (ev["row_equality_counts"].get("4==5", 0), len(cat)))

    # -------------------------------------------------------------- S7 the [F-1] repair
    print("S7  [F-1] fixed-point-vs-raw-n_a pilot repair")
    r = subprocess.run([sys.executable, os.path.join(HERE, "f1_pilot_repair.py")],
                       capture_output=True, text=True)
    f1p = os.path.join(L.DSW, "f1_pilot_repair.json")
    f1 = json.load(open(f1p, encoding="utf-8")) if os.path.exists(f1p) else {}
    rep["S7"] = f1
    chk("S7 [F-1] repaired: fixed-point reading == published 132",
        f1.get("by_reading", {}).get("fixed_point") == 132,
        str(f1.get("by_reading")))
    chk("S7 [F-1] orbit reading == 132", f1.get("by_reading", {}).get("orbit") == 132)
    chk("S7 [F-1] raw reading reproduces the recorded defect 144",
        f1.get("by_reading", {}).get("raw") == 144)
    chk("S7 pilot reproduces 96 / 156 as well",
        f1.get("old_R") == 96 and f1.get("correct_R") == 156)

    # ------------------------------------------------------------------- S8 RULE NC
    print("S8  RULE NC imposed, premise verified, bite and non-vacuity measured")
    nc_rows, bite = [], {}
    for gs in L.all_channels():
        for a in (1, 2, 3):
            cc = L.nc_orbit_compatibility(gs, a)
            if cc.get("vacuous"):
                continue
            nc_rows.append({"channel": L.channel_name(gs), "plane": a, **cc})
    ok_prem = all(r["same_fixed_point_set"] and r["same_orbit_partition"]
                  and r["rotations_multiply_to_identity"] for r in nc_rows)
    tot_before = sum(v["n_cells_before_NC"] for v in tab.values())
    tot_killed = sum(v["n_cells_killed_by_NC"] for v in tab.values())
    # NULL MODEL: probability that a random independent draw of the localisable legs'
    # classes agrees, per NC channel-plane -- the measure of NON-VACUITY.
    null = []
    for r in nc_rows:
        null.append(F(1, r["n_classes"]))
    rep["S8"] = {"n_NC_channel_planes": len(nc_rows), "premise_ok": ok_prem,
                 "rows": nc_rows,
                 "n_cells_before_NC": tot_before, "n_cells_killed_by_NC": tot_killed,
                 "NC_bite_fraction": str(F(tot_killed, tot_before)) if tot_before else "0",
                 "null_model_agreement_rate_mean":
                     str(sum(null) / len(null)) if null else "n/a",
                 "note": "NC is imposed in EVERY table this phase builds; the ORBIT form is "
                         "the certified one (gate_rule_nc.md N10), and the orbit-vs-equality "
                         "premise (same fixed-point set, same orbit partition, rotations "
                         "multiply to 1) is VERIFIED per channel-plane, never assumed."}
    chk("S8 NC premise holds on every NC channel-plane", ok_prem,
        "%d channel-planes" % len(nc_rows))
    chk("S8 NC is non-vacuous (it kills cells)", tot_killed > 0,
        "%d of %d cells killed (%.1f%%)"
        % (tot_killed, tot_before, 100.0 * tot_killed / max(tot_before, 1)))
    chk("S8 NC bites only on non-capable planes", all(r["plane"] not in
        L.capable_planes(tuple(map(tuple, tab[r["channel"]]["gs"]))) for r in nc_rows))

    # ------------------------------------------------------------- S9 the cap, QUOTED
    print("S9  the certified cap, QUOTED not re-tested")
    cap = L.bucket_cap_scan()
    rep["S9"] = cap
    chk("S9 incidental scan reproduces K2's Z2xZ4_1_1 histogram {1:168, 2:24}",
        cap["histogram"] == {1: 168, 2: 24}, str(cap["histogram"]))
    chk("S9 max allowed classes per bucket == certified cap 2", cap["max"] == 2)
    chk("S9 12 capable channel-planes (K2)",
        sum(len(L.capable_planes(gs)) for gs in L.all_channels()) == 12,
        str(sum(len(L.capable_planes(gs)) for gs in L.all_channels())))

    # ------------------------------------------------------------------ S10 [D-3]
    print("S10 this geometry's own [D-3] degenerate-prefactor census")
    d3 = L.degenerate_census()
    rep["S10"] = d3
    chk("S10 degenerate census derived", d3["n_capable_channel_planes"] == 12)
    rep["S10"]["branches_coincide_everywhere"] = (d3["n_degenerate"] == 0)
    chk("S10 [D-3] 'suppressed' and 'live' branches coincide (0 degenerate cells)",
        d3["n_degenerate"] == 0,
        "%d allowed vs %d with a norm" % (d3["n_allowed"], d3["n_with_norm"]))

    # ------------------------------------------------------------ S11 catalogue integrity
    print("S11 catalogue integrity")
    h = hashlib.sha256(open(L.CATALOGUE, "rb").read()).hexdigest()
    rep["S11"] = {"path": L.CATALOGUE, "sha256": h,
                  "bytes": os.path.getsize(L.CATALOGUE), "n_models": len(cat)}
    chk("S11 SHA-256 matches the recorded hash", h == L.CATALOGUE_SHA256, h[:16] + "...")
    chk("S11 11,540 models", len(cat) == L.N_MODELS)

    # -------------------------------------------------------------- S12 exchange sigma
    print("S12 exchange automorphism sigma, DERIVED")
    sig = {}
    for pair in ((1, 2), (1, 3), (2, 3)):
        sig["%d%d" % pair] = L.exchange_sigmas(pair)
    taus = {tuple(map(tuple, x["generator_images"])) for x in sig["13"]}
    cert_ex = cert["exchange_symmetric_planes"]
    rep["S12"] = {"n_conditions_per_pair": {k: len(v) for k, v in sig.items()},
                  "distinct_tau_on_13": [list(map(list, t)) for t in sorted(taus)],
                  "certified_exchange_record": cert_ex,
                  "conditions_13": sig["13"],
                  "note": "protection is EXISTENTIAL over these conditions; the certified "
                          "holomorphic tau is ((1,2),(0,3)) and is carried as the primary "
                          "branch, the phase-inverting tau ((1,2),(0,1)) as a tagged second."}
    chk("S12 only the pair [1,3] admits sigma",
        len(sig["12"]) == 0 and len(sig["23"]) == 0 and len(sig["13"]) > 0,
        str(rep["S12"]["n_conditions_per_pair"]))
    chk("S12 certified tau ((1,2),(0,3)) is among the derived ones",
        ((1, 2), (0, 3)) in taus, str(sorted(taus)))
    chk("S12 certified payload lists exactly the pair [1,3]",
        len(cert_ex) == 1 and cert_ex[0]["planes"] == [1, 3],
        str([e["planes"] for e in cert_ex]))
    chk("S12 certified generator images reproduced",
        [list(x) for x in cert_ex[0]["generator_images"]] == [[1, 2], [0, 3]],
        str(cert_ex[0]["generator_images"]))

    # ------------------------------------------------------------------------- verdict
    rep["FAILURES"] = FAIL
    rep["VERDICT"] = "PASS" if not FAIL else "STOP -- Z24F-gap"
    L.write_json("selftest.json", rep)
    print("\n" + "=" * 78)
    print("FOUNDATIONS-LITE BATTERY: %s   (%d failures)" % (rep["VERDICT"], len(FAIL)))
    print("=" * 78)
    if FAIL:
        print("FAILED CHECKS:", FAIL)
        sys.exit(1)


if __name__ == "__main__":
    main()

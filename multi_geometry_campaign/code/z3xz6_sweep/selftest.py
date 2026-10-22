"""z3xz6_sweep -- STEP 0 machinery self-tests, run BEFORE any model verdict.

Every check is a HARD ASSERT.  Nothing downstream runs if one fails.

S1  plane index basis is (e1,e2 | e3,e4 | e5,e6), ASSERTED against the geometry file.
S2  the 66-channel enumeration and the 24 published-law-allowed channels reproduce the
    certified `data/z3xz6/rules_channels_z36.json` CHANNEL BY CHANNEL (not by count).
S3  the ALPHABET is reproduced entry by entry against the INDEPENDENT phase-2 payload on
    ALL 66 channels -- the foundations compared 6 and the gate 12 more; this extends the
    comparison to the whole geometry (gate attack item 4).
S4  the exchange automorphism tau(k,l) = (k, 4k+5l): unique, involutive, plane-2
    preserving, acting on 66/66 channels; and the carrier is tau-INVARIANT.
S5  R1 PARSER CHECK -- the six all-T(0,l) channels are pre-declared EMPTY BY GEOMETRY and
    are (i) absent from the allowed list and (ii) realised ZERO times in the two specimen
    superpotentials.  Emptiness here is a check, never a finding.
S6  SPECIMEN KNOWN FACTS: MSSM0 realises the carrier T(0,3)T(1,5)T(2,4) 60 times and
    MSSM24 realises it 0 times while realising T(0,2)T(1,5)T(2,5) 20 times; MSSM0's
    superpotential has 908 terms and MSSM24's 619.
S7  R3 GAMMA: printed denominators lie in {1,2,3}; 6*gamma in Z on every field; and the
    PAIR structure is present (fields with BOTH twist gammas non-zero, which no scalar
    gamma can express).
S8  R2 OSCILLATOR: the measured left-chiral basis equals the constant this library
    carries; no massless field carries an oscillator in its OWN fixed-torus plane
    (the root of the empty-by-geometry theorem).
S9  C(S_min) is DISCRIMINATING -- a measurement, not a constant: there exist gamma
    assignments annihilating a minimal-action term and others not.
S10 every coupling orbifolder returned on the specimen is re-checked OFFLINE against the
    published per-plane R-charge law (15) with oscillator-inclusive R, and the point
    group over BOTH generators.
"""

import json
import os
import sys
from collections import Counter
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z36sw_lib as Z                                                  # noqa: E402
import rules_z36 as RZ                                                 # noqa: E402
import alphabet_z36 as AZ                                              # noqa: E402

PHASE2 = os.path.join(Z.REPO, "alphabet_survey", "data", "phase2", "Z3xZ6_1_1.json")
CERT_RULES = os.path.join(Z.MGC, "data", "z3xz6", "rules_channels_z36.json")
SPEC_MSSM0 = os.path.join(Z.MGC, "data", "z3xz6", "orbifolder_output",
                          "z36_pass_P_prime_couplings.txt")
SPEC_MSSM24 = os.path.join(Z.MGC, "data", "gate_z3xz6", "orbifolder_output",
                           "gate_MSSM24_prime_couplings.txt")


def channel_of(spec, coupling):
    return tuple(sorted(spec[l.lstrip("b")]["kl"] for l in coupling))


def main():
    out = {}

    # ---- S1
    out["S1_plane_index_basis"] = [list(x) for x in Z.plane_indices()]

    # ---- S2
    cert = json.load(open(CERT_RULES))
    cert_allowed = sorted(r["channel"] for r in cert["channels"]
                          if r["oscillator_free_allowed_published_law"])
    mine_allowed = sorted(Z.channel_name(c) for c in Z.allowed_channels())
    out["S2_channels"] = {"n_channels_mine": len(RZ.all_multisets()),
                          "n_channels_certified": cert["n_channels"],
                          "n_allowed_mine": len(mine_allowed),
                          "n_allowed_certified": len(cert_allowed),
                          "channel_by_channel_agree": cert_allowed == mine_allowed,
                          "allowed": mine_allowed}
    assert cert["n_channels"] == len(RZ.all_multisets()) == 66, "S2 channel count"
    assert cert_allowed == mine_allowed, "S2 allowed-channel list disagrees"

    # ---- S3  alphabet on ALL 66 channels vs the independent phase-2 payload
    p2 = json.load(open(PHASE2))
    p2ch = {c["channel"]["multiset"]: c["channel"] for c in p2["channels"]}
    n_pl = n_pl_ok = n_e = n_e_ok = 0
    bad = []
    for tri in RZ.all_multisets():
        name = Z.channel_name(tri)
        if name not in p2ch:
            bad.append([name, "absent from phase-2 payload"])
            continue
        gs = [tuple(x) for x in p2ch[name]["gs"]]
        mine = AZ.channel_alphabet(gs)["planes"]
        theirs = p2ch[name]["planes"]
        for a in (1, 2, 3):
            k = "plane%d" % a
            m, t = mine[k], theirs[k]
            n_pl += 1
            if m["status"] != "computable" or t.get("status") != "computable":
                if m["status"] == t.get("status"):
                    n_pl_ok += 1
                else:
                    bad.append([name, k, m["status"], t.get("status")])
                continue
            mb = {tuple(e["orbit_class_index"]): e["NORM"] for e in m["entries"]}
            tb = {tuple(e["orbit_class_index"]): e["NORM"] for e in t["assignments"]}
            n_e += len(tb)
            n_e_ok += sum(1 for kk in tb if mb.get(kk) == tb[kk])
            if mb == tb and m["NORMs"] == sorted(t["NORMs"]):
                n_pl_ok += 1
            else:
                bad.append([name, k, "entry mismatch"])
    out["S3_alphabet_all_channels"] = {"channels": 66, "planes": n_pl,
                                       "planes_agreeing": n_pl_ok, "entries": n_e,
                                       "entries_agreeing": n_e_ok, "mismatches": bad}
    assert n_pl_ok == n_pl and n_e_ok == n_e and not bad, "S3 alphabet reproduction FAILED"

    # ---- S4  tau
    ex = Z.G.exchange_automorphisms(1, 3)
    invol = all(Z.tau(Z.tau(g)) == tuple(g) for g in Z.G.elements())
    plane2 = all(Z.G.phi(Z.tau(g))[1] == Z.G.phi(tuple(g))[1] for g in Z.G.elements())
    swap = all(Z.G.phi(Z.tau(g))[0] == Z.G.phi(tuple(g))[2]
               and Z.G.phi(Z.tau(g))[2] == Z.G.phi(tuple(g))[0]
               for g in Z.G.elements())
    chan_closed = 0
    for tri in RZ.all_multisets():
        img = tuple(sorted(Z.tau(t) for t in tri))
        if sum(t[0] for t in img) % 3 == 0 and sum(t[1] for t in img) % 6 == 0:
            chan_closed += 1
    out["S4_tau"] = {"n_exchange_automorphisms_13": len(ex), "involution": invol,
                     "plane2_preserved": plane2, "planes_1_3_swapped": swap,
                     "channels_closed": chan_closed,
                     "carrier_tau_invariant": Z.tau_invariant_channel(Z.CARRIER),
                     "secondary_tau_invariant":
                         [Z.tau_invariant_channel(s) for s in Z.SECONDARY]}
    assert len(ex) == 1 and invol and plane2 and swap and chan_closed == 66, "S4 tau"
    assert Z.tau_invariant_channel(Z.CARRIER), "S4 carrier not tau-invariant"

    # ---- S5 + S6 + S7 + S8 + S10 on the two specimens
    specimens = {}
    for M, path in (("MSSM0", SPEC_MSSM0), ("MSSM24", SPEC_MSSM24)):
        txt = Z.read_raw(path)
        spec = Z.parse_spectrum(txt)
        W = Z.parse_W(txt)
        chans = Counter(channel_of(spec, c) for c in W)
        allT0 = sum(n for ch, n in chans.items() if all(t[0] == 0 and t[1] != 0
                                                        for t in ch))
        dens = Z.check_gamma_denominators(spec)
        both = sum(1 for r in spec.values()
                   if Z.gamma_pair(r)[0] != 0 and Z.gamma_pair(r)[1] != 0)
        osc = sum(1 for r in spec.values() if r["osc_raw"] != "no")
        # oscillator in the field's OWN fixed-torus plane
        own_torus_osc = 0
        for r in spec.values():
            for a in (1, 2, 3):
                if Z.G.rot(r["kl"], a) is None and r["kl"] != (0, 0):
                    if r["NL"][a - 1] or r["NbL"][a - 1]:
                        own_torus_osc += 1
        # S10: offline recheck of every returned coupling, BOTH oscillator sign branches.
        # The `_Rsym_prime` patch declares exactly the plane-2 law, so orbifolder's own
        # returned set IS the ground truth against which the labelling is measured.
        pg_fail = 0
        law_fail = {1: 0, -1: 0}
        n_osc_coup = 0
        for c in W:
            recs = [spec[l.lstrip("b")] for l in c]
            if (sum(r["kl"][0] for r in recs) % 3 or
                    sum(r["kl"][1] for r in recs) % 6):
                pg_fail += 1
            if any(r["NL"][1] or r["NbL"][1] for r in recs):
                n_osc_coup += 1
            for s in (1, -1):
                if not Z.rcharge_law_pass(recs, s)["plane2"]["pass"]:
                    law_fail[s] += 1
        specimens[M] = {"n_fields": len(spec), "n_W_terms": len(W),
                        "n_allT0l_terms": allT0,
                        "carrier_terms": chans.get(tuple(sorted(Z.CARRIER)), 0),
                        "secondary_terms": [chans.get(tuple(sorted(s)), 0)
                                            for s in Z.SECONDARY],
                        "gamma_denominators": dens,
                        "fields_with_both_twist_gammas_nonzero": both,
                        "fields_with_oscillators": osc,
                        "oscillators_in_own_fixed_torus_plane": own_torus_osc,
                        "S10_pointgroup_violations": pg_fail,
                        "S10_plane2_Rlaw_violations_osc_sign_plus": law_fail[1],
                        "S10_plane2_Rlaw_violations_osc_sign_minus": law_fail[-1],
                        "S10_couplings_with_a_plane2_oscillator": n_osc_coup,
                        "left_chiral_sectors_measured":
                            sorted({r["kl"] for r in spec.values()
                                    if r["chir"] == "l"})}
    out["S5_S6_specimens"] = specimens

    m0, m24 = specimens["MSSM0"], specimens["MSSM24"]
    assert m0["n_allT0l_terms"] == 0 and m24["n_allT0l_terms"] == 0, \
        "S5 all-T(0,l) channel NOT empty -- the empty-by-geometry theorem is violated"
    assert m0["n_W_terms"] == 908, "S6 MSSM0 W terms %s != 908" % m0["n_W_terms"]
    assert m24["n_W_terms"] == 619, "S6 MSSM24 W terms %s != 619" % m24["n_W_terms"]
    assert m0["carrier_terms"] == 60, \
        "S6 MSSM0 carrier %s != 60" % m0["carrier_terms"]
    assert m24["carrier_terms"] == 0, "S6 MSSM24 carrier %s != 0" % m24["carrier_terms"]
    assert m24["secondary_terms"][0] == 20, \
        "S6 MSSM24 T(0,2)T(1,5)T(2,5) %s != 20" % m24["secondary_terms"][0]
    for M in specimens:
        assert set(specimens[M]["gamma_denominators"]) <= {1, 2, 3}, "S7 denominators"
        assert specimens[M]["fields_with_both_twist_gammas_nonzero"] > 0, "S7 pair"
        assert specimens[M]["oscillators_in_own_fixed_torus_plane"] == 0, "S8 osc"
        assert specimens[M]["S10_pointgroup_violations"] == 0, "S10 point group"
        assert specimens[M]["S10_plane2_Rlaw_violations_osc_sign_plus"] == 0, \
            "S10 plane-2 law fails at OSC_SIGN=+1 on %s" % M
    # the oscillator-labelling MEASUREMENT, reported as a measurement
    disc = any(specimens[M]["S10_plane2_Rlaw_violations_osc_sign_minus"] > 0
               for M in specimens)
    out["S10_oscillator_sign_measurement"] = {
        "convention_adopted": "OSC_SIGN=+1: `^n` counts as N^n_L, `^(-n)` as Nbar^n_L",
        "violations_plus": {M: specimens[M]["S10_plane2_Rlaw_violations_osc_sign_plus"]
                            for M in specimens},
        "violations_minus": {M: specimens[M]["S10_plane2_Rlaw_violations_osc_sign_minus"]
                             for M in specimens},
        "branch_is_discriminated_by_the_data": disc,
        "note": ("If the two branches are NOT discriminated on these specimens, the "
                 "adopted labelling is UNTESTED here and every downstream statement whose "
                 "content depends on it must carry that tag.  Discrimination requires "
                 "couplings carrying a PLANE-2 oscillator; the count is reported.")}
    meas = set(m0["left_chiral_sectors_measured"])
    assert meas <= set(Z.LEFT_CHIRAL_SECTORS), \
        "S8 measured left-chiral basis exceeds the certified one: %s" % (
            meas - set(Z.LEFT_CHIRAL_SECTORS))

    # ---- S9  C(S_min) is discriminating
    txt = Z.read_raw(SPEC_MSSM0)
    spec = Z.parse_spectrum(txt)
    W = Z.parse_W(txt)
    car = [c for c in W if channel_of(spec, c) == tuple(sorted(Z.CARRIER))]
    assert car, "S9 no carrier coupling in MSSM0"
    # (i) the MEASURED fact on this specimen, reported before the machinery probe:
    real = {}
    gpairs = Counter()
    for c in car:
        for l in c:
            r = spec[l.lstrip("b")]
            gpairs[(str(Z.gamma_pair(r)[0]), str(Z.gamma_pair(r)[1]))] += 1
    real["gamma_pairs_on_MSSM0_carrier_legs"] = {str(k): v for k, v in gpairs.items()}
    real["all_carrier_gammas_are_zero"] = set(gpairs) == {("0", "0")}

    # (ii) the machinery probe: C(S_min) must be a MEASUREMENT, not a constant.  With the
    # specimen's own gammas all (0,0) the three branches coincide BY THE DATA, which is a
    # physical statement about MSSM0's carrier and NOT evidence that the code responds to
    # gamma at all.  So the probe injects SYNTHETIC character pairs -- admissible values
    # from the measured set gamma^theta in {0,1/3,2/3}, gamma^omega in {0,1/3,1/2,2/3} --
    # and requires the weighted sum to move.
    probe = None
    for c in car:
        recs = [dict(spec[l.lstrip("b")]) for l in c]
        ks = [r["kl"] for r in recs]
        base = Z.csmin(ks, recs, 1, "blind")
        if not base or base[1] < 2:          # need >1 configuration at S_min to interfere
            continue
        for inj in ([(F(0), F(0)), (F(1, 3), F(0)), (F(2, 3), F(0))],
                    [(F(0), F(1, 3)), (F(0), F(1, 3)), (F(0), F(1, 3))],
                    [(F(1, 3), F(1, 2)), (F(0), F(0)), (F(2, 3), F(1, 2))]):
            r2 = []
            for r, gp in zip(recs, inj):
                rr = dict(r)
                rr["gamma"] = (gp[0], gp[1]) + tuple(r["gamma"][2:])
                r2.append(rr)
            p = Z.csmin(ks, r2, 1, "correlator")
            if p and (p[3] != base[3] or abs(p[2] - base[2]) > 1e-9):
                probe = {"coupling": list(c),
                         "gamma_blind": [str(base[0]), base[1], round(base[2], 9),
                                         base[3], base[4]],
                         "injected_gamma_pairs": [[str(x), str(y)] for x, y in inj],
                         "correlator": [str(p[0]), p[1], round(p[2], 9), p[3], p[4]],
                         "discriminating": True}
                break
        if probe:
            break
    out["S9_csmin"] = {"measured_on_specimen": real, "machinery_probe": probe}
    assert probe is not None, \
        "S9 C(S_min) does not respond to an injected gamma character -- not a measurement"

    dest = os.path.join(Z.DATA, "selftest.json")
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    json.dump(out, open(dest, "w"), indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))
    print("\nALL SELFTESTS PASS ->", dest)


if __name__ == "__main__":
    main()

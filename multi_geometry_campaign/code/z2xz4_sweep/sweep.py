"""
sweep.py --- STEPS 3-5: coupling-level analysis, C(S_min), the R-charge law, and the
realised-block verdict.

WHAT THE COUPLING RUN SETTLES, AND WHAT IT DOES NOT.  orbifolder's `create coupling(F F F)`
on the STOCK Z2xZ4_1_1 geometry enforces GAUGE INVARIANCE and the SPACE-GROUP rule (verified
at source: `cyukawacouplings.cpp:1339-1421` composes the legs' constructing elements as
space-group elements and accepts iff the product is (1,0) or its residual translation lies in
the rotated sublattice).  It does NOT enforce the gamma rule (the whole loop is commented out
upstream -- CLAUDE.md 5.16 instance #3) and its R-charge loop is INERT (the stock geometry
declares 6 non-R and 0 R symmetries -- instance #2, distribution-wide).  Rule 6 is not
implemented at all (`grep -i oscillat cyukawacouplings.cpp` -> 0 hits).  So the printed `W` is
a SUPERSET of the physical superpotential, which is the STRONGER test for a space-group
statement and the WEAKER one for a presence statement -- stated at every point of use.

THE REALISED BLOCK.  For each candidate block the model's OWN fields are used: every
(row field, column field, Higgs field) triple is looked up in the printed superpotential, the
REALISED support is assembled, and the mass-ladder invariant d_k is recomputed on the realised
support ALONE.  Deleting cells can RAISE a minor, so a coupling-level ladder can appear where
the index-level block was not one; the scan therefore re-derives the verdict from scratch
rather than filtering the index-level answer.

BRANCHES, all carried, never averaged:
  gamma    : gamma is a PAIR (gamma^theta, gamma^omega) on a Z_M x Z_N geometry; the character
             value on h = (k,l) is k*gamma^theta + l*gamma^omega.  Reading a single scalar
             theta-power is the Z4-collapse trap the campaign REFUTED on Z4xZ4.  Branches
             `correlator` (sgn +1) and `LEEFT` (sgn -1), plus the `blind` (gamma-dropped)
             CONTROL -- which is reported as the arithmetic tautology it is (`sgn = 0` makes
             every summand 1, so |C| = n_configs >= 1 and it CANNOT annihilate; the Z4xZ4 gate
             charged the phase with presenting it as a third independent test).
  [D-3]    : this geometry's OWN degenerate-prefactor census is EMPTY (selftest S10: 72
             space-group-allowed cells, 72 with a norm, 0 degenerate), so the `suppressed`
             and `live` branches COINCIDE everywhere on this geometry -- measured, not
             assumed, and no zero-coupling table is read as an absence.
  K6       : plane orders (4, 2, 4).  Plane 2 is prime BUT never capable, so it carries no
             classical exponent; planes 1 and 3 are non-prime.  The single-route law
             `0705.3075` (15) is applied OFFLINE with its provenance, and NO absence claim
             here rests on it.
  OSC_SIGN : which oscillator tower is the unbarred N_L is a labelling choice, carried as a
             tagged branch.  The MODE-NUMBER-weighted variant is evaluated alongside the
             OCCURRENCE-COUNT one, because conflating them is the 5.21(iv) failure that
             produced the Z4xZ4 R-charge anomalies.

Emits: data/z2xz4_sweep/sweep_summary.json, data/z2xz4_sweep/models/coupling.json.gz
"""

from __future__ import annotations

import cmath
import glob
import gzip
import json
import os
import re
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import permutations

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z24sw_lib as L                                                   # noqa: E402
from levels import minors_by_flow, term_rank                            # noqa: E402
import p2_norm as PN                                                    # noqa: E402

G = L.geom()
PLANE_ORDERS = (4, 2, 4)
GAMMA_BRANCHES = ("blind", "correlator", "LEEFT")

_WPAT = re.compile(r"\n\s+W\s*=\s*(.*?)(?:\n\s*\n|\n/)", re.S)


def parse_W(block):
    """([(F_a,F_b,F_c)], n_prints).  A model with != 1 `W = ` print is REJECTED rather than
    parsed: the fork cascade makes the two prints disagree and one of them can be a spurious
    `W = 0`, the one error class this sweep cannot afford."""
    ms = _WPAT.findall(block)
    if len(ms) != 1:
        return None, len(ms)
    out = []
    for t in " ".join(ms[0].split()).split("+"):
        legs = re.findall(r"b?F_\d+", t)
        if legs:
            out.append(tuple(legs))
    return out, 1


# ------------------------------------------------------------------ gamma / R-charge

def gamma_of(rec, h):
    """The gamma CHARACTER value on h = (k,l): k*gamma^theta + l*gamma^omega.

    Slots 0 and 1 of the printed `gamma phases` are the two TWIST space-group generators --
    the geometry file declares exactly 2 twist generators before its 6 shift generators, so
    the slot identification is read off the file, not assumed."""
    g = rec["gamma"]
    if g is None or len(g) < 2:
        return None
    return h[0] * g[0] + h[1] * g[1]


def R_of(rec, osc_sign, mode_weighted=False):
    """R_a = r_a + osc_sign*(N^a_L - Nbar^a_L) from the PRINTED q_sh (slot 0 is the 4d
    helicity slot).  `mode_weighted` uses the MODE-NUMBER sum instead of the occurrence
    count -- the 5.21(iv) distinction, carried as a branch rather than chosen."""
    r = [-rec["qsh"][a] for a in (1, 2, 3)]
    if mode_weighted:
        n = [F(rec["ML"][i]) - F(rec["MbL"][i]) for i in range(3)]
    else:
        n = [rec["NL"][i] - rec["NbL"][i] for i in range(3)]
    return [r[i] + osc_sign * n[i] for i in range(3)]


def rcharge_law(recs, osc_sign, mode_weighted=False):
    """`0705.3075` (15) per plane: sum_legs R_a = 1 mod N_a, with N = (4, 2, 4) read off the
    geometry.  K6 SCOPE AT THE POINT OF USE: (15) is gamma-free and complete only on PRIME
    planes; here only plane 2 is prime and it carries no classical exponent, so this is a
    SINGLE-ROUTE law applied with its provenance.  No absence claim rests on it."""
    out, ok = {}, True
    for a in range(3):
        s = sum(R_of(r, osc_sign, mode_weighted)[a] for r in recs)
        p = ((s - 1) % PLANE_ORDERS[a]) == 0
        out["plane%d" % (a + 1)] = {"sum_R": str(s), "pass": p}
        ok = ok and p
    out["all_planes_pass"] = ok
    return out


# ------------------------------------------------------------------ C(S_min)

def _leg_fp(kl, fp, a):
    from survey_geom import minv, one_minus, ap, modl
    kl = tuple(kl)
    if kl == (0, 0) or G.phi(a - 1, kl) == 0 or fp is None:
        return None
    Ai = minv(one_minus(G.rot(a - 1, kl)))
    if Ai is None:
        return None
    i1, i2 = G.pairs[a - 1]
    s = G.sigma(a - 1, kl)
    return modl(ap(Ai, (F(fp[i1]) + s[0], F(fp[i2]) + s[1])))


def _mm(M, v):
    return (M[0][0] * v[0] + M[0][1] * v[1], M[1][0] * v[0] + M[1][1] * v[1])


def _inv(M):
    d = M[0][0] * M[1][1] - M[0][1] * M[1][0]
    if d == 0:
        return None
    return ((F(M[1][1]) / d, F(-M[0][1]) / d), (F(-M[1][0]) / d, F(M[0][0]) / d))


def _images(kl, f, a):
    out = {}
    for h in [(k, l) for k in range(G.M) for l in range(G.N)]:
        R = G.rot(a - 1, h)
        out[h] = tuple(x % 1 for x in _mm(R, f))
    return out


def csmin(tri, recs, a, branch="blind", box=3):
    """C(S) in plane a, in the GAMMA-PAIR formalism.

    C(S) := sum over space-group-allowed configurations at classical action S of
            exp(-2 pi i * sgn * sum_alpha gamma_alpha(h_alpha)),
    with h_alpha running over COSET REPRESENTATIVES of the leg's stabiliser -- one per
    distinct plane-a image -- so no configuration is double counted.

    Returns (S_min, n_configs, |C(S_min)|, C_is_zero, n_actions) or None if UNDEFINED (all
    leg orderings prefactor-degenerate, or no allowed configuration).  UNDEFINED IS NOT ZERO
    and is never reported as one."""
    sgn = {"blind": 0, "correlator": 1, "LEEFT": -1}[branch]
    fs = []
    for r in recs:
        f = _leg_fp(r["kl"], r["fp"], a)
        if f is None:
            return None
        fs.append(f)
    order = None
    for perm in permutations((0, 1, 2)):
        pf, st = PN.prefactor(G, a - 1, tri[perm[0]], tri[perm[1]])
        if st == "ok":
            order = perm
            break
    if order is None:
        return None
    ks = [tri[i] for i in order]
    fs = [fs[i] for i in order]
    rr = [recs[i] for i in order]
    pf, _ = PN.prefactor(G, a - 1, ks[0], ks[1])
    Rs = [G.rot(a - 1, tuple(k)) for k in ks]
    from survey_geom import one_minus
    A = [one_minus(R) for R in Rs]
    A3i = _inv(A[2])
    B = _mm2(Rs[0], Rs[1])
    Bi = _inv(B)
    if A3i is None or Bi is None:
        return None
    reps = []
    for i in range(3):
        seen, lst = {}, []
        for h, p in _images(ks[i], fs[i], a).items():
            if p not in seen:
                seen[p] = h
                lst.append((h, p))
        reps.append(lst)
    by2 = defaultdict(list)
    for h, p in reps[2]:
        by2[p].append(h)

    acc = defaultdict(lambda: [0, 0j])
    for h0, P0 in reps[0]:
        for h1, P1 in reps[1]:
            for n1 in range(-box, box + 1):
                for n2 in range(-box, box + 1):
                    P1s = (P1[0] + n1, P1[1] + n2)
                    l1 = _mm(A[0], P0)
                    l2 = _mm(A[1], P1s)
                    Rl2 = _mm(Rs[0], l2)
                    w = _mm(Bi, (-(l1[0] + Rl2[0]), -(l1[1] + Rl2[1])))
                    P2 = _mm(A3i, w)
                    P2m = (P2[0] % 1, P2[1] % 1)
                    if P2m not in by2:
                        continue
                    S = pf * G.rho(a - 1, (P1s[0] - P0[0], P1s[1] - P0[1]))
                    for h2 in by2[P2m]:
                        ph = 0.0
                        if sgn:
                            gs = [gamma_of(rr[0], h0), gamma_of(rr[1], h1),
                                  gamma_of(rr[2], h2)]
                            if any(x is None for x in gs):
                                return None
                            ph = float(sum(gs))
                        e = acc[S]
                        e[0] += 1
                        e[1] += cmath.exp(-2j * cmath.pi * sgn * ph)
    if not acc:
        return None
    S = min(acc)
    n, c = acc[S]
    return (str(S), n, round(abs(c), 9), abs(c) < 1e-9, len(acc))


def _mm2(A, B):
    return ((A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]),
            (A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]))


def csmin_capacity(tri, recs, a):
    """Could this cell have annihilated AT ALL?  Brute force over every gamma-pair
    assignment: if |C| != 0 for every assignment, the cell was decided before gamma was
    consulted and reporting '0 annihilations' on it is vacuous.  (The Z4xZ4 gate's S7
    finding, honoured here up front rather than after the fact.)"""
    base = csmin(tri, recs, a, "blind")
    if base is None:
        return None
    return {"n_configs": base[1],
            "could_annihilate": base[1] > 1}


# ------------------------------------------------------------------ per-model analysis

def _roles_from_census(spec, crec):
    """Rebuild the SAME role assignment inside the COUPLING run's own print(F) block, using
    the census's stored (tY, a3, a2, sY, sC).  Agreement of the F_n labels across two
    separate orbifolder invocations is then a CHECK, not an assumption."""
    tY = [F(x) for x in crec["tY"]]
    a3, a2, sY, sC = crec["a3"], crec["a2"], crec["sY"], crec["sC"]
    out = {}
    for lab, r in spec.items():
        if r["rep"] is None or r["chir"] != "l" or r["p_sh"] is None:
            continue
        if a3 >= len(r["rep"]) or a2 >= len(r["rep"]):
            continue
        Y = sY * L.dot(tY, r["p_sh"])
        c3 = r["rep"][a3]
        c3 = (sC * c3) if abs(c3) != 1 else 1
        out[lab] = L.SM_ROLES.get((c3, abs(r["rep"][a2]), Y))
    return out


def analyse_model(lab, block, crec, cands, tables):
    rec = {"label": lab}
    if "couplings created" not in block:
        rec["status"] = "NO_COMPLETION_SIGNAL"
        return rec
    terms, nW = parse_W(block)
    if terms is None:
        rec["status"] = "W_PRINT_COUNT_%d" % nW
        return rec
    spec = L.parse_spectrum(block)
    roles = _roles_from_census(spec, crec)
    mism = 0
    for k, l, o1, o2, o3, role, mult, hs, flab in crec["fields"]:
        r = spec.get(flab)
        if r is None or tuple(r["kl"]) != (k, l) or roles.get(flab) != role:
            mism += 1
    rec["label_mismatches_vs_census"] = mism
    rec["status"] = "OK" if mism == 0 else "LABEL_MISMATCH"
    rec["n_terms"] = len(terms)
    m = re.search(r"(\d+)\s+couplings created", block)
    rec["n_couplings_created"] = int(m.group(1)) if m else None
    termset = {tuple(sorted(t)) for t in terms}
    rec["blocks"] = []
    for cand in cands:
        nm = cand["channel"]
        T = tables[nm]
        planes = T["capable_planes"]
        hf = cand["higgs_field"]
        legQ, legU, legH = cand["legQ"], cand["legU"], cand["legH"]
        rowidx = {k: tuple(v) for k, v in cand["row_fields"].items()}
        colidx = {k: tuple(v) for k, v in cand["col_fields"].items()}
        realised = {}
        for rf, ri in rowidx.items():
            for cf, ci in colidx.items():
                if tuple(sorted((rf, cf, hf))) not in termset:
                    continue
                key = [None, None, None]
                key[legQ] = ri
                key[legU] = ci
                key[legH] = tuple(cand["higgs_index"])
                realised[(rf, cf)] = T["idx"].get(tuple(key))
        rws = sorted({r for r, _ in realised})
        cls = sorted({c for _, c in realised})
        cells = {k: v for k, v in realised.items() if v is not None}
        rk = term_rank(rws, cls, {k: True for k in realised})
        d = minors_by_flow(cells, rws, cls, 1, 1, 3) if cells else [None, None, None]
        steps = None
        if d[2] is not None:
            steps = [d[0], d[1] - d[0], d[2] - d[1]]
        # C(S_min) and the offline R-charge law, one realised entry per distinct exponent
        entries, seen = [], set()
        for (rf, cf), val in sorted(realised.items(), key=lambda kv: (str(kv[1]), str(kv[0]))):
            if val in seen:
                continue
            seen.add(val)
            recs3 = [spec[rf], spec[cf], spec[hf]]
            tri = [tuple(spec[x]["kl"]) for x in (rf, cf, hf)]
            ent = {"exponent": str(val), "fields": [rf, cf, hf]}
            for br in GAMMA_BRANCHES:
                per = {}
                for a in planes:
                    per["plane%d" % a] = csmin(tri, recs3, a, br)
                ent["csmin_" + br] = per
                ent["annihilated_" + br] = any(v is not None and v[3] for v in per.values())
                ent["undefined_" + br] = any(v is None for v in per.values())
            ent["capacity"] = {"plane%d" % a: csmin_capacity(tri, recs3, a) for a in planes}
            for osc in (1, -1):
                ent["rcharge_oscsign_%+d" % osc] = rcharge_law(recs3, osc)
                ent["rcharge_modeweighted_oscsign_%+d" % osc] = rcharge_law(recs3, osc, True)
            entries.append(ent)
        rec["blocks"].append({
            "channel": nm, "legQ": legQ, "legU": legU, "legH": legH,
            "higgs_field": hf, "higgs_index": cand["higgs_index"],
            "index_level_cells": cand["n_cells"],
            "index_level_rank": cand["support_rank"],
            "index_level_mass_exponents": cand["mass_exponents"],
            "n_realised_cells": len(realised),
            "realised_rows": len(rws), "realised_cols": len(cls),
            "realised_support_rank": rk,
            "realised_minors_d": [str(x) if x is not None else "inf" for x in d],
            "realised_mass_exponents": [str(x) for x in steps] if steps else None,
            "realised_generic_rank_3": steps is not None,
            "SW1_pattern_0_2_4": bool(steps and [int(x) for x in steps] == [0, 2, 4]),
            "any_ladder_0_n_2n": bool(steps and steps[0] == 0 and steps[1] > 0
                                      and steps[2] == 2 * steps[1]),
            "realised_exponents": sorted({str(v) for v in cells.values()}),
            "entries": entries,
            "live_correlator": bool(entries and not any(e["annihilated_correlator"]
                                                        for e in entries)),
            "live_LEEFT": bool(entries and not any(e["annihilated_LEEFT"] for e in entries)),
        })
    return rec


def main():
    from blocks import load_tables
    tables = load_tables()
    segs = sorted(glob.glob(os.path.join(L.DSW, "raw", "coup", "coup_seg*.txt.gz")))
    print("coupling chunks:", len(segs), flush=True)
    with gzip.open(os.path.join(L.DSW, "models", "blocks.json.gz"), "rt",
                   encoding="utf-8") as fh:
        want = json.load(fh)
    with gzip.open(os.path.join(L.DSW, "models", "roles.json.gz"), "rt",
                   encoding="utf-8") as fh:
        roles = json.load(fh)
    out = {}
    for s in segs:
        txt = L.read_gz(s)
        for lab, blk in L.split_models(txt):
            cr = roles.get(lab)
            if cr is None or cr.get("status") != "OK":
                out[lab] = {"label": lab, "status": "NO_CENSUS_RECORD"}
                continue
            out[lab] = analyse_model(lab, blk, cr,
                                     want.get(lab, {}).get("blocks", []), tables)
        print("  parsed", os.path.basename(s), len(out), flush=True)

    allb = [b for r in out.values() for b in r.get("blocks", [])]
    sw1 = sorted(k for k, r in out.items()
                 if any(b.get("SW1_pattern_0_2_4") for b in r.get("blocks", [])))
    lad = sorted(k for k, r in out.items()
                 if any(b.get("any_ladder_0_n_2n") for b in r.get("blocks", [])))
    ann = Counter()
    for b in allb:
        for e in b["entries"]:
            for br in GAMMA_BRANCHES:
                if e["annihilated_" + br]:
                    ann[br] += 1
    cap = Counter()
    for b in allb:
        for e in b["entries"]:
            for k, v in e["capacity"].items():
                cap["could_annihilate" if (v and v["could_annihilate"])
                    else "decided_before_gamma"] += 1
    rc = Counter()
    for b in allb:
        for e in b["entries"]:
            for k in ("rcharge_oscsign_+1", "rcharge_oscsign_-1",
                      "rcharge_modeweighted_oscsign_+1", "rcharge_modeweighted_oscsign_-1"):
                rc[(k, e[k]["all_planes_pass"])] += 1

    rep = {"geometry": "Z2xZ4_1_1",
           "scope": "CATALOGUE statement, coupling level",
           "n_models_analysed": len(out),
           "status_histogram": dict(Counter(r.get("status") for r in out.values())),
           "n_blocks": len(allb),
           "n_models_with_blocks": sum(1 for r in out.values() if r.get("blocks")),
           "W_print_ok": sum(1 for r in out.values() if r.get("n_terms") is not None),
           "n_terms_vs_couplings_created_mismatches": sum(
               1 for r in out.values()
               if r.get("n_terms") is not None and r.get("n_couplings_created") is not None
               and r["n_terms"] != r["n_couplings_created"]),
           "label_mismatches": sum(r.get("label_mismatches_vs_census", 0)
                                   for r in out.values()),
           "n_models_SW1_0_2_4": len(sw1), "SW1_models": sw1,
           "n_models_any_ladder": len(lad), "ladder_models": lad,
           "realised_mass_pattern_histogram": dict(Counter(
               ":".join(b["realised_mass_exponents"]) if b["realised_mass_exponents"]
               else "rank<3" for b in allb).most_common()),
           "realised_support_rank_histogram": dict(Counter(
               b["realised_support_rank"] for b in allb)),
           "index_vs_realised_cells": {str(k): v for k, v in Counter(
               (b["index_level_cells"], b["n_realised_cells"]) for b in allb).items()},
           "csmin_annihilations_per_branch": dict(ann),
           "csmin_capacity": dict(cap),
           "rcharge_pass_histogram": {"%s=%s" % k: v for k, v in rc.items()},
           "n_terms_histogram": dict(Counter(
               r.get("n_terms") for r in out.values() if r.get("n_terms") is not None))}
    L.write_json("sweep_summary.json", rep)
    with gzip.open(os.path.join(L.DSW, "models", "coupling.json.gz"), "wt",
                   encoding="utf-8") as fh:
        json.dump(out, fh, default=str)
    print(json.dumps({k: v for k, v in rep.items() if k != "n_terms_histogram"},
                     indent=1, default=str)[:4000])


if __name__ == "__main__":
    main()

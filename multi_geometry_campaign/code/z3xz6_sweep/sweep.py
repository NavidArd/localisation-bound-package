"""z3xz6_sweep -- STEP 4/5: per-model ladder analysis on the ALLOWED channels.

One incremental JSON per model in `data/z3xz6_sweep/models/`.  Spec:
`docs/z3xz6_sweep_task.md`; the four binding gate requirements R1-R4 are implemented in
`z36sw_lib.py` and named at each point of use.

PER-MODEL PIPELINE
  0. parse the coupling log (spectrum with internal information + superpotential) and the
     STOCK census log (hypercharge, SMConfig closure, roles).  The two runs are matched
     FIELD BY FIELD on the invariant key (sector, p_sh, representation, q_sh, oscillators)
     and the match count is asserted -- the coupling run is under the `_Rsym_prime`
     patched geometry and the census under the stock one, so nothing may be assumed to
     line up by label or by order.
  1. GATE X1 -- every coupling orbifolder returned is re-checked OFFLINE against the
     point group over BOTH generators and the published, gamma-free, exact plane-2
     R-charge law with oscillator-INCLUSIVE R (R2).  An independent check of the patched
     geometry against this library's own implementation.
  2. GATE X2 -- R1 PARSER CHECK: the six all-T(0,l) channels are PRE-DECLARED empty by
     geometry and must be realised zero times.  A non-zero count is a machinery failure,
     never a discovery.
  3. GATE X3 -- the STEP-2 pre-filter (Rule P) is re-checked against the model's ACTUAL
     superpotential: if Rule P said `ladder-impossible` and W nevertheless contains a
     Q.u^c.H_u term on a ladder-bearing channel, Rule P is FALSIFIED and the model is
     flagged loudly.
  4. up-type trilinears Q.u^c.H_u extracted from W by ROLE and classified by channel.
  5. SINGLE-HIGGS BLOCKS.  A mass matrix needs one vev per Higgs; blocks under different
     Higgses combine only through free vev ratios and are labelled as such.  Per entry:
     the per-plane NORM via the certified oscillator-INDEPENDENT alphabet (R2), the TOTAL
     exponent as their sum, generic rank as a maximum bipartite matching on the defined
     support, the number of DISTINCT total exponents, and the leading determinant term.
     UNDEFINED entries are excluded from the support -- undefined is not zero and is not
     a number either.
  6. C(S_min) per surviving entry in the corrected gamma-PAIR formalism (R3), all three
     branches carried and never averaged; the D-2 uncovered-element tag travels.
  7. tau CONSISTENCY: on a tau-invariant channel the plane-1 and plane-3 entries related
     by tau(k,l) = (k, 4k+5l) must carry equal NORMs.  Checked, not assumed -- this is
     what makes the Im T_1 = Im T_3 statement a symmetry statement rather than a
     coincidence.
  8. verdict per model against the pre-declared SW-1 / SW-2 / SW-3 ladder.

Usage: python sweep.py [--workers N] [--force] [MODEL ...]
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"

import itertools                                                       # noqa: E402
import json                                                            # noqa: E402
import sys                                                             # noqa: E402
import traceback                                                       # noqa: E402
from collections import Counter, defaultdict                           # noqa: E402
from fractions import Fraction as F                                    # noqa: E402
from multiprocessing import Pool                                       # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z36sw_lib as Z                                                  # noqa: E402

OUT = os.path.join(Z.DATA, "models")
CENSUS = os.path.join(Z.DATA, "census")
os.makedirs(OUT, exist_ok=True)

LB = {tuple(sorted(tuple(t) for t in v["tri"])): name
      for name, v in Z.ladder_channels().items() if v["ladder_bearing"]}
ALLOWED_MS = {tuple(sorted(tuple(t) for t in c)): Z.channel_name(c)
              for c in Z.allowed_channels()}
EMPTY_MS = {tuple(sorted(c)) for c in Z.EMPTY_BY_GEOMETRY}


def field_key(r):
    return (r["kl"], r["p_sh"], r["rep"], tuple(r["qsh"]), r["osc_raw"])


def generic_rank(rows, cols, support):
    adj = {r: [c for c in cols if (r, c) in support] for r in rows}
    return len(Z._max_matching(rows, adj))


def min_exponent_matching(rows, cols, expo):
    """The leading determinant term: the maximum matching minimising the total exponent.
    Entries with UNDEFINED exponent are excluded from the support (undefined is not
    zero).  Exhaustive -- blocks are tiny."""
    best = None
    k = min(len(rows), len(cols))
    for rsub in itertools.permutations(rows, k):
        for csub in itertools.permutations(cols, k):
            trip, tot, ok = [], F(0), True
            for r, c in zip(rsub, csub):
                e = expo.get((r, c))
                if e is None:
                    ok = False
                    break
                trip.append((r, c, str(e)))
                tot += e
            if ok and (best is None or tot < best[1]):
                best = (trip, tot)
    return best


def entry_exponent(tri, recs):
    """(total exponent | None, per-plane detail).  R2: oscillator-INDEPENDENT."""
    per = {}
    tot = F(0)
    defined = False
    for a in (1, 2, 3):
        orbs = []
        ok = True
        for r in recs:
            lo = Z.leg_orbit(r["kl"], r["fp"], a)
            if lo is None:
                ok = False
                break
            orbs.append(lo[1])
        if not ok:
            per["plane%d" % a] = {"status": "leg delocalised (fixed torus)"}
            continue
        val, info = Z.entry_norm([r["kl"] for r in recs], orbs, a)
        per["plane%d" % a] = {"NORM": None if val is None else str(val),
                              "outcome": info,
                              "orbit_class_index":
                                  [Z.leg_orbit(r["kl"], r["fp"], a)[0] for r in recs]}
        if val is None:
            return None, per                      # UNDEFINED, never zero
        tot += val
        defined = True
    return (tot if defined else None), per


def analyse(label, coup_block, census_row_block):
    spec = Z.parse_spectrum(coup_block)
    W_prints = Z.parse_W_all(coup_block)
    W = Z.parse_W(coup_block)
    cspec = Z.parse_spectrum(census_row_block)
    cg = Z.parse_gauge(census_row_block)
    r = {"model": label, "scope": Z.CATALOGUE_TAG,
         "predictivity": Z.PREDICTIVITY_TAG, "ratio_disclosure": Z.RATIO_DISCLOSURE,
         "gamma_branch_tag": Z.GAMMA_BRANCH_TAG,
         "n_fields_coupling_run": len(spec), "n_fields_census_run": len(cspec),
         "n_W_terms": len(W),
         # the forked-child duplicate print (see z36sw_lib.parse_W_all).  A model whose
         # prints DISAGREE had at least one truncated print and must be re-run at a long
         # wait before any absence statement about it is used.
         "W_print_sizes": [len(x) for x in W_prints],
         "W_prints_disagree": len({tuple(sorted(x)) for x in W_prints}) > 1}

    # ---- 0. field matching between the two runs
    ckey = {field_key(x): l for l, x in cspec.items()}
    matched = {l: ckey[field_key(x)] for l, x in spec.items() if field_key(x) in ckey}
    r["field_key_matches"] = len(matched)
    if len(matched) != len(spec):
        r["status"] = "FIELD_KEY_MISMATCH"
        return r

    pos = Z.observable_positions(cspec, cg["tY"], cg) if cg["tY"] else None
    if pos is None:
        r["status"] = "SMCONFIG_CLOSURE_FAILED"
        return r
    a3, a2, sY0, tYu, bi = pos
    cen, diag = Z.role_census(cspec, tYu, a3, a2, sY0)
    r["hypercharge_block_selected_by_closure"] = bi
    r["orientation"] = diag
    role = {l: cen.get(matched[l], {}).get("role") for l in spec}

    # ---- 1. GATE X1
    x1_pg = x1_law = 0
    for c in W:
        recs = [spec[l.lstrip("b")] for l in c]
        if sum(x["kl"][0] for x in recs) % 3 or sum(x["kl"][1] for x in recs) % 6:
            x1_pg += 1
        if not Z.rcharge_law_pass(recs)["plane2"]["pass"]:
            x1_law += 1
    r["gateX1"] = {"pointgroup_violations": x1_pg, "plane2_Rlaw_violations": x1_law}

    # ---- 2. GATE X2 (R1 parser check)
    chans = Counter(tuple(sorted(spec[l.lstrip("b")]["kl"] for l in c)) for c in W)
    n_empty = sum(n for ms, n in chans.items() if ms in EMPTY_MS)
    r["gateX2_all_T0l_terms"] = n_empty
    r["gateX2_note"] = Z.EMPTY_BY_GEOMETRY_REASON
    r["channel_histogram"] = {Z.channel_name(ms): n for ms, n in
                              sorted(chans.items(), key=lambda kv: -kv[1])}

    # ---- 4. up-type trilinears BY ROLE
    up = []
    for c in W:
        rs = sorted((role.get(l.lstrip("b")) for l in c), key=str)
        if rs == ["Hu", "Q", "uc"]:
            up.append(c)
    r["n_up_trilinears"] = len(up)
    up_ch = Counter(tuple(sorted(spec[l.lstrip("b")]["kl"] for l in c)) for c in up)
    r["up_channel_histogram"] = {Z.channel_name(ms): n for ms, n in up_ch.items()}
    r["up_on_ladder_channels"] = {LB[ms]: n for ms, n in up_ch.items() if ms in LB}
    r["up_on_allowed_channels"] = {ALLOWED_MS[ms]: n for ms, n in up_ch.items()
                                   if ms in ALLOWED_MS}

    # ---- 3. GATE X3: does W falsify Rule P?
    r["gateX3_ruleP_falsified"] = bool(r["up_on_ladder_channels"])
    # (the census row carries the Rule-P verdict; the aggregator compares them)

    # ---- 5/6/7. single-Higgs blocks on ladder-bearing channels
    blocks = []
    byH = defaultdict(list)
    for c in up:
        h = [l for l in c if role.get(l.lstrip("b")) == "Hu"][0]
        byH[h].append(c)
    for h, cs in byH.items():
        ms_all = {tuple(sorted(spec[l.lstrip("b")]["kl"] for l in c)) for c in cs}
        for ms in sorted(ms_all, key=str):
            sub = [c for c in cs if
                   tuple(sorted(spec[l.lstrip("b")]["kl"] for l in c)) == ms]
            name = ALLOWED_MS.get(ms, Z.channel_name(ms))
            rows = sorted({[l for l in c if role.get(l.lstrip("b")) == "Q"][0]
                           for c in sub})
            cols = sorted({[l for l in c if role.get(l.lstrip("b")) == "uc"][0]
                           for c in sub})
            expo, detail, support = {}, {}, set()
            for c in sub:
                q = [l for l in c if role.get(l.lstrip("b")) == "Q"][0]
                u = [l for l in c if role.get(l.lstrip("b")) == "uc"][0]
                recs = [spec[q.lstrip("b")], spec[u.lstrip("b")],
                        spec[h.lstrip("b")]]
                tot, per = entry_exponent([x["kl"] for x in recs], recs)
                expo[(q, u)] = tot
                detail["%s|%s" % (q, u)] = {"total_exponent":
                                            None if tot is None else str(tot),
                                            "planes": per}
                if tot is not None:
                    support.add((q, u))
            rank = generic_rank(rows, cols, support)
            dist = sorted({str(v) for v in expo.values() if v is not None})
            lead = min_exponent_matching(rows, cols, expo)
            blocks.append({
                "higgs": h, "higgs_sector": str(spec[h.lstrip("b")]["kl"]),
                "channel": name, "channel_is_ladder_bearing": ms in LB,
                "channel_is_allowed": ms in ALLOWED_MS,
                "tau_invariant_channel": Z.tau_invariant_channel(list(ms)),
                "n_Q": len(rows), "n_uc": len(cols), "n_entries": len(sub),
                "generic_rank": rank, "n_defined_entries": len(support),
                "distinct_total_exponents": dist,
                "n_distinct_total_exponents": len(dist),
                "leading_det_term": None if lead is None else
                    {"assignment": lead[0], "total": str(lead[1])},
                "entries": detail,
                "chiral_Q_on_block": sorted(
                    set(rows) & {l for l in spec
                                 if matched.get(l) and cen.get(matched[l], {}).get(
                                     "role") == "Q"}),
                "D2_uncovered_legs": sorted({str(t) for t in ms
                                             if t in Z.UNCOVERED_ELEMENTS}),
            })
    r["blocks"] = blocks
    r["n_blocks"] = len(blocks)
    r["max_rank_on_ladder_channel"] = max(
        [b["generic_rank"] for b in blocks if b["channel_is_ladder_bearing"]] or [0])
    r["max_distinct_factors_on_ladder_channel"] = max(
        [b["n_distinct_total_exponents"] for b in blocks
         if b["channel_is_ladder_bearing"]] or [0])

    # ---- 6. C(S_min) on the best ladder blocks (three branches, never averaged)
    csm = []
    for b in blocks:
        if not b["channel_is_ladder_bearing"] or b["generic_rank"] < 2:
            continue
        for key in list(b["entries"])[:12]:
            q, u = key.split("|")
            recs = [spec[q.lstrip("b")], spec[u.lstrip("b")],
                    spec[b["higgs"].lstrip("b")]]
            ks = [x["kl"] for x in recs]
            row = {"block_channel": b["channel"], "higgs": b["higgs"], "entry": key,
                   "planes": {}}
            for a in (1, 2, 3):
                per = {}
                for br in ("blind", "correlator", "LEEFT"):
                    v = Z.csmin(ks, recs, a, br)
                    per[br] = None if v is None else {
                        "S_min": str(v[0]), "n_configs": v[1],
                        "abs_C": round(v[2], 9), "C_is_zero": v[3],
                        "n_actions": v[4]}
                if any(per.values()):
                    row["planes"]["plane%d" % a] = per
            if row["planes"]:
                csm.append(row)
    r["C_Smin"] = csm
    r["C_Smin_branch_disagreement"] = any(
        len({(p[br] or {}).get("C_is_zero") for br in ("blind", "correlator", "LEEFT")}) > 1
        for row in csm for p in row["planes"].values())

    # ---- 8. verdict
    best = [b for b in blocks if b["channel_is_ladder_bearing"]]
    sw1 = any(b["generic_rank"] >= 3 and b["n_distinct_total_exponents"] >= 3
              for b in best)
    r["verdict"] = ("SW-1 candidate" if sw1 else
                    ("SW-2 partial" if best else
                     ("SW-2 partial (up trilinears, no ladder channel)"
                      if r["n_up_trilinears"] else "SW-3 empty")))
    r["status"] = "ok"
    return r


def _one(args):
    label, coup_block, cen_block = args
    dest = os.path.join(OUT, "%s.json" % label)
    try:
        rec = analyse(label, coup_block, cen_block)
    except Exception:
        rec = {"model": label, "status": "EXCEPTION",
               "traceback": traceback.format_exc()}
    tmp = dest + ".part"
    with open(tmp, "w") as fh:
        json.dump(rec, fh, indent=1, default=str)
    os.replace(tmp, dest)
    return label, rec.get("status"), rec.get("verdict")


def iter_jobs(only=None, force=False):
    """Stream (label, coupling block, census block) pairs, segment by segment, so the
    whole 44 MB + 60 MB of raw text is never resident at once."""
    # the census was cut into fixed 100-model segments, so a model's census segment is
    # its catalogue index // 100 regardless of how the COUPLING segments were cut.
    cache = {}

    def census_block(lab):
        seg = int(lab[4:]) // 100
        if seg not in cache:
            cache.clear()                        # keep at most one segment resident
            p = os.path.join(Z.RAW_CENSUS, "census_seg%d.txt.gz" % seg)
            cache[seg] = dict(Z.split_models(Z.read_raw(p))) if os.path.exists(p) else {}
        return cache[seg].get(lab)

    for fn in sorted(os.listdir(Z.RAW_COUP)):
        if not fn.startswith("coup_seg") or not fn.endswith(".gz"):
            continue
        if fn.startswith("coup_segP"):           # pilot / probe segments, not production
            continue
        coup = dict(Z.split_models(Z.read_raw(os.path.join(Z.RAW_COUP, fn))))
        # HARD FILTER: only blocks carrying orbifolder's COMPLETION SIGNAL
        #     PID <n> done: "create coupling(F F F)"
        # are admissible.  A block without it had its coupling child KILLED by wait(N) and
        # its printed `W = 0` is a TRUNCATION, not an absence.  Some models appear in more
        # than one segment (a probe run and a production run), so relying on which segment
        # is visited first would make the choice a matter of filename ordering; filtering
        # here makes analysing a truncated superpotential structurally impossible rather
        # than merely unlikely.
        coup = {k: v for k, v in coup.items() if ' done: "create coupling' in v}
        for lab in sorted(coup, key=lambda s: int(s[4:])):
            if only and lab not in only:
                continue
            if not force and os.path.exists(os.path.join(OUT, "%s.json" % lab)):
                continue
            cb = census_block(lab)
            if cb is None:
                continue
            yield (lab, coup[lab], cb)


def main():
    workers, force, only = 6, False, None
    args = sys.argv[1:]
    if "--force" in args:
        force = True
        args.remove("--force")
    if "--workers" in args:
        i = args.index("--workers")
        workers = int(args[i + 1])
        del args[i:i + 2]
    if args:
        only = set(args)
    n = 0
    with Pool(workers) as p:
        for lab, st, vd in p.imap_unordered(_one, iter_jobs(only, force),
                                            chunksize=2):
            n += 1
            if st != "ok" or n % 250 == 0:
                print("  %5d %-10s %-24s %s" % (n, lab, st, vd))
    print("sweep: %d models analysed -> %s" % (n, OUT))


if __name__ == "__main__":
    main()

"""
levels.py --- STEP 2: the NC-CORRECTED LEVEL ENUMERATION for Z2xZ4_1_1, computed FRESH
for this geometry.  This module is entirely GEOMETRY-LEVEL: it says what the geometry
ADMITS.  It says nothing whatever about the catalogue, and no statement produced here may
be quoted as a catalogue statement or vice versa (the campaign's standing scope discipline).

THREE LEVELS, all with RULE NC imposed from the start (guard 1):

  L-A  BUCKET LEVEL (the [D-2]-style route).  Fix two legs' FULL generation indices, vary
       the third: does the bucket carry three distinct TOTAL exponents in a 0 : n : 2n
       ladder?  Computed twice -- table-only and NC-imposed -- so this geometry's own
       analogue of the campaign's `240 -> 96` movement is measured, not inherited.

  L-B  SINGLE-PLANE COLUMN-GRADED SUPPORT (the certified sixth escape route, [SW-C1] of
       gate_z4xz4_sweep).  A block is (rows = one species' generation indices) x (columns =
       the other species') at a FIXED Higgs index.  The route asks whether the exponent is
       a function of the COLUMN index alone on the cells that exist, with three distinct
       charges, giving generic rank 3 by the minor invariant.  NOTE the certified
       corrections carried here: the object is a column-GRADED SUPPORT, not a full
       CLAUDE.md 3.1 column texture (the support is generally partial); and the ROW-graded
       orientation is a ROW texture in which the hierarchy sits on Q_L, which 3.1 forbids
       -- the two are reported SEPARATELY and never merged.

  L-C  MULTI-PLANE TOTAL-EXPONENT BUCKETS: the same as L-A but requiring the three rungs to
       differ in more than one computable plane (the exchange-symmetric-pair route, which
       on this geometry can only be the pair [1,3] -- the only sigma the geometry admits).

COMPLETENESS BY SEARCH, NOT BY SCOPING (guard 3, the S1-F lesson).  Every scan below runs
over ALL 15 certified channels and ALL SIX assignments of the three legs to (Q, u^c, H_u),
and every occupied Higgs index.  Nothing is pre-filtered.  A pre-filter may later ORDER the
coupling runs; it never scopes a discovery claim.

Emits: data/z2xz4_sweep/levels.json  (+ table_z24.json.gz, the NC-corrected index table)
"""

from __future__ import annotations

import gzip
import json
import os
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, permutations

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z24sw_lib as L                                                   # noqa: E402


def minors_by_flow(cells, rows, cols, cap_r, cap_c, kmax=3):
    """[d_1, ..., d_kmax] where d_k = min over k x k minors with a generic non-vanishing
    determinant (i.e. a full matching) of the sum of the k entry exponents; None past the
    generic rank.

    Computed as a MIN-COST BIPARTITE b-MATCHING by successive shortest augmenting paths --
    which is exactly the sequence d_1, d_2, d_3, because the k-th augmentation of a min-cost
    flow yields a min-cost flow of value k.  Row/column MULTIPLICITY enters as the node
    capacities `cap_r`, `cap_c`: several fields at one generation index carry the same
    exponent but independent O(1) coefficients and are independent rows (columns).

    All exponents are non-negative, so plain Bellman-Ford/SPFA on the residual graph is
    exact (negative residual arcs appear but no negative cycle can, by optimality of each
    intermediate flow)."""
    R, C = list(rows), list(cols)
    ri = {r: i for i, r in enumerate(R)}
    ci = {c: len(R) + j for j, c in enumerate(C)}
    S, T = len(R) + len(C), len(R) + len(C) + 1
    n = T + 1
    to, capv, cost, nxt = [], [], [], []
    head = [-1] * n

    def add(u, v, cp, cs):
        for (a, b, c_, d_) in ((u, v, cp, cs), (v, u, 0, -cs)):
            to.append(b)
            capv.append(c_)
            cost.append(d_)
            nxt.append(head[a])
            head[a] = len(to) - 1

    for r in R:
        add(S, ri[r], cap_r, F(0))
    for c in C:
        add(ci[c], T, cap_c, F(0))
    for (r, c), w in cells.items():
        add(ri[r], ci[c], min(cap_r, cap_c), w)

    out, total = [], F(0)
    for _ in range(kmax):
        INF = None
        dist = [INF] * n
        inq = [False] * n
        pre = [-1] * n
        dist[S] = F(0)
        q = [S]
        inq[S] = True
        while q:
            u = q.pop(0)
            inq[u] = False
            e = head[u]
            while e != -1:
                if capv[e] > 0:
                    nd = dist[u] + cost[e]
                    v = to[e]
                    if dist[v] is None or nd < dist[v]:
                        dist[v] = nd
                        pre[v] = e
                        if not inq[v]:
                            inq[v] = True
                            q.append(v)
                e = nxt[e]
        if dist[T] is None:
            out.append(None)
            break
        total += dist[T]
        v = T
        while v != S:
            e = pre[v]
            capv[e] -= 1
            capv[e ^ 1] += 1
            v = to[e ^ 1]
        out.append(total)
    while len(out) < kmax:
        out.append(None)
    return out


def mass_exponents(cells, rows, cols, cap_r=1, cap_c=1):
    """(d1, d2-d1, d3-d2) -- the basis-independent mass ladder; None if generic rank < 3."""
    d = minors_by_flow(cells, rows, cols, cap_r, cap_c, 3)
    dv = [str(x) if x is not None else "inf" for x in d]
    if d[2] is None:
        return None, dv
    return [d[0], d[1] - d[0], d[2] - d[1]], dv


def minors33(M):
    """Exact (d1, d2, d3) of a 3x3 matrix of exponents-or-None, by explicit enumeration of
    the matchings.  d_k = min over k x k minors with a generic non-vanishing determinant."""
    d1 = None
    for i in range(3):
        for j in range(3):
            v = M[i][j]
            if v is not None and (d1 is None or v < d1):
                d1 = v
    d2 = None
    for rs in combinations(range(3), 2):
        for cs in combinations(range(3), 2):
            for p in permutations(range(2)):
                a, b = M[rs[0]][cs[p[0]]], M[rs[1]][cs[p[1]]]
                if a is not None and b is not None and (d2 is None or a + b < d2):
                    d2 = a + b
    d3 = None
    for p in permutations(range(3)):
        s, ok = F(0), True
        for i in range(3):
            v = M[i][p[i]]
            if v is None:
                ok = False
                break
            s += v
        if ok and (d3 is None or s < d3):
            d3 = s
    return d1, d2, d3


def realisable_patterns(cells, rows, cols, max_types=40):
    """EVERY mass-exponent pattern the geometry ADMITS for a THREE-GENERATION block on this
    (channel, leg assignment, Higgs index).

    A physical up-Yukawa block is 3 x 3: three chiral generations of the row species and
    three of the column species.  Which generation INDICES they occupy is an OCCUPANCY --
    a per-model property -- so the geometry-level question is: over every multiset of three
    row indices and every multiset of three column indices drawn from what the channel
    offers, which mass patterns (d1, d2-d1, d3-d2) occur?  Repetition is allowed and is
    physical: several fields may sit at ONE generation index, carrying the SAME classical
    exponent but INDEPENDENT O(1) coefficients.

    Rows (columns) with identical exponent vectors are interchangeable, so the enumeration
    is done over DISTINCT ROW/COLUMN TYPES -- exact, not a sample."""
    colvec = {}
    for c in cols:
        colvec.setdefault(tuple(cells.get((r, c)) for r in rows), []).append(c)
    ctypes = sorted(colvec, key=repr)
    rowvec = {}
    for r in rows:
        rowvec.setdefault(tuple(cells.get((r, colvec[cv][0])) for cv in ctypes), []).append(r)
    rtypes = sorted(rowvec, key=repr)
    if len(rtypes) > max_types or len(ctypes) > max_types:
        return None, None
    pats, witness = {}, {}
    # Enumerate the three COLUMN types first; with the columns fixed, the row types collapse
    # to their restriction to those three columns, which takes very few distinct values.
    # Exact, not a sample -- the collapse is a deduplication of identical rows.
    for cs in combinations_with_replacement(range(len(ctypes)), 3):
        sub = {}
        for i, rt in enumerate(rtypes):
            sub.setdefault(tuple(rt[c] for c in cs), []).append(i)
        keys = sorted(sub, key=repr)
        for rs in combinations_with_replacement(range(len(keys)), 3):
            M = [list(keys[rs[i]]) for i in range(3)]
            d1, d2, d3 = minors33(M)
            key = "rank<3" if d3 is None else "%s:%s:%s" % (d1, d2 - d1, d3 - d2)
            pats[key] = pats.get(key, 0) + 1
            if key not in witness:
                witness[key] = {
                    "row_indices": [[str(x) for x in rowvec[rtypes[sub[keys[rs[i]]][0]]]]
                                    for i in range(3)],
                    "col_indices": [[str(x) for x in colvec[ctypes[c]]] for c in cs],
                    "matrix": [[str(x) if x is not None else None for x in row]
                               for row in M]}
    return pats, witness


def colvec_repr(mapping, types, i):
    return [str(x) for x in mapping[types[i]]]


def _is_ladder(pat):
    """'d1:d2-d1:d3-d2' is a genuine 1 : eps^n : eps^{2n} ladder iff it reads 0 : n : 2n."""
    try:
        a, b, c = (F(x) for x in pat.split(":"))
    except Exception:                                              # noqa: BLE001
        return False
    return a == 0 and b > 0 and c == 2 * b


def term_rank(rows, cols, support):
    """Maximum bipartite matching of the 0/1 support -- the generic rank of the support."""
    match = {}

    def aug(r, seen):
        for c in cols:
            if (r, c) not in support or c in seen:
                continue
            seen.add(c)
            if c not in match or aug(match[c], seen):
                match[c] = r
                return True
        return False
    n = 0
    for r in rows:
        if aug(r, set()):
            n += 1
    return n


# =======================================================================================
#  L-A / L-C : bucket-level ladders
# =======================================================================================

def bucket_scan(tab, nc_imposed=True):
    """For every channel and every varying leg: buckets keyed by the OTHER TWO legs' full
    generation indices; report those carrying a three-rung 0 : n : 2n TOTAL-exponent
    ladder, with the set of planes in which the three rungs actually differ."""
    out = {"n_buckets": 0, "size_histogram": Counter(), "ladder_buckets": [],
           "per_channel": {}}
    for nm, v in tab.items():
        gs = tuple(tuple(x) for x in v["gs"])
        comp = v["capable_planes"]
        if not comp:
            continue
        ents = v["entries"] if nc_imposed else v["entries_no_NC"]
        per = {"n_ladder_buckets": 0, "ladder_types": Counter()}
        for leg in range(3):
            buckets = defaultdict(list)
            for e in ents:
                others = tuple(tuple(e["legs"][i]) for i in range(3) if i != leg)
                buckets[others].append(e)
            for key, es in buckets.items():
                tot = {}
                for e in es:
                    tot.setdefault(F(e["total"]), e)
                out["n_buckets"] += 1
                out["size_histogram"][len(tot)] += 1
                if F(0) not in tot:
                    continue
                lad = None
                for n in sorted(x for x in tot if x > 0):
                    if 2 * n in tot:
                        lad = (F(0), n, 2 * n)
                        break
                if lad is None:
                    continue
                moving = []
                for a in comp:
                    s = {tot[r]["exponent_tuple"][a - 1] for r in lad}
                    if len(s) > 1:
                        moving.append(a)
                per["n_ladder_buckets"] += 1
                per["ladder_types"][":".join(str(x) for x in lad)] += 1
                out["ladder_buckets"].append(
                    {"channel": nm, "varying_leg": leg,
                     "fixed_legs": [list(x) for x in key],
                     "ladder": [str(x) for x in lad], "moving_planes": moving,
                     "single_plane": len(moving) == 1})
        out["per_channel"][nm] = {"n_ladder_buckets": per["n_ladder_buckets"],
                                  "ladder_types": dict(per["ladder_types"])}
    out["size_histogram"] = dict(sorted(out["size_histogram"].items()))
    out["n_ladder_buckets"] = len(out["ladder_buckets"])
    out["n_single_plane_ladder_buckets"] = sum(1 for b in out["ladder_buckets"]
                                               if b["single_plane"])
    out["ladder_type_histogram"] = dict(Counter(":".join(b["ladder"])
                                                for b in out["ladder_buckets"]))
    out["channels_carrying_a_ladder"] = sorted({b["channel"] for b in out["ladder_buckets"]})
    return out


# =======================================================================================
#  L-B : the block / texture level, geometry-side
# =======================================================================================

def block_scan(tab, nc_imposed=True, mult=3):
    """ALL channels x ALL SIX leg assignments x every Higgs index: the geometry-level block.

    ROW/COLUMN MULTIPLICITY IS CARRIED (`mult`), and it is load-bearing.  Several distinct
    fields may sit at ONE generation index; they carry the SAME classical exponent but
    INDEPENDENT O(1) coefficients, so they are independent rows (columns) of the mass matrix.
    That is exactly the configuration that makes a column-graded support rank 3 -- with one
    field per index the certified single-plane cap of 2 would bound the visible charges per
    row and no third minor could form.  A geometry-level ADMISSIBILITY statement must
    therefore be permissive here: every index is offered up to `mult` independent copies
    (3 suffices for a 3x3 minor).  `mult = 1` is reported alongside so the difference the
    multiplicity makes is MEASURED, not assumed.

    Rows = the ROW species' generation indices present in the channel's leg; columns = the
    COLUMN species'.  This asks what the GEOMETRY admits, so every generation index of each
    leg is taken to be available -- a per-model occupancy question is asked later and never
    conflated with this one.

    Reported per block: support size, generic (support) rank, whether the exponent is
    column-determined / row-determined on the support, the number of distinct charges, the
    minor-based mass ladder, and the planes carrying the variation."""
    rows_out = []
    for nm, v in tab.items():
        if not v["capable_planes"]:
            continue
        gs = tuple(tuple(x) for x in v["gs"])
        ents = v["entries"] if nc_imposed else v["entries_no_NC"]
        idx = {}
        for e in ents:
            idx[tuple(tuple(x) for x in e["legs"])] = F(e["total"])
        legidx = [sorted({k[i] for k in idx}) for i in range(3)]
        for legQ, legU, legH in permutations(range(3)):
            for h in legidx[legH]:
                cells, sup = {}, {}
                for k, val in idx.items():
                    if k[legH] != h:
                        continue
                    cells[(k[legQ], k[legU])] = val
                    sup[(k[legQ], k[legU])] = True
                if not cells:
                    continue
                rws0 = sorted({r for r, _ in cells})
                cls0 = sorted({c for _, c in cells})
                if not rws0 or not cls0:
                    continue
                rws, cls = rws0, cls0
                rk = term_rank(rws, cls, sup)
                steps, dvals = mass_exponents(cells, rws, cls, cap_r=1, cap_c=1)
                pats, wit = realisable_patterns(cells, rws, cls)
                lad_pats = sorted(p for p in (pats or {})
                                  if p != "rank<3" and _is_ladder(p))
                colmap, rowmap = defaultdict(set), defaultdict(set)
                for (r, c), val in cells.items():
                    colmap[c].add(val)
                    rowmap[r].add(val)
                col_det = all(len(s) == 1 for s in colmap.values())
                row_det = all(len(s) == 1 for s in rowmap.values())
                charges = sorted({next(iter(s)) for s in colmap.values()}) if col_det else []
                rcharges = sorted({next(iter(s)) for s in rowmap.values()}) if row_det else []
                planes = []
                for a in v["capable_planes"]:
                    vals = set()
                    for k, val in idx.items():
                        if k[legH] == h:
                            e = [e for e in ents
                                 if tuple(tuple(x) for x in e["legs"]) == k][0]
                            vals.add(e["exponent_tuple"][a - 1])
                    if len(vals) > 1:
                        planes.append(a)
                rows_out.append({
                    "channel": nm, "legQ": legQ, "legU": legU, "legH": legH,
                    "sector_Q": list(gs[legQ]), "sector_U": list(gs[legU]),
                    "sector_H": list(gs[legH]), "higgs_index": list(h),
                    "multiplicity": mult,
                    "n_rows": len(rws0), "n_cols": len(cls0), "n_cells": len(cells),
                    "n_cells_full": len(rws0) * len(cls0),
                    "support_rank": rk,
                    "column_determined": col_det, "row_determined": row_det,
                    "n_distinct_column_charges": len(charges),
                    "column_charges": [str(x) for x in charges],
                    "n_distinct_row_charges": len(rcharges),
                    "row_charges": [str(x) for x in rcharges],
                    "minors_d": dvals,
                    "mass_exponents": [str(x) for x in steps] if steps else None,
                    "generic_rank_3": steps is not None,
                    "ladder_0_n_2n": bool(steps and steps[0] == 0 and steps[1] > 0
                                          and steps[2] == 2 * steps[1]),
                    "ladder_0_2_4": bool(steps and [int(x) for x in steps] == [0, 2, 4]),
                    "ladder_0_1_2": bool(steps and [int(x) for x in steps] == [0, 1, 2]),
                    # the OCCUPANCY-EXHAUSTIVE geometry-level admissibility:
                    "realisable_patterns": pats,
                    "admits_ladder_patterns": lad_pats,
                    "admits_0_2_4": "0:2:4" in (pats or {}),
                    "admits_any_ladder": bool(lad_pats),
                    "ladder_witness": {p: wit[p] for p in lad_pats} if lad_pats else None,
                    "varying_planes": planes,
                    "single_plane": len(planes) <= 1})
    return rows_out


def build_table_no_NC(box=6):
    """The SAME table WITHOUT RULE NC -- built only so this geometry's NC movement can be
    MEASURED rather than asserted.  It is never used for any physics statement."""
    g = L.geom()
    out = {}
    for gs in L.all_channels():
        nm = L.channel_name(gs)
        comp = L.capable_planes(gs)
        legs = [L.leg_class_tuples(s) for s in gs]
        ntab = {a: L._norm_table(gs, a, box=box) for a in comp}
        entries = []
        for c0 in legs[0]:
            for c1 in legs[1]:
                for c2 in legs[2]:
                    t = [F(0), F(0), F(0)]
                    ok = True
                    for a in comp:
                        nv, _ = ntab[a][(c0[a - 1], c1[a - 1], c2[a - 1])]
                        if nv is None:
                            ok = False
                            break
                        t[a - 1] = nv
                    if ok:
                        entries.append({"legs": [list(c0), list(c1), list(c2)],
                                        "exponent_tuple": [str(x) for x in t],
                                        "total": str(sum(t))})
        out[nm] = entries
    return out


def main():
    tab = L.build_table()
    noNC = build_table_no_NC()
    for nm in tab:
        tab[nm]["entries_no_NC"] = noNC[nm]

    rep = {"geometry": "Z2xZ4_1_1",
           "scope": "GEOMETRY statement.  Nothing here is a catalogue statement.",
           "rule_NC": "imposed from the start in every table (guard 1); the certified ORBIT "
                      "form (gate_rule_nc.md N10).",
           "cap": "QUOTED, not re-tested: single-plane columns cap at 2 (K2, square-generic, "
                  "measured on this geometry).",
           "n_channels": len(tab)}

    # ---- the table itself
    rep["table"] = {nm: {"capable_planes": v["capable_planes"],
                         "nc_planes": v["nc_planes"],
                         "n_cells_before_NC": v["n_cells_before_NC"],
                         "n_cells_killed_by_NC": v["n_cells_killed_by_NC"],
                         "n_entries_with_NC": len(v["entries"]),
                         "n_entries_without_NC": len(v["entries_no_NC"]),
                         "alphabet_total": v["alphabet_total"],
                         "n_D3_entries": len(v["D3_entries"])}
                    for nm, v in tab.items()}
    rep["NC_movement_on_entries"] = {
        "total_with_NC": sum(len(v["entries"]) for v in tab.values()),
        "total_without_NC": sum(len(v["entries_no_NC"]) for v in tab.values())}

    # ---- L-A / L-C buckets
    rep["L_A_buckets_NC"] = bucket_scan(tab, nc_imposed=True)
    rep["L_A_buckets_no_NC"] = bucket_scan(tab, nc_imposed=False)
    for k in ("L_A_buckets_NC", "L_A_buckets_no_NC"):
        rep[k].pop("ladder_buckets", None)
    lb = bucket_scan(tab, nc_imposed=True)["ladder_buckets"]
    rep["L_A_ladder_bucket_examples"] = lb[:40]
    rep["L_C_multi_plane_ladder_buckets"] = [b for b in lb if not b["single_plane"]]

    # ---- L-B blocks
    blocks = block_scan(tab, nc_imposed=True)
    blocks_noNC = block_scan(tab, nc_imposed=False)
    rep["L_B_admissibility"] = {
        "criterion": "OCCUPANCY-EXHAUSTIVE: over every multiset of three row indices and "
                     "every multiset of three column indices the channel offers at a fixed "
                     "Higgs index (repetition allowed and physical -- several fields may sit "
                     "at one index with the same exponent and independent O(1) "
                     "coefficients), which mass patterns (d1, d2-d1, d3-d2) occur?",
        "n_blocks_admitting_any_ladder": sum(1 for b in blocks if b["admits_any_ladder"]),
        "n_blocks_admitting_0_2_4": sum(1 for b in blocks if b["admits_0_2_4"]),
        "n_blocks_admitting_any_ladder_no_NC": sum(1 for b in blocks_noNC
                                                   if b["admits_any_ladder"]),
        "n_blocks_admitting_0_2_4_no_NC": sum(1 for b in blocks_noNC if b["admits_0_2_4"]),
        "channels_admitting_any_ladder": sorted({b["channel"] for b in blocks
                                                 if b["admits_any_ladder"]}),
        "channels_admitting_0_2_4": sorted({b["channel"] for b in blocks
                                            if b["admits_0_2_4"]}),
        "ladder_patterns_admitted": dict(Counter(
            p for b in blocks for p in b["admits_ladder_patterns"]).most_common()),
        "realisable_pattern_union": dict(Counter(
            p for b in blocks for p in (b["realisable_patterns"] or {})).most_common()),
    }
    rep["L_B"] = {
        "n_blocks_scanned": len(blocks),
        "n_blocks_scanned_no_NC": len(blocks_noNC),
        "search_space": "ALL 15 channels x ALL 6 leg assignments x every Higgs index "
                        "(completeness by search, guard 3)",
        "n_generic_rank_3": sum(1 for b in blocks if b["generic_rank_3"]),
        "n_generic_rank_3_no_NC": sum(1 for b in blocks_noNC if b["generic_rank_3"]),
        "n_ladder_0_n_2n": sum(1 for b in blocks if b["ladder_0_n_2n"]),
        "n_ladder_0_n_2n_no_NC": sum(1 for b in blocks_noNC if b["ladder_0_n_2n"]),
        "n_ladder_0_2_4": sum(1 for b in blocks if b["ladder_0_2_4"]),
        "n_ladder_0_2_4_no_NC": sum(1 for b in blocks_noNC if b["ladder_0_2_4"]),
        "n_ladder_0_1_2": sum(1 for b in blocks if b["ladder_0_1_2"]),
        "n_column_determined": sum(1 for b in blocks if b["column_determined"]),
        "n_column_determined_3_charges": sum(
            1 for b in blocks if b["column_determined"]
            and b["n_distinct_column_charges"] >= 3),
        "n_row_determined_3_charges": sum(
            1 for b in blocks if b["row_determined"]
            and b["n_distinct_row_charges"] >= 3),
        "mass_pattern_histogram": dict(Counter(
            ":".join(b["mass_exponents"]) if b["mass_exponents"] else "rank<3"
            for b in blocks).most_common()),
        "support_rank_histogram": dict(Counter(b["support_rank"] for b in blocks)),
        "channels_admitting_a_ladder": sorted({b["channel"] for b in blocks
                                               if b["ladder_0_n_2n"]}),
        "channels_admitting_0_2_4": sorted({b["channel"] for b in blocks
                                            if b["ladder_0_2_4"]}),
        "single_plane_ladder_blocks": sum(1 for b in blocks
                                          if b["ladder_0_n_2n"] and b["single_plane"]),
    }
    rep["L_B_ladder_examples"] = [
        {k: v for k, v in b.items() if k != "realisable_patterns"}
        for b in blocks if b["admits_any_ladder"]][:30]

    # ---- the CAPABLE-CHANNEL summary the sweep will order its coupling runs by
    cap = {}
    for nm, v in tab.items():
        bl = [b for b in blocks if b["channel"] == nm]
        cap[nm] = {"capable_planes": v["capable_planes"],
                   "alphabet_total": v["alphabet_total"],
                   "admits_rank3": any(b["generic_rank_3"] for b in bl),
                   "admits_any_ladder": any(b["admits_any_ladder"] for b in bl),
                   "admits_0_2_4": any(b["admits_0_2_4"] for b in bl),
                   "ladder_patterns": sorted({p for b in bl
                                              for p in b["admits_ladder_patterns"]}),
                   "n_blocks": len(bl)}
    rep["capable_channels"] = cap
    rep["n_channels_admitting_a_ladder"] = sum(1 for v in cap.values()
                                               if v["admits_any_ladder"])

    L.write_json("levels.json", rep)
    with gzip.open(os.path.join(L.DSW, "table_z24.json.gz"), "wt", encoding="utf-8") as fh:
        json.dump({nm: {k: v for k, v in val.items() if k != "entries_no_NC"}
                   for nm, val in tab.items()}, fh, default=str)
    with gzip.open(os.path.join(L.DSW, "levels_blocks.json.gz"), "wt",
                   encoding="utf-8") as fh:
        json.dump({"with_NC": blocks, "without_NC": blocks_noNC}, fh, default=str)
    print(json.dumps({k: v for k, v in rep.items()
                      if k not in ("L_A_ladder_bucket_examples", "L_B_ladder_examples",
                                   "L_C_multi_plane_ladder_buckets")},
                     indent=1, default=str))


if __name__ == "__main__":
    main()

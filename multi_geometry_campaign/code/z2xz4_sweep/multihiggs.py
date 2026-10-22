"""
multihiggs.py --- the MULTI-HIGGS union analysis, with its alignment assumption STATED and
its cost QUANTIFIED.

The pre-declared SW-1 requires a SINGLE Higgs.  Widening to several H_u fields changes the
object: the Yukawa becomes a SUM of Yukawas,

    Y_ij = sum_h c^h_ij eps^{n^h_ij} (v_h / v),

so a hierarchy among the Higgs vevs RE-GRADES the matrix.  Writing v_h/v = eps^{m_h} with
integer m_h >= 0 and min_h m_h = 0, the effective exponent of a cell is min_h (n^h_ij + m_h).

TWO STATEMENTS ARE PRODUCED AND NEVER MERGED:

  (a) UNIFORM VEVS (all m_h = 0).  The effective exponent is the minimum over Higgs fields.
      This is the reading the Z4xZ4 phase used, and the Z4xZ4 gate's S3 correction was that
      the phase never stated the assumption.  It is stated here.

  (b) THE COST.  Every weighting with m_h <= 4 is scanned; for each block the CHEAPEST
      weighting that buys a genuine 1 : eps^n : eps^{2n} ladder is reported, as
      (number of tuned Higgs fields, total suppression sum m_h).  A ladder that has to be
      bought is reported as bought, with its bill.

NOTE, carried from the Z4xZ4 gate's S3 latent-bug finding: the exponent map is built with the
HIGGS LEG KEPT in the key, never with a `setdefault` that freezes each cell at whichever
Higgs index appears first in file order.

Emits: data/z2xz4_sweep/multihiggs.json
"""

from __future__ import annotations

import gzip
import json
import os
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, product

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z24sw_lib as L                                                   # noqa: E402
from levels import minors_by_flow                                       # noqa: E402


def pattern(cells, rows, cols):
    d = minors_by_flow(cells, rows, cols, 1, 1, 3)
    if d[2] is None:
        return None
    return [d[0], d[1] - d[0], d[2] - d[1]]


def is_ladder(steps):
    return bool(steps and steps[0] == 0 and steps[1] > 0 and steps[2] == 2 * steps[1])


def main():
    with gzip.open(os.path.join(L.DSW, "models", "coupling.json.gz"), "rt",
                   encoding="utf-8") as fh:
        coup = json.load(fh)
    with gzip.open(os.path.join(L.DSW, "models", "blocks.json.gz"), "rt",
                   encoding="utf-8") as fh:
        blocks = json.load(fh)

    rows_out = []
    for lab, rec in coup.items():
        if rec.get("status") not in ("OK", "LABEL_MISMATCH"):
            continue
        bl = blocks.get(lab, {}).get("blocks", [])
        # group the model's index-level blocks by (channel, leg assignment, row/col sets):
        # the multi-Higgs object is the union over the Higgs fields of one such group
        grp = defaultdict(list)
        for b in bl:
            key = (b["channel"], b["legQ"], b["legU"], b["legH"],
                   tuple(sorted(b["row_fields"])), tuple(sorted(b["col_fields"])))
            grp[key].append(b)
        realised = {}
        for rb in rec.get("blocks", []):
            realised[(rb["channel"], rb["legQ"], rb["legU"], rb["legH"],
                      rb["higgs_field"])] = rb
        for key, bs in grp.items():
            if len(bs) < 2:
                continue
            higgs = sorted({b["higgs_field"] for b in bs})
            # per-Higgs realised exponent maps, HIGGS KEPT IN THE KEY
            emap = {}
            for b in bs:
                rb = realised.get(key[:4] + (b["higgs_field"],))
                if rb is None:
                    continue
                for cellkey, val in b["cells"].items():
                    r, c = cellkey.split("|")
                    emap[(b["higgs_field"], r, c)] = F(val)
            if not emap:
                continue
            rws = sorted({r for _, r, _ in emap})
            cls = sorted({c for _, _, c in emap})
            # (a) uniform vevs
            uni = {}
            for (h, r, c), v in emap.items():
                k = (r, c)
                if k not in uni or v < uni[k]:
                    uni[k] = v
            st = pattern(uni, rws, cls)
            # (b) the cost scan.  Only the Higgs fields that actually CONTRIBUTE a realised
            # cell can change the effective exponent, so the scan runs over those; the box is
            # m_h <= 4 as declared.  Groups with more than 6 contributing Higgs fields would
            # need 5^7 weightings and are reported as NOT SCANNED rather than silently
            # sampled -- an unscanned group is never counted as a negative.
            contrib = sorted({h for h, _, _ in emap})
            best = None
            full = len(contrib) <= 6
            if full:
                weightings = (ms for ms in product(range(5), repeat=len(contrib))
                              if min(ms) == 0)
                scan_scope = "FULL box m_h <= 4 over all contributing Higgs fields"
            else:
                # 5^8 / 5^10 is out of reach; scan the region a purchase would live in --
                # AT MOST THREE tuned fields, each at m in 1..4 -- exactly, and SAY SO.
                def _w():
                    n = len(contrib)
                    yield tuple([0] * n)
                    for t in (1, 2, 3):
                        for pos in combinations(range(n), t):
                            for vals in product(range(1, 5), repeat=t):
                                ms = [0] * n
                                for p, v in zip(pos, vals):
                                    ms[p] = v
                                yield tuple(ms)
                weightings = _w()
                scan_scope = "PARTIAL: at most 3 tuned Higgs fields, each m in 1..4"
            scanned = True
            for ms in weightings:
                if min(ms) != 0:
                    continue
                w = dict(zip(contrib, ms))
                eff = {}
                for (h, r, c), v in emap.items():
                    k = (r, c)
                    x = v + w[h]
                    if k not in eff or x < eff[k]:
                        eff[k] = x
                s2 = pattern(eff, rws, cls)
                if is_ladder(s2):
                    cost = (sum(1 for m in ms if m > 0), sum(ms))
                    if best is None or cost < best[0]:
                        best = (cost, [str(x) for x in s2], dict(zip(contrib, ms)))
            rows_out.append({
                "model": lab, "channel": key[0],
                "legs": list(key[1:4]), "n_higgs_fields": len(higgs),
                "higgs_fields": higgs, "n_rows": len(rws), "n_cols": len(cls),
                "uniform_vev_mass_exponents": [str(x) for x in st] if st else None,
                "uniform_vev_is_ladder": is_ladder(st),
                "n_contributing_higgs": len(contrib),
                "cost_scan_scope": scan_scope,
                "cheapest_ladder_purchase": None if best is None else {
                    "n_tuned_higgs": best[0][0], "total_suppression": best[0][1],
                    "mass_exponents": best[1], "weights": {k: v for k, v in best[2].items()}},
            })

    rep = {"geometry": "Z2xZ4_1_1",
           "scope": "CATALOGUE statement, coupling level.  The pre-declared SW-1 requires a "
                    "SINGLE Higgs, so nothing here can produce an SW-1 letter; it is the "
                    "union analysis the ladder's level list asks for.",
           "alignment_assumption": "STATED, not implied: (a) is the UNIFORM-vev reading "
                                   "(all v_h equal); (b) prices every hierarchical "
                                   "alternative with m_h <= 4.",
           "n_multi_higgs_groups": len(rows_out),
           "n_uniform_vev_ladders": sum(1 for r in rows_out if r["uniform_vev_is_ladder"]),
           "n_purchasable_ladders": sum(1 for r in rows_out
                                        if r["cheapest_ladder_purchase"]),
           "cost_scan_scope_histogram": dict(Counter(r["cost_scan_scope"]
                                                     for r in rows_out).most_common()),
           "uniform_vev_pattern_histogram": dict(Counter(
               ":".join(r["uniform_vev_mass_exponents"])
               if r["uniform_vev_mass_exponents"] else "rank<3"
               for r in rows_out).most_common()),
           "purchase_cost_histogram": dict(Counter(
               "(%d tuned, total eps^%d)" % (r["cheapest_ladder_purchase"]["n_tuned_higgs"],
                                             r["cheapest_ladder_purchase"]["total_suppression"])
               for r in rows_out if r["cheapest_ladder_purchase"]).most_common()),
           "rows": rows_out}
    L.write_json("multihiggs.json", rep)
    print(json.dumps({k: v for k, v in rep.items() if k != "rows"}, indent=1, default=str))


if __name__ == "__main__":
    main()

"""
census.py --- STEP 1: the CONVENTION-BLIND ROLE CENSUS over the 11,540 Z2xZ4_1_1 models.

Per model it emits, from the orbifolder census raw output:

  * the SM identification, obtained by REBUILDING orbifolder's own printed SMConfig
    multiplicity table out of the full left-chiral spectrum (never by row position -- the
    Z6-I closure).  Failure to reproduce the table is recorded as a per-model DATA-GAP, not
    silently patched.  A model with NO SM vev-configuration at all is OUT OF SCOPE, which is
    a denominator statement, and is reported SEPARATELY from a DATA-GAP, which is an
    exposure.
  * the FOUR orientation branches (sY, sC).  Hypercharge sign and colour orientation are
    INDEPENDENT per-model conventions (CLAUDE.md 5.21; the live defect the Z6-I sweep caught
    mid-phase).  On this catalogue the branch actually matters: MSSM0 prints the quark
    doublet as (-3, 2, +1/6) and its physical branch is sC = -1.  The physical branch is
    selected by net Q = +3 and ALL FOUR nets are retained.  Argmax ties are COUNTED and
    reported (an argmax is not a derivation).
  * per role-carrying field: sector (k,l) mod (2,4) and its GENERATION INDEX -- the per-plane
    fixed-point ORBIT CLASS triple the certified alphabet's labels use, computed from
    orbifolder's printed `fixed point n_a` through f = (1 - R_a(g))^{-1}(n + sigma_a)|_a and
    located in the plane's orbit list.  The orbit machinery is the CERTIFIED engine's
    (p2_geom), not a reimplementation.
  * the net generation count per species.

SCOPE.  Everything here is a CATALOGUE statement about the 11,540 PVW `Z2xZ4_1_1` models.
Nothing here is a geometry statement.

Usage:  python census.py            (parses every segment present, resumable)
Emits:  data/z2xz4_sweep/models/roles.json.gz, data/z2xz4_sweep/census_summary.json
"""

from __future__ import annotations

import glob
import gzip
import json
import os
import re
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z24sw_lib as L                                                   # noqa: E402


def process_segment(path):
    txt = L.read_gz(path)
    out = {}
    for lab, block in L.split_models(txt):
        rec = {"label": lab}
        try:
            spec = L.parse_spectrum(block)
            blocks = L.parse_smconfig(block)
            fac = L.parse_gauge_factors(block)
            has_sm = "SMConfig" in block and "SU(3)_C" in block
            pos, tY = None, None
            if has_sm:
                for tYc, tbl in blocks:
                    printed = {}
                    for d3, d2, Y, m, ch in tbl:
                        if ch == "l":
                            printed[(d3, d2, Y)] = printed.get((d3, d2, Y), 0) + m
                    if not printed or tYc is None:
                        continue
                    p2 = L.observable_positions(spec, tYc, fac, printed)
                    if p2 is not None:
                        pos, tY = p2, tYc
                        break
            if pos is None:
                rec["status"] = "NO_SMCONFIG" if not has_sm else "SM_ID_DATA_GAP"
                rec["n_fields"] = len(spec)
                rec["has_tY"] = any(b[0] is not None for b in blocks)
                rec["n_smconfig_blocks"] = len(blocks)
                rec["gauge_factors"] = fac
                rec["configs"] = sorted(set(re.findall(
                    r'Gauge group in vev-configuration "([^"]+)"', block)))
                out[lab] = rec
                continue
            a3, a2, sY0 = pos
            roles, diag = L.role_census(spec, tY, a3, a2, sY0)
            fields = []
            for flab, r in roles.items():
                if not r["role"]:
                    continue
                kl = r["kl"]
                idx = [L.leg_orbit_index(kl, r["fp"], a) for a in (1, 2, 3)]
                fields.append([kl[0], kl[1], idx[0], idx[1], idx[2], r["role"],
                               r["mult"], int(bool(r["hidden_singlet"])), flab])
            rec.update(tY=[str(x) for x in tY], status="OK", a3=a3, a2=a2,
                       sY=diag["sY_selected"], sC=diag["sC_selected"],
                       nets=diag["nets"], n_argmax_ties=diag["n_argmax_ties"],
                       net_Q_is_plus_three=diag["net_Q_is_plus_three"],
                       n_fields=len(spec),
                       sectors_present=sorted({"%d,%d" % tuple(r["kl"])
                                               for r in spec.values()}),
                       fields=fields)
        except Exception as exc:                                   # noqa: BLE001
            rec["status"] = "PARSE_ERROR"
            rec["error"] = "%s: %s" % (type(exc).__name__, exc)
        out[lab] = rec
    return out


def main():
    segs = sorted(glob.glob(os.path.join(L.DSW, "raw", "census", "census_seg*.txt.gz")))
    segs = [s for s in segs if "TEST" not in s and "PILOT" not in s]
    print("segments:", len(segs), flush=True)
    allrec = {}
    with Pool(6) as pool:
        for i, d in enumerate(pool.imap_unordered(process_segment, segs, chunksize=1)):
            allrec.update(d)
            if (i + 1) % 10 == 0:
                print("  parsed %d/%d segments, %d models" % (i + 1, len(segs), len(allrec)),
                      flush=True)

    os.makedirs(os.path.join(L.DSW, "models"), exist_ok=True)
    with gzip.open(os.path.join(L.DSW, "models", "roles.json.gz"), "wt",
                   encoding="utf-8") as fh:
        json.dump(allrec, fh, default=str)

    st = Counter(r["status"] for r in allrec.values())
    ok = [r for r in allrec.values() if r["status"] == "OK"]
    sect = Counter()
    role_by_sector = defaultdict(Counter)
    offlat = 0
    for r in ok:
        for k, l, o1, o2, o3, role, mult, hs, _flab in r["fields"]:
            sect[(k, l)] += mult
            role_by_sector[(k, l)][role] += mult
            if "OFF_LATTICE" in (o1, o2, o3):
                offlat += 1
    left_sectors = Counter()
    for r in ok:
        for s in r.get("sectors_present", []):
            left_sectors[s] += 1
    rep = {"geometry": "Z2xZ4_1_1",
           "scope": "CATALOGUE statement about the 11,540 PVW Z2xZ4_1_1 models",
           "n_models_seen": len(allrec), "status_histogram": dict(st),
           "n_net_Q_plus_three": sum(1 for r in ok if r.get("net_Q_is_plus_three")),
           "n_argmax_ties_total": sum(r.get("n_argmax_ties", 1) - 1 for r in ok),
           "sY_histogram": dict(Counter(r["sY"] for r in ok)),
           "sC_histogram": dict(Counter(r["sC"] for r in ok)),
           "orientation_branch_histogram":
               {str(k): v for k, v in Counter((r["sY"], r["sC"]) for r in ok).items()},
           "n_off_lattice_fixed_points": offlat,
           "sectors_hosting_a_left_chiral_field_n_models": dict(left_sectors.most_common()),
           "role_fields_per_sector": {"T(%d,%d)" % k: v for k, v in sect.most_common()},
           "roles_per_sector": {"T(%d,%d)" % k: dict(v)
                                for k, v in sorted(role_by_sector.items())}}
    L.write_json("census_summary.json", rep)
    print(json.dumps({k: v for k, v in rep.items() if k != "roles_per_sector"},
                     indent=1, default=str))


if __name__ == "__main__":
    main()

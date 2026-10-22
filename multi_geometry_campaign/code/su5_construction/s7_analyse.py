"""S7 -- analyse the CB-1 candidate spectra: the SM battery and the two-10-role test.

Per model, using the CERTIFIED Z3xZ6 census machinery (z36sw_lib: closure-selected
hypercharge block, four-branch orientation, role census with role = None on hidden
non-singlets -- the G-9 convention of THIS geometry), we record

  * gauge group in TestConfig and SMConfig, k_Y = 2|t_Y|^2
  * net generations, convention-blind (all four (sY, sC) branches emitted; the physical
    one is the branch with net Q = +3, exactly as the certified census does)
  * every T(0,1)/T(0,5) cell: the 10-role RECORD multiplicity, and -- because the
    Z3xZ6 role census sets role = None on hidden non-singlets (gate D5 defect B) --
    a SECOND pass with that gate removed, which is where a 10-role in a hidden SU(2)
    DOUBLET would show up.  Both counts are reported; neither is called "the" count.
  * the STATE count per (cell, role): records expanded by their hidden multiplicity.

Output: data/su5_construction/s7_analyse.json
"""
from __future__ import annotations

import glob
import gzip
import json
import os
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F

import sc_lib as S

sys.path.insert(0, os.path.join(S.CAMP, "code", "z3xz6_sweep"))
import z36sw_lib as SW                                                    # noqa: E402

RAW = os.path.join(S.OUT, "raw")
# z36sw_lib._MARK hard-codes the MSSM label prefix of the catalogue; the CB models carry
# their own prefix, so the model splitter is re-expressed here (same cut points).
import re as _re                                                          # noqa: E402
_MARKCB = _re.compile(r"/>\s*/(CB\d+)>")


def split_models_cb(text):
    marks = [(m.start(), m.group(1)) for m in _MARKCB.finditer(text)]
    out = []
    for i, (pos, lab) in enumerate(marks):
        end = marks[i + 1][0] if i + 1 < len(marks) else len(text)
        out.append((lab, text[pos:end]))
    return out
SECT = {(0, 1), (0, 5)}
TEN = ("Q", "uc", "ec", "Qbar", "ucbar", "ecbar")


def hidden_mult(rep, a3, a2):
    m = 1
    for k, d in enumerate(rep):
        if k not in (a3, a2):
            m *= abs(d)
    return m


def main():
    out = dict(models=0, parsed=0, no_smconfig=0, kY=Counter(), groups=Counter(),
               net_gen=Counter(), gen3=[], ten_cells=Counter(),
               ten_record_mult=Counter(), ten_state_mult=Counter(),
               ten_hidden_nonsinglet=[], doubled=[], per_model=[])
    for pth in sorted(glob.glob(os.path.join(RAW, "*.txt.gz"))):
        txt = gzip.open(pth, "rt", errors="replace").read()
        for lab, blk in split_models_cb(txt):
            out["models"] += 1
            pre = blk.split("Analyzing")[0]
            spec = SW.parse_spectrum(pre)
            gauge = SW.parse_gauge(blk)
            if not spec or gauge.get("tY") is None:
                out["no_smconfig"] += 1
                continue
            pos = SW.observable_positions(spec, gauge["tY"], gauge)
            if pos is None:
                out["no_smconfig"] += 1
                continue
            a3, a2, sY, tv, bi = pos
            roles, diag = SW.role_census(spec, tv, a3, a2, sY)
            out["parsed"] += 1
            out["kY"][str(2 * S.dot(tv, tv))] += 1
            out["groups"][str(gauge.get("gauge_sm"))] += 1
            nets = diag.get("nets") if isinstance(diag, dict) else None
            netQ = diag.get("net", {}) if isinstance(diag, dict) else {}
            rec = dict(model=lab, kY=str(2 * S.dot(tv, tv)),
                       gauge_sm=gauge.get("gauge_sm"),
                       gauge_test=gauge.get("gauge_test"), nets=nets)
            # ---- per-cell 10-role counts at the loophole sectors -------------------
            cells = defaultdict(list)
            for l, r in spec.items():
                if r["p_sh"] is None or r["chir"] != "l" or not r["rep"]:
                    continue
                if tuple(r["kl"]) not in SECT:
                    continue
                cells[(tuple(r["kl"]), tuple(r["fp"]) if r["fp"] else None)].append(r)
            cellinfo = []
            for key, recs in cells.items():
                bycount = Counter()
                bystate = Counter()
                for r in recs:
                    ro = roles.get(r["label"], {}).get("role")
                    hm = hidden_mult(r["rep"], a3, a2)
                    if ro is None:
                        # gate removed: recompute the role ignoring the hidden-singlet gate
                        Y = sY * S.dot(tv, r["p_sh"])
                        c3, c2 = r["rep"][a3], abs(r["rep"][a2])
                        ro = SW_role_from_qn(c3, c2, Y, diag)
                    if ro in TEN:
                        bycount[ro] += 1
                        bystate[ro] += hm
                        if hm > 1:
                            out["ten_hidden_nonsinglet"].append(
                                dict(model=lab, kl=list(key[0]), fp=key[1], role=ro,
                                     hidden_mult=hm, p_sh=[str(x) for x in r["p_sh"]],
                                     p2=str(S.dot(r["p_sh"], r["p_sh"]))))
                if bycount:
                    out["ten_cells"][str(key[0])] += 1
                    for ro, c in bycount.items():
                        out["ten_record_mult"][c] += 1
                    for ro, c in bystate.items():
                        out["ten_state_mult"][c] += 1
                    if max(bystate.values()) >= 2:
                        out["doubled"].append(dict(model=lab, kl=list(key[0]),
                                                   fp=key[1],
                                                   record=dict(bycount),
                                                   state=dict(bystate)))
                    cellinfo.append(dict(kl=list(key[0]), fp=key[1],
                                         records=dict(bycount), states=dict(bystate)))
            rec["loophole_cells"] = cellinfo
            out["per_model"].append(rec)

    for k in ("kY", "groups", "net_gen", "ten_cells", "ten_record_mult",
              "ten_state_mult"):
        out[k] = {str(a): b for a, b in out[k].items()}
    with open(os.path.join(S.OUT, "s7_analyse.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({k: v for k, v in out.items()
                      if k not in ("per_model", "doubled", "ten_hidden_nonsinglet")},
                     indent=1))
    print("doubled cells:", len(out["doubled"]),
          " ten-role hidden non-singlets:", len(out["ten_hidden_nonsinglet"]))


def SW_role_from_qn(c3, c2, Y, diag):
    """role from (SU(3) dim signed, SU(2) dim, Y) with the SAME table the certified
    census uses -- reproduced here ONLY to see hidden non-singlets, which that census
    deliberately drops (gate D5 defect B).  Sign conventions come from `diag`."""
    a = abs(c3)
    if a == 3 and c2 == 2 and Y == F(1, 6):
        return "Q"
    if a == 3 and c2 == 1 and Y == F(-2, 3):
        return "uc"
    if a == 1 and c2 == 1 and Y == F(1):
        return "ec"
    if a == 3 and c2 == 2 and Y == F(-1, 6):
        return "Qbar"
    if a == 3 and c2 == 1 and Y == F(2, 3):
        return "ucbar"
    if a == 1 and c2 == 1 and Y == F(-1):
        return "ecbar"
    return None


if __name__ == "__main__":
    main()

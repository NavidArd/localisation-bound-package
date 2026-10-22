"""su5_ten_bound RECON -- premise-memo evidence (CLAUDE.md 5.20 (i)/(ii)/(iii)).

Answers, on real on-disk Z2xZ4_1_1 census output, BEFORE any derivation is attempted:
  P1  Is p_sh mod Lambda constant across all fields of a cell?  (If not, C2 is false and
      the whole local-shift route dies here.)
  P2  Does the local group computed from p_sh mod Lambda reproduce, when intersected over
      ALL cells of a model, orbifolder's OWN printed 4D gauge group?  (the known-invariant
      check -- 5.19 rule (i))
  P3  What ARE the local groups?  Do local SU(5)/SO(10)/E6 factors containing the SM
      actually occur?
  P4  Level-1 unitarity: are all SM-role states level-1 weights of every local factor?
  P5  k_Y = 2 |t_Y|^2 -- is the hypercharge GUT-normalised (5/3) in this catalogue?

Output: data/su5_ten_bound/recon.json
"""
from __future__ import annotations

import glob
import gzip
import json
import os
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
CAMP = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(CAMP, "code", "z2xz4_sweep"))
import s5lib as S                                                       # noqa: E402
import z24sw_lib as SW                                                  # noqa: E402

OUT = os.path.join(CAMP, "data", "su5_ten_bound")
CENSUS = os.path.join(CAMP, "data", "z2xz4_sweep", "raw", "census")
os.makedirs(OUT, exist_ok=True)

ROLES10 = ("Q", "uc", "ec")
ROLES5 = ("dc", "L", "Hu")


def gauge_group_from_shifts(V1, V2, W):
    """4D group = roots invariant under EVERY space-group generator's shift."""
    R = S.e8xe8_roots()
    for V in [V1, V2] + list(W):
        R = [a for a in R if S.dot(a, V).denominator == 1]
    return R


def printed_4d_group(block):
    m = None
    import re
    m = re.search(r"Gauge group:\s*(.*)", block)
    if not m:
        m = re.search(r'Gauge group in vev-configuration "TestConfig1":\s*(.*)', block)
    return m.group(1).strip() if m else None


def name_of_rootset(R):
    if not R:
        return "U(1)^16"
    Sp = S.simple_roots(R)
    comps, adj = S._components(Sp)
    facs = []
    for c in comps:
        kind, r = S._classify(c, adj, Sp)
        facs.append(S.group_name(kind, r))
    nm = " x ".join(sorted(facs))
    if len(Sp) < 16:
        nm += " x U(1)^%d" % (16 - len(Sp))
    return nm


def main(nmodels=12):
    rep = dict(models=[], P1_violations=0, P1_cells=0,
               local_group_names=Counter(), sm_in_factor=Counter(),
               level_violations=0, level_checked=0, kY=Counter(),
               role10_local_rep=Counter(), role5_local_rep=Counter(),
               gauge_match=Counter())
    done = 0
    for pth in sorted(glob.glob(os.path.join(CENSUS, "*.txt.gz"))):
        txt = gzip.open(pth, "rt", errors="replace").read()
        for lab, blk in SW.split_models(txt):
            if done >= nmodels:
                break
            pre = blk.split("Analyzing")[0]
            spec = SW.parse_spectrum(pre, M=2, N=4)
            factors = SW.parse_gauge_factors(blk)
            cfgs = SW.parse_smconfig(blk)
            if not cfgs or cfgs[0][0] is None or not spec:
                continue
            tY, tbl = cfgs[0]
            printed = defaultdict(int)
            for (c3, c2, Y, mult, ch) in tbl:
                printed[(c3, c2, Y)] += mult
            pos = SW.observable_positions(spec, tY, factors, dict(printed))
            if pos is None:
                continue
            roles, diag = SW.role_census(spec, tY, *pos)
            done += 1
            mrep = dict(label=lab, file=os.path.basename(pth))

            # ---- P5  k_Y
            kY = 2 * S.dot(tY, tY)
            mrep["k_Y"] = str(kY)
            rep["kY"][str(kY)] += 1

            # ---- P1  p_sh mod Lambda constant per cell
            cells = defaultdict(list)
            for l, r in spec.items():
                if r["p_sh"] is None or r["chir"] != "l":
                    continue
                cells[(tuple(r["kl"]), tuple(r["fp"]) if r["fp"] else None)].append(r)
            nbad = 0
            for key, recs in cells.items():
                ks = {S.reduce_mod_lattice(r["p_sh"]) for r in recs}
                rep["P1_cells"] += 1
                if len(ks) != 1:
                    nbad += 1
            rep["P1_violations"] += nbad
            mrep["P1_bad_cells"] = nbad

            # ---- P2  intersect local groups over all cells -> 4D group
            allR = S.e8xe8_roots()
            inter = allR
            for key, recs in cells.items():
                V = recs[0]["p_sh"]
                inter = [a for a in inter if S.dot(a, V).denominator == 1]
            mrep["gauge_from_cells"] = name_of_rootset(inter)
            mrep["gauge_printed"] = printed_4d_group(blk)
            rep["gauge_match"][mrep["gauge_from_cells"]] += 1

            # ---- P3/P4  local groups + local reps of SM-role states
            percell = []
            for key, recs in sorted(cells.items(), key=lambda kv: str(kv[0])):
                (kl, fp) = key
                if kl == (0, 0):
                    continue
                rrole = {r["label"]: roles.get(r["label"], {}).get("role") for r in recs}
                if not any(rrole.values()):
                    continue
                V = recs[0]["p_sh"]
                g = S.analyse_group(V)
                S.attach_orders(g)
                rep["local_group_names"][g["name"]] += 1
                # which factor carries SU(3)_c ?  the factor whose simple roots have
                # nonzero pairing with the state's colour weight -- identified below by
                # the factor under which a coloured role is non-singlet.
                ent = dict(kl=kl, fp=fp, group=g["name"], roles={})
                for r in recs:
                    ro = rrole.get(r["label"])
                    if ro is None:
                        continue
                    tags = []
                    for f in g["factors"]:
                        rr = S.rep_of(r["p_sh"], f, g["simple"])
                        rep["level_checked"] += 1
                        if not rr["ok"]:
                            tags.append((f["name"], "NONINT"))
                            continue
                        if not rr["level1_ok"]:
                            rep["level_violations"] += 1
                        if any(rr["hw"]):
                            tags.append((f["name"], rr.get("name") or str(rr["hw"])))
                    ent["roles"].setdefault(ro, []).append(tags)
                    key2 = " + ".join("%s:%s" % t for t in tags) or "all-singlet"
                    if ro in ROLES10:
                        rep["role10_local_rep"][key2] += 1
                    elif ro in ROLES5:
                        rep["role5_local_rep"][key2] += 1
                percell.append(ent)
            mrep["n_cells_with_roles"] = len(percell)
            mrep["sample_cells"] = percell[:6]
            rep["models"].append(mrep)
        if done >= nmodels:
            break
    for k in ("local_group_names", "sm_in_factor", "kY", "role10_local_rep",
              "role5_local_rep", "gauge_match"):
        rep[k] = dict(rep[k].most_common(25))
    with open(os.path.join(OUT, "recon.json"), "w") as fh:
        json.dump(rep, fh, indent=1, default=str)
    print(json.dumps({k: v for k, v in rep.items() if k != "models"},
                     indent=1, default=str)[:6000])
    print("--- first model ---")
    print(json.dumps(rep["models"][0], indent=1, default=str)[:4000])


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 12)

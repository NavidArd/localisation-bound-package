"""S11 -- the COMPLETE hypercharge battery (CLAUDE.md 5.21) over all CB models.

The 5.21 lesson exactly: an embedding can reproduce every Standard-Model value and still
be ANOMALOUS.  t_Y here is traceless inside the block-1 SU(5) and gives k_Y = 5/3 and the
right hypercharges on every state -- but the four-dimensional spectrum does NOT fill
complete SU(5) multiplets (that is what breaking SU(5) means), so

    Tr Y  =  sum over all left-chiral STATES of Y

is NOT automatically zero, and where it is non-zero the U(1)_Y is anomalous and is not a
hypercharge at all.  Every defining condition is therefore enumerated and checked:

  k_Y = 2|t_Y|^2 = 5/3                     (exact, by construction)
  t_Y . (SU(3), SU(2) simple roots) = 0    (exact, by construction)
  Tr Y = 0            <- the condition the construction does NOT get for free
  Tr Y^3 = 0          <- reported alongside (pure-U(1)_Y cubic anomaly)
  net 3 generations in Q, u^c, e^c         (convention-blind, over STATES)
  a doubled 10-role at a T(0,1)/T(0,5) cell

Output: data/su5_construction/s11_battery2.json
"""
import glob, gzip, json, os, sys
from collections import Counter
from fractions import Fraction as F
import sc_lib as S
sys.path.insert(0, os.path.join(S.CAMP, "code", "z3xz6_sweep"))
import z36sw_lib as SW
from s8_battery import P, DELTA, TY, A3, A2, split_models_cb, su3_dim, su2_dim, role_of, RAW, SECT


def main():
    out = dict(models=0, funnel=Counter(), trY_hist=Counter(), passing=[], rows=[])
    for pth in sorted(glob.glob(os.path.join(RAW, "*.txt.gz"))):
        txt = gzip.open(pth, "rt", errors="replace").read()
        for lab, blk in split_models_cb(txt):
            out["models"] += 1
            spec = SW.parse_spectrum(blk.split("Analyzing")[0])
            if not spec:
                continue
            left = [r for r in spec.values()
                    if r["chir"] == "l" and r["p_sh"] and r["rep"]]
            trY = trY3 = F(0)
            nets = {1: Counter(), -1: Counter()}
            hits = 0
            for r in left:
                d3, _ = su3_dim(r["p_sh"])
                d2 = su2_dim(r["p_sh"])
                tot = 1
                for d in r["rep"]:
                    tot *= abs(d)
                y = S.dot(TY, r["p_sh"])
                trY += tot * y
                trY3 += tot * y ** 3
                hid = F(tot, d3 * d2)
                for sY in (1, -1):
                    ro = role_of(d3, d2, y, sY, 1)
                    if ro is None:
                        continue
                    if ro.endswith("bar"):
                        nets[sY][ro[:-3]] -= int(hid)
                    else:
                        nets[sY][ro] += int(hid)
                if (tuple(r["kl"]) in SECT and S.dot(r["p_sh"], r["p_sh"]) == F(31, 18)
                        and abs(S.dot(r["p_sh"], DELTA)) == 1):
                    d3b, _ = su3_dim(r["p_sh"])
                    ro = (role_of(d3b, d2, y, 1, 1) or role_of(d3b, d2, y, -1, 1))
                    if ro in ("Q", "uc", "ec", "Qbar", "ucbar", "ecbar"):
                        hits += 1
            sYp = 1 if nets[1].get("Q") == 3 else (-1 if nets[-1].get("Q") == 3 else 0)
            gen3Q = sYp != 0
            gen3all = gen3Q and all(nets[sYp].get(k) == 3 for k in ("Q", "uc", "ec"))
            out["trY_hist"][str(trY)] += 1
            tag = "trY0" if trY == 0 else "trYx"
            out["funnel"][(tag, "genQ3" if gen3Q else "genQx",
                           "gen3all" if gen3all else "-",
                           "dbl" if hits else "-")] += 1
            out["rows"].append(dict(model=lab, trY=str(trY), trY3=str(trY3),
                                    netQ={str(k): dict(v) for k, v in nets.items()},
                                    doubled=hits))
            if trY == 0 and gen3Q and hits:
                out["passing"].append(lab)
    out["funnel"] = {"|".join(k): v for k, v in sorted(out["funnel"].items())}
    out["trY_hist"] = dict(out["trY_hist"])
    json.dump(out, open(os.path.join(S.OUT, "s11_battery2.json"), "w"), indent=1)
    print(json.dumps(out["funnel"], indent=1))
    print("Tr Y = 0 models:", out["trY_hist"].get("0", 0), "of", out["models"])
    print("Tr Y = 0 AND net Q = 3 AND doubled:", len(out["passing"]), out["passing"][:10])


if __name__ == "__main__":
    main()

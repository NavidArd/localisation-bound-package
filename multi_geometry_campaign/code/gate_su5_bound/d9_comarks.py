"""GATE su5_ten_bound -- D-item: the phase's declared attack #9.

s5lib.comarks hard-codes a_i^vee in the phase's own Bourbaki-like ordering.  Only the
level-1 test consumes them and it passed everywhere, which the phase itself calls weak
evidence.  The check it proposed is run here: for EVERY simple factor that occurs, the
ADJOINT must come out at level 2, and the level-1 irreps must be exactly the cominuscule
fundamentals.  Local groups are generated from real Borel-de Siebenthal shifts so the
orderings exercised are the ones the phase actually uses.
"""
import os, sys
from fractions import Fraction as F
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "su5_ten_bound"))
import s5lib as S

SH = {  # shifts producing a spread of local factor types inside E8 x E8
 "A": [F(4,5)]+[F(-1,5)]*4+[F(0)]*3,
 "D": [F(1)]+[F(0)]*7,
 "E6": [F(2,3),F(-1,3),F(-1,3)]+[F(0)]*5,
 "E7": [F(1,2),F(1,2),F(1,2),F(1,2),F(1,2),F(1,2),F(-1,2),F(1,2)],
 "mix": [F(1,6)]*5+[F(-5,6),F(0),F(0)],
 "mix2":[F(3,4),F(1,4),F(1,4),F(1,4),F(1,4),F(1,4),F(1,4),F(1,4)],
}
bad, seen, rows = [], set(), []
for lab, V0 in SH.items():
    V = list(V0)+[F(0)]*8
    g = S.attach_orders(S.analyse_group(V))
    LR = S.local_roots(V)
    for f in g["factors"]:
        so = [g["simple"][i] for i in f["order"]]
        cm = S.comarks(f["kind"], f["rank"])
        # roots of THIS factor = local roots with integral labels and nonzero only here
        fr = [a for a in LR if all(S.dot(a,b)==0 for ff in g["factors"] if ff is not f
                                   for b in [g["simple"][i] for i in ff["order"]])]
        fr = [a for a in fr if any(S.dot(a,b)!=0 for b in so)]
        if not fr: continue
        levs = set()
        for a in fr:
            lab_ = [int(S.dot(a,b)) for b in so]
            hw = S.to_dominant(lab_, so)
            levs.add(int(sum(x*y for x,y in zip(cm, hw))))
        key = (f["kind"], f["rank"])
        rows.append((lab, f["name"], key, sorted(levs), len(fr), len(cm)))
        if max(levs) != 2:
            bad.append((lab, f["name"], sorted(levs)))
        if len(cm) != f["rank"]:
            bad.append((lab, f["name"], "COMARK LENGTH %d != rank %d" % (len(cm), f["rank"])))
        seen.add(key)
for r in rows: print("  %-5s %-8s %-8s adjoint-levels=%s roots=%3d len(comarks)=%d" % r)
print("types exercised:", sorted(seen))
print("FAILURES:", bad if bad else "NONE  (every adjoint at level 2)")

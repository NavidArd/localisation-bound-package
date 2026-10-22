"""V6 -- why the certified catalogue never realises the loophole, structurally.

The gate derived that a doubled 10-role at T(0,1)/T(0,5) requires the class of V2 to be
Aut(Lambda)-equivalent to one of exactly TWO classes.  The size of the |p|^2 = 31/18 shell
of the coset V2 + Lambda is an Aut(Lambda)-INVARIANT of the class and needs no knowledge
of where the model's Standard Model sits.  If no catalogue model's T(0,1) coset has the
shell size of either target class, the catalogue cannot host the configuration -- an
SM-free explanation of the measured absence.
"""
import json, os, sys
from collections import Counter
from fractions import Fraction as F
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gs_lib as G
from g3_l1 import read_catalogue
CAT = None


def _init():
    global CAT
    CAT = read_catalogue()


def do(lab):
    rows = CAT[lab]
    V2 = rows[1]
    out = {}
    for l in (1, 5):
        Vg = G.smul(l, V2)
        sh = G.full_shell(Vg, F(31, 18))
        out[l] = len(sh)
    return lab, out


def main():
    _init()
    labs = sorted(CAT)
    with Pool(6, initializer=_init) as pl:
        res = pl.map(do, labs, chunksize=20)
    targets = {}
    g2 = json.load(open(os.path.join(G.OUT, "g2_classes.json")))
    for r in g2["representatives"]:
        p = [G.rat(x) for x in r["p"]]
        targets[r["case"]] = len(G.full_shell(p, F(31, 18)))
    print("target class shell sizes at |p|^2 = 31/18:", targets)
    h1 = Counter(o[1] for _, o in res)
    h5 = Counter(o[5] for _, o in res)
    print("catalogue T(0,1) shell-size histogram:", dict(sorted(h1.items())))
    print("catalogue T(0,5) shell-size histogram:", dict(sorted(h5.items())))
    matches = [lab for lab, o in res
               if o[1] in targets.values() or o[5] in targets.values()]
    print("catalogue models whose T(0,1) or T(0,5) shell size matches a target:",
          len(matches), matches[:20])
    json.dump(dict(targets=targets, hist1=dict(h1), hist5=dict(h5),
                   n_models=len(res), matches=matches),
              open(os.path.join(G.OUT, "g11_catalogue.json"), "w"), indent=1)


if __name__ == "__main__":
    main()

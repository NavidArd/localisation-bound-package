"""Referee #8 item 2.3: do lattice-representative changes of V2 (V2 -> V2 + nu, nu in
Lambda, complete MI set passing) produce spectra OUTSIDE the nine (a,b) brother classes
of b1/b5 (V1 -> V1+lam, W3=W4 -> W3+mu)?  Builds, for CB1534, every MI-valid V2+nu with
nu in {+-e_i +- e_j}, one per (nu.V1 mod 1, nu.W3 mod 1) value, plus the nine (a,b)
references, for b2_run_orbifolder.sh.
Output: review/data/paperA_referee8/b11_decks.txt, b11_decks.json"""
import json, os, sys
_PKG = os.environ.get("PAPERA_PKG") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b1_brothers as B1
S = B1.S
OUT = B1.OUT
def main():
    meta, decks = {}, []
    seen = {}
    for nu in B1.small_lattice_vectors():
        nv2 = B1.add(B1.V2, nu)
        ok, _ = S.mi_report(B1.V1, nv2, [B1.Z, B1.Z, B1.W3, B1.W3, B1.Z, B1.Z])
        if not ok:
            continue
        key = (str(B1.frac1(S.dot(nu, B1.V1))), str(B1.frac1(S.dot(nu, B1.W3))))
        if key in seen:
            continue
        seen[key] = nu
        lab = "NU%d" % len(seen)
        decks.append(B1.deck(lab, B1.V1, nv2, B1.W3))
        meta[lab] = dict(nu=[str(x) for x in nu], nuV1=key[0], nuW=key[1])
    # the nine (a,b) references, one representative each
    ref = {"00": ([F(0)]*16, [F(0)]*16)}
    for lam in B1.small_lattice_vectors():
        for mu in [[F(0)]*16] + B1.small_lattice_vectors():
            a = B1.frac1(S.dot(lam, B1.V2)) * 3; b = B1.frac1(S.dot(mu, B1.V2)) * 3
            if a.denominator != 1 or b.denominator != 1:
                continue
            k = "%d%d" % (a, b)
            if k in ref:
                continue
            nv1, nw = B1.add(B1.V1, lam), B1.add(B1.W3, mu)
            if S.mi_report(nv1, B1.V2, [B1.Z, B1.Z, nw, nw, B1.Z, B1.Z])[0]:
                ref[k] = (lam, mu)
        if len(ref) == 9:
            break
    for k, (lam, mu) in sorted(ref.items()):
        nw = B1.add(B1.W3, mu)
        decks.append(B1.deck("AB" + k, B1.add(B1.V1, lam), B1.V2, nw))
        meta["AB" + k] = dict(lam=[str(x) for x in lam], mu=[str(x) for x in mu])
    open(os.path.join(OUT, "b11_decks.txt"), "w", newline="\n").write("".join(decks))
    json.dump(meta, open(os.path.join(OUT, "b11_decks.json"), "w"), indent=1)
    print(len(decks), sorted(ref), {k: (v["nuV1"], v["nuW"]) for k, v in meta.items() if k.startswith("NU")})
main()

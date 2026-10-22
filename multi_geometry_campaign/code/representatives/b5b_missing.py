"""Referee #8 item 2.3: the (model, a, b) brother slots that b5's small search (vectors
+-e_i +- e_j) left unbuilt.  Searches a larger lattice set (all E8 vectors with entries in
{-1,0,1} and two or four non-zero entries, plus the 128 spinor-type vectors (+-1/2)^8 of
even minus count, in either block) for the PAIR (lam, mu) that fills the slot:

    V1 -> V1 + lam,  W3 = W4 -> W3 + mu,   3 lam.V2 = a (mod 3),  3 mu.V2 = b (mod 3),

with the complete orbifolder modular-invariance set (sc_lib.mi_report) on the combined data.
If none passes, the obstruction is recorded.

REPAIR (package v2.2, 2026-10-22).  The version first deposited searched lam only, with
lam.V2 = a/3, and wrote mu = 0 for every slot.  For the 232 slots with b = 0 that is
correct.  For the 80 slots with a = 0 and b != 0 (40 class-B pairs x (0,1), (0,2)) it
"found" a lam with lam.V2 = 0 and mu = 0 and labelled the deck (0, b): those decks realise
the (0,0) phase class and their spectra are record-identical to the pair's (0,0) spectrum.
Found by the independent v10 verification (paper15/review/paperA_v10_verification.md
sec. 1.3); the superseded decks, outputs and build records are kept, labelled, in
`brothers/superseded_v2.2/` (see the NOTICE there).  This version:
  * searches mu whenever b != 0 (and lam, starting from lam = 0, whenever a = 0);
  * for b = 0 keeps mu = 0 and the original lam order, so the 232 correct decks are
    reproduced vector for vector;
  * RE-COMPUTES, for every deck it writes, (3 lam.V2 mod 3, 3 mu.V2 mod 3) from the deck's
    own vectors minus the parent's and asserts that it equals the slot's (a, b) label
    (every generated deck carries a check of the invariant that defines its slot);
  * labels the decks whose slot needed a mu search `Y<parent>a<a>b<b>` (not `X...`), so that
    no label of a superseded deck is ever reused for a different object.  The search order
    for those slots (lam over [0] + candidates with lam.V2 = 0 mod 1, outer; mu over the
    candidates with 3 mu.V2 = b mod 3, inner; first MI-valid pair) is that of the
    verification's builder `v1b_missing_b_slots.py`, so the decks written are the 80 it
    built and ran (deposited verbatim as `decks_902.txt` / `out_902.txt.gz`).
The b = 0 decks go to brothers/decks_900.txt, decks_901.txt (200 per segment); the mu-search
decks to the segments that follow (decks_902.txt).  The merged build record that
b10_brother_analysis.py reads (`build_all.json`: b5's `build.json` with this script's decks
merged in) is written here too.
Output: brothers/build_missing.json, brothers/build_all.json, brothers/decks_9NN.txt
"""
import json
import os
_PKG = os.environ.get("PAPERA_PKG") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
import sys
from itertools import combinations, product
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b5_build_brother_census as B5  # noqa: E402

S = B5.S
Z = [F(0)] * 16


def big_cands():
    out = []
    for blk in (0, 8):
        for k in (2, 4):
            for idx in combinations(range(8), k):
                for sg in product((1, -1), repeat=k):
                    v = [F(0)] * 16
                    for i, s in zip(idx, sg):
                        v[blk + i] = F(s)
                    out.append(v)
        for sg in product((1, -1), repeat=8):
            if sg.count(-1) % 2 == 0:
                v = [F(0)] * 16
                for i, s in enumerate(sg):
                    v[blk + i] = F(s, 2)
                out.append(v)
    return out


def phase3(v, V2):
    """3 v.V2 mod 3 as an int in {0,1,2}, or None if 3 v.V2 is not an integer."""
    t = 3 * S.dot(v, V2)
    if t.denominator != 1:
        return None
    return int(t) % 3


def slot_label(vec_V1, vec_W3, parent_V1, parent_W3, V2):
    """(a, b) re-computed from a deck's own vectors: the invariant that defines its slot."""
    lam = [x - y for x, y in zip(vec_V1, parent_V1)]
    mu = [x - y for x, y in zip(vec_W3, parent_W3)]
    return phase3(lam, V2), phase3(mu, V2)


def main():
    bj = json.load(open(os.path.join(B5.OUT, "build.json")))
    miss = bj["missing"]
    A = dict(B5.read_decks(os.path.join(B5.CAMP, "su5_construction", "s5_models.txt")))
    Bd = dict(B5.read_decks(os.path.join(B5.CAMP, "gate_su5_construction",
                                         "caseB_models.txt")))
    C = big_cands()
    xdecks, ydecks, meta, still = [], [], {}, []
    for lab, a, b in miss:
        cls = "A" if lab in A else "B"
        rows = A.get(lab) or Bd.get(lab)
        V1, V2, W = rows[0], rows[1], rows[4]
        if b == 0:
            lams = [l for l in C if phase3(l, V2) == a]       # v2.2 order, unchanged
            if a == 0:
                lams = [Z] + lams
            mus = [Z]
        else:
            lams = [Z] + [l for l in C if phase3(l, V2) == a] if a == 0 else \
                [l for l in C if phase3(l, V2) == a]
            mus = [u for u in C if phase3(u, V2) == b]
        found = None
        for lam in lams:
            nv1 = B5.add(V1, lam)
            for mu in mus:
                nw = B5.add(W, mu)
                ok, _bad = S.mi_report(nv1, V2, [Z, Z, nw, nw, Z, Z])
                if ok:
                    found = (lam, mu, nv1, nw)
                    break
            if found:
                break
        if not found:
            fails = set()
            for lam in lams:
                for mu in mus:
                    nw = B5.add(W, mu)
                    fails |= set(S.mi_report(B5.add(V1, lam), V2, [Z, Z, nw, nw, Z, Z])[1])
            still.append(dict(model=lab, a=a, b=b, W_zero=all(x == 0 for x in W),
                              failing=sorted(fails)))
            continue
        lam, mu, nv1, nw = found
        # the slot's defining invariant, re-computed from the deck's own vectors
        assert slot_label(nv1, nw, V1, W, V2) == (a, b), (lab, a, b)
        mu_search = (b != 0)
        nl = "%s%sa%db%d" % ("Y" if mu_search else "X", lab, a, b)
        txt = ("begin model\nLabel:%s\nSpaceGroup:Geometry/Geometry_Z3xZ6_1_1.txt"
               "\nLattice:E8xE8\nShifts and Wilsonlines:\n%s\nend model\n"
               % (nl, "\n".join(B5.fmt(r) for r in [nv1, V2, Z, Z, nw, nw, Z, Z])))
        (ydecks if mu_search else xdecks).append(txt)
        meta[nl] = dict(parent=lab, cls=cls, a=a, b=b, lam=[str(x) for x in lam],
                        mu=[str(x) for x in mu])
    seg = 0
    for group in (xdecks, ydecks):
        for k in range(0, len(group), B5.SEG):
            with open(os.path.join(B5.OUT, "decks_9%02d.txt" % seg), "w",
                      newline="\n") as fh:
                fh.write("".join(group[k:k + B5.SEG]))
            seg += 1
    n_ok = sum(1 for m in meta.values())
    json.dump(dict(built=len(meta), built_lam_only=len(xdecks), built_with_mu=len(ydecks),
                   label_invariant_checked=n_ok, meta=meta, still_missing=still),
              open(os.path.join(B5.OUT, "build_missing.json"), "w"), indent=1)
    allmeta = dict(bj["meta"])
    allmeta.update(meta)
    json.dump(dict(n_models=bj["n_models"], n_decks=bj["n_decks"], missing=bj["missing"],
                   meta=allmeta, missing_after_wide_search=still),
              open(os.path.join(B5.OUT, "build_all.json"), "w"))
    print("missing", len(miss), "built", len(meta), "(lam only %d, with mu %d)"
          % (len(xdecks), len(ydecks)), "still", len(still), still[:3],
          "segments", seg, "total decks", len(allmeta))


if __name__ == "__main__":
    main()

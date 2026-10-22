"""GATE su5_ten_bound -- independent massless-shell enumerator.

Written fresh (the phase's `su5_ten_bound/shell.py` is NOT imported) so that the loophole
recount is not a re-run of the code under review.  The phase's own defect log records that
its first version applied the D8 even-sum condition to x rather than to x - V, returning an
EMPTY shell for every non-trivial coset; this implementation is checked against three
independent invariants in `selftest()`:

  (i)   V = 0 gives the E8 theta series 1, 240, 2160, 6720 at norm^2 0, 2, 4, 6;
  (ii)  E8 x E8 at V = 0 gives 1, 480, 61920 at norm^2 0, 2, 4;
  (iii) any vector x of the coset is found in its own shell.

Lambda_{E8} = D8 u (D8 + s), s = (1/2)^8, D8 = integer vectors of EVEN sum, so
    x in V + Lambda_{E8}  iff  x - V - h in Z^8 with even coordinate sum, h in {0, s}
-- the parity condition is RELATIVE to the coset representative V.
"""
from __future__ import annotations

import itertools
from collections import defaultdict
from fractions import Fraction as F


def e8_coset_shell(V8, cmax):
    """Every x in V8 + Lambda_{E8} with |x|^2 <= cmax.  Returns [(x, |x|^2)]."""
    sV = sum(V8)
    out = set()
    for h in (F(0), F(1, 2)):
        ranges = []
        for j in range(8):
            base = (V8[j] + h) % 1
            vals = [base + k for k in range(-8, 9) if (base + k) ** 2 <= cmax]
            if not vals:
                ranges = None
                break
            ranges.append(vals)
        if ranges is None:
            continue
        for cand in itertools.product(*ranges):
            n = sum(x * x for x in cand)
            if n > cmax:
                continue
            d = sum(cand) - sV
            if d.denominator != 1 or int(d) % 2 != 0:
                continue
            out.add((tuple(cand), n))
    return sorted(out)


def shell(V16, c):
    """Every p in V16 + Lambda_{E8xE8} with |p|^2 == c exactly."""
    A = e8_coset_shell(list(V16[:8]), c)
    B = e8_coset_shell(list(V16[8:]), c)
    byn = defaultdict(list)
    for v, n in B:
        byn[n].append(v)
    out = []
    for v, n in A:
        for wv in byn.get(c - n, ()):
            out.append(v + wv)
    return out


def selftest():
    z = [F(0)] * 8
    th = {}
    for c in (F(0), F(2), F(4), F(6)):
        th[str(c)] = len([1 for _v, n in e8_coset_shell(z, c) if n == c])
    assert th == {"0": 1, "2": 240, "4": 2160, "6": 6720}, th
    z16 = [F(0)] * 16
    th2 = {str(c): len(shell(z16, c)) for c in (F(0), F(2), F(4))}
    assert th2 == {"0": 1, "2": 480, "4": 61920}, th2
    # coset invariant: a point of the coset lies in its own shell
    V = [F(2, 3), F(2, 3), F(-1, 3), F(-1, 3), F(-1, 3), F(0), F(0), F(0),
         F(1, 2), F(1, 2), F(0), F(0), F(0), F(0), F(0), F(0)]
    c = sum(x * x for x in V)
    assert tuple(V) in set(shell(V, c)), "coset representative not in its own shell"
    return dict(e8_theta=th, e8xe8_theta=th2, coset_invariant=True)


if __name__ == "__main__":
    print(selftest())

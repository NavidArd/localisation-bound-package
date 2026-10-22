# The class boundary, re-enumerated from the mathematics

The construction's negative result is only worth anything if the search behind it was
**exhaustive**, and it is exhaustive because exactly **two** inequivalent shift classes
admit the configuration at all. That "two" is the load-bearing number in the whole
construction: it is what turns a search into an exhaustion.

A number that load-bearing should not rest on one implementation. This directory holds a
second one.

---

## What is here

| file | |
|---|---|
| `paperA_two_classes_enum.py` | the enumeration, written from the mathematics |
| `paperA_two_classes_enum.json` | its output |
| `paperA_two_classes_enum.log` | its full transcript, stage by stage |
| `paperA_two_classes_bfs_certificate.json` | the breadth-first orbit certificate, one row per class (added 2026-10-22) |
| `check_bfs_certificate.py` | an independent checker that reads only the certificate (added 2026-10-22) |

```
cd two_classes && python paperA_two_classes_enum.py     # about 20 seconds
```

It needs Python and `numpy`. Its only input is the primary geometry file
`orbifolder-1.2.1/Geometry/Geometry_Z3xZ6_1_1.txt`, from which it re-derives the twist
vectors rather than taking them from anywhere else.

---

## Why it is independent, concretely

It imports **nothing** from `multi_geometry_campaign/code/su5_construction/` or
`gate_su5_construction/`. Built from scratch inside this one file: the `E8` membership
test, the 240 roots, the Standard-Model Cartan span, the ten norm-6/5 weights, the search
box over the eight admissible budgets, both branches of the root's location, the lattice
class key, the equivalence group, and the orbit reduction.

**Its invariant checks run before any count is believed** — and they are the reason it
found what it found:

* the `E8` theta series reproduced from its own shell enumerator: `1, 240, 2160, 6720` at
  norms `0, 2, 4, 6`;
* the twists parsed from the primary geometry file and asserted, with `delta_c = 5/36` and
  `2(1 − delta_c) = 31/18` recomputed from them;
* the class key verified **invariant** under 200 random lattice translations, **and**
  verified *not* to identify `x` with `x + e_1` — a key that passes only the first test is
  a key that merges classes, which is exactly the defect it went on to find;
* the count of first-block roots orthogonal to `W_SM` verified to be 20, twelve of them
  supported on coordinates 6, 7, 8.

---

## What it found

**The endpoint is unchanged. Two intermediate counts were wrong.**

| stage | this enumeration | as published (`g2_classes.py`) | |
|---|---:|---:|---|
| admissible budgets `B_v` | 8 | 8 | agrees |
| norm-6/5 weights enumerated | 10 | 10 | agrees |
| candidate roots `Delta` | 260 | 260 | agrees |
| branch 2, `B_v = −2` | 2,400 | 2,400 | agrees |
| **branch 2, `B_v = −5/2`** | **9,600** | **0 — never enumerated** | **differs** |
| **total solutions** (10 branch) | **13,200** | 3,000 | **differs** |
| **distinct lattice classes** (both branches) | **1,260** | 119 | **differs** |
| **inequivalent classes** | **2** | **2** | **agrees** |

Two independent defects, each of the kind an invariant check catches and a re-run of the
same code never would:

1. **A strict inequality where the derivation gives a weak one.** The published branch-2
   enumeration imposed `rem = 1/2` where the condition is `rem >= 1/2`, so the entire
   budget `B_v = −5/2` was skipped — 9,600 solutions and 480 classes never looked at. All
   of them lie in an orbit already found, which is why the answer did not move.
2. **A class key that quotients by too large a lattice.** The key divided out
   `<Z^8, (1/2)^8>`, which strictly contains `E8 (+) E8`, so distinct classes were merged:
   119 where the correct key gives 300 on that same subset.

**Reconciliation is explicit, not asserted.** Restricted to the published enumeration's own
subset, this code returns 6,000 solutions and 238 coarse-key classes against that
enumeration's 3,000 and 119 — so the discrepancy is located, not merely observed.

---

## Why "two" survives both defects

Three certificates, in the log:

* **Closure.** Every one of the 1,260 classes, under each of the 264 generators, maps to a
  class **inside** the enumerated set — 332,640 applications, zero escapes. Without closure
  a union–find computes connected components of a subgraph rather than true orbits.
* **Union–find** over those generators leaves two components.
* **An explicit word per class.** A breadth-first search from the two representatives
  reaches **all 1,260**, recording a parent and a connecting generator for each; composing
  the chain gives an automorphism carrying that class to its representative. None is
  orphaned, so no third class hides in the set.

  **The rows are stored** (added 2026-10-22). Until then the search was run and its
  completeness asserted, but no row was written, so the certificate could not be checked
  without re-running the enumeration. `paperA_two_classes_bfs_certificate.json` now holds
  all 1,260 rows -- class representative, parent row, connecting generator (the root of
  the reflection, or the global `-1`) -- seeded at the two representatives the paper
  prints, `V2^A` (row 0, 1,200 classes) and `V2^B` (row 1, 60 classes), maximum depth 4.
  `python check_bfs_certificate.py` reads only that file and verifies, in exact arithmetic,
  that every momentum has `|p|^2 = 31/18`, `|P_W p|^2 = 6/5` and order six, that every
  generator is a root reflection of the declared kind or the global `-1`, that each
  generator carries the parent's representative into the row's class modulo `E8 + E8`,
  that every chain of parents ends at row 0 or row 1, and that the 1,260 classes are
  pairwise distinct. Any single row can be checked by hand from its parent. Writing the
  rows changed nothing else: `paperA_two_classes_enum.json` is byte-identical to the
  deposited copy, and the log gains only the line that reports the certificate.

The two representatives, both of order 6 with `|p|² = 31/18`, `|w|² = 6/5`, `MI2 = 10`:

```
branch 2   V2    = (-2/3,-2/3, 1/3, 1/3, 1/3, 0, 0, 0 | -1/2,-1/2, 0,0,0,0,0,0)
           Delta = ( 0, 0, 0, 0, 0, 0, 0, 0          |   -1,  -1, 0,0,0,0,0,0)
                   Delta lies OUTSIDE the Standard-Model factor

branch 1   V2    = (-2/3,-2/3, 1/3, 1/3, 1/3, 0, 1/2, -1/2 | 0^8)
           Delta = ( 0, 0, 0, 0, 0, 0,  -1,    1      | 0^8)
                   Delta lies INSIDE the Standard-Model factor
```

An automorphism fixing the Standard Model cannot move a vector between the two `E8`
factors, so the two cannot merge — which is a stronger statement than "the generators used
did not merge them."

---

## Bookkeeping, so the numbers are not misread

The counts above are of the `10` branch alone unless marked. Adjoining the `10bar` branch
doubles every count exactly, with no class self-conjugate.

| | 10 branch | both branches | stored as |
|---|---:|---:|---|
| solutions `(p, Delta)` | 13,200 | **26,400** | `n_solutions_10_branch_only`, `n_solutions` |
| distinct `p` | 12,600 | **25,200** | `n_distinct_p` |
| lattice classes | 630 | **1,260** | `n_classes_exact` |
| inequivalent classes | 2 | 2 | `n_inequivalent` |

The JSON stores the **both-branch** figures. The repair memo that prompted this work quotes
the 10-branch halves, `12,600` and `630`; those are the same numbers in a different unit,
not a disagreement. The orbit count is two if the global `−1` is in the equivalence group
and four if it is not — the four being the two classes read once each way — and the JSON
records both (`n_inequivalent`, `n_inequivalent_without_minus1`). Either bookkeeping gives
**two**.

Claims `CON-01a`…`CON-03c`, `CON-01r`, `CON-02r`. The superseded rows `CON-01` and `CON-02`
are kept and labelled rather than deleted, so the correction can be audited.

---

## What this does not establish

The two-class exhaustion holds **given** that the Standard Model sits inside the block-1
`SU(5)` — hypotheses `S1`–`S4`, declared at the top of the script as scope, not derived
there. A Standard Model embedded otherwise, with the state still carrying a norm-6/5 role,
is **untested**, and is the loophole's last open direction. The paper must say so.

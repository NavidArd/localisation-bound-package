# The hidden-`SU(2)` signature — data, and the scope it holds at

A referee asked the sharpest question anyone has asked of this paper: **a root of the
gauge group left unbroken at a fixed point need not be a root of the gauge group left
unbroken in four dimensions.** If it is only the former, then two states differing by it
are not the two weights of a four-dimensional `SU(2)`, and the signature the paper reads
off a doubled norm-6/5 role is not there.

The referee is right in general, and the files in this directory establish both halves of
that: **where the lift fails, and where it does not.**

---

## 1. What is deposited here

| file | what it holds |
|---|---|
| `questionA_centralizer.json` | the lift index at **every** twisted constructing element of all **seven** campaign geometries — 533 of them — computed exactly by Smith normal form, with the escaping generators named wherever the index exceeds 1 |
| `questionB_Z2xZ4_1_1.json`, `questionB_Z3xZ6_1_1.json` | the catalogue-wide data test: every pair of distinct printed records sharing a Standard-Model role, both catalogues, and for each whether its `Δ` is a root of the **unbroken four-dimensional** group |
| `questionB.log` | the transcript of that run |
| `questionC_samefp.json` | conjugacy-class behaviour of the local shift `V_g`, and the count of distinct fixed points sharing one shift |
| `guard1d_projection.json` | the projection formula, implemented from orbifolder's own source, run against orbifolder's **printed** spectrum on the two constructed models |

Code: `../../code/root_lift/` (`rl_lib.py` and `rl_a…rl_d`). Regenerate with

```
cd multi_geometry_campaign/code/root_lift
python rl_a_centralizer.py && python rl_b_data.py && python rl_c_samefp.py && python rl_d_projection.py
```

The scripts read the same deposited raw census output and the same primary geometry files
as everything else in this package; they need Python and `numpy` only.

---

## 2. The scope the signature is deposited at

**The claim is scoped to the sectors `T(0,1)` and `T(0,5)` of `Z3xZ6_1_1`, and to those
sectors only.** At that scope it is derived, not asserted:

* each of those two sectors has a **single** constructing space-group element, and it
  carries **no lattice translation**, so the elements commuting with it include the whole
  point group; and
* this geometry's only Wilson line the lattice permits to be non-trivial, `A_3 ≐ A_4`,
  lies in **plane 2** — exactly the plane the sector's twist leaves fixed
  (`allowed order of Wilson lines = 1 1 3 3 1 1`), so `e_3, e_4` are pure translations
  *inside* the centraliser.

Together these give `Δ · V_1`, `Δ · V_2` and `Δ · A_i ∈ ℤ` for **all eight** generators.
The lift index is exactly **1**, and `Δ` is a root of the unbroken four-dimensional gauge
group. Both facts are read off the primary geometry file; neither is fitted.

### The negative half, which is the honest one

| geometry | twisted constructing elements | of which the lift holds |
|---|---:|---:|
| `Z3xZ6_1_1` | 100 | **5** |
| `Z6-I_1_1` | 42 | 1 |
| `Z2xZ4_1_1` | 74 | **0** |
| `Z4xZ4_1_1` | 123 | 0 |
| `Z3xZ3_1_1` | 108 | 0 |
| `Z6-II_1_1` | 44 | 0 |
| `Z4_1_1` | 42 | 0 |
| **total** | **533** | **6** |

Source: `questionA_centralizer.json`.

**On `Z2xZ4_1_1` — this paper's other catalogue — the lift fails at every one of the 74
twisted constructing elements**, because that geometry's Wilson lines occupy all three
planes. At its `T(0,1)` the escaping generators are `A_1 ≐ A_2` and `A_5 ≐ A_6`. No
statement of the signature's form may be made there, and the paper must not make one.

Of the five lifting elements on `Z3xZ6_1_1`, `T(0,1)` and `T(0,5)` are the two sectors
where the lifting element is the **only** element. `T(0,2)`, `T(0,3)` and `T(0,4)` carry
further constructing elements at which the lift fails, so the sector-level statement holds
for two sectors, not five.

### The data test agrees, both directions

Every twisted cell of both catalogues, keyed by the full `(k, l, n_1…n_6)` tuple:

| | `Z2xZ4_1_1` | `Z3xZ6_1_1` | total |
|---|---:|---:|---:|
| distinct-record same-role pairs | 24,406 | 7,888 | **32,294** |
| `Δ ∈ Λ`, `\|Δ\|² = 2`, `Δ ⊥ W_SM` | 24,406 | 7,888 | 32,294 |
| **`Δ` a root of the unbroken 4d group** | **0** | **0** | **0** |
| route G vs route I disagreements | 0 | 0 | **0** |

Two readings. First, the general failure is **realised, not hypothetical** — not one of
the 32,294 pairs printed as two records has a `Δ` that lifts. Second, the derivation's own
prediction holds exactly: at an index-1 element, two same-role states sharing right-mover
and oscillator *must* differ by a 4d root, so they must be printed as **one** record — and
such a pair occurs **0 times** in either catalogue.

At the flagship cell it is verified directly. `CB1534`, `T(0,1)`, fixed point
`(0,0,0,0,0,0)`: 20 candidate weights in the shell at `|p|² = 31/18`, **12 survive**, they
form **one** Weyl orbit of size 12, and orbifolder prints **one** record `(1, −3, 2, 1, 2)`
of dimension `3 × 2 × 2 = 12`. The partner is `Δ = (0⁸ | 1, 1, 0⁶)` with `Δ · V_1 = 0`,
`Δ · V_2 = 1`, `Δ · A_i = 0` on all six — integral on every generator.
Source: `guard1d_projection.json`.

---

## 3. What travels with the claim, and what is not licensed

**State the mechanism, never the conclusion,** when the claim is carried anywhere new. The
mechanism is two conditions: (i) a translation-free constructing element, and (ii) every
non-trivial Wilson line inside the torus the sector's twist leaves fixed. `Z6-I_1_1`'s
`T(3,0)` satisfies both — the same mechanism on a different geometry — which is the
transferable content.

**Not licensed:** any sentence of the form "a local root is a four-dimensional root"; any
use of the signature on `Z2xZ4_1_1`; any extension to constructing elements carrying a
non-zero translation, including those of `T(0,2)`, `T(0,3)` and `T(0,4)` of this same
geometry.

---

## 4. Known limitations, carried rather than buried

1. **The premise set is four items, not two.** The two states must share the right-moving
   *weight* **and** its excitation, and the left-moving oscillator excitation, and sit at
   one constructing element. An earlier statement of this ("same right-mover `q_sh`")
   under-specified it.
2. **The oscillator premise is load-bearing, and the data show it failing.** Nine pairs do
   sit at index-1 elements — all nine at the signature's own two sectors — and all nine
   carry **different oscillator excitations**, which is exactly the case where the phase
   does not cancel. For all nine, `Δ` is not a 4d root. The signature survives because no
   *shared*-oscillator index-1 pair exists in either catalogue. **That is a measurement
   over 16,162 models, not a theorem**, and on any new geometry with a small oscillator
   mode it must be re-checked before the signature is quoted.
3. **The boxed converse is proved one way.** "Index 1 ⟺ translation-free ∧ Wilson lines in
   the fixed plane" is established as a **sufficient** condition; the reverse implication
   would need an existence construction and is not supplied.
4. **Right-mover uniqueness is imported, not re-derived here** (31,053 records, no
   exception). If it fails on a new geometry the derivation stops.
5. **The projection formula was exercised whole on two models, not on a catalogue** — 112
   `(cell, q_sh)` groups reproduced exactly, including the load-bearing one.
6. **`Z6-I_1_1`'s lifting element is derived, not measured** — no `Z6-I` catalogue was
   consulted.
7. **Index > 1 means "not forced", not "never".** A particular model can pair `Δ`
   integrally with the escaping Wilson lines by accident; the hidden multiplicities 2, 3, 4
   and 5 seen at index > 1 cells are exactly that. Nothing here claims those impossible —
   only that no theorem covers them.

---

## 5. On `31,822` versus `32,294`

This pass counts **32,294** same-role record pairs where the main chain counts **31,822**.
The difference is definitional and in the safe direction. `z36sw_lib`'s role census assigns
`role = None` to records that are not hidden singlets — which is exactly the configuration
under review — so this phase re-derived the roles carrying the hidden multiplicity instead.
The pair set here is a **superset**, and the count of `Δ` that lift is zero on the larger
set.

`32,294` is the homogeneous number, applying one definition to both catalogues; `31,822`
sums two different role conventions across the two. **`32,294` is the number the paper
prints.** See `../../../DEFINITIONS.md` §2.11.

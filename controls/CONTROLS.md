# The three positive controls

A bound that only ever explains an absence explains nothing. These three controls show the
bound working in both directions — forbidding where it should forbid, permitting where it
should permit, and capping at exactly the predicted number — on the same roles, in the
same catalogues.

Every row below is a machine-checked claim in `../claims.json` (`CTRL-*`).

---

## Control (a) — a norm-4/5 role BELOW its threshold, at multiplicity 1

**Why this is a control and not a restatement.** These are the *same roles* — `d^c`, `L`
and their conjugates — that reach multiplicity **four** at higher budget in the same two
catalogues. So the asymmetry the paper reports tracks `2B − |w|^2` and **not** the role
label, which is the obvious objection to the whole construction.

The threshold for a norm-4/5 role is `B >= 13/20`, equivalently `R >= 1/2`.

| geometry | `R` | `B` | fields | (cell, role) pairs | multiplicity |
|---|---|---|---:|---:|---|
| `Z2xZ4_1_1` | 13/40 | 9/16 | 3,868 | 3,673 | **1 in every one** |
| `Z2xZ4_1_1` | 3/40 | 7/16 | 170 | 170 | **1 in every one** |
| `Z3xZ6_1_1` | 11/30 | 7/12 | 18 | 18 | **1 in every one** |
| `Z3xZ6_1_1` | 23/90 | 19/36 | 40 | 40 | **1 in every one** |
| `Z3xZ6_1_1` | 4/45 | 4/9 | 586 | 586 | **1 in every one** |
| **total** | | | **4,682** | **4,487** | **zero exceptions** |

The 195-field difference between 4,682 and 4,487 is the `Z2xZ4_1_1` fields at `R = 13/40`
that sit at a cell whose role also carries a `33/40` partner; they are counted once, in
the `33/40` bucket. See `../DEFINITIONS.md` §4.

**Headline instance, predicted before it was looked at.** `Z2xZ4_1_1`, sector `T(1,3)`,
budget `B = 7/16` (the sector's `B_0 = 11/16` lowered by an oscillator of mode `1/4`):

```
|w|^2 = 4/5      R = 2B - |w|^2 = 7/8 - 4/5 = 3/40 = 0.075  <  1/2
ceiling  n <= 1/(1 - 3/40) = 40/37 = 1.081   =>   n <= 1
measured 170 fields, 170 (cell, role) groups, multiplicity 1 in every one
```

Data: `../results_deposited/su5_ten_bound/verify_Z2xZ4_1_1.json` ->
`B_by_sector["(1, 3)"]["7/16"] = 170` and `R_by_fam["5bar"]["3/40"] = 170`;
`../results_deposited/gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json` ->
`multR_rec["5bar|3/40"] = {"1": 170}`.

Regenerate: `python verify.py Z2xZ4_1_1` and `python d1_mechanism.py Z2xZ4_1_1`.

---

## Control (b) — a norm-6/5 role ABOVE threshold, doubled: model CB1534

The bound permits a repeat only where `B >= 17/20`. Two `Z3xZ6_1_1` sectors, `T(0,1)` and
`T(0,5)`, satisfy that. **CB1534 is an explicitly constructed model in which the repeat
occurs**, so the permission is not vacuous.

CB1534 was built by the enumeration, censused with orbifolder, and then **re-run from its
embedding data alone by an independent deck**, reproducing 147 of 147 records identically
along with the gauge group, the `U(1)` generators, the simple roots and the anomaly
traces.

| quantity | value |
|---|---|
| field / sector / fixed point | `F_1` / `T(0,1)` / `(0,0,0,0,0,0)` |
| printed representation | `(1, 3bar, 2, 1, 2)` |
| oscillator | none |
| **hidden multiplicity** | **2** |
| role (triality-aware) | `Q` — a genuine `(3, 2)_{+1/6}` |
| `\|p_sh\|^2` | 31/18 |
| `B` | 31/36 |
| `\|w\|^2` | 6/5 |
| `R` | 47/90 |
| `p . Delta` | −1 |
| gauge group | `SU(4) x SU(3) x SU(2)` and `E_6 x SU(2)` and `U(1)^3` |
| `k_Y` | 5/3 |
| `Tr Y`, `Tr Y^3` | 0, 0 |
| net `Q` (triality-aware) | **+3** |

**The arithmetic that makes this a control rather than an anecdote:**

```
threshold   B = 31/36 = 0.8611...   >=  17/20 = 0.85          met, by 1/90
            equivalently in R:  R = 47/90  >=  1/2            met, by 1/45
bound       2*sqrt(47/90) = 1.4453  >=  sqrt(2) = 1.4142      SATURATED from the ALLOWED side, by 2.2%
ceiling     1/(1 - 47/90) = 90/43 = 2.093  =>  n <= 2         and exactly 2 is realised
```

> **Corrected in v2.0.** Earlier text gave the `B` margin as `1/45`. It is `1/90`:
> `31/36 − 17/20 = 155/180 − 153/180 = 2/180`. `1/45` is the margin in `R`, where the
> threshold is `1/2` and `R = 47/90`, and the two differ by the factor of 2 in
> `R = 2B − |w|²`. Both are stated above so neither is mistaken for the other. No other
> number moves.

The theorem is therefore **tight at this cell in both directions**: it forbids the
configuration everywhere below `17/20`, permits it here, and caps it at two.

Data: `../results_deposited/gate_su5_construction/g1c_compare.json` ->
`dossiers.CB1534`, `compare.CB1534`. Input deck: `../construction/CB1534/`.

### CB1041 — DEMOTED

An earlier flagship, **CB1041**, is deposited alongside CB1534 and marked demoted. Its
doubled cells reproduce, and its record-level comparison is clean (214 / 214 identical
fresh vs stored). It was demoted because a **triality-blind role table** counted exotics
as quark doublets: `s8_battery.py::su3_dim` returns a `(dimension, sign)` pair and every
caller discards the sign, so `(3bar, 2)_{+1/6}` and `(3, 2)_{−1/6}` were both labelled
`Q`. Under the triality-aware recount CB1041 has `|net Q| = 4`, not 3, so it is not a
three-generation model. **Paper A should print CB1534 only.** See
`../construction/CB1041_DEMOTED/README.md`.

---

## Control (c) — norm-4/5 roles ABOVE threshold, doubled, AT the ceiling

Four independent saturations of `n <= floor(1/(1−R))`, at four different values of `R`,
on both geometries:

| geometry | `R` | `B` | ceiling | measured multiplicity histogram | status |
|---|---|---|---:|---|---|
| `Z2xZ4_1_1` | 23/40 | 11/16 | **2** | `{1: 41,129, 2: 88}` | **saturated** |
| `Z3xZ6_1_1` | 53/90 | 25/36 | **2** | `{1: 833, 2: 206}` | **saturated** |
| `Z2xZ4_1_1` | 7/10 | 3/4 | **3** | `{1: 130,464, 2: 20,188, 3: 883}` | **saturated** |
| `Z3xZ6_1_1` | 34/45 | 7/9 | **4** | `{1: 91,885, 2: 6,509, 3: 204, 4: 4}` | **saturated** |

In each row the measured maximum equals the ceiling exactly, and in no row is it
exceeded.

**Two non-saturated rows, printed alongside so the ceiling is not mistaken for a fit:**

| geometry | `R` | ceiling | measured max | status |
|---|---|---:|---:|---|
| `Z2xZ4_1_1` | 33/40 | 5 | 4 | one short |
| `Z3xZ6_1_1` | 83/90 | 12 | 2 | slack |

The complete ceiling-versus-measurement table over all realised `R`, both catalogues, is
in `claims.json` under the `GRAM-*` ids, and reproduced in `../claims.md`.

Data: `../results_deposited/gate_su5_bound/d1_mechanism_<geom>.json::multR_rec`.

---

## And the mechanism, measured

Not a control, but the strongest positive evidence in the package, and it belongs beside
them. All **31,822** repeated same-role pairs across both catalogues differ by exactly a
**root** of `E8 x E8` orthogonal to the Standard-Model Cartan span.

> **On the denominator.** 31,822 is the set this script tests; it sums two different role
> conventions across the two geometries, and the homogeneous number is **32,294** (see
> `../DEFINITIONS.md` §2.11). The larger set is a **superset**, and it was independently
> tested in full by `root_lift/rl_b_data.py`: on all 32,294, `Delta` lies in the lattice,
> has norm 2, and is orthogonal to the Standard-Model Cartan span. So the result below is
> not an artifact of the smaller denominator — it holds on the larger one too.

| test | `Z2xZ4_1_1` | `Z3xZ6_1_1` |
|---|---:|---:|
| pairs tested | 24,406 | 7,416 |
| `Delta` not in `Lambda` | 0 | 0 |
| `Delta = 0` | 0 | 0 |
| `P_W Delta != 0` | 0 | 0 |
| `\|Delta\|^2` odd or `< 2` | 0 | 0 |
| Cauchy–Schwarz failure | 0 | 0 |
| `\|Delta\|^2 = 2` and `Delta` is one of the 480 roots | **24,406 / 24,406** | **7,416 / 7,416** |
| Gram ceiling failures | 0 | 0 |

That is the corollary observed directly on 31,822 instances, not inferred.

Data: `../results_deposited/gate_su5_bound/d1_mechanism_<geom>.json`.
Regenerate: `python d1_mechanism.py Z2xZ4_1_1` and `python d1_mechanism.py Z3xZ6_1_1`.

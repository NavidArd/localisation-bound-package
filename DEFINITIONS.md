# Definitions — the counted objects, and which denominator each number uses

> **Dated note, 2026-09-10 — the paper's counts have moved to one uniform convention, and the
> tables in this file have not yet been rewritten to match.** The paper now counts under a
> single role convention (a record carries a role by its colour, isospin and hypercharge
> assignment alone, with no condition on the remaining gauge factors) and a single embedding
> rule (each catalogue model once, at the first printed embedding whose Standard-Model
> identification closes), with the banner parse repaired so that a model printing an `E_n`
> factor is read correctly. **The units this file defines are unchanged; the values attached to
> them are not.** For the values the paper prints, read the rows of `claims.json` whose
> `paper_section` is `verification (uniform convention)`, and their sources under
> `results_deposited/uniform_census/`. In particular the headline ratio is now **2.414**
> (`109,717 / 45,450`) rather than `2.413`, the same quantity in the same matched units; the
> same-role record-pair count is **33,743** rather than `32,294`, one convention on both
> geometries in both cases; and the census covers **16,138** distinct models rather than
> `15,535`. Every movement is upward except two per-geometry rows on `Z3xZ6_1_1`, where ten
> models the census output contains twice are now counted once. **No verdict moves:** zero
> violations of the bound in both units, the same 65 unprotected norm-6/5 fields, every
> norm-6/5 multiplicity and hidden dimension still one.
>
> **What below is stale, in full.** Not only §§2.3, 2.7, 2.9 and 2.11, as an earlier
> wording of this note said, but also **the preamble immediately below it, the §3
> quick-reference table and §4**. In particular the preamble's bolded sentence *"If you
> read only one thing here: the ratio to print is 2.413"* is the previous convention's
> value of that ratio and is **superseded**: the ratio to print is **2.414**
> (`109,717 / 45,450`), the same quantity in the same matched units, and `claims.json`
> row `UNI-30` carries it. The units this file defines are unchanged throughout; the
> values attached to them, wherever they appear below, are the previous convention's.
> Rewriting them is an open item.

Several numbers in the paper count different things over the same data, and the ratios
between them are large (up to a factor of 2.5). A reader who takes two of them to be the
same quantity will conclude that the paper contradicts itself. This file fixes each unit
once, says which numbers use it, and gives the file and key where each is stored.

Nothing here is new physics. It is a glossary, and it exists because the units were the
single most common source of confusion in earlier drafts.

**If you read only one thing here:** the ratio to print is **2.413** — norm-6/5
`(cell, role)` pairs per model-cell carrying any norm-6/5 role, `104,889 / 43,466`, over
both catalogues in matched units. `2.48` is the same ratio on one geometry, and
`104,889 / 42,932 = 2.443` is unit-mismatched and must not be printed at all. §4 lists
every ratio in circulation and what each is right about.

---

## 1. The chain of objects

```
catalogue file           one text file per orbifold geometry
   -> model              one entry in that file: a gauge embedding
      -> embedding       one vev-configuration printed for that model
         -> record       one printed left-chiral field entry
            -> field     one physical left-chiral multiplet
               -> state  one field counted with its hidden multiplicity
```

and, cutting across the last three:

```
sector (k,l)  +  full fixed-point tuple n_1..n_6   ==  cell   (a constructing element)
cell  +  one model                                 ==  model-cell
model-cell  +  role                                ==  (fixed point, role) pair
```

**Only 128 cells exist** across both geometries — 58 and 70. Every catalogue-wide total of
order 10⁵ or more counts **model-cells**, fields, or pairs. None of them counts fixed
points, and the single most common misreading of this work is to think one does.

---

## 2. The units, one at a time

### 2.1 `model`

One `begin model ... end model` block of the catalogue file: a label, a lattice, two
shift vectors and six Wilson lines. This is the unit in which the catalogues are sized.

| | `Z2xZ4_1_1` | `Z3xZ6_1_1` | total |
|---|---:|---:|---:|
| models in the catalogue file | **11,540** | **4,612** | 16,152 |
| model blocks in the raw census output | 11,540 | 4,622 | 16,162 |
| **census runs** | 10,925 | 4,620 | **15,545** |
| **distinct models censused** | 10,925 | **4,610** | **15,535** |
| excluded, with reason | 615 | 2 | 617 |

**Print 15,535, not 15,545.** "We scanned N models" is read as N *distinct* models, and
ten `Z3xZ6_1_1` models are censused twice (§2.2). The as-run figure 15,545 is a count of
**census runs** and is what the deposited tables store; both are given here so neither is
mistaken for the other. The de-duplication moves no conclusion: the ten duplicated models
contribute no unprotected `R = 47/90` field, no below-threshold pair, and no exception of
any kind, and the headline ratio is unmoved (2.4129 de-duplicated against 2.4131 as run).

Source: `results_deposited/gate_su5_bound/d5_excluded.json` -> `<geom>.loaded`,
`.ok`, `.skip.closure_failed`, `.n_distinct_labels`, `.duplicate_labels`.
Claim rows `A-01` (census runs) and `DEN-07` (distinct).

### 2.2 `record` (of a model) — and where the duplication actually comes from

One model block as it appears in the file. On `Z3xZ6_1_1` **ten labels appear twice in the
census output** (`MSSM0`–`MSSM9`), so 4,622 blocks were parsed under 4,612 distinct labels.

> **The duplication is NOT a defect of the published catalogue.** The published
> `Z3xZ6_1_1` file has 4,612 blocks, all labels distinct and all contents distinct. The
> repetition is a **pilot-segment overlap in this program's own census run**:
> `census_segPILOT.txt.gz` re-runs `MSSM0…MSSM9`, which also appear in
> `census_seg0.txt.gz`. Earlier versions of this file and of `EXCLUSIONS.md` attributed it
> upstream; that attribution was wrong and is corrected here.
>
> **Direction:** the effect is over-counting, so every absence claim in the paper is safe.
> `Z2xZ4_1_1` is clean — segments `census_seg0 … census_seg57`, no pilot,
> `duplicate_labels = 0`.

Both copies are censused, and the duplication is reported rather than silently
de-duplicated. See `exclusions/EXCLUSIONS.md` §2.

### 2.3 `embedding` (vev-configuration)

orbifolder can print more than one vev-configuration per model. **The two catalogues use
opposite selection conventions, and both are stated here:**

| | which printed hypercharge block is read |
|---|---|
| `Z2xZ4_1_1` | the **first** (`z24sw_lib` and the census take `cfgs[0]`) |
| `Z3xZ6_1_1` | the **last** (`z36sw_lib.parse_gauge` takes `blocks[-1]`) |

53 censused `Z3xZ6_1_1` models print more than one block and are affected in principle.
**Measured impact: none** — on all 55 two-block models, both choices give identical role
tallies. Earlier versions of this file asserted "the first" for everything; that was an
assertion rather than a measurement, and it was wrong for one of the two catalogues.

This matters because a model whose first configuration does not close can still close
under a later one. On `Z2xZ4_1_1`, **599 of the 615 excluded models do exactly that**;
they were recovered and re-censused separately, adding 4,910 further (fixed point,
`10`-role) pairs, **all at multiplicity 1**. The paper's absence claim therefore covers
11,524 of 11,540 `Z2xZ4_1_1` models, not 10,925. See `exclusions/EXCLUSIONS.md` §3.

### 2.4 `record` (of a field) — the printed spectrum entry

One `F_n` block in the orbifolder spectrum print-out: sector, fixed point, representation,
`U(1)` charges, left-moving `p_sh`, right-moving `q_sh`, oscillators, gamma phases.

**Only left-chiral (`_l`) records are counted.** The CPT partner is never counted
separately.

| | `Z2xZ4_1_1` | `Z3xZ6_1_1` | total |
|---|---:|---:|---:|
| left-chiral records **at twisted cells** | 2,106,082 | 1,006,547 | **3,112,629** |
| the same, **including the untwisted sector** | 2,247,793 | 1,050,282 | 3,298,075 |

**3,112,629 is the twisted count.** The untwisted sector `(0,0)` is dropped before this
number is taken, and calling it "left-chiral records" without qualification loses 185,446
of them. Everything downstream in this package is twisted-only, because the theorem is a
statement about twisted fixed points.

Source: `results_deposited/gate_su5_bound/d2_pshdup_<geom>.json`.

### 2.5 `field`

A record that carries a Standard-Model role (§2.9). This is the unit of the two large
aggregates.

| | `Z2xZ4_1_1` | `Z3xZ6_1_1` | total |
|---|---:|---:|---:|
| role-carrying fields | 312,050 | 169,678 | **481,728** |
| of which norm-6/5 (`10`-family) | 67,555 | 37,334 | **104,889** |
| of which norm-4/5 (`5bar`-family) | 244,495 | 132,344 | **376,839** |

Source: `results_deposited/su5_ten_bound/verify_<geom>.json` -> `fields`,
`w2_by_fam`.

> **104,889 is a FIELD count.** It is not a count of fixed points and not a count of
> models. It happens also to equal the number of (fixed point, `10`-role) pairs (§2.8),
> because that multiplicity is 1 everywhere — which is the result. The two units
> coincide *here* only because the theorem holds.

### 2.6 `state`

A field weighted by its **hidden multiplicity** — the product of its dimensions under the
non-Standard-Model gauge factors. A field printed once as a doublet of a hidden `SU(2)`
is one record but two states.

For the `10`-family the two units coincide: hidden multiplicity is **1 on all 104,889**
norm-6/5 fields, with no exception. For the `5bar`-family they differ materially:

| per-(cell, role) multiplicity | 1 | 2 | 3 | 4 |
|---|---:|---:|---:|---:|
| `Z2xZ4_1_1`, counting **records** | 198,429 | 21,664 | 910 | 2 |
| `Z2xZ4_1_1`, counting **states** | 189,003 | 29,723 | 2,198 | 81 |

Both are correct; they answer different questions. Quote which one is meant.
Source: `d1_mechanism_<geom>.json` -> `mult_rec_by_fam`, `mult_hid_by_fam`,
`hid_by_fam`.

### 2.7 `fixed point` / `cell`, and `constructing element`

A **cell** is the pair (sector `(k,l)`, full fixed-point tuple `n_a`), `a = 1..6` — the
key orbifolder prints. The full six-component tuple is used; truncating it merges
distinct fixed points and inflates multiplicity.

| | `Z2xZ4_1_1` | `Z3xZ6_1_1` | total |
|---|---:|---:|---:|
| **model-cells**, twisted | 554,922 | 296,309 | **851,231** |
| distinct twisted cells = **constructing elements** | 58 | 70 | **128** |
| model-cells carrying an **unbarred** norm-6/5 role | 26,796 | 16,136 | **42,932** |
| model-cells carrying a norm-6/5 role of **either chirality** | 27,216 | 16,250 | **43,466** |
| model-cells carrying a norm-4/5 role of either chirality | 122,977 | 68,532 | **191,509** |

> **851,231 counts MODEL-CELLS, not fixed points.** Across both catalogues only **128**
> twisted cells occur at all — 58 and 70, one per populated constructing element, matching
> an independent affine orbit count sector by sector. 851,231 is those 128 counted once
> per model. **Any total of order 10⁵ or larger in this paper is a count of model-cells,
> of fields, or of pairs — never of fixed points.**

> **Two distinct-cell counts exist and they are not interchangeable.** The pattern
> census counts model-cells carrying `Q`, `u^c` or `e^c` — 26,796 and 16,136, summing to
> **42,932**. Counting model-cells that carry *any* norm-6/5-family role, conjugates
> included, gives 27,216 and 16,250, summing to **43,466**. The difference is not simply
> the conjugate cells: 495 `Z2xZ4_1_1` cells carry a conjugate role, but 75 of them also
> carry an unbarred one, so `26,796 + 495 − 75 = 27,216`; on `Z3xZ6_1_1`,
> `16,136 + 164 − 50 = 16,250`.
>
> **Which one to use depends on the numerator, and the two are not swappable:**
>
> * **43,466** is the matched-unit denominator of **104,889**, because 104,889 counts all
>   six norm-6/5 roles, conjugates included. The ratio is **2.413**. *(claim `DEN-01`)*
> * **42,932** may be printed **only** as *"model-cells carrying at least one of `Q`,
>   `u^c`, `e^c`"* — the denominator the occupancy-pattern statements are made over. It is
>   **not** the denominator of 104,889; `104,889 / 42,932 = 2.443` is unit-mismatched,
>   numerator counting conjugate roles where the denominator does not, and must not be
>   printed. *(claim `A-23`)*
>
> Never write *"cells carrying a `10`-role"* — the phrase is ambiguous between the two.
> Source: `patterns_<geom>.json::cells10` and
> `d5_independence.json::task2_what_is_67555.distinct_cells_measured_by_this_gate`.

A **constructing element** is the space-group element
`g = (theta^k omega^l, sum_a n_a e_a)` — the point-group element *and* the translation.
It is not the same object as a geometric fixed point: one point of the covering torus can
be fixed by several inequivalent elements, and the local shift
`V_g = k V_1 + l V_2 + sum_a n_a W_a` depends on the translation through the Wilson lines.

**The theorem's hypothesis H3 is stated on the constructing element, not on the point.**
The printed cell key is what the census uses; the two agree here because the key was
certified bijective against an independent affine point-group orbit count — in every
populated sector, six on `Z2xZ4_1_1` and thirteen on `Z3xZ6_1_1`, the printed
fixed-point count equals the simultaneous orbit count exactly. So the key cannot split
one physical fixed point into two and hide a doubling. The sectors printing zero records
are precisely the CPT partners of populated ones.
Source: `results_deposited/gate_su5_bound/d5_cellkey.json`, `d5_keysplit.json`,
`d5_keycount.json`.

**The conclusion does not depend on the key convention.** The independent recount
(`d6_pattern_recount.py`) counts the occupancy patterns under three different keys, and
the two derived-impossible patterns stay at zero in all three:

| key convention | `111` | `011` | `100` | `110` | `101` | `010` | `001` |
|---|---:|---:|---:|---:|---:|---:|---:|
| full `(sector, n_a)` — the one used | 17,555 | 4,820 | 4,421 | **0** | **0** | **0** | **0** |
| space-group orbit key (coarser) | 17,543 | 4,798 | 4,389 | **0** | **0** | **0** | **0** |
| sector only (much coarser) | 8,571 | 782 | 673 | **0** | **0** | **0** | **0** |

Coarsening a key can only *merge* cells, which can only *raise* a multiplicity — the
unsafe direction for an absence claim. That the forbidden patterns stay empty even under
the coarsest key is therefore the strong direction of the test.
Source: `results_deposited/gate_su5_bound/d6_patterns_Z2xZ4_1_1.json`.

### 2.8 `(fixed point, role) pair`

One cell together with one role label. **This is the unit of the multiplicity
measurement** and the unit in which the ceiling is stated.

| | `Z2xZ4_1_1` | `Z3xZ6_1_1` | total |
|---|---:|---:|---:|
| (cell, norm-6/5 role) pairs | 67,555 | 37,334 | **104,889** |
| (cell, norm-4/5 role) pairs | 221,005 | 125,146 | **346,151** |

**346,151 is a real total and the record never printed it.** Placing the one-geometry
221,005 beside the both-catalogue 104,889 compares one geometry with two; print either
both per geometry or both as totals. *(claims `DEN-04`, `DEN-05`)*

The arithmetic that ties the pair count to the cell count on `Z2xZ4_1_1`:

```
17,555 cells x 3 roles + 4,421 x 1 + 4,820 x 2   = 66,726     (patterns 111, 100, 011)
    87 cells x 3       +   248 x 1 +   160 x 2   =    829     (conjugate roles)
                                                  -------
                                                    67,555
```

so 26,796 distinct cells carry 67,555 (cell, role) pairs — a ratio of 2.52. Both
numbers are right; they are not the same number.
Source: `results_deposited/su5_ten_bound/patterns_Z2xZ4_1_1.json` -> `pat`, `patbar`,
`cells10`.

> **26,134 of the 42,932 cells carry three `10`-role fields.** So "every fixed point
> carrying a `10` field carries exactly one" is FALSE as English. The true statement is
> per role: every one of the 104,889 (cell, role) pairs has multiplicity one.

### 2.9 `role`

The label a state's `(SU(3)_c, SU(2)_L, Y)` assignment gives it: `Q`, `u^c`, `e^c`,
`d^c`, `L`, `Lbar`, and their conjugates. Roles are taken verbatim from the certified
`role_census`; this package does not re-open the Standard-Model identification.

**Two norm classes, five classes up to conjugation:**

| norm class | roles | `\|w\|^2` |
|---|---|---|
| norm-6/5 (the `SU(5)` **10**) | `Q`, `u^c`, `e^c` + conjugates | 6/5 |
| norm-4/5 (the `5`/`5bar`) | `d^c`, `L` + conjugates | 4/5 |

`(1, 2)_{+1/2}` — which the census labels `Hu` and the construction parser labels `Lbar`
— is the **conjugate of `L`**, and lies in the **5**, not the `5bar`. `\|w\|^2` is
invariant under conjugation, so no number moves; the label matters because it is the kind
of thing a referee reads as a sign that conjugation is not being tracked.

### 2.11 `same-role record pair` — and why the paper prints 32,294, not 31,822

A pair of **distinct printed left-chiral records** at one model-cell carrying the **same**
Standard-Model role. This is the unit the proof is tested on.

| convention | `Z2xZ4_1_1` | `Z3xZ6_1_1` | total |
|---|---:|---:|---:|
| library role convention (`d1_mechanism.py`) | 24,406 | 7,416 | **31,822** |
| one uniform convention, both catalogues (`root_lift/rl_b_data.py`) | 24,406 | **7,888** | **32,294** |

The two differ on `Z3xZ6_1_1` only, and the reason is definitional: `z36sw_lib`'s role
census assigns `role = None` to records that are not hidden singlets, while `z24sw_lib`
does not. **So 31,822 sums two different role conventions across the two geometries**;
32,294 applies one definition to both.

> **Print 32,294**, defined as *"pairs of distinct printed left-chiral records at one
> model-cell carrying the same role."* The larger set is a **superset** of the smaller, and
> every conclusion drawn on it — `|Δ|² = 2` in every case, and zero pairs whose `Δ` lifts
> to the four-dimensional group — holds on the superset, which is the safe direction.
> *(claim `DEN-06`; the superseded row is `A-18`)*

A consequence worth naming: the `Z3xZ6_1_1` "state = record" equalities in `d5_states.json`
are artifacts of the same gate and **must not** be quoted as a hidden-multiplicity
measurement on that geometry. `104,889` is unaffected by the convention; `481,728` and
`376,839` are not.

### 2.10 `budget B`, `R`, `\|w\|^2`

* `B = \|p_sh\|^2 / 2 = 1 - dc(g) - N_osc` — the left-moving massless budget.
* `\|w\|^2 = \|P_W p\|^2` where `W` is the 4-dimensional Standard-Model Cartan span.
  Obtained by **explicit projection** on all 481,728 fields, never inferred from `k_Y`.
* `R = 2B - \|w\|^2` — the squared norm orthogonal to `W`.

`B` is computed from the **printed** `p_sh`, not from `dc`. `dc` enters only the
invariant check `I1`, which passed 481,728 / 481,728 exactly.

---

## 3. Quick reference — which unit each headline number uses

| number | unit | claim | where |
|---:|---|---|---|
| **15,535** | **distinct models** censused | `DEN-07` | `d5_excluded.json::<geom>.ok` less `.duplicate_labels` |
| 15,545 | **census runs** (10,925 + 4,620) | `A-01` | `verify_<geom>.json::models` |
| **851,231** | **model-cells**, twisted — 128 cells once per model | `A-02` | `verify_<geom>.json::cells` |
| 128 | twisted cells = **constructing elements** (58 + 70) | — | `d5_cellkey.json` |
| 3,112,629 | left-chiral **records at twisted cells** | `A-04` | `d2_pshdup_<geom>.json` |
| 3,298,075 | the same **including the untwisted sector** | — | `denominators` table |
| 481,728 | role-carrying **fields** | `A-03` | `verify_<geom>.json::fields` |
| **104,889** | norm-6/5 **fields** = (cell, role) **pairs** | `A-05`, `DEN-05` | `verify_<geom>.json::w2_by_fam`, `::mult_by_fam.10` |
| 376,839 | norm-4/5 **fields** | `A-06` | `verify_<geom>.json::w2_by_fam` |
| **346,151** | (cell, norm-4/5 role) **pairs** | `DEN-04` | `verify_<geom>.json::mult_by_fam.5bar` |
| **43,466** | **model-cells** with a norm-6/5 role of either chirality — the matched-unit denominator of 104,889 | `DEN-01` | `d5_independence.json::task2…distinct_cells…` |
| 42,932 | **model-cells** carrying `Q`, `u^c` or `e^c` — the occupancy denominator, and **not** the denominator of 104,889 | `A-23` | `patterns_<geom>.json::cells10` |
| 191,509 | model-cells with a norm-4/5 role of either chirality | `DEN-02` | `d5_independence.json::task2…` |
| **2.413** | the ratio to print: 104,889 / 43,466, matched units, both catalogues | — | derived |
| 2.482 | the same on `Z2xZ4_1_1` **alone** — this is the "2.48" the record quotes | `DEN-03` | derived |
| 26,134 | model-cells carrying **three** norm-6/5 fields | `A-24` | `patterns_<geom>.json::pat["111"]` |
| **32,294** | same-role record **pairs**, one convention both catalogues | `DEN-06` | `questionB_<geom>.json::pairs` |
| 31,822 | the same under **two** conventions — superseded | `A-18` | `d1_mechanism_<geom>.json::pairs` |
| **4,487** | below-threshold (cell, norm-4/5-role) pairs | `CTRL-a-total` | `d1_mechanism_<geom>.json::multR_rec` |
| 4,682 | the **fields** in those pairs | — | `verify_<geom>.json::R_by_fam` |
| **21,981** | `(cell, p_sh)` **GROUPS** of ≥2 records — *not* a count of cells | `A-25` | `d2_pshdup_<geom>.json::dup_groups` |
| 4,032 | SM-valid constructed **models** (2,912 + 1,120) | `CON-14` | `construction/` |
| **2,174** | of those, models carrying a **doubled** norm-6/5 role | `CON-24` | `g8c_doubled_smvalid.json` |
| 6 of 533 | constructing elements where the root lift is **forced** | `SIG-01` | `questionA_summary.json` |

---

## 4. The ratios most likely to be misread

**The one to print is 2.413** — norm-6/5 (cell, role) pairs per model-cell carrying any
norm-6/5 role, `104,889 / 43,466`, both catalogues, matched units. Everything below is a
different ratio over the same data, and each is right about something else.

1. **2.482** — the same ratio on `Z2xZ4_1_1` **alone** (67,555 / 27,216). This is the
   "2.48" that circulated as a both-catalogue figure; it is one geometry. On `Z3xZ6_1_1`
   alone it is 2.298.
2. **2.443** — `104,889 / 42,932`. **Unit-mismatched and must not be printed:** the
   numerator counts conjugate roles, the denominator does not.
3. **2.52** — (cell, role) pairs to *unbarred*-role cells on `Z2xZ4_1_1` (67,555 vs
   26,796). Quoting any pair count with the words "fixed points" is the error this file
   exists to prevent.
4. **1.807** — the norm-4/5 counterpart, `346,151 / 191,509` (`Z2xZ4` 1.797,
   `Z3xZ6` 1.826).
5. **1.012** — between the two distinct-cell conventions themselves (43,466 vs 42,932).
   Small, which is why the confusion survives; see §2.7.
6. **1.09** — the `4,682` fields to the `4,487` pairs below threshold. The gap is the
   195 `Z2xZ4_1_1` fields at `R = 13/40` sitting at a cell whose role also carries a
   `33/40` partner; those are counted once, in the `33/40` bucket.
7. **1.15** — record-counted vs state-counted norm-4/5 multiplicity-1 buckets on
   `Z2xZ4_1_1` (198,429 vs 189,003). Hidden multiplicity is the difference.

---

## 5. One value that names two unrelated things

**`4,910`.** It is the number of additional norm-6/5 (cell, role) pairs recovered from the
599 models of §2.3, **and** the number of generated `Q·u^c·H_u` trilinears in an unrelated
over-kill test elsewhere in the campaign. They are not the same object. Only the first
belongs to this paper; if both ever appear, they must be distinguished in words at each
occurrence.

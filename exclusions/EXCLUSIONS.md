# Excluded, failed, duplicate and recovered cases

Every model in both catalogues is accounted for here. Nothing was silently dropped; the
accounting closes to the declared catalogue sizes exactly, and that closure is itself a
machine check (`d5_excluded.json::<geom>.accounting_ok = true`).

Data: `results_deposited/gate_su5_bound/d5_excluded.json`,
`d5_independence.json` (task 5), `d5_states.json`.
Scripts: `multi_geometry_campaign/code/gate_su5_bound/d5_excluded.py`, `d5_states.py`.

---

## 1. The accounting, both catalogues

| | `Z2xZ4_1_1` | `Z3xZ6_1_1` |
|---|---:|---:|
| orbifolder-declared orbifolds loaded | 11,540 | 4,622 |
| independent regex block count | 11,540 | 4,622 |
| library `split_models` block count | 11,540 | 4,622 |
| **censused (used for every number)** | **10,925** | **4,620** |
| excluded: first configuration does not close | 615 | 2 |
| excluded: no SMConfig printed | 0 | 0 |
| excluded: no hypercharge generator (**DATA-GAP**) | **0** | **0** |
| excluded: parse error | 0 | 0 |
| excluded: SM basis failure | 0 | 0 |
| **silently dropped** | **0** | **0** |
| censused + excluded = loaded | 11,540 ✓ | 4,622 ✓ |

> **Neither Paper A catalogue contains a hypercharge DATA-GAP model.** The 12
> hypercharge DATA-GAP models that appear elsewhere in this research programme belong to
> the `Z3xZ3` 77-model catalogue and a different paper. They are not in scope here, and
> no number in Paper A is conditional on them.

---

## 2. The ten duplicate records (`Z3xZ6_1_1`) — a pilot-segment overlap, not a catalogue defect

The `Z3xZ6_1_1` **census output** contains **4,622 model blocks under 4,612 distinct
labels**: ten labels appear twice (`MSSM0`–`MSSM9`).

| | |
|---|---|
| blocks in the published catalogue file | 4,612, all labels distinct, all contents distinct |
| model blocks parsed from the census output | 4,622 |
| duplicate labels | **10** |
| examples | `d5_excluded.json::Z3xZ6_1_1.duplicate_label_examples` (all ten listed) |

> **Where the duplication actually comes from.** Earlier versions of this file, and of
> `DEFINITIONS.md` §2.2, attributed it to the published catalogue. **That was wrong.** The
> published file has 4,612 blocks with distinct labels *and* distinct content. The
> repetition is a **pilot segment of this program's own census run**:
> `census_segPILOT.txt.gz` re-runs `MSSM0…MSSM9`, which also appear in
> `census_seg0.txt.gz`. `Z2xZ4_1_1` has no pilot segment (`census_seg0 … census_seg57`)
> and no duplicates.

Both copies are parsed and censused, and the overlap is reported rather than
de-duplicated, so that the as-run figures in the deposited tables are what the scripts
actually produced. The distinct-model denominators are given in `DEFINITIONS.md` §2.1:
**15,535 distinct models over 15,545 census runs**, and that is what the paper prints.

**Direction of the exposure.** Counting a model twice can only *inflate* a multiplicity
census. The paper's `10`-role claim is an **absence** claim — that no (cell, role) pair
has multiplicity above one — and an over-count cannot manufacture an absence. The
duplication is therefore safe in the direction that matters, and this is stated rather
than assumed.

`Z2xZ4_1_1` has no duplicate labels (`duplicate_labels = 0`).

---

## 3. The 599 recovered models (`Z2xZ4_1_1`)

**The problem.** orbifolder can print more than one vev-configuration per model. Every
number in the paper is read from the **first** one (`cfgs[0]`). A model whose first
configuration does not close is excluded — but it may close under a later one.

| | |
|---|---:|
| models excluded as closure-failed | 615 |
| **of those, closing under a later printed configuration** | **599** |
| genuinely unrecoverable (only one configuration printed) | **16** |
| `k_Y` of the excluded set, read off every printed generator | `5/3` in **615 / 615** |

The 599 were recovered and censused separately:

| | |
|---|---:|
| additional (cell, `10`-role) pairs | **4,910** |
| their multiplicity, counting **records** | 1 in all 4,910 |
| their multiplicity, counting **states** | 1 in all 4,910 |
| distinct cells carrying a `10`-role | 2,013 |
| additional (cell, `4/5`-role) pairs | 8,755 at mult 1, 1,463 at mult 2 |
| **hits (any `10`-role doubling)** | **0** |

**Consequence for the paper.** The absence claim covers **11,524 of 11,540**
`Z2xZ4_1_1` models, not 10,925. If the paper quotes 10,925 as the scope, it understates
its own coverage; if it quotes 11,540, it overstates it. 11,524 is the correct number,
and the 16 residual models all print `k_Y = 5/3`, so none of them can carry a different
threshold.

> **REPORTED MISMATCH (not corrected here).** The gate document
> `gate_su5_ten_bound.md` §6.5 and defect register G-10 say **"18 genuinely
> unrecoverable"**. The arithmetic is `615 − 599 = 16`, and the on-disk breakdown agrees:
> `d5_excluded.json::Z2xZ4_1_1.excl_ncfg = {"2": 584, "1": 16, "4": 15}` — the 16 models
> printing a single configuration are exactly the unrecoverable ones. **The number in the
> prose is 18; the number on disk is 16.** An independent recount from the raw output
> confirms 16 directly: 10,941 models print a single embedding = 10,925 censused + 16
> excluded. Both the claim "11,524 of 11,540" and the
> claim "0 hits" are unaffected. This package reports the discrepancy and does not
> resolve it; the source documents were not edited.

`Z3xZ6_1_1`: `excl_any_cfg_closes = 0`. Its **2** excluded models print one configuration
each and are genuinely uncloseable. Both carry `k_Y = 5/3`.

---

## 4. The 4,309 gated-out role records (`Z3xZ6_1_1`)

The two parser libraries do not define "a Standard-Model role field" the same way:
`z36sw_lib.role_census` assigns `role = None` to every field that is not a hidden
singlet; `z24sw_lib.role_census` has no such gate. The two geometries' role counts are
therefore **not the same quantity**, and the paper must say so.

Ungating `Z3xZ6_1_1` recovers:

| | |
|---|---:|
| role records seen **with** the gate (the number the paper prints) | 169,678 |
| role records seen **without** the gate | 173,987 |
| **dropped** | **4,309** |
| of which in the norm-6/5 family | **0** |
| of which in the norm-4/5 family | 4,309 |
| by hidden multiplicity | 2: 4,196 · 3: 110 · 4: 2 · 5: 1 |
| new `10`-role doublings after ungating | **0** |

**The norm-6/5 result is untouched**, and closing the gap strengthens it: with the gate
removed, the state-level pair test still returns zero violations on both geometries.
What must travel with the numbers is the definitional note — the `4/5`-family
multiplicity figures for `Z3xZ6_1_1` are hidden-singlet-only, and reach 5 when ungated.

---

## 5. Embedding coverage (the two-count rule)

A separate exposure: orbifolder prints a spectrum table only for the last-selected
configuration, so a multi-configuration run can leave embeddings **tagged but never
censused**. The rule is to count both, per model.

| geometry | models | embeddings tagged | embeddings printed | uncensused | recovered | **hits** |
|---|---:|---:|---:|---:|---:|---:|
| `Z3xZ6_1_1` | 4,612 | 4,667 | 4,667 | **0** | — | **0** |
| `Z2xZ4_1_1` | 11,540 | 12,169 | 11,526 | 643 | 643 | **0** |

`Z3xZ6_1_1` is **clean by construction** — nothing was uncensused. `Z2xZ4_1_1`'s gap was
closed and every recovered embedding re-censused, with zero hits.

> **REPORTED FIGURE DIFFERENCE (not reconciled here).** The campaign records two
> different `Z2xZ4_1_1` recovery counts under different scopes: **643** for the sweep
> census (the population this package uses) and **1,211** for the earlier construct-gate
> closure. They are different populations, not a contradiction, but the package does not
> attempt to reconcile them, and Paper A should quote only the 643 sweep figure if it
> quotes one at all. The campaign-wide recovery total of **1,994 embeddings across all
> seven geometries, 0 hits** spans five geometries Paper A does not use; only the two
> rows above are in Paper A's scope.

Source documents: `multi_geometry_campaign/results/uncensused_cleanup_findings.md`,
`gate_uncensused_cleanup.md`. Those two documents are *not* deposited here — they cover
the whole campaign, and this package is scoped to Paper A — but the two rows above are
quoted from them verbatim.

---

## 6. Cases outside the theorem (not exclusions, but the same discipline)

**21,981 `(cell, gauge weight)` GROUPS** across the two catalogues carry two or more
printed left-chiral records sharing one `p_sh`, distinguished only by oscillator content.
For those the theorem's `Delta` vanishes and the inequality is vacuous — they are
**outside** hypothesis H5, not counterexamples to it.

> **The unit is groups, not cells.** `d2_pshdup.py::dup_groups` counts groups, and one
> model-cell can host several; the number of distinct model-cells hosting them is smaller
> and is **not** derived in this package. Earlier text here, and in `DEFINITIONS.md` §3,
> called 21,981 a cell count. The `H5` conclusion — zero on a Standard-Model role — is
> unaffected under either unit.

| | `Z2xZ4_1_1` | `Z3xZ6_1_1` | total |
|---|---:|---:|---:|
| `(cell, p_sh)` groups of two or more records | 17,202 | 4,779 | **21,981** |
| of those, carrying **any** Standard-Model role | 0 | 0 | **0** |
| of those, carrying a **norm-6/5** role | 0 | 0 | **0** |
| oscillator-carrying norm-6/5 fields anywhere | **0** | 11 | 11 |

On `Z2xZ4_1_1` the gap is closed by arithmetic rather than by measurement: a norm-6/5
role needs `B >= 3/5`, and the smallest oscillator mode leaves `B <= 13/16 - 1/4 = 9/16
< 3/5`, so no oscillator norm-6/5 field can exist there at all. On `Z3xZ6_1_1` it is
closed only by measurement — the 11 oscillator norm-6/5 fields sit at 11 distinct cells,
`I1` exact on all 11, none doubled.

**This is the residual most likely to be raised by a referee, and it has no theoretical
backstop.** On any new geometry with a small oscillator mode it must be re-checked before
the bound is quoted. It is recorded in `limitations`, not buried in a proof.

Data: `results_deposited/gate_su5_bound/d2_pshdup_<geom>.json`,
`d2_osc10_Z3xZ6_1_1.json`.

---

## 7. Cases excluded on the construction side

| | |
|---:|---|
| **312** | class-B lifts excluded as not Standard-Model models — the SM roots are present but merged into a larger gauge factor and `t_Y` is not orthogonal to the simple roots. **160 of them retain a four-dimensional `SU(5)`.** They are excluded from the 4,032 denominator and are not counted against the negative result. Source: `g8b_smvalid.json::pattern`. |
| **1,784** | class-A models carrying a doubled norm-6/5 role on which nothing beyond counting was done (1,812 − 28). |
| **32 / 23** | catalogue model-sectors with an empty shell at `T(0,1)` / `T(0,5)`. Source: `s1_liftspace.json::histogram`, `g11_catalogue.json::hist1/hist5`. |

# Reproduction record

Every regeneration script in this package was run **inside the package**, from the
deposited raw data, when the package was assembled. This file records what ran, how long it
took, and whether the regenerated table matched the deposited one key by key.

The comparison is `compare_to_deposited.py`; the claim map is re-resolved by
`claims_check.py`; both are deposited here so the whole pass can be repeated with
`bash regeneration/run_all.sh`.

**Two passes are on disk, and both are kept.** `regeneration/logs2/` and
`run_all_status.txt` are the transcripts and the raw record of **this** pass, the one this
file describes; `run_summary.json` is its machine-readable twin, regenerated from the raw
record by `write_reproduction_table.py` so the two cannot drift. `regeneration/logs/` and
`compare_output.txt` are the earlier pass made when the package was first assembled, kept
rather than deleted so that the two can be compared. Where they disagree on a wall time,
this pass is the current measurement; no number in any table differs between them.

---

## Addendum, 2026-10-22 (later) — pre-deposit repair of v2.2

Run inside the package after the repair recorded in `CHANGELOG.md` (v2.2 repair addendum)
and `multi_geometry_campaign/data/representatives/brothers/superseded_v2.2/NOTICE.md`.
Statuses appended to `multi_geometry_campaign/data/representatives/logs/status.txt`.

| step | exit | wall time | result |
|---|---:|---:|---|
| `b5b_missing.py` (repaired) | 0 | ~1 min | 312 slots built (232 `lambda` only, 80 with a `mu` search), 0 still missing; every deck's label asserted against its vectors; `decks_900/901` identical deck for deck to the 232 correct v2.2 decks; `decks_902.txt` byte-identical to the verification's `v1b_decks.txt` |
| `b6_run_brothers.sh` (orbifolder 1.2.1, stock geometry) | 0, 0 | 41 s, 8 s | segments 900 and 901 re-run; 232 / 232 per-model blocks text-identical to the superseded output. Segment 902's output is the verification's run, copied verbatim (exit 0, 14 s) |
| `b10_brother_analysis.py` (with `label_check`) | 0 | 2,048 s | `label_check` 36,288 pass, 0 fail; `differs_from_census_rep_by_ab` (0,1) 4,026, (0,2) 4,023 (sum 32,052); every other key unchanged; per-deck table differs in exactly the 80 replaced rows |
| `b7_census_rank.py Z3xZ6_1_1` (with the `E_6` patch) | 0 | 1,521 s | 4,612 models, 165,970 groups, 0 violations; identical to the verification's patched re-run |
| `b7_census_rank.py Z2xZ4_1_1` (with the `E_6` patch) | 0 | 3,825 s | 11,526 models, 303,817 groups, 0 violations; identical to the verification's patched re-run; five rows' counts change (e.g. `5bar|33/40` `n_groups` 25,000 -> 25,006), no claim row moves |

`claims_check.py --write`: **250 of 250 resolve and match, 0 mismatches, 0 unresolved**
(30 superseded); `claims.md` regenerated. `compare_to_deposited.py`: unchanged
(76 + 1 of 77 tables).

## Addendum, 2026-10-22 — lattice representatives (package v2.2)

The seven analysis steps of `multi_geometry_campaign/code/representatives/` were run inside
the package, serially, from the deposited orbifolder output (BLAS threads pinned to 1;
logs and statuses in `multi_geometry_campaign/data/representatives/logs/`):

| step | exit | wall time | output against the original computation |
|---|---:|---:|---|
| `b4_triality_nets.py` | 0 | 1 s | byte-identical |
| `b3_compare.py` | 0 | 2 s | byte-identical |
| `b11b_compare.py` | 0 | 1 s | byte-identical |
| `b8_allbins.py` | 0 | 283 s | identical on every key; one additive key (`printed_records_oscfree_twisted`) |
| `b10_brother_analysis.py` | 0 | 1,560 s | `b10_brothers.json` byte-identical; new per-deck table |
| `b7_census_rank.py Z2xZ4_1_1` | 0 | 3,424 s | byte-identical |
| `b7_census_rank.py Z3xZ6_1_1` | 0 | 1,264 s | byte-identical |

The orbifolder steps (`b2_run_orbifolder.sh`, `b6_run_brothers.sh`) were not repeated; their
raw output is deposited. `claims_check.py`: 248 of 248 resolve and match (30 superseded).

## Addendum, 2026-09-10 — the uniform-convention census (package v2.1)

**Four** scripts and five deposited tables were added after the pass described below, and
they were run **inside the package, from the deposited raw data**, exactly as the rest were.
The table below is five *runs* of three of those four; the fourth, `uniform_census.py`,
carries the census itself and was exercised through the `E_n`-aware wrapper rather than
invoked on its own.

| script | argument | exit | wall time |
|---|---|---:|---:|
| `uniform_census/uniform_census_e6.py` | `Z3xZ6_1_1` | 0 | 1,060 s |
| `uniform_census/uniform_census_e6.py` | `Z2xZ4_1_1` | 0 | 1,703 s |
| `uniform_census/uniform_totals.py` | — | 0 | < 1 s |
| `uniform_census/rl_b_uniform.py` | `Z2xZ4_1_1` | 0 | 1,244 s |
| `uniform_census/rl_b_uniform.py` | `Z3xZ6_1_1` | 0 | 347 s |

Re-run on 2026-09-10 with the two additive mixed-budget counters in place, plus the
signature-premise recount added the same day:

| script | argument | exit | wall time |
|---|---|---:|---:|
| `uniform_census/uniform_census_e6.py` | `Z3xZ6_1_1` | 0 | 641 s |
| `uniform_census/uniform_census_e6.py` | `Z2xZ4_1_1` | 0 | 1,024 s |
| `uniform_census/uniform_totals.py` | — | 0 | < 1 s |
| `uniform_census/sig_premise_counts.py` | — | 0 | 34 s |

The two re-run accumulators agree with the ones deposited earlier on **every leaf of every
accumulator** — `uniform`, `asrun` and `current`, both geometries — except `seconds` and the
new `mixed_groups` / `mixed_fields` keys.

Six workers, BLAS threads pinned to one. `regeneration/run_all.sh` now runs all five, before
the comparison and claim-map steps.

Both censuses were also run outside the package, from the same raw data, by the fix pass that
produced them; the two runs agree on **every key of every accumulator**, `uniform`, `asrun`
and `current` alike, on both geometries. The four-dimensional-root test uses a different
engine (`root_lift/rl_b_data.py`) with its own parse and independently returns the same model
and pair counts as the census: `11,526` models and `25,883` pairs on `Z2xZ4_1_1`, `4,612` and
`7,860` on `Z3xZ6_1_1`.

`claims_check.py` after the addendum: **223 of 223 rows resolve, 0 mismatches, 0 unresolved**,
of which 28 are marked superseded — the same quantities under the previous convention, still
resolved and still true of the files they came from.

---

## Summary

| | |
|---|---:|
| scripts run | **27** |
| exited cleanly | **27** |
| exited non-zero | **0** |
| wall time, serial, six workers | **2.9 h** |
| deposited tables **actually regenerated** by this pass | **28** |
| deposited tables compared | **76** |
| matching key for key | **75** |
| matching except for absolute paths and non-regenerated commentary | **1** |
| **differing substantively** | **0** |
| claim rows re-resolved against the deposited files | **178 / 178** |

```
75 matched exactly, 1 matched up to absolute paths / non-regenerated commentary,
0 differ substantively, 0 not re-run (of 76 deposited tables)
```

**Why "178 / 178" plus "0 substantive" is the whole reproducibility claim.**
`claims_check.py` resolves every claim against `results_deposited/`, the pristine copies;
`compare_to_deposited.py` shows the freshly regenerated workspace tables equal those
pristine copies on every key. The two together mean each of the 178 numbers the paper
prints is reproduced *by re-running the code on the raw data*, not merely read back out of
a file that shipped with it. Neither check alone would establish that.

**Read the two table counts together.** 76 tables are compared; 28 of them were rewritten
by this pass, and the rest were compared as they stand. The 28 are listed at the end of this
file. The remainder are produced either by a script needing the orbifolder binary or by a
step whose driver is an inline run -- see the reproduction gaps in
`../construction/README.md` section 3. `python regeneration/write_reproduction_table.py --which`
regenerates that list from file modification times rather than from prose.

---

## The runs

Environment: CPython 3.12.0, numpy 2.2.4, Windows 11 (64-bit), six worker processes, BLAS
threads pinned to 1, scripts run **serially**. Wall times differ on other hardware; nothing
about a result depends on the worker count, because the arithmetic is exact rationals and
every tally is order-independent.

| script | argument | exit | wall time | transcript |
|---|---|---:|---:|---|
| `su5_ten_bound/verify.py` | Z2xZ4_1_1 | 0 | 10 min | `regeneration/logs2/verify_Z2xZ4_1_1.log` |
| `su5_ten_bound/verify.py` | Z3xZ6_1_1 | 0 | 6 min | `regeneration/logs2/verify_Z3xZ6_1_1.log` |
| `su5_ten_bound/pattern_census.py` | Z2xZ4_1_1 | 0 | 7 min | `regeneration/logs2/pattern_census_Z2xZ4_1_1.log` |
| `su5_ten_bound/pattern_census.py` | Z3xZ6_1_1 | 0 | 9 min | `regeneration/logs2/pattern_census_Z3xZ6_1_1.log` |
| `gate_su5_bound/d1_mechanism.py` | Z2xZ4_1_1 | 0 | 28 min | `regeneration/logs2/d1_mechanism_Z2xZ4_1_1.log` |
| `gate_su5_bound/d1_mechanism.py` | Z3xZ6_1_1 | 0 | 10 min | `regeneration/logs2/d1_mechanism_Z3xZ6_1_1.log` |
| `gate_su5_bound/d2_pshdup.py` | Z2xZ4_1_1 | 0 | 10 min | `regeneration/logs2/d2_pshdup_Z2xZ4_1_1.log` |
| `gate_su5_bound/d2_pshdup.py` | Z3xZ6_1_1 | 0 | 8 min | `regeneration/logs2/d2_pshdup_Z3xZ6_1_1.log` |
| `gate_su5_bound/d2_osc10.py` | - | 0 | 5 min | `regeneration/logs2/d2_osc10.log` |
| `gate_su5_bound/d5_excluded.py` | - | 0 | 8 min | `regeneration/logs2/d5_excluded.log` |
| `gate_su5_bound/d5_keycount.py` | - | 0 | 45 s | `regeneration/logs2/d5_keycount.log` |
| `gate_su5_bound/d6_pattern_recount.py` | Z2xZ4_1_1 | 0 | 6 min | `regeneration/logs2/d6_pattern_recount_Z2xZ4_1_1.log` |
| `gate_su5_bound/d6_pattern_recount.py` | Z3xZ6_1_1 | 0 | 4 min | `regeneration/logs2/d6_pattern_recount_Z3xZ6_1_1.log` |
| `gate_su5_bound/d8_alltwists.py` | - | 0 | 1 s | `regeneration/logs2/d8_alltwists.log` |
| `gate_su5_bound/d8_classwide.py` | - | 0 | 4 s | `regeneration/logs2/d8_classwide.log` |
| `gate_su5_bound/d9_comarks.py` | - | 0 | 20 s | `regeneration/logs2/d9_comarks.log` |
| `gate_su5_construction/g2_classes.py` | - | 0 | 8 s | `regeneration/logs2/g2_classes.log` |
| `gate_su5_construction/g4_lift.py` | - | 0 | 22 s | `regeneration/logs2/g4_lift.log` |
| `su5_construction/s0_invariants.py` | - | 0 | 26 s | `regeneration/logs2/s0_invariants.log` |
| `su5_construction/s2_sizing.py` | - | 0 | 2 s | `regeneration/logs2/s2_sizing.log` |
| `two_classes/paperA_two_classes_enum.py` | - | 0 | 13 min | `regeneration/logs2/paperA_two_classes_enum.log` |
| `root_lift/rl_a_centralizer.py` | - | 0 | 46 s | `regeneration/logs2/rl_a_centralizer.log` |
| `root_lift/rl_b_data.py` | - | 0 | 18 min | `regeneration/logs2/rl_b_data.log` |
| `root_lift/rl_c_samefp.py` | - | 0 | 12 min | `regeneration/logs2/rl_c_samefp.log` |
| `root_lift/rl_d_projection.py` | - | 0 | 17 min | `regeneration/logs2/rl_d_projection.log` |
| `gate_su5_construction/g8c_doubled_smvalid.py` | - | 0 | 1 s | `regeneration/logs2/g8c_doubled_smvalid.log` |
| `root_lift/rl_e_summary.py` | - | 0 | 1 s | `regeneration/logs2/rl_e_summary.log` |

27 scripts, 27 exited cleanly, 0 exited non-zero, 2.9 h of wall time in total.

Serial execution is deliberate. The scripts build multiprocessing pools, and on Windows two
overlapping pools can collide on a concurrent open of a geometry file (see "One transient
failure" below). `run_all.sh` runs one at a time and the collision does not arise.

---

## Comparison against the deposited tables

**Every deposited table was compared, and none differs substantively.**

One table matches only up to two classes of difference that no re-run can remove, and the
comparison script labels them rather than normalising them away:

* **`PATH`** -- the deposited tables record the absolute path of the author's working tree
  in their provenance fields; a regenerated copy records the package root. The file
  **name**, the **line number** and the **quoted text** of every provenance row still match
  exactly, and those are what carry the twist vectors. 11 keys.
* **`NOTE`** -- keys that are commentary the gate's author added to the table after the
  fact. No script emits them, so no re-run can reproduce them. They are prose, not data, and
  no claim rests on them. 147 keys, in two blocks: `_gate_notes` (45) and
  `_gate_verdict_vs_phase_7_1` (102).

```
MATCH*     gate_su5_bound/d8_classwide.json  (147 NOTE, 11 PATH -- no substantive difference)
```

> **A defect in this comparison script, found by this pass and fixed.** Its `NOTE_PREFIXES`
> list named only `_gate_notes`, so the 102 keys of the `_gate_verdict_vs_phase_7_1` block
> were reported as **substantive** differences. They are not: every one is a key present in
> the deposited copy and absent from the regenerated one, which is the signature of
> hand-added commentary. The list now matches the `_gate_` prefix, and the 158 differences
> in that file resolve as 147 NOTE + 11 PATH + **0 substantive**. The budgets themselves --
> the per-geometry keys `claims.json` actually reads -- reproduced exactly both before and
> after the fix. Recorded here rather than silently corrected, because a comparison script
> that mis-labels its own output is exactly the kind of thing this file exists to catch.

---

## Two expected non-clean outcomes

**`d9_comarks.py` prints `FAILURES` and still exits 0.** The failures are a known, recorded
defect in a helper: `s5lib.comarks("D", r)` returns a list of length `r+1` instead of `r`.
It is harmless because the only consumer reads it through `zip`, which truncates, and the
first `r` entries are correct for every rank that occurs. The script's actual finding --
that **every adjoint comes out at level 2**, for every simple type exercised -- is what
validates the level-1 identification, and it holds. Reproducing the failure lines is the
correct behaviour.

**`verify.py Z3xZ6_1_1` reports 65 loophole cells.** These are the fields at `R = 47/90` in
sectors `T(0,1)` and `T(0,5)` that the bound does **not** protect. They are the paper's
stated non-universality, not an error; their multiplicity of 1 is established by exact shell
enumeration, not by the theorem.

---

## What this pass added over the previous one

* The **two-class re-enumeration** (`two_classes/paperA_two_classes_enum.py`, 13 min) was
  run inside the package, and its output is **identical** to the deposited copy, key for key.
* The **four root_lift scripts** were run inside the package -- the lift index at all 533
  twisted constructing elements, the catalogue-wide 32,294-pair test, the
  same-fixed-point analysis, and orbifolder's own projection formula checked against
  orbifolder's printed spectrum. All four reproduce their deposited tables exactly.
* Two scripts new in this version -- `g8c_doubled_smvalid.py` (the 1,482 to 2,174
  correction) and `rl_e_summary.py` (the 6-of-533 totals) -- were run and their outputs
  deposited.
* `claims_check.py` and `claims_md.py` were written and deposited, so the claim map is now
  resolved and rendered by scripts rather than by hand. All **178** rows resolve.
* `build_manifest.py` was written and deposited, so `manifest.json` and `CHECKSUMS.sha256`
  are built by a script and can be re-verified with `--check`.

---

## How to repeat this

```
pip install -r environment/requirements.txt
bash regeneration/run_all.sh                  # everything, serially, about 3 hours
```

or piecemeal, each script from its own directory, then

```
python regeneration/compare_to_deposited.py   # workspace vs pristine reference
python regeneration/claims_check.py           # every claim row re-resolved
python regeneration/claims_md.py              # rebuild claims.md from claims.json
python regeneration/build_manifest.py --check # integrity files vs the tree
```

Each analysis script writes into `multi_geometry_campaign/data/` (and `two_classes/` for the
re-enumeration). The pristine copies in `results_deposited/` are never written to, so the
comparison stays meaningful however many times you re-run.

---

## One transient failure, and what causes it

On an earlier pass, `d2_pshdup.py Z3xZ6_1_1` died with

```
PermissionError: [Errno 13] Permission denied:
    ...\orbifolder-1.2.1\Geometry\Geometry_Z3xZ6_1_1.txt
```

**This is a Windows file-sharing artifact, not a defect in the data or the code.**
`z36sw_lib` builds its geometry object at **import** time, so every worker process opens
`Geometry_Z3xZ6_1_1.txt` as the pool starts. When several pools overlap -- a previous
script's workers still tearing down as the next one's start -- Windows can refuse one of
those concurrent opens. The file itself is intact: it reads normally, and the same script
re-run on its own succeeded and produced a table matching the deposited copy on every key.
It did not recur in this pass, which ran serially throughout.

If you hit it, re-run the single script. It has no analogue on Linux or macOS, where
concurrent reads of one file are unrestricted.

---

## What is *not* re-run here

Anything needing the orbifolder binary: the census runs themselves, and the construction and
gate runs that produced the deposited `.txt.gz` files. Those are reproducible only with an
orbifolder 1.2.1 build -- see `../environment/ENVIRONMENT.md` sections 2-3 -- and the raw
output they produced is deposited so that every number downstream of them can be checked
without one.

The two construction reproduction gaps recorded in `../construction/README.md` section 3 are
also untouched by this pass: the class-B Standard-Model-validity split has no producing
script on disk, and the class-B builder sets no random seed.

---

## The 28 tables this pass rewrote

Regenerated from the raw data, then compared against the pristine copies:

* `multi_geometry_campaign/data/gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json`
* `multi_geometry_campaign/data/gate_su5_bound/d1_mechanism_Z3xZ6_1_1.json`
* `multi_geometry_campaign/data/gate_su5_bound/d2_osc10_Z3xZ6_1_1.json`
* `multi_geometry_campaign/data/gate_su5_bound/d2_pshdup_Z2xZ4_1_1.json`
* `multi_geometry_campaign/data/gate_su5_bound/d2_pshdup_Z3xZ6_1_1.json`
* `multi_geometry_campaign/data/gate_su5_bound/d5_excluded.json`
* `multi_geometry_campaign/data/gate_su5_bound/d5_keycount.json`
* `multi_geometry_campaign/data/gate_su5_bound/d6_patterns_Z2xZ4_1_1.json`
* `multi_geometry_campaign/data/gate_su5_bound/d6_patterns_Z3xZ6_1_1.json`
* `multi_geometry_campaign/data/gate_su5_bound/d8_alltwists.json`
* `multi_geometry_campaign/data/gate_su5_bound/d8_classwide.json`
* `multi_geometry_campaign/data/gate_su5_construction/g2_classes.json`
* `multi_geometry_campaign/data/gate_su5_construction/g4_lift_A.json`
* `multi_geometry_campaign/data/gate_su5_construction/g4_lift_B.json`
* `multi_geometry_campaign/data/gate_su5_construction/g8c_doubled_smvalid.json`
* `multi_geometry_campaign/data/root_lift/guard1d_projection.json`
* `multi_geometry_campaign/data/root_lift/questionA_centralizer.json`
* `multi_geometry_campaign/data/root_lift/questionA_summary.json`
* `multi_geometry_campaign/data/root_lift/questionB_Z2xZ4_1_1.json`
* `multi_geometry_campaign/data/root_lift/questionB_Z3xZ6_1_1.json`
* `multi_geometry_campaign/data/root_lift/questionC_samefp.json`
* `multi_geometry_campaign/data/su5_construction/s0_invariants.json`
* `multi_geometry_campaign/data/su5_construction/s2_sizing.json`
* `multi_geometry_campaign/data/su5_ten_bound/patterns_Z2xZ4_1_1.json`
* `multi_geometry_campaign/data/su5_ten_bound/patterns_Z3xZ6_1_1.json`
* `multi_geometry_campaign/data/su5_ten_bound/verify_Z2xZ4_1_1.json`
* `multi_geometry_campaign/data/su5_ten_bound/verify_Z3xZ6_1_1.json`
* `two_classes/paperA_two_classes_enum.json`

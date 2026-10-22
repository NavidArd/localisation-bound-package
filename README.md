# Reproduction package — a localisation bound for Standard-Model representations in heterotic orbifolds

**Version 2.2** · **Author:** Navid Ardakanian (independent researcher) —
<n.ardakanian@gmail.com>

This is the data-and-code deposit for the paper. It contains the catalogue inputs, the
raw output of the constructor, the parsed tables, the constructed model, and every script
that produces a number the paper prints. **Nothing here reads anything outside this
directory.**

The paper's claims are computational, and this deposit exists so that a referee can
re-derive them without re-implementing anything. Every number is listed in
[`claims.json`](claims.json) with the script that regenerates it and the file that stores
it; [§4](#4-what-reproduces-what) gives one entry point per claim.

> **New in v2.2 (2026-10-22): lattice representatives.** The spectrum depends on the
> vectors `V1`, `W3 = W4`, not only on their classes modulo the lattice; v2.1 built one
> representative per pair of classes. v2.2 adds all nine representatives of each of the
> 4,032 Standard-Model-intact pairs (36,288 orbifolder spectra, decks, raw output and
> analysis), the projection check run at every constructing element and right-mover, and
> the mean-budget and rank checks on both catalogue censuses —
> [`multi_geometry_campaign/code/representatives/README.md`](multi_geometry_campaign/code/representatives/README.md).
> The seven analysis steps were re-run inside the package from the deposited raw output and
> reproduce their outputs byte for byte (`b8` gains one additive key). Claim map after v2.2:
> **248 of 248** rows resolve and match, 30 marked superseded (`CHANGELOG.md`).
>
> **Repaired before deposit (2026-10-22).** An independent verification found that 80 of the
> 36,288 decks — the slots `(0,1)`, `(0,2)` of 40 class-B pairs — were copies of their pair's
> `(0,0)` deck under `(0, b)` labels (`b5b_missing.py` never searched the Wilson-line shift
> `mu`). The builder is fixed, the 80 genuine decks and their orbifolder output are in, every
> deck's label is now re-checked against its vectors (36,288 / 36,288), and `b7` parses the
> four `E_6`-first models it had missed. No printed number moves; `REP-06` becomes 32,052.
> Claim map: **250 of 250** resolve and match (30 superseded). The superseded files are kept in
> `multi_geometry_campaign/data/representatives/brothers/superseded_v2.2/` (`NOTICE.md`).
>
> **Reproduction status, as deposited (v2.1).** All **227** mapped claims resolve
> against the deposited files and match; 28 of them are marked superseded, being the
> same quantities under the counting convention the paper used before v2.1. The
> uniform-census scripts added in v2.1 --- four of them, run five times, together with
> the signature-premise recount added on 2026-09-10 --- were run inside the package from
> the raw data and all exited cleanly (`regeneration/REPRODUCTION.md`, addendum). The paragraph
> below records the v2.0 pass, which is unchanged.
>
> **Reproduction status, v2.0 pass.** All **178** mapped claims resolve against the
> deposited files and match. All **27** regeneration scripts were run end to end inside
> this package, from the raw data, in one serial pass of 2.9 hours; all exited cleanly, and
> all **76** deposited tables compare with **zero substantive differences** — including the
> 28 this pass rewrote from the raw output. Details, wall times, and the one table that
> differs only in absolute paths and hand-added commentary:
> [`regeneration/REPRODUCTION.md`](regeneration/REPRODUCTION.md).
>
> Four known gaps travel with the package and are listed in
> [`CHANGELOG.md`](CHANGELOG.md) — two on the construction side, one superseded file, one
> unit whose finer count is not derived here — together with two premises closed by
> measurement rather than by proof, and one reported mismatch between a source document's
> prose and its own data.

---

## 1. Start here

| if you want to | read |
|---|---|
| check one number | [`claims.md`](claims.md) — every claim, its script, its file, its status |
| know what a number counts | [`DEFINITIONS.md`](DEFINITIONS.md) — the units, and the ratios most often misread |
| see the bound working, not just forbidding | [`controls/CONTROLS.md`](controls/CONTROLS.md) |
| see the model where the repeat happens | [`construction/CB1534/`](construction/CB1534/) |
| check that the search was exhaustive | [`two_classes/`](two_classes/) — the boundary re-enumerated from scratch |
| know whether the hidden-`SU(2)` signature is real, and where | [`multi_geometry_campaign/data/root_lift/`](multi_geometry_campaign/data/root_lift/) and §8 below |
| know what was left out and why | [`exclusions/EXCLUSIONS.md`](exclusions/EXCLUSIONS.md) |
| rebuild the raw output from scratch | [`environment/ENVIRONMENT.md`](environment/ENVIRONMENT.md) |
| know what was re-run when this package was assembled | [`regeneration/REPRODUCTION.md`](regeneration/REPRODUCTION.md) |
| trace a number back to the catalogue | [`dataflow.md`](dataflow.md) |

---

## 2. Software required

| | |
|---|---|
| Python | 3.10 or later (3.12.0 used) |
| Python packages | **`numpy` only** (2.2.4 used). No scipy, no sympy |
| the constructor | **orbifolder 1.2.1** (Nilles, Ramos-Sánchez, Vaudrevange, Wingerter, [arXiv:1110.5229](https://arxiv.org/abs/1110.5229)), GPL v3 |

`orbifolder` is needed **only** to regenerate the raw output from scratch. Every script in
`multi_geometry_campaign/code/` runs against the deposited raw output and needs Python
alone.

The arithmetic that matters is exact — budgets, weight norms, `R` and every threshold are
`fractions.Fraction`, never floating point — so no result depends on the BLAS or on the
worker count.

```
pip install "numpy>=1.24"
```

Set `OMP_NUM_THREADS=MKL_NUM_THREADS=OPENBLAS_NUM_THREADS=1` before running anything;
every script does so at module top, and it affects wall time only.

---

## 3. Layout

```
claims.json  claims.md     every printed number -> script -> file -> status
DEFINITIONS.md            what each number counts; the multiple denominators, fixed
dataflow.md               catalogue -> raw output -> parsed table -> printed number
CHANGELOG.md              what is new in v2.0 and what is carried from v1.0
manifest.json             every file with its SHA-256, size and role
CHECKSUMS.sha256          the same hashes in sha256sum format

orbifolder-1.2.1/         redistributed upstream input: COPYING, AUTHORS, README,
                          configure.ac, Dockerfile, and all 138 Geometry files (GPL v3)
alphabet_survey/          the certified geometry/norm engines the parsers import,
                          and their two phase-2 tables
multi_geometry_campaign/  THE WORKSPACE - code and data at the paths the scripts expect
    code/su5_ten_bound/           the bound: verify, pattern census, shell enumeration
    code/gate_su5_bound/          the SECOND IMPLEMENTATION of the central counts
    code/su5_construction/        the class-A enumeration
    code/gate_su5_construction/   the second implementation, and the class-B enumeration
    code/root_lift/               the hidden-SU(2) signature: the lift index at all 533
                                  constructing elements, and the catalogue-wide test
    code/representatives/         v2.2: the nine lattice representatives of every
                                  Standard-Model-intact embedding, the projection check
                                  at every bin, the mean-budget and rank checks
    code/z2xz4_sweep/  z3xz6_sweep/  z3xz6/   the certified orbifolder-output parsers
    data/literature/parr2020_full/           the two catalogue slices consumed
    data/z2xz4_sweep/raw/census/             58 gzipped raw orbifolder segments
    data/z3xz6_sweep/raw/census/             48 gzipped raw orbifolder segments
    data/su5_construction/                   the 2,912-model deck, 30 raw segments, tables
    data/gate_su5_construction/              the 1,432-model deck, 3 raw runs, tables
    data/root_lift/                          the signature tables, and their own README
    data/representatives/                    v2.2: 36,288 decks in 183 segments, their
                                             raw output, the model-level decks and outputs
    data/su5_ten_bound/  data/gate_su5_bound/  the parsed tables
two_classes/              the referee-standard re-enumeration of the class boundary,
                          written from the mathematics and importing nothing from the
                          campaign's libraries: code, output, log, and the
                          breadth-first orbit certificate written out row by row,
                          with its independent checker
results_deposited/        pristine copies of every parsed table, for comparison
regeneration/             what was re-run here, its logs, the comparison, the claim
                          resolver and the claims.md renderer
construction/             CB1534 (flagship) and CB1041 (demoted), decks and raw output
controls/                 the three positive controls
exclusions/               every excluded, failed, duplicate and recovered case
environment/              orbifolder version and build, drivers, patches, requirements
```

**Why the workspace mirrors the campaign's directory names.** The deposited scripts
resolve their inputs relative to a repository root. Reproducing that layout means the
scripts run unmodified. The only edit made to any script is described in
[`environment/PATCHES_TO_PACKAGE_CODE.md`](environment/PATCHES_TO_PACKAGE_CODE.md): in
thirteen files, the single hard-coded literal `F:\physics` is replaced by the package
root. Every other script is byte-identical to the original, apart from five files that are
new in this package and have no original to differ from. The twelve scripts of
`code/representatives/` (v2.2) were written for this package and resolve every path from
their own location; their README lists the three additive changes made when they were
copied in.

**`results_deposited/` versus `multi_geometry_campaign/data/`.** The same parsed tables
appear in both. The workspace copies are what the scripts read and **overwrite** when you
re-run them; `results_deposited/` is the pristine reference, and is what `claims.json`
resolves against. Re-run freely — the reference copies are never touched.

---

## 4. What reproduces what

Run each from its own directory. Wall times are **measured**, from the serial pass
recorded in `regeneration/REPRODUCTION.md`, on six workers.

### The theorem's counts

| entry point | establishes | time |
|---|---|---:|
| `code/su5_ten_bound/verify.py Z2xZ4_1_1` | the full-catalogue test: models, cells, fields, `\|w\|^2` by role, `R` histograms, per-(cell, role) multiplicity, violations, `k_Y`, invariant `I1` | 10 min |
| `code/su5_ten_bound/verify.py Z3xZ6_1_1` | the same, and the 65 unprotected fields at `R = 47/90` | 6 min |
| `code/su5_ten_bound/pattern_census.py <geom>` | the `(Q, u^c, e^c)` occupancy patterns; the distinct-cell counts 26,796 and 16,136; the two derived-impossible patterns at zero | 7-9 min |
| `code/su5_ten_bound/shell.py` | exact massless-shell enumeration at the unprotected cells | |
| `code/su5_ten_bound/localgut.py`, `measure.py`, `recon.py`, `coloc_probe.py` | the local-group structure and the premise-memo evidence | |

### The second implementation of the central counts

The referee asked for an independent re-derivation. These scripts were written against the
same raw output with their own lattice-membership test, their own Weyl-group construction,
their own accumulator, their own regexes and their own role dictionary.

| entry point | establishes | time |
|---|---|---:|
| `code/gate_su5_bound/d1_mechanism.py <geom>` | **the proof tested step by step on all 31,822 doubled pairs**; `\|Delta\|^2 = 2` exactly in every one; hidden multiplicity; the Gram ceiling against measurement | 10-28 min |
| `code/gate_su5_bound/d6_pattern_recount.py <geom>` | an independent recount of the occupancy census — exact agreement on all fourteen numbers | 4-6 min |
| `code/gate_su5_bound/d5_spotcheck.py` | 30 cells re-extracted with independent regexes and an independent branch criterion: 161 fields compared, 160 agree (the one difference is the definitional gate of `exclusions/EXCLUSIONS.md` §4) | |
| `code/gate_su5_bound/d5_keycount.py`, `d5_cellkey.py`, `d5_keysplit.py` | the cell key certified **bijective** — it cannot split one fixed point and hide a doubling | |
| `code/gate_su5_bound/d5_shellkey.py` | key-convention independence, at shell level | |
| `code/gate_su5_bound/d7_shell.py`, `d8_shell_complete.py` | a **second, independent** shell enumerator; the two agree row by row | |
| `code/gate_su5_bound/d2_pshdup.py <geom>` | the 21,981 cells outside hypothesis H5, and that none carries a Standard-Model role | 8-10 min |
| `code/gate_su5_bound/d2_osc10.py` | the 11 oscillator norm-6/5 fields on `Z3xZ6_1_1`, none doubled | 5 min |
| `code/gate_su5_bound/d5_excluded.py` | the complete exclusion accounting | 8 min |
| `code/gate_su5_bound/d9_comarks.py` | every adjoint at level 2 — the level-1 identification is sound | 20 s |

### The budgets table

| entry point | establishes | time |
|---|---|---:|
| `code/gate_su5_bound/d8_alltwists.py` | the twist vectors of all seven geometries, read off the primary geometry files | 1 s |
| `code/gate_su5_bound/d8_classwide.py` | `delta_c` and `B` per sector, sector counts, `B_max` / `B_min`, the SUSY condition, and which sectors permit a repeat — with a provenance record (file and line) for every twist row | 4 s |
| `code/gate_su5_bound/d8_crosscheck.py`, `d8_verdict.py` | the recomputed budgets against the measured ones | |

### The construction

| entry point | establishes | time |
|---|---|---:|
| `code/gate_su5_construction/g2_classes.py` | the class boundary as originally published. **Its intermediate counts are superseded** — see the row below — but its answer, exactly two, is not | 8 s |
| `two_classes/paperA_two_classes_enum.py` | **the class boundary re-enumerated from the mathematics, importing nothing from the campaign's libraries.** 26,400 solutions in 1,260 lattice classes reducing to **exactly two**, with an orbit certificate reaching every class from the two representatives. Supersedes `g2_classes.py`'s 3,000 and 119, and reconciles with them explicitly | 13 min |
| `code/gate_su5_construction/g4_lift.py` | both lift spaces: 750,384 -> 2,912 and 696,600 -> 1,432 | 22 s |
| `code/su5_construction/s5_build.py` | the 2,912 class-A embeddings and the orbifolder deck | |
| `code/gate_su5_construction/g7_buildB.py` | the 1,432 class-B embeddings | |
| `code/su5_construction/s11_battery2.py` | the class-A funnel: `Tr Y = 0` in 481 | |
| `code/gate_su5_construction/g5_battery.py` | the triality-aware recount: net `Q = ±3` in 339, three generations in all three roles in **0** | |
| `code/gate_su5_construction/g8_batteryB.py` | the class-B funnel, three generations in all three roles in **0** | |
| `code/gate_su5_construction/g1c_compare.py` | CB1534 reproduced 147/147 from a fresh deck built from the embedding vectors alone | |
| `code/gate_su5_construction/g8c_doubled_smvalid.py` | how many of the 4,032 Standard-Model-valid models carry a doubled norm-6/5 role: **2,174**, correcting the 1,482 that was carried in prose | 1 s |
| `code/su5_construction/s0_invariants.py` | the catalogue's own modular-invariance validation, 4,612/4,612 | 26 s |

### Lattice representatives (v2.2)

| entry point | establishes | time |
|---|---|---:|
| `code/representatives/b5_build_brother_census.py`, `b5b_missing.py` | the nine `(a, b)` representatives of each of the 4,032 Standard-Model-intact pairs: 36,288 decks, all passing the complete modular-invariance set, each with its slot `(a, b)` re-computed from its own vectors (`b5b` asserts it; `b10` re-checks all 36,288). `b5b` repaired before deposit — see the note below this table | |
| `code/representatives/b6_run_brothers.sh` | orbifolder on all 36,288 (needs an orbifolder build) | 0.5 h on 6 workers |
| `code/representatives/b10_brother_analysis.py` | 36,288 spectra, every deck's label re-checked against its vectors (36,288 pass); three generations at once in **0**; two of three in 216; a doubled 10-like representation in 21,168 (supersedes the 2,174 of `g8c_doubled_smvalid.py`, which is the `(0,0)` slice); the `(0,0)` nets equal the certified ones on 4,032/4,032 | 26 min |
| `code/representatives/b1_brothers.py`, `b3_compare.py`, `b11_v2_brothers.py`, `b11b_compare.py` | on the model (7.1): nine distinct spectra; every lattice shift of `V1`, `W3`, `V2` tried reproduces one of them | 5 s |
| `code/representatives/b8_allbins.py` | the projection code of `rl_d_projection.py` at all 140 bins of CB1534 and CB1041, both chiralities: 0 mismatches, every empty bin empty | 5 min |
| `code/representatives/b7_census_rank.py <geom>` | the mean-budget ceiling (0 violations) and the rank restriction (tighter than the norm only at `R = 83/90`, 12 -> 9) on both catalogues, all 11,526 + 4,612 models (with the `E_6` banner parse of `uniform_census_e6.py`) | 57 + 21 min |

**Pre-deposit repair (2026-10-22).** As first assembled, `b5b_missing.py` searched the
shift of `V1` only and set `mu = 0`, so the 80 slots `(0,1)`, `(0,2)` of 40 class-B pairs held
copies of the `(0,0)` deck under `(0, b)` labels. The repaired builder searches `mu`; the 80
genuine decks and their orbifolder output are in `decks_902.txt` / `out_902.txt.gz`; the
superseded files are kept in `data/representatives/brothers/superseded_v2.2/` with a
`NOTICE.md`. No printed number moves (0 / 216 / 21,168 / 36,288 / 4,032); `REP-06` moves
from 31,972 to 32,052. `b7` now applies the `E_6` parse patch (it had covered 16,134 of the
16,138 models); its verdicts are unchanged. Details: `CHANGELOG.md`, v2.2 repair addendum.

Two construction numbers have reproduction gaps, reported in
[`construction/README.md`](construction/README.md) §3: the class-B
Standard-Model-validity split (**1,120** / 312) has no producing script on disk, and the
class-B builder sets no random seed.

### The signature, at its stated scope

| entry point | establishes | time |
|---|---|---:|
| `code/root_lift/rl_a_centralizer.py` | the lift index at **all 533** twisted constructing elements of the seven geometries — forced at 6, and at **none** of `Z2xZ4_1_1`'s 74 | 46 s |
| `code/root_lift/rl_b_data.py <geom>` | the catalogue-wide test: of 32,294 same-role record pairs, **zero** have a `Delta` that is a root of the unbroken four-dimensional group | 20 min |
| `code/root_lift/rl_c_samefp.py` | what "the same fixed point" must mean, and why a key built on the local shift would manufacture doublings | 4 min |
| `code/root_lift/rl_d_projection.py` | orbifolder's own projection formula, reimplemented from its source and run against its printed spectrum on CB1534 and CB1041 | 17 min |
| `code/root_lift/rl_e_summary.py` | the totals the paper quotes, 6 of 533 | 1 s |

### The claim map itself

| entry point | establishes | time |
|---|---|---:|
| `regeneration/claims_check.py` | every row of `claims.json` re-resolved from the deposited files and compared with what the paper prints — this is what makes the `status` column a measurement | 2 s |
| `regeneration/claims_md.py` | rebuilds `claims.md` from `claims.json`, so the two cannot drift | 1 s |
| `regeneration/compare_to_deposited.py` | the regenerated workspace tables against the pristine reference copies, key by key | 10 s |
| `regeneration/run_all.sh` | all of the above, serially, end to end | 2.9 h |

---

## 5. Redistributed material and licences

| what | licence |
|---|---|
| the code in `multi_geometry_campaign/code/` and `alphabet_survey/code/` | **MIT** — `LICENSE-CODE.txt` |
| the data and documentation in this package | **CC BY 4.0** — `LICENSE-DATA.txt` |
| `orbifolder-1.2.1/` and `environment/orbifolder_patches/patched_geometry_files/` | **GNU GPL v3** — not covered by either of the above; see the NOTICE in each directory |
| `multi_geometry_campaign/data/literature/parr2020_full/` | ancillary files of arXiv:2003.01732 (Parr, Vaudrevange, Wimmer); see the NOTICE there |

The orbifolder **program** is not redistributed — only its licence, its authorship files,
its build recipe, and the seven geometry input files the analysis reads.

---

## 6. Size and distribution

> **v2.2 (2026-10-22; after the pre-deposit repair): 1,171 files, 445 MB** (`manifest.json`). The
> representative census adds 183 gzipped orbifolder outputs and 183 deck segments under
> `multi_geometry_campaign/data/representatives/` (plus the eight superseded files of
> `brothers/superseded_v2.2/`, 11 MB); still no file above 10 MB, and still well
> inside GitHub's and Zenodo's limits. The breakdown below is the v2.0 one; `manifest.json`
> carries the current per-file roles and sizes.

**734 files, 203 MB.** No file exceeds 10 MB; the largest is a 9.4 MB gzipped
orbifolder run. Where it goes:

| | files | MB |
|---|---:|---:|
| raw orbifolder output, gzipped | 139 | 158.8 |
| parsed tables (workspace copies) | 100 | 14.5 |
| upstream catalogue slices and archives | 5 | 13.5 |
| pristine reference copies of the tables | 76 | 11.8 |
| engine input tables | 2 | 1.3 |
| upstream geometry files and licences (all 138 geometries) | 144 | 1.1 |
| analysis code | 145 | 1.1 |
| the constructed models: decks, vectors, raw output | 13 | 0.5 |
| everything else (docs, logs, drivers, patches, controls) | 110 | 0.8 |

`manifest.json` carries this breakdown per file, with a role for each;
`python regeneration/build_manifest.py` regenerates both it and `CHECKSUMS.sha256`.

**The package is deposited whole — code, data and raw output in one archive.** No split is
needed. 203 MB is comfortably inside a normal repository: well under GitHub's 1 GB soft
limit and its 5 GB hard one, and under Zenodo's 50 GB per record. The raw output is the
part a referee most needs, because it is what every downstream number is computed from.

If a git host's limits ever make the raw output awkward, the natural split is along a line
the layout already draws: everything except `multi_geometry_campaign/data/*/raw/` (the 139
gzipped files, 159 MB) is **44 MB** and stands on its own, since the parsed tables and the
pristine reference copies are both deposited. The raw output would then live on the archive
only. **That split is not applied here**, and applying it later would mean a referee could
no longer re-derive a table from the constructor's own output without a second download —
which is the whole point of depositing it.

Git LFS is not used and is not needed; no file is near the 100 MB blob limit. `.gitignore`
excludes only `__pycache__/` and `*.pyc`, which are the sole files not checksummed.

---

## 7. Integrity

`CHECKSUMS.sha256` covers every file; `.gitattributes` disables end-of-line conversion, so
`sha256sum -c CHECKSUMS.sha256` verifies on Windows clones as well as Unix ones.
`manifest.json` carries the same hashes with a size and a role for each file. Both are
built by a deposited script, so they can be regenerated and checked rather than trusted:

```
python regeneration/build_manifest.py --check     # verify against the tree, write nothing
```

The `Z2xZ4_1_1` catalogue file additionally has its SHA-256 and byte count **recorded
inside the analysis code**, as `z24sw_lib.CATALOGUE_SHA256` and `CATALOGUE_BYTES`. Note
what that does and does not do: the constants are **not** checked at parse time — they are
checked by `code/z2xz4_sweep/selftest.py` (check `S11`), which must be run deliberately.
Run it before trusting a rebuilt catalogue:

```
cd multi_geometry_campaign/code/z2xz4_sweep && python selftest.py
```

The deposited file matches: 7,560,329 bytes, SHA-256 `6c73de64…a6d2b7d`.

---

## 8. The hidden-`SU(2)` signature — deposited, at a stated scope

Whether the root `Delta` separating two states at a saturating fixed point is a root of the
**four-dimensional** gauge group, or only of the **local** group at that point, is the
sharpest question asked of this paper. It is now settled in both directions, and the files
are in `multi_geometry_campaign/data/root_lift/` with their own
[`README`](multi_geometry_campaign/data/root_lift/README.md).

* **In general the lift fails, and that is measured, not conceded.** Of the **533** twisted
  constructing elements across the seven campaign geometries, the lift is forced by the
  geometry at **6**. On `Z2xZ4_1_1` — this paper's other catalogue — it is forced at
  **none** of 74, and no statement of the signature's form is licensed there. Across both
  catalogues, of the **32,294** same-role pairs printed as two records, **zero** have a
  `Delta` that lifts.
* **At the scope the paper claims it, it is a theorem.** The claim is scoped to sectors
  `T(0,1)` and `T(0,5)` of `Z3xZ6_1_1`. Each has a **single, translation-free**
  constructing element, and this geometry's only permitted non-trivial Wilson line lies in
  the plane the twist leaves fixed — so `Delta` pairs integrally with all eight generators
  and the lift index is exactly **1**. Both facts are read off the primary geometry file.
* **Two premises are load-bearing and one of them is closed by measurement only.** The two
  states must share right-mover and oscillator excitation. Nine index-1 pairs exist, all at
  the signature's own two sectors, and all nine carry *different* oscillator excitations —
  the case where the cancellation fails. The signature survives because no *shared*-
  oscillator index-1 pair exists in either catalogue: a measurement over 16,162 models, not
  a theorem, and it must be re-checked on any new geometry.

**Not licensed:** "a local root is a four-dimensional root"; any use of the signature on
`Z2xZ4_1_1`; any extension to constructing elements carrying a lattice translation,
including `T(0,2)`, `T(0,3)`, `T(0,4)` of this same geometry.

Entry points: `code/root_lift/rl_a_centralizer.py` (the index at all 533 elements),
`rl_b_data.py` (the catalogue-wide test), `rl_c_samefp.py`, `rl_d_projection.py`
(orbifolder's own projection formula, reproduced), `rl_e_summary.py` (the totals).
Claims `SIG-01`…`SIG-06`.

---

## 9. What this package does not contain

* **The wider campaign.** Five further geometries, the coupling runs, the R-symmetry
  census and the arithmetic-alphabet work are outside this paper's scope and are not
  deposited. The budgets table's five arithmetic-only rows are the one place the paper
  reaches beyond the two measured catalogues, and they are labelled as arithmetic there
  and in `claims.json`. The one exception is `root_lift`, whose §8 index table covers all
  seven geometries, because the *negative* half of the signature statement needs them.
* **The host-cell count behind `21,981`.** That number counts `(cell, gauge weight)`
  **groups**; how many distinct model-cells host them is a smaller number and is not
  derived here. See `DEFINITIONS.md` §3.
* **A DOI.** See `CITATION.cff`: the v2.0 DOI **does not exist yet** and is marked as a
  placeholder. The DOI printed in earlier drafts of the manuscript belongs to v1.0 of a
  **different** paper's package and must not be reused.

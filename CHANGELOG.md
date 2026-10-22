# Changelog

## v2.2, repair addendum (before deposit) — 2026-10-22

**A labelling defect in the representative census, found by the independent verification
of the paper's v10 draft, and a parse gap in `b7`. Both repaired before deposit. No number
the paper prints moves (0 / 216 / 21,168 / 36,288 / 4,032); two unprinted rows do.**

* **The defect.** `b5b_missing.py`, which fills the 312 slots `b5` left unbuilt, searched
  the shift `lambda` of `V1` only and wrote `mu = 0` for every slot. For the 80 slots
  `(a, b) = (0,1), (0,2)` of 40 class-B pairs (`GB0980` ... `GB1203`) it therefore wrote
  decks with `lambda.V2 = 0` and `mu = 0` — realisations of the `(0,0)` phase class — under
  `(0, b)` labels. Their spectra are record-identical to the pair's `(0,0)` spectrum. No step
  re-computed a deck's slot from its vectors, and the claim map compares numbers, not labels,
  so `REP-01` ("nine per pair, one per `(a, b)`") and the README's "312 slots needed the wider
  lattice search of `b5b`" were false of the deposit at 80 of 36,288 slots.
* **The fix.** `b5b_missing.py` now searches `mu` when `b != 0`, keeps `mu = 0` and the
  original `lambda` order when `b = 0` (the 232 correct decks are reproduced vector for
  vector), and asserts for every deck it writes that `(3 lambda.V2 mod 3, 3 mu.V2 mod 3)`,
  re-computed from the deck's vectors, equals the slot label. The 80 genuine decks are
  labelled `Y<parent>a0b<b>` and sit in a new segment `decks_902.txt`; they and their
  orbifolder output `out_902.txt.gz` are **copied verbatim** from the verification
  (`v1b_decks.txt`, `v1c_out.txt.gz`), and the repaired builder regenerates `decks_902.txt`
  byte for byte. `decks_900.txt`, `decks_901.txt` now hold the 232 `b = 0` decks only and
  were re-run through orbifolder with `b6` (every per-model block text-identical to the
  previous output). The repaired `b5b` also writes `build_all.json`, which had no script.
* **The check that would have caught it.** `b10_brother_analysis.py` now reads every deck
  back and re-computes its `(a, b)` from its vectors against its label (key `label_check`):
  **36,288 pass, 0 fail**; on the superseded files it gives 36,208 / 80. New claim rows
  `REP-01c` (36,288) and `REP-01d` (0).
* **What moved.** `REP-06`, non-`(0,0)` spectra differing from their pair's `(0,0)`:
  **31,972 -> 32,052** (`(0,1)` 3,986 -> 4,026; `(0,2)` 3,983 -> 4,023). The coarse
  coincidences fall from 284 to 204. `REP-01`'s wording corrected. Unchanged: 36,288 spectra
  (`REP-01`), 0 missing, the `(0,0)` control 4,032 / 4,032, **0** three-generation
  (`REP-03`), **216** two-of-three (`REP-04`), **21,168** doubled (`REP-05`; the 80 genuine
  decks are doubled, as the 80 copies were). The per-deck table differs from the previous one
  in exactly the 80 replaced rows.
* **Kept and labelled.** The previous `b5b_missing.py` (hash-identical), `decks_900/901`,
  `out_900/901`, `build_missing.json` and `build_all.json` are in
  `multi_geometry_campaign/data/representatives/brothers/superseded_v2.2/`, with a
  `NOTICE.md` recording the defect, the provenance of the replacement files and their hashes.
  `environment/PATCHES_TO_PACKAGE_CODE.md` records every code edit of this repair.
* **`b7_census_rank.py` and the `E_6` parse.** `b7` imported `uniform_census` without the
  `E_6` banner-parse repair of `uniform_census_e6.py`, which the paper's tables use, and so
  covered 16,134 of the 16,138 models (`MSSM967`, `MSSM1409`, `MSSM825`, `MSSM2669`
  missing). Each worker now applies `uniform_census_e6._patch()`. Re-run on both catalogues:
  11,526 + 4,612 models; 0 violations of the mean-budget ceiling; mixed groups 197 and 9;
  the rank restriction tighter than the norm at `R = 83/90` only, at 238 / 238 groups; group
  counts now equal the paper's tables (5bar-like 231,323 and 128,747; 10-like 72,494 and
  37,223). The `b7` claim rows (`REP-09` ... `REP-10a`) are unchanged; the 33/40 row's
  `n_groups` moves 25,000 -> 25,006 (not a claim row). The "known residual" bullet of the v2.2
  entry below described this gap as a grouping difference; it was not one, and the bullet is
  corrected in place.
* `claims_check.py`: **250 of 250 resolve and match, 0 mismatches, 0 unresolved**
  (30 superseded). `claims.md`, `manifest.json`, `CHECKSUMS.sha256` regenerated;
  `build_manifest.py` labels the `superseded_v2.2/` files as superseded. `CITATION.cff`
  unchanged (version 2.2, DOI placeholder).

## v2.2 — 2026-10-22

**The spectrum depends on the lattice representatives of the shift vectors and Wilson lines,
not only on their classes. v2.1 built one representative per pair of classes; v2.2 builds
all nine. The paper's three-generation negative stands on all 36,288; two of its counts
change. No theorem moves.**

* **Why.** Shifting `V1` or `W3 = W4` by a lattice vector leaves the gauge group, the
  untwisted sector and every class-level condition of the paper's section 7.2 unchanged,
  but can change the projection phases of twisted states — the equivalence between lattice
  shifts and discrete torsion of Ploger, Ramos-Sanchez, Ratz, Vaudrevange, *Mirage Torsion*
  (hep-th/0702176, sections 3.1, 3.3, 4.1). On `Z3xZ6_1_1` the phases depend only on
  `(a, b) = 3 (lambda.V2, mu.V2) mod 3`, which modular invariance leaves free. The premise
  of the `su5_construction` enumeration that "only their classes mod 3 Lambda matter for
  the spectrum" is false at fixed (zero) discrete torsion.
* **New: `multi_geometry_campaign/code/representatives/`** (twelve scripts and a README),
  **`multi_geometry_campaign/data/representatives/`** (36,288 decks in 182 segments, 182
  gzipped orbifolder outputs, build records, run log; the model-level decks and outputs)
  and **`results_deposited/representatives/`** (nine output tables, including the per-deck
  table `b10_per_deck_nets.csv.gz` of net indices and the doubled-10-like flag).
  * the representative census: 36,288 spectra, 0 missing; the `(0,0)` nets equal the
    certified per-model nets on 4,032 / 4,032; three net generations of `Q`, `u^c`, `e^c`
    at once in **0**; two of the three in **216** (was 22); a doubled 10-like representation
    in **21,168** (was 2,174, which is the `(0,0)` slice); 31,972 of the 32,256 non-`(0,0)`
    spectra differ from their pair's `(0,0)` spectrum;
  * on the model (7.1): nine distinct spectra; 18 lattice shifts of `V1` and `W3` and nine
    of `V2` each reproduce one of them;
  * the projection code (`root_lift/rl_d_projection.py`, imported unchanged) at **every**
    one of the 140 (twisted constructing element, massless right-mover) bins of CB1534 and
    of CB1041, both chiralities: 0 mismatches; every empty bin empty; it reproduces the 88
    oscillator-free twisted multiplets of CB1534;
  * on both catalogue censuses: the mean-budget ceiling `n <= 1/(1 - Rbar)` with **0**
    violations; the mixed-budget groups (197 at `Rbar = 23/40`, 9 at `Rbar = 34/45`); the
    rank restriction (`n <= 9` in one `E8`), tighter than `floor(1/(1-R))` only at
    `R = 83/90` (12 -> 9), at all 238 groups there.
* **`claims.json`** gains the `REP-*` rows; `CON-24` (2,174) and `CON-25` (362) are marked
  superseded by `REP-05`. `claims_check.py`: **248 of 248 resolve and match, 0 mismatches, 0 unresolved** (30 superseded). `claims.md` regenerated.
* **Re-run inside the package** from the deposited raw output: `b4`, `b3`, `b11b` (outputs
  byte-identical to the originals), `b8` (identical on every original key, plus the new
  key `printed_records_oscfree_twisted`), `b10` and `b7` (both byte-identical to the originals; `b10` 26 min, `b7` 57 min and 21 min). The orbifolder steps
  (`b2`, `b6`) were not repeated.
* **`run_all.sh`** gains the seven analysis steps. `README.md`, `CITATION.cff` (version
  2.2), `manifest.json` and `CHECKSUMS.sha256` updated.
* **Known residual — corrected by the repair addendum above.** As first written, this bullet
  said that `b7_census_rank.py` groups states slightly differently from `uniform_census.py`'s
  `(cell, representation)` slots (for example 128,643 against 128,747 5bar-like groups on
  `Z3xZ6`). That was wrong: the groups are identical, and the gap was four models that `b7`
  did not parse (no `E_6` banner patch). With the patch the counts agree exactly. No `b7`
  count is printed in the paper, which states the mean-budget ceiling as a theorem.

## v2.1, second addendum — 2026-10-22

**Three package-side corrections found while checking the paper's Data-availability text
against this deposit. No number, theorem or verdict moves.**

* **The breadth-first certificate is now stored.** `two_classes/paperA_two_classes_enum.py`
  computed a parent and a connecting generator for every one of the 1,260 enumerated
  classes and asserted completeness, but wrote only the count (`CON-03b`). It now also
  writes `two_classes/paperA_two_classes_bfs_certificate.json` (pristine copy under
  `results_deposited/two_classes/`): one row per class, seeded at the two representatives
  the paper prints (1,200 rows from `V2^A`, 60 from `V2^B`, maximum depth 4). The new
  `two_classes/check_bfs_certificate.py` reads only that file and verifies every parent
  link, every generator, that every chain ends at a printed representative, and that the
  1,260 classes are distinct. The enumeration's main output is **byte-identical** to the
  deposited copy; its log gains one line. `claims.json` gains `CON-03d` (1,260);
  `claims_check.py`: **227 of 227 resolve and match, 0 mismatches, 0 unresolved** (28
  superseded); `claims.md` regenerated. `run_all.sh` runs the checker after the enumeration.
* **The paper no longer claims that this package holds the running-mass arithmetic behind
  its value of `r`.** It never did; the sentence is removed from the paper rather than the
  files added, since nothing in the paper rests on that arithmetic.
* **`CITATION.cff` and this README's title** now carry the paper's title and the totals it
  prints (16,138 models, 880,589 twisted model-cells, 502,396 fields in a Standard-Model
  representation; `results_deposited/uniform_census/uniform_totals.json`). They carried the
  title and the pre-v2.1 totals of an earlier draft (15,535 / 851,231 / 481,728).
* **`construction/CB1534/README.md`:** dated note on the row `t_Y . t_anom = 0`, which is
  vacuous — the model has no anomalous `U(1)`; every printed abelian anomaly trace is zero.
* `manifest.json` and `CHECKSUMS.sha256` regenerated. (`CITATION.cff` had been edited after
  the previous regeneration without one; that drift is closed by the same run.)

## v2.1, addendum — 2026-09-10

**Three numbers the paper prints were not in the claim map; they are now, at their
uniform-convention values.** The paper's signature section printed two record counts
(`31,053` and `12,406`) taken as-run, with ten `Z3 x Z6` models counted twice, and its
ceiling-annexe caption had dropped a mixed-budget numeral (`195`) that the v2.1
accumulators could not recompute. All three are recomputed here and the paper moves to
`30,922`, `12,350` and `197`. **No verdict, no theorem and no row of the paper's
verification table moves.**

* **New:** `multi_geometry_campaign/code/uniform_census/sig_premise_counts.py` and
  `results_deposited/uniform_census/sig_premise_counts.json`. The script recounts the
  `Z3 x Z6` catalogue at `T(0,1)` and `T(0,5)` under three model populations — as-run,
  each model once, and each model once minus the two `E_n`-first models — and combines the
  "each model once" figures with the constructed-model half read from
  `gate_su5_construction/g6_signature.json`. It reproduces the as-run `31,053` and
  `12,406` as well, so the decomposition is the one that produced the previous values.
  The `E_n` repair is shown to contribute nothing: a raw record count does not depend on
  the Standard-Model closure test, so the whole movement is the de-duplication.
* **Extended:** `uniform_census.py` gains two additive counters, `mixed_groups` and
  `mixed_fields` — `(cell, role)` groups holding records of two distinct `R`, keyed
  `family|R_max` and `family|R_max|R`. `uniform_totals.py` deposits them under
  `mixed_budget`. Both censuses were re-run inside the package; the new accumulators agree
  with the deposited ones on **every leaf of every accumulator, both geometries**, save
  `seconds` and the new keys. On the previous convention the counter returns `195`, the
  numeral the paper used to print, which is what identifies `197` as its successor.
* **`claims.json`** gains `SIG-07` (30,922), `SIG-08` (12,350) and `UNI-46` (197).
  `claims_check.py`: **226 of 226 resolve and match, 0 mismatches, 0 unresolved**, of which
  28 are superseded.
* **`run_all.sh`** gains `sig_premise_counts.py`. `manifest.json` and `CHECKSUMS.sha256`
  regenerated: **745 files, 203.5 MB** (`manifest.json` carries the exact byte total, which
  this file's own size enters), all 744 hashes verify.
* **Corrections to this version's own documentation.** `README.md` and
  `regeneration/REPRODUCTION.md` said "five uniform-census scripts"; there are **four**,
  run five times — `uniform_census.py` is the census itself and is exercised through the
  `E_n`-aware wrapper rather than invoked on its own. The "New in the package" list below
  omits `rl_b_uniform.py`, which is on disk, in the manifest, in `CHECKSUMS` and run twice
  by `run_all.sh`; and its `run_all.sh` bullet says "all three" where five runs are meant.
* **R-5 widened.** `DEFINITIONS.md`'s dated note listed only §§2.3, 2.7, 2.9, 2.11 as
  stale. Its preamble, its §3 quick-reference table and its §4 are stale too, and the
  preamble's bolded *"the ratio to print is 2.413"* contradicted both the note above it and
  the paper. The note now says so and names `2.414` (`UNI-30`) as the value.

## v2.1 — 2026-09-10

**One counting convention, and a repaired banner parse.** Nothing in the science moves; the
paper's verification table is now counted one way on both geometries, with every catalogue
model counted once, and the four models a defective regex had been dropping are back in.

### What changed

* **One role convention, both geometries.** A record carries a Standard-Model role if its
  colour, isospin and hypercharge assignment says so, with no condition on the representation
  it carries under the remaining gauge factors. Previously the field and multiplicity rows used
  a narrower rule on `Z3 x Z6` (hidden singlets only) while the pair row did not, and the
  difference — 472 pairs — was disclosed in the table's caption rather than removed. It is now
  removed.
* **One embedding rule, and each model counted once.** Each catalogue model is censused at the
  first printed embedding whose Standard-Model identification closes. On `Z2 x Z4` that brings
  in the 599 models that close only at a later printed configuration, which the previous
  aggregates excluded; on `Z3 x Z6` the ten models this program's own pilot segment re-runs are
  counted once rather than twice.
* **The `E_6` banner parse.** `z24sw_lib.parse_gauge_factors` and `z36sw_lib.parse_banner`
  matched `SU(n)`, `SO(n)` and `E(n)` but not `E_6`, so on a model printing an `E_6` factor
  every later factor was shifted by one position. Where the `E_6` is printed after the
  Standard-Model factors this changes no answer; where it is printed first — `MSSM967`,
  `MSSM1409`, `MSSM825`, `MSSM2669` — the Standard-Model identification cannot be rebuilt and
  the model was skipped. All four close under the repaired parse, all four at `k_Y = 5/3`, none
  with a violation. With them in, the `Z3 x Z6` census is complete at 4,612 of 4,612.
* **A state unit beside the record unit.** Every count the printed output resolves is now
  given twice: once per printed record, and once per state, a record of hidden dimension `d`
  standing for the `d` gauge momenta of equal `|p|^2` it collects. The bound is checked in
  both. The state test is the stronger one, because two states inside one record share `R`.

### Headline movements

| quantity | v2.0 | v2.1 |
|---|---:|---:|
| distinct models censused | 15,535 | **16,138** |
| twisted model-cells | 851,231 | **880,589** |
| role-carrying fields (records / states) | 481,728 | **502,396 / 518,373** |
| of norm 6/5 | 104,889 | **109,717** |
| of norm 4/5 (records / states) | 376,839 | **392,679 / 408,656** |
| (cell, norm-4/5 role) pairs | 346,151 | **360,070** |
| same-role record pairs | 32,294 | **33,743** |
| matched-unit ratio to print | 2.413 | **2.414** |
| **violations of the bound** | **0** | **0**, and **0** in the state unit |
| **unprotected norm-6/5 fields** | **65** | **65** |
| **(cell, gauge weight) groups carrying a role** | **0** | **0** |

Two per-geometry rows move **down** on `Z3 x Z6` — 296,309 model-cells to 295,734 and 37,334
norm-6/5 fields to 37,223 — because the de-duplication removes more there than the four
recovered models and the wider role convention add. Both both-catalogue totals still rise.

### New in the package

* `multi_geometry_campaign/code/uniform_census/` — `uniform_census.py` (the census),
  `uniform_census_e6.py` (the same census with the two banner parsers rebound at runtime; no
  library file is edited), `uniform_totals.py` (the both-catalogue reduction the claim map
  resolves against).
* `results_deposited/uniform_census/` — the two accumulators, their reduction, and the
  four-dimensional-root test re-run over the same pair population.
* `regeneration/run_all.sh` runs all three, before the comparison and claim-map steps.
* `claims.json` gains a `verification (uniform convention)` section, one row per number the
  paper's verification table prints. The rows carrying the previous convention's values are
  marked `superseded_by`; they are still resolved and still compared, because a superseded
  number must remain true of the file it came from.

### Known gaps added in this version

* **R-5** — `DEFINITIONS.md` §§2.3, 2.7, 2.9, 2.11 still tabulate the previous convention's
  values. A dated note at its head says so and points at the new rows; the tables themselves
  are not rewritten.
* **R-6** — `root_lift/rl_b_data.py`'s `notperp` counter is a dead branch: it builds its
  Standard-Model basis from integer factor positions inside a bare `try/except`, so the
  orthogonality test never executes and the counter cannot increment. No number in the package
  or the paper is taken from it. The orthogonality of all 33,743 same-role pairs is measured
  independently, against a basis recovered per model by explicit projection, by
  `uniform_census_e6.py`: zero failures.
* **R-7** — the two banner parsers are repaired at runtime only. The certified libraries and
  the `su5_ten_bound` census they feed are unchanged, so the files under
  `results_deposited/su5_ten_bound/` still carry the pre-repair values; that is why the rows
  resolved against them are marked superseded rather than rewritten.

---

## v2.0 — 2026-09-05

**A new package for a different paper.** Version 1.0 accompanied *A selection-rule and
decoupling audit of a heterotic `Z3 x Z3` vacuum* and covered the 77-model `Z3 x Z3`
catalogue. This package accompanies the localisation-bound paper and covers the
`Z2 x Z4` and `Z3 x Z6` catalogues. **No file is carried over.** They share an author, a
constructor, and a set of conventions — nothing else.

> **The DOI in earlier drafts of the manuscript pointed at v1.0 of the other paper's
> package.** That citation was wrong and must not be repeated. No DOI exists for this
> package yet; `CITATION.cff` marks the field as a placeholder.

### What is in v2.0

* **Two complete catalogues, raw.** 106 gzipped orbifolder segments — **15,535 distinct
  models** over 15,545 census runs, **851,231 model-cells** (the 128 twisted constructing
  elements counted once per model), 481,728 role-carrying fields — with the parsers, the
  parsed tables, and the scripts that regenerate every number.
* **A second independent implementation of the central counts**, as the referee asked for:
  its own lattice-membership test, its own Weyl-group construction and search, its own
  accumulator, its own regexes, its own role dictionary, its own chirality criterion, and a
  second shell enumerator. The proof is tested step by step on all 31,822 doubled pairs,
  not just at its conclusion — and the mechanism is independently re-tested on the larger,
  convention-homogeneous set of 32,294.
* **The hidden-`SU(2)` signature, deposited with its scope.** The lift index at all 533
  twisted constructing elements of the seven campaign geometries; the catalogue-wide test
  showing that of 32,294 same-role record pairs **zero** have a `Delta` that is a root of
  the unbroken four-dimensional group; and the derivation that at `Z3xZ6_1_1`'s `T(0,1)`
  and `T(0,5)` — and there only — the index is exactly 1 and the lift is forced.
* **The class boundary re-enumerated to referee standard**, from the mathematics, importing
  nothing from the campaign's own libraries: `two_classes/`.
* **The budgets table's inputs**: the seven geometry files, and the script that reads the
  twist vectors off them with a per-row record of the file and line each came from.
* **The construction**: CB1534's input deck, embedding vectors and complete raw output,
  plus the same output re-derived by an independent deck; CB1041 kept and marked
  **DEMOTED** with the reason; both two-class enumerations — the original and the
  independent re-derivation that supersedes its intermediate counts — with their complete
  output and certificates, and the 4,032-model funnel data.
* **The three positive controls**, each with its data file — the bound forbidding below
  threshold, permitting above it, and capping at exactly the predicted number.
* **A machine-readable claim map** (`claims.json`, 178 rows) in which every value is
  *resolved from the deposited files* rather than transcribed, so the status column is a
  measurement — **and the resolver itself, `regeneration/claims_check.py`, is now
  deposited**, together with `claims_md.py`, which renders `claims.md` from it so the two
  cannot drift. Every reduction a row applies (a sum, a maximum, a length, a ratio) is
  declared in the row rather than performed implicitly.
* **A definitions section** fixing every counted unit. The referee flagged that the paper
  uses several denominators over the same data with ratios up to 2.5 between them;
  `DEFINITIONS.md` names each one, gives the matched-unit ratio to print (**2.413**), and
  names the two that must not be printed.
* **A complete log of excluded, failed, duplicate and recovered cases** —
  including the two-count rule's tagged-versus-printed embedding numbers per model.
* **An environment specification**: the orbifolder release identified by the version
  string the binary itself printed into the raw output, the build recipe, the WSL
  invocation, the exact command block, pinned Python requirements, and SHA-256 for every
  catalogue file consumed and every raw output file.
* **A record of the runs made when this package was assembled** —
  `regeneration/REPRODUCTION.md` — with wall times, outcomes, and any number that did not
  reproduce.

### What changed in the underlying results since the last internal freeze

These are corrections to the science, folded in here so that the deposit and the paper
agree:

* **The `H_u` label.** `(1, 2)_{+1/2}` lies in the `5`, not the `5bar`; it is the
  conjugate of `L`. No number moves — `|w|^2` is invariant under conjugation — but the
  label was wrong and is corrected throughout. The package's primary vocabulary is now
  **norm-6/5** and **norm-4/5 roles**, with the `SU(5)` labels as declared shorthand.
* **The triality correction.** The construction's role table discarded a colour-conjugation
  sign, counting exotics as quark doublets. The negative result is unchanged (still zero
  models with three net generations in all three roles), the class-A doubled count is
  unchanged (1,812), and `Tr Y = 0` is unchanged (481); but net `Q = ±3` moves 326 → 339,
  and the flagship changes from CB1041 to CB1534. Both models are deposited.
* **Hypothesis H5 made explicit.** The proof's original justification for `Delta != 0` was
  false: 21,981 `(cell, gauge weight)` **groups** carry two records sharing one gauge
  momentum, distinguished only by oscillator content. None carries a Standard-Model role.
  The gap is closed by **measurement**, not by proof, and is now stated as a limitation
  rather than hidden inside a proof.
* **The Gram ceiling** `n <= 1/(1−R)` is stated with its domain condition `R < 1` as a
  hypothesis, not as an empirical remark, and is checked against measurement at every
  realised `R`.
* **"Independent pipeline" retracted.** The two chains share the Standard-Model
  identification byte for byte; only the tallying layer is independent. `dataflow.md` §3
  says exactly what is shared and what is not, and names the evidence that *is* independent.
* **The 599 recovered models.** The absence claim covers 11,524 of 11,540 `Z2xZ4_1_1`
  models, not 10,925.
* **The `Z3xZ6_1_1` duplication is ours, not the catalogue's.** Ten models are censused
  twice because a pilot segment of *this program's own run* re-runs `MSSM0…MSSM9`, which
  also appear in segment 0. The published file has 4,612 blocks, all labels and all
  contents distinct. v1.0-era text called it a catalogue defect; that was wrong.
  Consequence: the paper prints **15,535 distinct models** over 15,545 census runs. The
  effect is over-counting, so every absence claim is safe, and the headline ratio is
  unmoved (2.4129 de-duplicated against 2.4131 as run).
* **The two catalogues read opposite embedding blocks.** `Z2xZ4_1_1` takes the first
  printed hypercharge block, `Z3xZ6_1_1` takes the last. Earlier text asserted "the first"
  for both. Measured impact: none — on all 55 two-block models both choices give identical
  role tallies. Documentation only, but it was asserted rather than measured.
* **The denominators are separated.** `851,231` counts **model-cells**, not fixed points —
  only 128 twisted cells exist across both geometries. `21,981` counts `(cell, gauge
  weight)` **groups**, not cells. `3,112,629` counts left-chiral records **at twisted
  cells**; with the untwisted sector it is 3,298,075. And the matched-unit denominator of
  `104,889` is **43,466**, not `42,932`: the numerator counts conjugate roles and the
  latter does not. The ratio to print is **2.413**; `2.48` is `Z2xZ4_1_1` alone.
* **`31,822` → `32,294`.** The smaller figure sums two different role conventions across
  the two geometries; the larger applies one to both, and is a superset on which every
  conclusion still holds.
* **The doubled-model count: 1,482 → 2,174.** Over the 4,032 Standard-Model-valid models
  the correct count is 1,812 (class A) + 362 (class B, doubled **and** valid). The 400 in
  `g8_batteryB.json` is over all 1,432 *built* class-B lifts and includes 38 that are not
  Standard-Model models. `g8c_doubled_smvalid.py` now derives all three figures. **The
  negative result is untouched:** zero models reach three generations in all three roles.
* **The class-boundary intermediate counts: 3,000 → 26,400 solutions, 119 → 1,260 lattice
  classes.** The published enumeration imposed `rem = 1/2` where the derivation gives
  `rem >= 1/2`, skipping a whole budget, and its class key quotiented by a lattice strictly
  containing `E8 (+) E8`, merging classes. **The answer is unchanged at exactly two**, and
  the re-enumeration carries an orbit certificate reaching every class from the two
  representatives. The superseded rows are kept and labelled.
* **The hidden-`SU(2)` signature moved from "asserted" to "derived, at a stated scope"** —
  and, in general, to "refuted, and measured to be refuted." See `README.md` §8.
* **CB1534's threshold margin: `1/45` → `1/90`.** `B − 17/20 = 31/36 − 17/20 = 1/90`.
  `1/45` is the margin in `R` (`47/90 − 1/2`), and the two differ by the factor of 2 in
  `R = 2B − |w|²`. Both are now printed side by side in `controls/CONTROLS.md` and
  `construction/CB1534/README.md`. The threshold is met either way; only the stated margin
  was wrong. Two arithmetic slips beside it: `2*sqrt(47/90)` is `1.4453`, not `1.4456`
  (`construction/CB1041_DEMOTED/README.md`), and it exceeds `sqrt(2)` by **2.2%**, not 3%
  (`controls/CONTROLS.md`). Neither changes a verdict; both are the kind of number a
  referee checks on a calculator.

### Known gaps in this version

Listed so that nobody has to discover them.

* **R-1** — the class-B Standard-Model-validity split (1,120 / 312) has no producing
  script on disk. The 4,032 headline depends on it. Its inputs are deposited.
* **R-2** — the class-B builder sets no random seed, so its representative vectors are not
  bit-for-bit reproducible. Its counts are.
* **R-3** — `su5_ten_bound/classwide.json`, the phase's own budgets table, has no
  producing script on disk. It is superseded and is **not** the source of any claim here:
  the budgets table is taken from `gate_su5_bound/d8_classwide.json`, which
  `d8_classwide.py` regenerates from the primary geometry files with a per-row record of
  the file and line each twist came from. The superseded file is kept, and
  `d8_verdict.py` compares the two.
* **Reported mismatch** — the source documents say "18 genuinely unrecoverable"
  `Z2xZ4_1_1` models; the arithmetic and the on-disk breakdown both say **16**, and an
  independent recount from the raw output confirms 16 directly. Reported, not resolved in
  the source documents; no downstream number is affected.
* **R-4** — the number of distinct model-cells hosting the 21,981 `(cell, gauge weight)`
  groups is not derived in this package. The unit correction is made; the host-cell count
  would need a second pass over the raw records and `d2_pshdup.py` does not emit it. No
  claim rests on it.
* **The signature's oscillator premise is closed by measurement, not by proof.** All nine
  index-1 pairs in the two catalogues carry *different* oscillator excitations — the case
  where the cancellation fails. The signature survives because no shared-oscillator index-1
  pair exists in either catalogue. That is a measurement over 16,162 models, and on any new
  geometry with a small oscillator mode it must be re-checked. The same premise is what
  hypothesis H5 rests on, on the bound's side.
* **The signature's boxed converse is proved one way.** "Index 1 ⟺ translation-free ∧
  Wilson lines in the fixed plane" is established as a **sufficient** condition only.

### Conventions carried from v1.0

Layout and discipline, not content: one entry point per claim; raw constructor output
deposited gzipped with the command files beside it; `CHECKSUMS.sha256` over every file with
`.gitattributes` disabling end-of-line conversion so the checksums verify on Windows
clones; MIT for code, CC BY 4.0 for data, GPL v3 for redistributed orbifolder inputs, each
with its own NOTICE; superseded files kept and labelled rather than deleted, so that a
correction can be audited.

Extended in v2.0: the same "keep and label" discipline now applies to **numbers**, not
only files. A superseded claim keeps its row in `claims.json`, carries a `superseded_by`
field, and is still resolved and compared on every run — a superseded number has to remain
true of the file it came from — but it is reported separately and is not what the paper
prints. Four rows are in that state: `A-18`, `CON-01`, `CON-02`, `CON-19`.

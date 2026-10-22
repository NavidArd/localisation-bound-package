# Superseded: the 80 mislabelled `(0, b)` decks of the v2.2 representative census

**Read this first: nothing in this directory enters any number the package or the paper
prints.** These are the files of the v2.2 deposit that the pre-deposit repair of
2026-10-22 replaced. They are kept, byte-identical, because the package keeps and labels a
superseded file rather than deleting it, so that the correction can be audited.

## The defect

The representative census builds, over each of the 4,032 Standard-Model-intact pairs of
classes, one embedding for each of the nine slots `(a, b) = 3 (lambda.V2, mu.V2) mod 3`,
where `V1 -> V1 + lambda` and `W3 = W4 -> W3 + mu`. `b5_build_brother_census.py` left 312
slots unbuilt — `(1,0)` x 116, `(2,0)` x 116, `(0,1)` x 40, `(0,2)` x 40, all over class-B
pairs — and `b5b_missing.py` filled them from a wider lattice search.

As deposited, `b5b_missing.py` searched **`lambda` only**, accepted a `lambda` with
`3 lambda.V2 = a (mod 3)`, and wrote **`mu = 0` for every slot**. For the 232 slots with
`b = 0` that is correct. For the 80 slots with `a = 0, b != 0` (40 class-B pairs,
`GB0980` ... `GB1203`, times `(0,1)` and `(0,2)`) it accepted a `lambda` with
`lambda.V2 = 0 (mod 1)` and `mu = 0`, and labelled the deck `(0, b)`. Those 80 decks realise
the `(0,0)` phase class, not `(0, b)`: their `mu.V2` is `0`, not `b/3`. Their spectra are
record-identical to their pair's `(0,0)` spectrum (80 / 80, at the level of sector, fixed
point, representation, chirality, shifted momentum and oscillators).

No check re-computed a deck's slot from its vectors: the label was carried by the build
record, and `b10_brother_analysis.py` read `(a, b)` from that record. `claims_check.py`
compares numbers with outputs and cannot see deck labels. The defect was found by the
independent verification of the paper's v10 draft
(`paper15/review/paperA_v10_verification.md`, section 1.3, in the author's working tree),
which re-computed `(3 lambda.V2, 3 mu.V2) mod 3` for every deck: 36,208 pass, 80 fail.

## What it changed, and what it did not

With the 80 genuine `(0, b)` embeddings in place of the 80 copies:

| quantity | v2.2 as first deposited | repaired | printed in the paper |
|---|---:|---:|---:|
| spectra parsed | 36,288 | 36,288 | 36,288 |
| three generations at once (`REP-03`) | 0 | 0 | 0 |
| two of the three (`REP-04`) | 216 | 216 | 216 |
| doubled 10-like (`REP-05`) | 21,168 | 21,168 | 21,168 |
| non-`(0,0)` spectra differing from their `(0,0)` (`REP-06`) | 31,972 | **32,052** | not printed |
| non-`(0,0)` spectra coinciding with their `(0,0)` (coarse: record count and nets) | 284 | **204** | not printed |

**No number the paper prints moves.** The 80 genuine spectra contribute 0 three-generation
and 0 two-of-three cases and are all doubled, exactly as the 80 copies were. The 204
remaining coincidences are coarse only: none is record-identical to its `(0,0)`. What the
defect did falsify, at 80 of 36,288 slots, was the description of the deposit: "nine per
pair, one per `(a, b)`" (`REP-01`), and "312 slots needed the wider lattice search of
`b5b`", which presented the 312 as correctly filled.

## The files here

| file | what it is | SHA-256 |
|---|---|---|
| `b5b_missing_v2.2.py` | the builder as first deposited (searches `lambda` only, writes `mu = 0`) | `8709c065...cd85a3` |
| `decks_900.txt`, `decks_901.txt` | its 312 decks: the 232 correct `(a, 0)` decks and the 80 mislabelled `X<parent>a0b1`, `X<parent>a0b2` | `f19d0541...38f41`, `fb2433cd...dedda4` (full hashes in `manifest.json`) |
| `out_900.txt.gz`, `out_901.txt.gz` | orbifolder 1.2.1 on those decks (the `b6` run of the v2.2 deposit) | in `manifest.json` |
| `build_missing.json` | that builder's record (its `meta` says `(0, b)` for the 80) | in `manifest.json` |
| `build_all.json` | the merged build record `b10` read in v2.2 | in `manifest.json` |

## What replaced them (in `..`, the census directory)

* **`b5b_missing.py`, repaired.** It searches `mu` whenever `b != 0` (and `lambda`, from
  `lambda = 0`, whenever `a = 0`), keeps `mu = 0` and the original `lambda` order for
  `b = 0`, and asserts for every deck it writes that `(3 lambda.V2 mod 3, 3 mu.V2 mod 3)`,
  re-computed from the deck's own vectors, equals the slot's label. It writes the 232
  `b = 0` decks to `decks_900.txt` (200) and `decks_901.txt` (32) — **vector for vector the
  232 correct decks above** — and the 80 `mu`-search decks to `decks_902.txt`, labelled
  `Y<parent>a0b<b>` so that no superseded label is reused for a different object. It also
  writes `build_all.json`, which in v2.2 had no script.
* **`decks_902.txt` and `out_902.txt.gz` are copied verbatim from the verification.** They
  were built by `paper15/review/code_paperA_v10_verify/v1b_missing_b_slots.py` (SHA-256
  `d130f64f...4da19e`) and run through the same orbifolder 1.2.1 binary, command sequence
  and stock `Geometry_Z3xZ6_1_1.txt` as `b6` by `v1c_run_orbifolder.sh` (`876dadda...767c5e`;
  exit 0, 14 s). Source files, in the author's working tree:
  `paper15/review/data/paperA_v10/v1b_decks.txt` (`8f3de83a...35086`) and
  `.../v1c_out.txt.gz` (`c9dcbee5...322c0`); the copies here have the same hashes. The
  repaired `b5b_missing.py` regenerates `decks_902.txt` **byte for byte**, which is the
  check that the builder and the verification's decks are the same objects. All 80 have
  their parent's gauge group. That run is not in `../run.log`, which records only `b6`
  runs.
* **`out_900.txt.gz`, `out_901.txt.gz` were re-run** with `b6_run_brothers.sh` on the
  repaired `decks_900.txt`, `decks_901.txt` (same binary and stock geometry; `../run.log`,
  the lines after the second `ALLDONE`). Each of the 232 per-model blocks is
  **text-identical** to its block in the superseded outputs here; the only byte that differs
  is the session's closing "Do you really want to quit?" prompt, which now follows a
  different last model.
* **`b10_brother_analysis.py`** reads both label prefixes and now re-computes every deck's
  slot from its vectors (`b10_brothers.json::label_check`: 36,288 read, 36,288 pass, 0 fail;
  run against the files in this directory it returns 36,208 pass and 80 fail, the 80 above).
  `claims.json` rows `REP-01c`, `REP-01d` record the check; `REP-06` moves to 32,052;
  `REP-01`'s wording is corrected.

## Licence

Data files: CC BY 4.0 (`LICENSE-DATA.txt`); `b5b_missing_v2.2.py`: MIT (`LICENSE-CODE.txt`).
The `.txt.gz` files are orbifolder output, as everywhere else in this package.

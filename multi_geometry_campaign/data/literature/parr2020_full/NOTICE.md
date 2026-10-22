# Redistributed model catalogues — attribution and licence

The catalogue files in this directory are **not the work of the author of the accompanying
paper**. They are the ancillary data files of a published paper, redistributed here
because every number in this package is computed from them and a reproduction package that
did not carry its own inputs would not be one.

## Authorship

> E. Parr, P. K. S. Vaudrevange, M. Wimmer, *Predicting the orbifold origin of the
> MSSM*, [arXiv:2003.01732](https://arxiv.org/abs/2003.01732).

The arXiv ID was fetched directly and its title and author list read off the record.

## Contents

| file | what it is |
|---|---|
| `Z2xZ4_OrbifoldMSSMs.tar.xz` | the upstream archive, as distributed with the paper |
| `Z3xZ6_OrbifoldMSSMs.tar.xz` | the same, for `Z3 x Z6` |
| `extracted/Z2xZ4_1_1.txt` | the `Z2xZ4_1_1` space group's slice, extracted from the archive above: **11,540** models, each with a label, a lattice, two shift vectors and six Wilson lines |
| `extracted/Z3xZ6_1_1.txt` | the same for `Z3xZ6_1_1`: **4,622 records under 4,612 distinct labels** — the ten duplicated labels are upstream, and are reported rather than removed; see `../../../../exclusions/EXCLUSIONS.md` §2 |

Only these two space groups are deposited. The upstream archives cover seven point groups;
the rest are outside this paper's scope.

## What was and was not modified

The `.tar.xz` archives are byte-identical to the upstream ancillary files. The two
`extracted/*.txt` files are the corresponding members of those archives, extracted without
modification — no reformatting, no de-duplication, no filtering. The census drivers slice
them into segments with `awk` at run time; that slicing writes into a scratch directory
and never modifies these files.

## Integrity

Both extracted files are covered by the package's `CHECKSUMS.sha256`.

`Z2xZ4_1_1.txt` additionally has its SHA-256 and byte count recorded **inside the analysis
code**, as `z24sw_lib.CATALOGUE_SHA256` and `CATALOGUE_BYTES`. The constants are
**declared, not enforced at parse time**; they are checked by
`../../../code/z2xz4_sweep/selftest.py`, check `S11`, which must be run deliberately. Run
it before trusting a re-extracted catalogue. The deposited file matches: 7,560,329 bytes,
SHA-256 `6c73de64aeddb5688a50789d8afcf58ff4cb1abe598ad12d49a7c4a39a6d2b7d`.

## Licence

These files are the work of the authors above and are redistributed under the terms their
arXiv ancillary files carry. They are **excluded** from this package's CC BY 4.0 licence;
see `LICENSE-DATA.txt` at the package root.

If you use these catalogues, cite arXiv:2003.01732.

# The one edit made to the deposited code

Every analysis script in this package is the script that produced the paper's numbers.
**Thirteen** of them carried the absolute path of the author's working tree as a literal.
Those thirteen are edited; **every other file is byte-identical to the original.**

---

## What was changed

The literal `F:\physics` — and the two longer literals built from it — replaced by the
package root, resolved at import time:

```python
import os as _os
_PKG_ROOT = _os.environ.get('PAPERA_PKG') or _os.path.abspath(
    _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..'))
```

and then, at the original line:

| original | in this package |
|---|---|
| `REPO = r"F:\physics"` | `REPO = _PKG_ROOT` |
| `CAMP = r"F:\physics\multi_geometry_campaign"` | `CAMP = _os.path.join(_PKG_ROOT, 'multi_geometry_campaign')` |
| `GEODIR = r"F:\physics\orbifolder-1.2.1\Geometry"` | `GEODIR = _os.path.join(_PKG_ROOT, 'orbifolder-1.2.1', 'Geometry')` |
| `OUT = r"F:\physics\multi_geometry_campaign\data\gate_su5_bound\d8_alltwists.json"` | the same path rebuilt from `_PKG_ROOT` |

The shim block is inserted immediately after the module docstring, and after any
`from __future__ import ...` line — that import must remain the first statement.

Nothing else was touched: no logic, no constant, no threshold, no output path relative to
the root.

---

## The thirteen files

| file | which literal it carried |
|---|---|
| `multi_geometry_campaign/code/su5_ten_bound/s5lib.py` | `REPO` |
| `multi_geometry_campaign/code/gate_su5_bound/d5_keycount.py` | `REPO` |
| `multi_geometry_campaign/code/gate_su5_bound/d8_alltwists.py` | `GEODIR`, `OUT` |
| `multi_geometry_campaign/code/gate_su5_bound/d8_classwide.py` | `GEODIR`, `OUT` |
| `multi_geometry_campaign/code/gate_su5_bound/d8_crosscheck.py` | `CAMP` |
| `multi_geometry_campaign/code/gate_su5_bound/d8_shell_complete.py` | `CAMP` |
| `multi_geometry_campaign/code/gate_su5_bound/d8_verdict.py` | `CAMP` |
| `multi_geometry_campaign/code/su5_construction/sc_lib.py` | `REPO` |
| `multi_geometry_campaign/code/gate_su5_construction/gs_lib.py` | `REPO` |
| `multi_geometry_campaign/code/z2xz4_sweep/z24sw_lib.py` | `REPO` |
| `multi_geometry_campaign/code/z2xz4_sweep/protection.py` | `REPO` |
| `multi_geometry_campaign/code/z2xz4_sweep/f1_pilot_repair.py` | `REPO` |
| `two_classes/paperA_two_classes_enum.py` | `GEOM` (one geometry file) |

Search any of them for `reproduction-package root shim` to find the block.

**`two_classes/` sits one level below the package root, not three**, so its shim resolves
the root with a single `..`. Everything else about that file is unchanged.

**`multi_geometry_campaign/code/root_lift/` needed no edit at all.** Those five scripts
already derive every path from their own location, so they are byte-identical to the
originals — as are the three files in `two_classes/` apart from the `GEOM` line above.
Two scripts are **new in this package** and have no original to differ from:
`gate_su5_construction/g8c_doubled_smvalid.py` and `root_lift/rl_e_summary.py`, both of
which only sum figures already stored in deposited tables, plus `regeneration/`'s
`claims_check.py`, `claims_md.py` and `run_all.sh`.

---

## Why the package mirrors the campaign's directory names

The alternative — flattening the tree into a tidier layout — would have required rewriting
path expressions inside dozens of scripts. Each rewrite is a chance to point a script at
the wrong file, and a reader would then have to trust that none of them did. Reproducing
the original directory names means the scripts run unmodified, and the only thing a reader
has to check is a single substitution, listed above.

That is why `multi_geometry_campaign/`, `alphabet_survey/` and `orbifolder-1.2.1/` appear
at the package root with those names.

---

## Relocating the package

The shim derives the root from each script's own location, so the package works wherever
it is unpacked. If you move individual files out of the tree, set the root explicitly:

```
set PAPERA_PKG=D:\somewhere\paperA_repro_package        (Windows)
export PAPERA_PKG=/somewhere/paperA_repro_package       (Unix)
```

---

## Verifying that nothing else changed

Every non-shimmed `.py` file should hash identically to its counterpart in the source
tree. For the twelve above, the diff should be exactly the inserted block plus the
right-hand side of the listed assignments. All 115 Python files in the package compile
under Python 3.12 (`py_compile`, zero failures) — checked when the package was assembled.

---

## Pre-deposit repair of v2.2 (2026-10-22): what was changed

The twelve scripts of `multi_geometry_campaign/code/representatives/` were written for this
package (they have no campaign original). Three of them, and the manifest builder, were
edited after the v2.2 assembly and before deposit. This section is the record of those
edits; nothing else in the package's code changed.

| file | what changed | why |
|---|---|---|
| `code/representatives/b5b_missing.py` | searches the shift `mu` of `W3 = W4` whenever the slot's `b != 0` (and `lambda` from `lambda = 0` when `a = 0`); keeps `mu = 0` and the original `lambda` order for `b = 0`; asserts, for every deck written, that `(3 lambda.V2 mod 3, 3 mu.V2 mod 3)` re-computed from the deck's vectors equals the slot label; labels `mu`-search decks `Y<parent>a<a>b<b>` and writes them to the segments after the `lambda`-only ones; writes the merged `build_all.json` | **defect**: as deposited it searched `lambda` only and wrote `mu = 0`, so 80 decks labelled `(0,1)`/`(0,2)` (40 class-B pairs) realised `(0,0)`. The version as deposited is kept, hash-identical, as `data/representatives/brothers/superseded_v2.2/b5b_missing_v2.2.py`; the defect, the files it produced and their replacements are recorded in `superseded_v2.2/NOTICE.md` |
| `code/representatives/b10_brother_analysis.py` | new function `label_check`: reads every deck back and re-computes its `(a, b)` from its vectors against its label and build record (also `V2` = parent's, `W4 = W3`); output key `label_check`; the label pattern reads both `X` and `Y` prefixes | the check that would have caught the defect: `b10` had taken `(a, b)` from the build record. Run on the superseded files it returns 36,208 pass / 80 fail |
| `code/representatives/b7_census_rank.py` | each worker calls `uniform_census_e6._patch()` (imported unchanged) before `uniform_census.worker` | as deposited it imported `uniform_census` without the `E_6` banner-parse repair and so covered 16,134 of the 16,138 models (`MSSM967`, `MSSM1409`, `MSSM825`, `MSSM2669` missing). The deposited CHANGELOG had described the resulting count gap as a grouping difference; it was not one |
| `regeneration/build_manifest.py` | files under `superseded_v2.2/` get the role "superseded v2.2 file, kept and labelled" (checked before the `.txt.gz` rule) | so that the manifest does not list the superseded orbifolder output as live raw output |

Data files added or replaced by the repair: `brothers/decks_900.txt`, `decks_901.txt`
(regenerated: the 232 correct `b = 0` decks, vector for vector as before),
`brothers/decks_902.txt` and `out_902.txt.gz` (**copied verbatim** from the v10
verification's `v1b_decks.txt` / `v1c_out.txt.gz`; the repaired `b5b` regenerates
`decks_902.txt` byte for byte), `brothers/out_900.txt.gz`, `out_901.txt.gz` (re-run with
`b6_run_brothers.sh`; every per-model block text-identical to the superseded output),
`brothers/build_missing.json`, `build_all.json`, and the re-run outputs of `b10` and `b7`
with their pristine copies under `results_deposited/representatives/`.

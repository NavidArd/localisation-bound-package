# R-symmetry geometry patches — NOT on this paper's critical path

**Read this first: nothing in this directory was used for any number in the paper.**

Every orbifolder run in this package was made against the **stock**, unmodified geometry
files deposited in `../../orbifolder-1.2.1/Geometry/`. Both census drivers say so in their
own headers, and the reason is worth repeating:

> *"the STOCK `Geometry_Z<..>_1_1.txt`. The massless spectrum and the gauge group are
> independent of declared discrete symmetries, so the census must NOT be run under a
> patched geometry — a patch could otherwise shape its own evidence."*

The patches exist because of a separate finding in this research programme: all 138
geometry files shipped with orbifolder 1.2.1 declare **zero R symmetries**, so the tool's
R-charge selection-rule loop is inert on every stock run. That matters for **coupling**
calculations. This paper computes no couplings; it reads the massless spectrum and the
gauge group, neither of which depends on declared discrete symmetries.

The patches are deposited for completeness, and because a reader auditing the
environment will otherwise wonder whether one was silently applied. The answer is no.

## Contents

| file | what it is |
|---|---|
| `patch_geometry.py`, `patch_geometry_z36.py`, `patch_geometry_z4.py`, `patch_geometry_z44.py`, `patch_geometry_z6i.py` | the scripts that emit patched geometry files, each declaring which R symmetries are exactly expressible as an orbifolder discrete symmetry and which are not |
| `patched_geometry_files/*.txt` | the emitted files, for `Z6-II` (four space groups, two variants each), `Z3xZ6`, `Z4`, `Z4xZ4`, `Z6-I` |

A `_Rsym_prime` variant declares only the R symmetries on planes of **prime** rotation
order, where the law is exactly expressible. A `_Rsym_old` / `_Rsym_all3` variant also
declares the non-prime planes, where orbifolder would apply an older, too-strong rule —
those variants are deposited as the counterfactual, not as the correct patch.

## Licence

The files under `patched_geometry_files/` are **modified GPL v3 works** — orbifolder
geometry files with local changes. They remain under the GNU GPL v3
(`../../orbifolder-1.2.1/COPYING`) and are excluded from the package's CC BY 4.0 licence.
The `.py` patch scripts are the paper author's own work and are MIT-licensed like the rest
of the code.

## One assumption, recorded rather than verified

Both census drivers stage the geometry conditionally:

```bash
[ -s ~/rcharge_build/Geometry/"$(basename "$GEO")" ] || cp -f "$GEO" ~/rcharge_build/Geometry/
```

— the stock file is copied into the build tree only if a file of that name is **not
already there**, and `~/rcharge_build` is the same build tree other phases of this
research programme used with patched geometry files. No check that the build tree held the
stock file at run time was recorded. The exposure is small, for the reason given above,
but it is an assumption and not a verification. See `../ENVIRONMENT.md` §3.

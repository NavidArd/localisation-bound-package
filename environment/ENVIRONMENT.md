# Environment

Two environments are involved, and they do different jobs.

* **orbifolder** (C++, GPL v3, run under WSL2) produced the raw text output. It is needed
  **only** to regenerate the raw output from scratch. The raw output is deposited here,
  gzipped, so nobody needs orbifolder to check any number in the paper.
* **Python** turns that raw output into every number. It needs `numpy` and nothing else.

---

## 1. Python

| | |
|---|---|
| version used for the deposited runs and the reproduction runs | **CPython 3.12.0** (Windows 11, 64-bit) |
| minimum | 3.10 (uses `from __future__ import annotations`, no 3.11+ syntax) |
| third-party packages | **`numpy` only** — version used: **2.2.4** |
| everything else | standard library: `fractions`, `gzip`, `json`, `re`, `os`, `sys`, `math`, `cmath`, `itertools`, `collections`, `multiprocessing` |

Pinned requirements are in `requirements.txt`. `scipy` and `sympy` are *not* used
anywhere in this package.

The arithmetic that matters is exact: budgets, norms, `R`, `|w|^2` and every threshold
are `fractions.Fraction`, never floating point. `numpy` is used for integer/rational
vector bookkeeping over the 16-dimensional lattice, not for numerical linear algebra, so
the results do not depend on the BLAS.

### BLAS threads

Set these **before** `import numpy`, in the parent process, and explicitly (not with
`setdefault`):

```
OMP_NUM_THREADS=1   MKL_NUM_THREADS=1   OPENBLAS_NUM_THREADS=1
```

Every script in this package already does so at module top. On Windows with
`multiprocessing`, omitting them causes worker oversubscription; it does not change any
result, only the wall time.

### Parallelism

The heavy scripts use `multiprocessing.Pool` over the gzipped raw segments, one segment
per task. Six workers were used. The tally is order-independent, so the worker count does
not affect any output.

---

## 2. orbifolder — the exact release, read from the source

**The version string the runs actually carry.** Every deposited raw output file begins
with the tool's own banner:

```
  ###########################################################################
  #  The C++ Orbifolder                                                     #
  #  Version: 1.2.1                                                         #
  #  by H.P. Nilles, S. Ramos-Sanchez, P.K.S. Vaudrevange and A. Wingerter  #
  ###########################################################################
```

That banner is emitted by `src/orbifolder/main.cpp`, lines 212–213, as two literal
`cout` statements. It is the authoritative identification of the build, and it can be
checked in this package by decompressing any file under
`multi_geometry_campaign/data/*/raw/` and reading the first six lines.

> **Recorded discrepancy, so that nobody has to rediscover it.** The distribution
> directory is named `orbifolder-1.2.1` and the binary prints `Version: 1.2.1`, but the
> autoconf metadata inside the same tree was never bumped: `configure.ac` line 1 reads
> `AC_INIT([Orbifolder], [1.0], ...)`, so the generated `configure` carries
> `PACKAGE_VERSION='1.0'` and `PACKAGE_STRING='Orbifolder 1.0'`. A reader who identifies
> the release from `configure` will get `1.0`; the printed banner and the distribution
> name say `1.2.1`. **This package cites 1.2.1, on the authority of the banner in the raw
> output.** `orbifolder-1.2.1/configure.ac` is deposited so the discrepancy is checkable.
> No build hash or `--version` flag output was recorded at run time; that is a gap, and
> the banner is the finest-grained identification available.

**Publication.** H. P. Nilles, S. Ramos-Sánchez, P. K. S. Vaudrevange, A. Wingerter,
*The orbifolder: a tool to study the low energy effective theory of heterotic orbifolds*,
Comput. Phys. Commun. **183** (2012) 1363, arXiv:1110.5229.

**Licence.** GNU GPL v3. The orbifolder **program** is not redistributed here. Its
`COPYING`, `AUTHORS`, `README`, `configure.ac` and `Dockerfile`, and **all 138** geometry
input files (of which the analysis reads seven), are redistributed unmodified under
`orbifolder-1.2.1/`; see the NOTICE there.

> **One dangling constant, so it is not mistaken for a missing input.**
> `multi_geometry_campaign/code/geomlib.py` line 52 defines
> `ALL138 = orbifolder-1.2.1/All138SpaceGroups.txt`. That file is **not** deposited, and
> nothing needs it: the constant is assigned and never read anywhere in this package —
> every geometry is loaded from its own `Geometry/Geometry_<label>.txt`, all 138 of which
> are here. The line is left as it stands because the file is otherwise byte-identical to
> the original.

### Building it

The container the runs were made in, from the unmodified 1.2.1 source tree
(`orbifolder-1.2.1/Dockerfile`, deposited verbatim):

```dockerfile
FROM ubuntu:22.04
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential autoconf automake libtool libgsl-dev libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*
COPY . /orbifolder
WORKDIR /orbifolder
RUN autoreconf -fi || true
RUN ./configure || (cat config.log | tail -50 && false)
RUN make -j$(nproc)
CMD ["bash"]
```

Build it in the root of the unpacked orbifolder 1.2.1 distribution; the binary lands at
`/orbifolder/src/orbifolder/orbifolder`. Requires GSL >= 1.9 and Boost >= 1.0
(`configure.ac` lines 19, 21).

---

## 3. How the raw output was produced

Not Docker: **WSL2 (Ubuntu) on Windows**, with the repository visible at `/mnt/f/physics`.
The drivers are deposited verbatim in `orbifolder_drivers/`.

| | |
|---|---|
| invocation | `wsl -d Ubuntu -- bash <driver> <SEGID> <FIRST> <LAST>` |
| binary | `~/rcharge_build/src/orbifolder/orbifolder` |
| working directory | `~/rcharge_build` |
| repository base inside WSL | `/mnt/f/physics/multi_geometry_campaign` |
| concurrency | 6 segments at a time; each segment is resumable (skipped if its `.txt.gz` already exists) |
| timeout | 43,200 s per segment (census, construction); 20,000 s (gate runs) |
| output | piped to `out.txt`, then `gzip -c` into the data directory; `out.txt` deleted |

### The command sequence, identical in every driver

```
cd <LABEL>
cd spectrum ; print(F) with internal information ; cd ..
cd gauge group ; print gauge group ; print U1 generators ; print simple roots ; print anomaly info ; cd ..
cd vev-config ; analyze config print SU(5) simple roots ; use config(SMConfig) ; print gauge group ; cd ..
cd spectrum ; print summary ; cd ..
cd ..
...
exit
yes
```

The spectrum prints come **before** `analyze config` because that command installs and
selects a new vev-configuration, which would change the `U(1)` basis under which `p_sh`
is reported. This ordering is load-bearing.

### Geometry: STOCK, not patched

Both census drivers state it explicitly in their headers, and it matters:

> *"the STOCK `Geometry_Z<..>_1_1.txt`. The massless spectrum and the gauge group are
> independent of declared discrete symmetries, so the census must NOT be run under a
> patched geometry — a patch could otherwise shape its own evidence."*

**No R-symmetry patch is on Paper A's critical path.** The patch scripts and the patched
geometry files are deposited under `orbifolder_patches/` for completeness and because the
brief asked for them, but they belong to the campaign's *coupling* runs, which Paper A
does not use. See `orbifolder_patches/NOTICE.md`.

### One unverified assumption, recorded

Both drivers stage the geometry file conditionally:

```bash
[ -s ~/rcharge_build/Geometry/"$(basename "$GEO")" ] || cp -f "$GEO" ~/rcharge_build/Geometry/
```

The stock file is copied into the build tree **only if a file of that name is not already
there**. `~/rcharge_build` is the same build tree other campaign phases used with
R-symmetry-patched geometry files. No check that
`~/rcharge_build/Geometry/Geometry_Z3xZ6_1_1.txt` was absent or stock at run time was
recorded. The exposure is small — declared discrete symmetries do not affect the massless
spectrum or the gauge group, which is all this analysis reads — but it is an assumption,
not a verification, and it is stated here rather than left implicit.

---

## 4. Reproducing without orbifolder

Everything in the paper is regenerated from the deposited raw output by Python alone:

```
pip install "numpy>=1.24"
set OMP_NUM_THREADS=1 & set MKL_NUM_THREADS=1 & set OPENBLAS_NUM_THREADS=1
cd multi_geometry_campaign/code/su5_ten_bound
python verify.py Z2xZ4_1_1
```

See the root `README.md` §4 for the entry point behind each claim, and
`regeneration/REPRODUCTION.md` for the runs made when this package was assembled, with
their wall times and outcomes.

---

## 5. Integrity

`CHECKSUMS.sha256` covers every file. `manifest.json` carries the same hashes with sizes
and a role for each file. `.gitattributes` disables end-of-line conversion for all paths,
so `sha256sum -c CHECKSUMS.sha256` verifies on Windows clones as well as Unix ones.

The `Z2xZ4_1_1` catalogue file additionally has its SHA-256 and byte count recorded
**inside the analysis code**, as `z24sw_lib.CATALOGUE_SHA256` and `CATALOGUE_BYTES`
(lines 94–95 of the deposited copy).

**These constants are not enforced at parse time.** They are declared, and they are
checked by `code/z2xz4_sweep/selftest.py`, check `S11` — a script that has to be run
deliberately. A substituted catalogue would therefore *not* fail loudly during an
ordinary analysis run. Run the self-test before trusting a rebuilt or re-extracted
catalogue:

```
cd multi_geometry_campaign/code/z2xz4_sweep
python selftest.py
```

The deposited file matches the recorded values: 7,560,329 bytes, SHA-256
`6c73de64aeddb5688a50789d8afcf58ff4cb1abe598ad12d49a7c4a39a6d2b7d`. The `Z3xZ6_1_1`
catalogue file has no such constant in the code; its hash is
`7fc44dc01d84d0f3cb5476c5288e73dbe74ccb979eaa6725ecfbcd6135ceca23`, recorded in
`CHECKSUMS.sha256` like every other file.

**Rebuilding the integrity files.** `regeneration/build_manifest.py` walks the tree,
hashes every file and writes both `CHECKSUMS.sha256` and `manifest.json`, assigning each
file a role. It is deposited so that amending the package is a scripted step rather than a
hand edit, and so that a referee can confirm the two files are consistent with the tree:

```
python regeneration/build_manifest.py --check     # verify, write nothing
python regeneration/build_manifest.py             # rewrite both
```

`__pycache__/` and `*.pyc` are the only exclusions — byte-compiled Python depends on the
interpreter build and on absolute paths, so hashing it would make the checksums
unverifiable elsewhere. `.gitignore` lists the same exclusion. `manifest.json` cannot
cover its own hash and does not claim to; `CHECKSUMS.sha256` is written first, and the
manifest covers it.

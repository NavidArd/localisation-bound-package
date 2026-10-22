# Redistributed orbifolder files — attribution and licence

The files in this directory are **not the work of the author of the accompanying paper**.
They are redistributed verbatim, unmodified, from the **orbifolder 1.2.1** distribution.

| file | what it is |
|---|---|
| `Geometry/Geometry_Z2xZ4_1_1.txt` | space-group and geometry data for `T^6/(Z2 x Z4)`: point group, twist vectors, space-group generators, lattice, fixed-point data, discrete-symmetry declarations |
| `Geometry/Geometry_Z3xZ6_1_1.txt` | the same for `T^6/(Z3 x Z6)` |
| `Geometry/Geometry_{Z3xZ3,Z4xZ4,Z6-I,Z6-II,Z4}_1_1.txt` | the same for the five further geometries whose twist vectors enter the budgets table |
| `COPYING` | the GNU General Public License, version 3 |
| `AUTHORS` | the orbifolder authors, as shipped |
| `README` | the distribution README, as shipped |
| `configure.ac` | the autoconf input, deposited because it carries the version discrepancy recorded in `../environment/ENVIRONMENT.md` §2 |
| `Dockerfile` | the build recipe used for the runs |

## Authorship

> H. P. Nilles, S. Ramos-Sánchez, P. K. S. Vaudrevange, A. Wingerter,
> *The orbifolder: a tool to study the low energy effective theory of heterotic
> orbifolds*, Comput. Phys. Commun. **183** (2012) 1363,
> [arXiv:1110.5229](https://arxiv.org/abs/1110.5229).

## Licence

**These files are covered by the GNU GPL v3** (`COPYING` in this directory), not by the
CC BY 4.0 / MIT licences that cover the rest of this package. Both `LICENSE-DATA.txt` and
`LICENSE-CODE.txt` at the package root exclude this directory explicitly.

Redistribution here is permitted by the GPL: the files are unmodified, the licence text
travels with them, and authorship is stated above. The orbifolder **program** is not
redistributed — only these input and metadata files. Build the program from the upstream
distribution as described in `../environment/ENVIRONMENT.md` §2.

## Why they are here

With these files, an orbifolder 1.2.1 build, and the drivers in
`../environment/orbifolder_drivers/`, every raw output file in this package can be
reproduced from scratch. The geometry files are additionally read *directly* by the
analysis — the twist vectors of the budgets table are parsed off them, with a per-row
record of the file and line each row came from — so they are an input to the results, not
only to the runs.

**These are the STOCK files.** No patch was applied for any run in this package; see
`../environment/orbifolder_patches/NOTICE.md`.

## Integrity

These files are covered by the package's `CHECKSUMS.sha256`, so their bit-identity to the
upstream distribution can be checked by anyone who has it.

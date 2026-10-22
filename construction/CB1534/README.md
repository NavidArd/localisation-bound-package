# CB1534 — the flagship

An explicitly constructed heterotic `E8 x E8` model on `T^6/(Z3 x Z6)` in which one
twisted fixed point carries **two** states of the same norm-6/5 Standard-Model role. It
is the paper's positive control (b): the sector where the bound *permits* a repeat is not
a vacuous permission.

---

## The embedding

Standard basis; `V1` of order 3, `V2` of order 6, one order-3 Wilson line `W3 = W4`,
all other Wilson lines zero. Exact rationals, also in `embedding_vectors.json`:

```
V1 = ( 1/2,  3/2,  1/2,  1/6,  7/6,  1/2,  1/2,  1/2 |   0,    0,   0, 0, 0, 0, 0, 0 )
V2 = ( 2/3,  2/3, -1/3, -1/3, -1/3,    0,    0,    0 | 1/2,  1/2,   0, 0, 0, 0, 0, 0 )
W3 = ( 2/3, -4/3,  2/3,    0,    0,    0,    0,    0 | 1/3, -1/3, 2/3, 0, 0, 0, 0, 0 )
```

Geometry file: `Geometry_Z3xZ6_1_1.txt`, **stock, unpatched** — deposited at
`../../orbifolder-1.2.1/Geometry/`.

## The model

| | |
|---|---|
| gauge group | `SU(4) x SU(3) x SU(2)` and `E_6 x SU(2)` and `U(1)^3` |
| printed simple roots | 13 |
| printed records | 147 |
| `k_Y` | 5/3 |
| `Tr Y`, `Tr Y^3` | 0, 0 |
| `t_Y . t_anom` | 0 — **vacuous; see the dated note below** |
| net `Q` (triality-aware) | **+3** |
| net `u^c`, `e^c`, `d^c`, `L` (triality-aware) | −3, +1, −3, −8 |

> **Dated note, 2026-10-22.** The row `t_Y . t_anom = 0` is true only vacuously: this model
> has **no anomalous `U(1)`**. Every abelian anomaly trace orbifolder prints for it is zero
> (`raw_output_CB1534.txt`, the block after "U(1)_i - U(1)_i - U(1)_i": all twelve
> `tr Q_i`, `tr Q_i^3`, `tr Q_i Q_j^2` = 0/1, and all fifteen `tr l(rep) Q_j` = 0/1), so
> the anomalous generator is the zero vector and its product with `t_Y` is zero for any
> `t_Y`. The informative statement is the stronger one: every abelian anomaly trace
> vanishes, and since `t_Y = (1/2) Q_1 − (1/6) Q_2` in the printed generators, `Tr Y = Tr
> Y^3 = 0` follows from the printed traces. The row is kept, unchanged, so that the record
> can be audited; it should not be quoted as a Green–Schwarz check. "Triality-aware" means
> that the colour representation of each state is read from its gauge momentum rather than
> from orbifolder's printed Dynkin label.

## The doubled cell

| | |
|---|---|
| field | `F_1` |
| sector | `T(0,1)` |
| fixed point | `(0,0,0,0,0,0)` |
| printed representation | `(1, 3bar, 2, 1, 2)` |
| oscillator | **none** |
| **hidden multiplicity** | **2** |
| role (triality-aware) | `Q` — a genuine `(3, 2)_{+1/6}` |
| `\|p_sh\|^2` | 31/18 |
| `B` | 31/36 |
| `\|w\|^2` | 6/5 |
| `R` | 47/90 |
| `p . Delta` | −1 |
| partner | same coset, same role, `\|q\|^2 = 31/18`, `\|w\|^2 = 6/5` |

Every hypothesis re-checked at this cell: **H1** `Delta` in the lattice with
`|Delta|^2 = 2`; **H2** the Standard Model in the four-dimensional gauge group; **H3**
the sector's unique translation-free constructing element, so `V_g = V2` exactly, with no
Wilson-line contribution and no fixed-point dependence; **H4** equal `|w|^2`; **H5**
distinct gauge momenta.

## The arithmetic

```
threshold   B = 31/36 = 0.8611...  >=  17/20 = 0.85              met, by 1/90
            in R:  R = 47/90 >= 1/2                              met, by 1/45
bound       2*sqrt(47/90) = 1.4453  >=  sqrt(2) = 1.4142         saturated from the ALLOWED side
ceiling     floor(1/(1 - 47/90)) = floor(90/43) = 2              and exactly 2 is realised
```

The `B` margin is `1/90` (`155/180 − 153/180`); `1/45` is the margin in `R`. Earlier text
gave `1/45` for both — corrected in v2.0, no other number affected.

The theorem is tight at this cell in both directions: it forbids the configuration
everywhere below `17/20`, permits it here, and caps it at two.

---

## Independent reproduction

The model was censused once as part of the 2,912-model enumeration, and then **re-run
from its embedding vectors alone** by a second implementation that writes its own
orbifolder deck and shares no code with the enumeration.

| | |
|---|---|
| records identical, fresh vs stored | **147 / 147** |
| gauge group, `U(1)` generators, simple roots, anomaly traces | identical |
| model section length | 84,019 bytes |
| SHA-256 of the model section, from **either** source | `5c284fffda999f46…` |

The two sources are a 8.6 MB segment file and a 425 kB gate file, and the section sits at
a different offset in each; the byte-identity is a reproduction result, not an artefact of
reading one file twice.

Three further models were run in the same fresh deck and also reproduced exactly:
**CB1370** 154 / 154, **CB2307** 219 / 219, **CB1041** 214 / 214.

---

## Files

| file | what it is |
|---|---|
| `phase_deck_s5_models_block.txt` | the model block from the enumeration's own deck `s5_models.txt` |
| `gate_deck_gate_models_block.txt` | the same model as written independently from the embedding vectors by `g1b_fresh_run.py` |
| `embedding_vectors.json` | `V1`, `V2`, `W3` as exact rationals |
| `raw_output_CB1534.txt` | the complete orbifolder print-out, cut from `cb_seg15.txt.gz` |
| `raw_output_CB1534_FRESH_gate_run.txt` | the same, cut from the independent run `gate_v1.txt.gz` |

The complete gzipped sources are in
`../../multi_geometry_campaign/data/su5_construction/raw/` and
`../../multi_geometry_campaign/data/gate_su5_construction/`.

## To re-run it

```
wsl -d Ubuntu -- bash environment/orbifolder_drivers/run_gate_orb.sh v1 gate_models.txt
```

with an orbifolder 1.2.1 build; see `../../environment/ENVIRONMENT.md` §3. The deck
`gate_models.txt` (four models, including this one) is deposited in the workspace.

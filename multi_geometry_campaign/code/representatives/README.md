# Lattice representatives of the embeddings (package v2.2)

The massless spectrum of a `T^6/(Z3 x Z6)` model depends on the **vectors** `V1`, `V2`,
`W3 = W4`, not only on their classes modulo the `E8 x E8` lattice `Lambda`: shifting a
gauge shift or a Wilson line by a lattice vector can act as discrete torsion
(Ploger, Ramos-Sanchez, Ratz, Vaudrevange, *Mirage Torsion*, JHEP 04 (2007) 063,
hep-th/0702176, sections 3.1, 3.3, 4.1). The gauge group and every class-level condition of
section 7.2 of the paper are unchanged by such a shift; the twisted spectrum is not.

**Convention, stated once.** A *pair of classes* is `(V1 + Lambda, W3 + Lambda)` over a fixed
`V2` of eq. (7.4). An *embedding* is a pair of **vectors** `(V1, W3)`, `W4 = W3`, in those
classes that passes the complete orbifolder modular-invariance set (`sc_lib.mi_report`).
On this geometry the projection phases of two embeddings over one pair of classes differ by
amounts fixed by `(a, b) = 3 (lambda.V2, mu.V2) mod 3`, where `V1 -> V1 + lambda`,
`W3 -> W3 + mu`; modular invariance (MI3, MI5) allows all nine values. `(a, b) = (0, 0)` is
the representative from which the 4,032 Standard-Model-intact embeddings of package v2.1
were built (`su5_construction/s5_models.txt`, `gate_su5_construction/caseB_models.txt`).

| script | what it does | output (in `../../data/representatives/`) |
|---|---|---|
| `b1_brothers.py` | 18 lattice shifts of `V1` and of `W3` of the model (7.1) (CB1534), with their `lambda.V2`, `lambda.W3`, `mu.V2`, `mu.V1` mod 1, all MI-valid | `b1_decks.txt`, `b1_decks.json` |
| `b2_run_orbifolder.sh` | orbifolder 1.2.1 on a deck file, stock geometry | `b2_out.txt`, `b11_out.txt` |
| `b3_compare.py` | record-by-record comparison of the 19 spectra with CB1534 | `b3_compare.json` |
| `b4_triality_nets.py` | the net indices `I_X = sum_M d_hidden(M) [1(M ~ X) - 1(M ~ Xbar)]`, colour read from the momentum (library + script) | `b4_nets.json` |
| `b11_v2_brothers.py` | lattice shifts of `V2` (all nine MI-valid `(nu.V1, nu.W3)` mod 1) and the nine `(a, b)` reference decks | `b11_decks.txt`, `b11_decks.json` |
| `b11b_compare.py` | each `V2` shift reproduces one of the nine `(a, b)` spectra | `b11_compare.json` |
| `b5_build_brother_census.py` + `b5b_missing.py` | all nine `(a, b)` embeddings over each of the 4,032 SM-intact pairs: 36,288 MI-valid decks. `b5` fills 35,976 slots from the vectors `+-e_i +- e_j`; the other 312 need the wider search of `b5b` — 232 `(a, 0)` slots by a shift `lambda` of `V1`, 80 `(0, b)` slots by a shift `mu` of `W3` (labelled `Y...`, segment 902). `b5b` re-computes each deck's `(a, b)` from its vectors and asserts the label. Repaired 2026-10-22: as first deposited `b5b` never searched `mu`, and the 80 `(0, b)` slots held copies of the `(0,0)` deck (`brothers/superseded_v2.2/NOTICE.md`) | `brothers/decks_NNN.txt`, `brothers/build*.json` |
| `b6_run_brothers.sh` | orbifolder 1.2.1 on all 183 segments (6 workers; segment 902, the 80 repaired `(0, b)` decks, was run by the v10 verification with the same binary and commands and is deposited verbatim) | `brothers/out_NNN.txt.gz`, `brothers/run.log` |
| `b10_brother_analysis.py` | per deck: its slot `(a, b)` re-computed from its vectors against its label (`label_check`); net indices, the doubled 10-like test, the gauge group; aggregates; the `(0,0)` control against the certified per-model nets of `g5_battery.json` / `g8_batteryB.json` | `b10_brothers.json`, `b10_per_deck_nets.csv.gz` |
| `b8_allbins.py` | the second code for the orbifold projection (`root_lift/rl_d_projection.py`, imported unchanged) run at **every** (twisted constructing element, massless right-mover) bin of `Geometry_Z3xZ6_1_1`, both chiralities, not only the bins the spectrum populates | `b8_allbins.json` (BR0 = CB1534, BR4, BR7), `b8_allbins_deposited.json` (CB1534, CB1041 from the deposited deck and output) |
| `b7_census_rank.py` | on both catalogue censuses (all 11,526 + 4,612 models: each worker applies the `E_6` banner parse of `uniform_census_e6.py`, as the paper's tables do; repaired 2026-10-22, before which 16,134 of the 16,138 models were covered): the mean-budget ceiling `n_states <= 1/(1 - Rbar)` on every group, and the rank restriction (largest `A_{n-1}` in the roots of one `E8` factor orthogonal to `W`, and inside the local group) | `b7_Z2xZ4_1_1.json`, `b7_Z3xZ6_1_1.json` |

Pristine copies of every output are in `results_deposited/representatives/`; `claims.json`
rows `REP-*` resolve against them.

**Results.** 36,288 / 36,288 spectra parsed; 36,288 / 36,288 decks realise the `(a, b)` their label names; the `(0,0)` nets agree with the certified
per-model nets on 4,032 / 4,032; none of the 36,288 has net `Q = u^c = e^c = +-3` at once;
216 reach two of the three; 21,168 carry a doubled 10-like record at `T(0,1)`/`T(0,5)`
(2,174 of them at `(0,0)`, the v2.1 count `CON-24`, now superseded by `REP-05`); 32,052 of the
32,256 non-`(0,0)` spectra differ from their pair's `(0,0)` spectrum in left-chiral record count
or nets (`REP-06`; the other 204 agree with it on those two coarse measures). On the model
(7.1) the nine `(a, b)` spectra are distinct (147, 153, 119, 136, 142, 132, 154, 164 and 132
multiplets, the two 132s with different nets); `lambda.W3` and `mu.V1` change nothing, and
every `V2` shift reproduces one of the nine. The projection code reproduces the
oscillator-free twisted spectrum at all 140 bins of CB1534 and of CB1041 (0 mismatches;
every empty bin empty). The mean-budget ceiling has 0 violations on both catalogues; the
rank restriction lowers a row's ceiling at one realised value only, `R = 83/90`
(`floor(1/(1-R)) = 12 -> 9`).

**What is not proved.** That the nine `(a, b)` choices exhaust the zero-torsion
representatives follows from the structure of the projection phase on this geometry (the
paper, section 7.2) and is confirmed on the model (7.1) by every shift tried; it is not
re-checked by brute force on all 4,032 pairs.

**Runtime.** The orbifolder steps (`b2`, `b6`) need an orbifolder 1.2.1 build
(`environment/`); `b6` took about 3 CPU-hours, half an hour on six workers
(`brothers/run.log`). Every Python step reads the deposited
raw output and runs without orbifolder; `b7` re-reads both complete catalogue censuses.

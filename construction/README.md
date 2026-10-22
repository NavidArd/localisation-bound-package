# The construction

The bound permits a second norm-6/5 role only where the sector budget reaches `17/20`.
Across the seven geometries examined, exactly two sectors do: `T(0,1)` and `T(0,5)` of
`T^6/(Z3 x Z6)`, at `B = 31/36`. Nothing in the catalogues realises the repeat there, so
the permission would be vacuous unless a model realising it can be exhibited.

This directory contains that model, the enumeration that found it, and the enumeration's
own boundary — the proof that exactly two classes of embedding can admit the
configuration at all, so that the search is exhaustive rather than merely large.

---

## 1. Layout

```
CB1534/                  the flagship: input deck, embedding vectors, raw output,
                         and the same output re-derived by an independent deck
CB1041_DEMOTED/          the earlier flagship, kept with the reason it was demoted
README.md                this file
```

The enumeration code and its complete output are in the workspace, at the paths the
scripts expect:

```
multi_geometry_campaign/code/su5_construction/        the enumeration (class A)
multi_geometry_campaign/code/gate_su5_construction/   the second implementation, and class B
multi_geometry_campaign/data/su5_construction/        s5_models.txt (the 2,912-model deck),
                                                      raw/cb_seg0..29.txt.gz, s*.json
multi_geometry_campaign/data/gate_su5_construction/   caseB_models.txt (1,432), gate_models.txt,
                                                      gate_caseB.txt.gz, gate_v1.txt.gz,
                                                      gate_equiv.txt.gz, g*.json
results_deposited/                                    pristine copies of every JSON above
```

---

## 2. The two-class enumeration, and why it is exhaustive

The configuration needs a local shift `V_g` of order 6 at `B = 31/36` carrying a root
`Delta` orthogonal to the Standard-Model Cartan span with `p . Delta = 1`. Enumerating
the Diophantine solutions over 8 admissible budget triples, 10 weights and 260 candidate
roots gives **3,000 solutions**, falling into **119 distinct lattice classes**, which
reduce under the lattice automorphism group to **exactly two inequivalent classes**.

That "two" is the load-bearing number: it is what turns a search into an exhaustion.

| | class A | class B |
|---|---:|---:|
| `(V1, W)` pairs passing the screens | 750,384 | 696,600 |
| inequivalent pairs | **2,912** | **1,432** |
| literal embeddings built, all passing MI0–MI6 | 2,912 | 1,432 |
| build failures | 0 | 0 |
| censused with orbifolder | 2,912 | 1,432 |
| fully Standard-Model-valid | 2,912 | **1,120** |
| excluded (SM roots merged into a larger factor) | — | 312 |

**4,032 = 2,912 + 1,120** Standard-Model-valid models over the complete loophole space.

Of those, **2,174 carry a doubled norm-6/5 role**, and **zero** reach three net
generations in all three norm-6/5 roles simultaneously; the best is two of three.

> **Corrected in v2.0.** Earlier text put this at **1,482**. Recounted from the deposited
> tables it is **2,174** = 1,812 (class A, all of which are Standard-Model-valid) + 362
> (class B, doubled **and** Standard-Model-valid). The class-B figure that appears in
> `g8_batteryB.json::summary.doubled_10role`, **400**, is over all 1,432 *built* lifts and
> includes 38 that are not Standard-Model models — a different unit, and not the one the
> 4,032 headline is stated over. Over all 4,344 built models the figure is 2,212.
> `g8c_doubled_smvalid.py` derives all three; claims `CON-24`, `CON-25`, `CON-26`.
> **The negative result is untouched:** zero in all three roles either way.

> **Why the co-occurrence fails is UNDERIVED.** The enumeration measures that the
> doubling and the third generation never co-occur; it does not explain it. This is
> stated as an open question in the source documents and must be stated as one in the
> paper.

Entry points: `g2_classes.py` (the class boundary), `g4_lift.py` (both lift spaces),
`s5_build.py` (class A embeddings), `g7_buildB.py` (class B embeddings),
`s11_battery2.py` and `g5_battery.py` (the funnels).

---

## 3. Two reproduction gaps, reported not repaired

**R-1 — `g8b_smvalid.json` has no producing script.** The class-B split into 1,120
Standard-Model-valid and 312 excluded lifts is on disk, and the **4,032** headline
depends on it, but no script in `code/gate_su5_construction/` writes that file; it appears
to have been an inline run. Its inputs (`g8_batteryB.json`, `g4_lift_B.json`) are
deposited, so the split can be recomputed, but this package cannot claim it *was*
recomputed from a deposited script. **Anything Paper A prints that rests on 1,120 or
4,032 carries this gap.**

**R-2 — `g7_buildB.py` sets no random seed.** `s5_build.py` fixes
`random.seed(20260820)`; the class-B builder does not. Class-B representative selection is
therefore not bit-for-bit reproducible. The counts (1,432 built, 0 failed) are
representative-independent; the specific vectors in `caseB_models.txt` are not, so that
file is deposited as the record of the run rather than as something a rerun will
reproduce.

A third item, not a gap but a caution: **`s8_battery.py` is flagged do-not-reuse.** It
carries an unweighted hypercharge trace (`trY_nonzero = 2,728`, superseded by 2,431 from
`s11_battery2.json`) and the triality-blind role table described below. Its
`n_with_T0l_ten = 1,812` is unaffected and is reproduced by the second implementation.

---

## 4. The triality correction

`s8_battery.py::su3_dim` returns a `(dimension, sign)` pair and **every caller discards
the sign**; `role_of` accepts a colour-conjugation argument it never uses. The role table
therefore labels any `(3 or 3bar, 2)_{+1/6}` state a quark doublet, so
`(3bar, 2)_{+1/6}` and `(3, 2)_{−1/6}` exotics were counted as quark doublets.

Triality is not a convention. The only legitimate relabelling of a spectrum is the overall
conjugation `p -> −p`, which flips triality **and** hypercharge together; flipping
hypercharge alone on a triality-blind table is not a relabelling of the same spectrum.

What moves and what does not:

| quantity | triality-blind (superseded) | triality-aware (correct) |
|---|---:|---:|
| class-A models with net `Q = ±3` | 326 | **339** |
| models whose verdict changes | — | 275 |
| models passing `Tr Y = 0` and net `Q = 3` and doubled | 28 | **37** (only 14 shared) |
| class-A models with a doubled norm-6/5 role | 1,812 | **1,812** (unchanged) |
| class-A models with `Tr Y = 0` | 481 | **481** (unchanged) |
| three net generations in all three roles | 0 | **0** (unchanged) |

**The negative result is unchanged.** What changed is which individual model is the
flagship — see `CB1041_DEMOTED/`.

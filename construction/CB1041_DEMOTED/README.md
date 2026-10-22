# CB1041 — DEMOTED. Do not print this model.

CB1041 was the enumeration's original flagship. It is kept here in full because a
demotion that leaves no trace is not auditable, and because its record-level data is
sound and was independently confirmed.

**Paper A should print CB1534 only.** See `../CB1534/`.

---

## Why it was demoted

CB1041 was selected as a model with a doubled norm-6/5 role *and* three net quark
doublets. The second half is wrong.

`code/su5_construction/s8_battery.py::su3_dim` returns a `(dimension, sign)` pair, and
every caller discards the sign:

```python
d3, _ = su3_dim(...)
```

`role_of(c3, c2, Y, sY, sC)` additionally accepts a colour-conjugation argument `sC`
that it never uses. The role table therefore assigns `Q` to any `(3 or 3bar, 2)_{+1/6}`
state, so `(3bar, 2)_{+1/6}` and `(3, 2)_{−1/6}` **exotics were counted as quark
doublets**.

Under the triality-aware recount CB1041 has:

| | triality-blind (as published) | triality-aware (correct) |
|---|---:|---:|
| net `Q` | +3 | **−4** |
| net `u^c` | −2 | +2 |
| net `e^c` | −2 | +2 |

`|net Q| = 4`, not 3. **CB1041 is not a three-generation model.**

There is a second, independent sign that the original dossier was inconsistent: it
recorded `net Q = +3` on the `s_Y = −1` branch while simultaneously labelling its own
field `F_6` as `u^c`, which is the `s_Y = +1` label. The two halves of the dossier
disagree with each other.

`s11_battery2.py` imports `role_of` from `s8_battery.py` and inherits the same defect.

This is recorded in the campaign's error register as an instance of its most repeated
failure family — a conjugate object silently substituted for the object itself.

---

## What survives

Everything at record level. CB1041's spectrum is real and reproduces exactly:

| | |
|---|---|
| gauge group | `SU(3) x SU(2)` and `SO(8) x SU(2) x SU(2)` and `U(1)^7` |
| four-dimensional roots | 36 |
| printed simple roots | 9 |
| `k_Y` | 5/3 |
| `Tr Y`, `Tr Y^3` | 0, 0 |
| `t_Y . t_anom` | 0 |
| printed records | 214 |
| representation-dimension mismatches vs the lattice | 0 |
| level-1 unitarity violations | 0 |
| **records identical, fresh run vs stored run** | **214 / 214** |

Its doubled cells are genuine and satisfy the theorem from the allowed side:

| | |
|---|---|
| `F_6` | `T(0,1)`, `(0,0,0,0,0,0)`, rep `(3,1,1,2,1)`, role `u^c` |
| `F_57` | `T(0,5)`, rep `(1,1,1,2,1)`, role `e^c`-conjugate |
| oscillator | none, both |
| hidden multiplicity | 2, both |
| `\|p_sh\|^2`, `B`, `\|w\|^2`, `R` | 31/18, 31/36, 6/5, 47/90 |
| `Delta` | `(0^8 \| 1,1,0^6)`, `\|Delta\|^2 = 2` |
| `p . Delta` | −1 |
| bound | `2*sqrt(47/90) = 1.4453 >= sqrt(2) = 1.4142` — saturated from the allowed side |

So CB1041 remains a valid **existence** witness for the doubling. It is not a valid
three-generation model, and it was the three-generation claim that made it the flagship.

---

## A recorded internal inconsistency in the source documents

The gate document's printed modular-invariance value tuple for CB1041
(`MI1 = 4, MI2 = 10, MI3 = 2, MI4 = 6, MI5 = 0, MI6 = 12`) does **not** match its own
JSON, `g1_cb1041.json::A_MI`, which records
`MI1 = 10, MI2 = 10, MI3 = 0, MI4 = 0, MI5 = 0, MI6_33 = 10, MI6_34 = 10`.

The PASS verdict is identical in both, and no downstream number depends on the tuple.
Reported here; not resolved, and the source documents were not edited.

---

## Files

| file | what it is |
|---|---|
| `phase_deck_s5_models_block.txt` | the model block as it appears in the enumeration's deck `s5_models.txt` |
| `gate_deck_gate_models_block.txt` | the same model as written independently by the second implementation from the embedding vectors alone |
| `embedding_vectors.json` | `V1`, `V2`, `W3` as exact rationals |
| `raw_output_CB1041.txt` | the complete orbifolder print-out, cut from `cb_seg10.txt.gz` |
| `raw_output_CB1041_FRESH_gate_run.txt` | the same, cut from the independent run `gate_v1.txt.gz` |

The two raw extracts are **byte-identical**. The model section is 129,908 bytes and
hashes to SHA-256 `d286160656d35740…` when cut from either source — from
`cb_seg10.txt.gz` (an 8.6 MB file, at one offset) and from `gate_v1.txt.gz` (a 425 kB
file, at another). The complete gzipped sources they were cut from are deposited in the
workspace, so neither extract is the sole evidence for anything.

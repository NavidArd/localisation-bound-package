# Claim-to-script map

Every number Paper A prints, the script that regenerates it, the deposited file and
key that store it, and whether the stored value matches.

**This file is generated from [`claims.json`](claims.json) by
`regeneration/claims_md.py`, and the `value` column is read out of the deposited
files by `regeneration/claims_check.py` rather than transcribed.** `status` compares
that value with the one the paper will print. Rows that are a sum or another
reduction across the two catalogues show their component keys and the reduction, and
it is recomputed when this file is built.

| | |
|---|---:|
| claims mapped | **250** |
| resolved and matching (`VERIFIED`) | **250** |
| mismatched or unresolved | **0** |
| of the matching rows, superseded and **not** printed by the paper | **30** |

**All 250 rows resolved against the deposited files and matched.**

A row marked `SUPERSEDED` is still true of the file it came from -- superseded
numbers are kept and labelled rather than deleted, so that a correction can be
audited -- but it is not the number the paper prints. The row that replaces it is
named in the status column.

File paths below are relative to `results_deposited/`; the same tables live in the
workspace at `multi_geometry_campaign/data/`, where re-runs overwrite them.
Scripts are relative to `multi_geometry_campaign/code/`, except the two-class
re-enumeration, which lives in `two_classes/`.

---

## Scope of the verification

| id | claim | value | unit | script | file :: key | status |
|---|---|---:|---|---|---|---|
| `A-01` | Census runs, both catalogues. NOT distinct models: ten Z3xZ6_1_1 models are censused twice, so these 15,545 runs cover 15,535 distinct models (DEN-07). | **15,545** | models | `su5_ten_bound/verify.py` | `su5_ten_bound/verify_Z2xZ4_1_1.json::models` + `su5_ten_bound/verify_Z3xZ6_1_1.json::models` | VERIFIED (SUPERSEDED by UNI-01) |
| `A-02` | MODEL-CELLS, twisted -- the 128 twisted constructing elements counted once per model. NOT a count of fixed points: the two geometries have 58 + 70 = 128 twisted cells in total. | **851,231** | model-cells | `su5_ten_bound/verify.py` | `su5_ten_bound/verify_Z2xZ4_1_1.json::cells` + `su5_ten_bound/verify_Z3xZ6_1_1.json::cells` | VERIFIED (SUPERSEDED by UNI-02) |
| `A-03` | Standard-Model role-carrying fields. | **481,728** | fields | `su5_ten_bound/verify.py` | `su5_ten_bound/verify_Z2xZ4_1_1.json::fields` + `su5_ten_bound/verify_Z3xZ6_1_1.json::fields` | VERIFIED (SUPERSEDED by UNI-04) |
| `A-04` | Left-chiral records parsed AT TWISTED CELLS. The untwisted sector is dropped before this count; including it gives 3,298,075. | **3,112,629** | records | `gate_su5_bound/d2_pshdup.py` | `gate_su5_bound/d2_pshdup_Z2xZ4_1_1.json::recs` + `gate_su5_bound/d2_pshdup_Z3xZ6_1_1.json::recs` | VERIFIED |

## The projection norms

| id | claim | value | unit | script | file :: key | status |
|---|---|---:|---|---|---|---|
| `A-05` | Fields with \|w\|^2 = 6/5 (norm-6/5 roles). THE 104,889. <br/>*Also equals the number of (cell, role) pairs, because that multiplicity is 1 everywhere. A FIELD count, not a fixed-point count.* | **104,889** | fields | `su5_ten_bound/verify.py` | `su5_ten_bound/verify_Z2xZ4_1_1.json::w2_by_fam.10.6/5` + `su5_ten_bound/verify_Z3xZ6_1_1.json::w2_by_fam.10.6/5` | VERIFIED (SUPERSEDED by UNI-06) |
| `A-06` | Fields with \|w\|^2 = 4/5 (norm-4/5 roles). | **376,839** | fields | `su5_ten_bound/verify.py` | `su5_ten_bound/verify_Z2xZ4_1_1.json::w2_by_fam.5bar.4/5` + `su5_ten_bound/verify_Z3xZ6_1_1.json::w2_by_fam.5bar.4/5` | VERIFIED (SUPERSEDED by UNI-08) |
| `A-07` | \|w\|^2 exceptions: no field takes any value other than 6/5 or 4/5. | **0** | exceptions | `su5_ten_bound/verify.py` | `su5_ten_bound/verify_Z2xZ4_1_1.json::smbasis_fail` | VERIFIED |
| `A-08` | k_Y = 5/3 measured in every censused Z2xZ4_1_1 model. | **10,925** | models | `su5_ten_bound/verify.py` | `su5_ten_bound/verify_Z2xZ4_1_1.json::kY.5/3` | VERIFIED (SUPERSEDED by UNI-43) |
| `A-09` | k_Y = 5/3 measured in every censused Z3xZ6_1_1 model. | **4,620** | models | `su5_ten_bound/verify.py` | `su5_ten_bound/verify_Z3xZ6_1_1.json::kY.5/3` | VERIFIED (SUPERSEDED by UNI-39) |
| `A-10` | Invariant I1 (\|p_sh\|^2/2 = 1 - dc - N_osc) failures. | **0** | fields | `su5_ten_bound/verify.py` | `su5_ten_bound/verify_Z2xZ4_1_1.json::I1_bad` | VERIFIED |
| `A-11` | Invariant I1 checked on every Z2xZ4_1_1 field. | **312,050** | fields | `su5_ten_bound/verify.py` | `su5_ten_bound/verify_Z2xZ4_1_1.json::I1_checked` | VERIFIED (SUPERSEDED by UNI-11) |

## The theorem, tested on the full catalogues

| id | claim | value | unit | script | file :: key | status |
|---|---|---:|---|---|---|---|
| `A-12` | Theorem violations, Z2xZ4_1_1. | **0** | pairs | `su5_ten_bound/verify.py` | `su5_ten_bound/verify_Z2xZ4_1_1.json::viol_n` | VERIFIED |
| `A-13` | Theorem violations, Z3xZ6_1_1. | **0** | pairs | `su5_ten_bound/verify.py` | `su5_ten_bound/verify_Z3xZ6_1_1.json::viol_n` | VERIFIED |
| `A-14` | Per-(cell, role) multiplicity of every norm-6/5 role is 1, Z2xZ4_1_1. | **67,555** | (cell,role) pairs | `su5_ten_bound/verify.py` | `su5_ten_bound/verify_Z2xZ4_1_1.json::mult_by_fam.10.1` | VERIFIED (SUPERSEDED by UNI-33) |
| `A-15` | Per-(cell, role) multiplicity of every norm-6/5 role is 1, Z3xZ6_1_1. | **37,334** | (cell,role) pairs | `su5_ten_bound/verify.py` | `su5_ten_bound/verify_Z3xZ6_1_1.json::mult_by_fam.10.1` | VERIFIED (SUPERSEDED by UNI-34) |
| `A-16` | Cells where the bound does NOT protect a norm-6/5 role, Z2xZ4_1_1. | **0** | cells | `su5_ten_bound/verify.py` | `su5_ten_bound/verify_Z2xZ4_1_1.json::loophole_n` | VERIFIED |
| `A-17` | Cells where the bound does NOT protect a norm-6/5 role, Z3xZ6_1_1. <br/>*All at R = 47/90, B = 31/36, in sectors T(0,1)/T(0,5).* | **65** | fields | `su5_ten_bound/verify.py` | `su5_ten_bound/verify_Z3xZ6_1_1.json::R_by_fam.10.47/90` | VERIFIED |

## The mechanism, tested step by step

| id | claim | value | unit | script | file :: key | status |
|---|---|---:|---|---|---|---|
| `A-18` | Same-role doubled pairs tested by the proof. THE 31,822. This sums TWO role conventions -- its Z2xZ4_1_1 half counts hidden non-singlets, its Z3xZ6_1_1 half does not. The homogeneous number is DEN-06 = 32,294, and that is what the paper prints. | **31,822** | pairs | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json::pairs` + `gate_su5_bound/d1_mechanism_Z3xZ6_1_1.json::pairs` | VERIFIED (SUPERSEDED by DEN-06) |
| `A-19` | Pairs whose difference Delta is exactly a root (\|Delta\|^2 = 2), Z2xZ4. | **24,406** | pairs | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json::M4_norm.2` | VERIFIED |
| `A-20` | Pairs whose difference Delta is exactly a root (\|Delta\|^2 = 2), Z3xZ6. | **7,416** | pairs | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z3xZ6_1_1.json::M4_norm.2` | VERIFIED |
| `A-M1_fail` | Proof step M1 (equal \|w\|^2): failures on Z2xZ4. | **0** | pairs | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json::M1_fail` | VERIFIED |
| `A-M2_fail` | Proof step M2 (aligning sigma exists): failures on Z2xZ4. | **0** | pairs | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json::M2_fail` | VERIFIED |
| `A-M3_zero` | Proof step M3 (Delta non-zero): failures on Z2xZ4. | **0** | pairs | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json::M3_zero` | VERIFIED |
| `A-M3_notinlam` | Proof step M3 (Delta in Lambda): failures on Z2xZ4. | **0** | pairs | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json::M3_notinlam` | VERIFIED |
| `A-M3_notperp` | Proof step M3 (P_W Delta = 0): failures on Z2xZ4. | **0** | pairs | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json::M3_notperp` | VERIFIED |
| `A-M3_odd` | Proof step M3 (\|Delta\|^2 even): failures on Z2xZ4. | **0** | pairs | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json::M3_odd` | VERIFIED |
| `A-M3_cs_fail` | Proof step M3 (Cauchy-Schwarz): failures on Z2xZ4. | **0** | pairs | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json::M3_cs_fail` | VERIFIED |
| `A-M4_notroot` | Proof step M4 (Delta is a root): failures on Z2xZ4. | **0** | pairs | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json::M4_notroot` | VERIFIED |
| `A-gram_fail` | Proof step M6 (Gram ceiling): failures on Z2xZ4. | **0** | pairs | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json::gram_fail` | VERIFIED |
| `A-21` | Norm-6/5 fields carrying hidden multiplicity > 1, Z2xZ4_1_1. | **0** | fields | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json::hid10_n` | VERIFIED (SUPERSEDED by UNI-35) |
| `A-22` | Norm-6/5 fields carrying hidden multiplicity > 1, Z3xZ6_1_1. | **0** | fields | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z3xZ6_1_1.json::hid10_n` | VERIFIED (SUPERSEDED by UNI-36) |

## The Gram ceiling against measurement

| id | claim | value | unit | script | file :: key | status |
|---|---|---:|---|---|---|---|
| `GRAM-Z3xZ6-83_90` | Gram ceiling at R = 83/90 (Z3xZ6): floor(1/(1-R)) = 12, measured maximum same-(cell,role) multiplicity = 2 [slack]. | **2** | (cell,role) multiplicity | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z3xZ6_1_1.json::multR_rec.5bar|83/90` *(max_int_key)* | VERIFIED |
| `GRAM-Z2xZ4-33_40` | Gram ceiling at R = 33/40 (Z2xZ4): floor(1/(1-R)) = 5, measured maximum same-(cell,role) multiplicity = 4 [one short]. | **4** | (cell,role) multiplicity | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json::multR_rec.5bar|33/40` *(max_int_key)* | VERIFIED |
| `GRAM-Z3xZ6-34_45` | Gram ceiling at R = 34/45 (Z3xZ6): floor(1/(1-R)) = 4, measured maximum same-(cell,role) multiplicity = 4 [SATURATED]. | **4** | (cell,role) multiplicity | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z3xZ6_1_1.json::multR_rec.5bar|34/45` *(max_int_key)* | VERIFIED |
| `GRAM-Z2xZ4-7_10` | Gram ceiling at R = 7/10 (Z2xZ4): floor(1/(1-R)) = 3, measured maximum same-(cell,role) multiplicity = 3 [SATURATED]. | **3** | (cell,role) multiplicity | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json::multR_rec.5bar|7/10` *(max_int_key)* | VERIFIED |
| `GRAM-Z3xZ6-53_90` | Gram ceiling at R = 53/90 (Z3xZ6): floor(1/(1-R)) = 2, measured maximum same-(cell,role) multiplicity = 2 [SATURATED]. | **2** | (cell,role) multiplicity | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z3xZ6_1_1.json::multR_rec.5bar|53/90` *(max_int_key)* | VERIFIED |
| `GRAM-Z2xZ4-23_40` | Gram ceiling at R = 23/40 (Z2xZ4): floor(1/(1-R)) = 2, measured maximum same-(cell,role) multiplicity = 2 [SATURATED]. | **2** | (cell,role) multiplicity | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json::multR_rec.5bar|23/40` *(max_int_key)* | VERIFIED |
| `GRAM-Z3xZ6-8_15` | Gram ceiling at R = 8/15 (Z3xZ6): floor(1/(1-R)) = 2, measured maximum same-(cell,role) multiplicity = 1 [slack]. | **1** | (cell,role) multiplicity | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z3xZ6_1_1.json::multR_rec.5bar|8/15` *(max_int_key)* | VERIFIED |
| `GRAM-Z3xZ6-47_90` | Gram ceiling at R = 47/90 (Z3xZ6): floor(1/(1-R)) = 2, measured maximum same-(cell,role) multiplicity = 1 [1 in catalogue; 2 constructed]. | **1** | (cell,role) multiplicity | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z3xZ6_1_1.json::multR_rec.10|47/90` *(max_int_key)* | VERIFIED |
| `GRAM-Z2xZ4-17_40` | Gram ceiling at R = 17/40 (Z2xZ4): floor(1/(1-R)) = 1, measured maximum same-(cell,role) multiplicity = 1 [tight]. | **1** | (cell,role) multiplicity | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json::multR_rec.10|17/40` *(max_int_key)* | VERIFIED |
| `GRAM-Z3xZ6-16_45` | Gram ceiling at R = 16/45 (Z3xZ6): floor(1/(1-R)) = 1, measured maximum same-(cell,role) multiplicity = 1 [tight]. | **1** | (cell,role) multiplicity | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z3xZ6_1_1.json::multR_rec.10|16/45` *(max_int_key)* | VERIFIED |
| `GRAM-Z2xZ4-13_40` | Gram ceiling at R = 13/40 (Z2xZ4): floor(1/(1-R)) = 1, measured maximum same-(cell,role) multiplicity = 1 [tight]. | **1** | (cell,role) multiplicity | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json::multR_rec.5bar|13/40` *(max_int_key)* | VERIFIED |
| `GRAM-Z2xZ4-3_40` | Gram ceiling at R = 3/40 (Z2xZ4): floor(1/(1-R)) = 1, measured maximum same-(cell,role) multiplicity = 1 [tight]. | **1** | (cell,role) multiplicity | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json::multR_rec.5bar|3/40` *(max_int_key)* | VERIFIED |
| `GRAM-Z3xZ6-11_30` | Gram ceiling at R = 11/30 (Z3xZ6): floor(1/(1-R)) = 1, measured maximum same-(cell,role) multiplicity = 1 [tight]. | **1** | (cell,role) multiplicity | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z3xZ6_1_1.json::multR_rec.5bar|11/30` *(max_int_key)* | VERIFIED |
| `GRAM-Z3xZ6-23_90` | Gram ceiling at R = 23/90 (Z3xZ6): floor(1/(1-R)) = 1, measured maximum same-(cell,role) multiplicity = 1 [tight]. | **1** | (cell,role) multiplicity | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z3xZ6_1_1.json::multR_rec.5bar|23/90` *(max_int_key)* | VERIFIED |
| `GRAM-Z3xZ6-4_45` | Gram ceiling at R = 4/45 (Z3xZ6): floor(1/(1-R)) = 1, measured maximum same-(cell,role) multiplicity = 1 [tight]. | **1** | (cell,role) multiplicity | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z3xZ6_1_1.json::multR_rec.5bar|4/45` *(max_int_key)* | VERIFIED |
| `GRAM-Z3xZ6-17_90` | Gram ceiling at R = 17/90 (Z3xZ6): floor(1/(1-R)) = 1, measured maximum same-(cell,role) multiplicity = 1 [tight]. | **1** | (cell,role) multiplicity | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z3xZ6_1_1.json::multR_rec.10|17/90` *(max_int_key)* | VERIFIED |
| `GRAM-Z2xZ4-7_40` | Gram ceiling at R = 7/40 (Z2xZ4): floor(1/(1-R)) = 1, measured maximum same-(cell,role) multiplicity = 1 [tight]. | **1** | (cell,role) multiplicity | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json::multR_rec.10|7/40` *(max_int_key)* | VERIFIED |

## Occupancy patterns

| id | claim | value | unit | script | file :: key | status |
|---|---|---:|---|---|---|---|
| `A-23` | Model-cells carrying at least one of Q, u^c, e^c -- the occupancy-pattern denominator. THE 42,932. NOT the denominator of 104,889, which counts conjugate roles too; that is DEN-01 = 43,466. | **42,932** | model-cells | `su5_ten_bound/pattern_census.py` | `su5_ten_bound/patterns_Z2xZ4_1_1.json::cells10` + `su5_ten_bound/patterns_Z3xZ6_1_1.json::cells10` | VERIFIED |
| `A-24` | Cells carrying all three norm-6/5 roles (pattern 111). <br/>*Why 'every fixed point carries exactly one' is false as English.* | **26,134** | cells | `su5_ten_bound/pattern_census.py` | `su5_ten_bound/patterns_Z2xZ4_1_1.json::pat.111` + `su5_ten_bound/patterns_Z3xZ6_1_1.json::pat.111` | VERIFIED |
| `A-pat110` | Occupancy pattern 110 (Q + u^c, no e^c), derived impossible: occurrences over both catalogues. | **0** | cells | `su5_ten_bound/pattern_census.py` | `su5_ten_bound/patterns_Z2xZ4_1_1.json::DERIVED_IMPOSSIBLE.110` + `su5_ten_bound/patterns_Z3xZ6_1_1.json::DERIVED_IMPOSSIBLE.110` | VERIFIED |
| `A-pat101` | Occupancy pattern 101 (Q + e^c, no u^c), derived impossible: occurrences over both catalogues. | **0** | cells | `su5_ten_bound/pattern_census.py` | `su5_ten_bound/patterns_Z2xZ4_1_1.json::DERIVED_IMPOSSIBLE.101` + `su5_ten_bound/patterns_Z3xZ6_1_1.json::DERIVED_IMPOSSIBLE.101` | VERIFIED |

## Limitations

| id | claim | value | unit | script | file :: key | status |
|---|---|---:|---|---|---|---|
| `A-25` | (fixed point, gauge weight) GROUPS of two or more left-chiral records sharing one p_sh -- the H5 gap. THE 21,981. This is a count of GROUPS, not of cells; they are hosted by 17,236 model-cells. <br/>*Configurations outside the theorem: same gauge momentum, distinguished only by oscillator content.* | **21,981** | (cell, p_sh) groups | `gate_su5_bound/d2_pshdup.py` | `gate_su5_bound/d2_pshdup_Z2xZ4_1_1.json::dup_groups` + `gate_su5_bound/d2_pshdup_Z3xZ6_1_1.json::dup_groups` | VERIFIED (SUPERSEDED by UNI-27) |
| `A-26` | Of those, cells carrying any Standard-Model role. <br/>*The gap is closed by MEASUREMENT, not by proof.* | **0** | cells | `gate_su5_bound/d2_pshdup.py` | `gate_su5_bound/d2_pshdup_Z2xZ4_1_1.json::dup_groups_role` + `gate_su5_bound/d2_pshdup_Z3xZ6_1_1.json::dup_groups_role` | VERIFIED (SUPERSEDED by UNI-29) |

## Exclusions

| id | claim | value | unit | script | file :: key | status |
|---|---|---:|---|---|---|---|
| `A-27` | Z2xZ4_1_1 catalogue size (models loaded). | **11,540** | models | `gate_su5_bound/d5_excluded.py` | `gate_su5_bound/d5_excluded.json::Z2xZ4_1_1.loaded` | VERIFIED |
| `A-28` | Z2xZ4_1_1 models excluded (first vev-configuration does not close). | **615** | models | `gate_su5_bound/d5_excluded.py` | `gate_su5_bound/d5_excluded.json::Z2xZ4_1_1.skip.closure_failed` | VERIFIED (SUPERSEDED by UNI-42) |
| `A-29` | Of those, models that DO close under a later printed configuration. | **599** | models | `gate_su5_bound/d5_excluded.py` | `gate_su5_bound/d5_excluded.json::Z2xZ4_1_1.excl_any_cfg_closes` | VERIFIED |
| `A-30` | Z3xZ6_1_1 model labels appearing twice in the CENSUS OUTPUT. The published catalogue file has 4,612 blocks with distinct labels and distinct content; the duplication is a pilot-segment overlap in this program's own run (census_segPILOT re-runs MSSM0-MSSM9, which are also in census_seg0). Over-counting, so absence claims are safe. | **10** | labels | `gate_su5_bound/d5_excluded.py` | `gate_su5_bound/d5_excluded.json::Z3xZ6_1_1.duplicate_labels` | VERIFIED |
| `A-31` | Z3xZ6_1_1 models excluded (uncloseable). | **2** | models | `gate_su5_bound/d5_excluded.py` | `gate_su5_bound/d5_excluded.json::Z3xZ6_1_1.skip.closure_failed` | VERIFIED (SUPERSEDED by UNI-44) |
| `A-32` | Excluded models found to have no hypercharge generator (DATA-GAP). <br/>*Neither Paper A catalogue contains a hypercharge DATA-GAP model.* | **0** | models | `gate_su5_bound/d5_excluded.py` | `gate_su5_bound/d5_excluded.json::Z2xZ4_1_1.excl_cfg0_no_generator` | VERIFIED |
| `A-33` | Exclusion accounting closes (censused + excluded = catalogue). | **True** | boolean | `gate_su5_bound/d5_excluded.py` | `gate_su5_bound/d5_excluded.json::Z2xZ4_1_1.accounting_ok` | VERIFIED |

## The budgets table (seven geometries)

| id | claim | value | unit | script | file :: key | status |
|---|---|---:|---|---|---|---|
| `BUD-Z2xZ4_1_1-sectors` | Twisted sectors of Z2xZ4_1_1. | **7** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z2xZ4_1_1.n_twisted` | VERIFIED |
| `BUD-Z2xZ4_1_1-Bmax` | Maximum oscillator-free budget B_max of Z2xZ4_1_1. | **13/16** | B | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z2xZ4_1_1.B_max` | VERIFIED |
| `BUD-Z2xZ4_1_1-Bmin` | Minimum oscillator-free budget B_min of Z2xZ4_1_1. | **11/16** | B | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z2xZ4_1_1.B_min` | VERIFIED |
| `BUD-Z2xZ4_1_1-loophole` | Sectors of Z2xZ4_1_1 permitting a second norm-6/5 role (B >= 17/20). | **0** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z2xZ4_1_1.ten_doubling_permitted` *(len)* | VERIFIED |
| `BUD-Z2xZ4_1_1-susy` | SUSY condition sum_a eta_a in Z: violations in Z2xZ4_1_1. | **[]** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z2xZ4_1_1.susy_violations` | VERIFIED |
| `BUD-Z2xZ4_1_1-5bar` | Sectors of Z2xZ4_1_1 forbidding a second norm-4/5 role. <br/>*Oscillator-free statement; 4,682 oscillator fields below 13/20 exist and none doubles.* | **[]** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z2xZ4_1_1.fivebar_doubling_forbidden` | VERIFIED |
| `BUD-Z4xZ4_1_1-sectors` | Twisted sectors of Z4xZ4_1_1. | **15** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z4xZ4_1_1.n_twisted` | VERIFIED |
| `BUD-Z4xZ4_1_1-Bmax` | Maximum oscillator-free budget B_max of Z4xZ4_1_1. | **13/16** | B | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z4xZ4_1_1.B_max` | VERIFIED |
| `BUD-Z4xZ4_1_1-Bmin` | Minimum oscillator-free budget B_min of Z4xZ4_1_1. | **11/16** | B | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z4xZ4_1_1.B_min` | VERIFIED |
| `BUD-Z4xZ4_1_1-loophole` | Sectors of Z4xZ4_1_1 permitting a second norm-6/5 role (B >= 17/20). | **0** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z4xZ4_1_1.ten_doubling_permitted` *(len)* | VERIFIED |
| `BUD-Z4xZ4_1_1-susy` | SUSY condition sum_a eta_a in Z: violations in Z4xZ4_1_1. | **[]** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z4xZ4_1_1.susy_violations` | VERIFIED |
| `BUD-Z4xZ4_1_1-5bar` | Sectors of Z4xZ4_1_1 forbidding a second norm-4/5 role. <br/>*Oscillator-free statement; 4,682 oscillator fields below 13/20 exist and none doubles.* | **[]** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z4xZ4_1_1.fivebar_doubling_forbidden` | VERIFIED |
| `BUD-Z3xZ3_1_1-sectors` | Twisted sectors of Z3xZ3_1_1. | **8** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z3xZ3_1_1.n_twisted` | VERIFIED |
| `BUD-Z3xZ3_1_1-Bmax` | Maximum oscillator-free budget B_max of Z3xZ3_1_1. | **7/9** | B | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z3xZ3_1_1.B_max` | VERIFIED |
| `BUD-Z3xZ3_1_1-Bmin` | Minimum oscillator-free budget B_min of Z3xZ3_1_1. | **2/3** | B | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z3xZ3_1_1.B_min` | VERIFIED |
| `BUD-Z3xZ3_1_1-loophole` | Sectors of Z3xZ3_1_1 permitting a second norm-6/5 role (B >= 17/20). | **0** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z3xZ3_1_1.ten_doubling_permitted` *(len)* | VERIFIED |
| `BUD-Z3xZ3_1_1-susy` | SUSY condition sum_a eta_a in Z: violations in Z3xZ3_1_1. | **[]** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z3xZ3_1_1.susy_violations` | VERIFIED |
| `BUD-Z3xZ3_1_1-5bar` | Sectors of Z3xZ3_1_1 forbidding a second norm-4/5 role. <br/>*Oscillator-free statement; 4,682 oscillator fields below 13/20 exist and none doubles.* | **[]** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z3xZ3_1_1.fivebar_doubling_forbidden` | VERIFIED |
| `BUD-Z6-I_1_1-sectors` | Twisted sectors of Z6-I_1_1. | **5** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z6-I_1_1.n_twisted` | VERIFIED |
| `BUD-Z6-I_1_1-Bmax` | Maximum oscillator-free budget B_max of Z6-I_1_1. | **3/4** | B | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z6-I_1_1.B_max` | VERIFIED |
| `BUD-Z6-I_1_1-Bmin` | Minimum oscillator-free budget B_min of Z6-I_1_1. | **2/3** | B | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z6-I_1_1.B_min` | VERIFIED |
| `BUD-Z6-I_1_1-loophole` | Sectors of Z6-I_1_1 permitting a second norm-6/5 role (B >= 17/20). | **0** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z6-I_1_1.ten_doubling_permitted` *(len)* | VERIFIED |
| `BUD-Z6-I_1_1-susy` | SUSY condition sum_a eta_a in Z: violations in Z6-I_1_1. | **[]** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z6-I_1_1.susy_violations` | VERIFIED |
| `BUD-Z6-I_1_1-5bar` | Sectors of Z6-I_1_1 forbidding a second norm-4/5 role. <br/>*Oscillator-free statement; 4,682 oscillator fields below 13/20 exist and none doubles.* | **[]** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z6-I_1_1.fivebar_doubling_forbidden` | VERIFIED |
| `BUD-Z6-II_1_1-sectors` | Twisted sectors of Z6-II_1_1. | **5** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z6-II_1_1.n_twisted` | VERIFIED |
| `BUD-Z6-II_1_1-Bmax` | Maximum oscillator-free budget B_max of Z6-II_1_1. | **7/9** | B | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z6-II_1_1.B_max` | VERIFIED |
| `BUD-Z6-II_1_1-Bmin` | Minimum oscillator-free budget B_min of Z6-II_1_1. | **25/36** | B | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z6-II_1_1.B_min` | VERIFIED |
| `BUD-Z6-II_1_1-loophole` | Sectors of Z6-II_1_1 permitting a second norm-6/5 role (B >= 17/20). | **0** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z6-II_1_1.ten_doubling_permitted` *(len)* | VERIFIED |
| `BUD-Z6-II_1_1-susy` | SUSY condition sum_a eta_a in Z: violations in Z6-II_1_1. | **[]** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z6-II_1_1.susy_violations` | VERIFIED |
| `BUD-Z6-II_1_1-5bar` | Sectors of Z6-II_1_1 forbidding a second norm-4/5 role. <br/>*Oscillator-free statement; 4,682 oscillator fields below 13/20 exist and none doubles.* | **[]** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z6-II_1_1.fivebar_doubling_forbidden` | VERIFIED |
| `BUD-Z4_1_1-sectors` | Twisted sectors of Z4_1_1. | **3** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z4_1_1.n_twisted` | VERIFIED |
| `BUD-Z4_1_1-Bmax` | Maximum oscillator-free budget B_max of Z4_1_1. | **3/4** | B | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z4_1_1.B_max` | VERIFIED |
| `BUD-Z4_1_1-Bmin` | Minimum oscillator-free budget B_min of Z4_1_1. | **11/16** | B | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z4_1_1.B_min` | VERIFIED |
| `BUD-Z4_1_1-loophole` | Sectors of Z4_1_1 permitting a second norm-6/5 role (B >= 17/20). | **0** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z4_1_1.ten_doubling_permitted` *(len)* | VERIFIED |
| `BUD-Z4_1_1-susy` | SUSY condition sum_a eta_a in Z: violations in Z4_1_1. | **[]** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z4_1_1.susy_violations` | VERIFIED |
| `BUD-Z4_1_1-5bar` | Sectors of Z4_1_1 forbidding a second norm-4/5 role. <br/>*Oscillator-free statement; 4,682 oscillator fields below 13/20 exist and none doubles.* | **[]** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z4_1_1.fivebar_doubling_forbidden` | VERIFIED |
| `BUD-Z3xZ6_1_1-sectors` | Twisted sectors of Z3xZ6_1_1. | **17** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z3xZ6_1_1.n_twisted` | VERIFIED |
| `BUD-Z3xZ6_1_1-Bmax` | Maximum oscillator-free budget B_max of Z3xZ6_1_1. | **31/36** | B | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z3xZ6_1_1.B_max` | VERIFIED |
| `BUD-Z3xZ6_1_1-Bmin` | Minimum oscillator-free budget B_min of Z3xZ6_1_1. | **2/3** | B | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z3xZ6_1_1.B_min` | VERIFIED |
| `BUD-Z3xZ6_1_1-loophole` | Sectors of Z3xZ6_1_1 permitting a second norm-6/5 role (B >= 17/20). | **2** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z3xZ6_1_1.ten_doubling_permitted` *(len)* | VERIFIED |
| `BUD-Z3xZ6_1_1-susy` | SUSY condition sum_a eta_a in Z: violations in Z3xZ6_1_1. | **[]** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z3xZ6_1_1.susy_violations` | VERIFIED |
| `BUD-Z3xZ6_1_1-5bar` | Sectors of Z3xZ6_1_1 forbidding a second norm-4/5 role. <br/>*Oscillator-free statement; 4,682 oscillator fields below 13/20 exist and none doubles.* | **[]** | sectors | `gate_su5_bound/d8_classwide.py` | `gate_su5_bound/d8_classwide.json::Z3xZ6_1_1.fivebar_doubling_forbidden` | VERIFIED |

## The construction

| id | claim | value | unit | script | file :: key | status |
|---|---|---:|---|---|---|---|
| `CON-01` | Diophantine solutions over 8 budgets x 10 weights x 260 Delta, AS PUBLISHED BY g2_classes.py. SUPERSEDED: that enumeration imposed rem = 1/2 where the derivation gives rem >= 1/2, so the budget B_v = -5/2 was never enumerated. | **3,000** | solutions | `gate_su5_construction/g2_classes.py` | `gate_su5_construction/g2_classes.json::n_solutions` | VERIFIED (SUPERSEDED by CON-01a) |
| `CON-02` | Distinct Lambda-classes among them, AS PUBLISHED BY g2_classes.py. SUPERSEDED: the class key quotiented by a lattice strictly containing E8 (+) E8, merging classes. | **119** | classes | `gate_su5_construction/g2_classes.py` | `gate_su5_construction/g2_classes.json::n_distinct_classes` | VERIFIED (SUPERSEDED by CON-02a) |
| `CON-03` | Inequivalent V2 classes admitting the loophole: EXACTLY TWO. This is the class boundary that makes the enumeration exhaustive. | **2** | classes | `gate_su5_construction/g2_classes.py` | `gate_su5_construction/g2_classes.json::n_inequivalent` | VERIFIED |
| `CON-04` | Class A: (V1, W) pairs passing the X/Y and MI4 screens. | **750,384** | pairs | `gate_su5_construction/g4_lift.py` | `gate_su5_construction/g4_lift_A.json::n_pairs` | VERIFIED |
| `CON-05` | Class A: inequivalent (V1, W) pairs. | **2,912** | pairs | `gate_su5_construction/g4_lift.py` | `gate_su5_construction/g4_lift_A.json::n_inequivalent` | VERIFIED |
| `CON-06` | Class A: literal embeddings built, all passing MI0-MI6. | **2,912** | models | `su5_construction/s5_build.py` | `su5_construction/s5_build.json::n_built` | VERIFIED |
| `CON-07` | Class A: modular-invariance representative failures. | **0** | models | `su5_construction/s5_build.py` | `su5_construction/s5_build.json::n_MI_representative_not_found` | VERIFIED |
| `CON-08` | Class B: (V1, W) pairs passing the screens. | **696,600** | pairs | `gate_su5_construction/g4_lift.py` | `gate_su5_construction/g4_lift_B.json::n_pairs` | VERIFIED |
| `CON-09` | Class B: inequivalent (V1, W) pairs. | **1,432** | pairs | `gate_su5_construction/g4_lift.py` | `gate_su5_construction/g4_lift_B.json::n_inequivalent` | VERIFIED |
| `CON-10` | Class B: literal embeddings built. | **1,432** | models | `gate_su5_construction/g7_buildB.py` | `gate_su5_construction/g7_buildB.json::n_built` | VERIFIED |
| `CON-11` | Class B: build failures. | **0** | models | `gate_su5_construction/g7_buildB.py` | `gate_su5_construction/g7_buildB.json::n_failed` | VERIFIED |
| `CON-12` | Class B: fully Standard-Model-valid lifts. <br/>*g8b_smvalid.json has no generating script in code/gate_su5_construction/. The 4,032 headline depends on this split. REPRODUCTION GAP.* | **1,120** | models | `NO PRODUCING SCRIPT ON DISK -- see reproduction gap R-1` | `gate_su5_construction/g8b_smvalid.json::pattern.(True, True, True, True, True)` | VERIFIED |
| `CON-13` | Class B: lifts excluded as not Standard-Model models. | **312** | models | `NO PRODUCING SCRIPT ON DISK -- see reproduction gap R-1` | `gate_su5_construction/g8b_smvalid.json::pattern.(True, False, False, False, True)` | VERIFIED |
| `CON-14` | SM-valid models over the COMPLETE two-class loophole space: 2,912 + 1,120. | **4,032** | models | `derived from CON-05 and CON-12` | `gate_su5_construction/g4_lift_A.json::n_inequivalent` + `gate_su5_construction/g8b_smvalid.json::pattern.(True, True, True, True, True)` | VERIFIED |
| `CON-15` | Class A models reaching three net generations in ALL THREE norm-6/5 roles simultaneously. | **0** | models | `gate_su5_construction/g5_battery.py` | `gate_su5_construction/g5_battery.json::summary.gen3_all_three_aware` | VERIFIED |
| `CON-16` | Class B models reaching three net generations in all three roles. | **0** | models | `gate_su5_construction/g8_batteryB.py` | `gate_su5_construction/g8_batteryB.json::summary.NET_3_IN_ALL_THREE` | VERIFIED |
| `CON-17` | Class A models carrying a doubled norm-6/5 role at T(0,1)/T(0,5). | **1,812** | models | `su5_construction/s8_battery.py` | `su5_construction/s8_battery.json::n_with_T0l_ten` | VERIFIED |
| `CON-18` | Class A models with Tr Y = 0 (non-anomalous hypercharge). | **481** | models | `su5_construction/s11_battery2.py` | `su5_construction/s11_battery2.json::trY_hist.0` | VERIFIED |
| `CON-19` | Class A models with net Q = +-3, TRIALITY-BLIND (superseded). <br/>*SUPERSEDED by CON-20. s8_battery.py::su3_dim discards the triality sign, so exotics were counted as quark doublets. DO NOT PRINT.* | **326** | models | `su5_construction/s8_battery.py` | `gate_su5_construction/g5_battery.json::summary.netQ3_blind_PHASE` | VERIFIED (SUPERSEDED by CON-20) |
| `CON-20` | Class A models with net Q = +-3, TRIALITY-AWARE (correct). <br/>*275 models change verdict between the two readings.* | **339** | models | `gate_su5_construction/g5_battery.py` | `gate_su5_construction/g5_battery.json::summary.netQ3_aware_GATE` | VERIFIED |
| `CON-21` | Models whose net-Q verdict changes under the triality correction. | **275** | models | `gate_su5_construction/g5_battery.py` | `gate_su5_construction/g5_battery.json::summary.netQ3_disagreements` | VERIFIED |
| `CON-22` | Level-1 unitarity violations across all class-A left-chiral records. | **0** | records | `gate_su5_construction/g5_battery.py` | `gate_su5_construction/g5_battery.json::summary.level1_skipped_total` | VERIFIED |
| `CON-23` | Gauge-Green-Schwarz orthogonality checked model by model. | **2,912** | models | `gate_su5_construction/g5_battery.py` | `gate_su5_construction/g5_battery.json::summary.gs_orth_checked` | VERIFIED |
| `CON-repro-CB1534` | CB1534: records identical between the stored run and a FRESH orbifolder run driven from the embedding data alone. | **147** | records | `gate_su5_construction/g1b_fresh_run.py + run_gate_orb.sh + g1c_compare.py` | `gate_su5_construction/g1c_compare.json::compare.CB1534.fresh_records` | VERIFIED |
| `CON-ident-CB1534` | CB1534: fresh vs stored comparison verdict. | **True** | boolean | `gate_su5_construction/g1c_compare.py` | `gate_su5_construction/g1c_compare.json::compare.CB1534.identical` | VERIFIED |
| `CON-repro-CB1041` | CB1041: records identical between the stored run and a FRESH orbifolder run driven from the embedding data alone. | **214** | records | `gate_su5_construction/g1b_fresh_run.py + run_gate_orb.sh + g1c_compare.py` | `gate_su5_construction/g1c_compare.json::compare.CB1041.fresh_records` | VERIFIED |
| `CON-ident-CB1041` | CB1041: fresh vs stored comparison verdict. | **True** | boolean | `gate_su5_construction/g1c_compare.py` | `gate_su5_construction/g1c_compare.json::compare.CB1041.identical` | VERIFIED |
| `CON-repro-CB1370` | CB1370: records identical between the stored run and a FRESH orbifolder run driven from the embedding data alone. | **154** | records | `gate_su5_construction/g1b_fresh_run.py + run_gate_orb.sh + g1c_compare.py` | `gate_su5_construction/g1c_compare.json::compare.CB1370.fresh_records` | VERIFIED |
| `CON-ident-CB1370` | CB1370: fresh vs stored comparison verdict. | **True** | boolean | `gate_su5_construction/g1c_compare.py` | `gate_su5_construction/g1c_compare.json::compare.CB1370.identical` | VERIFIED |
| `CON-repro-CB2307` | CB2307: records identical between the stored run and a FRESH orbifolder run driven from the embedding data alone. | **219** | records | `gate_su5_construction/g1b_fresh_run.py + run_gate_orb.sh + g1c_compare.py` | `gate_su5_construction/g1c_compare.json::compare.CB2307.fresh_records` | VERIFIED |
| `CON-ident-CB2307` | CB2307: fresh vs stored comparison verdict. | **True** | boolean | `gate_su5_construction/g1c_compare.py` | `gate_su5_construction/g1c_compare.json::compare.CB2307.identical` | VERIFIED |
| `CON-24` | Standard-Model-valid loophole models carrying a doubled norm-6/5 role, over the complete 4,032. Supersedes the 1,482 that was carried in prose. <br/>*SUPERSEDED in v2.2: one representative per pair of classes; the paper now prints the nine-representative total REP-05 (21,168), whose (0,0) slice REP-05a reproduces this row.* | **2,174** | models | `gate_su5_construction/g8c_doubled_smvalid.py` | `gate_su5_construction/g8c_doubled_smvalid.json::doubled_over_sm_valid_total` | VERIFIED (SUPERSEDED by REP-05) |
| `CON-25` | Class B: doubled AND Standard-Model-valid. The unrestricted class-B figure, 400, counts 38 lifts that are not Standard-Model models. <br/>*SUPERSEDED in v2.2: the class-B part of CON-24 (one representative per pair) is no longer printed; see REP-05.* | **362** | models | `gate_su5_construction/g8c_doubled_smvalid.py` | `gate_su5_construction/g8c_doubled_smvalid.json::classB.doubled_over_sm_valid` | VERIFIED (SUPERSEDED by REP-05) |
| `CON-26` | The same count over all 4,344 BUILT models, unrestricted. A different unit, stated so the two are not confused. | **2,212** | models | `gate_su5_construction/g8c_doubled_smvalid.py` | `gate_su5_construction/g8c_doubled_smvalid.json::doubled_over_all_built_total` | VERIFIED |

## Control (a) - below threshold

| id | claim | value | unit | script | file :: key | status |
|---|---|---:|---|---|---|---|
| `CTRL-a-Z2xZ4-13_40` | Control (a): norm-4/5 role at R = 13/40 on Z2xZ4, BELOW the 13/20 threshold -- 3673 (cell, role) pairs, 3868 fields, multiplicity 1 in every one. | `{"1": 3673}` | (cell,role) pairs | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json::multR_rec.5bar|13/40` | VERIFIED |
| `CTRL-a-Z2xZ4-3_40` | Control (a): norm-4/5 role at R = 3/40 on Z2xZ4, BELOW the 13/20 threshold -- 170 (cell, role) pairs, 170 fields, multiplicity 1 in every one. | `{"1": 170}` | (cell,role) pairs | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json::multR_rec.5bar|3/40` | VERIFIED |
| `CTRL-a-Z3xZ6-11_30` | Control (a): norm-4/5 role at R = 11/30 on Z3xZ6, BELOW the 13/20 threshold -- 18 (cell, role) pairs, 18 fields, multiplicity 1 in every one. | `{"1": 18}` | (cell,role) pairs | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z3xZ6_1_1.json::multR_rec.5bar|11/30` | VERIFIED |
| `CTRL-a-Z3xZ6-23_90` | Control (a): norm-4/5 role at R = 23/90 on Z3xZ6, BELOW the 13/20 threshold -- 40 (cell, role) pairs, 40 fields, multiplicity 1 in every one. | `{"1": 40}` | (cell,role) pairs | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z3xZ6_1_1.json::multR_rec.5bar|23/90` | VERIFIED |
| `CTRL-a-Z3xZ6-4_45` | Control (a): norm-4/5 role at R = 4/45 on Z3xZ6, BELOW the 13/20 threshold -- 586 (cell, role) pairs, 586 fields, multiplicity 1 in every one. | `{"1": 586}` | (cell,role) pairs | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z3xZ6_1_1.json::multR_rec.5bar|4/45` | VERIFIED |
| `CTRL-a-total` | Control (a) total: (cell, role) pairs carrying a norm-4/5 role below its threshold, ALL at multiplicity 1. THE 4,487 (over 4,682 fields). <br/>*The same roles reach multiplicity 4 at higher budget in the same catalogues, so the asymmetry tracks 2B - \|w\|^2, not the role label. Field total: 4682.* | **4,487** | (cell,role) pairs | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json::multR_rec.5bar|13/40` + `gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json::multR_rec.5bar|3/40` + `gate_su5_bound/d1_mechanism_Z3xZ6_1_1.json::multR_rec.5bar|11/30` + `gate_su5_bound/d1_mechanism_Z3xZ6_1_1.json::multR_rec.5bar|23/90` + `gate_su5_bound/d1_mechanism_Z3xZ6_1_1.json::multR_rec.5bar|4/45` *(sum_hist_values)* | VERIFIED (SUPERSEDED by UNI-15) |

## Control (b) - above threshold, doubled (CB1534)

| id | claim | value | unit | script | file :: key | status |
|---|---|---:|---|---|---|---|
| `CTRL-b-CB1534` | Control (b): CB1534 carries TWO states of one norm-6/5 role (a genuine Q) at one constructing element, B = 31/36 >= 17/20, R = 47/90, 2*sqrt(47/90) = 1.4453 >= sqrt(2) = 1.4142 -- the bound saturated from the ALLOWED side; ceiling floor(1/(1-47/90)) = 2, and exactly 2 is realised. | **2** | hidden multiplicity | `gate_su5_construction/g1c_compare.py` | `gate_su5_construction/g1c_compare.json::dossiers.CB1534.hits.0.hid` | VERIFIED |
| `CTRL-b-CB1534-R` | Control (b): CB1534's doubled cell sits at R = 47/90. | **47/90** | R | `gate_su5_construction/g1c_compare.py` | `gate_su5_construction/g1c_compare.json::dossiers.CB1534.hits.0.R` | VERIFIED |
| `CTRL-b-CB1534-w2` | Control (b): CB1534's doubled field has \|w\|^2 = 6/5 (a norm-6/5 role). | **6/5** | |w|^2 | `gate_su5_construction/g1c_compare.py` | `gate_su5_construction/g1c_compare.json::dossiers.CB1534.hits.0.w2` | VERIFIED |
| `CTRL-b-CB1534-osc` | Control (b): CB1534's doubled field is oscillator-free. | **no** | oscillator | `gate_su5_construction/g1c_compare.py` | `gate_su5_construction/g1c_compare.json::dossiers.CB1534.hits.0.osc` | VERIFIED |
| `CTRL-b-CB1534-role` | Control (b): the doubled role is a genuine Q (triality-aware). | **Q** | role | `gate_su5_construction/g1c_compare.py` | `gate_su5_construction/g1c_compare.json::dossiers.CB1534.hits.0.role` | VERIFIED |
| `CTRL-b-CB1534-netQ` | Control (b): CB1534 net Q = +3, triality-aware. | **3** | generations | `gate_su5_construction/g1c_compare.py` | `gate_su5_construction/g1c_compare.json::dossiers.CB1534.net_aware.Q` | VERIFIED |
| `CTRL-b-CB1534-trY` | Control (b): CB1534 has Tr Y = 0. | **0** | trace | `gate_su5_construction/g1c_compare.py` | `gate_su5_construction/g1c_compare.json::dossiers.CB1534.trY` | VERIFIED |

## Control (c) - at the ceiling

| id | claim | value | unit | script | file :: key | status |
|---|---|---:|---|---|---|---|
| `CTRL-c-Z2xZ4-23_40` | Control (c): norm-4/5 role at R = 23/40 on Z2xZ4, ABOVE threshold -- ceiling floor(1/(1-R)) = 2, measured maximum multiplicity = 2. SATURATED. | `{"1": 41129, "2": 88}` | multiplicity histogram | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json::multR_rec.5bar|23/40` | VERIFIED |
| `CTRL-c-Z3xZ6-53_90` | Control (c): norm-4/5 role at R = 53/90 on Z3xZ6, ABOVE threshold -- ceiling floor(1/(1-R)) = 2, measured maximum multiplicity = 2. SATURATED. | `{"1": 833, "2": 206}` | multiplicity histogram | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z3xZ6_1_1.json::multR_rec.5bar|53/90` | VERIFIED |
| `CTRL-c-Z2xZ4-7_10` | Control (c): norm-4/5 role at R = 7/10 on Z2xZ4, ABOVE threshold -- ceiling floor(1/(1-R)) = 3, measured maximum multiplicity = 3. SATURATED. | `{"1": 130464, "2": 20188, "3": 883}` | multiplicity histogram | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z2xZ4_1_1.json::multR_rec.5bar|7/10` | VERIFIED |
| `CTRL-c-Z3xZ6-34_45` | Control (c): norm-4/5 role at R = 34/45 on Z3xZ6, ABOVE threshold -- ceiling floor(1/(1-R)) = 4, measured maximum multiplicity = 4. SATURATED. | `{"1": 91885, "2": 6509, "3": 204, "4": 4}` | multiplicity histogram | `gate_su5_bound/d1_mechanism.py` | `gate_su5_bound/d1_mechanism_Z3xZ6_1_1.json::multR_rec.5bar|34/45` | VERIFIED |

## The two-class boundary, re-enumerated independently

| id | claim | value | unit | script | file :: key | status |
|---|---|---:|---|---|---|---|
| `CON-01a` | Solutions (p, Delta) admitting a doubled norm-6/5 role at T(0,1), both the 10 and 10bar branches, from the independent re-enumeration. Supersedes CON-01. | **26,400** | (p, Delta) solutions | `two_classes/paperA_two_classes_enum.py` | `two_classes/paperA_two_classes_enum.json::n_solutions` | VERIFIED |
| `CON-01b` | The same, 10 branch alone (the 10bar branch doubles every count exactly). | **13,200** | (p, Delta) solutions | `two_classes/paperA_two_classes_enum.py` | `two_classes/paperA_two_classes_enum.json::n_solutions_10_branch_only` | VERIFIED |
| `CON-01c` | Distinct gauge momenta p among them, both branches. The 10-branch half of this, 12,600, is the figure the repair memo quotes. | **25,200** | distinct p | `two_classes/paperA_two_classes_enum.py` | `two_classes/paperA_two_classes_enum.json::n_distinct_p` | VERIFIED |
| `CON-02a` | Distinct Lambda-classes under the exact key, both branches. The 10-branch half, 630, is the figure the repair memo quotes. Supersedes CON-02. | **1,260** | lattice classes | `two_classes/paperA_two_classes_enum.py` | `two_classes/paperA_two_classes_enum.json::n_classes_exact` | VERIFIED |
| `CON-02b` | The same under the older, coarser key, which merges classes and is why the published count was low. | **974** | lattice classes | `two_classes/paperA_two_classes_enum.py` | `two_classes/paperA_two_classes_enum.json::n_classes_coarse` | VERIFIED |
| `CON-03a` | Inequivalent V2 classes admitting the loophole, re-derived independently: EXACTLY TWO. The endpoint is unchanged by the correction to CON-01/CON-02. | **2** | classes | `two_classes/paperA_two_classes_enum.py` | `two_classes/paperA_two_classes_enum.json::n_inequivalent` | VERIFIED |
| `CON-03b` | Orbit certificate: every enumerated class reached from the two representatives by an explicit word, so no third class hides in the set. | **1,260** | classes reached | `two_classes/paperA_two_classes_enum.py` | `two_classes/paperA_two_classes_enum.json::orbit_certificate_classes_reached` | VERIFIED |
| `CON-03d` | Breadth-first certificate, written out row by row: each of the 1,260 enumerated classes with its representative momentum, its parent and the connecting generator, seeded at the two printed representatives V2^A and V2^B (1,200 and 60 rows). two_classes/check_bfs_certificate.py verifies every parent link, that every chain ends at a printed representative, and that the 1,260 classes are distinct. Added 2026-10-22: before this, the certificate was computed and asserted (CON-03b) but no row was stored. | **1,260** | classes | `two_classes/paperA_two_classes_enum.py (writer) + two_classes/check_bfs_certificate.py (checker)` | `two_classes/paperA_two_classes_bfs_certificate.json::n_classes` | VERIFIED |
| `CON-03c` | Closure: generator applications checked to land inside the enumerated set, zero escapes. Without it a union-find computes subgraph components, not orbits. | **332,640** | generator applications | `two_classes/paperA_two_classes_enum.py` | `two_classes/paperA_two_classes_enum.json::closure_checks` | VERIFIED |
| `CON-01r` | Reconciliation: the re-enumeration restricted to the published enumeration's own subset gives 6,000 solutions where g2_classes.py gave 3,000. | **6,000** | solutions | `two_classes/paperA_two_classes_enum.py` | `two_classes/paperA_two_classes_enum.json::reconciliation.gate_style_solutions` | VERIFIED |
| `CON-02r` | Reconciliation: on that same subset the coarse key gives 238 classes where g2_classes.py gave 119; the exact key gives 300. | **238** | lattice classes | `two_classes/paperA_two_classes_enum.py` | `two_classes/paperA_two_classes_enum.json::reconciliation.coarse_classes_on_gate_subset` | VERIFIED |

## Denominators and units

| id | claim | value | unit | script | file :: key | status |
|---|---|---:|---|---|---|---|
| `DEN-01` | Model-cells carrying a norm-6/5 role of EITHER chirality. THE MATCHED-UNIT DENOMINATOR of 104,889. | **43,466** | model-cells | `gate_su5_bound/d5_independence.json (gate record)` | `gate_su5_bound/d5_independence.json::task2_what_is_67555.distinct_cells_measured_by_this_gate.Z2xZ4_1_1_10_family` + `gate_su5_bound/d5_independence.json::task2_what_is_67555.distinct_cells_measured_by_this_gate.Z3xZ6_1_1_10_family` | VERIFIED (SUPERSEDED by UNI-16) |
| `DEN-02` | Model-cells carrying a norm-4/5 role of either chirality. | **191,509** | model-cells | `gate_su5_bound/d5_independence.json (gate record)` | `gate_su5_bound/d5_independence.json::task2_what_is_67555.distinct_cells_measured_by_this_gate.Z2xZ4_1_1_5bar_family` + `gate_su5_bound/d5_independence.json::task2_what_is_67555.distinct_cells_measured_by_this_gate.Z3xZ6_1_1_5bar_family` | VERIFIED (SUPERSEDED by UNI-17) |
| `DEN-03` | Norm-6/5 (cell, role) pairs per model-cell carrying any norm-6/5 role, on Z2xZ4_1_1 alone. This is the 2.48 the record quotes: it is ONE geometry, not both. | **2.482** | pairs per model-cell | `derived from A-05 and DEN-01` | `su5_ten_bound/verify_Z2xZ4_1_1.json::w2_by_fam.10.6/5` + `gate_su5_bound/d5_independence.json::task2_what_is_67555.distinct_cells_measured_by_this_gate.Z2xZ4_1_1_10_family` *(ratio3)* | VERIFIED |
| `DEN-04` | (cell, norm-4/5 role) PAIRS, both catalogues. The pair unit beside the field unit of A-06; the record never printed this as a total. | **346,151** | (cell,role) pairs | `su5_ten_bound/verify.py` | `su5_ten_bound/verify_Z2xZ4_1_1.json::mult_by_fam.5bar` + `su5_ten_bound/verify_Z3xZ6_1_1.json::mult_by_fam.5bar` *(sum_hist_values)* | VERIFIED (SUPERSEDED by UNI-14) |
| `DEN-05` | (cell, norm-6/5 role) PAIRS, both catalogues, counted in the pair unit. Equals the field count 104,889 because that multiplicity is 1 everywhere: the result. | **104,889** | (cell,role) pairs | `su5_ten_bound/verify.py` | `su5_ten_bound/verify_Z2xZ4_1_1.json::mult_by_fam.10` + `su5_ten_bound/verify_Z3xZ6_1_1.json::mult_by_fam.10` *(sum_hist_values)* | VERIFIED (SUPERSEDED by UNI-13) |
| `DEN-06` | Same-role record pairs under ONE role convention applied to both catalogues. THE NUMBER TO PRINT, in place of 31,822. | **32,294** | record pairs | `root_lift/rl_b_data.py` | `root_lift/questionB_Z2xZ4_1_1.json::pairs` + `root_lift/questionB_Z3xZ6_1_1.json::pairs` | VERIFIED (SUPERSEDED by UNI-22) |
| `DEN-07` | Distinct Z3xZ6_1_1 models censused. Ten were censused twice because a pilot segment re-runs MSSM0-MSSM9, which are also in segment 0. <br/>*With Z2xZ4_1_1's 10,925 (duplicate_labels = 0) the total is 15,535 distinct models over 15,545 census runs.* | **4,610** | distinct models | `gate_su5_bound/d5_excluded.py` | `gate_su5_bound/d5_excluded.json::Z3xZ6_1_1.ok` + `gate_su5_bound/d5_excluded.json::Z3xZ6_1_1.duplicate_labels` *(difference)* | VERIFIED (SUPERSEDED by UNI-39) |

## The hidden-SU(2) signature, at its stated scope

| id | claim | value | unit | script | file :: key | status |
|---|---|---:|---|---|---|---|
| `SIG-02` | Same-role pairs printed as two records whose Delta IS a root of the unbroken four-dimensional group: ZERO of 32,294, both catalogues. The general failure the referee names is realised, not hypothetical. | **0** | record pairs | `root_lift/rl_b_data.py` | `root_lift/questionB_Z2xZ4_1_1.json::by.pairs` + `root_lift/questionB_Z3xZ6_1_1.json::by.pairs` *(sum_matching, key contains `4droot=True`)* | VERIFIED (SUPERSEDED by UNI-41) |
| `SIG-03` | Pairs at a lift-index-1 element sharing right-mover AND oscillator, the configuration the derivation says must merge into one record: ZERO, as required. | **0** | record pairs | `root_lift/rl_b_data.py` | `root_lift/questionB_Z2xZ4_1_1.json::liftsector_distinct_record_sameq` + `root_lift/questionB_Z3xZ6_1_1.json::liftsector_distinct_record_sameq` | VERIFIED |
| `SIG-04` | Route G versus route I disagreements in the four-dimensional test: zero. | **0** | disagreements | `root_lift/rl_b_data.py` | `root_lift/questionB_Z2xZ4_1_1.json::route_disagree` + `root_lift/questionB_Z3xZ6_1_1.json::route_disagree` | VERIFIED |
| `SIG-01` | Twisted constructing elements over the seven campaign geometries at which the local-to-four-dimensional root lift is FORCED by the geometry: 6 of 533. The referee's general objection is correct almost everywhere, and this is the honest half of the statement. | **6** | constructing elements | `root_lift/rl_e_summary.py` | `root_lift/questionA_summary.json::n_lift` | VERIFIED |
| `SIG-01b` | Twisted constructing elements surveyed, seven geometries. | **533** | constructing elements | `root_lift/rl_e_summary.py` | `root_lift/questionA_summary.json::n_ce` | VERIFIED |
| `SIG-05` | Z2xZ4_1_1 -- this paper's OTHER catalogue -- constructing elements at which the lift is forced: ZERO of 74. No statement of the signature's form is licensed there. | **0** | constructing elements | `root_lift/rl_e_summary.py` | `root_lift/questionA_summary.json::per_geometry.Z2xZ4_1_1.lift_forced_at` | VERIFIED |
| `SIG-06` | Z3xZ6_1_1 constructing elements at which the lift is forced: 5 of 100. The two sectors where the lifting element is the sector's ONLY element -- T(0,1) and T(0,5) -- are the whole scope of the signature claim. | **5** | constructing elements | `root_lift/rl_e_summary.py` | `root_lift/questionA_summary.json::per_geometry.Z3xZ6_1_1.lift_forced_at` | VERIFIED |
| `SIG-07` | Records examined for the right-moving-weight premise of the signature section: the 2,912 models built over the FIRST order-six shift class (3,801 + 3,680 at T(0,1) and T(0,5); the second class's 1,432 models are not in this half) together with the Z3xZ6_1_1 catalogue at the same two sectors, counted once per model as Table 1 counts them (23,441). One distinct right-moving q_sh per sector, no exception. The as-run catalogue figure, in which ten models are censused twice, is 31,053; the paper prints the de-duplicated total. | **30,922** | records | `uniform_census/sig_premise_counts.py` | `uniform_census/sig_premise_counts.json::printed_by_the_paper.right_mover_records` | VERIFIED |
| `SIG-08` | Catalogue half of the left-moving-oscillator premise: records at \|p_sh\|^2 = 31/18 in T(0,1) and T(0,5) of Z3xZ6_1_1, each model counted once, all with the oscillator flag absent. As-run, with the ten duplicated models counted twice, the figure is 12,406. | **12,350** | records | `uniform_census/sig_premise_counts.py` | `uniform_census/sig_premise_counts.json::printed_by_the_paper.osc_absent_catalogue` | VERIFIED |

## verification (uniform convention)

| id | claim | value | unit | script | file :: key | status |
|---|---|---:|---|---|---|---|
| `UNI-01` | Distinct models censused, both catalogues, under one role convention and one embedding rule, each model counted once. THE 16,138. | **16,138** | models | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.models` | VERIFIED |
| `UNI-02` | Twisted MODEL-CELLS -- the 128 twisted constructing elements counted once per model in which they are populated. Not a count of fixed points. | **880,589** | model-cells | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.twisted_model_cells` | VERIFIED |
| `UNI-03` | Distinct twisted constructing elements across the two space groups. A property of the geometries, not an aggregate. | **128** | constructing elements | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.constructing_elements_twisted` | VERIFIED |
| `UNI-04` | Standard-Model role-carrying fields, record unit. | **502,396** | fields | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.role_carrying_fields_rec` | VERIFIED |
| `UNI-05` | The same in the state unit: each record expanded into the d states of the hidden representation it carries. | **518,373** | states | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.role_carrying_fields_st` | VERIFIED |
| `UNI-06` | Fields of norm 6/5, record unit. | **109,717** | fields | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.fields_10_rec` | VERIFIED |
| `UNI-07` | Fields of norm 6/5, state unit. Equal to the record count because every norm-6/5 record has hidden dimension one. | **109,717** | states | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.fields_10_st` | VERIFIED |
| `UNI-08` | Fields of norm 4/5, record unit. | **392,679** | fields | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.fields_5bar_rec` | VERIFIED |
| `UNI-09` | Fields of norm 4/5, state unit. | **408,656** | states | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.fields_5bar_st` | VERIFIED |
| `UNI-10` | Models at k_Y = 5/3: every censused model. | **16,138** | models | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.kY_5_over_3` | VERIFIED |
| `UNI-11` | Fields on which the massless budget identity \|p\|^2/2 = 1 - dc - N_osc was checked. | **502,396** | fields | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.budget_identity_checked` | VERIFIED |
| `UNI-12` | Failures of that identity. | **0** | fields | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.budget_identity_bad` | VERIFIED |
| `UNI-13` | (cell, norm-6/5 role) PAIRS. Equal to the field count because that multiplicity is one everywhere. | **109,717** | (cell,role) pairs | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.slots_10` | VERIFIED |
| `UNI-14` | (cell, norm-4/5 role) PAIRS. | **360,070** | (cell,role) pairs | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.slots_5bar` | VERIFIED |
| `UNI-15` | Of those, the ones below the norm-4/5 threshold, R < 1/2. | **4,566** | (cell,role) pairs | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.below_slots_5bar` | VERIFIED |
| `UNI-16` | Model-cells carrying a norm-6/5 role of EITHER chirality -- the matched-unit denominator of 109,717. | **45,450** | model-cells | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.cells_with_a_10_role` | VERIFIED |
| `UNI-17` | Model-cells carrying a norm-4/5 role of either chirality. | **198,385** | model-cells | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.cells_with_a_5bar_role` | VERIFIED |
| `UNI-18` | Violations of sqrt(R) + sqrt(R') >= sqrt(2), record unit. | **0** | record pairs | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.violations_record` | VERIFIED |
| `UNI-19` | Violations in the STATE unit, the stronger test: two states inside one record share R, so a doubled record below threshold would violate the bound outright. | **0** | state pairs | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.violations_state` | VERIFIED |
| `UNI-20` | Norm-6/5 fields at a cell whose budget lets the threshold permit a repeat. | **65** | fields | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.unprotected_10_fields` | VERIFIED |
| `UNI-21` | The same count in the state unit. | **65** | states | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.unprotected_10_fields_state` | VERIFIED |
| `UNI-22` | Same-role pairs of distinct records. THE 33,743. | **33,743** | record pairs | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.pairs_rec` | VERIFIED |
| `UNI-23` | Same-role pairs of distinct STATES. | **51,879** | state pairs | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.pairs_st` | VERIFIED |
| `UNI-24` | Of those, the ones lying across distinct records. | **35,491** | state pairs | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.pairs_st_cross` | VERIFIED |
| `UNI-25` | Same-role pairs whose difference is not in Lambda. | **0** | record pairs | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.delta_not_in_lambda` | VERIFIED |
| `UNI-26` | Same-role pairs whose difference is not orthogonal to W, tested against a basis recovered per model by explicit projection. | **0** | record pairs | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.delta_not_perp_W` | VERIFIED |
| `UNI-27` | (cell, gauge weight) GROUPS of two or more records -- the H5 residual. | **23,619** | groups | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.dup_groups` | VERIFIED |
| `UNI-28` | Model-cells hosting at least one such group. | **18,498** | model-cells | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.dup_cells` | VERIFIED |
| `UNI-29` | Of those groups, the ones carrying a Standard-Model role. | **0** | groups | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.dup_groups_carrying_a_role` | VERIFIED |
| `UNI-30` | Norm-6/5 (cell, role) pairs per model-cell carrying a norm-6/5 role. THE RATIO TO PRINT. | **2.414** | pairs per model-cell | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.ratio_10_pairs_per_cell` | VERIFIED |
| `UNI-31` | The same ratio for the norm-4/5 family. | **1.815** | pairs per model-cell | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.ratio_5bar_pairs_per_cell` | VERIFIED |
| `UNI-32` | Single-budget groups exceeding the ceiling floor(1/(1-R)), state unit. | **0** | groups | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.ceiling_overruns` | VERIFIED |
| `UNI-33` | Largest per-(cell, role) multiplicity of a norm-6/5 role, Z2xZ4_1_1: one. | **1** | records per (cell,role) | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_e6_Z2xZ4_1_1.json::uniform.mult_rec_by_fam.10` *(max_int_key)* | VERIFIED |
| `UNI-34` | The same on Z3xZ6_1_1. | **1** | records per (cell,role) | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_e6_Z3xZ6_1_1.json::uniform.mult_rec_by_fam.10` *(max_int_key)* | VERIFIED |
| `UNI-35` | Largest hidden dimension of a norm-6/5 record, Z2xZ4_1_1: one. | **1** | hidden dimension | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_e6_Z2xZ4_1_1.json::uniform.hid_by_fam.10` *(max_int_key)* | VERIFIED |
| `UNI-36` | The same on Z3xZ6_1_1. | **1** | hidden dimension | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_e6_Z3xZ6_1_1.json::uniform.hid_by_fam.10` *(max_int_key)* | VERIFIED |
| `UNI-37` | Models entering at a later printed embedding than the first. | **599** | models | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.models_recovered_at_a_later_embedding` | VERIFIED |
| `UNI-38` | Z2xZ4_1_1 models censused. | **11,526** | models | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::per_geometry.models.Z2xZ4_1_1` | VERIFIED |
| `UNI-39` | Z3xZ6_1_1 models censused -- the whole catalogue. | **4,612** | models | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::per_geometry.models.Z3xZ6_1_1` | VERIFIED |
| `UNI-40` | Models on which the Standard-Model Cartan span could not be recovered. | **0** | models | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::totals.smbasis_fail` | VERIFIED |
| `UNI-43` | Z2xZ4_1_1 models closing at the FIRST printed embedding. | **10,927** | models | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_e6_Z2xZ4_1_1.json::asrun.models` | VERIFIED |
| `UNI-42` | Z2xZ4_1_1 models that print a single embedding and do not close: the residual the census does not reach. | **14** | models | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_e6_Z2xZ4_1_1.json::blocks` + `uniform_census/uniform_e6_Z2xZ4_1_1.json::uniform.models` *(difference)* | VERIFIED |
| `UNI-44` | Z3xZ6_1_1 models that do not close: none. | **0** | models | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_e6_Z3xZ6_1_1.json::distinct_labels` + `uniform_census/uniform_e6_Z3xZ6_1_1.json::uniform.models` *(difference)* | VERIFIED |
| `UNI-41` | Same-role pairs of distinct records whose difference IS a root of the unbroken four-dimensional gauge group: zero, over the uniform pair population. | **0** | record pairs | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/rootlift_uniform_Z2xZ4_1_1.json::pairs_4droot_true` + `uniform_census/rootlift_uniform_Z3xZ6_1_1.json::pairs_4droot_true` | VERIFIED |
| `UNI-45` | Same-role pairs of distinct records over which the four-dimensional-root test was run, uniform convention. | **33,743** | record pairs | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/rootlift_uniform_Z2xZ4_1_1.json::pairs` + `uniform_census/rootlift_uniform_Z3xZ6_1_1.json::pairs` | VERIFIED |
| `UNI-46` | Mixed-budget (cell, role) groups on Z2xZ4_1_1 at R_max = 33/40 -- norm-4/5 groups holding records of two distinct R. Each holds exactly one field at R = 13/40, so the group count and the field count coincide; the ceiling annexe's caption prints it as a field count. Under the previous convention this was 195. | **197** | (cell, role) groups | `uniform_census/uniform_census_e6.py + uniform_totals.py` | `uniform_census/uniform_totals.json::mixed_budget.Z2xZ4_1_1.groups.5bar|33/40` | VERIFIED |

## construction (representatives, v2.2)

| id | claim | value | unit | script | file :: key | status |
|---|---|---:|---|---|---|---|
| `REP-01` | Embeddings over the 4,032 Standard-Model-intact pairs of classes, nine per pair, one per slot (a, b) = 3(lambda.V2, mu.V2) mod 3, built MI-valid and run through orbifolder; spectra parsed. Every deck's (a, b) is re-computed from its own vectors and agrees with its label (REP-01c, REP-01d). As first assembled, the 80 slots (0,1), (0,2) of 40 class-B pairs held copies of the (0,0) deck (b5b_missing.py never searched mu); repaired before deposit, superseded files in data/representatives/brothers/superseded_v2.2/. | **36,288** | spectra | `representatives/b5_build_brother_census.py + b5b_missing.py + b6_run_brothers.sh + b10_brother_analysis.py` | `representatives/b10_brothers.json::decks_with_spectrum` | VERIFIED |
| `REP-01b` | Of those, decks with no orbifolder output. | **0** | spectra | `representatives/b10_brother_analysis.py` | `representatives/b10_brothers.json::n_missing_output` | VERIFIED |
| `REP-01c` | Label check: decks whose (a, b), re-computed from their own vectors as (3 lambda.V2 mod 3, 3 mu.V2 mod 3) against the parent embedding, equals the slot named by the label and the build record (also V2 = parent's, W4 = W3). | **36,288** | decks | `representatives/b10_brother_analysis.py` | `representatives/b10_brothers.json::label_check.passed` | VERIFIED |
| `REP-01d` | Label check: decks failing it (the superseded v2.2 files give 80, the mislabelled (0,1)/(0,2) decks). | **0** | decks | `representatives/b10_brother_analysis.py` | `representatives/b10_brothers.json::label_check.failed` | VERIFIED |
| `REP-02` | Control: the (a, b) = (0, 0) decks are the v2.1 census representatives; their net indices agree with the certified per-model nets of g5_battery.json / g8_batteryB.json. | **4,032** | models | `representatives/b10_brother_analysis.py` | `representatives/b10_brothers.json::control_vs_certified.agree` | VERIFIED |
| `REP-03` | Embeddings with net Q = u^c = e^c = +-3 at once (triality-aware, hidden-dimension weighted), over all nine representatives of every SM-intact pair. The paper's three-generation negative. | **0** | spectra | `representatives/b10_brother_analysis.py` | `representatives/b10_brothers.json::n_three_generations` | VERIFIED |
| `REP-04` | Embeddings reaching two of the three (net = 3e for two of Q, u^c, e^c, one sign e). | **216** | spectra | `representatives/b10_brother_analysis.py` | `representatives/b10_brothers.json::two_of_three_by_ab` *(sum_values)* | VERIFIED |
| `REP-05` | Embeddings carrying a doubled 10-like representation (a left-chiral 10-like record at T(0,1)/T(0,5) with \|p\|^2 = 31/18 and hidden dimension 2), all nine representatives. Supersedes CON-24 (2,174 = the (0,0) slice). | **21,168** | spectra | `representatives/b10_brother_analysis.py` | `representatives/b10_brothers.json::doubled_by_ab` *(sum_values)* | VERIFIED |
| `REP-05a` | The (0,0) slice of REP-05: reproduces CON-24. | **2,174** | spectra | `representatives/b10_brother_analysis.py` | `representatives/b10_brothers.json::doubled_by_ab.00` | VERIFIED |
| `REP-06` | Non-(0,0) embeddings whose spectrum (left-chiral record count or nets) differs from their pair's (0,0) spectrum. Was 31,972 before the repair of the 80 mislabelled (0,b) decks, which were copies of the (0,0) deck. | **32,052** | spectra | `representatives/b10_brother_analysis.py` | `representatives/b10_brothers.json::differs_from_census_rep_by_ab` *(sum_values)* | VERIFIED |
| `REP-07` | Model (7.1): each of the nine lattice shifts of V2 reproduces one of the nine (a, b) spectra (record count of the shift that reproduces (a, b) = (0, 0)). | **147** | multiplets | `representatives/b11_v2_brothers.py + b2_run_orbifolder.sh + b11b_compare.py` | `representatives/b11_compare.json::NU3.n_records` | VERIFIED |
| `REP-08` | Projection code at every bin, CB1534: (twisted constructing element, massless right-mover) bins of Geometry_Z3xZ6_1_1, both chiralities. | **140** | bins | `representatives/b8_allbins.py` | `representatives/b8_allbins_deposited.json::CB1534.bins` | VERIFIED |
| `REP-08a` | CB1534: left-chirality bins carrying printed oscillator-free states. | **44** | bins | `representatives/b8_allbins.py` | `representatives/b8_allbins_deposited.json::CB1534.populated_bins_tested` | VERIFIED |
| `REP-08b` | CB1534: oscillator-free twisted left-chiral multiplets printed (the multiplets the code reproduces). | **88** | multiplets | `representatives/b8_allbins.py` | `representatives/b8_allbins_deposited.json::CB1534.printed_records_oscfree_twisted` | VERIFIED |
| `REP-08c` | CB1534: survivors of the projection over all 140 bins = twice the 393 printed oscillator-free twisted states (both chiralities). | **786** | states | `representatives/b8_allbins.py` | `representatives/b8_allbins_deposited.json::CB1534.survivors_total` | VERIFIED |
| `REP-08d` | CB1534: agreement at every bin (no left mismatch, no right mismatch, no unprinted survivor). | **True** | boolean | `representatives/b8_allbins.py` | `representatives/b8_allbins_deposited.json::CB1534.agree` | VERIFIED |
| `REP-08e` | CB1041 (the second embedding of the Discussion): agreement at every one of its 140 bins. | **True** | boolean | `representatives/b8_allbins.py` | `representatives/b8_allbins_deposited.json::CB1041.agree` | VERIFIED |
| `REP-09` | Mean-budget ceiling n_states <= 1/(1 - Rbar): violations over every group of the Z2xZ4 catalogue. | **0** | groups | `representatives/b7_census_rank.py` | `representatives/b7_Z2xZ4_1_1.json::n_violations` | VERIFIED |
| `REP-09a` | The same, Z3xZ6 catalogue. | **0** | groups | `representatives/b7_census_rank.py` | `representatives/b7_Z3xZ6_1_1.json::n_violations` | VERIFIED |
| `REP-09b` | Mixed-budget groups of the R_max = 33/40 row (Z2xZ4): one state at 33/40 and one at 13/40 (Rbar = 23/40). | `{"[('33/40', 1), ('13/40', 1)]": 197}` | groups | `representatives/b7_census_rank.py` | `representatives/b7_Z2xZ4_1_1.json::rows.Z2xZ4_1_1|5bar|33/40.mixed` | VERIFIED |
| `REP-09c` | Mixed-budget groups of the R_max = 83/90 row (Z3xZ6): one state at 83/90 and one at 53/90 (Rbar = 34/45). | `{"[('83/90', 1), ('53/90', 1)]": 9}` | groups | `representatives/b7_census_rank.py` | `representatives/b7_Z3xZ6_1_1.json::rows.Z3xZ6_1_1|5bar|83/90.mixed` | VERIFIED |
| `REP-10` | Rank restriction at R = 83/90: min(floor(1/(1-R)), largest A_{n-1} in one E8 orthogonal to W) = 9 (norm 12). | **9** | multiplicity | `representatives/b7_census_rank.py` | `representatives/b7_Z3xZ6_1_1.json::rows.Z3xZ6_1_1|5bar|83/90.min_norm_lat` | VERIFIED |
| `REP-10a` | Groups at R = 83/90 where the rank restriction is tighter than the norm: all of them. | **238** | groups | `representatives/b7_census_rank.py` | `representatives/b7_Z3xZ6_1_1.json::rows.Z3xZ6_1_1|5bar|83/90.lat_tighter_cells` | VERIFIED |

---

## Regenerating any row

Run the script in the `script` column from its own directory; it overwrites the
table in `multi_geometry_campaign/data/`. Then

```
python regeneration/compare_to_deposited.py    # workspace vs pristine reference
python regeneration/claims_check.py            # every row re-resolved
python regeneration/claims_md.py               # rebuild this file
```

"""Paper A -- the UNIFORM-CONVENTION, STATE-LEVEL census (external referee, rounds 2-3).

The referee's standing request, charged in rounds 2 and 3 and carried to the PI as item 5:
recompute every catalogue-wide verification count of Table 1 under ONE role convention
applied to both geometries, with each DISTINCT model counted once, and with the STATE
rather than the printed record as the primary unit.  This script produces that census.
`reconcile.py` turns its JSON into the reconciliation table and the drop-in LaTeX.

THE CONVENTIONS, STATED ONCE (CLAUDE.md 5.19 (i)); every count below is under these.
----------------------------------------------------------------------------------------
U1  ROLE -- UNGATED.  A left-chiral record carries a Standard-Model role if and only if its
    (SU(3)_c, SU(2)_L, Y) assignment is one of the ten role signatures.  No condition of any
    kind is placed on the record's representation under the hidden gauge factors.  This is
    the convention `z24sw_lib.role_census` already uses; `z36sw_lib.role_census` instead
    sets role = None on every record that is not a singlet of every hidden factor, and that
    asymmetry is the 472-pair seam disclosed in the Table 1 caption of paperA_v5.  Here the
    ungated rule is applied to BOTH geometries, and the roles are assigned by this file --
    never read back from either library's `role_census`.

U2  ORIENTATION.  The two per-model sign conventions -- sY (overall sign of Y) and sC
    (which SU(3) weight is called the fundamental) -- are part of the Standard-Model
    IDENTIFICATION, not of the counting, and are taken from the certified per-model
    identification (`role_census`'s own four-branch selection by net Q = +3).  Invariant
    V-BR below re-derives the selection UNGATED and STATE-weighted and reports any model on
    which the two disagree; the count of such models is printed, not smoothed.

U3  EMBEDDING.  orbifolder may print several Standard-Model vev-configurations for one
    model.  The embedding used is the FIRST printed one at which the Standard-Model
    identification closes -- i.e. at which the full left-chiral spectrum reproduces
    orbifolder's own printed multiplicity table exactly.  On Z3xZ6_1_1 that is what
    `observable_positions` already does.  On Z2xZ4_1_1 the census took the first printed
    configuration and stopped, so 615 models failed; walking on to the next printed
    configuration recovers 599 of them.  Under U3 those 599 models are IN, counted once.

U4  MODEL.  One published catalogue entry, counted ONCE.  `census_segPILOT.txt.gz` re-runs
    MSSM0..MSSM9 of Z3xZ6_1_1, which census_seg0.txt.gz also covers; a label seen in more
    than one segment is kept from the first segment in sorted order and skipped thereafter.
    The duplicate map is MEASURED here (pass 1), not assumed, on both geometries.

U5  CELL KEY.  A cell is (sector (k,l), FULL six-component fixed-point tuple n_1..n_6)
    (CLAUDE.md 5.22 addendum (ii)).  `fp_len` records the observed tuple length so a
    truncated key cannot pass unnoticed.  Only twisted cells, (k,l) != (0,0), are counted.

U6  CHIRALITY.  Only left-chiral (`_l`) records are counted; the CPT partner is never
    counted separately.

U7  UNITS.  Two, reported side by side for every row.
      RECORD unit: one printed F_n entry = one count.  This is Table 1's present unit.
      STATE  unit: one record is expanded into d = prod_{k not in (a3,a2)} |rep_k| counts,
                   d being the dimension of the hidden-gauge representation it carries.
                   d is not an assumption: it is the same quantity orbifolder's own printed
                   Standard-Model multiplicity table is built from, and the closure of U3
                   is the check that this file reads it the same way orbifolder wrote it.
    All d states of one record share |p|^2, hence share B and R; they differ only in
    directions the printed output does not resolve.  Statements that need the individual
    hidden momenta (the Delta test) are therefore RECORD-level only and are labelled so.

U8  |w|^2 BY PROJECTION.  R = 2B - |w|^2 with |w|^2 = |P_W p|^2 obtained by explicit
    projection onto W = span(SU(3)_c simple roots, SU(2)_L simple root, t_Y), the basis
    being recovered per model from the unbroken four-dimensional gauge group.  Never
    assumed from k_Y.

WHAT IS IMPORTED AND WHAT IS WRITTEN HERE.  Imported as black boxes: the certified output
PARSERS (`split_models`, `parse_spectrum`, `parse_smconfig`/`parse_gauge`,
`parse_gauge_factors`), the certified closure `observable_positions`, `role_census` used
ONLY for its orientation diagnostic, and the exact lattice/projection helpers of `s5lib`.
Written fresh here: the embedding walk, the role assignment, the hidden-dimension
expansion, every key, every group, every histogram and every aggregate.

INVARIANTS CHECKED BEFORE ANY NEW NUMBER IS BELIEVED
----------------------------------------------------
V-FP   every fixed-point key is a full six-component tuple.
V-W2   |w|^2 takes exactly one value per norm family: 6/5 and 4/5.
V-104  in AS-RUN mode (U3 replaced by "first printed configuration", U4 replaced by "as
       run") the ungated norm-6/5 field count reproduces 67,555 (Z2xZ4_1_1) and 37,334
       (Z3xZ6_1_1), total 104,889 -- the certified figures.
V-322  in the same mode the same-role record-pair count reproduces 24,406 and 7,888, total
       32,294 -- the certified homogeneous (ungated) figures.
V-BR   the orientation branch of U2 agrees with the ungated state-weighted re-derivation.

Output: paper15/review/data/referee_round3/uniform_<geom>.json
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"

import glob                                                             # noqa: E402
import gzip                                                             # noqa: E402
import json                                                             # noqa: E402
import multiprocessing as mp                                            # noqa: E402
import re                                                               # noqa: E402
import sys                                                              # noqa: E402
import time                                                             # noqa: E402
from collections import Counter, defaultdict                            # noqa: E402
from fractions import Fraction as F                                     # noqa: E402

import numpy as np                                                      # noqa: E402

# Package copy: every path is derived from this file's own location, so nothing
# outside the package is read or written.
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "..", ".."))
CAMP = os.path.join(REPO, "multi_geometry_campaign")
OUT = os.path.join(REPO, "results_deposited", "uniform_census")

TEN_P = ("Q", "uc", "ec")
TEN_C = ("Qbar", "ucbar", "ecbar")
FIVE_P = ("dc", "L")
FIVE_C = ("dcbar", "Hu")            # Hu = (1,2)_{+1/2} = the conjugate of L
ROLES10 = set(TEN_P) | set(TEN_C)
ROLES5 = set(FIVE_P) | set(FIVE_C)

GEOM = {
    "Z2xZ4_1_1": dict(lib="z2xz4_sweep", mod="z24sw_lib", M=2, N=4,
                      census=os.path.join(CAMP, "data", "z2xz4_sweep", "raw", "census"),
                      v1=(F(0), F(1, 2), F(-1, 2)), v2=(F(1, 4), F(0), F(-1, 4))),
    "Z3xZ6_1_1": dict(lib="z3xz6_sweep", mod="z36sw_lib", M=3, N=6,
                      census=os.path.join(CAMP, "data", "z3xz6_sweep", "raw", "census"),
                      v1=(F(1, 3), F(0), F(-1, 3)), v2=(F(0), F(1, 6), F(-1, 6))),
}

# oscillator token; the MODE number is what the mass formula wants (CLAUDE.md 5.21 (iv))
_OSC = re.compile(r"(-?)a_\(\s*(-?\d+(?:/\d+)?)\s*\)\^(?:\((-\d)\)|(\d))")

_S5 = None
_R2 = None


def _init():
    global _S5, _R2
    sys.path.insert(0, os.path.join(CAMP, "code", "su5_ten_bound"))
    import s5lib                                                        # noqa: E402
    _S5 = s5lib
    _R2 = np.array([[int(2 * x) for x in a] for a in s5lib.e8xe8_roots()], dtype=np.int64)


def mask_of(V):
    """Boolean mask over e8xe8_roots(): alpha . V in Z.  Exact integer arithmetic."""
    D = 1
    for x in V:
        D = D * F(x).denominator // np.gcd(D, F(x).denominator)
    Vi = np.array([int(F(x) * D) for x in V], dtype=np.int64)
    return ((_R2 @ Vi) % (2 * D)) == 0


def mode_sum(rec):
    """N_osc as the sum of the oscillator MODE numbers, not of occurrences."""
    if "ML" in rec:
        return sum(F(x) for x in rec["ML"]) + sum(F(x) for x in rec["MbL"])
    s = (rec.get("osc") or rec.get("osc_raw") or "no").strip()
    if s in ("no", ""):
        return F(0)
    return sum(abs(F(m)) for _sg, m, _n, _p in _OSC.findall(s))


def pair_ok(a, b):
    """(sqrt(a) + sqrt(b))^2 >= 2, exactly over Q -- equation (mainbound)."""
    t = 2 - a - b
    if t <= 0:
        return True
    return 4 * a * b >= t * t


def ceiling_of(R):
    """floor(1/(1-R)); None outside the corollary's domain R < 1."""
    if R >= 1:
        return None
    return int(1 / (1 - R))


# ------------------------------------------------------------------ accumulator
def newacc():
    return dict(
        models=0, models_recovered=0, embed_index=Counter(), skipped=Counter(),
        kY=Counter(),
        recs_l_all=0, recs_l_tw=0, fp_len=Counter(),
        mcells_tw=0, geo_cells_tw=set(),
        I1_checked=0, I1_bad=0, I1_bad_ex=[],
        w2_by_fam=defaultdict(Counter),
        # fields
        f_rec_by_fam=Counter(), f_st_by_fam=Counter(),
        f_rec_by_role=Counter(), f_st_by_role=Counter(),
        hid_by_fam=defaultdict(Counter),
        # (cell, role) slots
        slots_by_fam=Counter(),
        mult_rec_by_fam=defaultdict(Counter), mult_st_by_fam=defaultdict(Counter),
        pairs_rec_by_fam=Counter(), pairs_st_by_fam=Counter(),
        pairs_st_cross_by_fam=Counter(),
        # occupancy of model-cells
        occ=Counter(),
        # budgets
        R_rec_by_fam=defaultdict(Counter), R_st_by_fam=defaultdict(Counter),
        below_f_rec=Counter(), below_f_st=Counter(),
        below_slots=Counter(), below_mult_rec=defaultdict(Counter),
        below_mult_st=defaultdict(Counter), below_hid_by_fam=defaultdict(Counter),
        loophole_hid=Counter(),
        # the theorem
        viol_rec=0, viol_st=0, viol_ex=[],
        loophole_f_rec=0, loophole_f_st=0, loophole_cells=set(),
        # ceiling annexe: (fam, Rmax) -> max multiplicities
        ceil_rec=defaultdict(int), ceil_st=defaultdict(int),
        ceil_rec_single=defaultdict(int), ceil_st_single=defaultdict(int),
        ceil_over=0,
        # mixed-budget (cell, role) groups: records of two or more distinct R in one
        # group.  Keyed "fam|Rmax" for the groups, "fam|Rmax|R" for their fields.
        # These are the numeral of the ceiling annexe's caption.
        mixed_groups=Counter(), mixed_fields=Counter(),
        # Delta test (RECORD level only)
        d_pairs=0, d_notlambda=0, d_d2=Counter(), d_notperp=0, d_bad_ex=[],
        # (cell, gauge weight) groups
        dup_groups=0, dup_cells=0, dup_groups_role=0, dup_groups_10=0,
        dup_size=Counter(), dup_groups_st=0,
        # V-BR
        br_disagree=0, br_disagree_ex=[],
        smbasis_fail=0,
    )


_INTKEYS = ("models", "models_recovered", "recs_l_all", "recs_l_tw", "mcells_tw",
            "I1_checked", "I1_bad", "viol_rec", "viol_st",
            "loophole_f_rec", "loophole_f_st", "ceil_over",
            "d_pairs", "d_notlambda", "d_notperp",
            "dup_groups", "dup_cells", "dup_groups_role", "dup_groups_10",
            "dup_groups_st", "br_disagree", "smbasis_fail")
_CNTKEYS = ("embed_index", "skipped", "kY", "fp_len", "f_rec_by_fam", "f_st_by_fam",
            "f_rec_by_role", "f_st_by_role", "slots_by_fam", "pairs_rec_by_fam",
            "pairs_st_by_fam", "pairs_st_cross_by_fam", "occ", "below_f_rec",
            "below_f_st", "below_slots", "d_d2", "dup_size", "loophole_hid",
            "mixed_groups", "mixed_fields")
_DDKEYS = ("w2_by_fam", "hid_by_fam", "mult_rec_by_fam", "mult_st_by_fam",
           "R_rec_by_fam", "R_st_by_fam", "below_mult_rec", "below_mult_st",
           "below_hid_by_fam")
_MAXKEYS = ("ceil_rec", "ceil_st", "ceil_rec_single", "ceil_st_single")
_SETKEYS = ("geo_cells_tw", "loophole_cells")
_EXKEYS = ("I1_bad_ex", "viol_ex", "d_bad_ex", "br_disagree_ex")


def merge_into(a, b):
    for k in _INTKEYS:
        a[k] += b[k]
    for k in _CNTKEYS:
        a[k].update(b[k])
    for k in _DDKEYS:
        for x, d in b[k].items():
            a[k][x].update(d)
    for k in _MAXKEYS:
        for x, v in b[k].items():
            if v > a[k][x]:
                a[k][x] = v
    for k in _SETKEYS:
        a[k] |= b[k]
    for k in _EXKEYS:
        a[k] = (a[k] + b[k])[:20]
    return a


# ------------------------------------------------------------------ per-model census
def census_model(geom, lab, spec, roles, hid, smbasis, acc, memo):
    """Accumulate one model.  `roles`/`hid`: {record label: role or None} / {label: d}.

    `memo` caches (B, |w|^2, N_osc) per record label across the two or three convention
    variants of the same model; the cached quantities do not depend on the variant."""
    G = GEOM[geom]
    S5 = _S5
    cells = defaultdict(list)
    for _l, r in spec.items():
        if r["chir"] != "l":
            continue
        acc["recs_l_all"] += 1
        if r["p_sh"] is None:
            continue
        fp = r["fp"]
        acc["fp_len"][len(fp) if fp else "none"] += 1
        if tuple(r["kl"]) == (0, 0):
            continue
        acc["recs_l_tw"] += 1
        cells[(tuple(r["kl"]), tuple(fp) if fp else None)].append(r)

    for key, recs in cells.items():
        kl, _fp = key
        acc["mcells_tw"] += 1
        acc["geo_cells_tw"].add(key)
        v = [(kl[0] * G["v1"][i] + kl[1] * G["v2"][i]) for i in range(3)]
        dc = sum((x % 1) * (1 - (x % 1)) for x in v) / 2

        # ---- (cell, gauge weight) groups of >= 2 records  (the H5 measurement)
        byp = defaultdict(list)
        for r in recs:
            byp[tuple(r["p_sh"])].append(r)
        hosted = 0
        for _psh, rs in byp.items():
            nst = sum(hid.get(r["label"], 1) for r in rs)
            if nst >= 2:
                acc["dup_groups_st"] += 1
            if len(rs) < 2:
                continue
            hosted += 1
            acc["dup_groups"] += 1
            acc["dup_size"][len(rs)] += 1
            ros = [roles.get(r["label"]) for r in rs]
            if any(x is not None for x in ros):
                acc["dup_groups_role"] += 1
            if any(x in ROLES10 for x in ros):
                acc["dup_groups_10"] += 1
        if hosted:
            acc["dup_cells"] += 1

        # ---- roles at this cell
        byrole = defaultdict(list)
        for r in recs:
            ro = roles.get(r["label"])
            if ro in ROLES10 or ro in ROLES5:
                byrole[ro].append(r)
        if not byrole:
            continue
        present = set(byrole)
        if present & set(TEN_P):
            acc["occ"]["ten_p"] += 1
        if present & set(TEN_C):
            acc["occ"]["ten_c"] += 1
        if present & ROLES10:
            acc["occ"]["ten_any"] += 1
        if (present & set(TEN_P)) and (present & set(TEN_C)):
            acc["occ"]["ten_both"] += 1
        if present & ROLES5:
            acc["occ"]["five_any"] += 1
        acc["occ"]["any_role"] += 1

        for ro, rs in byrole.items():
            fam = "10" if ro in ROLES10 else "5bar"
            Rs, Ds, m_st = [], [], 0
            for r in rs:
                mk = memo.get(r["label"])
                if mk is None:
                    B = S5.dot(r["p_sh"], r["p_sh"]) / 2
                    mk = (B, S5.proj_norm2(r["p_sh"], smbasis), mode_sum(r))
                    memo[r["label"]] = mk
                B, w2, Mosc = mk
                acc["I1_checked"] += 1
                if B != 1 - dc - Mosc:
                    acc["I1_bad"] += 1
                    if len(acc["I1_bad_ex"]) < 8:
                        acc["I1_bad_ex"].append(dict(model=lab, kl=list(kl), role=ro,
                                                     B=str(B), dc=str(dc)))
                Rv = 2 * B - w2
                d = hid.get(r["label"], 1)
                Rs.append(Rv)
                Ds.append(d)
                m_st += d
                acc["w2_by_fam"][fam][str(w2)] += 1
                acc["f_rec_by_fam"][fam] += 1
                acc["f_st_by_fam"][fam] += d
                acc["f_rec_by_role"][ro] += 1
                acc["f_st_by_role"][ro] += d
                acc["hid_by_fam"][fam][d] += 1
                acc["R_rec_by_fam"][fam][str(Rv)] += 1
                acc["R_st_by_fam"][fam][str(Rv)] += d
                if Rv < F(1, 2):
                    acc["below_f_rec"][fam] += 1
                    acc["below_f_st"][fam] += d
                    acc["below_hid_by_fam"][fam][d] += 1
                if fam == "10" and Rv >= F(1, 2):
                    acc["loophole_f_rec"] += 1
                    acc["loophole_f_st"] += d
                    acc["loophole_hid"][d] += 1
                    acc["loophole_cells"].add((lab, key))
            m_rec = len(rs)
            acc["slots_by_fam"][fam] += 1
            acc["mult_rec_by_fam"][fam][m_rec] += 1
            acc["mult_st_by_fam"][fam][m_st] += 1
            acc["pairs_rec_by_fam"][fam] += m_rec * (m_rec - 1) // 2
            acc["pairs_st_by_fam"][fam] += m_st * (m_st - 1) // 2
            acc["pairs_st_cross_by_fam"][fam] += (m_st * (m_st - 1) // 2
                                                  - sum(d * (d - 1) // 2 for d in Ds))
            Rmax = max(Rs)
            if Rmax < F(1, 2):
                acc["below_slots"][fam] += 1
                acc["below_mult_rec"][fam][m_rec] += 1
                acc["below_mult_st"][fam][m_st] += 1

            # ---- the ceiling annexe, both units
            ck = "%s|%s" % (fam, Rmax)
            if m_rec > acc["ceil_rec"][ck]:
                acc["ceil_rec"][ck] = m_rec
            if m_st > acc["ceil_st"][ck]:
                acc["ceil_st"][ck] = m_st
            if len(set(Rs)) == 1:
                if m_rec > acc["ceil_rec_single"][ck]:
                    acc["ceil_rec_single"][ck] = m_rec
                if m_st > acc["ceil_st_single"][ck]:
                    acc["ceil_st_single"][ck] = m_st
                cap = ceiling_of(Rmax)
                if cap is not None and m_st > cap:
                    acc["ceil_over"] += 1
            else:
                acc["mixed_groups"][ck] += 1
                for Rv in Rs:
                    acc["mixed_fields"]["%s|%s" % (ck, Rv)] += 1

            # ---- equation (mainbound), record unit and state unit
            if m_rec >= 2:
                srt = sorted(Rs, reverse=True)
                ok = any(pair_ok(srt[i], srt[j])
                         for i in range(len(srt)) for j in range(i + 1, len(srt)))
                if not ok:
                    acc["viol_rec"] += 1
                    if len(acc["viol_ex"]) < 10:
                        acc["viol_ex"].append(dict(unit="record", model=lab, kl=list(kl),
                                                   role=ro, R=[str(x) for x in Rs]))
            if m_st >= 2:
                # the state list repeats each record's R with multiplicity d, so a record
                # of d >= 2 contributes the within-record pair (R, R)
                stR = []
                for Rv, d in zip(Rs, Ds):
                    stR.extend([Rv] * d)
                srt = sorted(stR, reverse=True)
                ok = any(pair_ok(srt[i], srt[j])
                         for i in range(len(srt)) for j in range(i + 1, len(srt)))
                if not ok:
                    acc["viol_st"] += 1
                    if len(acc["viol_ex"]) < 10:
                        acc["viol_ex"].append(dict(unit="state", model=lab, kl=list(kl),
                                                   role=ro, R=[str(x) for x in Rs],
                                                   d=Ds))

            # ---- Delta = p - p'' for every same-role pair of distinct RECORDS
            if m_rec >= 2:
                for i in range(m_rec):
                    for j in range(i + 1, m_rec):
                        D = [F(x) - F(y) for x, y in zip(rs[i]["p_sh"], rs[j]["p_sh"])]
                        acc["d_pairs"] += 1
                        if not S5.in_lattice(D):
                            acc["d_notlambda"] += 1
                            if len(acc["d_bad_ex"]) < 10:
                                acc["d_bad_ex"].append(dict(kind="notlambda", model=lab,
                                                            kl=list(kl), role=ro))
                            continue
                        acc["d_d2"][str(S5.dot(D, D))] += 1
                        if S5.proj_norm2(D, smbasis) != 0:
                            acc["d_notperp"] += 1
                            if len(acc["d_bad_ex"]) < 10:
                                acc["d_bad_ex"].append(dict(kind="notperp", model=lab,
                                                            kl=list(kl), role=ro))


def build_smbasis(spec, roles, tY, allr, cells_first_psh):
    """W = span(SU(3)_c simple roots, SU(2)_L simple root, t_Y), recovered per model from
    the unbroken 4d gauge group (the intersection of the local groups over populated
    cells), the two factors being pinned by the requirement that the Q/d^c and Q/L records
    are non-trivial under them."""
    S5 = _S5
    inter = None
    for V in cells_first_psh:
        m = mask_of(V)
        inter = m if inter is None else (inter & m)
    if inter is None:
        return None
    R4 = [allr[i] for i in np.nonzero(inter)[0]]
    if not R4:
        return None
    Sp4 = S5.simple_roots(R4)
    if not Sp4:
        return None
    comps4, adj4 = S5._components(Sp4)
    qrec = [r for r in spec.values() if roles.get(r["label"]) in ("Q", "Qbar")]
    drec = [r for r in spec.values() if roles.get(r["label"]) in ("dc", "dcbar")]
    lrec = [r for r in spec.values() if roles.get(r["label"]) in ("L", "Hu")]
    f3, f2 = None, None
    for c in comps4:
        kind, rk = S5._classify(c, adj4, Sp4)
        so = [Sp4[i] for i in c]
        if kind == "A" and rk == 2 and (qrec or drec):
            if all(any(S5.dot(x["p_sh"], a) != 0 for a in so) for x in qrec + drec):
                f3 = so if f3 is None else "AMBIG"
        if kind == "A" and rk == 1 and (qrec or lrec):
            if all(any(S5.dot(x["p_sh"], a) != 0 for a in so) for x in qrec + lrec):
                f2 = so if f2 is None else "AMBIG"
    if isinstance(f3, list) and isinstance(f2, list):
        return list(f3) + list(f2) + [tuple(tY)]
    return None


SM_ROLES = {
    (3, 2, F(1, 6)): "Q", (-3, 1, F(-2, 3)): "uc", (-3, 1, F(1, 3)): "dc",
    (1, 2, F(-1, 2)): "L", (1, 1, F(1)): "ec", (1, 2, F(1, 2)): "Hu",
    (-3, 2, F(-1, 6)): "Qbar", (3, 1, F(2, 3)): "ucbar", (3, 1, F(-1, 3)): "dcbar",
    (1, 1, F(-1)): "ecbar",
}


def assign_roles(spec, tY, a3, a2, sY0, sY, sC):
    """U1 + U7: ungated role and hidden dimension for every left-chiral record."""
    S5 = _S5
    roles, hid = {}, {}
    for lab_f, r in spec.items():
        if r["rep"] is None or r["chir"] != "l" or r["p_sh"] is None:
            continue
        if a3 >= len(r["rep"]) or a2 >= len(r["rep"]):
            continue
        Y = sY * sY0 * S5.dot(tY, r["p_sh"])
        c3 = r["rep"][a3]
        c3 = (sC * c3) if abs(c3) != 1 else 1
        c2 = abs(r["rep"][a2])
        d = 1
        for kk, dd in enumerate(r["rep"]):
            if kk not in (a3, a2):
                d *= abs(dd)
        roles[lab_f] = SM_ROLES.get((c3, c2, Y))
        hid[lab_f] = d
    return roles, hid


def worker(args):
    path, geom, skip = args
    if _S5 is None:
        _init()
    G = GEOM[geom]
    sys.path.insert(0, os.path.join(CAMP, "code", G["lib"]))
    Z = __import__(G["mod"])
    S5 = _S5
    allr = S5.e8xe8_roots()
    skip = set(skip)
    # THREE convention variants of the same models, so that every difference between the
    # printed table and the uniform one is attributable to exactly one cause:
    #   current -- as run, first printed configuration, and the hidden-singlet GATE on
    #              Z3xZ6_1_1 only: this is the convention pair Table 1 of paperA_v5 prints
    #   asrun   -- as run, first printed configuration, UNGATED on both: isolates the gate
    #   uniform -- U1-U8 in full: ungated, first CLOSING configuration, distinct models
    curr, asrun, unif = newacc(), newacc(), newacc()
    txt = gzip.open(path, "rt", errors="replace").read()

    for lab, blk in Z.split_models(txt):
        pre = blk.split("Analyzing")[0]
        try:
            if geom == "Z2xZ4_1_1":
                spec = Z.parse_spectrum(pre, M=G["M"], N=G["N"])
                factors = Z.parse_gauge_factors(blk)
                cfgs = Z.parse_smconfig(blk)
                if not cfgs or not spec:
                    for _a in (curr, asrun, unif):
                        _a["skipped"]["no_smconfig"] += 1
                    continue
                found = None
                for ci, (tYc, tbl) in enumerate(cfgs):
                    if tYc is None:
                        continue
                    printed = defaultdict(int)
                    for row in tbl:
                        printed[(row[0], row[1], row[2])] += row[3]
                    pos = Z.observable_positions(spec, tYc, factors, dict(printed))
                    if pos is not None:
                        found = (ci, tYc, pos[0], pos[1], pos[2])
                        break
                if found is None:
                    for _a in (curr, asrun, unif):
                        _a["skipped"]["closure_failed"] += 1
                    continue
                ci, tY, a3, a2, sY0 = found
                _r, diag = Z.role_census(spec, tY, a3, a2, sY0)
            else:
                spec = Z.parse_spectrum(pre)
                gauge = Z.parse_gauge(blk)
                if gauge["tY"] is None or not spec:
                    for _a in (curr, asrun, unif):
                        _a["skipped"]["no_hypercharge"] += 1
                    continue
                pos = Z.observable_positions(spec, gauge["tY"], gauge)
                if pos is None:
                    for _a in (curr, asrun, unif):
                        _a["skipped"]["closure_failed"] += 1
                    continue
                a3, a2, sY0, tY, ci = pos
                _r, diag = Z.role_census(spec, tY, a3, a2, sY0)
        except Exception as e:                                          # noqa: BLE001
            for _a in (curr, asrun, unif):
                _a["skipped"]["parse_error:%s" % type(e).__name__] += 1
            continue

        sY, sC = diag["sY_selected"], diag["sC_selected"]
        roles, hid = assign_roles(spec, tY, a3, a2, sY0, sY, sC)

        # V-BR: re-derive the branch UNGATED and STATE-weighted
        best, bestkey = None, None
        for tsY in (1, -1):
            for tsC in (1, -1):
                ro2, hd2 = assign_roles(spec, tY, a3, a2, sY0, tsY, tsC)
                netQ = sum(hd2[k] for k, v in ro2.items() if v == "Q") - \
                    sum(hd2[k] for k, v in ro2.items() if v == "Qbar")
                netU = sum(hd2[k] for k, v in ro2.items() if v == "uc") - \
                    sum(hd2[k] for k, v in ro2.items() if v == "ucbar")
                if best is None or (netQ, netU) > best:
                    best, bestkey = (netQ, netU), (tsY, tsC)
        br_bad = (bestkey != (sY, sC))

        cells_first = {}
        for _l, r in spec.items():
            if r["p_sh"] is None or r["chir"] != "l":
                continue
            k = (tuple(r["kl"]), tuple(r["fp"]) if r["fp"] else None)
            cells_first.setdefault(k, r["p_sh"])
        smb = build_smbasis(spec, roles, tY, allr, list(cells_first.values()))
        if smb is None:
            for _a in (curr, asrun, unif):
                _a["smbasis_fail"] += 1
            continue

        # the gate is exactly `hidden dimension == 1`
        roles_g = {k: (v if hid.get(k, 1) == 1 else None) for k, v in roles.items()}

        # "as run" = the embedding the published census actually used: the FIRST printed
        # configuration on Z2xZ4_1_1 (whence its 615 failures), and the first CLOSING one
        # on Z3xZ6_1_1, whose `observable_positions` already walks the printed blocks.
        asrun_ok = (ci == 0) if geom == "Z2xZ4_1_1" else True
        targets = []
        if asrun_ok:
            targets.append((asrun, roles))
            if geom == "Z3xZ6_1_1":
                targets.append((curr, roles_g))
        if lab not in skip:
            targets.append((unif, roles))
        memo = {}
        for acc, rr in targets:
            acc["models"] += 1
            acc["embed_index"][ci] += 1
            if ci > 0:
                acc["models_recovered"] += 1
            acc["kY"][str(2 * S5.dot(tY, tY))] += 1
            if br_bad:
                acc["br_disagree"] += 1
                if len(acc["br_disagree_ex"]) < 10:
                    acc["br_disagree_ex"].append(dict(model=lab, lib=[sY, sC],
                                                      ungated_state=list(bestkey)))
            census_model(geom, lab, spec, rr, hid, smb, acc, memo)
    return curr, asrun, unif


def labels_of(args):
    path, geom = args
    sys.path.insert(0, os.path.join(CAMP, "code", GEOM[geom]["lib"]))
    Z = __import__(GEOM[geom]["mod"])
    txt = gzip.open(path, "rt", errors="replace").read()
    return path, [lab for lab, _b in Z.split_models(txt)]


def jsonable(a):
    o = {}
    for k, v in a.items():
        if isinstance(v, set):
            o[k] = len(v)
        elif isinstance(v, Counter):
            o[k] = {str(x): y for x, y in sorted(v.items(), key=lambda t: str(t[0]))}
        elif isinstance(v, defaultdict):
            if v and isinstance(next(iter(v.values())), Counter):
                o[k] = {str(x): {str(p): q for p, q in sorted(d.items(),
                                                             key=lambda t: str(t[0]))}
                        for x, d in v.items()}
            else:
                o[k] = {str(x): y for x, y in sorted(v.items(), key=lambda t: str(t[0]))}
        else:
            o[k] = v
    return o


def main():
    geom = sys.argv[1]
    files = sorted(glob.glob(os.path.join(GEOM[geom]["census"], "*.txt.gz")))
    if len(sys.argv) > 2:
        files = files[:int(sys.argv[2])]
    t0 = time.time()

    # ---- pass 1 (U4): measure the duplicate-label map; do not assume it
    with mp.Pool(6) as pool:
        labmap = dict(pool.map(labels_of, [(f, geom) for f in files]))
    seen, skipmap, dups = set(), {}, Counter()
    for f in files:
        sk = []
        for lab in labmap[f]:
            if lab in seen:
                sk.append(lab)
                dups[lab] += 1
            else:
                seen.add(lab)
        skipmap[f] = sk
    n_blocks = sum(len(v) for v in labmap.values())
    print("pass1: %d blocks, %d distinct labels, %d duplicate occurrences in %d labels"
          % (n_blocks, len(seen), sum(dups.values()), len(dups)))

    # ---- pass 2
    with mp.Pool(6) as pool:
        res = pool.map(worker, [(f, geom, skipmap[f]) for f in files])
    C, A, U = newacc(), newacc(), newacc()
    for c, a, u in res:
        merge_into(C, c)
        merge_into(A, a)
        merge_into(U, u)
    if geom == "Z2xZ4_1_1":
        # the published census is already ungated on this geometry, so `current` and
        # `asrun` are the same object and only one of them is computed
        C = A
    out = dict(geometry=geom, n_segments=len(files), seconds=round(time.time() - t0, 1),
               blocks=n_blocks, distinct_labels=len(seen),
               duplicate_occurrences=sum(dups.values()),
               duplicate_labels=sorted(dups)[:40],
               current=jsonable(C), asrun=jsonable(A), uniform=jsonable(U))
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "uniform_%s.json" % geom), "w") as fh:
        json.dump(out, fh, indent=1)
    for nm, acc in (("CURRENT", C), ("ASRUN", A), ("UNIFORM", U)):
        print("--- %s %s: models %d (recovered %d) fields10 rec %d st %d | "
              "pairs10 rec %d | fields5 rec %d st %d | pairs5 rec %d | viol %d/%d"
              % (nm, geom, acc["models"], acc["models_recovered"],
                 acc["f_rec_by_fam"]["10"], acc["f_st_by_fam"]["10"],
                 acc["pairs_rec_by_fam"]["10"],
                 acc["f_rec_by_fam"]["5bar"], acc["f_st_by_fam"]["5bar"],
                 acc["pairs_rec_by_fam"]["5bar"], acc["viol_rec"], acc["viol_st"]))
        print("    same-role record pairs total %d ; w2 %s ; fp_len %s ; I1_bad %d ; "
              "notperp %d ; notlambda %d ; d2 %s ; br_disagree %d"
              % (acc["pairs_rec_by_fam"]["10"] + acc["pairs_rec_by_fam"]["5bar"],
                 {k: dict(v) for k, v in acc["w2_by_fam"].items()},
                 dict(acc["fp_len"]), acc["I1_bad"], acc["d_notperp"],
                 acc["d_notlambda"], dict(acc["d_d2"]), acc["br_disagree"]))


if __name__ == "__main__":
    main()

"""z3xz6_sweep -- shared library for the 4,612-model `Z3xZ6_1_1` ladder sweep.

Spec: `docs/z3xz6_sweep_task.md`.  Authorisation: `results/gate_z3xz6_foundations.md`
(PASS-WITH-CORRECTIONS) under its FOUR BINDING REQUIREMENTS, which are implemented here
and named at the point of use:

  R1  CHANNEL RETARGETING.  The maximal 4x4 nine-factor block channel `T(0,1)T(0,2)T(0,3)`
      and all six all-`T(0,l)` channels are EMPTY BY GEOMETRY (the `R_2 = 0` theorem:
      every leg is an N=2 sector whose fixed torus is plane 2, integer oscillator moding
      against a massless budget < 1 in all 17 twisted sectors forces N^2_L = Nbar^2_L = 0,
      so `sum R_2 = 0 != 1 mod 3` for every assignment, in both conjugation orientations).
      They are PRE-DECLARED empty and used only as a PARSER CHECK (`EMPTY_BY_GEOMETRY`);
      finding them empty is never reported as a result.  The sweep runs on the ALLOWED
      channels, led by the gate-found carrier `T(0,3)T(1,5)T(2,4)`.
  R2  OSCILLATOR DISCIPLINE.  `S_cl` carries NO oscillator dependence (1301.2322
      (2.20)-(2.26): oscillators enter the prefactor and the power in (2.26), never the
      exponent), so the ALPHABET is oscillator-INDEPENDENT and needs no extension.  The
      R-CHARGE accounting is oscillator-INCLUSIVE, `R_a = r_a + N^a_L - Nbar^a_L`, and the
      per-field oscillator numbers are parsed and carried.
  R3  GAMMA IN THE CORRECTED PAIR FORMALISM.  `gamma = (gamma^theta, gamma^omega)`, a
      CHARACTER of the stabiliser (orbifolder `cstate.cpp:95-130`), read off the printed
      `gamma phases` slots 0 and 1 -- never modelled.  Denominators are asserted in
      {1,2,3} (the exponent of Z3xZ6 is 6, so 6*gamma in Z; ninths are structurally
      impossible).  The three-branch tag (blind / correlator `-gamma` / LEEFT `+gamma`)
      travels with every C(S_min).  The SIX UNCOVERED NON-GENERATOR ELEMENTS of the owned
      non-prime lift (gate DEFECT-D-2) are handled by `UNCOVERED_ELEMENTS` below.
  R4  HYPERCHARGE FIRST.  `hypercharge_z36.py` runs the dual route with EVERY defining
      condition (C1..C8) before any verdict row.

CONVENTIONS, STATED ONCE (CLAUDE.md 5.19(i)), each with its invariant asserted in
`selftest.py` BEFORE any count is printed:

  * sectors are PAIRS `(k, l)` in Z_3 x Z_6.  The `%6`-scalar sector arithmetic of the
    Z_N modules (`zsw_lib`, `gza_census`) is NOT portable here and is not imported.
  * orbifolder prints `right-moving q_sh` with slot 0 the 4d helicity slot and slots 1..3
    the three complex planes; the literature H-momentum is `r_a = -q_sh_a`.  MEASURED on
    this geometry, 243/243 fields, `r_a = 1 - frac(phi_a) (mod 1)` and `= 0` on a fixed
    torus (foundations addendum A2) -- so this is not a choice that happened to work.
  * a printed oscillator exponent `^(-a)` is the BARRED oscillator `Nbar^a_L`; `^a` is
    unbarred `N^a_L`.
  * lattice vectors 1,2 span plane 1; 3,4 plane 2; 5,6 plane 3 (ASSERTED against the
    geometry file, never assumed).
  * MEASURED left-chiral basis (foundations addendum A4, derived == measured exactly):
    the 14 `_l` sectors are U + the nine N=2 sectors + T(1,5), T(2,3), T(2,4), T(2,5);
    `T(1,1), T(1,2), T(1,3), T(2,1)` print only as CPT conjugates.  A zero count on one
    of those four labels is LABELLING, NOT ABSENCE.
  * hypercharge sign `sY` and colour orientation `sC` (3 vs 3bar) are INDEPENDENT
    per-model conventions; all four branches are computed and the physical one selected
    by net Q = +3.  (The Z6-I sweep found 7 of its first 13 models printing Q as
    (-3,2,+1/6); MSSM0 of THIS catalogue prints it the same way.)

Anti-claim discipline enforced in code: `CATALOGUE_TAG` vs `GEOMETRY_TAG` are stamped on
every emitted statement, and `RATIO_DISCLOSURE` ([G-3]) plus `PREDICTIVITY_TAG` ([G-1])
travel with every ladder statement.
"""

from __future__ import annotations

import cmath
import gzip
import json
import os
import re
import sys
from collections import defaultdict
from fractions import Fraction as F
from itertools import permutations, product

HERE = os.path.dirname(os.path.abspath(__file__))
MGC = os.path.abspath(os.path.join(HERE, "..", ".."))
REPO = os.path.abspath(os.path.join(MGC, ".."))
sys.path.insert(0, os.path.join(MGC, "code", "z3xz6"))

import z36lib as ZL                                                    # noqa: E402
import rules_z36 as RZ                                                 # noqa: E402
import alphabet_z36 as AZ                                              # noqa: E402

DATA = os.path.join(MGC, "data", "z3xz6_sweep")
RAW_CENSUS = os.path.join(DATA, "raw", "census")
RAW_COUP = os.path.join(DATA, "raw", "coup")

G = ZL.Geom()
M3, N6 = G.M, G.N                                     # 3, 6
PLANE_ORDER = [G.plane_order(a) for a in (1, 2, 3)]   # (6, 3, 6)

# ---------------------------------------------------------------- tags that TRAVEL
GEOMETRY_TAG = ("GEOMETRY statement: a property of the Z3xZ6_1_1 space group, true for "
                "every model built on it.")
CATALOGUE_TAG = ("CATALOGUE statement: a property of the 4,612 PVW `Z3xZ6_1_1` models; "
                 "NOT a statement about the geometry or about all admissible embeddings.")
RATIO_DISCLOSURE = ("[G-3] in-band index-level ratio exactly 2.0 and band-robust, but the "
                    "UNRESTRICTED index-level ratio set over the Z_M x Z_N class is "
                    "{3/2, 2, 3}, and 3/2 misses the pre-declared band floor 1.6 by 0.1.")
PREDICTIVITY_TAG = ("[G-1] Z3xZ6_1_1 is roto-translation free; its [1,3] exchange "
                    "automorphism tau(k,l) = (k, 4k+5l) acts on 66/66 channels and is an "
                    "involution preserving plane 2, so the equal-eps locus Im T_1 = Im T_3 "
                    "is a SYMMETRY POINT of the space group: a ladder realised here is a "
                    "PREDICTION, not a fit.")
GAMMA_BRANCH_TAG = ("gamma-law sign branch: `blind` (gamma dropped), `correlator` (-gamma) "
                    "and `LEEFT` (+gamma) are all carried; the printed gamma set is not "
                    "self-negative (1/3 != 2/3), so the Z4-style structural collapse of the "
                    "branch question does NOT occur here.  Never averaged across branches.")
RULEGAP_G1 = ("RULE-GAP G-1: the non-prime gamma correction in planes 1 and 3 has NO "
              "published two-generator form (1301.2322 fn.1 declines Z_N x Z_M).  The "
              "programme's lift is emitted PER GENERATOR while EIGHT group elements carry a "
              "non-prime plane and SIX of them are not powers of omega (gate DEFECT-D-2).")
RULEGAP_G2 = ("RULE-GAP G-2: 1301.2322 excludes discrete Wilson lines and 4,588 of the "
              "4,612 catalogue models carry one (order 3, plane 2).  Catalogue-wide "
              "exposure, larger than Z6-I's or Z4's.")

# ---- R1: the pre-declared empty-by-geometry channels (PARSER CHECK ONLY) ----------
EMPTY_BY_GEOMETRY = tuple(sorted(
    tuple(sorted(t)) for t in
    [((0, a), (0, b), (0, c)) for a in range(6) for b in range(6) for c in range(6)
     if (a + b + c) % 6 == 0 and 0 not in (a, b, c)]))
EMPTY_BY_GEOMETRY_REASON = (
    "R_2 = r_2 + N^2_L - Nbar^2_L = 0 identically for every massless field of a T(0,l) "
    "sector (all five are N=2 with fixed torus in PLANE 2; integer oscillator moding "
    "against a massless budget < 1 in all 17 twisted sectors forces N^2_L = Nbar^2_L = 0), "
    "so sum_legs R_2 = 0 != 1 mod 3 for EVERY assignment of massless fields, in BOTH "
    "conjugation orientations and at every oscillator level.  GEOMETRY-level theorem "
    "(gate T3.a/T3.b).  Pre-declared: emptiness here is a parser check, never a finding.")

# ---- R1: the carrier and its two secondaries -------------------------------------
CARRIER = ((0, 3), (1, 5), (2, 4))
SECONDARY = (((0, 2), (1, 5), (2, 5)), ((0, 4), (1, 5), (2, 3)))

# ---- R3: the six non-generator elements carrying a non-prime plane ---------------
NONPRIME_ELEMENTS = ((0, 1), (0, 5), (1, 1), (1, 3), (1, 5), (2, 1), (2, 3), (2, 5))
UNCOVERED_ELEMENTS = tuple(g for g in NONPRIME_ELEMENTS
                           if g not in tuple((0, l) for l in range(6)))
GAMMA_DENOMINATORS_ALLOWED = (1, 2, 3, 6)

# ---- the MEASURED left-chiral basis ----------------------------------------------
LEFT_CHIRAL_SECTORS = ((0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 0), (1, 4),
                       (1, 5), (2, 0), (2, 2), (2, 3), (2, 4), (2, 5))
RIGHT_CHIRAL_LABELS = ((1, 1), (1, 2), (1, 3), (2, 1))


def sector_label(g):
    return "U" if tuple(g) == (0, 0) else "T(%d,%d)" % (g[0] % 3, g[1] % 6)


def channel_name(tri):
    return "".join(sector_label(g) for g in tri)


# ===================================================================== parsing
def _rat(t):
    t = t.strip()
    return F(t) if t else F(0)


def _vec_pairs(s):
    out = []
    for grp in re.findall(r"\(([^)]*)\)", s):
        out.extend(_rat(x) for x in grp.split(","))
    return out


_MARK = re.compile(r"/>\s*/(MSSM\d+)>")
_OSC = re.compile(r"a\(\s*(-?\d+)\s*\)|a_\(\s*(-?\d+)\s*\)|\^\(?(-?\d+)\)?")


def split_models(text):
    """[(label, block)] -- one block per model, cut at the `/> /MSSMn>` prompt marks."""
    marks = [(m.start(), m.group(1)) for m in _MARK.finditer(text)]
    out = []
    for i, (pos, lab) in enumerate(marks):
        end = marks[i + 1][0] if i + 1 < len(marks) else len(text)
        out.append((lab, text[pos:end]))
    return out


def read_raw(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", errors="replace") as fh:
        return fh.read()


_OSC_TOK = re.compile(r"(-?)a_\(\s*(-?\d+(?:/\d+)?)\s*\)\^(?:\((-\d)\)|(\d))")


def parse_oscillators(s):
    """`no` or a list of `a_(+/-i)`-style tokens -> (N_L, Nbar_L) per plane (1-based).

    R2 / the CLAUDE.md 5.19 conjugation family, handled by MEASUREMENT not by choice.
    The token grammar is the gate's certified one (`gate_z3xz6/gate_settle.py::OSC`):
    orbifolder prints `-a_( mode )^(-n)` or `-a_( mode )^n`, the trailing index n being
    the plane and the parenthesised-negative form one of the two conjugate towers.
    WHICH tower is the unbarred `N^n_L` is a labelling choice; this library does not
    assume it.  `selftest.py` S10 fixes it by MEASUREMENT -- it requires the offline
    plane-2 R-charge law to reproduce orbifolder's OWN enforced selection on every one of
    the 908 + 619 couplings returned under the `_Rsym_prime` patch -- and reports both
    branches with their violation counts.  A plane index 0 or > 3 raises.
    """
    NL, NbL = {1: 0, 2: 0, 3: 0}, {1: 0, 2: 0, 3: 0}
    s = (s or "no").strip()
    if s in ("no", ""):
        return NL, NbL, s
    hits = _OSC_TOK.findall(s)
    if not hits:
        raise ValueError("unparsable oscillator string %r" % s)
    for _sgn, _mode, neg, pos in hits:
        a = abs(int(neg if neg else pos))
        if a == 0 or a > 3:
            raise ValueError("oscillator plane index %d in %r" % (a, s))
        (NbL if neg else NL)[a] += 1
    return NL, NbL, s


def parse_rep_token(x):
    x = x.strip()
    m = re.match(r"^(-?\d+)(?:_([vcs]))?$", x)
    if not m:
        raise ValueError("unparsable rep token %r" % x)
    return int(m.group(1)), m.group(2)


def parse_spectrum(block):
    """{F_n: record} from `print(F) with internal information`."""
    heads = [(m.start(), m.group(1))
             for m in re.finditer(r"(?m)^\s{2,8}(F_\d+)\s*$", block)]
    out = {}
    for i, (pos, lab) in enumerate(heads):
        end = heads[i + 1][0] if i + 1 < len(heads) else min(len(block), pos + 3000)
        b = block[pos:end]

        def grab(pat):
            m = re.search(pat, b)
            return m.group(1) if m else None

        sec = grab(r"sector \(k,l\)\s*:\s*\(([^)]*)\)")
        if sec is None:
            continue
        kl = [int(x) for x in sec.split(",")]
        kl = (kl[0] % 3, kl[1] % 6)
        fps = grab(r"fixed point n_a\s*:\s*\(([^)]*)\)")
        fp = tuple(int(x) for x in fps.split(",")) if fps else None
        qs = tuple(_rat(x) for x in
                   grab(r"right-moving q_sh\s*:\s*\(([^)]*)\)").split(","))
        gam = tuple(_rat(x) for x in
                    grab(r"gamma phases\s*:\s*\(([^)]*)\)").split(","))
        psh_s = grab(r"left-moving p_sh\s*:\s*(.*)")
        psh = tuple(_vec_pairs(psh_s)) if psh_s else None
        repm = re.search(
            r"representation\s*:\s*\(([^)]*)\)_([lr])\s+U\(1\):\s*\(([^)]*)\)", b)
        rep = tags = chir = q = None
        if repm:
            parsed = [parse_rep_token(x) for x in repm.group(1).split(",")]
            rep = tuple(p[0] for p in parsed)
            tags = tuple(p[1] for p in parsed)
            chir = repm.group(2)
            q = [_rat(x) for x in repm.group(3).split(",")]
        NL, NbL, osc_raw = parse_oscillators(grab(r"oscillators\s*:\s*(.*)"))
        out[lab] = dict(label=lab, kl=kl, fp=fp, qsh=qs, gamma=gam, p_sh=psh,
                        rep=rep, rep_tags=tags, chir=chir, q=q, osc_raw=osc_raw,
                        NL=(NL[1], NL[2], NL[3]), NbL=(NbL[1], NbL[2], NbL[3]))
    return out


_WPAT = re.compile(r"\n\s+W\s*=\s*(.*?)(?:\n\s*\n|\n/)", re.S)


def _terms(txt):
    out = []
    for t in " ".join(txt.split()).split("+"):
        legs = re.findall(r"b?F_\d+", t)
        if legs:
            out.append(tuple(legs))
    return out


def parse_W_all(block):
    """EVERY `W = ` print in the block, in order.

    `create coupling(F F F)` forks a child that inherits the redirected stdin, so in a
    MULTI-model orbifolder run the child re-executes the tail of the command stream and
    each model prints its superpotential TWICE.  The two prints are NOT always equal: a
    model whose coupling creation had not finished when the first `print superpotential`
    executed prints `W = 0`, and only the later print carries the real superpotential
    (measured: 1 of 28 models in one wait=12 segment).  Taking the FIRST print -- the
    obvious reading, and the one a SINGLE-model pilot cannot distinguish because no
    duplicate exists there -- would record a SPURIOUS ABSENCE, which is the one error
    class a sweep about absence cannot afford.
    """
    return [_terms(m.group(1)) for m in _WPAT.finditer(block)]


def parse_W(block):
    """The model's superpotential: the LONGEST of its prints (see `parse_W_all`).

    Longest rather than last, because either print may be the truncated one, and a
    truncated print is a PREFIX of the real superpotential (orbifolder appends terms as
    it finds them) rather than a different set.  `sweep.py` additionally flags every
    model whose prints disagree; those models are re-run at a long wait and the re-run
    supersedes.
    """
    alls = parse_W_all(block)
    return max(alls, key=len) if alls else []


_VEC16 = re.compile(r"\(([^()]*)\),\s*\(([^()]*)\)")
_SMROW = re.compile(r"^\s*(\d+)\s*\(\s*(-?\d+),\s*(-?\d+)\)_(\w+)\s+U\(1\)\s*:\s*"
                    r"\(\s*(-?[\d/]+)\s*\)")


def _vec16(line):
    m = _VEC16.search(line)
    if not m:
        return None
    a = [_rat(x) for x in m.group(1).split(",")]
    b = [_rat(x) for x in m.group(2).split(",")]
    if len(a) != 8 or len(b) != 8:
        return None
    return a + b


def parse_gauge(block):
    """t_Y, U(1) generators, simple roots, SU(5) origin, anomaly traces, SMConfig table.

    Segmented per `Hypercharge generator` occurrence and the LAST complete (tY, su5) pair
    used, so a model carrying several vev-configurations can never pair one config's
    hypercharge with another's SU(5) roots (the CLAUDE.md 5.19 swap shape).
    """
    d = {"u1": [], "roots": [], "su5": [], "tY": None, "sm": [], "trQ": {},
         "gauge_test": None, "gauge_sm": None, "n_hyper_blocks": 0, "configs": []}
    mode = None
    cur_tY, cur_su5, blocks = None, [], []
    for ln in block.split("\n"):
        s = ln.strip()
        if s.startswith('Gauge group in vev-configuration "TestConfig1"'):
            d["gauge_test"] = s.split(":", 1)[1].strip()
        if s.startswith('Gauge group in vev-configuration "SMConfig'):
            d["gauge_sm"] = s.split(":", 1)[1].strip()
            d["sm"] = []
            mode = "sm"
            continue
        if s.startswith("Gauge group in vev-configuration"):
            d["configs"].append(s.split(":", 1)[0].split('"')[1])
            if mode == "sm":
                mode = None
        m = re.match(r"^tr Q_(\d+)\s*=\s*(-?\d+)/(\d+)\s*$", s)
        if m:
            d["trQ"][int(m.group(1))] = F(int(m.group(2)), int(m.group(3)))
        if s.startswith("U(1) generators"):
            mode = "u1"
            continue
        if s.startswith("Simple roots"):
            mode = "roots"
            continue
        if s.startswith("SU(5) origin"):
            mode = "su5"
            continue
        if s.startswith("Hypercharge generator"):
            if cur_tY is not None:
                blocks.append((cur_tY, cur_su5))
            cur_tY, cur_su5 = None, []
            d["n_hyper_blocks"] += 1
            mode = "tY"
            continue
        if mode in ("u1", "roots", "su5", "tY"):
            v = _vec16(s)
            if v is None:
                if s and not s.startswith("("):
                    mode = None
                continue
            if mode == "tY":
                cur_tY = v
                mode = None
            elif mode == "su5":
                cur_su5.append(v)
            else:
                d[mode].append(v)
        elif mode == "sm":
            mm = _SMROW.match(ln)
            if mm:
                d["sm"].append({"mult": int(mm.group(1)), "d3": int(mm.group(2)),
                                "d2": int(mm.group(3)), "chir": mm.group(4),
                                "Y": _rat(mm.group(5))})
            elif s.startswith("/") or s.startswith("Do you"):
                mode = None
    if cur_tY is not None:
        blocks.append((cur_tY, cur_su5))
    if blocks:
        d["tY"], d["su5"] = blocks[-1]
    d["tY_blocks"] = blocks
    d["n_hyper_blocks"] = len(blocks)
    return d


def parse_banner(b):
    """`SU(3)_C x SU(2)_L x [SU(2)] and [SO(8)] and ...` -> [('SU',3), ...] in order."""
    out = []
    if not b:
        return out
    head = b.split(" and U(1)")[0]
    for tok in re.findall(r"(SU|SO|SP|E)[_]?\((\d+)\)", head):
        if tok[0] == "SU" and tok[1] == "1":
            continue
        out.append((tok[0], int(tok[1])))
    return out


# ============================================================ gamma (R3)
def gamma_pair(rec):
    """(gamma^theta, gamma^omega) -- slots 0 and 1 of the printed `gamma phases`.

    The 8 printed slots match the geometry file's 2 twist + 6 shift space-group
    generators (verified in the foundations pass), so slots 0,1 are the twist pair.
    """
    g = rec["gamma"]
    return (g[0], g[1])


def gamma_of(rec, h):
    """The character value gamma_alpha(h) = k*gamma^theta + l*gamma^omega for h=(k,l)."""
    gt, gw = gamma_pair(rec)
    return h[0] * gt + h[1] * gw


def check_gamma_denominators(spec):
    """Assert 6*gamma in Z on every printed twist gamma.

    The exponent of Z3 x Z6 is 6 and gamma is a CHARACTER of a subgroup, so every value
    satisfies 6*gamma in Z; denominator 9 is structurally impossible (gate T2.a, which
    is why the foundations' ten-value scalar set was retracted).  Raises on violation.
    """
    dens = set()
    for r in spec.values():
        for x in gamma_pair(r):
            if (6 * x).denominator != 1:
                raise ValueError("gamma %s violates 6*gamma in Z (field %s)"
                                 % (x, r["label"]))
            dens.add(x.denominator)
        if (3 * r["gamma"][0]).denominator != 1:
            raise ValueError("gamma^theta %s violates 3*gamma in Z" % r["gamma"][0])
    return sorted(dens)


# ============================================================ R-charge (R2)
OSC_SIGN = 1          # fixed by MEASUREMENT in selftest.py S10, not by assumption


def R_of(rec, osc_sign=None):
    """R_a = r_a + OSC_SIGN * (N^a_L - Nbar^a_L), oscillator-INCLUSIVE (R2), from the
    PRINTED q_sh.  r_a = -q_sh_a with slot 0 the 4d helicity slot (measured 243/243)."""
    s = OSC_SIGN if osc_sign is None else osc_sign
    r = [-rec["qsh"][a] for a in (1, 2, 3)]
    return [r[i] + s * (rec["NL"][i] - rec["NbL"][i]) for i in range(3)]


def rcharge_law_pass(recs, osc_sign=None):
    """0705.3075 (15) per plane: sum_legs R_a = 1 mod N_a, N = (6,3,6).  Oscillator-
    inclusive.  Returns {plane: (sum, pass)} plus the all-planes verdict."""
    out = {}
    ok_all = True
    for a in range(3):
        s = sum(R_of(r, osc_sign)[a] for r in recs)
        ok = ((s - 1) % PLANE_ORDER[a]) == 0
        out["plane%d" % (a + 1)] = {"sum_R": str(s), "modulus": PLANE_ORDER[a],
                                    "pass": ok}
        ok_all = ok_all and ok
    out["all_planes_pass"] = ok_all
    return out


# ============================================================ hypercharge / roles
def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


SM_ROLES = {
    (3, 2, F(1, 6)): "Q",
    (-3, 1, F(-2, 3)): "uc",
    (-3, 1, F(1, 3)): "dc",
    (1, 2, F(-1, 2)): "L",
    (1, 1, F(1)): "ec",
    (1, 2, F(1, 2)): "Hu",
    (-3, 2, F(-1, 6)): "Qbar",
    (3, 1, F(2, 3)): "ucbar",
    (3, 1, F(-1, 3)): "dcbar",
    (1, 1, F(-1)): "ecbar",
}


def hidden_singlet(rep, a3, a2):
    return all(abs(d) == 1 for k, d in enumerate(rep) if k not in (a3, a2))


def observable_positions(spec, tY, gauge):
    """(a3, a2, sY, tY_used) identified by REBUILDING the printed SMConfig multiplicity
    table from the full left-chiral spectrum -- never by row position (the Z6-I closure,
    ported to (k,l) sectors).

    WHICH hypercharge generator.  A model may print SEVERAL `Hypercharge generator`
    blocks (56 of the 4,612 print two), and taking the last one positionally pairs a
    hypercharge with a table it does not belong to -- a silent conjugation/duality-family
    swap of exactly the CLAUDE.md 5.19 shape, which shows up as a closure FAILURE rather
    than as a wrong number only because the closure is checked.  The block is therefore
    SELECTED BY THE SAME CLOSURE: every printed block is tried and the one that
    reproduces the table exactly is returned, together with its index.
    """
    printed = defaultdict(int)
    for st in gauge["sm"]:
        if st["chir"] != "l":
            continue
        printed[(st["d3"], st["d2"], st["Y"])] += st["mult"]
    if not printed:
        return None
    left = [r for r in spec.values() if r["rep"] and r["chir"] == "l" and r["p_sh"]]
    fac = parse_banner(gauge["gauge_test"])
    cands = [b[0] for b in gauge.get("tY_blocks") or []] or ([tY] if tY else [])
    for bi, tv in enumerate(cands):
        for a3 in range(len(fac)):
            if fac[a3] != ("SU", 3):
                continue
            for a2 in range(len(fac)):
                if a2 == a3 or fac[a2] != ("SU", 2):
                    continue
                for sY in (1, -1):
                    tbl = defaultdict(int)
                    for r in left:
                        if a3 >= len(r["rep"]) or a2 >= len(r["rep"]):
                            tbl = None
                            break
                        hid = 1
                        for k, dd in enumerate(r["rep"]):
                            if k not in (a3, a2):
                                hid *= abs(dd)
                        c3 = r["rep"][a3] if abs(r["rep"][a3]) != 1 else 1
                        tbl[(c3, abs(r["rep"][a2]), sY * dot(tv, r["p_sh"]))] += hid
                    if tbl is not None and dict(tbl) == dict(printed):
                        return (a3, a2, sY, tv, bi)
    return None


def role_census(spec, tY, a3, a2, sY0):
    """{label: dict(role, Y, kl)} plus orientation diagnostics.

    sY (overall sign of Y) and sC (which SU(3) weight is called the fundamental) are
    INDEPENDENT per-model conventions.  All FOUR branches are computed; the physical one
    is selected by net Q = +3.  Scanning only the combined flip leaves the
    (Y fixed, colour flipped) branch unreachable -- the live defect the Z6-I sweep caught
    mid-phase, and MSSM0 of THIS catalogue prints its quark doublet as (-3, 2, +1/6),
    i.e. exactly in that branch.
    """
    branches, nets = {}, {}
    for sY in (1, -1):
        for sC in (1, -1):
            out = {}
            for lab, r in spec.items():
                if r["rep"] is None or r["chir"] != "l" or r["p_sh"] is None:
                    continue
                if a3 >= len(r["rep"]) or a2 >= len(r["rep"]):
                    continue
                Y = sY * sY0 * dot(tY, r["p_sh"])
                if not hidden_singlet(r["rep"], a3, a2):
                    out[lab] = dict(label=lab, role=None, Y=Y, kl=r["kl"])
                    continue
                c3 = r["rep"][a3]
                c3 = (sC * c3) if abs(c3) != 1 else 1
                c2 = abs(r["rep"][a2])
                out[lab] = dict(label=lab, role=SM_ROLES.get((c3, c2, Y)), Y=Y,
                                kl=r["kl"])
            branches[(sY, sC)] = out

    def net(o, pos, neg):
        return (sum(1 for v in o.values() if v["role"] == pos)
                - sum(1 for v in o.values() if v["role"] == neg))
    for key, o in branches.items():
        nets[key] = {"Q": net(o, "Q", "Qbar"), "uc": net(o, "uc", "ucbar"),
                     "dc": net(o, "dc", "dcbar"), "ec": net(o, "ec", "ecbar")}
    best = max(nets, key=lambda k: (nets[k]["Q"], nets[k]["uc"]))
    diag = dict(sY_selected=best[0], sC_selected=best[1], sY_from_closure=sY0,
                nets={"sY=%d,sC=%d" % k: v for k, v in nets.items()},
                net_Q_is_plus_three=(nets[best]["Q"] == 3),
                net_all_three=all(nets[best][k] == 3 for k in ("Q", "uc", "dc", "ec")))
    return branches[best], diag


# ---------------------------------------------------------------- chirality
def _max_matching(left, adj):
    match = {}

    def try_aug(u, seen):
        for v in adj.get(u, ()):
            if v in seen:
                continue
            seen.add(v)
            if v not in match or try_aug(match[v], seen):
                match[v] = u
                return True
        return False
    for u in left:
        try_aug(u, set())
    return match


def _missable(left, right, adj):
    """Vertices some MAXIMUM matching leaves unmatched (alternating-path sweep)."""
    match = _max_matching(left, adj)
    inv = {u: v for v, u in match.items()}
    radj = {}
    for u, vs in adj.items():
        for v in vs:
            radj.setdefault(v, []).append(u)
    miss = set()
    stack = [u for u in left if u not in inv]
    seenL, seenR = set(stack), set()
    miss |= set(stack)
    while stack:
        u = stack.pop()
        for v in adj.get(u, ()):
            if inv.get(u) == v or v in seenR:
                continue
            seenR.add(v)
            w = match.get(v)
            if w is not None and w not in seenL:
                seenL.add(w)
                miss.add(w)
                stack.append(w)
    stack = [v for v in right if v not in match]
    seenR2, seenL2 = set(stack), set()
    miss |= set(stack)
    while stack:
        v = stack.pop()
        for u in radj.get(v, ()):
            if match.get(v) == u or u in seenL2:
                continue
            seenL2.add(u)
            w = inv.get(u)
            if w is not None and w not in seenR2:
                seenR2.add(w)
                miss.add(w)
                stack.append(w)
    return miss, len(match)


def pg_allowed_pair(kl1, kl2):
    """The Z3 x Z6 point-group rule for a BILINEAR: sum k = 0 mod 3 AND sum l = 0 mod 6.

    Two generators, so BOTH congruences.  (The `%6`-only rule of the Z_N modules would
    admit mass terms this geometry forbids -- the non-portability flagged in the
    foundations' 7.3.)
    """
    return (kl1[0] + kl2[0]) % 3 == 0 and (kl1[1] + kl2[1]) % 6 == 0


def chirality(cen, pos_role, neg_role):
    """A state is CHIRAL iff at least one MAXIMUM matching leaves it unmatched.

    Two left-chiral states pair into a string-scale mass term iff their observable reps
    are conjugate, hypercharges opposite, both hidden non-abelian singlets, and their
    BILINEAR is point-group allowed.  Hidden/anomalous U(1) charges are soaked by
    UNTWISTED singlet vevs (kl = (0,0), which do not move the point-group sums);
    pairings needing TWISTED singlet insertions are NOT counted as decoupling.
    Generous / necessary-condition semantics.
    """
    left = [l for l, r in cen.items() if r["role"] == pos_role]
    right = [l for l, r in cen.items() if r["role"] == neg_role]
    adj = {u: [v for v in right if pg_allowed_pair(cen[u]["kl"], cen[v]["kl"])]
           for u in left}
    miss, nm = _missable(left, right, adj)
    return miss, len(left), len(right), nm


# ============================================================== localization
_PLANE_IDX = None


def plane_indices():
    """(i1, i2) lattice-vector indices spanning plane a, ASSERTED against the geometry
    file rather than assumed."""
    global _PLANE_IDX
    if _PLANE_IDX is not None:
        return _PLANE_IDX
    lv = G.f.lattice_vectors()
    out = []
    for a in range(3):
        c0, c1 = 2 * a, 2 * a + 1
        idx = [i for i, w in enumerate(lv)
               if abs(w[c0]) + abs(w[c1]) > 1e-9
               and all(abs(w[c]) < 1e-9 for c in range(6) if c not in (c0, c1))]
        assert len(idx) == 2, "plane %d does not factorize: %s" % (a + 1, idx)
        assert idx == [2 * a, 2 * a + 1], \
            "lattice-vector index basis is not (e1,e2|e3,e4|e5,e6): %s" % idx
        out.append(tuple(idx))
    _PLANE_IDX = tuple(out)
    return _PLANE_IDX


def leg_fixed_point(kl, fp, a):
    """f = (1 - R_a(g))^{-1} n|_plane a, exact rationals mod 1; None on a fixed torus."""
    R = G.rot(tuple(kl), a)
    if R is None or fp is None:
        return None
    A = ZL.msub_from_id(R)
    Ai = ZL.minv_rat(A)
    i1, i2 = 2 * (a - 1), 2 * (a - 1) + 1
    n = (F(fp[i1]), F(fp[i2]))
    w = ZL.mapply_rat(Ai, n)
    return (w[0] % 1, w[1] % 1)


_ORB_CACHE = {}


def plane_orbits(kl, a):
    key = (tuple(kl), a)
    if key not in _ORB_CACHE:
        _ORB_CACHE[key] = G.orbits_in_plane(tuple(kl), a)
    return _ORB_CACHE[key]


def leg_orbit(kl, fp, a):
    """(orbit_index, orbit_members) of the leg's plane-a fixed point; None on a torus."""
    f = leg_fixed_point(kl, fp, a)
    if f is None:
        return None
    obs = plane_orbits(kl, a)
    if obs is None:
        return None
    for i, o in enumerate(obs):
        if f in o:
            return (i, o, f)
    raise ValueError("fixed point %s not in any plane-%d orbit of %s" % (f, a, kl))


# ============================================================ NORM (R2) + C(S_min) (R3)
def entry_norm(tri, orbs, a):
    """NORM of one block entry in plane a via the CERTIFIED foundations machinery.

    R2: the alphabet is oscillator-INDEPENDENT (S_cl carries no oscillator dependence,
    1301.2322 (2.20)-(2.26)), so nothing about the legs' oscillator numbers enters here.
    Returns (Fraction|None, outcome).  None means UNDEFINED, never zero.
    """
    val, info, _per = AZ.entry_norm(a, list(tri), list(orbs))
    return val, info


def _group_images(kl, f, a):
    """{h: image of f in plane a under h in Z3xZ6} for the leg's plane-a fixed point."""
    out = {}
    for h in G.elements():
        R = G.rot(h, a)
        if R is None:
            out[h] = f
        else:
            y = ZL.mapply(R, f)
            out[h] = (y[0] % 1, y[1] % 1)
    return out


def csmin(tri, recs, a, branch="blind", box=4):
    """C(S) in plane a, in the corrected GAMMA-PAIR formalism (R3).

    C(S) := sum over space-group-allowed configurations at classical action S of
            exp(-2 pi i * sgn * sum_alpha gamma_alpha(h_alpha)),
    with `gamma_alpha(h) = k gamma^theta_alpha + l gamma^omega_alpha` the CHARACTER value
    (never a scalar theta-power: on a Z_M x Z_N geometry gamma is a pair, measured), and
    h_alpha running over COSET REPRESENTATIVES of the leg's stabiliser -- one
    representative per distinct plane-a image, so no 6d configuration is counted twice
    within the plane.

    branch: `blind` (sgn=0, gamma dropped), `correlator` (sgn=+1, the -gamma law),
            `LEEFT` (sgn=-1, the +gamma law).  Never averaged (GAMMA_BRANCH_TAG).

    Returns (S_min, n_configs, |C(S_min)|, C_is_zero, n_actions) or None if undefined
    (all leg orderings degenerate / space-group forbidden -- UNDEFINED is not zero).
    """
    sgn = {"blind": 0, "correlator": 1, "LEEFT": -1}[branch]
    ks = list(tri)
    fs = []
    for r in recs:
        f = leg_fixed_point(r["kl"], r["fp"], a)
        if f is None:
            return None
        fs.append(f)
    # a non-degenerate leg ordering (guard: reject 1 - theta^{n g} = 0)
    order = None
    for perm in permutations((0, 1, 2)):
        pf, _ = AZ.prefactor(a, ks[perm[0]], ks[perm[1]])
        if pf is not None:
            order = perm
            break
    if order is None:
        return None
    ks = [tri[i] for i in order]
    fs = [fs[i] for i in order]
    rr = [recs[i] for i in order]
    pf, _ = AZ.prefactor(a, ks[0], ks[1])

    Rs = [G.rot(tuple(k), a) for k in ks]
    if any(R is None for R in Rs):
        return None
    A = [ZL.msub_from_id(R) for R in Rs]
    A3inv = ZL.minv_rat(A[2])
    Binv = ZL.minv_rat(ZL.mmul(Rs[0], Rs[1]))

    # coset representatives: one h per distinct plane-a image
    reps = []
    for i in range(3):
        img = _group_images(ks[i], fs[i], a)
        seen, lst = {}, []
        for h, p in img.items():
            if p not in seen:
                seen[p] = h
                lst.append((h, p))
        reps.append(lst)
    by_point2 = defaultdict(list)
    for h, p in reps[2]:
        by_point2[p].append(h)

    acc = defaultdict(lambda: [0, 0j])
    for h0, P0 in reps[0]:
        for h1, P1 in reps[1]:
            for n1 in range(-box, box + 1):
                for n2 in range(-box, box + 1):
                    P1s = (P1[0] + n1, P1[1] + n2)
                    l1 = ZL.mapply(A[0], P0)
                    l2 = ZL.mapply(A[1], P1s)
                    rhs = (-(l1[0] + ZL.mapply(Rs[0], l2)[0]),
                           -(l1[1] + ZL.mapply(Rs[0], l2)[1]))
                    w = ZL.mapply_rat(Binv, rhs)
                    P2 = ZL.mapply_rat(A3inv, w)
                    P2m = (P2[0] % 1, P2[1] % 1)
                    if P2m not in by_point2:
                        continue
                    S = pf * AZ.rho(a, (P1s[0] - P0[0], P1s[1] - P0[1]))
                    for h2 in by_point2[P2m]:
                        ph = 0.0
                        if sgn:
                            ph = float(gamma_of(rr[0], h0) + gamma_of(rr[1], h1)
                                       + gamma_of(rr[2], h2))
                        e = acc[S]
                        e[0] += 1
                        e[1] += cmath.exp(-2j * cmath.pi * sgn * ph)
    if not acc:
        return None
    S = min(acc)
    n, c = acc[S]
    return (S, n, abs(c), abs(c) < 1e-9, len(acc))


# ================================================================ channels
_ALLOWED = None


def allowed_channels():
    """The certified ALLOWED-channel list: 0705.3075 (15) in all three planes,
    oscillator-free, MINIMAL H-momentum reading (measured, catalogue-wide).

    R1: the six all-T(0,l) channels are pre-declared EMPTY BY GEOMETRY and are asserted
    ABSENT from this list rather than merely happening not to appear.
    """
    global _ALLOWED
    if _ALLOWED is not None:
        return _ALLOWED
    out = []
    for tri in RZ.all_multisets():
        n_u = sum(1 for g in tri if tuple(g) == (0, 0))
        ok = False
        for assign in product((1, 2, 3), repeat=n_u):
            legs, it = [], iter(assign)
            for g in tri:
                legs.append(("U", next(it)) if tuple(g) == (0, 0) else tuple(g))
            if RZ.law_published_perplane(legs)["all_planes_pass"]:
                ok = True
                break
        if ok:
            out.append(tuple(tuple(t) for t in tri))
    for ch in EMPTY_BY_GEOMETRY:
        assert tuple(ch) not in [tuple(sorted(c)) for c in out], \
            "an all-T(0,l) channel appeared in the allowed list: %s" % (ch,)
    _ALLOWED = tuple(out)
    return _ALLOWED


def allowed_channels_twisted():
    """The allowed channels with NO untwisted leg -- the ones that can carry a
    localisation-based suppression ladder at all (an untwisted leg is DELOCALISED, so its
    plane-a orbit is undefined and no NORM exists on that channel)."""
    return tuple(c for c in allowed_channels()
                 if not any(tuple(t) == (0, 0) for t in c))


_LADDER = None


def ladder_channels():
    """Allowed channels that are LADDER-BEARING, defined on the TOTAL exponent.

    A rank-3 ladder needs three DISTINCT classical factors, i.e. >= 3 distinct values of
    the TOTAL suppression exponent summed over the channel's computable planes -- NOT >= 3
    letters in a single plane.  The distinction is load-bearing and is exactly what makes
    this geometry interesting: the carrier's per-plane alphabets are {0,1,3,4} in planes 1
    and 3, and it is their SUM that realises the nine-factor set {0,...,8}.  A filter
    written per plane would have been a different (and weaker) object.

    Computed from the certified alphabet machinery, never asserted.
    """
    global _LADDER
    if _LADDER is not None:
        return _LADDER
    out = {}
    for tri in allowed_channels():
        name = channel_name(tri)
        if any(tuple(t) == (0, 0) for t in tri):
            out[name] = {"tri": [list(t) for t in tri], "computable_planes": [],
                         "letters": {}, "total_exponents": [], "n_total": 0,
                         "ladder_bearing": False,
                         "why": "untwisted leg is DELOCALISED: no NORM on any plane"}
            continue
        al = AZ.channel_alphabet(tri)["planes"]
        letters, planes = {}, []
        for a in (1, 2, 3):
            p = al.get("plane%d" % a, {})
            if p.get("status") != "computable":
                continue
            vals = sorted({F(e["NORM"]) for e in p["entries"] if e["NORM"] is not None})
            if not vals:
                continue
            planes.append(a)
            letters[a] = [str(v) for v in vals]
        tot = sorted({sum(c) for c in product(*[[F(x) for x in letters[a]]
                                               for a in planes])}) if planes else []
        out[name] = {"tri": [list(t) for t in tri], "computable_planes": planes,
                     "letters": letters, "total_exponents": [str(x) for x in tot],
                     "n_total": len(tot), "ladder_bearing": len(tot) >= 3,
                     "tau_invariant": tau_invariant_channel(tri)}
    _LADDER = out
    return out


# ================================================================ exchange (tau)
def tau(kl):
    """The unique [1,3] exchange automorphism tau(k,l) = (k, 4k+5l)."""
    k, l = kl[0] % 3, kl[1] % 6
    return (k % 3, (4 * k + 5 * l) % 6)


def tau_invariant_channel(tri):
    return sorted(tuple(tau(t)) for t in tri) == sorted(tuple(t) for t in tri)


def load_geom():
    return G

"""z3xz6_sweep -- STEP 2: the role census over all 4,612 `Z3xZ6_1_1` models.

Reads the STOCK-geometry census logs (`data/z3xz6_sweep/raw/census/census_seg*.txt.gz`)
and emits one incremental JSON per segment plus a catalogue-level roll-up.  Runs the
convention-blind machinery of `z36sw_lib`:

  * hypercharge t_Y from the model's OWN printed generator (the dual route having been
    verified in STEP 1 on 25 models, all eight defining conditions);
  * observable positions (a3, a2, sY) by the SMConfig CLOSURE -- the printed multiplicity
    table is REBUILT from the full left-chiral spectrum, never read by row position;
  * the role census over ALL FOUR (sY, sC) sign branches, the physical one selected by
    net Q = +3 (sY and sC are INDEPENDENT conventions; MSSM0 of this catalogue prints its
    quark doublet as (-3, 2, +1/6), i.e. in the branch a combined-flip-only scan cannot
    reach -- the live defect the Z6-I sweep caught mid-phase);
  * chirality of Q / u^c by maximum-matching semantics with the Z3 x Z6 point group over
    BOTH generators (sum k = 0 mod 3 AND sum l = 0 mod 6);
  * the D-2 uncovered-element tag: whether the model's candidate ladder legs sit on one
    of the SIX non-generator elements carrying a non-prime plane, for which the owned
    gamma lift is NOT defined (gate DEFECT-D-2).

============================ THE PRE-FILTER, STATED IN FULL ============================
The spec requires the pre-filter's rule to be written down so the gate can attack it.

  RULE P.  A model is classified `ladder-impossible` WITHOUT consulting its coupling run
  iff there is NO triple (Q, u, H) of its massless left-chiral fields with

     (i)   role(Q) = Q, role(u) = u^c, role(H) = H_u under the selected (sY, sC) branch;
     (ii)  Q and u are CHIRAL -- left unmatched by at least one MAXIMUM matching against
           their conjugate partners under the Z3 x Z6 bilinear point-group rule, with
           untwisted singlet vevs allowed to soak hidden/anomalous U(1) charges;
     (iii) the sector multiset {kl(Q), kl(u), kl(H)} is one of the ELEVEN LADDER-BEARING
           ALLOWED channels -- allowed by the published per-plane R-charge law
           0705.3075 (15), all three legs twisted (an untwisted leg is DELOCALISED and
           has no NORM), and admitting >= 3 distinct TOTAL suppression exponents.

  NECESSITY.  A rank-3 single-Higgs up-type block with three DISTINCT classical factors
  requires at least one `Q u^c H_u` trilinear that is (a) allowed -- otherwise orbifolder
  returns no such term at all -- and (b) on a channel admitting three distinct total
  exponents -- otherwise three distinct factors are impossible there.  Chirality of Q and
  u^c is required because a vector-like pair is removed at the string scale.

  DECLARED LOOSENESS.  Rule P asks for ONE triple, not nine, and says nothing about
  whether the triples assemble into a rank-3 matrix.  It is therefore strictly WEAKER
  than the ladder condition: it can only OVER-admit, never reject a model that could
  realise the ladder.  The total-exponent formulation in (iii) is deliberate -- a
  per-plane formulation would be a DIFFERENT and weaker object, and would have missed the
  carrier itself, whose per-plane alphabets are {0,1,3,4} twice and whose nine-factor set
  is their SUM.

  NOT LOAD-BEARING.  The coupling runs were executed on ALL 4,612 models regardless, so
  every Rule-P classification is re-checked against the model's actual superpotential
  (`sweep.py` gate X2).  A filter that gated nothing cannot have hidden a survivor.
========================================================================================

Usage: python census.py [--workers N] [--seg S]
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"                       # CLAUDE.md 4: set, not setdefault

import json                                                            # noqa: E402
import sys                                                             # noqa: E402
from collections import Counter, defaultdict                           # noqa: E402
from multiprocessing import Pool                                       # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z36sw_lib as Z                                                  # noqa: E402

OUT = os.path.join(Z.DATA, "census")
os.makedirs(OUT, exist_ok=True)


def ladder_bearing_multisets():
    """{frozen sector multiset: channel name} over the eleven ladder-bearing channels."""
    out = {}
    for name, v in Z.ladder_channels().items():
        if v["ladder_bearing"]:
            out[tuple(sorted(tuple(t) for t in v["tri"]))] = name
    return out


LB = ladder_bearing_multisets()
ALLOWED_MS = {tuple(sorted(tuple(t) for t in c)): Z.channel_name(c)
              for c in Z.allowed_channels()}


def analyse_model(label, block):
    spec = Z.parse_spectrum(block)
    gauge = Z.parse_gauge(block)
    r = {"model": label, "n_fields": len(spec),
         "gauge_TestConfig": gauge["gauge_test"], "gauge_SMConfig": gauge["gauge_sm"]}
    tY = gauge["tY"]
    if tY is None:
        r["status"] = "NO_HYPERCHARGE_PRINTED"
        return r
    pos = Z.observable_positions(spec, tY, gauge)
    if pos is None:
        r["status"] = "SMCONFIG_CLOSURE_FAILED"
        return r
    a3, a2, sY0, tY, bi = pos
    cen, diag = Z.role_census(spec, tY, a3, a2, sY0)
    r["observable_positions"] = {"a3": a3, "a2": a2, "sY_from_closure": sY0,
                                 "hypercharge_block_selected_by_closure": bi,
                                 "n_hypercharge_blocks": gauge["n_hyper_blocks"]}
    r["orientation"] = diag
    r["status"] = "ok"

    roles = Counter(v["role"] for v in cen.values() if v["role"])
    r["role_counts"] = dict(roles)

    chir = {}
    for pos_role, neg_role in (("Q", "Qbar"), ("uc", "ucbar"), ("dc", "dcbar"),
                               ("ec", "ecbar")):
        miss, nl, nr, nm = Z.chirality(cen, pos_role, neg_role)
        chir[pos_role] = {"n_left": nl, "n_right": nr, "matched": nm,
                          "chiral_labels": sorted(miss & {l for l, v in cen.items()
                                                          if v["role"] == pos_role}),
                          "n_chiral": len(miss & {l for l, v in cen.items()
                                                  if v["role"] == pos_role})}
    r["chirality"] = {k: {kk: vv for kk, vv in v.items() if kk != "chiral_labels"}
                      for k, v in chir.items()}

    Qs = chir["Q"]["chiral_labels"]
    Us = chir["uc"]["chiral_labels"]
    Hs = sorted(l for l, v in cen.items() if v["role"] == "Hu")
    r["n_chiral_Q"] = len(Qs)
    r["n_chiral_uc"] = len(Us)
    r["n_Hu"] = len(Hs)
    r["sectors"] = {
        "chiral_Q": dict(Counter(str(cen[l]["kl"]) for l in Qs)),
        "chiral_uc": dict(Counter(str(cen[l]["kl"]) for l in Us)),
        "Hu": dict(Counter(str(cen[l]["kl"]) for l in Hs)),
    }

    # ---- RULE P
    hits = defaultdict(int)
    for q in Qs:
        for u in Us:
            for h in Hs:
                ms = tuple(sorted([cen[q]["kl"], cen[u]["kl"], cen[h]["kl"]]))
                if ms in LB:
                    hits[LB[ms]] += 1
    r["ruleP_ladder_channel_triples"] = dict(hits)
    r["ruleP_pass"] = bool(hits)
    r["ruleP_carrier_triples"] = hits.get(Z.channel_name(Z.CARRIER), 0)

    # a strictly weaker diagnostic: triples on ANY allowed channel (ladder-bearing or not)
    anyc = defaultdict(int)
    for q in Qs:
        for u in Us:
            for h in Hs:
                ms = tuple(sorted([cen[q]["kl"], cen[u]["kl"], cen[h]["kl"]]))
                if ms in ALLOWED_MS:
                    anyc[ALLOWED_MS[ms]] += 1
    r["allowed_channel_triples"] = dict(anyc)

    # ---- D-2 uncovered-element tag (gate DEFECT-D-2)
    unc = sorted({str(cen[l]["kl"]) for l in Qs + Us + Hs
                  if tuple(cen[l]["kl"]) in Z.UNCOVERED_ELEMENTS})
    r["D2_uncovered_element_legs"] = unc
    r["D2_uncovered_tag"] = bool(unc)

    # ---- oscillator + gamma census (R2/R3 bookkeeping)
    r["n_fields_with_oscillators"] = sum(1 for v in spec.values()
                                         if v["osc_raw"] != "no")
    r["gamma_denominators"] = Z.check_gamma_denominators(spec)
    r["n_fields_both_twist_gammas_nonzero"] = sum(
        1 for v in spec.values()
        if Z.gamma_pair(v)[0] != 0 and Z.gamma_pair(v)[1] != 0)
    return r


def do_segment(fn):
    dest = os.path.join(OUT, fn.replace(".txt.gz", ".json"))
    if os.path.exists(dest):
        return dest
    txt = Z.read_raw(os.path.join(Z.RAW_CENSUS, fn))
    rows = []
    for lab, block in Z.split_models(txt):
        try:
            rows.append(analyse_model(lab, block))
        except Exception as e:                                    # loud, never silent
            rows.append({"model": lab, "status": "EXCEPTION", "error": repr(e)})
    tmp = dest + ".part"
    with open(tmp, "w") as fh:
        json.dump({"segment": fn, "n_models": len(rows), "models": rows}, fh,
                  indent=1, default=str)
    os.replace(tmp, dest)
    return dest


def main():
    workers = 6
    only = None
    for i, a in enumerate(sys.argv[1:]):
        if a == "--workers":
            workers = int(sys.argv[i + 2])
        if a == "--seg":
            only = sys.argv[i + 2]
    segs = sorted(f for f in os.listdir(Z.RAW_CENSUS)
                  if f.endswith(".txt.gz") and "PILOT" not in f)
    if only:
        segs = [f for f in segs if only in f]
    print("census: %d segments, %d workers" % (len(segs), workers))
    with Pool(workers) as p:
        for i, d in enumerate(p.imap_unordered(do_segment, segs)):
            print("  [%2d/%2d] %s" % (i + 1, len(segs), os.path.basename(d)))
    print("done")


if __name__ == "__main__":
    main()

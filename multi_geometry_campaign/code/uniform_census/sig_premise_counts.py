r"""The two record counts the signature section's Premises paragraph prints.

The paper writes, in the Premises paragraph of the signature section:

  * "the right-moving weight is the same for every one of the 30,922 records examined in
    the population relevant here, which is the models built over the first of the two
    shift classes together with the second catalogue at these two sectors, the catalogue
    side counted once per model as in Table 1";
  * "the left-moving oscillator flag reads absent in 2,239 of 2,239 and in 12,350 of
    12,350 records at the relevant norms".

Both are composite over two disjoint bodies of counting at T(0,1) and T(0,5) of
Z3xZ6_1_1:

  constructed-model half   `gate_su5_construction/g6_signature.json`, keys `cb_qsh`
                           (records) and `cb_osc_at_31_18` (records at |p_sh|^2 = 31/18).
                           Its population is `su5_construction/raw`, the 2,912 models
                           built over the FIRST shift class -- their common order-six
                           shift is (2/3,2/3,-1/3,-1/3,-1/3,0,0,0 | 1/2,1/2,0^6), the
                           class-A representative the paper prints.  The 1,432 models of
                           the second class (`gate_su5_construction/caseB_models.txt`)
                           are NOT in this half.  No duplicate labels, so this half is
                           unaffected by the catalogue de-duplication.
  catalogue half           recounted here off `z3xz6_sweep/raw/census`, each catalogue
                           model counted once -- the embedding-and-duplicate convention
                           of Table 1.  The census output contains ten of its models
                           twice, an overlap between two of its segments; the as-run
                           figures are emitted alongside so the difference is visible.

The de-duplication is the whole of the movement: `distinct` and `distinct_minus_e6`
below are equal, because a raw record count does not depend on the Standard-Model
closure test and the two models whose `E_n` factor is printed first were therefore
already in.

Run:  python sig_premise_counts.py
Out:  results_deposited/uniform_census/sig_premise_counts.json
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"

import glob                                                             # noqa: E402
import gzip                                                             # noqa: E402
import json                                                             # noqa: E402
import re                                                               # noqa: E402
import time                                                             # noqa: E402
from collections import Counter, defaultdict                            # noqa: E402
from fractions import Fraction as F                                     # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
CEN = os.path.join(REPO, "multi_geometry_campaign", "data", "z3xz6_sweep", "raw", "census")
G6 = os.path.join(REPO, "multi_geometry_campaign", "data", "gate_su5_construction",
                  "g6_signature.json")
OUT = os.path.join(REPO, "results_deposited", "uniform_census")

MARK = re.compile(r"/>\s*/(MSSM\d+)>")
HEAD = re.compile(r"(?m)^\s{2,8}(F_\d+)\s*$")
E6MODELS = {"MSSM825", "MSSM2669"}          # the two Z3xZ6 models printing E_6 first
SECT = ((0, 1), (0, 5))
P2 = F(31, 18)                              # the oscillator-free budget of these sectors


def rat(t):
    t = t.strip()
    return F(t) if t else F(0)


def catalogue_halves():
    """Records at T(0,1) / T(0,5), under three model populations."""
    pops = ("as_run", "distinct", "distinct_minus_e6")
    rec = {p: Counter() for p in pops}
    qsh = {p: defaultdict(Counter) for p in pops}
    n31 = {p: Counter() for p in pops}
    osc31 = {p: Counter() for p in pops}
    chir = {p: Counter() for p in pops}
    seen, dups = set(), Counter()
    for f in sorted(glob.glob(os.path.join(CEN, "*.txt.gz"))):
        txt = gzip.open(f, "rt", errors="replace").read()
        marks = [(m.start(), m.group(1)) for m in MARK.finditer(txt)]
        for i, (pos, lab) in enumerate(marks):
            end = marks[i + 1][0] if i + 1 < len(marks) else len(txt)
            pre = txt[pos:end].split("Analyzing")[0]
            which = ["as_run"]
            if lab in seen:
                dups[lab] += 1
            else:
                which.append("distinct")
                if lab not in E6MODELS:
                    which.append("distinct_minus_e6")
            seen.add(lab)
            heads = [(m.start(), m.group(1)) for m in HEAD.finditer(pre)]
            for j, (p0, _l) in enumerate(heads):
                e = heads[j + 1][0] if j + 1 < len(heads) else min(len(pre), p0 + 3000)
                b = pre[p0:e]
                ms = re.search(r"sector \(k,l\)\s*:\s*\(([^)]*)\)", b)
                if not ms:
                    continue
                kl = [int(x) for x in ms.group(1).split(",")]
                kl = (kl[0] % 3, kl[1] % 6)
                if kl not in SECT:
                    continue
                mq = re.search(r"right-moving q_sh\s*:\s*\(([^)]*)\)", b)
                mp = re.search(r"left-moving p_sh\s*:\s*(.*)", b)
                mr = re.search(r"representation\s*:\s*\(([^)]*)\)_([lr])", b)
                mo = re.search(r"oscillators\s*:\s*(.*)", b)
                comp = []
                if mp:
                    for grp in re.findall(r"\(([^)]*)\)", mp.group(1)):
                        comp.extend(rat(x) for x in grp.split(","))
                p2 = sum(x * x for x in comp) if len(comp) == 16 else None
                osc = mo.group(1).strip() if mo else "no"
                for w in which:
                    rec[w][kl] += 1
                    chir[w][mr.group(2) if mr else "?"] += 1
                    if mq:
                        qsh[w][kl][mq.group(1).replace(" ", "")] += 1
                    if p2 == P2:
                        n31[w][kl] += 1
                        osc31[w]["no" if osc in ("no", "") else "yes"] += 1
    out = {}
    for w in pops:
        out[w] = dict(
            records_T01=rec[w][(0, 1)], records_T05=rec[w][(0, 5)],
            records_total=rec[w][(0, 1)] + rec[w][(0, 5)],
            chirality=dict(chir[w]),
            distinct_qsh_T01=dict(qsh[w][(0, 1)]),
            distinct_qsh_T05=dict(qsh[w][(0, 5)]),
            n_distinct_qsh_T01=len(qsh[w][(0, 1)]),
            n_distinct_qsh_T05=len(qsh[w][(0, 5)]),
            n_at_31_18={str(k): v for k, v in n31[w].items()},
            n_at_31_18_total=sum(n31[w].values()),
            osc_flag_at_31_18=dict(osc31[w]))
    out["duplicate_labels"] = sorted(dups)
    out["duplicate_occurrences"] = sum(dups.values())
    return out


def main():
    t0 = time.time()
    cat = catalogue_halves()
    g6 = json.load(open(G6))
    cb_rec = g6["cb_qsh"]["(0, 1)"]["n_records"] + g6["cb_qsh"]["(0, 5)"]["n_records"]
    cb_osc = (g6["cb_osc_at_31_18"]["((0, 1), 'no')"]
              + g6["cb_osc_at_31_18"]["((0, 5), 'no')"])
    cb_qsh_distinct = [g6["cb_qsh"]["(0, 1)"]["n_distinct"],
                       g6["cb_qsh"]["(0, 5)"]["n_distinct"]]
    printed = dict(
        right_mover_records=cb_rec + cat["distinct"]["records_total"],
        osc_absent_constructed=cb_osc,
        osc_absent_catalogue=cat["distinct"]["n_at_31_18_total"],
        distinct_qsh_per_sector=sorted(set(
            cb_qsh_distinct
            + [cat["distinct"]["n_distinct_qsh_T01"],
               cat["distinct"]["n_distinct_qsh_T05"]])))
    as_run = dict(
        right_mover_records=cb_rec + cat["as_run"]["records_total"],
        osc_absent_catalogue=cat["as_run"]["n_at_31_18_total"])
    out = dict(
        geometry="Z3xZ6_1_1", sectors=[list(s) for s in SECT],
        seconds=round(time.time() - t0, 1),
        constructed_half=dict(models=g6["cb_models"], records=cb_rec,
                              records_at_31_18_osc_absent=cb_osc,
                              distinct_qsh_per_sector=cb_qsh_distinct),
        catalogue_half=cat,
        printed_by_the_paper=printed,
        as_run_for_comparison=as_run)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "sig_premise_counts.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("constructed half: %d records, %d at 31/18 with the oscillator flag absent"
          % (cb_rec, cb_osc))
    print("catalogue half  : %d records as-run, %d counted once"
          % (cat["as_run"]["records_total"], cat["distinct"]["records_total"]))
    print("                  %d at 31/18 as-run, %d counted once"
          % (cat["as_run"]["n_at_31_18_total"], cat["distinct"]["n_at_31_18_total"]))
    print("PRINTED: right-mover records %d ; oscillator flag absent %d of %d and %d of %d"
          % (printed["right_mover_records"], cb_osc, cb_osc,
             printed["osc_absent_catalogue"], printed["osc_absent_catalogue"]))
    print("distinct right-moving q_sh per sector: %s" % printed["distinct_qsh_per_sector"])


if __name__ == "__main__":
    main()

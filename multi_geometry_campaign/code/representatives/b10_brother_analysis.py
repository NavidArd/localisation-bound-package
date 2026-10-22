"""Referee #8, item 2.3: analyse the brother census (b5 build, b6 orbifolder run).

Per deck: triality-aware, hidden-dimension-weighted net indices (b4_triality_nets.nets);
the doubled 10-like test (a left-chiral 10-like record at T(0,1)/T(0,5) with |p|^2 = 31/18
and hidden dimension 2 -- the form section 8 of paperA_v9d shows the doubling takes);
the first-E8 gauge group line.  Criteria exactly as the certified gate
(gate_su5_construction g5_battery.py / g8_batteryB.py):
  three generations at once  <=>  net(Q) = net(u^c) = net(e^c) = 3e for one e = +-1;
  two of three               <=>  two of the three equal 3e for one e.
Control: the (a,b) = (0,0) decks are the census representatives themselves; their nets are
compared with the certified per-model nets of g5_battery.json / g8_batteryB.json.

Label check (package v2.2 repair, 2026-10-22).  Every deck in brothers/decks_*.txt is read
back and its slot is RE-COMPUTED from its own vectors: (a, b) = (3 lam.V2 mod 3, 3 mu.V2 mod 3)
with lam = V1 - V1(parent), mu = W3 - W3(parent); V2 must equal the parent's, W4 = W3, and
the label must be [XY]<parent>a<a>b<b>.  The result is recorded under `label_check`.  As
first deposited this script took (a, b) from the build record and never looked at the
vectors, which is how 80 decks labelled (0,1)/(0,2) that realise (0,0) went unseen
(b5b_missing.py; brothers/superseded_v2.2/NOTICE.md).  The decks the repaired builder writes
for slots that need a mu search are labelled Y<parent>a<a>b<b>; both prefixes are read.

Output: review/data/paperA_referee8/b10_brothers.json
"""
import glob
import gzip
import json
import os
_PKG = os.environ.get("PAPERA_PKG") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
import re
import sys
from collections import Counter
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b4_triality_nets as T  # noqa: E402
import b5_build_brother_census as B5  # noqa: E402

D = os.path.join(_PKG, "multi_geometry_campaign", "data", "representatives", "brothers")
GATE = os.path.join(_PKG, "multi_geometry_campaign", "data", "gate_su5_construction")
TEN = ("Q", "uc", "ec")


def doubled(spec):
    n = 0
    for r in spec.values():
        if r["chir"] != "l" or not r["p_sh"] or r["kl"] not in ((0, 1), (0, 5)):
            continue
        p = r["p_sh"]
        if T.S.dot(p, p) != F(31, 18):
            continue
        rl, _sg = T.role(p)
        if rl not in TEN:
            continue
        tot = 1
        for d in r["rep"]:
            tot *= abs(d)
        c3 = T.colour(p)
        dsm = (3 if c3 in (3, -3) else 1) * (2 if (p[3] - p[4]) != 0 else 1)
        if tot // dsm == 2:
            n += 1
    return n


def label_check(meta):
    """re-compute every deck's (a, b) from its vectors and compare with its label/record."""
    par = dict(B5.read_decks(os.path.join(B5.CAMP, "su5_construction", "s5_models.txt")))
    par.update(B5.read_decks(os.path.join(B5.CAMP, "gate_su5_construction",
                                          "caseB_models.txt")))

    def ph(v, V2):
        t = 3 * B5.S.dot(v, V2)
        return int(t) % 3 if t.denominator == 1 else None

    seen, fails = set(), []
    for pth in sorted(glob.glob(os.path.join(D, "decks_*.txt"))):
        for lab, rows in B5.read_decks(pth):
            seen.add(lab)
            m = meta.get(lab)
            if m is None:
                fails.append(dict(label=lab, why="no build record"))
                continue
            P = par[m["parent"]]
            lam = [x - y for x, y in zip(rows[0], P[0])]
            mu = [x - y for x, y in zip(rows[4], P[4])]
            ab = (ph(lam, rows[1]), ph(mu, rows[1]))
            why = []
            if ab != (m["a"], m["b"]):
                why.append("vectors give (a,b) = %s, label says (%d,%d)" % (ab, m["a"], m["b"]))
            if rows[1] != P[1]:
                why.append("V2 differs from the parent's")
            if rows[4] != rows[5]:
                why.append("W4 != W3")
            if not re.fullmatch(r"[XY]%sa%db%d" % (m["parent"], m["a"], m["b"]), lab):
                why.append("label does not name its slot")
            if why:
                fails.append(dict(label=lab, why="; ".join(why)))
    for lab in sorted(set(meta) - seen):
        fails.append(dict(label=lab, why="build record without a deck"))
    return dict(decks_read=len(seen), passed=len(seen) - sum(1 for f in fails
                                                            if f["label"] in seen),
                failed=len(fails), failures=fails[:20])


def main():
    meta = json.load(open(os.path.join(D, "build_all.json")))["meta"]
    lchk = label_check(meta)
    rows = {}
    for pth in sorted(glob.glob(os.path.join(D, "out_*.txt.gz"))):
        txt = gzip.open(pth, "rt", errors="replace").read()
        for lab, blk in T.split(txt, pat=r"/> /([XY]\w+?a\db\d)>"):
            spec = T.SW.parse_spectrum(blk.split("/gauge group>")[0])
            gg = re.search(r'Gauge group in vev-configuration "TestConfig1": (.*)', blk)
            rows[lab] = dict(nets=T.nets(spec.values()) if spec else None,
                             n_left=sum(1 for r in spec.values() if r["chir"] == "l"),
                             dbl=doubled(spec) if spec else 0,
                             gauge=gg.group(1).strip() if gg else None)
    # certified control
    cert = {}
    for fn in ("g5_battery.json", "g8_batteryB.json"):
        d = json.load(open(os.path.join(GATE, fn)))
        for r in d.get("rows", d.get("per_model", [])):
            if "net_aware" in r:
                cert[r["model"]] = r["net_aware"]
    ctrl = Counter()
    for lab, r in rows.items():
        m = meta[lab]
        if (m["a"], m["b"]) == (0, 0) and m["parent"] in cert and r["nets"]:
            c = cert[m["parent"]]
            ctrl["agree" if all(int(c[k]) == r["nets"][k] for k in TEN) else "DISAGREE"] += 1
    missing = sorted(set(meta) - set(rows))
    gen3, two3, dbl, byab, spec_differs = [], Counter(), Counter(), Counter(), Counter()
    base = {meta[l]["parent"]: r for l, r in rows.items()
            if (meta[l]["a"], meta[l]["b"]) == (0, 0)}
    for lab, r in rows.items():
        if not r["nets"]:
            continue
        m = meta[lab]
        ab = "%d%d" % (m["a"], m["b"])
        byab[ab] += 1
        n = r["nets"]
        hit = any(all(n[k] == 3 * e for k in TEN) for e in (1, -1))
        h2 = any(sum(n[k] == 3 * e for k in TEN) >= 2 for e in (1, -1))
        if hit:
            gen3.append(dict(label=lab, parent=m["parent"], a=m["a"], b=m["b"], nets=n,
                             doubled=r["dbl"]))
        if h2:
            two3[ab] += 1
        if r["dbl"]:
            dbl[ab] += 1
        b0 = base.get(m["parent"])
        if b0 and ab != "00" and (b0["n_left"] != r["n_left"] or b0["nets"] != n):
            spec_differs[ab] += 1
    out = dict(decks_built=len(meta), decks_with_spectrum=len(rows), label_check=lchk,
               decks_missing_output=missing[:20], n_missing_output=len(missing),
               control_vs_certified=dict(ctrl), by_ab=dict(byab),
               three_generations_at_once=gen3, n_three_generations=len(gen3),
               n_three_gen_parents=len({g["parent"] for g in gen3}),
               n_three_gen_with_doubling=sum(1 for g in gen3 if g["doubled"]),
               two_of_three_by_ab=dict(two3), doubled_by_ab=dict(dbl),
               differs_from_census_rep_by_ab=dict(spec_differs))
    import csv
    with gzip.open(os.path.join(os.path.dirname(D), "b10_per_deck_nets.csv.gz"), "wt",
                   newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["label", "parent", "class", "a", "b", "n_left_records", "Q", "uc",
                    "ec", "dc", "L", "doubled_10like_records_T0l"])
        for lab in sorted(rows):
            r, m = rows[lab], meta[lab]
            n = r["nets"] or {}
            w.writerow([lab, m["parent"], m["cls"], m["a"], m["b"], r["n_left"]]
                       + [n.get(k, "") for k in ("Q", "uc", "ec", "dc", "L")]
                       + [r["dbl"]])
    json.dump(out, open(os.path.join(os.path.dirname(D), "b10_brothers.json"), "w"),
              indent=1, default=str)
    print(json.dumps({k: v for k, v in out.items() if k != "three_generations_at_once"},
                     indent=1))
    for g in gen3[:20]:
        print(g)


if __name__ == "__main__":
    main()

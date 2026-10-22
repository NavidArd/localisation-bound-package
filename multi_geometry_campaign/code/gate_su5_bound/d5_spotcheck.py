"""GATE gate_su5_bound D5 -- FRESH 30-CELL SPOT-CHECK against the raw orbifolder text.

Everything here is parsed by THIS FILE's own regexes off the gzipped census text; nothing
about the Standard-Model identification is imported.  For each selected cell we extract

    sector (k,l), fixed point n_a, and for every left-chiral field at that cell:
    left-moving p_sh, representation, U(1) charges, oscillators,

plus the model's printed `Hypercharge generator` and its printed SM-configuration table.
We then assign SM roles OURSELVES:

    c3 = signed SU(3) entry of the representation tuple  (triality: +3 vs -3)
    c2 = |SU(2) entry|
    Y  = sY * (t_Y . p_sh)                       [t_Y read from the printed generator]
    k_Y = 2 |t_Y|^2                              [GUT normalisation, measured not assumed]

with (a3, a2) located by our own reconstruction of orbifolder's printed table and the
sign branch (sC, sY) fixed by OUR OWN criterion -- net chirality exactly (3,3,3,3,3) for
(Q, u^c, d^c, L, e^c) -- which is NOT the criterion z24sw_lib.role_census uses (it takes
argmax of (net Q, net u^c)).  The two answers are then compared field by field.

Cells are chosen to span both geometries and several sectors, with quotas: at least 5 cells
carrying all three 10-roles and at least 5 cells carrying a DOUBLED 5bar role.

Outputs data/gate_su5_bound/d5_spotcheck.json
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"

import glob                                                            # noqa: E402
import gzip                                                            # noqa: E402
import json                                                            # noqa: E402
import re                                                              # noqa: E402
import sys                                                             # noqa: E402
from collections import Counter, defaultdict                           # noqa: E402
from fractions import Fraction as F                                    # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CAMP = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(CAMP, "data", "gate_su5_bound")

GEOM = {
    "Z2xZ4_1_1": dict(lib="z2xz4_sweep", mod="z24sw_lib", M=2, N=4,
                      census=os.path.join(CAMP, "data", "z2xz4_sweep", "raw", "census")),
    "Z3xZ6_1_1": dict(lib="z3xz6_sweep", mod="z36sw_lib", M=3, N=6,
                      census=os.path.join(CAMP, "data", "z3xz6_sweep", "raw", "census")),
}

# ---------------------------------------------------------------- own parsing
_BLK = re.compile(r"/> /(MSSM\d+)> /")
_FLD = re.compile(r"(?m)^\s{2,8}(F_\d+)\s*$")
_ROW = re.compile(r"(?m)^\s*(\d+)\s+\(\s*(-?\d+),\s*(-?\d+)\)_([lr])\s+"
                  r"U\(1\)\s*:\s*\(\s*([^)]*)\)")
_GG = re.compile(r'Gauge group in vev-configuration "TestConfig1":\s*(.*)')


def rat(x):
    x = x.strip()
    return F(x) if x else F(0)


def vec(s):
    out = []
    for g in re.findall(r"\(([^)]*)\)", s):
        out += [rat(x) for x in g.split(",")]
    return out


def own_split(text):
    idx = [(m.start(), m.group(1)) for m in _BLK.finditer(text)]
    return [(lab, text[p:(idx[i + 1][0] if i + 1 < len(idx) else len(text))])
            for i, (p, lab) in enumerate(idx)]


def own_fields(pre, M, N):
    """every printed field record, verbatim, by our own regex."""
    heads = [(m.start(), m.group(1)) for m in _FLD.finditer(pre)]
    out = []
    for i, (p, lab) in enumerate(heads):
        b = pre[p:(heads[i + 1][0] if i + 1 < len(heads) else min(len(pre), p + 3000))]

        def g(pat):
            m = re.search(pat, b)
            return m.group(1) if m else None

        sec = g(r"sector \(k,l\)\s*:\s*\(([^)]*)\)")
        if sec is None:
            continue
        kl = [int(x) for x in sec.split(",")]
        kl = (kl[0] % M, kl[1] % N)
        fps = g(r"fixed point n_a\s*:\s*\(([^)]*)\)")
        fp = tuple(int(x) for x in fps.split(",")) if fps else None
        ps = g(r"left-moving p_sh\s*:\s*(.*)")
        psh = tuple(vec(ps)) if ps else None
        m = re.search(r"representation\s*:\s*\(([^)]*)\)_([lr])\s+U\(1\):\s*\(([^)]*)\)", b)
        rep = chir = u1 = None
        if m:
            rep = tuple(int(re.match(r"^\s*(-?\d+)", x).group(1))
                        for x in m.group(1).split(","))
            chir = m.group(2)
            u1 = tuple(rat(x) for x in m.group(3).split(","))
        osc = g(r"oscillators\s*:\s*(.*)")
        out.append(dict(label=lab, kl=kl, fp=fp, p_sh=psh, rep=rep, chir=chir,
                        u1=u1, osc=(osc or "").strip()))
    return out


def own_configs(blk):
    """[(t_Y, printed_table_rows)] in print order; rows = (c3, c2, Y, mult, chir)."""
    cuts = [m.start() for m in re.finditer(r"Hypercharge generator:", blk)]
    out = []
    for i, p in enumerate(cuts):
        seg = blk[p:(cuts[i + 1] if i + 1 < len(cuts) else len(blk))]
        m = re.search(r"Hypercharge generator:\s*\n(.*)", seg)
        tY = vec(m.group(1)) if m else None
        # the FIRST contiguous SM table in this segment.  orbifolder re-prints the same
        # table later in the block on some geometries; collecting every _ROW match in the
        # segment would double every multiplicity.
        h = re.search(r'Gauge group in vev-configuration "[^"]*":[^\n]*U\(1\)_Y[^\n]*\n',
                      seg)
        rows = []
        if h:
            for ln in seg[h.end():].split("\n"):
                mm = _ROW.match(ln)
                if mm:
                    rows.append((int(mm.group(2)), int(mm.group(3)), rat(mm.group(5)),
                                 int(mm.group(1)), mm.group(4)))
                elif ln.strip() == "":
                    continue
                else:
                    break
        out.append((tY, rows))
    return out


def own_factors(blk):
    m = _GG.search(blk)
    if not m:
        return []
    s = m.group(1).split("and U(1)")[0]
    return [(t, int(n)) for t, n in re.findall(r"(SU|SO|E)\((\d+)\)", s)]


SM = {(3, 2, F(1, 6)): "Q", (-3, 2, F(-1, 6)): "Qbar",
      (-3, 1, F(-2, 3)): "uc", (3, 1, F(2, 3)): "ucbar",
      (-3, 1, F(1, 3)): "dc", (3, 1, F(-1, 3)): "dcbar",
      (1, 2, F(-1, 2)): "L", (1, 2, F(1, 2)): "Hu",
      (1, 1, F(1)): "ec", (1, 1, F(-1)): "ecbar"}


def own_roles(fields, tY, rows, factors):
    """OUR assignment.  Returns (roles, meta) or (None, why)."""
    left = [r for r in fields if r["rep"] and r["chir"] == "l" and r["p_sh"]]
    printed = defaultdict(int)
    for c3, c2, Y, mult, ch in rows:
        printed[(c3, c2, Y)] += mult
    printed = dict(printed)
    cands = []
    for a3, f3 in enumerate(factors):
        if f3 != ("SU", 3):
            continue
        for a2, f2 in enumerate(factors):
            if a2 == a3 or f2 != ("SU", 2):
                continue
            for sYt in (1, -1):
                tbl = defaultdict(int)
                ok = True
                for r in left:
                    if max(a3, a2) >= len(r["rep"]):
                        ok = False
                        break
                    hid = 1
                    for k, d in enumerate(r["rep"]):
                        if k not in (a3, a2):
                            hid *= abs(d)
                    c3v = r["rep"][a3] if abs(r["rep"][a3]) != 1 else 1
                    tbl[(c3v, abs(r["rep"][a2]),
                         sYt * sum(a * b for a, b in zip(tY, r["p_sh"])))] += hid
                if ok and dict(tbl) == printed:
                    cands.append((a3, a2, sYt))
    if not cands:
        return None, "no (a3,a2,sY) reproduces the printed table"
    a3, a2, sYt = cands[0]
    # ---- our own branch criterion: net chirality exactly (3,3,3,3,3)
    winners = []
    for sC in (1, -1):
        for sY in (1, -1):
            roles, net = {}, Counter()
            for r in left:
                hid = 1
                for k, d in enumerate(r["rep"]):
                    if k not in (a3, a2):
                        hid *= abs(d)
                c3v = r["rep"][a3]
                c3v = (sC * c3v) if abs(c3v) != 1 else 1
                c2v = abs(r["rep"][a2])
                Y = sY * sYt * sum(a * b for a, b in zip(tY, r["p_sh"]))
                ro = SM.get((c3v, c2v, Y))
                roles[r["label"]] = dict(role=ro, Y=str(Y), mult=hid, c3=c3v, c2=c2v)
                if ro:
                    net[ro] += hid
            nets = {p: net[p] - net[m] for p, m in
                    (("Q", "Qbar"), ("uc", "ucbar"), ("dc", "dcbar"),
                     ("ec", "ecbar"), ("L", "Hu"))}
            if all(v == 3 for v in nets.values()):
                winners.append((sC, sY, roles, nets))
    if len(winners) != 1:
        return None, "net-(3,3,3,3,3) branch not unique: %d" % len(winners)
    sC, sY, roles, nets = winners[0]
    kY = 2 * sum(x * x for x in tY)
    return roles, dict(a3=a3, a2=a2, sY_table=sYt, sC=sC, sY=sY, nets=nets,
                       kY=str(kY), n_candidates=len(cands))


# ---------------------------------------------------------------- driver
def collect(geom, want=15, nfiles=2):
    G = GEOM[geom]
    sys.path.insert(0, os.path.join(CAMP, "code", G["lib"]))
    Z = __import__(G["mod"])
    R10 = {"Q", "uc", "ec", "Qbar", "ucbar", "ecbar"}
    R5 = {"dc", "L", "Hu", "dcbar"}
    picked, quota = [], Counter()
    seen_sec = Counter()
    files = sorted(glob.glob(os.path.join(G["census"], "*.txt.gz")))[:nfiles]
    for path in files:
        txt = gzip.open(path, "rt", errors="replace").read()
        for lab, blk in own_split(txt):
            if len(picked) >= want:
                break
            pre = blk.split("Analyzing")[0]
            fields = own_fields(pre, G["M"], G["N"])
            cfgs = own_configs(blk)
            factors = own_factors(blk)
            if not cfgs or cfgs[0][0] is None:
                continue
            tY, rows = cfgs[0]
            mine, meta = own_roles(fields, tY, rows, factors)
            if mine is None:
                continue
            # ---- the certified route, for comparison only
            try:
                spec = Z.parse_spectrum(pre, M=G["M"], N=G["N"]) if geom == "Z2xZ4_1_1" \
                    else Z.parse_spectrum(pre)
                if geom == "Z2xZ4_1_1":
                    pr = defaultdict(int)
                    for row in rows:
                        pr[(row[0], row[1], row[2])] += row[3]
                    pos = Z.observable_positions(spec, tY, Z.parse_gauge_factors(blk),
                                                 dict(pr))
                    if pos is None:
                        continue
                    theirs, _ = Z.role_census(spec, tY, *pos)
                else:
                    gauge = Z.parse_gauge(blk)
                    pos = Z.observable_positions(spec, gauge["tY"], gauge)
                    if pos is None:
                        continue
                    theirs, _ = Z.role_census(spec, gauge["tY"], pos[0], pos[1], pos[2])
            except Exception as e:                                     # noqa: BLE001
                continue
            cells = defaultdict(list)
            for r in fields:
                if r["chir"] != "l" or r["p_sh"] is None:
                    continue
                if tuple(r["kl"]) == (0, 0):
                    continue
                cells[(tuple(r["kl"]), r["fp"])].append(r)
            for key, recs in sorted(cells.items()):
                if len(picked) >= want:
                    break
                mult = Counter()
                for r in recs:
                    ro = mine[r["label"]]["role"]
                    if ro:
                        mult[ro] += 1
                ten = {r for r in mult if r in R10}
                five_dbl = any(mult[r] >= 2 for r in mult if r in R5)
                three_ten = len({x.replace("bar", "") for x in ten}) >= 3
                tag = None
                if three_ten and quota["ten3"] < 5:
                    tag = "ten3"
                elif five_dbl and quota["dbl5"] < 5:
                    tag = "dbl5"
                elif seen_sec[key[0]] < 1 and mult:
                    tag = "sector"
                if tag is None:
                    continue
                quota[tag] += 1
                seen_sec[key[0]] += 1
                cmp_rows, agree, disagree = [], 0, 0
                for r in recs:
                    a = mine[r["label"]]["role"]
                    b = (theirs.get(r["label"]) or {}).get("role")
                    if a == b:
                        agree += 1
                    else:
                        disagree += 1
                    cmp_rows.append(dict(
                        field=r["label"], rep=list(r["rep"]),
                        u1=[str(x) for x in r["u1"]],
                        p_sh=[str(x) for x in r["p_sh"]],
                        osc=r["osc"], Y_mine=mine[r["label"]]["Y"],
                        c3=mine[r["label"]]["c3"], c2=mine[r["label"]]["c2"],
                        hidden_mult=mine[r["label"]]["mult"],
                        role_mine=a, role_role_census=b, same=(a == b)))
                picked.append(dict(
                    geometry=geom, model=lab, file=os.path.basename(path), tag=tag,
                    sector=list(key[0]), fixed_point=list(key[1]) if key[1] else None,
                    hypercharge_generator=[str(x) for x in tY],
                    k_Y=meta["kY"], branch=dict(a3=meta["a3"], a2=meta["a2"],
                                                sY_table=meta["sY_table"],
                                                sC=meta["sC"], sY=meta["sY"]),
                    nets=meta["nets"],
                    printed_SM_table=[[r[0], r[1], str(r[2]), r[3], r[4]] for r in rows],
                    role_multiplicity_records={k: v for k, v in mult.items()},
                    role_multiplicity_states={
                        k: sum(mine[r["label"]]["mult"] for r in recs
                               if mine[r["label"]]["role"] == k) for k in mult},
                    n_fields=len(recs), agree=agree, disagree=disagree,
                    fields=cmp_rows))
        if len(picked) >= want:
            break
    return picked


def main():
    cells = collect("Z2xZ4_1_1", 15, 2) + collect("Z3xZ6_1_1", 15, 2)
    tot_a = sum(c["agree"] for c in cells)
    tot_d = sum(c["disagree"] for c in cells)
    summ = dict(n_cells=len(cells), fields_compared=tot_a + tot_d,
                agreements=tot_a, disagreements=tot_d,
                by_tag=dict(Counter(c["tag"] for c in cells)),
                by_geom=dict(Counter(c["geometry"] for c in cells)),
                sectors=dict(Counter("%s|%s" % (c["geometry"], c["sector"])
                                     for c in cells)),
                cells_with_ten3=sum(1 for c in cells if c["tag"] == "ten3"),
                cells_with_doubled_5bar=sum(1 for c in cells if c["tag"] == "dbl5"),
                kY_values=dict(Counter(c["k_Y"] for c in cells)),
                state_vs_record_mismatch=[
                    dict(model=c["model"], geometry=c["geometry"], sector=c["sector"],
                         fp=c["fixed_point"], rec=c["role_multiplicity_records"],
                         st=c["role_multiplicity_states"])
                    for c in cells
                    if c["role_multiplicity_records"] != c["role_multiplicity_states"]])
    print(json.dumps(summ, indent=1, default=str)[:6000])
    for c in cells:
        if c["disagree"]:
            print("DISAGREE", json.dumps(c, indent=1, default=str)[:4000])
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "d5_spotcheck.json"), "w") as fh:
        json.dump(dict(summary=summ, cells=cells), fh, indent=1, default=str)


if __name__ == "__main__":
    main()

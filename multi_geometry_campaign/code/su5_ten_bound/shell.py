"""su5_ten_bound STEP 5 -- EXACT MASSLESS-SHELL ENUMERATION at the loophole cells.

The bound is a NECESSARY condition.  Where it does not bite -- the Z3xZ6_1_1 sectors
T(0,1) and T(0,5), where R = 47/90 >= 1/2 -- the question of whether a fixed point can
carry TWO 10-role states is decided by an EXACT count, not by an inequality:

    enumerate every p in V_g + Lambda with |p|^2 = 2B,
    project each onto the SM Cartan span, and count the states carrying each SM role.

This is the complete left-moving massless level of the cell BEFORE the orbifold
projection.  If the shell already contains only one state per 10-role, the loophole is
closed at that cell by lattice geometry alone.  If it contains two or more, the census's
multiplicity 1 is due to the ORBIFOLD PROJECTION, and the loophole is open with a
precisely stated construction target.

Enumeration is exact: per E8 factor the coset representative fixes each coordinate's
fractional part up to the D8/(D8+s) branch, and |x|^2 <= 2B bounds each coordinate, so the
search is a small finite box; the even-sum condition selects Lambda_{E8}.

Output: data/su5_ten_bound/shell_<geom>.json
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"

import glob                                                             # noqa: E402
import gzip                                                             # noqa: E402
import itertools                                                        # noqa: E402
import json                                                             # noqa: E402
import sys                                                              # noqa: E402
from collections import Counter, defaultdict                            # noqa: E402
from fractions import Fraction as F                                     # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CAMP = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
import s5lib as S                                                       # noqa: E402
from localgut import lgroup                                             # noqa: E402
from verify import GEOM                                                 # noqa: E402

OUT = os.path.join(CAMP, "data", "su5_ten_bound")


def e8_coset_shell(V8, cmax):
    """All x in V8 + Lambda_{E8} with |x|^2 <= cmax.  V8: 8 Fractions.

    Lambda_{E8} = D8 u (D8 + s), s = (1/2)^8, D8 = integer vectors of EVEN sum.  In both
    branches x - V8 (resp. x - V8 - s) must have even coordinate sum, i.e.

        sum(x) - sum(V8)  in  2Z

    -- the parity condition is RELATIVE to the coset representative.  Testing sum(x) in 2Z
    is the V8 = 0 special case and silently returns an EMPTY shell for a generic local
    shift; it was caught by the invariant that a printed state must lie in its own shell.
    """
    sV = sum(V8)
    out = []
    for half in (F(0), F(1, 2)):
        ranges = []
        ok = True
        for j in range(8):
            base = (V8[j] + half) % 1
            lo, hi = -cmax, cmax
            vals = []
            k = -6
            while k <= 6:
                x = base + k
                if x * x <= cmax:
                    vals.append(x)
                k += 1
            if not vals:
                ok = False
                break
            ranges.append(vals)
        if not ok:
            continue
        for cand in itertools.product(*ranges):
            n = sum(x * x for x in cand)
            if n > cmax:
                continue
            if (sum(cand) - sV) % 2 != 0:
                continue
            out.append((tuple(cand), n))
    # de-duplicate (the two `half` branches can coincide when V8 is half-integral)
    return sorted(set(out))


def shell(V16, c):
    """All p in V16 + Lambda with |p|^2 == c, as a list of 16-tuples."""
    A = e8_coset_shell(V16[:8], c)
    B_ = e8_coset_shell(V16[8:], c)
    byn = defaultdict(list)
    for v, n in B_:
        byn[n].append(v)
    out = []
    for v, n in A:
        need = c - n
        for w in byn.get(need, ()):
            out.append(v + w)
    return out


ROLE_W2 = {"10": F(6, 5), "5bar": F(4, 5)}


def analyse_cell(V16, c, smbasis, tY, f3, f2):
    """Decompose the shell by SM quantum numbers.

    The SM role of a shell state is read from its SU(3)/SU(2) Dynkin labels and its
    hypercharge Y = tY . p -- the same three defining conditions the census uses (5.21).
    """
    sts = shell(V16, c)
    cnt = Counter()
    detail = defaultdict(list)
    for p in sts:
        Y = S.dot(tY, p)
        l3 = tuple(S.dot(p, a) for a in f3)
        l2 = tuple(S.dot(p, a) for a in f2)
        if any(x.denominator != 1 for x in l3 + l2):
            cnt["NONINT"] += 1
            continue
        w2 = S.proj_norm2(p, smbasis)
        key = (tuple(int(x) for x in l3), tuple(int(x) for x in l2), str(Y))
        cnt[str(key)] += 1
        detail[str(key)].append(str(w2))
    return len(sts), cnt, detail


def run(geom="Z3xZ6_1_1", sectors=((0, 1), (0, 5)), maxcells=25, only=None):
    G = GEOM[geom]
    sys.path.insert(0, os.path.join(CAMP, "code", G["lib"]))
    Z = __import__(G["mod"])
    rep = dict(geometry=geom, sectors=[list(s) for s in sectors], cells=[],
               summary=Counter())
    seen = 0
    files = sorted(glob.glob(os.path.join(G["census"], "*.txt.gz")))
    if only:
        files = [f for f in files if os.path.basename(f) in only_files]
    for pth in files:
        if seen >= maxcells:
            break
        txt = gzip.open(pth, "rt", errors="replace").read()
        for lab, blk in Z.split_models(txt):
            if seen >= maxcells:
                break
            if only and lab not in only:
                continue
            pre = blk.split("Analyzing")[0]
            try:
                spec = Z.parse_spectrum(pre)
                gauge = Z.parse_gauge(blk)
                if gauge["tY"] is None or not spec:
                    continue
                pos = Z.observable_positions(spec, gauge["tY"], gauge)
                if pos is None:
                    continue
                a3, a2, sY0, tY, bi = pos
                roles, _ = Z.role_census(spec, tY, a3, a2, sY0)
            except Exception:                                          # noqa: BLE001
                continue
            cells = defaultdict(list)
            for l, r in spec.items():
                if r["p_sh"] is None or r["chir"] != "l":
                    continue
                cells[(tuple(r["kl"]), tuple(r["fp"]) if r["fp"] else None)].append(r)
            # 4D group -> SM basis
            import numpy as np
            allr = S.e8xe8_roots()
            R2 = np.array([[int(2 * x) for x in a] for a in allr], dtype=np.int64)
            inter = None
            for key, recs in cells.items():
                V = recs[0]["p_sh"]
                D = 1
                for x in V:
                    D = D * x.denominator // np.gcd(D, x.denominator)
                m = ((R2 @ np.array([int(x * D) for x in V], dtype=np.int64))
                     % (2 * D)) == 0
                inter = m if inter is None else (inter & m)
            R4 = [allr[i] for i in np.nonzero(inter)[0]]
            Sp4 = S.simple_roots(R4)
            comps4, adj4 = S._components(Sp4)
            qrec = [r for r in spec.values()
                    if roles.get(r["label"], {}).get("role") in ("Q", "Qbar")]
            drec = [r for r in spec.values()
                    if roles.get(r["label"], {}).get("role") in ("dc", "dcbar")]
            lrec = [r for r in spec.values()
                    if roles.get(r["label"], {}).get("role") in ("L", "Hu")]
            f3 = f2 = None
            for c_ in comps4:
                kind, rk = S._classify(c_, adj4, Sp4)
                so = [Sp4[i] for i in c_]
                if kind == "A" and rk == 2 and (qrec or drec):
                    if all(any(S.dot(x["p_sh"], a) != 0 for a in so) for x in qrec + drec):
                        f3 = so if f3 is None else "AMBIG"
                if kind == "A" and rk == 1 and (qrec or lrec):
                    if all(any(S.dot(x["p_sh"], a) != 0 for a in so) for x in qrec + lrec):
                        f2 = so if f2 is None else "AMBIG"
            if not (isinstance(f3, list) and isinstance(f2, list)):
                continue
            smbasis = list(f3) + list(f2) + [tuple(tY)]

            for key, recs in cells.items():
                if tuple(key[0]) not in sectors:
                    continue
                tens = [r for r in recs
                        if roles.get(r["label"], {}).get("role") in ("Q", "uc", "ec")]
                if not tens:
                    continue
                V = recs[0]["p_sh"]
                c = S.dot(tens[0]["p_sh"], tens[0]["p_sh"])
                if any(r.get("osc_raw", "no") not in ("no", "") for r in tens):
                    continue
                sts = set(shell(list(V), c))
                inshell = all(tuple(r["p_sh"]) in sts for r in tens)
                n_shell, cnt, det = analyse_cell(list(V), c, smbasis, list(tY), f3, f2)
                g = lgroup(V)
                # which printed roles, and how many shell states share each key
                printed = Counter(roles.get(r["label"], {}).get("role") for r in tens)
                tenkeys = {}
                for r in tens:
                    ro = roles.get(r["label"], {}).get("role")
                    l3 = tuple(int(S.dot(r["p_sh"], a)) for a in f3)
                    l2 = tuple(int(S.dot(r["p_sh"], a)) for a in f2)
                    k2 = str((l3, l2, str(S.dot(tY, r["p_sh"]))))
                    tenkeys[ro] = (k2, cnt.get(k2, 0))
                ent = dict(model=lab, kl=list(key[0]), fp=list(key[1]) if key[1] else None,
                           group=g["name"], p2=str(c), n_shell=n_shell,
                           printed_states_found_in_shell=bool(inshell),
                           printed_10=dict(printed),
                           shell_multiplicity_of_printed_10={
                               k: v[1] for k, v in tenkeys.items()},
                           has_commutant_SU2=any(
                               f["name"] == "SU(2)" for f in g["factors"]))
                rep["cells"].append(ent)
                rep["summary"]["invariant_printed_in_shell_%s" % inshell] += 1
                for k, v in tenkeys.items():
                    rep["summary"]["shell_mult_%d" % v[1]] += 1
                seen += 1
                if seen >= maxcells:
                    break
    rep["summary"] = dict(rep["summary"])
    with open(os.path.join(OUT, "shell_%s.json" % geom), "w") as fh:
        json.dump(rep, fh, indent=1, default=str)
    print(json.dumps(rep, indent=1, default=str)[:6000])


if __name__ == "__main__":
    import json as _j
    d = _j.load(open(os.path.join(OUT, "verify_Z3xZ6_1_1.json")))
    only = {c["model"] for c in d["loophole_cells"]}
    only_files = {c["file"] for c in d["loophole_cells"]}
    globals()["only_files"] = only_files
    print("targeting %d loophole models in %d files" % (len(only), len(only_files)))
    run(only=only, maxcells=40)

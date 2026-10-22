"""root_lift GUARD 1(d) -- reproduce orbifolder's printed twisted spectrum from the
projection formula, to validate the centraliser and the conventions BEFORE any verdict.

Method.  For a model whose shift/Wilson data are on disk (`begin model` block) and whose
orbifolder output is on disk, take each twisted cell (k, l, n_1..n_6) that orbifolder
populated, and

  1. check  V_g = k V_1 + l V_2 + sum n_i A_i  ==  p_sh  (mod Lambda_{E8xE8})  for every
     printed field at that cell   -- this validates conventions C2 and C5 against
     orbifolder's OWN printed left-moving momenta, index by index;
  2. enumerate every p in V_g + Lambda with |p|^2 = 2B (oscillator-free budget);
  3. impose the projection (convention C7) with the printed right-mover q_sh, over the
     centraliser Z_S(g) computed in rl_lib;
  4. compare the survivor set with the weights orbifolder actually printed at that cell
     with that q_sh -- every printed weight must survive, and every survivor must lie in
     the Weyl orbit (under the UNBROKEN 4d gauge group) of a printed weight.

Step 4 is the sharp one: if the centraliser were too large, printed states would fail; if
too small, extra Weyl orbits would appear.

Output: data/root_lift/guard1d_projection.json
"""
from __future__ import annotations

import gzip
import json
import os
import re
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.dirname(HERE)
CAMP = os.path.dirname(CODE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(CODE, "su5_ten_bound"))

import rl_lib as R                                                      # noqa: E402
import s5lib as S5                                                      # noqa: E402
import shell as SH                                                      # noqa: E402

OUT = os.path.join(CAMP, "data", "root_lift")
os.makedirs(OUT, exist_ok=True)


# ------------------------------------------------------------------ model file parsing

def parse_models(path):
    out = {}
    txt = open(path).read()
    for blk in txt.split("begin model")[1:]:
        lab = re.search(r"Label:(\S+)", blk).group(1)
        sg = re.search(r"SpaceGroup:\S*Geometry_(\S+?)\.txt", blk).group(1)
        body = blk.split("Shifts and Wilsonlines:")[1].split("end model")[0]
        rows = []
        for ln in body.strip().splitlines():
            ln = ln.strip()
            if not ln:
                continue
            rows.append([F(x.strip()) for x in ln.split(",")])
        assert len(rows) == 8 and all(len(r) == 16 for r in rows), lab
        out[lab] = dict(geom=sg, V1=rows[0], V2=rows[1], A=rows[2:8])
    return out


# ------------------------------------------------------------------ spectrum parsing

_FLD = re.compile(r"\n\s+F_(\d+)\s*\n(.*?)(?=\n\s+F_\d+\s*\n|\Z)", re.S)


def parse_spectrum_block(txt):
    """Parse one model's printed spectrum block into records."""
    recs = []
    for m in _FLD.finditer(txt):
        b = m.group(2)
        kl = re.search(r"sector \(k,l\)\s*:\s*\((\d+),\s*(\d+)\)", b)
        fp = re.search(r"fixed point n_a\s*:\s*\(([^)]*)\)", b)
        ps = re.search(r"left-moving p_sh\s*:\s*\(([^)]*)\),\s*\(([^)]*)\)", b)
        qs = re.search(r"right-moving q_sh\s*:\s*\(([^)]*)\)", b)
        rep = re.search(r"representation\s*:\s*\(([^)]*)\)_(\w)", b)
        osc = re.search(r"oscillators\s*:\s*(.*)", b)
        if not (kl and fp and ps and qs):
            continue
        p = [F(x.strip()) for x in ps.group(1).split(",")] + \
            [F(x.strip()) for x in ps.group(2).split(",")]
        recs.append(dict(
            F="F_%s" % m.group(1),
            kl=(int(kl.group(1)), int(kl.group(2))),
            fp=tuple(int(x) for x in fp.group(1).split(",")),
            p=p,
            q=[F(x.strip()) for x in qs.group(1).split(",")],
            rep=rep.group(1).strip() if rep else None,
            chir=rep.group(2) if rep else None,
            osc=(osc.group(1).strip() if osc else "no")))
    return recs


_MARK = re.compile(r"/>\s*/([A-Za-z0-9_]+)>")


def split_models(txt):
    """One block per model, cut at the `/> /LABEL>` prompt marks (same convention as
    code/z3xz6_sweep/z36sw_lib.py::split_models, widened to non-MSSM labels)."""
    marks = [(m.start(), m.group(1)) for m in _MARK.finditer(txt)]
    out = {}
    for i, (pos, lab) in enumerate(marks):
        end = marks[i + 1][0] if i + 1 < len(marks) else len(txt)
        out[lab] = txt[pos:end]
    return out


# ------------------------------------------------------------------ the projection

def centralizer_elements(S, ce, nmax=3):
    """Explicit finite list of centraliser elements (k,l,n) of the constructing element
    `ce`, with the pure-translation part truncated to coefficients in [0, nmax) along the
    integer-kernel basis.  All Wilson lines have order <= nmax on these geometries, so the
    truncation loses no VALUE of V_h."""
    cd = S.centralizer_data(ce)
    kern = cd["kern"]
    els = []
    import itertools
    for (kh, lh) in cd["PG"]:
        n0 = cd["part"][(kh, lh)]
        for coef in itertools.product(range(nmax), repeat=len(kern)):
            n = n0.copy()
            for c, kv in zip(coef, kern):
                n = n + c * kv
            els.append((kh, lh, tuple(int(x) for x in n)))
    return els


def group_from_roots(roots):
    """Semi-simple name of the group generated by an explicit root set (s5lib pieces)."""
    if not roots:
        return "U(1)^16"
    S = S5.simple_roots(roots)
    comps, adj = S5._components(S)
    facs = [S5.group_name(*S5._classify(c, adj, S)) for c in comps]
    nm = " x ".join(sorted(facs)) or "-"
    if len(S) < 16:
        nm += " x U(1)^%d" % (16 - len(S))
    return nm


def Vh(mod, k, l, n):
    V1, V2, A = mod["V1"], mod["V2"], mod["A"]
    return [k * V1[i] + l * V2[i] + sum(n[j] * A[j][i] for j in range(6))
            for i in range(16)]


def vh(S, k, l):
    return tuple(k * S.v1[a] + l * S.v2[a] for a in range(3))


def survives(p, q, mod, S, ce, els):
    """Convention C7 -- exact over Q."""
    kg, lg = ce[0], ce[1]
    ng = ce[2:8]
    Vg = Vh(mod, kg, lg, ng)
    vg = vh(S, kg, lg)
    for (kh, lh, nh) in els:
        Vhh = Vh(mod, kh, lh, nh)
        vhh = vh(S, kh, lh)
        lm = sum(p[i] * Vhh[i] for i in range(16))
        lm -= F(1, 2) * (sum(Vg[i] * Vhh[i] for i in range(16))
                         - sum(F(vg[a]) * F(vhh[a]) for a in range(3)))
        # right-mover: q_sh is printed as (q_0, q_1, q_2, q_3) with q_0 the helicity slot
        rm = kh * sum(F(S.v1[a]) * q[a + 1] for a in range(3)) \
            + lh * sum(F(S.v2[a]) * q[a + 1] for a in range(3))
        if (lm - rm) % 1 != 0:
            return False
    return True


def main():
    R.selftest(verbose=False)
    print("[guard 1a-1c] PASS on all seven campaign geometries\n")

    models = {}
    models.update(parse_models(os.path.join(CAMP, "data", "gate_su5_construction",
                                            "gate_models.txt")))
    txt = gzip.open(os.path.join(CAMP, "data", "gate_su5_construction",
                                 "gate_v1.txt.gz"), "rt", errors="replace").read()
    blocks = split_models(txt)

    report = {}
    for lab in ["CB1534", "CB1041"]:
        if lab not in models or lab not in blocks:
            print("skip", lab, lab in models, lab in blocks)
            continue
        mod = models[lab]
        S = R.SpaceGroup(mod["geom"])
        recs = [r for r in parse_spectrum_block(blocks[lab]) if r["chir"] == "l"]
        cells = {}
        for r in recs:
            cells.setdefault((r["kl"], r["fp"]), []).append(r)

        # ---- step 1: V_g == p_sh (mod Lambda) on every printed field
        bad_p1 = []
        for (kl, fp), rs in cells.items():
            Vg = Vh(mod, kl[0], kl[1], fp)
            for r in rs:
                d = [r["p"][i] - Vg[i] for i in range(16)]
                if not R.in_lambda(d):
                    bad_p1.append((lab, kl, fp, r["F"]))
        print("%s: cells=%d printed left-chiral fields=%d   P1 (V_g == p_sh mod Lambda) "
              "violations=%d" % (lab, len(cells), len(recs), len(bad_p1)))

        # ---- unbroken 4d roots, from ALL EIGHT generators (guard 2 / CLAUDE 5.21)
        allr = S5.e8xe8_roots()
        gens = [mod["V1"], mod["V2"]] + mod["A"]
        R4 = [a for a in allr if all(S5.dot(a, g) % 1 == 0 for g in gens)]
        grp4name = group_from_roots(R4)

        # ---- steps 2-4, on the loophole sectors
        det = []
        for (kl, fp), rs in sorted(cells.items()):
            ce = (kl[0], kl[1], ) + tuple(fp)
            lifts = (S.lift_verdict(ce)["index"] == 1)
            els = centralizer_elements(S, ce)
            Vg = Vh(mod, kl[0], kl[1], fp)
            v = vh(S, kl[0], kl[1])
            dc = sum((F(x) % 1) * (1 - (F(x) % 1)) for x in v) / 2
            c = 2 * (1 - dc)
            cand = SH.shell(Vg, c)
            for q in {tuple(r["q"]) for r in rs if r["osc"] in ("no", "")}:
                printed = [r for r in rs if tuple(r["q"]) == q and r["osc"] in ("no", "")]
                surv = [p for p in cand if survives(p, list(q), mod, S, ce, els)]
                pset = {tuple(p) for p in surv}
                miss = [r["F"] for r in printed if tuple(r["p"]) not in pset]
                # Weyl orbits of the survivors under the unbroken 4d group
                orbits = []
                left = list(pset)
                while left:
                    seed = left[0]
                    orb = {seed}
                    frontier = [seed]
                    while frontier:
                        x = frontier.pop()
                        for a in R4:
                            n2 = S5.dot(a, a)
                            y = tuple(F(x[i]) - 2 * S5.dot(x, a) / n2 * F(a[i])
                                      for i in range(16))
                            if y in pset and y not in orb:
                                orb.add(y)
                                frontier.append(y)
                    orbits.append(orb)
                    left = [z for z in left if z not in orb]
                det.append(dict(kl=list(kl), fp=list(fp), q=[str(x) for x in q],
                                lifts=lifts,
                                n_candidates=len(cand), n_survivors=len(surv),
                                n_printed_records=len(printed),
                                n_weyl_orbits=len(orbits),
                                orbit_sizes=sorted(len(o) for o in orbits),
                                printed_missing=miss,
                                printed_reps=[r["rep"] for r in printed],
                                match=(not miss and len(orbits) == len(printed))))
                if (not det[-1]["match"]) or kl in [(0, 1), (0, 5)]:
                    print("   %-8s fp=%s lift=%-5s cand=%3d surv=%2d printed=%d orbits=%d "
                          "sizes=%s match=%s%s"
                          % (str(kl), fp, lifts, len(cand), len(surv),
                             len(printed), len(orbits), sorted(len(o) for o in orbits),
                             det[-1]["match"], "" if not miss else "  MISSING %s" % miss))
        report[lab] = dict(geom=mod["geom"], n_cells=len(cells), n_fields=len(recs),
                           P1_violations=len(bad_p1),
                           gauge4d=grp4name,
                           n_roots_4d=len(R4), loophole_cells=det,
                           all_match=all(d["match"] for d in det))
        nl = sum(1 for d in det if d["lifts"])
        multi = [d for d in det if d["n_printed_records"] > 1]
        print("   -> cells x q tested: %d   all_match = %s   (4d roots: %d, group %s)"
              % (len(det), report[lab]["all_match"], len(R4), grp4name))
        print("      lift-sector (cell,q) groups: %d   non-lift: %d" % (nl, len(det) - nl))
        print("      (cell,q) groups with >1 printed record: %d  [of which lift-sector: %d]"
              % (len(multi), sum(1 for d in multi if d["lifts"])))
        print()

    with open(os.path.join(OUT, "guard1d_projection.json"), "w") as fh:
        json.dump(report, fh, indent=1, default=str)
    print("wrote", os.path.join(OUT, "guard1d_projection.json"))


if __name__ == "__main__":
    main()

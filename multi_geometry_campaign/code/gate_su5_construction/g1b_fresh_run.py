"""V1 -- build a FRESH orbifolder input deck for CB1041 (and the gate's replacement
flagship candidates) from the EMBEDDING DATA ONLY.

The model block is written here from (V1, V2, W3) as recorded in the phase's flagship
dossier / build JSON; nothing is copied from data/su5_construction/s5_models.txt.  The
resulting spectrum is then compared record-by-record against the stored segment output.
"""
from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gs_lib as G

WANT = ["CB1041", "CB1534", "CB1370", "CB2307"]


def block(label, V1, V2, W3):
    Z = ", ".join(["0"] * 16)
    rows = [", ".join(V1), ", ".join(V2), Z, Z, ", ".join(W3), ", ".join(W3), Z, Z]
    return ("begin model\n"
            "Label:%s\n"
            "SpaceGroup:Geometry/Geometry_Z3xZ6_1_1.txt\n"
            "Lattice:E8xE8\n"
            "Shifts and Wilsonlines:\n" % label
            + "\n".join(rows) + "\nend model\n")


def main():
    src = json.load(open(os.path.join(G.PHASE_DATA, "s12_flagship.json")))
    dossiers = {}
    if isinstance(src.get("all"), list):
        for d in src["all"]:
            if isinstance(d, dict) and "model" in d:
                dossiers[d["model"]] = d
    dossiers[src["flagship"]["model"]] = src["flagship"]
    build = json.load(open(os.path.join(G.PHASE_DATA, "s5_build.json")))
    bym = {m.get("label", m.get("model")): m for m in build["models"]}
    out, used = [], {}
    for lab in WANT:
        if lab in dossiers:
            d = dossiers[lab]
            V1, V2, W3 = d["V1"], d["V2"], d["W3"]
        else:
            m = bym[lab]
            V1, W3 = m["V1"], m["W3"]
            V2 = build["V2_target"]
        used[lab] = dict(V1=V1, V2=V2, W3=W3)
        f1 = [G.rat(x) for x in V1]
        f2 = [G.rat(x) for x in V2]
        f3 = [G.rat(x) for x in W3]
        Z = [G.rat("0")] * 16
        ok, bad = G.mi_check(f1, f2, [Z, Z, f3, f3, Z, Z])
        print(lab, "MI:", ok, bad, "orders", G.order_of(f1), G.order_of(f2),
              G.order_of(f3))
        out.append(block(lab, V1, V2, W3))
    p = os.path.join(G.OUT, "gate_models.txt")
    open(p, "w", newline="\n").write("".join(out))
    json.dump(used, open(os.path.join(G.OUT, "gate_models.json"), "w"), indent=1)
    print("wrote", p)


if __name__ == "__main__":
    main()

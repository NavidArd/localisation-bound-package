"""Referee #8 item 2.3: compare the V2-representative brothers NU* of CB1534 with the nine
(a,b) reference spectra AB**, record by record (sector, fixed point, p_sh, oscillators,
right-mover, representation).  Output: review/data/paperA_referee8/b11_compare.json"""
import json, os, re, sys
_PKG = os.environ.get("PAPERA_PKG") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b4_triality_nets as T
D = T.D
txt = open(os.path.join(D, "b11_out.txt"), errors="replace").read()
blocks = T.split(txt, pat=r"/> /((?:NU|AB)\w+)>")
def keyset(blk):
    spec = T.SW.parse_spectrum(blk.split("print summary")[0])
    return Counter((r["kl"], r["fp"], tuple(map(str, r["p_sh"])), r["osc_raw"],
                    tuple(map(str, r["qsh"])), r["rep"]) for r in spec.values()
                   if r["chir"] == "l"), T.nets(spec.values())
ks = {lab: keyset(b) for lab, b in blocks}
out = {}
for lab, (k, n) in ks.items():
    if lab.startswith("NU"):
        match = [r for r, (k2, _n) in ks.items() if r.startswith("AB") and k2 == k]
        out[lab] = dict(matches=match, nets=n, n_records=sum(k.values()))
json.dump(out, open(os.path.join(D, "b11_compare.json"), "w"), indent=1)
for l, v in out.items(): print(l, v)

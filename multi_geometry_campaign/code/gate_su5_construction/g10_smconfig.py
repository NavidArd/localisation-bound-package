"""V5 attacks 5+6 -- cross-check the gate's role census against orbifolder's OWN
SMConfig identification, on every CB model where the finder actually installed one."""
import glob, gzip, json, os, re, sys
from collections import Counter
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gs_lib as G
RAW = os.path.join(G.PHASE_DATA, "raw")


def do(p):
    out = []
    txt = gzip.open(p, "rt", errors="replace").read()
    for lab, blk in G.split_models(txt):
        m = re.search(r'Gauge group in vev-configuration "SMConfig":\s*(.+)', blk)
        if not m:
            continue
        tail = blk.split("use config(SMConfig)")[-1]
        out.append((lab, m.group(1).strip(), tail[:4000]))
    return out


def main():
    paths = sorted(glob.glob(os.path.join(RAW, "*.txt.gz")))
    with Pool(6) as pl:
        res = [x for c in pl.map(do, paths) for x in c]
    print("CB models where orbifolder installed an SMConfig:", len(res))
    g5 = {r["model"]: r for r in
          json.load(open(os.path.join(G.OUT, "g5_battery.json")))["rows"]}
    rows = []
    for lab, gg, tail in res:
        r = g5[lab]
        nums = re.findall(r"#\(([^)]*)\)\s*=\s*(-?\d+)", tail)
        rows.append(dict(model=lab, gauge=gg, orb_summary=dict(nums),
                         gate_net=r["net_aware"], phase_blind_netQ=r["netQ_blind"],
                         trY=r["trY"], hits=r["hits_aware"]))
        print("  ", lab, "|", gg)
        print("     orbifolder summary:", dict(nums))
        print("     gate net (triality-aware):", r["net_aware"],
              " phase blind netQ:", r["netQ_blind"], " trY:", r["trY"])
        print("     tail sample:", " ".join(tail.split())[:300])
    json.dump(rows, open(os.path.join(G.OUT, "g10_smconfig.json"), "w"), indent=1)


if __name__ == "__main__":
    main()

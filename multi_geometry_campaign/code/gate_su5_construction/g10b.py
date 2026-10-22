import glob, gzip, os, sys
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gs_lib as G
RAW = os.path.join(G.PHASE_DATA, "raw")
def do(p):
    out = []
    txt = gzip.open(p, "rt", errors="replace").read()
    for lab, blk in G.split_models(txt):
        if "No vev-configuration identified" not in blk:
            i = blk.find("analyze config")
            out.append((lab, " ".join(blk[i:i + 2500].split())))
    return out
if __name__ == "__main__":
    paths = sorted(glob.glob(os.path.join(RAW, "*.txt.gz")))
    with Pool(6) as pl:
        res = [x for c in pl.map(do, paths) for x in c]
    print('CB models WITHOUT "No vev-configuration identified":', len(res))
    for lab, t in res[:8]:
        print("==", lab)
        print(t[:1100])
        print()

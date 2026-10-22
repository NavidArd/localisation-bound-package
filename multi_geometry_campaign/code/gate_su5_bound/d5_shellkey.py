"""GATE -- the cell-key question, answered at the level the theorem actually works at.

CLAUDE.md 5.22 warns that a cell key which SPLITS one physical fixed point into two would
HIDE a same-cell doubling -- the dangerous direction for an absence claim.  For THIS result
the worry is answerable structurally, because the theorem's hypothesis H3 is really "same
LOCAL SHIFT V_g", not "same printed key": V_g == p_sh (mod Lambda), so two states sharing a
V_g are excluded by the theorem whether or not they share a printed key, and two states with
different V_g are at genuinely different fixed points.

The sharp form: at a local shift V_g, the COMPLETE massless shell {p in V_g+Lambda :
|p|^2 = 2B} may contain at most ONE state of any given 10-role key when R < 1/2.  If that
holds, then no key convention -- merging, splitting, orbit-vs-representative -- can hide a
10-role doubling, because the union over ALL cells sharing that V_g still contains one state.

This script tests it directly on Z2xZ4_1_1 (where every 10-role budget is protected).

Output: data/gate_su5_bound/d5_shellkey_Z2xZ4.json
"""
import os
for _v in ("OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS"): os.environ[_v]="1"
import glob, gzip, json, sys, multiprocessing as mp
from collections import Counter, defaultdict
from fractions import Fraction as F
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); CAMP=os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(CAMP,"code","su5_ten_bound")); sys.path.insert(0,HERE)
OUT=os.path.join(CAMP,"data","gate_su5_bound")
R10={"Q","uc","ec","Qbar","ucbar","ecbar"}
CENSUS=os.path.join(CAMP,"data","z2xz4_sweep","raw","census")
_S=None;_R2=None
def _init():
    global _S,_R2
    import s5lib; _S=s5lib
    _R2=np.array([[int(2*x) for x in a] for a in _S.e8xe8_roots()],dtype=np.int64)
def worker(path):
    if _S is None: _init()
    import d7_shell as SH
    sys.path.insert(0, os.path.join(CAMP,"code","z2xz4_sweep")); Z=__import__("z24sw_lib")
    S=_S; allr=S.e8xe8_roots()
    acc=dict(models=0,cells_tested=0,max10=Counter(),n10=Counter(),inshell=Counter(),bad=[])
    txt=gzip.open(path,"rt",errors="replace").read(); nm=0
    for lab,blk in Z.split_models(txt):
        if nm>=12: break
        pre=blk.split("Analyzing")[0]
        try:
            spec=Z.parse_spectrum(pre,M=2,N=4); factors=Z.parse_gauge_factors(blk)
            cfgs=Z.parse_smconfig(blk)
            if not cfgs or cfgs[0][0] is None or not spec: continue
            tY,tbl=cfgs[0]; printed=defaultdict(int)
            for row in tbl: printed[(row[0],row[1],row[2])]+=row[3]
            pos=Z.observable_positions(spec,tY,factors,dict(printed))
            if pos is None: continue
            roles,_=Z.role_census(spec,tY,*pos)
        except Exception: continue
        nm+=1; acc["models"]+=1
        cells=defaultdict(list)
        for l,r in spec.items():
            if r["p_sh"] is None or r["chir"]!="l": continue
            cells[(tuple(r["kl"]),tuple(r["fp"]) if r["fp"] else None)].append(r)
        inter=None
        for key,recs in cells.items():
            V=recs[0]["p_sh"]; D=1
            for x in V: D=D*x.denominator//np.gcd(D,x.denominator)
            m=((_R2@np.array([int(x*D) for x in V],dtype=np.int64))%(2*D))==0
            inter=m if inter is None else (inter&m)
        R4=[allr[i] for i in np.nonzero(inter)[0]] if inter is not None else []
        if not R4: continue
        Sp4=S.simple_roots(R4); comps4,adj4=S._components(Sp4)
        qr=[r for r in spec.values() if roles.get(r["label"],{}).get("role") in ("Q","Qbar")]
        dr=[r for r in spec.values() if roles.get(r["label"],{}).get("role") in ("dc","dcbar")]
        lr=[r for r in spec.values() if roles.get(r["label"],{}).get("role") in ("L","Hu")]
        f3=f2=None
        for c_ in comps4:
            kind,rk=S._classify(c_,adj4,Sp4); so=[Sp4[i] for i in c_]
            if kind=="A" and rk==2 and (qr or dr):
                if all(any(S.dot(x["p_sh"],a)!=0 for a in so) for x in qr+dr): f3=so if f3 is None else "X"
            if kind=="A" and rk==1 and (qr or lr):
                if all(any(S.dot(x["p_sh"],a)!=0 for a in so) for x in qr+lr): f2=so if f2 is None else "X"
        if not (isinstance(f3,list) and isinstance(f2,list)): continue
        smb=list(f3)+list(f2)+[tuple(tY)]
        done=set()
        for key,recs in cells.items():
            kl,fp=key
            if tuple(kl)==(0,0): continue
            tens=[r for r in recs if roles.get(r["label"],{}).get("role") in R10]
            if not tens: continue
            if any((r.get("osc") or r.get("osc_raw") or "no").strip() not in ("no","") for r in tens): continue
            V=list(recs[0]["p_sh"]); c=S.dot(tens[0]["p_sh"],tens[0]["p_sh"])
            sig=(S.reduce_mod_lattice(V),str(c))
            if sig in done: continue
            done.add(sig)
            sts=SH.shell(V,c); ss=set(sts)
            acc["inshell"][all(tuple(r["p_sh"]) in ss for r in tens)]+=1
            c10=Counter()
            for x in sts:
                l3=tuple(S.dot(x,a) for a in f3); l2=tuple(S.dot(x,a) for a in f2)
                if any(t.denominator!=1 for t in l3+l2): continue
                if S.proj_norm2(x,smb)!=F(6,5): continue
                c10[str((tuple(int(t) for t in l3),tuple(int(t) for t in l2),str(S.dot(tY,x))))]+=1
            mx=max(c10.values()) if c10 else 0
            acc["cells_tested"]+=1; acc["max10"][mx]+=1; acc["n10"][sum(c10.values())]+=1
            if mx>1 and len(acc["bad"])<10:
                acc["bad"].append(dict(model=lab,kl=list(kl),fp=list(fp) if fp else None,max10=mx))
    for k in ("max10","n10","inshell"): acc[k]=dict(acc[k])
    return acc
if __name__=="__main__":
    files=sorted(glob.glob(os.path.join(CENSUS,"*.txt.gz")))[:10]
    with mp.Pool(6,initializer=_init) as p: accs=p.map(worker,files)
    out=dict(models=0,cells_tested=0); C={k:Counter() for k in ("max10","n10","inshell")}; bad=[]
    for a in accs:
        out["models"]+=a["models"]; out["cells_tested"]+=a["cells_tested"]
        for k in C: C[k].update(a[k])
        bad+=a["bad"]
    for k in C: out[k]=dict(C[k])
    out["bad"]=bad[:20]
    json.dump(out,open(os.path.join(OUT,"d5_shellkey_Z2xZ4.json"),"w"),indent=1,default=str)
    print(json.dumps(out,indent=1,default=str))

#!/usr/bin/env python3
import os,sys,math,subprocess,glob,re,json
import eddy_mk, post_eddy
root=eddy_mk.root
W=post_eddy.W; SIG=post_eddy.SIG
J0=1.6667e7   # source current density for 6.0 A leg current

VERS={"fullcu":(0.0,0.0),"topcut":(3.3e-3,0.0),"dualcut":(3.3e-3,3.3e-3)}
SCALES=[1.0,0.6,0.4]

def solve(ver,xct,xcb,scale):
    tag="%s_s%g"%(ver,scale)
    d=os.path.join(root,tag)
    d=eddy_mk.build(tag,xct,xcb,scale,J0)
    vt=glob.glob(os.path.join(d,"mesh","case_t*.vtu"))
    log=open(os.path.join(d,"run.log")).read()
    nrm=re.search(r"SS \(ITER=1\) \(NRM,RELC\): \( *([\d.E+-]+)",log)
    nrm=float(nrm.group(1)) if nrm else float('nan')
    hdr=[int(x) for x in open(os.path.join(d,"mesh","mesh.header")).readline().split()]
    nelem=hdr[1] if hdr else 0
    res={"tag":tag,"ver":ver,"scale":scale,"nrm":nrm,"nelem":nelem}
    if vt:
        o,g,Ar,Ai,pts,ids=post_eddy.analyze(vt[0])
        res["P_top"]=o["top"][0]; res["P_bot"]=o["bot"][0]
        res["Jmax_top"]=o["top"][1]; res["Jmax_bot"]=o["bot"][1]
        res["P_tot"]=o["top"][0]+o["bot"][0]
    return res

if __name__=="__main__":
    out=[]
    for ver,(xct,xcb) in VERS.items():
        for s in SCALES:
            r=solve(ver,xct,xcb,s)
            print(json.dumps(r)); out.append(r)
    json.dump(out,open(os.path.join(root,"results_raw.json"),"w"),indent=1)

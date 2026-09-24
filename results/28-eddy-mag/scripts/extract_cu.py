#!/usr/bin/env python3
"""Extract REAL copper pour polygons from epru for 3 versions; compute areas."""
import sys, math, json
sys.path.insert(0,"/mnt/raid10/sim-work/tps62933/fea3")
from parse_lib import load_pcbs, parse_path
F="/mnt/raid10/sim-work/tps62933/epro/pourSim.epru"
MIL=25.4e-6; L2=(245.0,-660.0)

def bbox(pts):
    xs=[p[0] for p in pts]; ys=[p[1] for p in pts]
    return (min(xs),max(xs),min(ys),max(ys))

def area_poly(pts):
    # shoelace (mil^2)
    a=0.0; n=len(pts)
    for i in range(n):
        x1,y1=pts[i]; x2,y2=pts[(i+1)%n]
        a+=x1*y2-x2*y1
    return abs(a)/2.0

pcbs=load_pcbs(F)
print("num PCB docs:",len(pcbs))
for i,d in enumerate(pcbs):
    board=None
    for t,pj in d["items"]:
        if t=="POLY" and pj.get("layerId")==11:
            board=bbox(parse_path(pj["path"]))
    # top-layer gnd pour = the copper (layerId 1 POUR named GND?) ; find pour polygons
    pours=[]
    for t,pj in d["items"]:
        if t=="POUR":
            pts=parse_path(pj["path"])
            pours.append((pj.get("layerId"),pj.get("name"),pj.get("netName"),len(pts),bbox(pts) if pts else None))
    regs=[]
    for t,pj in d["items"]:
        if t.upper().startswith("REGION") and "COPPER" in (pj.get("prohibitType") or []):
            pts=parse_path(pj["path"])
            regs.append((pj.get("layerId"),pj.get("name"),bbox(pts) if pts else None, area_poly(pts) if pts else 0))
    print("\n=== doc %d uuid=%s type=%s ==="%(i,d.get("uuid"),d.get("type")))
    print("  board bbox mil:",board)
    print("  pours(layer,name,net,npts,bbox):")
    for p in pours: print("    ",p)
    print("  REGION(COPPER):")
    for r in regs: print("    ",r)

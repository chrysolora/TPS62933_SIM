#!/usr/bin/env python3
"""geom_verify: read the three copper versions directly from pourSim.epru and
compare the ACTUAL cutout geometry against the rectangles ASSUMED by thermal_3x3.
Read-only. Writes only into geom_verify/.
"""
import json, math, os, sys
import numpy as np

EPRU = "/mnt/raid10/sim-work/tps62933/epro/pourSim.epru"
OUT  = "/mnt/raid10/sim-work/tps62933/geom_verify"
os.makedirs(OUT, exist_ok=True)
MIL2MM = 0.0254
BOARD_W_MM, BOARD_H_MM = 26.50, 50.50
BOARD_H_MIL = BOARD_H_MM / MIL2MM   # 1988.19 mil

# ---------------------------------------------------------------- load docs
def load_docs(f):
    lines = open(f, encoding="utf-8", errors="replace").read().split("\n")
    docs=[]; cur=None
    for l in lines:
        if not l.strip() or "||" not in l: continue
        h,p = l.split("||",1)
        if p.endswith("|"): p=p[:-1]
        try: hj=json.loads(h)
        except: continue
        t=hj.get("type")
        if t=="DOCHEAD":
            try: pj=json.loads(p)
            except: pj={}
            cur={"type":pj.get("docType"),"uuid":pj.get("uuid"),"items":[]}; docs.append(cur)
        elif cur is not None:
            try: pj=json.loads(p)
            except: pj={}
            cur["items"].append((t,pj,hj))
    return docs

# ---------------------------------------------------------------- geometry
def r_rect(x,y,w,h,rot):
    corners=[(0,0),(w,0),(w,-h),(0,-h)]
    th=math.radians(rot); ct,st=math.cos(th),math.sin(th)
    return [(x+dx*ct-dy*st, y+dx*st+dy*ct) for dx,dy in corners]

def parse_path(path):
    if not isinstance(path,list): return []
    if len(path)==1 and isinstance(path[0],list): path=path[0]
    if path and path[0]=="R":
        x,y,w,h=path[1],path[2],path[3],path[4]
        rot=path[5] if len(path)>5 else 0
        return r_rect(x,y,w,h,rot)
    pts=[]; i=0; n=len(path)
    while i<n:
        tok=path[i]
        if isinstance(tok,str):
            if tok=="ARC":
                i+=1
                if i<n and isinstance(path[i],(int,float)): i+=1   # skip angle
                continue
            i+=1; continue
        else:
            if i+1<n and isinstance(path[i+1],(int,float)):
                pts.append((float(path[i]),float(path[i+1]))); i+=2
            else: i+=1
    return pts

def poly_area_mil(pts):
    a=0.0
    for i in range(len(pts)):
        x1,y1=pts[i]; x2,y2=pts[(i+1)%len(pts)]
        a += x1*y2 - x2*y1
    return abs(a)/2.0

def bbox(pts):
    xs=[p[0] for p in pts]; ys=[p[1] for p in pts]
    return min(xs),min(ys),max(xs),max(ys)

def centroid_poly(pts):
    # area centroid (mil)
    a=0.0; cx=0.0; cy=0.0; n=len(pts)
    for i in range(n):
        x1,y1=pts[i]; x2,y2=pts[(i+1)%n]
        cross=x1*y2-x2*y1; a+=cross; cx+=(x1+x2)*cross; cy+=(y1+y2)*cross
    a*=0.5
    if abs(a)<1e-9: return (sum(p[0] for p in pts)/n, sum(p[1] for p in pts)/n)
    return (cx/(6*a), cy/(6*a))

docs=load_docs(EPRU)
pcbs=[d for d in docs if d["type"]=="PCB"]
# order: PCB8 (full), PCB8_1 (top), PCB8_2 (dual)
names={}
for i,d in enumerate(pcbs):
    title=[p for t,p,_ in d["items"] if t=="META"][0].get("title")
    names[i]=title
print("PCBs:", {i:names[i] for i in range(len(pcbs))})

LAYERS={}  # layerId->type
for t,p,hj in pcbs[0]["items"]:
    if t=="LAYER":
        lid=json.loads(hj["id"])[1] if isinstance(hj.get("id"),str) else None
        LAYERS[lid]=p.get("layerType")

def items_of(doc,typ): return [p for t,p,_ in doc["items"] if t==typ]

result={"note":"cutout geometry read directly from pourSim.epru; mil->mm x0.0254","versions":{}}

# pour outlines (definitions) per version per layer  (should be identical)
pour_defs={}
for i,d in enumerate(pcbs):
    per={}
    for p in items_of(d,"POUR"):
        lid=p.get("layerId"); nm=p.get("name"); net=p.get("netName")
        pts=parse_path(p.get("path"))
        rec={"name":nm,"net":net,"layer":lid,"layerType":LAYERS.get(lid),
             "bbox_mil":[round(v,3) for v in bbox(pts)],
             "area_mm2":round(poly_area_mil(pts)*MIL2MM**2,3)}
        per.setdefault(lid,[]).append(rec)
    pour_defs[i]=per

# prohibit regions per version
regions={}
for i,d in enumerate(pcbs):
    lst=[]
    for p in items_of(d,"REGION"):
        if p.get("regionType")!="PROHIBIT": continue
        lid=p.get("layerId"); pts=parse_path(p.get("path"))
        xs=[q[0] for q in pts]; ys=[q[1] for q in pts]
        bx=[min(xs),min(ys),max(xs),max(ys)]
        w=bx[2]-bx[0]; h=bx[3]-bx[1]
        cen=centroid_poly(pts)
        bb_cen=((bx[0]+bx[2])/2,(bx[1]+bx[3])/2)
        lst.append({
            "layer":lid,"layerType":LAYERS.get(lid),
            "name":p.get("name"),
            "n_vertices":len(pts),
            "bbox_mil":[round(v,3) for v in bx],
            "bbox_w_mm":round(w*MIL2MM,3),"bbox_h_mm":round(h*MIL2MM,3),
            "bbox_area_mm2":round(w*h*MIL2MM**2,3),
            "poly_area_mm2":round(poly_area_mil(pts)*MIL2MM**2,3),
            "centroid_mil":[round(cen[0],2),round(cen[1],2)],
            "bbox_center_mil":[round(bb_cen[0],2),round(bb_cen[1],2)],
            "vertices_mil":[[round(x,3),round(y,3)] for x,y in pts],
        })
    regions[i]=lst

print("\n=== PROHIBIT regions per version ===")
for i,lst in regions.items():
    print(f"\n-- PCB[{i}] {names[i]} : {len(lst)} prohibit regions")
    for r in lst:
        print("   layer=%s(%s) nv=%d bbox_mm=%.3fx%.3f poly_area=%.3f mm2 cen_mil=%s"%(
            r["layer"],r["layerType"],r["n_vertices"],r["bbox_w_mm"],r["bbox_h_mm"],
            r["poly_area_mm2"],r["centroid_mil"]))

# base regions (present in ALL versions) vs NEW (cutout)
base_ids=set()
for p,hj in [(p,hj) for t,p,hj in pcbs[0]["items"] if t=="REGION"]:
    pass
def region_ids(doc):
    out={}
    for t,p,hj in doc["items"]:
        if t=="REGION" and p.get("regionType")=="PROHIBIT":
            out[hj.get("id")]=p
    return out
ids=[region_ids(d) for d in pcbs]
base=set(ids[0].keys()) & set(ids[1].keys()) & set(ids[2].keys())
new1=set(ids[1].keys())-set(ids[0].keys())
new2=set(ids[2].keys())-set(ids[0].keys())
print("\nbase region ids (all versions):",len(base))
print("NEW in TopCutout (PCB8_1):",new1)
print("NEW in DualCutout(PCB8_2):",new2)

# ---- build the cutout records (the NEW regions), matched by id
def find_region(lst, target_layer):
    return [r for r in lst if r["layer"]==target_layer]

cut={}
for i,newset in [(1,new1),(2,new2)]:
    for r in regions[i]:
        # match by bbox matching a new id geometry (identify by centroid proximity)
        pass

# Simpler: the new regions in each version are the last entries (highest tickets).
# We identified them by id; recompute directly from the PCB docs preserving id.
def new_region_records(i, newids):
    recs=[]
    for t,p,hj in pcbs[i]["items"]:
        if t!="REGION" or p.get("regionType")!="PROHIBIT": continue
        if hj.get("id") not in newids: continue
        lid=p.get("layerId"); pts=parse_path(p.get("path"))
        xs=[q[0] for q in pts]; ys=[q[1] for q in pts]
        bx=[min(xs),min(ys),max(xs),max(ys)]; w=bx[2]-bx[0]; h=bx[3]-bx[1]
        cen=centroid_poly(pts); bb_cen=((bx[0]+bx[2])/2,(bx[1]+bx[3])/2)
        recs.append({"id":hj.get("id"),"layer":lid,"layerType":LAYERS.get(lid),
            "bbox_mil":[round(v,3) for v in bx],"bbox_w_mm":round(w*MIL2MM,3),
            "bbox_h_mm":round(h*MIL2MM,3),"bbox_area_mm2":round(w*h*MIL2MM**2,3),
            "poly_area_mm2":round(poly_area_mil(pts)*MIL2MM**2,3),
            "centroid_mil":[round(cen[0],2),round(cen[1],2)],
            "bbox_center_mil":[round(bb_cen[0],2),round(bb_cen[1],2)],
            "vertices_mil":[[round(x,3),round(y,3)] for x,y in pts]})
    return recs

cut_top = new_region_records(1,new1)
cut_dual= new_region_records(2,new2)

# --- assumed rectangles from thermal_3x3 REPORT (mm, in "flipped" coordinates)
# convert assumed mm center back to epru mil:  y_epru_mil = y_assumed_mil - BOARD_H_MIL
def assumed_to_epru(center_mm):
    xmil=center_mm[0]/MIL2MM; ymil=center_mm[1]/MIL2MM
    return (xmil, ymil-BOARD_H_MIL)
assumed=[
 {"label":"R1 L2","w_mm":6.60,"h_mm":7.00,"center_mm":(6.22,33.74),"center_mil_assumed":(245,-660)},
 {"label":"R2 L1","w_mm":4.35,"h_mm":4.56,"center_mm":(16.38,18.75)},
]
print("\nassumed R1 epru-center mil:",assumed_to_epru(assumed[0]["center_mm"]))
print("assumed R2 epru-center mil:",assumed_to_epru(assumed[1]["center_mm"]))

result["prohibit_regions_all"]={names[i]:regions[i] for i in regions}
result["cutouts_topcut"]=cut_top
result["cutouts_dualcut"]=cut_dual
result["assumed"]=assumed
result["pour_outlines"]={names[i]:{str(k):v for k,v in pour_defs[i].items()} for i in pour_defs}

json.dump(result, open(os.path.join(OUT,"copper_areas.json"),"w"), indent=1)

# ---------------------------------------------------------------- plotting
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MPoly
from matplotlib.path import Path

def to_pos(pts):
    # epru mil -> plotting mm with y flipped to positive-up like the report
    return [(x*MIL2MM, (BOARD_H_MIL-y)*MIL2MM) for x,y in pts]

fig,axes=plt.subplots(1,3,figsize=(16,7))
titles=["FullCopper (PCB8)","TopCutout (PCB8_1)","DualCutout (PCB8_2)"]
board=[(0,0),(BOARD_W_MM,0),(BOARD_W_MM,BOARD_H_MM),(0,BOARD_H_MM)]

# GND pour outline (top POUR4 = R(-5,70,1225,2080))
gnd_top=parse_path(["R",-5,70,1225,2080,0,0])
for ax,i in zip(axes,range(3)):
    ax.add_patch(MPoly(to_pos(board),closed=True,fill=False,ec="k",lw=1.5,label="board outline"))
    ax.add_patch(MPoly(to_pos(gnd_top),closed=True,fill=True,fc="#cfe8ff",ec="#3a7",lw=0.8,alpha=0.6,label="top GND pour (POUR4)"))
    # plot all prohibit regions (base, grey) and cutout (red)
    idset = ids[i]
    for t,p,hj in pcbs[i]["items"]:
        if t!="REGION" or p.get("regionType")!="PROHIBIT": continue
        pts=parse_path(p.get("path"))
        isnew = hj.get("id") in (new1|new2)
        lay=p.get("layerId")
        col = "red" if isnew else "grey"
        lw  = 1.8 if isnew else 0.8
        ax.add_patch(MPoly(to_pos(pts),closed=True,fill=isnew,fc="red",alpha=0.35,ec=col,lw=lw))
    # overlay ASSUMED rectangles (dashed) for L2 and L1
    for a in assumed:
        c=a["center_mm"]; w=a["w_mm"]; h=a["h_mm"]
        x0,y0=c[0]-w/2, c[1]-h/2
        ax.add_patch(MPoly([(x0,y0),(x0+w,y0),(x0+w,y0+h),(x0,y0+h)],closed=True,fill=False,ec="orange",lw=1.6,ls="--",label="assumed rect (thermal_3x3)"))
    ax.set_title(titles[i]); ax.set_aspect("equal"); ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
    ax.set_xlim(-2,BOARD_W_MM+2); ax.set_ylim(-2,BOARD_H_MM+2)
    ax.grid(alpha=0.2)
    if i==0:
        ax.legend(fontsize=7,loc="upper right")
fig.suptitle("Copper versions from pourSim.epru: real cutout (red) vs assumed rectangles (orange dashed)",fontsize=13)
plt.tight_layout(); plt.savefig(os.path.join(OUT,"copper_diff.png"),dpi=130)
print("wrote copper_diff.png")
print(json.dumps({"topcut":cut_top,"dualcut":cut_dual},indent=1))

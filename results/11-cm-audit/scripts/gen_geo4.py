#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gen_geo4.py -- known-good (gen_geo3) air geometry, with:
   - board core thickness tb (SW top face -> bottom GND) configurable;
   - optional top-layer coplanar GND strips (0V) around the SW bbox.
Single material/volume "air"; dielectric (FR4) applied in the .sif as eps_r(z).
"""
import os, sys
sys.path.insert(0, "/mnt/raid10/sim-work/tps62933/cm_pcb")
from geom_defs import SW_POLY_MIL, clip_halfplane, MIL
OUT = "/mnt/raid10/sim-work/tps62933/cm_audit"
BX0, BX1 = 0.0, 1043.3071*MIL
BY0, BY1 = -1958.189*MIL, 30.0*MIL
def sw_m(v):
    p = list(SW_POLY_MIL)
    if v in ("topcut","dualcut"): p = clip_halfplane(p, lambda q: q[1] <= -797.795)
    return [(x*MIL, y*MIL) for x, y in p]
def add_rect(A, base, x0,y0,x1,y1, z, lc):
    A("Point(%d) = {%.9g,%.9g,%.9g,%.9g};"%(base+0,x0,y0,z,lc))
    A("Point(%d) = {%.9g,%.9g,%.9g,%.9g};"%(base+1,x1,y0,z,lc))
    A("Point(%d) = {%.9g,%.9g,%.9g,%.9g};"%(base+2,x1,y1,z,lc))
    A("Point(%d) = {%.9g,%.9g,%.9g,%.9g};"%(base+3,x0,y1,z,lc))
    ls=[]
    for k in range(4):
        A("Line(%d) = {%d,%d};"%(base+10+k, base+k, base+(k+1)%4)); ls.append(base+10+k)
    return ls
def geo_text(v, d, tb, topgnd, clr=0.0003, lc_sw=0.0008, lc_far=0.03):
    sw = sw_m(v)
    swx0=min(p[0] for p in sw); swx1=max(p[0] for p in sw)
    swy0=min(p[1] for p in sw); swy1=max(p[1] for p in sw)
    L=[]; A=L.append
    A("SetFactory(\"OpenCASCADE\");"); A("Mesh.Algorithm3D = 1;")
    M=max(0.06,0.6*d); ax0,ax1=BX0-M,BX1+M; ay0,ay1=BY0-M,BY1+M; az0,az1=-d,0.03
    A("Box(1) = {%.9g,%.9g,%.9g,%.9g,%.9g,%.9g};"%(ax0,ay0,az0,ax1-ax0,ay1-ay0,az1-az0))
    pid=100
    for k,(x,y) in enumerate(sw): A("Point(%d) = {%.9g,%.9g,%.9g,%.9g};"%(pid+k,x,y,0.0,lc_sw))
    n=len(sw); sl=[]
    for k in range(n): A("Line(%d) = {%d,%d};"%(1000+k,pid+k,pid+(k+1)%n)); sl.append(1000+k)
    A("Curve Loop(5000) = {%s};"%",".join(map(str,sl))); A("Plane Surface(5001) = {5000};")
    lo=add_rect(A,10000,BX0,BY0,BX1,BY1,-tb,lc_far)
    A("Curve Loop(5001) = {%s};"%",".join(map(str,lo))); A("Plane Surface(5002) = {5001};")
    surfs=[5001,5002]
    if topgnd:
        hx0,hx1=swx0-clr,swx1+clr; hy0,hy1=swy0-clr,swy1+clr
        strips=[(BX0,BY0,BX1,hy0),(BX0,hy1,BX1,BY1),(BX0,hy0,hx0,hy1),(hx1,hy0,BX1,hy1)]
        for i,(x0,y0,x1,y1) in enumerate(strips):
            ll=add_rect(A,20000+i*100,x0,y0,x1,y1,0.0,lc_far)
            A("Curve Loop(%d) = {%s};"%(5010+i,",".join(map(str,ll)))); A("Plane Surface(%d) = {%d};"%(5011+i,5010+i))
            surfs.append(5011+i)
    A("BooleanFragments{ Volume{1}; Delete; }{ Surface{%s}; Delete; }"%",".join(map(str,surfs)))
    A("Physical Surface(\"sw\") = Surface In BoundingBox{%.9g,%.9g,%.9g,%.9g,%.9g,%.9g};"%(swx0-1e-5,swy0-1e-5,-1e-5,swx1+1e-5,swy1+1e-5,1e-5))
    A("Physical Surface(\"gnd\") = Surface In BoundingBox{%.9g,%.9g,%.9g,%.9g,%.9g,%.9g};"%(BX0-1e-5,BY0-1e-5,-tb-1e-5,BX1+1e-5,BY1+1e-5,-tb+1e-5))
    if topgnd:
        A("Physical Surface(\"topgnd\") = Surface In BoundingBox{%.9g,%.9g,%.9g,%.9g,%.9g,%.9g};"%(BX0-1e-5,BY0-1e-5,-1e-5,BX1+1e-5,BY1+1e-5,1e-5))
    A("Physical Surface(\"ref\") = Surface In BoundingBox{%.9g,%.9g,%.9g,%.9g,%.9g,%.9g};"%(ax0-1e-4,ay0-1e-4,az0-1e-4,ax1+1e-4,ay1+1e-4,az0+1e-4))
    A("Physical Volume(\"air\") = {1};")
    A("Field[1] = Box;"); A("Field[1].VIn = %.9g;"%lc_sw); A("Field[1].VOut = %.9g;"%lc_far)
    A("Field[1].XMin = %.9g;"%(swx0-0.004)); A("Field[1].XMax = %.9g;"%(swx1+0.004))
    A("Field[1].YMin = %.9g;"%(swy0-0.004)); A("Field[1].YMax = %.9g;"%(swy1+0.004))
    A("Field[1].ZMin = -0.006;"); A("Field[1].ZMax = 0.004;")
    A("Field[2] = Min;"); A("Field[2].FieldsList = {1};"); A("Background Field = 2;")
    A("Mesh.MeshSizeExtendFromBoundary = 0;"); A("Mesh.MeshSizeFromPoints = 0;"); A("Mesh.MeshSizeFromCurvature = 0;")
    return "\n".join(L)+"\n"
if __name__ == "__main__":
    import argparse
    ap=argparse.ArgumentParser()
    ap.add_argument("--var",default="fullcu"); ap.add_argument("--d",type=float,default=0.1)
    ap.add_argument("--tb",type=float,default=0.00143); ap.add_argument("--topgnd",action="store_true")
    ap.add_argument("--clr",type=float,default=0.0003); ap.add_argument("--tag",default="")
    a=ap.parse_args()
    fn=os.path.join(OUT,"geo4_%s_%s.geo"%(a.var,a.tag.lstrip("_")))
    open(fn,"w").write(geo_text(a.var,a.d,a.tb,a.topgnd,clr=a.clr)); print(fn)

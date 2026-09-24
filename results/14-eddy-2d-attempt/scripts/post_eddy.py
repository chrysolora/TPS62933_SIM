#!/usr/bin/env python3
"""Post-process Elmer 2D harmonic Az results -> P_eddy (W/m), |J| fields.
Loss: P = int 0.5 * rho * |J|^2 dA ; J = Js - i*omega*sigma*A  (sigma in Cu, Js in coil).
"""
import numpy as np, math, os, sys
import vtk
SIG=5.8e7; RHO=1.0/SIG; W=2*math.pi*805e3

def load(fn):
    r=vtk.vtkXMLUnstructuredGridReader(); r.SetFileName(fn); r.Update()
    return r.GetOutput()

def fields(g):
    n=g.GetNumberOfPoints()
    def arr(name):
        a=g.GetPointData().GetArray(name)
        if a is None: return np.zeros(n)
        if a.GetNumberOfComponents()==1:
            return np.array([a.GetTuple1(i) for i in range(n)])
        return np.array([a.GetTuple(i) for i in range(n)])
    Ar=arr("potential re"); Ai=arr("potential im")
    # try alternative names
    if not Ar.any():
        for nm in g.GetPointData().GetArrayName(i) if hasattr(g.GetPointData(),'GetArrayName') else []:
            pass
    pts=np.array([g.GetPoint(i) for i in range(n)])
    return Ar,Ai,pts

def tri_cells(g,pts):
    n=g.GetNumberOfCells()
    ids=[]
    for c in range(n):
        cell=g.GetCell(c)
        if cell.GetNumberOfPoints()==3:
            p=[cell.GetPointId(k) for k in range(3)]
            ids.append(p)
    return np.array(ids)

def region_centroid(pts,ids):
    return pts[ids].mean(axis=1)

def integrate_region(g,Ar,Ai,pts,ids,mask):
    # mask: boolean over cells
    sub=ids[mask]
    if len(sub)==0: return 0.0, 0.0, 0
    p0=pts[sub[:,0]]; p1=pts[sub[:,1]]; p2=pts[sub[:,2]]
    area=0.5*np.abs((p1[:,0]-p0[:,0])*(p2[:,1]-p0[:,1])-(p2[:,0]-p0[:,0])*(p1[:,1]-p0[:,1]))
    A=(Ar[sub[:,0]]+Ar[sub[:,1]]+Ar[sub[:,2]])/3.0
    Ai_=(Ai[sub[:,0]]+Ai[sub[:,1]]+Ai[sub[:,2]])/3.0
    Jm2=(W*SIG)**2*(A**2+Ai_**2)          # |J|^2 = (w sigma |A|)^2   (copper, Js=0)
    P=0.5*RHO*Jm2
    return float((P*area).sum()), float(np.sqrt(Jm2).max()), len(sub)

def analyze(fn):
    g=load(fn); Ar,Ai,pts=fields(g); ids=tri_cells(g,pts)
    cent=pts[ids].mean(axis=1)
    out={}
    regions={"top":(cent[:,1]>=-4e-5)&(cent[:,1]<=1e-6)&(np.abs(cent[:,0])<0.0062),
             "bot":(cent[:,1]>=-1.036e-3)&(cent[:,1]<=-0.999e-3)&(np.abs(cent[:,0])<0.0062)}
    for k,m in regions.items():
        P,Jmax,n=integrate_region(g,Ar,Ai,pts,ids,m)
        out[k]=(P,Jmax,n)
    return out, g,Ar,Ai,pts,ids

if __name__=="__main__":
    import glob
    fn=sys.argv[1]
    o,g,Ar,Ai,pts,ids=analyze(fn)
    print(fn)
    tot=0
    for k,(P,Jm,n) in o.items():
        print("  %s: P=%.6e W/m  |J|max=%.3e A/m2  ncells=%d"%(k,P,Jm,n)); tot+=P
    print("  TOTAL P=%.6e W/m"%tot)

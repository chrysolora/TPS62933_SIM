#!/usr/bin/env python3
"""Post-process axisymmetric Elmer harmonic Az -> P_eddy = int 0.5*rho*|J|^2 dA ; J=-i w sigma A."""
import numpy as np, math, sys
import vtk
SIG=5.8e7; RHO=1.0/SIG
def load(fn):
    r=vtk.vtkXMLUnstructuredGridReader(); r.SetFileName(fn); r.Update(); return r.GetOutput()
def analyze(fn,W):
    g=load(fn); n=g.GetNumberOfPoints()
    pd=g.GetPointData(); names=[pd.GetArrayName(i) for i in range(pd.GetNumberOfArrays())]
    def find(sub):
        for nm in names:
            if nm and sub in nm.lower(): 
                a=pd.GetArray(nm)
                if a.GetNumberOfComponents()==1: return np.array([a.GetTuple1(i) for i in range(n)])
        return None
    Ar=find("re"); Ai=find("im")
    if Ar is None or Ai is None:
        raise SystemExit("arrays not found: %s"%names)
    pts=np.array([g.GetPoint(i) for i in range(n)])
    cells=[]; 
    for c in range(g.GetNumberOfCells()):
        cell=g.GetCell(c)
        if cell.GetNumberOfPoints()==3:
            cells.append([cell.GetPointId(k) for k in range(3)])
    ids=np.array(cells); cent=pts[ids].mean(axis=1)
    TCU=3.5e-5; TB=1.0e-3
    def integ(mask):
        sub=ids[mask]
        if len(sub)==0: return 0.0,0.0,0
        p0,p1,p2=pts[sub[:,0]],pts[sub[:,1]],pts[sub[:,2]]
        area=0.5*np.abs((p1[:,0]-p0[:,0])*(p2[:,1]-p0[:,1])-(p2[:,0]-p0[:,0])*(p1[:,1]-p0[:,1]))
        A=(Ar[sub[:,0]]+Ar[sub[:,1]]+Ar[sub[:,2]])/3.0
        Ai_=(Ai[sub[:,0]]+Ai[sub[:,1]]+Ai[sub[:,2]])/3.0
        J2=(W*SIG)**2*(A**2+Ai_**2)
        P=0.5*RHO*J2
        return float((P*area).sum()), float(np.sqrt(J2).max()), len(sub)
    mtop=(cent[:,1]>=-1.2e-4)&(cent[:,1]<=1e-6)&(cent[:,0]<0.012)
    mbot=(cent[:,1]>=-1.05e-3)&(cent[:,1]<=-0.95e-3)&(cent[:,0]<0.012)
    Pt,Jt,nt=integ(mtop); Pb,Jb,nb=integ(mbot)
    return dict(P_top=Pt,P_bot=Pb,P_tot=Pt+Pb,Jmax_top=Jt,Jmax_bot=Jb,n_top=nt,n_bot=nb,names=names)
if __name__=="__main__":
    fn=sys.argv[1]; W=float(sys.argv[2]); print(analyze(fn,W))

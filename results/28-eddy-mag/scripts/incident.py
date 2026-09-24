#!/usr/bin/env python3
"""Run axisymmetric magnetostatics (Cu sigma=0) -> incident Bz(r,z); integrate thin-sheet eddy loss."""
import numpy as np, math, os, sys, subprocess
import vtk
SIG=5.8e7; RHO=1/SIG; MU0=4*math.pi*1e-7; TCU=3.5e-5; TB=1.0e-3
def load(fn):
    r=vtk.vtkXMLUnstructuredGridReader(); r.SetFileName(fn); r.Update(); return r.GetOutput()
def get_A(g):
    n=g.GetNumberOfPoints(); pd=g.GetPointData()
    names=[pd.GetArrayName(i) for i in range(pd.GetNumberOfArrays())]
    def f(s):
        for nm in names:
            if nm and s in nm.lower():
                a=pd.GetArray(nm)
                if a.GetNumberOfComponents()==1: return np.array([a.GetTuple1(i) for i in range(n)])
        return None
    return f("re"), f("im"), np.array([g.GetPoint(i) for i in range(n)])
def Bz_profile(fn, zplane):
    """incident Bz(r) at horizontal plane z=zplane, axisymmetric: Bz=(1/r) d(rA)/dr"""
    g=load(fn); Ar,Ai,pts=get_A(g)
    sel=np.abs(pts[:,1]-zplane)<8e-6
    P=pts[sel]; A=Ar[sel]
    o=np.argsort(P[:,0]); r=P[o,0]; a=A[o]
    # collapse duplicates by averaging on rounded r
    key=np.round(r,6); ur={}; 
    for k,ri,ai in zip(key,r,a): ur.setdefault(k,[]).append(ai)
    r=np.array(sorted(ur)); a=np.array([np.mean(ur[k]) for k in r])
    Bz=np.zeros_like(r)
    for i in range(len(r)):
        if r[i]<1e-7: Bz[i]=2*( (r[i+1]*a[i+1]-0)/(r[i+1]**2) ) if len(r)>1 else 0
        else:
            j=min(i+1,len(r)-1); k=max(i-1,0)
            Bz[i]=( (r[j]*a[j]-r[k]*a[k])/(r[j]-r[k]) )/r[i]
    return r,Bz
def integrate_thinsheet(r,Bz,rmax,rmin=0.0,W=2*math.pi*784e3,A_cu_mode="area"):
    # P = SIG*W^2*TCU^3/24 * int Bz^2 dA ; dA=2*pi*r*dr over r in [rmin,rmax]
    m=(r>=rmin)&(r<=rmax); rr=r[m]; b=Bz[m]**2
    if len(rr)<2: return 0.0,0.0
    # trapezoid integral of b * 2 pi r dr
    integ=np.trapz(b*2*math.pi*rr, rr)
    Pfac=SIG*W**2*TCU**3/24.0
    return Pfac*integ, Pfac
if __name__=="__main__":
    base=sys.argv[1]
    fn=os.path.join(base,"mesh","case_t0001.vtu")
    for zlab,zp in [("top",0.0),("bot",-TB)]:
        r,Bz=Bz_profile(fn,zp)
        P,_=integrate_thinsheet(r,Bz,0.030,0.0)
        print(zlab,"incident Bz max=%.4g T  Bz@1mm=%.4g  P_full_disc(r<6mm)=%.4g W/m"%(
            np.nanmax(np.abs(Bz)), np.interp(1e-3,r,Bz), integrate_thinsheet(r,Bz,0.006,0.0)[0]))

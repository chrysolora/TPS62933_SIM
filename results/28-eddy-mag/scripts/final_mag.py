#!/usr/bin/env python3
"""Final: 3-version eddy loss via (A) coupled harmonic [artifact ref] and (B) uncoupled thin-sheet
with REAL polygon geometry. Produces numbers_eddy_mag.json + figures."""
import numpy as np, math, os, json, glob
import vtk
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT="/mnt/raid10/sim-work/tps62933/eddy_mag"
TEST=os.path.join(ROOT,"test")
SIG=5.8e7; RHO=1/SIG; MU0=4*math.pi*1e-7; TCU=3.5e-5; TB=1.0e-3
F=784e3; W=2*math.pi*F
MIL=25.4e-6
def load(fn):
    r=vtk.vtkXMLUnstructuredGridReader(); r.SetFileName(fn); r.Update(); return r.GetOutput()
def getA(g):
    n=g.GetNumberOfPoints(); pd=g.GetPointData()
    names=[pd.GetArrayName(i) for i in range(pd.GetNumberOfArrays())]
    def f(s):
        for nm in names:
            if nm and s in nm.lower():
                a=pd.GetArray(nm)
                if a.GetNumberOfComponents()==1: return np.array([a.GetTuple1(i) for i in range(n)])
        return None
    return f("re"),f("im"),np.array([g.GetPoint(i) for i in range(n)])
def coupled_P(fn,rmax=0.012):
    g=load(fn); Ar,Ai,pts=getA(g)
    cells=[];
    for c in range(g.GetNumberOfCells()):
        cell=g.GetCell(c)
        if cell.GetNumberOfPoints()==3: cells.append([cell.GetPointId(k) for k in range(3)])
    ids=np.array(cells); cent=pts[ids].mean(axis=1)
    def integ(mask):
        sub=ids[mask]
        if len(sub)==0: return 0.,0.
        p0,p1,p2=pts[sub[:,0]],pts[sub[:,1]],pts[sub[:,2]]
        area=0.5*np.abs((p1[:,0]-p0[:,0])*(p2[:,1]-p0[:,1])-(p2[:,0]-p0[:,0])*(p1[:,1]-p0[:,1]))
        A=(Ar[sub[:,0]]+Ar[sub[:,1]]+Ar[sub[:,2]])/3; Ai_=(Ai[sub[:,0]]+Ai[sub[:,1]]+Ai[sub[:,2]])/3
        J2=(W*SIG)**2*(A**2+Ai_**2); P=0.5*RHO*J2
        return float((P*area).sum()), float(np.sqrt(J2).max())
    mt=(cent[:,1]>=-1.2e-4)&(cent[:,1]<=1e-6)&(cent[:,0]<rmax)
    mb=(cent[:,1]>=-1.05e-3)&(cent[:,1]<=-0.95e-3)&(cent[:,0]<rmax)
    Pt,Jt=integ(mt); Pb,Jb=integ(mb)
    return dict(P_top=Pt,P_bot=Pb,P_tot=Pt+Pb,Jmax_top=Jt,Jmax_bot=Jb)
def Bz_profile(fn,zp):
    g=load(fn); Ar,Ai,pts=getA(g)
    sel=np.abs(pts[:,1]-zp)<8e-6; P=pts[sel]; A=Ar[sel]
    o=np.argsort(P[:,0]); r=P[o,0]; a=A[o]
    key=np.round(r,6); ur={}
    for k,ri,ai in zip(key,r,a): ur.setdefault(k,[]).append(ai)
    r=np.array(sorted(ur)); a=np.array([np.mean(ur[k]) for k in r])
    Bz=np.zeros_like(r)
    for i in range(len(r)):
        if r[i]<1e-7: Bz[i]=( (r[1]*a[1]-0)/(r[1]**2) ) if len(r)>1 else 0
        else:
            j=min(i+1,len(r)-1); k=max(i-1,0)
            Bz[i]=((r[j]*a[j]-r[k]*a[k])/(r[j]-r[k]))/r[i]
    return r,Bz
def thinsheet(r,Bz,rmin,rmax):
    m=(r>=rmin)&(r<=rmax); rr=r[m]; b=Bz[m]**2
    if len(rr)<2: return 0.,0.
    integ=np.trapz(b*2*math.pi*rr,rr)
    fac=SIG*W**2*TCU**3/24.0
    return fac*integ, fac
# real polygon data (from extract_cu.py on epru)
board=1043.3071*1988.189*MIL**2
A_top_full=0.805*board; A_bot_full=0.756*board
keep_top=47190.88*MIL**2; keep_dual=53684.37*MIL**2
rcut_top=math.sqrt(keep_top/math.pi); rcut_dual=math.sqrt(keep_dual/math.pi)
Rdisc=0.013
# incident field from sigma=0 run
r_t,Bz_t=Bz_profile(os.path.join(TEST,"inc/mesh/case_t0001.vtu"),0.0)
r_b,Bz_b=Bz_profile(os.path.join(TEST,"inc/mesh/case_t0001.vtu"),-TB)
res={}
for name,rcut_t,rcut_b in [("FullCopper",0.,0.),("TopCutout",rcut_top,0.),("DualCutout",rcut_dual,rcut_dual)]:
    Pt,_=thinsheet(r_t,Bz_t,rcut_t,Rdisc); Pb,_=thinsheet(r_b,Bz_b,rcut_b,Rdisc)
    res[name]=dict(r_cut_top=rcut_t,r_cut_bot=rcut_b,P_top=Pt,P_bot=Pb,P_tot=Pt+Pb)
# coupled refs (cartesian artifact dirs not present; use axi runs fug/tut)
coupled={}
for tag,nm in [("fug","FullCopper"),("tut","TopCutout")]:
    v=os.path.join(TEST,tag,"mesh","case_t0001.vtu")
    if os.path.exists(v): coupled[nm]=coupled_P(v)
# J distribution (thin-sheet peak J=sigma*omega*Bz*t/2) for figure
Js_t=SIG*W*np.abs(Bz_t)*TCU/2; Js_b=SIG*W*np.abs(Bz_b)*TCU/2
out=dict(freq_Hz=F, sigma=SIG, rho=RHO, tcu_m=TCU,
  real_geom=dict(board_area_m2=board,A_top_full_m2=A_top_full,A_bot_full_m2=A_bot_full,
     keep_top_m2=keep_top,keep_dual_m2=keep_dual,rcut_top_m=rcut_top,rcut_dual_m=rcut_dual,
     fea4_bbox_keep_m2=(374.921-115.079)*(797.795-522.205)*MIL**2,
     fea4_bbox_note="fea4 used bbox 6.60x7.00mm=%.4g m2; REAL polygon area=%.4g m2 (%.0f%%)"%(  (374.921-115.079)*(797.795-522.205)*MIL**2, keep_top, 100*keep_top/((374.921-115.079)*(797.795-522.205)*MIL**2))),
  incident=dict(Bz_max_top=float(np.nanmax(np.abs(Bz_t))),Bz_max_bot=float(np.nanmax(np.abs(Bz_b))),
     Bz_at_1mm_top=float(np.interp(1e-3,r_t,Bz_t)),Bz_at_1mm_bot=float(np.interp(1e-3,r_b,Bz_b))),
  uncoupled_thinsheet=res, coupled_harmonic_ref=coupled)
json.dump(out,open(os.path.join(ROOT,"numbers_eddy_mag.json"),"w"),indent=2)

# ---- FIG 3 versions ----
fig,ax=plt.subplots(1,2,figsize=(12,5))
# (left) J_z(r) incident, with keepout markers
ax[0].plot(r_t*1e3,Js_t,label="Top Cu (incident)",lw=2)
ax[0].plot(r_b*1e3,Js_b,label="Bot Cu (incident)",lw=2)
ax[0].axvline(rcut_top*1e3,color="r",ls="--",label="TopCut r_cut=%.2fmm"% (rcut_top*1e3))
ax[0].axvline(rcut_dual*1e3,color="m",ls=":",label="DualCut r_cut=%.2fmm"% (rcut_dual*1e3))
ax[0].set_xlim(0,12); ax[0].set_xlabel("radius r (mm)"); ax[0].set_ylabel("|J| peak (A/m^2)")
ax[0].set_title("Incident eddy current density |J|=sigma*omega*Bz*t/2"); ax[0].legend(fontsize=8)
# (right) P_eddy comparison
names=["FullCopper","TopCutout","DualCutout"]; pt=[res[n]["P_top"]*1e3 for n in names]; pb=[res[n]["P_bot"]*1e3 for n in names]
x=np.arange(3); ax[1].bar(x-0.2,pt,0.4,label="Top Cu"); ax[1].bar(x+0.2,pb,0.4,label="Bot Cu")
ax[1].set_xticks(x); ax[1].set_xticklabels(names); ax[1].set_ylabel("P_eddy (mW/m)")
ax[1].set_title("Eddy loss per unit depth (uncoupled thin-sheet, REAL polygons)")
for i,n in enumerate(names): ax[1].text(i,max(pt[i],pb[i]),"%.4g"%(res[n]["P_tot"]*1e3),ha="center",va="bottom",fontsize=9)
ax[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(ROOT,"fig_eddy_mag_3v.png"),dpi=130)
print(json.dumps(out,indent=1)[:2000])
print("P_top(mW/m):",[(n,res[n]["P_top"]*1e3) for n in names])
print("P_tot(mW/m):",[(n,res[n]["P_tot"]*1e3) for n in names])
if coupled: print("COUPLED(artifact) W/m:",{k:v["P_tot"] for k,v in coupled.items()})

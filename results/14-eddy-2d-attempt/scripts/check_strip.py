#!/usr/bin/env python3
"""Analytic self-check: thin Cu plate 6mm x 35um in uniform AC tangential field B0=1mT.
Replicates loss_redo check2d. Expect ratio P_elmer/P_analytic ~0.8."""
import os,subprocess,math,re,glob,sys
import numpy as np
import post_eddy
W=post_eddy.W; SIG=post_eddy.SIG; RHO=post_eddy.RHO
d=os.path.join(post_eddy.__file__.rsplit("/",1)[0],"check_strip")
os.makedirs(d,exist_ok=True)
w=0.006; t=3.5e-5; L=0.02; B0=1e-3
geo='''Mesh.MshFileVersion=4.1; Mesh.SaveAll=1;
SetFactory("OpenCASCADE");
L=%g; w=%g; t=%g;
Rectangle(1)={-L,-L,0,2*L,2*L};
Rectangle(2)={-w/2,-t/2,0,w,t};
air[]=BooleanDifference{Surface{1};Delete;}{Surface{2};};
Physical Surface("Air")={air[]};
Physical Surface("Cu")={2};
Field[1]=Box; Field[1].VIn=2e-5; Field[1].VOut=4e-3;
Field[1].XMin=-w/2-2*t; Field[1].XMax=w/2+2*t; Field[1].YMin=-2*t; Field[1].YMax=2*t;
Background Field=1;
'''%(L,w,t)
open(os.path.join(d,"build.geo"),"w").write(geo)
subprocess.run(["gmsh","-2",os.path.join(d,"build.geo"),"-o",os.path.join(d,"mesh.msh")],cwd=d,capture_output=True)
subprocess.run(["ElmerGrid","14","2","mesh.msh","-out","mesh"],cwd=d,capture_output=True)
# find outer bcs by coords
nodes={}
for ln in open(os.path.join(d,"mesh","mesh.nodes")):
    p=ln.split()
    if len(p)>=5: nodes[int(p[0])]=(float(p[2]),float(p[3]))
o=set()
for ln in open(os.path.join(d,"mesh","mesh.boundary")):
    p=ln.split()
    if len(p)<7: continue
    a=nodes.get(int(p[-2])); b=nodes.get(int(p[-1]))
    if a and b and all((abs(abs(q[0])-L)<4e-3 or abs(abs(q[1])-L)<4e-3) for q in (a,b)): o.add(int(p[0]))
sif='''Header
  Mesh DB "." "mesh"
End
Simulation
  Coordinate System = Cartesian 2D
  Simulation Type = Steady State
  Steady State Max Iterations = 1
  Angular Frequency = %f
End
Body 1
  Target Bodies(1) = 1
  Equation = 1
  Material = 2
End
Body 2
  Target Bodies(1) = 2
  Equation = 1
  Material = 1
End
Equation 1
  Active Solvers(1) = 1
End
Solver 1
  Equation = "MagHarm"
  Procedure = "MagnetoDynamics2D" "MagnetoDynamics2DHarmonic"
  Linear System Solver = Direct
  Linear System Direct Method = UMFPACK
End
Solver 2
  Exec Solver = After All
  Procedure = "ResultOutputSolve" "ResultOutputSolver"
  Output File Name = case
  Output Format = vtu
  Binary Output = Logical False
End
Material 1
  Electric Conductivity = 5.8e7
  Relative Permeability = 1.0
End
Material 2
  Electric Conductivity = 0.0
  Relative Permeability = 1.0
End
Boundary Condition 1
  Target Boundaries(%d) = %s
  Potential = Variable Coordinate 2
    Real MATC "%g*tx"
End
'''%(W,len(o)," ".join(map(str,sorted(o))),B0)
open(os.path.join(d,"case.sif"),"w").write(sif)
r=subprocess.run(["ElmerSolver","case.sif"],cwd=d,capture_output=True,text=True,timeout=600)
open(os.path.join(d,"run.log"),"w").write(r.stdout+r.stderr)
vt=sorted(glob.glob(os.path.join(d,"mesh","case_t*.vtu")))
g=post_eddy.load(vt[0]); Ar,Ai,pts=post_eddy.fields(g); ids=post_eddy.tri_cells(g,pts)
cent=pts[ids].mean(axis=1)
mask=np.abs(cent[:,1])<t
p0=pts[ids[:,0]][mask]; p1=pts[ids[:,1]][mask]; p2=pts[ids[:,2]][mask]
area=0.5*np.abs((p1[:,0]-p0[:,0])*(p2[:,1]-p0[:,1])-(p2[:,0]-p0[:,0])*(p1[:,1]-p0[:,1]))
sub=ids[mask]
A=(Ar[sub[:,0]]+Ar[sub[:,1]]+Ar[sub[:,2]])/3; Ai_=(Ai[sub[:,0]]+Ai[sub[:,1]]+Ai[sub[:,2]])/3
Jm2=(W*SIG)**2*(A**2+Ai_**2)
P_elmer=float((0.5*RHO*Jm2*area).sum())
P_an=SIG*W**2*B0**2*w*t**3/24.0
print("SELF-CHECK strip: P_elmer=%.6e W/m  P_analytic=%.6e W/m  ratio=%.3f"%(P_elmer,P_an,P_elmer/P_an))
print("  Jmax=%.3e A/m2"%(np.sqrt(Jm2).max()))
open(os.path.join(d,"RESULT.txt"),"w").write("P_elmer=%.6e\nP_analytic=%.6e\nratio=%.4f\nJmax=%.4e\n"%(P_elmer,P_an,P_elmer/P_an,np.sqrt(Jm2).max()))

#!/usr/bin/env python3
"""Axisymmetric 2D harmonic eddy model (r,z) for TPS62933 L2 pour eddy loss.
Tests whether Axi Symmetric fixes the non-monotonic Cartesian-2D artifact.
"""
import os,sys,subprocess,math,glob,re,json
import numpy as np
SIG=5.8e7; RHO=1.0/SIG
MU0=4*math.pi*1e-7
ROOT=os.path.dirname(os.path.abspath(__file__))

def geo(Rtop,Rbot,rcut_top,rcut_bot,scale=1.0,
        HWc=6.0e-3,H_air=0.030,Z_air=0.02,TCU=3.5e-5,TB=1.0e-3,
        r_core=3.25e-3,z_core0=0.2e-3,z_core1=3.2e-3,
        r_coil0=2.0e-3,r_coil1=2.6e-3,z_coil0=0.9e-3,z_coil1=1.5e-3):
    L=[]
    L.append('Mesh.MshFileVersion=4.1; Mesh.SaveAll=1;')
    L.append('SetFactory("OpenCASCADE");')
    L.append('Rectangle(1)={0,%g,0,%g,%g};'%(-0.008,H_air,Z_air))   # air (r 0..0.03, z -0.008..0.012)
    # core
    L.append('Rectangle(11)={0,%g,0,%g,%g};'%(z_core0,r_core,z_core1-z_core0))
    L.append('Rectangle(12)={%g,%g,0,%g,%g};'%(r_coil0,z_coil0,r_coil1-r_coil0,z_coil1-z_coil0))  # coil window
    ct=[];cb=[]
    ct=['20']; L.append('Rectangle(20)={%g,%g,0,%g,%g};'%(rcut_top,-TCU,Rtop-rcut_top,TCU))
    cb=['21']
    zt=-TB-TCU
    L.append('Rectangle(21)={%g,%g,0,%g,%g};'%(rcut_bot,zt,Rbot-rcut_bot,TCU))
    L.append('cored[]=BooleanDifference{Surface{11};Delete;}{Surface{12};};')
    L.append('air[]=BooleanDifference{Surface{1};Delete;}{Surface{cored[],12,%s};};'%(','.join(ct+cb)))
    L.append('Physical Surface(1)={air[]};')
    L.append('Physical Surface(2)={cored[]};')
    L.append('Physical Surface(3)={12};')
    L.append('Physical Surface(4)={%s};'%(','.join(ct)))
    L.append('Physical Surface(5)={%s};'%(','.join(cb)))
    lc_cu=2.0e-5*scale; lc_co=1.5e-4*scale; lc_core=2.0e-4*scale
    L.append('Field[1]=Box; Field[1].VIn=%g; Field[1].VOut=2e-3;'%lc_co)
    L.append('Field[1].XMin=0; Field[1].XMax=0.008; Field[1].YMin=-1.2e-3; Field[1].YMax=3.4e-3;')
    L.append('Field[2]=Box; Field[2].VIn=%g; Field[2].VOut=2e-3;'%lc_core)
    L.append('Field[2].XMin=0; Field[2].XMax=0.004; Field[2].YMin=0; Field[2].YMax=3.4e-3;')
    L.append('Field[3]=Box; Field[3].VIn=%g; Field[3].VOut=1e-3;'%lc_cu)
    L.append('Field[3].XMin=0; Field[3].XMax=0.008; Field[3].YMin=-1.1e-3; Field[3].YMax=2*%g;'%TCU)
    L.append('Field[4]=Min; Field[4].FieldsList={1,2,3}; Background Field=4;')
    return "\n".join(L)+"\n"

def get_names(d):
    txt=open(os.path.join(d,"mesh.msh"),encoding="utf-8",errors="ignore").read()
    m=re.search(r"\$PhysicalNames\n(.*?)\$EndPhysicalNames",txt,re.S)
    b={}
    for ln in (m.group(1).splitlines() if m else []):
        mm=re.match(r'\s*(\d+)\s+(\d+)\s+"(.*)"',ln)
        if mm: b[mm.group(3)]=int(mm.group(2))
    return b

def outer_bcs(d,tol=3e-4):
    nodes={}
    for ln in open(os.path.join(d,"mesh","mesh.nodes")):
        p=ln.split()
        if len(p)>=5:
            try: nodes[int(p[0])]=(float(p[2]),float(p[3]))
            except: pass
    o=set()
    for ln in open(os.path.join(d,"mesh","mesh.boundary")):
        p=ln.split()
        if len(p)<7: continue
        try:
            bc=int(p[0]); a=nodes[int(p[-2])]; b=nodes[int(p[-1])]
        except: continue
        if all((abs(q[0]-0.03)<tol or abs(q[1]+0.008)<tol or abs(q[1]-0.012)<tol) for q in (a,b)):
            o.add(bc)
    return sorted(o)

def write_sif(d,b,outer,J0,W):
    s='Header\n  Mesh DB "." "mesh"\nEnd\n'
    s+="Simulation\n  Coordinate System = Axi Symmetric\n  Simulation Type = Steady State\n  Steady State Max Iterations = 1\n  Angular Frequency = %f\nEnd\n"%W
    def bd(n,body,mat,bf=None):
        t="Body %d\n  Target Bodies(1) = %d\n  Equation = 1\n  Material = %d\n"%(n,body,mat)
        if bf: t+="  Body Force = %d\n"%bf
        return t+"End\n"
    # Elmer renumbers bodies 1..5 in gmsh-tag order: Air=1,Core=2,Coil=3,CuTop=4,CuBot=5
    s+=bd(1,1,2); s+=bd(2,2,3); s+=bd(3,3,4,4)
    s+=bd(4,4,1); s+=bd(5,5,1)
    s+="Equation 1\n  Active Solvers(1) = 1\nEnd\n"
    s+='''Solver 1
  Equation = "MagHarm"
  Procedure = "MagnetoDynamics2D" "MagnetoDynamics2DHarmonic"
  Linear System Solver = Direct
  Linear System Direct Method = UMFPACK
  Nonlinear System Max Iterations = 1
End
Solver 2
  Exec Solver = After All
  Procedure = "ResultOutputSolve" "ResultOutputSolver"
  Output File Name = case
  Output Format = vtu
  Binary Output = Logical False
End
'''
    s+="Material 1\n  Name=\"Cu\"\n  Electric Conductivity = 5.8e7\n  Relative Permeability = 1.0\nEnd\n"
    s+="Material 2\n  Name=\"Air\"\n  Electric Conductivity = 0.0\n  Relative Permeability = 1.0\nEnd\n"
    s+="Material 3\n  Name=\"Core\"\n  Electric Conductivity = 0.0\n  Relative Permeability = 40.0\nEnd\n"
    s+="Material 4\n  Name=\"Coil\"\n  Electric Conductivity = 0.0\n  Relative Permeability = 1.0\nEnd\n"
    s+="Body Force 4\n  Current Density = %r\nEnd\n"%(J0)
    s+="Boundary Condition 1\n  Target Boundaries(%d) = %s\n  Potential = Real 0.0\nEnd\n"%(len(outer)," ".join(map(str,outer)))
    open(os.path.join(d,"case.sif"),"w").write(s)

def build(tag,Rtop,Rbot,rcut_top,rcut_bot,scale=1.0,J0=1.667e7,W=2*math.pi*784e3,run=True):
    d=os.path.join(ROOT,tag); os.makedirs(d,exist_ok=True)
    open(os.path.join(d,"build.geo"),"w").write(geo(Rtop,Rbot,rcut_top,rcut_bot,scale))
    r=subprocess.run(["gmsh","-2",os.path.join(d,"build.geo"),"-o",os.path.join(d,"mesh.msh")],cwd=d,capture_output=True,text=True)
    if r.returncode!=0: open(os.path.join(d,"gmsh.err"),"w").write(r.stdout+r.stderr); print("GMSH FAIL",tag)
    subprocess.run(["ElmerGrid","14","2","mesh.msh","-out","mesh"],cwd=d,capture_output=True)
    b=get_names(d); outer=outer_bcs(d)
    write_sif(d,b,outer,J0,W)
    if run:
        r=subprocess.run(["ElmerSolver","case.sif"],cwd=d,capture_output=True,text=True,timeout=900)
        open(os.path.join(d,"run.log"),"w").write(r.stdout+r.stderr)
    return d

if __name__=="__main__":
    tag=sys.argv[1]; Rt=float(sys.argv[2]); Rc=float(sys.argv[3]); rct=float(sys.argv[4]); rcb=float(sys.argv[5])
    scale=float(sys.argv[6]) if len(sys.argv)>6 else 1.0
    d=build(tag,Rt,Rt,rct,rcb,scale)
    log=open(os.path.join(d,"run.log")).read()
    print(tag,"done, P?")

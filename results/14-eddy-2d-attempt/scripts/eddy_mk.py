#!/usr/bin/env python3
# 2D planar Az harmonic eddy model for TPS62933 inductor copper-pour loss.
import os,sys,subprocess,math,glob,re
import numpy as np
SIG=5.8e7; RHO=1/SIG; W=2*math.pi*805e3
MU0=4*math.pi*1e-7
root=os.path.dirname(os.path.abspath(__file__))

def geo(xct,xcb,scale=1.0,HWc=0.006,TCU=3.5e-5,TB=1.0e-3):
    # copper strips limited to +/-HWc (6 mm) around inductor; fine local mesh required.
    L=[]
    L.append('Mesh.MshFileVersion=4.1; Mesh.SaveAll=1;')
    L.append('SetFactory("OpenCASCADE");')
    L.append('Rectangle(1)={-0.02,-0.008,0,0.04,0.018};')
    L.append('Rectangle(11)={-3.25e-3,0.2e-3,0,6.5e-3,3.0e-3};')    # core
    L.append('Rectangle(12)={-2.6e-3,0.9e-3,0,0.6e-3,0.6e-3};')      # coil L
    L.append('Rectangle(13)={2.0e-3,0.9e-3,0,0.6e-3,0.6e-3};')       # coil R
    ct=[];cb=[]
    if xct<=0:
        L.append('Rectangle(14)={-%g,-%g,0,%g,%g};'%(HWc,TCU,2*HWc,TCU)); ct=['14']
    else:
        L.append('Rectangle(14)={-%g,-%g,0,%g-%g,%g};'%(HWc,TCU,HWc,xct,TCU)); ct.append('14')
        L.append('Rectangle(15)={%g,-%g,0,%g-%g,%g};'%(xct,TCU,HWc,xct,TCU)); ct.append('15')
    if xcb<=0:
        L.append('Rectangle(16)={-%g,-%g-%g,0,%g,%g};'%(HWc,TB,TCU,2*HWc,TCU)); cb=['16']
    else:
        L.append('Rectangle(16)={-%g,-%g-%g,0,%g-%g,%g};'%(HWc,TB,TCU,HWc,xcb,TCU)); cb.append('16')
        L.append('Rectangle(17)={%g,-%g-%g,0,%g-%g,%g};'%(xcb,TB,TCU,HWc,xcb,TCU)); cb.append('17')
    tools=['12','13']+ct+cb
    L.append('cored[]=BooleanDifference{Surface{11};Delete;}{Surface{12,13};};')
    L.append('air[]=BooleanDifference{Surface{1};Delete;}{Surface{cored[],%s};};'%(','.join(tools)))
    L.append('Physical Surface("Air")={air[]};')
    L.append('Physical Surface("Core")={cored[]};')
    L.append('Physical Surface("CoilL")={12};')
    L.append('Physical Surface("CoilR")={13};')
    L.append('Physical Surface("CuTop")={%s};'%(','.join(ct)))
    L.append('Physical Surface("CuBot")={%s};'%(','.join(cb)))
    # mesh sizing: fine on copper (4e-6), graded around inductor, coarser far
    lc_cu=2.0e-5*scale; lc_co=1.2e-4*scale; lc_core=1.5e-4*scale
    L.append('Field[1]=Box; Field[1].VIn=%g; Field[1].VOut=2e-3;'%lc_co)
    L.append('Field[1].XMin=-0.0035; Field[1].XMax=0.0035; Field[1].YMin=-1.2e-3; Field[1].YMax=1.6e-3;')
    L.append('Field[2]=Box; Field[2].VIn=%g; Field[2].VOut=2e-3;'%lc_core)
    L.append('Field[2].XMin=-3.5e-3; Field[2].XMax=3.5e-3; Field[2].YMin=-1.1e-3; Field[2].YMax=3.4e-3;')
    L.append('Field[3]=Box; Field[3].VIn=%g; Field[3].VOut=1e-3;'%lc_cu)
    L.append('Field[3].XMin=%g; Field[3].XMax=%g; Field[3].YMin=%g; Field[3].YMax=%g;'%(-HWc-2e-4,HWc+2e-4,-TB-3*TCU,2*TCU))
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
        try: bc=int(p[0]); a=nodes[int(p[-2])]; b=nodes[int(p[-1])]
        except: continue
        if all((abs(abs(q[0])-0.02)<tol or abs(q[1]+0.008)<tol or abs(q[1]-0.01)<tol) for q in (a,b)):
            o.add(bc)
    return sorted(o)

def write_sif(d,b,outer,J0):
    s='Header\n  Mesh DB "." "mesh"\nEnd\n'
    s+="Simulation\n  Coordinate System = Cartesian 2D\n  Simulation Type = Steady State\n  Steady State Max Iterations = 1\n  Angular Frequency = %f\nEnd\n"%W
    def bd(n,body,mat,bf=None):
        t="Body %d\n  Target Bodies(1) = %d\n  Equation = 1\n  Material = %d\n"%(n,body,mat)
        if bf: t+="  Body Force = %d\n"%bf
        return t+"End\n"
    s+=bd(1,b["CuTop"],1); s+=bd(2,b["CuBot"],1); s+=bd(3,b["Core"],3)
    s+=bd(4,b["CoilL"],4,4); s+=bd(5,b["CoilR"],4,5); s+=bd(6,b["Air"],2)
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
    s+="Body Force 4\n  Current Density = %r\nEnd\n"%(-J0)
    s+="Body Force 5\n  Current Density = %r\nEnd\n"%(J0)
    s+="Boundary Condition 1\n  Target Boundaries(%d) = %s\n  Potential = Real 0.0\nEnd\n"%(len(outer)," ".join(map(str,outer)))
    open(os.path.join(d,"case.sif"),"w").write(s)

def build(tag,xct,xcb,scale=1.0,J0=1.6667e7,run=True):
    d=os.path.join(root,tag); os.makedirs(d,exist_ok=True)
    open(os.path.join(d,"build.geo"),"w").write(geo(xct,xcb,scale))
    r=subprocess.run(["gmsh","-2",os.path.join(d,"build.geo"),"-o",os.path.join(d,"mesh.msh")],cwd=d,capture_output=True,text=True)
    if r.returncode!=0: open(os.path.join(d,"gmsh.err"),"w").write(r.stdout+r.stderr)
    subprocess.run(["ElmerGrid","14","2","mesh.msh","-out","mesh"],cwd=d,capture_output=True)
    b=get_names(d); outer=outer_bcs(d)
    write_sif(d,b,outer,J0)
    if run:
        r=subprocess.run(["ElmerSolver","case.sif"],cwd=d,capture_output=True,text=True,timeout=1200)
        open(os.path.join(d,"run.log"),"w").write(r.stdout+r.stderr)
    return d

if __name__=="__main__":
    tag=sys.argv[1]; xct=float(sys.argv[2]); xcb=float(sys.argv[3])
    scale=float(sys.argv[4]) if len(sys.argv)>4 else 1.0
    d=build(tag,xct,xcb,scale)
    log=open(os.path.join(d,"run.log")).read()
    m=re.search(r"SS \(ITER=1\) \(NRM,RELC\): \( *([\d.E+-]+)",log)
    print(tag,"NRM=",m.group(1) if m else "?", "unused:", "yes" if "Unused keywords:" in log and log.split("Unused keywords:")[1].strip()[:20] else "none")

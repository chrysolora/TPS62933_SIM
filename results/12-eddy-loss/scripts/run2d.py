#!/usr/bin/env python3
# 2D planar (Az) harmonic eddy model of pour under L2. x=h-plane, y=vertical(z)
import os,sys,subprocess,math,numpy as np
SIG=5.8e7; RHO=1/SIG; FREQ=784e3; W=2*math.pi*FREQ
J0=6.0/(0.6e-3*0.6e-3)   # coil leg current density for 6A turn, A/m^2
TCU=3.5e-5; TB=1.0e-3; HW=0.015
root="/mnt/raid10/sim-work/tps62933/loss_redo/model2d"

def build_geo(xct,xcb):
    xc=xct
    return f'''SetFactory("OpenCASCADE");
Rectangle(1)={{-0.02,-0.008,0,0.04,0.018}};
Rectangle(2)={{-3.25e-3,0.2e-3,0,6.5e-3,3.0e-3}};
Rectangle(3)={{-2.6e-3,0.9e-3,0,0.6e-3,0.6e-3}};
Rectangle(4)={{2.0e-3,0.9e-3,0,0.6e-3,0.6e-3}};
Rectangle(5)={{-{HW},-{TCU},0,{HW}-{xct},{TCU}}};
Rectangle(6)={{{xct},-{TCU},0,{HW}-{xct},{TCU}}};
Rectangle(7)={{-{HW},-{TB}-{TCU},0,{HW}-{xcb},{TCU}}};
Rectangle(8)={{{xcb},-{TB}-{TCU},0,{HW}-{xcb},{TCU}}};
BooleanFragments{{Surface{{1,2,3,4,5,6,7,8}};Delete;}}{{}};
Physical Surface("CuTop")={{Surface In BoundingBox{{-0.01501,-3.6e-5,-1e-9,0.01501,1e-9,1e-9}}}};
Physical Surface("CuBot")={{Surface In BoundingBox{{-0.01501,-1.0351e-3,-1e-9,0.01501,-0.9999e-3,1e-9}}}};
Physical Surface("Core")={{Surface In BoundingBox{{-3.26e-3,0.199e-3,-1e-9,3.26e-3,3.201e-3,1e-9}}}};
Physical Surface("CoilL")={{Surface In BoundingBox{{-2.61e-3,0.899e-3,-1e-9,-1.999e-3,1.501e-3,1e-9}}}};
Physical Surface("CoilR")={{Surface In BoundingBox{{1.999e-3,0.899e-3,-1e-9,2.601e-3,1.501e-3,1e-9}}}};
Physical Surface("Air")={{Surface In BoundingBox{{-0.0201,-0.0081,-1e-9,0.0201,0.0101,1e-9}}}};
Field[1]=Box; Field[1].VIn=1.2e-4; Field[1].VOut=3e-3;
Field[1].XMin=-0.004; Field[1].XMax=0.004; Field[1].YMin=-1.2e-3; Field[1].YMax=1.2e-3;
Field[2]=Box; Field[2].VIn=2e-4; Field[2].VOut=3e-3;
Field[2].XMin=-0.016; Field[2].XMax=0.016; Field[2].YMin=-1.1e-3; Field[2].YMax=1e-3;
Field[3]=Min; Field[3].FieldsList={{1,2}}; Background Field=3;
Mesh.MshFileVersion=2.2; Mesh.SaveAll=1;
'''

def mesh_names(d):
    bodies={}; bcs={}
    cur=None
    for ln in open(os.path.join(d,"mesh","mesh.names")):
        p=ln.split()
        if not p: continue
        if p[0]=="$" and len(p)>=4 and p[2]=="=":
            bodies[p[1]]=int(p[3])
        if p[0]=="!" and "boundaries" in ln: cur="bc"
        if p[0]=="$" and len(p)>=4 and "line_bc" in p[1]:
            bcs[int(p[3])]=p[1]
    return bodies,bcs

def outer_bcs(d,names,extent=0.02,tol=2e-4):
    nodes={}
    for ln in open(os.path.join(d,"mesh","mesh.nodes")):
        p=ln.split()
        if len(p)>=5:
            try: nodes[int(p[0])]=(float(p[2]),float(p[3]))
            except: pass
    outer=set()
    for ln in open(os.path.join(d,"mesh","mesh.boundary")):
        p=ln.split()
        if len(p)<7: continue
        try: bc=int(p[1]); n1=int(p[-2]); n2=int(p[-1])
        except: continue
        a=nodes.get(n1); b=nodes.get(n2)
        if a and b and all(abs(q[0])>extent-tol or abs(q[1]-0.010)<tol or abs(q[1]+0.008)<tol for q in (a,b)):
            outer.add(bc)
    return sorted(outer)

def write_sif(d,bodies,outer):
    b=bodies
    def bodytxt(n,body,name,mat,bf=None):
        s=f"Body {n}\n  Target Bodies(1) = {body}\n  Name = \"{name}\"\n  Equation = 1\n  Material = {mat}\n"
        if bf: s+=f"  Body Force = {bf}\n"
        return s+"End\n"
    s='Header\n  Mesh DB "." "mesh"\nEnd\n'
    s+="Simulation\n  Coordinate System = Cartesian 2D\n  Simulation Type = Steady State\n  Steady State Max Iterations = 1\n  Angular Frequency = %f\nEnd\n"%W
    s+="Constants\n  Permeability of Vacuum = 1.25664e-6\nEnd\n"
    # bodies: CuTop=1,CuBot=2,Core=3,CoilL=4,CoilR=5,Air=6
    s+=bodytxt(1,b["CuTop"],"CuTop",1)
    s+=bodytxt(2,b["CuBot"],"CuBot",1)
    s+=bodytxt(3,b["Core"],"Core",3)
    s+=bodytxt(4,b["CoilL"],"CoilL",4,4)
    s+=bodytxt(5,b["CoilR"],"CoilR",4,5)
    s+=bodytxt(6,b["Air"],"Air",2)
    s+="Equation 1\n  Name = \"E\"\n  Active Solvers(1) = 1\nEnd\n"
    s+='''Solver 1
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
'''
    s+="Material 1\n  Name = \"Cu\"\n  Electric Conductivity = 5.8e7\n  Relative Permeability = 1.0\nEnd\n"
    s+="Material 2\n  Name = \"Air\"\n  Electric Conductivity = 0.0\n  Relative Permeability = 1.0\nEnd\n"
    s+="Material 3\n  Name = \"Core\"\n  Electric Conductivity = 0.0\n  Relative Permeability = 40.0\nEnd\n"
    s+="Material 4\n  Name = \"Coil\"\n  Electric Conductivity = 0.0\n  Relative Permeability = 1.0\nEnd\n"
    s+="Body Force 4\n  Name = \"srcL\"\n  Current Density = %r\nEnd\n"%(-J0)
    s+="Body Force 5\n  Name = \"srcR\"\n  Current Density = %r\nEnd\n"%(J0)
    s+="Boundary Condition 1\n  Target Boundaries(%d) = %s\n  Potential = Real 0.0\nEnd\n"%(len(outer)," ".join(map(str,outer)))
    open(os.path.join(d,"case.sif"),"w").write(s)

def run(xct,xcb,tag):
    d=os.path.join(root,tag); os.makedirs(d,exist_ok=True)
    open(os.path.join(d,"build.geo"),"w").write(build_geo(xct,xcb))
    subprocess.run(["gmsh","-2",os.path.join(d,"build.geo"),"-o",os.path.join(d,"mesh.msh")],cwd=d,capture_output=True)
    subprocess.run(["ElmerGrid","14","2","mesh.msh","-out","mesh"],cwd=d,capture_output=True)
    bodies,bcs=mesh_names(d); outer=outer_bcs(d,None)
    write_sif(d,bodies,outer)
    r=subprocess.run(["ElmerSolver","case.sif"],cwd=d,capture_output=True,text=True,timeout=600)
    open(os.path.join(d,"run.log"),"w").write(r.stdout+r.stderr)
    return d,bodies,outer

def loss(d):
    import vtk
    r=vtk.vtkXMLUnstructuredGridReader(); r.SetFileName(os.path.join(d,"mesh","case_t0001.vtu")); r.Update(); g=r.GetOutput()
    n=g.GetNumberOfPoints(); pts=np.array([g.GetPoint(i) for i in range(n)])
    pr=g.GetPointData().GetArray("potential re"); pi=g.GetPointData().GetArray("potential im")
    Ar=np.array([pr.GetTuple1(i) for i in range(n)]); Ai=np.array([pi.GetTuple1(i) for i in range(n)])
    out={}
    bands={"top":(-TCU-1e-6,1e-6),"bot":(-TB-TCU-1e-6,-TB+1e-6)}
    for band,(lo,hi) in bands.items():
        P=0.0;V=0.0
        for c in range(g.GetNumberOfCells()):
            if g.GetCellType(c)!=vtk.VTK_TRIANGLE: continue
            ids=g.GetCell(c).GetPointIds(); k=ids.GetNumberOfIds()
            p=np.array([pts[ids.GetId(q)] for q in range(k)]); cen=p.mean(axis=0)
            if lo<=cen[1]<=hi and abs(cen[0])<=0.01501:
                A=0.5*abs((p[1,0]-p[0,0])*(p[2,1]-p[0,1])-(p[2,0]-p[0,0])*(p[1,1]-p[0,1]))
                ar=np.mean(Ar[[ids.GetId(q) for q in range(k)]]); ai=np.mean(Ai[[ids.GetId(q) for q in range(k)]])
                Jm2=W*W*SIG*SIG*(ar*ar+ai*ai)
                P+=0.5*RHO*Jm2*A; V+=A
        out[band]=(P,V)
    return out

if __name__=="__main__":
    cases=[("fullcu",0.0,0.0),("topcut",0.0033,0.0),("dualcut",0.0033,0.0033)]
    for tag,xct,xcb in cases:
        d,bodies,outer=run(xct,xcb,tag)
        L=loss(d)
        print("=== %s (topcut=%.4f botcut=%.4f) ==="%(tag,xct,xcb), "bodies",bodies,"outer",outer)
        for band,(P,V) in L.items(): print("  %s P=%.5e W/m Vcu=%.4e"%(band,P,V))

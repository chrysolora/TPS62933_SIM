import vtk,numpy as np,math
SIGMA=5.8e7; RHO=1/SIGMA; W=2*math.pi*805e3
def load(fn):
    r=vtk.vtkXMLUnstructuredGridReader(); r.SetFileName(fn); r.Update(); return r.GetOutput()
def elems(g):
    n=g.GetNumberOfPoints(); pts=np.array([g.GetPoint(i) for i in range(n)])
    pr=g.GetPointData().GetArray("potential re"); pi=g.GetPointData().GetArray("potential im")
    Ar=np.array([pr.GetTuple1(i) for i in range(n)]); Ai=np.array([pi.GetTuple1(i) for i in range(n)])
    out=[]
    for c in range(g.GetNumberOfCells()):
        if g.GetCellType(c)!=vtk.VTK_TRIANGLE: continue
        ids=g.GetCell(c).GetPointIds(); k=ids.GetNumberOfIds()
        p=np.array([pts[ids.GetId(q)] for q in range(k)])
        A=0.5*abs((p[1,0]-p[0,0])*(p[2,1]-p[0,1])-(p[2,0]-p[0,0])*(p[1,1]-p[0,1]))
        ar=np.mean(Ar[[ids.GetId(q) for q in range(k)]]); ai=np.mean(Ai[[ids.GetId(q) for q in range(k)]])
        out.append((p.mean(axis=0),A,ar,ai))
    return out
def loss_in(fn, xmin,xmax,zmin,zmax, Js=0.0):
    g=load(fn); P=0.0; V=0.0
    for (cen,A,ar,ai) in elems(g):
        if xmin<=cen[0]<=xmax and zmin<=cen[1]<=zmax:
            Jr=Js - 0.0; # J = Js - i w sigma (Ar + i Ai) = Js + w*s*Ai  (real),  -i*w*s*Ar (imag neg)
            Jrr=Js + W*SIGMA*ai; Jii=-W*SIGMA*ar
            Jm2=Jrr*Jrr+Jii*Jii
            P+=0.5*RHO*Jm2*A; V+=A
    return P,V
if __name__=="__main__":
    import sys
    print(loss_in(*[sys.argv[1]]+[float(x) for x in sys.argv[2:]]))

import vtk, numpy as np, sys, os
FEA4="/mnt/raid10/sim-work/tps62933/fea4"
TB=0.001; TCU=0.000035
# copper z bands
BANDS={"top":(-TCU-2e-6, 2e-6), "bot":(-TB-TCU-2e-6, -TB+2e-6)}

def read(fn):
    r=vtk.vtkXMLUnstructuredGridReader(); r.SetFileName(fn); r.Update(); return r.GetOutput()

def tet_vol(p):
    a,b,c,d=p
    return abs(np.dot(np.cross(b-a,c-a),d-a))/6.0

def integrate(ver):
    g=read(os.path.join(FEA4,ver,"mesh","case_t0001.vtu"))
    npts=g.GetNumberOfPoints()
    pts=np.array([g.GetPoint(i) for i in range(npts)])
    cd=g.GetCellData()
    jh=cd.GetArray("joule heating e")
    jre=cd.GetArray("current density re e"); jim=cd.GetArray("current density im e")
    sigma=5.8e7
    res={}
    for band,(lo,hi) in BANDS.items():
        P_el=0.0; P_j2=0.0; Vcu=0.0; nc=0
        for c in range(g.GetNumberOfCells()):
            ids=g.GetCell(c).GetPointIds(); k=ids.GetNumberOfIds()
            p=[pts[ids.GetId(q)] for q in range(k)]
            zs=[x[2] for x in p]
            if min(zs)>=lo and max(zs)<=hi:
                if g.GetCellType(c)==vtk.VTK_TETRA:
                    v=tet_vol(p)
                elif g.GetCellType(c)==vtk.VTK_WEDGE:
                    # wedge volume
                    v=0.0
                    a,b,cc=p[0],p[1],p[2]; d,e,f=p[3],p[4],p[5]
                    v=abs(np.dot(np.cross(b-a,cc-a),d-a))/6.0+abs(np.dot(np.cross(e-d,f-d),a-d))/6.0+0
                    # approximate wedge as prism via 3 tets
                else:
                    v=0.0
                P_el+=jh.GetTuple1(c)*v
                jr=np.array(jre.GetTuple(c)); ji=np.array(jim.GetTuple(c))
                jm2=float(np.dot(jr,jr)+np.dot(ji,ji))
                P_j2+=0.5*(1.0/sigma)*jm2*v
                if jh.GetTuple1(c)>1e-6: Vcu+=v; nc+=1
        res[band]=(P_el,P_j2,Vcu,nc)
    return res

if __name__=="__main__":
    for ver in (sys.argv[1:] or ["fullcu","topcut","dualcut"]):
        r=integrate(ver)
        print("=== %s ==="%ver)
        for band,(P_el,P_j2,V,nc) in r.items():
            print("  %-4s  P_jouleElmer=%.4e W  P_fromJ=%.4e W  Vcu=%.4e m^3 cells=%d"%(band,P_el,P_j2,V,nc))

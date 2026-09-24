import vtk, numpy as np, sys, os
d = sys.argv[1]; TB = float(sys.argv[2]); tcu = float(sys.argv[3])
r = vtk.vtkXMLUnstructuredGridReader(); r.SetFileName(os.path.join(d, "mesh/case_t0001.vtu")); r.Update()
g = r.GetOutput()
n = g.GetNumberOfPoints()
pts = np.array([g.GetPoint(i) for i in range(n)])
nd = g.GetPointData()
Jre = np.array(nd.GetArray("current density re")); Jim = np.array(nd.GetArray("current density im"))
Jmag = np.sqrt(np.sum(Jre**2 + Jim**2, axis=1))

def slab(lo, hi, outname):
    cells = []; used = set()
    for c in range(g.GetNumberOfCells()):
        ids = g.GetCell(c).GetPointIds(); k = ids.GetNumberOfIds()
        zs = [pts[ids.GetId(q)][2] for q in range(k)]
        if min(zs) >= lo and max(zs) <= hi:
            cells.append((g.GetCellType(c), [ids.GetId(q) for q in range(k)]))
            used.update(ids.GetId(q) for q in range(k))
    remap = {o: i for i, o in enumerate(sorted(used))}
    ng = vtk.vtkUnstructuredGrid(); p = vtk.vtkPoints()
    for o in sorted(used):
        q = pts[o]; p.InsertNextPoint(q[0], q[1], q[2])
    ng.SetPoints(p)
    a = vtk.vtkDoubleArray(); a.SetName("Jmag"); a.SetNumberOfComponents(1); a.SetNumberOfTuples(len(remap))
    for o, i in remap.items(): a.SetValue(i, float(Jmag[o]))
    ng.GetPointData().AddArray(a)
    for ct, idl in cells:
        il = vtk.vtkIdList()
        for o in idl: il.InsertNextId(remap[o])
        ng.InsertNextCell(ct, il)
    w = vtk.vtkXMLUnstructuredGridWriter(); w.SetFileName(os.path.join(d, outname))
    w.SetDataModeToAscii(); w.SetInputData(ng); w.Write()
    print(outname, "cells", len(cells), "pts", len(used))

slab(-tcu - 1e-6, 1e-6, "mesh/copper_top.vtu")
slab(-TB - tcu - 1e-6, -TB + 1e-6, "mesh/copper_bot.vtu")

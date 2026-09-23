import vtk
import numpy as np

r = vtk.vtkXMLUnstructuredGridReader()
r.SetFileName("mesh/case_t0001.vtu")
r.Update()
g = r.GetOutput()
pts = np.array([g.GetPoint(i) for i in range(g.GetNumberOfPoints())])

ndata = g.GetPointData()
Jre = np.array(ndata.GetArray("current density re"))
Jim = np.array(ndata.GetArray("current density im"))
Jmag = np.sqrt(np.sum(Jre**2 + Jim**2, axis=1))

# select copper cells (z strictly within the 35um slab)
cells = []
used = set()
for c in range(g.GetNumberOfCells()):
    ids = g.GetCell(c).GetPointIds()
    n = ids.GetNumberOfIds()
    zs = [pts[ids.GetId(k)][2] for k in range(n)]
    if min(zs) >= -36e-6 and max(zs) <= 1e-6:
        cells.append((g.GetCellType(c), [ids.GetId(k) for k in range(n)]))
        used.update(ids.GetId(k) for k in range(n))
print("copper cells:", len(cells), "points:", len(used))

remap = {old: i for i, old in enumerate(sorted(used))}
new = vtk.vtkUnstructuredGrid()
np_ = vtk.vtkPoints()
for old in sorted(used):
    p = pts[old]; np_.InsertNextPoint(p[0], p[1], p[2])
new.SetPoints(np_)
jm = vtk.vtkDoubleArray(); jm.SetName("Jmag"); jm.SetNumberOfComponents(1)
jm.SetNumberOfTuples(len(remap))
for old, i in remap.items():
    jm.SetValue(i, float(Jmag[old]))
new.GetPointData().AddArray(jm)

for ctype, idl in cells:
    ids = vtk.vtkIdList()
    for o in idl: ids.InsertNextId(remap[o])
    new.InsertNextCell(ctype, ids)

w = vtk.vtkXMLUnstructuredGridWriter()
w.SetFileName("mesh/copper.vtu")
w.SetDataModeToAscii()
w.SetInputData(new)
w.Write()
print("wrote mesh/copper.vtu")

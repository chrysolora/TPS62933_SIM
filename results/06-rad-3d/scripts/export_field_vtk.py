#!/usr/bin/env python3
"""Export the near-field |E| on the V slice to a VTK structured grid (.vts)
so ParaView can display it together with board.step."""
import os
import numpy as np
import h5py
import vtk
from vtk.util import numpy_support as ns

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.join(HERE, 'prop', 'sim')
frame = int(os.environ.get('FRAME', '-1'))

f = h5py.File(os.path.join(SIM, 'E_sliceV.h5'), 'r')
td = f['FieldData/TD']
ks = sorted(td.keys(), key=lambda s: int(s))
kk = ks[frame]
x = np.array(f['Mesh/x']); y = np.array(f['Mesh/y']); z = np.array(f['Mesh/z'])
a = np.array(td[kk])
e = np.sqrt((a * a).sum(0))[:, 0, :]           # (nx, nz)
vref = float(e.max()) or 1.0
dB = 20.0 * np.log10(np.maximum(e, vref*1e-6)/vref)

nx, nz = len(x), len(z)
sg = vtk.vtkStructuredGrid()
sg.SetDimensions(nx, 1, nz)
pts = vtk.vtkPoints()
pts.SetNumberOfPoints(nx * nz)
idx = 0
for iz in range(nz):
    for ix in range(nx):
        pts.SetPoint(idx, x[ix]*1e3, y[0]*1e3, z[iz]*1e3)
        idx += 1
sg.SetPoints(pts)
arr = ns.numpy_to_vtk(dB.reshape(-1), deep=True, array_type=vtk.VTK_FLOAT)
arr.SetName('E_dB')
sg.GetPointData().SetScalars(arr)
w = vtk.vtkStructuredGridWriter()
w.SetFileName(os.path.join(HERE, 'field_slice.vtk'))
w.SetInputData(sg)
w.SetFileTypeToBinary()
w.Write()
print('wrote field_slice.vtk frame', kk, 't=%.3f ns' % (float(td[kk].attrs['time'])*1e9))

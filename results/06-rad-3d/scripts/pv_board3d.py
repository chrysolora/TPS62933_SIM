#!/usr/bin/env python3
"""ParaView render: board.stl (reconstructed board model) + near-field |E|
slice -> fig_board3d.png.  Run headless:
   LIBGL_ALWAYS_SOFTWARE=1 xvfb-run -a -s '-screen 0 1200x900x24' pvbatch pv_board3d.py
"""
import os
from paraview.simple import *

os.chdir(os.path.dirname(os.path.abspath(__file__)))
view = CreateView('RenderView')
view.ViewSize = [1200, 900]
view.Background = [0.05, 0.05, 0.08]
view.OrientationAxesVisibility = 1

# --- board model (STEP) ---
board = OpenDataFile('board.stl')
bd = Show(board, view)
bd.Representation = 'Surface With Edges'
bd.DiffuseColor = [0.75, 0.55, 0.25]
bd.LineWidth = 1.0
try:
    bd.EdgeColor = [0.1, 0.1, 0.1]
except Exception:
    pass

# --- near-field |E| slice ---
field = OpenDataFile('field_slice.vtk')
fd = Show(field, view)
fd.Representation = 'Surface'
ColorBy(fd, ('POINTS', 'E_dB'))
lut = GetColorTransferFunction('E_dB')
lut.ApplyPreset('Turbo', True)
lut.RescaleTransferFunction(-90.0, 0.0)
lut.NumberOfTableValues = 128
fb = GetScalarBar(lut, view)
fb.Title = '|E|  [dB rel. max]'
fb.ComponentTitle = ''
fb.Visibility = 1

# camera: oblique 3/4 view centred on the board
cam = GetActiveCamera()
cam.SetFocalPoint(3.0, -8.0, 10.0)
cam.SetPosition(120.0, -140.0, 95.0)
cam.SetViewUp(0.0, 0.0, 1.0)
cam.SetViewAngle(30)
ResetCamera()
cam.Zoom(1.15)

SaveScreenshot('fig_board3d.png', view, ImageResolution=[1200, 900])
print('saved fig_board3d.png')

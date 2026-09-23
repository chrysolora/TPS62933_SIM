import sys
from paraview.simple import *

src = OpenDataFile("case_t0001.vtu")
src.UpdatePipeline()

cd = src.CellData
pd = src.PointData
print("== CELL arrays ==")
for i in range(cd.GetNumberOfArrays()):
    a = cd.GetArray(i)
    if a and a.GetName():
        try:
            r = a.GetRange()
        except Exception:
            r = None
        print("  ", a.GetName(), a.GetNumberOfComponents(), r)
print("== POINT arrays ==")
for i in range(pd.GetNumberOfArrays()):
    a = pd.GetArray(i)
    if a and a.GetName():
        try:
            r = a.GetRange()
        except Exception:
            r = None
        print("  ", a.GetName(), a.GetNumberOfComponents(), r)

# choose field: joule heating (cell)
field = "joule heating e"

# ---- top view: horizontal slice through the copper plane ----
view = GetActiveViewOrCreate('RenderView')
view.ViewSize = [1100, 850]
view.Background = [1.0, 1.0, 1.0]

sl = Slice(Input=src)
sl.SliceType = 'Plane'
sl.SliceType.Origin = [0.0, 0.0, 0.0]
sl.SliceType.Normal = [0.0, 0.0, 1.0]
sl.UpdatePipeline()

disp = Show(sl, view)
ColorBy(disp, ('CELLS', field))
lut = GetColorTransferFunction(field)
lut.ApplyPreset('Rainbow Desaturated', True)
lut.ResetCamera()

# camera top-down
cam = view.GetActiveCamera()
cam.SetPosition(0, 0, 0.03)
cam.SetFocalPoint(0, 0, 0)
cam.SetViewUp(0, 1, 0)
view.ResetCamera()
cam.Zoom(1.3)
Render()

print("saving top slice 1")
SaveScreenshot("fea_field_top.png", view, ImageResolution=[1100, 850])

# ---- 3D-ish view with slice tilted ----
cam.SetPosition(0.020, -0.020, 0.018)
cam.SetFocalPoint(0, 0, 0.001)
cam.SetViewUp(0, 0, 1)
view.ResetCamera()
cam.Zoom(1.2)
Render()
print("saving 3d 2")
SaveScreenshot("fea_field_3d.png", view, ImageResolution=[1100, 850])
print("done")

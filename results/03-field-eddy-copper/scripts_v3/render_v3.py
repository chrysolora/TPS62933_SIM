from paraview.simple import *

src = OpenDataFile("mesh/copper.vtu")
src.UpdatePipeline()

v = GetActiveViewOrCreate("RenderView")
v.UseColorPaletteForBackground = 0
v.BackgroundColorMode = "Single Color"
v.Background = [1.0, 1.0, 1.0]
v.OrientationAxesVisibility = 0

d = Show(src, v)
d.Representation = "Surface"
d.Ambient = 0.35
d.Diffuse = 0.75
d.Specular = 0.1

ColorBy(d, ("POINTS", "Jmag"))
lut = GetColorTransferFunction("Jmag")
lut.ApplyPreset("Turbo", True)
lut.UseLogScale = 1
lut.RescaleTransferFunction(1e-1, 3e7)
sb = GetScalarBar(lut, v)
sb.Title = "|J| (A/m^2)"
sb.TitleFontSize = 16
sb.LabelFontSize = 13
sb.Orientation = "Vertical"

CX, CY = 0.006908, -0.007340   # copper bbox centre

# ---------- top view ----------
v.ViewSize = [1000, 1400]
cam = v.GetActiveCamera()
cam.ParallelProjectionOn()
cam.SetFocalPoint(CX, CY, -17e-6)
cam.SetPosition(CX, CY, 0.08)
cam.SetViewUp(0, 1, 0)
cam.SetParallelScale(0.0305)
v.Background = [1.0, 1.0, 1.0]
Render()
SaveScreenshot("fig_field_top_v3.png", v, ImageResolution=[1000, 1400])
print("saved top")

# ---------- 3d view ----------
v.ViewSize = [1200, 900]
cam.SetFocalPoint(CX, CY, 0.0)
cam.SetPosition(CX + 0.030, CY - 0.055, 0.030)
cam.SetViewUp(0, 0, 1)
cam.SetParallelScale(0.040)
Render()
SaveScreenshot("fig_field_3d_v3.png", v, ImageResolution=[1200, 900])
print("saved 3d")

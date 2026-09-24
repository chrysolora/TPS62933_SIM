#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate Gmsh .geo for the SW-copper <-> reference-plane capacitance model.
2D-shell electrodes (never 3D solids for the 35um copper)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from geom_defs import SW_POLY_MIL, KEEPOUT_A, KEEPOUT_B, clip_halfplane, MIL

OUT = os.path.dirname(os.path.abspath(__file__))
BX0, BX1 = 0.0, 1043.3071 * MIL
BY0, BY1 = -1958.189 * MIL, 30.0 * MIL

def sw_m(variant):
    poly = list(SW_POLY_MIL)
    if variant in ("topcut", "dualcut"):
        poly = clip_halfplane(poly, lambda p: p[1] <= -797.795)
    return [(x * MIL, y * MIL) for x, y in poly]

def rect_m(rect):
    x0, x1, y0, y1 = rect
    return x0 * MIL, x1 * MIL, y0 * MIL, y1 * MIL

def add_rect(A, base, x0, y0, x1, y1, z, lc):
    """emit 4 points + 4 lines for a rectangle; return (loop_tag, list_of_line_tags)"""
    A('Point(%d) = {%.9g,%.9g,%.9g,%.9g};' % (base+0, x0, y0, z, lc))
    A('Point(%d) = {%.9g,%.9g,%.9g,%.9g};' % (base+1, x1, y0, z, lc))
    A('Point(%d) = {%.9g,%.9g,%.9g,%.9g};' % (base+2, x1, y1, z, lc))
    A('Point(%d) = {%.9g,%.9g,%.9g,%.9g};' % (base+3, x0, y1, z, lc))
    lines = []
    for k in range(4):
        lt = base + 10 + k
        A('Line(%d) = {%d,%d};' % (lt, base+k, base+(k+1) % 4))
        lines.append(lt)
    return lines

def geo_text(variant, d, lc_sw=0.0008, lc_far=0.03):
    sw = sw_m(variant)
    L = []
    A = L.append
    A('SetFactory("OpenCASCADE");')
    A('Mesh.Algorithm3D = 1;')
    M = max(0.06, 0.6 * d)
    ax0, ax1 = BX0 - M, BX1 + M
    ay0, ay1 = BY0 - M, BY1 + M
    az0, az1 = -d, 0.03
    A('Box(1) = {%.9g,%.9g,%.9g,%.9g,%.9g,%.9g};' % (ax0, ay0, az0, ax1-ax0, ay1-ay0, az1-az0))
    swx0 = min(p[0] for p in sw); swx1 = max(p[0] for p in sw)
    swy0 = min(p[1] for p in sw); swy1 = max(p[1] for p in sw)
    # ---- SW polygon (2D shell) ----
    pid = 100
    for k, (x, y) in enumerate(sw):
        A('Point(%d) = {%.9g,%.9g,0,%.9g};' % (pid+k, x, y, lc_sw))
    n = len(sw)
    swlines = []
    for k in range(n):
        A('Line(%d) = {%d,%d};' % (1000+k, pid+k, pid+(k+1) % n))
        swlines.append(1000+k)
    A('Curve Loop(5000) = {%s};' % ",".join(str(t) for t in swlines))
    A('Plane Surface(5001) = {5000};')
    # ---- bottom GND plane (board rect, optionally with L2 carve holes) ----
    if variant == "dualcut":
        # bottom GND = board rect minus L2 keepout A (interior) -> 4 rectangles
        h = rect_m(KEEPOUT_A)
        cx0, cx1, cy0, cy1 = h[0], h[1], h[2], h[3]
        rects = [(BX0, BY0, cx0, BY1), (cx1, BY0, BX1, BY1),
                 (cx0, BY0, cx1, cy0), (cx0, cy1, cx1, BY1)]
        loop_tags = []
        for i, (x0, y0, x1, y1) in enumerate(rects):
            base = 10000 + i * 100
            ll = add_rect(A, base, x0, y0, x1, y1, -0.001, lc_far)
            lt = 5001 + i
            A('Curve Loop(%d) = {%s};' % (lt, ",".join(str(t) for t in ll)))
            loop_tags.append(lt)
        for i, lt in enumerate(loop_tags):
            A('Plane Surface(%d) = {%d};' % (5002 + i, lt))
        gndsurf = 5002  # first; gnd selected by bbox later
    else:
        lo = add_rect(A, 10000, BX0, BY0, BX1, BY1, -0.001, lc_far)
        A('Curve Loop(5001) = {%s};' % ",".join(str(t) for t in lo))
        A('Plane Surface(5002) = {5001};')
        gndsurf = 5002
    A('BooleanFragments{ Volume{1}; Delete; }{ Surface{5001,%d}; Delete; }' % gndsurf)
    A('Physical Surface("sw") = Surface In BoundingBox{%.9g,%.9g,%.9g,%.9g,%.9g,%.9g};'
      % (swx0-1e-5, swy0-1e-5, -1e-5, swx1+1e-5, swy1+1e-5, 1e-5))
    A('Physical Surface("gnd") = Surface In BoundingBox{%.9g,%.9g,%.9g,%.9g,%.9g,%.9g};'
      % (BX0-1e-5, BY0-1e-5, -0.001-1e-5, BX1+1e-5, BY1+1e-5, -0.001+1e-5))
    A('Physical Surface("ref") = Surface In BoundingBox{%.9g,%.9g,%.9g,%.9g,%.9g,%.9g};'
      % (ax0-1e-4, ay0-1e-4, az0-1e-4, ax1+1e-4, ay1+1e-4, az0+1e-4))
    A('Physical Volume("air") = {1};')
    A('Field[1] = Box;')
    A('Field[1].VIn = %.9g;' % lc_sw)
    A('Field[1].VOut = %.9g;' % lc_far)
    A('Field[1].XMin = %.9g;' % (swx0-0.004))
    A('Field[1].XMax = %.9g;' % (swx1+0.004))
    A('Field[1].YMin = %.9g;' % (swy0-0.004))
    A('Field[1].YMax = %.9g;' % (swy1+0.004))
    A('Field[1].ZMin = -0.004;')
    A('Field[1].ZMax = 0.004;')
    A('Field[2] = Min;')
    A('Field[2].FieldsList = {1};')
    A('Background Field = 2;')
    A('Mesh.MeshSizeExtendFromBoundary = 0;')
    A('Mesh.MeshSizeFromPoints = 0;')
    A('Mesh.MeshSizeFromCurvature = 0;')
    return "\n".join(L) + "\n"

if __name__ == "__main__":
    variant = sys.argv[1]; d = float(sys.argv[2])
    lc = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0008
    tag = sys.argv[4] if len(sys.argv) > 4 else ""
    fn = os.path.join(OUT, "geo_%s_d%03d%s.geo" % (variant, int(d*1000), tag))
    open(fn, "w").write(geo_text(variant, d, lc))
    print(fn)

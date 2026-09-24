#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Real-copper geometry for TPS62933 CM parasitic capacitance (task cm_pcb).
Units: coordinates defined in mil, converted to SI metres for the FEM.
Sources (verified, not guessed):
  board outline : rad/geom/board_outline.csv
  SW net pour   : epro/pourSim.epru POUR75 (netName $1N7)  [U21.SW pin5 -> L2]
  L2 keepouts   : fea4/<ver>/geo_meta.json (REG_PROHIBIT bboxes)
"""
MIL = 0.0254e-3
BOARD_MIL = (0.0, 1043.3071, -1958.189, 30.0)   # x0,x1,y0,y1

SW_POLY_MIL = [
    (315.0,-905.0),(285.0,-905.0),(255.0,-875.0),(255.0,-872.4),
    (255.0,-847.4),(242.6,-835.0),(195.0,-835.0),(175.0,-855.0),
    (175.0,-900.0),(170.0,-905.0),(150.0,-905.0),(142.5,-897.5),
    (142.5,-851.7),(152.5,-827.5),(160.0,-820.0),(180.0,-771.8),
    (180.0,-740.0),(190.0,-730.0),(305.6,-730.0),(315.0,-739.4),(315.0,-740.0),
]
# L2 copper-prohibit region A (x0,x1,y0,y1) mil  (6.60 x 7.00 mm over L2)
KEEPOUT_A = (115.079, 374.921, -797.795, -522.205)
# region B (present in topcut/dualcut) - far from SW pour, listed for completeness
KEEPOUT_B = (314.32, 485.68, -680.0, -500.43)

def area_mm2(poly_mil):
    n=len(poly_mil); A=0.0
    for i in range(n):
        x1,y1=poly_mil[i]; x2,y2=poly_mil[(i+1)%n]; A+=x1*y2-x2*y1
    return abs(A)/2*(MIL*1000)**2

def clip_halfplane(poly, inside):
    out=[]; n=len(poly)
    def cross(a,b):
        lo,hi=a,b
        for _ in range(80):
            m=((lo[0]+hi[0])/2,(lo[1]+hi[1])/2)
            if inside(m)==inside(lo): lo=m
            else: hi=m
        return ((lo[0]+hi[0])/2,(lo[1]+hi[1])/2)
    for i in range(n):
        a=poly[i]; b=poly[(i+1)%n]; ia=inside(a); ib=inside(b)
        if ia:
            out.append(a)
            if not ib: out.append(cross(a,b))
        elif ib: out.append(cross(a,b))
    return out

def sw_polygon(variant):
    """SW net copper polygon(s) in mil for a given pour variant."""
    poly=SW_POLY_MIL
    if variant in ("topcut","dualcut"):
        # remove copper inside keepout A (a Y-cut here: y must be <= -797.795)
        poly=clip_halfplane(poly, lambda p: p[1] <= -797.795)
    return [poly]

if __name__=="__main__":
    import json
    out={}
    print("SW pour area as designed = %.4f mm^2"%area_mm2(SW_POLY_MIL))
    for v in ("fullcu","topcut","dualcut"):
        polys=sw_polygon(v); A=sum(area_mm2(p) for p in polys)
        out[v]=A
        print("%-8s SW area = %.4f mm^2 (n=%d)"%(v,A,len(polys)))

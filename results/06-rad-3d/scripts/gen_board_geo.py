#!/usr/bin/env python3
"""Generate board.geo (OpenCASCADE) -> board.step for the TPS62933 board.

Geometry source: rad/geom/{board_outline.csv,components.csv}.
  * board solid   : exact outline extruded to FR4 thickness 1.51 mm (exact)
  * component boxes: package-code based footprint sizes (ESTIMATED, see table)
  * L2 (ZEMS0650) : footprint taken from the L2 keep-out region in
                    rad/b1/board_b1.py (=6.6 x 7.0 mm, exact w.r.t. that model),
                    height estimated.

All component heights & most footprints are ESTIMATES from the package code
(see COMP_SIZE below).  This STEP is *reconstructed from exported board
geometry*, NOT the vendor's original mechanical STEP.
"""
import csv, os, math

HERE = os.path.dirname(os.path.abspath(__file__))
GEOM = os.path.abspath(os.path.join(HERE, '..', 'rad', 'geom'))

H_FR4 = 1.51
z_bot, z_top = -H_FR4 / 2, H_FR4 / 2

# package-code -> (L, W, H) footprint size in mm.  ESTIMATES (EIA/industry std).
COMP_SIZE = {
    '0603WAF3572T5E': (1.6, 0.8, 0.6), '0603WAF3833T5E': (1.6, 0.8, 0.6),
    '0603WAF2672T5E': (1.6, 0.8, 0.6), '0603CG2R7C500NT': (1.6, 0.8, 0.9),
    '0603WAF1002T5E': (1.6, 0.8, 0.6), '0603WAF1403T5E': (1.6, 0.8, 0.6),
    '0603WAF4422T5E': (1.6, 0.8, 0.6),
    'CC0603KRX7R9BB104': (1.6, 0.8, 0.9),
    'CL31A106KBHNNNE': (3.2, 1.6, 1.0),          # 1206 MLCC
    '1206W4F100LT5E': (3.2, 1.6, 0.7),
    'HGC1210R5476M250NSVK': (3.2, 2.5, 2.5),     # 1210 polymer cap
    'XL-1608UGC-04': (1.6, 0.8, 0.6),            # 1608 LED
    'ZEMS0650-150M': (6.6, 7.0, 3.0),            # L2 power inductor (footprint=keepout)
    'TPS62933FDRLR': (1.6, 1.6, 0.6),            # SOT-563
    'SS36_C7420367': (4.3, 2.6, 2.2),            # SMA / DO-214AC
    'JFC2410-1200TS': (6.1, 2.5, 2.5),           # 2410 fuse
    'FXL0420-1R0-M': (4.2, 4.2, 1.8),
    'JVJ50v10M5x5': (5.3, 5.3, 5.4),             # 5x5 alu electrolytic
    'KF2EDGV-3.81-2P': (7.6, 7.5, 7.0),          # pluggable terminal
    'SMF24A_C19077514': (2.6, 1.5, 1.0),         # SMF / DO-219AB
    'Test-Point': (1.0, 1.0, 0.5),
}

def read_outline():
    with open(os.path.join(GEOM, 'board_outline.csv')) as f:
        xs, ys = [], []
        for r in csv.DictReader(f):
            xs.append(float(r['x_m'])); ys.append(float(r['y_m']))
    return min(xs) * 1e3, max(xs) * 1e3, min(ys) * 1e3, max(ys) * 1e3

def read_components():
    out = []
    with open(os.path.join(GEOM, 'components.csv')) as f:
        for r in csv.DictReader(f):
            out.append((r['device'], float(r['x_m']) * 1e3, float(r['y_m']) * 1e3,
                        float(r['angle_deg']), int(r['layerId'])))
    return out

def main():
    X0, X1, Y0, Y1 = read_outline()
    comps = read_components()
    geo = []
    geo.append('SetFactory("OpenCASCADE");\n')
    n = 0
    # board solid (exact outline)
    n += 1
    geo.append('Box(%d) = {%.6f,%.6f,%.6f, %.6f,%.6f,%.6f};\n' %
               (n, X0, Y0, z_bot, X1 - X0, Y1 - Y0, H_FR4))
    misses = []
    for dev, cx, cy, ang, layer in comps:
        size = COMP_SIZE.get(dev)
        if size is None:
            misses.append(dev); size = (2.0, 2.0, 1.0)
        L, W, H = size
        n += 1
        zz = z_top if layer == 1 else z_bot - H
        geo.append('Box(%d) = {%.6f,%.6f,%.6f, %.6f,%.6f,%.6f};\n' %
                   (n, cx - L / 2, cy - W / 2, zz, L, W, H))
        if abs(ang) > 1e-6:
            geo.append('Rotate { {0,0,1}, {%.6f,%.6f,%.6f}, %.6f } In Volume{%d};\n' %
                       (cx, cy, zz, math.radians(ang), n))
    ids = ','.join(str(i) for i in range(1, n + 1))
    geo.append('Compound Volume{%s};\n' % ids)
    with open(os.path.join(HERE, 'board.geo'), 'w') as f:
        f.write(''.join(geo))
    print('board outline mm: X[%.3f,%.3f] Y[%.3f,%.3f] size %.2fx%.2f' %
          (X0, X1, Y0, Y1, X1 - X0, Y1 - Y0))
    print('volumes: 1 board + %d components' % (n - 1))
    if misses:
        print('WARNING default size used for:', sorted(set(misses)))

if __name__ == '__main__':
    main()

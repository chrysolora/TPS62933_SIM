#!/usr/bin/env python3
"""Radiated-emission model for the TPS62933 buck -- three pour variants.

  --ver full | top | dual     (FullCopper / TopCutout / DualCutout)
  --fast                      coarse/short validation run

Only the *copper geometry* changes between variants (L2 keepout region on the
top layer, the bottom layer, or both). Source, loop, cable, mesh strategy,
boundaries and NF2FF are identical.

Copper is built as thin PEC sheets (no solid-body meshing):
  * top  = signal pours (from rad/geom/poured_polygons.json)
           + top GND pour  = board rect - signal-pour clearances
                             - hot-loop clearance - L2 keepout (top variants)
  * bottom = board rect - L2 keepout (dual variant)
Approximation: GND pours are reconstructed as board rectangle minus signal-pour
bounding boxes (inflated by a clearance) -- the exact GND polygons were not
exported.  Documented as a B1 limitation.
"""
import os, sys, csv, json, time, shutil
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS

FAST = '--fast' in sys.argv
VER  = 'full'
for a in sys.argv:
    if a.startswith('--ver'):
        VER = a.split('=')[1] if '=' in a else sys.argv[sys.argv.index(a) + 1]
CUT_TOP, CUT_BOT = {'full': (False, False), 'top': (True, False),
                    'dual': (True, True)}[VER]

HERE = os.path.dirname(os.path.abspath(__file__))
GEOM = os.path.abspath(os.path.join(HERE, '..', 'geom'))
OUT  = HERE                     # write directly into rad/b1/
os.makedirs(OUT, exist_ok=True)

# ---------------- board / stack-up (mm) ----------------
X0, X1 = -6.223, 20.277
Y0, Y1 = -32.974, 17.526
H_FR4  = 1.51; z_top = H_FR4/2; z_bot = -H_FR4/2
EPS_FR4 = 4.5
PX, PY = 4.0, -6.3
CX, CY = -1.22, -28.97
CABLE_L = 100.0
F_MIN, F_MAX = 30e6, 1e9

# L2 keepout region (mil -> mm): x[115.079,374.921] y[-797.795,-522.205]
KEEP = (115.079*0.0254 - 245*0.0254, 374.921*0.0254 - 245*0.0254,
        -797.795*0.0254 + 660*0.0254, -522.205*0.0254 + 660*0.0254)
LOOPCLEAR = (1.0, 7.0, -7.8, -4.8)          # copper-free slot for the source loop
CLR = 0.3                                   # pour clearance around signal pours
BOARD = (X0, X1, Y0, Y1)

def rects(name):
    out = []
    with open(os.path.join(GEOM, name)) as f:
        for r in csv.DictReader(f):
            out.append((r['name'], r['net'], float(r['x0_m'])*1e3, float(r['x1_m'])*1e3,
                        float(r['y0_m'])*1e3, float(r['y1_m'])*1e3))
    return out

top_pours = rects('pours_top.csv')
# json signal polygons (not the GND plane)
poly = json.load(open(os.path.join(GEOM, 'poured_polygons.json')))
def mil2mm(pts):
    return np.array([[(p[0]-245.0)*25.4e-3, (p[1]+660.0)*25.4e-3] for p in pts])
sig_polys = []
for p in poly:
    pts = mil2mm(p['pts'])
    if len(pts) < 3: continue
    wx, wy = pts[:, 0].ptp(), pts[:, 1].ptp()
    cx, cy = pts[:, 0].mean(), pts[:, 1].mean()
    if wx < 0.3 and wy < 0.3: continue
    if not (X0-2 < cx < X1+2 and Y0-2 < cy < Y1+2): continue   # drop off-board junk
    sig_polys.append(pts)

def xsub(a, k):
    ax0, ax1, ay0, ay1 = a; kx0, kx1, ky0, ky1 = k
    if ax1 <= kx0 or ax0 >= kx1 or ay1 <= ky0 or ay0 >= ky1:
        return [a]
    o = []
    if ay0 < ky0: o.append((ax0, ax1, ay0, min(ay1, ky0)))
    if ay1 > ky1: o.append((ax0, ax1, max(ay0, ky1), ay1))
    yy0, yy1 = max(ay0, ky0), min(ay1, ky1)
    if yy1 > yy0:
        if ax0 < kx0: o.append((ax0, min(ax1, kx0), yy0, yy1))
        if ax1 > kx1: o.append((max(ax0, kx1), ax1, yy0, yy1))
    return [r for r in o if r[1]-r[0] > 1e-4 and r[3]-r[2] > 1e-4]

def sub_many(r, ks):
    cur = [r]
    for k in ks:
        nxt = []
        for c in cur: nxt += xsub(c, k)
        cur = nxt
    return cur

# ---- top GND = board - signal clearances - loopclear - keepout(if cut_top) ----
gnd_excl = [(x0-CLR, x1+CLR, y0-CLR, y1+CLR) for (_, net, x0, x1, y0, y1) in top_pours
            if 'GND' not in (net or '')]
excl = gnd_excl + [LOOPCLEAR]
if CUT_TOP: excl += [KEEP]
top_gnd_rects = sub_many(BOARD, excl)
# signal pours kept only if they do not overlap the loop clearance
def overlap(a, b):
    return not (a[1] <= b[0] or a[0] >= b[1] or a[3] <= b[2] or a[2] >= b[3])
lc = (LOOPCLEAR[0], LOOPCLEAR[1], LOOPCLEAR[2], LOOPCLEAR[3])
sig_polys = [p for p in sig_polys
             if not overlap((p[:,0].min(),p[:,0].max(),p[:,1].min(),p[:,1].max()), lc)]

# ---- bottom GND = board - keepout(if cut_bot) ----
bot_rects = sub_many(BOARD, [KEEP]) if CUT_BOT else [BOARD]

print('ver=%s cut_top=%s cut_bot=%s | top_gnd pieces=%d sig_polys=%d bot pieces=%d'
      % (VER, CUT_TOP, CUT_BOT, len(top_gnd_rects), len(sig_polys), len(bot_rects)))

# ================= FDTD =================
CSX = ContinuousStructure()
NrTS = 40000 if FAST else 160000
FDTD = openEMS(NrTS=NrTS, EndCriteria=1e-4, CellConstantMaterial=True)
FDTD.SetCSX(CSX)
FDTD.SetGaussExcite(0.5*(F_MIN+F_MAX), 0.5*(F_MAX-F_MIN))
FDTD.SetBoundaryCond(['PML_8']*6)

sub = CSX.AddMaterial('FR4', epsilon=EPS_FR4, kappa=1e-3)
sub.AddBox(priority=0, start=[X0, Y0, z_bot], stop=[X1, Y1, z_top])

gnd_b = CSX.AddMetal('gnd_bot')
for (x0, x1, y0, y1) in bot_rects:
    gnd_b.AddBox(priority=10, start=[x0, y0, z_bot], stop=[x1, y1, z_bot])

cu_t = CSX.AddMetal('cu_top')
for (x0, x1, y0, y1) in top_gnd_rects:
    cu_t.AddBox(priority=10, start=[x0, y0, z_top], stop=[x1, y1, z_top])
for p in sig_polys:
    cu_t.AddPolygon(points=[p[:,0].tolist(), p[:,1].tolist()],
                    norm_dir=2, elevation=z_top, priority=10)

cable = CSX.AddMetal('cable')
cable.AddBox(priority=11, start=[CX-0.4, CY-0.4, z_top], stop=[CX+0.4, CY+0.4, z_top+CABLE_L])

# ---- hot loop + small source ----
VHX = 2.5; GAP = 0.25; Xa, Xb = PX-VHX, PX+VHX
loop = CSX.AddMetal('loop')
loop.AddBox(priority=12, start=[Xa, PY-0.15, z_top], stop=[PX, PY+0.15, z_top])
loop.AddBox(priority=12, start=[PX+GAP, PY-0.15, z_top], stop=[Xb, PY+0.15, z_top])
for xx in (Xa, Xb):
    loop.AddBox(priority=12, start=[xx-0.15, PY-0.15, z_bot], stop=[xx+0.15, PY+0.15, z_top])
port = FDTD.AddLumpedPort(1, 50.0, [PX, PY, z_top], [PX+GAP, PY, z_top], 'x',
                          excite=1.0, priority=5, edges2grid="yz")

# ---- mesh ----
mesh = CSX.GetGrid(); mesh.SetDeltaUnit(1e-3)
SNAP = 0.25
def addl(d, vals, snap=SNAP):
    v = np.round(np.asarray(vals, float)/snap)*snap
    mesh.AddLine(d, sorted(set(np.round(v, 6).tolist())))

resb = 0.8 if not FAST else 1.5
bx = list(np.arange(X0, X1+1e-9, resb)) + [X0, X1, PX, CX, PX+GAP, Xa-0.15, Xa+0.15, Xb-0.15, Xb+0.15]
by = list(np.arange(Y0, Y1+1e-9, resb)) + [Y0, Y1, PY, CY, PY-0.15, PY+0.15]
for (x0, x1, y0, y1) in top_gnd_rects + bot_rects:
    bx += [x0, x1]; by += [y0, y1]
for p in sig_polys:
    bx += [p[:,0].min(), p[:,0].max()]; by += [p[:,1].min(), p[:,1].max()]
addl(0, bx); addl(1, by)
bz = [z_bot, z_top, 0.0] + list(np.arange(z_bot, z_top+1e-9, 0.5)) \
     + list(np.arange(z_top, z_top+CABLE_L+1e-9, 4.0))
addl(2, bz, 1e-4)
AX, AY = (80, 90) if FAST else (130, 150)
ZLO, ZHI = (-80, 110) if FAST else (-130, 175)
addl(0, [-AX, AX]); addl(1, [-AY, AY]); addl(2, [ZLO, ZHI], 1e-4)
mesh.SmoothMeshLines('all', 15.0 if not FAST else 22.0, 1.5)
g = [np.array(mesh.GetLines(d)) for d in range(3)]
print('mesh x/y/z %d/%d/%d cells=%d' % (len(g[0]), len(g[1]), len(g[2]),
      (len(g[0])-1)*(len(g[1])-1)*(len(g[2])-1)))

nf2ff = FDTD.CreateNF2FFBox(start=[-35, -55, max(-45, ZLO+20)],
                            stop=[40, 55, min(125, ZHI-20)])
Sim = os.path.join(OUT, 'sim')
if os.path.exists(Sim):
    shutil.rmtree(Sim)
os.makedirs(Sim, exist_ok=True)
t0 = time.time()
FDTD.Run(Sim, cleanup=(not FAST), verbose=1, numthreads=12)
print('run done %.0f s' % (time.time()-t0))

freq = np.linspace(F_MIN, F_MAX, 401)
port.CalcPort(Sim, freq)
np.savez(os.path.join(OUT, 'port.npz'), freq=freq, uf_inc=port.uf_inc,
         if_tot=port.if_tot, uf_tot=port.uf_tot, P_acc=port.P_acc)
theta = np.arange(0, 181, 5.0); phi = np.arange(0, 360, 10.0)
nf = nf2ff.CalcNF2FF(Sim, freq, theta, phi, center=[0, -7.5, 0], verbose=0)
np.savez(os.path.join(OUT, 'nf2ff.npz'), freq=np.asarray(nf.freq),
         theta=np.asarray(nf.theta), phi=np.asarray(nf.phi),
         Dmax=np.asarray(nf.Dmax), Prad=np.asarray(nf.Prad),
         E_norm=np.asarray(nf.E_norm).reshape(len(freq), len(theta), len(phi)))
i100 = int(np.argmin(np.abs(freq-1e8)))
print('ver=%s uf_inc@100M=%.3e if_tot@100M=%.3e Dmax@100M=%.3f' % (
      VER, abs(port.uf_inc[i100]), abs(port.if_tot[i100]), nf.Dmax[i100]))
print('done', VER)

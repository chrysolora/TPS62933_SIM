#!/usr/bin/env python3
"""B2 radiated-emission model for the TPS62933 buck (openEMS).

Extends B1 with the three B1 gaps:

  1. upstream supply ripple as an *independent* source
       LM50-20B24 (金升阳), fsw = 65 kHz typ, ripple <=150 mVpp (datasheet
       worst) / ~100 mVpp (user).  CLASS-I module -> ripple reaches the
       outside world as a **common-mode** drive through the internal
       Y-capacitor <-> case/PE path.
  2. 10 cm / 0.75 mm^2 copper cable, series inductance realised physically
       (1.0 mm cross-section conductor) and given analytically
       L = 2e-7 * l * (ln(2l/a) - 0.75)  ~= 105 nH.
  3. reference ground plane (CISPR-like) + the CM return path.

>>> Key wiring fact (user, measured): the harness to the test board has ONLY
    +24 V and GND -- there is NO PE / earth wire at the board.  The CM loop
    therefore closes as:
        board GND - cable - Ycap(upstream) - case/PE -[ C_p ]- back to board
    i.e. C_p (board <-> earth parasitic capacitance) is the only return
    element and is modelled EXPLICITLY (sensitivity 2 / 10 / 50 pF).

Two excitations per variant (identical geometry / mesh):
  --excite dm   hot-loop lumped port (port 1)
  --excite cm   cable lumped port (port 2, source in series in the CM loop)

usage: build_rad_b2.py --ver=full|top|dual --excite=dm|cm [--cp=2|10|50] [--fast]
"""
import os, sys, csv, json, time, shutil
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS

FAST = '--fast' in sys.argv
VER  = 'full'; EXC = 'dm'; CP = 10.0
for a in sys.argv:
    if a.startswith('--ver='):    VER = a.split('=')[1]
    if a.startswith('--excite='): EXC = a.split('=')[1]
    if a.startswith('--cp='):     CP  = float(a.split('=')[1])
assert VER in ('full', 'top', 'dual'), VER
assert EXC in ('dm', 'cm'), EXC
CUT_TOP, CUT_BOT = {'full': (False, False), 'top': (True, False),
                    'dual': (True, True)}[VER]

HERE = os.path.dirname(os.path.abspath(__file__))
GEOM = '/mnt/raid10/sim-work/tps62933/rad/geom'
TAG  = '%s_%s%s' % (VER, EXC, '_pe' if '--pe' in sys.argv else '')
OUT  = '/mnt/raid10/sim-work/tps62933/rad_conv/ver_' + TAG + os.environ.get('SUF','')
if os.path.exists(OUT):
    shutil.rmtree(OUT)
os.makedirs(OUT, exist_ok=True)

# ---------------- board / stack-up (mm) ----------------
Q = 0.25   # FDTD base grid [mm]; ALL geometry quantized to this -> no geometry-induced micro-cells
def q(v): return float(round(float(v)/Q)*Q)
X0, X1 = q(-6.223), q(20.277)
Y0, Y1 = q(-32.974), q(17.526)
H_FR4  = 1.50; z_top = H_FR4/2; z_bot = -H_FR4/2   # 1.51->1.50mm rounded to grid (see report)
EPS_FR4 = 4.5
PX, PY = 4.0, q(-6.3)         # hot-loop centre
CX, CY = q(-1.22), q(-28.97)  # input terminal block (KF2EDGV-3.81-2P)
F_MIN, F_MAX = 30e6, 1e9

# ---- B2 additions ----
CABLE_L  = 100.0              # 10 cm cable
CABLE_W  = 1.0                # conductor cross-section [mm] (~0.98 mm eq. dia.)
YE       = q(CY - CABLE_L)       # far (supply) end of the cable
CM_GAP   = 0.5                # lumped-port gap (y direction) in the cable
Cd       = 0.5                # lumped-element gap height [mm]
Z_GP     = q(z_bot - 10.0)       # reference ground plane (10 mm below board)
GPX      = (q(-88.0), q(88.0)); GPY = (q(-140.0), q(30.0))   # plane (176 x 170 mm, reduced)
CY_YC    = 2.2e-9             # upstream internal Y-cap (assumed, 2.2 nF/line)
C_GEOM   = 1.19e-12           # board-pour <-> plane geometric cap (analytic)
CP_LUMP  = max(CP*1e-12 - C_GEOM, 0.1e-12)   # top-up so total C_p == CP

# L2 keepout (mil -> mm)
KEEP = (115.079*0.0254 - 245*0.0254, 374.921*0.0254 - 245*0.0254,
        -797.795*0.0254 + 660*0.0254, -522.205*0.0254 + 660*0.0254)
LOOPCLEAR = (1.0, 7.0, q(-7.8), q(-4.8))
CLR = 0.25
BOARD = (X0, X1, Y0, Y1)

def rects(name):
    out = []
    with open(os.path.join(GEOM, name)) as f:
        for r in csv.DictReader(f):
            out.append((r['name'], r['net'], q(float(r['x0_m'])*1e3), q(float(r['x1_m'])*1e3),
                        q(float(r['y0_m'])*1e3), q(float(r['y1_m'])*1e3)))
    return out

top_pours = rects('pours_top.csv')
poly = json.load(open(os.path.join(GEOM, 'poured_polygons.json')))
def mil2mm(pts):
    return np.array([[(p[0]-245.0)*25.4e-3, (p[1]+660.0)*25.4e-3] for p in pts])
sig_polys = []
for p in poly:
    pts = np.round(mil2mm(p['pts'])/Q)*Q
    if len(pts) < 3: continue
    wx, wy = pts[:, 0].ptp(), pts[:, 1].ptp()
    cx, cy = pts[:, 0].mean(), pts[:, 1].mean()
    if wx < 0.3 and wy < 0.3: continue
    if not (X0-2 < cx < X1+2 and Y0-2 < cy < Y1+2): continue
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

gnd_excl = [(x0-CLR, x1+CLR, y0-CLR, y1+CLR) for (_, net, x0, x1, y0, y1) in top_pours
            if 'GND' not in (net or '')]
excl = gnd_excl + [LOOPCLEAR]
if CUT_TOP: excl += [KEEP]
top_gnd_rects = sub_many(BOARD, excl)

def overlap(a, b):
    return not (a[1] <= b[0] or a[0] >= b[1] or a[3] <= b[2] or a[2] >= b[3])
lc = (LOOPCLEAR[0], LOOPCLEAR[1], LOOPCLEAR[2], LOOPCLEAR[3])
sig_polys = [p for p in sig_polys
             if not overlap((p[:,0].min(),p[:,0].max(),p[:,1].min(),p[:,1].max()), lc)]

bot_rects = sub_many(BOARD, [KEEP]) if CUT_BOT else [BOARD]

print('ver=%s excite=%s cp=%.1fpF cut_top=%s cut_bot=%s | top_gnd=%d sig=%d bot=%d'
      % (VER, EXC, CP, CUT_TOP, CUT_BOT, len(top_gnd_rects), len(sig_polys), len(bot_rects)))

# ================= FDTD =================
CSX = ContinuousStructure()
NrTS = int(os.environ.get('NRTS', 40000 if FAST else 120000))
FDTD = openEMS(NrTS=NrTS, EndCriteria=float(os.environ.get("ENDCRIT",1e-9)), CellConstantMaterial=True)
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

# ---- reference ground plane (PEC thin sheet) ----
gp = CSX.AddMetal('ref_gnd_plane')
gp.AddBox(priority=11, start=[GPX[0], GPY[0], Z_GP], stop=[GPX[1], GPY[1], Z_GP])

# ---- 10 cm cable (physical conductor) + CM return post to the plane ----
cbl = CSX.AddMetal('cable')
cbl.AddBox(priority=12, start=[CX-CABLE_W/2, YE+CM_GAP, z_top-CABLE_W/2],
                             stop=[CX+CABLE_W/2, CY,        z_top+CABLE_W/2])
# board-end via: cable lands on the 24VIN pour; tie it to the GND pour
# (represents the 3x10uF input caps acting as an RF short at HF).
cbl.AddBox(priority=12, start=[CX-CABLE_W/2, CY-CABLE_W/2, z_bot],
                             stop=[CX+CABLE_W/2, CY+CABLE_W/2, z_top+CABLE_W/2])
cbl.AddBox(priority=12, start=[CX-CABLE_W/2, YE-CABLE_W, Z_GP+Cd],
                             stop=[CX+CABLE_W/2, YE,        z_top+CABLE_W/2])

# ---- upstream Y-cap: series lumped C in the CM post (post bottom -> plane) ----
ycap = CSX.AddLumpedElement('Ycap', ny=2, caps=False, C=CY_YC, LEtype=0)
ycap.AddBox(priority=14, start=[CX-CABLE_W/2, YE-CABLE_W, Z_GP],
                         stop=[CX+CABLE_W/2, YE,        Z_GP+Cd])
# ---- board <-> earth parasitic C_p (lumped, at board centre) ----
XCP, YCP = 0.0, -5.0
cbl.AddBox(priority=12, start=[XCP-0.5, YCP-0.5, Z_GP+Cd],
                        stop=[XCP+0.5, YCP+0.5, z_bot])      # post up to board
cpel = CSX.AddLumpedElement('Cp', ny=2, caps=False, C=CP_LUMP, LEtype=0)
cpel.AddBox(priority=14, start=[XCP-0.5, YCP-0.5, Z_GP],
                         stop=[XCP+0.5, YCP+0.5, Z_GP+Cd])
if '--pe' in sys.argv:      # PE-comparison: board shorted to the plane (PEC)
    cbl.AddBox(priority=15, start=[XCP-0.5, YCP-0.5, Z_GP],
                            stop=[XCP+0.5, YCP+0.5, Z_GP+Cd])

# ---- hot loop + small DM source (identical to B1) ----
VHX = 2.5; GAP = 0.25; Xa, Xb = PX-VHX, PX+VHX
loop = CSX.AddMetal('loop')
loop.AddBox(priority=13, start=[Xa, PY-0.25, z_top], stop=[PX, PY+0.25, z_top])
loop.AddBox(priority=13, start=[PX+GAP, PY-0.25, z_top], stop=[Xb, PY+0.25, z_top])
for xx in (Xa, Xb):
    loop.AddBox(priority=13, start=[xx-0.25, PY-0.25, z_bot], stop=[xx+0.25, PY+0.25, z_top])
port_dm = FDTD.AddLumpedPort(1, 50.0, [PX, PY, z_top], [PX+GAP, PY, z_top], 'x',
                             excite=(1.0 if EXC == 'dm' else 0.0), priority=5,
                             edges2grid="yz")
port_cm = FDTD.AddLumpedPort(2, 50.0, [CX, YE, z_top], [CX, YE+CM_GAP, z_top], 'y',
                             excite=(1.0 if EXC == 'cm' else 0.0), priority=5,
                             edges2grid="xz")

# ---- mesh ----
mesh = CSX.GetGrid(); mesh.SetDeltaUnit(1e-3)
SNAP = 0.25
MINGAP = 0.10   # [mm] collapse mesh lines closer than this -> caps dt, kills micro-cells
def addl(d, vals, snap=SNAP, mingap=MINGAP):
    v = np.sort(np.unique(np.round(np.asarray(vals, float)/snap)*snap))
    if len(v) == 0: return
    out = [v[0]]
    for x in v[1:]:
        if x - out[-1] < mingap:
            out[-1] = 0.5*(out[-1] + x)
        else:
            out.append(x)
    mesh.AddLine(d, sorted(set(np.round(np.asarray(out), 6).tolist())))

resb = float(os.environ.get('RESB', 1.0 if not FAST else 1.5))
bx = list(np.arange(X0, X1+1e-9, resb)) + [X0, X1, PX, CX, PX+GAP, Xa-0.15, Xa+0.15,
      Xb-0.15, Xb+0.15, CX-CABLE_W/2, CX+CABLE_W/2, XCP-0.5, XCP+0.5, GPX[0], GPX[1]]
by = list(np.arange(Y0, Y1+1e-9, resb)) + [Y0, Y1, PY, CY, PY-0.25, PY+0.25, YE,
      YE+CM_GAP, YE-CABLE_W, YCP-0.5, YCP+0.5, CY-CABLE_W/2, CY+CABLE_W/2,
      GPY[0], GPY[1]]
for (x0, x1, y0, y1) in top_gnd_rects + bot_rects:
    bx += [x0, x1]; by += [y0, y1]
for p in sig_polys:
    bx += [p[:,0].min(), p[:,0].max()]; by += [p[:,1].min(), p[:,1].max()]
addl(0, bx); addl(1, by)
bz = [z_bot, z_top, 0.0, Z_GP, Z_GP+Cd, z_top-CABLE_W/2, z_top+CABLE_W/2] \
     + list(np.arange(z_bot, z_top+1e-9, 0.5))
addl(2, bz, 0.05)
AX, AY = (160, 210)
ZLO, ZHI = (-90, 170)
addl(0, [-AX, AX]); addl(1, [-AY, AY]); addl(2, [ZLO, ZHI], 0.05)
mesh.SmoothMeshLines('all', float(os.environ.get('SMOOTH', 22.0 if FAST else 18.0)), 1.5)
g = [np.array(mesh.GetLines(d)) for d in range(3)]
print('mesh x/y/z %d/%d/%d cells=%d' % (len(g[0]), len(g[1]), len(g[2]),
      (len(g[0])-1)*(len(g[1])-1)*(len(g[2])-1)))
for dd in range(3):
    df=np.diff(g[dd]); ii=int(np.argmin(df))
    print('dir%d min %.4f at %.4f  (n=%d)'%(dd,df.min(),g[dd][ii],len(g[dd])))
    j=np.argsort(df)[:6]
    print('   smallest gaps:', [(round(float(g[dd][k]),4), round(float(df[k]),4)) for k in j])
print('min cell dx/dy/dz [mm]: %.4f %.4f %.4f' % (np.diff(g[0]).min(),
      np.diff(g[1]).min(), np.diff(g[2]).min()))
if '--meshonly' in sys.argv:
    print('meshonly -> exit'); sys.exit(0)

nf2ff = FDTD.CreateNF2FFBox(start=[-92, -150, -45], stop=[92, 40, 110])
Sim = os.path.join(OUT, 'sim')
os.makedirs(Sim, exist_ok=True)
t0 = time.time()
FDTD.Run(Sim, cleanup=(not FAST), verbose=1, numthreads=2)
print('run done %.0f s' % (time.time()-t0))

freq = np.linspace(F_MIN, F_MAX, 101)
port_dm.CalcPort(Sim, freq)
port_cm.CalcPort(Sim, freq)
np.savez(os.path.join(OUT, 'port.npz'), freq=freq,
         dm_if_tot=port_dm.if_tot, dm_uf_tot=port_dm.uf_tot, dm_uf_inc=port_dm.uf_inc,
         cm_if_tot=port_cm.if_tot, cm_uf_tot=port_cm.uf_tot, cm_uf_inc=port_cm.uf_inc)
theta = np.arange(0, 181, 5.0); phi = np.arange(0, 360, 15.0)
nf = nf2ff.CalcNF2FF(Sim, freq, theta, phi, center=[0, -20, 0], verbose=0)
np.savez(os.path.join(OUT, 'nf2ff.npz'), freq=np.asarray(nf.freq),
         theta=np.asarray(nf.theta), phi=np.asarray(nf.phi),
         Dmax=np.asarray(nf.Dmax), Prad=np.asarray(nf.Prad),
         E_norm=np.asarray(nf.E_norm).reshape(len(freq), len(theta), len(phi)))
i100 = int(np.argmin(np.abs(freq-1e8)))
print('ver=%s exc=%s if_dm@100M=%.3e if_cm@100M=%.3e Dmax@100M=%.3f'
      % (VER, EXC, abs(port_dm.if_tot[i100]), abs(port_cm.if_tot[i100]), nf.Dmax[i100]))
print('done', TAG)

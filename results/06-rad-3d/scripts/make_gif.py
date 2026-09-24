#!/usr/bin/env python3
"""Render the near-field |E| time snapshots into propagation.gif.

Reads the openEMS TD E-field dumps produced by build_prop.py:
   prop/sim/E_sliceV.h5 : x-z plane through the source (side view)
   prop/sim/E_sliceH.h5 : x-y plane 3 mm above the board (top view)
and writes PNG frames + propagation.gif (ffmpeg).

Amplitude is the RAW (uncalibrated) FDTD |E| -> only the qualitative outward
propagation is meaningful (absolute level NOT reliable, see REPORT.md).
Colour maps use a fixed gamma-scaled linear range (bright base colour map),
NOT a log range.
"""
import os, glob, subprocess, shutil
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm, Normalize

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.join(HERE, 'prop', 'sim')
FR = os.path.join(HERE, 'prop', 'frames')
os.makedirs(FR, exist_ok=True)

X0, X1, Y0, Y1 = -6.223, 20.277, -32.974, 17.526
CX, CY = -1.22, -28.97
Z_BOARD = 0.755

def load(fn):
    f = h5py.File(fn, 'r')
    td = f['FieldData/TD']
    ks = sorted(td.keys(), key=lambda s: int(s))
    x = np.array(f['Mesh/x']); y = np.array(f['Mesh/y']); z = np.array(f['Mesh/z'])
    frames, times = [], []
    for k in ks:
        a = np.array(td[k])                       # (3, nx, ny, nz)
        e = np.sqrt((a * a).sum(0))
        frames.append(e)
        times.append(float(td[k].attrs.get('time', 0.0)))
    return (x * 1e3, y * 1e3, z * 1e3), frames, np.array(times)

(xV, yV, zV), fV, tV = load(os.path.join(SIM, 'E_sliceV.h5'))
(xH, yH, zH), fH, tH = load(os.path.join(SIM, 'E_sliceH.h5'))
n = min(len(fV), len(fH))
fV, fH, tV = fV[:n], fH[:n], tV[:n]
print('frames V/H', len(fV), len(fH), 'n=', n)

# squeeze to 2D for each slice
SV = [np.squeeze(a) for a in fV]                  # (nz, nx) or (nx, nz)
SH = [np.squeeze(a) for a in fH]
# determine orientation
def orient(fr, x, y, z, name):
    if fr.shape == (len(z), len(x)):
        return fr, x, z
    if fr.shape == (len(x), len(z)):
        return fr.T, x, z
    # horizontal slice: expects (len(y), len(x)) or (len(x),len(y))
    if fr.shape == (len(y), len(x)):
        return fr, x, y
    if fr.shape == (len(x), len(y)):
        return fr.T, x, y
    raise ValueError('shape %s %s' % (name, fr.shape))

# global colour scale (percentile of all |E|, gamma-scaled)
allv = []
for a, b in zip(fV, fH):
    allv.append(a.ravel())
    allv.append(b.ravel())
allv = np.concatenate(allv)
vref = float(allv.max())
if vref <= 0: vref = 1.0
DMIN = -85.0
XLIM = (-35.0, 45.0)
ZLIM = (-30.0, 55.0)
YLIM = (-55.0, 30.0)
print('vref(global max)=%.3e' % vref)
norm = plt.Normalize(vmin=DMIN, vmax=0.0)
cmap = plt.get_cmap('turbo')

def todB(a):
    return 20.0 * np.log10(np.maximum(a, vref * 10.0 ** (DMIN / 20.0)) / vref)

txt = []
for i in range(n):
    aV, xVp, zVp = orient(SV[i], xV, yV, zV, 'V')
    aH, xHp, yHp = orient(SH[i], xH, yH, zH, 'H')
    aV = todB(aV); aH = todB(aH)
    fig, axs = plt.subplots(1, 2, figsize=(12.6, 4.6))
    a = axs[0]
    m = a.pcolormesh(xVp, zVp, aV, shading='auto', cmap=cmap, norm=norm)
    a.add_patch(plt.Rectangle((X0, -Z_BOARD), X1 - X0, 2 * Z_BOARD, fill=False,
                              ec='w', lw=1.4))
    a.plot([CX, CX], [Z_BOARD, Z_BOARD + 100], color='w', lw=2.0, alpha=0.8)
    a.set_xlabel('x [mm]'); a.set_ylabel('z [mm]')
    a.set_title('side view  (y=%.1f mm cut)' % yV[0])
    a.set_aspect('equal'); a.set_xlim(*XLIM); a.set_ylim(*ZLIM)
    b = axs[1]
    m2 = b.pcolormesh(xHp, yHp, aH, shading='auto', cmap=cmap, norm=norm)
    b.add_patch(plt.Rectangle((X0, Y0), X1 - X0, Y1 - Y0, fill=False, ec='w', lw=1.4))
    b.plot([CX, CX], [CY, CY + 100], color='w', lw=2.0, alpha=0.8)
    b.set_xlabel('x [mm]'); b.set_ylabel('y [mm]')
    b.set_title('top view  (z=+3 mm cut)')
    b.set_aspect('equal'); b.set_xlim(*XLIM); b.set_ylim(*YLIM)
    fig.colorbar(m, ax=axs, fraction=0.02, pad=0.01, label='|E| [dB rel. global max]')
    fig.suptitle('TPS62933 near-field |E| propagation  —  t = %.2f ns (frame %d/%d)' %
                 (tV[i] * 1e9, i + 1, n), fontsize=12)
    fn = os.path.join(FR, 'frame_%03d.png' % i)
    fig.savefig(fn, dpi=90, bbox_inches='tight'); plt.close(fig)
    txt.append(fn)
print('rendered', len(txt), 'frames')

gif = os.path.join(HERE, 'propagation.gif')
mp4 = os.path.join(HERE, 'propagation.mp4')
subprocess.run(['ffmpeg', '-y', '-framerate', '7', '-i', os.path.join(FR, 'frame_%03d.png'),
                '-vf', 'scale=1000:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse',
                gif], check=True, capture_output=True)
subprocess.run(['ffmpeg', '-y', '-framerate', '7', '-i', os.path.join(FR, 'frame_%03d.png'),
                '-pix_fmt', 'yuv420p', '-vf', 'scale=1000:-2', mp4], check=True, capture_output=True)
print('wrote', gif, os.path.getsize(gif), 'bytes')
print('wrote', mp4, os.path.getsize(mp4), 'bytes')

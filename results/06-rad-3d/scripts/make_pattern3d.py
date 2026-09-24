#!/usr/bin/env python3
"""3D far-field pattern (30/100/300 MHz) from rad/b1 far-field data.

Normalisation kept identical to rad/b1/post_b1.py:
    H(f,th,ph) = E_raw / |if_tot(f)|          [Ohm/m]
    E_3m       = H * |I_loop(f)| / 3          [V/m]
    I_loop     = exact harmonic spectrum of the real input hot-loop current
                 (Ipk=4.3 A, fsw=805 kHz, D=0.5, tr=tf=5 ns)
Pattern radius r is the *normalised* linear field (r = E/E_max), so only the
shape is meaningful -- absolute amplitude is NOT reliable (see REPORT.md).
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm

HERE = os.path.dirname(os.path.abspath(__file__))
RADB1 = os.path.abspath(os.path.join(HERE, '..', 'rad', 'b1'))

IPK, FSW, DUTY, TR = 4.3, 805e3, 0.5, 5e-9

def harmon_amp(freqs):
    T = 1.0 / FSW; N = 200001
    t = np.linspace(0, T, N, endpoint=False)
    top = max(DUTY * T - 2 * TR, 0.0)
    x = np.interp(t, [0, TR, TR + top, 2 * TR + top, T], [0, IPK, IPK, 0, 0])
    return np.abs(np.fft.rfft(x) / N)

fd = np.load(os.path.join(RADB1, 'nf2ff_b1.npz'))
pd = np.load(os.path.join(RADB1, 'port_b1.npz'))
freq = fd['freq']; theta = fd['theta']; phi = fd['phi']
E_raw = fd['E_norm']
if_tot = np.abs(pd['if_tot'])
H = E_raw / if_tot[:, None, None]

_SPEC = harmon_amp(None)

def harmonic_amp_at(f0):
    """|I_n| at the nearest harmonic k*FSW (k = round(f0/FSW))."""
    k = max(1, int(round(f0 / FSW)))
    return _SPEC[k]

targets = [30e6, 100e6, 300e6]
TH, PH = np.meshgrid(np.radians(theta), np.radians(phi), indexing='ij')

fig = plt.figure(figsize=(16, 5.6))
for n, ftgt in enumerate(targets):
    idx = int(np.argmin(np.abs(freq - ftgt)))
    f_act = freq[idx]
    I = harmonic_amp_at(ftgt)
    E3 = H[idx] * I / 3.0
    dB = 20 * np.log10(np.maximum(E3, 1e-30) / 1e-6)
    r = 10.0 ** ((dB - dB.max()) / 20.0)
    r = np.clip(r, 0.0, 1.0)
    Xr = r * np.sin(TH) * np.cos(PH)
    Yr = r * np.sin(TH) * np.sin(PH)
    Zr = r * np.cos(TH)
    ax = fig.add_subplot(1, 3, n + 1, projection='3d')
    norm = plt.Normalize(dB.min(), dB.max())
    fc = cm.jet(norm(dB))
    ax.plot_surface(Xr, Yr, Zr, facecolors=fc, rstride=1, cstride=1,
                    linewidth=0, antialiased=True, shade=False)
    # schematic board model at origin (scaled, not to scale)
    bx, by = 26.5, 50.5
    s = 0.45 / max(bx, by)
    ax.plot([-bx/2*s, bx/2*s, bx/2*s, -bx/2*s, -bx/2*s],
            [-by/2*s, -by/2*s, by/2*s, by/2*s, -by/2*s],
            [0, 0, 0, 0, 0], color='k', lw=1.2)
    ax.set_title('%.0f MHz  (E$_{max}$=%.0f dBµV/m @3m)' % (f_act / 1e6, dB.max()))
    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z (board normal)')
    ax.set_box_aspect((1, 1, 0.9)); ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(-1, 1)
    ax.set_xticks([-1, 0, 1]); ax.set_yticks([-1, 0, 1]); ax.set_zticks([-1, 0, 1])
    ax.view_init(elev=22, azim=-60)

sm = cm.ScalarMappable(norm=plt.Normalize(-60, 0), cmap='jet')
sm.set_array([])
cb = fig.colorbar(sm, ax=fig.axes, fraction=0.02, pad=0.01)
cb.set_label('normalised pattern [dB rel. max]')
fig.suptitle('TPS62933 board far-field 3D pattern (r ∝ normalised E) — board model schematic at origin', fontsize=12)
fig.savefig(os.path.join(HERE, 'fig_pattern_3d.png'), dpi=130, bbox_inches='tight')
print('saved fig_pattern_3d.png')
for ftgt in targets:
    idx = int(np.argmin(np.abs(freq - ftgt)))
    E3 = H[idx] * harmonic_amp_at(ftgt) / 3.0
    dB = 20 * np.log10(np.maximum(E3, 1e-30) / 1e-6)
    print('%.0f MHz (f=%.2f MHz): Emax=%.1f dBuV/m  Dmax=%.3f dBi' %
          (ftgt / 1e6, freq[idx] / 1e6, dB.max(), fd['Dmax'][idx]))

#!/usr/bin/env python3
"""B1 post-processing: far field (30MHz-1GHz) vs CISPR 32 Class B @3m.

Absolute far field via the transfer function:
    H(f,th,ph) = E_raw(f,th,ph) / |if_tot(f)|          [Ohm/m]
    E_3m(f)    = H(f) * |I_loop(f)| / 3                [V/m]
where I_loop(f) is the exact harmonic spectrum of the real input hot-loop
switching current (trapezoidal pulse train: Ipk, fsw, duty, tr=tf=5ns).

Normalisation: X_true = X_raw/|uf_inc|; the ratio E_raw/if_tot cancels the
common dt-scaling, so multiplying by the true loop-current amplitude gives a
genuine field (method cross-checked against the validated patch smoke test).
An analytic small-loop (magnetic-dipole) estimate is printed as an independent
order-of-magnitude sanity check.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

B1 = os.path.dirname(os.path.abspath(__file__))

# ---- real hot-loop switching current (from results_v2/E_full.txt) ----
IPK  = 4.3          # A   peak loop current (max i(l2), steady state) from E_full.txt
FSW  = 805e3        # Hz
DUTY = 0.5
TR   = 5e-9         # s   rise = fall (user given)

def harmonic_amplitudes(freqs):
    """|I_n| [A] of a trapezoidal pulse train (one period DFT, peak amplitude)."""
    T = 1.0 / FSW
    N = 200001
    t = np.linspace(0.0, T, N, endpoint=False)
    top = max(DUTY * T - 2 * TR, 0.0)
    x = np.interp(t, [0, TR, TR + top, 2 * TR + top, T], [0, IPK, IPK, 0, 0])
    return np.abs(np.fft.rfft(x) / N)

def nearest_harmonics(fmin, fmax):
    k = np.arange(1, int(fmax / FSW) + 1)
    f = k * FSW
    m = (f >= fmin) & (f <= fmax)
    return f[m], k[m]

fd = np.load(os.path.join(B1, 'nf2ff_b1.npz'))
pd = np.load(os.path.join(B1, 'port_b1.npz'))
freq = fd['freq']; theta = fd['theta']; phi = fd['phi']
E_raw = fd['E_norm']                              # (nf, nth, nph)
if_tot = np.abs(pd['if_tot']); uf_inc = np.abs(pd['uf_inc'])
Prad = fd['Prad']; Dmax = fd['Dmax']

H = E_raw / if_tot[:, None, None]                 # Ohm/m

fh, kh = nearest_harmonics(30e6, 1e9)
Ih = harmonic_amplitudes(fh)[kh]
Hmax = H.max(axis=(1, 2))
Hmax_h = np.interp(fh, freq, Hmax)
E_3m = Hmax_h * Ih / 3.0
E_dBuVm = 20 * np.log10(np.maximum(E_3m, 1e-30) / 1e-6)

iw = int(np.argmax(E_dBuVm)); fw = fh[iw]
idx = int(np.argmin(np.abs(freq - fw)))
Ew_raw = E_raw[idx] / if_tot[idx] * Ih[iw] / 3.0
Ew_dB = 20 * np.log10(np.maximum(Ew_raw, 1e-30) / 1e-6)

# ---- analytic small-loop (magnetic dipole) reference at the worst harmonic ----
eta = 376.73; c = 2.99792458e8
A_loop = 7.5e-6                                  # m^2  (5mm x 1.51mm)
k = 2 * np.pi * fw / c
E1m_loop = eta * k**2 * A_loop * Ih[iw] / (4 * np.pi)   # V/m at 1 m
print('worst harmonic         : %.1f MHz  E(max)=%.1f dBuV/m  |I|=%.3e A' % (fw/1e6, E_dBuVm[iw], Ih[iw]))
print('analytic loop est @3m  : %.1f dBuV/m  (E1m=%.3e V/m)' % (20*np.log10(E1m_loop/3/1e-6), E1m_loop))
print('E range over band      : %.1f .. %.1f dBuV/m' % (E_dBuVm.min(), E_dBuVm.max()))
print('Dmax range             : %.3f .. %.3f dBi' % (Dmax.min(), Dmax.max()))
print('|if_tot| range (raw)   : %.2e .. %.2e' % (if_tot.min(), if_tot.max()))
# ports currents (true A) for reference
print('true sim loop current @100MHz = %.4g A' % (if_tot[idx]/uf_inc[idx]))

# ---- CISPR 32 Class B @3m (30-230MHz:30, 230M-1G:37 dBuV/m @10m ; +10 dB @3m) ----
lim_f = np.array([30e6, 230e6, 230e6, 1e9])
lim_v = np.array([40.0, 40.0, 47.0, 47.0])

# =============== figure 1: far-field pattern / hotspot ===============
TH, PH = np.meshgrid(np.radians(theta), np.radians(phi), indexing='ij')
R = np.maximum(Ew_dB - Ew_dB.max() + 60, 0.0)     # dB scale, 60 dB dynamic
Xr = R * np.sin(TH) * np.cos(PH); Yr = R * np.sin(TH) * np.sin(PH); Zr = R * np.cos(TH)

fig = plt.figure(figsize=(12, 5.2))
ax = fig.add_subplot(1, 2, 1, projection='3d')
ax.plot_surface(Xr, Yr, Zr, rstride=1, cstride=1, cmap='jet',
                linewidth=0, antialiased=True, alpha=0.95)
ax.set_title('Far-field |E| pattern @ %.0f MHz' % (fw/1e6))
ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
ax.set_box_aspect((1, 1, 0.6))
ax2 = fig.add_subplot(1, 2, 2)
pc = ax2.pcolormesh(phi, theta, Ew_dB, shading='auto', cmap='turbo')
ax2.set_xlabel('phi [deg]'); ax2.set_ylabel('theta [deg] (0=board normal +z)')
ax2.set_title('E dBµV/m @3m, %.0f MHz' % (fw/1e6))
fig.colorbar(pc, ax=ax2, label='dBµV/m')
fig.tight_layout()
fig.savefig(os.path.join(B1, 'fig_farfield_3d.png'), dpi=130)
print('saved fig_farfield_3d.png')

# =============== figure 2: vs CISPR 32 Class B @3m ===============
fig2, ax = plt.subplots(figsize=(9, 5.5))
ax.semilogx(fh/1e6, E_dBuVm, lw=0.9, color='C0', label='B1 model (max over angles)')
ax.semilogx(lim_f/1e6, lim_v, 'r--', lw=2, label='CISPR 32 Class B @3m (30/37 @10m +10dB)')
ax.set_xlabel('Frequency [MHz]'); ax.set_ylabel('E-field @3 m [dBµV/m]')
ax.set_title('TPS62933 B1 board-level far field vs CISPR 32 Class B (3 m)')
ax.grid(True, which='both', alpha=0.3); ax.set_ylim(0, 90)
ax.legend(loc='lower left')
fig2.tight_layout()
fig2.savefig(os.path.join(B1, 'fig_farfield_vs_cispr.png'), dpi=130)
print('saved fig_farfield_vs_cispr.png')

np.savez(os.path.join(B1, 'farfield_b1.npz'), fh=fh, Ih=Ih, E_dBuVm=E_dBuVm,
         fw=fw, theta=theta, phi=phi, Ew_dB=Ew_dB, Hmax=Hmax, freq=freq,
         if_tot=if_tot, uf_inc=uf_inc, Prad=Prad, Dmax=Dmax)
print('done post')

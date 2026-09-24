#!/usr/bin/env python3
"""Three-variant radiated-emission comparison (FullCopper / TopCutout / DualCutout)."""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
VERS = ['full', 'top', 'dual']
LABEL = {'full': 'FullCopper', 'top': 'TopCutout', 'dual': 'DualCutout'}
COL = {'full': 'C0', 'top': 'C1', 'dual': 'C3'}
IPK, FSW, DUTY, TR = 4.3, 805e3, 0.5, 5e-9

def harmonics(fmin, fmax):
    k = np.arange(1, int(fmax/FSW)+1); f = k*FSW
    m = (f >= fmin) & (f <= fmax); return f[m], k[m]

def wave():
    T = 1.0/FSW; N = 200001
    t = np.linspace(0, T, N, endpoint=False)
    top = max(DUTY*T - 2*TR, 0.0)
    x = np.interp(t, [0, TR, TR+top, 2*TR+top, T], [0, IPK, IPK, 0, 0])
    return np.abs(np.fft.rfft(x)/N)

fh, kh = harmonics(30e6, 1e9)
Ih = wave()[kh]

res = {}
for v in VERS:
    fd = np.load(os.path.join(HERE, 'ver_'+v, 'nf2ff.npz'))
    pd = np.load(os.path.join(HERE, 'ver_'+v, 'port.npz'))
    freq = fd['freq']; E = fd['E_norm']; ift = np.abs(pd['if_tot'])
    H = E / ift[:, None, None]
    Hmax = H.max(axis=(1, 2))
    dB = 20*np.log10(np.maximum(np.interp(fh, freq, Hmax)*Ih/3.0, 1e-30)/1e-6)
    res[v] = dict(freq=freq, th=fd['theta'], ph=fd['phi'], E=E, ift=ift, dB=dB)

def band_val(dB, t, half=2):
    """max over the +/-half harmonics (avoids the exact D=0.5 even-harmonic nulls)."""
    i = int(np.argmin(np.abs(fh - t)))
    return dB[max(0, i-half):i+half+1].max(), i

targets = [30e6, 100e6, 300e6]
print('%-11s %8s %8s %8s' % ('variant', *[('%.0fMHz' % (t/1e6)) for t in targets]))
vals = {}
for v in VERS:
    vals[v] = [band_val(res[v]['dB'], t) for t in targets]
    print('%-11s' % LABEL[v], ' '.join('%8.1f' % x[0] for x in vals[v]))
print('\nDelta vs FullCopper (dB):')
for v in ['top', 'dual']:
    print('%-11s' % LABEL[v], ' '.join('%8.2f' % (vals[v][k][0]-vals['full'][k][0])
                                       for k in range(len(targets))))

# =============== fig 1 ===============
fig, ax = plt.subplots(figsize=(9.5, 5.5))
for v in VERS:
    ax.semilogx(fh/1e6, res[v]['dB'], lw=0.9, color=COL[v], label=LABEL[v])
ax.semilogx([30, 230, 230, 1000], [40, 40, 47, 47], 'k--', lw=1.6, label='CISPR 32 B @3m')
ax.set_xlabel('Frequency [MHz]'); ax.set_ylabel('E @3m [dBµV/m]')
ax.set_title('Radiated far field — three pour variants (identical source/mesh)')
ax.grid(True, which='both', alpha=0.3); ax.set_ylim(0, 90); ax.legend(loc='lower left')
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'fig_rad_compare_freq.png'), dpi=130)
print('saved fig_rad_compare_freq.png')

# =============== fig 2: hotspot at max-spread harmonic (50-400 MHz) ===============
k = np.arange(len(fh))[(fh > 50e6) & (fh < 400e6)]
spread = np.array([max(res[v]['dB'][i] for v in VERS) - min(res[v]['dB'][i] for v in VERS) for i in k])
iw = int(k[int(np.argmax(spread))]); fw = fh[iw]
print('hotspot figure at %.1f MHz (max spread %.2f dB)' % (fw/1e6, spread.max()))

maps = {}
for v in VERS:
    j = int(np.argmin(np.abs(res[v]['freq']-fw)))
    maps[v] = 20*np.log10(np.maximum(res[v]['E'][j]/res[v]['ift'][j]*Ih[iw]/3.0, 1e-30)/1e-6)
vmin = min(m.min() for m in maps.values()); vmax = max(m.max() for m in maps.values())
fig2, axs = plt.subplots(1, 3, figsize=(15, 4.6))
for a, v in zip(axs, VERS):
    pc = a.pcolormesh(res[v]['ph'], res[v]['th'], maps[v], shading='auto',
                      cmap='turbo', vmin=vmin, vmax=vmax)
    a.set_title('%s  %.0f MHz' % (LABEL[v], fw/1e6))
    a.set_xlabel('phi [deg]'); a.set_ylabel('theta [deg]')
    fig2.colorbar(pc, ax=a, label='dBµV/m')
fig2.tight_layout(); fig2.savefig(os.path.join(HERE, 'fig_rad_compare_hotspot.png'), dpi=130)
print('saved fig_rad_compare_hotspot.png')

np.savez(os.path.join(HERE, 'compare.npz'), fh=fh, Ih=Ih, fw=fw,
         **{('dB_'+v): res[v]['dB'] for v in VERS})
print('done compare')

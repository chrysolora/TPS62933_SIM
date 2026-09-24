#!/usr/bin/env python3
"""B2 post-processing: three load cases x three pour variants, CM/DM split.

Physical calibration (all magnitudes flagged as estimates where noted):
  * DM  : hot-loop trapezoid current, fsw=805 kHz, D=0.5, tr=tf=5 ns.
          Ipk scaled with load: 3.0A->4.3A (B1), 1.5A->2.15A, 0A->0.43A
          (light-load/PFM estimate, not 0 -- switching still present).
  * CM  : upstream LM50-20B24 ripple trapezoid, fsw=65 kHz, D=0.5,
          tr=tf in {50,100} ns (ASSUMPTION, datasheet silent), Vpp in
          {150 mV (worst), 100 mV (user)}.
          I_cm(f) = Vharm(f) / |Z_loop(f)|,
          Z_loop = R_cable + jw L_cable + 1/(jw Cy) + 1/(jw Cp)
          L_cable=105 nH (0.75mm^2, 10cm free-space formula), R_cable~0.05 ohm,
          Cy=2.2 nF (assumed), Cp in {2,10,50} pF (assumed).
  * E_3m = H(f) * I(f) / 3   (B1 normalisation; dt factor cancels).
"""
import os, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
VERS = ['full', 'top', 'dual']
LABEL = {'full': 'FullCopper', 'top': 'TopCutout', 'dual': 'DualCutout'}
COL = {'full': 'C0', 'top': 'C1', 'dual': 'C3'}
LOADS = [('0A', 0.43), ('1.5A', 2.15), ('3.0A', 4.3)]

DM_FSW, DM_D, DM_TR = 805e3, 0.5, 5e-9
UP_FSW, UP_D = 65e3, 0.5
L_CABLE, R_CABLE, C_Y = 105e-9, 0.05, 2.2e-9

def harm(fmin, fmax, fsw, A, D, tr):
    N = 400001
    t = np.linspace(0, 1.0/fsw, N, endpoint=False)
    T = 1.0/fsw; top = max(D*T - 2*tr, 0.0)
    x = np.interp(t, [0, tr, tr+top, 2*tr+top, T], [0, A, A, 0, 0])
    X = np.abs(np.fft.rfft(x)/N)
    k = np.arange(1, int(fmax/fsw)+1); f = k*fsw
    m = (f >= fmin) & (f <= fmax)
    return f[m], X[k[m]]

def zloop(f, Cp):
    w = 2*np.pi*f
    return R_CABLE + 1j*w*L_CABLE + 1.0/(1j*w*C_Y) + 1.0/(1j*w*Cp)

def interp_cplx(H, fsrc, fdst):
    """H (nf,nth,nph) complex on fsrc -> fdst, vectorised over angles."""
    nf, a, b = H.shape
    Hf = H.reshape(nf, -1)
    out = np.empty((len(fdst), a*b), complex)
    for c in range(a*b):
        out[:, c] = np.interp(fdst, fsrc, Hf[:, c].real) + \
                    1j*np.interp(fdst, fsrc, Hf[:, c].imag)
    return out.reshape(len(fdst), a, b)

# ---------- load & precompute ----------
fh_dm, Ih_dm_unit = harm(30e6, 1e9, DM_FSW, 1.0, DM_D, DM_TR)
fh_cm, Vh_cm_unit = harm(30e6, 1e9, UP_FSW, 1.0, UP_D, 50e-9)
FAP = fh_dm                                  # common output grid (DM harmonics)
Vh_cm_on_dm = np.interp(FAP, fh_cm, Vh_cm_unit)

data = {}
for v in VERS:
    dd = np.load(os.path.join(HERE, 'ver_%s_dm' % v, 'nf2ff.npz'))
    pd = np.load(os.path.join(HERE, 'ver_%s_dm' % v, 'port.npz'))
    dc = np.load(os.path.join(HERE, 'ver_%s_cm' % v, 'nf2ff.npz'))
    pc = np.load(os.path.join(HERE, 'ver_%s_cm' % v, 'port.npz'))
    f0 = dd['freq']
    Hdm = dd['E_norm'] / np.abs(pd['dm_if_tot'])[:, None, None]
    Hcm = dc['E_norm'] / np.abs(pc['cm_if_tot'])[:, None, None]
    data[v] = dict(th=dd['theta'], ph=dd['phi'],
                   Hdm=interp_cplx(Hdm, f0, FAP),
                   Hcm=interp_cplx(interp_cplx(Hcm, dc['freq'], FAP), FAP, FAP))
print('precomputed %d variants on %d-pt grid' % (len(VERS), len(FAP)))

def build(v, Ipk, Vpp, tr, Cp):
    Ih = Ih_dm_unit*Ipk
    Vh = np.interp(FAP, *harm(30e6, 1e9, UP_FSW, Vpp/2.0, UP_D, tr))
    Icm = Vh/np.abs(zloop(FAP, Cp))
    Edm = data[v]['Hdm']*Ih[:, None, None]/3.0
    Ecm = data[v]['Hcm']*Icm[:, None, None]/3.0
    return Edm, Ecm

def dbmax(E):
    return 20*np.log10(np.maximum(np.abs(E).max(axis=(1, 2)), 1e-30)/1e-6)

def band(dB, t, half=2):
    i = int(np.argmin(np.abs(FAP - t)))
    return dB[max(0, i-half):i+half+1].max()

lines = []
def P(s=''):
    print(s); lines.append(s)
def fmt(x): return '%7.1f' % x if np.isfinite(x) else '   -inf'

P('=== B2 numbers ===')
P('upstream LM50-20B24: fsw=65kHz, ripple 150mVpp(worst)/100mVpp(user)')
P('cable 10cm/0.75mm^2: L=105nH, R~0.05ohm ; Cy=2.2nF ; Cp=2/10/50pF')
P('30MHz = 461st harmonic of 65kHz')
P('')
T = [30e6, 100e6, 300e6]

P('--- E@3m [dBuV/m]  (DM+CM, 150mV, tr=50ns, Cp=10pF) ---')
for t in T:
    P('  @%.0f MHz' % (t/1e6))
    for lab, Ipk in LOADS:
        row = []
        for v in VERS:
            Edm, Ecm = build(v, Ipk, 150e-3, 50e-9, 10e-12)
            row.append(fmt(band(dbmax(Edm+Ecm), t)))
        P('    %-5s Full/Top/Dual: %s' % (lab, ' '.join(row)))

P('')
P('--- DM vs CM vs total, FullCopper, 3.0A (Cp=10pF,150mV,50ns) ---')
Edm, Ecm = build('full', 4.3, 150e-3, 50e-9, 10e-12)
for t in T:
    P('  @%.0f MHz  DM=%s  CM=%s  tot=%s'
      % (t/1e6, fmt(band(dbmax(Edm), t)), fmt(band(dbmax(Ecm), t)),
         fmt(band(dbmax(Edm+Ecm), t))))

P('')
P('--- with / without upstream ripple (FullCopper,3.0A,150mV,50ns,10pF) ---')
for t in T:
    a, b = band(dbmax(Edm), t), band(dbmax(Edm+Ecm), t)
    P('  @%.0f MHz without=%s with=%s delta=%5.2f dB' % (t/1e6, fmt(a), fmt(b), b-a))

P('')
P('--- Cp sensitivity (FullCopper,3.0A,150mV,50ns) CM-only @3m ---')
for Cp in [2e-12, 10e-12, 50e-12]:
    _, Ecm2 = build('full', 4.3, 150e-3, 50e-9, Cp)
    P('  Cp=%4.0f pF : %s' % (Cp*1e12, ' '.join(fmt(band(dbmax(Ecm2), t)) for t in T)))

P('')
P('--- upstream tr sensitivity (FullCopper,3.0A,150mV,Cp=10pF) CM-only ---')
for tr in [50e-9, 100e-9]:
    _, Ecm2 = build('full', 4.3, 150e-3, tr, 10e-12)
    P('  tr=%3.0f ns : %s' % (tr*1e9, ' '.join(fmt(band(dbmax(Ecm2), t)) for t in T)))

P('')
P('--- ripple amplitude 100 vs 150 mV (FullCopper,3.0A,50ns,10pF) CM-only ---')
for Vpp in [100e-3, 150e-3]:
    _, Ecm2 = build('full', 4.3, Vpp, 50e-9, 10e-12)
    P('  %3.0f mV : %s' % (Vpp*1e3, ' '.join(fmt(band(dbmax(Ecm2), t)) for t in T)))

open(os.path.join(HERE, 'numbers_b2.txt'), 'w').write('\n'.join(lines) + '\n')

# ================= figures =================
LIM = lambda f: np.where(f < 230e6, 40.0, 47.0)

fig, axs = plt.subplots(1, 3, figsize=(16, 5), sharey=True)
for ax, (lab, Ipk) in zip(axs, LOADS):
    for v in VERS:
        Edm, Ecm = build(v, Ipk, 150e-3, 50e-9, 10e-12)
        ax.semilogx(FAP/1e6, dbmax(Edm+Ecm), lw=0.9, color=COL[v], label=LABEL[v])
    ax.semilogx(FAP/1e6, LIM(FAP), 'k--', lw=1.6, label='CISPR 32 B @3m')
    ax.set_title('load %s' % lab); ax.set_xlabel('Freq [MHz]')
    ax.grid(True, which='both', alpha=0.3); ax.set_ylim(0, 90)
axs[0].set_ylabel('E @3m [dBuV/m]'); axs[0].legend(loc='lower left', fontsize=8)
fig.suptitle('B2 radiated far field (DM + upstream CM) - three pour variants, three loads')
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'fig_rad_compare_freq_b2.png'), dpi=120)
print('saved fig_rad_compare_freq_b2.png')

fig2, ax = plt.subplots(figsize=(9.5, 5.5))
ax.semilogx(FAP/1e6, dbmax(Edm), 'C0', lw=1.1, label='DM (board switching)')
ax.semilogx(FAP/1e6, dbmax(Ecm), 'C3', lw=1.1, label='CM (upstream ripple via cable+Cp)')
ax.semilogx(FAP/1e6, dbmax(Edm+Ecm), 'k', lw=1.3, label='total')
ax.semilogx(FAP/1e6, LIM(FAP), 'k--', lw=1.4, label='CISPR 32 B @3m')
ax.set_xlabel('Freq [MHz]'); ax.set_ylabel('E @3m [dBuV/m]')
ax.set_title('CM vs DM breakdown (FullCopper, 3.0A, 150mV, Cp=10pF)')
ax.grid(True, which='both', alpha=0.3); ax.set_ylim(0, 90); ax.legend(loc='lower left')
fig2.tight_layout(); fig2.savefig(os.path.join(HERE, 'fig_rad_cm_dm_breakdown.png'), dpi=120)
print('saved fig_rad_cm_dm_breakdown.png')

# hotspot
tot = {v: dbmax(build(v, 4.3, 150e-3, 50e-9, 10e-12)[0] + build(v, 4.3, 150e-3, 50e-9, 10e-12)[1])
       for v in VERS}
kk = np.arange(len(FAP))[(FAP > 50e6) & (FAP < 400e6)]
spread = np.array([max(tot[v][i] for v in VERS) - min(tot[v][i] for v in VERS) for i in kk])
iw = int(kk[int(np.argmax(spread))]); fw = FAP[iw]
maps = {}
for v in VERS:
    Edm, Ecm = build(v, 4.3, 150e-3, 50e-9, 10e-12)
    maps[v] = 20*np.log10(np.maximum(np.abs(Edm+Ecm)[iw], 1e-30)/1e-6)
vmin = min(m.min() for m in maps.values()); vmax = max(m.max() for m in maps.values())
fig3, axs = plt.subplots(1, 3, figsize=(15, 4.6))
for a, v in zip(axs, VERS):
    pc = a.pcolormesh(data[v]['ph'], data[v]['th'], maps[v], shading='auto',
                      cmap='turbo', vmin=vmin, vmax=vmax)
    a.set_title('%s  %.0f MHz' % (LABEL[v], fw/1e6))
    a.set_xlabel('phi [deg]'); a.set_ylabel('theta [deg]')
    fig3.colorbar(pc, ax=a, label='dBuV/m')
fig3.suptitle('B2 far-field hotspot (3.0A, total), spread %.2f dB' % spread.max())
fig3.tight_layout(); fig3.savefig(os.path.join(HERE, 'fig_rad_compare_hotspot_b2.png'), dpi=120)
print('saved fig_rad_compare_hotspot_b2.png (%.0f MHz)' % (fw/1e6))

fig4, ax = plt.subplots(figsize=(9.5, 5.5))
ax.semilogx(FAP/1e6, dbmax(Edm), 'C0', lw=1.1, label='without upstream ripple (DM only)')
ax.semilogx(FAP/1e6, dbmax(Edm+Ecm), 'C3', lw=1.2, label='with upstream ripple (DM+CM)')
ax.semilogx(FAP/1e6, LIM(FAP), 'k--', lw=1.4, label='CISPR 32 B @3m')
ax.set_xlabel('Freq [MHz]'); ax.set_ylabel('E @3m [dBuV/m]')
ax.set_title('Effect of upstream ripple source (FullCopper, 3.0A, Cp=10pF)')
ax.grid(True, which='both', alpha=0.3); ax.set_ylim(0, 90); ax.legend(loc='lower left')
fig4.tight_layout(); fig4.savefig(os.path.join(HERE, 'fig_rad_with_vs_without_input_ripple.png'), dpi=120)
print('saved fig_rad_with_vs_without_input_ripple.png')

np.savez(os.path.join(HERE, 'compare_b2.npz'), fh_dm=FAP, fw=fw)
print('done post_b2')

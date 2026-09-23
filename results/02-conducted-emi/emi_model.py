#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
emi_model.py -- TPS62933 buck : frequency-domain DM conducted-EMI (CISPR 32 approx.)
                + analytic soft-start curve.

Run:  cd /mnt/raid10/sim-work/tps62933/emi && python3 emi_model.py
Out :  fig_emi_cispr.png , fig_startup.png , emi_numbers.txt

--------------------------------------------------------------------------------
WHY THIS METHOD (not a re-run of the transient)
--------------------------------------------------------------------------------
The board input filter is  3x10uF || (L1=1uH) || 20uF+100nF (+ RC damping) .
The available (unencrypted) transient model cannot start with >30uF at the VIN
pin, and results_v2/ transient is NOT at a settled operating point (Vout still
drifting 10.6->11.5 V over 80 us, duty not locked).  Therefore the source is
taken ANALYTICALLY (trapezoidal switch-side input current) and the filter is
solved in the FREQUENCY DOMAIN (linear nodal analysis).  The data-derived
spectrum (FFT of chopping i(L2) with v(SW)) is computed as an independent,
documented cross-check only.

DM source = pulsating current drawn at U21.VIN  (high-side switch current).
Transimpedance Z_t(f) = V_LISN / I_src  from a linear network:
    node A = LISN port (= PPHV, fuse F1 = ideal short)   [30uF bank, 50ohm/1uF, 50uH]
    node B = PPHV_OUT_FILTER = U21.VIN (source node)     [20uF, 100nF, R85+20uF]
    L1 = 1uH between A and B
V_lisn_harmonic = I_src(k*fsw) * |Z_t(k*fsw)| ;  dBuV = 20log10(V/1uV).

ALL parasitics (ESR/ESL/DCR/tr) are ENGINEERING ESTIMATES -- see REPORT_emi.md.
Unknowns (C1 value, TVS off-cap, fuse R/L, datasheet Iss) are NOT invented.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RES = '/mnt/raid10/sim-work/tps62933/results_v2/'
OUT = '/mnt/raid10/sim-work/tps62933/emi/'
FSW = 805e3          # measured switching frequency (Hz), confirmed from results_v2
VIN = 24.0
VOUT = 12.0

# ----------------------------------------------------------------------------
# helper: complex impedances
# ----------------------------------------------------------------------------
def Zc(f, C): return 1.0 / (1j * 2 * np.pi * f * C)
def Zl(f, L): return 1j * 2 * np.pi * f * L

# ----------------------------------------------------------------------------
# 1. SOURCE  : trapezoidal switch-side input current -> harmonic amplitudes
# ----------------------------------------------------------------------------
def trapezoid_harmonics(fsw, ipk, duty, tr, tf, kmax, nper=32, ppp=2048):
    """|single-sided amplitude| of harmonics k*fsw for an ideal trapezoidal
    pulse train (DC removed).  Numeric -> no formula-scale errors."""
    T = 1.0 / fsw
    N = nper * ppp
    dt = (nper * T) / N
    t = np.arange(N) * dt
    tt = t % T
    on = duty * T
    x = np.zeros_like(t)
    # rising
    r = tt < tr
    x[r] = ipk * tt[r] / tr
    # flat top
    fl = (tt >= tr) & (tt < on)
    x[fl] = ipk
    # falling
    fa = (tt >= on) & (tt < on + tf)
    x[fa] = ipk * (1 - (tt[fa] - on) / tf)
    x -= x.mean()
    X = np.fft.rfft(x * np.hanning(N))
    f = np.fft.rfftfreq(N, dt)
    amp = 2 * np.abs(X) / np.hanning(N).sum()
    ks = np.arange(1, kmax + 1)
    fh = ks * fsw
    idx = np.array([int(round(fr / f[1])) for fr in fh])
    return fh, amp[idx]

# ----------------------------------------------------------------------------
# 2. FILTER + LISN : |Z_t(f)| = |V_LISN / I_src| via 2-node nodal solve
# ----------------------------------------------------------------------------
def transimpedance(f, par):
    f = np.atleast_1d(np.asarray(f, float))
    Zt = np.zeros(len(f), dtype=complex)
    for i, fr in enumerate(f):
        def capZ(C, esr, esl):
            return esr + Zl(fr, esl) + Zc(fr, C)
        # ---- node A : LISN port / PPHV -----------------------------------
        yA = (
            3 / capZ(10e-6, par['esr_10u'], par['esl_10u'])      # C17||C62||C63
            + 1.0 / (Zc(fr, 1e-6) + par['r_lisn'])               # LISN 50ohm via 1uF
            + 1.0 / Zl(fr, 50e-6)                                # LISN 50uH -> mains
        )
        # ---- L1 branch A<->B ---------------------------------------------
        Zl1 = Zl(fr, par['L1']) + par['dcr_L1']
        yl1 = 1.0 / Zl1
        if par['cp_L1'] > 0:
            yl1 += 1.0 / Zc(fr, par['cp_L1'])
        yA += yl1
        # ---- node B : PPHV_OUT_FILTER = U21.VIN (source node) ------------
        yB = (
            2 / capZ(10e-6, par['esr_10u'], par['esl_10u'])      # C73||C74
            + 1.0 / capZ(100e-9, par['esr_100n'], par['esl_100n'])  # C72
            + 1.0 / (par['r85'] + par['esr_10u'] / 2 + Zc(fr, 20e-6))  # R85 + C16||C64
            + yl1
        )
        Y = np.array([[yA, -yl1], [-yl1, yB]], complex)
        V = np.linalg.solve(Y, np.array([0.0, 1.0], complex))   # unit current at B
        Zt[i] = V[0]
    return np.abs(Zt)

PAR = dict(L1=1e-6, dcr_L1=30e-3, cp_L1=1.0e-12,
           esr_10u=3e-3,  esl_10u=1.2e-9,
           esr_100n=20e-3, esl_100n=1.0e-9,
           r85=100e-3, r_lisn=50.0)

# ----------------------------------------------------------------------------
# 3. CISPR 32 limits (QP, dBuV)
# ----------------------------------------------------------------------------
def cispr_B(f):
    f = np.asarray(f, float)
    return np.where(f < 0.5e6,
                    66.0 - 10.0 * (np.log10(f) - np.log10(0.15e6)) /
                    (np.log10(0.5e6) - np.log10(0.15e6)),
                    np.where(f < 5e6, 56.0, 60.0))
def cispr_A(f):
    f = np.asarray(f, float)
    return np.where(f < 0.5e6,
                    79.0 - 6.0 * (np.log10(f) - np.log10(0.15e6)) /
                    (np.log10(0.5e6) - np.log10(0.15e6)),
                    np.where(f < 5e6, 73.0, 73.0))

# ----------------------------------------------------------------------------
# 4. cases
# ----------------------------------------------------------------------------
# Ipk = inductor DC current (= load in CCM); D = Vout/Vin in CCM (losses -> ~0.5)
cases = [
    ('noload', 0.02, 0.22, 0.15, '#1f77b4'),
    ('half',   1.50, 1.50, 0.50, '#2ca02c'),
    ('full',   3.00, 3.00, 0.50, '#d62728'),
]
TR = 15e-9     # switching rise/fall time (estimate)

def load_transient(fname):
    d = np.loadtxt(RES + fname)
    return dict(t=d[:, 0], vout=d[:, 1], vsw=d[:, 3], il2=d[:, 7])

def data_spectrum(fname, kmax):
    """Cross-check only: FFT of chopping i(L2) by v(SW)."""
    d = load_transient(fname)
    t, vsw, il2 = d['t'], d['vsw'], d['il2']
    T = 1.0 / FSW
    K, PPP = 16, 1024
    N = K * PPP
    dt = K * T / N
    t0 = t[-1] - K * T
    tu = t0 + np.arange(N) * dt
    iu = np.interp(tu, t, il2)
    vu = np.interp(tu, t, vsw)
    x = np.where(vu > VIN / 2, iu, 0.0)
    x -= x.mean()
    X = np.fft.rfft(x * np.hanning(N))
    f = np.fft.rfftfreq(N, dt)
    amp = 2 * np.abs(X) / np.hanning(N).sum()
    ks = np.arange(1, kmax + 1)
    idx = np.array([int(round(k * FSW / f[1])) for k in ks])
    return ks * FSW, amp[idx]

kmax = int(30e6 / FSW)
lines = []
def log(s=''):
    print(s); lines.append(s)

log('=' * 76)
log('TPS62933 DM conducted-EMI model  (fsw = %.1f kHz, engine = trapezoid)' % (FSW / 1e3))
log('=' * 76)

results = {}
data_chk = {}
for name, iload, ipk, duty, col in cases:
    fh, ah = trapezoid_harmonics(FSW, ipk, duty, TR, TR, kmax)
    Zt = transimpedance(fh, PAR)
    V = ah * Zt
    dbuv = 20 * np.log10(np.maximum(V, 1e-30) / 1e-6)
    lim = cispr_B(fh)
    results[name] = dict(f=fh, dbuv=dbuv, ah=ah, Zt=Zt, lim=lim,
                         marg=lim - dbuv, col=col, iload=iload, ipk=ipk, duty=duty)
data_chk['full'] = data_spectrum('E_full.txt', kmax)
data_chk['half'] = data_spectrum('E_half.txt', kmax)
data_chk['noload'] = data_spectrum('E_noload.txt', kmax)

log('')
log('%-7s %-6s %-5s %8s %10s %10s %9s' %
    ('case', 'Iload', 'D', 'I_src(1)', 'Zt(1)', 'V(1)dBuV', 'margin(1)'))
for name, iload, ipk, duty, col in cases:
    r = results[name]
    log('%-7s %-6.2f %-5.2f %8.3f %10.3g %10.1f %9.1f' %
        (name, iload, duty, r['ah'][0], r['Zt'][0], r['dbuv'][0], r['marg'][0]))

log('')
log('Worst-case margin vs CISPR 32 Class B (min over 150k-30M):')
for name, iload, ipk, duty, col in cases:
    r = results[name]
    j = np.argmin(r['marg'])
    log('  %-7s : margin_min = %6.1f dB @ %.3f MHz  (V=%.1f dBuV, limit=%.1f)'
        % (name, r['marg'][j], r['f'][j] / 1e6, r['dbuv'][j], r['lim'][j]))
    bad = r['marg'] < 0
    log('            violations = %d harmonics' % bad.sum())

# ---- resonance / transfer features ----------------------------------------
fg = np.logspace(np.log10(150e3), np.log10(30e6), 4000)
Ztg = transimpedance(fg, PAR)
log('')
log('Filter+LISN |Z_t(f)| :  value @0.805MHz=%.3g ohm, @5MHz=%.3g, @30MHz=%.3g ohm'
    % (np.interp(805e3, fg, Ztg), np.interp(5e6, fg, Ztg), np.interp(30e6, fg, Ztg)))
Cser = 1.0 / (1.0 / (20e-6 + 0.1e-6 + 20e-6) + 1.0 / 30e-6)
fres = 1.0 / (2 * np.pi * np.sqrt(PAR['L1'] * Cser))
log('L1 * C_series resonance :  C_ser=%.2f uF (30uF || 40.1uF) -> f=%.1f kHz  (below CISPR band)'
    % (Cser * 1e6, fres / 1e3))
for C, esl, nm in [(30e-6, 1.2e-9/3, '30uF bank SRF'), (40e-6, 1.2e-9/2, '40uF bank SRF'),
                   (100e-9, 1.0e-9, '100nF SRF')]:
    log('  %-16s : f_SRF ~ %.2f MHz' % (nm, 1 / (2 * np.pi * np.sqrt(esl * C)) / 1e6))

# ---- data cross-check ------------------------------------------------------
log('')
log('Cross-check (data-derived chopped-i(L2) FFT, NOT settled -> indicative only):')
for k in ['full', 'half']:
    fd, ad = data_chk[k]
    r = results[k]
    ratio = 20 * np.log10(np.maximum(ad, 1e-12) / np.maximum(r['ah'], 1e-12))
    log('  %-5s : median(data/model) = %+.1f dB, spread = %.1f..%.1f dB'
        % (k, np.median(ratio), np.percentile(ratio, 5), np.percentile(ratio, 95)))

# ---- rough CM estimate (order-of-magnitude, parasitic Cp ASSUMED few pF) -----
# Authorised rough estimate (task sect 4, trap#4).  CM loop: SW-node switching
# voltage couples to earth via stray capacitance Cp, returns through the LISN
# (two 50 ohm in parallel = 25 ohm for CM).  Cp value is an ASSUMPTION.
Cp = 5e-12
fsw_v, vsw_h = trapezoid_harmonics(FSW, VIN, 0.5, TR, TR, kmax)   # SW voltage harmonics
log('')
log('ROUGH CM order-of-magnitude  (Cp = %.0f pF  ASSUMED ; CM LISN Z ~ 25 ohm):' % (Cp * 1e12))
log('   %-4s %-9s %-11s %-13s %-8s %s' % ('h', 'f[MHz]', 'V_sw_h[V]', 'I_cm[uA]', 'V_cm[dBuV]', 'CISPR-B'))
for k in [1, 3, 5, 10, 20]:
    fr = k * FSW
    Icm = 2 * np.pi * fr * Cp * vsw_h[k - 1]
    dbv = 20 * np.log10(Icm * 25 / 1e-6)
    log('   %-4d %-9.2f %-11.2f %-13.1f %-8.1f %.1f'
        % (k, fr / 1e6, vsw_h[k - 1], Icm * 1e6, dbv, cispr_B(np.array([fr]))[0]))
log('   -> CM is strongly Cp-dependent; at Cp=5pF it exceeds Class B by ~20-45 dB,')
log('      i.e. CM(not DM) is expected to be the binding constraint on this design.')

# ----------------------------------------------------------------------------
# FIG 1 : EMI spectrum vs CISPR
# ----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(11.5, 6.8))
fgp = np.logspace(np.log10(1.5e5), np.log10(3e7), 800)
ax.plot(fgp / 1e6, cispr_B(fgp), 'k-', lw=2.2, label='CISPR 32 Class B (QP)')
ax.plot(fgp / 1e6, cispr_A(fgp), 'k--', lw=1.3, label='CISPR 32 Class A (QP)')
for name, iload, ipk, duty, col in cases:
    r = results[name]
    lab = '%s  (Iload=%.2f A)  margin_min=%+.1f dB' % (name, iload, r['marg'].min())
    ax.semilogx(r['f'] / 1e6, r['dbuv'], '-o', ms=3.5, color=col, label=lab)
    bad = r['marg'] < 0
    if bad.any():
        ax.semilogx(r['f'][bad] / 1e6, r['dbuv'][bad], 'x', ms=9, mew=2.2, color=col)
# data cross-check (dashed, same colour, no markers)
for name, col in [('full', '#d62728'), ('half', '#2ca02c'), ('noload', '#1f77b4')]:
    fd, ad = data_chk[name]
    r = results[name]
    Zt = r['Zt']
    dbv = 20 * np.log10(np.maximum(ad * Zt, 1e-30) / 1e-6)
    ax.semilogx(fd / 1e6, dbv, ':', lw=1.0, color=col, alpha=0.55)
ax.plot([], [], 'k:', lw=1.0, label='data-derived cross-check (dotted, indicative)')
ax.set_xlim(0.15, 30)
ax.set_ylim(-20, 95)
ax.set_xlabel('Frequency [MHz]   (log)')
ax.set_ylabel('DM noise on 50 $\\Omega$ LISN   [dB$\\mu$V]')
ax.set_title('TPS62933 conducted emission (DM, frequency-domain model) vs CISPR 32\n'
             'harmonics of fsw = 805 kHz ; input filter 3x10uF | 1uH | 20uF+100nF+damping')
ax.grid(True, which='both', alpha=0.3)
ax.legend(loc='upper right', fontsize=8.2)
fig.tight_layout()
fig.savefig(OUT + 'fig_emi_cispr.png', dpi=140)
log('')
log('Wrote ' + OUT + 'fig_emi_cispr.png')

# ----------------------------------------------------------------------------
# FIG 2 : soft-start (analytic)
# ----------------------------------------------------------------------------
Tss = 14e-3
Vref, Css = 0.8, 100e-9
Iss_impl = Css * Vref / Tss
t = np.linspace(0, 20e-3, 4000)
v = np.where(t < Tss, VOUT * t / Tss, VOUT)
fig2, ax2 = plt.subplots(figsize=(9, 5.5))
ax2.plot(t * 1e3, v, 'r-', lw=2.6, label='Vout analytical ramp  0 -> 12 V,  Tss = 14 ms')
ax2.axvline(Tss * 1e3, color='k', ls='--', lw=1.2)
ax2.annotate('Tss = 14 ms\n(user-measured, authoritative)',
             xy=(14, 12), xytext=(15.0, 6.5),
             arrowprops=dict(arrowstyle='->'), fontsize=10)
ax2.axhline(VOUT, color='gray', ls=':', lw=1)
ax2.text(0.35, 12.35, '12 V', color='gray', fontsize=9)
ax2.text(0.6, 4.0,
         'ANALYTIC, NOT A SIMULATION\n'
         r'$T_{ss}\approx C_{ss}\,V_{REF}/I_{ss}$'+'\n'
         'Css=100 nF, Vref=0.8 V, Tss=14 ms\n'
         '=> Iss_implied = %.2f uA\n'
         '(datasheet Iss unavailable;\n 14 ms used as authoritative)' % (Iss_impl * 1e6),
         fontsize=9, bbox=dict(boxstyle='round', fc='#fff7cc', ec='gray'))
ax2.set_xlabel('Time [ms]')
ax2.set_ylabel('Vout [V]')
ax2.set_ylim(0, 13.5)
ax2.set_xlim(0, 20)
ax2.set_title('TPS62933F soft-start, analytic ramp (Css = 100 nF)')
ax2.grid(alpha=0.3)
ax2.legend(loc='lower right')
fig2.tight_layout()
fig2.savefig(OUT + 'fig_startup.png', dpi=140)
log('Wrote ' + OUT + 'fig_startup.png')
log('Iss implied by 14ms (Css=100nF, Vref=0.8V): Iss = %.2f uA' % (Iss_impl * 1e6))

with open(OUT + 'emi_numbers.txt', 'w') as f:
    f.write('\n'.join(lines) + '\n')
log('')
log('wrote emi_numbers.txt')
